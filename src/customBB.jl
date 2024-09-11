mutable struct TreeNode
    id::Int64
    BMP::JuMP.Model
	BSPs::Vector{Vector{Vector{JuMP.Model}}} #TODO In theory all subpaths are valid for all branches. 2 options: Keep one set of BSPs with all columns added (one large one), or keep copies (memory intensive and we may need to generate the same column multiple times)
    Subpaths::Pool
    integer::Bool
    xSol,
	zSol,
    # branches::Dict{Tuple{Int, Int, Int,}} # Branches are already enforced in the BMP
    parent::Int64
end
TreeNode(i::Int64, m::JuMP.Model, p::Pool, s::Solution, v::Vector{Vector{Bool}}, n::Int64) = TreeNode(i, m, p, false, s, v, n)
TreeNode() = TreeNode(Jump.Model(),Pool(),false, Solution(),Vector{Vector{Bool}}(),0)

mutable struct Summary
    LB::Float64         # lower bound (updated when dequeueing)
    UB::Float64         # upper bound (updated at each node solved)
    time::Float64       # elapsed time (only updated after solving children or when time limit is reached or optimal converged)
    nSolved::Int64      # number of B&B tree nodes solved to LP optimality
    nQueue::Int64       # length of the queue (unexplored nodes)
    # bestSol::Solution   # best solution so far
    # CGIts::Int64        # total column generation iterations
    # cols::Int64         # total new generated columns (not counting initial dummy ones)
end
Summary() = Summary(zero(Float64), typemax(Float64), zero(Float64), zero(Int64), zero(Int64)) #, Solution(), zero(Int64), zero(Int64))

function findBestBranch(inst::Instance, sol::Solution)
    @unpack N, P, Pi = inst
    maxC = maximum(length.(Pi))
    solTimes = [Vector{Int64}() for n in 1:N, p in 1:maxC]
    solBerths = [Vector{Int64}() for n in 1:N, p in 1:maxC]
    for sch in sol.schedules
        i = sch.ship
        for (c, vis) in enumerate(sch.visits)
            # println("$i, $vis")
            p,b,t = vis
            push!(solTimes[i,c], t)
            push!(solBerths[i,c], b)
        end
    end
    devTimes = zeros(Float64, N, maxC)
    avgTimes = zeros(Float64, N, maxC)
    devBerths = zeros(Float64, N, maxC)
    avgBerths = zeros(Float64, N, maxC)
    # avoid NaN values (not compatible with argmax())
    for n in 1:N, (c, p) in enumerate(Pi[n])
        if length(solTimes[n,c]) > EPSVAL
            avgTimes[n,c] = mean(solTimes[n,c])
            avgBerths[n,c] = mean(solBerths[n,c])
            if length(solTimes[n,c]) > 1 + EPSVAL
                devTimes[n,c] = std(solTimes[n,c])
                devBerths[n,c] = std(solBerths[n,c])
            end
        end
    end
    maxDevT = maximum(devTimes)
    maxDevB = maximum(devBerths)
    branchOnTime = maxDevT > maxDevB ? true : false
    if branchOnTime
        idx = argmax(devTimes)
        val = avgTimes[idx]
        ship = idx[1]
        visit = idx[2]
        return branchOnTime, ship, visit, Int(round(val,digits=0))
    else
        idx = argmax(devBerths)
        val = avgBerths[idx]
        ship = idx[1]
        visit = idx[2]
        return branchOnTime, ship, visit, Int(round(val,digits=0))
    end
end

function priceNode!(inst::Instance,
	# sol::Solution,
	# rmp::JuMP.Model,
	# pool::Pool,
	tn::TreeNode,
	summary::Summary,
	queue::PriorityQueue{Int64, Float64, Base.Order.ForwardOrdering},
	# node::Int64, #can the dequeue be done inside this function instead?
	# bestSol::Solution,
	gs::Vector{MyGraph},
	totalNodes::Int64,
	orders::Vector{Vector{Int64}},
	Params::Vector{Tuple{Int64, Int64, Int64, Int64}},
	Duals::Vector{Vector{Vector{Int64}}},
	AffCons::Vector{Vector{Vector{Int64}}},
	Node::Array{Int64, 4})
	@unpack N = inst
	@unpack sol, rmp, pool, visit = tn
	# @unpack UB, bestSol, nSolved, nQueue = summary
	doneCG = false
	it = 0
	cols = 0 # new generated columns
	objVal = zero(Float64)
	solVal = Vector{Float64}()
	while !doneCG #&& it <50
		it += 1
		optimize!(rmp)
		objVal = objective_value(rmp)
		println(" CG IT ", it, "\t obj = ", (round(objVal, digits = 1)), "\t solve time = ", (round(solve_time(rmp), digits = 2)), "\t Cols = ", length(rmp[:lambda]))
		solVal = value.(rmp[:lambda])
		dPart = dual.(rmp[:part])
		dPack = dual.(rmp[:pack])
		nodeBonus = [zeros(Float64, totalNodes) for i in 1:N]
		doneCG = true
		for i in 1:N
			for n ∈ 2:totalNodes-1
				if dPack[n] < -EPSVAL
					nodeBonus[i][collect(Duals[i][n])] .+= dPack[n]
				end
			end
			path, cost = shortestPath(gs[i], orders[i], visit[i], nodeBonus[i], 1, totalNodes, totalNodes)
			rc = cost - dPart[i]
			# println("ship $i, path $path, rc $rc")
			doneSub = (rc > -EPSVAL)
			if !doneSub
				doneCG = false
				col = addColumn(i, path, AffCons[i], rmp, gs[i], Params)
				push!(pool.cols, col)
				cols += 1
				# println("Ship $i, visits: ", col.visits)
			end
		end
	end
	### solve MIP ###
	# set_binary.(rmp[:lambda])
	# optimize!(rmp)
	# objVal = objective_value(rmp)
	# solVal = value.(rmp[:lambda])
	#################
	clearSol(sol)
	sol.cost = objVal
	sol.solVal = solVal
	tn.int = true
	for (i, v) in enumerate(solVal)
		if v > EPSVAL && i > N
			push!(sol.schedules, pool.cols[i-N])
			tn.int = v < 1 - EPSVAL ? false : true
		end
	end
	# tn.int ? println(" Solution is INTEGER") : println(" Solution is FRACTIONAL")
	# println(" Sol cols: ", length(sol.schedules))
	# println(" OBJ = ",objVal)
	if tn.int
		if objVal < summary.UB
			# best solution so found
			summary.UB = objVal
			summary.bestSol = deepcopy(sol)
			# drawSolution(sol,inst)
		end
	end
	if !tn.int && objVal < summary.UB
		queue[summary.nSolved+1] = objVal # node to be explored
	end
	summary.nSolved += 1
	summary.nQueue = length(queue)
	summary.CGIts += it
	summary.cols += cols
	return summary
end

function priceChildren!(inst::Instance,
	branchOnTime::Bool,
	ship::Int64,
	visit::Int64,
	val::Int64,
	solTree::Vector{TreeNode},
	node::Int64,
	summary::Summary,
	queue::PriorityQueue{Int64, Float64, Base.Order.ForwardOrdering},
	# nodeCount::Int64, #TODO: can the dequeue doen inside this function instead?
	# bestSol::Solution,
	gs::Vector{MyGraph},
	totalNodes::Int64,
	orders::Vector{Vector{Int64}},
	Node::Array{Int64, 4},
	Params::Vector{Tuple{Int64, Int64, Int64, Int64}},
	Duals, #::Vector{Vector{Set{Int64}}},
	AffCons, #::Vector{Vector{Set{Int64}}}
	startTime,
	tl,
	bigIP,
	)
	tn = solTree[node]
	visits, rmps = updateVisitableNodesAndPools(inst, branchOnTime, ship, visit, val, tn, Node)
	for idx in 1:2
		# nodeCount += 1
		child = TreeNode(summary.nSolved + 1, rmps[idx], deepcopy(tn.pool), Solution(), visits[idx], node)
		# println("Solving child ", idx - 1)
		# summary = priceNode!(inst, child,
		#             summary, queue,
		#             gs, totalNodes, orders,
		#             Params, Duals, AffCons)
		summary = priceNodeLowMem!(inst, child,
			summary, queue,
			gs, totalNodes, orders,
			Params, Duals, AffCons, Node,
			startTime, tl,
			bigIP,
		)
		push!(solTree, child)
	end
	# we no longer need the parent node #TODO: checf if true
	clearNode!(solTree[node])
	return summary
end

"""
runs the algorithm on a given instance

### Keywords
(see function arguments)
### Returns
* lower bound,
* upper bound,
* number of iterations,
* time spent solving MP,
* time spent solving SP,
* number of Benders cuts generated,
* algorithm run time,
* number of sub-paths,
* sol object
"""
function runAlg(
    R::RouteSettingData,
    Inst::InstanceSettingData,
    Subpaths::Vector{Vector{Vector{Pool}}},                # pool of sub-paths (precomputed)
    all_subpath_graphs::Array{Vector{Vector{Vector{TEN}}},2},      # pool of time-expanded networks (used to generate sub-paths)
    all_load_expanded_graphs::Vector{Vector{TSLGraph}},              # pool of time-space-load networks (used to add sub-paths in the sub-problem)
    CG::Bool,      # solve second-stage with column generation
    # PO::Bool,      # add pareto-optimal cuts
    Heur::Bool,    # heuristic CG or not
    normalized::Bool,
    rootnode::Bool, # relax MP and solve rootnode with Benders
    num_focus::Int64=NUM_FOCUS,
    timelimit::Int=MP_TIME_LIMIT, # time limit for the algorithm execution
    )
    @unpack numL, numS, TDmp, TDsp, K, M, δ = Inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times = R
    # build first stage model
    MP = firstStageModel(R, Inst, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus)
    # we want to solve the root node first (set to false if not)
    local undo_relax
    if rootnode
        undo_relax = relax_integrality(MP)
    end
    # build second stage models (with all sub-paths in SPs)
    subproblems = [[[buildGenericSecondStage(R, Subpaths[l][t][s], Inst, l,t,s, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus) for s in 1:numS] for t in eachindex(R.Lines[l].freq)] for l in 1:numL]

    done = false
    LB = typemin(Float64)       # lower bound
    UB = typemax(Float64)       # upper bound
    # LBs = Vector{Float64}()     # lower bounds at each iteration
    # UBs = Vector{Float64}()     # upper bounds at each iteration
    cuts = ConstraintRef[]      # pool of Benders cuts
    it = 0                      # iteration counter
    RNtime = 0.               # time spent solving the root node
    solveTMP = 0.0              # time spent in MP
    # MPtimes = Vector{Float64}()
    solveTSP = 0.0              # time spent in SP
    FScosts = 0.0               # first-stage costs (0)
    SScosts = 0.0               # second stage costs
    sol = Sol(numS)   # initialize solution
    num_second_stage_vars = 0
    # numCols = Vector{Int}([sum(length(SPs[l,t][s].all) for l in 1:numL, t in 1:T, s in 1:numS)])
    # numCuts = Vector{Int}([length(cuts)])

    # settings for the PO cuts
    # # q0, x0, z0, t0 = core_point(R)  # core point of MP
    # local x0
    # local z0
    # ϵ = 0.
    # λ = 0.5
    # if PO
    #     ϵ = EPS                   # epsilon tolerance (for PO cuts). It needs to be ϵ <= MIP_GAP (to assume equivalent solutions of SP)
    # end
    local xVal
    local zVal
    MP_fractional = true
    startT = time()             # algorithm time tracker
    while !done                 # run until convergence
        it += 1
        optimize!(MP)                       # solve MP
        # check that MP is solved to optimality to retrieve valid Benders LB
        term_status = JuMP.termination_status(MP)
        objMP = term_status == MOI.OPTIMAL ? objective_value(MP) : objective_bound(MP)
        if term_status != MOI.OPTIMAL
            println(" Master problem not solved to optimality. STATUS: ", term_status)
        end
        FScosts = value(MP[:first_stage_costs])
        solvetime = solve_time(MP)          # solve time
        solveTMP += solvetime
        # push!(MPtimes, solvetime)
        # println("OBJ MP: ", objMP, " in ", round(solvetime, digits=2), " seconds.")
        xVal = value.(MP[:x])               # first stage X vars
        zVal = value.(MP[:z])               # first stage Z vars
        thetaVal = value.(MP[:theta])       # first stage Theta vars

        # core point update for PO cuts
        # if it == 1
        #     x0 = deepcopy(xVal)
        #     z0 = deepcopy(zVal)
        # else
        #     x0 = (λ * xVal) .+ ((1-λ) * x0)
        #     for l in 1:numL, t in 1:T, i in R.L[l].O, h in 0:R.H, s in 1:numS
        #         z0[l,t,i,h,s] = λ*zVal[l,t,i,h,s] + (1-λ)*z0[l,t,i,h,s]
        #     end
        # end

        # push!(LBs, objMP)
        LB = objMP > LB + 0.01 ? objMP : LB # update lower bound if better
        objSPS = [Vector{Float64}() for s in 1:numS]
        doneSPs = true
        for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS
            # println("$l, $t, $s")
            startCG = time()
            # update RHS with first stage solution
            set_normalized_rhs(subproblems[l][t][s][:c1], (xVal[l, t])) # + ϵ*x0[l,t]))
            set_normalized_rhs(subproblems[l][t][s][:c2], (xVal[l, t])) # + ϵ*x0[l,t]))
            for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
                if (l_,t_) == (l,t)
                    set_normalized_rhs(subproblems[l][t][s][:c4][p,(l,t)], (zVal[s,p,(l,t)])) # + ϵ*z0[l,t,k,h,s]))
                end
            end

            if CG
                # start CG procedure
                doneCG = false
                cgIt = 0        # column generation iterations
                t0 = Lines[l].ref_stop_time[Lines[l].ref_stop_id[1]][t]
                while !doneCG
                    cgIt += 1
                    # solve RMP (RMP is the restricted version of SP)
                    optimize!(subproblems[l][t][s])
                    solvetime = solve_time(subproblems[l][t][s])
                    # solveTSP += solvetime
                    if termination_status(subproblems[l][t][s]) != MOI.OPTIMAL
                        println("Term status of SP $l $t $s: ", termination_status(subproblems[l][t][s]))
                    end
                    objSP = objective_value(subproblems[l][t][s])
                    # println(" CG IT $cgIt, OBJ= ", objSP)
                    # obtain dual solution
                    # d1 = round(dual(subproblems[l][t][s][:c1]), digits=4)
                    # d2 = round(dual(subproblems[l][t][s][:c2]), digits=4)
                    # d3 = round.(dual.(subproblems[l][t][s][:c3]), digits=4)
                    d1 = dual(subproblems[l][t][s][:c1])
                    d2 = dual(subproblems[l][t][s][:c2])
                    d3 = dual.(subproblems[l][t][s][:c3])
                    d4 = dual.(subproblems[l][t][s][:c4])
                    # for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
                    #     if (l_,t_) == (l,t)
                    #         d4[p,(l,t)] = round(d4[p,(l,t)], digits=4)
                    #     end
                    # end

                    doneCG = true # we are done unless we find a column with negative reduced cost
                    for i in 1:length(Lines[l].ref_stop_id)-1, j in i+1:min(i + K, length(Lines[l].ref_stop_id))
                        g = all_subpath_graphs[l,s][t][i][j-i]
                        src = Lines[l].ref_stop_id[i]
                        snk = Lines[l].ref_stop_id[j]
                        labels, prev = Heur ? labelSettingAlgInCG(l,t, g, 1, g.order, Lines[l].capacity, d4) : fullLabelSettingAlgInCG(l,t,g, 1, g.order, Lines[l].capacity, d4)
                        # recreate paths
                        paths, labs = getPathsAndLabelsInCG(g.V, 1, labels, prev)
                        t1 = ceil(Int, Dates.value(Dates.Second(R.Lines[l].ref_stop_time[src][t]-t0))/TDsp)
                        t2 = ceil(Int, Dates.value(Dates.Second(R.Lines[l].ref_stop_time[snk][t]-t0))/TDsp)
                        for (idx, label) in enumerate(labs)
                            WT = sum(label[2]) # total walk + wait time
                            if WT < LARGE_NUMBER # otherwise it's not a valid path (NOTE: this check should not be necessary)
                                # compute number of passengers picked up
                                Q = 0
                                # total_weighted_delay = 0.
                                for (orig, time, num) in label[1]
                                    Q += num
                                    # arrival_dev = Dates.value(Dates.Second(R.Lines[l].arrival_time[t] - Demand[s][(orig,time)].request_dropoff_time))
                                    # if normalized
                                    #     arrival_dev /= Taxi_travel_times[orig]
                                    # end
                                    # if arrival_dev > -0.01
                                    #     total_weighted_delay += num*δ*arrival_dev
                                    # else
                                    #     total_weighted_delay += num*(δ/2)*abs(arrival_dev)
                                    # end
                                end
                                cost = WT #+ total_weighted_delay - M*Q # the cost of the sub-path
                                for c in 0:Int(Lines[l].capacity)-Q # add the sub-path for all combinations of "occupancy levels"
                                    count = length(Subpaths[l][t][s].all) # id of the sub-path (also helps keep count)
                                    node1 = all_load_expanded_graphs[l][t].p2n[i, t1+1, c+1]   # find tail node in SP
                                    node2 = all_load_expanded_graphs[l][t].p2n[j, t2+1, c+Q+1] # find head node in SP
                                    # compute reduced cost
                                    rc = label[3] - d3[node1] + d3[node2] # reduced cost fom lab setting alg + dual source and sink
                                    if rc < -0.001
                                        count += 1
                                        doneCG = false # if negative reduced cost, we have not found the LP optimal solution yet
                                        subP = SubPath(count, cost, Q, label[1], paths[idx], node1, node2)
                                        # add sub-path object to pool
                                        push!(Subpaths[l][t][s].all, subP)
                                        push!(Subpaths[l][t][s].in[node2], count)
                                        push!(Subpaths[l][t][s].out[node1], count)
                                        for (k, h, num) in label[1]
                                            push!(Subpaths[l][t][s].pax[k, h], count)
                                        end
                                        # add column to subproblem
                                        addSubPath!(R, l,t,s, subproblems[l][t][s], subP, node1, node2)
                                    end

                                end
                            end
                        end
                    end
                    # if we are done add objective value and generate Benders cut
                    if doneCG
                        # println("DONE CG")
                        push!(objSPS[s], objSP)
                        # add the cut only if the objective is better than the current first-stage theta value
                        if !Heur || (Heur && objSP > thetaVal[l,t,s] + EPS) # && xVal[l,t] > 0.01) || xVal[l,t] < 0.01
                            doneSPs = false # not done with Benders
                            # if optimal generate full benders cut with all the duals and add to RMP
                            # cRef = @constraint(MP, (d1 + d2) * MP[:x][l, t] + sum(d4[p,(l,t)] * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
                            cRef = @constraint(MP, round(d1 + d2, digits=6) * MP[:x][l, t] + sum(round(d4[p,(l,t)], digits=6) * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
                            push!(cuts, cRef)
                            # update RHS after adding first cut
                            if it == 1
                                delete(MP, MP[:recourse][l,t,s])
                            end
                        end
                    elseif time() - startT > timelimit
                        doneSPs = true
                    end
                end
            else
                # solve SP
                optimize!(subproblems[l][t][s])
                # solvetime = solve_time(SP[l][t][s])
                # solveTSP += solvetime
                objSP = objective_value(subproblems[l][t][s])
                push!(objSPS[s], objSP)
                if true #objSP > thetaVal[l, t, s] + 0.001
                    doneSPs = false
                    # d1 = round(dual(subproblems[l][t][s][:c1]), digits=4)
                    # d2 = round(dual(subproblems[l][t][s][:c2]), digits=4)
                    d1 = dual(subproblems[l][t][s][:c1])
                    d2 = dual(subproblems[l][t][s][:c2])
                    d4 = dual.(subproblems[l][t][s][:c4])
                    # for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
                    #     if (l_,t_) == (l,t)
                    #         d4[p,(l,t)] = round(d4[p,(l,t)], digits=4)
                    #     end
                    # end
                    # if optimal generate full benders cut with all the duals and add to RMP
                    # cRef = @constraint(MP, (d1 + d2) * MP[:x][l, t] + sum(d4[p,(l,t)] * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
                    cRef = @constraint(MP, round(d1 + d2, digits=6) * MP[:x][l, t] + sum(round(d4[p,(l,t)], digits=6) * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
                    push!(cuts, cRef)
                    # update RHS after adding first cut
                    if it == 1
                        delete(MP, MP[:recourse][l,t,s])
                    end
                end
            end
            solveTSP += time() - startCG
        end
        # push!(numCuts, length(cuts))
        # push!(numCols, sum(length(SPs[l,t][s].all) for l in 1:numL, t in 1:T, s in 1:numS) )

        # compute second-stage costs
        SScosts = sum(R.Pi[s] * sum(objSPS[s]) for s in 1:numS)
        ub = SScosts + FScosts
        # push!(UBs, ub)
        #   valid UB = first-stage costs + sum of all subproblem objs
        UB = ub < UB - 0.01 ? ub : UB # update upper bound if better
        println("IT $it, LB = $LB, UB = $UB, numCuts = ", length(cuts))
        # STOP if no more cuts can be added or if we have converged
        if rootnode && (doneSPs || abs((UB - LB) / LB) < MIP_GAP || time() - startT > timelimit)
            rootnode = false
            UB = typemax(Float64)
            RNtime = time() - startT
            println("ROOT NODE terminated in $(round(RNtime, digits=2))")
            undo_relax()
            if time() - startT > timelimit
                done = true
                MP_fractional = true
            end
        elseif doneSPs || abs((UB - LB) / LB) < MIP_GAP || time() - startT > timelimit
            done = true
            MP_fractional = false
        end
        if done
            ipTime = time()
            if MP_fractional
                optimize!(MP)                       # solve MP
                # check that MP is solved to optimality to retrieve valid Benders LB
                term_status = JuMP.termination_status(MP)
                objMP = term_status == MOI.OPTIMAL ? objective_value(MP) : objective_bound(MP)
                if term_status != MOI.OPTIMAL
                    println(" Master problem not solved to optimality. STATUS: ", term_status)
                end
                FScosts = value(MP[:first_stage_costs])
                solvetime = solve_time(MP)          # solve time
                solveTMP += solvetime
                # push!(MPtimes, solvetime)
                println("OBJ MP: ", objMP, " in ", round(solvetime, digits=2), " seconds.")
                xVal = value.(MP[:x])               # first stage X vars
                zVal = value.(MP[:z])               # first stage Z vars
                # thetaVal = value.(MP[:theta])       # first stage Theta vars
                LB = objMP > LB + 0.01 ? objMP : LB
            end
            # solve IP version of the SPs to check how far is the integer (second-stage) solution
            println("FS costs: ", FScosts)#, " SS costs: ", SScosts)
            sol.cost = FScosts
            # for l in 1:numL, t in 1:T, s in 1:numS, k in R.L[l].O, h in 0:R.H
            #     if haskey(R.D[s], (k,h))
            #         if zVal[l,t,k,h,s] > 0.01
            #         #     if k in [42429324, 42445917, 42448317, 42445910, 42439191]
            #         #         if h in [10,11,14,15]
            #         #             println("ZVal $l, $t, $k, $h > 0")
            #         #         end
            #         #     end
            #             sol.PaxAssigned += R.Pi[s]*R.D[s][k,h].num*zVal[l, t, k, h, s]
            #         end
            #     end
            # end
            for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS
                # set back RHS just in case
                set_normalized_rhs(subproblems[l][t][s][:c1], xVal[l, t])
                set_normalized_rhs(subproblems[l][t][s][:c2], xVal[l, t])
                num_second_stage_vars += length(Subpaths[l][t][s].all)
                for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
                    if (l_,t_) == (l,t)
                        set_normalized_rhs(subproblems[l][t][s][:c4][p,(l,t)], zVal[s,p,(l,t)])
                    end
                end
                set_binary.(subproblems[l][t][s][:y]) # second-stage variables as binary
                optimize!(subproblems[l][t][s])
                if has_values(subproblems[l][t][s])
                    # addSubproblemSol!(m, subproblems[l][t][s], sol, xVal[l,t], l,t,s, R, Inst, Subpaths[l][t][s], subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs)
                    objSP = objective_value(subproblems[l][t][s])
                    sol.cost += R.Pi[s]*objSP   # objective is second-stage costs
                else
                    objSP = typemax(Float64)
                    sol.cost += R.Pi[s]*objSP
                end
                # for (id, a) in enumerate(Subpaths[l][t][s].all)
                #     if value(subproblem[l][t][s][:y][a.id]) > 0.999
                #         @unpack c, Q, W = Subpaths[l][t][s].all[id]
                #         # sol.PaxServed += Q
                #         # sol.TotalTime += W
                #         # push!(sol.paths[l,t,s], a)
                #         # for (i,h) in a.O
                #         #     println(i,", ",h, " by ", l, " t ",t," in s ", s)
                #         # end
                #     end
                # end
                # numSPs += length(SPs[l,t][s].all)
            end
            # for l in 1:numL, t in eachindex(R.Lines[l].freq)
            #     sol.xVal[l,t] = xVal[l,t]
            # end
            # println("SOLVE IP: ", time() - ipTime)
        end
    end
    # plotAlgRun(LBs, UBs, MPtimes, numCuts, numCols)
    totT = time() - startT # stop algorithm time
    println("ALG TIME: ", totT)
    # println(" MP time: ", solveTMP)
    # println( "SP+CG time: ", solveTSP)
    println("IT $it, LB = $LB, UB = $UB, IP(UB) = $(sol.cost), numCuts = ", length(cuts), " TOTAL T: ", totT)
    return LB, sol.cost, totT, it, length(cuts), num_second_stage_vars, RNtime, solveTMP, solveTSP, xVal
end

function findBranch(Inst::InstanceSettingData, xVal, zVal)
	xFrac = abs.(xVal .- 0.5)
	zFrac = abs.(zVal .- 0.5)
	xIdx = argmin(xVal)
	zIdx = argmin(zVal)

	return xIdx, zIdx
end

function clearNode!(tn::TreeNode)
    tn.BMP = Model()
	tn.BSPs = Vector{Vector{Vector{JuMP.Model}}}()
    tn.Subpaths = Pool()
    tn.integer = false
    tn.xSol = []
	tn.xSol = []
    #? Do we need to keep the parent node?
end

function solveBSPs!(
	R::RouteSettingData,
    Inst::InstanceSettingData,
    Subpaths::Vector{Vector{Vector{Pool}}},                # pool of sub-paths (precomputed)
    all_subpath_graphs::Array{Vector{Vector{Vector{TEN}}},2},      # pool of time-expanded networks (used to generate sub-paths)
    all_load_expanded_graphs::Vector{Vector{TSLGraph}},              # pool of time-space-load networks (used to add sub-paths in the sub-problem)
	BMP::JuMP.Model,
	subproblems::Vector{Vector{Vector{JuMP.Model}}},
    CG::Bool,      # solve second-stage with column generation
    Heur::Bool,    # heuristic CG or not
    normalized::Bool,
	timelimit::Int=MP_TIME_LIMIT, # time limit for the algorithm execution
	)
	@unpack numL, numS, TDmp, TDsp, K, M, δ = Inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times = R

	xVal = value.(BMP[:x])               # first stage X vars
	zVal = value.(BMP[:z])               # first stage Z vars
	thetaVal = value.(BMP[:theta])       # first stage Theta vars
	objSPS = [Vector{Float64}() for s in 1:numS]
	doneSPs = true
	doneCG = true
	for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS
		# println("$l, $t, $s")
		startCG = time()
		# update RHS with first stage solution
		set_normalized_rhs(subproblems[l][t][s][:c1], (xVal[l, t])) # + ϵ*x0[l,t]))
		set_normalized_rhs(subproblems[l][t][s][:c2], (xVal[l, t])) # + ϵ*x0[l,t]))
		for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
			if (l_,t_) == (l,t)
				set_normalized_rhs(subproblems[l][t][s][:c4][p,(l,t)], (zVal[s,p,(l,t)])) # + ϵ*z0[l,t,k,h,s]))
			end
		end

		if CG
			# start CG procedure
			doneCG = false
			cgIt = 0        # column generation iterations
			t0 = Lines[l].ref_stop_time[Lines[l].ref_stop_id[1]][t]
			while !doneCG
				cgIt += 1
				# solve RMP (RMP is the restricted version of SP)
				optimize!(subproblems[l][t][s])
				solvetime = solve_time(subproblems[l][t][s])
				# solveTSP += solvetime
				if termination_status(subproblems[l][t][s]) != MOI.OPTIMAL
					println("Term status of SP $l $t $s: ", termination_status(subproblems[l][t][s]))
				end
				objSP = objective_value(subproblems[l][t][s])
				# println(" CG IT $cgIt, OBJ= ", objSP)
				# obtain dual solution
				d1 = dual(subproblems[l][t][s][:c1])
				d2 = dual(subproblems[l][t][s][:c2])
				d3 = dual.(subproblems[l][t][s][:c3])
				d4 = dual.(subproblems[l][t][s][:c4])

				doneCG = true # we are done unless we find a column with negative reduced cost
				for i in 1:length(Lines[l].ref_stop_id)-1, j in i+1:min(i + K, length(Lines[l].ref_stop_id))
					g = all_subpath_graphs[l,s][t][i][j-i]
					src = Lines[l].ref_stop_id[i]
					snk = Lines[l].ref_stop_id[j]
					labels, prev = Heur ? labelSettingAlgInCG(l,t, g, 1, g.order, Lines[l].capacity, d4) : fullLabelSettingAlgInCG(l,t,g, 1, g.order, Lines[l].capacity, d4)
					# recreate paths
					paths, labs = getPathsAndLabelsInCG(g.V, 1, labels, prev)
					t1 = ceil(Int, Dates.value(Dates.Second(R.Lines[l].ref_stop_time[src][t]-t0))/TDsp)
					t2 = ceil(Int, Dates.value(Dates.Second(R.Lines[l].ref_stop_time[snk][t]-t0))/TDsp)
					for (idx, label) in enumerate(labs)
						WT = sum(label[2]) # total walk + wait time
						if WT < LARGE_NUMBER # otherwise it's not a valid path (NOTE: this check should not be necessary)
							# compute number of passengers picked up
							Q = 0
							# total_weighted_delay = 0.
							for (orig, time, num) in label[1]
								Q += num
							end
							cost = WT #+ total_weighted_delay - M*Q # the cost of the sub-path
							for c in 0:Int(Lines[l].capacity)-Q # add the sub-path for all combinations of "occupancy levels"
								count = length(Subpaths[l][t][s].all) # id of the sub-path (also helps keep count)
								node1 = all_load_expanded_graphs[l][t].p2n[i, t1+1, c+1]   # find tail node in SP
								node2 = all_load_expanded_graphs[l][t].p2n[j, t2+1, c+Q+1] # find head node in SP
								# compute reduced cost
								rc = label[3] - d3[node1] + d3[node2] # reduced cost fom lab setting alg + dual source and sink
								if rc < -0.001
									count += 1
									doneCG = false # if negative reduced cost, we have not found the LP optimal solution yet
									subP = SubPath(count, cost, Q, label[1], paths[idx], node1, node2)
									# add sub-path object to pool
									push!(Subpaths[l][t][s].all, subP)
									push!(Subpaths[l][t][s].in[node2], count)
									push!(Subpaths[l][t][s].out[node1], count)
									for (k, h, num) in label[1]
										push!(Subpaths[l][t][s].pax[k, h], count)
									end
									# add column to subproblem
									addSubPath!(R, l,t,s, subproblems[l][t][s], subP, node1, node2)
								end

							end
						end
					end
				end
				# if we are done add objective value and generate Benders cut
				if doneCG
					# println("DONE CG")
					push!(objBSPs[s], objSP)
					# add the cut only if the objective is better than the current first-stage theta value
					if !Heur || (Heur && objSP > thetaVal[l,t,s] + EPS) # && xVal[l,t] > 0.01) || xVal[l,t] < 0.01
						doneSPs = false # not done with Benders
						# if optimal generate full benders cut with all the duals and add to RMP
						# cRef = @constraint(MP, (d1 + d2) * MP[:x][l, t] + sum(d4[p,(l,t)] * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
						cRef = @constraint(BMP, round(d1 + d2, digits=6) * BMP[:x][l, t] + sum(round(d4[p,(l,t)], digits=6) * BMP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= BMP[:theta][l, t, s])
						push!(cuts, cRef)
						# update RHS after adding first cut
						if it == 1
							delete(BMP, BMP[:recourse][l,t,s])
						end
					end
				elseif time() - startT > timelimit
					doneSPs = false
					doneCG = false
				end
			end
		else
			# solve SP
			optimize!(subproblems[l][t][s])
			objSP = objective_value(subproblems[l][t][s])
			push!(objBSPs[s], objSP)
			if true #objSP > thetaVal[l, t, s] + 0.001
				doneSPs = false
				d1 = dual(subproblems[l][t][s][:c1])
				d2 = dual(subproblems[l][t][s][:c2])
				d4 = dual.(subproblems[l][t][s][:c4])
				# cRef = @constraint(MP, (d1 + d2) * MP[:x][l, t] + sum(d4[p,(l,t)] * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
				cRef = @constraint(BMP, round(d1 + d2, digits=6) * BMP[:x][l, t] + sum(round(d4[p,(l,t)], digits=6) * BMP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= BMP[:theta][l, t, s])
				push!(cuts, cRef)
				# update RHS after adding first cut
				if it == 1
					delete(BMP, BMP[:recourse][l,t,s])
				end
			end
		end
		solveTSP += time() - startCG
	end
	return doneSPs, doneCG, BMP, Inst.Pi .* objBSPs
end

function solveIntBSPs!(
	R::RouteSettingData,
    Inst::InstanceSettingData,
	BMP::JuMP.Model,
	subproblems::Vector{Vector{Vector{JuMP.Model}}},
	)
	@unpack numL, numS, TDmp, TDsp, K, M, δ = Inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times = R
	xVal = value.(BMP[:x])               # first stage X vars
	zVal = value.(BMP[:z])               # first stage Z vars
	objIntBSPs = 0.
	for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS
		# set back RHS just in case
		set_normalized_rhs(subproblems[l][t][s][:c1], xVal[l, t])
		set_normalized_rhs(subproblems[l][t][s][:c2], xVal[l, t])
		for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
			if (l_,t_) == (l,t)
				set_normalized_rhs(subproblems[l][t][s][:c4][p,(l,t)], zVal[s,p,(l,t)])
			end
		end
		set_binary.(subproblems[l][t][s][:y]) # second-stage variables as binary
		optimize!(subproblems[l][t][s])
		if has_values(subproblems[l][t][s])
			# addSubproblemSol!(m, subproblems[l][t][s], sol, xVal[l,t], l,t,s, R, Inst, Subpaths[l][t][s], subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs)
			objSP = objective_value(subproblems[l][t][s])
			objIntBSPs += R.Pi[s]*objSP   # objective is second-stage costs
		else
			objSP = typemax(Float64)
			objIntBSPs += R.Pi[s]*objSP
		end
	end
	return objIntBSPs
end

function checkIntegrality(m::JuMP.Model)
	xVal = value.(m[:x])
	zVal = value.(m[:z])
	for (i,v) in enumerate(xVal)
		if EPS < v < 1- EPS
			return false
		end
	end
	for (i,v) in enumerate(zVal)
		if EPS < v < 1- EPS
			return false
		end
	end
	return true
end

function main(filename::String, time_limit=MP_TIME_LIMIT)
	CG = true
	HEUR = true
	TRANSIT = true
	TRAIN = true
	RELAX_MP = true
	trial_info = first(eachrow(filter(row -> row[:trial_id] == trial_id, trials)))
    
    num_lines = Int(trial_info[:num_lines])
    demand_hor = Int(trial_info[:demand_hor])
    num_scen = Int(trial_info[:num_scen])
    maxDev = convert(Float64, "max_dev" in names(trial_info) ? trial_info[:max_dev] : MAX_SCHEDULE_DEV)
    maxWalk = convert(Float64, "max_walk" in names(trial_info) ? trial_info[:max_walk] : MAX_WALKING_METERS)
    capacity = convert(Int, "capacity" in names(trial_info) ? trial_info[:capacity] : CAPACITY)
    fleet_size = Int("fleet_size" in names(trial_info) ? trial_info[:fleet_size] : FLEET_SIZE)
    K = Int("K" in names(trial_info) ? trial_info[:K] : 0) + 1
    max_sch_dev = Int("max_sch_dev" in names(trial_info) ? trial_info[:max_sch_dev] : MAX_SCHEDULE_DEV)

	inst = InstanceSettingData(numL, demT, numS, TRAIN, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,max_sch_dev,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    R, m = buildRouteSettingData(inst);
    
    (all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs, enumeration_time), precomp_time_sp = @timed generateSubPathSet(m, R, inst, CG, !TRANSIT);
    # LB, UB, alg_time, it, num_cuts, num_vars, rootnode_time, MP_time, SP_time,_ = runAlg(R,inst, all_subpaths, all_subpath_graphs, all_load_expanded_graphs, CG, HEUR, NORMALIZED, RELAX_MP)


	@unpack numL, numS, TDmp, TDsp, K, M, δ = inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times = R
    # build first stage model
    BMP = firstStageRelaxedModel(R, inst, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus)
    # build second stage models (with all sub-paths in SPs)
    subproblems = [[[buildGenericSecondStage(R, Subpaths[l][t][s], inst, l,t,s, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus) for s in 1:numS] for t in eachindex(R.Lines[l].freq)] for l in 1:numL]


	# initialize tree
		# priority queue
		# tree
	summary = Summary() # we keep track of valid LB and UB (to the full IP problem)
	queue = PriorityQueue{Int64, Float64}()
	intermediateSols = PriorityQueue{Int64, Float64}()
	solTree = Vector{TreeNode}()

	TreeNode = TreeNode(1, BMP, subproblems, all_subpaths, false, [], [], 0)
	# push!(solTree, rootNode)
	# solve root node
	startTime = time()
	# # start with root node
	# node, lb = dequeue_pair!(queue)
	# TreeNode = solTree[node]
	BendersDone = false
	while summary.UB - summary.LB > EPS && !isempty(queue) && time() - startTime < time_limit
		if BendersDone
			# select node to solve
			node, lb = dequeue_pair!(queue)
			summary.LB = lb > summary.LB + EPS ? lb : summary.LB
			TreeNode = solTree[node]
		end
		optimize!(TreeNode.BMP)
		# done = false
		# while !done
		if has_values(TreeNode.BMP)
			objBMP = objective_value(TreeNode.BMP)
			if objBMP < summary.UB - EPS
				intBMP = checkIntegrality(TreeNode.BMP)
				if intBMP
					# solve all BSPs with CG
					objBSPs = Vector{Float64}()
					# BendersDone = true
					BendersDone, doneCG, TreeNode.BMP, objBSPs = solveBSPs!(R,inst,TreeNode.Subpaths,all_subpath_graphs,all_load_expanded_graphs,TreeNode.BMP,TreeNode.BSPs,CG, HEUR,NORMALIZED)

					if BendersDone
						validLB = objBMP + sum(objBSPs) # objBSPs are already multiplied by inst.Pi
						if validLB > summary.UB - EPS
							# fathom node
							clearNode!(solTree[node])
							BendersDone = true
							continue
						else
							# add BMPsol and LB to solPool: enough with not clearing the node
							# add node info to sol queue
							intermediateSols[node] = lb
							# find integral (optimal?) solution to BSPs and update UB
								#? Two options: if solved to integer optimality, we don't need the postprocess but we may need to solve a lot
								#? 				if not solved any, we may need to explore the entire B&B tree because we never update the UB
							objIntBSPs = solveIntBSPs!(R,inst,TreeNode.BMP,TreeNode.BSPs)
							summary.UB = TreeNode.BMP[:first_stage_costs] + objIntBSPs < summary.UB - EPS ? TreeNode.BMP[:first_stage_costs] + objIntBSPs : summary.UB 

							# choose new variable to branch on
							# FIND BRANCHING (most fractional value)
							#! remember first stage solution is integer
							# We should be fine fathoming the node and continuing with another node already in the queue
							
						end
					else
						#Go back to solve BMP
						continue
					end
				else
					# add node to queue only if it is not fathomed
					queue[node] = validLB
					# choose new variable to branch on
					node, lb = dequeue_pair!(queue) # we branch on next node in queue
					summary.LB = lb > summary.LB + EPS ? lb : summary.LB
					TreeNode = solTree[node]
					# FIND BRANCHING (most fractional value)
					xIdx, zIdx = findBranch(inst, JuMP.value(TreeNode.BMP[:x]), JuMP.value(TreeNode.BMP[:z]))
					if EPS < value(TreeNode.BMP[:x])[xIdx] < 1 - EPS
						# branch on xVals first
						#TODO Copy models and everything 
						#TODO Add branching constraints in BMP
						#TODO add children to tree
					elseif EPS < value(TreeNode.BMP[:x])[xIdx] < 1 - EPS
						# if no fractional xVal, move to zVal
						#TODO Copy models and everything 
						#TODO Add branching constraints in BMP
						#TODO add children to tree
					end
					#? How do we enforce that the next two nodes to be solved are these children?
					
					# BendersDone = true #TODO Check if needed
					continue
				end
			else
				# fathom node
				clearNode!(solTree[node])
				BendersDone = true #TODO Check if needed
				continue
			end
		else
			term_status = JuMP.termination_status(tn.BMP)
			if term_status != MOI.OPTIMAL
				println(" Benders master problem not solved to optimality. STATUS: ", term_status)
			end
			# fathom node if infeasible
			clearNode!(solTree[node])
			#? What to do if time limit?
			BendersDone = true #TODO Check if needed
			continue
		end
		

		



		summary = solveBMP!(
			inst, R, 
			solTree[1],
			summary, queue,
			BMP, subproblems,
			all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs,
			startTime,
		)
		summary.LB = solTree[1].sol.cost
		summary.time = time() - startTime
	# enter while loop (while not empty queue or not converged)
		# dequeue
		# find branch (on what)
		# solve children
		# store sol if needed (if BMP integer)
		# update queue and bounds

	# WHAT DO I NEED TO KEEP TRACK OF:
		# - xSol (if integer)
		# parent
		# BMP (with cuts so far)
		# branching (visitable nodes)








	inst = readInstFromFile(filename) #seed, Nin, Nout, quay discretization
	els = split(filename, '_')
	seed = parse(Int64, els[3])
	Nout = parse(Int64, els[5])
	# inst.N = X
	@unpack N, P, Bp, qli = inst
	# println("Instance: ", N, " ships, ", P, " ports and berth discret. of ", qli, " meters.")
	totalNodes, Node, Params = nodeIDmx(inst)
	Duals, AffCons = dualIDmxLowMem(inst, totalNodes, Node)
	# Duals, AffCons = dualIDmx(inst,totalNodes,Node)
	gs, orders, arcIdx = createGraphs(inst, totalNodes, Node, Params)
	# return gs
	# orders = [collect(1:gs[i].V) for i in 1:N]
	visitable = [ones(Bool, totalNodes) for i in 1:N]

	pool = Pool()
	# pool = generatePoolOfPaths(inst, Node, totalNodes, gs)
	sol = Solution()
	solTree = Vector{TreeNode}()

	# Amin = fill(100., N)
	# Amax = fill(100., N)
	# Bmin = fill(100., totalNodes)
	# Bmax = fill(100., totalNodes)
	# dMin = fill(10000., N)
	# dMax = fill(10000., N)
	# eMin = fill(1000., totalNodes)
	# eMax = fill(1000., totalNodes)

	rmp = Model(CPLEX.Optimizer)
	set_silent(rmp)
	@variable(rmp, lambda[1:N] >= 0)
	# @variable(rmp, 0 <= αMin[i in 1:N] <= Amin[i])
	# @variable(rmp, 0 <= αMax[i in 1:N] <= Amax[i])
	# @variable(rmp, 0 <= βMin[i in 2:totalNodes-1] <= Bmin[i])
	# @variable(rmp, 0 <= βMax[i in 2:totalNodes-1] <= Bmax[i])
	@objective(rmp, Min, sum(200000 * lambda[i] for i in 1:N)) #+  dMax[i]*αMax[i] - dMin[i]*αMin[i] for i in 1:N) + sum(eMax[i]*βMax[i] -  eMin[i]*βMin[i] for i in 2:totalNodes-1))
	@constraint(rmp, part[i in 1:N], lambda[i] >= 1) #+ αMax[i] - αMin[i] == 1)
	@constraint(rmp, pack[i in 2:totalNodes-1], 0 <= 1) #βMax[i] - βMin[i] <= 1)
	set_optimizer_attribute(rmp, "CPXPARAM_Threads", NumberThreads)
	# set_optimizer_attribute(rmp, "CPX_PARAM_LPMETHOD", 2)
	bigIP = BigGSPP()

	for col in pool.cols
		nq = ceil(Int, inst.shipsIn[col.ship].l / inst.qli)
		colCost = 0.0 #col.cost
		path = col.path
		lambdaNew = @variable(rmp, base_name = "lambda", lower_bound = 0)
		touchedConstraints = ConstraintRef[]
		vals = Float64[]
		push!(touchedConstraints, rmp[:part][col.ship])
		push!(vals, 1.0)
		for id in 1:length(path)-1
			n1 = path[id]
			n2 = path[id+1]
			# println(id)
			for (src, dst, cost) in gs[col.ship].outArcs[n1]
				if dst == n2
					colCost += cost
					# println("adding ", cost)
					break
				end
			end
			if 1 < n1
				(p, c, b, t) = Params[n1]
				hand = ceil(Int, inst.h[col.ship][id-1][b])
				t_ = t + hand - 1
				for c_ in 1:inst.xP[col.ship, p], b_ in b:b+nq-1, t2 in t:t_
					n_ = Node[p, c_, b_, t2+1]
					push!(touchedConstraints, rmp[:pack][n_])
					push!(vals, 1.0)
				end
			end
		end
		println(colCost)
		set_objective_coefficient(rmp, lambdaNew, colCost)
		set_normalized_coefficient.(touchedConstraints, lambdaNew, vals)
		push!(rmp[:lambda], lambdaNew)
	end
	cols = length(pool.cols)
	summary = Summary()
	# @unpack LB, UB, time, nSolved, nQueue = summary
	# UB = typemax(Float64)
	# LB = zero(Float64)
	# bestSol = Solution()
	queue = PriorityQueue{Int64, Float64}()
	# nodeCount = 1
	rootNode = TreeNode(1, rmp, pool, sol, visitable, 0)
	push!(solTree, rootNode)
	startTime = time()
	# SOLVE ROOT NODE
	summary = priceNodeLowMem!(inst, solTree[1],
		summary, queue,
		gs, totalNodes, orders,
		Params, Duals, AffCons, Node,
		startTime,
		tl,
		bigIP,
	)
	summary.LB = solTree[1].sol.cost
	summary.time = time() - startTime
	# println("LB: ", round(summary.LB), "\tUB: ", round(summary.UB), "\tNodes: ", summary.nSolved, "\tQueue: ", summary.nQueue, "\tElapsed: ", summary.time)
	# println("$seed\t$N\t$Nout\t$qli\t$(round(summary.LB))\t$(round(summary.UB))\t$(round(summary.time))")
	# return round(summary.LB), round(summary.UB), round(summary.time)
	while summary.UB - summary.LB > EPSVAL && !isempty(queue) && time() - startTime < 0.9 * tl
		node, lb = dequeue_pair!(queue)
		summary.LB = lb > summary.LB ? lb : summary.LB
		tn = solTree[node]
		# FIND BEST BRANCHING
		branchOnTime, ship, visit, val = findBestBranch(inst, tn.sol)
		# println("Branch on time = $branchOnTime, Ship $ship, Port $port, Val $val")
		# PRICE CHILDREN
		summary = priceChildren!(inst, branchOnTime, ship, visit, val,
			solTree, node, summary, queue,
			gs, totalNodes, orders,
			Node, Params, Duals, AffCons,
			startTime, tl,
			bigIP)
		summary.time = time() - startTime
		# println("LB: ", round(summary.LB), "\tUB: ", round(summary.UB), "\tNodes: ", summary.nSolved, "\tQueue: ", summary.nQueue, "\tElapsed: ", summary.time)
	end
	if summary.UB - summary.LB > EPSVAL && !isempty(queue)
		objUB, solUB = solveBigGSPP(inst, bigIP, 0.1 * tl, Params)
		if objUB < summary.UB
			summary.UB = objUB
			summary.bestSol = solUB
		end
		# summary.UB = objUB < summary.UB ? objUB : summary.UB
	else
		summary.LB = summary.UB
	end
	# println.(summary.bestSol.schedules)
	computeCostsSeparately!(inst, summary.bestSol)
	@unpack W, D, H, F = summary.bestSol
	# println("COSTS: $W, $D, $H, $F and sum", W +D+ H+ F)
	# println("DONE!")
	# for sch in summary.bestSol.schedules
	#     println(sch)
	# end
	summary.time = time() - startTime
	# drawSolution(summary.bestSol, inst)
	# pen = isInfeasible(inst, summary.bestSol)
	# println("PENALTY: ", pen)
	# checkFeasibilitySol(summary.bestSol,inst)
	println("$seed\t$N\t$Nout\t$qli\t$(round(summary.LB))\t$(round(summary.UB))\t$(round(summary.time))\t$(round(W))\t$(round(D))\t$(round(H))\t$(round(F))")
	return summary #.bestSol #summary,inst #, gs
end

