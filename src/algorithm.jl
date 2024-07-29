"""
direct method implementation

### Keywords
* `numL` - number of lines (pool size)
* `demT` - demand horizon (in seconds)
* `numS` - number of scenarios (in-sample)
* `K` - potential length of subpath (in number of consecutive reference stops)
* `CG` - column generation enabled
* `TRANSIT` - fixed-transit benchmark enabled
* `maxDev` - max deviation (meters)
* `maxWalk` - max pax walking (meters)
* `capacity` - vehicle capacity
### Returns
* `sol` - Sol object with solution and level of service summary
* `R` - RouteSettingData object
"""
function test_direct(numL::Int, demand_horizon::Int, numS::Int, K::Int, trainData::Bool, TRANSIT::Bool, maxDev::Float64, maxWalk::Float64, capacity::Int, fleet_size::Int)
    CG = false # in direct models, we only do full enumeration
    inst = InstanceSettingData(numL, demand_horizon, numS, trainData, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,MAX_SCHEDULE_DEV,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    # inst = InstanceSettingData(numL, demT, numS, K)
    R, m = buildRouteSettingData(inst);
    all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs = generateSubPathSet(m, R, inst, CG, TRANSIT);
    model = directTwoStageModel(R, all_subpaths, inst);
    optimize!(model)
    sol = buildSol(m, model, R, inst, all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs)
    println(sol.cost)
    return m, inst, sol, R, model
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
function runAlgAndGetSol(
    m::MapData,
    R::RouteSettingData,
    Inst::InstanceSettingData,
    Subpaths::Vector{Vector{Vector{Pool}}},                # pool of sub-paths (precomputed)
    all_subpath_graphs::Array{Vector{Vector{Vector{TEN}}},2},      # pool of time-expanded networks (used to generate sub-paths)
    all_load_expanded_graphs::Vector{Vector{TSLGraph}},     # pool of time-space-load networks (used to add sub-paths in the sub-problem)
    subpath_road_networks,
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
                    addSubproblemSol!(m, subproblems[l][t][s], sol, xVal[l,t], l,t,s, R, Inst, Subpaths[l][t][s], subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs)
                    # objSP = objective_value(subproblems[l][t][s])
                    # sol.cost += R.Pi[s]*objSP   # objective is second-stage costs
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
    return LB, sol.cost, totT, it, length(cuts), num_second_stage_vars, RNtime, solveTMP, solveTSP, xVal, sol
end

"""
test Benders-based algorithm
"""
function test_benders(numL::Int, demand_horizon::Int, numS::Int, K::Int, trainData::Bool, CG::Bool, HEUR::Bool, maxDev::Float64, maxWalk::Float64, capacity::Int, fleet_size::Int, TRANSIT::Bool, num_focus::Int64=NUM_FOCUS, time_limit::Int=MP_TIME_LIMIT)
    RELAX_MP = true
    inst = InstanceSettingData(numL, demand_horizon, numS, trainData, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,MAX_SCHEDULE_DEV,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    R, m = buildRouteSettingData(inst);
    # return R, inst
    all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs = generateSubPathSet(m, R, inst, CG, TRANSIT);
    LB, UB, alg_time, it, num_cuts, num_second_stage_vars, RNtime, solveTMP, solveTSP, xVal = runAlg(R,inst, all_subpaths, all_subpath_graphs, all_load_expanded_graphs, CG, HEUR, NORMALIZED, RELAX_MP, num_focus, time_limit)
    return R, xVal
end

"""
runs the algorithm with a given fixed first-stage solution

### Keywords
(see function arguments)
### Returns
* lower bound,
* upper bound,
* sol object
"""
function solve_second_stage(
    m::MapData,
    R::RouteSettingData,
    Inst::InstanceSettingData,
    xVal,
    Subpaths::Vector{Vector{Vector{Pool}}},                # pool of sub-paths (precomputed)
    subpath_road_networks,
    all_subpath_graphs::Array{Vector{Vector{Vector{TEN}}},2},      # pool of time-expanded networks (used to generate sub-paths)
    all_load_expanded_graphs::Vector{Vector{TSLGraph}},              # pool of time-space-load networks (used to add sub-paths in the sub-problem)
    CG::Bool,      # solve second-stage with column generation
    Heur::Bool,    # heuristic CG or not
    num_focus::Int64=NUM_FOCUS,
    timelimit::Int=MP_TIME_LIMIT, # time limit for the algorithm execution
    )
    @unpack numL, numS, TDmp, TDsp, K, M, δ = Inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times = R
    # build first stage model
    MP = firstStageModel(R, Inst, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus)
    # fix x values
    for l in 1:numL, t in eachindex(R.Lines[l].freq)
        fix(MP[:x][l,t], xVal[l,t], force=true)
    end
    # we want to solve the root node first (set to false if not)
    rootnode = false
    # undo_relax = relax_integrality(MP)
    # build second stage models (with all sub-paths in SPs)
    subproblems = [[[buildGenericSecondStage(R, Subpaths[l][t][s], Inst, l,t,s, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus) for s in 1:numS] for t in eachindex(R.Lines[l].freq)] for l in 1:numL]

    done = false
    LB = typemin(Float64)       # lower bound
    UB = typemax(Float64)       # upper bound
    cuts = ConstraintRef[]      # pool of Benders cuts
    it = 0                      # iteration counter
    FScosts = 0.0               # first-stage costs (0)
    SScosts = 0.0               # second stage costs
    sol = Sol(numS)   # initialize solution
    num_subpaths = 0                  # number of sub-path

    local zVal
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

        println("OBJ MP: ", objMP, " in ", round(solvetime, digits=2), " seconds.")
        zVal = value.(MP[:z])               # first stage Z vars
        thetaVal = value.(MP[:theta])       # first stage Theta vars

        # push!(LBs, objMP)
        LB = objMP > LB + 0.01 ? objMP : LB # update lower bound if better
        objSPS = [Vector{Float64}() for s in 1:numS]
        doneSPs = true
        for l in 1:numL, t in eachindex(R.Lines[l].freq)
            for s in 1:numS
                if xVal[l,t] > EPS
                    # startCG = time()
                    # update RHS with first stage solution
                    set_normalized_rhs(subproblems[l][t][s][:c1], (xVal[l,t])) # + ϵ*x0[l,t]))
                    set_normalized_rhs(subproblems[l][t][s][:c2], (xVal[l,t])) # + ϵ*x0[l,t]))
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
                            if it == 1
                                delete(MP, MP[:recourse][l,t,s])
                            end
                        end
                    end
                    #  solveTSP += time() - startCG
                else
                    cRef = @constraint(MP, 0 <= MP[:theta][l, t, s])
                    push!(cuts, cRef)
                end
            end
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
            println("ROOT NODE SOLVED in $(round(RNtime, digits=2))")
            undo_relax()
        elseif doneSPs || abs((UB - LB) / LB) < MIP_GAP || time() - startT > timelimit
            done = true
            ipTime = time()
            # solve IP version of the SPs to check how far is the integer (second-stage) solution
            println("FS costs: ", FScosts)#, " SS costs: ", SScosts)
            sol.cost = FScosts
            for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS
                # set back RHS just in case
                # set_normalized_rhs(subproblems[l][t][s][:c1], xVal[l, t])
                # set_normalized_rhs(subproblems[l][t][s][:c2], xVal[l, t])
                for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
                    if (l_,t_) == (l,t)
                        set_normalized_rhs(subproblems[l][t][s][:c4][p,(l,t)], zVal[s,p,(l,t)])
                    end
                end
                set_binary.(subproblems[l][t][s][:y]) # second-stage variables as binary
                optimize!(subproblems[l][t][s])
                if has_values(subproblems[l][t][s])
                    # objSP = objective_value(subproblems[l][t][s])
                    # sol.cost += R.Pi[s]*objSP   # objective is second-stage costs
                    addSubproblemSol!(m, subproblems[l][t][s], sol, xVal[l,t], l,t,s, R, Inst, Subpaths[l][t][s], subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs)
                # else
                #     objSP = typemax(Float64)
                #     sol.cost += R.Pi[s]*objSP
                end
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
    return LB, UB, sol
    # return SScosts, sol
end

"""
Function to run out-of-sample tests.

### Keywords
* `R` - route-network data
* `Inst` - instance data including only the out-of-sample scenarios
* `xVal` - x values of first-stage solution
* `zVal` - z values of first-stage solution
* `s` - scenario id
* `CG` - if we want to solve the second-stage with column generation
* `HEUR` - if we want to use the heuristic CG method
* `TRANSIT` - if we want to generate only fixed-transit subpaths
### Returns
* `sol` - Solution
"""
function out_of_sample_test(inst::InstanceSettingData, xVal, CG::Bool, HEUR::Bool, TRANSIT::Bool, num_focus::Int=NUM_FOCUS)
    R, m = buildRouteSettingData(inst);
    all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs = generateRestrictedSubPathSet(m, R, inst, xVal, CG, TRANSIT);
    if CG
        LB, UB, sol = solve_second_stage(m, R, inst, xVal, all_subpaths, subpath_road_networks, all_subpath_graphs, all_load_expanded_graphs, CG, HEUR, num_focus)
        return sol, R, m
    end
    model = directTwoStageModelWithFixFirstStage(R, all_subpaths, inst,xVal)
    optimize!(model)
    sol = buildSol(m, model, R, inst, all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs)
    println(sol.cost)
    return sol, R, m
end

"""
Lazy cut callback implementation (single tree)
"""

function solveSubproblem(
    l,
    t,
    s, 
    inst, 
    R,
    CG,
    Heur,
    all_subpath_graphs,
    all_load_expanded_graphs,
    xVal, 
    zVal, 
    subproblem,
    Subpaths, 
    timelimit
    )

    @unpack numL, numS, TDmp, TDsp, K, M, δ, Fleet, Theta  = inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times, Trip_passengers, Freq_Overlap  = R

    # println("$l, $t, $s")
    startCG = time()
    # update RHS with first stage solution
    set_normalized_rhs(subproblem[:c1], (xVal[l, t])) # + ϵ*x0[l,t]))
    set_normalized_rhs(subproblem[:c2], (xVal[l, t])) # + ϵ*x0[l,t]))
    for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips
        if (l_,t_) == (l,t)
            set_normalized_rhs(subproblem[:c4][p,(l,t)], (zVal[s,p,(l,t)])) # + ϵ*z0[l,t,k,h,s]))
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
            optimize!(subproblem)
            solvetime = solve_time(subproblem)
            # solveTSP += solvetime
            if termination_status(subproblem) != MOI.OPTIMAL
                println("Term status of SP $l $t $s: ", termination_status(subproblem))
            end
            objSP = objective_value(subproblem)
            # println(" CG IT $cgIt, OBJ= ", objSP)
            # obtain dual solution
            # d1 = round(dual(subproblems[l][t][s][:c1]), digits=4)
            # d2 = round(dual(subproblems[l][t][s][:c2]), digits=4)
            # d3 = round.(dual.(subproblems[l][t][s][:c3]), digits=4)
            d1 = dual(subproblem[:c1])
            d2 = dual(subproblem[:c2])
            d3 = dual.(subproblem[:c3])
            d4 = dual.(subproblem[:c4])
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
                                addSubPath!(R, l,t,s, subproblem, subP, node1, node2)
                            end

                        end
                    end
                end
            end
            # if we are done add objective value and generate Benders cut
            if doneCG
                return objSP, d1, d2, d4
            elseif time() - startCG > timelimit
                return objSP, d1, d2, d4
            end
        end
    else
        # solve SP
        optimize!(subproblem)
        # solvetime = solve_time(SP[l][t][s])
        # solveTSP += solvetime
        objSP = objective_value(subproblem)
        push!(objSPS[s], objSP)
        if true #objSP > thetaVal[l, t, s] + 0.001
            return objSP, d1, d2, d4
        end
    end
    return -1
end
function runAlgWithLazyCuts(
    R::RouteSettingData,
    inst::InstanceSettingData,
    Subpaths::Vector{Vector{Vector{Pool}}},                # pool of sub-paths (precomputed)
    all_subpath_graphs::Array{Vector{Vector{Vector{TEN}}},2},      # pool of time-expanded networks (used to generate sub-paths)
    all_load_expanded_graphs::Vector{Vector{TSLGraph}},              # pool of time-space-load networks (used to add sub-paths in the sub-problem)
    CG::Bool,      # solve second-stage with column generation
    # PO::Bool,      # add pareto-optimal cuts
    Heur::Bool,    # heuristic CG or not
    normalized::Bool,
    num_focus::Int64=NUM_FOCUS,
    num_threads::Int=NUM_THREADS,
    mip_gap::Float64=MIP_GAP,
    timelimit::Int=MP_TIME_LIMIT, # time limit for the algorithm execution
    )
    @unpack numL, numS, TDmp, TDsp, K, M, δ, Fleet, Theta  = inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times, Trip_passengers, Freq_Overlap  = R

    MP = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "TimeLimit" => timelimit, "MIPGap" => mip_gap, "Threads" => num_threads, "OutputFlag" => 0, "NumericFocus" => num_focus #=, "Method" => 1 "FeasibilityTol" => 1e-8=#)) # "Cuts" => 3

    JuMP.@variables MP begin
        # if passenger p is assigned to trip (l,t) in scenario s
        z[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], Bin                             
        # if line l operates at frequency t
        x[l in 1:numL, eachindex(R.Lines[l].freq)], Bin
        # recourse for line l, freq t, and scenario s
        theta[l in 1:numL, eachindex(R.Lines[l].freq), 1:numS]
    end

    JuMP.@expression(MP,
        first_stage_costs,
        sum( WEIGHT_LINE_COST*R.Lines[l].cost * x[l, t] for l in 1:numL, t in eachindex(R.Lines[l].freq))
    )

    JuMP.@expression(MP,
        second_stage_costs,
        sum(R.Pi[s] * theta[l,t,s] for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS)
    )

    JuMP.@objective(MP, Min,
        first_stage_costs +
        second_stage_costs
    )

    JuMP.@constraints MP begin
        fleet[t in 1:Num_freq], sum(x[l, t_] for l in 1:numL, t_ in Freq_Overlap[l,t]) <= Fleet
        cover[s in 1:numS, p in eachindex(Demand[s])], sum(z[s,p,(l,t)] for (l,t) in Demand[s][p].candidateTrips) <= 1
        minLoad[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(Demand[s][p].num*z[s,p,(l,t)] for p in Trip_passengers[l,t,s]) >= (1 - Theta)R.Lines[l].capacity*x[l, t]
        maxLoad[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(Demand[s][p].num*z[s,p,(l,t)] for p in Trip_passengers[l,t,s]) <= (1 + Theta)R.Lines[l].capacity*x[l, t]

        recourse[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], theta[l,t,s] >= -LARGE_NUMBER
    end


    # build second stage models (with all sub-paths in SPs)
    subproblems = [[[buildGenericSecondStage(R, Subpaths[l][t][s], inst, l,t,s, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus) for s in 1:numS] for t in eachindex(R.Lines[l].freq)] for l in 1:numL]

    done = false
    LB = typemin(Float64)       # lower bound
    UB = typemax(Float64)       # upper bound
    # LBs = Vector{Float64}()     # lower bounds at each iteration
    # UBs = Vector{Float64}()     # upper bounds at each iteration
    # cuts = ConstraintRef[]      # pool of Benders cuts
    num_cuts = 0
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

    startT = time()             # algorithm time tracker

    
    function my_callback_function(cb_data)
        # @unpack S, B, L, Lines, ODlines, TransitODs, F, Trips,  D, Cb, Pi, Horizon, OD = inst
        # status = callback_node_status(cb_data, m)
        # if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        #     # `callback_value(cb_data, x)` is not integer (to some tolerance).
        #     # If, for example, your lazy constraint generator requires an
        #     # integer-feasible primal solution, you can add a `return` here.
        #     # return
        # elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
        #     # `callback_value(cb_data, x)` is integer (to some tolerance).
        # else
        #     @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        #     # `callback_value(cb_data, x)` might be fractional or integer.
        # end
        # add colgen code here
        x_k = callback_value.(cb_data, x)
        z_k = callback_value.(cb_data, z)
        theta_k = callback_value.(cb_data, theta)

        for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS

            objSP, d1, d2, d4 = solveSubproblem(l,t,s, inst, R, CG, Heur, all_subpath_graphs, all_load_expanded_graphs, x_k, z_k, subproblems[l][t][s], Subpaths, timelimit)
            
            if !Heur || (Heur && objSP > theta_k[l,t,s] + EPS)
                cRef = @build_constraint(round(d1 + d2, digits=6) * MP[:x][l, t] + sum(round(d4[p,(l,t)], digits=6) * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
                # push!(cuts, cRef)
                num_cuts += 1
                JuMP.MOI.submit(MP, MOI.LazyConstraint(cb_data), cRef)
            end
        end
        return
    end

    MOI.set(MP, MOI.LazyConstraintCallback(), my_callback_function)
    # return m

    optimize!(MP)
    # check if solution found (optimal or not)
    if termination_status(MP) != MOI.OPTIMAL
        println("Term status of MP: ", termination_status(MP))
    end
    LB = objective_bound(MP)
    # check if we have a solution 
    if has_values(MP)
        # UB = objective_value(m)
        sol.cost = value(MP[:first_stage_costs])
        #solve the integer subproblems
        xVal = value.(MP[:x])               # first stage X vars
        zVal = value.(MP[:z])               # first stage Z vars
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
    end
    # plotAlgRun(LBs, UBs, MPtimes, numCuts, numCols)
    totT = time() - startT # stop algorithm time
    println("ALG TIME: ", totT)
    # println(" MP time: ", solveTMP)
    # println( "SP+CG time: ", solveTSP)
    println("IT $it, LB = $LB, IP(UB) = $(sol.cost), numCuts = ", num_cuts, " TOTAL T: ", totT)
    return LB, sol.cost, totT, num_cuts #, num_second_stage_vars, RNtime, solveTMP, solveTSP, xVal

end

