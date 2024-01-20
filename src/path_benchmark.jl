############################
## PATH-BASED FORMULATION ##
############################
"""
Path object.
This object contains the main information of a path.
"""
struct Path
    id::Int
    cost::Float64                    # cost
    load::Int64                      # num of passengers picked up
    pax::Vector{Tuple{Int,Int,Int}} # set of  (origin, time instant (i,h), number of passengers) picked up
end
"""
Pool to define all the arcs (i.e., path), of a given sub-problem
"""
mutable struct PathPool
    all::Vector{Path}                               # set of all paths
    pax::Dict{Tuple{Int,Int},Vector{Int}}       # set of paths picking up pax at (k,h)
end
PathPool() = PathPool(Vector{Path}(), Dict{Tuple{Int,Int},Vector{Int}}())
"""
builds the path-based model

### Keywords
* `R` - route setting data
* `time_limit_sec` - time limit to solve the problem
* `mip_gap` - accepted optimality gap
### Returns
* JuMP model of the sub-problem
"""
function generatePaths(
    R::RouteSettingData,
    l::Int,
    t::Int,
    s::Int,
    Subpaths::Pool,
    )

    snk = length(Subpaths.in)
    candPool = Vector{Tuple{Path, Int}}()
    count = 0
    pathCount = 0
    pathPool = PathPool()
    for (k,h) in eachindex(R.Demand[s])
        pathPool.pax[(k,h)] = Vector{Path}()
    end
    for id in Subpaths.out[1]
        a = Subpaths.all[id]
        count += 1
        p = Path(count,a.cost, a.load, a.pax)
        push!(candPool, (p, a.j))
    end
    done = false
    limit = false
    while !done && !limit
        done = true
        for idx in eachindex(candPool)
            if candPool[idx][2] != snk
                el = popat!(candPool, idx)
                (p, i) = el
                done = false
                for id in Subpaths.out[i]
                    # add a to p as a new path
                    a = Subpaths.all[id]
                    count += 1
                    new_cost = p.cost + a.cost
                    new_load = p.load + a.load
                    if new_load > R.Lines[l].capacity + 0.01
                        println("Capacity exceeded")
                    end
                    new_pax = union(p.pax, a.pax)
                    new_path = Path(count, new_cost, new_load, new_pax)
                    push!(candPool, (new_path, a.j))
                    if a.j == snk
                        pathCount += 1
                        for (k,h,num) in new_pax
                            push!(pathPool.pax[k,h], pathCount)
                        end
                        path = Path(pathCount, new_cost, new_load, new_pax)
                        push!(pathPool.all, path)
                    end
                    if count > 1000000 # set limit of paths
                        limit = true
                        break
                    end
                end
            end
        end
    end
    if limit
        # add every path to pool (all are valid)
        for el in candPool
            (p,i) = el
            @unpack cost, load, pax = p
            pathCount += 1
            for (k,h,num) in pax
                push!(pathPool.pax[k,h], pathCount)
            end
            path = Path(pathCount, cost, load, pax)
            push!(pathPool.all, path)
        end
    end
    if length(pathPool.all) < 0.01
        println("No paths found, count = ", count)
        # add empty path
        pathCount += 1
        path = Path(pathCount, 0., 0, Vector{Tuple{Int,Int,Int}}())
        push!(pathPool.all, path)
    end
    return pathPool
end
"""
Build direct path-based two-stage stochastic model.

### Keywords
* `R` - RouteSettingData
* `Ps` - Pool of paths
* `Pax` - Set of paths covering each passenger
* `time_limit_sec` - maximum solve time in seconds
* `mip_gap` - epsilon tolerance
### Returns
* JuMP Model
"""
function directTwoStageModelPathBased(
    R::RouteSettingData,
    Paths::Vector{Vector{Vector{PathPool}}},
    Inst::InstanceSettingData,
    time_limit_sec::Int64=FULL_MIP_TIME_LIMIT,
    mip_gap::Float64=MIP_GAP,
    num_threads::Int64=NUM_THREADS,
    )

    @unpack numL, numS, TDmp, TDsp, Fleet, Theta = Inst # master and subproblem time discretization, and normalized (true/false) in-vehicle time and delay
    @unpack Num_freq, Demand, Trip_passengers, Freq_Overlap = R
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "TimeLimit" => time_limit_sec, "MIPGap" => mip_gap, "Threads" => num_threads, "OutputFlag" => 0))


    JuMP.@variables m begin
        # if passenger p is assigned to trip (l,t) in scenario s
        z[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], Bin                             
        # if line l operates at frequency t
        x[l in 1:numL, eachindex(R.Lines[l].freq)], Bin
        # 1 if subpath l,t,s,i is selected
        y[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, i in eachindex(Paths[l][t][s].all)], Bin
    end

    JuMP.@expression(m,
        first_stage_costs,
        sum( WEIGHT_LINE_COST*R.Lines[l].cost * x[l, t] for l in 1:numL, t in eachindex(R.Lines[l].freq))
    )

    JuMP.@expression(m,
        second_stage_costs,
        # generalized cost of travel of passengers served by paths
        sum( R.Pi[s] * (
            sum( p.cost * y[l, t, s, p.id] for p in Paths[l][t][s].all))
        for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS )
    )

    JuMP.@objective(m, Min,
        first_stage_costs +
        second_stage_costs
    )

    JuMP.@constraints m begin
        # first stage
        fleet[t in 1:Num_freq], sum(x[l, t_] for l in 1:numL, t_ in Freq_Overlap[l,t]) <= Fleet
        cover[s in 1:numS, p in eachindex(Demand[s])], sum(z[s,p,(l,t)] for (l,t) in Demand[s][p].candidateTrips) <= 1
        minLoad[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(Demand[s][p].num*z[s,p,(l,t)] for p in Trip_passengers[l,t,s]) >= (1 - Theta)R.Lines[l].capacity*x[l, t]
        maxLoad[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(Demand[s][p].num*z[s,p,(l,t)] for p in Trip_passengers[l,t,s]) <= (1 + Theta)R.Lines[l].capacity*x[l, t]

        # second-stage
        c1[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(y[l,t,s,p.id] for p in Paths[l][t][s].all) == x[l, t]
        c2[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], sum(y[l,t,s,id] for id in Paths[l][t][s].pax[p]) <= z[s,p,(l,t)]

    end

    return m
end

function testPathBasedModel(
    numL::Int, 
    demand_horizon::Int, 
    numS::Int, 
    K::Int, 
    trainData::Bool, 
    TRANSIT::Bool, 
    maxDev::Float64, 
    maxWalk::Float64, 
    capacity::Int, 
    fleet_size::Int,
    )

    CG = false # in direct models, we only do full enumeration
    inst = InstanceSettingData(numL, demand_horizon, numS, trainData, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,MAX_SCHEDULE_DEV,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    # inst = InstanceSettingData(numL, demT, numS, K)
    R, m = buildRouteSettingData(inst);
    all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs = generateSubPathSet(m, R, inst, CG, TRANSIT);
    all_paths = [[Vector{PathPool}() for t in eachindex(R.Lines[l].freq)] for l in eachindex(R.Lines)]
    for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS
        pool = generatePaths(R, l,t,s, all_subpaths[l][t][s])
        push!(all_paths[l][t], pool)
    end
    model = directTwoStageModelPathBased(R, all_paths, inst)
    optimize!(model)
    obj = objective_value(model)
    println("OBJ: ", obj, " in time: ", solve_time(model))

end

function testRunPathBased(case::Int)
    time = Dates.format(now(), "HH-MM")
    filename = string("testFull_$case-", time)
    setL = [2, 10, 50, 99]
    setT = [2, 8, 16]
    setS = [10, 100]#, 1000]
    setK = [1]#, 2]
    setTDmp = [180] #, 60, 180]
    totCases = length(setL)*length(setT)*length(setK)*length(setS)*length(setTDmp) #*length(setCG)*length(setPO)
    M = reshape(1:totCases, length(setL), length(setT), length(setS), length(setK), length(setTDmp)) #length(setCG), length(setPO)
    idx = findfirst(isequal(case), M)
    L = setL[idx[1]]
    T = setT[idx[2]]
    S = setS[idx[3]]
    K = setK[idx[4]]
    TDmp = setTDmp[idx[5]]
    filename = string( filename, "_P.txt")
    f = open(filename, "w")
    obj, solveT, preC, numPaths = testPathBasedModel(L, T, S, K, TDmp)
    write(f, "$L\t$T\t$S\t$K\t$TDmp\t$obj\t$solveT\t$preC\t$numPaths\n")
    flush(f)
    close(f)
end