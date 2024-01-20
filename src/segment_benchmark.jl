###########################
## SEGMENT-BASED FORMULATION ##
###########################
"""
Arc object.
This object contains the main information of an arc (arc-based model).
"""
struct Arc
    src::Int  #TEMP: source OSM id
    snk::Int  #TEMP: source OSM id
    id::Int
    cost::Float64                    # cost
    load::Int64                      # num of passengers picked up
    pax::Vector{Tuple{Int,Int,Int}} # set of  (origin, time instant (i,h), number of passengers) picked up
    checkpoints::Vector{Tuple{Int,Int}}
end
"""
Pool to define all the arcs, of a given sub-problem)
"""
mutable struct ArcPool
    all::Vector{Arc}                               # set of all arcs
    in::Vector{Vector{Int}}                     # set of arcs ids incoming to node i
    out::Vector{Vector{Int}}                    # set of arcs ids outgoing from node i
    pax::Dict{Tuple{Int,Int},Vector{Int}}       # set of arcs picking up pax at (k,h)
    ref_arcs::Vector{Vector{Int}}               # set of arcs finishing in a node corresponding to a reference stop, and reference visit time (schedule)
end
ArcPool(nn::Int) = ArcPool(Vector{Arc}(), [Vector{Int}() for n in 1:nn], [Vector{Int}() for n in 1:nn], Dict{Tuple{Int,Int},Vector{Int}}(), Vector{Vector{Int}}())
"""
generates the arc-based model

### Keywords
* `R` - route setting data
### Returns
* `As` set of arcs to build the sub-problems
"""
function generateArcSet(
    m::MapData,
    R::RouteSettingData,        # input
    Inst::InstanceSettingData,
    )
    @unpack Lines, All_stops, Num_freq, Walk_graph, Demand, Taxi_travel_times = R
    @unpack numL, numS, K, MaxDev, MaxWalk, MaxWait, MaxSchDev, routeSpeedFactor, TDmp, TDsp, M, λ, μ, σ, δ = Inst
    all_arcs = [[Vector{ArcPool}() for t in eachindex(Lines[l].freq)] for l in 1:numL]

    ref_nodes = [[[Vector{Int}() for i in 1:length(Lines[l].ref_stop_id)] for t in eachindex(Lines[l].freq)] for l in 1:numL]
    node_checkpoints = [Dict{Int, Vector{Tuple{Int,Int}}}() for l in 1:numL]

    for (l, line) in enumerate(Lines)
        @unpack all_stops, ref_stop_id, ref_stop_time, capacity, ref_route_nodes, all_route_nodes = line
        src = ref_stop_id[1]                            # start ref stop of sub-path
        idxI = findfirst(isequal(src), ref_route_nodes)
        snk = ref_stop_id[end]                            # end ref stop of sub-path
        idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
        idxJ += idxI
        line_network = buildLineNetwork(m, ref_route_nodes[idxI:idxJ], Set(all_route_nodes), all_stops, TDsp, routeSpeedFactor, MaxDev)
        # divide nodes per checkpoint pairs
        for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
            src_ = ref_stop_id[i]                            # start ref stop of sub-path
            idxI_ = findfirst(isequal(src_), ref_route_nodes)
            snk_ = ref_stop_id[j]                            # end ref stop of sub-path
            idxJ_ = findfirst(isequal(snk_), ref_route_nodes[idxI_+1:end])
            idxJ_ += idxI_
            for ref_route_node in ref_route_nodes[idxI_:idxJ_]
                ind = getLocRangeNodes(m, m.nodes[ref_route_node], MaxDev, Set(all_stops))
                for loc_node in ind
                    if haskey(node_checkpoints[l], loc_node)
                        push!(node_checkpoints[l][loc_node], (src_,snk_))
                    else
                        node_checkpoints[l][loc_node] = [(src_,snk_)]
                    end
                end
            end
        end
        for (k,v) in node_checkpoints[l]
            node_checkpoints[l][k] = unique!(v)
        end
        for t in eachindex(line.freq)
            ti = ref_stop_time[src][t]
            tj = ref_stop_time[snk][t]
            maxT = floor(Int, Dates.value(Dates.Second(tj - ti))/TDsp)
            
            n2p = Vector{Tuple{Int,Int,Int}}()
            push!(n2p, (0, 0, 0))
            p2n = zeros(Int, line_network.V, maxT + 1, capacity + 1)
            countN = 1
            for i in eachindex(line_network.B2A)
                for t in 0:maxT
                    for c in 0:capacity
                        countN += 1
                        p2n[i, t+1, c+1] = countN
                        push!(n2p, (i, t, c))
                    end
                end
            end
            countN += 1
            push!(n2p, (0, 0, 0))
            load_expanded_graph = TSLGraph(countN, n2p, p2n)
            # populate ref_arcs
            for (idx, n) in enumerate(ref_stop_id)
                ref_time = floor(Int, Dates.value(Dates.Second(ref_stop_time[n][t] - ti))/TDsp)
                x = line_network.A2B[n]
                for c in 0:capacity
                    node = load_expanded_graph.p2n[x, ref_time+1, c+1]
                    push!(ref_nodes[l][t][idx], node)
                end
            end
            generic_graph = buildGenericTEN(line_network, maxT)
            DFSSort!(generic_graph) # sort the DAG in topological order
            for s in 1:numS
                g = deepcopy(generic_graph)
                segment_demand = splitPassengersBySegments(m, l, line, t, Demand[s])
                pool = ArcPool(load_expanded_graph.V)
                count = 0
                for p in eachindex(Demand[s])
                    pool.pax[p] = Vector{Int}()
                end
                line_demand = Vector{Tuple{Int,Int}}()
                for idxR in 1:length(ref_stop_id)-1
                    append!(line_demand, segment_demand[idxR])
                end
                g.inA, g.outA = updateGenTENwithDemand(m, R, Inst, ti, tj, line_demand, g, line_network, l,t,s, NORMALIZED)
                for i in g.order
                    if i != 1 && i != g.V
                        (n1, t1) = g.n2p[i] # (x,t) in g (TDsp)
                        # get osm id of n1
                        src_osm = line_network.B2A[n1]
                        # get checkpoint pairs of origin node
                        src_checkpoints = node_checkpoints[l][src_osm]
                        for (j, cus, time) in g.outA[i]
                            if j != g.V
                                (n2, t2) = g.n2p[j]
                                # get osm id of n2
                                snk_osm = line_network.B2A[n2]
                                # get checkpoint pairs of dest node
                                snk_checkpoints = node_checkpoints[l][snk_osm]
                                # get arc checkoints
                                arc_chkp_osm = intersect(src_checkpoints, snk_checkpoints)
                                arc_chkp = Vector{Tuple{Int, Int}}()
                                for (src_, snk_) in arc_chkp_osm
                                    t1_ = floor(Int, Dates.value(Dates.Second(ref_stop_time[src_][t]  - ti))/TDsp)
                                    t2_ = floor(Int, Dates.value(Dates.Second(ref_stop_time[snk_][t]  - ti))/TDsp)
                                    if  t1 > t1_ - EPS && t2 < t2_ + EPS
                                        push!(arc_chkp, (src_,snk_))
                                    end
                                end
                                # check all combinations of customers
                                comb = collect(powerset(cus, 0, capacity))
                                for cc in comb
                                    Q = 0
                                    cost = 0.
                                    for c_ in cc
                                        ic = findfirst(isequal(c_), cus)
                                        Q += c_[3]
                                        cost += time[ic]
                                    end
                                    if Q <= capacity
                                        # add arcs for all C levels
                                        for c in 0:capacity-Q
                                            count += 1
                                            node1 = load_expanded_graph.p2n[n1, t1+1, c+1]
                                            node2 = load_expanded_graph.p2n[n2, t2+1, c+Q+1]
                                            arc = Arc(src_osm, snk_osm, count, cost, Q, cc, arc_chkp)
                                            push!(pool.all, arc)
                                            push!(pool.in[node2], count)
                                            push!(pool.out[node1], count)
                                            for (k, h, num) in cc
                                                push!(pool.pax[k, h], count)
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                # add sink arcs (never would have customers)
                node2 = load_expanded_graph.V
                for (i, cus, time) in g.inA[g.V]
                    (n1, t1) = g.n2p[i] # (x,t) in g (TDsp)
                    # add arcs for all C levels
                    for c in 0:capacity
                        count += 1
                        node1 = load_expanded_graph.p2n[n1, t1+1, c+1]
                        arc = Arc(0,0,count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Tuple{Int, Int}}())
                        push!(pool.all, arc)
                        push!(pool.in[node2], count)
                        push!(pool.out[node1], count)
                    end
                end
                # add source arcs (never would have customers)
                node1 = 1
                for (j, cus, time) in g.outA[1]
                    (n2, t2) = g.n2p[j]
                    # add arcs for only c=0
                    count += 1
                    node2 = load_expanded_graph.p2n[n2, t2+1, 1]
                    arc = Arc(0,0,count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Tuple{Int, Int}}())
                    push!(pool.all, arc)
                    push!(pool.in[node2], count)
                    push!(pool.out[node1], count)
                end
                for r_s in 1:length(ref_stop_id)
                    ref_stop_arcs = Vector{Int}()
                    for n in ref_nodes[l][t][r_s]
                        for arc_id in pool.in[n]
                            push!(ref_stop_arcs, arc_id)
                        end
                    end
                    push!(pool.ref_arcs, ref_stop_arcs)
                end
                push!(all_arcs[l][t], pool)
            end
        end
    end
    return all_arcs
end
"""
builds the arc-based model

### Keywords
* `R` - route setting data
* `time_limit_sec` - time limit to solve the problem
* `mip_gap` - accepted optimality gap
### Returns
* JuMP model of the sub-problem
"""
function directArcBasedModel(
    R::RouteSettingData,
    Arcs::Vector{Vector{Vector{ArcPool}}},
    Inst::InstanceSettingData,
    time_limit_sec::Int64=FULL_MIP_TIME_LIMIT,
    mip_gap::Float64=MIP_GAP,
    num_threads::Int64=NUM_THREADS,
    )

    @unpack numL, numS, TDmp, TDsp, Fleet, Theta, K = Inst # master and subproblem time discretization, and normalized (true/false) in-vehicle time and delay
    @unpack Num_freq, Demand, Trip_passengers, Freq_Overlap = R
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "TimeLimit" => time_limit_sec, "MIPGap" => mip_gap, "Threads" => num_threads, "OutputFlag" => 0))

    checkpointPairsIn = [Dict{Int, Vector{Int}}() for l in 1:numL]
    checkpointPairsOut = [Dict{Int, Vector{Int}}() for l in 1:numL]
    checkpointPairs = [Vector{Tuple{Int, Int}}() for l in 1:numL]

    for l in 1:numL
        for i in 1:length(R.Lines[l].ref_stop_id)-1, j in i+1:min(i+K,length(R.Lines[l].ref_stop_id))
            n1 = R.Lines[l].ref_stop_id[i]
            n2 = R.Lines[l].ref_stop_id[j]
            if haskey(checkpointPairsIn[l], n2)
                push!(checkpointPairsIn[l][n2], n1)
            else
                checkpointPairsIn[l][n2] = [n1]
            end
            if haskey(checkpointPairsOut[l], n1)
                push!(checkpointPairsOut[l][n1], n2)
            else
                checkpointPairsOut[l][n1] = [n2]
            end
            push!(checkpointPairs[l], (n1, n2))
        end
    end 

    JuMP.@variables m begin
        # if passenger p is assigned to trip (l,t) in scenario s
        z[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], Bin                             
        # if line l operates at frequency t
        x[l in 1:numL, eachindex(R.Lines[l].freq)], Bin
        # 1 if subpath l,t,s,i is selected
        y[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, i in eachindex(Arcs[l][t][s].all)], Bin
        # 1 if there is a subpath between checkpoints u and v
        b[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS,  c in checkpointPairs[l]], Bin
    end

    JuMP.@expression(m,
        first_stage_costs,
        sum( WEIGHT_LINE_COST*R.Lines[l].cost * x[l, t] for l in 1:numL, t in eachindex(R.Lines[l].freq))
    )

    JuMP.@expression(m,
        second_stage_costs,
        # generalized cost of travel of passengers served by arcs/segments
        sum( R.Pi[s] * (
            sum( a.cost * y[l, t, s, a.id] for a in Arcs[l][t][s].all))
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
        c1[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(y[l,t,s,id] for id in Arcs[l][t][s].out[1]) == x[l, t]
        c2[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(y[l,t,s,id] for id in Arcs[l][t][s].in[length(Arcs[l][t][s].in)]) == x[l, t]
        c3[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, i in 2:length(Arcs[l][t][s].in)-1], sum(y[l,t,s,id] for id in Arcs[l][t][s].out[i]) - sum(y[l,t,s,id] for id in Arcs[l][t][s].in[i]) == 0
        c4[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], sum(y[l,t,s,id] for id in Arcs[l][t][s].pax[p]) <= z[s,p,(l,t)]
        #
        c5[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(b[l,t,s,(R.Lines[l].ref_stop_id[1],c)] for c in checkpointPairsOut[l][R.Lines[l].ref_stop_id[1]]) == x[l, t] 
        c6[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(b[l,t,s,(c,R.Lines[l].ref_stop_id[end])] for c in checkpointPairsIn[l][R.Lines[l].ref_stop_id[end]]) == x[l, t]
        c7[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, u in R.Lines[l].ref_stop_id[2:end-1]], sum(b[l,t,s,(u,c)] for c in checkpointPairsOut[l][u]) - sum(b[l,t,s,(c,u)] for c in checkpointPairsIn[l][u]) == 0

        c8[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, (u,v) in checkpointPairs[l]], sum(y[l,t,s,id] for id in Arcs[l][t][s].ref_arcs[findfirst(isequal(v), R.Lines[l].ref_stop_id)]) >= b[l,t,s,(u,v)]
        c9[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, a in Arcs[l][t][s].all; a.src != 0 && a.snk != 0], y[l,t,s,a.id] <= sum(b[l,t,s,c] for c in a.checkpoints)
    end
    return m
end

function testArcBasedModel(numL::Int, 
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

    inst = InstanceSettingData(numL, demand_horizon, numS, trainData, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,MAX_SCHEDULE_DEV,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    R, m = buildRouteSettingData(inst);
    all_arcs = generateArcSet(m, R, inst);
    model = directArcBasedModel(R, all_arcs, inst);
    optimize!(model)
    println("OBJ: ", objective_value(model), " in solve time: ", solve_time(model))

    xVal = value.(model[:x])
    yVal = value.(model[:y])
    for l in 1:numL, t in eachindex(R.Lines[l].freq)
        if xVal[l,t] > EPS
            for s in 1:numS
                for i in eachindex(all_arcs[l][t][s].all)
                    if yVal[l,t,s,i] > EPS
                        if length(all_arcs[l][t][s].all[i].pax) > EPS
                            println("Line $l, time $t, scen $s, serves: ", all_arcs[l][t][s].all[i].pax)
                        end
                    end
                end
            end
        end
    end
    return model
end