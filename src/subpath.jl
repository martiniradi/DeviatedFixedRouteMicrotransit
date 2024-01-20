"""
Sub-path object.
This object contains the main information of a sub-path.
"""
struct SubPath
    id::Int                     # identifier
    cost::Float64               # cost (g_a)
    load::Int64                 # num of passengers picked up
    pax::Vector{Tuple{Int,Int,Int}}            # vector of AggregPassenger ids
    graph_path::Vector{Int}       # path in the time-expanded network
    i::Int                      # tail node in load-expanded network
    j::Int                      # head node in load-expanded network
end
"""
Pool to define all the arcs (i.e., sub-path), of a given sub-problem (TSLN)
"""
mutable struct Pool
    all::Vector{SubPath}                        # set of all sub-paths
    in::Vector{Vector{Int}}                     # set of sub-path ids incoming to node i
    out::Vector{Vector{Int}}                    # set of sub-path ids outgoing from node i
    pax::Dict{Tuple{Int,Int},Vector{Int}}       # set of sub-path ids picking up pax at (k,h)
end
Pool(nn::Int) = Pool(Vector{SubPath}(), [Vector{Int}() for n in 1:nn], [Vector{Int}() for n in 1:nn], Dict{Tuple{Int,Int},Vector{Int}}())

"""
run a label setting algorithm on graph g (TEN), with two labels: customers picked up, and stopped picked up (i.e., saved waiting time).

### Keywords
* `g` - TEN
* source` - source node
* `nodeOrder` - topological order of nodes
* `capacity` - maximum number of passengers to pick up
### Returns
* non-dominated labels at each node, and the parent node and label of each node
"""
function fullEnumAlg(
    g::TEN,
    source::Int64,
    nodeOrder::Vector{Int64},
    capacity::Int,
    )
    # take all combinations (no-dominance criteria)
    nn = g.V
    labels = [Vector{Pair{Vector{Tuple{Int, Int, Int}}, Vector{Float64}}}([Pair(Vector{Tuple{Int, Int, Int}}(), Vector{Int}([1000000.]))]) for n in 1:nn]
    labels[source] = [Pair(Vector{Tuple{Int, Int, Int}}(), Vector{Float64}())]
    # keep track of prev of each label
    prev = [Vector{Tuple{Int, Pair{Vector{Tuple{Int, Int, Int}}, Vector{Float64}}}}([(0, Pair(Vector{Tuple{Int, Int, Int}}(), Vector{Float64}()))]) for n in 1:nn]
    for i in nodeOrder
        if i != source
            for (src, cus, time) in g.inA[i]
                for label in labels[src]
                    # ignore nodes that cannot be reached
                    if sum(label.second) < LARGE_NUMBER
                        # compute num of pax at label
                        numCust = 0
                        for c_ in label.first
                            numCust += c_[3]
                        end
                        # initialize vector for new labels, which will be label + all comb of the customers in the arc
                        newLabels = Vector{Pair{Vector{Tuple{Int, Int, Int}}, Vector{Float64}}}()
                        comb = collect(powerset(cus, 0, capacity-numCust)) # we limit the number of elements to the capacity - passengers already in label (in worst case there is one pax per element)
                        for cus_ in comb
                            #NOTE there will be many of the combinations that will result in the same label
                            # compute labels for dst
                            C = deepcopy(label.first)
                            WT = deepcopy(label.second)
                            # check if they can be picked up
                            fit = true
                            Q_ = deepcopy(numCust)
                            for (ic, c) in enumerate(cus_)
                                # if customer not already picked up
                                if !(c in C)
                                    Q_ += cus_[ic][3]
                                    if Q_ > capacity + EPS
                                        fit = false
                                    end
                                end
                            end
                            if fit
                                # generate valid new label
                                for (ic, c) in enumerate(cus_)
                                    icC = findfirst(isequal(c), cus) # actual index of c in cus (to retreive the corresponding time)
                                    # if it is already picked up, check if t is better
                                    if c in C
                                        icS = findfirst(isequal(c), C)
                                        # we update the time if better (dual values do not depend on where the customer is picked up)
                                        if time[icC] <  WT[icS] - EPS
                                            WT[icS] = time[icC]
                                        end
                                    else
                                        push!(C, cus_[ic])
                                        push!(WT, time[icC])
                                    end
                                end
                                # there are probably many labels with same customers (but some of them may have different waiting times, depending on the predecessors)
                                # we can discard the ones with same pax but worse times directly (dominated)
                                push!(newLabels, Pair(C, WT))
                            end
                        end
                        # for the new labels created updated the ones in label[i]
                        # add new labels (only if they are not identical to exisiting ones)
                        for nLab in newLabels
                            C = nLab.first
                            WT = nLab.second
                            dominated = false
                            for (idx, labY) in enumerate(labels[i])
                                if labY.first == C && labY.second == WT
                                    dominated = true
                                end
                            end
                            if dominated == false
                                push!(labels[i], nLab)
                                push!(prev[i], (src, label))
                            end
                        end
                    end
                end
            end
        end
    end
    return labels, prev
end

function fullLabelSettingAlgInCG(
    l::Int,
    t::Int,
    g::TEN,
    source::Int64,
    nodeOrder::Vector{Int64},
    capacity::Int,
    d4,
    )
    # We need to keep track of all possible combinations
    nn = g.V
    labels = [Vector{Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}([(Vector{Tuple{Int, Int, Int}}(), Vector{Float64}(), typemax(Float64))]) for n in 1:nn]
    labels[source] = [(Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), 0.)]
    # keep track of prev of each undominated label
    prev = [Vector{Tuple{Int, Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}}([(0, (Vector{Tuple{Int, Int, Int}}(), Vector{Float64}(), 0.))]) for n in 1:nn]
    for i in nodeOrder
        if i != source
            for (src, cus, time) in g.inA[i] 
                for label in labels[src]
                    if label[3] < LARGE_NUMBER
                        numCust = 0
                        for c_ in label[1]
                            numCust += c_[3]
                        end
                        newLabels = Vector{Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}()
                        comb = collect(powerset(cus, 0, capacity-numCust)) 
                        existingCombs = Vector{Vector{Tuple{Int, Int, Int}}}()
                        for cus_ in comb
                            C = deepcopy(label[1])
                            WT = deepcopy(label[2])
                            RC = deepcopy(label[3])
                            # check if they can be picked up
                            fit = true
                            Q_ = deepcopy(numCust)
                            for (ic, c) in enumerate(cus_)
                                if !(c in C)
                                    Q_ += cus_[ic][3]
                                    if Q_ > capacity
                                        fit = false
                                    end
                                end
                            end
                            if fit
                                # generate valid new label
                                for (ic, c) in enumerate(cus_)
                                    icC = findfirst(isequal(c), cus) # actual index of c in cus (to retreive the corresponding time)
                                    # if it is already picked up, check if t is better
                                    if c in C
                                        icS = findfirst(isequal(c), C)
                                        # we update the time if better (dual values do not depend on where the customer is picked up)
                                        if time[icC] <  WT[icS] - EPS
                                            RC -= WT[icS] - time[icC]
                                            WT[icS] = time[icC]
                                        end
                                    else
                                        push!(C, cus_[ic])
                                        push!(WT, time[icC])
                                        # update RC with new customers cus_
                                        RC -= -time[icC] + d4[(cus_[ic][1], cus_[ic][2]),(l,t)]
                                    end
                                end
                                Csorted = sort(C, by = x -> x[1])
                                if !(Csorted in existingCombs)
                                    push!(newLabels, (C, WT, RC))
                                    push!(existingCombs, Csorted)
                                end
                            end
                        end

                        for nLab in newLabels
                            (C, WT, RC) = nLab
                            toDelete = Set{Int}()
                            dominated = false
                            for (idx, labY) in enumerate(labels[i])
                                cLC = issubset(labY[1], C)
                                cCL = issubset(C, labY[1])
                                if cLC && cCL
                                    if labY[3] < RC + EPS
                                        # keep LabY
                                        dominated = true
                                    else
                                        push!(toDelete, idx)
                                    end
                                end
                            end
                            if dominated == false
                                toDel = sort(collect(toDelete))
                                deleteat!(labels[i], toDel)
                                deleteat!(prev[i], toDel)
                                push!(labels[i], nLab)
                                push!(prev[i], (src, label))
                            end
                        end
                    end
                end
            end
        end
    end
    return labels, prev
end

function labelSettingAlgInCG(
    l::Int,
    t::Int,
    g::TEN,
    source::Int64,
    nodeOrder::Vector{Int64},
    capacity::Int,
    d4,
    )
    nn = g.V
    labels = [Vector{Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}([(Vector{Tuple{Int, Int, Int}}(), Vector{Float64}(), typemax(Float64))]) for n in 1:nn]
    labels[source] = [(Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), 0.)]
    # keep track of prev of each undominated label
    prev = [Vector{Tuple{Int, Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}}([(0, (Vector{Tuple{Int, Int, Int}}(), Vector{Float64}(), 0.))]) for n in 1:nn]
    for i in nodeOrder
        if i != source
            for (src, cus, time) in g.inA[i] # cus should be vector (there could be more than one customer related to a node)
                for label in labels[src]
                    # compute labels for dst
                    C = deepcopy(label[1])
                    WT = deepcopy(label[2])
                    RC = deepcopy(label[3])
                    Q = 0
                    for c_ in C
                        Q += c_[3]
                    end
                    # check customers in c
                    for (ic, c) in enumerate(cus)
                        # if it is already picked up, check if t is better
                        if c in C
                            icS = findfirst(isequal(c), C)
                            if time[ic] <  WT[icS] - 0.01
                                RC -= WT[icS] - time[ic]
                                WT[icS] = time[ic]
                            end
                        else
                            # if it is new: add it if the vehicle has capacity
                            # we always pickup all possible pax in a stop (it means they can make it)
                            if cus[ic][3] + Q <= capacity
                                redcost = -time[ic] + d4[(cus[ic][1], cus[ic][2]),(l,t)]
                                if redcost > 0.01
                                    push!(C, cus[ic])
                                    push!(WT, time[ic])
                                    RC -= redcost
                                end
                            end
                        end
                    end
                    # we keep the label with the lowest RC
                    toDelete = Set{Int}()
                    dominated = false
                    for (idx, labY) in enumerate(labels[i])
                        RC_i = labY[3]
                        if RC < RC_i - 0.001
                            # delete RC_i
                            push!(toDelete, idx)
                        else
                            # RC dominated
                            dominated = true
                        end
                    end
                    if dominated == false
                        toDel = sort(collect(toDelete))
                        newLab = (C, WT, RC)
                        deleteat!(labels[i], toDel)
                        deleteat!(prev[i], toDel)
                        push!(labels[i], newLab)
                        push!(prev[i], (src, label))
                    end
                end
            end
        end
    end
    return labels, prev
end

"""
generate the paths from the output of the label setting algorithm,
we do this by tracking back parent nodes starting from the sink node

### Keywords
* `node` - start node (sink node)
* `src` - source node
* `labels` - set of non-dominated labels at each node
* `prevNode` - for each label and node, the parent node and label
### Returns
* the set of paths from source to sink, with their corresponding final label
"""
function getPathsAndLabels(node::Int64, src::Int64,
                            labels::Array{Array{Pair{Array{Tuple{Int64,Int64,Int64},1},Array{Float64,1}},1},1},
                            prevNode::Array{Array{Tuple{Int64,Pair{Array{Tuple{Int64,Int64,Int64},1},Array{Float64,1}}},1},1},
                            ) #, g, startT::Int, TD::Int, TDs::Int) #::Vector{Vector{Int64}})
    paths = Vector{Vector{Int64}}()
    labs = Vector{Pair{Vector{Tuple{Int, Int, Int}}, Vector{Float64}}}()
    for idx in eachindex(prevNode[node])
        path = Vector{Int64}()
        push!(path, node)
        (pN, pLab) = prevNode[node][idx]
        while pN != 0 
            push!(path, pN)
            i = findfirst(isequal(pLab), labels[pN])
            (pN, pLab) = prevNode[pN][i]
        end
        if path[end] == src
            push!(paths, reverse(path))
            push!(labs, labels[node][idx])
        end
    end
    return paths, labs
end
"""
### Keywords
* `node` - start node (sink node)
* `src` - source node
* `labels` - set of non-dominated labels at each node
* `prevNode` - for each label and node, the parent node and label
### Returns
* the set of paths from source to sink, with their corresponding final label
"""
function getPathsAndLabelsInCG(node::Int64, src::Int64,
    labels::Vector{Vector{Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}},
    prevNode::Vector{Vector{Tuple{Int, Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}}},
    ) 
    paths = Vector{Vector{Int64}}()
    labs = Vector{Tuple{Vector{Tuple{Int, Int, Int}}, Vector{Float64}, Float64}}()
    for idx in eachindex(prevNode[node])
        path = Vector{Int64}()
        push!(path, node)
        (pN, pLab) = prevNode[node][idx]
        while pN != 0 
            push!(path, pN)
            i = findfirst(isequal(pLab), labels[pN])
            (pN, pLab) = prevNode[pN][i]
        end
        if path[end] == src
            push!(paths, reverse(path))
            push!(labs, labels[node][idx])
        end
    end
    return paths, labs
end
"""
function used to pre-compute the set of graphs used in the algorithm, and subpaths

### Keywords
(see function arguments)
### Returns
* `all_subpaths` - pool of sub-paths for each line, freq, and scenario
* `subpath_road_networks` - road network used for subpath routing
* `all_load_expanded_graphs` - second stage load-expanded networks
* `all_subpath_graphs` - time-expanded networks for finding subpaths
"""
function generateSubPathSet(m::MapData,
    R::RouteSettingData,        # model input
    Inst::InstanceSettingData,  # instance data
    CG::Bool,                   # solve second-stage with column generation (if true, we start with an empty set of sub-paths)
    TRANSIT::Bool,
    )
    @unpack Lines, All_stops, Num_freq, Walk_graph, Demand, Taxi_travel_times = R
    @unpack K, MaxDev, MaxWalk, MaxWait, MaxSchDev, routeSpeedFactor, TDmp, TDsp, M, λ, μ, σ, δ = Inst
    
    all_subpaths = [[Vector{Pool}() for t in eachindex(Lines[l].freq)] for l in eachindex(Lines)]
    subpath_road_networks = [[Vector{MyNetwork}() for i in 1:length(Lines[l].ref_stop_id)-1] for l in eachindex(Lines)]
    generic_subpath_graphs = [[[Vector{TEN}() for i in 1:length(Lines[l].ref_stop_id)-1] for t in eachindex(Lines[l].freq)] for l in eachindex(Lines)]
    all_subpath_graphs = [[[Vector{TEN}() for i in 1:length(Lines[l].ref_stop_id)-1] for t in eachindex(Lines[l].freq)] for l in eachindex(Lines), s in eachindex(Demand)]
    all_load_expanded_graphs = [Vector{TSLGraph}() for l in eachindex(Lines)]
    enumeration_time = 0.

    for (l, line) in enumerate(Lines)
        @unpack all_stops, ref_stop_id, ref_stop_time, capacity, ref_route_nodes, all_route_nodes = line
        for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
            src = ref_stop_id[i]                            # start ref stop of sub-path
            idxI = findfirst(isequal(src), ref_route_nodes)
            snk = ref_stop_id[j]                            # end ref stop of sub-path
            idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
            idxJ += idxI
            subpath_network = TRANSIT ? buildLineNetwork(m, ref_route_nodes[idxI:idxJ], Set(all_route_nodes), all_stops, TDsp, routeSpeedFactor, 0.1) : buildLineNetwork(m, ref_route_nodes[idxI:idxJ], Set(all_route_nodes), all_stops, TDsp, routeSpeedFactor, MaxDev)
            push!(subpath_road_networks[l][i], subpath_network)
        end
        for t in eachindex(line.freq)
            subproblemT = ceil(Int, Dates.value(Dates.Second(ref_stop_time[ref_stop_id[end]][t] - ref_stop_time[ref_stop_id[1]][t]))/TDsp)
            load_expanded_graph = initializeTSLgraph(line, subproblemT, t, TDsp, REFSTOP_TIME_GAP)
            push!(all_load_expanded_graphs[l], load_expanded_graph)
            for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
                src = ref_stop_id[i]                            # start ref stop of sub-path
                idxI = findfirst(isequal(src), ref_route_nodes)
                snk = ref_stop_id[j]                            # end ref stop of sub-path
                idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
                idxJ += idxI
                ti = ref_stop_time[src][t]
                tj = ref_stop_time[snk][t]
                subpath_graphT = floor(Int, Dates.value(Dates.Second(tj-ti))/TDsp)
                # build generic TEN
                g = buildGenericTEN(subpath_road_networks[l][i][j-i], subpath_graphT)
                DFSSort!(g) # sort the DAG in topological order
                push!(generic_subpath_graphs[l][t][i], g)
            end
        end
    end

    for (l, line) in enumerate(Lines)
        @unpack all_stops, ref_stop_id, ref_stop_time, capacity, ref_route_nodes, all_route_nodes = line
        for t in eachindex(line.freq)
            t0 = ref_stop_time[ref_stop_id[1]][t]
            subproblemT = ceil(Int, Dates.value(Dates.Second(ref_stop_time[ref_stop_id[end]][t] - t0))/TDsp)
            load_expanded_graph = all_load_expanded_graphs[l][t]
            # @showprogress "Demand update for line $l time $t "
            for s in eachindex(Demand)
                # passengers to be considered in the subproblem
                # divide scenario demand to its closest segment (a segment is the space betwen two consecutive ref stops). This is to avoid a solution picking up the same passenger multiple times.
                segment_demand = splitPassengersBySegments(m, l, line, t, Demand[s])
                pool = Pool(load_expanded_graph.V)
                count = 0
                for p in eachindex(Demand[s])
                    pool.pax[p] = Vector{Int}()
                end
                for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
                    src = ref_stop_id[i]                            # start ref stop of sub-path
                    idxI = findfirst(isequal(src), ref_route_nodes)
                    snk = ref_stop_id[j]                            # end ref stop of sub-path
                    idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
                    idxJ += idxI
                    ti = ref_stop_time[src][t]
                    tj = ref_stop_time[snk][t]
                    g = deepcopy(generic_subpath_graphs[l][t][i][j-i])
                    subpath_segment_demand = Vector{Tuple{Int,Int}}()
                    for idxR in i:j-1
                        append!(subpath_segment_demand, segment_demand[idxR])
                    end
                    g.inA, g.outA = TRANSIT ? updateGenTENwithDemandForTransit(m, R, Inst, ti, tj, subpath_segment_demand, g, subpath_road_networks[l][i][j-i], l,t,s, NORMALIZED) : updateGenTENwithDemand(m, R, Inst, ti, tj, subpath_segment_demand, g, subpath_road_networks[l][i][j-i], l,t,s, NORMALIZED)
                    push!(all_subpath_graphs[l,s][t][i], g)
                    t1 = ceil(Int, Dates.value(Dates.Second(ref_stop_time[src][t]-t0))/TDsp)
                    t2 = ceil(Int, Dates.value(Dates.Second(ref_stop_time[snk][t]-t0))/TDsp)
                    emptySubpath = false # we should always guarantee an "empty" path.
                    if !CG # precompute subset of "empty" sub-paths
                        start_enum_timer = time()
                        labels, prev = fullEnumAlg(g, 1, g.order, capacity)
                        # recreate paths
                        paths, labs = getPathsAndLabels(g.V, 1, labels, prev)
                        for (idx, label) in enumerate(labs)
                            WT = sum(label.second)
                            if WT < LARGE_NUMBER 
                                Q = 0
                                for (orig, time, num) in label.first
                                    Q += num
                                end
                                cost = WT 
                                if Q == 0
                                    emptySubpath = true
                                end
                                for c in 0:capacity-Q
                                    count += 1
                                    node1 = load_expanded_graph.p2n[i, t1+1, c+1]
                                    node2 = load_expanded_graph.p2n[j, t2+1, c+Q+1]
                                    subP = SubPath(count, cost, Q, label.first, paths[idx], node1, node2)
                                    push!(pool.all, subP)
                                    push!(pool.in[node2], count)
                                    push!(pool.out[node1], count)
                                    for (k, h, num) in label.first
                                        push!(pool.pax[k, h], count)
                                    end
                                end
                            end
                        end
                        enumeration_time += time() - start_enum_timer
                    end
                    if !emptySubpath
                        for c in 0:capacity
                            count += 1
                            node1 = load_expanded_graph.p2n[i, t1+1, c+1]
                            node2 = load_expanded_graph.p2n[j, t2+1, c+1]
                            subP = SubPath(count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), node1, node2)
                            push!(pool.all, subP)
                            push!(pool.in[node2], count)
                            push!(pool.out[node1], count)
                        end
                    end
                end
                # add set of arcs from the artificial source node from the first stop in the route (empty vehicle)
                count += 1
                node1 = 1
                node2 = load_expanded_graph.p2n[1, 1, 1]
                subP = SubPath(count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), node1, node2)
                push!(pool.all, subP)
                push!(pool.in[node2], count)
                push!(pool.out[node1], count)
                # add set of arcs to the artificial sink node from the last stop in the route
                for c in 0:capacity
                    count += 1
                    node1 = load_expanded_graph.p2n[length(ref_stop_id), subproblemT+1, c+1]
                    node2 = load_expanded_graph.V
                    subP = SubPath(count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), node1, node2)
                    push!(pool.all, subP)
                    push!(pool.in[node2], count)
                    push!(pool.out[node1], count)
                end
                push!(all_subpaths[l][t], pool)
            end
        end
    end
    # println("DONE!")
    return all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs, enumeration_time
end
"""
function used to pre-compute the set of subpaths given a first-stage solution

### Keywords
(see function arguments)
### Returns
* `all_subpaths` - pool of sub-paths for each line, freq, and scenario
* `subpath_road_networks` - road network used for subpath routing
* `all_load_expanded_graphs` - second stage load-expanded networks
* `all_subpath_graphs` - time-expanded networks for finding subpaths
"""
function generateRestrictedSubPathSet(m::MapData,
    R::RouteSettingData,        # model input
    Inst::InstanceSettingData,  # instance data
    xVal,
    CG::Bool,                   # solve second-stage with column generation (if true, we start with an empty set of sub-paths)
    TRANSIT::Bool,
    )
    @unpack Lines, All_stops, Num_freq, Walk_graph, Demand, Taxi_travel_times = R
    @unpack K, MaxDev, MaxWalk, MaxWait, MaxSchDev, routeSpeedFactor, TDmp, TDsp, M, λ, μ, σ, δ = Inst
    
    all_subpaths = [[Vector{Pool}() for t in eachindex(Lines[l].freq)] for l in eachindex(Lines)]
    subpath_road_networks = [[Vector{MyNetwork}() for i in 1:length(Lines[l].ref_stop_id)-1] for l in eachindex(Lines)]
    generic_subpath_graphs = [[[Vector{TEN}() for i in 1:length(Lines[l].ref_stop_id)-1] for t in eachindex(Lines[l].freq)] for l in eachindex(Lines)]
    all_subpath_graphs = [[[Vector{TEN}() for i in 1:length(Lines[l].ref_stop_id)-1] for t in eachindex(Lines[l].freq)] for l in eachindex(Lines), s in eachindex(Demand)]
    all_load_expanded_graphs = [Vector{TSLGraph}() for l in eachindex(Lines)]
    
    for (l, line) in enumerate(Lines)
        @unpack all_stops, ref_stop_id, ref_stop_time, capacity, ref_route_nodes, all_route_nodes = line
        for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
            src = ref_stop_id[i]                            # start ref stop of sub-path
            idxI = findfirst(isequal(src), ref_route_nodes)
            snk = ref_stop_id[j]                            # end ref stop of sub-path
            idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
            idxJ += idxI
            subpath_network = TRANSIT ? buildLineNetwork(m, ref_route_nodes[idxI:idxJ], Set(all_route_nodes), all_stops, TDsp, routeSpeedFactor, 0.1) : buildLineNetwork(m, ref_route_nodes[idxI:idxJ], Set(all_route_nodes), all_stops, TDsp, routeSpeedFactor, MaxDev)
            push!(subpath_road_networks[l][i], subpath_network)
        end
        for t in eachindex(line.freq)
            subproblemT = ceil(Int, Dates.value(Dates.Second(ref_stop_time[ref_stop_id[end]][t] - ref_stop_time[ref_stop_id[1]][t]))/TDsp)
            load_expanded_graph = initializeTSLgraph(line, subproblemT, t, TDsp, REFSTOP_TIME_GAP)
            push!(all_load_expanded_graphs[l], load_expanded_graph)
            if xVal[l,t] > EPS
                for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
                    src = ref_stop_id[i]                            # start ref stop of sub-path
                    idxI = findfirst(isequal(src), ref_route_nodes)
                    snk = ref_stop_id[j]                            # end ref stop of sub-path
                    idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
                    idxJ += idxI
                    ti = ref_stop_time[src][t]
                    tj = ref_stop_time[snk][t]
                    subpath_graphT = ceil(Int, Dates.value(Dates.Second(tj-ti))/TDsp)
                    # build generic TEN
                    g = buildGenericTEN(subpath_road_networks[l][i][j-i], subpath_graphT)
                    DFSSort!(g) # sort the DAG in topological order
                    push!(generic_subpath_graphs[l][t][i], g)
                end
            end
        end
    end

    for (l, line) in enumerate(Lines)
        @unpack all_stops, ref_stop_id, ref_stop_time, capacity, ref_route_nodes, all_route_nodes = line
        for t in eachindex(line.freq)
            t0 = ref_stop_time[ref_stop_id[1]][t]
            subproblemT = ceil(Int, Dates.value(Dates.Second(ref_stop_time[ref_stop_id[end]][t] - t0))/TDsp)
            load_expanded_graph = all_load_expanded_graphs[l][t]
            for s in eachindex(Demand)
                # passengers to be considered in the subproblem
                # divide scenario-assigned passengers to its closest segment (a segment is the space betwen two consecutive ref stops). This is to avoid a solution picking up the same passenger multiple times.
                segment_demand = splitPassengersBySegments(m, l, line, t, Demand[s])
                pool = Pool(load_expanded_graph.V)
                count = 0
                for p in eachindex(Demand[s])
                    pool.pax[p] = Vector{Int}()
                end
                if xVal[l,t] > EPS
                    for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
                        src = ref_stop_id[i]                            # start ref stop of sub-path
                        idxI = findfirst(isequal(src), ref_route_nodes)
                        snk = ref_stop_id[j]                            # end ref stop of sub-path
                        idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
                        idxJ += idxI
                        ti = ref_stop_time[src][t]
                        tj = ref_stop_time[snk][t]
                        g = deepcopy(generic_subpath_graphs[l][t][i][j-i])
                        subpath_segment_demand = Vector{Tuple{Int,Int}}()
                        for idxR in i:j-1
                            append!(subpath_segment_demand, segment_demand[idxR])
                        end
                        g.inA, g.outA = TRANSIT ? updateGenTENwithDemandForTransit(m, R, Inst, ti, tj, subpath_segment_demand, g, subpath_road_networks[l][i][j-i], l,t,s, NORMALIZED) : updateGenTENwithDemand(m, R, Inst, ti, tj, subpath_segment_demand, g, subpath_road_networks[l][i][j-i], l,t,s, NORMALIZED)
                        push!(all_subpath_graphs[l,s][t][i], g)
                        t1 = ceil(Int, Dates.value(Dates.Second(ref_stop_time[src][t]-t0))/TDsp)
                        t2 = ceil(Int, Dates.value(Dates.Second(ref_stop_time[snk][t]-t0))/TDsp)
                        emptySubpath = false # we should always guarantee an "empty" path.
                        if !CG
                            labels, prev = fullEnumAlg(g, 1, g.order, capacity)
                            # recreate paths
                            paths, labs = getPathsAndLabels(g.V, 1, labels, prev)
                            for (idx, label) in enumerate(labs)
                                WT = sum(label.second)
                                if WT < LARGE_NUMBER
                                    Q = 0
                                    for (orig, time, num) in label.first
                                        Q += num
                                    end
                                    cost = WT 
                                    if Q == 0
                                        emptySubpath = true
                                    end
                                    for c in 0:capacity-Q
                                        count += 1
                                        node1 = load_expanded_graph.p2n[i, t1+1, c+1]
                                        node2 = load_expanded_graph.p2n[j, t2+1, c+Q+1]
                                        subP = SubPath(count, cost, Q, label.first, paths[idx], node1, node2)
                                        push!(pool.all, subP)
                                        push!(pool.in[node2], count)
                                        push!(pool.out[node1], count)
                                        for (k, h, num) in label.first
                                            push!(pool.pax[k, h], count)
                                        end
                                    end
                                end
                            end
                        end
                        if !emptySubpath
                            for c in 0:capacity
                                count += 1
                                node1 = load_expanded_graph.p2n[i, t1+1, c+1]
                                node2 = load_expanded_graph.p2n[j, t2+1, c+1]
                                subP = SubPath(count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), node1, node2)
                                push!(pool.all, subP)
                                push!(pool.in[node2], count)
                                push!(pool.out[node1], count)
                            end
                        end
                    end
                    # add set of arcs from the artificial source node from the first stop in the route (empty vehicle)
                    count += 1
                    node1 = 1
                    node2 = load_expanded_graph.p2n[1, 1, 1]
                    subP = SubPath(count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), node1, node2)
                    push!(pool.all, subP)
                    push!(pool.in[node2], count)
                    push!(pool.out[node1], count)
                    # add set of arcs to the artificial sink node from the last stop in the route
                    for c in 0:capacity
                        count += 1
                        node1 = load_expanded_graph.p2n[length(ref_stop_id), subproblemT+1, c+1]
                        node2 = load_expanded_graph.V
                        subP = SubPath(count, 0., 0, Vector{Tuple{Int, Int, Int}}(), Vector{Int}(), node1, node2)
                        push!(pool.all, subP)
                        push!(pool.in[node2], count)
                        push!(pool.out[node1], count)
                    end
                end
                push!(all_subpaths[l][t], pool)
            end
        end
    end
    println("DONE!")
    return all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs
end
