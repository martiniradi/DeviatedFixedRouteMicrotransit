"""
AKA Flat network:
    - Each node represents a physical location
    - Each arc is a road segment connecting two nodes (with a weight equal to the travel time to traverse such arc)
    - I use this flat network to define the potential paths (deviations in space) within a line.
"""
mutable struct MyNetwork
    V::Int64                                    # num nodes. E.g., number of crossings.
    src::Int64                                  # source vertex (local id). E.g., starting stop of the line.
    snk::Int64                                  # sink vertex (local id). E.g., ending stop of the line.
    outA::Vector{Vector{Tuple{Int64,Float64}}}    # array of outgoing arcs (j, tt) from a node i. Each arc is (i,j) with tarvel time tt.
    inA::Vector{Vector{Tuple{Int64,Float64}}}     # array of incoming arcs (i, tt) to a node j. Each arc is (i,j) with tarvel time tt.
    B2A::Vector{Int64}                          # translate: local id -> OSM node id (global)
    A2B::Dict{Int64,Int64}                      # translate: OSM node id -> local id
end
MyNetwork(nn::Int64) = MyNetwork(nn, 0, 0, [Vector{Tuple{Int64,Float64}}() for n = 1:nn], [Vector{Tuple{Int64,Float64}}() for n = 1:nn], Vector{Int64}(), Dict{Int64,Int64}())
"""
Load-expanded graph:
    - Each node represents a triplet (x, t, c) location, time,  num of pax + artificial source (1 and sink nodes (V)
    - Each arc is a sub-path (with a given passenger set).
    - V: number of nodes
    - n2p: from node id to triplet (x,t,c)
    - p2n: from triplet (x,t,c) to node idx
"""
struct TSLGraph
    V::Int                                          # num of nodes
    n2p::Vector{Tuple{Int,Int,Int}}           # from node id to (x,t,c)
    p2n::Array{Int,3}                          # from (x,t+1,c+1) to node id
end

"""
Time-expanded network: Graph for finding sub-paths
    - Each node represents a physical location/stop (x) and time instant (t)
    - The weight of each arc is a tuple of the passengers that can be picked up (Vector of Tuple{orig, time, number of pax}), and the associated weighted sum of walking and waiting times.
    - inA: incoming arcs to node n
    - outA: outgoing arcs from node n
    - n2p: from node id to tuple (x, t)
    - p2n: from Tuple (x,t) to node id
    - order: topological order of the nodes in the graph
"""
mutable struct TEN
    V::Int
    inA::Vector{Vector{Tuple{Int,Vector{Tuple{Int,Int,Int}},Vector{Float64}}}}
    outA::Vector{Vector{Tuple{Int,Vector{Tuple{Int,Int,Int}},Vector{Float64}}}}
    n2p::Vector{Tuple{Int,Int}}                # from node id to (x,t)
    p2n::Array{Int,2}                          # from (x,t) to node id
    order::Vector{Int}                          # Depth first sorted order of nodes. Used for the DP algorithm
end
"""
compute the shortest path using Dijkstra's algorithm
- N:: network
- src: source node
- snk: sink node
RETURN:
- path (sequence of nodes), distance
"""
function dijkstra_shortest_path(g::MyNetwork, src::Int64) #, dst::Int64
    nvg = g.V
    dists = fill(typemax(Float64), nvg)
    parents = zeros(Int64, nvg)
    visited = zeros(Bool, nvg)
    H = PriorityQueue{Int64,Float64}()

    dists[src] = zero(Float64)
    visited[src] = true
    H[src] = zero(Float64)

    while !isempty(H)
        u = dequeue!(H)

        d = dists[u] # Cannot be typemax if `u` is in the queue
        for (v, tt) in g.outA[u]
            alt = d + tt

            if !visited[v]
                visited[v] = true
                dists[v] = alt
                parents[v] = u
                H[v] = alt
            elseif alt < dists[v] - 0.01
                dists[v] = alt
                parents[v] = u
                H[v] = alt
            end
        end
    end
    return parents, dists
end
function dijkstra_shortest_path(g::MyNetwork, src::Int64, dst::Int64)
    nvg = g.V
    dists = fill(typemax(Float64), nvg)
    parents = zeros(Int64, nvg)
    visited = zeros(Bool, nvg)
    H = PriorityQueue{Int64,Float64}()

    dists[src] = zero(Float64)
    visited[src] = true
    H[src] = zero(Float64)

    while !isempty(H)
        u = dequeue!(H)

        d = dists[u] # Cannot be typemax if `u` is in the queue
        for (v, tt) in g.outA[u]
            alt = d + tt

            if !visited[v]
                visited[v] = true
                dists[v] = alt
                parents[v] = u
                H[v] = alt
            elseif alt < dists[v] - 0.01
                dists[v] = alt
                parents[v] = u
                H[v] = alt
            end
        end
    end
    node = dst
    path = Vector{Int64}()
    push!(path, node)
    while (parents[node] != 0)
        node = parents[node]
        push!(path, node)
    end
    return reverse(path), dists[dst]
end
"""
create a sparse network from an actual mapdata subnetwork, where nodes only refer to stops and not other intersections.
INPUT:
- N: The MyNetwork to spardify. Here each node corresponds to a crossing and not necessarily a stop.
- ss: Vector of stops
- TD: time discretization (in seconds)
we sparsify by computing the shortest route between any pair of adjacent stops
RETURN:
- Sparse network (MyNetwork)
"""
function fromNetworkToSparseNetwork(m::MapData, N::MyNetwork, ls::Vector{Int}, TD::Int ) #, speedFactor::Float64)
    @unpack V, outA, inA, B2A, A2B = N
    newN = MyNetwork(length(ls))
    newN.B2A = ls
    for (ii, n) in enumerate(newN.B2A)
        newN.A2B[n] = ii
        if n == B2A[N.src]
            newN.src = ii
        end
        if n == B2A[N.snk]
            newN.snk = ii
        end
    end
    # travel times are in seconds (coming from the non-sparse network)
    # udpate those time units to match the ones from the formulation TD
    for ii in 1:length(ls)
        s1 = ls[ii]
        src = A2B[s1]
        # compute shortest path on N
        parents, dists = dijkstra_shortest_path(N, src)
        for jj in 1:length(ls)
            if ii != jj
                s2 = ls[jj]
                dst = A2B[s2]
                node = dst
                path = Vector{Int64}()
                push!(path, node)
                while (parents[node] != 0)
                    node = parents[node]
                    push!(path, node)
                end
                reverse!(path)
                tt = dists[dst]
                if tt < 9999
                    # if none of the nodes in path is a stop --> add arc
                    osmPath = B2A[path]
                    if length(intersect(osmPath[2:end-1], ls)) < 0.01
                        newTT = round(tt / TD) # keep unrounded
                        if newTT < 0.01
                            newTT = 1.
                        end
                        # add arc
                        ii_ = newN.A2B[s1]
                        jj_ = newN.A2B[s2]
                        push!(newN.inA[jj_], (ii_, newTT))
                        push!(newN.outA[ii_], (jj_, newTT))
                    end
                end
            end
        end
    end
    return newN
end
"""
Build the MyNetwork object for a given route segment of a line.

### Keywords
* `m` - MapData
* `r` - OSM nodes forming the route segment
* `nn` - all OSM nodes forming the route segment and the nodes within reach (given the maximum deviation)
* `ss` - all stops within reach (given the maximum deviation)
* `TD` - time discretization (in seconds)
* `speedFactor` - how slower the vehicles goes compared to the travel times given by `get_speeds()`
* `MaxDev` - maximum deviation from reference route (in meters from any crossing in `r`)
### Returns
* MyNetwork object (already sparse)
"""
function buildLineNetwork(m::MapData, r::Vector{Int}, nn::Set{Int}, ss::Vector{Int}, TD::Int64, speedFactor::Float64, MaxDev::Float64)
    B2A = Vector{Int}() # surrounding nodes
    for n in r
        ind = getLocRangeNodes(m, m.nodes[n], MaxDev, nn)
        for ii in ind
            push!(B2A, ii)
        end
    end
    sort!(B2A)
    unique!(B2A)
    subN = MyNetwork(length(B2A))
    subN.B2A = B2A
    ls = Vector{Int}() # surrounding stops
    tau = OpenStreetMapX.create_weights_matrix(m, OpenStreetMapX.network_travel_times(m, SPEEDS))
    for (ii, n) in enumerate(subN.B2A)
        if n in ss
            push!(ls, n)
        end
        subN.A2B[n] = ii
        if n == r[1]
            subN.src = ii
        end
        if n == r[end]
            subN.snk = ii
        end
    end
    for (ii, n) in enumerate(subN.B2A)
        for on in outneighbors(m.g, m.v[n])
            n_ = m.n[on]
            if haskey(subN.A2B, n_)
                jj = subN.A2B[n_]
                tt = speedFactor * tau[m.v[n], on]
                push!(subN.inA[jj], (ii, tt))
                push!(subN.outA[ii], (jj, tt))
            end
        end
    end
    subNsparse = fromNetworkToSparseNetwork(m, subN, ls, TD)
    return subNsparse
end
"""
Initialize the node set of the TSLN for a given line and timeslot

### Keywords
* `line` - Line object
* `T` - maximum time instant (given by the arrival to end of the line)
* `timeslot` - the timeslot id
* `gap` - (NOT USED) time deviation allowed at each ref stop visited
### Returns
* TSLN object (no arcs)
"""
function initializeTSLgraph(line::Line, T::Int, timeslot::Int, TDsp, gap::Int)
    @unpack ref_stop_id, ref_stop_time, capacity = line
    n2p = Vector{Tuple{Int,Int,Int}}()
    push!(n2p, (0, 0, 0))
    p2n = zeros(Int, length(ref_stop_id), T + 1, capacity + 1)
    count = 1
    t0 = ref_stop_time[ref_stop_id[1]][timeslot]
    for (i, n) in enumerate(ref_stop_id)
        time = ceil(Int, Dates.value(Dates.Second(ref_stop_time[n][timeslot] - t0))/TDsp)
        for t in max(0, time-gap):min(time+gap, T)
            for c in 0:capacity
                count += 1
                p2n[i, t+1, c+1] = count
                push!(n2p, (i, t, c))
            end
        end
    end
    count += 1
    push!(n2p, (0, 0, 0))
    TSLN = TSLGraph(count, n2p, p2n)
    return TSLN
end
"""
build a time-expanded network (TEN object) for a given route segment of a line, and a frequency

### Keywords
* `N` - route-network corresponding to the segment area
* `T` - discretized time-interval between stopping times of start and end reference stop in the segment
### Returns
* time-expanded network without any demand associated
"""
function buildGenericTEN(N::MyNetwork, T::Int)
    """
    build a time expanded network, for the entire line + max deviation in time and space
    #     N::MyNetwork # subset of flat network, only with locations valid for the TEN (arcs constrained to ensure no cycles)
    #     V::Int
    #     inA - incoming arcs to each node
    #     outA - outgoing arcs from each node
    #     n2p::Vector{Tuple{Int, Int}}        # from node id to (x,t)
    #     p2n::Array{Int, 2}                  # from (x,t) to node id
    """
    X = N.V
    p2n = zeros(Int64, X, T + 1)
    n2p = Tuple{Int64,Int64}[]
    count = 1 # artificial source node
    push!(n2p, (0, 0))
    for (n, src) in enumerate(N.B2A)
        for t in 0:T
            count += 1
            p2n[n, t+1] = count
            push!(n2p, (n, t))
        end
    end
    count += 1
    push!(n2p, (0,0)) # sink node
    inA = [Vector{Tuple{Int,Vector{Tuple{Int,Int, Int}},Vector{Float64}}}() for n in 1:count]
    outA = [Vector{Tuple{Int,Vector{Tuple{Int,Int, Int}},Vector{Float64}}}() for n in 1:count]
    # add artificial source and sink arcs
    src = N.src
    node = 1
    node2 = p2n[src, 1]
    push!(inA[node2], (node, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
    push!(outA[node], (node2, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
    earliestEndTime = T #max(0, T-gap)
    for t in earliestEndTime:T
        snk = N.snk
        node = p2n[snk, t+1]
        node2 = count
        push!(inA[node2], (node, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
        push!(outA[node], (node2, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
    end
    for (src, n) in enumerate(N.B2A)
        for (dst, tt) in N.outA[src]
            for t in 0:T-Int(tt)
                node = p2n[src, t+1]
                t2 = t + Int(tt) #+ rand(0:2)
                if t2 <= T
                    node2 = p2n[dst, t2+1]
                    push!(inA[node2], (node, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
                    push!(outA[node], (node2, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
                end
            end
        end
        for t in 0:T-1
            node = p2n[src, t+1]
            node2 = p2n[src, t+2]
            push!(inA[node2], (node, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
            push!(outA[node], (node2, Vector{Tuple{Int,Int,Int}}(), Vector{Float64}()))
        end
    end
    return TEN(count, inA, outA, n2p, p2n, collect(1:count))
end
"""
run a Depth-first search on a TEN (DAG) to sort the nodes in topological order
"""
function DFSSort!(g)
    vcolor = zeros(UInt8, g.V)
    verts = Vector{Int}()
    for v in 1:g.V
        vcolor[v] != 0 && continue
        S = Vector{Int}([v])
        vcolor[v] = 1
        while !isempty(S)
            u = S[end]
            w = 0
            for (n, cus, t) in g.outA[u]
                if vcolor[n] == 1
                    error("The input graph contains at least one loop.")
                elseif vcolor[n] == 0
                    w = n
                    break
                end
            end
            if w != 0
                vcolor[w] = 1
                push!(S, w)
            else
                vcolor[u] = 2
                push!(verts, u)
                pop!(S)
            end
        end
    end
    g.order = reverse(verts)
end
"""
Associate each AggregPassenger index (loc, time) to a segment between consecutive reference stop in a line

### Keywords
* `m` - MapData
* `l` - line id
* `line` - Line object
* `freq` - the freq id
* `D` - Dictionary: (orig, time) => AggregPassenger object
### Returns
* Vector of AggregPassenger indices divided by segments between consecutive reference stops in a line
"""
function splitPassengersBySegments(m::MapData, l::Int, line::Line, freq::Int, D::Dict{Tuple{Int64,Int64},AggregPassenger})
    # Assign each passenger to one set of sub-paths (that do not allow counting same passenger twice)
    @unpack ref_stop_id, ref_stop_time = line
    segment_demand = [Vector{Tuple{Int,Int}}() for i in 1:length(ref_stop_id)-1]
    for (id, p)  in collect(D)
        if (l, freq) in p.candidateTrips
            # find best segment
            stop_idx = nearest_node_index(m, p.origin_stop, ref_stop_id)
            if (stop_idx != 1 && stop_idx != length(ref_stop_id) && OpenStreetMapX.distance(m.nodes[p.origin_stop], m.nodes[ref_stop_id[stop_idx-1]]) < OpenStreetMapX.distance(m.nodes[p.origin_stop], m.nodes[ref_stop_id[stop_idx+1]])) || stop_idx == length(ref_stop_id)
                stop_idx -= 1
            end
            push!(segment_demand[stop_idx], id)
        end
    end
    return segment_demand
end

"""
update a time-expanded network with the potential demand of a given scenario.
We do this by updating the labels of the arcs

### Keywords
* `m` - MapData
* `R` - RouteSettingData
* `ti` - stopping time at the start reference stop of the time-expanded network
* `tj` - stopping time at the end stop of the time-expanded network
* `D` - vector of tuples (indices of AggregPassenger objects in R.Demand[s])
* `g` - time-expanded network
* `N` - road-network corresponding to the subpath routing area
* `l` - line id
* `t` - freq id
* `s` - scenario id
### Returns
* updated arc set of the time-expanded network with demand
"""
function updateGenTENwithDemand(m::MapData,
                                R::RouteSettingData,
                                Inst::InstanceSettingData,
                                ti::Time,
                                tj::Time,
                                D::Vector{Tuple{Int,Int}},
                                g::TEN,
                                N::MyNetwork,
                                l::Int,
                                t::Int,
                                s::Int,
                                normalized::Bool,
                                )
    @unpack Taxi_travel_times, Demand, Walk_graph = R
    @unpack MaxWait, M, λ, μ, σ, MaxWalk, TDmp, TDsp = Inst
    @unpack inA, outA, n2p, p2n = g
    for id in D
        (origin_stop, time_MP) = id
        p = Demand[s][id]
        @unpack num, request_dropoff_time = p
        be_ready_time = p.be_ready_times[(l,t)]
        trip_delay = p.trip_delays[(l,t)] # already weighted and normalized
        walk_stops = getLocRangeNodes(m, m.nodes[origin_stop], MaxWalk, Set(N.B2A))
        # compute shortest paths from passenger origin
        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(Walk_graph, m.v[origin_stop])
        for cand_stop in walk_stops
            # compute walking distance
            walkPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[cand_stop])
            walk_dist = 0.
            for (idx, node) in enumerate(walkPath[1:end-1])
                node2 = walkPath[idx+1]
                walk_dist += Walk_graph.weights[node, node2]
            end
            if walk_dist <= MaxWalk
                walk_time_seconds = ceil(Int, walk_dist/WALKSPEED_M_PER_SEC)
                # walk_time_minutes = walk_time_seconds/60
                src = N.A2B[cand_stop]
                earliest_pickup_SP = max(ceil(Int,  Dates.value(Dates.Second(be_ready_time + Dates.Second(walk_time_seconds)-ti))/TDsp),0)
                latest_pickup_SP = min(floor(Int, Dates.value(Dates.Second(be_ready_time + Dates.Second(walk_time_seconds + MaxWait)-ti))/TDsp), floor(Int, Dates.value(Dates.Second(tj-ti))/TDsp))
                if earliest_pickup_SP <= latest_pickup_SP
                    for pickup_time_SP in earliest_pickup_SP:latest_pickup_SP
                        node = p2n[src, pickup_time_SP+1]
                        pickup_time = ti + Dates.Second(pickup_time_SP*TDsp)
                        in_vehicle_time = Dates.value(Dates.Second(R.Lines[l].arrival_time[t] - pickup_time))
                        if normalized
                            in_vehicle_time /= Taxi_travel_times[origin_stop]
                        end
                        wait_time_seconds = Dates.value(Dates.Second(pickup_time - (be_ready_time + Dates.Second(walk_time_seconds))))
                        for (i, v) in enumerate(outA[node])
                            (node2, C, WT) = v
                            (xC, tC) = node2 != g.V ? n2p[node2] : (0,0)
                            if xC != src # not a holding arc
                                push!(C, (origin_stop, time_MP, num))
                                passengerService = round(num*(λ*walk_time_seconds + μ*wait_time_seconds + σ*in_vehicle_time + trip_delay - M), digits=6)
                                push!(WT, passengerService)
                                outA[node][i] = (node2, C, WT)
                                for (j, v2) in enumerate(inA[node2])
                                    (node_, C_, WT_) = v2
                                    (src_, t2_) = n2p[node_]
                                    if src_ == src
                                        push!(C_, (origin_stop, time_MP, num))
                                        push!(WT_, passengerService)
                                        inA[node2][j] = (node_, C_, WT_)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return inA, outA
end
"""
equivalent of updateGenTENwithDemand() for the fixed-transit version,
demand can only be assigned to nodes/arcs related to transit stops and scheduled times

### Keywords
* `m` - MapData
* `R` - RouteSettingData
* `ti` - stopping time at the start reference stop of the time-expanded network
* `tj` - stopping time at the end stop of the time-expanded network
* `D` - vector of tuples (indices of AggregPassenger objects in R.Demand[s])
* `g` - time-expanded network
* `N` - road-network corresponding to the subpath routing area
* `l` - line id
* `t` - freq id
* `s` - scenario id
### Returns
* updated arc set of the time-expanded network with demand
"""
function updateGenTENwithDemandForTransit(m::MapData, 
    R::RouteSettingData,
    Inst::InstanceSettingData,
    ti::Time,
    tj::Time,
    D::Vector{Tuple{Int,Int}},
    g::TEN,
    N::MyNetwork,
    l::Int,
    t::Int,
    s::Int,
    normalized::Bool,
    )
    @unpack Taxi_travel_times, Demand, Walk_graph = R
    @unpack MaxWait, M, λ, μ, σ, MaxWalk, TDmp, TDsp = Inst
    @unpack inA, outA, n2p, p2n = g

    for id in D
        (origin_stop, time_MP) = id
        p = Demand[s][id]
        @unpack num, request_dropoff_time = p
        be_ready_time = p.be_ready_times[(l,t)]
        trip_delay = p.trip_delays[(l,t)]
        walk_stops = getLocRangeNodes(m, m.nodes[origin_stop], MaxWalk, Set(R.Lines[l].transit_stop_id))
        # compute shortest paths from passenger origin
        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(Walk_graph, m.v[origin_stop])
            for cand_stop in walk_stops
                if haskey(N.A2B, cand_stop)
                    # compute walking distance
                    walkPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[cand_stop])
                    walk_dist = 0.
                    for (idx, node) in enumerate(walkPath[1:end-1])
                        node2 = walkPath[idx+1]
                        walk_dist += Walk_graph.weights[node, node2]
                    end
                    if walk_dist <= MaxWalk
                        walk_time_seconds = ceil(Int, walk_dist/WALKSPEED_M_PER_SEC )
                        src = N.A2B[cand_stop]
                        earliest_pickup = be_ready_time + Dates.Second(walk_time_seconds)
                        latest_pickup = min(be_ready_time + Dates.Second(walk_time_seconds + MaxWait), tj)
                        # check if it can be picked at the scheduled time
                        if earliest_pickup <= R.Lines[l].transit_stop_time[cand_stop][t] <= latest_pickup
                            pickup_time_SP = ceil(Int,  Dates.value(Dates.Second(R.Lines[l].transit_stop_time[cand_stop][t]-ti))/TDsp)
                            node = p2n[src, pickup_time_SP+1]
                            pickup_time = ti + Dates.Second(pickup_time_SP*TDsp)
                            in_vehicle_time = Dates.value(Dates.Second(R.Lines[l].arrival_time[t] - pickup_time))
                            if normalized
                                in_vehicle_time /= Taxi_travel_times[cand_stop]
                            end
                            wait_time = Dates.value(Dates.Second(pickup_time - (be_ready_time + Dates.Second(walk_time_seconds))))
                            # waiting time MUST be below the maximum allowed
                            for (i, v) in enumerate(outA[node])
                            (node2, C, WT) = v
                            (xC, tC) = node2 != g.V ? n2p[node2] : (0,0)
                            if xC != src # not a holding arc
                                push!(C, (origin_stop, time_MP, num))
                                passengerService = num*(λ*walk_time_seconds + μ*wait_time + σ*in_vehicle_time + trip_delay - M)
                                push!(WT, passengerService)
                                outA[node][i] = (node2, C, WT)
                                for (j, v2) in enumerate(inA[node2])
                                    (node_, C_, WT_) = v2
                                    (src_, t2_) = n2p[node_]
                                    if src_ == src
                                        push!(C_, (origin_stop, time_MP, num))
                                        push!(WT_, passengerService)
                                        inA[node2][j] = (node_, C_, WT_)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return inA, outA
end