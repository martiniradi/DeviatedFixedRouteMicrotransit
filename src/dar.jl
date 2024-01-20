"""
This files includes the main functions to run the dial-a-ride version of the Microtransit Network Design (MiND-DAR) problem.

"""


"""
Information about a specific location (Mainly pickup adn dropoff info)
    - coordinates
    - closest OSM node
    - closest stop
    - time
"""
struct Loc
    coord::Tuple{Float64, Float64}
    node::Int64     # closest OSM node to loc
    stop::Int64     # closest OSM node that is a stop
    time::Time      # time
end
"""
Information about a single passenger in the MiND-DAR case:
    - Index
    - Origin coordinates (lat, lon)
    - Destination coordinates (lat, lon)
    - (orig, dest) Node and stop associated to passenger
    - Taxi pickup and dropoff time
"""
struct PassengerDAR
    id::Int         # idx
    orig::Loc # origin info (taxi pickup time)
    dest::Loc # destination info (requested dropoff time)
end

function readStopsDAR(m::MapData)
    data_dir = get_data_path()
    df = CSV.read(joinpath(data_dir,"input", "MiND-DAR_setup", "Stops.csv"), DataFrame)
    stops = Vector{Int}()
    for r in eachrow(df)
        n = point_to_nodes((r.Lat, r.Long), m)
        push!(stops, n)
    end
    println("Number of stops: ",length(stops))
    return unique!(stops)
end

function readLinesDAR(m::MapData, selectedLines::Vector{Int}, operating_horizon_seconds::Int, numL::Int, refNum::Int, Capacity::Int, speedFactor::Float64, MaxDev::Float64, exitPoints::Vector{Int}, exitTimesSeconds::Vector{Float64}, frequency_period_seconds::Int, maximum_delay::Int, start_operations::Dates.Time)
    dwell_time = 30 # stop duration in transit stops ( in seconds)

    T = ceil(Int, operating_horizon_seconds/frequency_period_seconds) # number of frequencies
    horizon_end = start_operations + Dates.Second(operating_horizon_seconds + maximum_delay)
    data_dir = get_data_path()

    df = CSV.read(joinpath(data_dir,"input", "MiND-DAR_setup", "Lines.csv"), DataFrame)

    S = readStopsDAR(m) # generate all candidate stops

    lines = Vector{Line}()

    i = 0       # line index counter
    line_ref_stops = [Vector{Int}() for i in 1:numL]
    line_all_stops = [Vector{Int}() for i in 1:numL]
    for r in eachrow(df)
        if r.LineID in selectedLines
            if r.SeqId == 1
                i += 1
            end
            n = point_to_nodes((r.Lat, r.Long), m)
            if !(n in S)
                println("Stop ", n, " not in preset ", r.LineID, ", ", r.SeqId)
                nearest_n = nearest_node(m, r.Lat, r.Long, S)
                push!(line_all_stops[i], nearest_n)
                if (r.SeqId - 1) % refNum == 0
                    push!(line_ref_stops[i], nearest_n)
                end
            else
                push!(line_all_stops[i], n)
                if (r.SeqId - 1) % refNum == 0
                    push!(line_ref_stops[i], n)
                end
            end
        end
    end
    println("TOTAL NUMBER OF STOPS: ", length(S))
    print("Creating lines...")
    for idx in eachindex(selectedLines)
        ref_stop_id = line_ref_stops[idx]
        ref_stop_time = Dict{Int,Vector{Time}}()
        transit_stop_time = Dict{Int,Vector{Time}}()
        ref_route_nodes = Vector{Int}([ref_stop_id[1]])
        ref_times = Vector{Float64}([0.])
        visit_time = 0.
        ref_dist = 0.
        for i in 1:length(line_all_stops[idx])-1
            n1 = line_all_stops[idx][i]
            n2 = line_all_stops[idx][i+1]
            route, rDist, rTime = fastest_route(m, n1, n2; speeds=SPEEDS)
            ref_dist += rDist
            append!(ref_route_nodes, route[2:end])
            # check how many stops are in between reference stops (intermediate_stops)
            intermediate_stops = length(intersect(route[1:end-1], S))
            visit_time += speedFactor*rTime + intermediate_stops*dwell_time
            push!(ref_times, visit_time)

        end

        line_duration_seconds = ceil(Int, ref_times[end])
        # find line frequencies to run
        line_frequencies = collect(start_operations:Dates.Second(frequency_period_seconds):horizon_end-Dates.Second(line_duration_seconds))
        arrival_times = [start_operations + Dates.Second(line_duration_seconds + (t-1)*frequency_period_seconds) for t in eachindex(line_frequencies)]

        # define schedules for all the frequencies
        ref_stop_time[ref_stop_id[1]] = [ start_operations + (t-1)*Dates.Second(frequency_period_seconds) for t in eachindex(line_frequencies)]
        transit_stop_time[ref_stop_id[1]] = [ start_operations + (t-1)*Dates.Second(frequency_period_seconds) for t in eachindex(line_frequencies)]
        for i in 1:length(line_all_stops[idx])-1
            n1 = line_all_stops[idx][i]
            n2 = line_all_stops[idx][i+1]
            if n2 in ref_stop_id
                ref_stop_time[n2] = [start_operations + Dates.Second(ceil(Int, ref_times[i+1]) + (t-1)*frequency_period_seconds) for t in eachindex(line_frequencies)]
            end
            transit_stop_time[n2] = [start_operations + Dates.Second(ceil(Int, ref_times[i+1]) + (t-1)*frequency_period_seconds) for t in eachindex(line_frequencies)]
        end

        # surrounding nodes
        all_route_nodes = Vector{Int}()
        for n in ref_route_nodes
            ind = getLocRangeNodes(m, m.nodes[n], MaxDev, Set(m.n))
            for i_ in ind
                push!(all_route_nodes, i_)
            end
        end
        sort!(all_route_nodes)
        unique!(all_route_nodes)
        #surrounding stops
        all_stops = Vector{Int}()
        for i in all_route_nodes
            if i in S
                push!(all_stops, i)
            end
        end
        # the cost is proportional to the length (in seconds) of the line
        cost = ref_times[end]- ref_times[1]

        push!(lines, Line(all_stops, line_frequencies, ref_stop_id, ref_stop_time, arrival_times, Capacity, cost, ref_route_nodes, all_route_nodes, line_all_stops[idx], transit_stop_time, ref_dist, 0))
    end
    println(" DONE!")
    return lines, S, T
end

"""
read CSV file about a demand scenario and creates the vector of Request objects

### Keywords
* `m` - MapData
* `demand_horizon_seconds` - demand horizon (from `start_time`) in seconds
* `all_stops` - Vector of stops (OSM node ids)
* `scenario` - scenario idx
* `trainSet` - training (0) or test (1)
* `start_time` - starting time of the day
### Returns
* the vector of Passenger objects
"""
function readDemandDAR(
    m::MapData,
    demand_horizon_seconds::Int,
    all_stops::Vector{Int},
    scenario::Int,
    trainSet::Int,
    taxi_travel_times::Dict{Tuple{Int,Int}, Float64},
    start_time::Dates.Time
    )
    df = CSV.read(joinpath(get_data_path(),  "input", "MiND-DAR_setup", "weekday_data_midtown_DAR", "taxiDataDAR_manhattan_midtown_$(scenario)_$(trainSet).csv"), DataFrame)

    # build list of Passengers from DataFrame
    D = Vector{PassengerDAR}()
    for row in eachrow(df)
        dropoffTime = Dates.Time(round(Dates.DateTime(row.dtime, "yyyy-mm-dd HH:MM:SS"), Dates.Minute))
        # accept if pickup and drop-off times are within the planning horizon
        if dropoffTime <= start_time + Dates.Second(demand_horizon_seconds) 
            p_node = row.p_node
            p_stop = nearest_node(m, row.plat, row.plong, all_stops)
            d_node = row.d_node
            d_stop = nearest_node(m, row.dlat, row.dlong, all_stops)
            if haskey(taxi_travel_times, (p_stop, d_stop)) && p_stop != d_stop
                pickupTime = dropoffTime - Dates.Second(round(Int,taxi_travel_times[p_stop, d_stop]))
                push!(D, PassengerDAR(row.id, Loc((row.plat, row.plong), p_node, p_stop, pickupTime), Loc((row.dlat, row.dlong), d_node, d_stop, dropoffTime)))
            end
        end
    end
    return D
end

"""
Information about a single origin and time instant:
    * Id
    * Number of passengers
    * Origin stop
    * Requested dropoff time
    * Be-ready time depending on trip
    * Time id for first stage clustering
    * Set of potential trips that can serve this passenger
"""
mutable struct AggregPassengerDAR
    id::Int
    pax::Vector{Int}                                # set of clustered passenger ids
    num::Int                                        # number of passengers
    orig::Loc # origin info
    dest::Loc # destination info
    # taxi_trip_time_seconds::Int                      # taxi trip duration in seconds
    be_ready_times::Dict{Tuple{Int,Int}, Time}  # be ready times for each possible trip (line, frequency t)
    # trip_delays::Dict{Tuple{Int,Int}, Float64}  # weighted arrival delay by trip
    first_stage_time_instant::Int               # time instant used in master/first-stage problem
    candidateTrips::Vector{Tuple{Int,Int}}      # set of candidate trips (line, frequency) that can serve this passenger
end
"""
aggregate demand based on the first-stage time discretization and create AggregPassenger object

### Keywords
* `D` - set of passenger objects
* `TDmp` - first-stage time discretization (in seconds)
* `set` - training (0) or test (1)
* `start_time` - starting time of the day
### Returns
* the dictionary of AggregPassenger objects
"""
function aggregateDemandDAR(D::Vector{PassengerDAR}, TDmp::Int, DAY_START::Time)
    AggregDem = Dict{Tuple{Int, Int, Int}, AggregPassengerDAR}()
    counter = 0
    for p in D
        @unpack id, orig, dest = p
        first_stage_dest_time_instant = ceil(Int,Dates.value(Dates.Second(dest.time - DAY_START))/TDmp)
        if haskey(AggregDem, (orig.stop, first_stage_dest_time_instant, dest.stop))
            AggregDem[(orig.stop, first_stage_dest_time_instant, dest.stop)].num += 1
            push!(AggregDem[(orig.stop, first_stage_dest_time_instant, dest.stop)].pax, p.id)
        else
            counter += 1
            AggregDem[(orig.stop, first_stage_dest_time_instant, dest.stop)] = AggregPassengerDAR(counter, Int[id], 1, orig, dest, Dict{Tuple{Int,Int},Time}(), first_stage_dest_time_instant, Vector{Tuple{Int, Int}}())
        end
    end
    return AggregDem 
end
"""
find candidate trips (line, frequency) for each AggregPassengerDAR
"""
function findCandidateTripsDAR!(m::MapData, D::Dict{Tuple{Int,Int,Int}, AggregPassengerDAR}, L::Vector{Line}, MaxArrivalScheduleDev::Int, MaxWalkDist::Float64)
    # loop through demand
    for (key, d) in collect(D)
        # loop through lines
        for (l, line) in enumerate(L)
            # first check if the line is a candidate
            potential_line = false 
            if d.orig.stop in line.all_stops && d.dest.stop in line.all_stops
                potential_line = true
            else
                range_nodes_orig = getLocRangeNodes(m, m.nodes[d.orig.stop], MaxWalkDist, Set(line.all_stops))
                range_nodes_dest = getLocRangeNodes(m, m.nodes[d.dest.stop], MaxWalkDist, Set(line.all_stops))
                if isempty(range_nodes_orig) == false && isempty(range_nodes_dest) == false
                    potential_line = true
                end
            end
            if potential_line
                nearest_ref_stop_dest, walk_dist = nearest_node_and_dist(m, m.nodes[d.dest.stop], line.ref_stop_id)
                ref_stop_idx_dest = findfirst(isequal(nearest_ref_stop_dest), line.ref_stop_id)
                # find the right frequency
                for freq in eachindex(line.arrival_time)
                    ref_stop_time_dest = line.ref_stop_time[nearest_ref_stop_dest][freq]
                    
                    # estimated arrival time as ref stop time + walking time to destination
                    estimated_arrival_time = ref_stop_time_dest + Dates.Second(ceil(Int, walk_dist/1.4))
                    # if this time is okay then consider
                    if estimated_arrival_time - Dates.Second(MaxArrivalScheduleDev) <= d.dest.time <= estimated_arrival_time + Dates.Second(MaxArrivalScheduleDev)
                        #   find nearest reference stop
                        nearest_ref_stop_orig = nearest_node(m, m.nodes[d.orig.stop], line.ref_stop_id)
                        ref_stop_idx_orig = findfirst(isequal(nearest_ref_stop_orig), line.ref_stop_id)
                        # we need to ensure that the nearest ref stop dest is after the nearest ref stop orig
                        if ref_stop_idx_orig < ref_stop_idx_dest
                            ref_stop_time_orig = line.ref_stop_time[nearest_ref_stop_orig][freq]
                            # compute as well the be-ready time
                            # be ready time as the ref stop time - maximum walking time
                            be_ready_time = ref_stop_time_orig - Dates.Second(ceil(Int, MaxWalkDist/1.4))
                            d.be_ready_times[(l, freq)] = be_ready_time
                            push!(d.candidateTrips, (l, freq))
                        end
                    end
                end
            end
        end
    end
end
"""
compute overlapping frequencies for fleet size on each line in the MiND-DAR. The line ends at the last checkpoint, not in LaGuardia
"""
function computeOverlapFrequenciesDAR(lines::Vector{Line}, taxi_travel_times::Dict{Tuple{Int,Int}, Float64}, demand_horizon::Int, Num_freq::Int, maximum_delay::Int)
    start_operations = DAY_START - Dates.Second(OPERATIONAL_HORIZON)
    finish_operations = DAY_START + Dates.Second(demand_horizon + maximum_delay)
    # define set of frequencies
    set_frequencies = collect(start_operations:Dates.Second(TIME_PERIOD_SECONDS):finish_operations)
    freq_overlap = [Vector{Int}() for l in eachindex(lines), t in 1:Num_freq]
    for (l, line) in enumerate(lines)
        # calculate total line duration + travel back to start (the way back is taken as taxi_travel_time)
        for (t,line_start) in enumerate(line.freq)
            line_end = line.ref_stop_time[line.ref_stop_id[end]][t]
            back_at_start = line_end + Dates.Second(ceil(Int, taxi_travel_times[line.ref_stop_id[end], line.ref_stop_id[1]]))
            for (i, freq) in enumerate(set_frequencies)
                if line_start <= freq < back_at_start && freq in line.freq
                    push!(freq_overlap[l,t], i)
                end
            end
        end
    end
    return freq_overlap
end
"""
read `TravelTimesDAR.txt` file and return dictionary: (stop id, stop_id) -> travel time in seconds
"""
function readTravelTimesDAR()
    A = Dict{Tuple{Int,Int}, Float64}() # (OSM node id, OSM node id) -> time in seconds
    data_dir = get_data_path()
    f = open(joinpath(data_dir, "input", "MiND-DAR_setup", "TravelTimesDAR.txt"),"r")
    for l in readlines(f)
        vals = parse.(Float64, split(l,"\t"))
        A[Int(vals[1]), Int(vals[2])] = vals[3]
    end
    return A
end
"""
Model inputs

### Attributes
* `Lines` - set of Line objects
* `All_stops` - set of all potential stops
* `Num_freq` - number of frequencies in planning horizon
* `Freq_Overlap` - set of trip frequencies that overlap with a given line and time slot
* `Demand` - dictionary of demand (AggregPassengerDAR), one per scenario
* `Single_requests` - vector of single passenger requests (before aggregating)
* `Trip_passengers` - Set of passengers ids that can be served by a trip (l,t)
* `Pi` - Dict (scenario => probability of scenario)
* `Walk_graph` - undirected graph version from the one in m (MapData), m.g
* `Taxi_travel_times` - Direct travel time between any two stops (in seconds)
"""
struct RouteSettingDataDAR
    Lines::Vector{Line}
    All_stops::Vector{Int64}
    Num_freq::Int64
    Freq_Overlap::Array{Vector{Int},2}
    Demand::Vector{Dict{Tuple{Int, Int, Int},AggregPassengerDAR}}
    # Single_requests::Vector{Passenger}
    Trip_passengers::Array{Vector{Tuple{Int, Int, Int}},3}
    Pi::Vector{Float64}
    Walk_graph::SimpleWeightedGraph
    Taxi_travel_times::Dict{Tuple{Int,Int}, Float64}
end

"""
create the RouteSettingData object with the main input for the algorithm

### Keywords
* `Inst` - InstanceSettingData
### Returns
* the RouteSettingDataDAR object
"""
function buildRouteSettingDataDAR(Inst::InstanceSettingData)
    @unpack numL, demT, numS, trainData, scenarios, opT, K, RefStopDens, TDmp, TDsp, MaxDev, Capacity, MaxWalk, MaxWait, MaxSchDev, λ, μ, δ, routeSpeedFactor, refSpeedFactor = Inst
    m = generateNetwork()
    walk_graph = generateUndirectedGraph(m) 
    lines, all_stops, num_freq = readLinesDAR(m, collect(1:numL), demT+opT, numL, RefStopDens, Capacity, refSpeedFactor, MaxDev, EXIT_POINTS_OSM, EXIT_TIMES_SECONDS, TIME_PERIOD_SECONDS,MaxSchDev, DAY_START - Dates.Second(opT))
    taxi_travel_times = readTravelTimesDAR()

    all_demand = Vector{Dict{Tuple{Int,Int, Int},AggregPassengerDAR}}()
    trip_passengers = [Vector{Tuple{Int, Int, Int}}() for l in 1:numL, t in 1:num_freq, s in 1:numS]
    data_set = trainData ? 0 : 1 

    for (s, scenario) in enumerate(scenarios)
        passenger_demand = readDemandDAR(m, demT, all_stops, scenario, data_set, taxi_travel_times, DAY_START)
        aggregated_demand = aggregateDemandDAR(passenger_demand,TDmp, DAY_START )
        findCandidateTripsDAR!(m, aggregated_demand, lines, MaxSchDev, MaxWalk)
        push!(all_demand, aggregated_demand)
        for p in eachindex(aggregated_demand)
            for (l,t) in aggregated_demand[p].candidateTrips
                push!(trip_passengers[l,t,s], p)
            end
        end
    end
    # scenario probabilities
    prob = 1/numS
    Pi = fill(prob, numS)

    freq_overlap = computeOverlapFrequenciesDAR(lines, taxi_travel_times, demT, num_freq, MaxSchDev)

    return RouteSettingDataDAR(lines, all_stops, num_freq, freq_overlap, all_demand, trip_passengers, Pi, walk_graph, taxi_travel_times), m
end

"""
Sub-path object.
This object contains the main information of a sub-path.
"""
struct SubPathDAR
    id::Int                     # identifier
    cost::Float64               # cost (g_a)
    load::Int64                 # num of passengers picked up
    paxPU::Vector{Tuple{Int,Int,Int,Int}}            # vector of AggregPassenger ids picked up + number
    paxDO::Vector{Tuple{Int,Int,Int,Int}}            # vector of AggregPassenger ids dropped off + number (negative)
    graph_path::Vector{Int}       # path in the time-expanded network
    i::Int                      # tail node in load-expanded network
    j::Int                      # head node in load-expanded network
end
"""
Pool to define all the arcs (i.e., sub-path), of a given sub-problem (TSLN)
"""
mutable struct PoolDAR
    all::Vector{SubPathDAR}                        # set of all sub-paths
    in::Vector{Vector{Int}}                     # set of sub-path ids incoming to node i
    out::Vector{Vector{Int}}                    # set of sub-path ids outgoing from node i
    paxPU::Dict{Tuple{Int,Int,Int},Vector{Int}}       # set of sub-path ids picking up pax at (k,h)
    paxDO::Dict{Tuple{Int,Int,Int},Vector{Int}}       # set of sub-path ids dropping off pax at (k,h)
end
PoolDAR(nn::Int) = PoolDAR(Vector{SubPathDAR}(), [Vector{Int}() for n in 1:nn], [Vector{Int}() for n in 1:nn], Dict{Tuple{Int,Int,Int},Vector{Int}}(), Dict{Tuple{Int,Int,Int},Vector{Int}}())
"""
Time-expanded network for MiND-DAR: Graph for finding sub-paths
    - Each node represents a physical location/stop (x) and time instant (t)
    - The weight of each arc is a tuple of the passengers that can be picked up (Vector of Tuple{orig, time, number of pax}), and the associated weighted sum of walking and waiting times.
    - inA: incoming arcs to node n
    - outA: outgoing arcs from node n
    - n2p: from node id to tuple (x, t)
    - p2n: from Tuple (x,t) to node id
    - order: topological order of the nodes in the graph
"""
mutable struct TEN_DAR
    V::Int
    inA::Vector{Vector{Tuple{Int,Vector{Tuple{Int,Int,Int,Int}},Vector{Float64}}}}
    outA::Vector{Vector{Tuple{Int,Vector{Tuple{Int,Int,Int,Int}},Vector{Float64}}}}
    n2p::Vector{Tuple{Int,Int}}                # from node id to (x,t)
    p2n::Array{Int,2}                          # from (x,t) to node id
    order::Vector{Int}                          # Depth first sorted order of nodes. Used for the DP algorithm
end
"""
build a time-expanded network (TEN object) for a given route segment of a line, and a frequency

### Keywords
* `N` - route-network corresponding to the segment area
* `T` - discretized time-interval between stopping times of start and end reference stop in the segment
### Returns
* time-expanded network without any demand associated
"""
function buildGenericTEN_DAR(N::MyNetwork, T::Int)
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
    inA = [Vector{Tuple{Int,Vector{Tuple{Int,Int,Int,Int}},Vector{Float64}}}() for n in 1:count]
    outA = [Vector{Tuple{Int,Vector{Tuple{Int,Int,Int,Int}},Vector{Float64}}}() for n in 1:count]
    # add artificial source and sink arcs
    src = N.src
    node = 1
    node2 = p2n[src, 1]
    push!(inA[node2], (node, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
    push!(outA[node], (node2, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
    earliestEndTime = T #max(0, T-gap)
    for t in earliestEndTime:T
        snk = N.snk
        node = p2n[snk, t+1]
        node2 = count
        push!(inA[node2], (node, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
        push!(outA[node], (node2, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
    end
    for (src, n) in enumerate(N.B2A)
        for (dst, tt) in N.outA[src]
            for t in 0:T-Int(tt)
                node = p2n[src, t+1]
                t2 = t + Int(tt) #+ rand(0:2)
                if t2 <= T
                    node2 = p2n[dst, t2+1]
                    push!(inA[node2], (node, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
                    push!(outA[node], (node2, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
                end
            end
        end
        for t in 0:T-1
            node = p2n[src, t+1]
            node2 = p2n[src, t+2]
            push!(inA[node2], (node, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
            push!(outA[node], (node2, Vector{Tuple{Int,Int,Int,Int}}(), Vector{Float64}()))
        end
    end
    return TEN_DAR(count, inA, outA, n2p, p2n, collect(1:count))
end

"""
Associate each AggregPassengerDAR index (orig, dest, time) to a segment between consecutive reference stop in a line for pickup and dropoff

### Keywords
* `m` - MapData
* `l` - line id
* `line` - Line object
* `freq` - the freq id
* `D` - Dictionary: (orig, time) => AggregPassenger object
### Returns
* Vectors (pickup and dropoff) of AggregPassenger indices divided by segments between consecutive reference stops in a line
"""
function splitPassengersBySegmentsDAR(m::MapData, l::Int, line::Line, freq::Int, D::Dict{Tuple{Int,Int, Int},AggregPassengerDAR})
    # Assign each passenger to one set of sub-paths (that do not allow counting same passenger twice)
    @unpack ref_stop_id, ref_stop_time = line
    segment_demand_pickup = [Vector{Tuple{Int,Int,Int}}() for i in 1:length(ref_stop_id)-1]
    segment_demand_dropoff = [Vector{Tuple{Int,Int,Int}}() for i in 1:length(ref_stop_id)-1]
    for (id, p)  in collect(D)
        if (l, freq) in p.candidateTrips
            # find best segment for pickup
            stop_idx_orig = nearest_node_index(m, p.orig.stop, ref_stop_id)
            if (stop_idx_orig != 1 && stop_idx_orig != length(ref_stop_id) && OpenStreetMapX.distance(m.nodes[p.orig.stop], m.nodes[ref_stop_id[stop_idx_orig-1]]) < OpenStreetMapX.distance(m.nodes[p.orig.stop], m.nodes[ref_stop_id[stop_idx_orig+1]])) || stop_idx_orig == length(ref_stop_id)
                stop_idx_orig -= 1
            end
            push!(segment_demand_pickup[stop_idx_orig], id)
            # find best segment for dropoff
            stop_idx_dest = nearest_node_index(m, p.dest.stop, ref_stop_id)
            if (stop_idx_dest != 1 && stop_idx_dest != length(ref_stop_id) && OpenStreetMapX.distance(m.nodes[p.dest.stop], m.nodes[ref_stop_id[stop_idx_dest-1]]) < OpenStreetMapX.distance(m.nodes[p.dest.stop], m.nodes[ref_stop_id[stop_idx_dest+1]])) || stop_idx_dest == length(ref_stop_id)
                stop_idx_dest -= 1
            end
            push!(segment_demand_dropoff[stop_idx_dest], id)
        end
    end
    return segment_demand_pickup, segment_demand_dropoff
end

"""
update a time-expanded network with the potential demand of a given scenario.
We do this by updating the labels of the arcs

### Keywords
* `m` - MapData
* `R` - RouteSettingDataDAR
* `ti` - stopping time at the start reference stop of the time-expanded network
* `tj` - stopping time at the end stop of the time-expanded network
* `D` - vector of tuples (indices of AggregPassengerDAR objects in R.Demand[s])
* `g` - time-expanded network
* `N` - road-network corresponding to the subpath routing area
* `l` - line id
* `t` - freq id
* `s` - scenario id
### Returns
* updated arc set of the time-expanded network with demand
"""
function updateGenTENwithDemandDAR(m::MapData,
                                R::RouteSettingDataDAR,
                                Inst::InstanceSettingData,
                                ti::Time,
                                tj::Time,
                                D_pickup::Vector{Tuple{Int,Int,Int}},
                                D_dropoff::Vector{Tuple{Int,Int,Int}},
                                g::TEN_DAR,
                                N::MyNetwork,
                                l::Int,
                                t::Int,
                                s::Int,
                                normalized::Bool,
                                )
    @unpack Taxi_travel_times, Demand, Walk_graph = R
    @unpack MaxWait, M, λ, μ, σ, δ, MaxWalk, MaxSchDev, TDmp, TDsp = Inst
    @unpack inA, outA, n2p, p2n = g
    # check separately pickup and dropoff
    for id in D_pickup
        (orig_stop, time_MP, dest_stop) = id
        p = Demand[s][id]
        @unpack num, orig, dest = p
        be_ready_time = p.be_ready_times[(l,t)]
        walk_stops = getLocRangeNodes(m, m.nodes[orig_stop], MaxWalk, Set(N.B2A))
        # compute shortest paths from passenger origin
        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(Walk_graph, m.v[orig_stop])
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
                latest_pickup_SP = min(ceil(Int, Dates.value(Dates.Second(be_ready_time + Dates.Second(walk_time_seconds + MaxWait)-ti))/TDsp), ceil(Int, Dates.value(Dates.Second(tj-ti))/TDsp))
                if earliest_pickup_SP <= latest_pickup_SP
                    for pickup_time_SP in earliest_pickup_SP:latest_pickup_SP
                        node = p2n[src, pickup_time_SP+1]
                        pickup_time = ti + Dates.Second(pickup_time_SP*TDsp)
                        pickup_time_seconds = Dates.value(Dates.Second(pickup_time - (DAY_START - Dates.Second(OPERATIONAL_HORIZON))))
                        if normalized
                            pickup_time_seconds /= Taxi_travel_times[orig_stop, dest_stop]
                        end
                        wait_time_seconds = Dates.value(Dates.Second(pickup_time - (be_ready_time + Dates.Second(walk_time_seconds))))
                        # waiting time MUST be below the maximum allowed
                        for (i, v) in enumerate(outA[node])
                            (node2, C, WT) = v
                            (xC, tC) = node2 != g.V ? n2p[node2] : (0,0)
                            if xC != src # not a holding arc
                                push!(C, (orig_stop, time_MP, dest_stop, num))
                                passengerService = num*(λ*walk_time_seconds + μ*wait_time_seconds - σ*pickup_time_seconds - M)
                                push!(WT, passengerService)    
                                outA[node][i] = (node2, C, WT)
                                for (j, v2) in enumerate(inA[node2])
                                    (node_, C_, WT_) = v2
                                    (src_, t2_) = n2p[node_]
                                    if src_ == src
                                        push!(C_, (orig_stop, time_MP, dest_stop, num))
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
    for id in D_dropoff
        (orig_stop, time_MP, dest_stop) = id
        p = Demand[s][id]
        @unpack num, orig, dest = p
        walk_stops = getLocRangeNodes(m, m.nodes[dest_stop], MaxWalk, Set(N.B2A))
        # compute shortest paths from passenger dest (undirected graph for walking)
        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(Walk_graph, m.v[dest_stop])
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
                # earliest/latest dropoff are given by the requested arrival time and the maximum schedule deviation 
                latest_dropoff_SP = min(floor(Int,  Dates.value(Dates.Second(dest.time + Dates.Second(MaxSchDev - walk_time_seconds)-ti))/TDsp),ceil(Int, Dates.value(Dates.Second(tj-ti))/TDsp))
                earliest_dropoff_SP = max(ceil(Int, Dates.value(Dates.Second(dest.time - Dates.Second(MaxSchDev + walk_time_seconds)-ti))/TDsp), 0)
                if earliest_dropoff_SP <= latest_dropoff_SP
                    for dropoff_time_SP in earliest_dropoff_SP:latest_dropoff_SP
                        node = p2n[src, dropoff_time_SP+1]
                        # pickup_time = ti + Dates.Second(pickup_time_SP*TDsp)
                        dropoff_time_seconds = Dates.value( Dates.Second(ti + Dates.Second(dropoff_time_SP*TDsp) - (DAY_START - Dates.Second(OPERATIONAL_HORIZON))))
                        if normalized
                            dropoff_time_seconds /= Taxi_travel_times[orig_stop, dest_stop]
                        end
                        # compute weighted delay/earliness
                        arrival_time = ti + Dates.Second(dropoff_time_SP*TDsp + walk_time_seconds)
                        delay_seconds = Dates.value(Dates.Second(arrival_time - dest.time))
                        weighted_delay = delay_seconds < -EPS ? abs((δ/2)*delay_seconds) : abs(δ*delay_seconds)
                        if normalized
                            weighted_delay /= Taxi_travel_times[orig_stop, dest_stop]
                        end
                        for (i, v) in enumerate(outA[node])
                            (node2, C, WT) = v
                            (xC, tC) = node2 != g.V ? n2p[node2] : (0,0)
                            if xC != src # not a holding arc
                                push!(C, (orig_stop, time_MP, dest_stop, -num))
                                passengerService = num*(λ*walk_time_seconds + σ*dropoff_time_seconds + weighted_delay)
                                push!(WT, passengerService)
                                outA[node][i] = (node2, C, WT)
                                for (j, v2) in enumerate(inA[node2])
                                    (node_, C_, WT_) = v2
                                    (src_, t2_) = n2p[node_]
                                    if src_ == src
                                        push!(C_, (orig_stop, time_MP, dest_stop, -num))
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
equivalent of updateGenTENwithDemandDAR() for the fixed-transit version,
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
function updateGenTENwithDemandForTransitDAR(m::MapData, 
    R::RouteSettingDataDAR,
    Inst::InstanceSettingData,
    ti::Time,
    tj::Time,
    D_pickup::Vector{Tuple{Int,Int,Int}},
    D_dropoff::Vector{Tuple{Int,Int,Int}},
    g::TEN_DAR,
    N::MyNetwork,
    l::Int,
    t::Int,
    s::Int,
    normalized::Bool,
    )
    @unpack Taxi_travel_times, Demand, Walk_graph = R
    @unpack MaxWait, M, λ, μ, σ, δ, MaxWalk, MaxSchDev, TDmp, TDsp,  = Inst
    @unpack inA, outA, n2p, p2n = g
    # check separately pickup and dropoff
    for id in D_pickup
        (orig_stop, time_MP, dest_stop) = id
        p = Demand[s][id]
        @unpack num, orig, dest = p
        be_ready_time = p.be_ready_times[(l,t)]
        walk_stops = getLocRangeNodes(m, m.nodes[orig_stop], MaxWalk, Set(R.Lines[l].transit_stop_id))
        # compute shortest paths from passenger origin
        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(Walk_graph, m.v[orig_stop])
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
                        pickup_time_seconds = Dates.value(Dates.Second(pickup_time - (DAY_START - Dates.Second(OPERATIONAL_HORIZON))))
                        if normalized
                            pickup_time_seconds /= Taxi_travel_times[orig_stop, dest_stop]
                        end
                        wait_time_seconds = Dates.value(Dates.Second(pickup_time - (be_ready_time + Dates.Second(walk_time_seconds))))
                        # waiting time MUST be below the maximum allowed
                        for (i, v) in enumerate(outA[node])
                            (node2, C, WT) = v
                            (xC, tC) = node2 != g.V ? n2p[node2] : (0,0)
                            if xC != src # not a holding arc
                                push!(C, (orig_stop, time_MP, dest_stop, num))
                                passengerService = num*(λ*walk_time_seconds + μ*wait_time_seconds - σ*pickup_time_seconds - M)
                                push!(WT, passengerService)
                                outA[node][i] = (node2, C, WT)
                                for (j, v2) in enumerate(inA[node2])
                                    (node_, C_, WT_) = v2
                                    (src_, t2_) = n2p[node_]
                                    if src_ == src
                                        push!(C_, (orig_stop, time_MP, dest_stop, num))
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
    for id in D_dropoff
        (orig_stop, time_MP, dest_stop) = id
        p = Demand[s][id]
        @unpack num, orig, dest = p
        walk_stops = getLocRangeNodes(m, m.nodes[dest_stop], MaxWalk, Set(R.Lines[l].transit_stop_id))
        # compute shortest paths from passenger origin
        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(Walk_graph, m.v[dest_stop])
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
                    latest_dropoff = min(dest.time + Dates.Second(MaxSchDev - walk_time_seconds),tj)
                    earliest_dropoff = max(dest.time - Dates.Second(MaxSchDev + walk_time_seconds), ti)
                    # check if it can be dropped off at the scheduled time
                    if earliest_dropoff <= R.Lines[l].transit_stop_time[cand_stop][t] <= latest_dropoff
                        dropoff_time_SP = ceil(Int,  Dates.value(Dates.Second(R.Lines[l].transit_stop_time[cand_stop][t]-ti))/TDsp)
                        node = p2n[src, dropoff_time_SP+1]
                        dropoff_time = ti + Dates.Second(dropoff_time_SP*TDsp)
                        dropoff_time_seconds = Dates.value(Dates.Second(dropoff_time - (DAY_START - Dates.Second(OPERATIONAL_HORIZON))))
                        if normalized
                            dropoff_time_seconds /= Taxi_travel_times[orig_stop, dest_stop]
                        end
                        # compute weighted delay/earliness
                        arrival_time = dropoff_time + Dates.Second(walk_time_seconds)
                        delay_seconds = Dates.value(Dates.Second(arrival_time - dest.time))
                        weighted_delay = delay_seconds < -EPS ? abs((δ/2)*delay_seconds) : abs(δ*delay_seconds)
                        if normalized
                            weighted_delay /= Taxi_travel_times[orig_stop, dest_stop]
                        end
                        for (i, v) in enumerate(outA[node])
                            (node2, C, WT) = v
                            (xC, tC) = node2 != g.V ? n2p[node2] : (0,0)
                            if xC != src # not a holding arc
                                push!(C, (orig_stop, time_MP, dest_stop, -num))
                                passengerService = num*(λ*walk_time_seconds + σ*dropoff_time_seconds + weighted_delay)
                                push!(WT, passengerService)
                                outA[node][i] = (node2, C, WT)
                                for (j, v2) in enumerate(inA[node2])
                                    (node_, C_, WT_) = v2
                                    (src_, t2_) = n2p[node_]
                                    if src_ == src
                                        push!(C_, (orig_stop, time_MP, dest_stop, -num))
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
run a label setting algorithm on graph g (TEN), with two labels: customers picked up, and stopped picked up (i.e., saved waiting time).

### Keywords
* `g` - TEN_DAR
* source` - source node
* `nodeOrder` - topological order of nodes
* `capacity` - maximum number of passengers to pick up
### Returns
* non-dominated labels at each node, and the parent node and label of each node
"""
function fullEnumAlgDAR(
    g::TEN_DAR,
    source::Int64,
    nodeOrder::Vector{Int64},
    capacity::Int,
    )
    nn = g.V
    labels = [Vector{Pair{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}}}([Pair(Vector{Tuple{Int, Int, Int, Int}}(), Vector{Int}([LARGE_NUMBER+1]))]) for n in 1:nn]
    labels[source] = [Pair(Vector{Tuple{Int, Int, Int, Int}}(), Vector{Float64}())]
    # keep track of prev of each label
    prev = [Vector{Tuple{Int, Pair{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}}}}([(0, Pair(Vector{Tuple{Int, Int, Int, Int}}(), Vector{Float64}()))]) for n in 1:nn]
    for i in nodeOrder
        if i != source
            for (src, cus, time) in g.inA[i]
                for label in labels[src]
                    # ignore nodes that cannot be reached
                    if sum(label.second) < LARGE_NUMBER
                        # compute num of pax at label
                        numCust = 0
                        for c_ in label.first
                            numCust += c_[4] # negative if dropoff
                        end
                        # initialize vector for new labels, which will be label + all comb of the customers in the arc
                        newLabels = Vector{Pair{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}}}()
                        comb = collect(powerset(cus,0,length(cus)))
                        for cus_ in comb
                            # compute labels for dst
                            C = deepcopy(label.first)
                            WT = deepcopy(label.second)
                            # check if they can be picked up
                            fit = true
                            Q_ = deepcopy(numCust)
                            for (ic, c) in enumerate(cus_)
                                # if customer not already picked up
                                if !(c in C)
                                    Q_ += cus_[ic][4]
                                end
                            end
                            if (Q_ > capacity + EPS) || (Q_ < -capacity - EPS)
                                fit = false
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

function fullLabelSettingAlgInCG_DAR(
    l::Int,
    t::Int,
    g::TEN_DAR,
    source::Int64,
    nodeOrder::Vector{Int64},
    capacity::Int,
    d4,
    d5,
    )
    nn = g.V
    labels = [Vector{Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}([(Vector{Tuple{Int, Int, Int, Int}}(), Vector{Float64}(), typemax(Float64))]) for n in 1:nn]
    labels[source] = [(Vector{Tuple{Int, Int, Int, Int}}(), Vector{Int}(), 0.)]
    # keep track of prev of each undominated label
    prev = [Vector{Tuple{Int, Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}}([(0, (Vector{Tuple{Int, Int, Int, Int}}(), Vector{Float64}(), 0.))]) for n in 1:nn]
    for i in nodeOrder
        if i != source
            for (src, cus, time) in g.inA[i] # cus should be vector (there could be more than one customer related to a node)
                for label in labels[src]
                    if label[3] < LARGE_NUMBER
                        numCust = 0
                        for c_ in label[1]
                            numCust += c_[4]
                        end
                        newLabels = Vector{Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}()
                        comb = collect(powerset(cus, 0, length(cus))) # at most we want to inspect X (X = capacty - passengers already in label) elements (worst case 1 pax in each el)
                        existingCombs = Vector{Vector{Tuple{Int, Int, Int, Int}}}()
                        for cus_ in comb
                            C = deepcopy(label[1])
                            WT = deepcopy(label[2])
                            RC = deepcopy(label[3])
                            # check if they can be picked up
                            fit = true
                            Q_ = deepcopy(numCust)
                            for (ic, c) in enumerate(cus_)
                                if !(c in C)
                                    Q_ += cus_[ic][4]
                                end
                            end
                            if (Q_ > capacity + EPS) || (Q_ < -capacity - EPS)
                                fit = false
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
                                        RC -= -time[icC]
                                        
                                        if cus_[ic][4] > EPS
                                            RC -=  d4[(cus_[ic][1], cus_[ic][2], cus_[ic][3]),(l,t)] + d5[(cus_[ic][1], cus_[ic][2], cus_[ic][3]),(l,t)]
                                        else
                                            RC -=  -d5[(cus_[ic][1], cus_[ic][2], cus_[ic][3]),(l,t)]
                                        end
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
                            # if any label in labels[dst] is changed
                            # add dst to Q
                            # go through all labels and see which are dominated (to delete), before adding the new one
                            toDelete = Set{Int}()
                            dominated = false
                            for (idx, labY) in enumerate(labels[i])
                                # check if labY is dominated by (C, WT)
                                # check coinciding customers (extra in labY that are not in C)
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

function labelSettingAlgInCG_DAR(
    l::Int,
    t::Int,
    g::TEN_DAR,
    source::Int64,
    nodeOrder::Vector{Int64},
    capacity::Int,
    d4,
    d5,
    )
    nn = g.V
    labels = [Vector{Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}([(Vector{Tuple{Int, Int, Int, Int}}(), Vector{Float64}(), typemax(Float64))]) for n in 1:nn]
    labels[source] = [(Vector{Tuple{Int, Int, Int, Int}}(), Vector{Int}(), 0.)]
    # keep track of prev of each undominated label
    prev = [Vector{Tuple{Int, Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}}([(0, (Vector{Tuple{Int, Int, Int, Int}}(), Vector{Float64}(), 0.))]) for n in 1:nn]
    for i in nodeOrder
        if i != source
            for (src, cus, time) in g.inA[i] # cus should be vector (there could be more than one customer related to a node)
                for label in labels[src]
                    if label[3] < LARGE_NUMBER
                        # compute labels for dst
                        C = deepcopy(label[1])
                        WT = deepcopy(label[2])
                        RC = deepcopy(label[3])
                        Q = 0 # current load
                        for c_ in C
                            Q += c_[4]
                        end
                        # check customers in c
                        for (ic, c) in enumerate(cus)
                            # if it is already picked up, check if t is better
                            if c in C
                                icS = findfirst(isequal(c), C)
                                if time[ic] <  WT[icS] - EPS
                                    RC -= WT[icS] - time[ic]
                                    WT[icS] = time[ic]
                                end
                            else
                                # if it is new: add it if the vehicle has capacity
                                # we always pickup all possible pax in a stop (it means they can make it)
                                if cus[ic][4] + Q <= capacity && cus[ic][4] + Q >= -capacity
                                    redcost = time[ic]
                                    if cus[ic][4] > EPS
                                        redcost -=  d4[(cus[ic][1], cus[ic][2], cus[ic][3]),(l,t)] + d5[(cus[ic][1], cus[ic][2], cus[ic][3]),(l,t)]
                                    else
                                        redcost -=  -d5[(cus[ic][1], cus[ic][2], cus[ic][3]),(l,t)]
                                    end
            
                                    if redcost < -EPS
                                        push!(C, cus[ic])
                                        push!(WT, time[ic])
                                        RC += redcost
                                        Q += cus[ic][4]
                                    end
                                end
                            end
                        end
                        # keep the label with the lowest RC
                        toDelete = Set{Int}()
                        dominated = false
                        for (idx, labY) in enumerate(labels[i])
                            RC_i = labY[3]
                            if RC < RC_i - EPS
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
function getPathsAndLabelsDAR(node::Int64, src::Int64,
                            labels::Array{Array{Pair{Array{Tuple{Int,Int,Int,Int},1},Array{Float64,1}},1},1},
                            prevNode::Array{Array{Tuple{Int,Pair{Array{Tuple{Int,Int,Int,Int},1},Array{Float64,1}}},1},1},
                            ) 
    paths = Vector{Vector{Int64}}()
    labs = Vector{Pair{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}}}()
    for idx in eachindex(prevNode[node])
        path = Vector{Int64}()
        push!(path, node)
        (pN, pLab) = prevNode[node][idx]
        while pN != 0 
            push!(path, pN)
            i = findfirst(isequal(pLab), labels[pN])
            (pN, pLab) = prevNode[pN][i]
        end
        # we only want paths that start from our "source" node
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
function getPathsAndLabelsInCG_DAR(node::Int64, src::Int64,
    labels::Vector{Vector{Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}},
    prevNode::Vector{Vector{Tuple{Int, Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}}},
    ) 
    paths = Vector{Vector{Int64}}()
    labs = Vector{Tuple{Vector{Tuple{Int, Int, Int, Int}}, Vector{Float64}, Float64}}()
    for idx in eachindex(prevNode[node])
        path = Vector{Int64}()
        push!(path, node)
        (pN, pLab) = prevNode[node][idx]
        while pN != 0 
            push!(path, pN)
            i = findfirst(isequal(pLab), labels[pN])
            (pN, pLab) = prevNode[pN][i]
        end
        # we only want paths that start from our "source" node
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
function generateSubPathSetDAR(m::MapData,
    R::RouteSettingDataDAR,        # model input
    Inst::InstanceSettingData,  # instance data
    CG::Bool,                   # solve second-stage with column generation (if true, we start with an empty set of sub-paths)
    TRANSIT::Bool,
    )
    @unpack Lines, All_stops, Num_freq, Walk_graph, Demand, Taxi_travel_times = R
    @unpack K, MaxDev, MaxWalk, MaxWait, MaxSchDev, routeSpeedFactor, TDmp, TDsp, M, λ, μ, σ, δ = Inst
    
    all_subpaths = [[Vector{PoolDAR}() for t in eachindex(Lines[l].freq)] for l in eachindex(Lines)]
    subpath_road_networks = [[Vector{MyNetwork}() for i in 1:length(Lines[l].ref_stop_id)-1] for l in eachindex(Lines)]
    generic_subpath_graphs = [[[Vector{TEN_DAR}() for i in 1:length(Lines[l].ref_stop_id)-1] for t in eachindex(Lines[l].freq)] for l in eachindex(Lines)]
    all_subpath_graphs = [[[Vector{TEN_DAR}() for i in 1:length(Lines[l].ref_stop_id)-1] for t in eachindex(Lines[l].freq)] for l in eachindex(Lines), s in eachindex(Demand)]
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
            for i in 1:length(ref_stop_id)-1, j in i+1:min(i + K, length(ref_stop_id))
                src = ref_stop_id[i]                            # start ref stop of sub-path
                idxI = findfirst(isequal(src), ref_route_nodes)
                snk = ref_stop_id[j]                            # end ref stop of sub-path
                idxJ = findfirst(isequal(snk), ref_route_nodes[idxI+1:end])
                idxJ += idxI
                ti = ref_stop_time[src][t]
                tj = ref_stop_time[snk][t]
                subpath_graphT = ceil(Int, Dates.value(Dates.Second(tj-ti))/TDsp)
                # build generic TEN_DAR
                g = buildGenericTEN_DAR(subpath_road_networks[l][i][j-i], subpath_graphT)
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
            for s in eachindex(Demand)
                # passengers to be considered in the subproblem
                # divide scenario demand to its closest segment (a segment is the space betwen two consecutive ref stops). This is to avoid a solution picking up the same passenger multiple times.
                segment_demand_pickup, segment_demand_dropoff = splitPassengersBySegmentsDAR(m, l, line, t, Demand[s])
                pool = PoolDAR(load_expanded_graph.V)
                count = 0
                for p in eachindex(Demand[s])
                    pool.paxPU[p] = Vector{Int}()
                    pool.paxDO[p] = Vector{Int}()
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
                    subpath_segment_demand_pickup = Vector{Tuple{Int,Int,Int}}()
                    subpath_segment_demand_dropoff = Vector{Tuple{Int,Int,Int}}()
                    for idxR in i:j-1
                        append!(subpath_segment_demand_pickup, segment_demand_pickup[idxR])
                        append!(subpath_segment_demand_dropoff, segment_demand_dropoff[idxR])
                    end
                    g.inA, g.outA = TRANSIT ? updateGenTENwithDemandForTransitDAR(m, R, Inst, ti, tj, subpath_segment_demand_pickup, subpath_segment_demand_dropoff, g, subpath_road_networks[l][i][j-i], l,t,s, NORMALIZED) : updateGenTENwithDemandDAR(m, R, Inst, ti, tj, subpath_segment_demand_pickup, subpath_segment_demand_dropoff, g, subpath_road_networks[l][i][j-i], l,t,s, NORMALIZED)
                    push!(all_subpath_graphs[l,s][t][i], g)
                    t1 = ceil(Int, Dates.value(Dates.Second(ref_stop_time[src][t]-t0))/TDsp)
                    t2 = ceil(Int, Dates.value(Dates.Second(ref_stop_time[snk][t]-t0))/TDsp)
                    emptySubpath = false # we should always guarantee an "empty" path.
                    if !CG # precompute subset of sub-paths
                        labels, prev = fullEnumAlgDAR(g, 1, g.order, capacity)
                        # recreate paths
                        paths, labs = getPathsAndLabelsDAR(g.V, 1, labels, prev)
                        for (idx, label) in enumerate(labs)
                            WT = sum(label.second)
                            if WT < LARGE_NUMBER
                                Q = 0
                                paxPU = Vector{Tuple{Int,Int,Int,Int}}()
                                paxDO = Vector{Tuple{Int,Int,Int,Int}}()
                                for (orig, time, dest, num) in label.first
                                    Q += num
                                    if num > EPS
                                        push!(paxPU, (orig, time, dest, num))
                                    else
                                        push!(paxDO, (orig, time, dest, num))
                                    end
                                end
                                cost = WT 
                                if isempty(label.first)
                                    emptySubpath = true
                                end
                                for c in max(0,-Q):min(capacity,capacity-Q)
                                    count += 1
                                    node1 = load_expanded_graph.p2n[i, t1+1, c+1]
                                    node2 = load_expanded_graph.p2n[j, t2+1, c+Q+1]
                                    
                                    subP = SubPathDAR(count, cost, Q, paxPU, paxDO, paths[idx], node1, node2)
                                    push!(pool.all, subP)
                                    push!(pool.in[node2], count)
                                    push!(pool.out[node1], count)
                                    for (orig, time, dest, num) in label.first
                                        if num > EPS
                                            push!(pool.paxPU[orig, time, dest], count)
                                        else
                                            push!(pool.paxDO[orig, time, dest], count)
                                        end
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
                            subP = SubPathDAR(count, 0., 0, Vector{Tuple{Int, Int, Int, Int}}(), Vector{Tuple{Int, Int, Int, Int}}(), Vector{Int}(), node1, node2)
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
                subP = SubPathDAR(count, 0., 0, Vector{Tuple{Int, Int, Int, Int}}(), Vector{Tuple{Int, Int, Int, Int}}(), Vector{Int}(), node1, node2)
                push!(pool.all, subP)
                push!(pool.in[node2], count)
                push!(pool.out[node1], count)
                # add set of arcs to the artificial sink node from the last stop in the route
                for c in 0:capacity
                    count += 1
                    node1 = load_expanded_graph.p2n[length(ref_stop_id), subproblemT+1, c+1]
                    node2 = load_expanded_graph.V
                    subP = SubPathDAR(count, 0., 0, Vector{Tuple{Int, Int, Int, Int}}(), Vector{Tuple{Int, Int, Int, Int}}(), Vector{Int}(), node1, node2)
                    push!(pool.all, subP)
                    push!(pool.in[node2], count)
                    push!(pool.out[node1], count)
                end
                push!(all_subpaths[l][t], pool)
            end
        end
    end
    println("DONE!")
    return all_subpaths, subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs
end

"""
Build entire two-stage stochastic model.

### Keywords
* `R` - RouteSettingData
* `Subpaths` - Pool of sub-paths
* `Inst` - instance data
* `time_limit_sec` - maximum solve time in seconds
* `mip_gap` - epsilon tolerance
### Returns
* JuMP Model
"""
function directTwoStageModelDAR(
    R::RouteSettingDataDAR,
    Subpaths::Vector{Vector{Vector{PoolDAR}}},
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
        y[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, i in eachindex(Subpaths[l][t][s].all)], Bin
    end

    JuMP.@expression(m,
        first_stage_costs,
        sum( WEIGHT_LINE_COST*R.Lines[l].cost * x[l, t] for l in 1:numL, t in eachindex(R.Lines[l].freq))
    )
    JuMP.@expression(m,
        second_stage_costs,
        # generalized cost of travel of passengers served by sub-paths
        sum( R.Pi[s] * (
            sum( a.cost * y[l, t, s, a.id] for a in Subpaths[l][t][s].all))
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
        c1[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(y[l,t,s,id] for id in Subpaths[l][t][s].out[1]) == x[l, t]
        c2[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(y[l,t,s,id] for id in Subpaths[l][t][s].in[length(Subpaths[l][t][s].in)]) == x[l, t]
        c3[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS, i in 2:length(Subpaths[l][t][s].in)-1], sum(y[l,t,s,id] for id in Subpaths[l][t][s].out[i]) - sum(y[l,t,s,id] for id in Subpaths[l][t][s].in[i]) == 0
        c4[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], sum(y[l,t,s,id] for id in Subpaths[l][t][s].paxPU[p]) <= z[s,p,(l,t)]
        c5[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], sum(y[l,t,s,id] for id in Subpaths[l][t][s].paxPU[p]) - sum(y[l,t,s,id] for id in Subpaths[l][t][s].paxDO[p]) == 0
    end

    return m
end

"""
function used in the column generation process to add the new column to the sub-problem,
we do that by adding the column in the rows with non-zero coefficients and in the obj. function

### Keywords
* `m` - JuMP.Model of sub-problem
* `sp` - SubPath object
* `i` - tail node in TSLN
* `j` - head node in TSLN
* `δ` - weight for arrival delay
### Returns
* model m gets the new variable added (no return)
"""
function addSubPathDAR!(
    l::Int,
    t::Int,
    s::Int,
    m::JuMP.Model,
    sp::SubPathDAR,
    i::Int,
    j::Int,
    PUs::Vector{Tuple{Int,Int,Int}},
    DOs::Vector{Tuple{Int,Int,Int}},
    name::String="y",
    )
    yNew = @variable(m, base_name=name, lower_bound=0)
    touchedConstraints = ConstraintRef[]
    vals = Float64[]
    push!(touchedConstraints, m[:c3][i])
    push!(vals, 1.)
    push!(touchedConstraints, m[:c3][j])
    push!(vals, -1.)

    for (orig,time,dest,num) in sp.paxPU
        push!(touchedConstraints, m[:c4][(orig,time,dest),(l,t)])
        push!(vals, 1.)
    end
    for (orig,time,dest) in setdiff(PUs,DOs)
        push!(touchedConstraints, m[:c5][(orig,time,dest),(l,t)])
        push!(vals, 1.)
    end
    for (orig,time,dest) in setdiff(DOs,PUs)
        push!(touchedConstraints, m[:c5][(orig,time,dest),(l,t)])
        push!(vals, -1.)
    end
    set_objective_coefficient(m, yNew, sp.cost)
    set_normalized_coefficient.(touchedConstraints, yNew, vals)
    push!(m[:y],yNew)
end

"""
build a sub-problem given a pool of sub-paths
"generic" means that all RHS is = 1 (all active), except for flow conservation constraints.

### Keywords
* `SPs` - pool of sub-paths (arcs) for the sub-problem
* `Pax` - pool of sub-paths that cover a specific passenger demand: tuple (loc, time)
* `time_limit_sec` - time limit to solve the problem
* `mip_gap` - accepted optimality gap
### Returns
* JuMP model of the sub-problem
"""
function buildGenericSecondStageDAR(
    R::RouteSettingDataDAR,
    Subpaths::PoolDAR,
    l::Int,
    t::Int,
    s::Int,
    time_limit_sec::Int64=MP_TIME_LIMIT,
    mip_gap::Float64=MIP_GAP,
    num_threads::Int64=NUM_THREADS,
    num_focus::Int64=NUM_FOCUS,
    )

    @unpack Num_freq, Demand, Trip_passengers, Freq_Overlap = R
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "TimeLimit" => time_limit_sec, "MIPGap" => mip_gap, "Threads" => num_threads, "OutputFlag" => 0, "NumericFocus" => num_focus))
    # set_silent(m)

    @variable(m, y[i in eachindex(Subpaths.all)], lower_bound=0)

    @objective(m, Min, sum( a.cost * y[a.id] for a in Subpaths.all))

    @constraint(m, c1, sum(y[id] for id in Subpaths.out[1]) == 1)
    @constraint(m, c2, sum(y[id] for id in Subpaths.in[length(Subpaths.in)]) == 1)
    @constraint(m, c3[i in 2:length(Subpaths.in)-1], sum(y[id] for id in Subpaths.out[i]) - sum(y[id] for id in Subpaths.in[i]) == 0)
    @constraint(m, c4[p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips; (l_,t_) == (l,t)], sum(y[id] for id in Subpaths.paxPU[p]) <= 1)
    @constraint(m, c5[p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips; (l_,t_) == (l,t)], sum(y[id] for id in Subpaths.paxPU[p]) - sum(y[id] for id in Subpaths.paxDO[p]) == 0)

    return m
end
"""
Build first stage model.

### Keywords
* `R` - RouteSettingDataDAR
* `Inst` - InstanceSettingData
* `time_limit_sec` - maximum solve time in seconds
* `mip_gap` - epsilon tolerance
### Returns
* JuMP Model
"""
function firstStageModelDAR(
    R::RouteSettingDataDAR,
    Inst::InstanceSettingData,
    time_limit_sec::Int64=MP_TIME_LIMIT,
    mip_gap::Float64=MIP_GAP,
    num_threads::Int64=NUM_THREADS,
    num_focus::Int64=NUM_FOCUS,
    )

    @unpack numL, numS, TDmp, TDsp, Fleet, Theta = Inst 
    @unpack Num_freq, Demand, Trip_passengers, Freq_Overlap = R
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "TimeLimit" => time_limit_sec, "MIPGap" => mip_gap, "Threads" => num_threads, "OutputFlag" => 0, "NumericFocus" => num_focus #=, "Method" => 1 "FeasibilityTol" => 1e-8=#)) # "Cuts" => 3


    JuMP.@variables m begin
        # if passenger p is assigned to trip (l,t) in scenario s
        z[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], Bin                             
        # if line l operates at frequency t
        x[l in 1:numL, eachindex(R.Lines[l].freq)], Bin
        # recourse for line l, freq t, and scenario s
        theta[l in 1:numL, eachindex(R.Lines[l].freq), 1:numS]
    end

    JuMP.@expression(m,
        first_stage_costs,
        sum( WEIGHT_LINE_COST*R.Lines[l].cost * x[l, t] for l in 1:numL, t in eachindex(R.Lines[l].freq))
    )

    JuMP.@expression(m,
        second_stage_costs,
        sum(R.Pi[s] * theta[l,t,s] for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS)
    )

    JuMP.@objective(m, Min,
        first_stage_costs +
        second_stage_costs
    )

    JuMP.@constraints m begin
        fleet[t in 1:Num_freq], sum(x[l, t_] for l in 1:numL, t_ in Freq_Overlap[l,t]) <= Fleet
        cover[s in 1:numS, p in eachindex(Demand[s])], sum(z[s,p,(l,t)] for (l,t) in Demand[s][p].candidateTrips) <= 1
        minLoad[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(Demand[s][p].num*z[s,p,(l,t)] for p in Trip_passengers[l,t,s]) >= (1 - Theta)R.Lines[l].capacity*x[l, t]
        maxLoad[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], sum(Demand[s][p].num*z[s,p,(l,t)] for p in Trip_passengers[l,t,s]) <= (1 + Theta)R.Lines[l].capacity*x[l, t]

        recourse[l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS], theta[l,t,s] >= -LARGE_NUMBER
    end

    return m
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
function runAlgDAR(
    m::MapData,
    R::RouteSettingDataDAR,
    Inst::InstanceSettingData,
    Subpaths::Vector{Vector{Vector{PoolDAR}}},                          
    all_subpath_graphs::Array{Vector{Vector{Vector{TEN_DAR}}},2},      
    all_load_expanded_graphs::Vector{Vector{TSLGraph}},                 
    subpath_road_networks,
    CG::Bool,       # solve second-stage with column generation
    Heur::Bool,     # heuristic CG or not
    rootnode::Bool, # relax MP and solve rootnode with Benders
    num_focus::Int64=NUM_FOCUS,
    timelimit::Int=MP_TIME_LIMIT, # time limit for the algorithm execution
    )
    @unpack numL, numS, TDmp, TDsp, K, M, δ = Inst
    @unpack Lines, Num_freq, Demand, Taxi_travel_times = R
    # build first stage model
    MP = firstStageModelDAR(R, Inst, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus)
    # we want to solve the root node first (set to false if not)
    local undo_relax
    if rootnode
        undo_relax = relax_integrality(MP)
    end
    subproblems = [[[buildGenericSecondStageDAR(R, Subpaths[l][t][s], l,t,s, MP_TIME_LIMIT, MIP_GAP, NUM_THREADS, num_focus) for s in 1:numS] for t in eachindex(R.Lines[l].freq)] for l in 1:numL]

    done = false
    LB = typemin(Float64)       # lower bound
    UB = typemax(Float64)       # upper bound
    cuts = ConstraintRef[]      # pool of Benders cuts
    it = 0                      # iteration counter
    RNtime = 0.                 # time spent solving the root node
    solveTMP = 0.0              # time spent in MP
    solveTSP = 0.0              # time spent in SP
    FScosts = 0.0               # first-stage costs (0)
    SScosts = 0.0               # second stage costs
    sol = SolDAR(numS)          # initialize solution
    num_second_stage_vars = 0

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
        println("OBJ MP: ", objMP, " in ", round(solvetime, digits=2), " seconds.")
        xVal = value.(MP[:x])               # first stage X vars
        zVal = value.(MP[:z])               # first stage Z vars
        thetaVal = value.(MP[:theta])       # first stage Theta vars

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
                    if termination_status(subproblems[l][t][s]) != MOI.OPTIMAL
                        println("Term status of SP $l $t $s: ", termination_status(subproblems[l][t][s]))
                    end
                    objSP = objective_value(subproblems[l][t][s])
                    # obtain dual solution
                    d1 = dual(subproblems[l][t][s][:c1])
                    d2 = dual(subproblems[l][t][s][:c2])
                    d3 = dual.(subproblems[l][t][s][:c3])
                    d4 = dual.(subproblems[l][t][s][:c4])
                    d5 = dual.(subproblems[l][t][s][:c5])

                    doneCG = true # we are done unless we find a column with negative reduced cost
                    for i in 1:length(Lines[l].ref_stop_id)-1, j in i+1:min(i + K, length(Lines[l].ref_stop_id))
                        g = all_subpath_graphs[l,s][t][i][j-i]
                        src = Lines[l].ref_stop_id[i]
                        snk = Lines[l].ref_stop_id[j]
                        labels, prev = Heur ? labelSettingAlgInCG_DAR(l,t, g, 1, g.order, Lines[l].capacity, d4, d5) : fullLabelSettingAlgInCG_DAR(l,t,g, 1, g.order, Lines[l].capacity, d4, d5)
                        # recreate paths
                        paths, labs = getPathsAndLabelsInCG_DAR(g.V, 1, labels, prev)
                        t1 = ceil(Int, Dates.value(Dates.Second(R.Lines[l].ref_stop_time[src][t]-t0))/TDsp)
                        t2 = ceil(Int, Dates.value(Dates.Second(R.Lines[l].ref_stop_time[snk][t]-t0))/TDsp)
                        for (idx, label) in enumerate(labs)
                            WT = sum(label[2]) # total walk + wait time
                            if WT < LARGE_NUMBER 
                                # compute number of passengers picked up
                                Q = 0
                                paxPU = Vector{Tuple{Int,Int,Int,Int}}()
                                paxDO = Vector{Tuple{Int,Int,Int,Int}}()
                                PUs = Vector{Tuple{Int,Int,Int}}()
                                DOs = Vector{Tuple{Int,Int,Int}}()
                                for (orig, time, dest, num) in label[1]
                                    Q += num
                                    if num > EPS
                                        push!(paxPU, (orig, time, dest, num))
                                        push!(PUs, (orig, time, dest))
                                    else
                                        push!(paxDO, (orig, time, dest, num))
                                        push!(DOs, (orig, time, dest))
                                    end
                                end
                                cost = WT 
                                for c in max(0,-Q):min(Lines[l].capacity,Lines[l].capacity-Q) # add the sub-path for all combinations of "occupancy levels"
                                    count = length(Subpaths[l][t][s].all) # id of the sub-path (also helps keep count)
                                    node1 = all_load_expanded_graphs[l][t].p2n[i, t1+1, c+1]   # find tail node in SP
                                    node2 = all_load_expanded_graphs[l][t].p2n[j, t2+1, c+Q+1] # find head node in SP
                                    # compute reduced cost
                                    rc = label[3] - d3[node1] + d3[node2] # reduced cost fom lab setting alg + dual source and sink
                                    
                                    if rc < -0.001
                                        count += 1
                                        doneCG = false # if negative reduced cost, we have not found the LP optimal solution yet
                                        subP = SubPathDAR(count, cost, Q, paxPU, paxDO, paths[idx], node1, node2)
                                        push!(Subpaths[l][t][s].all, subP)
                                        push!(Subpaths[l][t][s].in[node2], count)
                                        push!(Subpaths[l][t][s].out[node1], count)
                                        for (orig, time, dest, num) in label[1]
                                            if num > EPS
                                                push!(Subpaths[l][t][s].paxPU[orig, time, dest], count)
                                            else
                                                push!(Subpaths[l][t][s].paxDO[orig, time, dest], count)
                                            end
                                        end
                                        # add column to subproblem
                                        addSubPathDAR!(l,t,s, subproblems[l][t][s], subP, node1, node2, PUs, DOs)                                  
                                    end
                                end
                                
                            end
                        end
                    end
                    # if we are done add objective value and generate Benders cut
                    if doneCG
                        push!(objSPS[s], objSP)
                        if !Heur || (Heur && objSP > thetaVal[l,t,s] + EPS) 
                            doneSPs = false # not done with Benders
                            # if optimal generate full benders cut with all the duals and add to RMP
                            cRef = @constraint(MP, (d1 + d2) * MP[:x][l, t] + sum(d4[p,(l, t)] * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
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
                objSP = objective_value(subproblems[l][t][s])
                push!(objSPS[s], objSP)
                if true #objSP > thetaVal[l, t, s] + 0.001
                    doneSPs = false
                    d1 = dual(subproblems[l][t][s][:c1])
                    d2 = dual(subproblems[l][t][s][:c2])
                    d4 = dual.(subproblems[l][t][s][:c4])
                    cRef = @constraint(MP, (d1 + d2) * MP[:x][l, t] + sum(d4[p,(l, t)] * MP[:z][s,p,(l, t)] for p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips if (l_,t_) == (l,t)) <= MP[:theta][l, t, s])
                    push!(cuts, cRef)
                    # update RHS after adding first cut
                    if it == 1
                        delete(MP, MP[:recourse][l,t,s])
                    end
                end
            end
            solveTSP += time() - startCG
        end
        # compute second-stage costs
        SScosts = sum(R.Pi[s] * sum(objSPS[s]) for s in 1:numS)
        ub = SScosts + FScosts
        UB = ub < UB - 0.01 ? ub : UB # update upper bound if better
        println("IT $it, LB = $LB, UB = $UB, numCuts = ", length(cuts))
        # STOP if no more cuts can be added or if we have converged
        if rootnode && (doneSPs || abs((UB - LB) / LB) < MIP_GAP || time() - startT > timelimit || abs(UB - LB) < MIP_GAP || UB < LB)
            rootnode = false
            UB = typemax(Float64)
            LB = -LARGE_NUMBER
            RNtime = time() - startT
            println("ROOT NODE terminated in $(round(RNtime, digits=2))")
            undo_relax()
            if time() - startT > timelimit
                done = true
                MP_fractional = true
            end
        elseif doneSPs || abs((UB - LB) / LB) < MIP_GAP || time() - startT > timelimit || abs(UB - LB) < MIP_GAP || UB < LB
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
                println("OBJ MP: ", objMP, " in ", round(solvetime, digits=2), " seconds.")
                xVal = value.(MP[:x])               # first stage X vars
                zVal = value.(MP[:z])               # first stage Z vars
                LB = objMP > LB + 0.01 ? objMP : LB
            end
            sol.cost = FScosts
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
                    addSubproblemSolDAR!(m, subproblems[l][t][s], sol, xVal[l,t], l,t,s, R, Inst, Subpaths[l][t][s], subpath_road_networks, all_load_expanded_graphs, all_subpath_graphs)
                    objSP = objective_value(subproblems[l][t][s])
                else
                    objSP = typemax(Float64)
                    sol.cost += R.Pi[s]*objSP
                end
            end
        end
    end
    totT = time() - startT # stop algorithm time
    println("ALG TIME: ", totT)
    println("IT $it, LB = $LB, UB = $UB, IP(UB) = $(sol.cost), numCuts = ", length(cuts), " TOTAL T: ", totT)
    return LB, sol.cost, totT, it, length(cuts), num_second_stage_vars, RNtime, solveTMP, solveTSP, xVal, sol, Subpaths
end

"""
Trip object

### Attributes
* `line` - line id
* `freq` - freq id
* `scenario` - scenario id
* `subpaths` - set of subpaths forming the trip
* `actual_schedule` - Vector of tuples (OSM id, Time)
* `actual_route` - Vector of OSM id nodes forming the actual route
"""
mutable struct TripDAR
    line::Int
    freq::Int
    scenario::Int
    subpaths::Vector{SubPathDAR}
    num_passengers::Int
    actual_schedule::Vector{Tuple{Int, Time}}
    actual_route::Vector{Int}
end
TripDAR(l::Int, t::Int, s::Int) = TripDAR(l,t,s,Vector{SubPathDAR}(), 0, Vector{Tuple{Int, Time}}(), Vector{Int}())
"""
Object that captures the level of service of passengers served
"""
mutable struct ServedPassengerDAR
    aggregPax_index::Tuple{Int, Int, Int} # index of the corresponding AggregPassengerDAR object
    servedTrip::Tuple{Int,Int} # (l, t)
    originStopOSMid::Int
    pickupTime::Time
    pickupStopOSMid::Int
    destStopOSMid::Int
    dropoffTime::Time
    dropoffStopOSMid::Int
    numberPassengers::Int
    walkTimePickupSeconds::Int
    walkTimeDropoffSeconds::Int
    waitTimeSeconds::Int
    arrivalDelaySeconds::Int
    walkPathPickup::Vector{Int} # set of OSM node ids forming the path
    walkPathDropoff::Vector{Int} # set of OSM node ids forming the path
    inVehicleTimeSeconds::Int
    taxiInVehicleTimeSeconds::Int # direct travel time
end
"""
Solution object

### Attributes
* `cost` - objective value
* `solve_time` - solving time
* `trips` - set of operating trips
* `num_passengers` - number of passengers served
* `served_passengers` - level of service of passengers served
* `transitDistanceMeters` - distance traveled by reference routes
* `actualDistanceMeters` - distance traveled by deviation schedule
"""
mutable struct SolDAR
    cost::Float64
    solve_time::Float64
    trips::Vector{TripDAR}
    num_passengers::Int
    served_passengers::Vector{Vector{ServedPassengerDAR}} # one vector per demand scenario
    transitDistanceMeters::Float64
    actualDistanceMeters::Float64
end
SolDAR() = SolDAR(0.,0., Vector{TripDAR}(),0,Vector{Vector{ServedPassengerDAR}}(),0.,0.)
SolDAR(numS::Int) = SolDAR(0.,0., Vector{TripDAR}(),0,[Vector{ServedPassengerDAR}() for s in 1:numS],0.,0.)


"""
Build solution and compute level of service
"""
function buildSolDAR(
        m::MapData,
        model::JuMP.Model,
        R::RouteSettingDataDAR,
        Inst::InstanceSettingData,
        Subpaths::Vector{Vector{Vector{PoolDAR}}},
        subpath_road_networks,
        all_load_expanded_graphs,
        all_subpath_graphs,
        normalized::Bool=NORMALIZED,
    )
    @unpack numL, numS, TDsp, MaxWalk, λ, μ, σ, δ = Inst
    sol = SolDAR(numS)
    sol.cost = objective_value(model)
    sol.solve_time = solve_time(model)
    xVal = value.(model[:x])
    yVal = value.(model[:y])
    for l in 1:numL, t in eachindex(R.Lines[l].freq)
        if xVal[l,t] > EPS
            line_start_time = R.Lines[l].ref_stop_time[R.Lines[l].ref_stop_id[1]][t]
            for s in 1:numS
                trip = TripDAR(l,t,s)
                actualLineSchedule = Vector{Tuple{Int,Time}}()
                subpathRoutes = [Vector{Int}() for i in 1:length(R.Lines[l].ref_stop_id)-1]
                # store pickup info for the matching dropoff subpath
                pickup_times = Dict{Tuple{Int, Int, Int}, Dates.Time}()
                dropoff_times = Dict{Tuple{Int, Int, Int}, Dates.Time}()
                walk_times_pickup = Dict{Tuple{Int, Int, Int}, Int}()
                walk_times_dropoff = Dict{Tuple{Int, Int, Int}, Int}()
                wait_times = Dict{Tuple{Int, Int, Int}, Int}()
                walk_paths_pickup = Dict{Tuple{Int, Int, Int}, Vector{Int}}()
                walk_paths_dropoff = Dict{Tuple{Int, Int, Int}, Vector{Int}}()
                pickup_nodes = Dict{Tuple{Int, Int, Int}, Int}()
                dropoff_nodes = Dict{Tuple{Int, Int, Int}, Int}()
                arrival_delays = Dict{Tuple{Int, Int, Int}, Int}()
                num_passengers = Dict{Tuple{Int, Int, Int}, Int}()
                for i in eachindex(Subpaths[l][t][s].all)
                    if yVal[l,t,s,i] > EPS
                        trip.num_passengers += length(Subpaths[l][t][s].all[i].paxPU) #Subpaths[l][t][s].all[i].load
                        push!(trip.subpaths, Subpaths[l][t][s].all[i])

                        # calculate level of service of served passengers
                        subpathSetPassengers = union(Subpaths[l][t][s].all[i].paxPU, Subpaths[l][t][s].all[i].paxDO)
                        subpathTENpath = Subpaths[l][t][s].all[i].graph_path
                        subpathTSLNsrcNode = Subpaths[l][t][s].all[i].i
                        subpathTSLNsnkNode = Subpaths[l][t][s].all[i].j
                        (x1,t1,c1) = all_load_expanded_graphs[l][t].n2p[subpathTSLNsrcNode]
                        (x2,t2,c2) = all_load_expanded_graphs[l][t].n2p[subpathTSLNsnkNode]
                        # build actual route
                        if x1 != 0 && x2 != 0
                            if length(subpathTENpath) > EPS && length(subpathSetPassengers) > EPS
                                for n in subpathTENpath[2:end-2]
                                    (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                                    osm_node = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                    if length(subpathRoutes[x1]) > EPS
                                        if osm_node != subpathRoutes[x1][end]
                                            push!(subpathRoutes[x1], osm_node)
                                        end
                                    else
                                        push!(subpathRoutes[x1], osm_node)
                                    end
                                end
                            else
                                push!(subpathRoutes[x1], R.Lines[l].ref_stop_id[x1])
                            end
                        end
                        # compute total number of customers picked up and their saved time
                        for (idxPassenger, subPathPassenger) in enumerate(subpathSetPassengers)
                            (passengerOriginStop,timeRequestTDmp, passengerDestStop, numberOfPassengers) = subPathPassenger
                            pax = R.Demand[s][passengerOriginStop,timeRequestTDmp, passengerDestStop] #NOTE: remember is an AggregPassengerDAR object
                            # split between pickup and dropoff element
                            if numberOfPassengers > EPS
                                # compute walk to pikup spot, wait and pickup time
                                be_ready_time = pax.be_ready_times[l,t]
                                walkTimeSeconds = typemax(Int)
                                waitTimeSeconds = typemax(Int)
                                walkingPath = Int[]
                                pickupStopOSMid = -1
                                pickupTime = Time(0)
                                bestCombo = typemax(Int)
                                dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(R.Walk_graph, m.v[passengerOriginStop])
                                for n in subpathTENpath[2:end-1]
                                    (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                                    tempPickupTime = Dates.Second((time + t1)*TDsp) + line_start_time
                                    # compute walking path from passengerOriginStop to loc
                                    tempWalkingPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[subpath_road_networks[l][x1][x2-x1].B2A[loc]])
                                    currentWalkingDistance = 0.
                                    for (idx, node) in enumerate(tempWalkingPath[1:end-1])
                                        node2 = tempWalkingPath[idx+1]
                                        currentWalkingDistance += R.Walk_graph.weights[node, node2]
                                    end
                                    if currentWalkingDistance <= MaxWalk
                                        tempWalkingTimeSeconds = ceil(Int, (currentWalkingDistance / WALKSPEED_M_PER_SEC))
                                        tempWaitingTimeSeconds = ceil(Int, Dates.value(Dates.Second(tempPickupTime - be_ready_time))+tempWalkingTimeSeconds)
                                        normPickupTime = normalized ? Dates.value(Dates.Second(tempPickupTime))/(R.Taxi_travel_times[passengerOriginStop, passengerDestStop]) : Dates.value(Dates.Second(tempPickupTime))
                                        if  tempWaitingTimeSeconds > -0.01 && λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds - σ*normPickupTime < bestCombo - 0.01
                                            walkTimeSeconds = tempWalkingTimeSeconds
                                            waitTimeSeconds = tempWaitingTimeSeconds
                                            walkingPath = [m.n[ind] for ind in tempWalkingPath]
                                            pickupStopOSMid = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                            pickupTime = tempPickupTime
                                            bestCombo = λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds - σ*normPickupTime
                                        end
                                    end
                                end
                                pickup_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = pickupTime
                                walk_times_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkTimeSeconds
                                wait_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = waitTimeSeconds
                                walk_paths_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkingPath
                                pickup_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = pickupStopOSMid
                                num_passengers[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = numberOfPassengers
                            else
                                # compute walk to destination, delay
                                walkTimeSeconds = typemax(Int)
                                arrivalDelaySeconds = typemax(Int)
                                walkingPath = Int[]
                                dropoffStopOSMid = -1
                                dropoffTime = Time(0)
                                bestCombo = typemax(Int)
                                dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(R.Walk_graph, m.v[passengerDestStop])
                                for n in subpathTENpath[2:end-1]
                                    (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                                    tempDropoffTime = Dates.Second((time + t1)*TDsp) + line_start_time
                                    # compute walking path from passengerDestStop to loc (undirected)
                                    tempWalkingPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[subpath_road_networks[l][x1][x2-x1].B2A[loc]])
                                    currentWalkingDistance = 0.
                                    for (idx, node) in enumerate(tempWalkingPath[1:end-1])
                                        node2 = tempWalkingPath[idx+1]
                                        currentWalkingDistance += R.Walk_graph.weights[node, node2]
                                    end
                                    if currentWalkingDistance <= MaxWalk
                                        tempWalkingTimeSeconds = ceil(Int, (currentWalkingDistance / WALKSPEED_M_PER_SEC))
                                        tempArrivalTime = tempDropoffTime + Dates.Second(tempWalkingTimeSeconds) 
                                        tempDelaySeconds = Dates.value(Dates.Second(pax.dest.time - tempArrivalTime))
                                        normDropoffTime = normalized ? Dates.value(Dates.Second(tempDropoffTime))/(R.Taxi_travel_times[passengerOriginStop, passengerDestStop]) : Dates.value(Dates.Second(tempDropoffTime))
                                        tempWeightedDelay = tempDelaySeconds > EPS ? δ*tempDelaySeconds : abs((δ/2)*tempDelaySeconds)
                                        normWeightedDelay = normalized ? tempWeightedDelay/(R.Taxi_travel_times[passengerOriginStop, passengerDestStop]) : tempWeightedDelay

                                        if  λ*tempWalkingTimeSeconds + σ*normDropoffTime + normWeightedDelay < bestCombo - 0.01
                                            walkTimeSeconds = tempWalkingTimeSeconds
                                            walkingPath = [m.n[ind] for ind in tempWalkingPath]
                                            dropoffStopOSMid = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                            arrivalDelaySeconds = tempDelaySeconds
                                            dropoffTime = tempDropoffTime
                                            bestCombo = λ*tempWalkingTimeSeconds + σ*normDropoffTime + normWeightedDelay
                                        end
                                    end
                                end
                                dropoff_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = dropoffTime
                                walk_times_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkTimeSeconds
                                walk_paths_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkingPath
                                dropoff_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = dropoffStopOSMid
                                arrival_delays[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = arrivalDelaySeconds
                            end

                        end
                        # build line schedule (stop, time)
                        if length(subpathTENpath)> 1.1
                            for n in subpathTENpath[2:end-1]
                                (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                                x = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                push!(actualLineSchedule, (x, Dates.Second((time+t1)*TDsp)+line_start_time))
                            end
                        end
                    end
                end
                # create the served passengers info
                for (passengerOriginStop, timeRequestTDmp, passengerDestStop) in keys(pickup_times)
                
                    taxiInVehicleTimeSeconds = ceil(Int, R.Taxi_travel_times[passengerOriginStop, passengerDestStop])
                    pickupTime = pickup_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    dropoffTime = dropoff_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    inVehicleTimeSeconds = Dates.value(Dates.Second(dropoffTime - pickupTime))
                    pickupStopOSMid = pickup_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    dropoffStopOSMid = dropoff_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    walkTimePickupSeconds = walk_times_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    walkTimeDropoffSeconds = walk_times_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    walkPathPickup = walk_paths_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    walkPathDropoff = walk_paths_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    waitTimeSeconds = wait_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    arrivalDelaySeconds = arrival_delays[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    numberOfPassengers = num_passengers[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
                    newPassenger = ServedPassengerDAR((passengerOriginStop, timeRequestTDmp, passengerDestStop), (l, t), passengerOriginStop, pickupTime, pickupStopOSMid, passengerDestStop, dropoffTime, dropoffStopOSMid, numberOfPassengers, walkTimePickupSeconds, walkTimeDropoffSeconds, waitTimeSeconds, arrivalDelaySeconds, walkPathPickup, walkPathDropoff, inVehicleTimeSeconds, taxiInVehicleTimeSeconds)
                    push!(sol.served_passengers[s], newPassenger)
                end


                push!(subpathRoutes[end], R.Lines[l].ref_stop_id[end])
                actual_route = Vector{Int}()
                for subpath_route in subpathRoutes
                    append!(actual_route, subpath_route)
                end
                trip.actual_route = actual_route
                trip.actual_schedule = sort(actualLineSchedule, by= x->x[2])
                actual_route = [first(x) for x in trip.actual_schedule]
                push!(sol.trips, trip)
                sol.num_passengers += trip.num_passengers

                # compute travel distance
                for i in 1:length(R.Lines[l].transit_stop_id)-1
                    n1 = R.Lines[l].transit_stop_id[i]
                    n2 = R.Lines[l].transit_stop_id[i+1]
                    sol.transitDistanceMeters += OpenStreetMapX.distance(m.nodes[n1], m.nodes[n2])
                end
                for i in 1:length(actual_route)-1
                    n1 = actual_route[i]
                    n2 = actual_route[i+1]
                    sol.actualDistanceMeters += OpenStreetMapX.distance(m.nodes[n1], m.nodes[n2])
                end
            end
        end
    end
    return sol
end
function addSubproblemSolDAR!(
    m::MapData,
    model::JuMP.Model,
    sol::SolDAR,
    xVal,
    l::Int,
    t::Int,
    s::Int,
    R::RouteSettingDataDAR,
    Inst::InstanceSettingData,
    Subpaths::PoolDAR,
    subpath_road_networks,
    all_load_expanded_graphs,
    all_subpath_graphs,
    normalized::Bool=NORMALIZED,
    )
    @unpack numL, numS, TDsp, MaxWalk, λ, μ, σ, δ = Inst
    sol.cost += R.Pi[s]*objective_value(model)
    sol.solve_time += solve_time(model)
    yVal = value.(model[:y])
    if xVal > EPS
        line_start_time = R.Lines[l].ref_stop_time[R.Lines[l].ref_stop_id[1]][t]
        trip = TripDAR(l,t,s)
        actualLineSchedule = Vector{Tuple{Int,Time}}()
        subpathRoutes = [Vector{Int}() for i in 1:length(R.Lines[l].ref_stop_id)-1]
        # store pickup info for the matching dropoff subpath
        pickup_times = Dict{Tuple{Int, Int, Int}, Dates.Time}()
        dropoff_times = Dict{Tuple{Int, Int, Int}, Dates.Time}()
        walk_times_pickup = Dict{Tuple{Int, Int, Int}, Int}()
        walk_times_dropoff = Dict{Tuple{Int, Int, Int}, Int}()
        wait_times = Dict{Tuple{Int, Int, Int}, Int}()
        walk_paths_pickup = Dict{Tuple{Int, Int, Int}, Vector{Int}}()
        walk_paths_dropoff = Dict{Tuple{Int, Int, Int}, Vector{Int}}()
        pickup_nodes = Dict{Tuple{Int, Int, Int}, Int}()
        dropoff_nodes = Dict{Tuple{Int, Int, Int}, Int}()
        arrival_delays = Dict{Tuple{Int, Int, Int}, Int}()
        num_passengers = Dict{Tuple{Int, Int, Int}, Int}()
        for i in eachindex(Subpaths.all)
            if yVal[i] > EPS
                trip.num_passengers += length(Subpaths.all[i].paxPU)
                push!(trip.subpaths, Subpaths.all[i])

                # calculate level of service of served passengers
                subpathSetPassengers = union(Subpaths.all[i].paxPU, Subpaths.all[i].paxDO)
                subpathTENpath = Subpaths.all[i].graph_path
                subpathTSLNsrcNode = Subpaths.all[i].i
                subpathTSLNsnkNode = Subpaths.all[i].j
                (x1,t1,c1) = all_load_expanded_graphs[l][t].n2p[subpathTSLNsrcNode]
                (x2,t2,c2) = all_load_expanded_graphs[l][t].n2p[subpathTSLNsnkNode]
                # build actual route
                if x1 != 0 && x2 != 0
                    if length(subpathTENpath) > EPS && length(subpathSetPassengers) > EPS
                        for n in subpathTENpath[2:end-2]
                            (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                            osm_node = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                            if length(subpathRoutes[x1]) > EPS
                                if osm_node != subpathRoutes[x1][end]
                                    push!(subpathRoutes[x1], osm_node)
                                end
                            else
                                push!(subpathRoutes[x1], osm_node)
                            end
                        end
                    else
                        push!(subpathRoutes[x1], R.Lines[l].ref_stop_id[x1])
                    end
                end
                # compute total number of customers picked up and their saved time
                for (idxPassenger, subPathPassenger) in enumerate(subpathSetPassengers)
                    (passengerOriginStop,timeRequestTDmp, passengerDestStop, numberOfPassengers) = subPathPassenger
                    pax = R.Demand[s][passengerOriginStop,timeRequestTDmp, passengerDestStop] #NOTE: remember is an AggregPassengerDAR object
                    # split between pickup and dropoff element
                    if numberOfPassengers > EPS
                        # compute walk to pikup spot, wait and pickup time
                        be_ready_time = pax.be_ready_times[l,t]
                        walkTimeSeconds = typemax(Int)
                        waitTimeSeconds = typemax(Int)
                        walkingPath = Int[]
                        pickupStopOSMid = -1
                        pickupTime = Time(0)
                        bestCombo = typemax(Int)
                        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(R.Walk_graph, m.v[passengerOriginStop])
                        for n in subpathTENpath[2:end-1]
                            (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                            tempPickupTime = Dates.Second((time + t1)*TDsp) + line_start_time
                            # compute walking path from passengerOriginStop to loc
                            tempWalkingPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[subpath_road_networks[l][x1][x2-x1].B2A[loc]])
                            currentWalkingDistance = 0.
                            for (idx, node) in enumerate(tempWalkingPath[1:end-1])
                                node2 = tempWalkingPath[idx+1]
                                currentWalkingDistance += R.Walk_graph.weights[node, node2]
                            end
                            if currentWalkingDistance <= MaxWalk
                                tempWalkingTimeSeconds = ceil(Int, (currentWalkingDistance / WALKSPEED_M_PER_SEC))
                                tempWaitingTimeSeconds = ceil(Int, Dates.value(Dates.Second(tempPickupTime - be_ready_time))+tempWalkingTimeSeconds)
                                normPickupTime = normalized ? Dates.value(Dates.Second(tempPickupTime))/(R.Taxi_travel_times[passengerOriginStop, passengerDestStop]) : Dates.value(Dates.Second(tempPickupTime))
                                if  tempWaitingTimeSeconds > -0.01 && λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds - σ*normPickupTime < bestCombo - 0.01
                                    walkTimeSeconds = tempWalkingTimeSeconds
                                    waitTimeSeconds = tempWaitingTimeSeconds
                                    walkingPath = [m.n[ind] for ind in tempWalkingPath]
                                    pickupStopOSMid = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                    pickupTime = tempPickupTime
                                    bestCombo = λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds - σ*normPickupTime
                                end
                            end
                        end
                        pickup_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = pickupTime
                        walk_times_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkTimeSeconds
                        wait_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = waitTimeSeconds
                        walk_paths_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkingPath
                        pickup_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = pickupStopOSMid
                        num_passengers[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = numberOfPassengers
                    else
                        # compute walk to destination, delay
                        walkTimeSeconds = typemax(Int)
                        arrivalDelaySeconds = typemax(Int)
                        walkingPath = Int[]
                        dropoffStopOSMid = -1
                        dropoffTime = Time(0)
                        bestCombo = typemax(Int)
                        dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(R.Walk_graph, m.v[passengerDestStop])
                        for n in subpathTENpath[2:end-1]
                            (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                            tempDropoffTime = Dates.Second((time + t1)*TDsp) + line_start_time
                            # compute walking path from passengerDestStop to loc (undirected)
                            tempWalkingPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[subpath_road_networks[l][x1][x2-x1].B2A[loc]])
                            currentWalkingDistance = 0.
                            for (idx, node) in enumerate(tempWalkingPath[1:end-1])
                                node2 = tempWalkingPath[idx+1]
                                currentWalkingDistance += R.Walk_graph.weights[node, node2]
                            end
                            if currentWalkingDistance <= MaxWalk
                                tempWalkingTimeSeconds = ceil(Int, (currentWalkingDistance / WALKSPEED_M_PER_SEC))
                                tempArrivalTime = tempDropoffTime + Dates.Second(tempWalkingTimeSeconds) 
                                tempDelaySeconds = Dates.value(Dates.Second(pax.dest.time - tempArrivalTime))
                                normDropoffTime = normalized ? Dates.value(Dates.Second(tempDropoffTime))/(R.Taxi_travel_times[passengerOriginStop, passengerDestStop]) : Dates.value(Dates.Second(tempDropoffTime))
                                tempWeightedDelay = tempDelaySeconds > EPS ? δ*tempDelaySeconds : abs((δ/2)*tempDelaySeconds)
                                normWeightedDelay = normalized ? tempWeightedDelay/(R.Taxi_travel_times[passengerOriginStop, passengerDestStop]) : tempWeightedDelay

                                if  λ*tempWalkingTimeSeconds + σ*normDropoffTime + normWeightedDelay < bestCombo - 0.01
                                    walkTimeSeconds = tempWalkingTimeSeconds
                                    walkingPath = [m.n[ind] for ind in tempWalkingPath]
                                    dropoffStopOSMid = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                    arrivalDelaySeconds = tempDelaySeconds
                                    dropoffTime = tempDropoffTime
                                    bestCombo = λ*tempWalkingTimeSeconds + σ*normDropoffTime + normWeightedDelay
                                end
                            end
                        end
                        dropoff_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = dropoffTime
                        walk_times_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkTimeSeconds
                        walk_paths_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = walkingPath
                        dropoff_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = dropoffStopOSMid
                        arrival_delays[(passengerOriginStop, timeRequestTDmp, passengerDestStop)] = arrivalDelaySeconds
                    end

                end
                # build line schedule (stop, time)
                if length(subpathTENpath)> 1.1
                    for n in subpathTENpath[2:end-1]
                        (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                        x = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                        push!(actualLineSchedule, (x, Dates.Second((time+t1)*TDsp)+line_start_time))
                    end
                end
            end
        end

        # create the served passengers info
        for (passengerOriginStop, timeRequestTDmp, passengerDestStop) in keys(pickup_times)
                
            taxiInVehicleTimeSeconds = ceil(Int, R.Taxi_travel_times[passengerOriginStop, passengerDestStop])
            pickupTime = pickup_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            dropoffTime = dropoff_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            inVehicleTimeSeconds = Dates.value(Dates.Second(dropoffTime - pickupTime))
            pickupStopOSMid = pickup_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            dropoffStopOSMid = dropoff_nodes[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            walkTimePickupSeconds = walk_times_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            walkTimeDropoffSeconds = walk_times_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            walkPathPickup = walk_paths_pickup[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            walkPathDropoff = walk_paths_dropoff[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            waitTimeSeconds = wait_times[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            arrivalDelaySeconds = arrival_delays[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            numberOfPassengers = num_passengers[(passengerOriginStop, timeRequestTDmp, passengerDestStop)]
            newPassenger = ServedPassengerDAR((passengerOriginStop, timeRequestTDmp, passengerDestStop), (l, t), passengerOriginStop, pickupTime, pickupStopOSMid, passengerDestStop, dropoffTime, dropoffStopOSMid, numberOfPassengers, walkTimePickupSeconds, walkTimeDropoffSeconds, waitTimeSeconds, arrivalDelaySeconds, walkPathPickup, walkPathDropoff, inVehicleTimeSeconds, taxiInVehicleTimeSeconds)
            push!(sol.served_passengers[s], newPassenger)
        end

        push!(subpathRoutes[end], R.Lines[l].ref_stop_id[end])
        actual_route = Vector{Int}()
        for subpath_route in subpathRoutes
            append!(actual_route, subpath_route)
        end
        trip.actual_route = actual_route
        trip.actual_schedule = sort(actualLineSchedule, by= x->x[2])
        actual_route = [first(x) for x in trip.actual_schedule]
        push!(sol.trips, trip)
        sol.num_passengers += trip.num_passengers

        # compute travel distance
        for i in 1:length(R.Lines[l].transit_stop_id)-1
            n1 = R.Lines[l].transit_stop_id[i]
            n2 = R.Lines[l].transit_stop_id[i+1]
            sol.transitDistanceMeters += OpenStreetMapX.distance(m.nodes[n1], m.nodes[n2])
        end
        for i in 1:length(actual_route)-1
            n1 = actual_route[i]
            n2 = actual_route[i+1]
            sol.actualDistanceMeters += OpenStreetMapX.distance(m.nodes[n1], m.nodes[n2])
        end
    end
    return nothing
end

"""
Summarise solution. Create file if DNE

### Keywords
* `m` - MapData (Manhattan)
* `inst` - InstanceSettingData
* `R` - RouteSettingDataDAR
* `sol` - Sol
* `directDist` - Dictionary (node ID => direct taxi trip distance in meters)
* `los_fn` - filename to record level of service
* `trial_id` - trial ID
* `lineDict` - dictionary of line index to line ID
### Writes
* one line of level of service metrics to `los_fn` per (l, t, s)
"""
function summarise_losDAR(
        m::MapData,
        inst::InstanceSettingData,
        R::RouteSettingDataDAR,
        sol::SolDAR,
        directDist::Dict{Tuple{Int,Int},Float64},
        los_fn::String,
        trial_id::Int,
        lineDict::Dict{Int,Int}
    )
    #--- sort served passengers and trips
    selectedTrips = unique([(trip.line, trip.freq, trip.scenario) for trip in sol.trips])
    servedPassengers = Dict((l, t, s) => ServedPassengerDAR[] for (l, t, s) in selectedTrips)
    trips = Dict((trip.line, trip.freq, trip.scenario) => trip for trip in sol.trips)
    for s in 1:inst.numS, pax in sol.served_passengers[s]
        (l, t) = pax.servedTrip
        push!(servedPassengers[(l, t, s)], pax)
    end

    #--- precompute taxi distances of rejected passengers
    # all individual passenger pickups
    allPax = Dict(s => [(pax, aggPax.orig.stop, aggPax.dest.stop) for aggPax in values(R.Demand[s]) for pax in aggPax.pax] for s in 1:inst.numS)
    passengerPickups = Dict(s => [] for s in 1:inst.numS)
    for trip in sol.trips
        for subpath in trip.subpaths
            for pax in subpath.paxPU
                aggPax = R.Demand[trip.scenario][(pax[1], pax[2], pax[3])]
                for p in aggPax.pax

                    push!(passengerPickups[trip.scenario], p)
                end
            end
        end
    end

    # compute taxi distances of passengers who weren't picked up by scenario
    taxiDist = Dict(s => 0. for s in 1:inst.numS)
    for (s, paxList) in allPax
        for (pax, originStop, destinStop) in paxList
            # served passenger
            if pax in passengerPickups[s]
                continue
            end
            # direct distance
            taxiDist[s] += directDist[(originStop, destinStop)]
        end
    end

    # write the first line
    f =  open(los_fn, "w")
    write(f, "trial_id,l,t,s,num_pickups,avg_detour,max_detour,avg_delay_sec,avg_delay_norm,max_delay_sec,max_delay_norm,avg_walkPU_sec,max_walkPU_sec,avg_walkDO_sec,max_walkDO_sec,avg_wait_sec,max_wait_sec,dist_m,taxi_dist_m\n")

    # aggregate level of service metrics
    for ((l, t, s), setPax) in servedPassengers
        # coverage
        if isempty(setPax)
            continue
        end
        num_pickups = sum(pax.numberPassengers for pax in setPax)

        # detour
        avg_detour = sum(pax.numberPassengers * pax.inVehicleTimeSeconds / pax.taxiInVehicleTimeSeconds for pax in setPax) / num_pickups
        max_detour = maximum([pax.inVehicleTimeSeconds / pax.taxiInVehicleTimeSeconds for pax in setPax])

        # displacement
        avg_delay_sec = sum(pax.numberPassengers * pax.arrivalDelaySeconds for pax in setPax) / num_pickups
        avg_delay_norm = sum(pax.numberPassengers * pax.arrivalDelaySeconds / pax.taxiInVehicleTimeSeconds for pax in setPax) / num_pickups
        max_delay_sec = maximum([pax.arrivalDelaySeconds for pax in setPax])
        max_delay_norm = maximum([pax.arrivalDelaySeconds / pax.taxiInVehicleTimeSeconds for pax in setPax])

        # walk PU
        avg_walkPU_sec = sum(pax.numberPassengers * pax.walkTimePickupSeconds for pax in setPax) / num_pickups
        max_walkPU_sec = maximum([pax.walkTimePickupSeconds for pax in setPax])
        # walk DO
        avg_walkDO_sec = sum(pax.numberPassengers * pax.walkTimeDropoffSeconds for pax in setPax) / num_pickups
        max_walkDO_sec = maximum([pax.walkTimeDropoffSeconds for pax in setPax])

        # wait
        avg_wait_sec = sum(pax.numberPassengers * pax.waitTimeSeconds for pax in setPax) / num_pickups
        max_wait_sec = maximum([pax.waitTimeSeconds for pax in setPax])

        # distance
        dist_m = 0. 
        trip = trips[(l, t, s)]
        dist_m = 0.
        for i in 1:length(trip.actual_route)-1
            n1 = trip.actual_route[i]
            n2 = trip.actual_route[i+1]
            dist_m += OpenStreetMapX.distance(m.nodes[n1], m.nodes[n2])
        end
        metrics = [trial_id, lineDict[l], t, s, num_pickups, avg_detour, max_detour, avg_delay_sec, avg_delay_norm, max_delay_sec, max_delay_norm, avg_walkPU_sec, max_walkPU_sec, avg_walkDO_sec, max_walkDO_sec, avg_wait_sec, max_wait_sec, dist_m, taxiDist[s]]
        write(f, join(metrics, ",") * "\n")
    end
    close(f)
    return nothing
end


"""
Custom study case for the MiND-DAR Model

    * 130 stops (13 streets W2E between 34th and 58th street, and 10 avenues)
    * 13 lines (one per street)
    
    TODO:
    - Read stops
    - Plot stops
    - Read lines
    - Plot lines
    - Generate demand
    - Plot demand
    - Run single-line test
    - Run whole algorithm
    - Compare transit with MT

"""

