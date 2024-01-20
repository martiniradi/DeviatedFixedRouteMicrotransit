"""
First-stage reference route

### Attributes
* `all_stops` - set of MapData nodes (stops) that the line could serve demand for
* `ref_stop_id` - set of reference stops (start/end of subpath)
* `freq` - set of frequencies for the line (times that each trip starts operating)
* `ref_stop_time` - Dict (stop ID => vector of allowed stopping times)
* `arrival_time` - time of arrival of the line to Laguardia at each frequency
* `capacity` - capacity of line
* `cost` - cost to operate line
* `ref_route_nodes` - vector of all OSM nodes forming the reference route (all intersections, stops and non-stops)
* `all_route_nodes` - vector of all OSM nodes within the maximum deviation from the reference route (all intersections, stops and non-stops)
* `transit_stop_id` - set of stops of the reference route when running as transit line
* `transit_stop_time` - Dict (stop ID => vector of scheduled stopping times)
* `ref_dist` - travel distance of reference route (in meters)
* `exit_point_idx` - idx corresponding to the line exit point (from EXIT_POINTS_OSM)
"""
struct Line
    all_stops::Vector{Int}
    freq::Vector{Time}
    ref_stop_id::Vector{Int}
    ref_stop_time::Dict{Int, Vector{Time}}
    arrival_time::Vector{Time}
    capacity::Int64
    cost::Float64
    ref_route_nodes::Vector{Int}
    all_route_nodes::Vector{Int}
    transit_stop_id::Vector{Int}
    transit_stop_time::Dict{Int, Vector{Time}}
    ref_dist::Float64
    exit_point_idx::Int
end
"""
read CSV files about lines and creates the vector of Line objects

### Keywords
* `m` - MapData
* `operating_horizon_seconds` - time between first and last operating vehicle operating (starting) a line
* `numL` - number of lines
* `refNum` - one reference stop every `refNum` line stops
* `Capacity` - maximum number of passengers per vehicle
* `speedFactor` - how slower our ref vehicle travels compared to the travel times in MapData (e.g., 2 = travel time is twice (half the speed))
* `MaxDev` - maximum vehicle deviation from reference route (in meters)
* `exitPoints` - OSM id of four exit points
* `exitTimesSeconds` - time from each exit point to Laguardia (in seconds)
* `frequency_period_seconds` - duration between timeslots (in seconds)
* `TD` - time discretization (in seconds)
* `start_operations` - time of the day when the first vehicle departs
### Returns
* the vector of Line objects
"""
function readLines(m::MapData, operating_horizon_seconds::Int, numL::Int, refNum::Int, Capacity::Int, speedFactor::Float64, MaxDev::Float64, exitPoints::Vector{Int}, exitTimesSeconds::Vector{Float64}, frequency_period_seconds::Int, maximum_delay::Int, start_operations::Dates.Time)

    dwell_time = 30 # stop duration in transit stops ( in seconds)

    T = ceil(Int, operating_horizon_seconds/frequency_period_seconds) # number of frequencies
    horizon_end = start_operations + Dates.Second(operating_horizon_seconds + maximum_delay)
    data_dir = get_data_path()
    SetOfValidLines = Vector{Int}()
    dfLineSets = CSV.read(joinpath(data_dir,"input", "MiND-VRP_setup", "refroutes", "line_sets.csv"), DataFrame)
    for lineID in eachrow(dfLineSets)
        if lineID.set_id == "popular"
            push!(SetOfValidLines, lineID.line_id )
        end
    end
    # println("numL ", numL)
    selectedLines = SetOfValidLines[1:numL] # we only want numL lines
    # println("Set of valid lines: ", selectedLines)

    S = generateStops(m) # generate all candidate stops
    for ex in exitPoints # add exit point nodes as stops
        if !(ex in S)
            push!(S, ex)
        end
    end

    df = CSV.read(joinpath(data_dir, "input", "MiND-VRP_setup",  "refroutes", "reflines.csv"), DataFrame)
    lines = Vector{Line}()

    i = 0       # line index counter
    line_ref_stops = [Vector{Int}() for i in 1:numL]
    line_all_stops = [Vector{Int}() for i in 1:numL]
    for r in eachrow(df)
        if r.line_id in selectedLines
            if r.seq == 1
                i += 1
            end
            if !(r.osm_id in S)
                # it is an exit point, find the closest
                dist = typemax(Float64)
                exit_idx = -1
                for idx in eachindex(exitPoints)
                    if OpenStreetMapX.distance(m.nodes[r.osm_id], m.nodes[exitPoints[idx]]) < dist -0.01
                        dist = OpenStreetMapX.distance(m.nodes[r.osm_id], m.nodes[exitPoints[idx]])
                        exit_idx = idx
                    end
                end
                # push!(line_all_stops[i], exitPoints[exit_idx])
                # if (r.seq - 1) % refNum == 0
                #     push!(line_ref_stops[i], exitPoints[exit_idx])
                # end
            else
                push!(line_all_stops[i], r.osm_id)
                if (r.seq - 1) % refNum == 0
                    push!(line_ref_stops[i], r.osm_id)
                end
            end
        end
    end
    # println("TOTAL NUMBER OF STOPS: ", length(S))
    # print("Creating lines...")
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

        # find exit point of line and compute arrival times to laguardia
        exitIdx = nearest_node_index(m, line_all_stops[idx][end], exitPoints)
        line_duration_seconds = ceil(Int, ref_times[end] + exitTimesSeconds[exitIdx])
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

        push!(lines, Line(all_stops, line_frequencies, ref_stop_id, ref_stop_time, arrival_times, Capacity, cost, ref_route_nodes, all_route_nodes, line_all_stops[idx], transit_stop_time, ref_dist, exitIdx))
    end
    # println(" DONE")
    return lines, S, T
end

"""
Information about a single passenger:
    - Index
    - Origin coordinates (lat, lon)
    - Node and stop associated to passenger
    - Request dropoff time (DO by taxi)
"""
struct Passenger
    id::Int         # idx
    lat::Float64    # origin latitude
    lon::Float64    # origin longitude,
    node::Int64     # closest OSM node to origin
    stop::Int64     # closest OSM node that is a stop
    time::Time      # request drop-off time
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
function readDemand(
        m::MapData,
        demand_horizon_seconds::Int,
        all_stops::Vector{Int},
        scenario::Int,
        trainSet::Int,
        start_time::Dates.Time
    )
    df = CSV.read(joinpath(get_data_path(),  "input", "MiND-VRP_setup", "weekday_data_all", "taxiData_manhattan_all_$(scenario)_$(trainSet).csv"), DataFrame)

    # build list of Passengers from DataFrame
    D = Vector{Passenger}()
    for row in eachrow(df)
        # drop-off time (rounded to the nearest minute)
        dropoffTime = Dates.Time(round(Dates.DateTime(row.dtime, "yyyy-mm-dd HH:MM:SS"), Dates.Minute))
        # accept if drop-off time is within the planning horizon
        if dropoffTime <= start_time + Dates.Second(demand_horizon_seconds)
            node = point_to_nodes((row.plat, row.plong), m)
            stop = nearest_node(m, row.plat, row.plong, all_stops)
            push!(D, Passenger(row.id, row.plat, row.plong, node, stop, dropoffTime))
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
    * Taxi pick up time
    * Be-ready time depending on trip
    * Time id for first stage clustering
    * Set of potential trips that can serve this passenger
"""
mutable struct AggregPassenger
    id::Int
    pax::Vector{Int}                                # set of clustered passenger ids
    num::Int                                        # number of passengers
    origin_stop::Int                                # stop of aggregated demand
    request_dropoff_time::Time                  # requested dropoff time
    taxi_pickup_time::Time                      # taxi pickup time from origin
    be_ready_times::Dict{Tuple{Int,Int}, Time}  # be ready times for each possible trip (line, frequency t)
    trip_delays::Dict{Tuple{Int,Int}, Float64}  # weighted arrival delay by trip
    first_stage_time_instant::Int               # time instant used in master/first-stage problem
    candidateTrips::Vector{Tuple{Int,Int}}      # set of candidate trips (line, frequency) that can serve this passenger
end
"""
aggregate demand based on the first-stage time discretization and create AggregPassenger object

### Keywords
* `D` - set of passenger objects
* `TDmp` - first-stage time discretization (in seconds)
* `taxi_travel_times` - dictionary of taxi trips durations (in seconds) from any stop in Manhattan
* `set` - training (0) or test (1)
* `start_time` - starting time of the day
### Returns
* the dictionary of AggregPassenger objects
"""
function aggregateDemand(D::Vector{Passenger}, TDmp::Int, taxi_travel_times::Dict{Int, Float64})
    AggregDem = Dict{Tuple{Int, Int}, AggregPassenger}()
    counter = 0
    for p in D
        time_instant = DAY_START + ceil(Dates.Second(p.time - DAY_START), Dates.Second(TDmp))
        first_stage_time_instant = ceil(Int,Dates.value(Dates.Second(p.time - DAY_START))/TDmp)
        if haskey(AggregDem, (p.stop, first_stage_time_instant))
            AggregDem[(p.stop, first_stage_time_instant)].num += 1
            push!(AggregDem[(p.stop, first_stage_time_instant)].pax, p.id)
        else
            counter += 1
            taxi_pickup_time = p.time - Dates.Second(round(Int, taxi_travel_times[p.stop]))
            AggregDem[(p.stop, first_stage_time_instant)] = AggregPassenger(counter, Int[p.id], 1, p.stop, time_instant, taxi_pickup_time, Dict{Tuple{Int,Int},Time}(), Dict{Tuple{Int,Int},Float64}(), first_stage_time_instant, Vector{Tuple{Int, Int}}())
        end
    end
    return AggregDem 
end
"""
find candidate trips (line, frequency) for each AggregPassenger
"""
function findCandidateTrips!(m::MapData, D::Dict{Tuple{Int,Int}, AggregPassenger}, L::Vector{Line}, MaxArrivalScheduleDev::Int, MaxWalkDist::Float64, Taxi_travel_times::Dict{Int, Float64}, δ::Float64, normalized::Bool)
    # loop through demand
    for (key, d) in collect(D)
        # loop through lines
        for (l, line) in enumerate(L)
            # first check if the line is a candidate
            potential_stop = false 
            if d.origin_stop in line.all_stops
                potential_stop = true
            else
                range_nodes = getLocRangeNodes(m, m.nodes[d.origin_stop], MaxWalkDist, Set(line.all_stops))
                if isempty(range_nodes) == false
                    potential_stop = true
                end
            end
            if potential_stop
                # second check which frequencies of the line are candidate
                for freq in eachindex(line.arrival_time)
                    if line.arrival_time[freq] - Dates.Second(MaxArrivalScheduleDev) <= d.request_dropoff_time <= line.arrival_time[freq] + Dates.Second(MaxArrivalScheduleDev)
                        push!(d.candidateTrips, (l, freq))
                        # compute as well the be-ready time
                        #   find nearest reference stop
                        nearest_ref_stop = nearest_node(m, m.nodes[d.origin_stop], line.ref_stop_id)
                        ref_stop_time = line.ref_stop_time[nearest_ref_stop][freq]
                        # be ready time as the ref stop time - maximum walking time
                        be_ready_time = ref_stop_time - Dates.Second(ceil(Int, MaxWalkDist/1.4))
                        d.be_ready_times[(l, freq)] = be_ready_time
                        # compute arrival delay
                        arrival_dev = Dates.value(Dates.Second(line.arrival_time[freq] - d.request_dropoff_time))
                        if normalized
                            arrival_dev /= Taxi_travel_times[d.origin_stop]
                        end
                        total_weighted_delay = 0.
                        if arrival_dev > -0.01
                            total_weighted_delay += d.num*δ*arrival_dev
                        else
                            total_weighted_delay += d.num*(δ/2)*abs(arrival_dev)
                        end
                        d.trip_delays[(l,freq)] = total_weighted_delay
                    end
                end
            end
        end
    end
end

"""
compute overlapping frequencies for fleet size on each line
"""
function computeOverlapFrequencies(lines::Vector{Line}, taxi_travel_times::Dict{Int, Float64}, demand_horizon::Int, Num_freq::Int, maximum_delay::Int)
    start_operations = DAY_START - Dates.Second(OPERATIONAL_HORIZON)
    finish_operations = DAY_START + Dates.Second(demand_horizon + maximum_delay)
    # define set of frequencies
    set_frequencies = collect(start_operations:Dates.Second(TIME_PERIOD_SECONDS):finish_operations)
    freq_overlap = [Vector{Int}() for l in eachindex(lines), t in 1:Num_freq]
    for (l, line) in enumerate(lines)
        # calculate total line duration + travel back to start (the way back is taken as taxi_travel_time)
        for (t,line_start) in enumerate(line.freq)
            line_end = line.arrival_time[t]
            back_at_start = line_end + Dates.Second(ceil(Int, taxi_travel_times[line.ref_stop_id[1]]))
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
Instance parameter inputs

### Attributes
* `numL` - number of lines
* `numT` - number of frequencies
* `numS` - number of scenarios
* `K` - number of reference stops that sub-path can skip + 1
* `opT` - How earlier than the demand horizon, the first vehicle starts operating
* `TDmp` - time discretization at MP level
* `TDsp` - time discretization at SP level
* `MaxDev` - max deviation in meters from the reference route
* `Capacity` - vehicle capacity in number of passengers
* `MaxWalk` - maximum walking distance for passengers to a stop (in meters)
* `MaxWait` - maximum waiting time for passengers (in seconds)
* `MaxSchDev` - maximum schedule deviation fir passengers (used to compute maximum delay/earliness at destination)
* `Fleet` - number of vehicles available
* `λ` - importance of walking time for passengers
* `μ` - importance of waiting time for passengers
* `σ` - importance of in-vehicle time for passengers
* `δ` - importance of arrival delay for passengers (half for arrival earliness)
* `routeSpeedFactor` - relative speed of vehicle respetive to ones given by getSpeeds() (e.g., 1.2 is 20% slower)
* `refSpeedFactor` - relative speed of vehicle for reference schedule by getSpeeds() (e.g., 1.2 is 20% slower)
"""
struct InstanceSettingData
    numL::Int                   # number of lines
    demT::Int                   # demand horizon (in seconds)
    numS::Int                   # number of scenarios
    trainData::Bool             # type of test (in-sample = true, out-of-sample = false)
    scenarios::Vector{Int}      # set of scenario indices
    K::Int                      # number of reference stops that sub-path can skip + 1
    opT::Int                    # operational horizon (the first vehicle starts operating X seconds earlier 6 am (start of demT))
    RefStopDens::Int            # density of reference stops in the line. One every X transit stops
    TDmp::Int                   # time discretization at MP level
    TDsp::Int                   # time discretization at SP level
    MaxDev::Float64             # max deviation in meters from the reference route
    Capacity::Int               # vehicle capacity in number of passengers
    MaxWalk::Float64            # maximum walking distance for passengers to a stop (in meters)
    MaxWait::Int                # maximum waiting time for passengers (in seconds)
    MaxSchDev::Int              # maximum schedule deviation fir passengers (TODO: better define this term)
    Fleet::Int                  # number of vehicles available
    Theta::Float64              # first-stage load factor
    M::Float64                  # importance of serving demand
    λ::Float64                  # importance of walking time for passengers
    μ::Float64                  # importance of waiting time for passengers
    σ::Float64                  # importance of in-vehicle time for passengers
    δ::Float64                  # importance of arrival delay for passengers
    routeSpeedFactor::Float64   # relative speed of vehicle respetive to ones given by getSpeeds() (e.g., 1.2 is 20% slower)
    refSpeedFactor::Float64     # relative speed of vehicle for reference schedule by getSpeeds() (e.g., 1.2 is 20% slower)
end
InstanceSettingData(l::Int, t::Int, s::Int, k::Int, trainData::Bool) = InstanceSettingData(l,t,s,trainData,collect(1:s),k,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,MAX_DEVIATION_METERS,CAPACITY,MAX_WALKING_METERS,MAX_WAITING_SECONDS,MAX_SCHEDULE_DEV,FLEET_SIZE, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
InstanceSettingData(l::Int, t::Int, S::Vector{Int}, k::Int, trainData::Bool) = InstanceSettingData(l,t,length(S),trainData,S,k,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,MAX_DEVIATION_METERS,CAPACITY,MAX_WALKING_METERS,MAX_WAITING_SECONDS,MAX_SCHEDULE_DEV,FLEET_SIZE, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
"""
Model inputs

### Attributes
* `Lines` - set of Line objects
* `All_stops` - set of all potential stops
* `Num_freq` - number of frequencies in planning horizon
* `Freq_Overlap` - set of trip frequencies that overlap with a given line and time slot
* `Demand` - dictionary of demand (AggregPassenger), one per scenario
* `Single_requests` - vector of single passenger requests (before aggregating)
* `Trip_passengers` - Set of passengers ids that can be served by a trip (l,t)
* `Pi` - Dict (scenario => probability of scenario)
* `Walk_graph` - undirected graph version from the one in m (MapData), m.g
* `Taxi_travel_times` - Direct travel time from any stop (in seconds)
"""
struct RouteSettingData
    Lines::Vector{Line}
    All_stops::Vector{Int64}
    Num_freq::Int64
    Freq_Overlap::Array{Vector{Int},2}
    Demand::Vector{Dict{Tuple{Int64,Int64},AggregPassenger}}
    Trip_passengers::Array{Vector{Tuple{Int, Int}},3}
    Pi::Vector{Float64}
    Walk_graph::SimpleWeightedGraph
    Taxi_travel_times::Dict{Int, Float64}
end
"""
create the RouteSettingData object with the main input for the algorithm

### Keywords
* `Inst` - InstanceSettingData
### Returns
* the RouteSettingData object
"""
function buildRouteSettingData(
        Inst::InstanceSettingData,
    )
    @unpack numL, demT, numS, trainData, scenarios, opT, K, RefStopDens, TDmp, TDsp, MaxDev, Capacity, MaxWalk, MaxWait, MaxSchDev, λ, μ, δ, routeSpeedFactor, refSpeedFactor = Inst

    m = generateNetwork()                           # read OSM data
    walk_graph = generateUndirectedGraph(m)         # create undirected version of m.g (for computing walking paths)
    lines, all_stops, num_freq = readLines(m, demT+opT, numL, RefStopDens, Capacity, refSpeedFactor, MaxDev, EXIT_POINTS_OSM, EXIT_TIMES_SECONDS, TIME_PERIOD_SECONDS,MaxSchDev, DAY_START - Dates.Second(opT))
    taxi_travel_times = readTravelTimes()
    all_demand = Vector{Dict{Tuple{Int64,Int64},AggregPassenger}}()
    trip_passengers = [Vector{Tuple{Int, Int}}() for l in 1:numL, t in 1:num_freq, s in 1:numS]
    data_set = trainData ? 0 : 1
    for (s, scenario) in enumerate(scenarios)
        passenger_demand = readDemand(m, demT, all_stops, scenario, data_set, DAY_START)
        aggregated_demand = aggregateDemand(passenger_demand,TDmp, taxi_travel_times)
        findCandidateTrips!(m, aggregated_demand, lines, MaxSchDev, MaxWalk, taxi_travel_times, δ, NORMALIZED)
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

    freq_overlap = computeOverlapFrequencies(lines, taxi_travel_times, demT, num_freq, MaxSchDev)

    return RouteSettingData(lines, all_stops, num_freq, freq_overlap, all_demand, trip_passengers, Pi, walk_graph, taxi_travel_times), m
end
