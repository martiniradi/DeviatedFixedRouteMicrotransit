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
mutable struct Trip
    line::Int
    freq::Int
    scenario::Int
    subpaths::Vector{SubPath}
    num_passengers::Int
    actual_schedule::Vector{Tuple{Int, Time}}
    actual_route::Vector{Int}
end
Trip(l::Int, t::Int, s::Int) = Trip(l,t,s,Vector{SubPath}(), 0, Vector{Tuple{Int, Time}}(), Vector{Int}())
"""
Object that captures the level of service of passengers served
"""
mutable struct ServedPassenger
    aggregPax_index::Tuple{Int, Int} # index of the corresponding AggregPassenger object
    originStopOSMid::Int
    pickupTrip::Tuple{Int,Int} # (l, t)
    pickupTime::Time
    pickupStopOSMid::Int
    numberPassengers::Int
    walkTimeSeconds::Int
    waitTimeSeconds::Int
    arrivalDelaySeconds::Int
    inVehicleTimeSeconds::Int
    walkPath::Vector{Int} # set of OSM node ids forming the path
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
mutable struct Sol
    cost::Float64
    solve_time::Float64
    trips::Vector{Trip}
    num_passengers::Int
    served_passengers::Vector{Vector{ServedPassenger}} # one vector per demand scenario
    transitDistanceMeters::Float64
    actualDistanceMeters::Float64
end
Sol() = Sol(0.,0., Vector{Trip}(),0,Vector{Vector{ServedPassenger}}(),0.,0.)
Sol(numS::Int) = Sol(0.,0., Vector{Trip}(),0,[Vector{ServedPassenger}() for s in 1:numS],0.,0.)

"""
Build solution and compute level of service
"""
function buildSol(
        m::MapData,
        model::JuMP.Model,
        R::RouteSettingData,
        Inst::InstanceSettingData,
        Subpaths::Vector{Vector{Vector{Pool}}},
        subpath_road_networks,
        all_load_expanded_graphs,
        all_subpath_graphs,
        normalized::Bool=NORMALIZED,
    )
    @unpack numL, numS, TDsp, MaxWalk, λ, μ, σ, δ = Inst
    sol = Sol(numS)
    sol.cost = objective_value(model)
    sol.solve_time = solve_time(model)
    xVal = value.(model[:x])
    yVal = value.(model[:y])
    for l in 1:numL, t in eachindex(R.Lines[l].freq)
        if xVal[l,t] > EPS
            line_start_time = R.Lines[l].ref_stop_time[R.Lines[l].ref_stop_id[1]][t]
            for s in 1:numS
                trip = Trip(l,t,s)
                actualLineSchedule = Vector{Tuple{Int,Time}}()
                subpathRoutes = [Vector{Int}() for i in 1:length(R.Lines[l].ref_stop_id)-1]
                for i in eachindex(Subpaths[l][t][s].all)
                    if yVal[l,t,s,i] > EPS
                        trip.num_passengers += Subpaths[l][t][s].all[i].load
                        push!(trip.subpaths, Subpaths[l][t][s].all[i])

                        # calculate level of service of served passengers
                        subpathSetPassengers = Subpaths[l][t][s].all[i].pax
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
                                        # if (x1 > 1 && osm_node != subpathRoutes[x1-1][end]) || x1 == 1
                                            push!(subpathRoutes[x1], osm_node)
                                        # end
                                    end
                                end
                            else
                                push!(subpathRoutes[x1], R.Lines[l].ref_stop_id[x1])
                            end
                        end
                        # compute total number of customers picked up and their saved time
                        for (idxPassenger, subPathPassenger) in enumerate(subpathSetPassengers)
                            (passengerOriginStop,timeRequestTDmp, numberOfPassengers) = subPathPassenger
                            pax = R.Demand[s][passengerOriginStop,timeRequestTDmp] #NOTE: remember is an AggregPassenger object
                            directInVehicleTimeSeconds = convert(Dates.Second, pax.request_dropoff_time - pax.taxi_pickup_time).value
                            be_ready_time = pax.be_ready_times[l,t]
                            # compute arrival delay (negative if early)
                            arrivalDelaySeconds = ceil(Int, Dates.value(Dates.Second(R.Lines[l].arrival_time[t] - pax.request_dropoff_time)))
                            # compute walking and waiting time
                            walkTimeSeconds = typemax(Int)
                            waitTimeSeconds = typemax(Int)
                            inVehicleTimeSeconds = typemax(Int)
                            walkingPath = Int[]
                            pickupStopOSMid = -1
                            pickupTime = Time(0)
                            bestCombo = typemax(Int)
                            dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(R.Walk_graph, m.v[passengerOriginStop])
                            for n in subpathTENpath[2:end-1]
                                (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                                pickupTime = Dates.Second((time + t1)*TDsp) + line_start_time
                                # compute walking path from passengerOriginStop to loc
                                tempWalkingPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[subpath_road_networks[l][x1][x2-x1].B2A[loc]])
                                currentWalkingDistance = 0.
                                for (idx, node) in enumerate(tempWalkingPath[1:end-1])
                                    node2 = tempWalkingPath[idx+1]
                                    currentWalkingDistance += R.Walk_graph.weights[node, node2]
                                end
                                # calculate in-vehicle time
                                if currentWalkingDistance <= MaxWalk
                                    tempWalkingTimeSeconds = ceil(Int, (currentWalkingDistance / WALKSPEED_M_PER_SEC))
                                    tempWaitingTimeSeconds = ceil(Int, Dates.value(Dates.Second(pickupTime - be_ready_time))+tempWalkingTimeSeconds)
                                    tempInVehicleTimeSeconds = ceil(Int, Dates.value(Dates.Second(R.Lines[l].arrival_time[t] - pickupTime)))
                                    tempInVehicleTimeNorm = normalized ? tempInVehicleTimeSeconds/(R.Taxi_travel_times[subpath_road_networks[l][x1][x2-x1].B2A[loc]]) : tempInVehicleTimeSeconds
                                    if  tempWaitingTimeSeconds > -0.01 && λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds + σ*tempInVehicleTimeNorm < bestCombo - 0.01
                                        walkTimeSeconds = tempWalkingTimeSeconds
                                        waitTimeSeconds = tempWaitingTimeSeconds
                                        walkingPath = [m.n[ind] for ind in tempWalkingPath]
                                        inVehicleTimeSeconds = tempInVehicleTimeSeconds
                                        pickupStopOSMid = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                        bestCombo = λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds + σ*tempInVehicleTimeNorm
                                    end
                                end
                            end
                            newPassenger = ServedPassenger((passengerOriginStop,timeRequestTDmp), passengerOriginStop, (l, t), pickupTime, pickupStopOSMid, numberOfPassengers, walkTimeSeconds, waitTimeSeconds, arrivalDelaySeconds, inVehicleTimeSeconds, walkingPath, directInVehicleTimeSeconds)
                            push!(sol.served_passengers[s], newPassenger)
                        end
                        # build line schedule (stop, time)
                        if length(subpathTENpath)> 1.1
                            for n in subpathTENpath[2:end-1]
                                (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                                #  we need to add t1 to time
                                x = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                push!(actualLineSchedule, (x, Dates.Second((time+t1)*TDsp)+line_start_time))
                            end
                        end
                    end
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
function addSubproblemSol!(
    m::MapData,
    model::JuMP.Model,
    sol::Sol,
    xVal,
    l::Int,
    t::Int,
    s::Int,
    R::RouteSettingData,
    Inst::InstanceSettingData,
    Subpaths::Pool,
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
        trip = Trip(l,t,s)
        actualLineSchedule = Vector{Tuple{Int,Time}}()
        subpathRoutes = [Vector{Int}() for i in 1:length(R.Lines[l].ref_stop_id)-1]
        for i in eachindex(Subpaths.all)
            if yVal[i] > EPS
                trip.num_passengers += Subpaths.all[i].load
                push!(trip.subpaths, Subpaths.all[i])

                # calculate level of service of served passengers
                subpathSetPassengers = Subpaths.all[i].pax
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
                                # if (x1 > 1 && osm_node != subpathRoutes[x1-1][end]) || x1 == 1
                                    push!(subpathRoutes[x1], osm_node)
                                # end
                            end
                        end
                    else
                        push!(subpathRoutes[x1], R.Lines[l].ref_stop_id[x1])
                    end
                end
                # compute total number of customers picked up and their saved time
                for (idxPassenger, subPathPassenger) in enumerate(subpathSetPassengers)
                    (passengerOriginStop,timeRequestTDmp, numberOfPassengers) = subPathPassenger
                    pax = R.Demand[s][passengerOriginStop,timeRequestTDmp] #NOTE: remember is an AggregPassenger object
                    directInVehicleTimeSeconds = convert(Dates.Second, pax.request_dropoff_time - pax.taxi_pickup_time).value
                    be_ready_time = pax.be_ready_times[l,t]
                    # compute arrival delay (negative if early)
                    arrivalDelaySeconds = ceil(Int, Dates.value(Dates.Second(R.Lines[l].arrival_time[t] - pax.request_dropoff_time)))
                    # compute walking and waiting time
                    walkTimeSeconds = typemax(Int)
                    waitTimeSeconds = typemax(Int)
                    inVehicleTimeSeconds = typemax(Int)
                    walkingPath = Int[]
                    pickupStopOSMid = -1
                    pickupTime = Time(0)
                    bestCombo = typemax(Int)
                    dijkstraState = SimpleWeightedGraphs.dijkstra_shortest_paths(R.Walk_graph, m.v[passengerOriginStop])
                    for n in subpathTENpath[2:end-1]
                        (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                        pickupTime = Dates.Second((time + t1)*TDsp) + line_start_time
                        # compute walking path from passengerOriginStop to loc
                        tempWalkingPath = SimpleWeightedGraphs.enumerate_paths(dijkstraState, m.v[subpath_road_networks[l][x1][x2-x1].B2A[loc]])
                        currentWalkingDistance = 0.
                        for (idx, node) in enumerate(tempWalkingPath[1:end-1])
                            node2 = tempWalkingPath[idx+1]
                            currentWalkingDistance += R.Walk_graph.weights[node, node2]
                        end
                        # calculate in-vehicle time
                        if currentWalkingDistance <= MaxWalk
                            tempWalkingTimeSeconds = ceil(Int, (currentWalkingDistance / WALKSPEED_M_PER_SEC))
                            tempWaitingTimeSeconds = ceil(Int, Dates.value(Dates.Second(pickupTime - be_ready_time))+tempWalkingTimeSeconds)
                            tempInVehicleTimeSeconds = ceil(Int, Dates.value(Dates.Second(R.Lines[l].arrival_time[t] - pickupTime)))
                            tempInVehicleTimeNorm = normalized ? tempInVehicleTimeSeconds/(R.Taxi_travel_times[subpath_road_networks[l][x1][x2-x1].B2A[loc]]) : tempInVehicleTimeSeconds
                            if  tempWaitingTimeSeconds > -0.01 && λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds + σ*tempInVehicleTimeNorm < bestCombo - 0.01
                                walkTimeSeconds = tempWalkingTimeSeconds
                                waitTimeSeconds = tempWaitingTimeSeconds
                                walkingPath = [m.n[ind] for ind in tempWalkingPath]
                                inVehicleTimeSeconds = tempInVehicleTimeSeconds
                                pickupStopOSMid = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                                bestCombo = λ*tempWalkingTimeSeconds + μ*tempWaitingTimeSeconds + σ*tempInVehicleTimeNorm
                            end
                        end
                    end
                    newPassenger = ServedPassenger((passengerOriginStop,timeRequestTDmp), passengerOriginStop, (l, t), pickupTime, pickupStopOSMid, numberOfPassengers, walkTimeSeconds, waitTimeSeconds, arrivalDelaySeconds, inVehicleTimeSeconds, walkingPath, directInVehicleTimeSeconds)
                    push!(sol.served_passengers[s], newPassenger)
                end
                # build line schedule (stop, time)
                if length(subpathTENpath)> 1.1
                    for n in subpathTENpath[2:end-1]
                        (loc, time) = all_subpath_graphs[l,s][t][x1][x2-x1].n2p[n]
                        #  we need to add t1 to time
                        x = subpath_road_networks[l][x1][x2-x1].B2A[loc]
                        push!(actualLineSchedule, (x, Dates.Second((time+t1)*TDsp)+line_start_time))
                    end
                end
            end
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
* `R` - RouteSettingData
* `sol` - Sol
* `directDist` - Dictionary (node ID => direct taxi trip distance in meters)
* `los_fn` - filename to record level of service
* `trial_id` - trial ID
* `lineDict` - dictionary of line index to line ID
### Writes
* one line of level of service metrics to `los_fn` per (l, t, s)
"""
function summarise_los(
        m::MapData,
        inst::InstanceSettingData,
        R::RouteSettingData,
        sol::Sol,
        directDist::Dict{Int,Float64},
        los_fn::String,
        trial_id::Int,
        lineDict::Dict{Int,Int}
    )
    #--- sort served passengers and trips
    selectedTrips = unique([(trip.line, trip.freq, trip.scenario) for trip in sol.trips])
    servedPassengers = Dict((l, t, s) => ServedPassenger[] for (l, t, s) in selectedTrips)
    trips = Dict((trip.line, trip.freq, trip.scenario) => trip for trip in sol.trips)
    for s in 1:inst.numS, pax in sol.served_passengers[s]
        (l, t) = pax.pickupTrip
        push!(servedPassengers[(l, t, s)], pax)
    end

    #--- precompute taxi distances of rejected passengers
    # all individual passenger pickups
    allPax = Dict(s => [(pax, aggPax.origin_stop) for aggPax in values(R.Demand[s]) for pax in aggPax.pax] for s in 1:inst.numS)
    passengerPickups = Dict(s => [] for s in 1:inst.numS)
    for trip in sol.trips
        for subpath in trip.subpaths
            for pax in subpath.pax
                aggPax = R.Demand[trip.scenario][(pax[1], pax[2])]
                for p in aggPax.pax
                    push!(passengerPickups[trip.scenario], p)
                end
            end
        end
    end

    # compute taxi distances of passengers who weren't picked up by scenario
    taxiDist = Dict(s => 0. for s in 1:inst.numS)
    for (s, paxList) in allPax
        for (pax, originStop) in paxList
            # served passenger
            if pax in passengerPickups[s]
                continue
            end
            # direct distance
            if originStop in keys(directDist)
                taxiDist[s] += directDist[originStop]
            else
                # compute new dist
                shortest_dist = typemax(Float64)
                for i in eachindex(EXIT_POINTS_OSM)
                    _, rDist, _ = fastest_route(m, originStop, EXIT_POINTS_OSM[i]; speeds=SPEEDS)
                    if rDist + EXIT_DIST_METERS[i] < shortest_dist - 0.01
                        shortest_dist = rDist + EXIT_DIST_METERS[i]
                    end
                end
                if shortest_dist < 999999
                    taxiDist[s] += shortest_dist
                else
                    taxiDist[s] += mean(values(directDist))
                end
            end
        end
    end

    # write the first line
    f =  open(los_fn, "w")
    write(f, "trial_id,l,t,s,num_pickups,avg_detour,max_detour,avg_delay_sec,avg_delay_norm,max_delay_sec,max_delay_norm,avg_walk_sec,max_walk_sec,avg_wait_sec,max_wait_sec,dist_m,taxi_dist_m\n")

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

        # walk
        avg_walk_sec = sum(pax.numberPassengers * pax.walkTimeSeconds for pax in setPax) / num_pickups
        max_walk_sec = maximum([pax.walkTimeSeconds for pax in setPax])

        # wait
        avg_wait_sec = sum(pax.numberPassengers * pax.waitTimeSeconds for pax in setPax) / num_pickups
        max_wait_sec = maximum([pax.waitTimeSeconds for pax in setPax])

        # distance
        dist_m = EXIT_DIST_METERS[R.Lines[l].exit_point_idx]
        trip = trips[(l, t, s)]
        for i in 1:length(trip.actual_schedule)-1
            n1, _ = trip.actual_schedule[i]
            n2, _ = trip.actual_schedule[i+1]
            dist_m += OpenStreetMapX.distance(m.nodes[n1], m.nodes[n2])
        end

        metrics = [trial_id, lineDict[l], t, s, num_pickups, avg_detour, max_detour, avg_delay_sec, avg_delay_norm, max_delay_sec, max_delay_norm, avg_walk_sec, max_walk_sec, avg_wait_sec, max_wait_sec, dist_m, taxiDist[s]]
        write(f, join(metrics, ",") * "\n")
    end
    close(f)
    return nothing
end

"""
Write the first-stage solution (x and z values) to a file `sol_fn`
"""
function writeFirstStageSol(R::RouteSettingData, Inst::InstanceSettingData, sol_fn::String, xVal)
    f = open(sol_fn, "w")
    for l in 1:Inst.numL
        for t in eachindex(R.Lines[l].freq)
            write(f, string(xVal[l,t],","))
        end
        write(f, "\n")
    end
    close(f)
    return nothing
end
"""
Read the first-stage solution (x and z values) from a file `sol_fn`
"""
function readFirstStageSol(sol_fn::String, R::RouteSettingData, Inst::InstanceSettingData)
    xVal = Dict{Tuple{Int, Int}, Float64}()
    f = readlines(sol_fn)
    counter = 0
    for l in 1:Inst.numL
        counter += 1
        if length(f[counter]) > EPS
            vals = parse.(Float64, split(strip(f[counter], ','), ","))
            for t in eachindex(R.Lines[l].freq)
                xVal[l,t] = vals[t]
            end
        end
    end
    return xVal
end
