"""
Generates the MapData object of the geographical area of study
### Returns
* Manhattan MapData
"""
function generateNetwork()
    # MapData object
    data_dir = get_data_path()
    m = get_map_data(joinpath(data_dir, "input", "maps", "manhattan.osm"), use_cache=false)
    return m
end

"""
Generates a undirected graph from the directed graph from MapData

### Keywords
* `m` - MapData
### Returns
* Graph - undirected graph
"""
function generateUndirectedGraph(
        m::MapData
    )
    g = SimpleWeightedGraph(length(m.n))
    for (n1,n2) in m.e
        i = m.v[n1]
        j = m.v[n2]
        w = m.w[i,j]
        SimpleWeightedGraphs.add_edge!(g, i, j, w)
    end
    return g
end

"""
Skip coordinates that fall away from actual roadways (sea, other neighborhoods or in the middle of central square)

### Keywords
* `x` - longitude
* `y` - latitude
### Returns
* Bool - whether (`x`,`y`) is in Manhattan
"""
function isstop(
        x,
        y
    )
    return inpolygon((x, y), MANHATTAN_POLY) == 1 && inpolygon((x, y), CS1) == 0 && inpolygon((x, y), CS2) == 0 && inpolygon((x, y), CS3) == 0 && inpolygon((x, y), CS4) == 0 && inpolygon((x, y), CS5) == 0
end

"""
Generate the set of potential stops on the network. Each station is one node from `m`.
### Keywords
* `m` - Manhattan MapData
### Returns
* Vector of stops (Int)
"""
function generateStops(
        m::MapData,
        numX::Int=NUMX,
        numY::Int=NUMY
    )
    # Divide Manhattan in a grid of 100x100m elements that are potential stops
    stops = Vector{Int}()
    @unpack min_y, max_y, min_x, max_x = m.bounds

    step_x = (max_x - min_x) / numX
    step_y = (max_y - min_y) / numY
    for x in min_x:step_x:max_x, y in min_y:step_y:max_y
        if isstop(x, y)
            n = point_to_nodes((y, x), m)
            push!(stops, n)
        end
    end
    # NOTE: we may have multiple stops linked to the same m.n node
    return unique!(stops)
end
"""
Given a location (lat, lon), find the nodes in the map `m` that are within `radius` meters from `loc` (within the `setnodes` subset)
distance is measured as a beeline

### Keywords
* `m` - MapData
* `loc` - ENU location
* `radius` - in meters
* `setnodes` - set of nodes
"""
function getLocRangeNodes(
        m::MapData,
        loc::ENU,
        radius::Float64,
        setnodes::Set{Int64}
    )
    indices = Int[]
    for ind in collect(setnodes)
        dist = OpenStreetMapX.distance(m.nodes[ind], loc)
        if dist < radius
            push!(indices, ind)
        end
    end
    return indices
end
"""
Given an OSM node, find the nearest node in the map `m` that is within the `nodes` subset)

In other words, find closest node to another point

### Keywords
* `m` - MapData
* `n` - OSM node
* `nodes` - set of nodes
"""
function nearest_node_index(m::MapData, n::Int, nodes::Vector{Int64})
    min_dist = typemax(Float64)
    best_index = 0
    loc = m.nodes[n]
    for (idx, i) in enumerate(nodes)
        dist = OpenStreetMapX.distance(m.nodes[i], loc)
        if dist < min_dist - 0.01
            min_dist = dist
            best_index = idx
        end
    end
    return best_index
end
function nearest_node(m::MapData, lat::Float64, lon::Float64, nodes::Vector{Int64})
    min_dist = Inf
    best_node = 0
    pointLLA = LLA(lat,lon)
    loc = OpenStreetMapX.ENU(pointLLA, m.bounds)
    for n in nodes
        dist = OpenStreetMapX.distance(m.nodes[n], loc)
        if dist < min_dist
            min_dist = dist
            best_node = n
        end
    end
    return best_node
end
function nearest_node(m::MapData, loc, nodes::Vector{Int64})
    min_dist = Inf
    best_node = 0
    for n in nodes
        dist = OpenStreetMapX.distance(m.nodes[n], loc)
        if dist < min_dist
            min_dist = dist
            best_node = n
        end
    end
    return best_node
end
function nearest_node_and_dist(m::MapData, loc, nodes::Vector{Int64})
    min_dist = Inf
    best_node = 0
    for n in nodes
        dist = OpenStreetMapX.distance(m.nodes[n], loc)
        if dist < min_dist
            min_dist = dist
            best_node = n
        end
    end
    return best_node, min_dist
end
"""
Compute the travel time by taxi for each of the origin stops and write them to file `TravelTimes.txt`

### Keywords
* `O` - set of origin stops
* `m` - MapData
* `exitPoints` - OSM nodes corresponding to the four exit points
* `exitTimesSeconds` - travel time in seconds from each of the exit points to Laguardia
"""
function computeTravelTimes(O::Vector{Int}, m::MapData, exitPoints::Vector{Int}, exitTimesSeconds::Vector{Float64})
    A = Dict{Int, Float64}() # OSM node -> time in seconds
    @showprogress "Computing arrival times " for o in O
        shortest_time = typemax(Float64)
        for i in eachindex(exitPoints)
            route, rDist, rTime = fastest_route(m, o, exitPoints[i]; speeds=SPEEDS)
            if rTime + exitTimesSeconds[i] < shortest_time - 0.01
                shortest_time = rTime + exitTimesSeconds[i]
            end
        end
        if shortest_time > 999999
            println("THIS TRIP CANNOT BE MADE")
        end
        A[o] = shortest_time
    end
    data_dir = get_data_path()
    f = open(joinpath(data_dir, "input", "MiND-VRP_setup", "TravelTimes.txt"), "w")
    for (k,v) in collect(A)
        write(f, "$k\t$v\n")
    end
    close(f)
    # return A
end
"""
read `TravelTimes.txt` file and return dictionary: (stop id) -> travel time in seconds
"""
function readTravelTimes()
    A = Dict{Int, Float64}() # OSM node -> time in seconds
    data_dir = get_data_path()
    f = open(joinpath(data_dir, "input", "MiND-VRP_setup","TravelTimes.txt"),"r")
    for l in readlines(f)
        vals = parse.(Float64, split(l,"\t"))
        A[Int(vals[1])] = vals[2]
    end
    return A
end
"""
Compute the travel time by taxi between all candidate stops and write them to file `TravelTimesDAR.txt`

### Keywords
* `O` - set of origin stops
* `m` - MapData
"""
function computeTravelTimesDAR(O::Vector{Int}, m::MapData)
    A = Dict{Tuple{Int,Int}, Float64}() # OSM node -> time in seconds
    @showprogress "Computing arrival times " for o1 in O, o2 in O
        route, rDist, rTime = fastest_route(m, o1, o2; speeds=SPEEDS)
        if rTime < LARGE_NUMBER
            A[o1,o2] = rTime
        else
            println("THIS TRIP CANNOT BE MADE")
        end
    end
    data_dir = get_data_path()
    f = open(joinpath(data_dir, "input", "MiND-DAR_setup", "TravelTimesDAR.txt"), "w")
    for (k,v) in collect(A)
        (o1,o2) = k
        write(f, "$o1\t$o2\t$v\n")
    end
    close(f)
    # return A
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
