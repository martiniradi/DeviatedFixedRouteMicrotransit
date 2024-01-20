"""
Retrieve data path
"""
function get_data_path()
    filepath = joinpath(dirname(@__DIR__), "data-path.txt")
    if isfile(filepath)
        open(filepath) do f
            return strip(readline(f))
        end
    else
        @warn "Please write path to data directory as the only line in a txt file called `data-path.txt` and store in the root directory."
    end
end

"""
Retrieve speeds dict
"""
function get_speeds()
    speeds_fn = joinpath(get_data_path(), "input", "speeds", "speeds_dict.csv")
    speed_dict = deepcopy(OpenStreetMapX.SPEED_ROADS_URBAN)
    speeds = CSV.read(speeds_fn, DataFrame)

    # for way types represented in Manhattan OSM data
    scales = []
    for row in eachrow(speeds)
        speed_dict[row[:way_type]] = row[:mean_speed_kph]
        push!(scales, row[:mean_speed_kph] / OpenStreetMapX.SPEED_ROADS_URBAN[row[:way_type]])
    end
    scale = sum(scales) / length(scales)

    # for other types, scale down proportionally
    for (k, v) in speed_dict
        if !(k in speeds[:, :way_type])
            speed_dict[k] = scale * v
        end
    end
    return speed_dict
end
