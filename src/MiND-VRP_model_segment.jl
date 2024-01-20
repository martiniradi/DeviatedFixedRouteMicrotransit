using ArgParse
using Combinatorics
using CSV
using DataFrames
using DataStructures
using Dates
using Gurobi
using JuMP
using LightGraphs
using OpenStreetMapX
using PolygonOps
using ProgressMeter
using SimpleWeightedGraphs
using UnPack

include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "networks.jl"))
include(joinpath(@__DIR__, "input.jl"))
include(joinpath(@__DIR__, "graphs.jl"))
include(joinpath(@__DIR__, "subpath.jl"))
include(joinpath(@__DIR__, "solution.jl"))
include(joinpath(@__DIR__, "segment_benchmark.jl"))

"""
In-sample solutions for the segment-based direct model

File outline:

* Instantiate with parameters:
    - number of lines
    - demand horizon
    - number of scenarios
    - capacity
    - max vehicle deviation
    - max passenger walking
* obtain:
    - objective value
    - algorithm run information
* record:
    - 
"""

CG = true
HEUR = true
TRANSIT = true
TRAIN = true


"""
Get trial info, initialize IO
"""
function parse_trial_info()
    # args
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--trial_id", "-t"
            help = "trial ID"
            arg_type = Int
            required = true
        "--output_dir", "-o"
            help = "output directory"
            arg_type = String
            required = true
        "--trial_fn", "-f"
            help = "trials filename"
            arg_type = String
            required = true
    end
    args = parse_args(s)

    # trial attributes
    trial_id = args["trial_id"]
    output_dir = joinpath(get_data_path(), "output", args["output_dir"])
    trials = CSV.read(joinpath(output_dir, args["trial_fn"]), DataFrame)
    trial_info = first(eachrow(filter(row -> row[:trial_id] == trial_id, trials)))
    
    num_lines = Int(trial_info[:num_lines])
    demand_hor = Int(trial_info[:demand_hor])
    num_scen = Int(trial_info[:num_scen])
    maxDev = convert(Float64, "max_dev" in names(trial_info) ? trial_info[:max_dev] : MAX_SCHEDULE_DEV)
    maxWalk = convert(Float64, "max_walk" in names(trial_info) ? trial_info[:max_walk] : MAX_WALKING_METERS)
    capacity = convert(Int, "capacity" in names(trial_info) ? trial_info[:capacity] : CAPACITY)
    fleet_size = Int("fleet_size" in names(trial_info) ? trial_info[:fleet_size] : FLEET_SIZE)
    K = Int("K" in names(trial_info) ? trial_info[:K] : 0) + 1
    max_sch_dev = Int("max_sch_dev" in names(trial_info) ? trial_info[:max_sch_dev] : MAX_SCHEDULE_DEV)

    return output_dir, num_lines, demand_hor, num_scen, capacity, maxDev, maxWalk, trial_id, fleet_size, K, max_sch_dev
end

function main()
    # instantiate
    output_dir, numL, demT, numS, capacity, maxDev, maxWalk, trial_id, fleet_size, K, max_sch_dev = parse_trial_info()
    ioStream = open(joinpath(output_dir, string(trial_id) * "_Mod_segment.csv"), "w")
    write(ioStream, join(["trial_id", "model_type", "LB", "UB", "solve_time", "num_second_stage_vars", "precomp_time", "build_model_time"], ",") * "\n")

    # dummy compilation run
    inst = InstanceSettingData(2, 1800, 2, TRAIN, collect(1:2), 1,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,200.,10,210.,MAX_WAITING_SECONDS,max_sch_dev,10, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    R, m = buildRouteSettingData(inst);    
    all_arcs = generateArcSet(m, R, inst);
    model_segment = directArcBasedModel(R, all_arcs, inst);
    optimize!(model_segment)

    
    
    inst = InstanceSettingData(numL, demT, numS, TRAIN, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,max_sch_dev,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    R, m = buildRouteSettingData(inst);

    #--- segment-based model
    (all_arcs), precomp_time_s = @timed generateArcSet(m, R, inst);
    (model_segment), build_model_time = @timed directArcBasedModel(R, all_arcs, inst);
    optimize!(model_segment)
    UB_s = has_values(model_segment) ? objective_value(model_segment) : typemax(Float64)
    LB_s = objective_bound(model_segment)
    solve_time_s = solve_time(model_segment)
    num_segments = 0
    for (l, line) in enumerate(R.Lines), t in eachindex(line.freq), s in 1:numS
        num_segments += length(all_arcs[l][t][s].all)
    end
    for (l, line) in enumerate(R.Lines), t in eachindex(line.freq), s in 1:numS
        num_segments += (length(line.ref_stop_id)-1)*K
    end
    write(ioStream, join([trial_id, "segment", LB_s, UB_s, solve_time_s, num_segments, precomp_time_s, build_model_time], ",") * "\n")

    close(ioStream)

    return nothing
end

####################################################
########################## PIPELINE
####################################################

main()
