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
include(joinpath(@__DIR__, "path_benchmark.jl"))

"""
In-sample solutions for the path-based direct model

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
    ioStream = open(joinpath(output_dir, string(trial_id) * "_Mod_path.csv"), "w")
    write(ioStream, join(["trial_id", "model_type", "LB", "UB", "solve_time", "num_second_stage_vars", "precomp_time", "enumeration_time"], ",") * "\n")

    # dummy compilation run
    inst = InstanceSettingData(2, 1800, 2, TRAIN, collect(1:2), 1,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,500.,10,210.,MAX_WAITING_SECONDS,max_sch_dev,10, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    R, m = buildRouteSettingData(inst); 
       
    all_subpaths, _, _, _, enumeration_time = generateSubPathSet(m, R, inst, !CG, !TRANSIT);
    all_paths = [[Vector{PathPool}() for t in eachindex(R.Lines[l].freq)] for l in eachindex(R.Lines)]
    for l in 1:2, t in eachindex(R.Lines[l].freq), s in 1:2
        pool = generatePaths(R, l,t,s, all_subpaths[l][t][s])
        push!(all_paths[l][t], pool)
    end
    model_path = directTwoStageModelPathBased(R, all_paths, inst)
    optimize!(model_path)

    
    
    inst = InstanceSettingData(numL, demT, numS, TRAIN, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,max_sch_dev,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR)
    R, m = buildRouteSettingData(inst);
    total_enumeration_time = 0.

    #--- path-based model
    (all_subpaths, _, _, _, enumeration_time), precomp_time_sp = @timed generateSubPathSet(m, R, inst, !CG, !TRANSIT);
    total_enumeration_time += enumeration_time
    all_paths = [[Vector{PathPool}() for t in eachindex(R.Lines[l].freq)] for l in eachindex(R.Lines)]
    num_paths = 0
    precomp_time_p = 0.
    for l in 1:numL, t in eachindex(R.Lines[l].freq), s in 1:numS
        (pool), path_enum_time = @timed generatePaths(R, l,t,s, all_subpaths[l][t][s])
        push!(all_paths[l][t], pool)
        num_paths += length(pool.all)
        total_enumeration_time += path_enum_time
        precomp_time_p += path_enum_time
    end
    model_path = directTwoStageModelPathBased(R, all_paths, inst)
    optimize!(model_path)
    UB_p = has_values(model_path) ? objective_value(model_path) : typemax(Float64)
    LB_p = objective_bound(model_path)
    solve_time_p = solve_time(model_path)
    write(ioStream, join([trial_id, "path", LB_p, UB_p, solve_time_p, num_paths, precomp_time_sp+precomp_time_p, total_enumeration_time], ",") * "\n")
    
    close(ioStream)

    return nothing
end

####################################################
########################## PIPELINE
####################################################

main()
