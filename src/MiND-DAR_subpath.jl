using ArgParse
using Colors
using Combinatorics
using CSV
using DataFrames
using DataStructures
using Dates
using Distances
using Gurobi
using JuMP
using LightGraphs
using OpenStreetMapX
using OpenStreetMapXPlot
using Plots
using PolygonOps
using ProgressMeter
using SimpleWeightedGraphs
using UnPack

include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "networks.jl"))
include(joinpath(@__DIR__, "input.jl"))
include(joinpath(@__DIR__, "graphs.jl"))
include(joinpath(@__DIR__, "dar.jl"));

"""
In-sample solutions for transit and microtransit
Out-of-sample solutions for transit and microtransit

File outline:

* Instantiate with parameters:
    - number of lines
    - demand horizon
    - number of scenarios
    - capacity
    - max vehicle deviation
    - max passenger walking
* obtain:
    - microtransit solution
    - oos level of service for microtransit solution
    - transit solution
    - oos level of service for transit solution
* record:
    - first-stage decision variable values for microtransit
    - first-stage decision variable values for transit
"""

CG = true
HEUR = true
TRANSIT = true
TRAIN = true
RELAX_MP = true
TIME_LIMIT_TESTS = 10 * 3600 # 10 hours

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
    demand_hor = trial_info[:demand_hor]
    num_scen = trial_info[:num_scen]
    maxDev = convert(Float64, "max_dev" in names(trial_info) ? trial_info[:max_dev] : MAX_SCHEDULE_DEV)
    maxWalk = convert(Float64, "max_walk" in names(trial_info) ? trial_info[:max_walk] : MAX_WALKING_METERS)
    capacity = convert(Int, "capacity" in names(trial_info) ? trial_info[:capacity] : CAPACITY)
    num_oos = "num_oos" in names(trial_info) ? trial_info[:num_oos] : 5
    fleet_size = Int("fleet_size" in names(trial_info) ?
    trial_info[:fleet_size] : FLEET_SIZE)
    scen_OOS = [1:num_oos;]
    K = Int("K" in names(trial_info) ? trial_info[:K] : 0) + 1
    num_focus = Int("num_focus" in names(trial_info) ? trial_info[:num_focus] : NUM_FOCUS)
    time_limit = Int("time_limit" in names(trial_info) ? trial_info[:time_limit] : TIME_LIMIT_TESTS)

    return output_dir, num_lines, demand_hor, num_scen, capacity, maxDev, maxWalk, trial_id, scen_OOS, fleet_size, K, num_focus, time_limit
end

function main()
    # instantiate
    output_dir, numL, demT, numS, capacity, maxDev, maxWalk, trial_id, scen_OOS, fleet_size, K, num_focus, time_limit = parse_trial_info()
    ioStream = open(joinpath(output_dir, string(trial_id) * "_xDAR.csv"), "w")
    write(ioStream, join(["trial_id", "sol_type", "l", "t"], ",") * "\n")

    inst = InstanceSettingData(numL, demT, numS, TRAIN, collect(1:numS), K,OPERATIONAL_HORIZON,REFNUM,TIME_DISC_MP_DAR,TIME_DISC_SP,maxDev,capacity,maxWalk,MAX_WAITING_SECONDS,MAX_SCHEDULE_DEV,fleet_size, THETA, WEIGHT_COVERAGE,WEIGHT_WALK,WEIGHT_WAIT,WEIGHT_INVEHICLE,WEIGHT_DELAY,SPEED_FACTOR,REF_SPEED_FACTOR_DAR)
    R, m = buildRouteSettingDataDAR(inst);

    #--- microtransit solution
    all_subpathsMT, subpath_road_networksMT, all_load_expanded_graphsMT, all_subpath_graphsMT = generateSubPathSetDAR(m, R, inst, CG, !TRANSIT);
    LB, UB, alg_time, it, num_cuts, num_second_stage_vars, RNtime, solveTMP, solveTSP, xValMT, sol, all_subpathsMT = runAlgDAR(m, R, inst, all_subpathsMT, all_subpath_graphsMT, all_load_expanded_graphsMT, subpath_road_networksMT,CG, HEUR, RELAX_MP, NUM_FOCUS, time_limit);
    #--- write network design solution
    for l in 1:numL, t in eachindex(R.Lines[l].freq)
        if xValMT[l, t] > 0.5
            write(ioStream, join([trial_id, "microtransitDAR", l, t], ",") * "\n")
        end
    end
    close(ioStream)

    return nothing
end

####################################################
########################## PIPELINE
####################################################

main()
