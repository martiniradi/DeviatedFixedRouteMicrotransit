module Microtransit

    # required packages
    using ArgParse
    using Colors
    using Combinatorics
    using CSV
    using DataFrames
    using DataStructures
    using Dates
    using DelimitedFiles
    using Gurobi
    using JuMP # v0.21.4
    using LightGraphs # v1.3.5
    using OpenStreetMapX # v0.2.3
    using OpenStreetMapXPlot
    using Parameters
    using Plots
    using PolygonOps
    using ProgressMeter
    using Random
    using Shapefile
    using SimpleWeightedGraphs
    using Statistics
    using StatsPlots
    using Plots.PlotMeasures

    include("utils.jl")
    include("parameters.jl")        # this file has all the constants
    include("networks.jl")          # this file has all the road-network related functions/structs

    include("input.jl")             # this file has everything input (pre-processing) related
    include("graphs.jl")            # this file has everything related to the different graphs/networks used in the model (time/load-expanded)
    include("subpath.jl")           # this file has everything related to subpaths and their generation
    include("direct.jl")            # this file has the direct implementation of the problem
    include("decomposition.jl")     # this file has functions related to decomposition methods (Benders and Column generation)
    include("solution.jl")          # this file has everything related to post-processing
    include("algorithm.jl")         # this file contains the algorithm (+variants) to run
    include("path_benchmark.jl")    # this file contains everything related to the path-based model
    include("segment_benchmark.jl") # this file contains everything related to the segment-based model
    include("dar.jl")               # this file contains all the functions adapted to the MiND-DAR model
end
