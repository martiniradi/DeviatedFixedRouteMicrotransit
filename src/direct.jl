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
function directTwoStageModel(
    R::RouteSettingData,
    Subpaths::Vector{Vector{Vector{Pool}}},
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
        c4[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], sum(y[l,t,s,id] for id in Subpaths[l][t][s].pax[p]) <= z[s,p,(l,t)]
    end

    return m
end
"""
Build entire two-stage stochastic model and fix first-stage solution.

### Keywords
* `R` - RouteSettingData
* `Subpaths` - Pool of sub-paths
* `Inst` - instance data
* `xVal` - first-stage x values
* `time_limit_sec` - maximum solve time in seconds
* `mip_gap` - epsilon tolerance
### Returns
* JuMP Model
"""
function directTwoStageModelWithFixFirstStage(
    R::RouteSettingData,
    Subpaths::Vector{Vector{Vector{Pool}}},
    Inst::InstanceSettingData,
    xVal,
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
        sum( R.Lines[l].cost * x[l, t] for l in 1:numL, t in eachindex(R.Lines[l].freq))
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
        c4[s in 1:numS, p in eachindex(Demand[s]), (l,t) in Demand[s][p].candidateTrips], sum(y[l,t,s,id] for id in Subpaths[l][t][s].pax[p]) <= z[s,p,(l,t)]
    end

    # fix first-stage x decisions                    
    for l in 1:numL, t in eachindex(R.Lines[l].freq)
        fix(x[l,t], xVal[l,t], force=true)
    end

    return m
end

