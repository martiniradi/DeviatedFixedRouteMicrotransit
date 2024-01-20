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
function addSubPath!(
    R::RouteSettingData,
    l::Int,
    t::Int,
    s::Int,
    m::JuMP.Model,
    sp::SubPath,
    i::Int,
    j::Int,
    )
    yNew = @variable(m, base_name="y", lower_bound=0)
    touchedConstraints = ConstraintRef[]
    vals = Float64[]
    push!(touchedConstraints, m[:c3][i])
    push!(vals, 1.)
    push!(touchedConstraints, m[:c3][j])
    push!(vals, -1.)
    for (k,h,num) in sp.pax
        push!(touchedConstraints, m[:c4][(k,h),(l,t)])
        push!(vals, 1.)
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
function buildGenericSecondStage(
    R::RouteSettingData,
    Subpaths::Pool,
    Inst::InstanceSettingData,
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
    @constraint(m, c4[p in eachindex(Demand[s]), (l_,t_) in Demand[s][p].candidateTrips; (l_,t_) == (l,t)], sum(y[id] for id in Subpaths.pax[p]) <= 1)

    return m
end
"""
Build first stage model.

### Keywords
* `R` - RouteSettingData
* `Inst` - InstanceSettingData
* `time_limit_sec` - maximum solve time in seconds
* `mip_gap` - epsilon tolerance
### Returns
* JuMP Model
"""
function firstStageModel(
    R::RouteSettingData,
    Inst::InstanceSettingData,
    time_limit_sec::Int64=MP_TIME_LIMIT,
    mip_gap::Float64=MIP_GAP,
    num_threads::Int64=NUM_THREADS,
    num_focus::Int64=NUM_FOCUS,
    )

    @unpack numL, numS, TDmp, TDsp, Fleet, Theta = Inst # master and subproblem time discretization, and normalized (true/false) in-vehicle time and delay
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