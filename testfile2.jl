n, m = 6, 6

p = zeros(n*m)
u1 = zeros(n*m)
for i = 1:n*m
    p[i] = i
    u1 = u[:,i]
    u1
end
#################################

using DynamicalSystems, CairoMakie
include("CW_4.jl")

x, y = zeros(n), zeros(m)
u = zeros(Float64, 2, n * m)
k = 0


for i = 1:n
    x[i] = -3 + (6 / n) * i
    for j = 1:m
        y[j] = -3 + (6 / m) * j
        k += 1
        u[:, k] = [x[i], y[j]]
    end

end
u

params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
total_time = 100.0
sampling_time = 0.01


n, m = 6, 6
p = zeros(n*m)
u1 = zeros(n*m)
for i = 1:n*m
    p[i] = i
    u1 = u[:,i]
    println(u1)

    for u0 in u1
    neuron = CoupledODEs(neuron_system!, u0, params)
    N, t = trajectory(neuron, total_time; Δt=sampling_time)
    v, w = N[:, 1], N[:, 2]
    #println(v, w)
time = collect(t)

    lines!(ax, v, w)
        #v_end, w_end = N[end, 1], N[end, 2]
        #scatter!(ax, v_end, w_end)
    end

end


############################
total_time = 100.0
sampling_time = 0.01
t = collect(0:sampling_time:total_time)
println(t)
#k = sin(t)
#=
fig = Figure()
ax = Axis(fig[:,:], xlabel = "x", ylabel = "y")
lines!(ax, t, k) 
fig
=#

######################################################

using CairoMakie, DynamicalSystems

include("CW_4.jl")

n, m = 6, 6 # n and m don't necessarily have to be the same size
x, y = zeros(n), zeros(m)
u = [Vector{Float64}(undef, 2) for _ in 1:n*m]
k = 0


fig = Figure()
ax = Axis(fig[:, :]; title="many trajectories", xlabel="v", ylabel="w")

increment_size_x = 6 / n
increment_size_y = 6 / m
start_x = -2.5
start_y = -2.5

for i = 1:n
    multiplier = i - 1
    x[i] = start_x + increment_size_x * multiplier
    #x[i] = -3 + (6 / n) * (i - 1)
    for j = 1:m
        multiplier = j - 1
        y[j] = start_y + increment_size_y * multiplier
        #y[j] = -2.5 + (6 / m) * (j - 1) 
        k += 1
        u[k] = [x[i], y[j]]
    end
end

params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
total_time = 100.0
sampling_time = 0.01
start_time = 0.0

v_plot = Vector{Vector{Float64}}()
#time = Vector
t_plot = Vector{Float64}()

time = zeros(n*m)
for i = 1:n*m
    neuron = CoupledODEs(neuron_system!, u[i], params)
    N, t = trajectory(neuron, total_time; Δt=sampling_time) # Discarding 2.5s of transient time
    v, w = N[:, 1], N[:, 2]
    v_end, w_end = N[end, 1], N[end, 2]
    lines!(ax, v, w, color=(:gray, 0.25))
    scatter!(ax, v_end, w_end, color="red", markersize=10)#, marker = :rect)
    #type = typeof(v)
    #print(type)
    push!(v_plot, v)
     
    
    #type = typeof(t)
     #println(type)

     #println(t[1])
#t_plot = t[i]
     #push!(t_plot, t)
end
#fig
#println(v_plot[1])
#println(t_plot[1])

#typeof(t_plot)
time = collect(start_time:sampling_time:total_time)
typeof(time) # Vector{Float64}
typeof(v_plot) # Vector{Vector{Float64}}
typeof(v_plot[1]) # Vector{Float64}
#println(time)

length(v_plot[1]) # same length
length(time) # "

println(time)
#println(v_plot[1])

#=
fig_test =  Figure()
ax_test = Axis(fig_test[:,]; xlabel = "v", ylabel = "y")
line!(ax_test, time, v_plot[1])
fig_test

=#

#std(w)
#mean(x)
############################################################################################################################################################################


using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")

# Initialise relevant functions, parameters and initial conditions

values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# threshold value
r = 0.2

# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0.0, 200.0, 400.0, 0.5

# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions, AllVMax, AllWMax, AllStdW = [], [], [], []

#prange = [-1.0, -0.9, -0.8]

for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end


# --- Matching 
MatchingThreshold = 0.15  # threshold for matching std(w) to an existing group  
# somewhere between 0.6-0.7 it links the broken attractor together

ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonValues = copy(AllStdW[1])                     # representative std(w) for each group
counts = fill(1, length(ComparisonValues))              # how many members assigned to each comparison vector
FirstRowIndices = collect(1:length(ComparisonValues)) # [1,2,3,...] (index for each group)
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 

# Process subsequent rows
for row in AllStdW[2:end] 
    row_assignment = Int[]  # will hold group indices for current row
    for CurrentValue in row
        # compute distances to current representatives
        Difference = abs.(ComparisonValues .- CurrentValue) # compare current row values to previous ones
        MinimumDifference, MinimumIndex = findmin(Difference)

        if MinimumDifference < MatchingThreshold
            # assign to the closest existing group and update its representative (running mean)
            push!(row_assignment, MinimumIndex)
            counts[MinimumIndex] += 1
            #PreviousValue
            ComparisonValues[MinimumIndex] = CurrentValue # value to compare to on next stwp 
            #ComparisonMeans[MinimumIndex] = (ComparisonMeans[MinimumIndex]*(counts[MinimumIndex]-1) + CurrentValue) / counts[MinimumIndex] # mean value of all values
        else
            # create a new group
            NewIndex = length(ComparisonValues) + 1
            push!(ComparisonValues, CurrentValue)
            push!(counts, 1)
            push!(row_assignment, NewIndex)
        end
    end
    push!(ColourIndex, row_assignment)
end

# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Labels
# Latex label
#Label(GlobalContinuationFigure[0, :], text = L"\\textbf{Effects of varying the parameter p}", fontsize = 30)
# Rich label
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Normal label
#Label(GlobalContinuationFigure[0, :], text = "Effects of varying the parameter p", fontsize = 30)
# Legend
# Build legend entries
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))] #length(Colours)]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))] # loop over number of colours used #length(Colours)]

# Add the legend to the figure
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")
# Show Figure
GlobalContinuationFigure

########################################################################################################################################################################################################


#using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
#include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")

# Initialise relevant functions, parameters and initial conditions

values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# threshold value
r = 0.2

# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0.0, 200.0, 400.0, 0.5

# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions, AllVMax, AllWMax, AllStdW = [], [], [], []

for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end


## Matching Step ## 
MatchingThreshold = 0.2   # threshold for matching std(w) to an existing group

ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonMeans = copy(AllStdW[1])                     # representative std(w) for each group (a current mean to compare to)
counts = fill(1, length(ComparisonMeans))              # how many members assigned to each rep
FirstRowIndices = collect(1:length(ComparisonMeans)) # [1,2,3,...] (index for each group)
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 

# Process subsequent rows
for row in AllStdW[2:end] 
    row_assignment = Int[]  # will hold group indices for current row
    for CurrentValue in row
        # compute distances to current representatives
        Difference = abs.(ComparisonMeans .- CurrentValue) # compare current row values to previous ones
        MinimumDifference, MinimumIndex = findmin(Difference)

        if MinimumDifference < MatchingThreshold
            # assign to the closest existing group and update its representative (running mean)
            push!(row_assignment, MinimumIndex)
            counts[MinimumIndex] += 1
            # Compare Sequenetially
            #PreviousValue[MinimumIndex] = val # value to compare to on next stwp 
            #ComparisonMeans[MinimumIndex] = CurrentValue # value to compare to on next stwp 
            # Compare to a mean
            ComparisonMeans[MinimumIndex] = (ComparisonMeans[MinimumIndex]*(counts[MinimumIndex]-1) + CurrentValue) / counts[MinimumIndex] # adjust means for new values
        else
            # create a new group
            NewIndex = length(ComparisonMeans) + 1
            push!(ComparisonMeans, CurrentValue)
            push!(counts, 1)
            push!(row_assignment, NewIndex)
        end
    end
    push!(ColourIndex, row_assignment)
end


# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Labels
# Latex label
#Label(GlobalContinuationFigure[0, :], text = L"\\textbf{Effects of varying the parameter p}", fontsize = 30)
# Rich label
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Normal label
#Label(GlobalContinuationFigure[0, :], text = "Effects of varying the parameter p", fontsize = 30)
# Legend
# Build legend entries
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))] #length(Colours)]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))] # loop over number of colours used #length(Colours)]

# Add the legend to the figure
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")
# Show Figure
GlobalContinuationFigure

##########################################################################################################################################################################################################


#using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
#include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")

# Initialise relevant functions, parameters and initial conditions

values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# threshold value
r = 0.2

# Time Values
StartTime = 0.0
TransientTime = 200.0
TotalTime = 400.0
SamplingTime = 0.5



# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions = []
AllVMax = []
AllWMax = []
AllStdW = []
#prange = [-1.0, -0.9, -0.8]

for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end

# Matching
MatchingThreshold = 0.2   # threshold for matching std(w) to an existing group
# Initial Step
PreviousValues = AllStdW[1]
ColourIndex = [collect(1:length(PreviousValues))] # index for each initial value
# Subsequent Steps
for CurrentRow in AllStdW[2:end]
    RowAssignment = Int[]
    for val in CurrentRow
        Differences = abs.(PreviousValues .- val)
        MinimumDifference, Index = findmin(Differences)
        if MinimumDifference < MatchingThreshold
            push!(RowAssignment, Index)
        else
            push!(RowAssignment, length(PreviousValues) + 1)
        end
    end
    push!(ColourIndex, RowAssignment)
    # Replace previous row — we forget earlier history
    PreviousValues = CurrentRow
end

# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Title Label (using rich)
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Legend - build legend entries
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))] 
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))] # loop over number of colours used #length(Colours)]
# Add the legend to the figure
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")

# Show Figure
GlobalContinuationFigure

###########################################################################################################################################################################################################

# This is the newer version
#using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
#include("CW_13.jl")
#include("CW_12.jl")
#include("fg3_FORALLOFTHEM.jl")

# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 100)#1000)
ics = [[x, y] for x in xg for y in yg ]

# Initialise relevant functions, parameters and initial conditions
r = 1.5
Parameters = [2.9, 0.66]
StartTime = 0
TransientTime = 200
TotalTime = 50
SamplingTime = 1

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
push!(ds, ODE)
end

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

# Basins of Attraction Plot
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:,1], limits = (-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
for i in 1:length(ics)
    x, y = ics[i] 
    if GroupIndex[i] == 1 
        scatter!(BasinAxis, x, y, color = "LightBlue", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, x, y, color = "PeachPuff", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, x, y, color = "PaleGreen", marker = :rect, markersize = 10)
    #elseif GroupIndex[i] == 4
     #   scatter!(BasinAxis, x, y, color = "LightPink", marker = :rect, markersize = 10)
    end
end
elem_1 = [PolyElement(color = "LightBlue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "PeachPuff", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "PaleGreen", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,2], [elem_1, elem_2, elem_3], ["0.1286", "0.1286", "0.7428"], "Basin of Attraction Fraction",  patchsize = (35, 35), rowgap = 10, framevisible = false)
BasinsFigure 

################################################################################################################################################################################################

#=
ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonMeans = copy(AllStdW[1])                     # representative std(w) for each group (a current mean to compare to)
counts = fill(1, length(ComparisonMeans))              # how many members assigned to each rep
FirstRowIndices = collect(1:length(ComparisonMeans)) # [1,2,3,...] (index for each group)
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 

# Process subsequent rows
for row in AllStdW[2:end] 
    row_assignment = Int[]  # will hold group indices for current row
    for CurrentValue in row
        # compute distances to current representatives
        Difference = abs.(ComparisonMeans .- CurrentValue) # compare current row values to previous ones
        MinimumDifference, MinimumIndex = findmin(Difference)

        if MinimumDifference < MatchingThreshold
            # assign to the closest existing group and update its representative (running mean)
            push!(row_assignment, MinimumIndex)
            counts[MinimumIndex] += 1
            # Compare Sequenetially
            #PreviousValue[MinimumIndex] = val # value to compare to on next stwp 
            #ComparisonMeans[MinimumIndex] = CurrentValue # value to compare to on next stwp 
            # Compare to a mean
            ComparisonMeans[MinimumIndex] = (ComparisonMeans[MinimumIndex]*(counts[MinimumIndex]-1) + CurrentValue) / counts[MinimumIndex] # adjust means for new values
        else
            # create a new group
            NewIndex = length(ComparisonMeans) + 1
            push!(ComparisonMeans, CurrentValue)
            push!(counts, 1)
            push!(row_assignment, NewIndex)
        end
    end
    push!(ColourIndex, row_assignment)
end=#
 
##################################################################################################################################################################################################

# 3D System Attractors Plot

# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0.0, 50.0, 100.0, 0.5 

# Threshold Value and Parameters
r = 0.2
Parameters = 0.16

# Initial Conditions
yg = xg = zg = (range(-3, 3; length = 20))
ics = [[x, y, z + 0.1] for x in xg for y in yg for z in zg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(DS26, InitialConditions, Parameters)
push!(ds, ODE)
end

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer2, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)

# Plot the system attractors - set it so that the labels are the FractionsBoA values
ThreeDimensionaSystemFigure = Figure(size = (800, 500))
TDSAxis = Axis3(ThreeDimensionaSystemFigure[:,1]; azimuth = 2.2, elevation = 0.35  , title = "System Attractors", xlabel = "x", ylabel = "y", zlabel = "z")
colours = ["blue", "red", "green", "pink"]
labels = ["$(FractionsBoA[1])", "$(FractionsBoA[2])", "$(FractionsBoA[3])"]
for i in 1:maximum(maximum(GroupIndex)) # change to a length at some point
    x, y, z  = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2], GroupCentresTrajectories[i][3]
    lines!(TDSAxis, x, y, z, label = labels[i], color = colours[i], linewidth = 10)
    #scatter!(TDSAxis, x[end], y[end], z[end], label = labels[i], color = colours[i], marker = :star6, markersize = 20)
end
axislegend("Basin of Attraction Fractions", position = :rt)

ThreeDimensionaSystemFigure # show figure

##########################################################################################################################################################################################################

# Featurise and Group 3
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)

# Show the Basins of Attraction Fractions
println(FractionsBoA) 
# next step - make this a dictionary which pairs the fractions with the group index
# check why it takes so long with more initial conditions
BoAFractionsDictionary = Dict("Attractor 1" => FractionsBoA[1], "Attractor 2" => FractionsBoA[2], "Attractor 3" => FractionsBoA[3])
BoAFractionsDictionary


##########################################################################################################################################################################################################


## Global continuation with a new system

# parameter range
μrange = vcat(2.8:0.01:3.6)
# initial conditions
values = 8
xg = range(-2.5, 2.5; length=values)
yg = range(-2.5, 2.5; length=values)
ics = [[x, y] for x in xg for y in yg]
# threshold value
r = 1.5
# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0.0, 200.0, 400.0, 0.5
# Initialise empty vectors to store features for different parameter values
AllBoAFractions, AllMeanX, AllMeanY, AllPeriods = [], [], [], []
# Find features at different parameter values
for μ in μrange
    # Change relevant parameters and intial conditions with θ 
    r = 1.5
    Parameters = [2.9, 0.66]
    StartTime, TransientTime, TotalTime, SamplingTime = 0, 200, 500, 1
    # Dynamical System
    ds = Vector{Any}()
    for InitialConditions in ics
        ODE = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
        push!(ds, ODE)
    end
    # Group and featurise for that parameter
    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)
    # Extract components
    x, y = Vector{Vector{Any}}(), Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(x, GroupCentresTrajectories[i][1]), push!(y, GroupCentresTrajectories[i][2])
    end
    # Calculate some features to compare as parameters change
    MeanX, MeanY, Period = Vector{Float64}(), Vector{Float64}(), Vector{Float64}() # Max v, Max w and Std w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MeanX, mean(x[i])), push!(MeanY, maximum(y[i])), push!(Period, period(GroupCentresTrajectories[i]))
    end
    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA), push!(AllMeanX, MeanX), push!(AllMeanY, MeanY), push!(AllPeriods, Period)
end

## Matching Step ## 
MatchingThreshold = 0.1   # threshold for matching std(w) to an existing group
ColourIndex = simple_matching(MatchingThreshold, AllPeriods)

#########################

# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; limits = (-1.1, 0.05, -0.05, 1.05), title = "BoA fractions against μ", xlabel = "μ", ylabel = "BoAFractions")
MeanXAxis = Axis(GlobalContinuationFigure[1,2]; title = "mean(x)) against μ", xlabel = "μ", ylabel = "Mean x")
MeanYxis = Axis(GlobalContinuationFigure[2,1]; title = "mean(y)) against μ", xlabel = "μ", ylabel = "Mean y")
PeriodAxis = Axis(GlobalContinuationFigure[2,2]; title = "Period against μ", xlabel = "μ", ylabel = "Period of attractor trajectory")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(μrange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, μrange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MeanXAxis, μrange[i], AllMeanX[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MeanYxis, μrange[i], AllMeanX[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(PeriodAxis, μrange[i], AllPeriods[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Title
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter μ"; fontweight = "bold", underline = true), fontsize = 30)
# Legend
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))]
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")
# Show Figure
GlobalContinuationFigure

####################################################################################################################################################################################################################