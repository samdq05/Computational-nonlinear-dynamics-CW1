
### matching test file ###


using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_4.jl"), include("CW_6.jl"), include("fg3_adjustedbasins.jl")#include("fg3_FORALLOFTHEM.jl")

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

    # perhaps change basin calcs to be findall(x -> x == 1, etc)

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


# --- Matching (replacement for your ColourIndex building code) ---
MatchingThreshold = 0.15  # threshold for matching std(w) to an existing group  
# somewhere between 0.6-0.7 it links the broken attractor together

ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonMeans = copy(AllStdW[1])                     # representative std(w) for each group
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
            #PreviousValue
            ComparisonMeans[MinimumIndex] = CurrentValue # value to compare to on next stwp 
            #ComparisonMeans[MinimumIndex] = (ComparisonMeans[MinimumIndex]*(counts[MinimumIndex]-1) + CurrentValue) / counts[MinimumIndex] # mean value of all values
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

############################

#=
# --- One-to-one stepwise matching ---
match_thresh = 0.4

ColourIndex = Vector{Vector{Int}}()

prev_std = AllStdW[1]
prev_groups = collect(1:length(prev_std))
push!(ColourIndex, prev_groups)
next_group_id = length(prev_std) + 1

for i in 2:length(AllStdW)
    current_std = AllStdW[i]
    current_groups = Int[]

    used_prev = falses(length(prev_std))

    for val in current_std
        dists = abs.(prev_std .- val)
        min_dist, idx = findmin(dists)

        if min_dist < match_thresh && !used_prev[idx]
            push!(current_groups, prev_groups[idx])
            used_prev[idx] = true
        else
            push!(current_groups, next_group_id)
            next_group_id += 1
        end
    end

    push!(ColourIndex, current_groups)
    prev_std = current_std
    prev_groups = current_groups
end=#
#=
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")

using Colors

n_groups = maximum(vcat(ColourIndex...))
palette = distinguishable_colors(n_groups)

# Each attractor group gets a consistent color across all parameter values
for (i_p, pval) in enumerate(prange)
    for (j, group_id) in enumerate(ColourIndex[i_p])
        scatter!(BoAFractionsAxis, pval, AllBoAFractions[i_p][j], color = palette[group_id])
        scatter!(MaxVAxis, pval, AllVMax[i_p][j], color = palette[group_id])
        scatter!(MaxWAxis, pval, AllWMax[i_p][j], color = palette[group_id])
        scatter!(StdAxis, pval, AllStdW[i_p][j], color = palette[group_id])
    end
end

Label(GlobalContinuationFigure[0, :],
      text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true),
      fontsize = 30)

      GlobalContinuationFigure =#

      ##################################################

      # Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; limits = ((0.0, 0.05), nothing), title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        lines!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], linewidth = 2)
        scatter!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        lines!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], linewidth = 2)
        scatter!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        lines!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], linewidth = 2)
        scatter!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        lines!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], linewidth = 2)
    end
end
# For the basins of attraction, I want to separate the points into separate lines

# Labels
# Latex label
#Label(GlobalContinuationFigure[0, :], text = L"\\textbf{Effects of varying the parameter p}", fontsize = 30)
# Rich label
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Normal label
#Label(GlobalContinuationFigure[0, :], text = "Effects of varying the parameter p", fontsize = 30)

# Add the legend to the figure
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))] #length(Colours)]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))] # loop over number of colours used #length(Colours)]
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")

# Show Figure
GlobalContinuationFigure

println(maximum(maximum(ColourIndex)))

# interesting to note that two green dots exist at the end, come back to later





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
# this works fine - this should do 3.4 fine as well

# How can we check for vanishing attractors?
# make note of when an attractor dissapears and record the value at that point
# --- Matching 
MatchingThreshold = 0.15  # threshold for matching std(w) to an existing group  
# somewhere between 0.6-0.7 it links the broken attractor together

ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)
RowLength = []
# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonValues = copy(AllStdW[1])                     # representative std(w) for each group
counts = fill(1, length(ComparisonValues))              # how many members assigned to each comparison vector
FirstRowIndices = collect(1:length(ComparisonValues)) # [1,2,3,...] (index for each group)
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 
# record row length
#i = 1
#push!(RowLength, length(ComparisonValues))

# Process subsequent rows
for row in AllStdW[2:end]
    RowAssignment = Int[]  # will hold group indices for current row
    #println(length(row)) # check the row length each time for debug
    #i += 1
    for CurrentValue in row
    #=    push!(RowLength, length(ComparisonValues))
        if RowLength[i-1] > RowLength[i] # check if row length is smaller, ie. an attractor has vanished
            # thing
        elseif RowLength[i-1] < RowLength[i] # check if row length is longer, ie. attractor has returned
            # other thing
        end=#
        # compute distances to current representatives
        Difference = abs.(ComparisonValues .- CurrentValue) # compare current row values to previous ones
        MinimumDifference, MinimumIndex = findmin(Difference)
        println(CurrentValue)
        println(ComparisonValues)

        if MinimumDifference < MatchingThreshold
            # assign to the closest existing group and update its representative (running mean)
            push!(RowAssignment, MinimumIndex)
            counts[MinimumIndex] += 1
            #PreviousValue
            ComparisonValues[MinimumIndex] = CurrentValue # value to compare to on next step (if close update that group) 
        #ComparisonMeans[MinimumIndex] = (ComparisonMeans[MinimumIndex]*(counts[MinimumIndex]-1) + CurrentValue) / counts[MinimumIndex] # mean value of all values
        else
            # create a new group
            NewIndex = length(ComparisonValues) + 1
            push!(ComparisonValues, CurrentValue) # if not close add a new comparison value
            push!(counts, 1)
            push!(RowAssignment, NewIndex)
        end
    end
    # record row length
    #push!(RowLength, length(row))
    # add to colour index
    push!(ColourIndex, RowAssignment)
end


##################################################################################################


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
# this works fine - this should do 3.4 fine as well

# simple matching
previous_values = AllStdW[1]
ColourIndex = [1:length(previous_values)]

for current_row in AllStdW[2:end]
    row_assignment = Int[]
    for val in current_row
        differences = abs.(previous_values .- val)
        min_diff, idx = findmin(differences)
        if min_diff < MatchingThreshold
            push!(row_assignment, idx)
        else
            push!(row_assignment, length(previous_values) + 1)
        end
    end
    push!(ColourIndex, row_assignment)
    # Replace previous row — we forget earlier history
    previous_values = current_row
end
## 2 ##
# simple matching
PreviousValues = AllStdW[1]
ColourIndex = [1:length(PreviousValues)]

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

