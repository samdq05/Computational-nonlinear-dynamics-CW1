### This one works - YAY! ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl")

# General featurise and group function 
function fg1(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

    ## Featurise ##

    # initialise a features vector (contains feature vectors for each intial condition)
    features = Vector{Vector{Float64}}(undef, length(ics))

    for i in 1:length(ics)
        # trajectory
        neuron = CoupledODEs(ds, ics[i], params)
        N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
        # find the features
        features[i] = featurizer(N)
    end

    ## Group ##

    GroupIndex = Vector{Int16}(undef, length(ics))

    k = round(Int, 50 / 0.5) # number of steps to discard due to transience 

    VecGroupCentreIndex = Vector{Int16}(undef, 3) # vector of the indices of group centres


    # do the first group
    GroupIndex[1] = 1
    GroupCentres = [[features[1][1], features[1][2]]] # index 1 is group 1 , index 2 is group 2 etc  etc etc

    # loop starts from 2
    for i in 2:length(ics)

        # Assign groups based on features
        std_v = features[i][1]
        std_w = features[i][2]
        ComparisonVector = [std_v, std_w]


        putInGroup = false
        # go through all groups. If not within any group, maake a new group, put it in there, and continue
        for GroupPointer in 1:length(GroupCentres)

            # check if comparison vector is within group center
            FeatureDistance = norm(GroupCentres[GroupPointer] - ComparisonVector)

            if FeatureDistance < r
                putInGroup = true
                GroupIndex[i] = GroupPointer
                break
            end
        end

        if putInGroup == false
            push!(GroupCentres, ComparisonVector) # Create new centre if the features aren't in an existing group
            GroupIndex[i] = length(GroupCentres)
        end
    end

    return GroupIndex
end

#Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

l1 = length(v)
l2 = length(w)
K = 0

ics = Vector{Vector{Float64}}(undef, l1 * l2)

for i in 1:l1
    for j in 1:l2
        K += 1
        ics[K] = [v[i], w[j]]
    end
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

# want a plot of the ics colour-coded based on their index
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:,:]; title = "Basins of Attraction", xlabel = "v", ylabel = "w")
for i in 1:length(ics)
    v, w = ics[i] 
    if GroupIndex[i] == 1 
        scatter!(BasinAxis, v, w, color = "blue", marker = :circle, markersize = 8)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, v, w, color = "red", marker = :circle, markersize = 8)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, v, w, color = "green", marker = :circle, markersize = 8)
    end
end
BasinsFigure


# Tells me how many times each index occurs
GroupDict = Dict{Int16, Int16}()

# count each index
for i in 1:length(GroupIndex)
    if haskey(GroupDict, GroupIndex[i])
        GroupDict[GroupIndex[i]] += 1
    else
        GroupDict[GroupIndex[i]] = 1
    end
end
println(GroupDict) # print how many 1s are in group index

# first value automatically becomes first group
# if its within THRESHOLD DIFFERENCE (calculated) of first group put it in the first group
# if it is not (and there are not any other groups) make a second group