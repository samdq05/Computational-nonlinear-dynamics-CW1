# General featurise and group function 
function fg3(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)


    ## Featurise ##

    # initialise a features vector (contains feature vectors for each intial condition)
    Features = Vector{Vector{Float64}}(undef, length(ics))

    for i in 1:length(ics)
        # trajectory
        ODES = CoupledODEs(ds, ics[i], Parameters)
        Trajectory, t = trajectory(ODES, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
        # find the features
        Features[i] = featurizer(Trajectory)
    end


    ## Group ##

    GroupIndex = Vector{Int16}(undef, length(ics))

    #k = round(Int, 50 / 0.5) # number of steps to discard due to transience 

    VecGroupCentreIndex = Vector{Int16}(undef, 3) # vector of the indices of group centres

    # do the first non-divergent group
    ##
    GroupIndex[1] = 1
    GroupCentres = Vector{Vector{Float64}}()
    FirstGroupCentre = Vector{Float64}()
    for j = 1:length(Features[1])
        push!(FirstGroupCentre, Features[1][j])
    end
    push!(GroupCentres, FirstGroupCentre) # add the first group centre
    GroupCentresIndex = [1]

    # Assume divergence is false unless a feature returns otherwise
    #divergence = false

    # loop starts from 2
    for i in 2:length(ics)

        ######## Focus on THISS AGHGAAGHGHG
    # First Check if the trajectories diverge based on the divergence feature
    if divergence == true
        # now i need to create a group for this condition and break the loop 
        put 
        break
    end
    ################################

        # Assign groups based on features
        ComparisonVector = Vector{Float64}()
        for n in 1:length(Features[i])

            push!(ComparisonVector, Features[i][n])
        end

        #ComparisonVector = [Features[i][1], Features[i][2], Features[i][3]]
        


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
            push!(GroupCentresIndex, i)
        end
    end

    # add all the divergent values to a group at the end 

    # fix this section to work with an arbitrary featurizer and trajectory
    GroupCentresTrajectories = Vector{Vector{Vector{Float64}}}()
    TransientTime = 0.0
    #k2 = round(Int, 50 ÷ 0.2)
    for i in GroupCentresIndex        # trajectory
        VecTemp = Vector{Vector{Float64}}()
        ODEs = CoupledODEs(ds, ics[i], Parameters)
        Trajectory, t = trajectory(ODEs, T; Ttr=TransientTime, Δt=0.2) # Discarding Ttr seconds of transient time
        for j in 1:length(Trajectory[1])
            push!(VecTemp, Trajectory[:, j])
        end
        push!(GroupCentresTrajectories, VecTemp)
    end

    # Calculate the Basin of Attraction fractions
    GroupDict = Dict{Int16,Int16}()
    FractionsBoA = Vector{Float64}()
    # count each index
    for i in 1:length(GroupIndex)
        if haskey(GroupDict, GroupIndex[i])
            GroupDict[GroupIndex[i]] += 1
        else
            GroupDict[GroupIndex[i]] = 1
        end
    end

    # Find the BoA Fractions
    for i in 1:length(GroupDict)
        TempFraction = GroupDict[i] / length(ics)
        push!(FractionsBoA, TempFraction)
    end

    return GroupIndex, GroupCentresTrajectories, FractionsBoA
end

#############################################################################################

#=
## Test ##

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("Featurizer2.7.jl"), include("DiscreteTimeDynamicalSystem.jl")

# Initialise relevant functions, parameters and initial conditions
ds = DTDS
r = 0.2
Parameters = [2.9, 0.66]
# Time Values
StartTime = 0.0
TransientTime = 0.0
TotalTime = 50
SamplingTime = 1.0

ics = Vector{Vector{Float64}}()
xg = yg = collect(range(-2.5, 2.5; length = 10))#1000))
for x in xg
    for y in yg
            ics_temp = [x, y]
            push!(ics, ics_temp)
    end
end

fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)
=#