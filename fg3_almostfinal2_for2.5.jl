# General featurise and group function 
function fg3(ds, featurizer, r, ics; Ttr, Dt, T)#Ttr=50.0, Dt=0.5, T=100.0)


    ## Featurise ##

    # initialise a features vector (contains feature vectors for each intial condition)
    Features = Vector{Vector{Float64}}(undef, length(ics))
    i = 1
    for InitialConditions in ics#i in 1:length(ics)
        # trajectory
        #ODES = CoupledODEs(ds, ics[i], Parameters) ## double check if it takes ode function or something else
        ODEs = CoupledODEs(ds, InitialConditions, Parameters)
        Trajectory, t = trajectory(ODEs, T; Ttr, Dt, t0 = 0)#; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
        # find the features
        Features[i] = featurizer(Trajectory)
        i += 1
    end


    ## Group ##

     # not sure if there is any point having the first bit separate - since it's not making the first group from the first initial conditions 
     # - I think the error may arise from it not thinking of the divergence as a bool but as a 1 or 0 instead of === for programmatically the same 
     # I've changed it to == to just check value comparison which seems to work for now, in the unlikely event that a values mean is exactly 1.0, it will
        # be considered a divergent trajectory, consider how to fix that   
    GroupIndex = Vector{Int16}(undef, length(ics))

    # do the first non-divergent group
    ##
    # Check first feature vector for Bool-true
    BoolTrueGroup = -1 # will be set to a valid group number later
    GroupCentres = Vector{Vector{Float64}}()
    GroupCentresIndex = Int[]
    if any(feature -> feature == true, Features[1])
        # If the first feature vector contains a true Bool, mark for BoolTrueGroup assignment
        GroupIndex[1] = 0 # placeholder, will set BoolTrueGroup after loop
    else
        # Otherwise, create the first group centre from the first feature vector
        FirstGroupCentre = Vector{Float64}()
        for j = 1:length(Features[1])
            push!(FirstGroupCentre, Features[1][j])
        end
        push!(GroupCentres, FirstGroupCentre) # add the first group centre
        GroupCentresIndex = [1]
        GroupIndex[1] = 1
    end

    # loop starts from 2
    for i in 2:length(ics)

        ######## Focus on THISS AGHGAAGHGHG
        # First Check if the trajectories diverge based on the divergence feature
        # Check for any true Bool feature
        if any(feature -> feature == true, Features[i])
            GroupIndex[i] = 0 # placeholder, will set BoolTrueGroup after loop
            continue
        end

        ################################

        # Assign groups based on features
        ComparisonVector = Vector{Float64}()
        for n in 1:length(Features[i])
            push!(ComparisonVector, Features[i][n])
        end


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

    # Now assign a unique group index for all Bool-true cases
    BoolTrueGroup = length(GroupCentres) + 1
    for i in 1:length(GroupIndex)
        if GroupIndex[i] == 0
            GroupIndex[i] = BoolTrueGroup
        end
    end

    # add all the divergent values to a group at the end 

    # fix this section to work with an arbitrary featurizer and trajectory
    GroupCentresTrajectories = Vector{Vector{Vector{Float64}}}()
    TransientTime = 0.0
    #k2 = round(Int, 50 ÷ 0.2)
    for i in GroupCentresIndex        # trajectory
        VecTemp = Vector{Vector{Float64}}()
        #ODEs = CoupledODEs(ds, ics[i], Parameters)
        ODEs = CoupledODEs(ds, ics[i], Parameters)
        Trajectory, t = trajectory(ODEs, T; Ttr, Dt, t0 = 0)#; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
        for j in 1:length(Trajectory[1, :])
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
