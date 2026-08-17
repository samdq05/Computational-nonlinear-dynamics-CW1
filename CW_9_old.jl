# General featurise and group function 
function fg3(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

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
    GroupCentresIndex = [1]

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
            push!(GroupCentresIndex, i)
        end
    end

# fix this section to work with an arbitrary featurizer and trajectory
    GroupCentresTrajectories = Vector{Vector{Vector{Float64}}}()
    TransientTime = 0.0
    #k2 = round(Int, 50 ÷ 0.2)
    for i in GroupCentresIndex        # trajectory
        neuron = CoupledODEs(ds, ics[i], params)
        N, t = trajectory(neuron, T; Ttr=TransientTime, Δt=0.2) # Discarding Ttr seconds of transient time
        v_temp, w_temp = N[:, 1], N[:, 2] # separate by component
        #v, w = v_temp[k2:end], w_temp[k2:end] # ignore section before transient time
        u_temp = [v_temp, w_temp]
        push!(GroupCentresTrajectories, u_temp)
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