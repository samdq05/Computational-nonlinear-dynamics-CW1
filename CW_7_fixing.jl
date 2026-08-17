# General featurise and group function 
function fg1(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

    ## Featurise ##

    # initialise a features vector (contains feature vectors for each intial condition)
    Features = Vector{Vector{Float64}}(undef, length(ics))
    i = 1
    for DS in ds
        # trajectory
        Trajectory, t = trajectory(DS, T; Ttr, Dt, t0 = 0.0)         
        # find the features
        Features[i] = featurizer(Trajectory)
        i += 1
    end

    ## Group ##

    GroupIndex = Vector{Int16}(undef, length(ics))
    GroupCentres = Vector{Vector{Float64}}()
        # Create the first group centre from the first feature vector
        FirstGroupCentre = Vector{Float64}()
        for j = 1:length(Features[1])
            push!(FirstGroupCentre, Features[1][j])
        end
        push!(GroupCentres, FirstGroupCentre) # add the first group centre
        GroupIndex[1] = 1

    # loop starts from 2
    for i in 2:length(ics)

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
        end
    end

    return GroupIndex
end