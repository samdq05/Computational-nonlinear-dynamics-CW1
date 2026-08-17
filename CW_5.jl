#########################################

function simple_grouping(features, r)

    GroupIndex = Vector{Int16}(undef, length(features))
    GroupCentres = Vector{Vector{Float64}}()
    # Create the first group centre from the first feature vector
    FirstGroupCentre = Vector{Float64}()
    for j = 1:length(features[1])
        push!(FirstGroupCentre, features[1][j])
    end
    push!(GroupCentres, FirstGroupCentre) # add the first group centre
    GroupIndex[1] = 1

    for i in 2:length(features)

        # Assign groups based on features
        ComparisonVector = Vector{Float64}()
        for n in 1:length(features[i])
            push!(ComparisonVector, features[i][n])
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