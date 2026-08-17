## Functions file

######################################################################################################################

# lorenz63 function for section 1.4
function lorenz63_rule!(du, u, p, t)
    x, y, z = u
    σ, β, ρ = p
    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
    return nothing # always `return nothing` for in-place form!
end

######################################################################################################################

# Poincaré map function for 1.6
function psos1(trajectory, index, value)

    n = size(trajectory, 1) # number of points in the trajectory
    x_store, y_store, z_store = zeros(n), zeros(n), zeros(n) # empty vectors potentially as large as every point
    j = 0 # initialise counter

    # trajecotry components
    x, y, z = columns(trajectory)
    
    # For a given index (x,y,z = 1,2,3 ) and value -- gives a point on the plane and a normal vector to the plane 
    if index == 1
        P = [value, 0.0, 0.0] # point
        v = [1.0, 0.0, 0.0] # normal
        coord = x

    elseif index == 2
        P = [0.0, value, 0.0] # point
        v = [0.0, 1.0, 0.0] # normal   
        coord = y

    elseif index == 3
        P = [0.0, 0.0, value] # point 
        v = [0.0, 0.0, 1.0] # normal   
        coord = z

    else
        error("index must be 1 (x), 2 (y), or 3 (z)") 
    end

    for i = 2:n

        if coord[i] >  value && coord[i-1] < value        
            j += 1 # add one to count
            B = trajectory[i,:] # point after
            A = trajectory[i-1,:] # point before
            r = dot(v, (P - A)) / dot(v, (B - A)) # interpolation fraction
            C = A + r * (B - A) # interpolated point
            x_store[j], y_store[j], z_store[j] = C 
        end

    end
    
    x_final, y_final, z_final = x_store[1:j], y_store[1:j], z_store[1:j]
    return x_final, y_final, z_final
end

######################################################################################################################

# Poincaré map function for 1.7
function psos2(trajectory, tvec, index, value)

    n = size(trajectory, 1) # number of points in the trajectory
    x_final = []
    y_final = []
    z_final = []
    t_final = []

   x, y, z = columns(trajectory)
    
    # For a given index (x,y,z = 1,2,3 ) and value - gives a point on the plane and a normal vector to the plane 
    if index == 1
        P = [value, 0.0, 0.0] # point
        v = [1.0, 0.0, 0.0] # normal
        coord = x

    elseif index == 2
        P = [0.0, value, 0.0] # point
        v = [0.0, 1.0, 0.0] # normal   
        coord = y

    elseif index == 3
        P = [0.0, 0.0, value] # point 
        v = [0.0, 0.0, 1.0] # normal   
        coord = z

    else
        error("index must be 1 (x), 2 (y), or 3 (z)") 
    end

    for i = 2:n

        if coord[i] >  value && coord[i-1] < value
            B = trajectory[i,:] # point after
            A = trajectory[i-1,:] # point before
            r = dot(v, (P - A)) / dot(v, (B - A)) # interpolation fraction
            C = A + r * (B - A) # interpolated point
            push!(x_final, C[1])
            push!(y_final, C[2])
            push!(z_final, C[3])
            push!(t_final, tvec[i])

        end

    end

    return x_final, y_final, z_final, t_final
end

######################################################################################################################

# Neuron 
function neuron_system!(du, initial_state, Parameters, time)
    v, w = initial_state
    a1, a2, a3, a4, c, ϵ, d, I, p = Parameters
    du[1] = p * v + a1 * v^3 - a2 * v^5 + a3 * v^7 - a4 * v^9 - w + I
    du[2] = ϵ * (1 + c * v + d * v^2 - w)
    return nothing
end

#######################################################################################################################

# Simple Grouping
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

#######################################################################################################################

# Featurise a given trajectory 
function featurizer(trajectory)

    # Extract components
    v, w = columns(trajectory)
    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_v = std(v)
    std_w = std(w)

    return [std_v, std_w] # feature vector

end

#######################################################################################################################

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

#############################################################################################################################

# General featurise and group function 
function fg2(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

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
        GroupCentresIndex = Int[]
        # Create the first group centre from the first feature vector
        FirstGroupCentre = Vector{Float64}()
        for j = 1:length(Features[1])
            push!(FirstGroupCentre, Features[1][j])
        end
        push!(GroupCentres, FirstGroupCentre) # add the first group centre
                GroupCentresIndex = [1]
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
            push!(GroupCentresIndex, i)
        end
    end

        # fix this section to work with an arbitrary featurizer and trajectory
    GroupCentresTrajectories = Vector{Vector{Vector{Float64}}}()
    TransientTime = 0.0
    for i in GroupCentresIndex        # trajectory
        VecTemp = Vector{Vector{Float64}}()
        DS = ds[i]
        Trajectory, t = trajectory(DS, T; Ttr, Dt, t0 = 0)
        for j in 1:length(Trajectory[1, :])
            push!(VecTemp, Trajectory[:, j])
        end
        push!(GroupCentresTrajectories, VecTemp)
    end

    return GroupIndex, GroupCentresTrajectories
end

#############################################################################################################################

# General featurise and group function 
function fg3(ds, featurizer, r, ics; Ttr, Dt, T)#Ttr=50.0, Dt=0.5, T=100.0)

    ## Featurise ##
    
    # initialise a features vector (contains feature vectors for each intial condition)
    Features = Vector{Vector{Float64}}(undef, length(ics))
    i = 1
    for DS in ds
        # trajectory
        Trajectory, t = trajectory(DS, T; Ttr, Dt, t0=0.0)
        # find the features
        Features[i] = featurizer(Trajectory)
        i += 1
    end


    ## Group ##
 
    GroupIndex = Vector{Int64}(undef, length(ics))

    # do the first non-divergent group
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

        # First Check if the trajectories diverge based on the divergence feature
        # Check for any true Bool feature
        if any(feature -> feature == true, Features[i])
            GroupIndex[i] = 0 # placeholder, will set the BoolTrueGroup after loop
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
            push!(GroupCentresIndex, i) # keep track of the group centres
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
    for i in GroupCentresIndex        # trajectory
        VecTemp = Vector{Vector{Float64}}()
        DS = ds[i]
        Trajectory, t = trajectory(DS, T; Ttr, Dt, t0=0)
        for j in 1:length(Trajectory[1, :])
            push!(VecTemp, Trajectory[:, j])
        end
        push!(GroupCentresTrajectories, VecTemp)
    end

    # Calculate Basins of Attraction Fractions
    FractionsBoA = zeros(maximum(maximum(GroupIndex)))
    # Bains of Attraction fractions
    for i in 1:maximum(maximum(GroupIndex)) #length(GroupCentres) # loop over number of attractors# find the number of elements with index i
        FractionsBoA[i] = length(findall(x -> x == i, GroupIndex)) / length(GroupIndex)  # Number of elements with index i over total number of elements
    end

    return GroupIndex, GroupCentresTrajectories, FractionsBoA
end

###############################################################################################################################################################

# Featurise a given trajectory 
function featurizer2(trajectory)

    # Extract Components
    x, y, z = columns(trajectory)
    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_x, std_y, std_z = std(x), std(y), std(z)

    return [std_x, std_y, std_z] # feature vector

end

###########################################################################################################################################################

# lorenz63 function for section 1.4
function DS26(du, u, p, t)
    x, y, z = u
    b = p
    du[1] = sin(y) - b * x 
    du[2] = sin(z) - b * y
    du[3] = sin(x) - b * z
    return nothing # always `return nothing` for in-place form!
end

##########################################################################################################################################################

# Featurizer for Task 2.7
function featurizer3(trajectory)

    divergence = false # assumme a trajectory doesn't diverge 

    #x, y = trajectory[:, 1], trajectory[:, 2]
    x, y = columns(trajectory)
    # Features - mean and standard deviations seem appropriate distinguishing characteristic
    mean_x = mean(x)
    mean_y = mean(y)
    #std_x = std(x)
    #std_y = std(y)
    # Special feature -- for diverging trajectories
    if abs(mean_x) > 10_000 || abs(mean_y) > 10_000 # check if trajectories diverge in any direction
        divergence = true
    end

    return [mean_x, mean_y, divergence] # std_x, std_y, divergence] # feature vector
end

#########################################################################################################################################################

# Task 2.7 -- Discrete time dynamical system
function DTDS(u, p, n)
    x,y = u
    μ, j = p
    xn = y
    yn = μ * y - y^3 - j * x
    return SVector(xn, yn)
end

#########################################################################################################################################################

## Simple Matching ## 
function simple_matching(MatchingThreshold, Features)
    #MatchingThreshold = 0.15  # threshold for matching std(w) to an existing group  

    # Initial row
    ComparisonValues = copy(Features[1])  # representative std(w) for each group
    counts = fill(1, length(ComparisonValues))  # start each with 1
    ColourIndex = [collect(1:length(ComparisonValues))]  # first row: unique colours

    # Loop over the rest of the rows
    for i in 2:length(Features)
        row = Features[i]
        row_assignment = zeros(Int, length(row))  # group colours for this row

        for k in 1:length(row)
            CurrentValue = row[k]

            # Find closest comparison value
            Difference = abs.(ComparisonValues .- CurrentValue)
            MinimumDifference, MinimumIndex = findmin(Difference)

            # Check threshold
            if MinimumDifference < MatchingThreshold
                # Match the colour of the corresponding comparison value - use the previous row's colour index at that comparison position
                if MinimumIndex <= length(ColourIndex[i-1])
                    MatchedIndex = ColourIndex[i-1][MinimumIndex]
                end

                row_assignment[k] = MatchedIndex
                counts[MatchedIndex] += 1

            else
                # Create a new group
                NewIndex = maximum(vcat(ColourIndex...)) + 1
                row_assignment[k] = NewIndex
                push!(counts, 1)
            end
        end

        # Update for next iteration
        push!(ColourIndex, row_assignment)
        ComparisonValues = row
    end
    return ColourIndex
end

#########################################################################################################################################################

## Vanish Matching ## 
function vanish_matching(MatchingThreshold, Features)
    # -Matching 

    ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)
    # Initialise Comparison Values and counts from the first row of Features (AllStdW)
    ComparisonValues = copy(Features[1])                     # representative std(w) for each group
    counts = fill(1, length(ComparisonValues))              # how many members assigned to each comparison vector
    FirstRowIndices = collect(1:length(ComparisonValues)) # [1,2,3,...] (index for each group)
    push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 

    # Process subsequent rows
    for row in Features[2:end]
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
                ComparisonValues[MinimumIndex] = CurrentValue # value to compare to on next step 
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
    return ColourIndex
end

#########################################################################################################################################################

# Estimate period of trajectory A
function period(A)
B = [round.(u; digits = 1) for u in A]
C = unique(B)
return length(C)
end

#########################################################################################################################################################