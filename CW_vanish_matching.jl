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