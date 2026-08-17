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
                else
                    # Safety fallback if lengths differ
                    MatchedIndex = maximum(vcat(ColourIndex...)) + 1
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