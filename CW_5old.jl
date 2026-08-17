function simple_grouping(features, r)

    # separate each tuple into x and y values
    x = Vector{Float64}(undef, 1000)
    y = Vector{Float64}(undef, 1000)
    for i = 1:1000
        x[i], y[i] = features[i][1], features[i][2]
    end

    Group_index = Vector{Int16}(undef, 1000) # empty vector of integers
    # Base case -- set the first group as defined by the first point
    Group_index[1] = 1
    for i in 2:1000
        distance = norm([x[i], y[i]] - [x[1], y[1]])
        if distance < r
            Group_index[i] = 1
        else
            Group_index[i] = 2
        end
    end
    return Group_index # tells us if a given point is in group one or two
end