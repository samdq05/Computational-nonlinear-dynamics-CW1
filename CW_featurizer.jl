# Featurise a given trajectory 
function featurizer(trajectory)

    # Extract components
    v, w = columns(trajectory)
    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_v = std(v)
    std_w = std(w)

    return [std_v, std_w] # feature vector

end

