# Featurise a given trajectory 
function featurizer2(trajectory)

    # Extract Components
    x, y, z = columns(trajectory)
    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_x, std_y, std_z = std(x), std(y), std(z)

    return [std_x, std_y, std_z] # feature vector

end

