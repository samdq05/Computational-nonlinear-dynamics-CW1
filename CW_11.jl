# Featurise a given trajectory 
function featurizer2(trajectory)

    x, y, z = trajectory[:, 1], trajectory[:, 2], trajectory[:, 3]   
    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_x = std(x)
    std_y = std(y)
    std_z = std(z)

    return [std_x, std_y, std_z] # feature vector

end

