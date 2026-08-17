# Featurise a given trajectory 
function featurizer2(trajectory)

    #k = round(Int, Ttr ÷ Dt) # number of steps to discard due to transience -- should be Ttr / Dt

    # extract relevant section of the trajecotry for each component
    #x_temp, y_temp, z_temp = trajectory[:, 1], trajectory[:, 2] # separate by component
    #x, y, z = x_temp[k:end], y_temp[k:end], z_temp[k:end] # ignore section before transient time -- work out if base trajectory ignores this anyway

    x, y, z = trajectory[:, 1], trajectory[:, 2], trajectory[:, 3]   
    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_x = std(x)
    std_y = std(y)
    std_z = std(z)

    return [std_x, std_y, std_z] # feature vector

end

