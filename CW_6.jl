# Featurise a given trajectory 
function featurizer(trajectory)

    k = round(Int, 50 ÷ 0.5) # number of steps to discard due to transience -- should be Ttr / Dt

    # extract relevant section of the trajecotry for each component
    v_temp, w_temp = trajectory[:, 1], trajectory[:, 2] # separate by component
    v, w = v_temp[k:end], w_temp[k:end] # ignore section before transient time

    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_v = std(v)
    std_w = std(w)

    return [std_v, std_w] # feature vector

end

