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

