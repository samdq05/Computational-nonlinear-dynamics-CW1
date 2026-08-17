# Featurizer for Task 2.7
function featurizer3(trajectory)

    divergence = false # assumme a trajectory doesn't diverge 

    x, y = trajectory[:, 1], trajectory[:, 2]
    # Features - mean seems an appropriate distinguishing characteristic
    mean_x = mean(x)
    mean_y = mean(y)
    # Special feature -- for diverging trajectories
    if abs(mean_x) > 10_000 || abs(mean_y) > 10_000 # check if trajectories diverge in any direction
        divergence = true
    end

    return [mean_x, mean_y, divergence] # feature vector
end

