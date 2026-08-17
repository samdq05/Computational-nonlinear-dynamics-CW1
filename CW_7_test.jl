# General featurise and group 
function fg1(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

    # initialise a features vector (contains feature vectors for each intial condition)
    features = Vector{Vector{Float64}}(undef, length(ics))

    for i in 1:length(ics)
        # trajectory
        neuron = CoupledODEs(ds, ics[i], params)
        N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
        # find the features
        features[i] = featurizer(N)

    end

    return features
end
