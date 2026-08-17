# Neuron 
function neuron_system!(du, initial_state, Parameters, time)
    v, w = initial_state
    a1, a2, a3, a4, c, ϵ, d, I, p = Parameters
    du[1] = p * v + a1 * v^3 - a2 * v^5 + a3 * v^7 - a4 * v^9 - w + I
    du[2] = ϵ * (1 + c * v + d * v^2 - w)
    return nothing
end