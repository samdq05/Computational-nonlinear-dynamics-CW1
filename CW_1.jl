# lorenz63 function for section 1.4
function lorenz63_rule!(du, u, p, t)
    x, y, z = u
    σ, β, ρ = p
    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
    return nothing # always `return nothing` for in-place form!
end