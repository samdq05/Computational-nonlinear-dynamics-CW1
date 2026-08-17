# lorenz63 function for section 1.4
function DS26(du, u, p, t)
    x, y, z = u
    b = p
    du[1] = sin(y) - b * x 
    du[2] = sin(z) - b * y
    du[3] = sin(x) - b * z
    return nothing # always `return nothing` for in-place form!
end