# Task 2.7 -- Discrete time dynamical system
function DTDS(u, p, n)
    x,y = u
    μ, j = p
    xn = y
    yn = μ * y - y^3 - j * x
    return SVector(xn, yn)
end
