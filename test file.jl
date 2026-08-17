import Pkg
using OrdinaryDiffEqTsit5

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob, Tsit5())

X, t = trajectory(sol, tspan)
##########################################################################

using DynamicalSystems

function henon_rule(u, p, n) # here `n` is "time", but we don't use it.
    x, y = u # system state
    a, b = p # system parameters
    xn = 1.0 - a*x^2 + y
    yn = b*x
    return SVector(xn, yn)
end

u0 = [0.2, 0.3]
p0 = [1.4, 0.3]

henon = DeterministicIteratedMap(henon_rule, u0, p0)

total_time = 10_000
X, t = trajectory(henon, total_time)
X

using CairoMakie
scatter(X)

###########################################################################################################

function 

###############################################################################################################
function lorenz63_rule!(du, u, p, t)
    x, y, z = u
    σ, β, ρ = p
    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
    return nothing # always `return nothing` for in-place form!
end

# 
total_time = 100.0
sampling_time = 0.02
# paramters and initial conditions
p0 = [10.0, 8/3, 28.0]

u1 = [10.0, 2.0, 20.0]
u2 = [-10.0, -2.0, 15.0]
u3 = [4.0, -1.0, 30.0]

# figure
using CairoMakie
using DynamicalSystems

fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "time", ylabel = "variable")
ax.azimuth = 5.3

for u0 in (u1, u2, u3)
# ODE trajectories
lorenz63 = CoupledODEs(lorenz63_rule!, u0, p0)
X, t = trajectory(lorenz63, total_time; Δt = sampling_time)
# test
println(X[1]) # shows initial conditions
println()

x, y, z = X[:,1], X[:,2], X[:,3]
lines!(ax, x, y, z)

end
fig


###############################################

function lorenz63_rule!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = -u[1] * u[3] + p[3] * u[1] - u[2]
    du[3] = u[1] * u[2] - p[2] * u[3]
    return nothing # always `return nothing` for in-place form!
end

# initial state and parameters
u0_1 = [10.0, 2.0, 20.0] # x,y,x
u0_2 = [-10.0, -2.0, 15.0]
u0_3 = [4.0, -1.0, 30.0]
p0 = [10.0, 8/3, 20.0] # σ,β,ρ

# calculate for each initial condition
lorenz63_1 = CoupledODEs(lorenz63_rule!, u0_1, p0)
lorenz63_2 = CoupledODEs(lorenz63_rule!, u0_2, p0)
lorenz63_3 = CoupledODEs(lorenz63_rule!, u0_3, p0)

total_time = 100.0
sampling_time = 0.02

X, t = trajectory(lorenz63_1, total_time; Ttr = 2.2, Δt = sampling_time) # Ttr signifies the transient time to discard (i.e. let the system settle)
Y, t = trajectory(lorenz63_2, total_time; Ttr = 2.2, Δt = sampling_time)
Z, t = trajectory(lorenz63_3, total_time; Ttr = 2.2, Δt = sampling_time)


fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "time", ylabel = "variable")
for var in columns(Y)
    lines!(ax, t, var)
end
fig


# Scatter plot separate
using Makie
figScat = Figure()
scatter(figScat[1,1], X)
scatter(figScat[1,2], Y)
scatter(figScat[2,1:2], Z)

# Scatter plot combined
fig, ax = Figure()
scatter(fig, ax, X, color := red, markersize = 1)
scatter!(fig, ax, Y, color := blue, markersize = 1)
scatter!(fig, ax, Z, color := green, markersize = 1 )
fig

fig = 


## axis 3 example
import Pkg
using CairoMakie
fig = Figure()

Axis3(fig[1, 1], aspect = (1, 1, 1), title = "aspect = (1, 1, 1)")
Axis3(fig[1, 2], aspect = (2, 1, 1), title = "aspect = (2, 1, 1)")
Axis3(fig[2, 1], aspect = (1, 2, 1), title = "aspect = (1, 2, 1)")
Axis3(fig[2, 2], aspect = (1, 1, 2), title = "aspect = (1, 1, 2)")

fig

## same axes
x = range(0, 10; length = 20)
lines(x, sin; label = "sin")  # makes a figure
scatter!(x, cos; label = "cos", color = "purple") # uses last-used figure
current_figure() # return (and hence display) last-used figure
##
## different axes
fig = Figure()
x = range(0, 10; length = 20)
ax1, li1 = lines(fig[1, 1], x, sin)
ax2, sc2 = scatter(fig[2, 1], x, cos; color = "purple")
Legend(fig[1:2, 2], [li1, sc2], ["ημίτονο", "cosinus"])
fig

##
fig = Figure()
ax, sc1 = scatter(fig[1,1], X)
ax, sc2 = scatter!(fig[1,1], Y)
ax, sc3 = scatter!(fig[1,1], Z)
fig

##
fig = Figure()
ax, l1 = lines(X)
ax, l2 = lines!(Y)
ax, l3 = lines!(Z)
fig




################################
# scatter test
x1 = [1,2,3]; y1 = [3,2,1]; z1 = [1,2,1]
x2 = [1,0,3]; y2 = [3,0,1]; z2 = [2,1,2]

using Makie
fig = Figure()
kwargs = (; marker = Rect)
scatter(fig[1,1], x1,y1,z1; kwargs..., markersize = [(10,20), (20,30), (40,30)]) 
scatter(fig[2,1], x2,y2,z2; kwargs..., markersize = 30) 




##############################################################################################

using CairoMakie
using AbstractPlotting.ColorSchemes

function demo()
    xs = LinRange(0.5, 6, 50)
    ys = LinRange(0.5, 6, 50)
    data1 = [sin(x^1.5) * cos(y^0.5) for x in xs, y in ys] .+ 0.1 .* randn.()
    data2 = [sin(x^0.8) * cos(y^1.5) for x in xs, y in ys] .+ 0.1 .* randn.()

    f = Figure(resolution = (1400, 1100), font = "Avenir Light")

    ax, co = contourf(f[1, 1][1, 1], xs, ys, data1, levels = 6,
        axis = (title = "Pyramidal Cells", ylabel = "Coronal Section"))
    contour!(xs, ys, data1, color = :black)

    contourf(f[1, 1][1, 2], xs, ys, data2, levels = 6, axis = (title = "Layer IV Neurons",))
    contour!(xs, ys, data2, color = :black)
    hidedecorations!.([ax, current_axis()])
    ax.ylabelvisible = true

    f[1, 1][1, 1:2, Bottom()] = Label(f, "Sagittal Section", padding = (0, 0, 0, 10))
    f[1, 1][2, 1:2] = Colorbar(f, height = 20, vertical = false, label = "Spike Rate",
        flipaxisposition = false, ticklabelalign = (:center, :top))

    _, sc1 = scatter(f[1, 2][1, 1], randn(100, 2) * [1 5; 3 1],
        color = :red, colorrange = (1, 10),
        axis = (title = "Particle Simulation", xlabel = "Velocity [m/s]",
            ylabel = "Acceleration [m/s²]"))
    sc2 = scatter!(randn(100, 2) * [0.5 -5; 1 0.1],
        color = :blue, colorrange = (1, 10), marker = 'x')

    f[1, 2][1, 2] = Legend(f, [sc1, sc2], ["Gamma", "Beta"], "Particles", framevisible = false)

    f[2, :] = Box(f, color = :gray90)
    f[2, :] = Label(f, "Group Measurements", padding = (0, 0, 5, 5))

    for group in 1:3
        f[3, :][1, group, Top()] = Box(f, color = :gray90)
        f[3, :][1, group, Top()] = Label(f, "Group $group", padding = (0, 0, 5, 5))
        for i in 1:3, j in 1:3
            f[3, :][1, group][i, j] = Axis(f, xticks = LinearTicks(4))
            i < 3 && hidexdecorations!(current_axis(), grid = false)
            j > 1 && hideydecorations!(current_axis(), grid = false)
            for n in 1:3
                lines!(0..10, cumsum(randn(1000)), color = ColorSchemes.Set1_4[n])
            end
        end
    end

    f[0, :] = Label(f, "Makie Complex Plot Demo", textsize = 30)
    
    f
end

################################################################################################



diffeq = (alg = Tsit5(), abstol = 1e-9, reltol = 1e-9)
lorenz63 = ContinuousDynamicalSystem(lorenz63_rule!, u0, p0; diffeq)






#=

# Basic Forward Euler solver to begin testing with
function FE(x0,y0,z0,t0,T,Δt)
    # initialise
    x = zeros(Int(round(T/Δt)+1));
    y = zeros(Int(round(T/Δt)+1));
    z = zeros(Int(round(T/Δt)+1));
    t = zeros(Int(round(T/Δt)+1));
    #initial conditions
    x[1]=x0; y[1]=y0; z[1]=z0; t[1]=t0; i=1;
    while t[i]<T
        x[i+1] = σ*(y[i]-x[i]);
        y[i+1] = -x[i]*y[i] + ρ*x[i] - y[i];
        z[i+1] = x[i]*y[i] - β*z[i];
        t[i+1] = t[i] + Δt; 
        i+=1
    end 
    return x, y, z, t
end
=#


#=

function f(pos_current, t, func)
    func(pos_current..., t)
    return k1
end

function RK4(ini_pos, T, Δt, k1, k2, k3, k4)

    # initial values
    t=zeros(Int(round(T/Δt))); # want it to round down as last step is before T since it stops when t>T

    # L-63 equations
    xdot(x,y,z) = σ*(y-x)
    ydot(x,y,z) = -x*z + ρ*x - y
    zdot(x,y,z) = x*y - β*z

    while t<T

        k1()

        i += 1
        t[i] += Δt
    end 
    return 

end

=#


expression(z) = z == 27.0
A = findall(expression, X)
x,y,z = A[:,1],A[:,2],A[:,2]
scatter!(ax, x, y, z, color = "red")
#


# consider points either side of 27
d = z .- 27.0 
p = [zeros(5001), zeros(5001)]
for i = 2:5001
    if d[i]>d[i-1]
        if d[i]*d[i-1]<0
            p[i] = d[i]
            p[i] = [d[i-1], d[i]] 
        end
    else p[i] = [0.0,0.0]
    end
end


#=
findall(x->x==2, X

x, y, z = ind2sub(X, find(x->x == 2,X))
idx = [x  y  z]
=#
    
###############
#=
fig = Figure()
ax1 = Axis3(fig[1, 1]; xlabel = "time", ylabel = "y", zlabel = "z")
ax.azimuth = 5.3
ax2 = axis(fig[1,2]); xlabel = "x", ylabel = "y") 
=#


#=
# paramters and initial conditions
p0 = [10.0, 8/3, 28.0]
u0 = [10.0,2.0,20.0] # initial condition
T = 1000 # total integration time

# figure
using CairoMakie, Makie, DynamicalSystems 

fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "time", ylabel = "y", zlabel = "z")
ax.azimuth = 5.3

# ODE trajectories
lorenz63 = CoupledODEs(lorenz63_rule!, u0, p0)
X, t = trajectory(lorenz63, total_time; Δt = sampling_time)
x, y, z = X[:,1], X[:,2], X[:,3]
plane_body = Rect3f(Point3f(-30,-30,27), Vec3f(70, 70, 0));
lines!(ax, x, y, z)
mesh!(ax, plane_body, color = (:blue, 0.25))

fig
=#
#################



fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "time", ylabel = "y", zlabel = "z")
lorenz63 = CoupledODEs(lorenz63_rule!, u0, p0)
X, t = trajectory(lorenz63, total_time; Δt = sampling_time)
x, y, z = X[:,1], X[:,2], X[:,3]
fig



function psos1(A, index, value) # trajectory A, 
end

using CairoMakie, DynamicalSystems

# paramters and initial conditions
p0 = [10.0, 8/3, 28.0]
u0 = [10.0,2.0,20.0] # initial condition
T = 1000 # total integration time

lorenz63 = CoupledODEs(lorenz63_rule!, u0, p0)
X, t = trajectory(lorenz63, total_time; Δt = sampling_time)
x, y, z = X[:,1], X[:,2], X[:,3]

# this works
fig = Figure(size = (1000,400))
# 3D axis
ax1 = Axis3(fig[1, 1:3]; title = "Poincaré Section for a trajectory with initial conditions (10,2,20)", titlealign = :center, xlabel = "x", ylabel = "y", zlabel = "z")
ax1.azimuth = 5.3
ax1.elevation = 0.1
# 2D axis
ax2 = Axis(fig[1,4]; title = "Poincaré Section in the z = 27 plane", titlealign = :right, xlabel = "x", ylabel = "y")
# trajectory plot
lines!(ax1, x, y, z, color = (:blue, 0.5))
# key values
n = length(X)
P = [0.0, 0.0, 27.0] # point on z=27 plane
v = [0.0, 0.0, 1.0] # normal vector to plane
# check each point
for i = 2:n
    if z[i] > 27.0 && z[i-1] < 27.0
        B = X[i,:] 
        A = X[i-1,:]
        r = v' * (P - A) / (v' * (B - A)) # return time
        C = A + r * (B - A) # interpolated point
        # plots
        scatter!(ax1, C[1], C[2], C[3], color = "red") # 3D axis
        scatter!(ax2, C[1], C[2], color = "red") # 2D axis
    end
end
fig


# seems like A::float64 will make it read inputs as floats



# Lorenz-63 function
function lorenz63_rule!(du, u, p, t)
    σ, β, ρ = p
    du[1] = p[1] * (u[2] - u[1])
    du[2] = -u[1] * u[3] + p[3] * u[1] - u[2]
    du[3] = u[1] * u[2] - p[2] * u[3]
    return nothing # always `return nothing` for in-place form!
end


##############################
# packages -- NOT SURE THIS IS NECESSARY
#using Makie, CairoMakie, DynamicalSystems, OrdinaryDiffEq, BenchmarkTools
#Pkg.status(["DynamicalSystems", "Makie", "CairoMakie", "OrdinaryDiffEq", "BenchmarkTools"]; mode = Pkg.PKGMODE_MANIFEST)
##############################


############## Poincare Section
# specifically for z=27 plane

using CairoMakie, DynamicalSystems

# paramters and initial conditions
p0 = [10.0, 8/3, 28.0]
u0 = [10.0,2.0,20.0] # initial condition
T = 1000 # total integration time

lorenz63 = CoupledODEs(lorenz63_rule!, u0, p0)
X, t = trajectory(lorenz63, total_time; Δt = sampling_time)
x, y, z = X[:,1], X[:,2], X[:,3]
n = length(X) # index

# this works
fig = Figure(size = (1000,400))
# 3D axis
ax1 = Axis3(fig[1, 1:3]; title = "Poincaré Section for a trajectory with initial conditions (10,2,20)", titlealign = :center, xlabel = "x", ylabel = "y", zlabel = "z")
ax1.azimuth = 5.3
ax1.elevation = 0.1
# 2D axis
ax2 = Axis(fig[1,4]; title = "Poincaré Section in the z = 27 plane", titlealign = :right, xlabel = "x", ylabel = "y")
# trajectory plot
lines!(ax1, x, y, z, color = (:blue, 0.5))
# key values
P = [0.0, 0.0, 27.0] # point on z=27 plane
v = [0.0, 0.0, 1.0] # normal vector to plane
# check each point
for i = 2:n
    if z[i] > 27.0 && z[i-1] < 27.0
        B = X[i,:] 
        A = X[i-1,:]
        r = v' * (P - A) / (v' * (B - A)) # return time
        C = A + r * (B - A) # interpolated point
        # plots
        scatter!(ax1, C[1], C[2], C[3], color = "red") # 3D axis
        scatter!(ax2, C[1], C[2], color = "red") # 2D axis
    end
end
fig

######################################################

#=
function psos1(trajectory, index, value)

    n = length(trajectory) # number of points in the trajectory
    x_store, y_store, z_store = zeros(n), zeros(n), zeros(n) # empty vectors potentially as large as every point
    j = 0 # initialise counter

    # trajecotry components
    x, y, z = trajectory[:,1], trajectory[:,2], trajectory[:,3]
    
    # For a given index (x,y,z = 1,2,3 ) and value -- gives a point on the plane and a normal vector to the plane 
    if index == 1
        # point on plane and normal vector
        P = [value, 0.0, 0.0]
        v = [1.0, 0.0, 0.0]

        for i = 2:n
            if x[i] > value && x[i-1] < value
                j += 1 # add one to count
                B = trajectory[i,:] # point before
                A = trajectory[i-1,:] # point after
                r = dot(v, (P - A)) / dot(v, (B - A)) # return time
                C = A + r * (B - A) # interpolated point
                x_store[j], y_store[j], z_store[j] = C                
            end
        end
    
    elseif index == 2
        # point on plane and normal vector
        P = [0.0, value, 0.0]
        v = [0.0, 1.0, 0.0]

        for i = 2:n
            if y[i] > value && y[i-1] < value
                j += 1 # add one to count
                B = trajectory[i,:] # point before
                A = trajectory[i-1,:] # point after
                r = dot(v, (P - A)) / dot(v, (B - A)) # return time
                C = A + r * (B - A) # interpolated point
                x_store[j], y_store[j], z_store[j] = C
            end
        end

    elseif index == 3
        # point on plane and normal vector
        P = [0.0, 0.0, value]
        v = [0.0, 0.0, 1.0]   

        for i = 2:n
            if z[i] > value && z[i-1] < value
                j += 1 # add one to count
                B = trajectory[i,:] # point before
                A = trajectory[i-1,:] # point after
                r = dot(v, (P - A)) / dot(v, (B - A)) # return time
                C = A + r * (B - A) # interpolated point
                x_store[j], y_store[j], z_store[j] = C
            end
        end
    else
        error("index must be 1, 2, or 3")
    end
    
    x_final, y_final, z_final = x_store[1:j], y_store[1:j], z_store[1:j]
    return x_final, y_final, z_final
end
=#

########################################################################

using CairoMakie, DynamicalSystems

n, m = 6, 6 # n and m don't necessarily have to be the same size
#u = [zeros(n), zeros(m)]
#x, y = zeros(n), zeros(m)
u = zeros(Float64, n*m, 2)

#fig = Figure()
#ax = Axis(fig[:, :]; title="many trajectories", xlabel="v", ylabel="w")

for i = 1:n
    for j = 1:m
        u[i,j] = -3 + (6 / m) * j
        #y[j] = -3 + (6 / m) * j
    end
    #x[i] = -3 + (6 / n) * i
    u[i,j] = -3 + (6 / n) * i
end

u

#u[1] = [x[1],y[1]]
#u[2] = [x[1], y[2]]


# u = collect(zip(x,y))
#u = collect([x,y])
#u = zip(x, y)


#=    
    for u0 in u
    neuron = CoupledODEs(neuron_system!, u0, params)
    N, t = trajectory(neuron, total_time; Δt = sampling_time)
    v, w = N[:,1], N[:,1]
    v_end, w_end = N[end,1], N[end,2]
    lines!(ax, v, w)
    scatter!(ax, v_end, w_end)
    end

    fig
   =#


   ##################################################
   using CairoMakie, DynamicalSystems

n, m = 6, 6 # n and m don't necessarily have to be the same size
#u = [zeros(n), zeros(m)]
x, y = zeros(n), zeros(m)
u = zeros(Float64, 2, n * m)

#fig = Figure()
#ax = Axis(fig[:, :]; title="many trajectories", xlabel="v", ylabel="w")

for i = 1:n
    for j = 1:m
        y[j] = -3 + (6 / m) * j
    end
    x[i] = -3 + (6 / n) * i
end

i = 0
#for i = 1:n*m
    for k = 1:6
        for j = 1:6
            i += 1
            u[:, i] = [x[k], y[j]]
        end
    end
#end


u


#=
u[:, 1] = [x[1], y[1]]
u[:, 2] = [x[1], y[2]]
=#
#u[1] = [x[1],y[1]]
#u[2] = [x[1], y[2]]


# u = collect(zip(x,y))
#u = collect([x,y])
#u = zip(x, y)


#=    
    for u0 in u
    neuron = CoupledODEs(neuron_system!, u0, params)
    N, t = trajectory(neuron, total_time; Δt = sampling_time)
    v, w = N[:,1], N[:,1]
    v_end, w_end = N[end,1], N[end,2]
    lines!(ax, v, w)
    scatter!(ax, v_end, w_end)
    end

    fig
   =#


   for i = 1:n
    x[i] = -3 + (6 / n) * i
    for j = 1:m
        y[j] = -3 + (6 / m) * j
        k += 1
        u[:, k] = [x[i], y[j]]
    end
    
end
u
#=
i = 0
    for k = 1:6
        for j = 1:6
            
        end
    end


u
=#


for i = 1:n*m
    #u1 = u[:,i]
    for u0 in u1#u[:,1]#[u[:,1], u[:,2], u[:,3]]
        neuron = CoupledODEs(neuron_system!, u0, params)
        N, t = trajectory(neuron, total_time; Δt=sampling_time)
        v, w = N[:, 1], N[:, 2]
        v_end, w_end = N[end, 1], N[end, 2]
        lines!(ax, v, w)
        scatter!(ax, v_end, w_end)
    end
    
end
fig



#=
p = zeros(n*m)
u1 = zeros(n*m)
for i = 1:n*m
    p[i] = i
    u1 = u[:,i]
    #println(u)

    for u0 in u1
#println(u0)
    #neuron = CoupledODEs(neuron_system!, u0, params)
    #N, t = trajectory(neuron, total_time; Δt=sampling_time)
    #v, w = N[:, 1], N[:, 2]
    #println(v, w)


    #lines!(ax, v, w)
        #v_end, w_end = N[end, 1], N[end, 2]
        #scatter!(ax, v_end, w_end)
    end

end
=#

typeof(u)
println(u)







using DynamicalSystems
include("CW_4.jl")

#function fg1(ds, featurizer, r, ics; Ttr = 50.0, Dt = 0.5, T = 100.0)

# Temp 
T = 100.0
Dt = 0.5
Ttr = 50.0
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
r = 0.05 # threshold value/ distance

n = length(ics)
Group_index = Vector{Int16}(undef, n)

#=
v_plot = Vector{Vector{Float64}}()
w_plot = Vector{Vector{Float64}}()
#time = Vector
t_plot = Vector{Float64}()
=#

k = round(Int, 50 ÷ 0.5); # number of steps to discard due to transience 


for i = 1:n
    neuron = CoupledODEs(neuron_system!, ics[i], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    v_temp, w_temp = N[:, 1], N[:, 2]
    v, w = v_temp[k:end], w_temp[k:end]
    v_end, w_end = N[end, 1], N[end, 2]

    push!(v_plot, v)
    push!(w_plot, w)

end

#println(v_plot[2])
#typeof(v_plot)    

# key things

# find standard deviation for each trajectory 
std_v = Vector{Float64}(undef, n) # empty vector for standard deviations of v for each trajectory
std_w = Vector{Float64}(undef, n) # " of w "

# First state 
vec_temp_1 = 0
p1 = 0;
distance_2 = 0
p2 = 0;
distance_3 = 0
p3 = 0;

for i = 1:n
    std_v[i] = std(v_plot[i])
    std_w[i] = std(w_plot[i])

    vec_temp = [std_v[i], std_w[i]]

    # First feature
    if p1 == 0
        vec_temp_1 = vec_temp
    end

    # Check if it's in the first group
    distance_1 = norm(vec_temp - vec_temp_1)

    if distance_1 < r # if it's below the threshold add it to that feature, otherwise make a new feature
        Group_index[i] = 1
        if p1 == 0
            vec_temp_1 = [std_v[i], std_w[i]]
            p1 += 1
        end

        # Second feature
        if p2 == 0
            vec_temp_2 = vec_temp
        end

        # If it isn't in the first group, check the second
        distance_2 = norm(vec_temp - vec_temp_2)

    elseif distance_2 < r
        Group_index[i] = 2
        if p2 == 0
            vec_temp_2 = [std_v[i], std_w[i]]
            p2 += 1
        end

        # Third feature
        if p3 == 0
            vec_temp_3 = vec_temp
        end

        # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong)
        distance_3 = norm(vec_temp - vec_temp_3)

    elseif distance_3 < r
        Group_index[i] = 3
        if p3 == 0
            vec_temp_3 = [std_v[i], std_w[i]]
            p3 += 1
        end

    end

end

println(Group_index)
#return Group_index # returns the a vector group indices for the initial conditions
#end



using DynamicalSystems
include("CW_4.jl")

#function fg1(ds, featurizer, r, ics; Ttr = 50.0, Dt = 0.5, T = 100.0)

# Temp 
T = 100.0
Dt = 0.5
Ttr = 50.0
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
r = 0.05 # threshold value/ distance

n = length(ics)
Group_index = Vector{Int16}(undef, n)

#=
v_plot = Vector{Vector{Float64}}()
w_plot = Vector{Vector{Float64}}()
#time = Vector
t_plot = Vector{Float64}()
=#

k = round(Int, 50 ÷ 0.5); # number of steps to discard due to transience 


for i = 1:n
    neuron = CoupledODEs(neuron_system!, ics[i], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    v_temp, w_temp = N[:, 1], N[:, 2]
    v, w = v_temp[k:end], w_temp[k:end]
    #v_end, w_end = N[end, 1], N[end, 2]

    push!(v_plot, v)
    push!(w_plot, w)

end

#println(v_plot[2])
#typeof(v_plot)    

# key things

# find standard deviation for each trajectory 
std_v = Vector{Float64}(undef, n) # empty vector for standard deviations of v for each trajectory
std_w = Vector{Float64}(undef, n) # " of w "

# First state 
vec_temp_1 = 0
p1 = 0;
distance_2 = 0
p2 = 0;
distance_3 = 0
p3 = 0;

for i = 1:n
    std_v[i] = std(v_plot[i])
    std_w[i] = std(w_plot[i])

    vec_temp = [std_v[i], std_w[i]]

    # First feature
    if p1 == 0
        vec_temp_1 = vec_temp
    end

    # Check if it's in the first group
    distance_1 = norm(vec_temp - vec_temp_1)

    if distance_1 < r # if it's below the threshold add it to that feature, otherwise make a new feature
        Group_index[i] = 1
        if p1 == 0
            vec_temp_1 = [std_v[i], std_w[i]]
            p1 += 1
        end

        # Second feature
        if p2 == 0
            vec_temp_2 = vec_temp
        end

        # If it isn't in the first group, check the second
        distance_2 = norm(vec_temp - vec_temp_2)

    elseif distance_2 < r
        Group_index[i] = 2
        if p2 == 0
            vec_temp_2 = [std_v[i], std_w[i]]
            p2 += 1
        end

        # Third feature
        if p3 == 0
            vec_temp_3 = vec_temp
        end

        # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong)
        distance_3 = norm(vec_temp - vec_temp_3)

    elseif distance_3 < r
        Group_index[i] = 3
        if p3 == 0
            vec_temp_3 = [std_v[i], std_w[i]]
            p3 += 1
        end

    end

end

println(Group_index) # seems to be working, now i just need to make this a function
#return Group_index # returns the a vector group indices for the initial conditions
#end

using DynamicalSystems, Statistics, LinearAlgebra
include("CW_4.jl")

#function fg1(ds, featurizer, r, ics; Ttr = 50.0, Dt = 0.5, T = 100.0)

# Temp 
T = 100.0
Dt = 0.5
Ttr = 50.0
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
r = 0.05 # threshold value/ distance

n = length(ics)
Group_index = Vector{Int16}(undef, n)


v_plot = Vector{Vector{Float64}}()
w_plot = Vector{Vector{Float64}}()
#time = Vector
t_plot = Vector{Float64}()


k = round(Int, 50 ÷ 0.5); # number of steps to discard due to transience 


for i = 1:n
    neuron = CoupledODEs(neuron_system!, ics[i], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    v_temp, w_temp = N[:, 1], N[:, 2]
    v, w = v_temp[k:end], w_temp[k:end]
    #v_end, w_end = N[end, 1], N[end, 2]

    push!(v_plot, v)
    push!(w_plot, w)

end

#println(v_plot[2])
#typeof(v_plot)    

# key things

# find standard deviation for each trajectory 
std_v = Vector{Float64}(undef, n) # empty vector for standard deviations of v for each trajectory
std_w = Vector{Float64}(undef, n) # " of w "

# First state 
vec_temp_1 = 0
p1 = 0;
distance_2 = 0
p2 = 0;
distance_3 = 0
p3 = 0;

for i = 1:n
    std_v[i] = std(v_plot[i])
    std_w[i] = std(w_plot[i])

    vec_temp = [std_v[i], std_w[i]]

    # First feature
    if p1 == 0
        vec_temp_1 = vec_temp
    end

    # Check if it's in the first group
    distance_1 = norm(vec_temp - vec_temp_1)

    if distance_1 < r # if it's below the threshold add it to that feature, otherwise make a new feature
        Group_index[i] = 1
        if p1 == 0
            vec_temp_1 = [std_v[i], std_w[i]]
            p1 += 1
        end

        # Second feature
        if p2 == 0
            vec_temp_2 = vec_temp
        end

        # If it isn't in the first group, check the second
        distance_2 = norm(vec_temp - vec_temp_2)

    elseif distance_2 < r
        Group_index[i] = 2
        if p2 == 0
            vec_temp_2 = [std_v[i], std_w[i]]
            p2 += 1
        end

        # Third feature
        if p3 == 0
            vec_temp_3 = vec_temp
        end

        # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong)
        distance_3 = norm(vec_temp - vec_temp_3)

    elseif distance_3 < r
        Group_index[i] = 3
        if p3 == 0
            vec_temp_3 = [std_v[i], std_w[i]]
            p3 += 1
        end

    end

end

println(Group_index) # seems to be working, now i just need to make this a function
#return Group_index # returns the a vector group indices for the initial conditions
#end


using DynamicalSystems, Statistics
include("CW_4.jl")

#function fg1(ds, featurizer, r, ics; Ttr = 50.0, Dt = 0.5, T = 100.0)

# Temp 
T = 100.0
Dt = 0.5
Ttr = 50.0
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
r = 0.05 # threshold value/ distance

n = length(ics)
Group_index = Vector{Int16}(undef, n)


#v_plot = Vector{Vector{Float64}}()
#w_plot = Vector{Vector{Float64}}()

k = round(Int, 50 ÷ 0.5); # number of steps to discard due to transience 

# find standard deviation for each trajectory 
std_v = Vector{Float64}(undef, n) # empty vector for standard deviations of v for each trajectory
std_w = Vector{Float64}(undef, n) # " of w "

# First state 
vec_temp_1 = 0
p1 = 0;
distance_2 = 0
p2 = 0;
distance_3 = 0
p3 = 0;

for i = 1:n
    neuron = CoupledODEs(neuron_system!, ics[i], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    v_temp, w_temp = N[:, 1], N[:, 2] # separate by component
    v, w = v_temp[k:end], w_temp[k:end] # ignore section before transient time

    # prob remove
    #v_end, w_end = N[end, 1], N[end, 2]

    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_v[i] = std(v)
    std_w[i] = std(w)

    # can probably remove later
    # push!(v_plot, v)
    # push!(w_plot, w)

    # Assign groups based on features
    vec_temp = [std_v[i], std_w[i]]

    # First feature group #
    if p1 == 0 # if this group is empty, automatically fill it an set as group centre
        vec_temp_1 = vec_temp
    end

    # Check if it's in the first group
    distance_1 = norm(vec_temp - vec_temp_1)

    if distance_1 < r # if it's below the threshold add it to that feature, otherwise make a new feature
        Group_index[i] = 1
        if p1 == 0
            vec_temp_1 = [std_v[i], std_w[i]]
            p1 += 1
        end

        # Second feature group #
        if p2 == 0 # if this group is empty, automatically fill it an set as group centre
            vec_temp_2 = vec_temp
        end

        distance_2 = norm(vec_temp - vec_temp_2) # If it isn't in the first group, check the second
    elseif distance_2 < r
        Group_index[i] = 2
        if p2 == 0
            vec_temp_2 = [std_v[i], std_w[i]]
            p2 += 1
        end

        # Third feature group #
        if p3 == 0 # if this group is empty, automatically fill it an set as group centre
            vec_temp_3 = vec_temp
        end

        distance_3 = norm(vec_temp - vec_temp_3) # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong - there should only be three groups)

    elseif distance_3 < r
        Group_index[i] = 3
        if p3 == 0
            vec_temp_3 = [std_v[i], std_w[i]]
            p3 += 1
        end

    end
end

println(Group_index)

# featurizer

using DynamicalSystems, Statistics
include("CW_4.jl")

function featurizer(trajectory)

    #std_v = Vector{Float64}(undef, n) # empty vector for standard deviations of v for each trajectory
    #std_w = Vector{Float64}(undef, n) # " of w "

    # extract relevant section of the trajecotry for each component
    v_temp, w_temp = trajectory[:, 1], trajectory[:, 2] # separate by component
    v, w = v_temp[k:end], w_temp[k:end] # ignore section before transient time

    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_v = std(v)
    std_w = std(w)

    return feature_vec = [std_v, std_w]

end


# Test

# params and initial conditions
T = 100.0
Dt = 0.5
Ttr = 50.0
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
r = 0.05 # threshold value/ distance
n = length(ics)

# initialise a features vector (contains feature vectors for each intial condition)
features = Vector{Vector{Float64}}(undef, n)

for i in 1:n
    # trajectory
    neuron = CoupledODEs(neuron_system!, ics[i], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    # find the features
    features[i] = featurizer(N)
end

println(features)



using DynamicalSystems, Statistics
include("CW_6.jl")

function fg1(ds, featurizer, r, ics; Ttr = 50.0, Dt = 0.5, T = 100.0)

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


# Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.05
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]] # begin with this for testing - the other one could cause problems unless it's all working
ds = neuron_system!

#=
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

n = length(v)
m = length(w)
k = 0
=#

#= use when confident the res is working
ics = Vector{Vector{Float64}}(undef, n*m)

for i in 1:n
    for j in 1:m
        k += 1 
        ics[k] = [v[i], w[j]]
    end
end
=#
#println(ics[1], ics[2], ics[end])
#length(ics)

fg1(ds, featurizer, r, ics; Ttr = 50, Dt = 0.5, T = 100)



v = collect(-2.5:0.5:2.5)
println(v)
w = collect(-2.5:0.5:2.5)
# start by testing for a small number, at the minute, I think I nedd to move the std into first for loop, or combine whole thing

# want to use a total of 100 values in each range to give a 100x100 area so 10000 initial conditions

#########################################################

using DynamicalSystems, Statistics
include("CW_4.jl")

#function fg1(ds, featurizer, r, ics; Ttr = 50.0, Dt = 0.5, T = 100.0)

n = length(ics)
Group_index = Vector{Int16}(undef, n)

k = round(Int, 50 ÷ 0.5); # number of steps to discard due to transience 

# find standard deviation for each trajectory 
std_v = Vector{Float64}(undef, n) # empty vector for standard deviations of v for each trajectory
std_w = Vector{Float64}(undef, n) # " of w "

# First state 
vec_temp_1 = 0
p1 = 0;
distance_2 = 0
p2 = 0;
distance_3 = 0
p3 = 0;

for i = 1:n
    neuron = CoupledODEs(neuron_system!, ics[i], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    v_temp, w_temp = N[:, 1], N[:, 2] # separate by component
    v, w = v_temp[k:end], w_temp[k:end] # ignore section before transient time

    # Features -- standard deviations seems an appropriate distinguishing characteristic
    std_v[i] = std(v)
    std_w[i] = std(w)

    # Assign groups based on features
    vec_temp = [std_v[i], std_w[i]]

    # First feature group #
    if p1 == 0 # if this group is empty, automatically fill it an set as group centre
        vec_temp_1 = vec_temp
    end

    # Check if it's in the first group
    distance_1 = norm(vec_temp - vec_temp_1)

    if distance_1 < r # if it's below the threshold add it to that feature, otherwise make a new feature
        Group_index[i] = 1
        if p1 == 0
            vec_temp_1 = [std_v[i], std_w[i]]
            p1 += 1
        end

        # Second feature group #
        if p2 == 0 # if this group is empty, automatically fill it an set as group centre
            vec_temp_2 = vec_temp
        end

        distance_2 = norm(vec_temp - vec_temp_2) # If it isn't in the first group, check the second
    elseif distance_2 < r
        Group_index[i] = 2
        if p2 == 0
            vec_temp_2 = [std_v[i], std_w[i]]
            p2 += 1
        end

        # Third feature group #
        if p3 == 0 # if this group is empty, automatically fill it an set as group centre
            vec_temp_3 = vec_temp
        end

        distance_3 = norm(vec_temp - vec_temp_3) # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong - there should only be three groups)

    elseif distance_3 < r
        Group_index[i] = 3
        if p3 == 0
            vec_temp_3 = [std_v[i], std_w[i]]
            p3 += 1
        end

    end
end

return Group_index

#end

#
# Temp 
T = 100.0
Dt = 0.5
Ttr = 50.0
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
r = 0.05 # threshold value/ distance

#####################################

using DynamicalSystems, Statistics
include("CW_4.jl"), include("CW_6.jl"), include("CW_7.jl")

# Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.05
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
#ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]] # begin with this for testing - the other one could cause problems unless it's all working
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

n = length(v)
m = length(w)
k = 0

ics = Vector{Vector{Float64}}(undef, n * m)

for i in 1:n
    for j in 1:m
        k += 1
        ics[k] = [v[i], w[j]]
    end
end
#println(ics[1], ics[2], ics[end])
#length(ics)

features = fg1(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

println(features)

#########################################


# featurizer

using DynamicalSystems, Statistics
include("CW_4.jl")
include("CW_6.jl")

# Test

# params and initial conditions
T = 100.0
Dt = 0.5
Ttr = 50.0
ics = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
r = 0.05 # threshold value/ distance
n = length(ics)

# initialise a features vector (contains feature vectors for each intial condition)
features = Vector{Vector{Float64}}(undef, n)

for i in 1:n
    # trajectory
    neuron = CoupledODEs(neuron_system!, ics[i], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    # find the features
    features[i] = featurizer(N)
end

println(features)
#typeof(features)




#=
    # initial trajectory
    neuron = CoupledODEs(ds, ics[1], params)
    N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
    v_temp, w_temp = N[:, 1], N[:, 2] # separate by component
    v, w = v_temp[k:end], w_temp[k:end] # ignore section before transient time
    u1 = [v, w]
    =#


    # General featurise and group function 
function fg1(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

    ## Featurise ##

    # initialise a features vector (contains feature vectors for each intial condition)
    features = Vector{Vector{Float64}}(undef, length(ics))

    for i in 1:length(ics)
        # trajectory
        neuron = CoupledODEs(ds, ics[i], params)
        N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
        # find the features
        features[i] = featurizer(N)
    end

    ## Group ## 

    Group_index = Vector{Int16}(undef, length(ics))

    k = round(Int, 50 ÷ 0.5) # number of steps to discard due to transience 

    # First state 
    vec_temp_1 = 0 # check if these can be removed
    distance_1 = 0
    p1 = 0
    vec_temp_2 = 0
    distance_2 = 0
    p2 = 0
    vec_temp_3 = 0
    distance_3 = 0
    p3 = 0

    for i = 1:length(ics)

        # Assign groups based on features
        std_v = features[i][1]
        std_w = features[i][2]
        vec_temp = [std_v, std_w]

        # First feature group #
        if p1 == 0 # if this group is empty, automatically fill it an set as group centre
            vec_temp_1 = vec_temp
        end

        # Check if it's in the first group
        distance_1 = norm(vec_temp - vec_temp_1)

        if distance_1 < r # if it's below the threshold add it to that feature, otherwise make a new feature
            Group_index[i] = 1
            if p1 == 0
                vec_temp_1 = [std_v, std_w]
                p1 += 1
            end

            # Second feature group #
            if p2 == 0 # if this group is empty, automatically fill it an set as group centre
                vec_temp_2 = vec_temp
            end

            distance_2 = norm(vec_temp - vec_temp_2) # If it isn't in the first group, check the second

        elseif distance_2 < r
            Group_index[i] = 2
            if p2 == 0
                vec_temp_2 = [std_v, std_w]
                p2 += 1
            end

            # Third feature group #
            if p3 == 0 # if this group is empty, automatically fill it an set as group centre
                vec_temp_3 = vec_temp
            end

            distance_3 = norm(vec_temp - vec_temp_3) # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong - there should only be three groups)

        elseif distance_3 < r
            Group_index[i] = 3
            if p3 == 0
                vec_temp_3 = [std_v, std_w]
                p3 += 1
            end

        end
    end

    return Group_index

end

using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl")

# General featurise and group function 
function fg1(ds, featurizer, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

    ## Featurise ##
#=
    # initialise a features vector (contains feature vectors for each intial condition)
    features = Vector{Vector{Float64}}(undef, length(ics))

    for i in 1:length(ics)
        # trajectory
        neuron = CoupledODEs(ds, ics[i], params)
        N, t = trajectory(neuron, T; Ttr, Δt=Dt) # Discarding Ttr seconds of transient time
        # find the features
        features[i] = featurizer(N)
    end
=#
# temp --
features = [[0.0,0.0], [0.1,0.0], [1.0, 1.0], [1.0,1.1], [2.0,2.0], [2.0,2.1]]
ics=[1,2,3,4,5,6]

    ## Group ## 

    Group_index = Vector{Int16}(undef, length(ics))

    k = round(Int, 50 ÷ 0.5) # number of steps to discard due to transience 

    # First state 
    Group_centre_1 = 0 # check if these can be removed
    distance_1 = 0
    p1 = 0
    Group_centre_2 = 0
    distance_2 = 0
    p2 = 0
    Group_centre_3 = 0
    distance_3 = 0
    p3 = 0

    failed = 0
    debug1 = 0
    debug2 = 0
    debug3 = 0
    debug4 = 0

    J0 = Vector{Int16}(undef, 3) # vector of the indices of group centres

    for i = 1:length(ics)

        # Assign groups based on features
        std_v = features[i][1]
        std_w = features[i][2]
        Comparison_vector = [std_v, std_w]

        # First feature group #
        if p1 == 0 # if this group is empty, automatically fill it an set as group centre
            Group_centre_1 = Comparison_vector
        end

        # Check if it's in the first group
        distance_1 = norm(Comparison_vector - Group_centre_1)

        if distance_1 < r # if it's below the threshold add it to that feature, otherwise make a new feature
            Group_index[i] = 1
            if p1 == 0
                Group_centre_1 = [std_v, std_w]
                p1 += 1
                J1 = i # index of the group centre
                J0[1] = J1
            end

            # Second feature group #
            if p2 == 0 # if this group is empty, automatically fill it an set as group centre
                Group_centre_2 = Comparison_vector
            end

            distance_2 = norm(Comparison_vector - Group_centre_2) # If it isn't in the first group, check the second

        elseif distance_2 < r
            Group_index[i] = 2
            if p2 == 0
                Group_centre_2 = [std_v, std_w]
                p2 += 1
                J2 = i # index of the group centre
                J0[2] = J2
            end

            if p2 == 1
                debug1 = distance_2
            end

            # Third feature group #
            if p3 == 0 # if this group is empty, automatically fill it an set as group centre
                Group_centre_3 = Comparison_vector
            end

            distance_3 = norm(Comparison_vector - Group_centre_3) # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong - there should only be three groups)

        elseif distance_3 < r
            Group_index[i] = 3
            if p3 == 0
                Group_centre_3 = [std_v, std_w]
                p3 += 1
                J3 = i # index of the group centre
                J0[3] = J3
            end

                        if p2 == 1
                debug2 = distance_1
                debug3 = distance_2
                debug4 = distance_3
            end

        else 
            failed += 1

        end
    end

    return Group_index, J0, failed, debug1, debug2, debug3, debug4

end

# Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

l1 = length(v)
l2 = length(w)
K = 0

ics = Vector{Vector{Float64}}(undef, l1 * l2)

for i in 1:l1
    for j in 1:l2
        K += 1
        ics[K] = [v[i], w[j]]
    end
end

indices, J0, failed, d1, d2, d3, d4 = fg1(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

#length(indices) # temp way to indicate it's working
#println(J0)
#=
# Plot the trajectories representing the group centres
fig = Figure()
ax = Axis(fig[:, :]; title="", xlabel="", ylabel="")

v, w = Vector{Vector{Float64}}(undef, 3), Vector{Vector}(undef, 3) # initialise empty vectors of vectors
count = 0
T = 100
Ttr = 50
Dt = 0.5

println(J0)
println(failed)


for J in J0
    count += 1
    println(J)
    println(ics[J])
    neuron = CoupledODEs(neuron_system!, ics[J], params)
    N, t = trajectory(neuron, T; Ttr=50, Δt=Dt) # Discarding Ttr seconds of transient time
    v_temp, w_temp = N[:, 1], N[:, 2] # separate by component
    v[count], w[count] = v_temp, w_temp   
end

lines!(ax, v[1], w[1], color="blue")
lines!(ax, v[2], w[2], color="red")
lines!(ax, v[3], w[3], color="green")

fig
=#
#=

#println(ics[5101])
#println(indices[5101])
println(J0)
println(d1)
println(d2)
println(d3)
println(d4)
=#

println(indices)


####################################################################################################################################################

####################################################################################################################################################

####################################################################################################################################################

##### This section involves useful sections that previously worked very well so could be worth referring back to in the event that something in Task 2 breaks

####################################################################################################################################################

### This one works - YAY! ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_7.jl")

#Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

l1 = length(v)
l2 = length(w)
K = 0

ics = Vector{Vector{Float64}}(undef, l1 * l2)

for i in 1:l1
    for j in 1:l2
        K += 1
        ics[K] = [v[i], w[j]]
    end
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

# want a plot of the ics colour-coded based on their index
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:,:], limits = (-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
for i in 1:length(ics)
    v, w = ics[i] 
    if GroupIndex[i] == 1 
        scatter!(BasinAxis, v, w, color = "LightBlue", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, v, w, color = "PeachPuff", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, v, w, color = "PaleGreen", marker = :rect, markersize = 10)
    end
end
BasinsFigure

#=
# Tells me how many times each index occurs
GroupDict = Dict{Int16, Int16}()

# count each index
for i in 1:length(GroupIndex)
    if haskey(GroupDict, GroupIndex[i])
        GroupDict[GroupIndex[i]] += 1
    else
        GroupDict[GroupIndex[i]] = 1
    end
end
println(GroupDict) # print how many 1s are in group index
=#

####################################################################################################################################################

include("fg3_almostfinal.jl") # causes a failure at some point
using DynamicalSystems, LinearAlgebra, Statistics, CairoMakie
include("CW_4.jl"), include("CW_6.jl")#, include("CW_9_old.jl")

## Test fg3 ##

#Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

l1 = length(v)
l2 = length(w)
K = 0

ics = Vector{Vector{Float64}}(undef, l1 * l2)

for i in 1:l1
    for j in 1:l2
        K += 1
        ics[K] = [v[i], w[j]]
    end
end

# Call the function
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

println(FractionsBoA) 
println(GroupIndex)

####################################################################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
#include("CW_10.jl"), include("Featurizer2.6.jl"), include("T2.6_for_CW_9.jl") # Once finalised, make this CW_9
include("CW_10.jl"), include("CW_11.jl"), include("fg3_almostfinal.jl") #include("T2.7_for_CW_9._2.jl")

# Function Rule
ds = DS26

# parameters and intial conditions
Parameters = 0.16
# Fixed set of initial conditions
ics = Vector{Vector{Float64}}()
# temp set for testing
#yg = xg = zg = collect(range(-3, 3; length=4))
#zg .+= 0.1

yg = xg = zg = collect(range(-3, 3; length = 20))
zg .+= 0.1

for x in xg
    for y in yg
        for z in zg
            ics_temp = [x, y, z]
            push!(ics, ics_temp)
        end
    end
end

# Time Values
StartTime = 0.0
TransientTime = 0.0
TotalTime = 100.0
SamplingTime = 0.2

# threshold value
r = 0.2

GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer2, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

#println(FractionsBoA)

### this newer version of fg3 agrees with the previous one

# Plot the system attractors
ThreeDimensionaSystemFigure = Figure()
TDSAxis = Axis3(ThreeDimensionaSystemFigure[:,:]; title = "System Attractors", xlabel = "x", ylabel = "y", zlabel = "z")
TDSAxis.azimuth = 5.3 # axis perspective
colours = ["blue", "red", "green", "pink"]
for i in 1:4 # change to a length at some point
    x, y, z  = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2], GroupCentresTrajectories[i][3]
    scatter!(TDSAxis, x[end], y[end], z[end], color = colours[i], marker = :star6, markersize = 20)
end
ThreeDimensionaSystemFigure

####################################################################################################################################################
#=
### 2.7 -- Discrete Dynamical System

So I am adding a special feature for divergence. But the divergence condition can't be consider in the norm calcs, so if it is true, it need to make some assumptions, perhaps I could add a divergence vector which tracks if a given initial condition diverges and compare subsequent conditions


Idea 1:
1) I make a vector of bools that tracks the index for every single function, whether or not it's divergent -- true or false

2) To generalise, I will need to have to somehow check if we have features which are bools or only floats, if it is only floats, then I can ignore this new section and proceed as usual, otherwise move on to the next step

3) If its a function with a divergence feature, then I need to consider a couple things: 
        -firstly if the first element is divergent, as before it needs to be the first group, but making it the group centre is meaningless, as it's interaction with the norm function will always produce norm, but I want it grouped with other divergent terms, so if I keep a Vector{Bool} I can check if the current Bool is true and put it into the divergence group given by some indicator

Idea 2:
1) Scrap the divergent feature, since it's divergence will be infinite if its mean is +-infinity so we could just keep that as a feature to the check (again, problems arise if the feature isn' something like mean as std gives NaN instead of Inf, then again, both would give false if compared to being smaller than a value as one is Inf and the other is NaN)

2) Continues similarly to the one above

Idea:
1) Explain the process
=#

############################################################################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("CW_10.jl"), include("Featurizer2.6.jl"), include("T2.6_for_CW_9.jl") # Once finalised, make this CW_9

# Function Rule
ds = DS26

# parameters and intial conditions
Parameters = 0.16
# Fixed set of initial conditions
ics = Vector{Vector{Float64}}()
# temp set for testing
#yg = xg = zg = collect(range(-3, 3; length=4))
#zg .+= 0.1

yg = xg = zg = collect(range(-3, 3; length = 20))
zg .+= 0.1

for x in xg
    for y in yg
        for z in zg
            ics_temp = [x, y, z]
            push!(ics, ics_temp)
        end
    end
end

# Time Values
StartTime = 0.0
TransientTime = 0.0
TotalTime = 100.0
SamplingTime = 0.2

# threshold value
r = 0.2

GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer2, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

println(FractionsBoA)

###################################################################################################################################################################

using CairoMakie, DynamicalSystems, Statistics

include("CW_4.jl")

n, m = 6, 6 # n and m don't necessarily have to be the same size
x, y = zeros(n), zeros(m)
u = [Vector{Float64}(undef, 2) for _ in 1:n*m]
k = 0


fig = Figure()
ax = Axis(fig[:, :]; title="many trajectories", xlabel="v", ylabel="w")

increment_size_x = 6 / n
increment_size_y = 6 / m
start_x = -2.5
start_y = -2.5

for i = 1:n
    multiplier = i - 1
    x[i] = start_x + increment_size_x * multiplier
    #x[i] = -3 + (6 / n) * (i - 1)
    for j = 1:m
        multiplier = j - 1
        y[j] = start_y + increment_size_y * multiplier
        #y[j] = -2.5 + (6 / m) * (j - 1) 
        k += 1
        u[k] = [x[i], y[j]]
    end
end

params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
total_time = 100.0
sampling_time = 0.01
start_time = 0.0

v_plot = Vector{Vector{Float64}}()
w_plot = Vector{Vector{Float64}}()
#time = Vector
t_plot = Vector{Float64}()

time = zeros(n*m)
for i = 1:n*m
    neuron = CoupledODEs(neuron_system!, u[i], params)
    N, t = trajectory(neuron, total_time; Δt=sampling_time) # Discarding 2.5s of transient time
    v, w = N[:, 1], N[:, 2]
    v_end, w_end = N[end, 1], N[end, 2]
    lines!(ax, v, w, color=(:gray, 0.25))
    scatter!(ax, v_end, w_end, color="red", markersize=10)#, marker = :rect)
    #type = typeof(v)
    #print(type)
    push!(v_plot, v)
    push!(w_plot, w)
     
    
    #type = typeof(t)
     #println(type)

     #println(t[1])
#t_plot = t[i]
     #push!(t_plot, t)
end
#fig
#println(v_plot[1])
#println(t_plot[1])

#typeof(t_plot)
time = collect(start_time:sampling_time:total_time)
typeof(time) # Vector{Float64}
typeof(v_plot) # Vector{Vector{Float64}} -- v_plot contains all 36 trajectories
typeof(v_plot[1]) # Vector{Float64}
#println(time)

length(v_plot[1]) # same length
length(time) # "

println(time)
println(v_plot[1])
println(v_plot[36])

fig_test =  Figure()
ax_test = Axis(fig_test[1,:]; limits=((0, 100), nothing), title = "v time series for trajectories", xlabel = "t", ylabel = "v")
ax2 = Axis(fig_test[2,:]; limits=((0, 100), nothing), title = "w time series for trajectories", xlabel = "t", ylabel = "w")
# top
lines!(ax_test, time, v_plot[1], color = "green")
lines!(ax_test, time, v_plot[22], color = "red")
lines!(ax_test, time, v_plot[20], color = "blue")
# bottom
lines!(ax2, time, w_plot[1], color = "green")
lines!(ax2, time, w_plot[22], color = "red")
lines!(ax2, time, w_plot[20], color = "blue")
fig_test

# Find Features
std_v1 = std(v_plot[1])
println(std_v1)
#mean(x)

# key things
r = 0.05 # threshold value/ distance
Num = 1 # number of groups 
# find standard deviation for each trajectory 
std_v = Vector{Float64}(undef, n*m) # empty vector for standard deviations of v for each trajectory
std_w = Vector{Float64}(undef, n*m)# " of w "
vec_group = Vector{Vector{Float64}}(undef, 3)
# base case
std_v[1] = std(v_plot[1])
std_w[1] = std(w_plot[1])
for i = 2:n*m
    std_v[i] = std(v_plot[i])
    std_w[i] = std(w_plot[i])
    vec_group[i] = [std_v[1], std_w[1]]
    vec_temp_current = [std_v[i], std_w[i]]
    #distance = norm(vec_temp_2 - vec_temp_1)
    #if distance < r # if it's below the threshold add it to that feature, otherwise make a new feature
    #    label =      
    #else 
     #   Num += 1
      #  vec_
end
# the next step is to start putting the values into groups and colour code them 

println(std_v)
typeof(std_v)

fig_feature = Figure()
ax = Axis(fig_feature[:,:]; title = "Feature Space", xlabel = "std_v", ylabel = "std_w")
scatter!(ax, std_v, std_w, color = " red")
fig_feature

##########################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("fg3_almostfinal2_for2.7.jl"), include("CW_12.jl"), include("CW_13.jl")
#include("Featurizer2.7.jl"), include("DiscreteTimeDynamicalSystem.jl"), include("T2.7_for_CW_9.jl")

# Initialise relevant functions, parameters and initial conditions
ds = DTDS
r = 0.2
Parameters = [2.9, 0.66]
# Time Values
StartTime = 0.0
TransientTime = 200.0
TotalTime = 50.0
SamplingTime = 1.0
#=
ics = Vector{Vector{Float64}}()
xg = yg = collect(range(-2.5, 2.5; length = 100))#1000))
for x in xg
    for y in yg
            ics_temp = [x, y]
            push!(ics, ics_temp)
    end
end=#
xg = yg = range(-2.5, 2.5; length = 100)#1000)
ics = [[x, y] for x in xg for y in yg ]


GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

#println(FractionsBoA)
#println(GroupIndex)


### Plot basin
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:,:], limits = (-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
for i in 1:length(ics)
    x, y = ics[i] 
    if GroupIndex[i] == 1 
        s1 = scatter!(BasinAxis, x, y, color = "LightBlue", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 2
        s2 = scatter!(BasinAxis, x, y, color = "PeachPuff", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 3
        s3 = scatter!(BasinAxis, x, y, color = "PaleGreen", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 4
        s4 = scatter!(BasinAxis, x, y, color = "LightPink", marker = :rect, markersize = 10)
    end
end
Legend
BasinsFigure 

##########################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("fg3_almostfinal.jl"), include("CW_12.jl"), include("CW_13.jl")
#include("Featurizer2.7.jl"), include("DiscreteTimeDynamicalSystem.jl"), include("T2.7_for_CW_9.jl")

# Initialise relevant functions, parameters and initial conditions
ds = DTDS
r = 0.2
Parameters = [2.9, 0.66]
# Time Values
StartTime = 0.0
TransientTime = 200.0
TotalTime = 50.0
SamplingTime = 1.0

xg = yg = range(-2.5, 2.5; length = 100)#1000)
ics = [[x, y] for x in xg for y in yg ]

Trajectories = Vector{StateSpaceSet{2, Float64, SVector{2, Float64}, Nothing}}()
for InitialConditions in ics
DynamicalSystem = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
Trajectory, t = trajectory(DynamicalSystem, TotalTime)
push!(Trajectories, Trajectory)
end
#println(Trajectories[1])
#println(ics[10000])
#ics[10000]
l1 = featurizer3(Trajectories[10000]) 
l2 = featurizer3(Trajectories[1])
l1[3] == l2[3]

#=
TestFigure2 = Figure()
TestAxis2 = Axis(TestFigure2[:,:]; limits = (-10.0, 10.0, -10.0, 10.0), title = "2.7 Test Trajectory", xlabel = "x", ylabel = "y")
for i = 1: length(Trajectories)
scatter!(TestAxis2, Trajectories[i][:,1][end], Trajectories[i][:,2][end], color = "blue")
end
TestFigure2=#
#GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

##########################################################################################################################################################

using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")


# Initialise relevant functions, parameters and initial conditions

values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0]

# Time Values
StartTime = 0.0
TransientTime = 400.0
TotalTime = 100.0
SamplingTime = 0.5

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
push!(ds, ODE)
end

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

GroupCentresTrajectories

####################################################################################################################################

v = []
w = []

for i in 1:3 # change to a length at some point
v, w = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2]
end

typeof(v)
v[1]
print(v)
max(v)

MaxV = max(v)
MaxW = max(w)
StdW = std(w)

#####################################################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("CW_13.jl")
include("CW_12.jl")
include("fg3_FORALLOFTHEM.jl")

# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 100)#1000)
ics = [[x, y] for x in xg for y in yg ]

# Initialise relevant functions, parameters and initial conditions
r = 1.5
Parameters = [2.9, 0.66]
StartTime = 0
TransientTime = 200
TotalTime = 50
SamplingTime = 1

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
push!(ds, ODE)
end

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

#################################################################################################################################

yg = xg = zg = collect(range(-3, 3; length = 20))
zg .+= 0.1
ics = [[x, y, z] for x in xg for y in yg for z in zg]

#########################################################################################################################################

yg = xg = zg = (range(-3, 3; length = 20))
ics = [[x, y, z + 0.1] for x in xg for y in yg for z in zg]

######################################################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("CW_10.jl"), include("CW_11.jl"), include("fg3_FORALLOFTHEM.jl")

# Time Values
StartTime = 0.0
TransientTime = 0.0
TotalTime = 100.0
SamplingTime = 0.5
# Threshold Value and Parameters
r = 0.2
Parameters = 0.16
# Initial Conditions
yg = xg = zg = (range(-3, 3; length = 20))
ics = [[x, y, z + 0.1] for x in xg for y in yg for z in zg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(DS26, InitialConditions, Parameters)
push!(ds, ODE)
end

#=
fig = Figure()
ax = Axis3(fig[:,:])
for DS in ds
            Trajectory, t = trajectory(DS, T=TotalTime; Ttr=TransientTime, Dt=SamplingTime, t0 = StartTime)  
            x, y, z = Trajectory[:,1], Trajectory[:,2], Trajectory[:,3]
scatter!(ax, x[end], y[end], z[end])
end
fig
=#
#=
figFeat = Figure()
AxFeat = Axis3(figFeat[:,:])
for DS in ds
            Trajectory, t = trajectory(DS, TotalTime; Ttr=TransientTime, Dt=SamplingTime, t0 = StartTime)  
            x, y, z = Trajectory[:,1], Trajectory[:,2], Trajectory[:,3]
    scatter!(AxFeat, std(x), std(y), std(z))
end
    figFeat=#

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer2, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)

println(GroupIndex)

# Plot the system attractors -- set it so that the labels are the FractionsBoA values
ThreeDimensionaSystemFigure = Figure(size = (800, 500))
TDSAxis = Axis3(ThreeDimensionaSystemFigure[:,1]; azimuth = 2.2, elevation = 0.35  , title = "System Attractors", xlabel = "x", ylabel = "y", zlabel = "z")
colours = ["blue", "red", "green", "pink"]
labels = ["$(FractionsBoA[1])", "$(FractionsBoA[2])", "$(FractionsBoA[3])"]#[L"0.0025", L"0.3325", L"0.3325"]
for i in 1:3# change to a length at some point
    x, y, z  = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2], GroupCentresTrajectories[i][3]
    scatter!(TDSAxis, x[end], y[end], z[end], label = labels[i], color = colours[i], marker = :star6, markersize = 20)
end
axislegend("Basin of Attraction Fractions", position = :rt)
ThreeDimensionaSystemFigure

#######################################################################################################################

### Test for 2.3 ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl")
include("CW_7_fixing.jl")

  ## Test ##

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length = 100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
push!(ds, ODE)
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)
#GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)
println(GroupIndex)

# want a plot of the ics colour-coded based on their index
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:,:], limits = (-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
i = 1
for InitialConditions in ics #i in 1:length(ics)
    v, w = InitialConditions
    if GroupIndex[i] == 1 
        scatter!(BasinAxis, v, w, color = "LightBlue", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, v, w, color = "PeachPuff", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, v, w, color = "PaleGreen", marker = :rect, markersize = 10)
    end
    i += 1
end
BasinsFigure

###############################################################################################################################################

using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_8_fixing.jl")


## Test ##

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length = 100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
push!(ds, ODE)
end

GroupIndex, GroupCentresTrajectories = fg2(ds, featurizer, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)

println(GroupIndex)
println(GroupCentresTrajectories)
length(GroupCentresTrajectories[3][1])

# Double check if it wants full group centre trajectory or just attractor trajectories
colours = ["blue", "red", "green"]
for i in 1:3 # change to a length at some point
    x, y = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2]
    lines!(BasinAxis, x, y, color = colours[i])
end
scatter!(BasinAxis, 0, 0, color = "green", marker = :star6, markersize = 20) # since the trajectory to the group centre is not clear
BasinsFigure

########################################################################################################################################

using DynamicalSystems, Statistics, LinearAlgebra
include("CW_4.jl"), include("CW_6.jl"), include("CW_7.jl")

# Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

l1 = length(v)
l2 = length(w)
K = 0

ics = Vector{Vector{Float64}}(undef, l1 * l2)

for i in 1:l1
    for j in 1:l2
        K += 1
        ics[K] = [v[i], w[j]]
    end
end

indices = fg1(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

#length(indices)
#println(indices)

##################################################################################################################################

### This one works - YAY! ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl")#, include("CW_7.jl") 
include("fg3_almostfinal.jl")

#Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p 
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0]
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

l1 = length(v)
l2 = length(w)
K = 0

ics = Vector{Vector{Float64}}(undef, l1 * l2)

for i in 1:l1
    for j in 1:l2
        K += 1
        ics[K] = [v[i], w[j]]
    end
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)
#GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)
#println(GroupIndex)

# want a plot of the ics colour-coded based on their index
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:,:], limits = (-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
for i in 1:length(ics)
    v, w = ics[i] 
    if GroupIndex[i] == 1 
        scatter!(BasinAxis, v, w, color = "LightBlue", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, v, w, color = "PeachPuff", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, v, w, color = "PaleGreen", marker = :rect, markersize = 10)
    end
end
BasinsFigure

##########################################################################################################################################

using DynamicalSystems, LinearAlgebra, Statistics, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), #include("CW_9_old.jl")
include("fg3_almostfinal.jl")
## Test fg3 ##

#Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ds = neuron_system!

# Grid of initial conditions
v = collect(-2.5:0.05:2.5)
w = collect(-2.5:0.05:2.5)

l1 = length(v)
l2 = length(w)
K = 0

ics = Vector{Vector{Float64}}(undef, l1 * l2)

for i in 1:l1
    for j in 1:l2
        K += 1
        ics[K] = [v[i], w[j]]
    end
end

# Call the function
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

println(FractionsBoA) 
println(GroupIndex)

# next step - make this a dictionary which pairs the fractions with the group index

############################################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
#include("CW_10.jl"), include("Featurizer2.6.jl"), include("T2.6_for_CW_9.jl") # Once finalised, make this CW_9
include("CW_10.jl"), include("CW_11.jl"), include("fg3_almostfinal.jl") #include("T2.7_for_CW_9._2.jl")

# Function Rule
ds = DS26

# parameters and intial conditions
Parameters = 0.16
# Fixed set of initial conditions
ics = Vector{Vector{Float64}}()
# temp set for testing
#yg = xg = zg = collect(range(-3, 3; length=4))
#zg .+= 0.1
yg = xg = zg = range(-3, 3; length = 20)
zg .+= 0.1
ics = [[x, y, z] for x in xg for y in yg for z in zg]

yg = xg = zg = collect(range(-3, 3; length = 20))
zg .+= 0.1

for x in xg
    for y in yg
        for z in zg
            ics_temp = [x, y, z]
            push!(ics, ics_temp)
        end
    end
end

# Time Values
StartTime = 0.0
TransientTime = 0.0
TotalTime = 100.0
SamplingTime = 0.2

# threshold value
r = 0.2

GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer2, r, ics; Ttr=50.0, Dt=0.5, T=100.0)

println(FractionsBoA)

### this newer version of fg3 agrees with the previous one

# Plot the system attractors
ThreeDimensionaSystemFigure = Figure(size = (800, 500))
TDSAxis = Axis3(ThreeDimensionaSystemFigure[:,1]; azimuth = 2.2, elevation = 0.35  , title = "System Attractors", xlabel = "x", ylabel = "y", zlabel = "z")
colours = ["blue", "red", "green", "pink"]
labels = [L"0.0025", L"0.3325", L"0.3325", L"0.3325"]
for i in 1:4# change to a length at some point
    x, y, z  = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2], GroupCentresTrajectories[i][3]
    scatter!(TDSAxis, x[end], y[end], z[end], label = labels[i], color = colours[i], marker = :star6, markersize = 20)
end
axislegend("Basin of Attraction Fractions", position = :rt)
ThreeDimensionaSystemFigure

###############################################################################################################################################################

using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("fg3_almostfinal2_for2.7.jl")
#include("Featurizer2.7.jl")
include("CW_13.jl")
include("CW_12.jl")

# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 100)#1000)
ics = [[x, y] for x in xg for y in yg ]

# Initialise relevant functions, parameters and initial conditions
ds = DTDS
r = 1.5
Parameters = [2.9, 0.66]
# Time Values
StartTime = 0
TransientTime = 200
TotalTime = 50
SamplingTime = 1
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

println(FractionsBoA)
#println(GroupIndex)
# this means i should expect three groups

BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:,1], limits = (-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
for i in 1:length(ics)
    x, y = ics[i] 
    if GroupIndex[i] == 1 
        scatter!(BasinAxis, x, y, color = "LightBlue", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, x, y, color = "PeachPuff", marker = :rect, markersize = 10)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, x, y, color = "PaleGreen", marker = :rect, markersize = 10)
    #elseif GroupIndex[i] == 4
     #   scatter!(BasinAxis, x, y, color = "LightPink", marker = :rect, markersize = 10)
    end
end
elem_1 = [PolyElement(color = "LightBlue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "PeachPuff", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "PaleGreen", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,2], [elem_1, elem_2, elem_3], ["0.1286", "0.1286", "0.7428"], "Basin of Attraction Fraction",  patchsize = (35, 35), rowgap = 10, framevisible = false)
BasinsFigure 

###########################################################################################################################################################

using LinearAlgebra, CairoMakie

include("CW_5.jl")
# input arguments
randompoint() = (rand([-1, 1]) + 0.1randn(), 0.1randn())
features = [randompoint() for _ in 1:1000]
r = 0.8 # threshold distance

GroupIndex = simple_grouping(features, r) # outputs a vector telling us if the point belongs in group 1 or 2

println(GroupIndex)

SimpleGroupingFigure = Figure()
SimpleGroupingAxis = Axis(SimpleGroupingFigure[:, :]; title = "Simple Grouping Plot", xlabel = "x", ylabel = "y")
for i in 1:1000
    if GroupIndex[i] == 1
        scatter!(SimpleGroupingAxis, features[i][1], features[i][2], color="blue")
    else
        scatter!(SimpleGroupingAxis, features[i][1], features[i][2], color="red")
    end
end
SimpleGroupingFigure

#######################################################################################################################################
### THE OLD ONE ====
using LinearAlgebra, CairoMakie

include("CW_5old.jl")

# input arguments
randompoint() = (rand([-1, 1]) + 0.1randn(), 0.1randn())
features = [randompoint() for _ in 1:1000]
r = 0.8 # threshold distance

I = simple_grouping(features, r) # outputs a vector telling us if the point belongs in group 1 or 2

fig = Figure()
ax = Axis(fig[:, :];)
    # separate each tuple into x and y values
    x = Vector{Float64}(undef, 1000)
    y = Vector{Float64}(undef, 1000)
for i in 1:1000
    x[i], y[i] = features[i][1], features[i][2]
    if I[i] == 1
        scatter!(ax, x[i], y[i], color="blue")
    else
        scatter!(ax, x[i], y[i], color="red")
    end
end
fig

################################################################################################################################################

#=
The grouping algorithm takes the first point and makes it the group centre. It then looks at a points' coordinates and considers the distance between 
it and the group centre. If it is less than the threshold value, they are put in the same group (assigned same group index). 
Otherwise they are put in the second group. The plot shows each point colour-coded to their group.



look at a coordinate test its distance from previous coords, if its sufficiently close, then add it to the group


1. start by populating state space with initial conditions - as it has an x and y value just make a for loop over both - want to test from -3 to 3 in each direction
2. next make a vector containing each initial condition
3. find the trajectories for each initial condition and mark the end points
4. extract features 

I created an algorithm that provided a vector of vectors and a real number (the threshold), will return the vector of group indices that each feature belongs in. 
It does this by looking at the distance (Euclidean Norm) between each vector and if it is less than the threshold it goes into Group 1, 
otherwise it goes into Group 2.=#

########################################################################################################################################################################

### Test for 2.3 ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_7_fixing.jl")

## Test ##

r = 0.4
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 100.0
TotalTime = 50.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length=100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
    ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
    push!(ds, ODE)
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)
#GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)
#println(GroupIndex)

    GroupDict = Dict{Int16,Int16}()
    # count each index
    for i in 1:length(GroupIndex)
        if haskey(GroupDict, GroupIndex[i])
            GroupDict[GroupIndex[i]] += 1
        else
            GroupDict[GroupIndex[i]] = 1
        end
    end



##############################################
vg = wg = range(-2.5, 2.5; length=100)
ics = [[v, w] for v in vg for w in wg]

# Fake output for example
result = [sin(v) * cos(w) for (v, w) in ics]

# Reshape back to grid
Z = reshape(result, (length(wg), length(vg)))'

using Plots

heatmap(wg, vg, Z)
###############################################


# not currently working
vplot = [v for v in vg]  
wplot = [w for w in wg]
vg = wg = (range(-2.5, 2.5; length=100))
HeatmapPlot = Matrix{Float64}(undef, 100, 100)
k = 1
for i in 1:100
    for j in 1:100
        HeatmapPlot[i,j] = GroupIndex[k]
        k += 1
    end
end
heatmap(vplot, wplot, HeatmapPlot)
heatmap(vg, wg, HeatmapPlot)

Z = reshape(GroupIndex, (length(wg), length(vg)))'
heatmap(vg, wg, Z)





#=
# want a plot of the ics colour-coded based on their index
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:, 1:3], limits=(-2.5, 2.5, -2.5, 2.5), title="Basins of Attraction", xlabel="v", ylabel="w")
i = 1
for InitialConditions in ics #i in 1:length(ics)
    v, w = InitialConditions
    if GroupIndex[i] == 1
        scatter!(BasinAxis, v, w, color="LightBlue", marker=:rect, markersize=10)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, v, w, color="PeachPuff", marker=:rect, markersize=10)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, v, w, color="PaleGreen", marker=:rect, markersize=10)
    end
    i += 1
end
BasinsFigure=#

########################################################################################################################################################

### Test for 2.3 ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_7_fixing.jl")

## Test ##

r = 0.4
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 100.0
TotalTime = 400.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length=100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
    ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
    push!(ds, ODE)
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)

    GroupDict = Dict{Int16,Int16}()
    # count each index
    for i in 1:length(GroupIndex)
        if haskey(GroupDict, GroupIndex[i])
            GroupDict[GroupIndex[i]] += 1
        else
            GroupDict[GroupIndex[i]] = 1
        end
    end

    Z = reshape(GroupIndex, (length(wg), length(vg)))'
heatmap(vg, wg, Z, colormap = Reverse(:deep))
heatmap(vg, wg, Z, colormap = ["LightBlue", "PeachPuff", "PaleGreen"])

########################################################################################################################################################

### Test for 2.3 ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_7_fixing.jl")

## Test ##

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 100.0
TotalTime = 400.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length=100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
    ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
    push!(ds, ODE)
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)
#GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)
#println(GroupIndex)

# want a plot of the ics colour-coded based on their index
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:, 1:3], limits=(-2.5, 2.5, -2.5, 2.5), title="Basins of Attraction", xlabel="v", ylabel="w")
i = 1
for InitialConditions in ics #i in 1:length(ics)
    v, w = InitialConditions
    if GroupIndex[i] == 1
        scatter!(BasinAxis, v, w, color="LightBlue", marker=:rect, markersize=10)
    elseif GroupIndex[i] == 2
        scatter!(BasinAxis, v, w, color="PeachPuff", marker=:rect, markersize=10)
    elseif GroupIndex[i] == 3
        scatter!(BasinAxis, v, w, color="PaleGreen", marker=:rect, markersize=10)
    end
    i += 1
end
BasinsFigure

#########################################################################################################################

### Test for 2.3 ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_7_fixing.jl")

## Test ##

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 200.0
TotalTime = 200.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length=100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
    ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
    push!(ds, ODE)
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)
GroupIndexMatrix = reshape(GroupIndex, (length(wg), length(vg)))'

# Figure
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:, 1:3], limits=(-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
heatmap!(BasinAxis, vg, wg, GroupIndexMatrix, colormap = ["LightBlue", "PeachPuff", "PaleGreen"])
BasinsFigure

# only took 17mins to run - Yay!!

using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_8_fixing.jl")

## Test ##

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 200.0
TotalTime = 200.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length = 100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
push!(ds, ODE)
end

GroupIndex, GroupCentresTrajectories = fg2(ds, featurizer, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)

BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:, 1:3], limits=(-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
# Double check if it wants full group centre trajectory or just attractor trajectories
colours = ["blue", "red", "green"]
for i in 1:3 # change to a length at some point
    x, y = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2]
    lines!(BasinAxis, x, y, color = colours[i])
end
scatter!(BasinAxis, 0, 0, color = "green", marker = :star6, markersize = 20) # since the trajectory to the group centre is not clear
elem_1 = [PolyElement(color = "blue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "red", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "green", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,4], [elem_1, elem_2, elem_3], ["Outer limit cycle", "Inner limit cycle", "Stable node"],
 "Attractors",  patchsize = (35, 35), rowgap = 10, framevisible = false)
BasinsFigure

###### Test with squares shapes instead 
colours = [(:blue, 0.3), (:red, 0.3), (:green, 0.3)]
markers = [:circle, :rect, :star6]
for i in 1:3 # change to a length at some point
    x, y = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2]
    scatter!(BasinAxis, x, y, color = colours[i], marker = markers[i], markersize = 7)
end
#scatter!(BasinAxis, 0, 0, color = "green", marker = :star6, markersize = 20) # since the trajectory to the group centre is not clear
elem_1 = [PolyElement(color = "blue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "red", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "green", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,4], [elem_1, elem_2, elem_3], ["Outer limit cycle", "Inner limit cycle", "Stable node"],
 "Attractors",  patchsize = (35, 35), rowgap = 10, framevisible = false)
BasinsFigure

#############################################################################################################################################################

# This is the newer version
using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("CW_13.jl")
include("CW_12.jl")
include("fg3_FORALLOFTHEM.jl")

# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 100)#1000)
ics = [[x, y] for x in xg for y in yg ]

# Initialise relevant functions, parameters and initial conditions
r = 1.5
Parameters = [2.9, 0.66]
StartTime = 0
TransientTime = 200
TotalTime = 50
SamplingTime = 1

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
    ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
    push!(ds, ODE)
end

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)
GroupIndexMatrix = reshape(GroupIndex, (length(yg), length(xg)))'

# Figure - see if I can get this to work with a title, this seems more proper and clear
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:, 1:3], limits=(-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "x", ylabel = "y")
Heatmap = heatmap!(BasinAxis, xg, yg, GroupIndexMatrix, colormap = ["LightBlue", "PeachPuff", "PaleGreen"])
Colorbar(BasinsFigure[:, end+1], Heatmap; title = "Group Index and Basin of Attraction Fractions")
BasinsFigure
## ==IMPORTANT== ##
###############################################################################################################################################################

using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")

# Initialise relevant functions, parameters and initial conditions

values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# threshold value
r = 0.2

# Time Values
StartTime = 0.0
TransientTime = 200.0
TotalTime = 400.0
SamplingTime = 0.5



# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions = []
AllVMax = []
AllWMax = []
AllStdW = []
#prange = [-1.0, -0.9, -0.8]

for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end


using LaTeXStrings

# Figure
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# The number of groups changes with the parameters, this can be indicated by length(GroupCentresTrajectories). 
# From this I then need to plot each of the values in that row against prange
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j])
        scatter!(MaxVAxis, prange[i], AllVMax[i][j])
        scatter!(MaxWAxis, prange[i], AllWMax[i][j])
        scatter!(StdAxis, prange[i], AllStdW[i][j])
    end
end
# Latex label
#Label(GlobalContinuationFigure[0, :], text = L"\\textbf{Effects of varying the parameter p}", fontsize = 30)
# Rich label
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Normal label
#Label(GlobalContinuationFigure[0, :], text = "Effects of varying the parameter p", fontsize = 30)
# Show Figure
GlobalContinuationFigure

# I now need to group the same attractors by colour
# the max value seems to be consistently in one group or another, this could be a useful way to group the data

## Mathcing ##
# So I have several vectors, lets focus first on the basins of attraction vector. Depending on the number of attractors, the elements in each row changes 
# (technincally, it's a vector of vectors, so more accurate to say size of vectors at each index). 
# The coursework sheet recommends I use the absolute difference of the standard deviation of the w coordinate of each attractor


####
#
####

ColourIndex = Vector{Vector{Int16}}() # gives an index to group each feature by for plotting
# needs to be a vector containing vectors of integers
# set first indinces
a_old, b_old = [], []
# b needs an index for each stdw in first row
p = 1
for stdw in AllStdW[1]
    push!(a_old, stdw) # intial one is to compare
    push!(b_old, p) # each time the index increases add it as an element
    p += 1
end
push!(ColourIndex, b_old) # add first indices as first vector in colour index
# so a_old now contains first stdw values for each attractor (in this case there's only one at this point)
for a in AllStdW[2:end]
    #println(a)
    # as we can see from the print, a just goes sequentially through each row.
    # next we need to start asigning each one to a group 
    #Differences = []
    ColourIndexTemp = [] # this is basically b
    b = []
    loopcompleted = 1
    for i in 1:length(a) 
        #println(i)
        for j in 1:length(a_old)
            println(j)
            println(b)
            println(b_old)
            DifferenceTemp = abs(a[i] - a_old[j]) # compare a values to a_old values, if any difference is small (below some threshold), then match it
            if DifferenceTemp < 0.01
                # match a[i] and a_old[j]
                b_temp = b_old[j] # store the relevant index 
                push!(b, b_temp)  # add it to the list (as we are testing in order of the a[i], the position in the list should match)
                #ColourIndexTemp = 
            else
                # if an a value can't be matched to the old group, then create a new one
                groups = length(a_old) #  check how many groups there were
                newgroup = groups + 1 # add one
                push!(b, newgroup) # add the new group index
            end

            #push!(Differences, DifferencesTemp)
        end
    end
    push!(ColourIndex, b)  # now add the indices for the next values
a_old = a # store current a to compare back to
println(loopcompleted)
loopcompleted += 1
end

####

############################################################################################################################################################################################################################################################

# --- Matching (replacement for your ColourIndex building code) ---
MatchingThreshold = 0.1   # threshold for matching std(w) to an existing group

ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise reps and counts from the first row of AllStdW
reps = copy(AllStdW[1])                     # representative std(w) for each group
counts = fill(1, length(reps))              # how many members assigned to each rep
FirstRowIndices = collect(1:length(reps)) # [1,2,3,...]
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 

# Process subsequent rows
for row in AllStdW[2:end]
    row_assignment = Int[]  # will hold group indices for this row
    for val in row
        # compute distances to current representatives
        dists = abs.(reps .- val)
        MinimumDistance, MinimumIndex = findmin(dists)

        if MinimumDistance < MatchingThreshold
            # assign to the closest existing group and update its representative (running mean)
            push!(row_assignment, MinimumIndex)
            counts[MinimumIndex] += 1
            reps[MinimumIndex] = (reps[MinimumIndex]*(counts[MinimumIndex]-1) + val) / counts[MinimumIndex]
        else
            # create a new group
            NewIndex = length(reps) + 1
            push!(reps, val)
            push!(counts, 1)
            push!(row_assignment, NewIndex)
        end
    end
    push!(ColourIndex, row_assignment)
end

########################################

# --- Stepwise matching based on AllStdW ---
match_thresh = 0.05  # threshold for std(w) similarity

ColourIndex = Vector{Vector{Int}}()

# ---- initialise first row ----
prev_std = AllStdW[1]
prev_groups = collect(1:length(prev_std))  # [1, 2, 3, ...]
push!(ColourIndex, prev_groups)
next_group_id = length(prev_std) + 1       # next unused group index

# ---- iterate over remaining rows ----
for i in 2:length(AllStdW)
    current_std = AllStdW[i]
    current_groups = Int[]

    for val in current_std
        # find closest attractor in previous step
        dists = abs.(prev_std .- val)
        min_dist, idx = findmin(dists)

        if min_dist < match_thresh
            # assign same group as closest attractor in previous step
            push!(current_groups, prev_groups[idx])
        else
            # create new group
            push!(current_groups, next_group_id)
            next_group_id += 1
        end
    end

    push!(ColourIndex, current_groups)
    # update "previous" state for next iteration
    prev_std = current_std
    prev_groups = current_groups
end


##############################################################################################################################################################

# --- Stepwise one-to-one matching for attractor groups ---
match_thresh = 0.05  # tolerance for std(w) similarity

ColourIndex = Vector{Vector{Int}}()

# Initialise with first set of attractors
prev_std = AllStdW[1]
prev_groups = collect(1:length(prev_std))
push!(ColourIndex, prev_groups)
next_group_id = length(prev_std) + 1

# Step through parameter values
for i in 2:length(AllStdW)
    current_std = AllStdW[i]
    current_groups = Int[]

    used_prev = falses(length(prev_std))  # to ensure one-to-one matches

    for val in current_std
        # Compute distances to previous attractors
        dists = abs.(prev_std .- val)
        min_dist, idx = findmin(dists)

        if min_dist < match_thresh && !used_prev[idx]
            # Match to previous attractor
            push!(current_groups, prev_groups[idx])
            used_prev[idx] = true
        else
            # Start a new group
            push!(current_groups, next_group_id)
            next_group_id += 1
        end
    end

    push!(ColourIndex, current_groups)
    prev_std = current_std
    prev_groups = current_groups
end

##############################################################################################################################################################################################################

# --- One-to-one stepwise matching ---
match_thresh = 0.4

ColourIndex = Vector{Vector{Int}}()

prev_std = AllStdW[1]
prev_groups = collect(1:length(prev_std))
push!(ColourIndex, prev_groups)
next_group_id = length(prev_std) + 1

for i in 2:length(AllStdW)
    current_std = AllStdW[i]
    current_groups = Int[]

    used_prev = falses(length(prev_std))

    for val in current_std
        dists = abs.(prev_std .- val)
        min_dist, idx = findmin(dists)

        if min_dist < match_thresh && !used_prev[idx]
            push!(current_groups, prev_groups[idx])
            used_prev[idx] = true
        else
            push!(current_groups, next_group_id)
            next_group_id += 1
        end
    end

    push!(ColourIndex, current_groups)
    prev_std = current_std
    prev_groups = current_groups
end

GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")

using Colors

n_groups = maximum(vcat(ColourIndex...))
palette = distinguishable_colors(n_groups)

# Each attractor group gets a consistent color across all parameter values
for (i_p, pval) in enumerate(prange)
    for (j, group_id) in enumerate(ColourIndex[i_p])
        scatter!(BoAFractionsAxis, pval, AllBoAFractions[i_p][j], color = palette[group_id])
        scatter!(MaxVAxis, pval, AllVMax[i_p][j], color = palette[group_id])
        scatter!(MaxWAxis, pval, AllWMax[i_p][j], color = palette[group_id])
        scatter!(StdAxis, pval, AllStdW[i_p][j], color = palette[group_id])
    end
end

Label(GlobalContinuationFigure[0, :],
      text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true),
      fontsize = 30)

      GlobalContinuationFigure


      #######################################################################################################################

      print(GroupCentresTrajectories)

v = Vector{Vector{Any}}()
w = Vector{Vector{Any}}()

for i in 1:3 # change to a length at some point
push!(v, GroupCentresTrajectories[i][1])
push!(w, GroupCentresTrajectories[i][2])
end

MaxV = Vector{Vector{Any}}()
MaxW = Vector{Vector{Any}}()
StdW = Vector{Vector{Any}}()
for i = 1:3
MaxV = maximum(v[i])
push!(MaxV, maximum(v[i]))
MaxW = maximum(w[i])
push!(MaxW, maximum(w[i]))
StdW = std(w[i])
push!(StdW, Maximum)
for



##################################################################################################################################################


# trying to get a stacked line plot
import Pkg
Pkg.add("PlotlyJS")

using PlotlyJS
x = collect(0:10)
y1 = x
y2 = ones(size(x))
t1 = scatter(x=x, y=y1, stackgroup="one")
t2 = scatter(x=x, y=y2, stackgroup="one")
plot([t1, t2])

#############

# Sample data
x = 1:12  # Represents months
y1 = [10, 12, 15, 20, 18, 25, 30, 28, 35, 40, 38, 42]  # Category 1
y2 = [5, 8, 10, 12, 14, 15, 18, 20, 23, 25, 28, 30]   # Category 2
y3 = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]     # Category 3

# Stack the data
y_stack = [y1 y2 y3]
y_cumulative = cumsum(y_stack, dims=2)

# Plotting
fig = Figure()
ax = Axis(fig[1, 1], title="Stacked Area Plot")

for i in 1:size(y_stack, 2)
    if i == 1
        fill_between!(ax, x, 0, y_cumulative[:, i], label="y$i")
    else
        fill_between!(ax, x, y_cumulative[:, i-1], y_cumulative[:, i], label="y$i")
    end
end

fig

# not a great solution - seems online user use chatgpt plus it isn't makie

#######################

using CairoMakie
f = Figure()
Axis(f[1, 1])

xs = 1:0.2:10
ys_low = -0.2 .* sin.(xs) .- 0.25
ys_high = 0.2 .* sin.(xs) .+ 0.25

band!(xs, ys_low, ys_high)
band!(xs, ys_low .- 1, ys_high .-1, color = :red)

f

#########################

function stacked_area_plot!(ax, x, ys)
    y_stack = reduce(hcat, ys)
    y_cumulative = cumsum(y_stack, dims=2)

    for i in 1:size(y_stack, 2)
        if i == 1
            band!(ax, x, zeros(length(x)), y_cumulative[:, i], label="Category $i")
        else
            band!(ax, x, y_cumulative[:, i-1], y_cumulative[:, i], label="Category $i")
        end
    end
    
    nothing
end
 
f = Figure()
ax = Axis(f[1, 1])

x = 1:12  
y1 = [10, 12, 15, 20, 18, 25, 30, 28, 35, 40, 38, 42] 
y2 = [5, 8, 10, 12, 14, 15, 18, 20, 23, 25, 28, 30]   
y3 = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25] 

ys = [y1, y2, y3]

stacked_area_plot!(ax, x, ys)
f

######
using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")
# Initialise relevant functions, parameters and initial conditions
values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]
# threshold value
r = 0.2
# Time Values
StartTime = 0.0
TransientTime = 200.0
TotalTime = 400.0
SamplingTime = 0.5
# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions = []
AllVMax = []
AllWMax = []
AllStdW = []
# 
for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end

println(AllBoAFractions)

#######

# y1, y2, y3, y4 up to yn, represent y values for each 
 
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value

    end
end
########################################################################################################

using CairoMakie
f = Figure()
Axis(f[1, 1])

xs = 1:0.2:10
ys_low = -0.2 .* sin.(xs) .- 0.5
ys_high = 0.2 .* sin.(xs) .+ 0.5

band!(xs, ys_low, ys_high)
band!(xs, ys_low .- 1, ys_high .-1, color = :red)

f

############################################################################################################################################


using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")

# Initialise relevant functions, parameters and initial conditions

# Parameter curve
θrange = range(0, 2π; length = 100)

values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
#ics = [[v, w] for w in wg for v in vg]
ics = [[v, w + I] for w in wg for v in vg]
# threshold value
r = 0.2

# Time Values
StartTime = 0.0
TransientTime = 200.0
TotalTime = 400.0
SamplingTime = 0.5

# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions = []
AllVMax = []
AllWMax = []
AllStdW = []
#prange = [-1.0, -0.9, -0.8]

for θ in θrange # change for prange once working
    p = -1.025 + 0.005 * cos(θ)
    I = 7.5 * sin(θ) 
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, I, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end


## Matching Step ## 
MatchingThreshold = 0.2   # threshold for matching std(w) to an existing group

ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonMeans = copy(AllStdW[1])                     # representative std(w) for each group (a current mean to compare to)
counts = fill(1, length(ComparisonMeans))              # how many members assigned to each rep
FirstRowIndices = collect(1:length(ComparisonMeans)) # [1,2,3,...] (index for each group)
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 

# Process subsequent rows
for row in AllStdW[2:end] 
    row_assignment = Int[]  # will hold group indices for current row
    for CurrentValue in row
        # compute distances to current representatives
        Difference = abs.(ComparisonMeans .- CurrentValue) # compare current row values to previous ones
        MinimumDifference, MinimumIndex = findmin(Difference)

        if MinimumDifference < MatchingThreshold
            # assign to the closest existing group and update its representative (running mean)
            push!(row_assignment, MinimumIndex)
            counts[MinimumIndex] += 1
            # Compare Sequenetially
            #PreviousValue[MinimumIndex] = val # value to compare to on next stwp 
            #ComparisonMeans[MinimumIndex] = CurrentValue # value to compare to on next stwp 
            # Compare to a mean
            ComparisonMeans[MinimumIndex] = (ComparisonMeans[MinimumIndex]*(counts[MinimumIndex]-1) + CurrentValue) / counts[MinimumIndex] # adjust means for new values
        else
            # create a new group
            NewIndex = length(ComparisonMeans) + 1
            push!(ComparisonMeans, CurrentValue)
            push!(counts, 1)
            push!(row_assignment, NewIndex)
        end
    end
    push!(ColourIndex, row_assignment)
end


# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against θ", xlabel = "θ", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against θ", xlabel = "θ", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against θ", xlabel = "θ", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against θ", xlabel = "θ", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(θrange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, θrange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxVAxis, θrange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxWAxis, θrange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(StdAxis, θrange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Labels
# Latex label
#Label(GlobalContinuationFigure[0, :], text = L"\\textbf{Effects of varying the parameter p}", fontsize = 30)
# Rich label
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter θ"; fontweight = "bold", underline = true), fontsize = 30)
# Normal label
#Label(GlobalContinuationFigure[0, :], text = "Effects of varying the parameter p", fontsize = 30)
# Legend
# Build legend entries
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))] #length(Colours)]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))] # loop over number of colours used #length(Colours)]

# Add the legend to the figure
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")
# Show Figure
GlobalContinuationFigure

#####
#=
y = Array{Float64, 2}(undef, length(AllBoAFractions), maximum(maximum(ColourIndex))) # 2D array rows for each parameter, width for all attractors
# each column in y represents a different attractor as indicated by the ColourIndex, this will either be 0 or a fraction of the BoA
for i = 1:length(θrange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
       y[i, j] = AllBoAFractions[i][j] 
    end
end=#
# doesn't currently work since it always 'fills up the y' instead of putting in approriate space, so make it depend on Colour Index to work
y = zeros(length(AllBoAFractions), maximum(maximum(ColourIndex)))
#y = Array{Float64, 2}(undef, length(AllBoAFractions), maximum(maximum(ColourIndex))) # 2D array rows for each parameter, width for all attractors
# each column in y represents a different attractor as indicated by the ColourIndex, this will either be 0 or a fraction of the BoA
for i = 1:length(θrange) # loop over however many p values are tested
    k = 1
    for j in ColourIndex[i] # = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        #println(j)
       y[i, j] = AllBoAFractions[i][k] 
       k += 1 # since AllBoAFractions doesn't have zeros for non-existent attractors
    end
end

# check there are no abnormally large values
findall(x -> x > 1, y)

f = Figure()
ax = Axis(f[1, 1])
ys = Vector{Vector{Float64}}()
for i in 1:size(y, 2) # number of columns in y
    push!(ys, y[:, i])
end
# now ys contains a column vector for each attractor
x = collect(θrange)
stacked_area_plot!(ax, x, ys)
f

#####################################################

# Initial Conditions
r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 200.0
TotalTime = 200.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length=1000))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
    ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
    push!(ds, ODE)
end

# Example
# define a state space grid to compute the basins on:
xg = yg = range(-2, 2; length = 201)
# find attractors using recurrences in state space:
mapper = AttractorsViaRecurrences(henon, (xg, yg); sparse = false)
# compute the full basins of attraction:
basins, attractors = basins_of_attraction(mapper; show_progress = false)

####################
# My Attempt
include("CW_4.jl")
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
InitialConditions = [1.0, 1.0]
neuron = CoupledODEs(neuron_system!, InitialConditions, Parameters)
vg = wg = range(-2.5, 2.5; length = 1000)
mapper = AttractorsViaRecurrences(neuron, (vg, wg); sparse = false,
    Δt = 0.01,                             # integration step
    Ttr = 400,                            # transient time
)
basins, attractors = basins_of_attraction(mapper; show_progress = false)
# plot
heatmap_basins_attractors((vg, wg), basins, attractors)
# 2
using CairoMakie
plot_attractors(attractors)

# using second mapper
include("CW_4.jl")
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
InitialConditions = [1.0, 1.0]
neuron = CoupledODEs(neuron_system!, InitialConditions, Parameters)
vg = wg = range(-2.5, 2.5; length = 100)
# featurizer 
mapper = AttractorsViaFeaturizing(neuron, (vg, wg); sparse = false,
    Δt = 0.01,                             # integration step
    Ttr = 400,                            # transient time
    )
basins, attractors = basins_of_attraction(mapper; show_progress = false)

###################################################################
values = 100 #21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ds = neuron
grid = (
    range(-2.5, 2.5; length = values), # x
    range(-2.5, 2.5; length = values), # y
)
sampler, = statespace_sampler(grid)
prange = -1.1:0.01:0.05

params(θ) = [1 => 5 + 0.5cos(θ), 2 => 0.1 + 0.01sin(θ)]
θs = range(0, 2π; length = 101)
pcurve = params.(θs)

#prange = 4.5:0.01:6
pidx = 9 # index of the parameter
mapper = AttractorsViaRecurrences(
    neuron,                                # your CoupledODEs system
    (range(-2.5, 2.5; length = values),    # x range
     range(-2.5, 2.5; length = values));   # y range
    Δt = 0.01,                             # integration step
    Ttr = 2000,                            # transient time
)
ascm = AttractorSeedContinueMatch(mapper)
fractions_cont, attractors_cont = global_continuation(
    ascm, prange, pidx, sampler;
    samples_per_parameter = 1_000,
    show_progress = true
)
animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, prange, pidx;
);


fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange,
)

fig = animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, prange, pidx;
    figure = (size = (600, 700),),
    axis = (ylabel = "y", xlabel = "x"),
    savename = "attracont_extra.mp4",
    add_legend = false,
    a2rs = [A -> minimum(A[:, 1]), A -> maximum(A[:, 2])],
    a2rs_ylabels = ["x-min", "y-max"],
    a2rs_ratio = 0.33,
    vline_kwargs = (linestyle = :dash, linewidth = 3, color = "red"),
)





######

include("CW_12.jl")
Parameters = [2.9, 0.66]
InitialConditions = [1.0, 1.0]
ds = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
xg = yg = range(-2.5, 2.5; length = 100)
mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse = false,
    Δt = 1,                             # integration step
    Ttr = 10000,                            # transient time
)
basins, attractors = basins_of_attraction(mapper; show_progress = false)
# plot
heatmap_basins_attractors((xg, yg), basins, attractors)


####
# this is useful to note
x, y = columns(X)
summary.((x, y))

##########################################

using DynamicalSystems, LinearAlgebra, Statistics, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), #include("CW_9_old.jl")
#include("fg3_FORALLOFTHEM.jl")
include("fg3_adjustedbasins.jl")
## Test fg3 ##

#Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, 0.01] # a1, a2, a3, a4, c, ϵ, d, I, p
#For example, we can probe an elipsoid defined as



StartTime = 0.0
TransientTime = 200.0
TotalTime = 200.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length = 21))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
for InitialConditions in ics
ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
push!(ds, ODE)
end

# Call the function
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)
GroupIndexMatrix = reshape(GroupIndex, (length(wg), length(vg)))'

# Figure
BasinsFigure = Figure(size = (1000, 800))
BasinAxis = Axis(BasinsFigure[:, 1:3], limits=(-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
heatmap!(BasinAxis, vg, wg, GroupIndexMatrix, colormap = ["LightBlue", "PeachPuff", "PaleGreen"])

# Double check if it wants full group centre trajectory or just attractor trajectories
colours = ["blue", "red", "green"]
markers = [:circle, :rect, :star6]
for i in 1:maximum(maximum(GroupIndex)) # change to a length at some point
    x, y = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2]
    #lines!(BasinAxis, x, y, color = colours[i])
    scatter!(BasinAxis, x, y, color = colours[i], marker = markers[i], markersize = 7)
end
#scatter!(BasinAxis, 0, 0, color = "green", marker = :star6, markersize = 20) # since the trajectory to the group centre is not clear
elem_1 = [PolyElement(color = "blue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "red", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "green", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,4], [elem_1, elem_2, elem_3], ["Outer limit cycle", "Inner limit cycle", "Stable node"],
 "Attractors",  patchsize = (35, 35), rowgap = 10, framevisible = false)
BasinsFigure

# note to self:
x, y = columns(Trajectory)

#################################################################################

values = 100 #21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ds = neuron
grid = (
    range(-2.5, 2.5; length = values), # x
    range(-2.5, 2.5; length = values), # y
)
sampler, = statespace_sampler(grid)
prange = -1.1:0.01:0.05

mapper = AttractorsViaRecurrences(
    neuron,                                # your CoupledODEs system
    (range(-2.5, 2.5; length = values),    # x range
     range(-2.5, 2.5; length = values));   # y range
    Δt = 0.01,                             # integration step
    Ttr = 4000,                            # transient time
    sparse = true,                  # slightly faster & less memory
    force_non_adaptive = true
)

params(θ) = [9 => -1.025 + 0.005 * cos(θ), 8 => 7.5 * sin(θ)]
θs = range(0, 2π; length = 100)
pcurve = params.(θs)

matcher = MatchBySSSetDistance(use_vanished = true)

ascm = AttractorSeedContinueMatch(mapper, matcher)

fractions_cont, attractors_cont = global_continuation(
	ascm, pcurve, sampler; samples_per_parameter = 1_000
)

animate_attractors_continuation(
    neuron, attractors_cont, fractions_cont, pcurve;
    savename = "curvecont.mp4"
);

#prange = 4.5:0.01:6
pidx = 8 & 9 # index of the parameters, come back to
mapper = AttractorsViaRecurrences(
    neuron,                                # your CoupledODEs system
    (range(-2.5, 2.5; length = values),    # x range
     range(-2.5, 2.5; length = values));   # y range
    Δt = 0.01,                             # integration step
    Ttr = 2000,                            # transient time
)
ascm = AttractorSeedContinueMatch(mapper)
fractions_cont, attractors_cont = global_continuation(
    ascm, prange, pidx, sampler;
    samples_per_parameter = 1_000,
    show_progress = true
)
animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, prange, pidx;
);


fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange,
)

fig = animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, prange, pidx;
    figure = (size = (600, 700),),
    axis = (ylabel = "y", xlabel = "x"),
    savename = "attracont_extra.mp4",
    add_legend = false,
    a2rs = [A -> minimum(A[:, 1]), A -> maximum(A[:, 2])],
    a2rs_ylabels = ["x-min", "y-max"],
    a2rs_ratio = 0.33,
    vline_kwargs = (linestyle = :dash, linewidth = 3, color = "red"),
)

#################################################
values = 100
grid = (
    range(-10.0, 10.0; length = values),
    range(-10.0, 10.0; length = values)
)

mapper = AttractorsViaRecurrences(
    neuron,
    grid;
    Δt = 0.01,
    Ttr = 3000
)

# Define 2D parameter curve (p8, p9)
params(θ) = [9 => -1.025 + 0.005 * cos(θ), 8 => 7.5 * sin(θ)]
θs = range(0, 2π; length = 100)
pcurve = params.(θs)

# Match attractors between steps
matcher = MatchBySSSetDistance(use_vanished = true)
ascm = AttractorSeedContinueMatch(mapper, matcher)

# Run continuation
fractions_cont, attractors_cont = global_continuation(
    ascm, pcurve, sampler;
    samples_per_parameter = 1_000
)

# Animate continuation
animate_attractors_continuation(
    neuron, attractors_cont, fractions_cont, pcurve;
    savename = "curvecont.mp4"
)


########################################################################################################################




using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_4.jl"), include("CW_6.jl"), include("fg3_FORALLOFTHEM.jl")

# Initialise relevant functions, parameters and initial conditions

values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# threshold value
r = 0.2

# Time Values
StartTime = 0.0
TransientTime = 200.0
TotalTime = 400.0
SamplingTime = 0.5



# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions = []
AllVMax = []
AllWMax = []
AllStdW = []
#prange = [-1.0, -0.9, -0.8]

for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end

# Matching
MatchingThreshold = 0.2  # threshold for matching std(w) to an existing group  
#AllValues = [] # test
# Initial Step
PreviousValues = AllStdW[1]
ColourIndex = [collect(1:length(PreviousValues))] # index for each initial value
#push!(AllValues, PreviousValues) # test
# Subsequent Steps
for CurrentRow in AllStdW[2:end]
    RowAssignment = Int[]
    for val in CurrentRow
        Differences = abs.(PreviousValues .- val) # compare values in each step
        MinimumDifference, Index = findmin(Differences) # add the minimum
        1
        println(Differences)
        println(MinimumDifference)
        if MinimumDifference < MatchingThreshold
            push!(RowAssignment, Index)
        else
            push!(RowAssignment, length(PreviousValues) + 1)
            #push!(PreviousValues, val)  # remember the new group for next comparison
        end
    end
    push!(ColourIndex, RowAssignment)
    # Replace previous row — we forget earlier history
    PreviousValues = CurrentRow
   # push!(AllValues, PreviousValues) # test
end

# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Title Label (using rich)
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Legend - build legend entries
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))] 
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))] # loop over number of colours used #length(Colours)]
# Add the legend to the figure
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")

# Show Figure
GlobalContinuationFigure


#########################################################################

using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_4.jl"), include("CW_6.jl"), include("CW_8_fixing.jl"), include("fg3_FORALLOFTHEM.jl")

## Test ##

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.01] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 200.0
TotalTime = 200.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length = 100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
push!(ds, ODE)
end

GroupIndex, GroupCentresTrajectories, BoAFractionsAxis = fg3(ds, featurizer, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)

BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[1,1])
# Double check if it wants full group centre trajectory or just attractor trajectories
colours = ["blue", "red", "green"]
markers = [:circle, :rect, :star6]
for i in 1:maximum(maximum(GroupIndex)) # change to a length at some point
    x, y = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2]
    #lines!(BasinAxis, x, y, color = colours[i])
    scatter!(BasinAxis, x, y, color = colours[i], marker = markers[i], markersize = 7)
end
#scatter!(BasinAxis, 0, 0, color = "green", marker = :star6, markersize = 20) # since the trajectory to the group centre is not clear
elem_1 = [PolyElement(color = "blue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "red", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "green", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,4], [elem_1, elem_2, elem_3], ["Outer limit cycle", "Inner limit cycle", "Stable node"],
 "Attractors",  patchsize = (35, 35), rowgap = 10, framevisible = false)
BasinsFigure


    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

println(std(w[3]))
# double check AllStdW is correct by manually calculating StdW
# at -1.2 there is one attractor
# stdw = 2.08354e-7
# at -1.03 attractor at origin 
# 3.7977952988092773e-7
# at -1.02 limit cycle and attractor - bifurcation
# 0.9050703394707829, 3.634608862795079e-7
# at -1.01 second limit cycle plus above - bifurcation
# 1.362577749153058, 0.9413708580021818, 3.493561378947487e-7


#############################################################################################
n_groups = maximum(maximum(ColourIndex))
palette = distinguishable_colors(n_groups)  # from ColorSchemes.jl or Colors.jl
Markers = repeat([:rect, :circle, :diamond, :cross, :star5, :star4], ceil(Int, n_groups/6))[1:n_groups]
Colours = [palette[i] for i in 1:n_groups]


for group_id in 1:n_groups
    xs = Float64[]
    ys_stdw = Float64[]
    for (k, p) in enumerate(prange)
        idxs = findall(==(group_id), ColourIndex[k])
        if !isempty(idxs)
            push!(xs, p)
            push!(ys_stdw, AllStdW[k][idxs[1]])
        end
    end
    lines!(StdAxis, xs, ys_stdw; color = Colours[group_id], linewidth = 1.5)
end 

using CairoMakie, Colors

GlobalContinuationFigure = Figure(size = (1000,1000))
axes = (
    boa = Axis(GlobalContinuationFigure[1,1], title="BoA fractions against p", xlabel="p", ylabel="BoA Fractions"),
    maxv = Axis(GlobalContinuationFigure[1,2], title="max(v) against p", xlabel="p", ylabel="Max v"),
    maxw = Axis(GlobalContinuationFigure[2,1], title="max(w) against p", xlabel="p", ylabel="Max w"),
    stdw = Axis(GlobalContinuationFigure[2,2], title="std(w) against p", xlabel="p", ylabel="Standard deviation of w")
)

# Dynamic colours and markers
n_groups = maximum(maximum(ColourIndex))
palette = distinguishable_colors(n_groups)
markers = repeat([:rect, :circle, :diamond, :cross, :star5, :star4], ceil(Int, n_groups/6))[1:n_groups]

for (i, p) in enumerate(prange)
    n_attr = length(AllStdW[i])
    for j in 1:n_attr
        c = ColourIndex[i][j]
        scatter!(axes.boa, p, AllBoAFractions[i][j], color=palette[c], marker=markers[c])
        scatter!(axes.maxv, p, AllVMax[i][j], color=palette[c], marker=markers[c])
        scatter!(axes.maxw, p, AllWMax[i][j], color=palette[c], marker=markers[c])
        scatter!(axes.stdw, p, AllStdW[i][j], color=palette[c], marker=markers[c])
    end
end

Label(GlobalContinuationFigure[0, :], text=rich("Effects of varying the parameter p"; fontweight="bold", underline=true), fontsize=30)

# Legend
legend_labels = ["Attractor $i" for i in 1:n_groups]
legend_elems = [MarkerElement(color=palette[i], marker=markers[i], markersize=15) for i in 1:n_groups]
Legend(GlobalContinuationFigure[3,1:2], legend_elems, legend_labels; orientation=:horizontal, framevisible=false, title="Attractor Groups")

GlobalContinuationFigure

################################################################################################################################################

# Return Times

# Paramters and initial conditions
p0 = [10.0, 8 / 3, 28.0]
u0 = [10.0, 2.0, 20.0] # initial condition
total_time = 1000.0 # total integration time
sampling_time = 0.02 # step time

# Set index and value
index = 3 # x,y,z = 1,2,3
value = 27.0

# Dynamical system
lorenz63 = CoupledODEs(lorenz63_rule!, u0, p0)
X, t = trajectory(lorenz63, total_time; Δt=sampling_time)
x, y, z = columns(X)

# Poincaré Surface of the Section
x_ps, y_ps, z_ps, t_ps = psos2(X, t, index, value) # index 

ReturnTimes = []
# return times vector
for i = 2:length(t_ps)
    push!(ReturnTimes, t_ps[i] - t_ps[i-1])
end

println(ReturnTimes) # return times vectors


##################################################################################################################################################################


### Test for 2.3 ###
using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_neuron.jl"), include("CW_featurizer.jl"), include("CW_7_fixing.jl")

## Test ##

r = 0.2
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
StartTime = 0.0
TransientTime = 200.0
TotalTime = 200.0
SamplingTime = 0.5

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length=100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
    ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
    push!(ds, ODE)
end

GroupIndex = fg1(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)
GroupIndexMatrix = reshape(GroupIndex, (length(wg), length(vg)))'

# Figure
BasinsFigure = Figure(size = (1000, 800))
BasinAxis = Axis(BasinsFigure[:, 1:3], limits=(-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "v", ylabel = "w")
heatmap!(BasinAxis, vg, wg, GroupIndexMatrix, colormap = ["LightBlue", "PeachPuff", "PaleGreen"])
BasinsFigure


## Test ##
include("CW_8_fixing.jl")

# Initial Conditions
vg = wg = (range(-2.5, 2.5; length = 100))
ics = [[v, w] for v in vg for w in wg]

# Dynamical System
ds = Vector{Any}()
i = 1
for InitialConditions in ics
ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
push!(ds, ODE)
end

GroupIndex, GroupCentresTrajectories = fg2(ds, featurizer, r, ics; Ttr = TransientTime, Dt = SamplingTime, T = TotalTime)

# Double check if it wants full group centre trajectory or just attractor trajectories
colours = ["blue", "red", "green"]
markers = [:circle, :rect, :star6]
for i in 1:3 # change to a length at some point
    x, y = GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2]
    #lines!(BasinAxis, x, y, color = colours[i])
    scatter!(BasinAxis, x, y, color = colours[i], marker = markers[i], markersize = 10)
end
#scatter!(BasinAxis, 0, 0, color = "green", marker = :star6, markersize = 20) # since the trajectory to the group centre is not clear
elem_1 = [PolyElement(color = "blue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "red", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "green", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,4], [elem_1, elem_2, elem_3], ["Outer limit cycle", "Inner limit cycle", "Stable node"],
 "Attractors",  patchsize = (35, 35), rowgap = 10, framevisible = false)
BasinsFigure

#####################################################################################################################

using DynamicalSystems, Statistics, LinearAlgebra, CairoMakie
include("CW_DTDS.jl"), include("CW_featurizer3.jl")
# Basins of Attraction and Basin fractions

# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 1000)
ics = [[x, y] for x in xg for y in yg ]

# Set key values
r = 1.5
Parameters = [2.9, 0.66]
StartTime, TransientTime, TotalTime, SamplingTime = 0, 200, 500, 1

# Dynamical System
ds = Vector{Any}()
Features = []
for InitialConditions in ics
ODE = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
push!(ds, ODE)
end

for DS in ds
traj, t = trajectory(DS, TotalTime; Δt = SamplingTime)
features_temp = featurizer3(traj)
push!(Features, features_temp)
end

Features

fig = Figure()
ax = Axis3(fig[1,1], limits = (-5,5,-5,5,-5,5))
ax.azimuth = π/4
xplot, yplot, zplot, wplot = [], [], [], []
for i in 1:length(Features)
xtemp ,ytemp, ztemp, wtemp = Features[i][1], Features[i][2], Features[i][3], Features[i][4]
push!(xplot, xtemp)
push!(yplot, ytemp)
push!(zplot, ztemp)
push!(wplot, wtemp)
end
scatter!(ax, zplot, wplot, yplot)
fig

fig = Figure()
ax = Axis(fig[1,1], limits = (-5,5,-5,5))
xplot, yplot, zplot, wplot = [], [], [], []
for i in 1:length(Features)
xtemp ,ytemp, ztemp, wtemp = Features[i][1], Features[i][2], Features[i][3], Features[i][4]
push!(xplot, xtemp)
push!(yplot, ytemp)
push!(zplot, ztemp)
push!(wplot, wtemp)
end
scatter!(ax, xplot, yplot)#zeros(length(xplot)))
fig

#####################################################################################################################

using DynamicalSystems, LinearAlgebra, Statistics, CairoMakie
include("CW_featurizer3.jl"), include("CW_fg3.jl"), include("CW_DTDS.jl")
# Basins of Attraction and Basin fractions

# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 1000)
ics = [[x, y] for x in xg for y in yg ]

# Set key values
r = 1.5
Parameters = [2.9, 0.66]
StartTime, TransientTime, TotalTime, SamplingTime = 0, 200, 500, 1

# Dynamical System
ds = Vector{Any}()
for InitialConditions in ics
ODE = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
push!(ds, ODE)
end

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)
GroupIndexMatrix = reshape(GroupIndex, (length(yg), length(xg)))'

# How many groups are there and what are their fractions
println(maximum(maximum(GroupIndex)))
println(FractionsBoA) # change fractions calc to workout divergent group also

# Figure
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:, :], limits=(-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "x", ylabel = "y")
colours1 = ["Crimson", "MidnightBlue", "LightYellow"]
colours2, markers = ["Green", "Yellow"], [:rect, :star5]
heatmap!(BasinAxis, xg, yg, GroupIndexMatrix, colormap = colours1)
for i in 1;length(GroupCentresTrajectories)
    scatter!(BasinAxis, GroupCentresTrajectories[i][1], GroupCentresTrajectories[i][2], color = colours2[i], marker = markers[i], markersize = 10)
end
elem_1 = [PolyElement(color = "LightBlue", strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = "PeachPuff", strokecolor = :black, strokewidth = 1)]
elem_3 = [PolyElement(color = "PaleGreen", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,end+1], [elem_1, elem_2, elem_3], ["0.1286", "0.1286", "0.7428"], "Basin of Attraction Fraction",  patchsize = (35, 35), rowgap = 10, framevisible = false)

BasinsFigure # show figure

#######################################################################################################################################################################################


include("CW_featurizer3.jl")
# Basins of Attraction and Basin fractions

# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 1000)
ics = [[x, y] for x in xg for y in yg ]

# Set key values
r = 1.5
Parameters = [2.9, 0.66]
StartTime, TransientTime, TotalTime, SamplingTime = 0, 200, 500, 1

# Dynamical System
ds = Vector{Any}()
for InitialConditions in ics
ODE = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
push!(ds, ODE)
end

# Featurize and Group
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)
GroupIndexMatrix = reshape(GroupIndex, (length(yg), length(xg)))'

# How many groups are there and what are their fractions
println(maximum(maximum(GroupIndex)))
println(FractionsBoA) # change fractions calc to workout divergent group also

# Attractor positions
println(GroupCentresTrajectories[1][1][end], "|",  GroupCentresTrajectories[1][2][end], "|", GroupCentresTrajectories[2][1][end], "|", GroupCentresTrajectories[2][2][end])

# Figure
BasinsFigure = Figure()
BasinAxis = Axis(BasinsFigure[:, :], limits=(-2.5, 2.5, -2.5, 2.5), title = "Basins of Attraction", xlabel = "x", ylabel = "y")
heatmap!(BasinAxis, xg, yg, GroupIndexMatrix, colormap = ["Crimson", "MidnightBlue", "LightYellow"])
colours, markers = ["Green", "Yellow"], [:rect, :star5]
for i in 1:length(GroupCentresTrajectories)
scatter!(BasinAxis, GroupCentresTrajectories[i][1][end], GroupCentresTrajectories[i][2][end], color = colours[i], marker = markers[i], markersize = 20)
end
elem_1 = [PolyElement(color = "Crimson", strokecolor = :black, strokewidth = 1), MarkerElement(color = "Green", marker = :rect, markersize = 15)]
elem_2 = [PolyElement(color = "MidnightBlue", strokecolor = :black, strokewidth = 1), MarkerElement(color = "Yellow", marker = :star5, markersize = 15)]
elem_3 = [PolyElement(color = "LightYellow", strokecolor = :black, strokewidth = 1)]
Legend(BasinsFigure[:,end+1], [elem_1, elem_2, elem_3], ["0.129052", "0.129052", "0.741896"], "Basin of Attraction Fraction",  patchsize = (35, 35), rowgap = 10, framevisible = false)

BasinsFigure # show figure

##############################################################################################################################################################################################################

using DynamicalSystems, LinearAlgebra, Statistics, CairoMakie
include("CW_fg3.jl"), include("CW_neuron.jl"), include("CW_featurizer.jl")

# Test - Global continuation

# Initialise relevant functions, parameters and initial conditions
values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# threshold value
r = 0.2
# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0.0, 200.0, 400.0, 0.5
# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions, AllVMax, AllWMax, AllStdW = [], [], [], []

for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging

    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end

    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)

    v = Vector{Vector{Any}}()
    w = Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1])
        push!(w, GroupCentresTrajectories[i][2])
    end

    # Calculate some values to compare as parameters change
    MaxV = Vector{Float64}() # Max v for each attractor
    MaxW = Vector{Float64}() # Max w for each attractor
    StdW = Vector{Float64}() # standard deviation of w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i]))
        push!(MaxW, maximum(w[i]))
        push!(StdW, std(w[i]))
    end

    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA) 
    push!(AllVMax, MaxV)
    push!(AllWMax, MaxW)
    push!(AllStdW, StdW)
end


# --- Matching 
MatchingThreshold = 0.15  # threshold for matching std(w) to an existing group  
# somewhere between 0.6-0.7 it links the broken attractor together

ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonValues = copy(AllStdW[1])                     # representative std(w) for each group
counts = fill(1, length(ComparisonValues))              # how many members assigned to each comparison vector
FirstRowIndices = collect(1:length(ComparisonValues)) # [1,2,3,...] (index for each group)
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 

# Process subsequent rows
for row in AllStdW[2:end] 
    row_assignment = Int[]  # will hold group indices for current row
    for CurrentValue in row
        # compute distances to current representatives
        Difference = abs.(ComparisonValues .- CurrentValue) # compare current row values to previous ones
        MinimumDifference, MinimumIndex = findmin(Difference)

        if MinimumDifference < MatchingThreshold
            # assign to the closest existing group and update its representative (running mean)
            push!(row_assignment, MinimumIndex)
            counts[MinimumIndex] += 1
            #PreviousValue
            ComparisonValues[MinimumIndex] = CurrentValue # value to compare to on next stwp 
            #ComparisonMeans[MinimumIndex] = (ComparisonMeans[MinimumIndex]*(counts[MinimumIndex]-1) + CurrentValue) / counts[MinimumIndex] # mean value of all values
        else
            # create a new group
            NewIndex = length(ComparisonValues) + 1
            push!(ComparisonValues, CurrentValue)
            push!(counts, 1)
            push!(row_assignment, NewIndex)
        end
    end
    push!(ColourIndex, row_assignment)
end
###################################################
ColourIndex = Vector{Vector{Int}}()  # indices for each row (p value)

# Initialise ComparisonMeans and counts from the first row of AllStdW
ComparisonValues = copy(AllStdW[1])                     # representative std(w) for each group
counts = fill(0, 100)              # how many members assigned to each comparison vector plus redundant values beyond that in case there are more groups (not most elegant solution)
for j in 1:length(ComparisonValues)
counts[j] = 1
end
FirstRowIndices = collect(1:length(ComparisonValues)) # [1,2,3,...] (index for each group)
push!(ColourIndex, FirstRowIndices)       # first row: one group per rep 
i = 1 # index the row number
println("row complete $i")

for row in AllStdW
    i += 1 # new row number
    row_assignment = zeros(Int, length(row)) # vector of zeros the same length as the row
    k = 1 # column/ value in the row are we on
    for CurrentValue in row
        Difference = abs.(ComparisonValues .- CurrentValue) # difference between all values in current row and current group representative
        MinimumDifference, MinimumIndex = findmin(Difference) # which is closet (most likely to be in the same group as the rep)
        # On the 9th iteration, the 2nd index matches the 1st index on the 8th, so the ColourIndex of the first index in colour index is 1 so I need to match it to 1
        MatchedIndex = ColourIndex[i-1][] # check the colour index corresponding to min index at this step -  next I need to tell it the 2nd value in the row should be a 1 in ColourIndex

        if MinimumDifference < MatchingThreshold 
            row_assignment[MinimumIndex] = MatchedIndex # set relevant index in the row to matched colour index
            #push!(row_assignment, MatchedIndex) # which index in comparison values does it match
            counts[MatchedIndex] += 1 

        else
            # create a new group
            NewIndex = maximum(maximum(ColourIndex)) + 1 # new index one higher than current highest
           # push!(ComparisonValues, CurrentValue)
            push!(counts, 1) # start the count at 1
            row_assignment[MinimumIndex] = NewIndex 
            #push!(row_assignment, NewIndex) 
        end
            k += 1
    end
        ComparisonValues = row # the current row will be compared to on the subsequent step
    #ColourIndex[i][]
    push!(ColourIndex, row_assignment)
    println("row complete $i")
end
##################################################################################################
ComparisonValues = copy(AllStdW[1])
counts = fill(0, 100)
for j in 1:length(ComparisonValues)
    counts[j] = 1
end

FirstRowIndices = collect(1:length(ComparisonValues))
ColourIndex = [FirstRowIndices]
println("row complete 1")

for i in 2:length(AllStdW)
    row = AllStdW[i]
    row_assignment = zeros(Int, length(row))
    k = 1

    for CurrentValue in row
        Difference = abs.(ComparisonValues .- CurrentValue)
        MinimumDifference, MinimumIndex = findmin(Difference)
        MatchedIndex = k <= length(ColourIndex[i-1]) ? ColourIndex[i-1][k] : maximum(ColourIndex[i-1]) + 1

        if MinimumDifference < MatchingThreshold
            row_assignment[MinimumIndex] = MatchedIndex
            counts[MatchedIndex] += 1
        else
            NewIndex = maximum(maximum(ColourIndex)) + 1
            push!(counts, 1)
            row_assignment[MinimumIndex] = NewIndex
        end

        k += 1
    end

    push!(ColourIndex, row_assignment)
    ComparisonValues = row
    println("row complete $i")
end
##################################################################################################################
# --- Setup ---
ComparisonValues = copy(AllStdW[1])  # representative std(w) for each group
counts = fill(1, length(ComparisonValues))  # start each with 1
ColourIndex = [collect(1:length(ComparisonValues))]  # first row: unique colours

println("row complete 1")

# --- Loop over remaining rows ---
for i in 2:length(AllStdW)
    row = AllStdW[i]
    row_assignment = zeros(Int, length(row))  # group colours for this row

    for k in 1:length(row)
        CurrentValue = row[k]

        # Step 1: find closest comparison value
        Difference = abs.(ComparisonValues .- CurrentValue)
        MinimumDifference, MinimumIndex = findmin(Difference)

        # Step 2: check threshold
        if MinimumDifference < MatchingThreshold
            # Step 3: match the colour of the corresponding comparison value
            # Use the previous row's colour index at that comparison position
            if MinimumIndex <= length(ColourIndex[i-1])
                MatchedIndex = ColourIndex[i-1][MinimumIndex]
            else
                # Safety fallback if lengths differ
                MatchedIndex = maximum(vcat(ColourIndex...)) + 1
            end

            row_assignment[k] = MatchedIndex
            counts[MatchedIndex] += 1

        else
            # Step 4: create a new group
            NewIndex = maximum(vcat(ColourIndex...)) + 1
            row_assignment[k] = NewIndex
            push!(counts, 1)
        end
    end

    # Step 5: update for next iteration
    push!(ColourIndex, row_assignment)
    ComparisonValues = row
    println("row complete $i")
end



# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Labels
# Latex label
#Label(GlobalContinuationFigure[0, :], text = L"\\textbf{Effects of varying the parameter p}", fontsize = 30)
# Rich label
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Normal label
#Label(GlobalContinuationFigure[0, :], text = "Effects of varying the parameter p", fontsize = 30)
# Legend
# Build legend entries
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))] #length(Colours)]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))] # loop over number of colours used #length(Colours)]

# Add the legend to the figure
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")
# Show Figure
GlobalContinuationFigure

########################################################################################################################################################

#  Finding the three attractors for the neuron system

# Initial state and parameters
InitialStates = [0.0, 0.0], [1.0, 1.0], [2.0, 0.0] # v, w
points = [L"(0.75,0.75)", L"(1.0,1.0)", L"(2.0,2.0)"]
Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, 0.07] # a1, a2, a3, a4, c, ϵ, d, I, p
TotalTime = 100.0
SamplingTime = 0.01

TrajectoriesFigure = Figure(size = (700, 500))
TrajectoriesAxis = Axis(TrajectoriesFigure[:, :]; title="Trajectories for neuron system", xlabel = "v", ylabel = "w")

for (u0, label) in zip(InitialStates, points)
    # ODE trajectories
    neuronODE = CoupledODEs(neuron_system!, u0, Parameters)
    N, t = trajectory(neuronODE, TotalTime; Δt = SamplingTime)
    v, w = columns(N)
    lines!(TrajectoriesAxis, v, w; label = label)
end
axislegend(TrajectoriesAxis, "Initial state", position = :rt)

TrajectoriesFigure # show figure

# between 0.06 and 0.07 it looks like the attractor at the centre vanishes and we are left with just the outer limit cycle


# input arguments
randompoint() = (rand([-1, 1]) + 0.1randn(), 0.1randn())
features = [randompoint() for _ in 1:1000]
r = 0.8 # threshold distance

gc1 = features[1]
comp_vec = features[2]
d = norm(gc1 .- comp_vec)

########################################################################################################################



using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_neuron.jl"), include("CW_featurizer.jl"), include("CW_fg3.jl")

# Initialise relevant functions, parameters and initial conditions
values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# threshold value
r = 0.2

# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0.0, 200.0, 400.0, 0.5
# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions, AllVMax, AllWMax, AllStdW = [], [], [], []


for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging
    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end
    # Group and featurise for that parameter
    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)
    # Extract components
    v, w = Vector{Vector{Any}}(), Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1]), push!(w, GroupCentresTrajectories[i][2])
    end
    # Calculate some features to compare as parameters change
    MaxV, MaxW, StdW = Vector{Float64}(), Vector{Float64}(), Vector{Float64}() # Max v, Max w and Std w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i])), push!(MaxW, maximum(w[i])), push!(StdW, std(w[i]))
    end
    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA), push!(AllVMax, MaxV), push!(AllWMax, MaxW), push!(AllStdW, StdW)
end

## Matching ## 
MatchingThreshold = 0.15  # threshold for matching std(w) to an existing group  

# Initial row
ComparisonValues = copy(AllStdW[1])  # representative std(w) for each group
counts = fill(1, length(ComparisonValues))  # start each with 1
ColourIndex = [collect(1:length(ComparisonValues))]  # first row: unique colours

# Loop over the rest of the rows
for i in 2:length(AllStdW)
    row = AllStdW[i]
    row_assignment = zeros(Int, length(row))  # group colours for this row

    for k in 1:length(row)
        CurrentValue = row[k]

        # Find closest comparison value
        Difference = abs.(ComparisonValues .- CurrentValue)
        MinimumDifference, MinimumIndex = findmin(Difference)

        # Check threshold
        if MinimumDifference < MatchingThreshold
            # Match the colour of the corresponding comparison value - use the previous row's colour index at that comparison position
            if MinimumIndex <= length(ColourIndex[i-1])
                MatchedIndex = ColourIndex[i-1][MinimumIndex]
            else
                # Safety fallback if lengths differ -- double check this later
                MatchedIndex = maximum(vcat(ColourIndex...)) + 1
            end

            row_assignment[k] = MatchedIndex
            counts[MatchedIndex] += 1

        else
            # Create a new group
            NewIndex = maximum(vcat(ColourIndex...)) + 1
            row_assignment[k] = NewIndex
            push!(counts, 1)
        end
    end

    # Update for next iteration
    push!(ColourIndex, row_assignment)
    ComparisonValues = row
end

# Figure 
BasinsFigure = Figure()
BasinsAxis = Axis(BasinsFIgure[1,1]; title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        lines!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]

# finish making it work as a line or stacked line area plot later
#####################################################################################################################################


using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems
include("CW_neuron.jl"), include("CW_featurizer.jl"), include("CW_fg3.jl"), include("CW_simple_matching.jl")


# Initialise relevant functions, parameters and initial conditions
values = 8 
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]

# Threshold value
r = 0.2

# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0.0, 200.0, 400.0, 0.5
# test different values of p
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
# these will contain values to compare to the p parameter
AllBoAFractions, AllVMax, AllWMax, AllStdW = [], [], [], []


for p in prange # change for prange once working
    Parameters = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, p] # a1, a2, a3, a4, c, ϵ, d, I, p 
    #println(p) # this helps indicate current loop whilst debugging
    # Dynamical System
    ds = Vector{Any}()
    i = 1
    for InitialConditions in ics
        ODE = CoupledODEs(neuron_system!, InitialConditions, Parameters)
        push!(ds, ODE)
    end
    # Group and featurise for that parameter
    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)
    # Extract components
    v, w = Vector{Vector{Any}}(), Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres
        push!(v, GroupCentresTrajectories[i][1]), push!(w, GroupCentresTrajectories[i][2])
    end
    # Calculate some features to compare as parameters change
    MaxV, MaxW, StdW = Vector{Float64}(), Vector{Float64}(), Vector{Float64}() # Max v, Max w and Std w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MaxV, maximum(v[i])), push!(MaxW, maximum(w[i])), push!(StdW, std(w[i]))
    end
    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA), push!(AllVMax, MaxV), push!(AllWMax, MaxW), push!(AllStdW, StdW)
end

## Matching ## 
MatchingThreshold = 0.15  # threshold for matching std(w) to an existing group  

# Initial row
ComparisonValues = copy(AllStdW[1])  # representative std(w) for each group
counts = fill(1, length(ComparisonValues))  # start each with 1
ColourIndex = [collect(1:length(ComparisonValues))]  # first row: unique colours

# Loop over the rest of the rows
for i in 2:length(AllStdW)
    row = AllStdW[i]
    row_assignment = zeros(Int, length(row))  # group colours for this row

    for k in 1:length(row)
        CurrentValue = row[k]

        # Find closest comparison value
        Difference = abs.(ComparisonValues .- CurrentValue)
        MinimumDifference, MinimumIndex = findmin(Difference)

        # Check threshold
        if MinimumDifference < MatchingThreshold
            # Match the colour of the corresponding comparison value - use the previous row's colour index at that comparison position
            if MinimumIndex <= length(ColourIndex[i-1])
                MatchedIndex = ColourIndex[i-1][MinimumIndex]
            else
                # Safety fallback if lengths differ
                MatchedIndex = maximum(vcat(ColourIndex...)) + 1
            end

            row_assignment[k] = MatchedIndex
            counts[MatchedIndex] += 1

        else
            # Create a new group
            NewIndex = maximum(vcat(ColourIndex...)) + 1
            row_assignment[k] = NewIndex
            push!(counts, 1)
        end
    end

    # Update for next iteration
    push!(ColourIndex, row_assignment)
    ComparisonValues = row
end


ColourIndex = simple_matching(MatchingThreshold, AllStdW)

# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; limits = (-1.1, 0.05, -0.05, 1.05), title = "BoA fractions against p", xlabel = "p", ylabel = "BoAFractions")
MaxVAxis = Axis(GlobalContinuationFigure[1,2]; title = "max(v) against p", xlabel = "p", ylabel = "Max v")
MaxWAxis = Axis(GlobalContinuationFigure[2,1]; title = "max(w) against p", xlabel = "p", ylabel = "Max w")
StdAxis = Axis(GlobalContinuationFigure[2,2]; title = "std(w) against p", xlabel = "p", ylabel = "Standard Deviation of w")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
        scatter!(BoAFractionsAxis, prange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxVAxis, prange[i], AllVMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MaxWAxis, prange[i], AllWMax[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(StdAxis, prange[i], AllStdW[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Rich label
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter p"; fontweight = "bold", underline = true), fontsize = 30)
# Legend
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))]
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")
# Show Figure
GlobalContinuationFigure

#######################################################################################################################################################

## Quick check for roots for task 3.3

import Pkg
Pkg.add("Roots")
##
using Roots

f(v) = 2.7778*v^3 - 2.3333*v^5 + 0.7619*v^7 - 0.0847*v^9 - (1 + 100v + 5v^2)

v_star = find_zero(f, 0.0)   # starting guess near 0
println(v_star)
v_star = find_zero(f, (-1, 1))  # bisection/secant hybrid

####
using Pkg
Pkg.add("NLsolve")
using NLsolve

function system!(F, x)
    v, w = x
    F[1] = 2.7778*v^3 - 2.3333*v^5 + 0.7619*v^7 - 0.0847*v^9 - w
    F[2] = 0.01*(1 + 100v + 5v^2 - w)
end

sol = nlsolve(system!, [0.0, 0.0])
println(sol.zero)


################################################################################################################################################################################################################
using LinearAlgebra, Statistics, CairoMakie, DynamicalSystems
include("CW_DTDS.jl") # DTDS function rule
include("CW_featurizer3.jl") # featurizer for task 2.7
include("CW_simple_matching.jl") # simple matching for task 3.1
include("CW_vanish_matching.jl") # vanish matching for taske 3.5
include("CW_period.jl") # find period of a given trajectory  
include("CW_fg3.jl") # featurise and group function 3

## Global continuation with a new system

# parameter range
μrange = vcat(2.8:0.01:3.6)
# Initial Conditions
xg = yg = range(-2.5, 2.5; length = 100)
ics = [[x, y] for x in xg for y in yg ]
# threshold value
r = 1.5
# Time Values
StartTime, TransientTime, TotalTime, SamplingTime = 0, 200, 500, 1
# Initialise empty vectors to store features for different parameter values
AllBoAFractions, AllMeanX, AllMeanY, AllPeriods = [], [], [], []
# Find features at different parameter values
for μ in μrange
    Parameters = [μ, 0.66]
    # Dynamical System
    ds = Vector{Any}()
    for InitialConditions in ics
    ODE = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
    push!(ds, ODE)
    end

    # Featurize and Group
    GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr = SamplingTime, Dt = SamplingTime, T = TotalTime)
    # test
    #println(maximum(GroupIndex)) 
    #println(length(GroupCentresTrajectories)) # if there is always a divergent value, then this value should be one less than the above
    # last group index is for he divergent group - assuming this always exists, I can tell it to ignore this in subsequent calcs and set it to some temp value like 0 or just omit it
    # Extract components
    x, y = Vector{Vector{Any}}(), Vector{Vector{Any}}()
    for i in 1:length(GroupCentresTrajectories) # loop for the number of group centres (divergent group not included in this)
        push!(x, GroupCentresTrajectories[i][1]), push!(y, GroupCentresTrajectories[i][2])
    end
    # Calculate some features to compare as parameters change
    MeanX, MeanY, Period = Vector{Float64}(), Vector{Float64}(), Vector{Float64}() # Max v, Max w and Std w for each attractor
    for i = 1:length(GroupCentresTrajectories)
        push!(MeanX, mean(x[i])), push!(MeanY, maximum(y[i])), push!(Period, period(GroupCentresTrajectories[i]))
    end
    # Add a zero at end of each feature vector for the divergent value (otherwise the grouping gets messed up later)
    #println(MeanX)
    push!(MeanX, 0), push!(MeanY, 0), push!(Period, 0)
    #println(MeanX) # these prints just check it's adding zeros correctly
    # Now make vectors which contain the values above for each p in prange
    push!(AllBoAFractions, FractionsBoA), push!(AllMeanX, MeanX), push!(AllMeanY, MeanY), push!(AllPeriods, Period)
end

## Matching Step ## 
MatchingThreshold = 0.2   # threshold for matching std(w) to an existing group
ColourIndex = simple_matching(MatchingThreshold, AllMeanX) # I believe I'm using the top one
ColourIndex = vanish_matching(MatchingThreshold, AllMeanX)

#########################

# Plot global continuation
GlobalContinuationFigure = Figure(size = (1000,1000))
# Axes
BoAFractionsAxis = Axis(GlobalContinuationFigure[1,1]; limits = (2.75, 3.65, -0.05, 1.05), title = "BoA fractions against μ", xlabel = "μ", ylabel = "BoAFractions")
MeanXAxis = Axis(GlobalContinuationFigure[2,1]; title = "mean(x) against μ", xlabel = "μ", ylabel = "Mean x")
MeanYxis = Axis(GlobalContinuationFigure[2,2]; title = "mean(y) against μ", xlabel = "μ", ylabel = "Mean y")
PeriodAxis = Axis(GlobalContinuationFigure[1,2]; title = "Period against μ", xlabel = "μ", ylabel = "Period of attractor trajectory")
# Markers and Colours for each attractor
Markers = [:rect, :circle, :diamond, :cross, :star5, :star4]
Colours = ["blue", "red", "green", "purple", "orange", "yellow" ]       
for i = 1:length(μrange) # loop over however many μ values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that μ value (0 take away one due to divergent group)
        scatter!(BoAFractionsAxis, μrange[i], AllBoAFractions[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MeanXAxis, μrange[i], AllMeanX[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(MeanYxis, μrange[i], AllMeanY[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
        scatter!(PeriodAxis, μrange[i], AllPeriods[i][j], color = Colours[ColourIndex[i][j]], marker = Markers[ColourIndex[i][j]], markersize = 7)
    end
end
# Title
Label(GlobalContinuationFigure[0, :], text = rich("Effects of varying the parameter μ"; fontweight = "bold", underline = true), fontsize = 30)
# Legend
GroupLabels = ["Attractor $i" for i in 1:maximum(maximum(ColourIndex))]
Attractors = [MarkerElement(color = Colours[i], marker = Markers[i], markersize = 15) for i in 1:maximum(maximum(ColourIndex))]
Legend(GlobalContinuationFigure[3, 1:2], Attractors, GroupLabels; orientation = :horizontal, framevisible = false, title = "Attractor Groups")
# Show Figure
GlobalContinuationFigure 

###################################################################################################################################################################

#=

The plots above show 4 bifurcations. The local bifurcations are when a bifurcation takes place at a fixed point in state space. Information about the flow around a bifurcation point (given by the Jacobian matrix) is sufficient to characterise the bifurcation. These are when
The global bifurcations are when extended invariant sets (eg. limit cycles) collide with other sets and the resulting change in topology of the trajectories are not restricted to a local neighbourhood, but have a global impact. These are when

I think they're all local, maybe the red limit cycle is global


There are 4 bifurcations in the continuation. The first local bifurcation is at $p = -1.02$, creating an additional attractor. Likely a transcritical bifurcation since an attractor has appeared.
 Then another local bifurcation at $p = -1.01$, creating anther additional attractor. There is what looks like a Hopf bifurcation (global bifurcation) at around . 
 There is a gradual weakening of the inner limit cycle (seen in the basin plot) and eventually there's just the outer attractor (which is a limit cycle) and the inner attractor suggesting a Hopf bifurcation. 
 A similar Hopf bifurcation happens to the attractor at the centre just around $p=-0.01$.
=#
