# Bains stack area plot

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

# y1, y2, y3, y4 up to yn, represent y values for each colour
# if ColourIndex == i then yi
y = Array{Float64, 2}(undef, length(AllBoAFractions), maximum(maximum(ColourIndex))) # 2D array rows for each parameter, width for all attractors
# each column in y represents a different attractor as indicated by the ColourIndex, this will either be 0 or a fraction of the BoA
for i = 1:length(prange) # loop over however many p values are tested
    for j = 1:length(AllBoAFractions[i]) # loop over however many groups there are for that p value
       y[i, j] = AllBoAFractions[i][j] 
    end
end

####

f = Figure()
ax = Axis(f[1, 1])
x = collect(prange)
y1, y2, y3, y4 = y[:,1], y[:,2], y[:,3], y[:,4]
ys = [y1, y2, y3, y4]
stacked_area_plot!(ax, x, ys)
f


f = Figure()
ax = Axis(f[1, 1])
ys = Vector{Vector{Float64}}()
for i in 1:size(y, 2) # number of columns in y
    push!(ys, y[:, i])
end
# now ys contains a column vector for each attractor
x = collect(prange)
stacked_area_plot!(ax, x, ys)
f

#=
StackedAreaPlot = Figure(size = (500, 500))
StackedAreaAxis = Axis(StackedAreaPlot[1,1]; limits = (-1.1, 0.1, 0.0, 1.0) title = "Basins of Attraction Fractions", xlabel = "p", ylabel = "fraction") # change to percentage

xPlot = prange
yPlot =#