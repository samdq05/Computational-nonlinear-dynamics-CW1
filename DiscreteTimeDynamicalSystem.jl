# Task 2.7 -- Discrete time dynamical system
function DTDS(u, p, n)
    x,y = u
    μ, j = p
    xn = y
    yn = μ * y - y^3 - j * x
    return SVector(xn, yn)
end

## Test - 1 ## 
using CairoMakie, LinearAlgebra, Statistics, DynamicalSystems

#InitialConditions = [1.0,1.0]
xg = yg = range(-2.5, 2.5; length = 100)#00)
ics = [[x, y] for x in xg for y in yg ]
Parameters = [2.9, 0.66]

#Trajectories = [], push!(Trajectories, Trajectory) # not right at the moment
#=
FigTest = Figure()
AxTest = Axis(FigTest[:,:], limits = (-5, 5, -5, 5), title = "end points with limits", xlabel = "x", ylabel = "y")
for InitialConditions in ics
#DynamicalSystemRule = DTDS(InitialConditions, Parameters, 0) -- unecessary step
DynamicalSystem = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)

TotalTime = 50#10_000
Trajectory, t = trajectory(DynamicalSystem, TotalTime)
x, y = Trajectory[:, 1], Trajectory[:, 2]
#Trajectory
scatter!(AxTest, x[end], y[end])
end
FigTest =#

include("fg3_almostfinal2_for2.7.jl")
include("Featurizer2.7.jl")

# Initialise relevant functions, parameters and initial conditions
ds = DTDS
r = 1.5#2.0
Parameters = [2.9, 0.66]
# Time Values
StartTime = 0
TransientTime = 200
TotalTime = 50
SamplingTime = 1
GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=SamplingTime, Dt=SamplingTime, T=TotalTime)

#println(FractionsBoA)
println(GroupIndex)
println(FractionsBoA)
# this means i should expect three groups
#=
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
BasinsFigure =#


## Test - 2 ##
#=
Parameters = [2.9, 0.66]
TotalTime = 50
ics = Vector{Vector{Float64}}()
xg = yg = collect(range(-2.5, 2.5; length = 10))#1000))
for x in xg
    for y in yg
            ics_temp = [x, y]
            push!(ics, ics_temp)
    end
end
Trajectories = Vector{StateSpaceSet{2, Float64, SVector{2, Float64}, Nothing}}()
for InitialConditions in ics
DynamicalSystem = DeterministicIteratedMap(DTDS, InitialConditions, Parameters)
Trajectory, t = trajectory(DynamicalSystem, TotalTime)
push!(Trajectories, Trajectory)
end

#println(Trajectories[1])

#TestFigure = Figure()
#TestAxis = Axis(TestFigure[:,:]; limits = (-100.0, 100.0, -100.0, 100.0), title = "2.7 Test Trajectory", xlabel = "x", ylabel = "y")
#scatter!(TestAxis, Trajectories[1][:,1], Trajectories[1][:,2], color = "blue")
#TestFigure

TestFigure2 = Figure()
TestAxis2 = Axis(TestFigure2[:,:]; limits = (-10.0, 10.0, -10.0, 10.0), title = "2.7 Test Trajectory", xlabel = "x", ylabel = "y")
for i = 1: length(Trajectories)
scatter!(TestAxis2, Trajectories[i][:,1][end], Trajectories[i][:,2][end], color = "blue")
end
TestFigure2

# From the plot above, it's clear that the features to consider should be the means 

# notice how many trajectories become enormous and tend to +-∞, these need there own special features that add them to a basin of divergence
=#
