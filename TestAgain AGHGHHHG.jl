using DynamicalSystems, CairoMakie, LinearAlgebra, Statistics
include("fg3_almostfinal2_for2.7.jl"), include("CW_12.jl"), include("CW_13.jl")
#include("Featurizer2.7.jl"), include("DiscreteTimeDynamicalSystem.jl"), include("T2.7_for_CW_9.jl")

# Initialise relevant functions, parameters and initial conditions
ds = DTDS
r = 1000
Parameters = [2.9, 0.66]
# Time Values
StartTime = 0
TransientTime = 200
TotalTime = 50
SamplingTime = 1
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


GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer3, r, ics; Ttr=TransientTime, Dt=SamplingTime, T=TotalTime)

#println(FractionsBoA)
println(GroupIndex)

### THIS IS DRIVING ME INSANE, NEED TO WORK OUT WHY IT KEEPS BREAKING AGHGHGHGHGHGH

#=
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
BasinsFigure =#