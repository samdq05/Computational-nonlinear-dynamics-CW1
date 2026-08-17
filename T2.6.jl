# Test for CW_10 function with CW_9 feature and grouping -- this script will be helpful for getting exactly 10000 points to plot
using DynamicalSystems, LinearAlgebra, Statistics, CairoMakie
include("CW_9.jl"), include("CW_10.jl"), include("Featurizer2.6.jl")

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

#=
for i in ics
    println(i)
end
length(ics)
=#

# I need to adjust fg3 so that it can take dynamical systems with different numbers of ODEs  
#fig = Figure()
#ax = Axis3(fig[:, :], title="Trajectory", xlabel="x", ylabel="y", zlabel="z")

#InitialConditions = [ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[8000], ics[2500], ics[3000]]#, ics[4], ics[5], ics[6]]

# Function Rule
ds = DS26

# Time Values
StartTime = 0.0
TransientTime = 0.0
TotalTime = 100.0
SamplingTime = 0.2

# Colours Vector for plot
#Colours = ["red", "blue", "green", "orange", "pink", "purple", "cyan", "aquamarine", "CadetBlue", "LightGray", "Tan", "Maroon"]
#i = 1

# Component Vectors
Trajectories = Vector{StateSpaceSet{3,Float64,SVector{3,Float64},Nothing}}()
Time = Vector
for ics in ics#InitialConditions
    RuleODE = CoupledODEs(ds, ics, Parameters)
    Trajectory, t = trajectory(RuleODE, TotalTime; Ttr=TransientTime, Δt=SamplingTime) # Discarding Ttr seconds of transient time
    #lines!(ax, Trajectory[:,1], Trajectory[:,2], Trajectory[:,3], color = Colours[i])
    #i += 1

    # Components for comparison
    push!(Trajectories, Trajectory) # add to trajectories vector

end
#fig


#typeof(Trajectories)

# looks like we have 2 different limit cycles,
# lets take a trajectory from each and plot there time series 

#######

# lets do a few time series to find which features will be useful 
#=
TimeSeriesFig = Figure()
TimeSeriesAxisX = Axis(TimeSeriesFig[1, 1]; title="time series of x", xlabel="time", ylabel="x")
TimeSeriesAxisY = Axis(TimeSeriesFig[2, 1]; title="time series of y", xlabel="time", ylabel="y")
TimeSeriesAxisZ = Axis(TimeSeriesFig[3, 1]; title="time series of z", xlabel="time", ylabel="z")

Time = collect(StartTime:SamplingTime:(TotalTime+TransientTime)) # set the time array for plotting

for traj in Trajectories[1:10]

    lines!(TimeSeriesAxisX, Time, traj[:, 1], color="red")
    lines!(TimeSeriesAxisY, Time, traj[:, 2], color="blue")
    lines!(TimeSeriesAxisZ, Time, traj[:, 3], color="green")
end

TimeSeriesFig
=#

# not particular useful with 64 values

# Plot the Feature space, let's try to see which features do a good job separating and grouping the data
#=
FeatureSpaceFigure = Figure()
FeatureSpaceAxis = Axis3(FeatureSpaceFigure[:, :], title="Feature Space", xlabel="std(x)", ylabel="std(y)", zlabel="std(z)")
for traj in Trajectories
    x, y, z = traj[:, 1], traj[:, 2], traj[:, 3]
    scatter!(FeatureSpaceAxis, std(x), std(y), std(z))
end
FeatureSpaceFigure
=#

# this shows all the features in state space, let's now make a featurizer for it
# this is almost perfect, but that special exception in the corner could ruin it

Features = Vector{Vector{Float64}}()
for traj in Trajectories
TrajectoryFeatures = featurizer2(traj)
push!(Features, TrajectoryFeatures) # adds the individual trajectory features to the overall features vector
#println(Features) #-- this tests it works
end

# now that I have a featurizer, let's group them.
##ComparisonVector = [Features[i]] # assuming I'm grouping the ith element of Features, it will look like this

# ok now lets call the modified fg3 function
include("T2.6_for_CW_9.jl")
# threshold value
r = 0.2

GroupIndex, GroupCentresTrajectories, FractionsBoA = fg3(ds, featurizer2, r, ics; Ttr=50.0, Dt=0.5, T=100.0)