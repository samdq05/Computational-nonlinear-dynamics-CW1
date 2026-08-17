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
features = [[0.0,0.0], [0.1,0.0], [1.0, 1.0], [1.0,1.1], [2.0,2.0], [2.0,2.1], [300.0,300.0]]
ics=[1,2,3,4,5,6, 7]

    ## Group ## 

    Group_index = Vector{Int16}(undef, length(ics))

    k = round(Int, 50 ÷ 0.5) # number of steps to discard due to transience 

    # First state 
    Group_centre_1 = [0.0, 0.0] # check if these can be removed
    distance_1 = 0.0
    p1 = 0
    Group_centre_2 = [0.0, 0.0]
    distance_2 = 0.0
    p2 = 0
    Group_centre_3 = [0.0, 0.0]
    distance_3 = 0.0
    p3 = 0

    Group_centre = [Group_centre_1, Group_centre_2, Group_centre_3]
    Distance = [distance_1, distance_2, distance_3]

    failed = 0
    debug1 = Vector{Float64}(undef, 7)
    debug2 = 0
    debug3 = 0
    debug4 = 0
    check_distance_1 = Vector{Float64}(undef, )

    J0 = Vector{Int16}(undef, 3) # vector of the indices of group centres

    for i = 1:length(ics)

        # Assign groups based on features
        std_v = features[i][1]
        std_w = features[i][2]
        Comparison_vector = [std_v, std_w]

        # First feature group #
        if p1 == 0 # if this group is empty, automatically fill it an set as group centre
            Group_centre[1] = Comparison_vector
        end

        # Check if it's in the first group
        Distance[1] = norm(Comparison_vector - Group_centre[1])
check_distance_1[i] = Distance[1]

        if Distance[1] < r # if it's below the threshold add it to that feature, otherwise make a new feature
            Group_index[i] = 1
            if p1 == 0
                #Group_centre[1] = [std_v, std_w]
                p1 += 1
                J1 = i # index of the group centre
                J0[1] = J1
            end

            # Second feature group #
            if p2 == 0 # if this group is empty, automatically fill it an set as group centre
                Group_centre[2] = Comparison_vector
            end
    

            Distance[2] = norm(Comparison_vector - Group_centre[2]) # If it isn't in the first group, check the second
            println("d2")
            println(Distance[2])

        elseif Distance[2] < r
            Group_index[i] = 2
            if p2 == 0
                #Group_centre[2] = [std_v, std_w]
                p2 += 1
                J2 = i # index of the group centre
                J0[2] = J2
            end

            if p2 >= 1
                debug1[i] = Distance[2]
                println(Comparison_vector)
            end

            # Third feature group #
            if p3 == 0 # if this group is empty, automatically fill it an set as group centre
                Group_centre[3] = Comparison_vector
            end

            Distance[3] = norm(Comparison_vector - Group_centre[3]) # If it isn't in the first or second group, check the third (it should be in this one at least, otherwise something has gone wrong - there should only be three groups)

            
        elseif Distance[3] < r 
            println("HERE")
            Group_index[i] = 3
            if p3 == 0
                #Group_centre[3] = [std_v, std_w]
                p3 += 1
                J3 = i # index of the group centre
                J0[3] = J3
            end

                        if p2 == 1
                debug2 = Distance[1]
                debug3 = Distance[2]
                debug4 = Distance[3]
            end

        else 
            failed += 1

        end
    end

    return Group_index, J0, failed, debug1, debug2, debug3, debug4, Distance, check_distance_1

end

# Test -- technically i want it 100 by 100, so come back and adjust later
r = 0.2
params = [2.7778, 2.3333, 0.7619, 0.0847, 100.0, 0.01, 5.0, 0.0, -1.0] # a1, a2, a3, a4, c, ϵ, d, I, p
ds = neuron_system!
#=
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
=#
ics=[1,2,3,4,5,6, 7]
indices, J0, failed, d1, d2, d3, d4, Dis, cd1 = fg1(ds, featurizer, r, ics; Ttr=50, Dt=0.5, T=100)

println(indices)
println(d1)
println(d2)
println(d3)
println(d4)
println(Dis)
println(cd1)
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

