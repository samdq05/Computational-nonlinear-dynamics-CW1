using Statistics, CairoMakie, LinearAlgebra

randompoint() = (rand([-1, 1]) + 0.1randn(), 0.1randn())
features = [randompoint() for _ in 1:1000]
#=
typeof(features)
f = features[1]
f[2]
=#
x = Vector{Float64}(undef, 1000)
y = Vector{Float64}(undef, 1000)
for i = 1:1000
x[i], y[i] = features[i][1], features[i][2]
end
x
typeof(x)

fig = Figure()
ax = Axis(fig[:,:]; title = "", xlabel = "", ylabel = "")
scatter!(ax, x, y)
fig

# from the scatter plot (also consider how points are generated) it's clear a threshold of around 0.5 is appropriate (approx.)

#features[1:1000][1] # not working as intended

maximum(x)
minimum(x)
maximum(x)-minimum(x)
Range = maximum(x)-minimum(x)

#=
Group_1_temp = Vector{Vector{Float64}}(undef, 1000)
Group_2_temp = Vector{Vector{Float64}}(undef, 1000)
k1, k2 = 1, 0
# Base case -- set the first group as defined by the intiial point
Group_1_temp[1] = [x[1], y[1]]
r = 0.4 # threshold distance
for i in 2:1000
    distance = norm([x[i], y[i]] - [x[i-1], y[i-1]])
if distance < r
    k1 += 1
    Group_1_temp[k1] = [x[i], y[i]]
else
    k2 += 1
    Group_2_temp[k2] = [x[i], y[i]]
end
end
# extract relevant info
Group_1 = Group_1_temp[1:k1]
Group_2 = Group_2_temp[1:k2]
=#
#################################
Group_1_temp = Matrix{Float64}(undef, 1000, 2)
Group_2_temp = Matrix{Float64}(undef, 1000, 2)
Group_index = Vector{Int16}(undef, 1000)# tells us what group a vector in the list belongs to
k1, k2 = 1, 0
# Base case -- set the first group as defined by the intiial point
Group_1_temp[1,:] = [x[1], y[1]]
Group_index[1] = 1
r = 0.5 # threshold distance
for i in 2:1000
    distance = norm([x[i], y[i]] - [x[1], y[1]])
if distance < r
    k1 += 1
    Group_1_temp[k1,:] = [x[i], y[i]]
    Group_index[i] = 1
else
    k2 += 1
    Group_2_temp[k2,:] = [x[i], y[i]]
    Group_index[i] = 2
end
end
# extract relevant info
Group_1 = Group_1_temp[1:k1,:]
Group_2 = Group_2_temp[1:k2,:]

Group_2[1,:] # first coordinate in Group 2

fig_new = Figure()
ax_new = Axis(fig_new[:,:]; title = "", xlabel = "", ylabel = "")
scatter!(ax_new, Group_1[:,1], Group_1[:,2], color = "blue")
scatter!(ax_new, Group_2[:,1], Group_2[:,2], color = "red")
fig_new

print(Group_index) # tells us if a given point is in group one or two


#=
vec = Vector{Float64}(undef, 1000)
for i = 1:1000
vec[i] = features[i]
end
=#


using Statistics, CairoMakie, LinearAlgebra

randompoint() = (rand([-1, 1]) + 0.1randn(), 0.1randn())
features = [randompoint() for _ in 1:1000]
#=
typeof(features)
f = features[1]
f[2]
=#
x = Vector{Float64}(undef, 1000)
y = Vector{Float64}(undef, 1000)
for i = 1:1000
x[i], y[i] = features[i][1], features[i][2]
end
x
typeof(x)

fig = Figure()
ax = Axis(fig[:,:]; title = "", xlabel = "", ylabel = "")
scatter!(ax, x, y)
fig

# from the scatter plot (also consider how points are generated) it's clear a threshold of around 0.5 is appropriate (approx.)

#features[1:1000][1] # not working as intended

maximum(x)
minimum(x)
maximum(x)-minimum(x)
Range = maximum(x)-minimum(x)

#=
Group_1_temp = Vector{Vector{Float64}}(undef, 1000)
Group_2_temp = Vector{Vector{Float64}}(undef, 1000)
k1, k2 = 1, 0
# Base case -- set the first group as defined by the intiial point
Group_1_temp[1] = [x[1], y[1]]
r = 0.4 # threshold distance
for i in 2:1000
    distance = norm([x[i], y[i]] - [x[i-1], y[i-1]])
if distance < r
    k1 += 1
    Group_1_temp[k1] = [x[i], y[i]]
else
    k2 += 1
    Group_2_temp[k2] = [x[i], y[i]]
end
end
# extract relevant info
Group_1 = Group_1_temp[1:k1]
Group_2 = Group_2_temp[1:k2]
=#
#################################
Group_1_temp = Matrix{Float64}(undef, 1000, 2)
Group_2_temp = Matrix{Float64}(undef, 1000, 2)
Group_index = Vector{Int16}(undef, 1000)# tells us what group a vector in the list belongs to
k1, k2 = 1, 0
# Base case -- set the first group as defined by the intiial point
Group_1_temp[1,:] = [x[1], y[1]]
Group_index[1] = 1
r = 0.5 # threshold distance
for i in 2:1000
    distance = norm([x[i], y[i]] - [x[1], y[1]])
if distance < r
    k1 += 1
    Group_1_temp[k1,:] = [x[i], y[i]]
    Group_index[i] = 1
else
    k2 += 1
    Group_2_temp[k2,:] = [x[i], y[i]]
    Group_index[i] = 2
end
end
# extract relevant info
Group_1 = Group_1_temp[1:k1,:]
Group_2 = Group_2_temp[1:k2,:]

Group_2[1,:] # first coordinate in Group 2

fig_new = Figure()
ax_new = Axis(fig_new[:,:]; title = "", xlabel = "", ylabel = "")
scatter!(ax_new, Group_1[:,1], Group_1[:,2], color = "blue")
scatter!(ax_new, Group_2[:,1], Group_2[:,2], color = "red")
fig_new

print(Group_index) # tells us if a given point is in group one or two


#=
vec = Vector{Float64}(undef, 1000)
for i = 1:1000
vec[i] = features[i]
end
=#