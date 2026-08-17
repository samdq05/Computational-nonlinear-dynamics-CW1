# Poincaré map function for 1.7
function psos2(trajectory, tvec, index, value)

    n = size(trajectory, 1) # number of points in the trajectory
    x_final = []
    y_final = []
    z_final = []
    t_final = []

   x, y, z = columns(trajectory)
    
    # For a given index (x,y,z = 1,2,3 ) and value - gives a point on the plane and a normal vector to the plane 
    if index == 1
        P = [value, 0.0, 0.0] # point
        v = [1.0, 0.0, 0.0] # normal
        coord = x

    elseif index == 2
        P = [0.0, value, 0.0] # point
        v = [0.0, 1.0, 0.0] # normal   
        coord = y

    elseif index == 3
        P = [0.0, 0.0, value] # point 
        v = [0.0, 0.0, 1.0] # normal   
        coord = z

    else
        error("index must be 1 (x), 2 (y), or 3 (z)") 
    end

    for i = 2:n

        if coord[i] >  value && coord[i-1] < value
            B = trajectory[i,:] # point after
            A = trajectory[i-1,:] # point before
            r = dot(v, (P - A)) / dot(v, (B - A)) # interpolation fraction
            C = A + r * (B - A) # interpolated point
            push!(x_final, C[1])
            push!(y_final, C[2])
            push!(z_final, C[3])
            push!(t_final, tvec[i])

        end

    end

    return x_final, y_final, z_final, t_final
end