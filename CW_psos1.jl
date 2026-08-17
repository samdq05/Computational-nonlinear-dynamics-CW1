# Poincaré map function for 1.6
function psos1(trajectory, index, value)

    n = size(trajectory, 1) # number of points in the trajectory
    x_store, y_store, z_store = zeros(n), zeros(n), zeros(n) # empty vectors potentially as large as every point
    j = 0 # initialise counter

    # trajecotry components
    x, y, z = columns(trajectory)
    
    # For a given index (x,y,z = 1,2,3 ) and value -- gives a point on the plane and a normal vector to the plane 
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
            j += 1 # add one to count
            B = trajectory[i,:] # point after
            A = trajectory[i-1,:] # point before
            r = dot(v, (P - A)) / dot(v, (B - A)) # interpolation fraction
            C = A + r * (B - A) # interpolated point
            x_store[j], y_store[j], z_store[j] = C 
        end

    end
    
    x_final, y_final, z_final = x_store[1:j], y_store[1:j], z_store[1:j]
    return x_final, y_final, z_final
end