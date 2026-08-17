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