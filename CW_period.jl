# Estimate period of trajectory A
function period(A)
B = [round.(u; digits = 1) for u in A]
C = unique(B)
return length(C)
end