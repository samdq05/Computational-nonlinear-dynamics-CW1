# initialise
prange = vcat(-1.1:0.01:-0.8, -0.15:0.01:0.05)
values = 21
vg = range(-2.5, 2.5; length = values)
wg = range(-2.5, 2.5; length = values)
ics = [[v, w] for w in wg for v in vg]