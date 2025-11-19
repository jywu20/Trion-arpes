using CairoMakie

f = Figure(size=(400, 300))
ax = Axis(f[1, 1])
xs = LinRange(0, 10, 1000)
lines!(ax, xs, sin.(xs))

y_bottom = ax.elements[:xgridlines][1][][1][2]
y_top = ax.elements[:xoppositeline][1][][1][2]
x_left = ax.elements[:xoppositeline][1][][1][1]
x_right = ax.elements[:xoppositeline][1][][2][1]

lines!(f.scene, [x_left - 10, x_left + 10], [y_top - 5, y_top + 5])

f
