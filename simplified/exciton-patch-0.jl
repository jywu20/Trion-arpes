# Plot the ARPES heatmap;

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("exciton.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
include("overlap.jl")
using CairoMakie
using LaTeXStrings
using Colors

caption_padding = 35

broaden = gaussian_broadening(20fs)
electron_color = colorant"deepskyblue2"
hole_color = colorant"coral2"

m_h = 0.21
m_e = 0.37
E_g = 2.84
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = m_h,
    m_e = m_e,
    w = w,
    E_g = E_g,
    E_B = 0.75, # Binding energy for the ehh trion mode
    a = 10.3,
    b = 25.2
)
exciton_direct = IntraValley2DExciton(
    m_h = m_h,
    m_e = m_e,
    E_g = E_g,
    E_B = 0.71,
    a = 10.4
)

# The trion momentum is set to be w
P_ratio = 1.0
P = P_ratio * w

# The exciton momentum is set to zero
Q = SA[0.0, 0.0]

rk, Avck = read_ex_wfc("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5", SVector{3, Float64}(0.333333333333333, 0.3333333333333, 0))

for iS in 1:30
    let f = Figure(size=(3000, 500))
        ax = Axis(f[1, 1], title="S = $iS")
        scatter!(ax, rk[1, :], rk[2, :], color=abs.(Avck[1, 1, :, iS]), colormap=reverse(cgrad(:grayC)), colorrange=(0, 1),)
        colsize!(f.layout, 1, Aspect(1, 1))
        hidedecorations!(ax, ticklabels = false, ticks = false)
        
        ax = Axis(f[1, 2])
        Avck_norm = map(1 : size(Avck)[3]) do ik
            norm(Avck[:, :, ik, iS])^2
        end
        scatter!(ax, rk[1, :], rk[2, :], color=abs.(Avck[1, 2, :, iS]), colormap=reverse(cgrad(:grayC)), colorrange=(0, 1),)
        colsize!(f.layout, 2, Aspect(1, 1))
        hidedecorations!(ax, ticklabels = false, ticks = false)
        
        Avck_iS_first_two_bands = Avck[1, 1, :, iS] + Avck[1, 2, :, iS]
        Avck_1s = Avck[1, 1, :, 2] + Avck[1, 2, :, 2]
        phase_factors = Avck_1s ./ abs.(Avck_1s)
        Avck_iS_first_two_bands = Avck_iS_first_two_bands ./ phase_factors
        Avck_iS_first_two_bands ./= Avck_iS_first_two_bands[1] / abs(Avck_iS_first_two_bands[1])

        ax = Axis(f[1, 3], title="norm = $(norm(Avck_iS_first_two_bands))")
        scatter!(ax, rk[1, :], rk[2, :], color=abs.(Avck_iS_first_two_bands), 
            colormap=reverse(cgrad(:grayC)),
            colorrange=(0, 1),
        )
        colsize!(f.layout, 3, Aspect(1, 1))
        hidedecorations!(ax, ticklabels = false, ticks = false)
        
        #Colorbar(f[1, 4], colormap=reverse(cgrad(:grayC)), )

        ax = Axis(f[1, 4], title="norm = $(norm(real.(Avck_iS_first_two_bands)))")
        intensity_limit = maximum(abs.(Avck_iS_first_two_bands))
        scatter!(ax, rk[1, :], rk[2, :], color=real.(Avck_iS_first_two_bands), 
            colormap=cgrad(:balance),
            colorrange=(-intensity_limit, intensity_limit),
        )
        colsize!(f.layout, 4, Aspect(1, 1))
        hidedecorations!(ax, ticklabels = false, ticks = false)
        
        ax = Axis(f[1, 5], title="norm = $(norm(imag.(Avck_iS_first_two_bands)))")
        intensity_limit = maximum(abs.(Avck_iS_first_two_bands))
        scatter!(ax, rk[1, :], rk[2, :], color=imag.(Avck_iS_first_two_bands), 
            colormap=cgrad(:balance),
            colorrange=(-intensity_limit, intensity_limit),
        )
        colsize!(f.layout, 5, Aspect(1, 1))
        hidedecorations!(ax, ticklabels = false, ticks = false)
        
        ax = Axis(f[1, 6])
        scatter!(ax, map(ik -> norm(rk[1:2, ik]), 1 : size(rk)[2]), real.(Avck_iS_first_two_bands))
        
        save("exciton-0-patch-$iS.png", f)

    end
end
