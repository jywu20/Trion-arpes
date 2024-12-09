include("model-data.jl")
include("units.jl")

function E_c1(trion::Intervalley2DChandraTrion, k_e::SVector{2, Float64})
    inv_eV * k_e' * k_e / 2trion.m_e + trion.E_g
end

function E_c2(trion::Intervalley2DChandraTrion, k_e::SVector{2, Float64})
    inv_eV * (k_e - trion.w)' * (k_e - trion.w) / 2trion.m_e + trion.E_g
end

function E_v1(trion::Intervalley2DChandraTrion, k_h::SVector{2, Float64})
    inv_eV * k_h' * k_h / 2trion.m_h
end

function E_v2(trion::Intervalley2DChandraTrion, k_h::SVector{2, Float64})
    inv_eV * (k_h - trion.w)' * (k_h - trion.w) / 2trion.m_h
end