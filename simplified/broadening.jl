"""
The main purpose of this function is to avoid the broadening function 
depending on a global variable σ,
whose type is variable, which makes it slow.
"""
function gaussian_broadening(σ::Float64)
    broadening(x::Float64) = @fastmath exp(- σ^2 * x^2)
    broadening
end