"""
    function CustomLinearRoquetSeawaterPolynomial(Θ, S, FT=Float64)
Return a linear `RoquetSeawaterPolynomial` with custom values for α and β which are computed
using GibbsSeaWater.jl from `Θ` and `S` (`S` is practical salinity).
"""
function CustomLinearRoquetSeawaterPolynomial(Θ, S, FT)

    α = gsw_alpha(S, Θ, 0)
    β = gsw_beta(S, Θ, 0)

    return SecondOrderSeawaterPolynomial{FT}(R₀₁₀ = α,
                                             R₁₀₀ = β)
end

"Extend `RoquetSeawaterPolynomial` for `:CustomLinear` coefficent set"
SecondOrderSeawaterPolynomials.RoquetSeawaterPolynomial(Θ, S, FT=Float64, coefficient_set=:CustomLinear) =
    eval(Symbol(coefficient_set, :RoquetSeawaterPolynomial))(Θ, S, FT)

"Constructor for `CustomLienarEquationOfState`. The `Θ` and `S` inputs are the values that
α and β (the coefficients in the equation of state) will be calculated at."
CustomLinearEquationOfState(Θ, S; reference_density = 1024.6) =
    BoussinesqEquationOfState(RoquetSeawaterPolynomial(Θ, S), reference_density)
