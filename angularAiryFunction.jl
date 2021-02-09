using HCubature, Plots

# Parameters of the simulation
ω = 50E-6
mfd = 10E-6
m = 5
R1 = .97
R2 = .97
N = 1.5
h = 100E-6

# Angular spectrum of a Gaussian bea,
gauss(θ,ϕ,λ) = exp(-π^2 / λ^2 * ω^2 * sin(θ)^2) / λ

# Airy function in reflection
r_airy(θ,λ) = √R1 + (1-R1) * √R2 * exp(im * 4π / λ * N * h * cos(θ)) / (1 + √R1 * √R2 * exp(im * 4π / λ * N * h * cos(θ)))

# Airy function in transmission
t_airy(θ,λ) = √(1 - R1) * √(1 - R2) * exp(im * 2π / λ * N * h * cos(θ)) / (1 + √R1 * √R2 * exp(im * 4π / λ * N * h * cos(θ)))

# Reflection Angular airy function
r_angularairy(θ,ϕ,λ) = gauss(θ,ϕ,λ) * r_airy(θ,λ)

# transmission Angular airy function
t_angularairy(θ,ϕ,λ) = gauss(θ,ϕ,λ) * t_airy(θ,λ)

# Intensity measured with an infinite detector
i_infinite(Efp, λ) = begin
	aux(θϕ) = abs2(Efp(θϕ[1],θϕ[2], λ))
	hcubature(aux, [0.; 0.], [π/2; 2π])[1][1]
end

# Intensity measured with an fibre based readout scheme
i_fibre(Efp, λ) = begin
	aux(θϕ) = Efp(θϕ[1],θϕ[2], λ) * exp(-mfd^2 * m^2 * π^2 / 4λ^2 * sin(θϕ[1])^2)
	return abs2(hcubature(aux, [0.; 0.], [π/2; 2π])[1][1])
end

# Calculation of an ITF. The data is not normalized
λ = 1555E-9:.1E-9:1560E-9
i_t = [i_infinite(t_angularairy, λi) for λi in λ]
i_r = [i_fibre(r_angularairy, λi) for λi in λ]

plot(λ, i_t)
plot(λ, i_r)
