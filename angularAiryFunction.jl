using HCubature, Plots, Jolab

# This script compares the results from the equations of [1] and compare it
# with the results computed from Jolab

# Parameters of the simulation
ω = 50E-6 # Gaussian spot size
R1 = .97 # First mirror reflectivity
R2 = .97 # Second mirror reflectivity
N = 1.5 # refractive index in the cavity
h = 100E-6 # cavity thickness

# Angular spectrum of a Gaussian beam
gauss(θ,ϕ,λ) = exp(-(2π)^2 / 16 / λ^2 * ω^2 * N^2 * sin(θ)^2) / λ

# Airy function in reflection (Equation 1 of [1])
r_airy(θ,λ) = √R1 + (1-R1) * √R2 * exp(im * 4π / λ * N * h * cos(θ)) / (1 + √R1 * √R2 * exp(im * 4π / λ * N * h * cos(θ)))

# Airy function in transmission (Equation 1 of [1])
t_airy(θ,λ) = √(1 - R1) * √(1 - R2) * exp(im * 2π / λ * N * h * cos(θ)) / (1 + √R1 * √R2 * exp(im * 4π / λ * N * h * cos(θ)))

# Reflection Angular airy function (Equation 2 of [1])
r_angularairy(θ,ϕ,λ) = gauss(θ,ϕ,λ) * r_airy(θ,λ)

# transmission Angular airy function (Equation 2 if [1])
t_angularairy(θ,ϕ,λ) = gauss(θ,ϕ,λ) * t_airy(θ,λ)

# Intensity measured with big detector
function i_infinite(Efp, λ)
	aux(θϕ) = abs2(Efp(θϕ[1],θϕ[2], λ)) * sin(θϕ[1])

	# hcubature performs numerical integration of the function
	return hcubature(aux, [0.; 0.], [π/2; 2π])[1][1]
end

# Calculation of an ITF. The data is not normalized
λ = range(1557E-9, 1559E-9, length = 100)
i_t = [i_infinite(t_angularairy, λi) for λi in λ]
i_r = [i_infinite(r_angularairy, λi) for λi in λ]

ref_1 = ReferenceFrame(0,0,0.)
ref_2 = ReferenceFrame(0.,0.,h)
mirror1_airy = Mirror(R1, 1, N, ref_1)
mirror2_airy = Mirror(R2, N, 1, ref_2)
fp_airy = [mirror1_airy, mirror2_airy]

# Sampling of the plane wave directions of propagation
# sx = sin(θ)cos(ϕ)
# sy = sin(θ)sin(ϕ)
nsx = range(-0.05, 0.05, length = 200)
nsy = range(-0.05, 0.05, length = 200)
jolab_r = zeros(length(λ))
jolab_t = zeros(length(λ))

for i in eachindex(λ)
	field = FieldAngularSpectrumScalar_gaussian(nsx, nsy, ω, λ[i], 1., 1, ReferenceFrame(0,0,0.))
	(rfield, tfield) = lightinteraction(fp_airy, field)
	jolab_r[i] = intensity(rfield)
	jolab_t[i] = intensity(tfield)
end

plot(λ, i_t ./ maximum(i_r))
plot!(λ, i_r ./ maximum(i_r))
plot!(λ, jolab_r)
plot!(λ, jolab_t)
