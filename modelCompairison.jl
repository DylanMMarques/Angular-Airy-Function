using Jolab

function data(λ, ω)
	# Sampling of the plane wave directions of propagation
	# sx = sin(θ)cos(ϕ)
	# sy = sin(θ)sin(ϕ)
	nsx = range(-0.1, 0.1, length = 200)
	nsy = range(-0.1, 0.1, length = 200)

	# Direction of the field in Jolab
	dir = 1

	# Initialization of refractive indexes to model the FP etalon
	n_air = Jolab.refractiveindex_air(printBool = false)
	n_fusedsilica = Jolab.refractiveindex_fusedsilica(printBool = false)
	n_cryo = Jolab.refractiveindex_cryolite(printBool = false)
	n_zns = Jolab.refractiveindex_zincsulfide(printBool = false)

	# Cavity thickness
	h = 102E-6

	# Design wavelength of the dielectric mirror
	λ_des = 1402E-9

	# Creation of the FP etalon as a multilayer structure
	ref_1 = ReferenceFrame(0.,0.,0.)
	n_M = cat(n_air, repeat([n_cryo, n_zns], 6), n_fusedsilica, repeat([n_zns, n_cryo], 6), dims = 1)
	h_M = λ_des ./ 4 ./ [n_M[i](λ_des) for i in 2:length(n_M)-1]
	h_M[13] = h
	fp_mls = MultilayerStructure(n_M, h_M, ref_1)

	# Creation of the FP etalon as two ideal mirrors
	# Performing simulations to calculate the reflectivity of the mirrors
	field = FieldAngularSpectrum_gaussian(nsx, nsy, 200E-6, λ[1], n_air(λ[1]), dir, ref)
	(R,tmp) = lightinteraction(mirror1, field)
	R = intensity(R)

	# Creating the FP etalon based on mirrors with reflectivity R
	ref_2 = ReferenceFrame(0.,0.,h)
	mirror1_airy = Mirror(R, n_air, n_fusedsilica, ref_1)
	mirror2_airy = Mirror(R, n_air, n_fusedsilica, ref_1)
	fp_airy = [mirror1_airy, mirror2_airy]

	# Initialization of the arrays to store the results
	itf_t_mls = zeros(length(λ))
	itf_r_mls = zeros(length(λ))
	itf_r_airy = zeros(length(λ))
	itf_t_airy = zeros(length(λ))

	fibre = SingleModeFibre(ω, n_air(λ[1]), 1, ref)

	Threads.@threads for i in eachindex(λ)
    	field = FieldAngularSpectrum_gaussian(nsx, nsy, ω, λ[i], n_air(λ[i]), dir, ref)
    	(rfield, tfield) = lightinteraction(fp_mls, field)
    	itf_r_mls[i] = signal(fibre, rfield)
    	itf_t_mls[i] = intensity(tfield)

    	(rfield, tfield) = lightinteraction(fp_airy, field)
    	itf_r_airy[i] = signal(fibre, rfield)
    	itf_t_airy[i] = intensity(tfield)
	end
	return (itf_r_airy, itf_t_airy, itf_r_mls, itf_t_mls)
end

# Number of wavelengths per ITF
sizeλ = 100

# Initialization of the arrays to store results
itf_t_mls = zeros(4, sizeλ)
itf_r_mls = zeros(4, sizeλ)
itf_r_airy = zeros(4, sizeλ)
itf_t_airy = zeros(4, sizeλ)
λ = zeros(4, sizeλ)

# Compute the ITF on a specific wavelength range and spot size
λ[1,:] = range(1400E-9, 1407E-9, length = sizeλ)
(itf_r_airy[1,:], itf_t_airy[1,:], itf_r_mls[1,:], itf_t_mls[1,:]) = data(λ, 25E-6)

# Compute the ITF on a specific wavelength range and spot size
λ[2,:] = range(1500E-9, 1507E-9, length = sizeλ)
(itf_r_airy[2,:], itf_t_airy[2,:], itf_r_mls[2,:], itf_t_mls[2,:]) = data(λ, 50e-6)

# Compute the ITF on a specific wavelength range and spot size
λ[3,:] = range(1549.5E-9, 1556.5E-9, length = sizeλ)
(itf_r_airy[3,:], itf_t_airy[3,:], itf_r_mls[3,:], itf_t_mls[3,:]) = data(λ, 75e-6)

# Compute the ITF on a specific wavelength range and spot size
λ[4,:] = range(1600E-9, 1607E-9, length = sizeλ)
(itf_r_airy[4,:], itf_t_airy[4,:], itf_r_mls[4,:], itf_t_mls[4,:]) = data(λ, 100e-6)

# Ploting the data for visualization
plot(λ[1,:], itf_r_airy[1,:])
plot!(λ[1,:], itf_t_airy[1,:])
plot!(λ[1,:], itf_r_mls[1,:])
plot!(λ[1,:], itf_t_mls[1,:])
