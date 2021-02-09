using Jolab

function data(λ, m)
	# Initialization of refractive indexes to model the FP etalon
	n_fusedsilica = Jolab.refractiveindex_fusedsilica(printBool = false)
	n_cryo = Jolab.refractiveindex_cryolite(printBool = false)
	n_zns = Jolab.refractiveindex_zincsulfide(printBool = false)
	n_air = Jolab.refractiveindex_air(printBool = false)

	## Defining the optical system in Jolab. The optical system and FP etalon is based on https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-28-5-7691&id=427957
	f_col = 10E-3
	fibre = SingleModeFibre(10E-6, n_air, 1, ReferenceFrame(0,0,0.))
	collimator = Lens(f_col, 1, ReferenceFrame(0,0,f_col))
	f_obj = f_col * m
	objective = Lens(f_obj, 1, ReferenceFrame(0,0,2f_col + f_obj))

	# Sampling of the plane wave directions of propagation
	# sx = sin(θ)cos(ϕ)
	# sy = sin(θ)sin(ϕ)
	nsx = range(-0.15, 0.15, length = 150)
	nsy = range(-0.15, 0.15, length = 150)

	# Direction of the field in Jolab
	dir = 1

	# Cavity thickness
	h = 102E-6

	# Design wavelength of the dielectric mirror
	λ_des = 1402E-9

	# Creation of the FP etalon as a multilayer structure
	ref_1 = ReferenceFrame(0.,0.,2f_obj + 2f_col)
	n_M = cat(n_air, repeat([n_cryo, n_zns], 6), n_fusedsilica, repeat([n_zns, n_cryo], 6), dims = 1)
	h_M = λ_des ./ 4 ./ [n_M[i](λ_des) for i in 2:length(n_M)-1]
	h_M[13] = h
	fp_mls = MultilayerStructure(n_M, h_M, ref_1)

	# Creation of the FP etalon as two ideal mirrors
	# Performing simulations to calculate the reflectivity of the mirrors
	nsx_col = range(-0.05, 0.05, length = 100)
	nsy_col = range(-0.05, 0.05, length = 100)
	field = FieldAngularSpectrum_gaussian(nsx_col, nsy_col, 200E-6, λ[1], n_fusedsilica(λ[1]), dir, ref_1)
	mirror1 = MultilayerStructure(n_M[14:end], h_M[14:end], ref_1)
	(R,tmp) = lightinteraction(mirror1, field)
	@show R = intensity(R) / intensity(field)

	# Creating the FP etalon based on mirrors with reflectivity R
	ref_2 = ReferenceFrame(0.,0.,h)
	mirror1_airy = Mirror(R, n_air, n_fusedsilica, ref_1)
	mirror2_airy = Mirror(R, n_fusedsilica, n_air, ref_1 + ref_2)
	fp_airy = [mirror1_airy, mirror2_airy]

	# Initialization of the arrays to store the results
	itf_t_mls = zeros(length(λ))
	itf_r_mls = zeros(length(λ))
	itf_r_airy = zeros(length(λ))
	itf_t_airy = zeros(length(λ))

	for i in eachindex(λ)
		## Propagation of the fields by the optical system in Jolab
		field = FieldAngularSpectrum_fromfibre(fibre, nsx, nsy, λ[i])
    	(~, field) = lightinteraction(collimator, field)
    	(~, field) = lightinteraction(objective, field)
    	(rfield, tfield) = lightinteraction(fp_mls, field)
    	itf_t_mls[i] = intensity(tfield)

    	(rfield, ~) = lightinteraction(objective, rfield)
    	(rfield, ~) = lightinteraction(collimator, rfield)
    	itf_r_mls[i] = signal(fibre, rfield)

    	(rfield, tfield) = lightinteraction(fp_airy, field)
    	itf_t_airy[i] = intensity(tfield)

    	(rfield, ~) = lightinteraction(objective, rfield)
    	(rfield, ~) = lightinteraction(collimator, rfield)
    	itf_r_airy[i] = signal(fibre, rfield)
	end
	return (itf_r_airy, itf_t_airy, itf_r_mls, itf_t_mls)
end

# Number of wavelengths per ITF
sizeλ = 200

# Initialization of the arrays to store results
itf_t_mls = zeros(4, sizeλ)
itf_r_mls = zeros(4, sizeλ)
itf_r_airy = zeros(4, sizeλ)
itf_t_airy = zeros(4, sizeλ)
λ = zeros(4, sizeλ)

# Compute the ITF on a specific wavelength range and spot size
λ[1,:] = range(1400E-9, 1407E-9, length = sizeλ)
(itf_r_airy[1,:], itf_t_airy[1,:], itf_r_mls[1,:], itf_t_mls[1,:]) = data(λ[1,:], 2.5)

# Compute the ITF on a specific wavelength range and spot size
λ[2,:] = range(1500E-9, 1507E-9, length = sizeλ)
(itf_r_airy[2,:], itf_t_airy[2,:], itf_r_mls[2,:], itf_t_mls[2,:]) = data(λ[2,:], 5)

# Compute the ITF on a specific wavelength range and spot size
λ[3,:] = range(1549.5E-9, 1556.5E-9, length = sizeλ)
(itf_r_airy[3,:], itf_t_airy[3,:], itf_r_mls[3,:], itf_t_mls[3,:]) = data(λ[3,:], 7.5)

# Compute the ITF on a specific wavelength range and spot size
λ[4,:] = range(1600E-9, 1607E-9, length = sizeλ)
(itf_r_airy[4,:], itf_t_airy[4,:], itf_r_mls[4,:], itf_t_mls[4,:]) = data(λ[4,:], 10)

# Ploting the data for visualization
plot(λ[1,:], itf_r_airy[1,:])
plot!(λ[1,:], itf_t_airy[1,:])
plot!(λ[1,:], itf_r_mls[1,:])
plot!(λ[1,:], itf_t_mls[1,:])
