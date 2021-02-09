using Jolab
using StaticArrays, HCubature
include("/home/dylan/Documents/julia/DataTreatment.jl")

function rangularAiryFunction(ω, R1, R2, n, α, h, λ)
	norm = 	ω * √(1 / 32 / π^3);
	k = 2π / λ
	kmax = √(-log(10E-5) * 16 / ω^2)
    @inline gauss(knsr²) = norm * exp((-ω^2 / 16) * knsr²[1])
	@inline knsz(knsr²) = √(complex(k^2 * (n + im * α)^2 - knsr²[1]))
    @inline phaseTerm(knsr²) = exp(2im * knsz(knsr²[1]) * h)
    @inbounds rgauss(knsr) = (knsr[1] * gauss(knsr[1]^2)^2 * (√R1 + (1 - R1) * √R2 * phaseTerm(knsr[1]^2) / (1 + √R1 * √R2 * phaseTerm(knsr[1]^2))))

	return 64π^6 * abs2(hcubature(rgauss, SVector(0.), SVector(kmax))[1][1])
end

function tangularAiryFunction(ω, R1, R2, n, α, h, λ)
	norm = 	ω * √(1 / 32 / π^3);
	k = 2π / λ
	kmax = √(-log(10E-5) * 16 / ω^2)
    @inline gauss(knsr²) = norm * exp((-ω^2 / 16) * knsr²[1])
	@inline knsz(knsr²) = √(complex(k^2 * (n + im * α)^2 - knsr²[1]))
    @inline phaseTerm(knsr²) = exp(2im * knsz(knsr²[1]) * h)
    @inbounds tgauss(knsr) = (knsr[1] * abs2(gauss(knsr[1]^2) * (√(1 - R1) * √(1 - R2) * √(phaseTerm(knsr[1]^2)) / (1 + √R1 * √R2 * phaseTerm(knsr[1]^2)))))

	return 8π^3 * abs(hcubature(tgauss, SVector(0.), SVector(kmax))[1][1])
end


function data(λ, ω)

	nsx = range(-0.1, 0.1, length = 200)
	dir = 1
	n_air = Jolab.refractiveindex_air(printBool = false)
	n_glass = Jolab.refractiveindex_fusedsilica(printBool = false)

	ref = ReferenceFrame(0.,0.,0.)
	n_cryo = Jolab.refractiveindex_cryolite(printBool = false)
	n_zns = Jolab.refractiveindex_zincsulfide(printBool = false)
	n_M = cat(n_air, repeat([n_cryo, n_zns], 6), n_glass, dims = 1)
	λ_des = 1402E-9
	h_M = λ_des ./ 4 ./ [n_M[i](λ_des) for i in 2:length(n_M)-1]
	mirror1 = MultilayerStructure(n_M, h_M, ref)
	# mirror1 = MultilayerStructure([1, .29359 + im * 10.252, 1], [100E-9], ref)
	h = 102E-6
	ref_2 = ReferenceFrame(0.,0.,h)
	mirror2 = MultilayerStructure(n_M[end:-1:1], h_M[end:-1:1], ref_2)
	# mirror2 = MultilayerStructure([1, .29359 + im * 10.252, 1], [100E-9], ref_2)
	# mirror2 = MultilayerStructure(n_M, h_M, ref_2)

	fp = [mirror1, mirror2]
	field = FieldAngularSpectrum_gaussian(nsx, nsx, ω, λ[1], n_air(λ[1]), dir, ref)
	fibre = SingleModeFibre(ω, n_air(λ[1]), 1, ref)
	(R,tmp) = lightinteraction(mirror1, field)
	R = intensity(R)

	itf_t_mls = zeros(length(λ))
	itf_r_mls = zeros(length(λ))
	itf_r_airy = zeros(length(λ))
	itf_t_airy = zeros(length(λ))

	Threads.@threads for i in eachindex(λ)
    	field = FieldAngularSpectrum_gaussian(nsx, nsx, ω, λ[i], n_air(λ[i]), dir, ref)
    	(rfield, tfield) = lightinteraction(fp, field)
    	itf_r_mls[i] = signal(fibre, rfield)
    	itf_t_mls[i] = intensity(tfield)
		itf_r_airy[i] = rangularAiryFunction(ω, R, R, n_glass(λ[i]), 0., h, λ[i])
		itf_t_airy[i] = tangularAiryFunction(ω, R, R, n_glass(λ[i]), 0., h, λ[i])
	end
	return (itf_r_airy, itf_t_airy, itf_r_mls, itf_t_mls)
end

ω = 50E-6
sizeM = 100
itf_t_mls = zeros(4, sizeM)
itf_r_mls = zeros(4, sizeM)
itf_r_airy = zeros(4, sizeM)
itf_t_airy = zeros(4, sizeM)

λ = range(1400E-9, 1407E-9, length = sizeM)
(itf_r_airy[1,:], itf_t_airy[1,:], itf_r_mls[1,:], itf_t_mls[1,:]) = data(λ, 25E-6)

λ = range(1500E-9, 1507E-9, length = sizeM)
(itf_r_airy[2,:], itf_t_airy[2,:], itf_r_mls[2,:], itf_t_mls[2,:]) = data(λ, 50e-6)

λ = range(1549.5E-9, 1556.5E-9, length = sizeM)
(itf_r_airy[3,:], itf_t_airy[3,:], itf_r_mls[3,:], itf_t_mls[3,:]) = data(λ, 75e-6)

λ = range(1600E-9, 1607E-9, length = sizeM)
(itf_r_airy[4,:], itf_t_airy[4,:], itf_r_mls[4,:], itf_t_mls[4,:]) = data(λ, 100e-6)


itf_r_airy = [alignxdata(λ, itf_r_mls[i,:], itf_r_airy[i,:], txlim = 7e-9) for i in 1:4]
itf_t_airy = [alignxdata(λ, itf_t_mls[i,:], itf_t_airy[i,:], txlim = 7E-9) for i in 1:4]
itf_r_mls = itf_r_mls'
itf_t_mls = itf_t_mls'
