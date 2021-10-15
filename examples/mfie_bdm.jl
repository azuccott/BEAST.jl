using CompScienceMeshes, BEAST

Γ = readmesh(joinpath(dirname(@__FILE__),"sphere2.in"))
# X, Y = raviartthomas(Γ), buffachristiansen(Γ)
X = brezzidouglasmarini(Γ)
Y = brezzidouglasmarini(Γ)
# X = raviartthomas(Γ)
# Y = raviartthomas(Γ)

ϵ, μ, ω = 1.0, 1.0, 1.0; κ = ω * √(ϵ*μ)
NK, Id = BEAST.DoubleLayerRotatedMW3D(im*κ), Identity()
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*μ*ω)*curl(E)
h = n × H

@hilbertspace j
@hilbertspace m
mfie = @discretise (NK-0.5Id)[m,j] == h[m]  j∈X m∈Y
u = gmres(mfie)

include("utils/postproc.jl")
include("utils/plotresults.jl")
