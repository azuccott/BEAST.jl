using CompScienceMeshes, BEAST, LinearAlgebra
#Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
#CompScienceMeshes.meshsphere(1.0,0.4)

Γ = CompScienceMeshes.meshcuboid(1.0,1.0,1.0,1.0/6)

X = lagrangecxd0(Γ)

vec1=zeros(numfunctions(X))
vec1[200]=1.0
vec1[300]=1.0


par=2
Δt, Nt = 0.1039049/par, 350*par
T = timebasisshiftedlagrange(Δt, Nt, 0)
#dT=derive(T)

U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U


width, delay, scaling =  80*Δt, 24.0, 1.0*80*Δt/8
gaussian = creategaussian(width, delay, scaling)
e = BEAST.planewave(point(0,0,1), 1.0, gaussian)

@hilbertspace j
@hilbertspace j′

SL = TDAcustic3D.acusticsinglelayer(speedofsound=1.0, numdiffs=0)
#BEAST.@defaultquadstrat (SL, W, V) #nothing #BEAST.AllAnalyticalQStrat
#BEAST.OuterNumInnerAnalyticQStrat(7)

tdacusticsl = @discretise SL[j′,j] == -1.0e[j′]   j∈V  j′∈W
xacusticsl = solve(tdacusticsl)


#@time Z = assemble(SL,W,V)

#corregere da qui 
import Plots

xacusticsl

Plots.plot(xtdhhsl[1,:],label=["Current_wilton_rule2"])#800,1400


Plots.plot(xacusticsl[1,:],label="Current_exact")

Plots.plot([xacusticsl[1,:],xtdhhsl[1,:]],label=["Solution_exact" "Solution_semi-analytical"],xlim=(200*par,300*par))#800,1400

Plots.xlabel!("t")
Plots.savefig("c:/Users/Administrator/Pictures/Plots_mot/new_mot_exact_vs_wilton_7.png")
#Plots.savefig("mot_exact_wilton_4.png") 

xacusticsl

pval=ConvolutionOperators.polyvals(Z)
findmax(norm(pval[i]) for i in 1:size(pval,1))

import Plotly
#fcr, geo = facecurrents(xefie[:,125], X)
#Plotly.plot(patch(geo, norm.(fcr)))

import BEAST.ConvolutionOperators
ConvolutionOperators.timeslice(Zs,1)
ConvolutionOperators.timeslice(Z,1)


for a in 1:1
    za=ConvolutionOperators.timeslice(Z,a)
    zawilton=ConvolutionOperators.timeslice(Zs,a)
    for i in 1:10 #numfunctions(X)
        #for j in 1:numfunctions(X)
            if  norm(za[i,i]-zawilton[i,i])>=10^(-6)
                println(a," ",i," ",za[i,i]," ",zawilton[i,i])
            end
        #end
    end
end



name=readline()
parse(Float64, name)

Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)

fcr, geo = facecurrents(ue, X)
Plotly.plot(patch(geo, norm.(fcr)))

