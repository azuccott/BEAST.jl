#=function minmax1d(vertex,edge) #perche questa e' gia inclusa in edgevertexinteraction?
        T = eltype(vertex)
        m = norm(vertex-edge[1])
        M = m
        s=edge[1]-edge[2]
        s/=norm(s)
        ev1=edge[1]-vertex
        x0=(edge[1]-dot(ev1,s)*s)
        
        a=dot((edge[2]-x0),s)
        b=dot((edge[1]-x0),s)
        if a<=0 && b>=0
           m=norm(vertex-x0)
           if abs(a)<abs(b)
                M=norm(vertex-edge[2])
           end   
            
        else
            for j in 1:length(edge)
                q = edge[j]
                d = norm(vertex-q)
                d < m && (m=d)
                d > M && (M=d)
            end
        end
        return m, M
end =# #ce errore in calcolo di d max ma sembrerebbe corretto invece in TimeDomainBEMInt (edgevvertrexgeo)

function rings1d(τ, σ, ΔR)
	m, M = minmax1d(τ, σ)
	r0 = floor(Int, m/ΔR+1) #non include anello piu interno perche poi si calcola: anello i meno anello i-1
	r1 = ceil(Int, M/ΔR) #include anello piu esterno
	return r0:r1
end


#TODO risolvere problema delta R non include nell'original quaddata (4d quaddata)#forse vanno corretti i tipi
#function quaddata1D(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
 #   testels::Vector{Simplex{3,0,3,1,T}}, trialels::Vector{Simplex{3,1,1,2,T}}, timeels, quadstrat::AllAnalyticalQStrat,ΔR) where T
 function quaddata1D(
    testels::Vector{SVector{3,T}}, trialels::Vector{CompScienceMeshes.Simplex{3,1,2,2,T}},ΔR) where T    
    #dimU=dimension(testels[1])
    #dimV=dimension(trialels[1])
    #rigerenerare delta R
    #quaddata 1D
    #@assert dimU+dimV==1
    #testelsboundary=skeleton(testels,dimU-1)
   # trialelsboundary=skeleton(trialels,dimV-1) 
    
        numnodes=length(testels)
        numedges=length(trialels)

        datavertexedge=Array{TimeDomainBEMInt.edgevertexgeo, 2}(undef, numnodes, numedges)
        rings=Array{UnitRange{Int},2}(undef, numnodes, numedges)
        datarings=Array{Vector{Tuple{Int,Vector}},2}(undef,numnodes,numedges)#il type va bene

        #fill datarings with zeross !

        for p in 1:numnodes
            τ = testels[p]#testels[p]
            for q in 1:numedges
                σ = trialels[q]
                edgevertgeo=TimeDomainBEMInt.edgevertexinteraction(τ,σ[1],σ[2])
                datavertexedge[p,q]=edgevertgeo
               # a,b=edgevertgeo.extint0[1],edgevertgeo.extint0[2]
                 # a cosa servono?
                mind=edgevertgeo.dmin
                maxd=edgevertgeo.dmax
                r0 = floor(Int, mind/ΔR+1) #non include anello piu interno perche poi si calcola: anello i meno anello i-1
	            r1 = ceil(Int, maxd/ΔR) 
                rngs=r0:r1
                rings[p,q]=rngs
                datarings[p,q]=[(0,[0.0,0.0])]
                for r in rngs
                    #r > numfunctions(timebasisfunction) && continue #serve? era in quaddata originale
                    #ι = ring(r,ΔR)#ma serve? se poi prendo solo il num 2

                    rp=τ #se e un simplex ok se no va messo chart(τ,1).vertices credo
                    #t2=ι[2]#needs a check
                    extint=TimeDomainBEMInt.edgevertexinteraction(r,edgevertgeo)
                    push!(datarings[p,q],extint)
    
                    # qr = quadrule(op, U, V, W, p, τ, q, σ, r, ι, qd, quadstrat)
                    #momintegrals!(z, op, U, V, W, τ, σ, ι, qr)
                end
            end
        end

        return datavertexedge, rings, datarings
    #else
     #   return "devo ancora scrivere"
    #end
end

function quaddata2D_edg_edg(
    testels::Vector{CompScienceMeshes.Simplex{3,1,2,2,T}}, trialels::Vector{CompScienceMeshes.Simplex{3,1,2,2,T}},ΔR,quaddata1D,connectivity) where T
    #(op::AcusticSingleLayerTDIO, testrefs, trialrefs, timerefs,
    #testels::Vector{CompScienceMeshes.Simplex{3,1,1,2,T}}, trialels::Vector{CompScienceMeshes.Simplex{3,1,1,2,T}}, timeels, quadstrat::AllAnalyticalQStrat) where T
    
    #nnodes=length(nodes) #i nodes sono salvati?

    numedges=length(testels)
    totrings=Array{UnitRange{Int},2}(undef, numedges, numedges)
    datavalues=Array{Vector{Tuple{Tuple,Tuple}},2}(undef,numedges,numedges)
    cnnct=connectivity#connectivity(edges,nodes)
    #edgevertexgeo,rings,datarings vanno specificati ancora. suppongo(quasi certamente) sia sufficiente prenderli come gli output di quaddata1D
    nnodes=size(cnnct,1)
    edgevertexgeo,rings,datarings=quaddata1D[1],quaddata1D[2],quaddata1D[3]
    for p in 2:numedges
        for q in 1:(p-1) #escluso p=q perche non ci serve da calcolare
            edge1=testels[p] #chart(edges,p)
            edge2=trialels[q] #chart(edges,q)
            
            vertind1=cnnct[1:nnodes,p].nzind
            vertsgn1=cnnct[1:nnodes,p].nzval
            vertind2=cnnct[1:nnodes,q].nzind
            vertsgn2=cnnct[1:nnodes,q].nzval

            if vertsgn1[1]==1
                a1,a2=edge1[1],edge1[2]
            else
                a2,a1=edge1[1],edge1[2]
            end

            if vertsgn2[1]==1
                b1,b2=edge2[1],edge2[2]
            else
                b2,b1=edge2[1],edge2[2]
            end

            geo1,rings1,datarings1=edgevertexgeo[vertind1[1],q],rings[vertind1[1],q],datarings[vertind1[1],q]
            geo2,rings2,datarings2=edgevertexgeo[vertind1[2],q],rings[vertind1[2],q],datarings[vertind1[2],q]
            geo3,rings3,datarings3=edgevertexgeo[vertind2[1],p],rings[vertind2[1],p],datarings[vertind2[1],p]
            geo4,rings4,datarings4=edgevertexgeo[vertind2[2],p],rings[vertind2[2],p],datarings[vertind2[2],p]

            geo=[geo1,geo2,geo3,geo4]
            rings_edgedg=[rings1,rings2,rings3,rings4]
            datarings_edgedg=[datarings1,datarings2,datarings3,datarings4]
            print([p,q]," sto per chiamare linelineglobal ")
            totrings[p,q],datavalues[p,q]=intlinelineglobal(a1,a2,b1,b2,geo,rings_edgedg,datarings_edgedg,[10^6,10^6,10^6],ΔR,Val{0}) 
        end
    end
    return totrings,datavalues
end


function intlinelineglobal(a1,a2,b1,b2,geo,rings,datarings,parcontrol,ΔR,UB::Type{Val{N}}) where N
        
    #nedges=length(edges)
   
    
     #vertices=[a1,a2,a1′,a2′]
    vertices=[a1,a2,b1,b2]
    l12=norm(vertices[1]-vertices[2])
    l12′=norm(vertices[3]-vertices[4])
    

   #geo1=edgevertexinteraction(a1,a1′,a2′) 
    #geo2=edgevertexinteraction(a2,a1′,a2′)
    #geo3=edgevertexinteraction(a1′,a1,a2)
    #geo4=edgevertexinteraction(a2′,a1,a2)
    
    #datatime=Array{Tuple}(undef,2,4)?
    I = WiltonInts84.maketuple(eltype(a1), UB)
    K = WiltonInts84.maketuple(typeof(a1), UB)
    
    x=geo[3].tangent

    #z=cross(a12′,x)
    #J=norm(z)
    #if J blablabla
    #z /= J
    #h=dot(r22,z)
    hdir=cross(geo[1].tangent,x)
    n=hdir/norm(hdir)
        sgnn=[+1,-1,-1,+1]
        h=dot(a2-b2,n)
        sgnh=[+1,-1,+1,-1]
        angletot=0.0
        dminv=Vector{eltype(a1)}(undef, 4)
        dmaxv=Vector{eltype(a1)}(undef, 4)
        ξ=Vector{typeof(a1)}(undef, 4)
        for j in 1:4
            dminv[j]=geo[j].dmin
            dmaxv[j]=geo[j].dmax 
            v=vertices[j]
            ξ[j]=v-n*h*sgnh[j]*sgnn[j] 
            angletot+=TimeDomainBEMInt.anglecontribution(ξ[j],sgnn[j]*n,geo[j])
        end
    if abs(angletot-2π)<100*eps(eltype(a1))
        dmin=abs(h)
    else
        dmin=min(dminv[1],dminv[2],dminv[3],dminv[4])
    end

    dmax=max(dmaxv[1],dmaxv[2],dmaxv[3],dmaxv[4])

    r0 = floor(Int, dmin/ΔR+1) 
	r1 = ceil(Int, dmax/ΔR)#+1) #recuperare deltaR
	ringtot = r0 : r1
    allint=Vector{typeof((I,K))}(undef,r1-r0+1)
    fill!(allint,(I,K))
    if norm(hdir) < (parcontrol[1])*eps(eltype(a1))
        
        I=intparallelsegment(a1,a2,b1,b2,temp1,temp2)[1] #TODO attenzione qui non compatibile con quello che stiamo scrivendo
    else
        n=hdir/norm(hdir)
        sgnn=[+1,-1,-1,+1]
        h=dot(a2-b2,n)
        sgnh=[+1,-1,+1,-1]
        for j in 1:4  
            for i in ringtot[1]:(rings[j][1]-1)
                print(" sono nei ring totali ")
                        v=vertices[j]
                        indx=i-ringtot[1]+1 
                        print(" ",i," ",h)
                        P, Q = TimeDomainBEMInt.arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],0,[0,0],i*ΔR,UB)
                        P=TimeDomainBEMInt.multiply(P,1/(l12′*l12*norm(hdir)))
                        Q=TimeDomainBEMInt.multiply(Q,1/(l12′*l12*norm(hdir)))

                        P=TimeDomainBEMInt.add(allint[indx][1],P)
                        Q=TimeDomainBEMInt.add(allint[indx][2],Q)
                        allint[indx]=(P,Q)
                        P, Q = TimeDomainBEMInt.arcsegcontribution(v,ξ[j],-sgnn[j]*n,sgnh[j]*h,geo[j],0,[0,0],(i-1)*ΔR,UB)
                        P=TimeDomainBEMInt.multiply(P,1/(l12′*l12*norm(hdir)))
                        Q=TimeDomainBEMInt.multiply(Q,1/(l12′*l12*norm(hdir)))
                        P=TimeDomainBEMInt.add(allint[indx][1],P)
                        Q=TimeDomainBEMInt.add(allint[indx][2],Q)
                        allint[indx]=(P,Q)
            end
            
            for i in rings[j]
                print(" sono nei ring locali ")
                    #shall I put some check like i*deltaR > h
                    v=vertices[j]
                    indx=i-ringtot[1]+1 #attenzione: il +2 invece di +1 e' perche gli anelli 1D sono sempre uno in piu e il primo indice invece e' sempre (0,[0,0])
                    indx_rings=i-rings[j][1]+1
                    #print(" indx_rings= ",indx_rings," indx= ",indx," ringtot=",ringtot," ringsj=",rings[j]," dataringsj=",datarings[j]," ")  
                    P, Q = TimeDomainBEMInt.arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],datarings[j][indx_rings+1][1],datarings[j][indx_rings+1][2],i*ΔR,UB)  
                        #allint[indx]=
                        P=TimeDomainBEMInt.multiply(P,1/(l12′*l12*norm(hdir)))
                        Q=TimeDomainBEMInt.multiply(Q,1/(l12′*l12*norm(hdir)))
                        P=TimeDomainBEMInt.add(allint[indx][1],P)
                        Q=TimeDomainBEMInt.add(allint[indx][2],Q)
                        allint[indx]=(P,Q)
                        P, Q = TimeDomainBEMInt.arcsegcontribution(v,ξ[j],-sgnn[j]*n,sgnh[j]*h,geo[j],datarings[j][indx_rings][1],datarings[j][indx_rings][2],(i-1)*ΔR,UB)
                        P=TimeDomainBEMInt.multiply(P,1/(l12′*l12*norm(hdir)))
                        Q=TimeDomainBEMInt.multiply(Q,1/(l12′*l12*norm(hdir)))
                        P=TimeDomainBEMInt.add(allint[indx][1],P)
                        Q=TimeDomainBEMInt.add(allint[indx][2],Q)
                        allint[indx]=(P,Q)
                        
            end 
            #probabilmente va bene cosi anche con ceil(int,frac) invece di ceil(int,frac+1)
            
            #=   for i in (relrings[j][2]+1):ringtot[j][2]
        
                    #shall I put some check like i*deltaR > h
                        saveP[i],saveQ[i]  = arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],1,[a,b],i*ΔR,UB) #dadefinire a e b 
                        save,and,subtract = arcsegcontribution(v,ξ[j],sgnn[j]*n,sgnh[j]*h,geo[j],1,[a,b],(i-1)*ΔR,UB)
            end =# #sembra che non serva a causa di ceil(int,frac+1)!
        end
            #I+=(1/abs(r12′[2]*r12[1]))*(3*(temp2^2-temp1^2)*d[1]-2*(temp2^3-temp1^3)*d[2])
    end
    return ringtot,allint #missing buidgrad since it is not yet adapted for int line line ci serve??
end