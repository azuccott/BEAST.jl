#TODO: rename to reorder_dof and make it depending on RefSpace 
function reorder_dof(space::NDLCDRefSpace,I)
    n = length(I)
    K = zeros(MVector{4,Int64})
    for i in 1:n
        for j in 1:n
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    return SVector(K),SVector{4,Int64}(1,1,1,1)
end


function reorder_dof(space::NDLCCRefSpace, I )

    n = length(I)
    J = zeros(MVector{4,Int64})
    for i in 1:n
        for j in 1:n
            if I[j] == i
                J[i] = j
                break
            end
        end
    end

    edges = collect(combinations(J,2))
    ref_edges  = collect(combinations(SVector{4,Int64}(1,2,3,4),2))
    n = length(edges)
    K = zeros(MVector{6,Int64})
    O = zeros(MVector{6,Int64})

    for i in 1:n
        for j in 1:n
            if edges[i] == ref_edges[j]
                K[i] = j
                O[i] = 1
                break
            elseif reverse(edges[i]) == ref_edges[j]
                K[i] = j
                O[i] = -1
                break
            end
        end
    end

    return SVector(K),SVector(O)
end

function reorder_dof(space::RTRefSpace,I)
    n = length(I)
    K =  zeros(MVector{3,Int64})
    for i in 1:n
        for j in 1:n
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    return SVector(K),SVector{3,Int64}(1,1,1)
end

function reorder_dof(space::LagrangeRefSpace{T,0,3,1},I) where T
   
    return SVector{1,Int64}(1),SVector{1,Int64}(1)
end

function reorder_dof(space::LagrangeRefSpace{T,1,3,3},I) where T
    n = length(I)
    K = zeros(MVector{3,Int64})
    for i in 1:n
        for j in 1:n
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    return SVector(K),SVector{3,Int64}(1,1,1)
end

function reorder_dof(space::LagrangeRefSpace{T,1,4,4},I) where T
    n = length(I)
    K = zeros(MVector{4,Int64})
    for i in 1:n
        for j in 1:n
            if I[j] == i
                K[i] = j
                break
            end
        end
    end

    return SVector(K),SVector{4,Int64}(1,1,1,1)
end

function integrate!(out, op::VIEOperator,
    test_functions::Space, test_ptr, test_chart,
    trial_functions::Space, trial_ptr, trial_chart,
    strat::SauterSchwab3DStrategy)

    local_test_space = refspace(test_functions)
    local_trial_space = refspace(trial_functions)

    #Find permutation of vertices to match location of singularity to SauterSchwab
    J, I= SauterSchwab3D.reorder(strat.sing)
      
    #Get permutation and rel. orientatio of DoFs 
    K,O1 = reorder_dof(local_test_space, I)
    L,O2 = reorder_dof(local_trial_space, J)
    #Apply permuation to elements
 
    if length(I) == 4
        test_tetrahedron_element  = simplex(
            test_tetrahedron_element.vertices[I[1]],
            test_tetrahedron_element.vertices[I[2]],
            test_tetrahedron_element.vertices[I[3]],
            test_tetrahedron_element.vertices[I[4]])
    elseif  length(I) == 3
        test_tetrahedron_element  = simplex(
            test_tetrahedron_element.vertices[I[1]],
            test_tetrahedron_element.vertices[I[2]],
            test_tetrahedron_element.vertices[I[3]])
    end

    #test_tetrahedron_element  = simplex(test_tetrahedron_element.vertices[I]...)

    if length(J) == 4
    trial_tetrahedron_element  = simplex(
        trial_tetrahedron_element.vertices[J[1]],
        trial_tetrahedron_element.vertices[J[2]],
        trial_tetrahedron_element.vertices[J[3]],
        trial_tetrahedron_element.vertices[J[4]])
    elseif  length(J) == 3
        trial_tetrahedron_element  = simplex(
        trial_tetrahedron_element.vertices[J[1]],
        trial_tetrahedron_element.vertices[J[2]],
        trial_tetrahedron_element.vertices[J[3]])
    end

    #trial_tetrahedron_element = simplex(trial_tetrahedron_element.vertices[J]...)

    #Define integral (returns a function that only needs barycentric coordinates)
    igd = VIEIntegrand(test_tetrahedron_element, trial_tetrahedron_element,
        op, local_test_space, local_trial_space)

    #Evaluate integral
    Q = SauterSchwab3D.sauterschwab_parameterized(igd, strat)
  
    #Undo permuation on DoFs
    for j in 1 : length(L)
        for i  in 1 : length(K)
            out[i,j] += Q[K[i],L[j]]*O1[i]*O2[j]
        end
    end
    nothing
end

function momintegrals!(z, biop::VIEOperator,
    test_functions::Space, tptr, tcell,
    trial_functions::Space, bptr, bcell,
    strat::DoubleQuadRule)

    # memory allocation here is a result from the type instability on strat
    # which is on purpose, i.e. the momintegrals! method is chosen based
    # on dynamic polymorphism.
    womps = strat.outer_quad_points
    wimps = strat.inner_quad_points

    M, N = size(z)

    for womp in womps
        tgeo = womp.point
        tvals = womp.value
        jx = womp.weight

        for wimp in wimps
            bgeo = wimp.point
            bvals = wimp.value
            jy = wimp.weight

            j = jx * jy
            kernel = kernelvals(biop, tgeo, bgeo)


            igd = integrand(biop, kernel, tvals, tgeo, bvals, bgeo)
        
            for m in 1 : length(tvals)
               
                for n in 1 : length(bvals)
                  
                    z[m,n] += j * igd[m,n]
                end
            end
        end
    end

    return z
end

function integrate!(out, op::VIEOperator,
    test_functions::Space, test_ptr, test_chart,
    trial_functions::Space, trial_ptr, trial_chart,
    strat::SauterSchwab3DStrategy, qbuffer)

    return integrate!(
        out, op,
        test_functions, test_ptr, test_chart,
        trial_functions, trial_ptr, trial_chart,
        strat,
    )
end



#TODO: use trial_ptr to get the cell material using material_array[trial_ptr] ...

# function integrate!(z, op::VIEOperator,
#     test_functions::Space, test_cellptr, test_chart,
#     trial_functions::Space, trial_cellptr, trial_chart,
#     strat::DoubleQuadRule)

#     tshs = refspace(test_functions)
#     bshs = refspace(trial_functions) 

#     igd = Integrand(op, tshs, bshs, test_chart, trial_chart) 

#     womps = strat.outer_quad_points
#     wimps = strat.inner_quad_points
    
#     for womp in womps
#         tgeo = womp.point
#         tvals = womp.value
#         M = length(tvals)
#         jx = womp.weight
        
#         for wimp in wimps
#             bgeo = wimp.point
#             bvals = wimp.value
#             N = length(bvals)
#             jy = wimp.weight

#             j = jx * jy

#             z1 = j * igd(tgeo, bgeo, tvals, bvals)
#             for n in 1:N
#                 for m in 1:M
#                     z[m,n] += z1[m,n]
#             end end
#         end
#     end

#     return z
# end
