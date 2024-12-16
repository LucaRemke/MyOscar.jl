########################################
# DEFINITION OF TORIC VARIETIES
########################################
# Functions that constructs toric varieties with our ordering of the rays and maximal cones
# the optional parameter non_redundant = true is mandatory for this purpose
########################################

#---------------------------------------
# Defines the n-dimensional projective space
#---------------------------------------
function define_projective_space(n::Int)
    ray_generators = Vector{Int64}[]
    gen_0 = [-1 for i in 1:n]
    push!(ray_generators, gen_0)
    
    for i in 1:n
        gen = [j == i ? 1 : 0 for j in 1:n]
        push!(ray_generators, gen)
    end
    
    max_cones = Vector{Int64}[]
    
    for i in 1:n+1
        mcone = [j for j in 1:n+1 if j != i]
        push!(max_cones, mcone)
    end
    
    P = normal_toric_variety(max_cones, ray_generators, non_redundant = true)
    
    return P
end;

#---------------------------------------
# The projective surface
# The order of the rays differs from define_projective_space(n::Int)
#---------------------------------------
function define_projective_surface()
    ray_generators = [[0,1], [-1,-1], [1,0]];
    max_cones = [[1,2], [2,3], [3,1]];
    p2 = normal_toric_variety(max_cones, ray_generators, non_redundant = true)
    
    return p2
end;

#---------------------------------------
# Defines the r-th Hirzebruch surface
#---------------------------------------
function define_hirzebruch_surface(r::Int)    
    ray_generators = [[0,1], [-r,-1], [1,0], [-1,0]]
    max_cones = [[1,4], [4,2], [2,3], [3,1]]
    hr = normal_toric_variety(max_cones, ray_generators, non_redundant = true)

    return hr
end;

#---------------------------------------
# Defines the product P1*P1
#---------------------------------------
function define_pp()
    ray_generators = [[1,0], [-1,0], [0,1], [0,-1]];
    max_cones = [[1,3], [1,4], [2,4], [2,3]];
    pp = normal_toric_variety(max_cones, ray_generators, non_redundant = true)
    
    return pp
end;

#---------------------------------------
# Defines the product P1*P1*P1
#---------------------------------------
function define_ppp()
    ray_generators = [[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]];
    max_cones = [[1,3,5], [1,3,6], [1,4,5], [1,4,6], [2,3,5], [2,3,6], [2,4,5], [2,4,6]];
    ppp = normal_toric_variety(max_cones, ray_generators, non_redundant = true)
    
    return ppp
end;

#---------------------------------------
# Defines the pentagon
#---------------------------------------
function define_pentagon()
    ray_generators = [[0,1], [-1,-1], [1,0], [-1,0], [0,-1]];
    max_cones = [[1,4], [4,2], [2,5], [5,3], [3,1]];
    pentagon = normal_toric_variety(max_cones, ray_generators, non_redundant = true)
    
    return pentagon
end;
    
#---------------------------------------
# Defines the hexagon
#---------------------------------------
function define_hexagon()
    ray_generators = [[0,1], [-1,-1], [1,0], [-1,0], [0,-1], [1,1]];
    max_cones = [[1,4], [4,2], [2,5], [5,3], [3,6], [6,1]];
    hexagon = normal_toric_variety(max_cones, ray_generators, non_redundant = true)
    
    return hexagon
end;


########################################
# DIVISORS
########################################
# Basic properties of divisors on toric varieties
########################################

#---------------------------------------
# representation of the divisor class of a toric divisor
#---------------------------------------
function get_representative_of_divisor(td::ToricDivisor)   
    return divisor_class(toric_divisor_class(td)) 
end;

#---------------------------------------
# shows the generators and the relations of the prime divisors in the class-group
# returns the map pi: Div_T -> Cl as matrix
#---------------------------------------
function show_generators_and_relations_of_classgroup(v::NormalToricVariety; print_output::Bool=true)

    fan_v = polyhedral_fan(v)
    nrays_v = n_rays(fan_v)        
    cl_group = class_group(v)
    cl_rank = torsion_free_rank(cl_group)
    all_reps = []
    
    generators = [([i == j ? 1 : 0 for j in 1:cl_rank]) for i in 1:cl_rank]
    
    if print_output
        println("The class group is ", cl_group)  
        println("The Picard group is ", picard_group(v)) 
        println() 
    end
    
    # identify the index of the base elements
    base = []
    
    for c in 1:nrays_v
        vec = zeros(Int, nrays_v)
        vec[c] = 1
        div_c = toric_divisor(v, vec)
        cl_c = get_representative_of_divisor(div_c)
        rep_c = Int64[]
              
        for j in 1:cl_rank
            push!(rep_c, cl_c[j])
        end
        
        if rep_c in generators
            push!(base, c)
            deleteat!(generators, findall(x->x==rep_c,generators)) 
        end
        
        push!(all_reps,rep_c)
    end
    
    # create the output strings
    for c in 1:nrays_v
        if c in base
            print_output && println("Divisorclass of D$c is: ", all_reps[c], " => generator of class-group")     
        else
            str = ""
            for k in 1:cl_rank
                fac = all_reps[c][k]
                
                if fac == 0
                    #nothing
                elseif fac == 1
                    str = str * "D$(base[k]) + "
                elseif fac == -1
                    str = str * "-D$(base[k]) + "
                else
                    str = str * "$fac*D$(base[k]) + "
                end
                
                if k == cl_rank
                    str = chop(str, tail = 2)
                end
            end
            
            print_output && println("Divisorclass of D$c is: ", all_reps[c], " => D$c ~ ", str)
        end    
    end
    
    mat_pi = hcat(all_reps...)
    
    if print_output
        println()
        println("The map pi: Div_T(X) --> Cl(X) is given by the output matrix")
    end
     
    return mat_pi
    
end;

#---------------------------------------
# Calculates the Cartier-Data of a given Cartier-Divisor
#---------------------------------------
function calculate_cartier_data(td::ToricDivisor, cones::Bool=true)

    # implement a check if divisor is not cartier
    coeffs_td_int = transform_rayvector(coefficients(td))  
    var = toric_variety(td)
    rays_var = rays(var)
    maxcones_var = maximal_cones(polyhedral_fan(var))
    
    rays_var_int = []
    cartier_data = []
    
    for r in rays_var
        push!(rays_var_int, transform_rayvector(r))    
    end
        
    rays_matrix = transpose(hcat(rays_var_int...))    
    
    for c in maxcones_var
        rays_c = rays(c)        
        ind = Int64[]
        
        for r in rays_c
            i = findfirst(item -> item == r, rays_var)
            push!(ind,i)
        end
        
        A = rays_matrix[ind,:]
        b = -coeffs_td_int[ind]
        sol = A\b
        
        if cones == true
            out = c
        else
            out = rays(c)
        end
        
        push!(cartier_data,(out, sol))
        
    end
    
    return cartier_data
end;

#---------------------------------------
# Checks if the input function characterize the nef divisors with coefficients in range
#---------------------------------------
function check_nef_conditions(
        var::NormalToricVariety, 
        range::UnitRange{Int64}, 
        condition_func; 
        elements::Vector{Int64}=ones(Int, n_rays(var))
    )

    check = true
    
    all_coeffs = generate_vectors(n_rays(polyhedral_fan(var)), range)
    
    if elements != ones(Int, n_rays(var))
        all_coeffs = [v .* elements for v in all_coeffs]
        all_coeffs = collect(Set(all_coeffs)) 
    end
    
    for coeff in all_coeffs
        td = toric_divisor(var, coeff)
        cond = condition_func(coeff...)
        if is_nef(td) != cond
            check = false
            break
        end
    end

    return check
end;

#---------------------------------------
# Calculates a decomposition of a divisor into two nef divisors
# Caution when using the input basepoint_choice, an index other than 1 may not exist.
#---------------------------------------
function calculate_nef_pair(
        tv::NormalToricVariety, 
        td::ToricDivisor; 
        show_coeffs::Bool=true, 
        basepoint_choice::Int=1
    )

    n = torsion_free_rank(class_group(tv))     
    coeffs = Int64[]
    for i in 1:n
        index_i = divisor_class(toric_divisor_class(td))[i]
        push!(coeffs, index_i)
    end
    
    nefFacets = []
    for f in facets(nef_cone(tv)) 
        push!(nefFacets, transpose(transform_rayvector(f.a))) 
    end
    nefFacets = vcat(nefFacets...)
    r = nrows(nefFacets)
    c = ncols(nefFacets)
    
    B = hcat(I(n), -I(n)) 
    equalities = (B, coeffs) 
    
    A = vcat(nefFacets, zeros(Int64,r,c))
    AA = vcat(zeros(Int64,r,c), nefFacets)
    AAA = hcat(A, AA)
    inequalitites = (AAA, zeros(Int64,2*r)) 
    
    P = polyhedron(inequalitites, equalities)
    # we make a choice here!
    pvv = minimal_faces(P).base_points[basepoint_choice]
    pv1 = transform_rayvector(pvv[1:n])
    pv2 = transform_rayvector(pvv[n+1:2*n])
    
    td1 = toric_divisor(toric_divisor_class(tv, pv1))
    td2 = toric_divisor(toric_divisor_class(tv, pv2)) 
    
    ctd1 = transform_rayvector(coefficients(td1)) 
    ctd2 = transform_rayvector(coefficients(td2)) 
    ctd = transform_rayvector(coefficients(toric_divisor(toric_divisor_class(td)))) 
    println("The input divisor is linear equivalent to: ", ctd)
    
    #plausibility check
    if !is_nef(td1) || !is_nef(td2)
        println("td1 is nef: ", is_nef(td1)) 
        println("td2 is nef: ", is_nef(td2))
        throw(ErrorException("At least one calculated divisor is not nef!")) 
    end 
    
    if ctd1-ctd2 != ctd
        println(ctd1, " - ", ctd2, " = ", ctd1-ctd2, " != ", ctd) 
        throw(ErrorException("The decomposition is wrong!")) 
    end
    
    if show_coeffs == true
        return (ctd1, ctd2) 
    else
        return (td1, td2)
    end 
end;


########################################
# BLOW-UPS
########################################

#---------------------------------------
# Gives the exceptional divisor of a blowup as Vector{Int64}
#---------------------------------------
function get_exceptional_divisor(morphism::Oscar.ToricBlowupMorphism{NormalToricVariety, NormalToricVariety})
    exc_div = exceptional_prime_divisor(morphism)
    blow_up = domain(morphism)
    cl_rank = torsion_free_rank(class_group(blow_up))
    E = transform_rayvector([get_representative_of_divisor(exc_div)[i] for i in 1:cl_rank])
    return E
end;