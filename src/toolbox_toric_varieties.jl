########################################
# DEFINITION OF CERTAIN TORIC VARIETIES
########################################

#---------------------------------------
# Construction of important toric varieties such that the labeling of the rays is fixed
# The optional parameter non_redundant = true is mandatory for this purpose
#---------------------------------------


@doc raw"""
    define_projective_space(n::Int)

Construct the n-dimensional projective space with fixed labeling of the rays.

# Examples
```jldoctest
julia> P2 = define_projective_space(2)
Normal toric variety

```
"""
function define_projective_space(n::Int)
    # Ray generators
    ray_generators = Vector{Vector{Int}}()
    # e_0
    push!(ray_generators, fill(-1, n))
    
    # e_1,...,e_n
    append!(ray_generators, [zeros(Int, n) for _ in 1:n])
    for i in 1:n
        ray_generators[i+1][i] = 1
    end
    
    #Maximal cones
    max_cones = [filter(j -> j != i, 1:n+1) for i in 1:n+1]
    
    return normal_toric_variety(max_cones, ray_generators, non_redundant = true)
end


@doc raw"""
    define_projective_surface()

Construct the projective surface with fixed labeling of the rays.
The labeling of the rays differs from define_projective_space(n::Int)

# Examples
```jldoctest
julia> P2 = define_projective_surface()
Normal toric variety

```
"""
function define_projective_surface()
    ray_generators = [[0,1], [-1,-1], [1,0]]
    max_cones = [[1,2], [2,3], [3,1]]
    return normal_toric_variety(max_cones, ray_generators, non_redundant = true)
end


@doc raw"""
    define_hirzebruch_surface(r::Int)

Construct the r-th Hirzebruch surface with fixed labeling of the rays.

# Examples
```jldoctest
julia> H2 = define_hirzebruch_surface(2)
Normal toric variety

```
"""
function define_hirzebruch_surface(r::Int)
    ray_generators = [[0,1], [-r,-1], [1,0], [-1,0]]
    max_cones = [[1,4], [4,2], [2,3], [3,1]]
    return normal_toric_variety(max_cones, ray_generators, non_redundant = true)
end


@doc raw"""
    define_pp()

Construct (P^1)^2 with fixed labeling of the rays.

# Examples
```jldoctest
julia> pp = define_pp()
Normal toric variety

```
"""
function define_pp()
    ray_generators = [[1,0], [-1,0], [0,1], [0,-1]]
    max_cones = [[1,3], [1,4], [2,4], [2,3]]
    return normal_toric_variety(max_cones, ray_generators, non_redundant = true)
end


@doc raw"""
    define_ppp()

Construct (P^1)^3 with fixed labeling of the rays.

# Examples
```jldoctest
julia> ppp = define_ppp()
Normal toric variety

```
"""
function define_ppp()
    ray_generators = [[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]]
    max_cones = [[1,3,5], [1,3,6], [1,4,5], [1,4,6], [2,3,5], [2,3,6], [2,4,5], [2,4,6]]
    return normal_toric_variety(max_cones, ray_generators, non_redundant = true)
end


@doc raw"""
    define_pentagon()

Construct the pentagon with fixed labeling of the rays.

# Examples
```jldoctest
julia> pentagon = define_pentagon()
Normal toric variety

```
"""
function define_pentagon()
    ray_generators = [[0,1], [-1,-1], [1,0], [-1,0], [0,-1]]
    max_cones = [[1,4], [4,2], [2,5], [5,3], [3,1]]
    return normal_toric_variety(max_cones, ray_generators, non_redundant = true)
end
    

@doc raw"""
    define_hexagon()

Construct the pentagon with fixed labeling of the rays.

# Examples
```jldoctest
julia> hexagon = define_hexagon()
Normal toric variety


```
"""
function define_hexagon()
    ray_generators = [[0,1], [-1,-1], [1,0], [-1,0], [0,-1], [1,1]]
    max_cones = [[1,4], [4,2], [2,5], [5,3], [3,6], [6,1]]
    return normal_toric_variety(max_cones, ray_generators, non_redundant = true)
end


########################################
# DIVISORS
########################################

#---------------------------------------
# Calculate properties of toric divisors and their classes
#---------------------------------------


@doc raw"""
    get_representative_of_divisor(td::ToricDivisor)

Return the divisor class of a toric divisor in Cl = Z^n

# Examples
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> td = toric_divisor(h1, [1,1,1,1]);

julia> get_representative_of_divisor(td)
Abelian group element [1, 2]

```
"""
function get_representative_of_divisor(td::ToricDivisor) 
    return divisor_class(toric_divisor_class(td)) 
end


function divisor_class_int(td::ToricDivisor)
    return transform_abelian_group_element(get_representative_of_divisor(td))
end


function compute_divisor_classes_and_base(variety)
    # 1. Pre-calculate fixed properties
    fan = polyhedral_fan(variety)
    nrays = n_rays(fan)
    cl = class_group(variety)
    cl_rank = torsion_free_rank(cl)

    # 2. Initialize containers
    # Pre-allocate the vector of all divisor classes for efficiency
    all_reps = Vector{Int}[]
    sizehint!(all_reps, nrays) # Optional: hint the expected size

    # Use a mutable set of generator representations for fast lookups and deletion
    # Ensure the generator creation is correct for the vector type from divisor_class_int
    generators_set = Set{Vector{Int}}()
    for i in 1:cl_rank
        # Create the unit vector: [0, ..., 1, ..., 0]
        gen = zeros(Int, cl_rank)
        gen[i] = 1
        push!(generators_set, gen)
    end
    
    base = Int[] # The indices of the generators found

    # 3. Iterate over the rays efficiently
    # Use a pre-allocated vector that is modified in place inside the loop (optional but fast)
    ray_vector = zeros(Int, nrays) 
    
    for c in 1:nrays
        # Reset the previous 1 and set the new one (if c > 1)
        if c > 1
            ray_vector[c-1] = 0 
        end
        ray_vector[c] = 1
        
        # Calculate the divisor class
        div_c = toric_divisor(variety, ray_vector)
        cl_c = divisor_class_int(div_c)

        # Use the Set for quick check and removal
        if cl_c in generators_set
            push!(base, c)
            delete!(generators_set, cl_c) # O(1) average time complexity for Set deletion
        end
        
        push!(all_reps, cl_c)
    end
    
    # Optional: Check if all generators were found
    if !isempty(generators_set)
        @warn "Only found $(cl_rank - length(generators_set)) out of $(cl_rank) class group generators."
    end

    return base, all_reps
end


@doc raw"""
    show_generators_and_relations_of_classgroup(v::NormalToricVariety; print_output::Bool=true)

Gives the map pi: Div_T(X) -> Cl(X).
Prints generators and relations of the prime divisors in the class group.

# Arguments
- `v::NormalToricVariety`: The toric variety of interest
- `print_output::Bool=true`: (optional) Determine if generators and relations are printed

# Examples
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> show_generators_and_relations_of_classgroup(h1)
The class group is Z^2
The Picard group is Z^2

Divisorclass of D1 is: [1, 0] => generator of class-group
Divisorclass of D2 is: [1, 0] => D2 ~ D1
Divisorclass of D3 is: [0, 1] => generator of class-group
Divisorclass of D4 is: [-1, 1] => D4 ~ -D1 + D3

The map pi: Div_T(X) --> Cl(X) is given by the matrix
1 1 0 -1
0 0 1 1

The anticanonical divisor class is [1, 2]
The toric variety is Fano: true
2×4 Matrix{Int64}:
 1  1  0  -1
 0  0  1   1

```
"""
function show_generators_and_relations_of_classgroup(v::NormalToricVariety; print_output::Bool=true)
    fan_v = polyhedral_fan(v)
    nrays_v = n_rays(fan_v)        
    cl_group = class_group(v)
    cl_rank = torsion_free_rank(cl_group)

    # print info
    if print_output
        println("The class group is ", cl_group)  
        println("The Picard group is ", picard_group(v)) 
        println() 
    end

    # explicit typed containers
    all_reps = Vector{Vector{Int}}()
    base = Int[]
    
    # generators are unit vectors
    generators = [([i == j ? 1 : 0 for j in 1:cl_rank]) for i in 1:cl_rank]
    
    
    for c in 1:nrays_v
        vec = zeros(Int, nrays_v)
        vec[c] = 1
        div_c = toric_divisor(v, vec)
        cl_c = get_representative_of_divisor(div_c)

        rep_c = Int64[]            
        for j in 1:cl_rank
            push!(rep_c, cl_c[j])
        end
        
        # check if it's a generator
        idx = findfirst(==(rep_c), generators)
        if idx !== nothing
            push!(base, c)
            deleteat!(generators, idx)
        end
        
        push!(all_reps,rep_c)
    end
    
    # create the output strings
    for c in 1:nrays_v
        if c in base
            print_output && println("Divisorclass of D$c is: ", all_reps[c], " => generator of class-group")     
        else
            buf = IOBuffer()
            for k in 1:cl_rank
                fac = all_reps[c][k]
                if fac == 0
                    continue
                elseif fac == 1
                    print(buf, "D$(base[k]) + ")
                elseif fac == -1
                    print(buf, "-D$(base[k]) + ")
                else
                    print(buf, "$fac*D$(base[k]) + ")
                end
            end
            str = String(take!(buf))
            str = isempty(str) ? "" : chop(str, tail=2)
            print_output && println("Divisorclass of D$c is: ", all_reps[c], " => D$c ~ ", str)
        end
    end
    
    # matrix form of the map
    mat_pi = hcat(all_reps...)
    
    if print_output
        println()
        println("The map pi: Div_T(X) --> Cl(X) is given by the matrix")
        for row in eachrow(mat_pi)
            println(join(row, " "))
        end
    end
    
    if print_output
        anti_coefs = ones(Int64, nrays_v)
        anti = mat_pi * anti_coefs
        ample = is_ample(toric_divisor(v, anti_coefs))
        println()
        println("The anticanonical divisor class is $anti")
        println("The toric variety is Fano: $ample")
    end
     
    return mat_pi
end


@doc raw"""
    calculate_cartier_data(td::ToricDivisor, cones::Bool=true)

Calculate the Cartier data of a Cartier divisor together with the corresponding cones.

# Arguments
- `td::ToricDivisor`: The toric divisor of interest
- `cones::Bool=true`: (optional) Determines if the corresponding cones or their rays are returned

# Examples
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> td = toric_divisor(h1, [1,1,1,1]);

julia> calculate_cartier_data(td)
4-element Vector{Tuple{Any, Vector{Float64}}}:
 (Polyhedral cone in ambient dimension 2, [1.0, -1.0])
 ⋮
 (Polyhedral cone in ambient dimension 2, [-1.0, -1.0])

```
"""
function calculate_cartier_data(td::ToricDivisor, cones::Bool=true)
    coeffs_td_int = transform_rayvector(coefficients(td))  
    var = toric_variety(td)
    rays_var = rays(var)
    maxcones_var = maximal_cones(polyhedral_fan(var))
    
    rays_var_int  = [transform_rayvector(r) for r in rays_var]
    cartier_data = Vector{Tuple{Any, Vector{Float64}}}()
    
    # compute ray matrix
    rays_matrix = transpose(hcat(rays_var_int...))    
    
    for c in maxcones_var
        rays_c = rays(c)
        
        # find indices of rays of the current cone
        ind = [findfirst(==(r), rays_var) for r in rays_c]
        
        A = rays_matrix[ind, :]
        b = -coeffs_td_int[ind]
        sol = A\b
        
        out = cones ? c : rays_c
        push!(cartier_data, (out, sol))
    end
    
    return cartier_data
end


@doc raw"""
    check_nef_conditions(
        var::NormalToricVariety, 
        range::UnitRange{Int64}, 
        condition_func; 
        elements::Vector{Int64}=ones(Int, n_rays(var))
    )

Check if the input function characterizes all nef divisors with coefficients in a range correctly.

# Arguments
- `var::NormalToricVariety`: The toric variety of interest
- `range::UnitRange{Int64}`: The coefficients of the divisors in the class group
- `condition_func`: Function with certain conditions on the coefficients of a toric divisor
- `elements::Vector{Int64}=ones(Int, n_rays(var))`: (optional) Applies a permutation to the coefficients (if another basis is required)

# Examples
```jldoctest
julia> p2 = define_projective_surface();

julia> coeff_range = -2:2;

julia> cond = (a, b, c) -> a + b + c >= 0;

julia> check_nef_conditions(p2, coeff_range, cond)
true

```
"""
function check_nef_conditions(
        var::NormalToricVariety, 
        range::UnitRange{Int64}, 
        condition_func; 
        elements::Union{Nothing, Vector{Int}}=nothing
    )

    nrays = n_rays(polyhedral_fan(var))
    all_coeffs = generate_vectors(nrays, range)
    
    if elements !== nothing
        all_coeffs = [v .* elements for v in all_coeffs]
        unique!(all_coeffs)  # in-place deduplication
    end
    
    for coeff in all_coeffs
        td = toric_divisor(var, coeff)
        cond = condition_func(coeff...)
        if is_nef(td) != cond
            return false
        end
    end

    return true
end


@doc raw"""
    calculate_nef_pair(
        tv::NormalToricVariety, 
        td::ToricDivisor; 
        show_coeffs::Bool=true, 
        basepoint_choice::Int=1
    )

Calculate a nef pair of a Q-Cartier divisor on a semiprojective toric variety.
This decomposition is in general not unique.

# Arguments
- `tv::NormalToricVariety`: The semiprojective toric variety
- `td::ToricDivisor`: The toric divisor for which a decomposition into nef divisors should be found
- `show_coeffs::Bool=true`: (optional) Determines if the coefficients of the divisors or the divisors are returned
- `basepoint_choice::Int=1`: (optional) Pick one of the possible decompositions

# Examples
```jldoctest
julia> hexagon = define_hexagon();

julia> td = toric_divisor(hexagon, [-4,-2,-1,-3,0,0]);

julia> calculate_nef_pair(hexagon, td)
The input divisor is linear equivalent to: [-4, -2, -1, -3, 0, 0]
([1, 1, -1, 2, 0, 0], [5, 3, 0, 5, 0, 0])

julia> calculate_nef_pair(hexagon, td; basepoint_choice=2)
The input divisor is linear equivalent to: [-4, -2, -1, -3, 0, 0]
([1, 2, -1, 2, 0, 0], [5, 4, 0, 5, 0, 0])

```
"""
function calculate_nef_pair(
        tv::NormalToricVariety, 
        td::ToricDivisor; 
        show_coeffs::Bool=true, 
        basepoint_choice::Int=1
    )

    n = torsion_free_rank(class_group(tv))

    # coefficients of the divisor class
    coeffs = Int64[]
    for i in 1:n
        index_i = divisor_class(toric_divisor_class(td))[i]
        push!(coeffs, index_i)
    end

    # collect facet inequalities
    nefFacets = vcat([transpose(transform_rayvector(f.a)) for f in facets(nef_cone(tv))]...)
    r, c = size(nefFacets)
    
    # equalities
    B = hcat(I(n), -I(n)) 
    equalities = (B, coeffs) 
    
    # inequalities
    A = vcat(nefFacets, zeros(Int64,r,c))
    AA = vcat(zeros(Int64,r,c), nefFacets)
    AAA = hcat(A, AA)
    inequalitites = (AAA, zeros(Int64,2*r)) 
    
    P = polyhedron(inequalitites, equalities)

    # choose basepoint
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
    
    return show_coeffs ? (ctd1, ctd2) : (td1, td2)
end


@doc raw"""
    get_exceptional_divisor(morphism::Oscar.ToricBlowupMorphism{NormalToricVariety, NormalToricVariety})

Return the exceptional divisor of a toric blowup morphism as a vector of integers.

# Arguments
- `morphism::Oscar.ToricBlowupMorphism{NormalToricVariety, NormalToricVariety}`: The toric blow up morphism

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> new_ray = [-1, -1];

julia> blowup_morph = blow_up(pentagon, new_ray);

julia> get_exceptional_divisor(blowup_morph)
3-element Vector{Int64}:
 0
 1
 0


```
"""
function get_exceptional_divisor(morphism::Oscar.ToricBlowupMorphism{NormalToricVariety, NormalToricVariety})
    exc_div = exceptional_prime_divisor(morphism)
    blow_up = domain(morphism)
    cl_rank = torsion_free_rank(class_group(blow_up))
    E = transform_rayvector([get_representative_of_divisor(exc_div)[i] for i in 1:cl_rank])
    return E
end


########################################
# SECONDARY FAN
########################################

function is_query_block_contained(
    query::Vector{Vector{Vector{T}}},
    origin::Vector{Vector{Vector{Vector{T}}}}
) where T

    # Convert matrices to hashable types (tuple-of-tuples)
    hashable_query = [Tuple.(mat) for mat in query]
    query_counts = countmap(hashable_query)

    for block in origin
        hashable_block = [Tuple.(mat) for mat in block]
        block_counts = countmap(hashable_block)

        # Check if all query matrices are contained in this block
        if all(p -> get(block_counts, p.first, 0) ≥ p.second, pairs(query_counts))
            return true  # found one block that contains the query
        end
    end

    return false  # no block contained the query
end


function calculate_secondary_fan(beta::Vector{Vector{T}}) where T
    #STEP 1 - Decide for a finite set of vectors in Z^d if they span R^d
    d = length(beta[1])
    combs = [collect(c) for c in Combinatorics.combinations(beta, d)] #find all subsets of right size
    independent_sets = [comb for comb in combs if rank(hcat(comb...)) == d] #check if they are linear independant    
    max_k = length(independent_sets)   
    
    # STEP 2 - Compute the intersection of rational cones
    maximal_sets = Vector{Vector{Vector{Vector{T}}}}()

    for k in max_k:-1:1
        all_subsets_k = collect(Combinatorics.combinations(independent_sets, k))      

        for s in all_subsets_k
            # Skip if subset of already accepted, in this case, A is not maximal!
            if !is_empty(maximal_sets)
                if is_query_block_contained(s, maximal_sets)
                    continue
                end
            end
            
            matrix = [vcat((x' for x in group)...) for group in s]
            cone = positive_hull(matrix[1])
            
            for m in matrix[2:end]
                mcone = positive_hull(m)
                cone = intersect(cone, mcone)
            end
            
            if dim(cone) == d
                push!(maximal_sets, s)
            end
        end
    end
    
    #STEP 3 - Compute the faces of a rational cone
    tau = Vector{Cone{QQFieldElem}}()
    
    for A in maximal_sets
        matrix = [vcat((x' for x in group)...) for group in A]
        cone = positive_hull(matrix[1])
            
        for m in matrix[2:end]
            mcone = positive_hull(m)
            cone = intersect(cone, mcone)
        end
        
        pf = polyhedral_fan(cone)
        for i in 0:dim(pf)
            pf_cones = cones(pf, i)
            append!(tau, pf_cones)
        end        
    end
    
    secondary_fan = polyhedral_fan(tau)
    return secondary_fan 
end


function find_adjacent_cones(cones::Vector{Vector{Int}})
    n = length(cones)
    sets = [Set(c) for c in cones]

    adj_pairs = Vector{Tuple{Tuple{Int,Int}, Vector{Int}}}()
    for i in 1:n, j in i+1:n
        common = intersect(sets[i], sets[j])
        if length(common) == 3
            push!(adj_pairs, ((i, j), sort(collect(common))))  # tuple of indices + vector of common elements
        end
    end
    return adj_pairs
end


function oriented_cocircuit(spanvecs::Vector{Vector{Int}}, vecs::Vector{Vector{Int}})
    d = length(spanvecs[1])
    M = [row[i] for row in spanvecs, i in 1:length(spanvecs[1])]

    if rank(M) != d-1
        error("Expected $(d-1) spanning vectors for a hyperplane in ℝ^$d")
    end

    # Calculate exact kernel
    k = kernel(matrix(ZZ, M'))

    # Kernel is generated by one element, transform ZZMatrix to Int-vector
    gen = Int.(k[1,:])

    # Canonical orientation: first nonzero entry must be positive
    for x in gen
        if x < 0
            gen = [-yi for yi in gen]
            break
        elseif x > 0
            break
        end
    end

    cocircuit = [
        sum(gen[i]*v[i] for i in 1:d) > 0 ? 1 :
        sum(gen[i]*v[i] for i in 1:d) < 0 ? -1 : 0 
        for v in vecs
    ]

    return cocircuit
end


function classify_vector(v::Vector{Int})
    npos = count(==(1), v)   # number of 1s
    nneg = count(==(-1), v)  # number of -1s

    # error for invalid or trivial cases
    if npos == 0 || nneg == 0
        error("Invalid vector: must contain at least one 1 and one -1")
    end

    if npos == 1 && nneg == 1
        return "degenerate"
    elseif npos == 1 || nneg == 1
        return "divisorial"
    elseif npos > 1 && nneg > 1
        return "flipping"
    else
        # This should not happen if input is valid
        error("Unexpected case in classify_vector")
    end
end


function extremal_position(beta::Vector{Vector{Int64}})
    position = Int64[]
    l = length(beta)
    for i in 1:l
        beta_i = beta[i]
        beta_hat = [beta[j] for j in 1:l if j != i]
        cone_beta_hat = positive_hull(beta_hat)
        contained = in(beta_i, cone_beta_hat)
        if contained == false
            push!(position, i)
        end
    end      
    return position
end


