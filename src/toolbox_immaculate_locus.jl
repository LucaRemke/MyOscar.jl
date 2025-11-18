########################################
# IMMACULATE LINE BUNDLES & COHOMOLOGY
########################################

# This is the improved version of get_immaculate_lb:
function ilb_range(var::NormalToricVariety, range::UnitRange{Int})
    ilb = immaculate_line_bundles(var)
    
    rk = rank(picard_group(var))
    classes = generate_vectors(rk, range)
    
    ilb_in_range = Vector{Int}[]
    
    for c in classes
        lb = toric_line_bundle(var, c)
        if contains(ilb, lb)
            push!(ilb_in_range, c)
        end
    end
    
    return ilb_in_range
end 


@doc raw"""
    get_immaculate_lb(
        variety::NormalToricVariety, 
        range::UnitRange{<:Integer};
        perm::Union{Nothing,AbstractVector{<:Integer}}=nothing, 
        output_td::Bool=false, 
        status::Bool=false
    )

Check for each divisor with coefficients in a specific range if it has vanishing cohomologies.
Return the immaculate line bundles.

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest
- `range::UnitRange{<:Integer}`: (optional) Determines if the corresponding cones or their rays are returned
- `perm::Union{Nothing,AbstractVector{<:Integer}}`: (optional) Permutation on the vectors
- `output_td::Bool=false`: (optional) Determines if linebundles or their coefficients are returned
- `status::Bool=falses`: (optional) Prints the current linebundle coefficients

# Examples
```jldoctest
julia> h0 = define_hirzebruch_surface(0);

julia> get_immaculate_lb(h0, -2:2)
4-element Vector{Vector{Int64}}:
 [-2, 1]
 ⋮
 [1, -2]


```
"""
function get_immaculate_lb(
        variety::NormalToricVariety, 
        range::UnitRange{<:Integer};
        perm::Union{Nothing,AbstractVector{<:Integer}}=nothing, 
        output_td::Bool=false, 
        status::Bool=false
    )

    nrays = n_rays(variety)
    n_pic = torsion_free_rank(picard_group(variety)) 
    n_zeros = nrays - n_pic

    immaculate = output_td ? ToricLineBundle[] : Vector{Int}[]

    coeffs = generate_vectors(n_pic, range) 
        
    for coeff in coeffs
        coeff_complete = vcat(coeff, fill(0, n_zeros)) 
        
        if !isnothing(perm)
            coeff_complete = coeff_complete[perm]
        end
        
        td = toric_divisor(variety, coeff_complete)
        lb = toric_line_bundle(td)
        cohom = all_cohomologies(lb)
        
        status && println(coeff_complete)
        
        if all(==(0), cohom)
            if output_td
                push!(immaculate, lb)
            else
                push!(immaculate, coeff)
            end
        end
    end

    return immaculate
end


########################################
# TEMPTING SUBSETS
########################################


@doc raw"""
    calculate_tempting_subsets(var::NormalToricVariety; zero_ray::Bool=false)

Return the tempting subsets of a fan.
For that the function iterates over all subsets of Sigma(1) and checks the cohomologies of -sum(D_rho) for rho in the current subset.

# Arguments
- `var::NormalToricVariety`: The toric divisor of interest
- `zero_ray::Bool`: (optional) Determines if the ray labeling starts with 0 or 1

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> calculate_tempting_subsets(pentagon)
12-element Vector{Vector{Int64}}:
 []
 ⋮
 [1, 2, 3, 4, 5]


```
"""
function calculate_tempting_subsets(var::NormalToricVariety; zero_ray::Bool=false)
    v_fan = polyhedral_fan(var)
    nrays = n_rays(v_fan)
    temptings = Vector{Int}[]
    
    for subset in powerset(1:nrays)
        coeffs = zeros(Int, nrays)
        coeffs[subset] .= -1
        divisor = toric_divisor(var, coeffs)
        lb = toric_line_bundle(divisor)
        cohom = all_cohomologies(lb)
        
        if any(!=(0), cohom)
            push!(temptings, subset)
        end
    end
    
    if zero_ray
        temptings = [vec .- 1 for vec in temptings]
    end
    
    return temptings   
end


@doc raw"""
    search_nontemptings_via_boundary_condition(
        variety::NormalToricVariety; 
        output_indizes::Bool=true, 
        output_nontemptings::Bool=true
    )

Eliminate the non temptings via the boundary condition
This approach works if the toric variety is smooth since rho is surjective in this case

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest
- `output_indizes::Bool`: (optional) Determines if the rays or their indizes are returned
- `output_nontemptings::Bool`: (optional) Determines if the temptings or the non-temptings are returned

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> search_nontemptings_via_boundary_condition(pentagon)
14-element Vector{Vector{Int64}}:
 [1]
 ⋮
 [2, 3, 4, 5]


```
"""
function search_nontemptings_via_boundary_condition(
        variety::NormalToricVariety; 
        output_indizes::Bool=true, 
        output_nontemptings::Bool=true
    )

    nrays_variety = n_rays(variety)
    ray_matrix = Matrix{Int}(matrix(ZZ, rays(variety)))
    ray_vector = [ray_matrix[row, :] for row in 1:size(ray_matrix, 1)]

    zeros_im = zeros(Int, nrays_variety)
    zeros_pim = zeros(Int, 1, dim(variety))

    output = Vector{Int}[]
            
    for subset in powerset(1:nrays_variety)  
        B = -ray_matrix     
        B[subset, :] .= B[subset, :] .* -1  
        P = polyhedron(B, zeros_im)
        
        is_nontempting = (P != convex_hull(zeros_pim))

        if output_nontemptings == is_nontempting
            push!(output, output_indizes ? subset : ray_vector[subset])
        end   
    end    

    return output
end


@doc raw"""
    eliminate_faces_as_temptings(
        variety::NormalToricVariety; 
        output_indizes::Bool=true, 
        output_nontemptings::Bool=true
    )

Return the non-temptings via faces of the fan

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest
- `output_indizes::Bool`: (optional) Determines if the rays or their indizes are returned
- `output_nontemptings::Bool`: (optional) Determines if the temptings or the non-temptings are returned

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> eliminate_faces_as_temptings(pentagon)
12-element Vector{Vector{Int64}}:
 []
 ⋮
 [1, 2, 3, 4, 5]


```
"""
function eliminate_faces_as_temptings(
        variety::NormalToricVariety; 
        output_indizes::Bool=true, 
        output_nontemptings::Bool=true
    )

    fan_variety = polyhedral_fan(variety)   
    nrays_variety = n_rays(variety)
    pset_nrays = collect(powerset(1:nrays_variety))

    ray_matrix = Matrix{Int}(matrix(ZZ, rays(variety)))
    ray_vector = [ray_matrix[row, :] for row in 1:size(ray_matrix, 1)]

    face_subsets = Vector{Int}[]

    for i in 1:dim(fan_variety)        
        for j in cones(fan_variety, i)
            rays_matrix_ij = Matrix{Int}(matrix(ZZ, rays(j)))
            rays_vector_ij = [rays_matrix_ij[row, :] for row in 1:size(rays_matrix_ij, 1)]
            
            idx_inside  = filter(!isnothing, indexin(rays_vector_ij, ray_vector))
            idx_outside = filter(!isnothing, indexin(setdiff(ray_vector, rays_vector_ij), ray_vector))

            push!(face_subsets, idx_inside)
            push!(face_subsets, idx_outside)
        end
    end

    remaining = setdiff(pset_nrays, face_subsets)

    if !output_nontemptings
        remaining = setdiff(pset_nrays, remaining)
    end

    if !output_indizes
        return [ray_vector[subset] for subset in remaining]
    else
        return remaining
    end
end


@doc raw"""
    get_temptings_as_primitive_collections(
        variety::NormalToricVariety;
        output_indizes::Bool=true, 
        output_nontemptings::Bool=true, 
        include_emptyset::Bool=true
    )

Return temptings by determining the primitive collections

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest
- `output_indizes::Bool`: (optional) Determines if the rays or their indizes are returned
- `output_nontemptings::Bool`: (optional) Determines if the temptings or the non-temptings are returned
- `include_emptyset::Bool`: (optional) Determines if the emptyset should be added to the output

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> get_temptings_as_primitive_collections(pentagon)
20-element Vector{Vector{Int64}}:
 [1]
 ⋮
 [2, 3, 4, 5]


```
"""
function get_temptings_as_primitive_collections(
        variety::NormalToricVariety;
        output_indizes::Bool=true, 
        output_nontemptings::Bool=true, 
        include_emptyset::Bool=true
    )

    nrays_variety = n_rays(variety)
    set_nrays = Set(1:nrays_variety)
    all_sets = [Set(v) for v in collect(powerset(1:nrays_variety))] 

    ray_matrix = Matrix{Int}(matrix(ZZ, rays(variety)))
    ray_vector = [ray_matrix[row, :] for row in 1:size(ray_matrix, 1)]

    output = Set{Set{Int}}()
    for pc in primitive_collections(variety)
        push!(output, pc)
        push!(output, setdiff(set_nrays, pc)) 
    end

    if include_emptyset
        push!(output, Set()) 
        push!(output, set_nrays)
    end

    if output_nontemptings
        output = setdiff(Set(all_sets), output)
    end

    # convert to sorted vectors, sort first by length then lexicographically
    out_vecs = sort([sort(collect(s)) for s in output], by = x -> (length(x), x))

    if !output_indizes
        out_vecs = [[ray_vector[j] for j in subset] for subset in out_vecs]
    end
        
    return out_vecs
end


@doc raw"""
    get_images_of_cube_vertices(variety::NormalToricVariety; print_output::Bool=true)

Calculate the images of the cube under pi

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest
- `print_output::Bool`: (optional) Determines if information is printed

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> get_images_of_cube_vertices(pentagon; print_output=false)
24-element Vector{Tuple{Vector{Int64}, Vector{Vector{Int64}}}}:
 ([0, 0, -2], [[0, -1, -1, -1, 0]])
 ⋮
 ([0, 0, -1], [[0, -1, 0, -1, 0], [0, 0, -1, 0, 0]])


```
"""
function get_images_of_cube_vertices(variety::NormalToricVariety; print_output::Bool=true)
    pi_v = show_generators_and_relations_of_classgroup(variety; print_output=false)
    nrays_v = n_rays(variety)
    cube_v = cube(nrays_v,-1,0) 
    vertices_v = vertices(cube_v)
    
    preimages = [transform_rayvector(v) for v in vertices_v]
    images    = [pi_v * p for p in preimages]
    
    # Group by image
    dict = Dict{Vector{Int}, Vector{Vector{Int}}}()
    for (im, pim) in zip(images, preimages)
        dict[im] = get(dict, im, Vector{Vector{Int}}())
        push!(dict[im], pim)
    end

    results = [(im, pims) for (im, pims) in dict]

    if print_output
        for (im, pims) in results
            im_str = "(" * join(im, ",") * ")"
            str = im_str
            for pim in pims
                pim_str = "(" * join((-1) .* pim, ",") * ")"
                str *= " = -pi$pim_str"
            end
            println(str)
        end
    end

    return results
end


@doc raw"""
    get_polytope_of_cube_vertices(variety::NormalToricVariety; print_output::Bool=true)

Calculate the immaculate line bundles comming from the cube

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest
- `print_output::Bool`: (optional) Determines if information is printed

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> get_polytope_of_cube_vertices(pentagon; print_output=true)
We have 32 subsets, where 12 are tempting and 20 are not.

The polytope P has 24 lattice points, of which 18 are vertices.
- 12 vertices come from the 12 tempting subsets.
- 6 vertices come from 6 non-tempting subsets.
The remaining 6 lattice points comes from the remaining 14 non-temptings.
- 2 are interior lattice points.
- 4 are boundary points (but no vertices).

In particular, the cube produces 12 immaculate line bundles:
[-1, 0, 0]
[0, -1, 0]
[0, 0, -1]
[0, 1, -1]
[-1, 1, 0]
[-1, 0, -1]
[-1, 1, -1]
[-1, 1, -2]
[-2, 1, -1]
[-1, 0, -2]
[-2, 0, -1]
[-2, 2, -2]
12-element Vector{Vector{Int64}}:
 [-1, 0, 0]
 ⋮
 [-2, 2, -2]


```
"""
function get_polytope_of_cube_vertices(variety::NormalToricVariety; print_output::Bool=true)
    pi_v = show_generators_and_relations_of_classgroup(variety; print_output=false)
    nrays_v = n_rays(variety)

    # Cube and its polytope image
    cube_v = cube(nrays_v,-1,0) 
    vertices_v = vertices(cube_v)   
    images = [pi_v * vec for vec in vertices_v]
    mat = transpose(hcat(images...)) 
    polytope = convex_hull(mat)
    
    # Determine the vertices corresponding to (non)tempting sets
    pset_nrays_v = collect(powerset(1:nrays_v)) 
    temptings = calculate_tempting_subsets(variety)
    
    # Images of tempting vertices
    temp_vertices = [begin
        v = zeros(Int, nrays_v); v[subset] .= -1; v
    end for subset in temptings]
    images_temp_vertices = unique([pi_v * vec for vec in temp_vertices])
    
     # Images of non-tempting vertices
    nontemptings = setdiff(pset_nrays_v, temptings)
    ntemp_vertices = [begin
        v = zeros(Int, nrays_v); v[subset] .= -1; v
    end for subset in nontemptings]
    images_ntemp_vertices = unique([pi_v * vec for vec in ntemp_vertices])
    
    ntemp_vertices = setdiff(transform_rayvector.(vertices(polytope)), images_temp_vertices)
    
    if print_output
        l_sub = length(pset_nrays_v)
        l_sub_t = length(temptings)
        l_im_t = length(images_temp_vertices)
        l_sub_nt = length(nontemptings)
        l_im_nt = length(images_ntemp_vertices)
        l_vert = n_vertices(polytope)
        l_lp = length(lattice_points(polytope)) 
        l_vert_nt = length(ntemp_vertices)
        l_lp_int = length(interior_lattice_points(polytope)) 
        l_lp_bound = length(boundary_lattice_points(polytope))

        println("We have $l_sub subsets, where $l_sub_t are tempting and $l_sub_nt are not.")
        println()
        println("The polytope P has $l_lp lattice points, of which $l_vert are vertices.")
        println("- $l_im_t vertices come from the $l_sub_t tempting subsets.")
        println("- $l_vert_nt vertices come from $l_vert_nt non-tempting subsets.")
        println("The remaining $(l_lp-l_vert) lattice points comes from the remaining $(l_sub_nt-l_vert_nt) non-temptings.")
        println("- $l_lp_int are interior lattice points.")
        println("- $(l_lp_bound-l_vert) are boundary points (but no vertices).")
        println()
        println("In particular, the cube produces $l_im_nt immaculate line bundles:")
        for i in images_ntemp_vertices
            println(i)
        end
    end
    
    return images_ntemp_vertices  
end


########################################
# MACULATE REGIONS
########################################


@doc raw"""
    get_maculate_regions(variety::NormalToricVariety)

Calculate maculate regions of a line bundle.
Return a vector of polyhedrons.

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest

# Examples
```jldoctest
julia> pentagon = define_pentagon();

julia> get_maculate_regions(pentagon)
12-element Vector{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 3
 ⋮
 Polyhedron in ambient dimension 3


```
"""
function get_maculate_regions(variety::NormalToricVariety)
    nrays_var = n_rays(variety)
    dim_pic = torsion_free_rank(picard_group(variety))
    pi = show_generators_and_relations_of_classgroup(variety; print_output=false)
    temptings = calculate_tempting_subsets(variety)
    
    regions = Polyhedron{QQFieldElem}[]
    
    for temp in temptings
        base_vec_t = [in(i, temp) ? -1 : 0 for i in 1:nrays_var]
        vertice_t = reshape(pi*base_vec_t,1,dim_pic)
        
        gens_t = Matrix{Int}(undef, 0, dim_pic)
        for i in 1:nrays_var
            v_i = (base_vec_t[i] == -1) ? 
                [j == i ? -1 : 0 for j in 1:nrays_var] : [j == i ? 1 : 0 for j in 1:nrays_var]
            gens_t = vcat(gens_t, transpose(pi*v_i)) 
        end
        
        polyhedron_t = convex_hull(vertice_t, gens_t)               
        push!(regions, polyhedron_t)
    end
    
    return regions
end


@doc raw"""
    cutout_maculate_regions(variety::NormalToricVariety, cutout::Polyhedron{QQFieldElem})

Calculate maculate regions, the H^i cones of these regions and the immaculate line bundles in a spefic area of coefficients.
If error "Denominator must be 1" occurs, enlarge the cutout.

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest
- `cutout::Polyhedron{QQFieldElem}`: Determines which part of the regions are shown

# Examples
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> cutout = convex_hull([-6 -6; -6 4; 4 -6; 4 4]);

julia> cutout_maculate_regions(h1, cutout);


```
"""
function cutout_maculate_regions(variety::NormalToricVariety, cutout::Polyhedron{QQFieldElem})
    lp_cutout = transform_rayvector.(lattice_points(cutout))
    temptations = calculate_tempting_subsets(variety)
    maculate_regions = get_maculate_regions(variety)

    lp_regions = Vector{Vector{DataFrame}}()
    lp_nonregions = lp_cutout
    
    # to visualize the regions we cutout a part of it
    # we obtain finitely many lattice point which can be plotted 
    for region in maculate_regions
        isec_region = intersect(cutout, region)
        if isempty(lattice_points(isec_region))
            continue
        end
        
        lp_region = transform_rayvector.(lattice_points(isec_region))
        lp_nonregions = [v1 for v1 in lp_nonregions if !(v1 in lp_region)]
        
        vertices_region = transform_rayvector.(vertices(isec_region))
        sorted_vertices_region = sort_points_clockwise(vertices_region)
        
        df_lp_region = convert_to_df(lp_region)
        df_sorted_vertices_region = convert_to_df(sorted_vertices_region)     
        push!(lp_regions, [df_lp_region, df_sorted_vertices_region])
     end
    
    # determine the H^i cone
    hcone = Int[]
    for temp in temptations
        dr_coeff = [in(i, temp) ? -1 : 0 for i in 1:n_rays(variety)]
        cohom = transform_rayvector(all_cohomologies(toric_line_bundle(toric_divisor(variety, dr_coeff)))) 
        i = findfirst(x -> x != 0, cohom) - 1
        push!(hcone, i)
    end
   
    # determine the immaculate linebundles in cutout
    df_lp_nonregions = convert_to_df(lp_nonregions)
    
    return (regions = lp_regions, hcone = hcone, immaculate = df_lp_nonregions)
end


@doc raw"""
    print_maculate_region_info(variety::NormalToricVariety)

Print a description of all maculate regions.

# Arguments
- `variety::NormalToricVariety`: The toric divisor of interest

# Examples
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> print_maculate_region_info(h1)
The maculate region of R = Int64[] has:
	 Basepoint [0, 0]
	 Generators [[1, 0], [-1, 1]]

The maculate region of R = [1, 2] has:
	 Basepoint [-2, 0]
	 Generators [[-1, 0], [0, 1]]

The maculate region of R = [3, 4] has:
	 Basepoint [1, -2]
	 Generators [[1, 0], [0, -1]]

The maculate region of R = [1, 2, 3, 4] has:
	 Basepoint [-1, -2]
	 Generators [[-1, 0], [1, -1]]


```
"""
function print_maculate_region_info(variety::NormalToricVariety)
    temptings = calculate_tempting_subsets(variety)
    regions = get_maculate_regions(variety)

    for i in eachindex(regions)
        region = regions[i]
        tempting = temptings[i]

        basepoint = vcat(transform_rayvector.(minimal_faces(region)[1])...)
        generators = transform_rayvector.(rays_modulo_lineality(region)[1])

        println("The maculate region of R = $tempting has:")
        println("\t Basepoint $basepoint")
        println("\t Generators $generators")
        i != length(regions) ? println("") : nothing
    end
end


########################################
# SPECIAL FUNCTIONS FOR TORIC VARIETIES OF PICARD RANK 3
########################################


@doc raw"""
    get_seed_and_chull_of_toricpic3(
        l::NTuple{3, Int}, 
        c::Tuple,
        range::UnitRange{Int}
    )

Calculate the seed and the c-hull of a toric variety of Picard rank 3.

# Arguments
- `l::NTuple{3, Int}`: Size of the three blocks of the map pi
- `c::Tuple`: Values inside the three blocks
- `range::UnitRange{Int}`: Coefficient range

# Examples
```jldoctest
julia> l = (2,2,2);

julia> c = ((0,0),(0,0),(0,0));

julia> get_seed_and_chull_of_toricpic3(l,c,-2:2);


```
"""
function get_seed_and_chull_of_toricpic3(
        l::NTuple{3, Int}, 
        c::Tuple,
        range::UnitRange{Int}
    )    
    l1, l2, l3 = l[1], l[2], l[3]
    c12, c13, c23 = c[1], c[2], c[3]
    
    # Validate length
    if length(c12) != l2 || length(c13) != l3 || length(c23) != l3
        return "Error: the entries of c must have the correct lengths"
    end
     
    c12_s = sum(c12)
    c13_s = sum(c13)
    c23_s = sum(c23)
    
    chull = DataFrame(x1 = Int[], x2 = Int[], x3 = Int[], x4 = String[])
    
    # First block
    a = [i for i in -(l1-1):-1]
    for index in a
        push!(chull, [index 0 0 "seed"])
        push!(chull, [index-c12_s -l2 0 "chull"])
        push!(chull, [index-c13_s -c23_s -l3 "chull"])
        push!(chull, [index-c12_s-c13_s -l2-c23_s -l3 "chull"])
    end
    
    # Second block
    b = [i for i in -(l2-1):-1]
    for index in b
        for r in range
            push!(chull, [r index 0 "seed"])
            push!(chull, [r index-c23_s -l3 "chull"])
        end
    end
    
    # Third block
    c = [i for i in -(l3-1):-1]
    all_tuples = generate_vectors(2, range)        
    for index in c
        for t in all_tuples
            push!(chull, [t[1] t[2] c "seed"])
        end
    end
      
    return  chull
end


@doc raw"""
    calculate_planary_parallelograms(coeffs::Vector{Int})

Calculate the planary parallelograms P1,P2 of a toric variety of Picard rank 3.

# Arguments
- `coeffs::Vector{Int}`: A coefficient vector of length 4.

# Examples
```jldoctest
julia> calculate_planary_parallelograms([1,2,3,4])
(Polyhedron in ambient dimension 2, Polyhedron in ambient dimension 2)


```
"""
function calculate_planary_parallelograms(coeffs::Vector{Int})
    if length(coeffs) != 4
        throw(ArgumentError("coeffs must be a vector of length 4")) 
    end
        
    a,b,c,d = coeffs[1], coeffs[2], coeffs[3], coeffs[4]
    p1 = convex_hull([-a-b-c+2 a-1; -a a-1; -b+d -c-d+1; c+d-2 -c-d+1])
    p2 = convex_hull([-a-b+1 a+b-2; d-1 -d; -a-b+1 a-c; d-1 -b-c-d+2])
    
    return (p1, p2)   
end