########################################
# IMMACULATE LINE BUNDLES
########################################
# Determining toric divisors/line bundles with vanishing cohomology
########################################

#---------------------------------------
# Determine the immaculate line bundles in range
#---------------------------------------
function get_immaculate_lb(variety::NormalToricVariety, range::UnitRange{Int64}; 
    perm=nothing, output_td::Bool=false, status::Bool=false)

nrays = n_rays(variety)
n_pic = torsion_free_rank(picard_group(variety)) 
n_zeros = nrays - n_pic

immaculate = Vector{Int64}[];

coeffs = generate_vectors(n_pic, range) 
    
for coeff in coeffs
    coeff_complete = vcat(coeff, fill(0, n_zeros)) 
    
    if !isnothing(perm)
        coeff_complete = coeff_complete[perm]
    end
    
    td = toric_divisor(variety, coeff_complete)
    lb = toric_line_bundle(td)
    cohom = all_cohomologies(lb)
    
    if status
        println(coeff_complete) 
    end
    
    if cohom == repeat([0], length(cohom)) 
        if output_td
            push!(immaculate, lb)
        else
            push!(immaculate, coeff)
        end
    end
end

return immaculate
end;


########################################
# TEMPTING SUBSETS
########################################
# Determining the tempting subsets of a toric variety
########################################

#---------------------------------------
# Iterates over all subsets of Sigma(1) and checks the cohomologies of -sum(D_rho) rho in subset
#---------------------------------------
function calculate_tempting_subsets(var::NormalToricVariety; zero_ray::Bool=false)
    v_fan = polyhedral_fan(var)
    nrays = n_rays(v_fan)
    pset_rays = collect(powerset(1:nrays)) 
    temptings = Vector{Int64}[]
    
    for subset in pset_rays
        coeffs = zeros(Int, nrays)
        coeffs[subset] .= -1
        divisor = toric_divisor(var, coeffs)
        lb = toric_line_bundle(divisor)
        cohom = all_cohomologies(lb)
        
        if cohom != zeros(Int, dim(var)+1)
            push!(temptings, subset)
        end
    end
    
    if zero_ray
        # We cant use powerset(0:nrays) since 0 is not a index
        temptings = [vec .- 1 for vec in temptings]
    end
    
    return temptings
    
end;

#---------------------------------------
# Eliminates non temptings via the boundary condition
# This approach works if the toric variety is smooth since rho is surjective in this case
#---------------------------------------
function search_nontemptings_via_boundary_condition(
    variety::NormalToricVariety; output_indizes::Bool=true, output_nontemptings::Bool=true)

nrays_variety = n_rays(variety)
pset_nrays = collect(powerset(1:nrays_variety)) 
ray_matrix = Matrix{Int64}(matrix(ZZ, rays(variety)))
ray_vector = [ray_matrix[row, :] for row in 1:size(ray_matrix, 1)]

zeros_im = zeros(Int64, nrays_variety)
zeros_pim = Matrix{Int64}(transpose(hcat(zeros(Int64, dim(variety)))))

output = Vector{Int64}[]
        
for subset in pset_nrays
    
    B = -ray_matrix
    
    B[subset, :] .= B[subset, :] .* -1  
    P = polyhedron(B, zeros_im)
    
    if output_nontemptings            
        if P != convex_hull(zeros_pim)
            if output_indizes
                push!(output, subset)
            else
                push!(output, ray_vector[subset])
            end                
        end
    else
        if P == convex_hull(zeros_pim)
            if output_indizes
                push!(output, subset)
            else
                push!(output, ray_vector[subset])
            end                
        end
    end        
   
end    

return output
end;

#---------------------------------------
# Give non-temptings via faces of the fan
#---------------------------------------
function eliminate_faces_as_temptings(
    variety::NormalToricVariety; output_indizes::Bool=true, output_nontemptings::Bool=true)

fan_variety = polyhedral_fan(variety)   
nrays_variety = n_rays(variety)
pset_nrays = collect(powerset(1:nrays_variety)) 
ray_matrix = Matrix{Int64}(matrix(ZZ, rays(variety)))
ray_vector = [ray_matrix[row, :] for row in 1:size(ray_matrix, 1)]

output = Vector{Int64}[]

for i in 1:dim(fan_variety) 
    cones_i = cones(fan_variety, i)
    
    for j in cones_i
        rays_matrix_ij = Matrix{Int64}(matrix(ZZ, rays(j)))
        rays_vector_ij = [rays_matrix_ij[row, :] for row in 1:size(rays_matrix_ij, 1)]
        
        push!(output, indexin(rays_vector_ij, ray_vector)) 
        push!(output, indexin(setdiff(ray_vector, rays_vector_ij), ray_vector))          
    end
end

output = setdiff(pset_nrays, output)

if output_nontemptings
    output = setdiff(pset_nrays, output) #eliminate duplicates
end

if !output_indizes
    output = [ray_vector[output[i]] for i in 1:length(output)]
end

return output
end;

#---------------------------------------
# Give temptings by determining the primitive collections
#---------------------------------------
function get_temptings_as_primitive_collections(variety::NormalToricVariety;
    output_indizes::Bool=true, output_nontemptings::Bool=true, include_emptyset::Bool=true)

nrays_variety = n_rays(variety)
pset_nrays = collect(powerset(1:nrays_variety)) 
set_nrays = Set(1:nrays_variety)
all_sets = [Set(v) for v in pset_nrays]    
ray_matrix = Matrix{Int64}(matrix(ZZ, rays(variety)))
ray_vector = [ray_matrix[row, :] for row in 1:size(ray_matrix, 1)]

output = Set{Int64}[]

for pc in primitive_collections(variety)
    diff = setdiff(set_nrays, pc) 
    push!(output, pc)
    push!(output, diff) 
end

#remove duplicates
output = Set(output)

if include_emptyset
    push!(output, Set()) 
    push!(output, setdiff(set_nrays, Set()))
end

if output_nontemptings
    output = setdiff(all_sets, output)
end

output = collect.(output)
output = sort.(output)
output = sort(output)
output = sort(output, by=vec -> length(vec)) 

if !output_indizes
    output = [ray_vector[output[i]] for i in 1:length(output)]
end
    
return output

end;

#---------------------------------------
# Calculate the images of the cube under pi
#---------------------------------------
function get_images_of_cube_vertices(variety::NormalToricVariety; print_output::Bool=true)
    pi_v = show_generators_and_relations_of_classgroup(variety; print_output=false)
    nrays_v = n_rays(variety)
    cube_v = cube(nrays_v,-1,0) 
    vertices_v = vertices(cube_v)
    
    preimages = Vector{Int64}[]
    images = Vector{Int64}[]
    
    for vert in vertices_v
        vert_t = transform_rayvector(vert)
        push!(preimages, vert_t) 
        push!(images, pi_v * vert_t) 
    end
    
    results = Tuple{Vector{Int64}, Vector{Vector{Int64}}}[]
    
    for (im, pim) in zip(images, preimages)
        found = false
        for (existing_im, existing_pim) in results
            if existing_im == im
                push!(existing_pim, pim)
                found = true
                break
            end
        end
        if !found
            push!(results, (im, [pim])) 
        end
    end
    
    if print_output == true    
        for res in results
            res1 = join(res[1], ",")
            str = "($res1)"

            for i in 1:length(res[2])
                res2 = join(-1*res[2][i], "")
                str = str * " = -pi($res2)"
            end

            println(str)
            str=nothing
        end
    end
    
    return results
end;

#---------------------------------------
# Calculate the immaculate line bundles comming from the cube
#---------------------------------------
function get_polytope_of_cube_vertices(variety::NormalToricVariety; print_output::Bool=true)
    pi_v = show_generators_and_relations_of_classgroup(variety; print_output=false)
    nrays_v = n_rays(variety)
    cube_v = cube(nrays_v,-1,0) 
    vertices_v = vertices(cube_v)
    
    images = [pi_v * vec for vec in vertices_v]
    mat = transpose(hcat(images...)) 
    polytope = convex_hull(mat)
    
    # Determine the vertices corresponding to (non)tempting sets
    pset_nrays_v = collect(powerset(1:nrays_v)) 
    temptings = calculate_tempting_subsets(variety)
    
    temp_vertices = []
    for subset in temptings
        temp = zeros(Int64, nrays_v) 
        temp[subset] .= -1
        push!(temp_vertices, temp)
    end
    images_temp_vertices = unique([pi_v * vec for vec in temp_vertices])
    
    nontemptings = setdiff(pset_nrays_v, temptings)
    ntemp_vertices = []
    for subset in nontemptings
        ntemp = zeros(Int64, nrays_v) 
        ntemp[subset] .= -1
        push!(ntemp_vertices, ntemp)
    end
    images_ntemp_vertices = unique([pi_v * vec for vec in ntemp_vertices])
    
    ntemp_vertices = setdiff(transform_rayvector.(vertices(polytope)), images_temp_vertices)
    
    # Print output
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
    
    return images_ntemp_vertices
    
end;


########################################
# MACULATE REGIONS
########################################
# Determining the maculate regions coming from the tempting subsets
########################################

#---------------------------------------
# Calculate maculate regions of a line bundle
# Return value is a vector of polyhedrons
#---------------------------------------
function get_maculate_regions(variety::NormalToricVariety)
    rays_var = rays(variety)
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
end;

#---------------------------------------
# Calculate maculate regions, the H^i cones of these regions and the 
# immaculate line bundle in a spefic area of coefficients
# if error "Denominator must be 1" occurs, enlarge the cutout
#---------------------------------------
function cutout_maculate_regions(variety::NormalToricVariety, cutout::Polyhedron{QQFieldElem})
    lp_cutout = transform_rayvector.(lattice_points(cutout))
    temptations = calculate_tempting_subsets(variety)
    maculate_regions = get_maculate_regions(variety)
    lp_regions = Vector{DataFrame}[]
    lp_nonregions = lp_cutout
    
    # to visualize the regions we cutout a part of it
    # we obtain finitely many lattice point which can be plotted
    
    for region in maculate_regions
        isec_region = intersect(cutout, region)
        
        lp_region = transform_rayvector.(lattice_points(isec_region))
        lp_nonregions = [v1 for v1 in lp_nonregions if !(v1 in lp_region)]
        
        vertices_region = transform_rayvector.(vertices(isec_region))
        sorted_vertices_region = Vector{Vector{Int64}}(sort_points_clockwise(vertices_region))
        
        df_lp_region = convert_to_df(lp_region)
        df_sorted_vertices_region = convert_to_df(sorted_vertices_region)    
        df_lp_region = [df_lp_region, df_sorted_vertices_region]
        push!(lp_regions, df_lp_region)
     end
    
    # determine the H^i cone
    hcone = Int64[]
    for temp in temptations
        dr_coeff = [in(i, temp) ? -1 : 0 for i in 1:n_rays(variety)]
        cohom = transform_rayvector(all_cohomologies(toric_line_bundle(toric_divisor(variety, dr_coeff)))) 
        i = findfirst(x -> x != 0, cohom) - 1
        push!(hcone, i)
    end

    
    # determine the immaculate linebundles in cutout
    df_lp_nonregions = convert_to_df(lp_nonregions)
    
    return [lp_regions, hcone, df_lp_nonregions]
end;

#---------------------------------------
# Prints a complete description of all maculate regions
#---------------------------------------
function print_maculate_region_info(variety::NormalToricVariety)
    temptings = calculate_tempting_subsets(variety)
    regions = get_maculate_regions(variety)

    for i in 1:length(regions)
        region = regions[i]
        tempting = temptings[i]

        basepoint = vcat(transform_rayvector.(minimal_faces(region)[1])...)
        generators = transform_rayvector.(rays_modulo_lineality(region)[1])

        println("The maculate region of R = $tempting has:")
        println("\t Basepoint $basepoint")
        println("\t Generators $generators")
        i != length(regions) ? println("") : nothing

    end
end;


########################################
# GENERATING SEED
########################################

#---------------------------------------
# Calculate the seed and the chull of a toric variety of Picard rank 3
#---------------------------------------
function get_seed_and_chull_of_toricpic3(l::Tuple, c::Tuple, range::UnitRange{Int64})
    n = length(l)
    
    if n != 3
        return "Error: parameter l must have exactly three entries"
    end
    
    l1, l2, l3 = l[1], l[2], l[3]
    c12, c13, c23 = c[1], c[2], c[3]
        
    if length(c12) != l2 || length(c13) != l3 || length(c23) != l3
        return "Error: the entries of c must have the correct lengths"
    end
    
    a = [i for i in -(l1-1):-1]
    b = [i for i in -(l2-1):-1]
    c = [i for i in -(l3-1):-1]
    c12_s = sum(c12)
    c13_s = sum(c13)
    c23_s = sum(c23)
    
    chull = DataFrame(x1 = Int64[], x2 = Int64[], x3 = Int64[], x4 = String[])
    
    for index in a
        push!(chull, [index 0 0 "seed"])
        push!(chull, [index-c12_s -l2 0 "chull"])
        push!(chull, [index-c13_s -c23_s -l3 "chull"])
        push!(chull, [index-c12_s-c13_s -l2-c23_s -l3 "chull"])
    end
      
    for index in b
        for r in range
            push!(chull, [r index 0 "seed"])
            push!(chull, [r index-c23_s -l3 "chull"])
        end
    end
    
    all_tuples = generate_vectors(2, range)        
    
    for index in c
        for t in all_tuples
            push!(chull, [t[1] t[2] c "seed"])
        end
    end
      
    return  chull
end


########################################
# TORIC VARIETIES OF PICARD RANK 3
########################################
# Special treatment for toric varieties of Picard rank 3
########################################

#---------------------------------------
# Calculate the planary parallelograms
#---------------------------------------
function calculate_planary_parallelograms(coeffs::Vector{Int64})
    if length(coeffs) != 4
        throw(ArgumentError("coeffs must be a vector of length 4")) 
    end
        
    a,b,c,d = coeffs[1], coeffs[2], coeffs[3], coeffs[4]
    p1 = convex_hull([-a-b-c+2 a-1; -a a-1; -b+d -c-d+1; c+d-2 -c-d+1])
    p2 = convex_hull([-a-b+1 a+b-2; d-1 -d; -a-b+1 a-c; d-1 -b-c-d+2])
    
    return (p1, p2)
    
end