# This toolbox contains self-definied functions which are used in several other scripts
# If you made any changes updated the toolbox via
# "File" -> "Download as" -> toolbox.jl
# Move the file to the desired folder

# load required packages

using Pkg
using Oscar
using Plots
gr()
using DataFrames
using CSV
using LinearAlgebra
using Combinatorics

# Functions that constructs toric varieties with our ordering of the rays and maximal cones:
# the optional parameter non_redundant = true is mandatory for this purpose
####################
# The projective space

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
    
    P = normal_toric_variety(ray_generators, max_cones, non_redundant = true)
    
    return P
end;

####################
# The Hirzebruch surface

function define_hirzebruch_surface(r::Int; kleinschmidt::Bool=false)    
    
    if kleinschmidt
        ray_generators = [[0,1], [-r,-1], [1,0], [-1,0]]
        max_cones = [[1,4], [4,2], [2,3], [3, 1]]
    else
        ray_generators = [[-1,r], [0,1], [1,0], [0,-1]]
        max_cones = [[1,2], [2,3], [3,4], [4,1]]
    end
        
    H = normal_toric_variety(ray_generators, max_cones, non_redundant = true)
    
    return H
end;

####################
# The product P*P*P

function define_ppp()
    ray_generators = [[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]];
    max_cones = [[1,3,5], [1,3,6], [1,4,5], [1,4,6], [2,3,5], [2,3,6], [2,4,5], [2,4,6]];
    ppp = normal_toric_variety(ray_generators, max_cones, non_redundant = true)
    
    return ppp
end;

####################
# The blowup of P*P along cone(e1,e2)

function define_blowup_pp()
    ray_generators = [[1,0], [0,-1], [-1,0], [0,1], [1,1]];
    max_cones = [[4,5], [5,1], [1,2], [2,3], [3,4]];
    blowup_pp = normal_toric_variety(ray_generators, max_cones, non_redundant = true)
    
    return blowup_pp
end;
    
####################
# The hexagon

function define_hexagon()
    ray_generators = [[1,0], [0, 1], [-1, 1], [-1, 0], [0, -1], [1, -1]];
    max_cones = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
    hexagon = normal_toric_variety(ray_generators, max_cones, non_redundant = true)
    
    return hexagon
end;

# Auxiliary functions that are used in other self-definied functions
####################
# chaining of the functions divisor_class and toric_divisor_class

function get_representative_of_divisor(td::ToricDivisor)
    
    return divisor_class(toric_divisor_class(td)) 
end;

####################
# Transforms a RayVector or Vector with ringelements into an integer vector
function transform_rayvector(R::Any)
    vect = Int64[]
    for j in 1:length(R)
        entry = Int64(R[j])
        push!(vect,entry)
    end
    
    return vect
end;

####################
# Generates all vectors of length n with entries in range r
function generate_vectors(n, r)
    if n == 1
        return [[x] for x in r]
    else
        subvectors = generate_vectors(n - 1, r)
        vectors = []
        for x in r, subvector in subvectors
            push!(vectors, [x, subvector...])
        end
        return vectors
    end
end;

####################
# Convert a vector of vectors into a dataframe
function convert_to_df(vec)
    matrix = transpose(hcat(vec...)) 
    df = DataFrame(matrix, :auto)
    
    return df
end;

####################
# Sorts a vector of point clockwise around its centroid
function sort_points_clockwise(points)
    centroid = [sum(p[1] for p in points) / length(points), sum(p[2] for p in points) / length(points)]
    angles = [atan(p[2] - centroid[2], p[1] - centroid[1]) for p in points]
    sorted_indices = sortperm(angles)
    sorted_points = [points[i] for i in sorted_indices]
    
    return sorted_points
end;

# shows the generators and the relations of the prime divisors in the class-group
# returns the map pi: Div_T -> Cl as matrix
####################
function show_generators_and_relations_of_classgroup(v::NormalToricVariety; print_output::Bool=true)
    fan_v = fan(v)
    nrays_v = nrays(fan_v)        
    cl_group = class_group(v)
    cl_rank = rank(cl_group)
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

# Functions for calculating the Cartier data
####################
# calculates the Cartier-Data of a specific Cartier-Divisor
function calculate_cartier_data(td::ToricDivisor, cones::Bool=true)
    # implement a check if divisor is not cartier
    coeffs_td_int = transform_rayvector(coefficients(td))  
    var = toric_variety(td)
    rays_var = rays(var)
    maxcones_var = maximal_cones(fan(var))
    
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

# Function that checks if the input conditions characterizes the nef divisors in range
####################
function check_nef_conditions(var::NormalToricVariety, range::UnitRange{Int64}, condition_func)
    check = true
    
    all_coeffs = generate_vectors(nrays(fan(var)), range)
    
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

# Function that calculates a decomposition of a divisor into two nef divisors
# Caution when using the input basepoint_choice, an index other than 1 may not exist.
####################
function calculate_nef_pair(tv::NormalToricVariety, td::ToricDivisor; show_coeffs::Bool=true, basepoint_choice::Int=1)
    n = rank(class_group(tv))     
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

# Function that are used for calculating the immaculate locus
####################
# determines the immaculate line bundles in range  via brute-force
function get_immaculate_lb(variety::NormalToricVariety, range::UnitRange{Int64}; 
        perm=nothing, output_td::Bool=false, status::Bool=false)
    
    n_rays = nrays(variety)
    n_pic = rank(picard_group(variety)) 
    n_zeros = n_rays - n_pic
    
    immaculate = [];
    
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

####################
# determines the immaculate regions of a line bundle
function get_immaculate_regions(variety::NormalToricVariety)
    rays_var = rays(variety)
    nrays_var = nrays(variety)
    dim_pic = rank(picard_group(variety))
    pi = show_generators_and_relations_of_classgroup(variety; print_output=false)
    temptings = get_temptings_via_bruteforce(variety)
    
    regions = []
    
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

# Functions for plotting
####################
# plot a dataframe as 1/2/3D plot
# plot a dataframe as 1/2/3D plot
function plot_data(df::DataFrame; df_title=nothing, ax_labels::Bool=true, style::Symbol=:origin, camerapoint=nothing)
    
    n_cols = size(df,2)
    name_cols = names(df)
    
    x_min, x_max = extrema(df[!, name_cols[1]])
        
    if n_cols == 1        
        plot_df = plot(df[!, name_cols[1]], zeros(Int64, size(df,1)),
            xlims=(x_min-1.5,x_max+1.5),
            xticks=Int(x_min-1):1:Int(x_max+1), 
            yaxis=false, ylims=(-1,10), #ylims neu
            framestyle=:zerolines,
            seriestype=:scatter, size=(600,200) #size neu
        )
        
    elseif n_cols in [2,3]
        y_min, y_max = extrema(df[!, name_cols[2]])
        
        if n_cols == 2
            plot_df = plot(df[!, name_cols[1]], df[!, name_cols[2]],
                xticks=Int(x_min):1:Int(x_max),
                yticks=Int(y_min):1:Int(y_max),
                framestyle=style,
                seriestype=:scatter
            )
            
        else
            z_min, z_max = extrema(df[!, name_cols[3]])
            
            plot_df = plot(df[!, name_cols[1]], df[!, name_cols[2]], df[!, name_cols[3]],
                xticks=x_min:1:x_max,
                yticks=y_min:1:y_max,
                zticks=z_min:1:z_max,
                framestyle=style,
                seriestype=:scatter
            )
        
            if !isnothing(camerapoint)
                plot!(camera=camerapoint)
            end
            
            if ax_labels
                plot!(xlabel="x", ylabel="y", zlabel="z")
            end
        end

    else
        throw(ErrorException("Tuples of size $n_cols can not be plotted!")) 
    end
    
    plot!(legend=nothing)
    
    if !isnothing(df_title)
        plot!(title=df_title)
    end
    
    return plot_df
end;

####################
# creates a gif from a vector of plots
function create_gif(plots::Vector{Any}, path::String, dim::Tuple{Int64, Int64}; 
        frames_per_sec::Int64=20, records::UnitRange{Int64}=0:360, angle::Int64=165)
    
    all_plots = plot(plots..., size=dim)
    anim = Animation()
    
    for i in records
        plot!(camera = (i, angle))
        frame(anim)
    end
    
    gif_plot = gif(anim, path; fps=frames_per_sec)
    
    return gif_plot
    
end;

####################
# creates a gif from a dataframe
function create_gif(df::DataFrame, path::String; 
        frames_per_sec::Int64=20, records::UnitRange{Int64}=0:360, angle::Int64=165)
    
    name_cols = names(df)    
    plots = []
    anim = Animation()
    for i in records
        plot_i = plot(df[!, name_cols[1]], df[!, name_cols[2]], df[!, name_cols[3]],
            seriestype=:scatter,
            framestyle=:origin,
            label = nothing,
            camera = (i, angle)
        )
        push!(plots, plot_i)
        frame(anim)
    end
    
    gif_plot = gif(anim, path; fps=frames_per_sec)
    
    return gif_plot
    
end;

# Functions for calculating temptings of a toric variety
####################
# iterates over all subsets of Sigma(1) and checks the cohomologies of -sum(D_rho) rho in subset
function get_temptings_via_bruteforce(var::NormalToricVariety; zero_ray::Bool=false)
    v_fan = fan(var)
    n_rays = nrays(v_fan)
    pset_rays = collect(powerset(1:n_rays)) 
    temptings = Vector[]
    
    for subset in pset_rays
        coeffs = zeros(Int, n_rays)
        coeffs[subset] .= -1
        divisor = toric_divisor(var, coeffs)
        lb = toric_line_bundle(divisor)
        cohom = all_cohomologies(lb)
        
        if cohom != zeros(Int, dim(var)+1)
            push!(temptings, subset)
        end
    end
    
    if zero_ray
        # We cant use powerset(0:n_rays) since 0 is not a index
        temptings = [vec .- 1 for vec in temptings]
    end
    
    return temptings
    
end;

####################
# eliminates non temptings via the boundary condition
# this approach works if the toric variety is smooth since rho is surjective in this case
function search_nontemptings_via_boundary_condition(
        variety::NormalToricVariety; output_indizes::Bool=true, output_nontemptings::Bool=true)
    
    nrays_variety = nrays(variety)
    pset_nrays = collect(powerset(1:nrays_variety)) 
    ray_matrix = Matrix{Int64}(matrix(ZZ, rays(variety)))
    ray_vector = [ray_matrix[row, :] for row in 1:size(ray_matrix, 1)]
    
    zeros_im = zeros(Int64, nrays_variety)
    zeros_pim = Matrix{Int64}(transpose(hcat(zeros(Int64, dim(variety)))))
    
    output =[]
            
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

####################
# give non temptings via faces of the fan
function eliminate_faces_as_temptings(
        variety::NormalToricVariety; output_indizes::Bool=true, output_nontemptings::Bool=true)
    
    fan_variety = fan(variety)   
    nrays_variety = nrays(variety)
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

####################
# give temptings by determining the primitive collections
function get_temptings_as_primitive_collections(variety::NormalToricVariety;
        output_indizes::Bool=true, output_nontemptings::Bool=true, include_emptyset::Bool=true)
    
    nrays_variety = nrays(variety)
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


