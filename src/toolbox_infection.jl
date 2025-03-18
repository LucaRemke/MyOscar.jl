########################################
# INFECTION PROCESS FOR DIFFERENT TORIC VARIETIES
########################################
# PENTAGON
########################################

#---------------------------------------
# Calculates the infection spreaded by the input vectors with the rules of the pentagon
#---------------------------------------
function calculate_infectionstep_of_pentagon(
        seq::Vector{Vector{Int64}}; 
        bound::Vector{Vector{Int64}}=Vector{Int64}[], 
        print_info::Bool=true
    )
    seq_set = Set(seq)
    infected_points = Vector{Int64}[]
    
    # Define rules as tuples: (offsets to check, offset to apply if rule matches)
    rules = [
        ([ [0,0,0], [1,1,-1], [1,0,0] ], [2,1,-1]),
        ([ [0,0,0], [-1,-1,1], [-1,0,0] ], [-2,-1,1]),
        ([ [0,0,0], [0,-1,1], [1,0,0] ], [1,-1,1]),
        ([ [0,0,0], [0,1,-1], [-1,0,0] ], [-1,1,-1]),
        ([ [0,0,0], [0,1,0], [1,1,-1] ], [1,2,-1]),
        ([ [0,0,0], [0,-1,0], [-1,-1,1] ], [-1,-2,1]),
        ([ [0,0,0], [-1,0,1], [0,1,0] ], [-1,1,1]),
        ([ [0,0,0], [1,0,-1], [0,-1,0] ], [1,-1,-1]),
        ([ [0,0,0], [0,-1,1], [-1,0,1] ], [-1,-1,2]),
        ([ [0,0,0], [0,1,-1], [1,0,-1] ], [1,1,-2])
    ]
    
    # Iterate over each vector in the sequence
    for vec in seq
        rule_number = 1
        for (targets, new_point_offset) in rules
            target_vectors = [vec + offset for offset in targets]
            infection = vec + new_point_offset
            
            # Check if all the target vectors are in the set
            # Check if the new infection is indeed new and not already contained (prevent duplication)
            # Check if the infection was not already found by another rule
            if issubset(target_vectors, seq_set) && !(in(infection, seq_set)) && 
                !(infection in infected_points)

                if is_empty(bound) || (infection in bound)
                    if print_info
                        v2 = target_vectors[2]
                        v3 = target_vectors[3]
                        superscript = rule_number % 2 == 1 ? "+" : "-"
                        subscript = ceil(Int, rule_number / 2)

                        #println("R$(subscript)$(superscript): $vec, $v2, $v3 \t --> $infection")
                        println("$(Tuple(infection)) \\xleftarrow{R_$(subscript)^$(superscript)} $(Tuple(vec)), $(Tuple(v2)), $(Tuple(v3))")
                    end
                    push!(infected_points, infection)
                end                
            end
            rule_number += 1
        end
    end
    
    if print_info
        println("\n$(length(infected_points)) points were infected.")
    end
    
    # Combine original sequence with new infected points
    output = vcat(seq, infected_points)
    
    return output
end;

#---------------------------------------
# Calculates alle infections after "steps" interior_lattice_points
# Checks if the infected points contain a target set, this can be considered as a termination rule for the algorithm
#---------------------------------------
function search_infections_of_pentagon(
        seq::Vector{Vector{Int64}}, 
        steps::Int64, 
        check_points::Set{Vector{Int64}};
        bound::Vector{Vector{Int64}}=Vector{Int64}[],
        print_info::Bool=true
    )
    all_infections = seq
    
    if print_info
        println("$(length(intersect(Set(seq), check_points))) of $(length(check_points)) points are already infected\n") 
    end 
    
    for i in 1:steps
        if print_info
            println("We reached step $i:")
            println("-------------------")
        end

        old_infection = copy(all_infections)
        all_infections = calculate_infectionstep_of_pentagon(all_infections; bound=bound, print_info=print_info)
        new_infections = setdiff(all_infections, old_infection)
        
        if issubset(check_points, Set(all_infections)) 
            if print_info
                println("All required points are infected")
            end
            return all_infections
        else
            required = intersect(new_infections, check_points)
            if print_info
                println("$(length(required)) of these points are in the target set.")
            end

            missing = length(setdiff(check_points, Set(all_infections)))
            if print_info 
                println("$missing points still to be infected.")
            end
        end
        
        if print_info
            println("###################\n") 
        end
            
    end
    
    return all_infections
end;

#---------------------------------------
# Checks if the searching for infection with input spreader is successful for the pentagon
#---------------------------------------
function search_infections_of_pentagon_bool(seq::Vector{Vector{Int64}}, steps::Int64, check_points::Set{Vector{Int64}})
    all_infections = seq  
    
    for i in 1:steps
        all_infections = calculate_infectionstep_of_pentagon(all_infections; print_info=false)
        
        if issubset(check_points, Set(all_infections))
            return [true, i]
        end            
    end
    
    return [false, nothing]
end;

#---------------------------------------
# Visualzie the infection process on the pentagon
#---------------------------------------
function visualize_infection_of_pentagon(
        seq::Vector{Vector{T}},
        plot_area::Vector{UnitRange{Int64}},
        color_vec::Vector{Int64},
        perspective::Tuple{Int64, Int64} = (10,20),
        plot_size::Tuple{Int64, Int64} = (300,300),
        point_size::Int64 = 5
    ) where T

    df = convert_to_df(seq)
    df.x4 = color_vec

    x_lim_min = first(plot_area[1])
    x_lim_max = last(plot_area[1])
    y_lim_min = first(plot_area[2])
    y_lim_max = last(plot_area[2])
    z_lim_min = first(plot_area[3])
    z_lim_max = last(plot_area[3])

    layer_plot = plot(
        xlims = (x_lim_min,x_lim_max), xticks = x_lim_min:1:x_lim_max,
        ylims = (y_lim_min,y_lim_max), yticks = y_lim_min:1:y_lim_max,
        zlims = (z_lim_min,z_lim_max), zticks = z_lim_min:1:z_lim_max,
        grid = :true, gridalpha = 0.6,
        framestyle = :box,
        camera = perspective,
        size = plot_size
    )

    for layer in plot_area[3]
        # Add the xy-plane of height layer
        z_plane = ones(length(plot_area[1]), length(plot_area[2])) * layer

        surface!(
            plot_area[1], plot_area[2], z_plane, 
            alpha = 0.4, color = :gray, colorbar=false
        )

        # Determine the elements of the exceptional sequence in the current layer
        df_seq_layer = filter(row -> row.x3 == layer, df)

        # Plot the elemetns in the current layer and label them
        colors = @. ifelse(df_seq_layer.x4 == 0, "green",
        @. ifelse(df_seq_layer.x4 == 1, "black", 
        ifelse(df_seq_layer.x4 == 3, "blue", "red")))

        plot!(
            df_seq_layer.x1, df_seq_layer.x2, df_seq_layer.x3, 
            seriestype = :scatter, framestyle = :axes, legend = :false, color = colors, markersize = point_size
        )
    end

    return layer_plot
end;

########################################
# HEXAGON
########################################

#---------------------------------------
# Calculates the infection spreaded by the input vectors with the rules of the hexagon 
#---------------------------------------
function calculate_infectionstep_of_hexagon(
    seq::Vector{Vector{Int64}}; 
    bound::Vector{Vector{Int64}}=Vector{Int64}[], 
    print_info::Bool=true
)
seq_set = Set(seq)
infected_points = Vector{Int64}[]

# Define rules as tuples: (offsets to check, offset to apply if rule matches)
rules = [
    ([ [0,0,0,0], [0,1,0,0], [1,0,0,0] ], [1,1,0,0]),
    ([ [0,0,0,0], [0,-1,0,0], [-1,0,0,0] ], [-1,-1,0,0]),
    ([ [0,0,0,0], [0,0,1,0], [1,0,0,0] ], [1,0,1,0]),
    ([ [0,0,0,0], [0,0,-1,0], [-1,0,0,0] ], [-1,0,-1,0]),
    ([ [0,0,0,0], [1,0,-1,1], [1,0,0,0] ], [2,0,-1,1]),
    ([ [0,0,0,0], [-1,0,1,-1], [-1,0,0,0] ], [-2,0,1,-1]),
    ([ [0,0,0,0], [0,0,1,0], [0,1,0,0] ], [0,1,1,0]),
    ([ [0,0,0,0], [0,0,-1,0], [0,-1,0,0] ], [0,-1,-1,0]),
    ([ [0,0,0,0], [0,1,-1,1], [0,1,0,0] ], [0,2,-1,1]),
    ([ [0,0,0,0], [0,-1,1,-1], [0,-1,0,0] ], [0,-2,1,-1]),
    ([ [0,0,0,0], [0,0,0,1], [0,0,1,0] ], [0,0,1,1]),
    ([ [0,0,0,0], [0,0,0,-1], [0,0,-1,0] ], [0,0,-1,-1]),
    ([ [0,0,0,0], [1,0,-1,1], [0,0,0,1] ], [1,0,-1,2]),
    ([ [0,0,0,0], [-1,0,1,-1], [0,0,0,-1] ], [-1,0,1,-2]),
    ([ [0,0,0,0], [0,1,-1,1], [0,0,0,1] ], [0,1,-1,2]),
    ([ [0,0,0,0], [0,-1,1,-1], [0,0,0,-1] ], [0,-1,1,-2]),
    ([ [0,0,0,0], [0,1,-1,1], [1,0,-1,1] ], [1,1,-2,2]),
    ([ [0,0,0,0], [0,-1,1,-1], [-1,0,1,-1] ], [-1,-1,2,-2])
]

# Iterate over each vector in the sequence
for vec in seq
    rule_number = 1
    for (targets, new_point_offset) in rules
        target_vectors = [vec + offset for offset in targets]
        infection = vec + new_point_offset
        
        # Check if all the target vectors are in the set
        # Check if the new infection is indeed new and not already contained (prevent duplication)
        # Check if the infection was not already found by another rule
        if issubset(target_vectors, seq_set) && !(in(infection, seq_set)) && 
            !(infection in infected_points)

            if is_empty(bound) || (infection in bound)
                if print_info
                    v2 = target_vectors[2]
                    v3 = target_vectors[3]
                    superscript = rule_number % 2 == 1 ? "+" : "-"
                    subscript = ceil(Int, rule_number / 2)

                    target_set = Set([-5, -4, -3, -2, 3, 4, 5])

                    if !any(x -> x in target_set, infection)

                    #println("R$(subscript)$(superscript): $vec, $v2, $v3 \t --> $infection")
                    println("$(Tuple(infection)) &\\xleftarrow{R_$(subscript)^$(superscript)} $(Tuple(vec)), $(Tuple(v2)), $(Tuple(v3)) \\\\[-0.25cm]")
                    end
                end
                push!(infected_points, infection)
            end                
        end
        rule_number += 1
    end
end

if print_info
    println("\n$(length(infected_points)) points were infected.")
end

# Combine original sequence with new infected points
output = vcat(seq, infected_points)

return output
end;

#---------------------------------------
# Calculates alle infections after "steps" interior_lattice_points
# Checks if the infected points contain a target set, this can be considered as a termination rule for the algorithm
#---------------------------------------
function search_infections_of_hexagon(
    seq::Vector{Vector{Int64}}, 
    steps::Int64, 
    check_points::Set{Vector{Int64}};
    bound::Vector{Vector{Int64}}=Vector{Int64}[],
    print_info::Bool=true
)
all_infections = seq

if print_info
    println("$(length(intersect(Set(seq), check_points))) of $(length(check_points)) points are already infected\n") 
end 

for i in 1:steps
    if print_info
        println("We reached step $i:")
        println("-------------------")
    end

    old_infection = copy(all_infections)
    all_infections = calculate_infectionstep_of_hexagon(all_infections; bound=bound, print_info=print_info)
    new_infections = setdiff(all_infections, old_infection)
    
    if issubset(check_points, Set(all_infections)) 
        if print_info
            println("All required points are infected")
        end
        return all_infections
    else
        required = intersect(new_infections, check_points)
        if print_info
            println("$(length(required)) of these points are in the target set.")
        end

        missing = length(setdiff(check_points, Set(all_infections)))
        if print_info 
            println("$missing points still to be infected.")
        end
    end
    
    if print_info
        println("###################\n") 
    end
        
end

return all_infections
end;

#---------------------------------------
# Checks if the searching for infection with input spreader is successful for the hexagon
#---------------------------------------
function search_infections_of_hexagon_bool(seq::Vector{Vector{Int64}}, steps::Int64, check_points::Set{Vector{Int64}})
all_infections = seq  

for i in 1:steps
    all_infections = calculate_infectionstep_of_hexagon(all_infections; print_info=false)
    
    if issubset(check_points, Set(all_infections))
        return [true, i]
    end            
end

return [false, nothing]
end;
