########################################
# INFECTION PROCESS FOR AN ARBITRARY VARIETY
########################################

# Generates the integral tuples
function generate_integral_tuples(n::Int, S::Set{Int})::Vector{Vector{Int}}
    # 1. Generate all subsets of S and convert to a list of sets
    all_subsets_S = collect(subsets(S))

    # 2. SORT the subsets by their size (cardinality), which ensures the order
    # from emptyset to full set (e.g., {}, {1}, {2}, {1,2})
    sort!(all_subsets_S, by=length, rev=true)

    result_vectors = Vector{Vector{Int}}()

    # --- Process each sorted subset SS ---
    for SS in all_subsets_S
        T = zeros(Int, n) # Initialize vector T of length n with zeros

        # Apply rules for indices k from 1 to n
        for k = 1:n
            if k in S
                # k is in S: T[k] = NOT (k in SS)
                if k in SS
                    # k in S AND k in SS => T[k] = 0 (Rule 2)
                    # T[k] = 0 (already set)
                else
                    # k in S AND k NOT in SS => T[k] = 1 (Rule 3)
                    T[k] = 1
                end
            # else: k not in S => T[k] = 0 (Rule 1, already set)
            end
        end

        push!(result_vectors, T)
    end
    
    return result_vectors
end


# Maps each integral tuple via pi into the Picard group, returns the rule of infection for a primitive collection
function map_intergral_tuples(P::AbstractMatrix, V_in::Vector{Vector{Int}})
    # V_in is the ordered output from generate_binary_vectors_ordered.

    # 1. Handle the empty case
    if isempty(V_in)
        # Determine the element type of the result based on P
        ResultType = eltype(P)
        return (Vector{Vector{ResultType}}(), zeros(ResultType, size(P, 1)))
    end

    # 2. Get all vectors EXCEPT the last one, preserving the sequential order
    other_vectors = V_in[1:end-1]
    
    # 3. Multiply the 'other_vectors' by P, preserving the sequential order
    # The map function iterates over 'other_vectors' in order.
    V_out = map(v -> P * v, other_vectors)
    
    # 4. Multiply the 'last_vector' by P
    last_result = P * V_in[end]

    # 5. Return the result as a tuple (Vector of results, last result)
    return (V_out, last_result)
end

# For a given rule, create the negative rule
function negate_rule(input_tuple::Tuple{Vector{Vector{T}}, Vector{T}}) where T <: Number
    V_out, last_result = input_tuple
    
    # 1. Negate all vectors in the first element (V_out)
    # The map function applies element-wise multiplication by -1 to each inner vector
    V_negated = map(v -> v * -1, V_out)
    
    # 2. Negate the last result vector
    last_negated = last_result * -1

    # 3. Return the result in the original tuple structure
    return (V_negated, last_negated)
end

# Apply a base change to the rules if you want to use another basis than the one Oscar chooses 
function apply_base_change(A::AbstractMatrix, input_rule::Tuple{Vector{Vector{T}}, Vector{T}}) where T <: Number
    V_out, last_result = input_rule
    
    # 1. Apply matrix A to all vectors in the first element (V_out)
    # The map function applies A*v to each inner vector in sequence, preserving order.
    V_transformed = map(v -> A * v, V_out)
    
    # 2. Apply matrix A to the last result vector
    last_transformed = A * last_result

    # 3. Return the newly transformed result in the original tuple structure
    return (V_transformed, last_transformed)
end

# Returns all infection rules 
function collect_infection_rules(
        nrays::Int, 
        all_pc::Vector{<:Any}, 
        pi::AbstractMatrix,
        A::Union{AbstractMatrix, Nothing}=nothing # Optional transition matrix A
    )::Vector{Tuple{Vector{Vector{eltype(pi)}}, Vector{eltype(pi)}}}
    
    # Initialize a vector to hold all resulting rule-tuples.
    # The type is specified to match the output of process_and_structure_vectors:
    # Tuple{Vector{Vector{ResultType}}, Vector{ResultType}}
    ResultType = eltype(pi)
    all_rules = Vector{Tuple{Vector{Vector{ResultType}}, Vector{ResultType}}}()
    
    # Determine which base change function to use: either identity or actual multiplication
    # We use a ternary operator to define the transformation function.
    transform_func = if A === nothing
        # If A is not provided, the transformation function is the identity function (f(x) = x)
        x -> x 
    else
        # If A IS provided, the transformation function is our matrix application function
        rule -> apply_base_change(A, rule) 
    end

    for pc in all_pc
        # 1. Generate Binary Vectors (Input for the next step)
        pc_integraltuples = generate_integral_tuples(nrays, pc)

        # 2. Process and Structure Vectors (Your 'map_intergral_tuples')
        pc_rule = map_intergral_tuples(pi, pc_integraltuples)
        
        # 3. Negate the Rule (Your 'negate_rule')
        pc_rule_neg = negate_rule(pc_rule)
        
        # Apply transformation and push the positive rule
        push!(all_rules, transform_func(pc_rule))
        
        # Apply transformation and push the negative rule
        push!(all_rules, transform_func(pc_rule_neg))
    end

    return all_rules
end

# Returns all infection rules 
function collect_infection_rules(tv::NormalToricVariety, A::Union{AbstractMatrix, Nothing}=nothing)
    # Collect necessary data
    pi = show_generators_and_relations_of_classgroup(tv; print_output=false)
    all_pc = primitive_collections(tv) 
    nrays = Oscar.n_rays(polyhedral_fan(tv))

    return collect_infection_rules(nrays, all_pc, pi, A)
end


# One infection step of a MES on a toric variety
function infection_step(
        variety::NormalToricVariety, 
        mes_set::Set{Vector{Int}}; 
        base_change::Union{AbstractMatrix, Nothing}=nothing,
        print_info::Bool=false
    )::Set{Vector{Int}}

    # 1. Initialize Sets for efficient lookups
    initial_mes_set = deepcopy(mes_set) # Original infected points
    becoming_ill_set = Set{Vector{Int}}() # New infections (for O(1) checks)

    # 2. Get Rules (Assuming collect_infection_rules uses variety to derive required info)
    # The output type of rules must be consistent: Tuple{Vector{Vector{Int}}, Vector{Int}}
    rules = collect_infection_rules(variety, base_change)

    # 3. Double Loop (O(N_mes * N_rules))
    for vec in initial_mes_set
        for (idx, (first, second)) in enumerate(rules)
            
            # Calculate origin points (points that must be infected for rule to fire)
            origins = [vec + shift for shift in first]
            target = vec + second

            # Check 1: All origin points are in the initial infected set (O(1) lookups)
            check1 = issubset(origins, initial_mes_set)

            # Check 2: The target point is not infected from the beginning (O(1) lookup)
            check2 = !(target in initial_mes_set)

            # Check 3: Has not been already infected in this step (O(1) lookup)
            check3 = !(target in becoming_ill_set)

            if check1 && check2 && check3
                push!(becoming_ill_set, target)

                # Print logic (using enumeration index idx)
                if print_info
                    superscript = idx % 2 == 1 ? "+" : "-"
                    subscript = ceil(Int, idx / 2)
                    
                    # Note: Using string interpolation for complex print statements
                    out_string = "$(Tuple(target)) &\\xleftarrow{R_$(subscript)^$(superscript)}"
                    out_string *= join([
                        " $(Tuple(el))" * (i == length(origins) ? " \\\\[-0.25cm]" : ",") 
                        for (i, el) in enumerate(origins)
                    ])
                    println(out_string)
                end
            end
        end
    end

    if print_info
        println("\n$(length(becoming_ill_set)) points were infected.")
    end

    # Combine the new infections with the old infections
    infected = union(initial_mes_set, becoming_ill_set)    
    return infected
end

# Infection procedure
function infection_procedure(
    variety::NormalToricVariety, 
    mes_set::Set{Vector{Int}},
    max_steps::Int,
    target_set::Set{Vector{Int}};
    base_change::Union{AbstractMatrix, Nothing}=nothing,
    print_info::Bool=false
)::Tuple{Set{Vector{Int}}, Union{Int, Nothing}}
    
    # 1. Initialization and Performance Setup
    # CRITICAL: Start with a Set for fast lookups.
    currently_infected = deepcopy(mes_set) # Start with a copy to avoid modifying the input set
    
    if print_info
        initial_intersec = intersect(currently_infected, target_set)
        println("$(length(initial_intersec)) of $(length(target_set)) points are initially infected\n")
    end
    
    # 2. Iteration Loop
    for current_step in 1:max_steps
        if print_info
            println("Infection step $current_step:")
            println("-------------------")
        end
        
        old_infection_size = length(currently_infected)
        
        # Assume infection_step (the previous function) is also improved 
        # to take and return a Set{Vector{Int}}.
        newly_infected = infection_step(variety, currently_infected; base_change, print_info)
        
        # Calculate new infections (fast set difference)
        new_infections_only = setdiff(newly_infected, currently_infected)
        
        # Update the currently infected set efficiently
        union!(currently_infected, newly_infected)
        
        # 3. Stagnation Check (New Early Exit)
        if length(currently_infected) == old_infection_size
            if print_info
                println("Infection stabilized: No new points were infected.")
                println("###################\n")
            end
            break # Exit the loop if no progress is made
        end

        # 4. Completion Check
        if issubset(target_set, currently_infected) 
            if print_info
                println("All required points are infected in $(current_step) steps.")
            end
            return (currently_infected, current_step)
        end
        
        # 5. Reporting (Only if print_info is true)
        if print_info
            # Required points gained in this step
            required_gained = intersect(new_infections_only, target_set)
            
            # Points still missing from the target
            missing_count = length(setdiff(target_set, currently_infected))
            
            println("$(length(required_gained)) target points gained in this step.")
            println("$missing_count points still missing from the target set.")
            println("###################\n")
        end
    end
    
    # 6. Failure to Complete
    if print_info
        println("The maximal number of $(max_steps) steps has been exceeded without infecting all targets.")
    end
    
    # Return the final set of infected points and nothing for the step count.
    return (currently_infected, nothing) 
end






########################################
# INFECTION PROCESS FOR THE PENTAGON
########################################


@doc raw"""
    calculate_infectionstep_of_pentagon(
        seq::Vector{Vector{Int}}; 
        bound::Vector{Vector{Int}}=Vector{Int64}[], 
        print_info::Bool=true,
        print_latex::Bool=false
    )

Calculate all points that are infected by the input points using the rules of the pentagon in exactly one step.

# Arguments
- `seq::Vector{Vector{Int}}`: Start sequence
- `bound::Vector{Vector{Int}}`: (optional) Determines the search range
- `print_info::Bool`: (optional) Prints additional information about the process
- `print_latex::Bool`: (optional) Prints latex code that describes the infection process

# Examples
```jldoctest
julia> mes = [[0,0,0], [0,-1,1], [-1,0,1], [0,0,1], [0,0,2]];

julia> calculate_infectionstep_of_pentagon(mes; print_info=false)
6-element Vector{Vector{Int64}}:
 [0, 0, 0]
 ⋮
 [-1, -1, 2]

```
"""
function calculate_infectionstep_of_pentagon(
        seq::Vector{Vector{Int}}; 
        bound::Vector{Vector{Int}}=Vector{Int64}[], 
        print_info::Bool=true,
        print_latex::Bool=false
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

                        if print_latex
                            println("$(Tuple(infection)) &\\xleftarrow{R_$(subscript)^$(superscript)} $(Tuple(vec)), $(Tuple(v2)), $(Tuple(v3)) \\\\[-0.25cm]")
                        else
                            println("$(Tuple(infection)) using R_$(subscript)^$(superscript) for $(Tuple(vec)), $(Tuple(v2)), $(Tuple(v3))")
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
end


@doc raw"""
    search_infections_of_pentagon(
        seq::Vector{Vector{Int}}, 
        steps::Int, 
        check_points::Set{Vector{Int}};
        bound::Vector{Vector{Int}}=Vector{Int64}[],
        print_info::Bool=true,
        print_latex::Bool=false
    )

Calculate all infections after `steps`-iterations.
After each step, it is checked if these points are a superset of `check_points`.
In this case the algorithm terminates.


# Arguments
- `seq::Vector{Vector{Int}}`: Start sequence
- `steps::Int`: Number of (maximal) iterations of the infection process
- `check_points::Set{Vector{Int}}`: Determines a stop criterion
- `bound::Vector{Vector{Int}}`: (optional) Determines the search range
- `print_info::Bool`: (optional) Prints additional information about the process
- `print_latex::Bool`: (optional) Prints latex code that describes the infection process

# Examples
```jldoctest
julia> mes = [[0,0,0], [0,-1,1], [-1,0,1], [0,0,1], [0,0,2]];

julia> target_set = Set([[x, y, z] for x in -1:1, y in -1:1, z in -1:1]);

julia> search_infections_of_pentagon(mes, 3, target_set; print_info=false)
15-element Vector{Vector{Int64}}:
 [0, 0, 0]
 ⋮
 [0, -2, 1]

```
"""
function search_infections_of_pentagon(
        seq::Vector{Vector{Int}}, 
        steps::Int, 
        check_points::Set{Vector{Int}};
        bound::Vector{Vector{Int}}=Vector{Int64}[],
        print_info::Bool=true,
        print_latex::Bool=false
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
        all_infections = calculate_infectionstep_of_pentagon(all_infections; bound=bound, print_info=print_info, print_latex=print_latex)
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
end


@doc raw"""
    search_infections_of_pentagon_bool(
        seq::Vector{Vector{Int}},
        steps::Int,
        check_points::Set{Vector{Int}}
    )

Calculate if the starting sequence infects all points in the target set.
Also returns the number of iterations which is needed for this.

# Arguments
- `seq::Vector{Vector{Int}}`: Start sequence
- `steps::Int`: Number of (maximal) iterations of the infection process
- `check_points::Set{Vector{Int}}`: Determines a stop criterion

# Examples
```jldoctest
julia> mes = [[0,0,0], [0,-1,1], [-1,0,1], [0,0,1], [0,0,2]];

julia> target_set = Set([[x, y, z] for x in -1:1, y in -1:1, z in -1:1]);

julia> search_infections_of_pentagon_bool(mes, 13, target_set)
2-element Vector{Int64}:
  1
 12

```
"""
function search_infections_of_pentagon_bool(
        seq::Vector{Vector{Int}},
        steps::Int,
        check_points::Set{Vector{Int}}
    )
    all_infections = seq  
    
    for i in 1:steps
        all_infections = calculate_infectionstep_of_pentagon(all_infections; print_info=false)
        
        if issubset(check_points, Set(all_infections))
            return [true, i]
        end            
    end
    
    return [false, nothing]
end


########################################
# INFECTION PROCESS FOR THE HEXAGON
########################################


@doc raw"""
    calculate_infectionstep_of_hexagon(
        seq::Vector{Vector{Int}}; 
        bound::Vector{Vector{Int}}=Vector{Int64}[], 
        print_info::Bool=true,
        print_latex::Bool=false
    )

Calculate all points that are infected by the input points using the rules of the hexagon in exactly one step.

# Arguments
- `seq::Vector{Vector{Int}}`: Start sequence
- `bound::Vector{Vector{Int}}`: (optional) Determines the search range
- `print_info::Bool`: (optional) Prints additional information about the process
- `print_latex::Bool`: (optional) Prints latex code that describes the infection process

# Examples
```jldoctest
julia> mes = [[0,0,0,0], [0,1,-1,1], [0,1,0,0], [0,1,0,1], [1,0,-1,0], [1,1,-1,1]];

julia> calculate_infectionstep_of_hexagon(mes; print_info=false)
11-element Vector{Vector{Int64}}:
 [0, 0, 0, 0]
 ⋮
 [-1, 1, 0, 0]


```
"""
function calculate_infectionstep_of_hexagon(
        seq::Vector{Vector{Int}}; 
        bound::Vector{Vector{Int}}=Vector{Int64}[], 
        print_info::Bool=true,
        print_latex::Bool=false
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
                            if print_latex
                                println("$(Tuple(infection)) &\\xleftarrow{R_$(subscript)^$(superscript)} $(Tuple(vec)), $(Tuple(v2)), $(Tuple(v3)) \\\\[-0.25cm]")
                            else
                                println("$(Tuple(infection)) using R_$(subscript)^$(superscript) for $(Tuple(vec)), $(Tuple(v2)), $(Tuple(v3))")
                            end                        
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
end


@doc raw"""
    search_infections_of_hexagon(
        seq::Vector{Vector{Int}}, 
        steps::Int, 
        check_points::Set{Vector{Int}};
        bound::Vector{Vector{Int}}=Vector{Int64}[],
        print_info::Bool=true,
        print_latex::Bool=false
    )

Calculate all infections after `steps`-iterations.
After each step, it is checked if these points are a superset of `check_points`.
In this case the algorithm terminates.


# Arguments
- `seq::Vector{Vector{Int}}`: Start sequence
- `steps::Int`: Number of (maximal) iterations of the infection process
- `check_points::Set{Vector{Int}}`: Determines a stop criterion
- `bound::Vector{Vector{Int}}`: (optional) Determines the search range
- `print_info::Bool`: (optional) Prints additional information about the process
- `print_latex::Bool`: (optional) Prints latex code that describes the infection process

# Examples
```jldoctest
julia> mes = [[0,0,0,0], [0,1,-1,1], [0,1,0,0], [0,1,0,1], [1,0,-1,0], [1,1,-1,1]];

julia> target_set = Set([[w, x, y, z] for w in -1:1, x in -1:1, y in -1:1, z in -1:1]);

julia> search_infections_of_hexagon(mes, 3, target_set; print_info=false)
44-element Vector{Vector{Int64}}:
 [0, 0, 0, 0]
 ⋮
 [1, 0, 0, 0]

```
"""
function search_infections_of_hexagon(
        seq::Vector{Vector{Int}}, 
        steps::Int, 
        check_points::Set{Vector{Int}};
        bound::Vector{Vector{Int}}=Vector{Int64}[],
        print_info::Bool=true,
        print_latex::Bool=false
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
        all_infections = calculate_infectionstep_of_hexagon(all_infections; bound=bound, print_info=print_info, print_latex=print_latex)
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
end

@doc raw"""
    search_infections_of_hexagon_bool(
        seq::Vector{Vector{Int}},
        steps::Int,
        check_points::Set{Vector{Int}}
    )

Calculate if the starting sequence infects all points in the target set.
Also returns the number of iterations which is needed for this.

# Arguments
- `seq::Vector{Vector{Int}}`: Start sequence
- `steps::Int`: Number of (maximal) iterations of the infection process
- `check_points::Set{Vector{Int}}`: Determines a stop criterion

# Examples
```jldoctest
julia> mes = [[0,0,0,0], [0,1,-1,1], [0,1,0,0], [0,1,0,1], [1,0,-1,0], [1,1,-1,1]];

julia> target_set = Set([[w, x, y, z] for w in -1:1, x in -1:1, y in -1:1, z in -1:1]);

julia> search_infections_of_hexagon_bool(mes, 13, target_set)
2-element Vector{Int64}:
  1
 10

```
"""
function search_infections_of_hexagon_bool(
        seq::Vector{Vector{Int}}, 
        steps::Int, 
        check_points::Set{Vector{Int}}
    )
    all_infections = seq  

    for i in 1:steps
        new_infections = calculate_infectionstep_of_hexagon(all_infections; print_info=false)

        # Check if all check_points are now infected
        if issubset(check_points, Set(new_infections))
            return [true, i]
        else
            # Stop if no change compared to previous step
            if new_infections == all_infections            
                return [false, nothing]
            end
        end

        all_infections = new_infections
    end

    return [false, nothing]
end