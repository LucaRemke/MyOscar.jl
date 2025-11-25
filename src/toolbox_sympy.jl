# add function if it works when pushing on github
# This is for all functions using sympy
# This file will be replaced by toolbox_sympy.jl if the push on github works

########################################
# load required packages
########################################
using PyCall
pyimport("sympy")
using SymPy

########################################
# CARTIER DATA OF TORIC DIVISORS
########################################

#---------------------------------------
# Calculates the Cartier-Data of general Cartier-Divisor
#---------------------------------------
function calculate_cartier_data(coeffs::Vector{Sym{PyObject}}, var::NormalToricVariety, cones::Bool=true)
    
    # implement a check if divisor is not cartier
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
        ind = Int[]

        for r in rays_c
            i = findfirst(item -> item == r, rays_var)
            push!(ind,i)
        end

        A = rays_matrix[ind,:]
        b = -coeffs[ind]
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
# Evaluate Cartier-Data of general Cartier-Divisor for a specific choice of coordinates
#---------------------------------------
function evaluate_cartier_data(data::Vector{Any}, symbols::Vector{Sym}, coeffs::Vector{Int})
    
    if length(symbols) != length(coeffs)
        throw(ArgumentError("Length of the input vectors does not match.")) 
    end
    
    all_data = []
    
    for i in 1:length(data)
        cdata = data[i][2]
        
        for j in 1:length(symbols)
            cdata = cdata.subs(symbols[j],coeffs[j])
        end
        
        cdata = vec(convert(Matrix{Int}, N(cdata))) 
        
        push!(all_data, cdata)
        
    end
    
    return all_data
end;


########################################
# BLOW-UPS
########################################

#---------------------------------------
# Pullback phiD of a divisor D via a blowup phi
#---------------------------------------
function pullback_divisor_of_blowup(
        morphism::Oscar.ToricBlowupMorphism{NormalToricVariety, NormalToricVariety},
        divisor::Vector{Sym{PyObject}}
    )
    origin = codomain(morphism)
    target = domain(morphism)

    cartierdata_origin = calculate_cartier_data(divisor, origin, true)
    cartierdata_mapped = []

    for data in cartierdata_origin
        cone = data[1]
        dim_cone = dim(cone)
        for c in cones(target, dim_cone)
            if issubset(c, cone) == true
                new_data = [rays(c), data[2]]
                push!(cartierdata_mapped, new_data)
                break
            end
        end
    end

    divisor_mapped = Sym{PyObject}[]

    for r in rays(target)
        for data in cartierdata_mapped
            cone = data[1]
            if r in cone
                r = transform_rayvector(r)
                lambda = -transpose(r)*data[2]
                push!(divisor_mapped, lambda)
                break
            end
        end
    end

    return divisor_mapped
end;


########################################
# EXCEPTIONAL SEQUENCES
########################################

#---------------------------------------
# Substitute the Symbol p such that the first occurence of p in a term is always p 
#---------------------------------------
function unify_sequence(V::Vector{Vector{Sym{PyObject}}}; param_sym::Sym{PyObject}=Sym("p"))
    for (i, vector) in enumerate(V)
        for (j, element) in enumerate(vector)
            # Try to convert the current element to an Integer
            convert = N(element)
            # The first time the convertion gives not an integer is the desired value
            if typeof(convert) == Sym{PyObject}
                # depending on the signs, the substitution differs
                sgn_p = N(convert.coeff(param_sym)) 
                if sgn_p == -1
                    subst = convert                
                else
                    subst = param_sym + (-1*N(convert)(param_sym=>0))
                end

                seq = [[subs(expr, param_sym, subst) for expr in vec] for vec in V]

                return seq
            end
        end
    end
    # if no symbol was found, return the input vector
    return V
end;

#---------------------------------------
# Calculate the meses of a given list V modulo helixing
#---------------------------------------
function calculate_exceptional_sequences_modulo_helixing(
        V::Vector{Vector{Vector{T}}},
        cd::Vector{Int}
    ) where T
        
    configuration = Vector{Vector{T}}[]

    for v in V
        set = Set(configuration)
        hex = calculate_all_helixing(v, cd)
        check = false
        for h in hex
            if unify_sequence(h) in set
                check = true
                break
            end
        end
        if check == false
            push!(configuration, v)
        end
    end

    return configuration
end;

#---------------------------------------
# Calculate the meses of a given list V modulo dualizing
#---------------------------------------
function calculate_exceptional_sequences_modulo_dualizing(
        V::Vector{Vector{Vector{T}}},
        cd::Vector{Int}
    ) where T
        
        configuration = Vector{Vector{T}}[]

        for v in V
            set = Set(configuration)
            dual = calculate_dualizing(v, cd)

            if !(unify_sequence(dual) in set)
                push!(configuration, v)
            end
        end

    return configuration
end;

#---------------------------------------
# Checks which flippings gives exceptional sequences
#---------------------------------------
function are_flippings_exceptional(
        var::Symbol, 
        V::Vector{Vector{Vector{T}}};
        par::Union{Int, Nothing}=nothing,
        s::Int=1
    ) where T

    l = length(V[1][1])
    perm = collect(permutations(1:l)) 

    flipping_all = Any[] #I got problems with push! when defining a precise type here

    for p in perm
        flipping_p = Vector{Vector{T}}[]
        cond = true
        for v in V
            flipping_p_v = calculate_flipping(v,p)
            push!(flipping_p, flipping_p_v)
            if !(is_exceptional(var, flipping_p_v; param=par, sgn=s) == true)
                cond = false
            end
        end
        insert = (p, flipping_p)
        push!(flipping_all, insert) 
        println("For permutation $p the sequences are exceptional: $cond")
    end           

    return flipping_all
end;

#---------------------------------------
# Finds the first occurence of a parameter in a sequence
#---------------------------------------
function find_parameter(seq::Vector{Vector{T}}) where T
    for (i, vector) in enumerate(seq)
        for (j, element) in enumerate(vector)
            convert = N(element)
            if typeof(convert) == Sym{PyObject}
                return (i,j)
            end
        end
    end
    return nothing
end


#---------------------------------------
# Substitute a sequence with a parameter with the value
#---------------------------------------
function substitute_sequence(
        seq::Vector{Vector{Sym{PyObject}}},
        param::Sym{PyObject},
        value::Int
    )
    seq_sub = [[subs(expr, param, value) for expr in vec] for vec in seq]
    params = find_parameter(seq_sub)
    if isnothing(params)
        return transform_rayvector.(seq_sub)
    else
        return seq_sub
    end
end

#---------------------------------------
# Substitute a sequence with a parameter with the value at the same index of the integer sequence
#---------------------------------------
function substitute_configuration(
        seq_int::Vector{Vector{Int}}, 
        seq_p::Vector{Vector{Sym{PyObject}}}; 
        param_sym::Sym{PyObject}=Sym("p")
    )
    vec, elem = find_parameter(seq_p)
    subst_var = seq_int[vec][elem]
    subst_seq_p = substitute_sequence(seq_p, param_sym, subst_var)
    
    return subst_seq_p
end;

#---------------------------------------
# Look if a concrete sequence corresponds to a sequence with general description
#---------------------------------------
function find_sequence_in_configurations(
        seqs_int::Vector{Vector{Int}},
        seqs_p::Vector{Vector{Vector{Sym{PyObject}}}};
        return_sequence::Bool=false
    )
    
    found = Vector{Vector{Sym{PyObject}}}[]

    for seq_p in seqs_p
        seq_p_subst = substitute_configuration(seqs_int, seq_p)
        if seq_p_subst == seqs_int
            if return_sequence == true
                push!(found, seq_p)
            else
                return true
            end
        end
    end

    if return_sequence == true
        return found
    else
        return false
    end;
end;

#---------------------------------------
# Checks if a vector of sequences can be described by a general rule
#---------------------------------------
function find_configurations(
        D::Dict{Vector{Int}, Vector{Vector{Vector{Int}}}};
        check_result::Bool=true, 
        param_sym::Sym{PyObject}=Sym("p")
    )
    general_description = Dict{Vector{Int}, Vector{Vector{Sym{PyObject}}}}()
    
    for conf in D
        parameter = nothing
        conf_key, conf_seq = conf[1], conf[2]
        
        # We can not find a general description for a single sequence
        if length(conf_seq) == 1
            continue
        end
        
        # We have at least 2 meses if we are here
        conf_seq_1, conf_seq_2 = conf_seq[1], conf_seq[2]
        
        # Build the difference to the the non-trivial elements
        conf_seq_diff = conf_seq_1 - conf_seq_2
        
        # Now we create a general desription of the sequences of the current type
        ## first we initialise it
        conf_seq_general = Vector{Sym{PyObject}}[]
        
        ## loop over the differences to define the vectors with general description
        for (pos, vec) in enumerate(conf_seq_diff)
            ### initialise vector for current iteration
            conf_seq_general_vec = Sym{PyObject}[]
            
            ### if there are no differences the vectors is fixed
            if vec == zeros(length(vec))
                conf_seq_general_vec = conf_seq_1[pos]
                push!(conf_seq_general, conf_seq_general_vec)
            ### if there are differences, things become more complicated
            else
                for (ind, elem) in enumerate(vec)
                    conf_seq_1_elem = conf_seq_1[pos][ind]
                    ### if element of difference vector is zero use element of origin vector
                    if elem == 0
                        push!(conf_seq_general_vec, conf_seq_1_elem)
                    ### otherwise we have to express it in terms of p
                    else                        
                        ### check if we already set the reference integer for p
                        #### if it is empty, we give it a value, this is our starting point for p
                        if parameter == nothing
                            parameter = conf_seq_1_elem
                            push!(conf_seq_general_vec, param_sym)
                        else
                            if sign(parameter) == sign(conf_seq_1_elem)
                                diff_elem = conf_seq_1_elem - parameter
                                vec_elem = param_sym + diff_elem
                            else
                                diff_elem = conf_seq_1_elem + parameter
                                vec_elem = -param_sym + diff_elem
                            end
                            push!(conf_seq_general_vec, vec_elem)
                        end
                    end
                end
                push!(conf_seq_general, conf_seq_general_vec) 
            end
        end
        general_description[conf_key] = conf_seq_general

        # Check if the general condition fits for all sequences
        if check_result == true
            for seq in conf_seq
                check = check_single_configurations_vector(conf_seq_general, seq)
                if check == false
                    throw(ErrorException("The defined general description does not fit for all sequences of this type!")) 
                end                
            end
        end
    end
        
    return general_description
end;

#---------------------------------------
# Checks if each single encoding fits into a general encoding
#---------------------------------------
function check_single_configurations(
        D_general::Dict{Vector{Int}, Vector{Vector{Sym{PyObject}}}},
        D_single::Dict{Vector{Int}, Vector{Vector{Vector{Int}}}}
    )

    values_general = collect(values(D_general))
    values_single = map(v -> v[1], values(D_single)) 

    for i in values_single
        found = find_sequence_in_configurations(i, values_general)
        if found == false
            return false
        end
    end

    return true
end;

#---------------------------------------
# Checks if each single encoding fits into a general encoding where inputs are vectors
#---------------------------------------
function check_single_configurations_vector(
        V_general::Vector{Vector{Sym{PyObject}}},
        V_single::Vector{Vector{Int}}
    )
    found = find_sequence_in_configurations(V_single, [V_general])
    if found == false
        return false
    else
        return true
    end
end;


########################################
# AUGMENTATIONS
########################################

"""
    create_symbolic_vector(l::Int, base::Vector{Int}; prefix="s")

Creates a vector of length l where elements at indices specified by base 
are unique SymPy symbols (e.g., s1, s2, ...) and all other elements are 0.

# Arguments
- l: The total length of the resulting vector.
- base: A vector of 1-based indices (subset of 1:l) where symbols should be placed.
- prefix: The string prefix for the unique symbols created (e.g., 's' -> s1, s2, ...).
"""
function create_symbolic_vector(l::Int, base::Vector{Int}; prefix="s")
    
    # This checks if the indices in base are within the bounds 1:l.
    if !isempty(base)
        checkbounds(1:l, base)
    end
    
    # Create unique SymPy symbols based on the length of BASE
    num_symbols = length(base)
    unique_symbols = [Sym("$(prefix)$(i)") for i in 1:num_symbols]
    
    # Initialize the result vector. 
    result_vector = fill(Sym(0), l) # Sym(0) is a Sym{PyObject} representing zero
    
    # Map symbols to the BASE indices using a single assignment
    result_vector[base] = unique_symbols
    
    return result_vector, unique_symbols
end



"""
    create_mapper(
        symbolic_vector::Vector, 
        symbols_to_map::Vector{Sym{PyObject}}
    )

Creates and returns a function that maps a numerical input vector (vec) 
to a numerical output vector based on the linear combinations defined 
by the symbolic_vector.

# Arguments
- symbolic_vector: The vector containing linear symbolic expressions 
                   (e.g., [a, -b, a+b]).
- symbols_to_map: The vector of unique symbols defining the order of 
                  input components (e.g., [a, b]).

# Returns
- A function `(vec::Vector{T}) -> output_vec::Vector{T}`.
"""
function create_mapper(
    symbolic_vector::Vector, 
    symbols_to_map::Vector{Sym{PyObject}}
)
    input_length = length(symbols_to_map)
    output_length = length(symbolic_vector)

    # 1. Create a map from symbol to input index (a -> 1, b -> 2, ...)
    symbol_to_index = Dict(s => i for (i, s) in enumerate(symbols_to_map))

    # 2. Pre-calculate the structure (coefficients and indices) for each output element
    structure = Vector{Vector{Tuple{Int, Int}}}() # [ [(coeff, index), ...], ...]

    for expr in symbolic_vector
        term_list = Vector{Tuple{Int, Int}}()
        
        # SymPy method to get the constant term (coefficient of 1)
        constant_term = expr.subs(Dict(s => 0 for s in symbols_to_map))

        # Handle the constant term (which corresponds to vec[0] conceptually, but is just a constant)
        if !iszero(constant_term)
            push!(term_list, (Int(N(constant_term)), 0)) # Use index 0 for constant
        end

        for symbol in symbols_to_map
            # Get the coefficient of the symbol (e.g., coeff of 'a' in 'a+b' is 1)
            coeff_sym = expr.coeff(symbol)
            
            # Use N() to convert the coefficient from Sym{PyObject} to a Julia number
            coeff = Int(N(coeff_sym)) 
            
            if !iszero(coeff)
                index = symbol_to_index[symbol]
                push!(term_list, (coeff, index))
            end
        end
        push!(structure, term_list)
    end

    # 3. Return the numerical closure
    return function (vec::Vector{T}) where T <: Number
        
        if length(vec) != input_length
            throw(ArgumentError("Input vector must have length $input_length, but got $(length(vec))."))
        end

        output_vec = Vector{T}(undef, output_length)
        
        # Iterate over the pre-calculated structure to compute the result
        for (i, terms) in enumerate(structure)
            result = zero(T) # Start the calculation
            
            for (coeff, index) in terms
                if index == 0
                    # Add the constant term
                    result += coeff 
                else
                    # Add coeff * vec[index]
                    result += coeff * vec[index]
                end
            end
            output_vec[i] = result
        end

        return output_vec
    end
end


function calculate_exceptional_sequences_after_blowup(
        variety::NormalToricVariety,
        blowup_ray::Vector{Int},
        sequences::Union{Vector{Vector{T}},Vector{Vector{Vector{T}}}}
    ) where T
    
    # Setup Variety Properties
    fan = polyhedral_fan(variety)
    number_rays = n_rays(fan)
    
    # Define Blowup + Exceptional Divisor
    blowup_morphism = blow_up(variety, blowup_ray; coordinate_name="r$(number_rays+1)")
    blowup = domain(blowup_morphism)
    pi_blowup = show_generators_and_relations_of_classgroup(blowup; print_output=false)
    exceptional_divisor = get_exceptional_divisor(blowup_morphism)
    
    # Define Pullback Function
    base, classes = compute_divisor_classes_and_base(variety)
    variety_base, variety_symbols = create_symbolic_vector(number_rays, base)

    pullback_base = pullback_divisor_of_blowup(blowup_morphism, variety_base)
    pullback_base_class = pi_blowup * pullback_base

    pullback_fct = create_mapper(pullback_base_class, variety_symbols)
    
    # Execute Augmentation
    sequence_depth = nesting_depth(sequences)
    
    if sequence_depth == 2
        return calculate_augmentation_on_surface(sequences, exceptional_divisor, pullback_fct)
    elseif sequence_depth == 3
        return calculate_all_augmentations_on_surface(sequences, exceptional_divisor, pullback_fct)
    end
end    