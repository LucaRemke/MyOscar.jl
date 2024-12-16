# This is for all functions using sympy
# This file will be replaced by toolbox_sympy.jl if the push on github works

########################################
# load required packages
########################################
using PyCall
pyimport("sympy")
using SymPy;


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
        ind = Int64[]

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
function evaluate_cartier_data(data::Vector{Any}, symbols::Vector{Sym}, coeffs::Vector{Int64})
    
    if length(symbols) != length(coeffs)
        throw(ArgumentError("Length of the input vectors does not match.")) 
    end
    
    all_data = []
    
    for i in 1:length(data)
        cdata = data[i][2]
        
        for j in 1:length(symbols)
            cdata = cdata.subs(symbols[j],coeffs[j])
        end
        
        cdata = vec(convert(Matrix{Int64}, N(cdata))) 
        
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
function unify_sequence(V::Vector{Vector{Sym{PyObject}}})
    for (i, vector) in enumerate(V)
        for (j, element) in enumerate(vector)
            # Try to convert the current element to an Integer
            convert = N(element)
            # The first time the convertion gives not an integer is the desired value
            if typeof(convert) == Sym{PyObject}
                # depending on the signs, the substitution differs
                sgn_p = N(convert.coeff(p)) 
                if sgn_p == -1
                    subst = convert                
                else
                    subst = p + (-1*N(convert)(p=>0))
                end

                seq = [[subs(expr, p, subst) for expr in vec] for vec in V]

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
        cd::Vector{Int64}
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
        cd::Vector{Int64}
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
        s::Int64=1
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
    return "No parameter was found"
end;


#---------------------------------------
# Substitute a sequence with a parameter with the value
#---------------------------------------
function substitute_sequence(
        seq::Vector{Vector{Sym{PyObject}}},
        param::Sym{PyObject},
        value::Int64
    )
    seq_sub = [[subs(expr, param, value) for expr in vec] for vec in seq]
    seq_sub = transform_rayvector.(seq_sub)
    return seq_sub
end;

#---------------------------------------
# Substitute a sequence with a parameter with the value at the same index of the integer sequence
#---------------------------------------
function substitute_configuration(
        seq_int::Vector{Vector{Int64}}, 
        seq_p::Vector{Vector{Sym{PyObject}}}
    )
    vec, elem = find_parameter(seq_p)
    subst_var = seq_int[vec][elem]
    subst_seq_p = substitute_sequence(seq_p, p, subst_var)
    
    return subst_seq_p
end;

#---------------------------------------
# Look if a concrete sequence corresponds to a sequence with general description
#---------------------------------------
function find_sequence_in_configurations(
        seqs_int::Vector{Vector{Int64}},
        seqs_p::Vector{Vector{Vector{Sym{PyObject}}}};
        return_sequence::Bool=false
    )   

    for seq_p in seqs_p
        seq_p_subst = substitute_configuration(seqs_int, seq_p)
        if seq_p_subst == seqs_int
            if return_sequence == true
                return seq_p
            else
                return true
            end
        end
    end

    return false
end;

#---------------------------------------
# Checks if a vector of sequences can be described by a general rule
#---------------------------------------
function find_configurations(
        D::Dict{Vector{Int64}, Vector{Vector{Vector{Int64}}}};
        check_result::Bool=true
    )
    general_description = Dict{Vector{Int64}, Vector{Vector{Sym{PyObject}}}}()
    p = Sym("p") 
    
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
                            push!(conf_seq_general_vec, p)
                        else
                            if sign(parameter) == sign(conf_seq_1_elem)
                                diff_elem = conf_seq_1_elem - parameter
                                vec_elem = p + diff_elem
                            else
                                diff_elem = conf_seq_1_elem + parameter
                                vec_elem = -p + diff_elem
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
        D_general::Dict{Vector{Int64}, Vector{Vector{Sym{PyObject}}}},
        D_single::Dict{Vector{Int64}, Vector{Vector{Vector{Int64}}}}
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
        V_single::Vector{Vector{Int64}}
    )
    found = find_sequence_in_configurations(V_single, [V_general])
    if found == false
        return false
    else
        return true
    end
end;
