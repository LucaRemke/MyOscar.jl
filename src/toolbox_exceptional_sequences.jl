########################################
# DEFINITION OF (NEGATIVE) IMMACULATE LOCUS
########################################

#---------------------------------------
# Define a structure which stores the values for the immaculate locus
# The following bases are used:
# P^n - D1 
# H_r - D1,D3 
# Pentagon - D1,D3,D2+D4+D5
# Hexagon - D1,D2,D3,D4 
# P1xP1xP1 - D1,D3,D5
#***************************************
# The parameters are as follows:
# difference: check if difference is in (negative) immaculate locus
# param: used for varieties like P^n or H_r to specify n or reached
# sgn: specifies if negative or positive immaculate locus is used
#---------------------------------------
struct NegativeImmaculateLocus
    name::Symbol
    func::Function
end;

# Initialize dictionary
nimmloc = Dict{Symbol, NegativeImmaculateLocus}()

#nimmloc[:p2] = NegativeImmaculateLocus(:p2, 
#    (difference, param) -> all(diff -> diff == [1] || diff == [2], difference)
#) 

nimmloc[:pn] = NegativeImmaculateLocus(:pn, 
    (difference, param, sgn) -> all(diff -> diff in [[i] for i in 1:param], sign(sgn)*difference)
) 

nimmloc[:h0] = NegativeImmaculateLocus(:h0, 
    (difference, param, sgn) -> all(diff -> diff[1] == 1 || diff[2] == 1, sign(sgn)*difference)
)

nimmloc[:hr] = NegativeImmaculateLocus(:hr,
    (difference, param, sgn) -> all(diff -> diff == [1,0] || diff == [1-param,2] || diff[2] == 1, sign(sgn)*difference)
)

nimmloc[:twisted540122] = NegativeImmaculateLocus(:twisted540122,
    (difference, param, sgn) -> all(diff -> diff[2] == 1 || diff[2] == 2 || diff[2] == 3 || 
        diff == [1,0] || diff == [2,0] || diff == [3,0] || diff == [4,0] || diff == [3,-1] || diff == [4,-1] ||
        diff == [-1,4] || diff == [-2,4] || diff == [-3,4] || diff == [-4,4] || diff == [-3,5] || diff == [-4,5],
        sign(sgn)*difference)
)

nimmloc[:pentagon] = NegativeImmaculateLocus(:pentagon,
    (difference, param, sgn) -> all(diff -> diff == [1,1,-1] || diff == [0,0,2] || 
                (diff[3] == 0 && (diff[1] == 1 || diff[2] == 1)) ||
                (diff[3] == 1 && (diff[1] == 0 || diff[2] == 0)), sign(sgn)*difference)
) 

nimmloc[:hexagon] = NegativeImmaculateLocus(:hexagon,
    (difference, param, sgn) -> all(diff -> 
        diff == [0,0,-1,1] || 
        diff == [0,0,1,-1] ||  
        diff == [2,2,0,2] ||
        diff == [2,2,-2,4] ||
        (diff[1] == diff[4] && diff[2] == 1 && diff[3] == -1) ||
        (diff[1] == diff[4] && diff[2] == 1 && diff[3] == 0) ||
        (diff[1] == 1 && diff[2] == diff[4] && diff[3] == 0) ||
        (diff[1] == 1 && diff[2] == diff[4] && diff[3] == -1) ||
        #
        (diff[1]+1 == diff[4] && diff[2] == 1 && diff[3] == -1) ||
        (diff[1]+1 == diff[4] && diff[2] == 1 && diff[3] == 0) ||
        (diff[1] == diff[2] && diff[1] == -diff[3] && diff[1]+1 == diff[4]) ||
        (diff[1]+1 == diff[2] && diff[1] == -diff[3] && diff[1]+1 == diff[4]) ||
        #
        (diff[1] == diff[2] && diff[1] == -diff[3]+1 && diff[1] == diff[4]) ||
        (diff[1]-1 == diff[2] && diff[1] == -diff[3]+1 && diff[1] == diff[4]) ||
        (diff[1] == 1 && diff[3] == 0 && diff[2]+1 == diff[4]) ||
        (diff[1] == 1 && diff[3] == -1 && diff[2]+1 == diff[4]),
        sign(sgn)*difference)
) 

nimmloc[:p1p1p1] = NegativeImmaculateLocus(:p1p1p1,
    (difference, param, sgn) -> all(diff -> diff[1] == 1 || diff[2] == 1 || diff[3] == 1, sign(sgn)*difference)
) 


########################################
# DICTIONARY OF CONFIGURATIONS
########################################

#---------------------------------------
# Define a dictionary which stores the functions for the configurations
# The following bases are used:
# P^n - D1 
# H_r - D1,D3 
# Pentagon - D1,D3,D2+D4+D5
# Hexagon - D1,D2,D3,D4 
# P1xP1xP1 - D1,D3,D5
#***************************************
# The parameters are as follows:
# difference: check if difference is in (negative) immaculate locus
# param: used for varieties like P^n or H_r to specify n or reached
# sgn: specifies if negative or positive immaculate locus is used
#---------------------------------------


#---------------------------------------
# Configurations for P^1xP^1
#---------------------------------------
function configurate_exceptional_sequences_on_h0(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    type_sequences_dict = Dict{Vector{Int64}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int64}()
        for vec in seq
            i = -1*sign(sgn)*vec #we want the negative immaculate locus
            
            if i == [-1,-1]
                push!(type_seq, 0)
            elseif i[2] == -1
                push!(type_seq, 1)
            elseif i[1] == -1
                push!(type_seq, 2)
            end
        end
        if haskey(type_sequences_dict, type_seq)
            push!(type_sequences_dict[type_seq], seq)
        else
            type_sequences_dict[type_seq] = [seq]
        end
    end

    return type_sequences_dict
end;

#---------------------------------------
# Configurations for H_r
#---------------------------------------
function configurate_exceptional_sequences_on_hr(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    type_sequences_dict = Dict{Vector{Int64}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int64}()
        for vec in seq
            i = -1*sign(sgn)*vec #we want the negative immaculate locus
            
            if i == [-1,0]
                push!(type_seq, 0)
            elseif i == [param-1,-2]
                push!(type_seq, 1)
            elseif i[2] == -1
                push!(type_seq, 2)
            end
        end
        if haskey(type_sequences_dict, type_seq)
            push!(type_sequences_dict[type_seq], seq)
        else
            type_sequences_dict[type_seq] = [seq]
        end
    end

    return type_sequences_dict
end;

#---------------------------------------
# Configurations for the pentagon
#---------------------------------------
function configurate_exceptional_sequences_on_pentagon(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    type_sequences_dict = Dict{Vector{Int64}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int64}()
        for vec in seq
            i = -1*sign(sgn)*vec #we want the negative immaculate locus
            
            if i == [-1,-1,1]
                push!(type_seq, 0)
            elseif i == [0,0,-2]
                push!(type_seq, 1)
            elseif i == [-1,-1,0]
                push!(type_seq, 2)
            elseif i == [0,0,-1]
                push!(type_seq, 3)
            elseif i[3] == 0 && i[1] == -1
                push!(type_seq, 4)
            elseif i[3] == 0 && i[2] == -1
                push!(type_seq, 5)
            elseif i[3] == -1 && i[1] == 0
                push!(type_seq, 6)
            elseif i[3] == -1 && i[2] == 0
                push!(type_seq, 7)
            end
        end
        if haskey(type_sequences_dict, type_seq)
            push!(type_sequences_dict[type_seq], seq)
        else
            type_sequences_dict[type_seq] = [seq]
        end
    end

    return type_sequences_dict
end;

#---------------------------------------
# Configurations for the hexagon
#---------------------------------------
function configurate_exceptional_sequences_on_hexagon(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    type_sequences_dict = Dict{Vector{Int64}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int64}()
        for vec in seq
            i = -1*sign(sgn)*vec #we want the negative immaculate locus
            
            if i == [0,0,1,-1]
                push!(type_seq, 0)
            elseif i == [0,0,-1,1]
                push!(type_seq, 1)
            elseif i == [-2,-2,0,-2]
                push!(type_seq, 2)
            elseif i == [-2,-2,2,-4]
                push!(type_seq, 3)
            # The points with multiple codes get an own type
            elseif i == [-1, 0, 0, -1]
                push!(type_seq, -1)
            elseif i == [-1, -2, 1, -2]
                push!(type_seq, -2)
            elseif i == [-1, -1, 1, -1]
                push!(type_seq, -3)
            elseif i == [-1, -1, 0, -1]
                push!(type_seq, -4)
            elseif i == [0, -1, 0, -1]
                push!(type_seq, -5)
            elseif i == [-1, -1, 0, -2]
                push!(type_seq, -6)
            elseif i == [-2, -1, 1, -2]
                push!(type_seq, -7)
            elseif i == [-1, -1, 1, -2]
                push!(type_seq, -8)
            elseif (i[1] == i[4] && i[2] == -1 && i[3] == 1)
                push!(type_seq, 4)
            elseif (i[1] == i[4] && i[2] == -1 && i[3] == 0)
                push!(type_seq, 5)
            elseif (i[2] == i[4] && i[1] == -1 && i[3] == 0)
                push!(type_seq, 6)
            elseif (i[2] == i[4] && i[1] == -1 && i[3] == 1)
                push!(type_seq, 7)
            elseif (i[1] == i[2] && i[1] == -i[3]-1 && i[1] == i[4])
                push!(type_seq, 8)
            elseif (i[1] == i[2]-1 && i[1] == -i[3]-1 && i[1] == i[4])
                push!(type_seq, 9)
            elseif (i[1] == i[2]+1 && i[1] == -i[3] && i[1] == i[4]+1)
                push!(type_seq, 10)
            elseif (i[1] == i[2] && i[1] == -i[3] && i[1] == i[4]+1)
                push!(type_seq, 11)
            elseif (i[1] == i[4]+1 && i[2] == -1 && i[3] == 0)
                push!(type_seq, 12)
            elseif (i[1] == i[4]+1 && i[2] == -1 && i[3] == 1)
                push!(type_seq, 13)
            elseif (i[1] == -1 && i[2] == i[4]+1 && i[3] == 0)
                push!(type_seq, 14)
            elseif (i[1] == -1 && i[2] == i[4]+1 && i[3] == 1)
                push!(type_seq, 15)
            end
        end
        if haskey(type_sequences_dict, type_seq)
            push!(type_sequences_dict[type_seq], seq)
        else
            type_sequences_dict[type_seq] = [seq]
        end
    end

    return type_sequences_dict
end;

#---------------------------------------
# Configurations for the product P^1 x P^1 x P^1
#---------------------------------------
function configurate_exceptional_sequences_on_p1p1p1(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    type_sequences_dict = Dict{Vector{Int64}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int64}()
        for vec in seq
            i = -1*sign(sgn)*vec #we want the negative immaculate locus
            
            if i == [-1,-1,-1]
                push!(type_seq, 0)
            elseif (i[1]==-1 && i[2]==-1)
                push!(type_seq, 1)
            elseif (i[1]==-1 && i[3]==-1)
                push!(type_seq, 2)
            elseif (i[2]==-1 && i[3]==-1)
                push!(type_seq, 3)
            elseif i[1] == -1
                push!(type_seq, 4)
            elseif i[2] == -1
                push!(type_seq, 5)
            elseif i[3] == -1 
                push!(type_seq, 6)
            end
        end
        if haskey(type_sequences_dict, type_seq)
            push!(type_sequences_dict[type_seq], seq)
        else
            type_sequences_dict[type_seq] = [seq]
        end
    end

    return type_sequences_dict
end;

#---------------------------------------
# Definition of the dictionary
#---------------------------------------
configurations = Dict{Symbol, Function}()
configurations[:h0] = configurate_exceptional_sequences_on_h0
configurations[:hr] = configurate_exceptional_sequences_on_hr
configurations[:pentagon] = configurate_exceptional_sequences_on_pentagon
configurations[:hexagon] = configurate_exceptional_sequences_on_hexagon
configurations[:p1p1p1] = configurate_exceptional_sequences_on_p1p1p1

########################################
# BRUTEFORCE COMPUTATION OF MESES
########################################

#---------------------------------------
# Checks if the difference of v with each vector in V is in the
# negative immaculate locus
#---------------------------------------
function is_in_nimmloc(
        v::Vector{T}, 
        V::Vector{Vector{T}}, 
        var::Symbol, 
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    
    difference = [v - w for w in V]
    cond = nimmloc[var]
    return cond.func(difference, param, sgn)
end;

#---------------------------------------
# Extend a given exceptional sequence by one element if it fullfills the exceptionality condition
#---------------------------------------
function extend_exceptional_sequences(
        E::Vector{Vector{Vector{T}}}, 
        V::Vector{Vector{T}}, 
        var::Symbol, 
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    
    new_sequences = Vector{Vector{Vector{T}}}()
    
    for seq in E
        for vec in V
            # Check if the difference between vec and all vectors in seq is in P
            if is_in_nimmloc(vec, seq, var, param, sgn)
                push!(new_sequences, push!(copy(seq), vec)) 
            end
        end
    end
    
    return new_sequences
end;

#---------------------------------------
# Generate all exceptional sequences of length search_depth-1,
# where the first element is the trivial element in the Picard group
#---------------------------------------
function generate_exceptional_sequences(
        var::Symbol,
        search_range::Vector{UnitRange{T}},
        search_depth::Int;
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1
    ) where T
    
    l = length(search_range)
    
    if l == 1
        considered_points_old = generate_vectors_from_ranges(search_range)
        considered_points = [[t] for t in considered_points_old]
    else
        considered_points = generate_vectors_from_ranges(search_range)
    end
    
    # Starting with sequences of length 1
    # Since we can shift exceptional sequences we cn set the origin to zero
    starting_points = [zeros(Int, l)]
    initial_sequences = [[v] for v in starting_points]

    # Perform iterations to extend sequences
    for i in 1:search_depth
        println("We reached step $i") 
        initial_sequences = extend_exceptional_sequences(initial_sequences, considered_points, var, param, sgn)
        if initial_sequences == Vector{Vector{Int64}}[]
            println("There are no sequences after $i steps")
            break
        else
            println("There are $(length(initial_sequences)) sequences of length $(i+1) \n")
        end
    end

    return initial_sequences
end; 

########################################
# CONDITIONS FOR (MAXIMAL) EXCEPTIONAL SEQUENCES
########################################

#---------------------------------------
# Check if a given sequence is exceptional (possible for all varieties in nimmloc)
#---------------------------------------
function is_exceptional(
        var::Symbol,
        V::Vector{Vector{T}};
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1,
        M::Matrix{Int64} = Matrix{Int}(I,length(V[1]),length(V[1]))
    ) where T

    # base change via M
    V = [M*v for v in V]

    l = length(V)

    for i in 2:l
        V_sub = V[1:i-1]
        V_i = V[i]
        if is_in_nimmloc(V_i, V_sub, var, param, sgn) == false
            return false
            break
        end
    end

    return true
end;


########################################
# AUGMENTATION + HELIXING
########################################


# Calculate the augmentions of a maximal exceptional sequence on a toric surface
# mes - the maximal exceptional sequence whose augmentations should be calculated
# E - the exceptional divisor of the blow up
# p - the pullback of divisors
#---------------------------------------
function calculate_augmentation_on_surface(mes::Vector{Vector{T}}, E::Vector{Int64}, pullback) where T
    
    augmentations = Vector{Vector{T}}[]
    l = length(mes)
    
    pmes = pullback.(mes)
    
    for i in 1:l
        aug = copy(pmes)
        insert!(aug, i, aug[i])
               
        add_excdiv = [j for j in 1:i-1 if j != i]
        push!(add_excdiv, i+1)
        
        for k in add_excdiv
            aug[k] = aug[k] + E
        end
        
        aug = [vec.-aug[1] for vec in aug]
        push!(augmentations, aug)
    end
    
    return augmentations
end;

#---------------------------------------
# Calculate the helixing of a maximal exceptional sequence
# mes - the maximal exceptional sequence that is helexed
# K - the anticanonical divisor -K_X in Pic(X) of a toric variety X 
#---------------------------------------
function calculate_helixing(mes::Vector{Vector{T}}, K::Vector{Int64}) where T
    l = length(mes)
    helexed_mes = mes[2:l]
    push!(helexed_mes, mes[1]+K)
    helexed_mes = [vec.-helexed_mes[1] for vec in helexed_mes]
    
    return helexed_mes
end;

#---------------------------------------
# Calculate all helixings of a maximal exceptional sequence
# mes - the maximal exceptional sequence that is helexed
# K - the anticanonical divisor -K_X in Pic(X) of a toric variety X 
#---------------------------------------
function calculate_all_helixing(mes::Vector{Vector{T}}, K::Vector{Int64}) where T
    period = length(mes)  
    
    all_helexed_mes = Vector{Vector{T}}[]
    current_mes = mes
    push!(all_helexed_mes, current_mes)    
    
    for i in 1:period-1
        current_helex = calculate_helixing(current_mes, K)
        push!(all_helexed_mes, current_helex)
        current_mes = current_helex
    end
    
    return all_helexed_mes
    
end;

#---------------------------------------
# Calculate the dualized sequence of a mes
# mes - the maximal exceptional sequence that is helexed
# K - the canonical divisor K_X in Pic(X) of a toric variety X 
#---------------------------------------

function calculate_dualizing(mes::Vector{Vector{T}}, K::Vector{Int64}) where T
    l = length(mes)
    
    inverted_mes = [-K .- mes[l+1-i] for i in 1:l]
    inverted_mes = [vec.-inverted_mes[1] for vec in inverted_mes]
    
    return inverted_mes
end;


########################################
# FLIPPING
########################################

#---------------------------------------
# Calculate the flipping for a given permutation of a maximal exceptional sequence on a toric surface
#---------------------------------------
function calculate_flipping(V::Vector{Vector{T}}, perm::Vector{Int64}) where T
    
    flipping = Vector{T}[]
    
    for v in V
        push!(flipping, v[perm])
    end
  
    return flipping
end;


########################################
# FULLNESS CONDITIONS
########################################

#---------------------------------------
# Checks if ph_ac(E1,...,En) > -dim X, compare Krah
#---------------------------------------
function is_not_full_by_krah(seq::Vector{Vector{Int64}}, variety::NormalToricVariety)
    
    l = length(seq)    
    if l != n_maximal_cones(polyhedral_fan(variety))
        throw(ErrorException("The sequence is not maximal and therefore not full")) 
    end
    
    for i in 2:l
        seq_sub = seq[1:i-1]
        seq_i = seq[i]
        diff = [seq_i - sub for sub in seq_sub]
        for d in diff
            tdc = toric_divisor_class(variety, d)
            tl = toric_line_bundle(toric_divisor(tdc)) 
            h0 = cohomology(tl, 0)
            if h0 != 0
                return false
            end
        end
    end
    
    return true
end;

#---------------------------------------
# Calculates the anticanonical pseudoheight
#---------------------------------------
function calculate_anticanonical_pseudoheight(seq::Vector{Vector{Int64}}, variety::NormalToricVariety)
    
    l = length(seq)    
    if l != n_maximal_cones(polyhedral_fan(variety))
        throw(ErrorException("The sequence is not maximal and therefore not full")) 
    end
    
    # calculate all possible e(Ei,Ej)
    relative_height = Dict{Vector{Int64}, Float64}()
    for i in 1:l
        for j in i+1:l
            seq_i, seq_j = seq[i], seq[j]
            diff_ij = seq_j - seq_i
            tdc = toric_divisor_class(variety, diff_ij)
            tl = toric_line_bundle(toric_divisor(tdc)) 
            cohom = all_cohomologies(tl)
            k = findfirst(x -> x != 0, cohom)
            if isnothing(k)
                k = Inf
            end
            
            relative_height[[i,j]] = k-1
        end
    end
    
    # calculate all possile values for which the infimum is build
    comb_all = filter(x -> length(x)>1, collect(combinations(1:l))) 
    
    # calculate the divisor class of the anticanonical line bundle
    acdc = divisor_class(anticanonical_divisor_class(variety)) 
    var_rank = torsion_free_rank(class_group(variety))   
    acdc_int = [Int64(acdc[i]) for i in 1:var_rank]
        
    # calculate the minimum
    # for length 0 we have e(Ei,Ei+wX^-1)=min H^l(wX^-1)
    wx_cohom = all_cohomologies(anticanonical_bundle(variety)) 
    wx_k = findfirst(x -> x != 0, wx_cohom)
    if isnothing(wx_k)
        wx_k = Inf
    end
    ph_ac = wx_k-1
    
    for comb in comb_all
        l_comb = length(comb)
        max_comb, min_comb = maximum(comb), minimum(comb)
        td_wx = toric_divisor_class(variety, -seq[max_comb] + seq[min_comb] + acdc_int)
        tl_wx = toric_line_bundle(toric_divisor(td_wx)) 
        cohom_wx = all_cohomologies(tl_wx)
        k_wx = findfirst(x -> x != 0, cohom_wx)
        if isnothing(k_wx)
            k_wx = Inf
        end
        ph_ac_comb = k_wx-1
        
        sub_comb = [[comb[a], comb[b]] for a in 1:l_comb for b in a+1:l_comb]
        for s in sub_comb
            val_s = relative_height[s]
            ph_ac_comb += val_s
        end
        
        ph_ac_comb += -l_comb
        
        if ph_ac_comb < ph_ac
            ph_ac = ph_ac_comb
        end
    end
    
    ph_ac = Int64(ph_ac)
    
    return ph_ac
end;













# old functions
#=
# function that checks if there are duplicates in vector of vectors of vectors
function has_duplicates(list::Vector{Vector{Vector{T}}}; print::Bool=false) where T
    set = Set(list)

    if length(list) == length(set)
        return false
    else
        print == true ? println("Duplicates found at indices:") : nothing
        for (i, vec) in enumerate(list)
            if vec in set
                delete!(set, vec)
            else
                print == true ? println(i) : nothing
            end
        end
        return true
    end
end

# function for calculating the inverting sequence of a mes
# mes - the maximal exceptional sequence that is helexed
# K - the anticanonical divisor -K_X in Pic(X) of a toric variety X 
function calculate_inverted_mes(mes::Vector{Vector{T}}, K::Vector{Int64}) where T
    l = length(mes)
    
    inverted_mes = [i==1 ? mes[i] : K .- mes[l+2-i] for i in 1:l]
    
    return inverted_mes
end

# function that looks for duplicates in the two inputs
function search_duplicates(list_source::Vector{Vector{Vector{T}}}, list_target::Vector{Vector{Vector{T}}}) where T
    l = length(list_source[1])
    
    # this makes the function much faster!
    source_dict = Dict(v => i for (i, v) in enumerate(list_source))
    matches = Vector{Int64}[]
    
    for(j, seq_j) in enumerate(list_target)
        if haskey(source_dict, seq_j)
            i = source_dict[seq_j]

            push!(matches, [i,j])
        end
    end

    
    return matches
    
end

# function for calculating the augmentions of the meses on the pentagon
function create_augmentations(meses::Vector{Vector{Vector{Sym{PyObject}}}}, E::Vector{Int64})
    augmentations = Vector{Vector{Sym{PyObject}}}[]
    for seq in meses
        mseq = [[vec[1], vec[1]+vec[2]+vec[3], -vec[1], vec[1]+vec[3]] for vec in seq]
        
        aug_1 = [mseq[1], mseq[1]+E, mseq[2], mseq[3], mseq[4], mseq[5]]
        aug_2 = [mseq[1], mseq[2]-E, mseq[2], mseq[3]-E, mseq[4]-E, mseq[5]-E]
        aug_3 = [mseq[1], mseq[2], mseq[3]-E, mseq[3], mseq[4]-E, mseq[5]-E]
        aug_4 = [mseq[1], mseq[2], mseq[3], mseq[4]-E, mseq[4], mseq[5]-E]
        aug_5 = [mseq[1], mseq[2], mseq[3], mseq[4], mseq[5]-E, mseq[5]]
        push!(augmentations, aug_1)
        push!(augmentations, aug_2)
        push!(augmentations, aug_3)
        push!(augmentations, aug_4)
        push!(augmentations, aug_5)
    end
    
    return augmentations
end


# function for calculating the helexed meses on the hexagon
function calculate_helexing(meses::Vector{Vector{Vector{Sym{PyObject}}}}, K::Vector{Int64})
    l = length(meses[1])  
    
    seq_helexed = Vector{Any}()
    
    for seq in meses
        helexed_part = [vec.+K for vec in seq]
        
        helexed_seq = vcat(seq, helexed_part)
        
        all_helexed_seq = Any[]
        for i in 1:l-1        
            helexed_subseq = helexed_seq[1+i:6+i]
            helexed_subseq_zero = [vec.-helexed_subseq[1] for vec in helexed_subseq]
            push!(all_helexed_seq, helexed_subseq_zero)
        end
        push!(seq_helexed, [seq, all_helexed_seq])
    end
    
    return seq_helexed
    
end
# calculate for each of the 100 augmentations their helexed sequences
# the period is 6, that is after 6 steps we obtain the starting sequence.
# therefore each of the 100 augmentations gives 5 new meses leading to
# 100 + 5*100 = 600 meses
K = [2,2,-1,3]
A = augmentations_on_hexagon
helexings = calculate_helexing(A, K);

# function for calculating the inverting sequence of a mes on the hexagon
function inverting_sequence(meses::Vector{Vector{Vector{Sym{PyObject}}}}, K::Vector{Int64})
    inverted_meses = Vector{Vector{Sym{PyObject}}}[]
    
    for seq in meses
        new_seq = [seq[1], K-seq[6], K-seq[5], K-seq[4], K-seq[3], K-seq[2]]
    
        push!(inverted_meses, new_seq)      
    end
    
    return inverted_meses
end

# calculating the inverted sequences for the 100 augmentations
# leading to 100 inverted sequences
K = [2,2,-1,3]
A = augmentations_on_hexagon
invertings = inverting_sequence(A, K)

# write all aug and their helexes in one vector (of length 600)
# preparation for the next step
helexing_vec = Vector{Vector{Sym{PyObject}}}[]
for hex in helexings
    push!(helexing_vec, hex[1])
    for i in 1:5
        push!(helexing_vec, hex[2][i])
    end
end
helexing_vec

# calculate the invertings for all 600 meses and substitute
K = [2,2,-1,3]
A = helexing_vec
invertings_all = inverting_sequence(A, K)

invertings_sub_all = Vector{Vector{Sym{PyObject}}}[]

for inv in invertings_all
    inv_sub = [[subs(expr, n, -n) for expr in vec] for vec in inv]
    push!(invertings_sub_all, inv_sub)
end

invertings_sub_all;

matches = 5
coincidences_all = Any[]

for inv in invertings_sub_all
    coincidences_inv = Any[]
    
    for (index, seq) in enumerate(helexing_vec)
        count = 0
        for i in 2:length(seq) # the zero vector always coincides         
            if inv[i] == seq[i]
                count += 1
            end
        end
        
        if count >= matches
            quotient, remainder = divrem(index, 6)
            if remainder == 0
                upper_index = 5
                lower_index = quotient
            else
                upper_index = remainder - 1
                lower_index = quotient + 1
            end

            push!(coincidences_inv, (seq, [lower_index, upper_index])) 
        end
    end
    
    println(length(coincidences_inv))    
    push!(coincidences_all, coincidences_inv)
end
coincidences_all;

# compare the results
for (index, c) in enumerate(coincidences_all)
    quotient, remainder = divrem(index, 6)
    if remainder == 0
        upper_index = 5
        lower_index = quotient
    else
        upper_index = remainder - 1
        lower_index = quotient + 1
    end    
    
    println([lower_index, upper_index], "\t origin sequence")
    
    for i in 1:length(c)
        println(c[i][2], "\t corresponding inverted sequence")
    end
    
    println("---------------------------------------------------------------")
end

matches = 5
coincidences_all_new = Any[]

for inv in invertings_all
    coincidences_inv = Any[]
    
    for (index, seq) in enumerate(helexing_vec)
        count = 0
        for i in 2:length(seq) # the zero vector always coincides         
            if inv[i] == seq[i]
                count += 1
            end
        end
        
        if count >= matches
            quotient, remainder = divrem(index, 6)
            if remainder == 0
                upper_index = 5
                lower_index = quotient
            else
                upper_index = remainder - 1
                lower_index = quotient + 1
            end

            push!(coincidences_inv, (seq, [lower_index, upper_index])) 
        end
    end
    
    println(length(coincidences_inv))    
    push!(coincidences_all_new, coincidences_inv)
end
coincidences_all_new;

# compare the results
for (index, c) in enumerate(coincidences_all_new)
    quotient, remainder = divrem(index, 6)
    if remainder == 0
        upper_index = 5
        lower_index = quotient
    else
        upper_index = remainder - 1
        lower_index = quotient + 1
    end    
    
    println([lower_index, upper_index], "\t origin sequence")
    
    for i in 1:length(c)
        println(c[i][2], "\t corresponding inverted sequence")
    end
    
    println("---------------------------------------------------------------")
end

=#