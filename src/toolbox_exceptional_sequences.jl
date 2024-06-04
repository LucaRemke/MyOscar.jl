# function for calculating the augmentions of a maximal exceptional sequence
# mes - the maximal exceptional sequence whose augmentations should be calculated
# E - the exceptional divisor of the blow up
# p - the pullback of divisors
function calculate_augmentations_of_mes(mes::Vector{Vector{T}}, E::Vector{Int64}, pullback) where T
    
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
end

# function for calculating the helexing of a maximal exceptional sequence
# mes - the maximal exceptional sequence that is helexed
# K - the anticanonical divisor -K_X in Pic(X) of a toric variety X 
function calculate_helex_of_mes(mes::Vector{Vector{T}}, K::Vector{Int64}) where T
    l = length(mes)
    helexed_mes = mes[2:l]
    push!(helexed_mes, mes[1]+K)
    helexed_mes = [vec.-helexed_mes[1] for vec in helexed_mes]
    
    return helexed_mes
    
end


# function for calculating all helexings of a maximal exceptional sequence
# mes - the maximal exceptional sequence that is helexed
# K - the anticanonical divisor -K_X in Pic(X) of a toric variety X 
function calculate_all_helex_of_mes(mes::Vector{Vector{T}}, K::Vector{Int64}) where T
    period = length(mes)  
    
    all_helexed_mes = Vector{Vector{T}}[]
    current_mes = mes
    push!(all_helexed_mes, current_mes)    
    
    for i in 1:period-1
        current_helex = calculate_helex_of_mes(current_mes, K)
        push!(all_helexed_mes, current_helex)
        current_mes = current_helex
    end
    
    return all_helexed_mes
    
end

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





# old functions
#=
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