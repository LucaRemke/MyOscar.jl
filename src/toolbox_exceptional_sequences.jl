########################################
# IMPLEMENTED TORIC VARIETIES AND THEIR BASES
#---------------------------------------
# The projective space P^n (D1) 
# The Hirzebruch surface H_r (D1,D3) 
# The Pentagon (D1,D3,D2+D4+D5)
# The Hexagon (D1,D2,D3,D4) 
# The product P^1xP^1xP^1 (D1,D3,D5)
########################################


########################################
# STRUCTURE ENCODING THE (NEGATIVE) IMMACULATE LOCI OF CERTAIN TORIC VARIETIES
########################################


"""
    NegativeImmaculateLocus

Structure to store the classes of line bundles in the negative immaculate locus, consisting of:

- `name::Symbol`: a symbolic identifier for the toric variety
- `func::Function`: a function that checks if a difference of two classes is in the immaculate locus

The function has signature `(difference, param, sgn) -> Bool` with:

- `difference`: The difference of two vectors
- `param`: A parameter for varieties like P^n or H_r
- `sgn`: Specifies if the negative or positive immaculate locus is of interest

Available toric varieties with bases:
- The projective space P^n (D1) 
- The Hirzebruch surface H_r (D1,D3) 
- The Pentagon (D1,D3,D2+D4+D5)
- The Hexagon (D1,D2,D3,D4) 
- The product P^1xP^1xP^1 (D1,D3,D5)
"""
struct NegativeImmaculateLocus
    name::Symbol
    func::Function
end

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


@doc raw"""
    configurate_exceptional_sequences_on_h0(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Associate to each sequences a labeling by assign a special value to each element of the sequence.

# Arguments
- `V::Vector{Vector{Vector{T}}}`: The sequences which should be encoded
- `param::Union{Int, Nothing}`: (optional) Not required here
- `sgn::Int`: (optional) Determines if the positive or negative immaculate locus should be considered
"""
function configurate_exceptional_sequences_on_h0(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    ) where T
    type_sequences_dict = Dict{Vector{Int}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int}()
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
end


@doc raw"""
    configurate_exceptional_sequences_on_hr(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Associate to each sequences a labeling by assign a special value to each element of the sequence.

# Arguments
- `V::Vector{Vector{Vector{T}}}`: The sequences which should be encoded
- `param::Union{Int, Nothing}`: (optional) Defines r for H_r
- `sgn::Int`: (optional) Determines if the positive or negative immaculate locus should be considered
"""
function configurate_exceptional_sequences_on_hr(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    ) where T
    type_sequences_dict = Dict{Vector{Int}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int}()
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
end


@doc raw"""
    configurate_exceptional_sequences_on_pentagon(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Associate to each sequences a labeling by assign a special value to each element of the sequence.

# Arguments
- `V::Vector{Vector{Vector{T}}}`: The sequences which should be encoded
- `param::Union{Int, Nothing}`: (optional) Not required here
- `sgn::Int`: (optional) Determines if the positive or negative immaculate locus should be considered
"""
function configurate_exceptional_sequences_on_pentagon(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    ) where T
    type_sequences_dict = Dict{Vector{Int}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int}()
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
end


@doc raw"""
    configurate_exceptional_sequences_on_hexagon(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Associate to each sequences a labeling by assign a special value to each element of the sequence.

# Arguments
- `V::Vector{Vector{Vector{T}}}`: The sequences which should be encoded
- `param::Union{Int, Nothing}`: (optional) Not required here
- `sgn::Int`: (optional) Determines if the positive or negative immaculate locus should be considered
"""
function configurate_exceptional_sequences_on_hexagon(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    ) where T
    type_sequences_dict = Dict{Vector{Int}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int}()
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
end


@doc raw"""
    configurate_exceptional_sequences_on_pentagon(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Associate to each sequences a labeling by assign a special value to each element of the sequence.

# Arguments
- `V::Vector{Vector{Vector{T}}}`: The sequences which should be encoded
- `param::Union{Int, Nothing}`: (optional) Not required here
- `sgn::Int`: (optional) Determines if the positive or negative immaculate locus should be considered
"""
function configurate_exceptional_sequences_on_p1p1p1(
        V::Vector{Vector{Vector{T}}},
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    ) where T
    type_sequences_dict = Dict{Vector{Int}, Vector{Vector{Vector{T}}}}()

    for seq in V
        type_seq = Vector{Int}()
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
end


# Definition of the dictionary
configurations = Dict{Symbol, Function}()
configurations[:h0] = configurate_exceptional_sequences_on_h0
configurations[:hr] = configurate_exceptional_sequences_on_hr
configurations[:pentagon] = configurate_exceptional_sequences_on_pentagon
configurations[:hexagon] = configurate_exceptional_sequences_on_hexagon
configurations[:p1p1p1] = configurate_exceptional_sequences_on_p1p1p1


########################################
# BRUTEFORCE COMPUTATION OF MAXIMAL EXCEPTIONAL SEQUENCES
########################################


@doc raw"""
    is_in_nimmloc(
        v::Vector{T}, 
        V::Vector{Vector{T}}, 
        var::Symbol, 
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Check if the difference of v and each vector in V is in the (negative) immaculate locus of the given toric variety.

# Arguments
- `v::Vector{T}`: The vector which should be checked
- `V::Vector{Vector{T}}`: The vectors to build the differences 
- `var::Symbol`: Symbol for the toric variety of interest 
- `param::Union{Int, Nothing}`: A parameter mandatory for P^n and H_r 
- `sgn::Int`: Determines if the neagative of positive immaculate locus should be considered

# Examples
```jldoctest
julia> is_in_nimmloc([1,1], [[0, 0], [1, 0], [0, 1]], :h0)
true

julia> is_in_nimmloc([-1,-1], [[0, 0], [1, 0], [0, 1]], :h0)
false


```
"""
function is_in_nimmloc(
        v::Vector{T}, 
        V::Vector{Vector{T}}, 
        var::Symbol, 
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    ) where T
    
    difference = [v - w for w in V]
    cond = nimmloc[var]
    return cond.func(difference, param, sgn)
end


@doc raw"""
    extend_exceptional_sequences(
        E::Vector{Vector{Vector{T}}}, 
        V::Vector{Vector{T}}, 
        var::Symbol, 
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Extend all given exceptional sequences E by elements in V if the exceptionality condition is fulfilled.
All sequences in E which can be extended by an element in V are returned.

# Arguments
- `E::Vector{Vector{Vector{T}}}`: The vector of exceptional sequences which should be checked
- `V::Vector{Vector{T}}`: The vectors to build the differences 
- `var::Symbol`: Symbol for the toric variety of interest 
- `param::Union{Int, Nothing}`: A parameter mandatory for P^n and H_r 
- `sgn::Int`: Determines if the neagative of positive immaculate locus should be considered

# Examples
```jldoctest
julia> E = [[[0, 0], [1, 0], [0, 1]], [[0, 0], [1, 1], [2, 1]]];

julia> V = [[1,1], [-1,-1], [1,2], [-1,-2]];

julia> extend_exceptional_sequences(E,V,:h0)
2-element Vector{Vector{Vector{Int64}}}:
 [[0, 0], [1, 0], [0, 1], [1, 1]]
 [[0, 0], [1, 1], [2, 1], [1, 2]]


```
"""
function extend_exceptional_sequences(
        E::Vector{Vector{Vector{T}}}, 
        V::Vector{Vector{T}}, 
        var::Symbol, 
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
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
end


@doc raw"""
    generate_exceptional_sequences(
        var::Symbol,
        search_range::Vector{UnitRange{T}},
        search_depth::Int;
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
    )

Generate all exceptional sequences of length `search_depth-1`, where the first element is the trivial element in the Picard group.

# Arguments
- `var::Symbol`: Symbol for the toric variety of interest
- `search_range::Vector{UnitRange{T}}`: Determines the used coefficients in the Picard group
- `search_depth::Int`: Determines the length of the sequences
- `param::Union{Int, Nothing}`: A parameter mandatory for P^n and H_r 
- `sgn::Int`: Determines if the neagative of positive immaculate locus should be considered

# Examples
```jldoctest
julia> generate_exceptional_sequences(:h0, [-1:1, -1:1], 3)
We reached step 1
There are 5 sequences of length 2

We reached step 2
There are 8 sequences of length 3

We reached step 3
There are 4 sequences of length 4

4-element Vector{Vector{Vector{Int64}}}:
 [[0, 0], [0, 1], [1, -1], [1, 0]]
 ⋮
 [[0, 0], [1, 0], [0, 1], [1, 1]]

```
"""
function generate_exceptional_sequences(
        var::Symbol,
        search_range::Vector{UnitRange{T}},
        search_depth::Int;
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1
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
end


@doc raw"""
    is_exceptional(
        var::Symbol,
        V::Vector{Vector{T}};
        param::Union{Int, Nothing}=nothing,
        sgn::Int64=1,
        M::Matrix{Int64} = Matrix{Int}(I,length(V[1]),length(V[1]))
    )

Check if a given sequence is exceptional (possible for all varieties in nimmloc).

# Arguments
- `var::Symbol`: Symbol for the toric variety of interest
- `V::Vector{Vector{T}}`: The sequence which should be checked for exceptionality
- `param::Union{Int, Nothing}`: A parameter mandatory for P^n and H_r 
- `sgn::Int`: Determines if the neagative of positive immaculate locus should be considered
- `M::Matrix{Int64}`: Applies a base change of the class group

# Examples
```jldoctest
julia> is_exceptional(:h0, [[0, 0], [1, 0], [0, 1], [1, 1]])
true

julia> is_exceptional(:h0, [[0, 0], [1, 0], [0, 1], [2, 1]])
false


```
"""
function is_exceptional(
        var::Symbol,
        V::Vector{Vector{T}};
        param::Union{Int, Nothing}=nothing,
        sgn::Int=1,
        M::Matrix{Int} = Matrix{Int}(I,length(V[1]),length(V[1]))
    ) where T

    # base change via M
    V = [M*v for v in V]
 
    for i in 2:length(V)
        if !is_in_nimmloc(V[i], V[1:i-1], var, param, sgn)
            return false
        end
    end

    return true
end

########################################
# # The above algorithm is fast but unflexible
# # This is the general approach
# range = -100:100
# range_vec = [range, range]
# h1 = define_hirzebruch_surface(1)
# h1_all, h1_coimm = co_immaculate_points_in_range(h1, range_vec)
# is_difference_coimmaculate(h1_coimm, [1,1], [0,0])
# seq = [[0, 0], [0, 1], [1, 1], [0, 2]]
# is_exceptional(h1_coimm, seq)
# seq = [[0, 0], [0, 1], [1, 1]]
# pt = [0, 2]
# is_exceptional_extendable(h1_coimm, seq, pt)
# variety = define_hirzebruch_surface(1)
# range = -100:100
# range_vec = [range, range]
# depth = 3
# generate_exceptional_sequences(variety, range_vec)
########################################

function get_picard_points(ranges::Vector{UnitRange{Int}})::Set{Vector{Int}}
    return Set(collect(coords) for coords in Iterators.product(ranges...))
end


function get_hypercube_vertices(ranges::Vector{UnitRange{Int}})::Set{Vector{Int}}
    
    # 1. Concise creation of boundary sets (replacing the first loop)
    # The map function applies the inner logic to every element (range) of the 'ranges' vector.
    # Logic: For each range R, create a Set containing only its first (min) and last (max) elements.
    boundaries = map(ranges) do R
        # Use an array literal [first(R), last(R)] and convert to a Set
        # If the range is single-point (e.g., 0:0), the Set will correctly contain only one element.
        Set([first(R), last(R)])
    end
    
    # 2. Use Iterators.product to generate all 2^d combinations
    all_vertex_coords = Iterators.product(boundaries...)
    
    # 3. Efficient collection into a Set (replacing the second loop)
    # Use a Set comprehension with a generator, similar to the optimization for Picard points.
    return Set(collect(coords) for coords in all_vertex_coords)
end

function vertices_to_matrix(vertices::Set{Vector{Int}})::Matrix{Int}
    
    x = length(vertices) # Number of rows (number of vertices)
    
    if x == 0
        return Matrix{Int}(undef, 0, 0) # Handle empty case
    end
    
    # 1. Determine dimension d based on the first element (Fast, since Set is not ordered)
    d = length(first(vertices)) 
    
    # 2. Initialize the Matrix (The key performance gain is here: one allocation)
    M = Matrix{Int}(undef, x, d)
    
    # 3. Populate the Matrix Row by Row (The highly efficient memory write loop)
    i = 1
    for vertex in vertices
        # Assign the vertex vector as the i-th row of the matrix M
        # This is a highly optimized memory copy in Julia's system.
        M[i, :] = vertex
        i += 1
    end
    
    return M
end

function maculate_points_in_range(immloc::ToricVanishingSet, range_vec::Vector{UnitRange{Int}})
    cube_vertices = get_hypercube_vertices(range_vec)
    matrix_vertices = vertices_to_matrix(cube_vertices)
    poly_picpoints = convex_hull(matrix_vertices)
    
    imm_poly = polyhedra(immloc)
    maculate_points = Set{Vector{Int}}()

    for poly in imm_poly
        intersection = intersect(poly, poly_picpoints)
        lp_intersection = collect(lattice_points(intersection))
        lp_intersection_t = transform_rayvector.(lp_intersection)
        union!(maculate_points, lp_intersection_t)
    end

    return maculate_points
end    


function co_immaculate_points_in_range(
        variety::NormalToricVariety, 
        range_vec::Vector{UnitRange{Int}}
    )
    
    # Considered points
    pic_points = get_picard_points(range_vec)
    
    # Immaculate locus
    imm_loc = immaculate_line_bundles(variety)
    imm_poly = polyhedra(imm_loc)
    
    # Maculate points in range
    mac_points = maculate_points_in_range(imm_loc, range_vec)
    
    # Immaculate points
    imm_points = setdiff(pic_points, mac_points)
    
    # Coimmaculate points
    coimm_points = -1 .* imm_points

    return (pic_points, Set(coimm_points))
end


function is_difference_coimmaculate(
        coimmaculate_set::Set{Vector{Int}},
        sj::Vector{T}, 
        si::Vector{T}
    )::Bool where T  
    # Calculate the difference vector, which represents the class of L_i tensor L_j^{-1}
    difference_vector = sj - si
    
    # Since 'immaculate_set' is a Set, membership checking (in) is highly efficient 
    # (average O(1) time complexity).
    return difference_vector in coimmaculate_set
end

function is_exceptional(
    coimmaculate_set::Set{Vector{Int}}, 
    seq::Vector{Vector{T}}
)::Bool where T
    
    l = length(seq)
    
    if l < 2
        return true
    end
        
    return all(2:l) do i
        all(1:i-1) do j
            is_difference_coimmaculate(coimmaculate_set, seq[i], seq[j])
        end
    end
end

function is_exceptional_extendable(
        coimmaculate_set::Set{Vector{Int}}, 
        seq::Vector{Vector{T}}, 
        pt_new::Vector{T}
    ) where T
    
    return all(seq) do v_j
        # Check the difference defined by the new pair (p_new, v_j)
        is_difference_coimmaculate(coimmaculate_set, pt_new, v_j)
    end
end

function generate_exceptional_sequences(
        variety::NormalToricVariety,
        coefficient_range::Vector{UnitRange{Int}},
        search_depth::Int=n_maximal_cones(polyhedral_fan(variety))-1;
        print_status::Bool=false
    )
    
    # Determine all considered points and coimmaculate points
    considered_points, coimm_points = co_immaculate_points_in_range(variety, coefficient_range)
    
    # Initialize the sequence with only one single element which is the zero in Pic(variety)
    cl_rank = torsion_free_rank(class_group(variety))
    current_sequences = Vector{Vector{Vector{Int}}}([ [zeros(Int, cl_rank)] ])
    
    for i in 1:search_depth
        extended_sequences = Vector{Vector{Vector{Int}}}()
        
        # Outer Loop: Iterate over the potential next points (P_new)
        for p_new in considered_points
            
            # --- CRITICAL STEP 1: INITIALIZE CACHE FOR P_new ---
            # This cache stores all vectors Q that are known to fail the check against this specific P_new.
            point_specific_forbidden_predecessors = Set{Vector{Int}}()
            
            # Inner Loop: Iterate over all existing valid sequences
            for seq in current_sequences 
                
                # --- CRITICAL STEP 2: PRUNING CHECK ---
                # Check if the sequence contains any vector Q already known to be forbidden
                # as a predecessor for this specific p_new.
                if !isempty(point_specific_forbidden_predecessors) && any(q -> q in point_specific_forbidden_predecessors, seq)
                    # If the sequence contains a forbidden predecessor Q, skip the check entirely.
                    continue
                end

                # --- CRITICAL STEP 3: INCREMENTAL CHECK & CACHE POPULATION ---
                # Check the new pairs (p_new, v_j). We must be able to identify which v_j failed.
                
                v_j_failed = false
                
                for v_j in seq
                    # Perform the expensive check: is (p_new - v_j) NOT in the coimmaculate locus?
                    if !is_difference_coimmaculate(coimm_points, p_new, v_j)
                        # Failure: v_j is a forbidden predecessor for p_new.
                        push!(point_specific_forbidden_predecessors, v_j)
                        v_j_failed = true
                        break # Stop checking this sequence; it failed.
                    end
                end
                
                # --- CRITICAL STEP 4: EXTEND ---
                if !v_j_failed
                    extended = push!(copy(seq), p_new)
                    push!(extended_sequences, extended)
                end
            end
        end
        
        print_status && println("ES after Step $(i): $(length(extended_sequences))")
        
        current_sequences = extended_sequences
    end
    
    return current_sequences 
end

########################################
# Functions to find configurations of mes
########################################

function get_canonical_base(
    s::Vector{Vector{Int}}, 
    R::Vector{Vector{Int}}
)::Tuple{Vector{Vector{Int}}, Int}
    
    N = length(s)
    if N != length(R)
        throw(DimensionMismatch("Sequences s and R must have the same length."))
    end
    
    # --- 1. Define the Objective Function (L1-norm) ---
    function score_sequence(base_seq::Vector{Vector{Int}})::Int
        # Calculate the L1-norm of the entire sequence: Sum of absolute values of all components
        return sum(sum(abs, v) for v in base_seq)
    end

    # --- 2. Determine Optimal Search Range for λ (Convexity) ---
    
    lambda_avg_numerator = 0.0
    lambda_avg_denominator = 0.0
    
    for i in 1:N
        lambda_avg_numerator += dot(s[i], R[i])
        lambda_avg_denominator += dot(R[i], R[i])
    end

    # Handle the trivial case where the shift sequence R is all zeros
    if iszero(lambda_avg_denominator)
        return s, 0
    end
    
    # The minimum is likely near the integer closest to this real-valued lambda
    lambda_real = lambda_avg_numerator / lambda_avg_denominator
    
    # Search range centered on the nearest integer solution
    lambda_test_start = floor(Int, lambda_real) - 2
    lambda_test_end = ceil(Int, lambda_real) + 2
    
    # --- 3. Search and Find Minimum L1-Norm ---
    
    min_score = typemax(Int)
    # Stores (lambda, base_seq) for all candidates at the minimum score
    candidate_bases = Vector{Tuple{Int, Vector{Vector{Int}}}}() 
    
    # NOTE: No need to initialize best_s_base here; it's handled by finding the first candidate
    
    for lambda_test in lambda_test_start:lambda_test_end
        
        s_base_candidate = Vector{Vector{Int}}(undef, N)
        current_score = 0
        
        for i in 1:N
            s_base_candidate[i] = s[i] - lambda_test * R[i]
            current_score += sum(abs, s_base_candidate[i])
        end
        
        if current_score < min_score
            min_score = current_score
            empty!(candidate_bases) # Reset candidates if a strictly better score is found
            push!(candidate_bases, (lambda_test, s_base_candidate))
        elseif current_score == min_score
            push!(candidate_bases, (lambda_test, s_base_candidate))
        end
    end
    
    # --- 4. Lexicographical Tie-Breaker ---
    
    if isone(length(candidate_bases))
        return candidate_bases[1][2], candidate_bases[1][1] # Return (s_can, lambda_can)
    end
    
    # If there are ties, select the lexicographically smallest sequence
    best_s_base = candidate_bases[1][2]
    best_lambda = candidate_bases[1][1]
    
    # Flatten the first sequence for comparison initialization
    best_flat = reduce(vcat, best_s_base)
    
    for k in 2:length(candidate_bases)
        lambda_candidate = candidate_bases[k][1]
        s_base_candidate = candidate_bases[k][2]

        current_flat = reduce(vcat, s_base_candidate)

        # Lexicographical comparison: If the new sequence is smaller
        if current_flat < best_flat 
            best_s_base = s_base_candidate
            best_lambda = lambda_candidate
            best_flat = current_flat
        end
    end
    
    return best_s_base, best_lambda
end

function get_remainder(s::Vector{Vector{Int}}, t::Vector{Vector{Int}})
    # Check dimensions early
    if length(s) != length(t)
        throw(DimensionMismatch("Sequences s and t must have the same length."))
    end
    
    delta_sequence = s .- t
    
    if all(iszero, delta_sequence)
        # Sequence is a vector of vectors of integers. Find the dimension D.
        D = length(s[1]) 
        zero_remainder = [zeros(Int, D) for _ in 1:length(s)]
        return (zero_remainder, 0)
    end
    
    # Find the Scalar GCD of all non-zero components
    non_zero_components = Int[]
    for delta_i in delta_sequence
        append!(non_zero_components, delta_i)
    end
    filter!(x -> x != 0, non_zero_components)

    lambda = reduce(gcd, non_zero_components)
    remainder = [div.(v, lambda) for v in delta_sequence]
    
    # Ensure remainder is uniquely represented by making its first non-zero component positive.
    first_non_zero_value = 0
    
    for R_i in remainder
        for r_ij in R_i
            if r_ij != 0
                first_non_zero_value = r_ij
                break # Found the first non-zero component
            end
        end
        if first_non_zero_value != 0
            break
        end
    end
    
    # If the first non-zero component is negative, flip the sign of the entire R_basis
    if first_non_zero_value < 0
        remainder = remainder .* -1
        lambda = -lambda # The scalar factor must also flip sign
    end
    
    return (remainder, lambda)
end

function get_pairwise_canonical_signatures(seq::Vector{Vector{Vector{Int}}})
    max_idx = length(seq)
    
    # Store the results
    signature_list = Vector{Tuple{Tuple{Int, Int}, Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}}}()

    # Calculate the base/remainder pairs
    for i in 1:max_idx-1
        seq_i = seq[i]

        for j in i+1:max_idx
            seq_j = seq[j]

            rem_ij, factor_ij = get_remainder(seq_i,seq_j)

            base_i, coeff_i = get_canonical_base(seq_i, rem_ij)
            base_j, coeff_j = get_canonical_base(seq_j, rem_ij) #Not necessary

            push!(signature_list, ((i, j), (base_i, rem_ij)))

        end
    end

    return signature_list
end

function count_canonical_key_occurrences(signature_list::Vector)
    
    # Initialize a dictionary to store the counts.
    # Key: (Canonical Base, Canonical Remainder)
    # Value: Count of occurrences
    key_counts = Dict{Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}, Int}()
    
    # Iterate through the signature list
    for signature_entry in signature_list
        # The canonical key is the second element of the signature_entry tuple:
        # signature_entry[1] is the (i, j) index tuple
        # signature_entry[2] is the (base, remainder) key tuple
        canonical_key = signature_entry[2]
        
        # Increment the count for the key. get(key_counts, key, 0) retrieves the 
        # current count or 0 if the key is new.
        key_counts[canonical_key] = get(key_counts, canonical_key, 0) + 1
    end
    
    return key_counts
end

function filter_canonical_keys_by_occurrence(
    key_counts::Dict{Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}, Int}, 
    X::Int
)::Vector{Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}}
    
    # Use a list comprehension for concise filtering based on the value (count)
    confirmed_keys = [
        key 
        for (key, count) in key_counts 
        if count >= X
    ]
    
    return confirmed_keys
end


function map_keys_to_cluster_indices(
    confirmed_keys::Vector{Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}},
    signature_list::Vector
)::Dict{Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}, Set{Int}}
    
    # Initialize the dictionary to store the mapping
    # The value is now a Set{Int} to store unique indices
    key_to_indices_set = Dict{
        Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}, 
        Set{Int}
    }()

    # Pre-populate the dictionary with all confirmed keys, initializing with an empty Set
    for key in confirmed_keys
        key_to_indices_set[key] = Set{Int}()
    end

    # Iterate through the full list of pairwise signatures
    for entry in signature_list
        indices_tuple = entry[1]     # The (i, j) index pair
        canonical_key = entry[2]     # The (base, remainder) key

        # Check if this canonical_key is one of the confirmed, highly-supported keys
        if haskey(key_to_indices_set, canonical_key)
            
            # If it is, add both indices i and j to the corresponding Set
            i, j = indices_tuple
            
            push!(key_to_indices_set[canonical_key], i)
            push!(key_to_indices_set[canonical_key], j)
        end
    end
    
    return key_to_indices_set
end

function get_union_of_all_index_sets(
    key_to_indices_set::Dict{Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}, Set{Int}}
)::Set{Int}
    
    # Initialize an empty Set to accumulate all indices
    total_union = Set{Int}()
    
    # Iterate through the values (which are the index Sets) of the dictionary
    for index_set in values(key_to_indices_set)
        # Use the union! function (mutating union) to efficiently add all 
        # elements of index_set into total_union
        union!(total_union, index_set)
    end
    
    return total_union
end

function generate_symbolic_general_classes(
        key_to_indices_set::Dict{Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}, Set{Int}}
    )
    
    p = Sym("p")
    symbolic_classes = Vector{Vector{Vector{Sym{PyObject}}}}()
    
    # Iterate through all unique canonical keys (which define the general patterns)
    for key in keys(key_to_indices_set)       
        S_can = key[1]  # The Canonical Base Sequence
        R_can = key[2]  # The Canonical Remainder Sequence (Shift Vector R)
        
        # First occurence of a parameter should be p
        conf = unify_sequence(S_can + p * R_can)       
        push!(symbolic_classes, conf)       
    end
    
    return symbolic_classes
end

function find_configurations_of_mes(seq_vec::Vector{Vector{Vector{Int}}}; threshold::Int=2)
    signatures = get_pairwise_canonical_signatures(seq_vec)
    signatures_occurences = count_canonical_key_occurrences(signatures)
    
    # Only consider the canonical_keys with multiple appearances
    signatures_relevant = filter_canonical_keys_by_occurrence(signatures_occurences, threshold)
    
    # Assign each signature to one of the relevant signatures
    classes = map_keys_to_cluster_indices(signatures_relevant, signatures)
    
    # Check if each class was mapped to a relevant class
    mapped = get_union_of_all_index_sets(classes)
    
    if length(mapped) != length(seq_vec)
        @warn "Not all sequences fit into a general description!"
    end
    
    conf_vec = generate_symbolic_general_classes(classes) 
    
    return conf_vec
end


########################################
# OPERATIONS ON EXCEPTIONAL SEQUENCES
########################################


@doc raw"""
    calculate_augmentation_on_surface(
        mes::Vector{Vector{T}},
        E::Vector{Int},
        pullback
    )

Calculate the augmentions of a maximal exceptional sequence on a toric surface.

# Arguments
- `mes::Vector{Vector{T}}`: the maximal exceptional sequence whose augmentations should be calculated
- `E::Vector{Int}`: the exceptional divisor of the blow up
- `pullback`: the pullback function of divisors

# Examples
```jldoctest
julia> p2_S = [[0], [1], [2]];

julia> h1_E = [-1,1];

julia> h1_pullback = (vec) -> [0, vec[1]];

julia> calculate_augmentation_on_surface(p2_S, h1_E, h1_pullback)
3-element Vector{Vector{Vector{Int64}}}:
 [[0, 0], [-1, 1], [0, 1], [0, 2]]
 [[0, 0], [1, 0], [0, 1], [1, 1]]
 [[0, 0], [0, 1], [1, 1], [0, 2]]

```
"""
function calculate_augmentation_on_surface(
        mes::Vector{Vector{T}},
        E::Vector{Int},
        pullback
    ) where T
    
    augmentations = Vector{Vector{T}}[]    
    pmes = pullback.(mes)
    
    for i in 1:length(mes)
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


function calculate_all_augmentations_on_surface(
    all_mes::Vector{Vector{Vector{T}}},
    E::Vector{Int},
    pullback
) where T
    
    # Initialize a single, flat vector to hold all resulting augmented sequences
    all_augmentations = Vector{Vector{T}}[]
    
    # Iterate over each "mes" sequence in the input container
    for mes_sequence in all_mes
        # Calculate augmentations for the current single sequence
        results = calculate_augmentation_on_surface(mes_sequence, E, pullback)
        
        # Append the results of the current sequence to the main container
        append!(all_augmentations, results)
    end
    
    return all_augmentations
end




@doc raw"""
    calculate_helixing(
        mes::Vector{Vector{T}},
        K::Vector{Int}
    )

Calculate the helixing of a maximal exceptional sequence.

# Arguments
- `mes::Vector{Vector{T}}`: the maximal exceptional sequence that is helexed
- `K::Vector{Int}`: the anticanonical divisor -K_X in Pic(X) of a toric variety X

# Examples
```jldoctest
julia> mes = [[0,0,0], [1,1,-1], [1,0,0], [0,1,0], [1,1,0]];

julia> calculate_helixing(mes, [1,1,1])
5-element Vector{Vector{Int64}}:
 [0, 0, 0]
 ⋮
 [0, 0, 2]

```
"""
function calculate_helixing(
        mes::Vector{Vector{T}},
        K::Vector{Int}
    ) where T
    l = length(mes)
    l < 2 && throw(ArgumentError("mes must contain at least 2 vectors"))

    helexed_mes = mes[2:l]
    push!(helexed_mes, mes[1]+K)
    helexed_mes = [vec.-helexed_mes[1] for vec in helexed_mes]

    if typeof(mes) == Vector{Vector{Sym{PyObject}}}
        return unify_sequence(helexed_mes)
    end
    
    return helexed_mes
end


@doc raw"""
    calculate_all_helixing(
        mes::Vector{Vector{T}},
        K::Vector{Int}
    )

Calculate all helixings of a maximal exceptional sequence.

# Arguments
- `mes::Vector{Vector{T}}`: the maximal exceptional sequence that is helexed
- `K::Vector{Int}`: the anticanonical divisor -K_X in Pic(X) of a toric variety X

# Examples
```jldoctest
julia> mes = [[0,0,0], [1,1,-1], [1,0,0], [0,1,0], [1,1,0]];

julia> calculate_all_helixing(mes, [1,1,1])
5-element Vector{Vector{Vector{Int64}}}:
 [[0, 0, 0], [1, 1, -1], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
 ⋮
 [[0, 0, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]]

```
"""
function calculate_all_helixing(
        mes::Vector{Vector{T}},
        K::Vector{Int}
    ) where T

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
end


@doc raw"""
    calculate_dualizing(
        mes::Vector{Vector{T}},
        K::Vector{Int}
    )

Calculate the dualized sequence of a mes.

# Arguments
- `mes::Vector{Vector{T}}`: the maximal exceptional sequence that is dualized
- `K::Vector{Int}`: the anticanonical divisor -K_X in Pic(X) of a toric variety X

# Examples
```jldoctest
julia> mes = [[0,0,0], [1,1,-1], [1,0,0], [0,1,0], [1,1,0]];

julia> calculate_dualizing(mes, [1,1,1])
5-element Vector{Vector{Int64}}:
 [0, 0, 0]
 ⋮
 [1, 1, 0]

```
"""
function calculate_dualizing(
        mes::Vector{Vector{T}},
        K::Vector{Int}
    ) where T
    l = length(mes)
    
    inverted_mes = [-K .- mes[l+1-i] for i in 1:l]
    inverted_mes = [vec.-inverted_mes[1] for vec in inverted_mes]

    if typeof(mes) == Vector{Vector{Sym{PyObject}}}
        return unify_sequence(inverted_mes)
    end
    
    return inverted_mes
end


@doc raw"""
    calculate_flipping(
        V::Vector{Vector{T}},
        perm::Vector{Int}
    )

Calculate the flipping for a given permutation of a maximal exceptional sequence on a toric surface.

# Arguments
- `mes::Vector{Vector{T}}`: the maximal exceptional sequence that is flipped
- `K::Vector{Int}`: the anticanonical divisor -K_X in Pic(X) of a toric variety X

# Examples
```jldoctest
julia> mes = [[0,0,0], [1,1,-1], [1,0,0], [0,1,0], [1,1,0]];

julia> calculate_flipping(mes, [2,3,1])
5-element Vector{Vector{Int64}}:
 [0, 0, 0]
 ⋮
 [1, 0, 1]

```
"""
function calculate_flipping(
        V::Vector{Vector{T}},
        perm::Vector{Int}
    ) where T
    
    flipping = Vector{T}[]
    
    for v in V
        push!(flipping, v[perm])
    end
  
    return flipping
end


########################################
# FULLNESS CONDITIONS
# cf. https://arxiv.org/abs/2304.01269
########################################

function relative_height(variety::NormalToricVariety, si::Vector{Int64}, sj::Vector{Int64})
    diff = sj - si
    tdc = toric_divisor_class(variety, diff)
    tl = toric_line_bundle(toric_divisor(tdc))
    cohom = all_cohomologies(tl)
    k = findfirst(x -> x != 0, cohom)
    
    if isnothing(k)
        return Inf
    end
            
    return k-1
end

function anticanonical_pseudoheight(variety::NormalToricVariety, seq::Vector{Vector{Int}})
    n = length(seq)
    KX = -transform_abelian_group_element(picard_class(canonical_bundle(variety)))
    best = Inf

    # helper to compute relative heights
    function e(i::Int, j::Int)
        relative_height(variety, seq[i], seq[j])
    end

    # helper: e(D_{a_p}, D_{a₀} + K_X)
    function e_anticanonical(i::Int, j::Int)
        target = seq[j] + KX  # shift by K_X
        relative_height(variety, seq[i], target)
    end

    # iterate over all strictly increasing chains of indices
    for p in 1:n
        for chain in Combinatorics.combinations(1:n, p)
            s = 0
            bad = false
            # sum over e(D_ai, D_{ai+1})
            for t in 1:(p-1)
                val = e(chain[t], chain[t+1])
                if val == Inf
                    bad = true
                    break
                end
                s += val
            end
            if bad
                continue
            end
            # add the last term e(D_{a_p}, D_{a_0} + K_X)
            lastval = e_anticanonical(chain[end], chain[1])
            if lastval == Inf
                continue
            end
            total = s + lastval - (p - 1)
            best = min(best, total)
            
            #println([total,best])
        end
    end

    return Int(best)
end

function kuznetsov_fullness_condition(variety::NormalToricVariety, seq::Vector{Vector{Int}})
    ph_ac = anticanonical_pseudoheight(variety, seq)  
    return ph_ac <= -dim(variety)
end

# OUTDATED/WRONG since it does not treat p correctly (it has to be reduced by -1)

@doc raw"""
    is_not_full_by_krah(
        seq::Vector{Vector{Int}},
        variety::NormalToricVariety
    )

Check if ph_ac(E1,...,En) > -dim X

# Arguments
- `seq::Vector{Vector{Int}}`: a maximal exceptional sequence
- `variety::NormalToricVariety`: the toric variety for which the sequence comes from

# Examples
```jldoctest
julia> mes = [[0,0,0], [1,1,-1], [1,0,0], [0,1,0], [1,1,0]];

julia> is_not_full_by_krah(mes, define_pentagon())
false

```
"""
function is_not_full_by_krah(
        seq::Vector{Vector{Int}},
        variety::NormalToricVariety
    )
    
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
end


@doc raw"""
    calculate_anticanonical_pseudoheight(
        seq::Vector{Vector{Int}},
        variety::NormalToricVariety
    )

Calculate the anticanonical pseudoheight.

# Arguments
- `seq::Vector{Vector{Int}}`: a maximal exceptional sequence
- `variety::NormalToricVariety`: the toric variety for which the sequence comes from

# Examples
```jldoctest
julia> mes = [[0,0,0], [1,1,-1], [1,0,0], [0,1,0], [1,1,0]];

julia> calculate_anticanonical_pseudoheight(mes, define_pentagon())
-3

```
"""
function calculate_anticanonical_pseudoheight(
        seq::Vector{Vector{Int}},
        variety::NormalToricVariety
    )
    
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
    comb_all = filter(x -> length(x)>1, collect(Combinatorics.combinations(1:l))) 
    
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
end






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