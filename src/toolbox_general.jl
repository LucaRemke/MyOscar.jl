########################################
# CREATE OR MANIPULATE VECTORS
########################################


@doc raw"""
    generate_vectors(n::Int, r::UnitRange{Int})

Generates all vectors of length n with entries in range r.

# Arguments
- `n::Int`: Length of the vectors
- `r::UnitRange{Int}`: Range of the coefficients in the vectors

# Examples
```jldoctest
julia> generate_vectors(2, 0:1)
4-element Vector{Vector{Int64}}:
 [0, 0]
 ⋮
 [1, 1]


```
"""
function generate_vectors(n::Int, r::UnitRange{Int})
    vectors = [Int[]]
    for _ in 1:n
        new_vectors = Vector{Vector{Int}}(undef, length(r) * length(vectors))
        idx = 1
        for v in vectors, x in r
            new_vectors[idx] = [v... , x]
            idx += 1
        end
        vectors = new_vectors
    end
    return vectors
end


@doc raw"""
    generate_vectors_from_ranges(r::Vector{UnitRange{T}})

Generates all vectors of with entries in a unit range.
The first range defiens the range for the first entries, etc.

# Arguments
- `r::Vector{UnitRange{T}}`: Ranges corresponding to the entries of the generated vector

# Examples
```jldoctest
julia> generate_vectors_from_ranges([0:1, 2:3])
4-element Vector{Vector{Int64}}:
 [0, 2]
 ⋮
 [1, 3]


```
"""
function generate_vectors_from_ranges(r::Vector{UnitRange{T}}) where {T<:Integer}
    vectors = [Vector{T}([])]
    for range in r
        new_vectors = Vector{Vector{T}}(undef, length(vectors) * length(range))
        idx = 1
        for vec in vectors, elem in range
            new_vectors[idx] = vcat(vec, elem)
            idx += 1
        end
        vectors = new_vectors
    end
    return vectors
end


@doc raw"""
    transform_rayvector(R::Any)

Transform a indexable object of any type into integer type.

# Arguments
- `R::Any`: Object which entries are transformed into integers

# Examples
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> td = toric_divisor(h1, [1,2,3,4]);

julia> transform_rayvector(coefficients(td))
4-element Vector{Int64}:
 1
 ⋮
 4


```
"""
function transform_rayvector(R::Any)
    vect = Int[]
    for j in eachindex(R)
        entry = Int(R[j])
        push!(vect,entry)
    end
    
    return vect
end


@doc raw"""
    transform_raysequence(seq::Vector{Vector{T}})

Transforms a Vector of RayVector or Vector with ring elements into a vector of integer vectors.

# Arguments
- `seq::Vector{Vector{T}}`: The vector of vectors which is transformed into a vector of integer vectors

# Examples
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> td1 = toric_divisor(h1, [1,2,3,4]);

julia> td2 = toric_divisor(h1, [4,3,2,1]);

julia> transform_raysequence([coefficients(td1), coefficients(td2)])
2-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [4, 3, 2, 1]


```
"""
transform_raysequence(seq::Vector{Vector{T}}) where T = [transform_rayvector(vec) for vec in seq]


"""
    transform_abelian_group_element(g::FinGenAbGroupElem) -> Vector{Int64}

Convert a finitely generated abelian group element `g` into a plain vector
of integers by iterating over its coordinates.

Each coordinate of `g` (accessed via `g[j]` for `j = 1:rank(parent(g))`)
is converted to `Int64` and collected into a standard Julia vector.

# Arguments
- `g::FinGenAbGroupElem`: an element of a finitely generated abelian group, 
  as provided by the Oscar/Hecke algebraic structures.

# Returns
- `Vector{Int64}`: a vector containing the integer coordinates of `g`
  with respect to the chosen generators of its parent group.

# Example
```jldoctest
julia> h1 = define_hirzebruch_surface(1);

julia> Kh1 = picard_class(canonical_bundle(h1));

julia> transform_abelian_group_element(Kh1)


```
"""
function transform_abelian_group_element(g::FinGenAbGroupElem)
    vect = Int[]
    for j in 1:rank(parent(g))
        entry = Int(g[j])
        push!(vect,entry)
    end
    
    return vect
end


function multiply_nested_vectors_map(M::Matrix{Int}, L::Vector{Vector{Vector{Int}}})
    # Map over the outermost layer (L)
    return map(L) do outer_vector_of_vectors
        
        # Map over the intermediate layer (Vector{Vector{Int}})
        return map(outer_vector_of_vectors) do inner_vector_of_vectors
            
            # Map over the innermost layer (Vector{Int}) and apply the multiplication (M * v)
            # Broadcasting (M * .) is used here to apply the matrix multiplication 
            # to every vector v inside the innermost collection.
            return M * inner_vector_of_vectors
        end
    end
end


@doc raw"""
    convert_to_df(vec::Vector{Vector{T}})

Convert a vector of vector into a dataframe.

# Arguments
- `vec::Vector{Vector{T}}`: The vector of vectors which is transformed into a dataframe

# Examples
```jldoctest
julia> convert_to_df([[1,2], [3,4]])
2×2 DataFrame
 Row │ x1     x2
     │ Int64  Int64
─────┼──────────────
  ⋮  │   ⋮      ⋮
      2 rows omitted


```
"""
function convert_to_df(vec::Vector{Vector{T}}) where T
    matrix = transpose(hcat(vec...)) 
    df = DataFrame(matrix, :auto)
    return df
end


@doc raw"""
    sort_points_clockwise(points::Vector{<:AbstractVector{<:Real}})

Sort a vector of point clockwise around its centroid.

# Arguments
- `points::Vector{<:AbstractVector{<:Real}}`: Points which are sorted clockwise
# Examples
```jldoctest
julia> sort_points_clockwise([[1,0], [-1,0], [0,1], [0,-1]])
4-element Vector{Vector{Int64}}:
 [0, -1]
 ⋮
 [-1, 0]


```
"""
function sort_points_clockwise(points::Vector{<:AbstractVector{<:Real}})
    n = length(points)
    cx = sum(p[1] for p in points) / n
    cy = sum(p[2] for p in points) / n

    angles = Vector{Float}(undef, n)
    for i in 1:n
        p = points[i]
        angles[i] = atan(p[2] - cy, p[1] - cx)
    end

    sorted_indices = sortperm(angles)
    sorted_points = [points[i] for i in sorted_indices]
    
    return sorted_points
end


########################################
# MISCELLANEOUS
########################################

function nesting_depth(v)
    depth = 0
    while v isa AbstractVector && !isempty(v)
        depth += 1
        v = first(v)
    end
    return depth
end

function transition_matrix(B::AbstractVector, C::AbstractVector)
    B_matrix = hcat(B...)
    C_matrix = hcat(C...)
    
    # 2. Calculation using left-division (C \ B)
    # This solves the system C * P = B for P, which is P_C <- B.
    # This is equivalent to inv(C) * B, but is numerically more stable and faster.
    P_C_to_B = C_matrix \ B_matrix
    
    return Int.(P_C_to_B)
end



########################################
# PRINT LATEX STRINGS
########################################

@doc raw"""
    print_latex_code(
        V::Vector{Vector{Vector{Vector{T}}}},
        path::String,
        name::String
    )

Transform a vector of vector of vectors into latex code.

# Arguments
- `V::Vector{Vector{Vector{Vector{T}}}}`: Vector which is transformed into latex code
- `path::String`: Define the path where the latex code is saved
- `name::String`: Define the name of the vector
"""
function print_latex_code(
        V::Vector{Vector{Vector{Vector{T}}}},
        path::String,
        name::String
    ) where T

    n_V = length(V)
    n_seq = length(V[1])
    n_vec = length(V[1][1])
    n_ele = length(V[1][1][1])

    cols_string = "\\begin{array}{|c|" * repeat("c|", n_vec-1) * "}\n"

    buf = IOBuffer()

    for index_hex in 1:n_V
        println(buf, "\\resizebox{\\textwidth}{!}{\$\\displaystyle")
        println(buf, "\\renewcommand{\\arraystretch}{2.5}")
        print(buf, cols_string)
        println(buf, "\\hline")
        println(buf, "$(name)_{$index_hex}^{0} &")
        println(buf, "\\begin{psmallmatrix}")

        # First sequence (helex = 0)
        seq0 = V[index_hex][1]
        for row in 1:n_ele
            for col in 1:n_vec
                el = seq0[col][row]
                if col != n_vec
                    print(buf, "$el & ")
                else
                    print(buf, el)
                end
            end
            println(buf, " \\\\")
        end
        println(buf, "\\end{psmallmatrix}")
        println(buf, "&")

        # Remaining sequences
        for helex in 1:n_seq-1
            println(buf, "$(name)_{$index_hex}^{$helex} &")
            println(buf, "\\begin{psmallmatrix}")
            seqh = V[index_hex][helex+1]
            for row in 1:n_ele
                for col in 1:n_vec
                    el = seqh[col][row]
                    if col != n_vec
                        print(buf, "$el & ")
                    else
                        print(buf, el)
                    end
                end
                println(buf, " \\\\")
            end
            println(buf, "\\end{psmallmatrix}")

            if (helex == 2 || helex == n_seq-1)
                println(buf, "\\\\[1.5ex]")
                println(buf, "\\hline")
            else
                println(buf, "&")
            end
        end

        println(buf, "\\end{array}\$")
        println(buf, "}")
        println(buf)
    end

    open(path, "w") do file
        write(file, String(take!(buf)))
    end

    return "File was saved as $path"
end