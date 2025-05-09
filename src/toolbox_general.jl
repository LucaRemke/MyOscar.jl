########################################
# CREATE OR MANIPULATE VECTORS
########################################

#---------------------------------------
# Generates all vectors of length n with entries in range r
#---------------------------------------
function generate_vectors(n::Int, r::UnitRange{Int})
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

#---------------------------------------
# Returns all possible vectors with values in the input ranges
#---------------------------------------

function generate_vectors_from_ranges(r::Vector{UnitRange{T}}) where T
    # Base case: If there's only one range, return the vectors in that range
    if length(r) == 1
        return [elem for elem in r[1]]
    end
    
    # Recursive case: Generate vectors for the first range
    vectors = generate_vectors_from_ranges(r[1:end-1])
    new_range = r[end]
    
    # Extend each vector with elements from the last range
    new_vectors = Vector{Int}[]
    for vec in vectors
        for elem in new_range
            push!(new_vectors, [vec..., elem])
        end
    end
    
    return new_vectors
end;

#---------------------------------------
# Transforms a RayVector or Vector with ring elements into an integer vector
#---------------------------------------
function transform_rayvector(R::Any)
    vect = Int64[]
    for j in eachindex(R)
        entry = Int64(R[j])
        push!(vect,entry)
    end
    
    return vect
end;

#---------------------------------------
# Transforms a Sequence of Type T into a integer sequence
#---------------------------------------
function transform_raysequence(seq::Vector{Vector{T}}) where T
    int_seq = Vector{Int64}[]
    for vec in seq
        int_vec = transform_rayvector(vec)
        push!(int_seq, int_vec)  
    end
    return int_seq
end;

#---------------------------------------
# Convert a vector of vectors into a dataframe
#---------------------------------------
function convert_to_df(vec::Vector{Vector{T}}) where T
    matrix = transpose(hcat(vec...)) 
    df = DataFrame(matrix, :auto)
    return df
end;

#---------------------------------------
# Sorts a vector of point clockwise around its centroid
#---------------------------------------
function sort_points_clockwise(points)
    centroid = [sum(p[1] for p in points) / length(points), sum(p[2] for p in points) / length(points)]
    angles = [atan(p[2] - centroid[2], p[1] - centroid[1]) for p in points]
    sorted_indices = sortperm(angles)
    sorted_points = [points[i] for i in sorted_indices]
    
    return sorted_points
end;


########################################
# PRINT LATEX STRINGS
########################################

#---------------------------------------
# Prints a vector of vector of vectors as a latex table
#---------------------------------------
function print_latex_code(
        V::Vector{Vector{Vector{Vector{T}}}},
        path::String,
        name::String
    ) where T

    n_V = length(V)
    n_seq = length(V[1])
    n_vec = length(V[1][1])
    n_ele = length(V[1][1][1])

    cols_string = "\\begin{array}{|c|"
    cols_string *= repeat("c|", n_vec-1)
    cols_string *= "}\n"

    string = ""

    for index_hex in 1:n_V
        string *= "\\resizebox{\\textwidth}{!}{\$\\displaystyle\n"
        string *= "\\renewcommand{\\arraystretch}{2.5}\n"    
        string *= cols_string
        string *= "\\hline\n"
        string *= "$(name)_{$(index_hex)}^{0} &\n"
        string *= "\\begin{psmallmatrix}\n"
        
        for row in 1:n_ele
            for col in 1:n_vec
                el = V[index_hex][1][col][row]
                col != n_vec ? string *= "$el & " : string *= "$el "
            end
            string *= "\\\\\n"
        end
        
        string *= "\\end{psmallmatrix}\n"
        string *= "&\n"

        for helex in 1:n_seq-1
            string *= "$(name)_{$(index_hex)}^{$(helex)} &\n"
            string *= "\\begin{psmallmatrix}\n"

            for row in 1:n_ele
                for col in 1:n_vec
                    el = V[index_hex][helex+1][col][row]
                    col != n_vec ? string *= "$el & " : string *= "$el "
                end
            string *= "\\\\\n"
            end

            string *= "\\end{psmallmatrix}\n"

            if (helex == 2 || helex == n_seq-1)
                string *= "\\\\[1.5ex]\n"
                string *= "\\hline\n"
            else
                string *= "&\n"
            end
        end

        string *= "\\end{array}\$\n"
        string *= "}\n"
        #string *= "\\vspace*{0.35cm}\n"
        string *= "\n"
                    
    end

    open(path, "w") do file
        write(file, string)
    end

    return "File was saved as $path"
end    