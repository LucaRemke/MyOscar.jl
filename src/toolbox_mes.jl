# load additionally required packages

using Plots.PlotMeasures


########################################
########################################
# THE pentagon
########################################
########################################

####################
# Returns all possible vectors with values in the input ranges

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
end

####################
# Checks if the difference of two vectors belongs to the negative immaculate locus of the pentagon

function is_in_locus_of_pentagon(v::Vector{T}, V::Vector{Vector{T}}) where T
    difference = [v - w for w in V]
    
    eval = all(diff -> diff == [1,1,-1] || diff == [0,0,2] || 
                (diff[3] == 0 && (diff[1] == 1 || diff[2] == 1)) ||
                (diff[3] == 1 && (diff[1] == 0 || diff[2] == 0)), difference)
    
    
    #eval = all(diff -> diff == [0,-1,0] || diff == [-2,2,-2] || 
    #            (diff[1] == -1 && (diff[2] == 0 || diff[2] == 1)) ||
    #            (diff[3] == -1 && (diff[2] == 0 || diff[2] == 1)), difference)
    
    return eval
end

####################
# Main part to build exceptional sequences, by checking the condition s-t in -Imm for all points in the search range

function advanced_extend_sequences(E::Vector{Vector{Vector{T}}}, V::Vector{Vector{T}}) where T
    new_sequences = Vector{Vector{Vector{T}}}()
    
    for seq in E
        for vec in V
            # Check if the difference between vec and all vectors in seq is in P
            if is_in_locus_of_pentagon(vec, seq)
                push!(new_sequences, push!(copy(seq), vec)) 
            end
        end
    end
    
    return new_sequences
end



####################
# Searches all exceptional sequences of length depth in the search_range

function advanced_brute_force_exseq_for_pentagon_in_zero(
    search_range::Vector{UnitRange{T}},
    search_depth::Int        
) where T

considered_points = generate_vectors_from_ranges(search_range)
    
# Starting with sequences of length 1
# Since we can shift exceptional sequences we cn set the origin to zero
starting_points = [[0,0,0]]
initial_sequences = [[v] for v in starting_points]

# Perform iterations to extend sequences
for i in 1:search_depth
    println("We reached step $i")
    initial_sequences = advanced_extend_sequences(initial_sequences, considered_points)
    if initial_sequences == Vector{Vector{Int64}}[]
        println("There are no sequences after $i steps")
        break
    else
        println("There are $(length(initial_sequences)) sequences of length $(i+1) \n")
    end
end

return initial_sequences
end 

####################
# Function that assigns each vector in an exceptional sequence a number which encodes the part of the negative immaculate locus the vector belongs to
# The idea is to find patterns in the encoded sequences

function group_by_type_sequence(V::Vector{Vector{Vector{T}}}) where T
    type_sequences_dict = Dict{Vector{Int64}, Vector{Vector{Vector{T}}}}()
    
    for seq in V
        type_seq = Vector{Int64}()
        for vec in seq
            if vec == [1,1,-1]
                push!(type_seq, 0)
            elseif vec == [0,0,2]
                push!(type_seq, 1)
            elseif vec == [1,1,0]
                push!(type_seq, 2)
            elseif vec == [0,0,1]
                push!(type_seq, 3)
            elseif vec[3] == 0 && vec[1] == 1
                push!(type_seq, 4)
            elseif vec[3] == 0 && vec[2] == 1
                push!(type_seq, 5)
            elseif vec[3] == 1 && vec[1] == 0
                push!(type_seq, 6)
            elseif vec[3] == 1 && vec[2] == 0
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


####################
# Plots an exceptional sequence in different layers

function plot_type(
    seq::Vector{Vector{Int64}};
    plot_area::Vector{UnitRange{Int64}} = fill(-3:3, 3),
    order::Bool = true
)

# Define empty vector for the plots
plot_seq = Any[]

# Define a dataframe for each contamination step
df_seq = convert_to_df(seq)

# Determine the layers of the plot
x_lim_min = first(plot_area[1])
x_lim_max = last(plot_area[1])
y_lim_min = first(plot_area[2])
y_lim_max = last(plot_area[2])
z_lim_min = first(plot_area[3])
z_lim_max = last(plot_area[3])

for layer in plot_area[3]

    # Define an empty plot
    layer_plot = plot(
        xlims = (x_lim_min,x_lim_max), xticks = x_lim_min:1:x_lim_max,
        ylims = (y_lim_min,y_lim_max), yticks = y_lim_min:1:y_lim_max,
        zlims = (z_lim_min,z_lim_max), zticks = z_lim_min:1:z_lim_max,
        grid = :true, gridalpha = 0.5,
        framestyle = :box,
        camera = (10,20),
        margin = -2mm
    )

    # Add the xy-plane of height layer
    z_plane = ones(length(plot_area[1]), length(plot_area[2])) * layer

    surface!(
        plot_area[1], plot_area[2], z_plane, 
        alpha = 0.5, color = :gray, colorbar=false
    )
    
    # Determine the elements of the exceptional sequence in the current layer
    df_seq_layer = filter(row -> row.x3 == layer, df_seq)
    
    # Plot the elemetns in the current layer and label them
    plot!(
        df_seq_layer.x1, df_seq_layer.x2, df_seq_layer.x3, 
        seriestype = :scatter, framestyle = :axes, legend = :false, color = :red, markersize = 3
    )

    labeling = findall(row -> row.x3 == layer, eachrow(df_seq)) 

    for i in labeling                
        annotate!([(
            df_seq.x1[i], df_seq.x2[i], df_seq.x3[i], 
            text("$(i-1)", 11, :green, :left, :bottom)
        )])
    end 

    
    
    push!(plot_seq, layer_plot)
    
    
end

return plot_seq

end

####################
# Plots vector of exceptional sequences

# function for plotting the contamination process
function plot_type_allseq(
    seq::Vector{Vector{Vector{Int64}}};
    plot_area::Vector{UnitRange{Int64}} = fill(-3:3, 3),
    order::Bool = true
)

# Define empty vector for the plots
plot_seq = Any[]

# Define a dataframe for each contamination step
df_seq = convert_to_df.(seq)
df1 = df_seq[1]
df2 = df_seq[2]

# Determine duplicates, as they should have another color
result_vector = Vector{Int}(undef, nrow(df1)) 
for i in 1:nrow(df1)
    result_vector[i] = isequal(df1[i, :], df2[i, :]) ? 1 : 0
end

for df in df_seq
    df.x4 = result_vector
end    

# Determine the layers of the plot
x_lim_min = first(plot_area[1])
x_lim_max = last(plot_area[1])
y_lim_min = first(plot_area[2])
y_lim_max = last(plot_area[2])
z_lim_min = first(plot_area[3])
z_lim_max = last(plot_area[3])

for df in df_seq
    
    # Define an empty plot
    layer_plot = plot(
        xlims = (x_lim_min,x_lim_max), xticks = x_lim_min:1:x_lim_max,
        ylims = (y_lim_min,y_lim_max), yticks = y_lim_min:1:y_lim_max,
        zlims = (z_lim_min,z_lim_max), zticks = z_lim_min:1:z_lim_max,
        grid = :true, gridalpha = 0.5,
        framestyle = :box,
        camera = (10,20)
    )

    for layer in plot_area[3]
        # Add the xy-plane of height layer
        z_plane = ones(length(plot_area[1]), length(plot_area[2])) * layer

        surface!(
            plot_area[1], plot_area[2], z_plane, 
            alpha = 0.5, color = :gray, colorbar=false
        )

        # Determine the elements of the exceptional sequence in the current layer
        df_seq_layer = filter(row -> row.x3 == layer, df)

        # Plot the elemetns in the current layer and label them
        colors = ifelse.(df_seq_layer.x4 .== 0, "blue", "red")
        
        plot!(
            df_seq_layer.x1, df_seq_layer.x2, df_seq_layer.x3, 
            seriestype = :scatter, framestyle = :axes, legend = :false, color = colors, markersize = 3
        )

        labeling = findall(row -> row.x3 == layer, eachrow(df)) 

        for i in labeling                
            annotate!([(
                df.x1[i], df.x2[i], df.x3[i], 
                text("$(i-1)", 11, :green, :left, :bottom)
            )])
        end
    end

    push!(plot_seq, layer_plot)
end
    

return plot_seq

end