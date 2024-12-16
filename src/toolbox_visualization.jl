using Plots.PlotMeasures

########################################
# PLOTING GENERAL DATA
########################################

#---------------------------------------
# Plot data of a Dataframe as 1/2/3D plot
#---------------------------------------
function plot_data(
        df::DataFrame; 
        df_title=nothing, 
        ax_labels::Bool=true, 
        style::Symbol=:origin, 
        camerapoint=nothing
    )
    
    n_cols = size(df,2)
    name_cols = names(df)
    
    x_min, x_max = extrema(df[!, name_cols[1]])
        
    if n_cols == 1        
        plot_df = plot(df[!, name_cols[1]], zeros(Int64, size(df,1)),
            xlims=(x_min-1.5,x_max+1.5),
            xticks=Int(x_min-1):1:Int(x_max+1), 
            yaxis=false, ylims=(-1,10), #ylims neu
            framestyle=:zerolines,
            seriestype=:scatter, size=(600,200) #size neu
        )
        
    elseif n_cols in [2,3]
        y_min, y_max = extrema(df[!, name_cols[2]])
        
        if n_cols == 2
            plot_df = plot(df[!, name_cols[1]], df[!, name_cols[2]],
                xticks=Int(x_min):1:Int(x_max),
                yticks=Int(y_min):1:Int(y_max),
                framestyle=style,
                seriestype=:scatter
            )
            
        else
            z_min, z_max = extrema(df[!, name_cols[3]])
            
            plot_df = plot(df[!, name_cols[1]], df[!, name_cols[2]], df[!, name_cols[3]],
                xticks=x_min:1:x_max,
                yticks=y_min:1:y_max,
                zticks=z_min:1:z_max,
                framestyle=style,
                seriestype=:scatter
            )
        
            if !isnothing(camerapoint)
                plot!(camera=camerapoint)
            end
            
            if ax_labels
                plot!(xlabel="x", ylabel="y", zlabel="z")
            end
        end

    else
        throw(ErrorException("Tuples of size $n_cols can not be plotted!")) 
    end
    
    plot!(legend=nothing)
    
    if !isnothing(df_title)
        plot!(title=df_title)
    end
    
    return plot_df
end;


########################################
# CREATING GIFS
########################################

#---------------------------------------
# Creates a gif from a dataframe
#---------------------------------------
function create_gif_from_df(
        df::DataFrame, 
        path::String; 
        frames_per_sec::Int64=20, 
        records::UnitRange{Int64}=0:360, 
        angle::Int64=165
    )

    name_cols = names(df)    
    plots = []
    anim = Animation()
    for i in records
        plot_i = plot(df[!, name_cols[1]], df[!, name_cols[2]], df[!, name_cols[3]],
            seriestype=:scatter,
            framestyle=:origin,
            label = nothing,
            camera = (i, angle)
        )
        push!(plots, plot_i)
        frame(anim)
    end

    gif_plot = gif(anim, path; fps=frames_per_sec)

    return gif_plot
end;

#---------------------------------------
# Creates a gif from a vector of plots
#---------------------------------------
function create_gif_from_vector(
        plots::Vector{Any}, 
        path::String, 
        dim::Tuple{Int64, Int64}; 
        frames_per_sec::Int64=20, 
        records::UnitRange{Int64}=0:360, 
        angle::Int64=165
    )

    all_plots = plot(plots..., size=dim)
    anim = Animation()

    for i in records
        plot!(camera = (i, angle))
        frame(anim)
    end

    gif_plot = gif(anim, path; fps=frames_per_sec)

    return gif_plot
end;


########################################
# PLOTTING IMMACULATE LOCUS + MACULATE REGIONS
########################################

#---------------------------------------
# Plot the maculate regions for a specific set of coordinates of a toric variety with Picard rank 2
#---------------------------------------
function plot_maculate_regions(regions::Vector{Any}; plotsize=(400,400), plottitle=nothing)

    maculate_lb = regions[1]
    hcone = regions[2]
    immaculate_lb = regions[3]
    plot_regions = plot(framestyle=:origin, size=plotsize, ticks=([], false), legend=:outertopright) 
    
    for i in eachindex(maculate_lb)
        region = maculate_lb[i]
             
        x_lp, y_lp = region[1].x1, region[1].x2
        x_bound, y_bound = region[2].x1, region[2].x2
        
        color_regions = palette(:devon10)[1+2*hcone[i]]
        
        plot!(x_bound,y_bound, seriestype=:shape, alpha=0.5, color=color_regions, label="H$(hcone[i])-cone")
        plot!(x_lp,y_lp, seriestype=:scatter, color=color_regions, label=:none)
    end
    
    # plot immaculate line bundles
    color_lb =  palette(:Oranges_9)[6]
    x_imm, y_imm = immaculate_lb.x1, immaculate_lb.x2
    plot!(x_imm,y_imm, seriestype=:scatter, color=color_lb, label="ImmZ(X)")
    
    if !isnothing(plottitle)
        title!(plottitle)
    end
    
    return plot_regions    
end;

#---------------------------------------
# Plot the maculate regions for a specific set of coordinates of a toric variety with Picard rank 3
#---------------------------------------
function plot3d_maculate_regions(
        regions::Vector{Any}; 
        plotall::Bool = false, 
        plotsize=(400,400), 
        plottitle="",
        plottitlesize = 12,
        plotlims=nothing,
        plotlegend=:outertopright
    )

    maculate_lb = regions[1]
    hcone = regions[2]
    immaculate_lb = regions[3]

    if plotall == false
        all_plots = []
    else
        plot_regions = plot(
            framestyle=:origin, size=plotsize, legend=plotlegend, title=plottitle, titlefontsize=plottitlesize
        ) 
    end
        
    for i in eachindex(maculate_lb)
        
        if plotall == false
            plot_regions = plot(
                framestyle=:origin, size=plotsize, legend=plotlegend, title=plottitle, titlefontsize=plottitlesize
            )
        end 
        
        region = maculate_lb[i]
            
        x_lp, y_lp, z_lp = region[1].x1, region[1].x2, region[1].x3
        x_bound, y_bound, z_bound = region[2].x1, region[2].x2, region[2].x3
        
        color_regions = palette(:devon10)[1+2*hcone[i]]
        plot!(x_lp,y_lp,z_lp, seriestype=:scatter, color=color_regions, label="H$(hcone[i])-cone")
        
        if !isnothing(plotlims)
            min, max = plotlims[1], plotlims[2]
            xlims!(min,max)
            ylims!(min,max)
            zlims!(min,max)           
        end
            
        plotall == false ? push!(all_plots, plot_regions) : nothing
    end

    # plot immaculate line bundles
    color_lb =  palette(:Oranges_9)[6]
    x_imm, y_imm, z_imm = immaculate_lb.x1, immaculate_lb.x2, immaculate_lb.x3

    if plotall == false
        plot_regions = plot(
            framestyle=:origin, size=plotsize, legend=plotlegend, title=plottitle, titlefontsize=plottitlesize
        )
    end

    plot!(x_imm,y_imm,z_imm, seriestype=:scatter, color=color_lb, label="ImmZ(X)")
        
    if plotall == false
        push!(all_plots, plot_regions)
        return all_plots
    else
        return plot_regions
    end
end;

#---------------------------------------
# Plot the projection of immaculate line bundles for toric varieties of Picard rank 3
#---------------------------------------
function visualize_projection_of_immaculate_linebundles_for_pic3(coeffs::Vector{Int64}; mysize::Tuple{Int64, Int64}=(800,800))
    
    # We start with para since calculate_planary_parallelograms checks for length
    para = calculate_planary_parallelograms(coeffs)
    p1, p2 = para[1], para[2]

    a,b,c,d = coeffs[1], coeffs[2], coeffs[3], coeffs[4]
        
    # Get the lattice points
    lp_1, lp_2 = lattice_points(p1), lattice_points(p2)
    tlp_1, tlp_2 = transform_rayvector.(lp_1), transform_rayvector.(lp_2)    
    df_1, df_2 = convert_to_df(tlp_1), convert_to_df(tlp_2)
    df_F = unique(vcat(df_1, df_2))
    
    # Get the vertices
    v1, v2 = vertices(p1), vertices(p2)
    tv1 = sort_points_clockwise(transform_rayvector.(v1))
    tv2 = sort_points_clockwise(transform_rayvector.(v2))    
    df_tv1 = convert_to_df(tv1)
    df_tv2 = convert_to_df(tv2)
    
    # Get the canonical divisor
    canonical = (-a-b+d, a-c-d)
    
    # Calculate the line segments of type (A)
    y = -c-d+1:a+b-1
    
    df_A = DataFrame(x1 = Int64[], x2 = Int64[])
    df_a = DataFrame(x1 = Int64[], x2 = Int64[])
    
    for i in y
        push!(df_A, [-i i])
        push!(df_a, [canonical[1]+i canonical[2]-i])
    end
    
    # Calculate the line segments of type (B)
    df_B = DataFrame(x1 = Int64[], x2 = Int64[])
    df_b = DataFrame(x1 = Int64[], x2 = Int64[])
    
    if b >= 2 && c >= 2
        push!(df_B, [-a-b a])
        push!(df_b, [canonical[1]+a+b canonical[2]-a])
    elseif b == 1 && c == 1
        # no immaculate linebundles in this case
    elseif b == 1
        for i in 0:c-1
            push!(df_B, [-a-b-i a])
            push!(df_b, [canonical[1]+a+b+i canonical[2]-a])
        end        
    elseif c == 1
        for i in 0:b-1
            push!(df_B, [-a-b a+i])
            push!(df_b, [canonical[1]+a+b canonical[2]-a-i])
        end        
    end    
    
    # Calculate the maximum of the dataframe for axis definition
    x_max = max(maximum(df_tv1.x1), maximum(df_tv2.x1))
    x_min = min(minimum(df_tv1.x1), minimum(df_tv2.x1))
    
    y_max = max(maximum(df_tv1.x2), maximum(df_tv2.x2))
    y_min = min(minimum(df_tv1.x2), minimum(df_tv2.x2))
    
    expand = 2
        
    plot_all = plot(
        grid=true, gridlinewidth=1, gridcolor=:gray, gridalpha=0.5,
        xticks=x_min-expand:1:x_max+expand, xlims=(x_min-expand, x_max+expand),
        yticks=y_min-expand:1:y_max+expand, ylims=(y_min-expand, y_max+expand),
        seriestype=:scatter,
        framestyle=:origin, size=mysize    
    )
    
    plot!(df_tv1.x1, df_tv1.x2, 
        seriestype=:shape, linecolor=palette(:tab10)[1], fillalpha=0.5, color = palette(:tab10)[1], label="P1"
    )
    plot!(df_tv2.x1, df_tv2.x2,
        seriestype=:shape, linecolor=palette(:tab10)[4], fillalpha=0.5, color = palette(:tab10)[4], label="P2"
    )
    
    # Type (F) --------------------------------
    scatter!(df_F.x1, df_F.x2, markersize=3, color = :black, label="(F)")
    
    # Type (A) --------------------------------
    scatter!(df_A.x1, df_A.x2, 
        markersize=3, markercolor = palette(:tab10)[2], markerstrokecolor = palette(:tab10)[2], label="(A)"
    )
    annotate!(df_A.x1, df_A.x2, text(string('A'), 11, palette(:tab10)[2], :bottom, :middle))
    
    scatter!(df_a.x1, df_a.x2, 
        markersize=3, markercolor = palette(:tab10)[2], markerstrokecolor = palette(:tab10)[2], label=:none
    )
    annotate!(df_a.x1, df_a.x2, text(string('a'), 11, palette(:tab10)[2], :top, :middle))
    
    # Type (B) --------------------------------
    
    if b != 1 || c != 1
    scatter!(df_B.x1, df_B.x2, 
        markersize=3, markercolor = palette(:tab10)[3], markerstrokecolor = palette(:tab10)[3], label="(B)"
    )
    annotate!(df_B.x1, df_B.x2, text(string('B'), 11, palette(:tab10)[3], :bottom, :middle))
    
    scatter!(df_b.x1, df_b.x2, 
        markersize=3, markercolor = palette(:tab10)[3], markerstrokecolor = palette(:tab10)[3], label=:none
    )
    annotate!(df_b.x1, df_b.x2, text(string('b'), 11, palette(:tab10)[3], :top, :middle)) 
    end
   
    return plot_all  
end;


########################################
# PLOTTING EXCEPTIONAL SEQUENCES
########################################

#---------------------------------------
# Plot a sequence in layers 
#---------------------------------------
function visualize3d_sequence_by_layers(
        seq::Vector{Vector{T}},
        plot_area::Vector{UnitRange{Int64}};
        x_shift::Float64=0.0,
        y_shift::Float64=0.0,
        cam::Tuple{Int64, Int64}=(10,20)
    ) where T

    # Define empty vector for the plots
    plot_seq = Any[]

    # Convert sequence into dataframe
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
            camera = cam,
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
                df_seq.x1[i] + x_shift, df_seq.x2[i] + y_shift, df_seq.x3[i], 
                text("$(i-1)", 11, :green, :right)
            )])
        end 

        push!(plot_seq, layer_plot)
    end

    return plot_seq
end;

#---------------------------------------
# Plot a sequence with indicating moving lines
#---------------------------------------
function visualize3d_sequence_moving_lines(
        seq::Vector{Vector{T}},
        plot_area::Vector{UnitRange{Int64}},
        color_vector::Vector{Int64},
        x_shift::Float64=0.0,
        y_shift::Float64=0.0;
        cam::Tuple{Int64, Int64}=(10,20),
        siz::Tuple{Int64, Int64}=(300,300),
        skip::Union{Int, Nothing}=nothing,
        transp::Float64=0.5
    ) where T

    df = convert_to_df(seq)
    df.x4 = color_vector

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
        grid = :true, gridalpha = 0.5,
        framestyle = :box,
        camera = cam,
        size = siz
    )

    for layer in plot_area[3]
        # Add the xy-plane of height layer
        z_plane = ones(length(plot_area[1]), length(plot_area[2])) * layer

        surface!(
            plot_area[1], plot_area[2], z_plane, 
            alpha = transp, color = :gray, colorbar=false
        )

        # Determine the elements of the exceptional sequence in the current layer
        df_seq_layer = filter(row -> row.x3 == layer, df)

        # Plot the elemetns in the current layer and label them
        colors = ifelse.(df_seq_layer.x4 .== 0, "blue", "red")

        plot!(
            df_seq_layer.x1, df_seq_layer.x2, df_seq_layer.x3, 
            seriestype = :scatter, framestyle = :axes, legend = :false, color = colors, markersize = 3
        )

        if layer != skip        
            labeling = findall(row -> row.x3 == layer, eachrow(df)) 
            for i in labeling                
                annotate!([(
                    df.x1[i] + x_shift, df.x2[i] + y_shift, df.x3[i], 
                    text("$(i-1)", 11, :green, :right)
                )])
            end
        end
    end

    return layer_plot
end;

#---------------------------------------
# Plot a sequence in layers for several p
#---------------------------------------
function visualize3d_sequence_process_by_layers(
        seq::Vector{Vector{Vector{T}}},
        plot_area::Vector{UnitRange{Int64}};
        x_shift::Float64=0.0,
        y_shift::Float64=0.0,
        cam::Tuple{Int64, Int64}=(10,20)
    ) where T

    # Define empty vector for the plots
    plot_seq = Any[]

    # Determine duplicates, as they should have another color
    df_seq = convert_to_df.(seq)
    df1 = df_seq[1]
    df2 = df_seq[2]

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
            camera = cam
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
                    df.x1[i] + x_shift, df.x2[i] + y_shift, df.x3[i], 
                    text("$(i-1)", 11, :green, :right)
                )])
            end
        end

        push!(plot_seq, layer_plot)
    end

    return plot_seq

end;