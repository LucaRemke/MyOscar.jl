module MyOscar

########################################
# load required packages
########################################
using Pkg
using Oscar
using Plots
gr()
using DataFrames
using CSV
using LinearAlgebra
using Combinatorics


########################################
########################################
include("toolbox_general.jl")

export generate_vectors
export transform_rayvector
export convert_to_df
export sort_points_clockwise

export print_latex_code


########################################
########################################
include("toolbox_toric_varieties.jl")

export define_projective_space
export define_projective_surface
export define_hirzebruch_surface
export define_pp
export define_ppp
export define_pentagon
export define_hexagon

export get_representative_of_divisor
export show_generators_and_relations_of_classgroup
export calculate_cartier_data
export check_nef_conditions
export calculate_nef_pair

export get_exceptional_divisor


########################################
########################################
include("toolbox_immaculate_locus.jl")

export get_immaculate_lb

export calculate_tempting_subsets
export search_nontemptings_via_boundary_condition
export eliminate_faces_as_temptings
export get_temptings_as_primitive_collections
export get_images_of_cube_vertices
export get_polytope_of_cube_vertices

export get_maculate_regions
export cutout_maculate_regions
export print_maculate_region_info

export get_seed_and_chull_of_toricpic3

export calculate_planary_parallelograms


########################################
########################################
include("toolbox_exceptional_sequences.jl")

export nimmloc

export configurate_exceptional_sequences_on_h0
export configurate_exceptional_sequences_on_hr
export configurate_exceptional_sequences_on_pentagon
export configurate_exceptional_sequences_on_hexagon
export configurate_exceptional_sequences_on_p1p1p1
export configurations

export is_in_nimmloc
export extend_exceptional_sequences
export generate_exceptional_sequences

export is_exceptional

export calculate_augmentation_on_surface
export calculate_helixing
export calculate_all_helixing
export calculate_dualizing

export calculate_flipping


########################################
########################################
include("toolbox_visualization.jl")

export plot_data

export create_gif_from_df
export create_gif_from_vector

export plot_maculate_regions
export plot3d_maculate_regions
export visualize_projection_of_immaculate_linebundles_for_pic3

export visualize3d_sequence_by_layers
export visualize3d_sequence_moving_lines
export visualize3d_sequence_process_by_layers


########################################
########################################
include("toolbox_sympy.jl")
# add function if it works when pushing on github


#=
These function might be needed for older code
##########
#include("toolbox_mes.jl")

#export generate_vectors_from_ranges
#export is_in_locus_of_pentagon
#export advanced_extend_sequences
#export advanced_brute_force_exseq_for_pentagon_in_zero
#export group_by_type_sequence
#export plot_type
#export plot_type_allseq


##########
#include("toolbox_exceptional_sequences.jl")
#
#export calculate_augmentations_of_mes
#export calculate_helex_of_mes
#export calculate_all_helex_of_mes
#export has_duplicates
#export calculate_inverted_mes
#export search_duplicates
=#

end;