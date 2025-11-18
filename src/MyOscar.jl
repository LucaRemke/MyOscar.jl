module MyOscar

########################################
# load required packages
########################################
using Pkg
using Revise
using Test 
using Documenter
using Plots
gr()
using DataFrames
using CSV
using LinearAlgebra
using Combinatorics
using Serialization
using StatsBase
using Oscar

# See https://github.com/oscar-system/Oscar.jl/blob/master/src/imports.jl
for i in names(Oscar)
    @eval import Oscar: $i
    @eval export $i
  end


########################################
########################################
include("toolbox_general.jl")

export generate_vectors
export generate_vectors_from_ranges
export transform_rayvector
export transform_raysequence
export transform_abelian_group_element
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

export is_query_block_contained
export calculate_secondary_fan 
export find_adjacent_cones
export oriented_cocircuit
export extremal_position
export classify_vector


########################################
########################################
include("toolbox_immaculate_locus.jl")

export ilb_range
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
include("toolbox_infection.jl")

export generate_integral_tuples
export map_intergral_tuples
export negate_rule
export apply_base_change
export collect_infection_rules
export infection_step
export infection_procedure

export calculate_infectionstep_of_pentagon
export search_infections_of_pentagon
export search_infections_of_pentagon_bool

export calculate_infectionstep_of_hexagon
export search_infections_of_hexagon
export search_infections_of_hexagon_bool


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

export relative_height
export anticanonical_pseudoheight
export kuznetsov_fullness_condition

#export is_not_full_by_krah
#export calculate_anticanonical_pseudoheight


########################################
########################################
include("toolbox_visualization.jl")

export visualize_tikz
export visualize_mode

export plot_data

export create_gif_from_df
export create_gif_from_vector
export create_gif_2D

export plot2d_maculate_regions
export plot3d_maculate_regions
export visualize_projection_of_immaculate_linebundles_for_pic3

export visualize3d_sequence_by_layers
export visualize3d_sequence_moving_lines
export visualize3d_sequence_process_by_layers

export visualize_infection_of_pentagon
export plot_infection_step_2D
export plot_infection_procedure


########################################
########################################
include("toolbox_sympy.jl")
# add function if it works when pushing on github


end