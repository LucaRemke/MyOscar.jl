module MyOscar

##########
include("toolbox.jl")

export define_projective_space
export define_projective_surface
export define_hirzebruch_surface
export define_ppp
export define_pentagon
export define_hexagon

export get_representative_of_divisor
export transform_rayvector
export generate_vectors
export convert_to_df
export sort_points_clockwise

export show_generators_and_relations_of_classgroup
export calculate_cartier_data
export check_nef_conditions
export calculate_nef_pair
export get_immaculate_lb
export get_maculate_regions
export cutout_maculate_regions
export plot_maculate_regions
export plot3d_maculate_regions
export print_maculate_region_info

export plot_data
export create_gif

export get_temptings_via_bruteforce
export search_nontemptings_via_boundary_condition
export eliminate_faces_as_temptings
export get_temptings_as_primitive_collections

export get_images_of_cube_vertices
export get_polytope_of_cube_vertices

export get_seed_and_chull_of_toricpic3


##########
include("toolbox_mes.jl")

export generate_vectors_from_ranges
export is_in_locus_of_pentagon
export advanced_extend_sequences
export advanced_brute_force_exseq_for_pentagon_in_zero
export group_by_type_sequence
export plot_type
export plot_type_allseq


##########
include("toolbox_exceptional_sequences.jl")

export calculate_augmentations_of_mes
export calculate_helex_of_mes
export calculate_all_helex_of_mes
export has_duplicates
export calculate_inverted_mes
export search_duplicates




#include("sympy.jl")

#export evaluate_cartier_data

end
