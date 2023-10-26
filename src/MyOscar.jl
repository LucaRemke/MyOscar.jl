module MyOscar

include("toolbox.jl")

export define_projective_space
export define_hirzebruch_surface
export define_ppp
export define_blowup_pp
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
export get_immaculate_regions

export plot_data
export create_gif

export get_temptings_via_bruteforce
export search_nontemptings_via_boundary_condition
export eliminate_faces_as_temptings
export get_temptings_as_primitive_collections

#include("sympy.jl")

#export evaluate_cartier_data

end
