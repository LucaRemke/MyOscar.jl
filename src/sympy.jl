# This is for all functions using sympy

using PyCall
pyimport("sympy")
using SymPy;


####################
# calculates the Cartier-Data of general Cartier-Divisor
function calculate_cartier_data(coeffs::Vector{Sym{PyObject}}, var::NormalToricVariety, cones::Bool=true)
    # implement a check if divisor is not cartier
    rays_var = rays(var)
    maxcones_var = maximal_cones(polyhedral_fan(var))

    rays_var_int = []
    cartier_data = []

    for r in rays_var
        push!(rays_var_int, transform_rayvector(r))    
    end

    rays_matrix = transpose(hcat(rays_var_int...))    

    for c in maxcones_var
        rays_c = rays(c)        
        ind = Int64[]

        for r in rays_c
            i = findfirst(item -> item == r, rays_var)
            push!(ind,i)
        end

        A = rays_matrix[ind,:]
        b = -coeffs[ind]
        sol = A\b

        if cones == true
            out = c
        else
            out = rays(c)
        end

        push!(cartier_data,(out, sol)) 

    end
    
    return cartier_data
end;

####################
# evaluates the Cartier-Data of general Cartier-Divisor for a specific choice of coordinates
function evaluate_cartier_data(data::Vector{Any}, symbols::Vector{Sym}, coeffs::Vector{Int64})
    if length(symbols) != length(coeffs)
        throw(ArgumentError("Length of the input vectors does not match.")) 
    end
    
    all_data = []
    
    for i in 1:length(data)
        cdata = data[i][2]
        
        for j in 1:length(symbols)
            cdata = cdata.subs(symbols[j],coeffs[j])
        end
        
        cdata = vec(convert(Matrix{Int64}, N(cdata))) 
        
        push!(all_data, cdata)
        
    end
    
    return all_data

end;