module SectionProperties


using CUFSM, LinearAlgebra, Statistics

# using ..Geometry

struct PlasticSectionProperties

    neutral_axis_location::Float64
    Z::Float64

end



function open_thin_walled(center, t)

    X_c = [center[i][1] for i in eachindex(center)]
    Y_c = [center[i][2] for i in eachindex(center)]

    #define elements for section property calc
    num_nodes = length(X_c)
    num_elem = num_nodes - 1
    coord = [X_c Y_c]
    ends = [1:num_nodes-1 2:num_nodes t]

    #calculate section properties
    section_properties = CUFSM.cutwp_prop2(coord,ends)

    return section_properties

end


function closed_thin_walled(center, t)

    X_c = [center[i][1] for i in eachindex(center)]
    Y_c = [center[i][2] for i in eachindex(center)]

    num_nodes = length(X_c)
    coord = [X_c Y_c]
    start_nodes = 1:num_nodes
    end_nodes = [2:num_nodes; 1]
    ends = [start_nodes end_nodes t]

    section_properties = CUFSM.cutwp_prop2(coord,ends)

    return section_properties

end



function calculate_plastic_neutral_axis_location(element_connectivity, element_thicknesses, node_geometry, about_axis)

    if about_axis == "x"

        axis_iterator = sort(unique(node_geometry[:, 2]))

    elseif about_axis == "y"

        axis_iterator = sort(unique(node_geometry[:, 1]))

    end

    num_axis_points = length(axis_iterator)

    sum_A_elements = Array{Float64}(undef, num_axis_points)

    for i=1:num_axis_points

        axis_location = axis_iterator[i]

        #Define an array of all the element areas.
        A_elements = calculate_axis_area(element_connectivity, element_thicknesses, node_geometry, axis_location, about_axis)

        sum_A_elements[i] = sum(A_elements)

    end

    #Find sign switch on area sum.
    plastic_neutral_axis_index = findfirst(x->x<0, sum_A_elements)

    #Decide if last positive or first negative is a better solution.
    if abs(sum_A_elements[plastic_neutral_axis_index]) > abs(sum_A_elements[plastic_neutral_axis_index-1])
        plastic_neutral_axis_index = plastic_neutral_axis_index - 1
    end

    plastic_neutral_axis_location = axis_iterator[plastic_neutral_axis_index]

    return plastic_neutral_axis_location

end


function calculate_plastic_section_modulus(node_geometry, element_connectivity, element_thicknesses, plastic_neutral_axis_location, about_axis)

    num_elem = size(element_connectivity)[1]

    Z_i = zeros(Float64, num_elem)

    for i=1:num_elem

        node_i = trunc(Int, element_connectivity[i, 1])
        node_j = trunc(Int, element_connectivity[i, 2])

        x_i = node_geometry[node_i, 1]
        x_j = node_geometry[node_j, 1]
        y_i = node_geometry[node_i, 2]
        y_j = node_geometry[node_j, 2]

        length_i = norm([x_j, y_j] - [x_i, y_i])
        A_i = length_i * element_thicknesses[i]

        if about_axis == "y"

            z_i = mean([x_i, x_j])

        elseif about_axis == "x"

            z_i = mean([y_i, y_j])

        end

        Z_i[i] = A_i * abs(z_i - plastic_neutral_axis_location)

    end

    Z = sum(Z_i)

    return Z

end

function calculate_plastic_section_properties(node_geometry, element_definitions, about_axis)

    #Assign cross-section discretization info.
    element_connectivity = element_definitions[:, 1:2]
    element_thicknesses = element_definitions[:, 3]

    #Approximate the plastic neutral axis location about the defined axis.
    plastic_neutral_axis_location = calculate_plastic_neutral_axis_location(element_connectivity, element_thicknesses, node_geometry, about_axis)

    #Calculate the plastic modulus.
    Z = calculate_plastic_section_modulus(node_geometry, element_connectivity, element_thicknesses, plastic_neutral_axis_location, about_axis)

    #Add these properties to an object.
    plastic_section_properties = PlasticSectionProperties(plastic_neutral_axis_location, Z)

    return plastic_section_properties

end

function circular_tube(outside_diameter, t, num_radial_segments)

    r = outside_diameter/2
    α = range(0.0, 2π, num_radial_segments)[1:end-1]
    x = (r - t/2) * sin.(α)  
    y = (r - t/2) * cos.(α) 

    y = y .- minimum(y) .+ t/2
    x = x .- minimum(x) .+ t/2

    center = [[x[i], y[i]] for i in eachindex(x)]

    properties = closed_thin_walled(center, ones(Float64, length(x)) .* t)

    return properties

end



function calculate_axis_area(element_connectivity, element_thicknesses, node_geometry, axis_location, about_axis)

    num_elem = size(element_connectivity)[1]

    A_elements = zeros(Float64, num_elem)

    for i=1:num_elem

        node_i = trunc(Int, element_connectivity[i, 1])
        node_j = trunc(Int, element_connectivity[i, 2])

        x_i = node_geometry[node_i, 1]
        x_j = node_geometry[node_j, 1]
        y_i = node_geometry[node_i, 2]
        y_j = node_geometry[node_j, 2]

        length_i = norm([x_j, y_j] - [x_i, y_i])

        A_i = length_i * element_thicknesses[i]

        if about_axis == "y"

            z_i = mean([x_i, x_j])

        elseif about_axis == "x"

            z_i = mean([y_i, y_j])

        end

        if (z_i - axis_location) > 0  #Find which side of the axis A_i is on.

            A_elements[i] = A_i

        elseif (z_i - axis_location) < 0

            A_elements[i] = -A_i
            
        end

    end

    return A_elements

end



end #module



