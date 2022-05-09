module Tools

using Statistics
using LinearAlgebra


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


end