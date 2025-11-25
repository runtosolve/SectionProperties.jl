using  RackSections, SectionProperties

E = 29500.0
ν = 0.30

#' ### Inputs

member_type = "column"
section_type = "closed_tube"
section_info = "with_holes_in_vertical"
H = 3.0 
D = 3.0
L = 0.75
R = 0.125 + 0.100
t = 0.100
dh_H = 0.710
dh_D = 0.531
de_H  = H/2 - 0.700
de_D = 0.875
hole_pitch_H = 2.0
hole_pitch_D = 2.0
hole_length_H = 1.086
hole_length_D = 0.531






section_inputs = RackSections.Columns.RectangularTubeInput(H, D, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

section = RackSections.Columns.rectangular_tube(section_inputs)


# coord = [section.geometry.x section.geometry.y]

# coord = vcat(coord, coord[1, :]')

# ends = [section.properties.element_info[:, 1] section.properties.element_info[:, 2] section.properties.element_info[:, 3]]
# ends[end, 2] = ends[end, 1] + 1



# properties = SectionProperties.cutwp_prop2(coord,ends)

