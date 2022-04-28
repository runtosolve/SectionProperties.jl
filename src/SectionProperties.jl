module SectionProperties

using Statistics
using LinearAlgebra
using CSV
using DataFrames
# using TriangleMesh
using StaticArrays
using Unitful, UnitfulUS

using Geometry

export AISC, wshape_nodes,
       assemble, Feature, Deck, surface_normals, avg_node_normals, xycoords_along_normal, create_CUFSM_node_elem, feature_geometry, get_xy_coordinates, area_from_cells, centroid_from_cells, moment_of_inertia_from_cells,
       triangular_mesh_properties, mesh, SectionProperties, rectangular_tube_geometry, define_w_shape_centerline_model, discretize_w_shape_centerline_model


struct WShape

    d::Float64
    tw::Float64
    bf::Float64
    tf::Float64
    kdes::Float64
    k1::Float64
    A::Float64
    Ix::Float64
    Iy::Float64
    J::Float64
    Cw::Float64
    Zx::Float64
    Zy::Float64
    Wno::Float64


end

struct CrossSectionBranch

    anchor::Tuple{Float64, Float64}
    direction::Float64
    magnitude::Float64
    n::Int64

end

#primitive, line element defined as vector, with n segments
struct Feature
    ΔL::Union{Array{Float64,1}, Vector{<:Unitful.Length}}
    θ::Array{Float64,1}
    n::Array{Int,1}
    radius::Union{Array{Float64,1}, Vector{<:Unitful.Length}}
    n_radius::Array{Int,1}
    closed_or_open::Int

end

struct Point

    x::Union{Float64, <:Unitful.Length}
    y::Union{Float64, <:Unitful.Length}

end

#deck cross-section definition
#this could be a model for other cross-section types
struct Deck

    features::Tuple{Feature,Feature,Feature}
    feature_map::Array{Int,1}

end


struct Open

    features::Feature
    feature_map::Array{Int,1}

end

struct PlasticSectionProperties

    neutral_axis_location::Float64
    Z::Float64

end

#Get the node and element properties just for the flange and lip of a C or Z section.
#This code grabs the bottom flange and lip.
# function CZflange_template(CorZ,H,Bc,Bt,Dc,Dt,r1,r2,r3,r4,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)

#     prop,node,elem,lengths,springs,constraints,geom,cz = CrossSection.CUFSMtemplate(CorZ,H,Bc,Bt,Dc,Dt,r1,r2,r3,r4,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)

#     if CorZ == 2
#         node[:, 2] = -node[:, 2]
#     end

#     index_xo = findall(x-> x==0.0, node[:, 2])
#     index_yo = findall(y-> y==0.0, node[:, 3])
#     index_o = intersect(index_xo, index_yo)
#     index = 1:(index_o[1]+1)

# 	nodeflange = node[index,:]
#     elemflange = elem[index[1:end-1],:]

#     return nodeflange, elemflange

# end


"""
    AISC(shape_name)

Accepts `shape_name` as a String and returns cross-section dimensions and section properties in the Struct `shape_info`.

For now only W-shapes have been tested, eg., `AISC("W14x90")` where the struct contains `d, tw, bf, tf, kdes, k1, A, Ix, Iy, J, Cw, Zx, Wno`.

"""


function AISC(shape_name)

    filename = string(@__DIR__, "/assets/aisc-shapes-database-v15.0.csv")

    data = CSV.File(filename)

    shape_row = findfirst(==(shape_name), data.AISC_Manual_Label)
    
    section_type = data.Type[shape_row]
    
    if section_type == "W"
    
        d = parse(Float64, data.d[shape_row])
        tw = parse(Float64, data.tw[shape_row])
        bf = parse(Float64, data.bf[shape_row])
        tf = parse(Float64, data.tf[shape_row])
        kdes = parse(Float64, data.kdes[shape_row])
    
        #get k1 from AISC table, it is in fraction format
        k1 = data.k1[shape_row]
        index = findfirst("/", k1)
    
        if length(k1) == 7 
            whole = Int(data.k1[shape_row][1] - '0')
        else
            whole = 0.0
        end
    
        if isempty(index) == false
            top_fraction = parse(Float64, data.k1[shape_row][index[1]-2:index[1]-1])
            bottom_fraction = parse(Float64, data.k1[shape_row][index[1]+1:index[1]+2])
        else
            top_fraction = 0.0
            bottom_fraction = 0.0
        end
    
        k1 = whole + top_fraction/bottom_fraction
    
        A = data.A[shape_row]
        Ix = data.Ix[shape_row]
        Iy = data.Iy[shape_row]
        J = parse(Float64, data.J[shape_row])
        Cw = parse(Float64, data.Cw[shape_row])
        Zx = data.Zx[shape_row]
        Zy = data.Zy[shape_row]
        Wno = parse(Float64, data.Wno[shape_row])
    
        shape_info = WShape(d, tw, bf, tf, kdes,k1, A, Ix, Iy, J, Cw, Zx, Zy, Wno)
    
    end

    return shape_info

end

"""
    wshape_nodes(shape_info, n)

Accepts the Struct `shape_info` generated using CrossSection.AISC and the discretization Vector `n` and outputs the outline x-y coordinates of a W shape 'xcoords' and 'ycoords'.

The Vector 'n' describes the number of segments in a quarter cross-section, i.e., `n = [half of outside flange face, flange thickness, half of inside flange face, flange-web radius, half of web]`.

"""


function wshape_nodes(shape_info, n)

    #from bottom of bottom flange, web centerline to left edge
    xcoords = zeros(n[1]+1)
    ycoords = zeros(n[1]+1)

    flange_range = 0.0 : -shape_info.bf / 2 / n[1] : -shape_info.bf / 2
    [xcoords[i] =  flange_range[i] for i in eachindex(flange_range)]
    ycoords .= 0.0

    #up along bottom flange thickness
    flange_thickness_range = shape_info.tf/n[2]:shape_info.tf/n[2]:shape_info.tf
    xcoords = [xcoords; ones(n[2])*xcoords[end]]
    ycoords = [ycoords; flange_thickness_range]

    #over to fillet radius at bottom flange - web intersection

    # flange_flat = shape_info.bf/2 - shape_info.k1
    flange_flat = shape_info.bf/2 - shape_info.tw/2 - (shape_info.kdes - shape_info.tf)

    inside_flange_range = (xcoords[end] + flange_flat/n[3]) : flange_flat/n[3] : (xcoords[end] + flange_flat)

    xcoords = [xcoords; inside_flange_range]
    ycoords = [ycoords; ones(n[3])*ycoords[end]]

    #go around the fillet
    radius = -xcoords[end] - shape_info.tw/2
    θ = (-π/2 + π/2/n[4]):π/2/n[4]: 0.0

    xo = xcoords[end]
    yo = ycoords[end] + radius

    x_radius = xo .+ radius .* cos.(θ)
    y_radius = yo .+ radius .* sin.(θ)

    # plot(x_radius, y_radius, markershape = :o, linetype = :scatter)

    xcoords = [xcoords; x_radius]
    ycoords = [ycoords; y_radius]

    #add web flat
    web_flat = shape_info.d/2 - shape_info.tf - radius

    web_flat_range = LinRange(ycoords[end] + web_flat/n[5], (ycoords[end] + web_flat), n[5])
    # web_flat_range = (ycoords[end] + web_flat/n[5]): web_flat/n[5]: (ycoords[end] + web_flat)
    xcoords = [xcoords; ones(n[5])*xcoords[end]]
    ycoords = [ycoords; web_flat_range]

    #mirror about horizontal axis
    ycoords_horz_flip = ycoords .- ycoords[end]
    ycoords_horz_flip = -ycoords_horz_flip
    ycoords_horz_flip = ycoords_horz_flip .+ ycoords[end]

    xcoords = [xcoords; reverse(xcoords)[2:end]]
    ycoords = [ycoords; reverse(ycoords_horz_flip)[2:end]]

    #mirror about vertical axis
    xcoords_vert_flip = reverse(-xcoords)[2:end-1]

    xcoords = [xcoords; xcoords_vert_flip]
    ycoords = [ycoords; reverse(ycoords)[2:end-1]]


    return xcoords, ycoords

end


# #define nodal coordinates within a feature
# function discretize_feature(feature)

#     dx = feature.Δx ./ feature.n
#     dy = feature.Δy ./ feature.n

#     return dx, dy

# end

#define feature geometry, typically a feature is repetitive
function feature_geometry(feature, dx, dy)

    xcoords = []
    ycoords = []

    num_lines = length(feature.Δx)

    for i = 1:num_lines

        if i==1
            xcoords = range(0.0, feature.Δx[i], length = feature.n[i] + 1)
            ycoords = range(0.0, feature.Δy[i], length = feature.n[i] + 1)
        else
            xcoords = [xcoords; xcoords[end] .+ range(dx[i], feature.Δx[i], length = feature.n[i])]
            ycoords = [ycoords; ycoords[end] .+ range(dy[i], feature.Δy[i], length = feature.n[i])]
        end

    end

    return xcoords, ycoords

end


#assemble an open cross-section as a series of line element features
function assemble(OpenSection)

    xcoords = []
    ycoords = []
    num_features = length(OpenSection.feature_map)

    for i =1:num_features

        feature_type = OpenSection.feature_map[i]

        # dx, dy = discretize_feature(OpenSection.features[feature_type])

        if i==1
            xcoords, ycoords = get_xy_coordinates(OpenSection.features[feature_type])
        else
            feature_xcoords, feature_ycoords = get_xy_coordinates(OpenSection.features[feature_type])
            xcoords = [xcoords; xcoords[end] .+ feature_xcoords[2:end]]
            ycoords = [ycoords; ycoords[end] .+ feature_ycoords[2:end]]
        end

    end

    return xcoords, ycoords

end








function create_CUFSM_node_elem(xcoords, ycoords, connectivity, t)

    num_elem = size(connectivity)[1]
    num_nodes = length(xcoords)

    node = zeros((num_nodes, 8))
    elem = zeros((num_elem, 5))

    for i=1:num_nodes
        node[i,:] = [i, xcoords[i], ycoords[i], 1, 1, 1, 1, 1.0]
    end

    for i=1:num_elem
        elem[i, :] = [i, connectivity[i,1], connectivity[i,2], t[i], 100]
    end

    return node, elem

end


function get_xy_coordinates(feature)

    xcoords = []
    ycoords = []

    #convert feature vectors to xy components
    Δxy = Geometry.vector_components(feature.ΔL, feature.θ)

    #number of straight line segments in the feature
    num_lines = size(Δxy)[1]

    #get xy line coordinates of feature
    xcoords_straight, ycoords_straight = Geometry.line_coordinates(Δxy, feature.closed_or_open)

    #calculate feature surface normals
    unitnormals = surface_normals(xcoords_straight, ycoords_straight, feature.closed_or_open)

    #calculate feature node normals
    nodenormals = avg_node_normals(unitnormals, feature.closed_or_open)

    #calculate interior corner angle magnitudes
    num_radius = length(feature.radius)
    interior_radius_angle = zeros(Float64, num_radius)
    for i=1:num_radius

        θ1 = feature.θ[i]
        θ2 = feature.θ[i+1]

        if sign(θ1) != sign(θ2)   #don't think this is fully general
            interior_radius_angle[i] = abs(180.0 - (abs(θ1) + abs(θ2)))
        else
            interior_radius_angle[i] = θ2 - θ1
        end

    end

    #calculate distance from curve PI to start and end of curve, along tangents
    # Tc = zeros(Float64, length(feature.radius))
    Tc = []

    for i = 1:length(feature.radius)

        radius = feature.radius[i]

        if feature.radius[i] > 0.0 * feature.radius[i]

            xy_PI = [xcoords_straight[i+1], ycoords_straight[i+1]]
            n = feature.n_radius[i]
            γ = interior_radius_angle[i]
            θ1 = feature.θ[i]
            θ2 = feature.θ[i+1]

            if (sign(θ1) <=0) & (sign(θ2) >= 0)
                PI_unit_normal = -nodenormals[i+1,:]
            else
                PI_unit_normal = nodenormals[i+1,:]
            end

            xy_curve, Δ, E, Tc_i, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(feature.radius[i], γ, xy_PI, PI_unit_normal, n)

        elseif feature.radius[i] == 0.0 * feature.radius[i]

            Tc_i = 0.0 * feature.radius[i]
        
        end

        if i == 1

            Tc = Tc_i

        else

            Tc = [Tc; Tc_i]

        end

    end

    #shorten vector lengths to include curves
    # ΔLc = zeros(Float64, num_lines)
    ΔLc = []
    for i = 1:length(feature.ΔL)

        if i == 1  #first
            ΔLc = feature.ΔL[i] - Tc[i]   #j end
        elseif i == length(feature.ΔL)  #last
            ΔLc = [ΔLc; (feature.ΔL[i] - Tc[end])]  #j end
        else  #others
            ΔLc = [ΔLc; (feature.ΔL[i] - Tc[i-1] - Tc[i])] #i and j ends
        end

    end

    Δxyc = Geometry.vector_components(ΔLc, feature.θ)

    #discretize feature
    dx = Δxyc[:,1] ./ feature.n
    dy = Δxyc[:,2] ./ feature.n

    #assemble feature
    for i = 1:num_lines

        if i==1   #first vector
            xcoords = range(0.0 * Δxyc[i,1], Δxyc[i,1], length = feature.n[i] + 1)
            ycoords = range(0.0 * Δxyc[i,2], Δxyc[i,2], length = feature.n[i] + 1)
        elseif i != num_lines
            xcoords = [xcoords; xcoords[end] .+ range(dx[i], Δxyc[i,1], length = feature.n[i])]
            ycoords = [ycoords; ycoords[end] .+ range(dy[i], Δxyc[i,2], length = feature.n[i])]
        end

        #add radius to end of vector
        if i < num_lines
            if feature.radius[i] > 0.0 * feature.radius[i]
                xy_PI = [xcoords_straight[i+1], ycoords_straight[i+1]]
                n = feature.n_radius[i]
                radius = feature.radius[i]
                γ = interior_radius_angle[i]
                θ1 = feature.θ[i]
                θ2 = feature.θ[i+1]

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    PI_unit_normal = -nodenormals[i+1,:]  #concave up
                
                # elseif (θ2 > θ1) | (θ1 > θ2)

                #     PI_unit_normal = -nodenormals[i+1,:] #added for Cee

                else
                    PI_unit_normal = nodenormals[i+1,:]  #concave down
                end

                xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    xy_curve = reverse(xy_curve, dims=1)    #concave up
                end

                xcoords = [xcoords; xy_curve[2:end,1]]
                ycoords = [ycoords; xy_curve[2:end,2]]

            end
        end

        if i== num_lines  #last vector
            xcoords = [xcoords; xcoords[end] .+ range(dx[i], Δxyc[i,1], length = feature.n[i])]
            ycoords = [ycoords; ycoords[end] .+ range(dy[i], Δxyc[i,2], length = feature.n[i])]
        end

    end

    return xcoords, ycoords

end



#calculate cross-sectional area from cells

function area_from_cells(Ai)

    A = sum(Ai)

end

#calculate cross-section centroid from cells

function centroid_from_cells(Ai, cxi, cyi)

    A = area_from_cells(Ai)

    cx = sum(cxi .* Ai) / A
    cy = sum(cyi .* Ai) / A

    return cx, cy

end


function moment_of_inertia_from_cells(Ai, ci, c)

    I = sum(((c .- ci) .^2 .* Ai))

    return I

end

#removed because of problems installing TriangleMesh.jl on Windows computers.
# #discretize cross-section with a triangular mesh
# function triangular_mesh(xcoords, ycoords, mesh_size)

#     num_nodes = length(xcoords)
#     num_segments = num_nodes

#     # n_point, n_point_marker, n_point_attribute, n_segment, n_holes
#     poly = TriangleMesh.Polygon_pslg(num_nodes, 1, 0, num_segments, 0)

#     node = [xcoords ycoords]
#     set_polygon_point!(poly, node)

#     node_marker = ones(Int, num_nodes, 1)
#     set_polygon_point_marker!(poly, node_marker)

#     segments = zeros(Int, num_segments, 2)
#     for i=1:num_segments

#         if i == num_segments
#             segments[i, 1:2] = [i, 1]
#         else
#             segments[i, 1:2] = [i, i+1]
#         end

#     end

#     set_polygon_segment!(poly, segments)

#     segment_markers = ones(Int, num_segments)
#     set_polygon_segment_marker!(poly, segment_markers)

#     #switches from https://www.cs.cmu.edu/~quake/triangle.html
#     switches = "penvVa" * string(mesh_size) * "D"

#     mesh = create_mesh(poly, switches)

#     return mesh

# end


# function triangular_mesh_properties(mesh)

#     #calculate cell area and centroid
#     Ai = zeros(Float64, mesh.n_cell)
#     cxi = zeros(Float64, mesh.n_cell)
#     cyi = zeros(Float64, mesh.n_cell)

#     for i = 1:mesh.n_cell

#         p1 = mesh.cell[1, i]
#         p2 = mesh.cell[2, i]
#         p3 = mesh.cell[3, i]

#         x1 = mesh.point[1, p1]
#         y1 = mesh.point[2,p1]
#         x2 = mesh.point[1, p2]
#         y2 = mesh.point[2,p2]
#         x3 = mesh.point[1, p3]
#         y3 = mesh.point[2,p3]

#     Ai[i] = Geometry.triangle_area(x1, y1, x2, y2, x3, y3)

#     cxi[i], cyi[i] = Geometry.triangle_centroid(x1, y1, x2, y2, x3, y3)

#     end

#     return Ai, cxi, cyi

# end



# function mesh(xcoords, ycoords, mesh_size)

#     section_mesh = CrossSection.triangular_mesh(xcoords, ycoords, mesh_size)
#     Ai, cxi, cyi = CrossSection.triangular_mesh_properties(section_mesh)

#     return Ai, cxi, cyi

#  end


 function rectangular_tube_geometry(feature)


    xcoords = []
    ycoords = []

    #convert feature vectors to xy components
    Δxy = Geometry.vector_components(feature.ΔL, feature.θ)

    #number of straight line segments in the feature
    num_lines = size(Δxy)[1]

    #get xy line coordinates of feature
    xcoords_straight, ycoords_straight = Geometry.line_coordinates(Δxy, feature.closed_or_open)

    #calculate feature surface normals
    unitnormals = surface_normals(xcoords_straight, ycoords_straight, feature.closed_or_open)

    #calculate feature node normals
    nodenormals = avg_node_normals(unitnormals, feature.closed_or_open)

    #calculate interior corner angle magnitudes
    num_radius = length(feature.radius)
    interior_radius_angle = zeros(Float64, num_radius)

    for i=1:num_radius

        θ1 = feature.θ[i] - 180

        if (feature.closed_or_open == 0) & (i==num_radius)  #for closed section, return to beginning node
            θ2 = feature.θ[1]
        else
            θ2 = feature.θ[i+1]
        end

        unit_vector_i = [cos(deg2rad(θ1)), sin(deg2rad(θ1))]
        unit_vector_j = [cos(deg2rad(θ2)), sin(deg2rad(θ2))]

        interior_radius_angle[i] = rad2deg(acos(dot(unit_vector_i, unit_vector_j)))

    end

    #calculate distance from curve PI to start and end of curve, along tangents
    Tc = zeros(Float64, length(feature.radius))

    for i = 1:length(feature.radius)

        radius = feature.radius[i]

        if radius > 0.0


            if feature.closed_or_open == 0
                xy_PI = [xcoords_straight[i], ycoords_straight[i]]
            elseif feature.closed_or_open == 1
                xy_PI = [xcoords_straight[i+1], ycoords_straight[i+1]]
            end

            n = feature.n_radius[i]
            γ = interior_radius_angle[i]


            θ1 = feature.θ[i]


            if (feature.closed_or_open == 0) & (i == length(feature.radius))  #for closed section, return to beginning node
                θ2 = feature.θ[1]
            else
                θ2 = feature.θ[i+1]
            end


            if (feature.closed_or_open == 0) & (i == length(feature.radius))

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    PI_unit_normal = -nodenormals[1,:]
                else
                    PI_unit_normal = nodenormals[1,:]
                end

            else

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    PI_unit_normal = -nodenormals[i+1,:]
                else
                    PI_unit_normal = nodenormals[i+1,:]
                end

            end


            xy_curve, Δ, E, Tc[i], xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

        elseif radius == 0.0
            Tc[i] = 0.0
        end

    end

    #Shorten vector lengths to include curves.
    ΔLc = zeros(Float64, num_lines)

    for i = 1:length(feature.ΔL)

        if i == 1  #first
            if feature.closed_or_open == 1
                ΔLc[i] = feature.ΔL[i] - Tc[i]   #j end
            elseif feature.closed_or_open == 0
                ΔLc[i] = feature.ΔL[i] - Tc[i] - Tc[end]
            end
        elseif i == length(feature.ΔL)  #last
            if feature.closed_or_open == 1
                ΔLc[i] = feature.ΔL[i] - Tc[end]  #j end
            elseif feature.closed_or_open == 0
                ΔLc[i] = feature.ΔL[i] - Tc[end] - Tc[i - 1]
            end
        else  #others
            ΔLc[i] = feature.ΔL[i] - Tc[i-1] - Tc[i] #i and j ends
        end

    end

    #Get xy components of flat lengths around tube.
    Δxyc = Geometry.vector_components(ΔLc, feature.θ)


    #Discretize feature.
    dx = Δxyc[:,1] ./ feature.n
    dy = Δxyc[:,2] ./ feature.n

    #Define vertical tube segment at y=0.

    ycoords = Tc[1]:dy[1]:(Δxyc[1,2] + Tc[1])
	xcoords = zeros(Float64, length(ycoords))

    #Define top left corner.

    n = feature.n_radius[1]
    radius = feature.radius[1]
    γ = interior_radius_angle[1]
    PI_unit_normal = nodenormals[2,:]
    xy_PI = [xcoords_straight[2], ycoords_straight[2]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:end, 1]]
    ycoords = [ycoords; xy_curve[2:end, 2]]

    #Define top flat.

    x_range = xcoords[end]:dx[2]:xcoords[end]+Δxyc[2,1]
    y_range = ycoords[end] * ones(length(x_range))
    xcoords = [xcoords; x_range[2:end]]
    ycoords = [ycoords; y_range[2:end]]

    #Define top right corner.

    n = feature.n_radius[2]
    radius = feature.radius[2]
    γ = interior_radius_angle[2]
    PI_unit_normal = nodenormals[3,:]
    xy_PI = [xcoords_straight[3], ycoords_straight[3]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:end, 1]]
    ycoords = [ycoords; xy_curve[2:end, 2]]

    #Define right vertical flat.

    y_range = ycoords[end]:dy[3]:(ycoords[end] + Δxyc[3,2])
    x_range = xcoords[end] * ones(length(y_range))
    xcoords = [xcoords; x_range[2:end]]
    ycoords = [ycoords; y_range[2:end]]

    #Define bottom right corner.

    n = feature.n_radius[3]
    radius = feature.radius[3]
    γ = interior_radius_angle[3]
    PI_unit_normal = nodenormals[4,:]
    xy_PI = [xcoords_straight[4], ycoords_straight[4]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:end, 1]]
    ycoords = [ycoords; xy_curve[2:end, 2]]

    #Define bottom horizontal flat.

    x_range = xcoords[end]:dx[4]:xcoords[end]+Δxyc[4,1]
    y_range = ycoords[end] * ones(length(x_range))
    xcoords = [xcoords; x_range[2:end]]
    ycoords = [ycoords; y_range[2:end]]

    #Define lower left corner.

    n = feature.n_radius[4]
    radius = feature.radius[4]
    γ = interior_radius_angle[4]
    PI_unit_normal = nodenormals[1,:]
    xy_PI = [xcoords_straight[1], ycoords_straight[1]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:(end-1), 1]]
    ycoords = [ycoords; xy_curve[2:(end-1), 2]]

    return xcoords, ycoords

end



function discretize_w_shape_centerline_model(shape, cross_section)

    num_branches = size(shape)[1]

    xcoords = []
    ycoords = []

    for i = 1:num_branches

        ΔL = shape[i].magnitude
        θ = shape[i].direction
        n = shape[i].n
        Δxy = Geometry.vector_components(ΔL, θ) 
        anchor = shape[i].anchor

        if i == 1
            xcoords = range(0.0, Δxy[1], length = n + 1) .+ anchor[1]
            ycoords = range(0.0, Δxy[2], length = n + 1) .+ anchor[2]

        else
            xcoords = [xcoords; range(0.0, Δxy[1], length = n + 1) .+ anchor[1]]
            ycoords = [ycoords; range(0.0, Δxy[2], length = n + 1) .+ anchor[2]]

        end

    end

    #Round here to help unique function.
    xycoords = [(round(xcoords[i], digits = 3), round(ycoords[i], digits = 3)) for i = 1:length(xcoords)]

    xycoords = unique(xycoords)

    coord = [y[i] for y in xycoords, i in 1:2]

    #Shift coordinates so that web is centered on x=0.
    coord[:, 1] = coord[:, 1] .- cross_section.bf/2

    #Shift coordinates so that bottom fiber is at y=0.
    coord[:, 2] = coord[:, 2] .+ cross_section.tf/2


    #Define element connectivity.

    num_elem = sum([shape[i].n for i=1:num_branches])
    
    node_start = 1
    node_end = shape[1].n

    node_i = node_start:node_end
    node_j = node_i .+ 1

    node_start = floor(Int, shape[1].n/2)+1
    node_end = node_end + 2

    node_i = [node_i; node_start]
    node_j = [node_j; node_end]

    node_start = node_end
    node_end = node_end + shape[2].n - 2

    node_i = [node_i; node_start:node_end]
    node_j = [node_j; (node_start:node_end) .+ 1]

    node_start = shape[1].n + shape[2].n + 2
    node_end = node_start + floor(Int, shape[2].n/2) - 1
    node_i_range = range(node_start, node_end-1)
    node_j_range = node_i_range .+ 1

    node_i = [node_i; node_i_range]
    node_j = [node_j; node_j_range]

    node_start = node_i[end] + 1
    node_end = shape[1].n + shape[2].n + 1

    node_i = [node_i; node_start]
    node_j = [node_j; node_end]

    node_start = node_j[end]
    node_end = node_i[end] + 1

    node_i = [node_i; node_start]
    node_j = [node_j; node_end]

    node_start = shape[1].n + shape[2].n + 2 + floor(Int, shape[3].n/2)
    node_end = node_start + floor(Int, shape[3].n/2) - 1
    node_i_range = range(node_start, node_end-1)
    node_j_range = node_i_range .+ 1

    node_i = [node_i; node_i_range]
    node_j = [node_j; node_j_range]

    t = [ones(Float64, shape[1].n)*cross_section.tf[1]; ones(Float64, shape[2].n)*cross_section.tw[1]; ones(Float64, shape[2].n)*cross_section.tf[1]]

    ends = [node_i node_j t]

    return coord, ends

end


function define_w_shape_centerline_model(bf, tf, d, n)

    num_branches = 3
    w_shape = Vector{CrossSectionBranch}(undef, num_branches)

    #first branch, bottom flange
    anchor = (0.0, 0.0)
    direction = 0.0
    magnitude = bf

    w_shape[1] = CrossSectionBranch(anchor, direction, magnitude, n[1])

    #second branch, web
    anchor = (bf/2, 0.0)
    direction = 90.0
    magnitude = d - tf

    w_shape[2] = CrossSectionBranch(anchor, direction, magnitude, n[2])

    #third branch, top flange
    anchor = (0.0, d - tf)
    direction = 0.0
    magnitude = bf

    w_shape[3] = CrossSectionBranch(anchor, direction, magnitude, n[3])

    return w_shape

end




function insert_cross_section_node(node_geometry, element_connectivity, element_thicknesses, new_node_geometry)

    num_elem = size(element_connectivity)[1]
    num_nodes = size(node_geometry)[1]

    element_position_marker = Array{Float64}(undef, num_elem)

    #Define the new node.
    new_node = [new_node_geometry[1], new_node_geometry[2]]

    for i = 1:num_elem
        
        node_i = Int(element_connectivity[i, 1])
        node_j = Int(element_connectivity[i, 2])

        #direction vector from new node to node i
        point_i = node_geometry[node_i, :]
        vector_i = point_i - new_node

        #direction vector from new node to node j
        point_j = node_geometry[node_j, :]
        vector_j = point_j - new_node

        # Point on the element line, element line vector
        element_i = Geometry.Line(SA[node_geometry[node_i,1], node_geometry[node_i,2]], SA[node_geometry[node_j,1] - node_geometry[node_i,1], node_geometry[node_j,2] - node_geometry[node_i,2]])

        #If dot product is positive, than means that both nodes are ahead of the new node.  If dot product is negative then position vectors are pointing in opposite directions.
        element_position_marker[i] = vector_i   ⋅  vector_j  

    end

    #This is the element to split.
    split_element_index = findfirst(x->x<0, element_position_marker)

    if split_element_index==nothing   #for cases where node does not fall within the line segment of an element

        updated_node_geometry = deepcopy(node_geometry)
        updated_element_connectivity = deepcopy(element_connectivity)
        updated_element_thicknesses = deepcopy(element_thicknesses)

    else

        #Add new node to the end of the node geometry array.
        updated_node_geometry = [node_geometry; new_node_geometry]

        #Update the element definitions to include the new node.
        new_node_number = num_nodes + 1
        new_element = [new_node_number  element_connectivity[split_element_index, 2]]  #This is the second element in the split.
        updated_element_connectivity = deepcopy(element_connectivity)
        updated_element_connectivity[split_element_index, 2] = new_node_number  #Update the first element in the split.
        updated_element_connectivity = [updated_element_connectivity[1:split_element_index, :]; new_element; updated_element_connectivity[split_element_index+1:end, :]]
        
        #Update the element thickness array.
        new_element_thickness = element_thicknesses[split_element_index]
        updated_element_thicknesses = deepcopy(element_thicknesses)
        updated_element_thicknesses = [updated_element_thicknesses[1:split_element_index, :]; new_element_thickness; updated_element_thicknesses[split_element_index+1:end, :]]

    end

    return updated_node_geometry, updated_element_connectivity, updated_element_thicknesses

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


function generate_cross_section_outer_inner_mid_nodes(outer_surface_object, t, closed_or_open)

    if closed_or_open == 1

        #Calculate the outward (up) facing surface coordinates.
        xcoords_top, ycoords_top = get_xy_coordinates(outer_surface_object)

    elseif closed_or_open == 0   #only works for rectangular closed tubes right now

        #Calculate the out-to-out tube surface coordinates.
        xcoords_top, ycoords_top = rectangular_tube_geometry(outer_surface_object)

    end

    #Shift coordinates so that bottom of section is at y=0.
    ycoords_top = ycoords_top .- minimum(ycoords_top) .+ t


    #Calculate element normals.
    unitnormals = surface_normals(xcoords_top, ycoords_top, closed_or_open)
    nodenormals = avg_node_normals(unitnormals, closed_or_open)

    #Calculate the midline coordinates.
    xcoords_mid, ycoords_mid = xycoords_along_normal(xcoords_top, ycoords_top, nodenormals, -t/2)

    #Calculate inward (down) facing surface coordinates.
    xcoords_bottom, ycoords_bottom = xycoords_along_normal(xcoords_top, ycoords_top, nodenormals, -t)

    coords_top = (x=xcoords_top, y=ycoords_top)

    coords_mid = (x=xcoords_mid, y=ycoords_mid)

    coords_bottom = (x=xcoords_bottom, y=ycoords_bottom)

    return coords_top, coords_mid, coords_bottom

end






end #module


