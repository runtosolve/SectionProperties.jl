module SectionProperties


using LinearAlgebra, Statistics

# using ..Geometry

struct PlasticSectionProperties

    neutral_axis_location::Float64
    Z::Float64

end


mutable struct SectionPropertiesObject

    node_geometry::Array{Float64, 2}
    element_info::Array{Float64, 2}
    A::Float64
    xc::Float64
    yc::Float64
    Ixx::Float64
    Iyy::Float64
    Ixy::Float64
    θ::Float64
    I1::Float64
    I2::Float64
    J::Float64
    xs::Float64
    ys::Float64
    Cw::Float64
    B1::Float64
    B2::Float64
    wn::Array{Float64, 1}

end



function cutwp_prop2(coord,ends)

    #Translated from Matlab to Julia on October 5, 2020.
    
    #### Matlab notes below.
    
    #Function modified for use in CUFSM by Ben Schafer in 2004 with permission
    #of Sarawit. removed elastic buckling calcs & kept only section
    #properties
    #
    #August 2005: additional modifications; program only handles
    #singly-branched open sections; | single cell closed sections; arbitrary
    #section designation added for other types.
    #
    #December 2006 bug fixes to B1 B2
    #
    #December 2015 extended to not crash on disconnected & arbitrary sections
    #
    #  Compute cross section properties
    #----------------------------------------------------------------------------
    #  Written by:
    #       Andrew T. Sarawit
    #       last revised:   Wed 10/25/01
    #
    #  Function purpose:
    #       This function computes the cross section properties: area; centroid;
    #       moment of inertia; torsional constant; shear center; warping constant;
    #       B1; B2; elastic critical buckling load & the deformed buckling shape
    #
    #  Dictionary of Variables
    #     Input Information:
    #       coord[i,1:2]   .==  node i's coordinates
    #                            coord[i,1] = X coordinate
    #                            coord[i,2] = Y coordinate
    #       ends[i,1:2]    .==  subelement i's nodal information
    #                            ends[i,1] = start node #
    #                            ends[i,2] = finish node #
    #                            ends[i,3] = element's thickness
    #       KL1            .==  effective unbraced length for bending about the 1-axis()
    #       KL2            .==  effective unbraced length for bending about the 2-axis()
    #       KL3            .==  effective unbraced length for twisting about the 3-axis()
    #       force          .==  type of force applied
    #                      .== "Pe"  : elastic critical axial force
    #                      .== "Me1" : elastic critical moment about the 1-axis()
    #                      .== "Me2" : elastic critical moment about the 2-axis()
    #       exy[1:2]     .==  Pe eccentricities coordinates
    #                            exy[1] = ex
    #                            exy[2] = ey
    #  Output Information:
    #       A              .==  cross section area()
    #       xc             .==  X coordinate of the centroid from orgin
    #       yc             .==  Y coordinate of the centroid from origin
    #       Ix             .==  moment of inertia about centroid X axes()
    #       Iy             .==  moment of inertia about centroid Y axes()
    #       Ixy            .==  product of inertia about centroid
    #       Iz             .==  polar moment of inertia about centroid
    #       theta          .==  rotation angle for the principal axes()
    #       I1             .==  principal moment of inertia about centroid 1 axes()
    #       I2             .==  principal moment of inertia about centroid 2 axes()
    #       J              .==  torsional constant
    #       xs             .==  X coordinate of the shear center from origin
    #       ys             .==  Y coordinate of the shear center from origin
    #       Cw             .==  warping constant
    #       B1             .==  int[y*(x^2+y^2),s,0,L]   *BWS, x,y=prin. crd.
    #       B2             .==  int[x*(x^2+y^2),s,0,L]
    #                          where: x = x1+s/L*(x2-x1)
    #                                 y = y1+s/L*(y2-y1)
    #                                 L = lenght of the element
    #       Pe[i]          .==  buckling mode i's elastic critical buckling load()
    #       dcoord         .==  node i's coordinates of the deformed buckling shape
    #                            coord[i,1,mode] = X coordinate
    #                            coord[i,2,mode] = Y coordinate
    #                          where: mode = buckling mode number
    #
    #  Note:
    #     J;xs;ys;Cw;B1;B2;Pe;dcoord is not computed for close-section
    #
    #----------------------------------------------------------------------------
    #
    # find nele  .== total number of elements
    #      nnode .== total number of nodes
    #      j     .== total number of 2 element joints
    
    
        nele = size(ends,1);
        node = ends[:,1:2]; node = node[:]
        nnode = 0
        j = 0
    
        while isempty(node) == false
            i = findall(x-> x==node[1], node)
            deleteat!(node, i)
            # node[i] = []
            if size(i,1)==2
                j = j+1
            end
            nnode = nnode+1
        end
    
        # classify the section type()
        #This section modified in April 2006 by BWS to create an "arbitrary" category
        if j .== nele
            section = "close"; #single cell()
        elseif j .== nele-1
            section = "open";
        else
            section = "arbitrary"; #arbitrary section unidentified
                #in the future it would be good to handle multiple cross-sections in
                #one model; etc.; for now the code will bomb if more than a single()
                #section is used - due to inability to calculate section properties.
                #2015 decided to treat the secton as a fully composite section &
        end
    
    
        # if the section is close re-order the element
        if section == "close"
            xnele = (nele-1)
            for i = 1:xnele
                en = deepcopy(ends); en[i,2] = 0
                index = findall(ends[i,2] .== en[:,1:2])  #update this for Julia, need to deal with CartesianIndex
                m=index[1][1]
                n=index[1][2]
                if n==1
                    ends[i+1,:] = en[m,:]
                    ends[m,:] = en[i+1,:]
                elseif n == 2
                    ends[i+1,:] = en[m,[2 1 3]]
                    ends[m,:] = en[i+1,[2 1 3]]
                end
            end
        end
    
        # t = zeros(Float64, nele)
        # xm = zeros(Float64, nele)
        # ym = zeros(Float64, nele)
        # xd = zeros(Float64, nele)
        # yd = zeros(Float64, nele)
        # L = zeros(Float64, nele)
    
        t = Vector{Float64}(undef, nele)
        xm = Vector{Float64}(undef, nele)
        ym = Vector{Float64}(undef, nele)
        xd = Vector{Float64}(undef, nele)
        yd = Vector{Float64}(undef, nele)
        L = Vector{Float64}(undef, nele)
    
    
        # find the element properties
        for i = 1:nele
            sn = ends[i,1]; fn = ends[i,2];
    
            sn = convert(Int, sn)
            fn = convert(Int, fn)
    
            # thickness of the element
            t[i] = ends[i,3]
            # compute the coordinate of the mid point of the element
            xm[i] = mean(coord[[sn fn],1])
            ym[i] = mean(coord[[sn fn],2])
            # compute the dimension of the element
            xd[i] = diff(vec(coord[[sn fn],1]))[1]
            yd[i] = diff(vec(coord[[sn fn],2]))[1]
            # compute the length of the element
            L[i] = norm([xd[i] yd[i]])
        end
    
        # compute the cross section area()
        A = sum(L.*t)
        # compute the centroid
        xc = sum(L.*t.*xm)/A
        yc = sum(L.*t.*ym)/A
    
        if abs(xc/sqrt(A)) .< 1e-12
            xc = 0   
        end
        if abs(yc/sqrt(A)) .< 1e-12
            yc = 0   
        end
    
        # compute the moment of inertia
        Ix = sum((yd.^2/12 .+(ym .-yc).^2).*L.*t)
        Iy = sum((xd.^2/12 .+(xm .-xc).^2).*L.*t)
        Ixy = sum((xd.*yd/12 .+(xm .-xc).*(ym .-yc)).*L.*t)
    
        if abs(Ixy/A^2) .< 1e-12
            Ixy = 0 
        end
    
        # compute the rotation angle for the principal axes()
        theta = angle(Ix-Iy-2*Ixy*1im)/2
    
        # coord12 = zeros(Float64, size(coord))
        coord12 = Matrix{Float64}(undef, size(coord))
    
        # transfer the section coordinates to the centroid principal coordinates
        coord12[:,1] = coord[:,1] .-xc
        coord12[:,2] = coord[:,2] .-yc
        coord12 = [cos(theta) sin(theta); -sin(theta) cos(theta)]*transpose(coord12)
        coord12 = transpose(coord12)
    
        # find the element properties
        for i = 1:nele
            sn = ends[i,1]
            fn = ends[i,2]
    
            sn = convert(Int, sn)
            fn = convert(Int, fn)
    
            # compute the coordinate of the mid point of the element
            xm[i] = mean(coord12[[sn fn],1])
            ym[i] = mean(coord12[[sn fn],2])
            # compute the dimension of the element
            xd[i] = diff(vec(coord12[[sn fn],1]))[1]
            yd[i] = diff(vec(coord12[[sn fn],2]))[1]
        end
    
        # compute the principal moment of inertia
        I1 = sum((yd.^2/12+ym.^2).*L.*t)
        I2 = sum((xd.^2/12+xm.^2).*L.*t)
    
        if section == "close"
    
            p = zeros(Float64, nele)
            # compute the torsional constant for close-section
            for i = 1:nele
                sn = ends[i,1]
                fn = ends[i,2]
    
                sn = Int(sn)
                fn = Int(fn)
    
                p[i] = ((coord[sn,1]-xc)*(coord[fn,2]-yc)-(coord[fn,1]-xc)*(coord[sn,2]-yc))/L[i]
            end
            J = 4*sum(p.*L/2)^2/sum(L./t)
            xs = NaN; ys = NaN; Cw = NaN; B1 = NaN; B2 = NaN; Pe = NaN; dcoord = NaN; wn=[NaN; NaN]
        elseif section == "open"
            # compute the torsional constant for open-section
            J = sum(L.*t.^3)/3
    
            # compute the shear center & initialize variables
            nnode = size(coord,1)
            # w = zeros((nnode,2))
            w = Matrix{Float64}(undef, (nnode,2))
            w[:,1] .= 0.0
            w[:,2] .= 0.0 
            w[convert(Int, ends[1,1]),1] = ends[1,1]
            # wo = zeros((nnode,2))
            wo = Matrix{Float64}(undef, (nnode,2))
            wo[:,1] .= 0.0
            wo[:,2] .= 0.0 
            wo[convert(Int, ends[1,1]),1] = ends[1,1]
            Iwx = 0.0 ; Iwy = 0.0 ; wno = 0.0 ; Cw = 0.0 
    
            for j = 1:nele
                i = 1
                while ((ends[i,1] in w[:,1]) & (ends[i,2] in w[:,1]))|(!(ends[i,1] in w[:,1])&(ends[i,2] ∉ w[:,1]))
                    i = i+1
                end
                sn = ends[i,1]
                fn = ends[i,2]
    
                sn = convert(Int, sn)
                fn = convert(Int, fn)
    
                p = ((coord[sn,1]-xc)*(coord[fn,2]-yc)-(coord[fn,1]-xc)*(coord[sn,2]-yc))/L[i]
                if w[sn,1]==0
                    w[sn,1] = sn;
                    w[sn,2] = w[fn,2]-p*L[i]
                elseif w[fn,1]==0
                    w[fn,1] = fn;
                    w[fn,2] = w[sn,2]+p*L[i]
                end
                Iwx = Iwx+(1/3*(w[sn,2]*(coord[sn,1]-xc)+w[fn,2]*(coord[fn,1]-xc))+1/6*(w[sn,2]*(coord[fn,1]-xc)+w[fn,2]*(coord[sn,1]-xc)))*t[i]* L[i];
                Iwy = Iwy+(1/3*(w[sn,2]*(coord[sn,2]-yc)+w[fn,2]*(coord[fn,2]-yc))+1/6*(w[sn,2]*(coord[fn,2]-yc)+w[fn,2]*(coord[sn,2]-yc)))*t[i]* L[i];
            end
    
            if (Ix*Iy-Ixy^2)!=0.0
                xs = (Iy*Iwy-Ixy*Iwx)/(Ix*Iy-Ixy^2)+xc
                ys = -(Ix*Iwx-Ixy*Iwy)/(Ix*Iy-Ixy^2)+yc
            else
                xs = xc; ys = yc
            end
    
            if abs(xs/sqrt(A)) .< 1e-12
                xs = 0
            end
            if abs(ys/sqrt(A)) .< 1e-12
                ys = 0
            end
    
            # compute the unit warping
            for j = 1:nele
                i = 1
                while ((ends[i,1] in wo[:,1]) & (ends[i,2] in wo[:,1]))|(!(ends[i,1] in wo[:,1])&(ends[i,2] ∉ wo[:,1]))
                    i = i+1
                end
                sn = ends[i,1]
                fn = ends[i,2]
    
                sn = convert(Int, sn)
                fn = convert(Int, fn)
    
                po = ((coord[sn,1]-xs)*(coord[fn,2]-ys)-(coord[fn,1]-xs)*(coord[sn,2]-ys))/L[i]
                if wo[sn,1]==0
                    wo[sn,1] = sn;
                    wo[sn,2] = wo[fn,2]-po*L[i]
                elseif wo[convert(Int, ends[i,2]),1]==0
                    wo[fn,1] = fn;
                    wo[fn,2] = wo[sn,2]+po*L[i]
                end
                wno = wno+1/(2*A)*(wo[sn,2]+wo[fn,2])*t[i]* L[i];
            end
            wn = wno .-wo[:,2]
    
            # compute the warping constant
            for i = 1:nele
                sn = ends[i,1]; fn = ends[i,2]
    
                sn = convert(Int, sn)
                fn = convert(Int, fn)
    
                Cw = Cw+1/3*(wn[sn]^2+wn[sn]*wn[fn]+wn[fn]^2)*t[i]* L[i]
            end
    
            # transfer the shear center coordinates to the centroid principal coordinates
            s12 = [cos(theta) sin(theta); -sin(theta) cos(theta)]* transpose([xs-xc ys-yc]);
            # compute the polar radius of gyration of cross section about shear center
            ro = sqrt((I1+I2)/A+s12[1]^2+s12[2]^2)
    
            # compute B1 & B2
            B1 = 0 ; B2 = B1
            for i = 1:nele
                sn = ends[i,1]
                fn = ends[i,2];
    
                sn = convert(Int, sn)
                fn = convert(Int, fn)
    
                x1 = coord12[sn,1]; y1 = coord12[sn,2]
                x2 = coord12[fn,1]; y2 = coord12[fn,2]
                B1 = B1+((y1+y2)*(y1^2+y2^2)/4+(y1*(2*x1^2+(x1+x2)^2)+y2*(2*x2^2+(x1+x2)^2))/12)*L[i]*t[i]
                B2 = B2+((x1+x2)*(x1^2+x2^2)/4+(x1*(2*y1^2+(y1+y2)^2)+x2*(2*y2^2+(y1+y2)^2))/12)*L[i]*t[i]
            end
            B1 = B1/I1-2*s12[2]
            B2 = B2/I2-2*s12[1];
    
            if abs(B1/sqrt(A)) .< 1e-12
                B1 = 0
            end
            if abs(B2/sqrt(A)) .< 1e-12
                B2 = 0
            end
        elseif section == "arbitrary"
            J = sum(L.*t.^3)/3
            xs = NaN; ys = NaN; Cw = NaN; B1 = NaN; B2 = NaN; Pe = NaN; dcoord = NaN; wn=NaN
    
            #use the open section algorithm; modified to handle multiple parts; but
            #not completely generalized *that is a work for a future day* (17 Dec
            #2015) the primary goal here is to not crash the section property
            #calculators & applied stress generators when users have done a built
            #up section...
    
            # compute the torsional constant for open-section
            J = sum(L.*t.^3)/3
    
            # compute the shear center & initialize variables
            nnode = size(coord,1)
            w = zeros(nnode,2); w[Int(ends[1,1]),1] = ends[1,1]
            wo = zeros(nnode,2); wo[Int(ends[1,1]),1] = ends[1,1]
            Iwx = 0; Iwy = 0; wno = 0; Cw = 0
    
            for j = 1:nele
                i = 1
                while ((ends[i,1] in w[:,1]) & (ends[i,2] in w[:,1]))|(!(ends[i,1] in w[:,1])&(ends[i,2] ∉ w[:,1]))
                        i = i+1
                        if i>nele #brute force catch to continue calculation for multi part
                            i=nele
                            break
                        end
                end
                sn = ends[i,1]
                fn = ends[i,2]
    
                sn = Int(sn)
                fn = Int(fn)
    
                p = ((coord[sn,1]-xc)*(coord[fn,2]-yc)-(coord[fn,1]-xc)*(coord[sn,2]-yc))/L[i]
                if w[sn,1]==0
                    w[sn,1] = sn;
                    w[sn,2] = w[fn,2]-p*L[i]
                elseif w[fn,1]==0
                    w[fn,1] = fn;
                    w[fn,2] = w[sn,2]+p*L[i]
                end
                Iwx = Iwx+(1/3*(w[sn,2]*(coord[sn,1]-xc)+w[fn,2]*(coord[fn,1]-xc))+1/6*(w[sn,2]*(coord[fn,1]-xc)+w[fn,2]*(coord[sn,1]-xc)))*t[i]* L[i];
                Iwy = Iwy+(1/3*(w[sn,2]*(coord[sn,2]-yc)+w[fn,2]*(coord[fn,2]-yc))+1/6*(w[sn,2]*(coord[fn,2]-yc)+w[fn,2]*(coord[sn,2]-yc)))*t[i]* L[i];
            end
    
            if (Ix*Iy-Ixy^2) != 0
                xs = (Iy*Iwy-Ixy*Iwx)/(Ix*Iy-Ixy^2)+xc
                ys = -(Ix*Iwx-Ixy*Iwy)/(Ix*Iy-Ixy^2)+yc
            else
                xs = xc; ys = yc
            end
    
            if abs(xs/sqrt(A)) .< 1e-12
                xs = 0
            end
            if abs(ys/sqrt(A)) .< 1e-12
                ys = 0
            end
    
            # compute the unit warping
            for j = 1:nele
                i = 1
                while ((ends[i,1] in wo[:,1]) & (ends[i,2] in wo[:,1]))|(!(ends[i,1] in wo[:,1])&(ends[i,2] ∉ wo[:,1]))
                        i = i+1
                        if i>nele #brute force catch to continue calculation for multi part
                            i=nele
                            break
                        end
                end
                sn = ends[i,1]; fn = ends[i,2]
    
                sn = Int(sn)
                fn = Int(fn)
    
                po = ((coord[sn,1]-xs)*(coord[fn,2]-ys)-(coord[fn,1]-xs)*(coord[sn,2]-ys))/L[i]
                if wo[sn,1]==0
                    wo[sn,1] = sn;
                    wo[sn,2] = wo[fn,2]-po*L[i]
                elseif wo[Int(ends[i,2]),1]==0
                    wo[fn,1] = fn;
                    wo[fn,2] = wo[sn,2]+po*L[i]
                end
                wno = wno+1/(2*A)*(wo[sn,2]+wo[fn,2])*t[i]* L[i];
            end
            wn = wno .- wo[:,2]
    
            # compute the warping constant
            for i = 1:nele
                sn = ends[i,1]; fn = ends[i,2]
    
                sn = Int(sn)
                fn = Int(fn)
    
                Cw = Cw+1/3*(wn[sn]^2+wn[sn]*wn[fn]+wn[fn]^2)*t[i]* L[i]
            end
    
            # transfer the shear center coordinates to the centroid principal coordinates
            s12 = [cos(theta) sin(theta); -sin(theta) cos(theta)]*transpose([xs-xc ys-yc]);
            # compute the polar radius of gyration of cross section about shear center
            ro = sqrt((I1+I2)/A+s12[1]^2+s12[2]^2)
    
            # compute B1 & B2
            B1 = 0; B2 = B1
            for i = 1:nele
                sn = ends[i,1]; fn = ends[i,2];
    
                sn = Int(sn)
                fn = Int(fn)
    
                x1 = coord12[sn,1]; y1 = coord12[sn,2]
                x2 = coord12[fn,1]; y2 = coord12[fn,2]
                B1 = B1+((y1+y2)*(y1^2+y2^2)/4+(y1*(2*x1^2+(x1+x2)^2)+y2*(2*x2^2+(x1+x2)^2))/12)*L[i]*t[i]
                B2 = B2+((x1+x2)*(x1^2+x2^2)/4+(x1*(2*y1^2+(y1+y2)^2)+x2*(2*y2^2+(y1+y2)^2))/12)*L[i]*t[i]
            end
            B1 = B1/I1-2*s12[2]
            B2 = B2/I2-2*s12[1];
    
            if abs(B1/sqrt(A)) .< 1e-12
                B1 = 0
            end
            if abs(B2/sqrt(A)) .< 1e-12
                B2 = 0
            end
    
        end
    
        section_properties = SectionPropertiesObject(coord,ends,A,xc,yc,Ix,Iy,Ixy,theta,I1,I2,J,xs,ys,Cw,B1,B2,wn)
    
        return section_properties
    
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
    section_properties = cutwp_prop2(coord,ends)

    return section_properties

end

function open_thin_walled(X_c, Y_c, t)

    # X_c = [center[i][1] for i in eachindex(center)]
    # Y_c = [center[i][2] for i in eachindex(center)]

    #define elements for section property calc
    num_nodes = length(X_c)
    num_elem = num_nodes - 1
    coord = [X_c Y_c]
    ends = [1:num_nodes-1 2:num_nodes t]

    #calculate section properties
    section_properties = cutwp_prop2(coord,ends)

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

    section_properties = cutwp_prop2(coord,ends)

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



