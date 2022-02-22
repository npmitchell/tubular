function [faces, vertices, detailedOutput] = ...
    meshPlaneIntersect(faces, vertices, xplane, plane_normal, preview)
%meshPlaneIntersect(faces, vertices, xplane, plane_normal)
% Determine which faces intersect the plane. For the moment, this uses only
% x=const planes. ToDo: generalize to any plane by rotating mesh so that
% plane normal is along x=const, then rotating back.
%         |            |
% -.---.--|---.----.---|--.-----.
%         |            |
%         |            |
%        intersectf   
%
% Parameters
% ----------
% faces : F x 3 int array
%   face connectivity array for input mesh. faces(i, :) gives the indices
%   into vertices of face i
% vertices : N x 3 float array
%   vertices of input mesh. vertices(i, :) gives the 3d coordinates of mesh
%   vertex i
% xplane : float
%   position in plane which together with plane_normal specifies where to
%   bisect the mesh 
% plane_normal : optional, length 3 numeric array
%   normal to the plane bisecting the mesh
%
% Returns
% -------
% faces : (F-I+Q) x 3  int array
%   new face list, without intersected faces but with new faces that
%   subsection the intersected faces
% vertices : N+P x 3 float array
%   new vertices, with additional vertices for subsection faces spanning
%   each intersected face
% detailedOutput : struct with fields
%   intersectf : I x 1 int array
%       the indices into input faces of faces which are intersected and thus
%       removed from the mesh.
%     
% 
% NPMitchell 2020

% ToDo: generalize to any plane by rotating mesh so that
% plane normal is along x=const, then rotating back.
if nargin > 3
    error('todo: Add functionality to support arbitrary intersecting planes')
else
    preview = false ;
end

addpts = [] ;       % will append new points for bisected faces
addtri = [] ;       % will append new triangles for subsectioned faces

% Which side of the plane is each vertex on?
triRight = reshape(vertices(faces, 1) > xplane, size(faces)) ;
% xright --> the faces are right of the plane=3, left=0
xright = sum(triRight, 2) ;
% Intersecting faces are those whose xright row is not 000 or 111
intersectf = find(xright < 3 & xright > 0) ;

% For each face that intersects, make new triangles
% mesh_surface.faces = mesh.f ;
% mesh_surface.vertices = mesh.v ;
% plane1.faces = [1 2 3; 3, 4, 1] ;
% plane1.vertices = [xlower, -ww, -ww; ...
%                    xlower, -ww, ww; ...
%                    xlower, ww, ww; ...
%                    xlower, ww, -ww] ;
% [intMatrix, intSurface] = SurfaceIntersection(mesh_surface, plane1)

% new line segments subsectioning faces
lnsegs = zeros(length(intersectf), 2) ; 

nvtx = size(vertices, 1) ;
for pp = 1:length(intersectf) 
    fidR = intersectf(pp) ;
    vIds = faces(fidR, :)' ;
    vtcs = vertices(vIds, :) ;
    vtx = vtcs(:, 1) ;

    % Sort so that the left point(s) is(are) first
    [vtx, sortId] = sort(vtx) ;
    vIds = vIds(sortId) ;
    vtcs = vtcs(sortId, :) ;
    ds = abs(vtx - xplane) ;
    side = (vtx - xplane > 0) ;
    
    if all(sortId' == [1 2 3]) || all(sortId' == [2 3 1]) || all(sortId' == [3 1 2])
        orient = false ;
    else
        orient = true ;
    end

    % For each segment that traverses the plane, get new point on
    % the plane
    a1 = nvtx + size(addpts, 1) + 1 ;   % index for added point #1
    a2 = nvtx + size(addpts, 1) + 2 ;   % index for added point #2

    % There are two cases for side: [001] or [011], since the first
    % point is on the left and the last point is on the right.
    if all(side == [0; 0; 1])
        % seg 13
        v13 = vtcs(1, :) - vtcs(3, :) ;
        addpt1 = vtcs(3, :) + ds(3) / (ds(1) + ds(3)) * v13 ;
        % seg 23
        v23 = vtcs(2, :) - vtcs(3, :) ;
        addpt2 = vtcs(3, :) + ds(3) / (ds(2) + ds(3)) * v23 ; 

        % Check if addpt1 and addpt2 already are added points
        if ~isempty(addpts) && any(vecnorm(addpt1 - addpts, 2, 2) < eps)
            a1 = nvtx + find(vecnorm(addpt1 - addpts, 2, 2) < eps) ;
            a2 = a2 - 1;
        else
            % This is a new point. Add it to addpts
            addpts = cat(1, addpts, addpt1) ; 
        end
        
        % Check if addpt2 already has been added
        if ~isempty(addpts) && any(vecnorm(addpt2 - addpts, 2, 2) < eps) 
            a2 = nvtx + find(vecnorm(addpt2 - addpts, 2, 2) < eps) ;
        else
            % This is a new point. Add it to addpts
            addpts = cat(1, addpts, addpt2) ; 
        end
        % Add triangles -- note ordering matters! chirality preserved
        %             2  . 
        %               /|\
        %              / | \ 
        %             /  \  * a2
        %            /    | |\
        %           /     | | \
        %          /      | |  \
        %         /       | |   \
        %        /        | |    \
        %       /         | |     \
        %      /           \|      \
        % 1   .-------------*-------. 3
        %                   a1
        %
        if orient
            tri1 = [vIds(2), vIds(1), a1] ;
            tri2 = [vIds(2), a1, a2] ;
            tri3 = [a2, a1, vIds(3)] ;
        else
            tri1 = [vIds(1), vIds(2), a1] ;
            tri2 = [a1, vIds(2), a2] ;
            tri3 = [a1, a2, vIds(3)] ;
        end
        % Which side are these triangles on
        % tside==1 is on right, previous mesh segment=0
        % tside = [tside 1 0 0 ]; 
        addtri = [addtri; tri1; tri2; tri3] ;

        % line segment vertex indices
        if orient
            lnsegs(pp, :) = [a2, a1] ;
        else
            lnsegs(pp, :) = [a1, a2] ;
        end
        
        % Check it
        if preview
            vtxCheck = cat(1, vertices, addpts) ;
            newtri = [tri1; tri2; tri3] ;
            trisurf(newtri, vtxCheck(:, 1), vtxCheck(:, 2), ...
                vtxCheck(:, 3)) ;
            hold on;
            plot3(vertices(vIds(1), 1), vertices(vIds(1), 2), ...
                vertices(vIds(1), 3), 'ro')
            plot3(vertices(vIds(2), 1), vertices(vIds(2), 2), ...
                vertices(vIds(2), 3), 'go')
            plot3(vertices(vIds(3), 1), vertices(vIds(3), 2), ...
                vertices(vIds(3), 3), 'bo')
            plot3(addpt1(1), addpt1(2), addpt1(3), 'ko')
            plot3(addpt2(1), addpt2(2), addpt2(3), 'ko')
            pause(0.01)
        end

    elseif all(side == [0; 1; 1])
        % seg 21
        v21 = vtcs(2, :) - vtcs(1, :) ;
        addpt1 = vtcs(1, :) + ds(1) / (ds(2) + ds(1)) * v21 ; 
        % seg 13
        v13 = vtcs(1, :) - vtcs(3, :) ;
        addpt2 = vtcs(3, :) + ds(3) / (ds(1) + ds(3)) * v13 ;

        % Check if addpt1 and addpt2 already are added points
        if ~isempty(addpts) && any(vecnorm(addpt1 - addpts, 2, 2) < eps)
            a1 = nvtx + find(vecnorm(addpt1 - addpts, 2, 2) < eps) ;
            a2 = a2 - 1;
        else
            % This is a new point. Add it to addpts
            addpts = cat(1, addpts, addpt1) ; 
        end
        if ~isempty(addpts) && any(vecnorm(addpt2 - addpts, 2, 2) < eps) 
            a2 = nvtx + find(vecnorm(addpt2 - addpts, 2, 2) < eps) ;
        else
            % This is a new point. Add it to addpts
            addpts = cat(1, addpts, addpt2) ; 
        end

        % Add triangles -- note ordering matters! chirality preserved
        %         2 .  
        %         /  |
        %       /    |
        %   a1 .     | 
        %    / |\    |
        % 1 .  | |   |
        %    \ | \   |
        %      *  |  |
        %   a2  \ \  |
        %         \\ |
        %           *  3
        %
        if orient
            tri1 = [a1, vIds(1), a2] ;
            tri2 = [a1, a2, vIds(3)] ;
            tri3 = [vIds(2), a1, vIds(3)] ;
        else
            tri1 = [vIds(1), a1, a2] ;
            tri2 = [a2, a1, vIds(3)] ;
            tri3 = [a1, vIds(2), vIds(3)] ;
        end
        addtri = [addtri; tri1; tri2; tri3] ;

        % line segment vertex indices
        if orient
            lnsegs(pp, :) = [a1, a2] ;
        else
            lnsegs(pp, :) = [a2, a1] ;
        end
        
        % Check it
        if preview
            vtxCheck = cat(1, vertices, addpts) ;

            % Inspect this face alone
            newtri = [tri1; tri2; tri3] ;
            trisurf(newtri, vtxCheck(:, 1), vtxCheck(:, 2), ...
                vtxCheck(:, 3)) ;
            hold on;
            plot3(vertices(vIds(1), 1), vertices(vIds(1), 2), ...
                vertices(vIds(1), 3), 'ro')
            plot3(vertices(vIds(2), 1), vertices(vIds(2), 2), ...
                vertices(vIds(2), 3), 'go')
            plot3(vertices(vIds(3), 1), vertices(vIds(3), 2), ...
                vertices(vIds(3), 3), 'bo')
            plot3(addpt1(1), addpt1(2), addpt1(3), 'ko')
            plot3(addpt2(1), addpt2(2), addpt2(3), 'ko')
            pause(0.01)
        end
        
    end
    clearvars tri1 tri2 tri3 a1 a2
end

% Add the new vertices and faces
ntri = size(faces, 1) ;
vertices = cat(1, vertices, addpts) ;
faces = cat(1, faces, addtri) ;

% Prune triangles that were intersected. Note all vertices are kept.
faces(intersectf, :) = [] ;

if nargout > 2
    detailedOutput = struct() ;
    detailedOutput.added_vertices = addpts ;
    detailedOutput.added_faces = addtri ;
    detailedOutput.intersected_faces = intersectf ;
    detailedOutput.line_segments = lnsegs ; 
    detailedOutput.daughters = ntri:(ntri + size(addtri, 1)) ;
    mothers = repmat(intersectf', 3, 1) ;
    detailedOutput.mothers = mothers(:) ;
end
