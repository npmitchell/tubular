function net = networkFromSegmentation(skel, threefold, very_far)
% todo: finish this code to create network from segmentation
% Not needed currently, where we just take the raw shapes of cells as
% polygons
%
% NPMitchell, much adapted from Nick Noll's tissueAnalysisSuite

    net = struct() ;
    branchPoints = zeros(size(skel));
    rows = 2:(size(skel,1)-1);
    cols = 2:(size(skel,2)-1);
    branchPoints(rows,cols) = skel(rows+1,cols) + skel(rows-1,cols) + ...
                              skel(rows,cols+1) + skel(rows,cols-1);
    branchPoints = branchPoints .* skel;
    branchPoints = branchPoints >= 3;
    
    % Find vertices
    cc = bwconncomp(~skel, 4) ;
    LL = labelmatrix(cc) ;
    
    Vdat = seg.find_vertices(LL,threefold, very_far, 4);
    Cdat = seg.find_cells(Vdat,LL);
    if (mode==1)
        [Vdat]=seg.find_bonds(Vdat,L);
    end 
    net.Vdat = Vdat;
    net.Cdat = Cdat;

    for ii=1:length(Struct)
        for jj=1:length(Struct(ii).Cdat) 
            Struct(ii).Cdat(jj).ncells=setdiff(unique([Struct(ii).Vdat(Struct(ii).Cdat(jj).nverts).ncells]),jj);
        end
    end
 for v = 1:length(Struct(ii).Vdat)
            Struct(ii).Vdat(v).vertxcoord = double(Struct(ii).Vdat(v).vertxcoord);
            Struct(ii).Vdat(v).vertycoord = double(Struct(ii).Vdat(v).vertycoord);
        end
    end 