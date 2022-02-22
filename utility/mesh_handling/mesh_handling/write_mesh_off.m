function write_mesh_off( fullfilename, tr )
%WRITE_MESH_OFF Saves triangulation as an .off file
%   fullfilename: 
%
%   tr: MATLAB triangulation
% 

%--------------------------------------------------------------------------
% Input processing
%--------------------------------------------------------------------------
V = tr.Points;
if size(V,2) < 3
    V = [V zeros(size(V,1),1)]; % ONLY WRITES 3D MESHES
end

F = tr.ConnectivityList; % ONLY WRITES TRIANGULATIONS
E = tr.edges;


if(exist('fullfilename','var') == 0)
    [filename, filefolder] = uiputfile('*.off', 'Write OFF-file');
    fullfilename = [filefolder filename];
end
[filefolder, filename] = fileparts( fullfilename );

fid = fopen(fullfilename, 'wt');

%--------------------------------------------------------------------------
% Write preamble
%--------------------------------------------------------------------------
fprintf(fid,'OFF\n');
fprintf(fid,['# ',fullfilename,'\n']);
fprintf(fid, '\n');

fprintf(fid, '%d %d %d\n', size(V,1), size(F,1), size(E,1) );

%--------------------------------------------------------------------------
% Write vertices 
%--------------------------------------------------------------------------
formatSpec =  '%0.6f %0.6f %0.6f\n';
for i=1:size(V,1)
   fprintf( fid, formatSpec, V(i,1), V(i,2), V(i,3) ); 
end

%--------------------------------------------------------------------------
% Write face connectivity list - NOTE: .off FILES USE 0-INDEXING
%--------------------------------------------------------------------------
formatSpec = '%d %d %d %d\n';
for i=1:size(F,1)
    fprintf( fid, formatSpec, 3, ...
        F(i,1)-1, F(i,2)-1, F(i,3)-1 );
end

fclose(fid);


end

