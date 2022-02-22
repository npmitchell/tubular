%Commented by Noah Mitchell 2015
% Calculate the Free energy F and force on each particle dF 
% Note: NL2 is a version of NL with zeros replaced by ones (taken care of by
% KL. In other words, we create artificial coupling of all nodes except 
% node 1 with node 1.
%vk made a slight modification of NL2 (9 Oct 12)
% NPM added force input argument, which provides a force (N x 3 array) on
% each of the N particles (could be zero for particle not experiencing
% forcing).

function [F,dF] = wFree_Energy4_dF_forcing(xyz,NL,BL,KL,TRISm,BLm,KLm,force,Fconst)

%calc Free energy
temp = (xyz(TRISm(:,1),:)-xyz(TRISm(:,2),:));
F = 1/2*sum(KLm.*(sqrt(dot(temp,temp,2))-BLm).^2) ;

%For each node where there is a boundary force, subtracting the work done
%by the boundary forces from the energy.
matrix = force .* xyz ;
F = F - sum(matrix(:)) + Fconst; 

% vk added this part
NL2=NL;
zeros_first_row = NL2(1,:)==0; %vk added to find the column numbers of the first row 
NL2(1, zeros_first_row)=2; %vk added to create artificial coupling of node 1 with node 2 (and not node 1!)
NL2(NL2==0)=1; %vk create artificial coupling of all nodes except node 1 with node 1


%Each row of dF is the net force on each particle in x,y,z.
%Preallocate dF for speed.
dF = zeros(size(xyz,1),3);

%Summing over nearest neighbors, sum the force on each particle
for kk=1:size(NL2,2)
    %Cartesian vectors of NNeighbors
    vec=(xyz-xyz(NL2(:,kk),:)); 
    
    %dis is the scalar distance between each particle and its kk neighbor.
    dis = sqrt(dot(vec,vec,2));
    %dis(1) = 1; %vk is no longer needed
    
    %Debugging
    %fprintf(['size(dis)=',num2str(length(dis)),'\n']) 
    %fprintf(['size(BL(:,kk))=',num2str(length(BL(:,kk))),'\n']) 
    
    %The Lagrangian strain is: (dis-BL(:,kk))./dis  
    % dF is K*displacement*normalvector, where vec/dis is the normal vector
    dF = dF+((KL(:,kk).*(dis-BL(:,kk))./dis)*[1,1,1]).*vec;  
    
end

% dF is the gradient of F, so it points opposite to the force.
% Apply extra force (given by traction) on nd nodes.
dF = dF - force; 
