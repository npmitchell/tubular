% wFRnd5(tol,xyz,NL,KL,BL,TRISm,KLm,BLm,nd,projectfcn,varargin)
%
% wtmi 2011 - Chicago
% Fletcher-Reeves minimization of lattice energy
% Last changes 15 mar 2011 - Early code was based on Lucks' fr code, rewritten for 
% passing projection functions, protecting bonds, speedup and exit
% commented npm 2014-2015
% NPM added gravity 2018

% nmaxiter : int
%   maximum number of iterations before energy minimization exits

function [xyz,val,iter,serie] = wFRnd6_pressure(tol,xyz,NL,BL,KL,TRI,...
    TRISm,BLm,KLm,pressure,nd,nmaxiter,projectfcn,varargin)

    serie = [];

    %NL(find(NL==0))=1;
    % Find the normals to all triangles in the network. Define the pressure to
    % be pointing in the normal direction
    normals = TRI2normals(TRI, xyz, 1, false);
    % pressure
    % normals(1)
    if normals(1) == NaN
        exit('getting NaNs')
    end

    forces = pressure * normals ;
    % forces(1:3, :)
    Fconst = sum(dot(xyz, normals));


    [F,dF] = wFree_Energy4_dF_forcing(xyz,NL,BL,KL,TRISm,BLm,KLm,forces,Fconst);
    % H is the negative gradient of the free energy
    H=-dF; 
    H(isnan(H))=0;
    H(nd,:) = 0;

    for i=1:nmaxiter
        FOld = F;
        HOld = H;

        % % Debugging:
        % fprintf(['size(xyz)=' ,num2str(size(xyz)),'\n'])
        % fprintf(['size(-dF)=' ,num2str(size(dF )),'\n'])

        % optimize amount of displacement in dir along dF to min FreeEnergy
        % Minimize the Free energy by moving along the negative gradient in
        % Free energy (which is the force, called H here)
        eps = fminbnd(@wFmin4_forcing,0,1,[],xyz,H,TRISm,BLm,KLm,forces,Fconst,projectfcn,varargin{:});
        % Note that the above could return [eps, fval, exitflag]
        % eps
        % fval
        % exitflag

        xyz = xyz + eps*H;
        [x,y,z] = projectfcn(xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});
        xyz = [x,y,z];

        % Obtain new normal vectors, update pressure (forcing) and Fconst
        normals = TRI2normals(TRI, xyz, 1, false);
        forces = pressure * normals;
        Fconst = sum(dot(xyz, normals));

        iter=i;

        [F,dF] = wFree_Energy4_dF_forcing(xyz,NL,BL,KL,TRISm,BLm,KLm,forces,Fconst);
        dF(isnan(dF))=0;
        dF(nd,:) = 0;
        H=-dF;

        serie = cat(1,serie,F);
        val = F;
        % F-FOld;

        if (2*abs(F-FOld) <= tol*(abs(F)+abs(FOld)))
            iter=i;
            val=F;
            % return the xyz positions, the free energy, iteration number,
            % and series of free energies (should be descending)
            return;
        else
            h_new = [H(:,1); H(:,2); H(:,3)];
            h_old = [HOld(:,1); HOld(:,2); HOld(:,3)];
            gamma = dot(h_new,h_new)/dot(h_old,h_old);

            if(dot(h_new,h_new) <= 0.)
              % 'Zero Gradient Exit'
              iter=i;
              val=F;
              return;
            else
              H = H+gamma*HOld; 
              H(nd,:) = 0;

            end
        end

    end
return;
