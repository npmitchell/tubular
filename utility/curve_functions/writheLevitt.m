function [Wr, wr] = writheLevitt(xyz, closed)
%WRITHELEVITT(xyz, s) Compute writhe by computing solid angle of segments
%   Wr = sum_i sum_j Omega_ij / 4pi = 2 sum_i=2..N sum_j<i Omega_ij / 4pi
%   Note that if the curve is closed, the first and last point in xyz need
%   not be identical. Note that the writhe density is zero for
%   adjacent segments, since the solid angle of two segs sharing a pt is 
%   zero (there is only a 1d curve of the unit sphere in which view they
%   cross, with a solid angle of zero). 
%   Numerically, these would jump around from -pi, -pi/2, 0, pi/2, and pi.
%   Therefore the writhe density ignores adjacent segments.
%   
% Parameters
% ----------
% xyz : N x 3 float array
%   The 3D curve whose writhe we compute
% ss : N x 1 float array, optional
%   The cumulative pathlength from the endpoint of the curve xyz(1, :) to a
%   given curve point
% 
% Returns
% -------
% Wr : float 
%   The writhe of the curve as computed using the Gauss Integral
% wr : N x 1 float array 
%   The writhe "density" obtained by integrating only over ds', not over ds
% 
% NPMitchell 2019
lsegs = linesegments(xyz, closed) ;

% Consider each linesegment and compare to all others
if nargout == 1
    Omega = zeros((size(lsegs, 1) - 1) * (size(lsegs, 1) * 0.5), 1) ;
    % wr = zeros(size(lsegs, 1), size(lsegs, 1)) ;
    dmyk = 1 ;
    for ii = 2:size(lsegs, 1)
        s1 = lsegs(ii, 1:3) ;
        s2 = lsegs(ii, 4:6) ;

        for jj = 1:(ii - 1)
            if ii == jj
                error('Should not enter this situation -> double-counting.')
            elseif abs(ii - jj) > 1
                % Note that the writhe density is very noisy and ill 
                % defined for adjacent segments so we skip those.
                s3 = lsegs(jj, 1:3) ;
                s4 = lsegs(jj, 4:6) ;
                r12 = s2 - s1 ;
                r13 = s3 - s1 ;
                r14 = s4 - s1 ;
                r23 = s3 - s2 ;
                r24 = s4 - s2 ;
                r34 = s4 - s3 ;

                % Define normalized cross products
                n1 = cross(r13, r14) ;
                if vecnorm(n1) > 0
                    n1 = n1 / vecnorm(n1) ;
                end
                n2 = cross(r14, r24) ;
                if vecnorm(n2) > 0
                    n2 = n2 / vecnorm(n2) ;
                end
                n3 = cross(r24, r23) ;
                if vecnorm(n3) > 0
                    n3 = n3 / vecnorm(n3) ;
                end
                n4 = cross(r23, r13) ;
                if vecnorm(n4) > 0
                    n4 = n4 / vecnorm(n4) ;
                end
                
                % Sum 
                ostar1 = asin(dot(n1, n2)) ;
                ostar2 = asin(dot(n2, n3)) ;
                ostar3 = asin(dot(n3, n4)) ;
                ostar4 = asin(dot(n1, n4)) ;
                
                % Put into range --> should already be in range!
                % ostars = [ostar1, ostar2, ostar3, ostar4] ;
                % for qq = 1:4
                %     oqq = ostars(qq) ;
                %     if oqq > 0
                %         if oqq > pi*0.5
                %             error('break')
                %              plot3([s1(:, 1), s2(:, 1)], ...
                %                    [s1(:, 2), s2(:, 2)], ...
                %                    [s1(:, 3), s2(:, 3)])
                %              hold on;
                % 
                %              plot3([s3(:, 1), s4(:, 1)], ...
                %                    [s3(:, 2), s4(:, 2)], ...
                %                    [s3(:, 3), s4(:, 3)])
                %              hold on;
                %              plot3(s1(:, 1), s1(:, 2), s1(:, 3), 'ro')
                %              plot3(s2(:, 1), s2(:, 2), s2(:, 3), 'gs')
                %              plot3(s3(:, 1), s3(:, 2), s3(:, 3), 'b^')
                %              plot3(s4(:, 1), s4(:, 2), s4(:, 3), 'k.')
                %         end
                %     else
                %         if oqq < -pi * 0.5
                %             error('break')
                %         end
                %     end
                % end
                % modulo_into_range(ostar1, -pi*0.5, pi*0.5) ;

                ostar = ostar1 + ostar2 + ostar3 + ostar4 ;
                r34xr12dr13 = dot(cross(r34, r12), r13) ;
                Omega(dmyk) = ostar * sign(r34xr12dr13) ;
                % wr(ii, jj) = Omega(dmyk) ;
                dmyk = dmyk + 1 ;
            end
        end
    end
    Omega = Omega(1:dmyk-1) ;
    Wr = nansum(Omega) / (2 * pi) ; % Note that there is a factor of 2 / 4pi
else
    % Return the writhe density as well
    Omega = zeros(size(lsegs, 1), 1) ;
    wr = zeros(size(lsegs, 1), 1) ;
    for ii = 1:size(lsegs, 1)
        s1 = lsegs(ii, 1:3) ;
        s2 = lsegs(ii, 4:6) ;
        for jj = 1:size(lsegs, 1)
            % Skip double-counting
            % Note that the writhe density is very noisy and ill defined
            % for adjacent segments.
            if abs(ii - jj) > 1
                s3 = lsegs(jj, 1:3) ;
                s4 = lsegs(jj, 4:6) ;
                r12 = s2 - s1 ;
                r13 = s3 - s1 ;
                r14 = s4 - s1 ;
                r23 = s3 - s2 ;
                r24 = s4 - s2 ;
                r34 = s4 - s3 ;

                % Define normalized cross products
                n1 = cross(r13, r14) ;
                if vecnorm(n1) > 0
                    n1 = n1 / vecnorm(n1) ;
                end
                n2 = cross(r14, r24) ;
                if vecnorm(n2) > 0
                    n2 = n2 / vecnorm(n2) ;
                end
                n3 = cross(r24, r23) ;
                if vecnorm(n3) > 0
                    n3 = n3 / vecnorm(n3) ;
                end
                n4 = cross(r23, r13) ;
                if vecnorm(n4) > 0
                    n4 = n4 / vecnorm(n4) ;
                end

                % Sum 
                ostar1 = asin(dot(n1, n2)) ;
                ostar2 = asin(dot(n2, n3)) ;
                ostar3 = asin(dot(n3, n4)) ;
                ostar4 = asin(dot(n1, n4)) ;
                ostar = real(ostar1 + ostar2 + ostar3 + ostar4) ;
                r34xr12dr13 = dot(cross(r34, r12), r13) ;
                Omega(ii) = Omega(ii) + ostar * sign(r34xr12dr13) ;
                
            end
        end
    end
    Wr = nansum(Omega) / (4 * pi) ; % Note that we div 4pi for solid sphere
    
    % Get density by dividing by ds
    if closed
        ds = vecnorm(diff([xyz; xyz(1,:)]), 2, 2) ;
    else
        ds = vecnorm(diff(xyz), 2, 2) ;
    end
    wr = ( Omega ./ ds ) / (4 * pi) ;
end

end

