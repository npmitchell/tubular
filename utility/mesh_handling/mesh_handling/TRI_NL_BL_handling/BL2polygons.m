function BL2polygons(xy, BL, NL, KL, PVx, PVy, PVxydict, viewmethod, check)
% Extract polygons from a lattice of points.
% Note that dangling bonds are removed, but no points are removed.
% This allows correct indexing for PVxydict keys, if supplied.
% This code fails if a site is its own NNN.
% adapted from extract_polygons() from NPMitchell's python lepm library.
% UNFINISHED: used Cdat2polygons() for now!
%
% Parameters
% ----------
% xy : NP x 2 float array
%     points living on vertices of dual to triangulation
% BL : Nbonds x 2 int array
%     Each row is a bond and contains indices of connected points
% NL : NP x NN int array (optional, speeds up calc if it is known there are no dangling bonds)
%     Neighbor list. The ith row has neighbors of the ith particle, padded with zeros
% KL : NP x NN int array (optional, speeds up calc if it is known there are no dangling bonds)
%     Connectivity list. The ith row has ones where ith particle is connected to NL[i,j]
% PVx : NP x NN float array (optional, for periodic lattices and speed)
%     ijth element of PVx is the x-component of the vector taking NL[i,j] to its image as seen by particle i
%     If PVx and PVy are specified, PVxydict need not be specified.
% PVy : NP x NN float array (optional, for periodic lattices and speed)
%     ijth element of PVy is the y-component of the vector taking NL[i,j] to its image as seen by particle i
%     If PVx and PVy are specified, PVxydict need not be specified.
% PVxydict : dict (optional, for periodic lattices)
%     dictionary of periodic bonds (keys) to periodic vectors (values)
%     If key = (i,j) and val = np.array([ 5.0,2.0]), then particle i sees particle j at xy[j]+val
%     --> transforms into:  ijth element of PVx is the x-component of the vector taking NL[i,j] to its image as seen
%     by particle i
% viewmethod : bool
%     View the results of many intermediate steps
% check: bool
%     Check the initial and final result
% eps : float
%     minimum value to consider nonzero in KL
%
% Returns
% ----------
% polygons : list of lists of ints
%     list of lists of indices of each polygon
% 
viewmethod = true
if size(xy, 1) < size(xy, 2) && size(xy, 2) == 2
    disp('Assuming that xy is given as transposed coord list. Reshaping...')
    xy = xy' ;
end
NP = size(xy, 1) ;

if isempty(KL) || isempty(NL)
    NL, KL = BL2NLandKL(BL, NP, [])
    if any(BL < 0)
        if len(PVxydict) > 0
            PVx, PVy = PVxydict2PVxPVy(PVxydict, NL, KL)
        else
            error('Must specify either PVxydict or KL and NL in extract_polygons_lattice()' +
                               ' when periodic bonds exist!')
        end
    end
elseif any(BL < 0)
    if isempty(PVx) || isempty(PVy)
        if isempty(PVxydict)
            error('Must specify either PVxydict or PVx and PVy in extract_polygons_lattice()' +
                               ' when periodic bonds exist!')
        else
            PVx, PVy = PVxydict2PVxPVy(PVxydict, NL, KL)
        end
    end 
end
NN = np.size(KL, 1)
% Remove dangling bonds
% dangling bonds have one particle with only one neighbor
finished_dangles = false
while ~finished_dangles
    dangles = np.where([np.count_nonzero(row) == 1 for row in KL])[0]
    if len(dangles) > 0
        % Check if need to build PVxy dictionary from PVx and PVy before changing NL and KL
        if any(BL < 0) && isempty(PVxydict)
            PVxydict = PVxy2PVxydict(PVx, PVy, NL, KL=KL)
        end
        
        % Make sorted bond list of dangling bonds
        dpair = np.sort(np.array([[d0, NL[d0, np.where(KL[d0] != 0)[0]]] for d0 in dangles]), axis=1)
        % Remove those bonds from BL
        BL = dh.setdiff2d(BL, dpair.astype(BL.dtype))
        % print 'dpair = ', dpair
        % print 'ending BL = ', BL
        NL, KL = BL2NLandKL(BL, NP=NP, NN=NN)

        % Now that NL and KL rebuilt (changed), (re)build PVx and PVy if periodic bcs
        if (BL < 0).any()
            if len(PVxydict) > 0
                PVx, PVy = PVxydict2PVxPVy(PVxydict, NL, KL)
    else
        finished_dangles = true

if viewmethod or check:
    print 'Plotting result after chopped dangles, if applicable...'
    display_lattice_2D(xy, BL, NL=NL, KL=KL, PVx=PVx, PVy=PVy, PVxydict=PVxydict,
                       title='Result after chopping dangling bonds', close=false)
    for i in range(len(xy))
        plt.text(xy[i, 0] + 0.2, xy[i, 1], str(i))
    plt.show()

% bond markers for counterclockwise, clockwise
used = np.zeros((len(BL), 2), dtype=bool)
polygons = []
finished = false
if viewmethod
    f, (ax1, ax2) = plt.subplots(1, 2)

% For periodicity, remember which bonds span periodic boundary
periB = np.array([(row < 0).any() for row in BL])

if periB.any() and PVxydict is None and (PVx is None or PVy is None)
    raise RuntimeError('Periodic boundaries have been detected, but no periodic vectors supplied to ' +
                       'extract_polygons_lattice()')

if ~periB.any()
    print 'no PBCs, calculating polygons...'
    while ~finished
        % Check if all bond markers are used in order A-->B
        % print 'Checking AB (A-->B): '
        todoAB = np.where(~used[:, 0])[0]
        % print 'len(todoAB) = ', len(todoAB)
        % print 'used = ', used
        % print 'todoAB = ', todoAB
        % print polygons
        if len(todoAB) > 0
            bond = BL[todoAB[0]]
            % if (bond == [21, 22]).all()
            %     for todoab in todoAB:
            %         ax1.plot([xy[BL[todoab, 0], 0], xy[BL[todoab, 1], 0]],
            %                  [xy[BL[todoab, 0], 1], xy[BL[todoab, 1], 1]], 'b-', lw=3)
            %     todoBA = np.where(~used[:, 1])[0]
            %     for todoba in todoBA:
            %         ax1.plot([xy[BL[todoba, 0], 0], xy[BL[todoba, 1], 0]],
            %                  [xy[BL[todoba, 0], 1], xy[BL[todoba, 1], 1]], 'g--')
            %     print 'bond = ', bond
            %     plt.pause(40)
            %     sys.exit()

            % bb will be list of polygon indices
            % Start with orientation going from bond[0] to bond[1]
            nxt = bond[1]
            bb = [bond[0], nxt]
            dmyi = 1

            % Now mark the new bond that has now been added to bb as used
            % Get index of used matching thisbond
            mark_used = np.where((np.logical_or(BL == bb[0], BL == bb[1])).all(axis=1))
            % print 'marking bond [', thisbond, '] as used'
            used[mark_used, 0] = true

            %%%%%%%%%%%%%%%
            % check
            if viewmethod
                ax1.plot(xy[:, 0], xy[:, 1], 'k.')
                ax1.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                             xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                             arrowprops=dict(arrowstyle="->",
                                             color="r",
                                             shrinkA=5, shrinkB=5,
                                             patchA=None,
                                             patchB=None,
                                             connectionstyle="arc3,rad=0.2", ), )
                for i in range(len(xy))
                    ax1.text(xy[i, 0] + 0.2, xy[i, 1], str(i))
                end
                ax2.imshow(used)
                ax1.set_aspect('equal')
            end
            %%%%%%%%%%%%%%%

            % as long as we haven't completed the full outer polygon, add next index
            while nxt ~= bond[0]
                n_tmp = NL[nxt, np.argwhere(KL[nxt]).ravel()]
                % Exclude previous boundary particle from the neighbors array, unless its the only one
                % (It cannot be the only one, if we removed dangling bonds)
                if len(n_tmp) == 1:
                    '''The bond is a lone bond, not part of a triangle.'''
                    neighbors = n_tmp
                else
                    neighbors = np.delete(n_tmp, np.where(n_tmp == bb[dmyi - 1])[0])
                end
                
                angles = np.mod(np.arctan2(xy[neighbors, 1] - xy[nxt, 1], xy[neighbors, 0] - xy[nxt, 0]).ravel() \
                                - np.arctan2(xy[bb[dmyi - 1], 1] - xy[nxt, 1],
                                             xy[bb[dmyi - 1], 0] - xy[nxt, 0]).ravel(), 2 * np.pi)
                nxt = neighbors[angles == max(angles)][0]
                bb.append(nxt)

                %%%%%%%%%%%%%%%
                % % Check
                % if viewmethod
                %     plt.annotate("", xy=(xy[bb[dmyi],0],xy[bb[dmyi],1] ), xycoords='data',
                %             xytext=(xy[nxt,0], xy[nxt,1]), textcoords='data',
                %             arrowprops=dict(arrowstyle="->",
                %                             color="r",
                %                             shrinkA=5, shrinkB=5,
                %                             patchA=None,
                %                             patchB=None,
                %                             connectionstyle="arc3,rad=0.2",),  )
                %
                %%%%%%%%%%%%%%%

                % Now mark the new bond that has now been extended (added) as used
                thisbond = [bb[dmyi], bb[dmyi + 1]]
                % Get index of used matching thisbond
                mark_used = np.where((np.logical_or(BL == bb[dmyi], BL == bb[dmyi + 1])).all(axis=1))

                % mark_used = np.where((BL == thisbond).all(axis=1))
                if ~used[mark_used, 0]:
                    % print 'marking bond [', thisbond, '] as used'
                    used[mark_used, 0] = true
                else
                    % Get index of used matching reversed thisbond (this list boolean is directional)
                    % mark_used = np.where((BL == thisbond[::-1]).all(axis=1))
                    % Used this bond in reverse order
                    used[mark_used, 1] = true
                end
                % print 'used = ', used
                dmyi += 1
            end

            polygons.append(bb)
            %%%%%%%%%%%%%%%
            % Check new polygon
            if viewmethod
                ax1.plot(xy[:, 0], xy[:, 1], 'k.')
                for i in range(len(xy))
                    ax1.text(xy[i, 0] + 0.2, xy[i, 1], str(i))
                end
                for dmyi in range(len(bb))
                    nxt = bb[np.mod(dmyi + 1, len(bb))]
                    ax1.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                                 xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                                 arrowprops=dict(arrowstyle="->",
                                                 color="r",
                                                 shrinkA=5, shrinkB=5,
                                                 patchA=None,
                                                 patchB=None,
                                                 connectionstyle="arc3,rad=0.2", ), )
                end
                ax2.cla()
                ax2.imshow(used)
                plt.pause(0.00001)
                %%%%%%%%%%%%%%%
            end
        else
            % Check for remaining bonds unused in reverse order (B-->A)
            % print 'CHECKING REVERSE (B-->A): '
            todoBA = np.where(~used[:, 1])[0]
            if len(todoBA) > 0
                bond = BL[todoBA[0]]

                %%%%%%%%%%%%%%%
                % % check
                % if viewmethod
                %     plt.annotate("", xy=(xy[bb[dmyi],0],xy[bb[dmyi],1] ), xycoords='data',
                %             xytext=(xy[nxt,0], xy[nxt,1]), textcoords='data',
                %             arrowprops=dict(arrowstyle="->",
                %                         color="b",
                %                         shrinkA=5, shrinkB=5,
                %                         patchA=None,
                %                         patchB=None,
                %                         connectionstyle="arc3,rad=0.6",),  )
                % %%%%%%%%%%%%%%%

                % bb will be list of polygon indices
                % Start with orientation going from bond[0] to bond[1]
                nxt = bond[0]
                bb = [bond[1], nxt]
                dmyi = 1

                % Now mark the new bond that has now been added to bb as used
                % Get index of used matching thisbond
                thisbond = [bb[dmyi], bb[dmyi - 1]]
                mark_used = np.where((BL == thisbond).all(axis=1))
                % print 'marking bond [', thisbond, '] as used'
                used[mark_used, 1] = true

                % as long as we haven't completed the full outer polygon, add nextIND
                while nxt != bond[1]:
                    n_tmp = NL[nxt, np.argwhere(KL[nxt]).ravel()]
                    % Exclude previous boundary particle from the neighbors array, unless its the only one
                    % (It cannot be the only one, if we removed dangling bonds)
                    if len(n_tmp) == 1:
                        '''The bond is a lone bond, not part of a triangle.'''
                        neighbors = n_tmp
                    else
                        neighbors = np.delete(n_tmp, np.where(n_tmp == bb[dmyi - 1])[0])
                    end
                    angles = np.mod(np.arctan2(xy[neighbors, 1] - xy[nxt, 1], xy[neighbors, 0] - xy[nxt, 0]).ravel() \
                                    - np.arctan2(xy[bb[dmyi - 1], 1] - xy[nxt, 1],
                                                 xy[bb[dmyi - 1], 0] - xy[nxt, 0]).ravel(), 2 * np.pi)
                    nxt = neighbors[angles == max(angles)][0]
                    bb.append(nxt)

                    %%%%%%%%%%%%%%%
                    % Check
                    % if viewmethod
                    %     plt.annotate("", xy=(xy[bb[dmyi],0],xy[bb[dmyi],1] ), xycoords='data',
                    %         xytext=(xy[nxt,0], xy[nxt,1]), textcoords='data',
                    %         arrowprops=dict(arrowstyle="->",
                    %                     color="b",
                    %                     shrinkA=5, shrinkB=5,
                    %                     patchA=None,
                    %                     patchB=None,
                    %                     connectionstyle="arc3,rad=0.6", %connectionstyle,
                    %                     ),  )
                    %%%%%%%%%%%%%%%

                    % Now mark the current bond as used --> note the inversion of the bond order to match BL
                    thisbond = [bb[dmyi + 1], bb[dmyi]]
                    % Get index of used matching [bb[dmyi-1],nxt]
                    mark_used = np.where((BL == thisbond).all(axis=1))
                    if len(mark_used) > 0
                        used[mark_used, 1] = true
                    else
                        raise RuntimeError('Cannot mark polygon bond as used: this bond was already used '
                                           'in its attempted orientation. (All bonds in first column '
                                           'should already be marked as used.)')
                    end
                    dmyi += 1

                polygons.append(bb)

                % Check new polygon
                if viewmethod
                    ax1.plot(xy[:, 0], xy[:, 1], 'k.')
                    for i in range(len(xy))
                        ax1.text(xy[i, 0] + 0.2, xy[i, 1], str(i))
                    end
                    for dmyi in range(len(bb))
                        nxt = bb[np.mod(dmyi + 1, len(bb))]
                        ax1.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                                     xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                                     arrowprops=dict(arrowstyle="->",
                                                     color="b",
                                                     shrinkA=5, shrinkB=5,
                                                     patchA=None,
                                                     patchB=None,
                                                     connectionstyle="arc3,rad=0.6", ), )
                    end
                    ax2.cla()
                    ax2.imshow(used)
                    plt.pause(0.00001)
                    %%%%%%%%%%%%%%%
                end
            else
                % All bonds have been accounted for
                finished = true
            end
else
    print 'detected periodicity...'
    % get particles on the finite (non-periodic) system's boundary. This allows massive speedup.
    KLfin = np.zeros_like(KL)
    KLfin[KL > 0] = 1
    % Create BLfin to pass to extract_boundary()
    prows = np.where(BL < 0)[0]
    nprows = np.setdiff1d(np.arange(len(BL)), prows)
    if check:
        print 'rows of BL that are periodic: ', prows
        print 'BL[prows] = ', BL[prows]
    BLfin = BL[nprows]
    finbd = extract_boundary(xy, NL, KLfin, BLfin, check=check)

    % If there were dangling points in the non-periodic representation, then we need to add those to finbd because
    % they will have periodic bonds attached to them.
    dangles = np.where(~KLfin.any(axis=1))[0]
    print 'dangles = ', dangles
    if len(dangles) > 0
        print 'Found dangling points in the finite/non-periodic representation. Adding to finbd...'
        finbd = np.hstack((finbd, np.array(dangles)))

    if check:
        print 'finite boundary: finbd = ', finbd
        plt.clf()
        display_lattice_2D(xy, BL, NL=NL, KL=KLfin, PVx=PVx, PVy=PVy, PVxydict=PVxydict,
                           title='Identified finite boundary', close=false)
        for i in range(len(xy))
            plt.text(xy[i, 0] + 0.2, xy[i, 1], str(i))
        plt.plot(xy[finbd, 0], xy[finbd, 1], 'ro')
        plt.show()
    first_check = true

    % Then erase periodicity in BL
    BL = np.abs(BL)

    while ~finished:
        if len(polygons) % 20 == 0:
            print 'constructed ', len(polygons), ' polygons...'
        % Check if all bond markers are used in order A-->B
        % print 'Checking AB (A-->B): '
        todoAB = np.where(~used[:, 0])[0]
        % print 'len(todoAB) = ', len(todoAB)
        % print 'used = ', used
        % print 'todoAB = ', todoAB
        if len(todoAB) > 0
            bond = BL[todoAB[0]]

            % bb will be list of polygon indices
            % Start with orientation going from bond[0] to bond[1]
            nxt = bond[1]
            bb = [bond[0], nxt]
            dmyi = 1

            % define 'previous angle' as backwards of current angle -- ie angle(prev-current_pos)
            % Must include effect of PV on this angle -- do in ref frame of nxt particle
            PVind = np.argwhere(NL[nxt] == bond[0])[0][0]
            addx = PVx[nxt, PVind]
            addy = PVy[nxt, PVind]
            xyb0 = xy[bond[0], :] + np.array([addx, addy])
            prev_angle = np.arctan2(xyb0[1] - xy[nxt, 1], xyb0[0] - xy[nxt, 0]).ravel()

            %%%%%%%%%%%%%%%
            % check
            if viewmethod
                if first_check:
                    ax1.plot(xy[:, 0], xy[:, 1], 'k.')
                    for i in range(len(xy))
                        ax1.text(xy[i, 0] + 0.2, xy[i, 1], str(i))
                    first_check = false

                ax1.annotate("", xy=(xy[bb[dmyi - 1], 0], xy[bb[dmyi - 1], 1]), xycoords='data',
                             xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                             arrowprops=dict(arrowstyle="->",
                                             color="r",
                                             shrinkA=5, shrinkB=5,
                                             patchA=None,
                                             patchB=None,
                                             connectionstyle="arc3,rad=0.2", ), )
                ax2.imshow(used, aspect=1. / len(used), interpolation='none')
                ax1.set_aspect('equal')
            %%%%%%%%%%%%%%%
            % define the displacment from the starting point that we have moved so far
            displ = xy[nxt] - xyb0

            % as long as we haven't completed the full outer polygon, add next index
            while nxt != bond[0] or abs(displ[0]**2 + displ[1]**2) > eps:
                % print nxt
                %            o     o neighbors
                %             \   /
                %              \ /
                %               o nxt
                %             /
                %           /
                %         o  bb[dmyi-1]
                %
                n_tmp = NL[nxt, np.argwhere(KL[nxt]).ravel()]
                % Exclude previous boundary particle from the neighbors array, unless its the only one
                % (It cannot be the only one, if we removed dangling bonds)
                if len(n_tmp) == 1:
                    '''The bond is a lone bond, not part of a triangle/polygon.'''
                    neighbors = n_tmp
                else
                    % Remove the current particle from the list of its next nearest neighbors
                    % Note that we may add this particle back later if bb[dmyi - 1] is its own NNN
                    neighbors = np.delete(n_tmp, np.where(n_tmp == bb[dmyi - 1])[0])
                    % Here, handle the case where a periodic bond links the neighbor back to the original particle,
                    % as in the bond linkage of 0-1-0.
                    if len(neighbors) == 0:
                        neighbors = n_tmp

                % check if neighbors CAN be connected across periodic bc--
                %  ie if particle on finite boundary (finbd)
                if nxt in finbd:
                    % Since on finite system boundary, particle could have periodic bonds
                    % Find x values to add to neighbors, by first getting indices of row of
                    % PV (same as of NL) matching neighbors
                    % PVinds = [np.argwhere(NL[nxt] == nnn)[0][0] for nnn in neighbors] <--- this assumed no 0-1-0
                    PVinds = []
                    for nnn in dh.unique_nosort(neighbors)
                        okinds = np.ravel(np.argwhere(np.logical_and(NL[nxt] == nnn, np.abs(KL[nxt]) > eps)))
                        % print 'neighbors = ', neighbors
                        % print 'okinds = ', okinds
                        % print 'NL = ', NL
                        % print 'KL = ', KL
                        % print NL[nxt] == nnn, np.abs(KL[nxt]) > eps
                        % print np.argwhere(np.logical_and(NL[nxt] == nnn, np.abs(KL[nxt]) > eps))
                        for okind in okinds:
                            PVinds.append(okind)

                    addx = PVx[nxt, PVinds]
                    addy = PVy[nxt, PVinds]

                    % print 'nxt = ', nxt
                    % print 'PVinds', PVinds
                    % print 'xy[neighbors, :] = ', xy[neighbors, :]
                    % print 'np.dstack([addx, addy])[0] = ', np.dstack([addx, addy])[0]

                    xynb = xy[neighbors, :] + np.dstack([addx, addy])[0]
                    xynxt = xy[nxt, :]
                    current_angles = np.arctan2(xynb[:, 1] - xynxt[1], xynb[:, 0] - xynxt[0]).ravel()
                    angles = np.mod(current_angles - prev_angle, 2 * np.pi)

                    if check:
                        print '\n'
                        print 'particle ', nxt, ' is on finbd'
                        print 'nxt = ', nxt
                        print 'neighbors = ', neighbors
                        print 'xy[neighbors,:] =', xy[neighbors, :]
                        print 'addxy = ', np.dstack([addx, addy])[0]
                        print 'xynb = ', xynb
                        print 'xynxt = ', xynxt
                        print 'current_angles = ', current_angles
                        print 'prev_angle = ', prev_angle
                        print 'angles = ', angles
                        print 'redefining nxt = ', neighbors[angles == max(angles)][0]

                    % redefine previous angle as backwards of current angle -- ie angle(prev-current_pos)
                    prev_angletmp = np.arctan2(xynxt[1] - xynb[:, 1], xynxt[0] - xynb[:, 0]).ravel()
                    prev_angle = prev_angletmp[angles == max(angles)][0]

                    % CHECK
                    % ax1 = plt.gca()
                    % ax1.plot(xy[:,0],xy[:,1],'k.')
                    % for i in range(len(xy))
                    %    ax1.text(xy[i,0]+0.2,xy[i,1],str(i))
                    % plt.show()

                else
                    current_angles = np.arctan2(xy[neighbors, 1] - xy[nxt, 1],
                                                xy[neighbors, 0] - xy[nxt, 0]).ravel()
                    angles = np.mod(current_angles - prev_angle, 2 * np.pi)
                    % redefine previous angle as backwards of current angle -- ie angle(prev-current_pos)
                    % prev_angle = np.arctan2(xy[bb[dmyi-1],1] - xynxt[1], xy[bb[dmyi-1],0] - xynxt[0] ).ravel()
                    xynxt = xy[nxt, :]
                    xynb = xy[neighbors, :]
                    prev_angletmp = np.arctan2(xynxt[1] - xy[neighbors, 1], xynxt[0] - xy[neighbors, 0]).ravel()
                    prev_angle = prev_angletmp[angles == max(angles)][0]

                nxt = neighbors[angles == max(angles)][0]
                bb.append(nxt)
                % update displacement
                displ += xynb[angles == max(angles)][0] - xynxt

                %%%%%%%%%%%%%%%
                % Check bond
                if viewmethod
                    % Check individually
                    % ax1 = plt.gca()
                    % ax1.plot(xy[:,0],xy[:,1],'k.')
                    if first_check:
                        for i in range(len(xy))
                            ax1.text(xy[i, 0] + 0.2, xy[i, 1], str(i))

                    plt.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                                 xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                                 arrowprops=dict(arrowstyle="->",
                                                 color="r",
                                                 shrinkA=5, shrinkB=5,
                                                 patchA=None,
                                                 patchB=None,
                                                 connectionstyle="arc3,rad=0.2", ), )

                %%%%%%%%%%%%%%%

                % Now mark the current bond as used
                % thisbond = [bb[dmyi-1], bb[dmyi]]
                % Get index of used matching thisbond
                mark_used = np.where((np.logical_or(BL == bb[dmyi - 1], BL == bb[dmyi])).all(axis=1))[0]
                % mark_used = np.where((BL == thisbond).all(axis=1))
                % print 'mark_used = ', mark_used
                % I adjusted the line below to allow multiple entries in mark_used (2018-04-26)'
                if ~(used[mark_used, 0]).all()
                    % print 'marking bond [', thisbond, '] as used'
                    marking, kk = true, 0
                    while marking:
                        if ~used[mark_used[kk], 0]:
                            used[mark_used[kk], 0] = true
                            marking = false
                        kk += 1
                else
                    % Get index of used matching reversed thisbond (this list boolean is directional)
                    % mark_used = np.where((BL == thisbond[::-1]).all(axis=1))
                    % Used this bond in reverse order
                    marking, kk = true, 0
                    while marking:
                        print 'mark_used = ', mark_used
                        print 'mark_used[kk] = ', mark_used[kk]
                        print 'used[mark_used[kk]] = ', used[mark_used[kk]]
                        print '--------------------------'
                        if ~used[mark_used[kk], 1]:
                            used[mark_used[kk], 1] = true
                            marking = false
                        % except IndexError:
                        %     print 'mark_used = ', mark_used
                        %     print 'used[mark_used] = ', used[mark_used[kk]]
                        %     print 'marking bond ', BL[mark_used[kk]]
                        %     print 'kk = ', kk
                        %     print 'bb = ', bb
                        %     print 'Encountered index error in marking bond used'
                        %     plt.show()
                        %     sys.exit()
                        kk += 1
                        if kk == len(mark_used)
                            marking = false

                % print 'used = ', used
                dmyi += 1
                if check:
                    print 'bb = ', bb

            polygons.append(bb)
            %%%%%%%%%%%%%%%
            % Check new polygon
            if viewmethod
                if first_check:
                    ax1.plot(xy[:, 0], xy[:, 1], 'k.')
                    for i in range(len(xy))
                        ax1.text(xy[i, 0] + 0.2, xy[i, 1], str(i))

                for dmyi in range(len(bb))
                    nxt = bb[np.mod(dmyi + 1, len(bb))]
                    ax1.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                                 xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                                 arrowprops=dict(arrowstyle="->",
                                                 color="r",
                                                 shrinkA=5, shrinkB=5,
                                                 patchA=None,
                                                 patchB=None,
                                                 connectionstyle="arc3,rad=0.2", ), )
                ax2.cla()
                ax2.imshow(used, aspect=1. / len(used), interpolation='none')
                print 'polygons = ', polygons
                % plt.show()
                plt.pause(0.00001)
                %%%%%%%%%%%%%%%

        else
            % Check for remaining bonds unused in reverse order (B-->A)
            % print 'CHECKING REVERSE (B-->A): '
            todoBA = np.where(~used[:, 1])[0]
            % print 'len(todoBA) = ', len(todoBA)
            if len(todoBA) > 0
                bond = BL[todoBA[0]]

                %%%%%%%%%%%%%%%
                % % check
                if viewmethod
                    plt.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                                 xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                                 arrowprops=dict(arrowstyle="->",
                                                 color="b",
                                                 shrinkA=5, shrinkB=5,
                                                 patchA=None,
                                                 patchB=None,
                                                 connectionstyle="arc3,rad=0.6", ), )
                % %%%%%%%%%%%%%%%

                % bb will be list of polygon indices
                % Start with orientation going from bond[0] to bond[1]
                nxt = bond[0]
                bb = [bond[1], nxt]
                dmyi = 1

                % define 'previous angle' as backwards of current angle -- ie angle(prev-current_pos)
                % Must include effect of PV on this angle -- do in ref frame of nxt particle
                PVind = np.argwhere(NL[nxt] == bond[1])[0][0]
                addx = PVx[nxt, PVind]
                addy = PVy[nxt, PVind]
                xyb0 = xy[bond[1], :] + np.array([addx, addy])
                prev_angle = np.arctan2(xyb0[1] - xy[nxt, 1], xyb0[0] - xy[nxt, 0])  % .ravel()

                % as long as we haven't completed the full outer polygon, add nextIND
                % define the displacment from the starting point that we have moved so far
                displ = xy[nxt] - xyb0

                % as long as we haven't completed the full outer polygon, add next index
                while nxt != bond[1] or abs(displ[0] ** 2 + displ[1] ** 2) > eps:
                    n_tmp = NL[nxt, np.argwhere(KL[nxt]).ravel()]
                    % Exclude previous boundary particle from the neighbors array, unless its the only one
                    % (It cannot be the only one, if we removed dangling bonds)
                    if len(n_tmp) == 1:
                        '''The bond is a lone bond, not part of a triangle.'''
                        neighbors = n_tmp
                    else
                        neighbors = np.delete(n_tmp, np.where(n_tmp == bb[dmyi - 1])[0])
                        % Add neighbors back in if this bond is not dangling but we have a NNN structure of 0-1-0
                        if len(neighbors) == 0:
                            neighbors = n_tmp

                    %%%%%%%%
                    % check if neighbors CAN be connected across periodic bc-- ie if particle is
                    % on the finite boundary (finbd)
                    if nxt in finbd:
                        % Since on finite system boundary, particle could have periodic bonds
                        % Find x values to add to neighbors, by first getting indices of row of PV
                        % (same as of NL) matching neighbors
                        % ALL CALCS in frame of reference of NXT particle
                        % PVinds = [np.argwhere(NL[nxt] == nnn)[0][0] for nnn in neighbors]
                        PVinds = []
                        for nnn in dh.unique_nosort(neighbors)
                            okinds = np.ravel(np.argwhere(np.logical_and(NL[nxt] == nnn, np.abs(KL[nxt]) > eps)))
                            for okind in okinds:
                                PVinds.append(okind)

                        addx = PVx[nxt, PVinds]
                        addy = PVy[nxt, PVinds]

                        xynb = xy[neighbors, :] + np.dstack([addx, addy])[0]
                        xynxt = xy[nxt, :]
                        % print '\n'
                        % print 'nxt = ', nxt
                        % print 'neighbors = ', neighbors
                        % print 'xy[neighbors,:] =', xy[neighbors,:]
                        % print 'addxy = ', np.dstack([addx, addy])[0]
                        % print 'xynb = ', xynb
                        % print 'xynxt = ', xynxt
                        current_angles = np.arctan2(xynb[:, 1] - xynxt[1], xynb[:, 0] - xynxt[0]).ravel()
                        angles = np.mod(current_angles - prev_angle, 2 * np.pi)
                        selectIND = np.where(angles == max(angles))[0][0]
                        % print 'selectIND = ', selectIND
                        % print 'current_angles = ', current_angles/np.pi
                        % print 'prev_angle = ', prev_angle/np.pi
                        % print 'angles = ', angles/np.pi

                        % redefine previous angle as backwards of current angle -- ie angle(nxt - neighbor )
                        prev_angletmp = np.arctan2(xynxt[1] - xynb[:, 1], xynxt[0] - xynb[:, 0]).ravel()
                        prev_angle = prev_angletmp[selectIND]

                        % print 'new prev_angle = ', prev_angle/np.pi
                        % print 'NL[nxt] = ', NL[nxt]
                        % print 'bb = ', bb
                        % % CHECK
                        % ax1 = plt.gca()
                        % ax1.plot(xy[:,0],xy[:,1],'k.')
                        % for i in range(len(xy))
                        %   ax1.text(xy[i,0]+0.2,xy[i,1],str(i))
                        % plt.arrow(xynxt[0], xynxt[1], np.cos(angles[selectIND]),
                        %           np.sin(angles[selectIND]),fc='r', ec='r')
                        % plt.arrow(xynb[selectIND,0], xynb[selectIND,1],
                        %           np.cos(prev_angle), np.sin(prev_angle),fc='b', ec='b')
                        % plt.show()

                    else
                        current_angles = np.arctan2(xy[neighbors, 1] - xy[nxt, 1],
                                                    xy[neighbors, 0] - xy[nxt, 0]).ravel()
                        angles = np.mod(current_angles - prev_angle, 2 * np.pi)
                        % redefine previous angle as backwards of current angle -- ie angle(prev-current_pos)
                        xynxt = xy[nxt, :]
                        xynb = xy[neighbors, :]
                        prev_angletmp = np.arctan2(xynxt[1] - xynb[:, 1], xynxt[0] - xynb[:, 0]).ravel()
                        selectIND = np.where(angles == max(angles))[0][0]
                        % print '\n'
                        % print 'nxt = ', nxt
                        % print 'bb = ', bb
                        % print 'neighbors = ', neighbors
                        % print 'current_angles = ', current_angles/np.pi
                        % print 'prev_angle = ', prev_angle/np.pi
                        % print 'angles = ', angles/np.pi
                        % print 'selectIND = ', selectIND
                        % print('xynxt[1] - xynb[:,1], xynxt[0] - xynb[:,0] = ', xynxt[1] - xynb[:,1],
                        %       xynxt[0] - xynb[:,0])
                        % print('np.arctan2(xynxt[1] - xynb[:,1], xynxt[0] - xynb[:,0]) = ',
                        %       np.arctan2(xynxt[1] - xynb[:,1], xynxt[0] - xynb[:,0]))
                        % print 'prev_angletmp = ', prev_angletmp/np.pi

                        prev_angle = prev_angletmp[selectIND]
                        % print 'new prev_angle = ', prev_angle/np.pi

                    %%%%%%%%%%%%%%%
                    nxt = neighbors[angles == max(angles)][0]
                    bb.append(nxt)
                    % update displacement of particle at nxt from first site (keeping track of periodic bonds)
                    displ += xynb[angles == max(angles)][0] - xynxt

                    %%%%%%%%%%%%%%%
                    % Check
                    if viewmethod
                        % If checking individual bonds
                        % ax1 = plt.gca()
                        % ax1.plot(xy[:,0],xy[:,1],'k.')
                        % for i in range(len(xy))
                        %    ax1.text(xy[i,0]+0.2,xy[i,1],str(i))

                        plt.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                                     xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                                     arrowprops=dict(arrowstyle="->",
                                                     color="b",
                                                     shrinkA=5, shrinkB=5,
                                                     patchA=None,
                                                     patchB=None,
                                                     connectionstyle="arc3,rad=0.6",
                                                     ), )
                        % plt.show()
                    %%%%%%%%%%%%%%%

                    % Now mark the current bond as used --> note the inversion of the bond order to match BL
                    thisbond = [bb[dmyi], bb[dmyi - 1]]
                    % Get index of used matching [bb[dmyi-1],nxt]
                    mark_used = np.where((BL == thisbond).all(axis=1))
                    if len(mark_used) > 0
                        used[mark_used, 1] = true
                    else
                        messg = 'Cannot mark polygon bond as used: this bond was already used in its attempted' + \
                                ' orientation. (All bonds in first column should already be marked as used.)'
                        raise RuntimeError(messg)

                    dmyi += 1

                polygons.append(bb)
                % print 'added polygon = ', bb

                % Check new polygon
                if viewmethod
                    if first_check:
                        ax1.plot(xy[:, 0], xy[:, 1], 'k.')
                        for i in range(len(xy))
                            ax1.text(xy[i, 0] + 0.2, xy[i, 1], str(i))

                    for dmyi in range(len(bb))
                        nxt = bb[np.mod(dmyi + 1, len(bb))]
                        ax1.annotate("", xy=(xy[bb[dmyi], 0], xy[bb[dmyi], 1]), xycoords='data',
                                     xytext=(xy[nxt, 0], xy[nxt, 1]), textcoords='data',
                                     arrowprops=dict(arrowstyle="->",
                                                     color="b",
                                                     shrinkA=5, shrinkB=5,
                                                     patchA=None,
                                                     patchB=None,
                                                     connectionstyle="arc3,rad=0.6", ), )
                    ax2.cla()
                    ax2.imshow(used)
                    % plt.show()
                    plt.pause(0.0001)
                    %%%%%%%%%%%%%%%

            else
                % All bonds have been accounted for
                print 'all finished with finding polygons...'
                finished = true
% check
if viewmethod
    plt.show()

% Check for duplicates (up to cyclic permutations and inversions) in polygons
% Note that we need to ignore the last element of each polygon (which is also starting pt)
keep = np.ones(len(polygons), dtype=bool)
for ii in range(len(polygons))
    print 'ii = ', ii
    polyg = polygons[ii]
    for p2 in polygons[ii + 1:]:
        if is_cyclic_permutation(polyg[:-1], p2[:-1])
            keep[ii] = false

polygons = [polygons[i] for i in np.where(keep)[0]]

% Remove duplicates via inversion (maybe not necessary?)

% Remove the polygon which is the entire lattice boundary, except dangling bonds
if ~periB.any()
    print 'le.extract_polygons_lattice: Removing entire lattice boundary from list of polygons...'
    boundary = extract_boundary(xy, NL, KL, BL)
    % print 'boundary = ', boundary
    keep = np.ones(len(polygons), dtype=bool)
    for ii in range(len(polygons))
        polyg = polygons[ii]
        if is_cyclic_permutation(polyg[:-1], boundary.tolist())
            keep[ii] = false
        elif is_cyclic_permutation(polyg[:-1], boundary[::-1].tolist())
            keep[ii] = false

    polygons = [polygons[i] for i in np.where(keep)[0]]

% Check order of each polygon so that it is oriented counterclockwise
% for polys in polygons:
%     angle_poly = 0
%     % Make sure that oriented counterclockwise
%     print 'polys = ', polys
%     for i in range(len(polys))
%         p0 = polys[ np.mod(i-1, len(polys)-1)]
%         p1 = polys[i]
%         p2 = polys[ np.mod(i+1,len(polys)-1) ]
%         print 'p0,p1,p2 = ', p0, p1, p2
%         angle_tmp = np.mod(np.arctan2(xy[p2,1]-xy[p1,1], xy[p2,0]-xy[p1,0]) - np.arctan2( xy[p1,1]-xy[p0,1],
%                            xy[p1,0]-xy[p0,0] ), 2*np.pi)
%         print 'angle_tmp = ', angle_tmp
%         angle_poly += angle_tmp
%
%     print 'angle = ', angle_poly/6.
print 'le: polygons = ', polygons
if check
    polygons2PPC(xy, polygons, BL=BL, PVxydict=PVxydict, check=true)
end

