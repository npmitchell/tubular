# spring_networks_matlab
Codes developed by William Irvine (UChicago) and I for spring network energetics and evolution


Common Variables
================
% pts : N x 3 float array
%   3D positions of points. Row i is the position of the ith particle.
% NL : array of dimension #pts x max(#neighbors)
%   The NL(i,:) contains indices for the neighbors for the ith point, 
%   buffered by zeros if a particle does not have the maximum # nearest 
%   neighbors.
% KL : NP x max(#neighbors) int array
%   spring connection/constant list, where 1 corresponds to a true connection,
%   0 signifies that there is not a connection, -1 signifies periodic bond
% TRI : #faces x 3 int array
%   triangle array, with TRI(i, :) being indices into pts of ith face
% TRISm : 
BM : array of length #pts x max(#neighbors)
    The (i,j)th element is the bond length of the bond connecting the ith particle to its jth neighbor (the particle with index NL[i,j]).
BL : array of dimension #bonds x 2
    Each row is a bond and contains indices of connected points. Negative values denote particles connected through periodic bonds.
bL : array of length #bonds
    The ith element is the length of of the ith bond in BL
LL : tuple of 2 floats
    Horizontal and vertical extent of the bounding box (a rectangle) through which there are periodic boundaries.
    These give the dimensions of the network in x and y, for S(k) measurements and other periodic things.
BBox : #vertices x 2 float array
    bounding polygon for the network, usually a rectangle
lp : dict
    The lattice parameters dictionary, with all keys needed for specifying path (these params vary depending on the
    value of the LatticeTop key).
polygons : list of int lists
    indices of xy points defining polygons.
NLNNN : array of length #pts x max(#next-nearest-neighbors)
    Next-nearest-neighbor array: The ith row contains indices for the next nearest neighbors for the ith point.
KLNNN : array of length #pts x max(#next-nearest-neighbors)
    Next-nearest-neighbor connectivity/orientation array:
    The ith row states whether a next nearest neighbors is counterclockwise (1) or clockwise (-1)
PVxydict : dict
    dictionary of periodic bonds (keys) to periodic vectors (values)
    If key = (i,j) and val = np.array([ 5.0, 2.0]), then particle i sees particle j at xy[j]+val
    --> transforms into:  ijth element of PVx is the x-component of the vector taking NL[i,j] to its image as 
    seen by particle i
    If network is a small periodic unit like a unit cell, such that one particle i sees particle j more than once,
    then PVxydict[(i, j)] has multiple vectors specifying each image, like
    ex. PVxydict = {(i, j): np.array([first_pv, second_pv, ...])
    In this case, denote that this is a unitcell with lp['unitcell'] == True to catch exceptions
PVx : NP x NN float array (optional, for periodic lattices)
    ijth element of PVx is the x-component of the vector taking NL[i,j] to its image as seen by particle i 
    Note that shape(PVx) == shape(NL)
PVy : NP x NN float array (optional, for periodic lattices)
    ijth element of PVy is the y-component of the vector taking NL[i,j] to its image as seen by particle i
    Note that shape(PVy) == shape(NL)
PVxnp : NP x NP float array
    ijth element of PVy is the x-component of the vector taking j = NL[i,k] to its image as seen by particle i
    Note that shape(PVxnp) == (len(xy), len(xy))
PVynp : NP x NP float array
    ijth element of PVy is the y-component of the vector taking j = NL[i,k] to its image as seen by particle i
    Note that shape(PVynp) == (len(xy), len(xy))
nljnnn : #pts x max(#NNN) int array
    nearest neighbor array matching NLNNN and KLNNN. nljnnn[i, j] gives the neighbor of i such that NLNNN[i, j] is
    the next nearest neighbor of i through the particle nljnnn[i, j]
kljnnn : #pts x max(#NNN) int array
    bond array describing periodicity of bonds matching NLNNN and KLNNN. kljnnn[i, j] describes the bond type
    (bulk -> +1, periodic --> -1) of bond connecting i to nljnnn[i, j]
klknnn : #pts x max(#NNN) int array
    bond array describing periodicity of bonds matching NLNNN and KLNNN. klknnn[i, j] describes the bond type
    (bulk -> +1, periodic --> -1) of bond connecting nljnnn[i, j] to NLNNN[i, j]
pvxnnn : NP x max(#NNN) float array
    ijth element of pvxnnn is the x-component of the vector taking NLNNN[i,j] to its image as seen by particle i
pvynnn : NP x max(#NNN) float array
    ijth element of pvynnn is the y-component of the vector taking NLNNN[i,j] to its image as seen by particle i
PV : 2 x 2 float array
    periodic lattice vectors, with x-dominant vector first, y-dominant vector second.
gxy : tuple of NX x NY arrays
    two-point correlation function in the positions of particles as function of vector distance x,y
gr :
    two-point correlation function in the positions of particles as function of scalar distance

