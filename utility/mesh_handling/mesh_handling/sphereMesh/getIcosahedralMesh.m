function [ P, tri ] = getIcosahedralMesh()
	
	phi = ( 1 + sqrt( 5 ) ) / 2;
	P = [ ...
		-1 phi 0 ; ...
		1 phi 0 ; ...
		-1 -phi 0 ; ...
		1 -phi 0 ; ...
		0 -1 phi ; ...
		0 1 phi ; ...
		0 -1 -phi ; ...
		0 1 -phi ; ...
		phi 0 -1 ; ...
		phi 0 1 ; ...
		-phi 0 -1 ; ...
		-phi 0 1   ...
	]';
	P = P ./ repmat( sqrt( sum( P.^2 ) ), 3, 1 );
	P = P - repmat( mean( P, 2 ), 1, size( P, 2 ) );
	P = randomlyRotate( P );
	tri = [ ...
		1 12 6 ; ...
		1 6 2 ; ...
		1 2 8 ; ...
		1 8 11 ; ...
		1 11 12 ; ...
		2 6 10 ; ...
		6 12 5 ; ...
		12 11 3 ; ...
		11 8 7 ; ...
		8 2 9 ; ...
		4 10 5 ; ...
		4 5 3 ; ...
		4 3 7 ; ...
		4 7 9 ; ...
		4 9 10 ; ...
		5 10 6 ; ...
		3 5 12 ; ...
		7 3 11 ; ...
		9 7 8 ; ...
		10 9 2  ...
	];

end