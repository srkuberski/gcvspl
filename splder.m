function y = splder( x, c, m, t, n )
% compute spline values
%
% y = splder( x, c, m, t, n )
%
% INPUT
% x : independent variables (numeric row [1, N])
% c : spline coefficients (numeric matrix [K, N])
% m : spline half order (numeric scalar)
% t : evaluation points (numeric row [1, T])
% n : order of derivative (numeric scalar)
%
% OUTPUT
% y : spline values (numeric matrix [K, T])
%
% REMARKS
% - array of independent variables x must be sorted
% - spline half orders m = 1, 2, 3, 4 correspond to linear, cubic, quintic, heptic splines (etc.)
% - physically, orders of derivate n = 0, 1, 2 correspond to position, velocity, acceleration values (etc.)
%
% REQUIREMENTS
% - the binary spldermex must be accessible from Matlab's search path

		% safeguard
	if nargin < 1 || ~isnumeric( x ) || ~isrow( x )
		error( 'invalid argument: x (numeric row [1, N])' );
	end
	N = numel( x );

	if nargin < 2 || ~ismatrix( c ) || size( c, 2 ) ~= N
		error( 'invalid argument: c (numeric matrix [K, N])' );
	end
	K = size( c, 1 );

	if nargin < 3 || ~isnumeric( m ) || ~isscalar( m )
		error( 'invalid argument: m (numeric scalar)' );
	end

	if nargin < 4 || ~isnumeric( t ) || ~isrow( t )
		error( 'invalid argument: t (numeric row [1, T])' );
	end

	if nargin < 5 || ~isnumeric( n ) || ~isscalar( n )
		error( 'invalid argument: n (numeric scalar' );
	end

		% compute spline values
	x = transpose( double( x ) );
	c = transpose( double( c ) );
	m = double( m );
	t = transpose( double( t ) );
	n = double( n );

	for k = 1:K
		y(:, k) = spldermex( x, c(:, k), m, t, n );
	end

	y = transpose( y );

end

