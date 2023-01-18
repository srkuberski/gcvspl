function [z, s] = splzer( x, c, m, n, xmin, xmax, ofs )
% spline zeros
%
% [z, s] = splzer( x, c, m, n, xmin = min( x ), xmax = max( x ), ofs = 0 )
%
% INPUT
% x : independent variables [sorted] (numeric row [1, N])
% c : spline coefficients (numeric row [1, N])
% m : spline half order [1-4: linear, cubic, quintic, heptic] (numeric scalar)
% n : order of derivative [0 for spline values] (numeric scalar)
% xmin, xmax : search interval (numeric scalar)
% ofs : derivative offset (numeric scalar)
%
% OUTPUT
% z : zeros of the spline (numeric row [1, Z])
% s : signs of the zeros (numeric row [1, Z])

		% safeguard
	if nargin < 1 || ~isnumeric( x ) || ~isrow( x )
		error( 'invalid argument: x (numeric row [1, N])' );
	end
	N = numel( x );

	if nargin < 2 || ~isnumeric( c ) || ~isrow( c ) || numel( c ) ~= N
		error( 'invalid argument: c (numeric row [1, N])' );
	end

	if nargin < 3 || ~isnumeric( m ) || ~isscalar( m ) || ~ismember( m, 1:4 )
		error( 'invalid argument: m (numeric scalar)' );
	end

	if nargin < 4 || ~isnumeric( n ) || ~isscalar( n ) || ~ismember( n, 0:2*m-1 )
		error( 'invalid argument: n (numeric scalar)' );
	end

	if nargin < 5
		xmin = min( x );
	end
	if ~isnumeric( xmin ) || ~isscalar( xmin )
		error( 'invalid argument: xmin (numeric scalar)' );
	end

	if nargin < 6
		xmax = max( x );
	end
	if ~isnumeric( xmax ) || ~isscalar( xmax )
		error( 'invalid argument: xmax (numeric scalar)' );
	end

	if nargin < 7
		ofs = 0;
	end
	if ~isnumeric( ofs ) || ~isscalar( ofs )
		error( 'invalid argument: ofs (numeric scalar)' );
	end

		% adjust search interval
	imin = x(1+2*m-2);
	imax = x(end-2*m-2);
	xmin = aux.clamp( xmin, imin, imax );
	xmax = aux.clamp( xmax, imin, imax );

		% zeros
	sp = fnder( spmak( augknt( x, m+1 ), c ), n );
	sp.coefs = sp.coefs+ofs;

	z = fnzeros( sp, [xmin, xmax] );
	z = transpose( unique( z(:) ) );

		% signs
	sp = fnder( sp, 1 );
	s = sign( fnval( sp, z ) );

end % function

