function [z, s] = splzer( x, c, m, n, xmin, xmax, ofs )
% compute spline zeros
%
% [z, s] = splzer( x, c, m, n, xmin = min( x ), xmax = max( x ), ofs = 0 )
%
% INPUT
% x : independent variables (numeric row [1, N])
% c : spline coefficients (numeric row [1, N])
% m : spline half order (numeric scalar)
% n : order of derivative (numeric scalar)
% xmin, xmax : search interval (numeric scalar)
% ofs : spline offset (numeric scalar)
%
% OUTPUT
% z : zero locations (numeric row [1, Z])
% s : zero signs (numeric row [1, Z])
%
% REMARKS
% - array of independent variables x must be sorted
% - spline half orders m = 1, 2, 3, 4 correspond to linear, cubic, quintic, heptic splines (etc.)
% - physically, orders of derivate n = 0, 1, 2 correspond to position, velocity, acceleration values (etc.)
% - this function requires the presence of the Curve Fitting Toolbox
%
% REQUIREMENTS
% - Matlab's Curve Fitting Toolbox must be installed on your computer

		% safeguard
	if nargin < 1 || ~isnumeric( x ) || ~isrow( x )
		error( 'invalid argument: x (numeric row [1, N])' );
	end
	N = numel( x );

	if nargin < 2 || ~isnumeric( c ) || ~isrow( c ) || numel( c ) ~= N
		error( 'invalid argument: c (numeric row [1, N])' );
	end

	if nargin < 3 || ~isnumeric( m ) || ~isscalar( m )
		error( 'invalid argument: m (numeric scalar)' );
	end

	if nargin < 4 || ~isnumeric( n ) || ~isscalar( n )
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
	xmin = max( xmin, x(1+2*m-2) );
	xmax = min( xmax, x(end-2*m-2) );

		% compute zero locations
	sp = fnder( spmak( augknt( x, m+1 ), c ), n );
	sp.coefs = sp.coefs+ofs;

	z = fnzeros( sp, [xmin, xmax] );
	z = transpose( unique( z(:) ) );

		% compute zero signs
	sp = fnder( sp, 1 );
	s = sign( fnval( sp, z ) );

end % function

