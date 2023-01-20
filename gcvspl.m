function [c, wk] = gcvspl( x, y, m, v, w )
% compute spline coefficients
%
% c = gcvspl( x, y, m, v, w = ones )
% 
% INPUT
% x : independent variables (numeric row [1, N])
% y : data to be smoothed (numeric matrix [K, N])
% m : spline half order (numeric scalar)
% v : prior given variance [negative: GCV regularization] (numeric scalar)
% w : weight factors (numeric row [1, N])
%
% OUTPUT
% c : spline coefficients (numeric matrix [K, N])
% wk : internal work vector (numeric row [1, 6])
%
% REMARKS
% - array of independent variables x must be sorted
% - spline half orders m = 1, 2, 3, 4 correspond to linear, cubic, quintic, heptic splines (etc.)
% - negative prior given variance v results in generalized, cross-validatory splines
%
% REQUIREMENTS
% - the binary gcvsplmex must be accessible from Matlab's search path

		% safeguard
	if nargin < 1 || ~isnumeric( x ) || ~isrow( x )
		error( 'invalid argument: x (numeric row [1, N])' );
	end
	N = numel( x );

	if nargin < 2 || ~isnumeric( y ) || size( y, 2 ) ~= N
		error( 'invalid argument: y (numeric matrix [K, N])' );
	end
	K = size( y, 1 );

	if nargin < 3 || ~isnumeric( m ) || ~isscalar( m )
		error( 'invalid argument: m (numeric scalar)' );
	end

	if nargin < 4 || ~isnumeric( v ) || ~isscalar( v )
		error( 'invalid argument: v (numeric scalar)' );
	end

	if nargin < 5
		w = ones( size( x ) );
	end
	if ~isnumeric( w ) || ~isrow( w ) || numel( w ) ~= N
		error( 'invalid argument: w (numeric row [1, N])' );
	end

		% compute spline coefficients
	x = transpose( double( x ) );
	y = transpose( double( y ) );
	m = double( m );
	v = double( v );
	w = transpose( double( w ) );

	[c, wk] = gcvsplmex( x, y, m, v, N*w/sum( w ) );

	c = transpose( c );

end % function

