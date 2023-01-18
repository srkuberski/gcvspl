function [c, wk] = gcvspl( x, y, m, v, w )
% c = gcvspl( x, y, m, v, w=ones )
% 
% INPUT
% x : independent variables [sorted] (double row [1, N])
% y : data to be smoothed (double matrix [K, N])
% m : spline half order [1-4: linear, cubic, quintic, heptic] (double scalar)
% v : prior given variance [negative: GCV regularization] (double scalar)
% w : weight factors (double row [1, N])
%
% OUTPUT
% c : spline coefficients (double matrix [K, N])
% wk : internal work vector (numeric row [1, 6])

		% safeguard
	if nargin < 1 || ~isnumeric( x ) || ~isrow( x )
		error( 'invalid argument: x (double row [1, N])' );
	end
	N = numel( x );

	if nargin < 2 || ~isnumeric( y ) || size( y, 2 ) ~= N
		error( 'invalid argument: y (double matrix [K, N])' );
	end
	K = size( y, 1 );

	if nargin < 3 || ~isnumeric( m ) || ~isscalar( m )
		error( 'invalid argument: m (double scalar)' );
	end

	if nargin < 4 || ~isnumeric( v ) || ~isscalar( v )
		error( 'invalid argument: v (double scalar)' );
	end

	if nargin < 5
		w = ones( size( x ) );
	end
	if ~isnumeric( w ) || ~isrow( w ) || numel( w ) ~= N
		error( 'invalid argument: w (double row [1, N])' );
	end

		% compute spline coefficients
	[c, wk] = gcvsplmex( transpose( x ), transpose( y ), m, v, transpose( N*w/sum( w ) ) );
	c = transpose( c );

end % function

