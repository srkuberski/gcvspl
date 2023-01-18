function y = splder( x, c, m, t, n )
% y = splder( x, c, m, t, n )
%
% INPUT
% x : independent variables [sorted] (numeric row [1, N])
% c : spline coefficients (numeric matrix [K, N])
% m : spline half order [1-4: linear, cubic, quintic, heptic] (numeric scalar)
% t : evaluation points (numeric row [1, T])
% n : order of derivative [0 for spline values] (numeric scalar)
%
% OUTPUT
% y : spline values or derivatives (numeric matrix [K, T])

	x = transpose( x );
	c = transpose( c );
	t = transpose( t );
	for k = 1:size( c, 2 )
		y(:, k) = spldermex( x, c(:, k), m, t, n );
	end
	y = transpose( y );

end

