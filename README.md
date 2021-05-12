# GCV spline smoothing

Stephan R. Kuberski, Department of Linguistics, University of Potsdam, November, 2019

This is a version of Woltring's classic generalized, cross-validatory (GCV) spline smoothing and differentiation code [1]. The original Fortran 77 code [2] was converted to C using the f2c converter [3]. The resulting C code `gcvspl.c` (as well as its Fortran 77 source `gcvspl.f`) is stored here, along with this package, for archiving purposes. Please adhere to their original copyright.

To make the code accessible to a wider range of Matlab users the following two MEX wrappers have been implemented in C. Both files have been developed and tested in Matlab 9.6.0 (R2019a) but should properly run on older releases too.

## gcvsplmex.c

The function computes a natural B-spline using the generalized cross-validation and mean-squared prediction error criteria of Craven & Wahba [4]. The model assumes uncorrelated, additive noise and essentially smooth, underlying functions. The independent coordinates may be spaced non-equidistantly.

In order to use this function, you need to compile it first (and only once).  This usually requires the installation of a C/C++ compiler (e.g., the GNU Compiler Collection gcc) using your operating system's software center. Once a C compiler is available on your system, start Matlab, change to the path where the files of this package are located and run:

~~~matlab
mex gcvsplmex.c gcvspl.c
~~~

The resulting binary `gcvsplmex.mexa64` (or similar on non-Unix systems) is a
Matlab executable with the following calling convention:

~~~matlab
% c = gcvsplmex( x, y, m, v )
% 
% INPUT
% x : independent variables [sorted] (double column [N, 1])
% y : data to be smoothed (double matrix [N, K])
% m : spline half order [1-4: linear, cubic, quintic, heptic] (double scalar)
% v : prior given variance [negative: GCV regularization] (double scalar)
%
% OUTPUT
% c : spline coefficients (double matrix [N, K])
~~~

## spldermex.c

The function computes the values of a B-spline or its derivatives at selected evaluation points. Woltring's underlying code is an adoption of Lyche *et al.* [5].

In order to use this function, you need to compile it first (and only once).  This usually requires the installation of a C/C++ compiler (e.g., the GNU Compiler Collection gcc) using your operating system's software center. Once a C compiler is available on your system, start Matlab, change to the path where the files of this package are located and run:

~~~matlab
mex spldermex.c gcvspl.c
~~~

The resulting binary `spldermex.mexa64` (or similar on non-Unix systems) is a
Matlab executable with the following calling convention:

~~~matlab
% y = spldermex( x, c, m, t, n )
%
% INPUT
% x : independent variables [sorted] (double column [N, 1])
% c : spline coefficients (double column [N, 1])
% m : spline half order [1-4: linear, cubic, quintic, heptic] (double scalar)
% t : evaluation points (double column [T, 1])
% n : order of derivative [0 for spline values] (double scalar)
%
% OUTPUT
% y : spline values or derivatives (double column [T, 1])
~~~

## Example

As an example of usage, consider the following lines of code, which compute the (fifth-order) quintic spline values and their first two derivatives (velocity and acceleration) for a given set of data points `(x,y)`:

~~~matlab
c = gcvsplmex( x, y, 3 ); % compute GCV spline coefficients
y0 = spldermex( x, c, 3, x, 0 ); % compute spline values at x
y1 = spldermex( x, c, 3, x, 1 ); % compute first derivatives at x (velocity)
y2 = spldermex( x, c, 3, x, 2 ); % compute second derivatives at x (acceleration)
~~~

In order to run the code above, make sure that the two binaries mentioned earlier (`gcvsplmex.mexa64` and `spldermex.mexa64`) are placed in a path which is known to Matlab (use `addpath` if needed). Also, note that by specification of more (or less) than `x` evaluation points, re-sampling of the data is easily possible.

## License

The content of this package is released to the public domain. Woltring's original software is available for unrestricted non-commercial use [2].

## References

[1] Woltring, Herman J. (1986). A Fortran package for generalized, cross-validatory spline smoothing and differentiation. *Advances in Engineering Software* 8(2). doi: [10.1016/0141-1195(86)90098-7](https://doi.org/10.1016/0141-1195(86)90098-7).

[2] International Society of Biomechanics. Signal processing software. url: [https\://isbweb.org/software/sigproc.html](https://isbweb.org/software/sigproc.html).  (retrieved in November, 2019).

[3] AT&T Bell Laboratories. A Fortran to C converter. url: [http\://www.netlib.org/f2c/](http://www.netlib.org/f2c/). (retrieved in November, 2019).

[4] Craven, Peter and Grace Wahba (1979). Smoothing noisy data with spline functions. *Numerische Mathematik* 31(4). doi: [10.1007/BF01404567](https://doi.org/10.1007/BF01404567)

[5] Lyche, Tom, Larry L. Schumaker, and Kamy Sephehrnoori (1983). Fortran subroutines for computing smoothing and interpolating natural splines.  *Advances in Engineering Software* 5(1). doi: [10.1016/0141-1195(83)90073-6](https://doi.org/10.1016/0141-1195(83)90073-6)

