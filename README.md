# GCV spline smoothing in Matlab

This is a Matlab interface to Woltring's classic generalized, cross-validatory (GCV) spline smoothing and differentiation code [1]. Woltring's original Fortran 77 code [2] was converted to C using the f2c converter [3]. Two Matlab MEX wrappers have been implemented to make the code available to a wider range of users.

- [Theoretical minimum](#theoretical-minimum)
- [Matlab interface](#matlab-interface)
- [Matlab example](#matlab-example)
- [License](#license)
- [References](#references)

## Theoretical minimum

- let $t_i$ ( $i=1,\ldots,n\ge2m$ ) be a set of strictly increasing (not necessarily equidistant) abscissa values with corresponding ordinates $x_i$ and positive weight factors $w_i$
- a natural spline $s_p(t)$ is a piecewise polynomial function on knot positions $t_i$ which minimizes the criterion function $C_p$ for suitably selected regularization parameter $p\ge0$

$$
C_p=\sum_{i=1}^nw_i[x_i-s_p(t_i)]^2+p\int_{t=-\infty}^{+\infty}|s_p^{(m)}(t)|^2dt
$$

- through the regularization parameter $p$, a trade-off can be effectuated between the smootheness of the spline (its $m$th derivative) and the goodness-of-fit to the given data
- for $p=0$, an exactly interpolating spline is obtained; in the limiting case $p\rightarrow\infty$, the spline becomes an $m$ th order (i.e., an $(m-1)$ th degree) polynomial which fits the data in a weighted least-squares sense
- Woltring[1] implemented ways to estimate the regularization parameter $p$ by generalized cross-validation (GCV) or from a given, predicted mean squared error (MSE)
- the spline polynomials have order $\le2m$ between the knots and order $m$ outside the knot range
- the polynomials are continuous at the knots up to and including the $(2m-2)$th derivative with vanishing $m$th and higher derivatives at the terminal knots $t_1$ and $t_n$

|                                         | linear | cubic | quintic | heptic |
|----------------------------------------:|:------:|:-----:|:-------:|:------:|
|                          Half order $m$ |    1   |   2   |    3    |    4   |
|                   Polynomial order $2m$ |    2   |   4   |    6    |    8   |
|                Polynomial degree $2m-1$ |    1   |   3   |    5    |    7   |
| Number of continuous derivatives $2m-2$ |    0   |   2   |    4    |    6   |
|            Minimum number of knots $2m$ |    2   |   4   |    6    |    8   |
|  Terminal knots w/o full support $2m-2$ |    0   |   2   |    4    |    6   |

*Table*: Key figures of frequently used spline representations

## Matlab interface

For the construction and evaluation of splines, Woltring's original Fortran code provided two separate functions, `gcvspl` and `splder`. Access to these functions are made available here by the use of two corresponding MEX wrappers. In order to use these wrappers with your Matlab installation, you need to first (and only once) compile them on your computer. The compilation requires the presence of a C/C++ compiler software on your computer (e.g., the GNU Compiler Collection gcc). If you do not know if this is the case or not, please refer to the Matlab documentation on how to install, setup, and test a proper MEX building environment before. Once a C compiler is available on your computer, start Matlab, change to the path where the files of this package are located and run:

~~~matlab
mex gcvsplmex.c gcvspl.c
mex spldermex.c splder.c
~~~

Alternatively, you can use the `Makefile` provided here. The resulting binaries `gcvsplmex.mexa64` and `spldermex.mex64` (or similar on non-Unix platforms) are Matlab executables which can be used as is. However, the package also provides (yet another and recommend) set of wrappers with additional error handling. Below you will find short documentation on how to use these.

### Spline construction with `gcvspl`

The `gcvspl` function computes a natural B-spline using the generalized cross-validation and mean-squared prediction error criteria of Craven & Wahba [4]. The model assumes uncorrelated, additive noise and essentially smooth, underlying functions. The independent coordinates may be spaced non-equidistantly.

~~~matlab
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
% wk : GCVSPL work vector (numeric row [1, 6])
~~~
*Code*: Type `help gcvspl` to see this in your Matlab installation

### Spline evaluation with `splder`

The `splder` function computes the values of a B-spline or its derivatives at selected evaluation points. Woltring's underlying code is an adoption of Lyche et al. [5].

~~~matlab
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
~~~
*Code*: Type `help splder` to see this in your Matlab installation

### Zero crossings with `splzer`

The `splzer` function computes the location (and sign) of zero crossings of a B-spline or its derivatives.

~~~matlab
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
~~~
*Code*: Type `help splzer` to see this in your Matlab installation

## Matlab example

As an example of usage, consider the following lines of code, which compute the (fifth-order) quintic spline values and their first two derivatives (velocity and acceleration) for a given set of data points `(x,y)`:

~~~matlab
c = gcvspl( x, y, 3 ); % compute GCV spline coefficients
y0 = splder( x, c, 3, x, 0 ); % compute spline values
y1 = splder( x, c, 3, x, 1 ); % compute first derivatives (velocity)
y2 = splder( x, c, 3, x, 2 ); % compute second derivatives (acceleration)
z = splzer( x, c, 3, 1 ); % compute zeros in the first derivative (velocity)
~~~
*Code*: A simple example of using the splines interface

Note that by specification of more (or less) than `x` evaluation points in the fourth argument to `splder`, re-sampling of the original data is possible.

## License

For archiving and documentation purposes, the original Fortran 77 source file (`gcvspl.f`) is duplicated here, as well as the intermediary C code files (`gcvspl.c`, `gcvspl.h`). Please adhere to their original copyright (unrestricted non-commercial use [2]). Everything else provided in this repository is released to the public domain.

## References

[1] Woltring, Herman J. (1986). A Fortran package for generalized, cross-validatory spline smoothing and differentiation. *Advances in Engineering Software* 8(2). doi: [10.1016/0141-1195(86)90098-7](https://doi.org/10.1016/0141-1195(86)90098-7).

[2] International Society of Biomechanics. Signal processing software. url: [https\://isbweb.org/software/sigproc.html](https://isbweb.org/software/sigproc.html).  (retrieved in November, 2019).

[3] AT&T Bell Laboratories. A Fortran to C converter. url: [http\://www.netlib.org/f2c/](http://www.netlib.org/f2c/). (retrieved in November, 2019).

[4] Craven, Peter and Grace Wahba (1979). Smoothing noisy data with spline functions. *Numerische Mathematik* 31(4). doi: [10.1007/BF01404567](https://doi.org/10.1007/BF01404567)

[5] Lyche, Tom, Larry L. Schumaker, and Kamy Sephehrnoori (1983). Fortran subroutines for computing smoothing and interpolating natural splines.  *Advances in Engineering Software* 5(1). doi: [10.1016/0141-1195(83)90073-6](https://doi.org/10.1016/0141-1195(83)90073-6)

