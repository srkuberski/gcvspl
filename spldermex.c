/* --------------------------------------------------------------------------
 * MEX wrapper for Woltring's SPLDER function
 *
 * y = spldermex( x, c, m, t, n )
 *
 * INPUT
 * x : independent variables [sorted] (double column [N, 1])
 * c : spline coefficients (double column [N, 1])
 * m : spline half order [1-4: linear, cubic, quintic, heptic] (double scalar)
 * t : evaluation points (double column [T, 1])
 * n : order of derivative (double scalar)
 *
 * OUTPUT
 * y : spline values or derivatives (double column [T, 1])
 *
 * REMARKS
 * Compile with 'mex spldermex.c gcvspl.c'.
 * ------------------------------------------------------------------------ */

#include <math.h>
#include "mex.h"
#include "gcvspl.h"

#define XIN pr[0]
#define CIN pr[1]
#define MIN pr[2]
#define TIN pr[3]
#define NIN pr[4]

#define YOUT pl[0]

void
mexFunction( int nl, mxArray ** pl, int nr, const mxArray ** pr )
{

		/* safeguard */
	if (nr < 1 || !mxIsDouble( XIN ) || mxGetN( XIN ) != 1)
		mexErrMsgTxt( "invalid argument: x (double column [nx, 1])" );
	long int nx = mxGetM( XIN );
	double * x = mxGetPr( XIN );

	if (nr < 2 || !mxIsDouble( CIN ) || mxGetM( CIN ) != nx || mxGetN( CIN ) != 1)
		mexErrMsgTxt( "invalid argument: c (double column [nx, 1])" );
	double * c = mxGetPr( CIN );

	if (nr < 3 || !mxIsDouble( MIN ) || !mxIsScalar( MIN ))
		mexErrMsgTxt( "invalid argument: m (double scalar)" );
	long int m = mxGetScalar( MIN );

	if (nr < 4 || !mxIsDouble( TIN ) || mxGetN( TIN ) != 1)
		mexErrMsgTxt( "invalid argument: t (double column [nx, 1])" );
	long int nt = mxGetM( TIN );
	double * t = mxGetPr( TIN );

	if (nr < 5 || !mxIsDouble( NIN ) || !mxIsScalar( NIN ))
		mexErrMsgTxt( "invalid argument: nx (double scalar)" );
	long int ider = mxGetScalar( NIN );

		/* prepare buffers */
	long int l = nx;
	double * q = (double *) mxCalloc( 2*m, sizeof( double ) );

		/* call Woltring's SPLDER */
	YOUT = mxCreateDoubleMatrix( nt, 1, 0 );
	double * y = mxGetPr( YOUT );

	for (long int ti = 0; ti < nt; ++ti)
		y[ti] = splder_( &ider, &m, &nx, &t[ti], x, c, &l, q );

		/* release buffers */
	mxFree( q );

}

