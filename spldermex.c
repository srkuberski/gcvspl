/* --------------------------------------------------------------------------
 * MEX wrapper for Woltring's SPLDER function
 *
 * y = spldermex( x, c, m, t, n )
 *
 * INPUT
 * x : independent variables [sorted] (double vector)
 * c : spline coefficients (double vector)
 * m : spline half order [1-4: linear, cubic, quintic, heptic] (double scalar)
 * t : evaluation points (double vector)
 * n : order of derivative (double scalar)
 *
 * OUTPUT
 * y : spline values or derivatives (double vector)
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
	if (nr < 1 || !mxIsDouble( XIN ) || (mxGetM( XIN ) != 1 && mxGetN( XIN ) != 1))
		mexErrMsgTxt( "invalid argument: x (double vector)" );
	if (nr < 2 || !mxIsDouble( CIN ) || (mxGetM( CIN ) != 1 && mxGetN( CIN ) != 1))
		mexErrMsgTxt( "invalid argument: c (double vector)" );
	if (nr < 3 || !mxIsDouble( MIN ) || !mxIsScalar( MIN ))
		mexErrMsgTxt( "invalid argument: m (double scalar)" );
	if (nr < 4 || !mxIsDouble( TIN ) || (mxGetM( TIN ) != 1 && mxGetN( TIN ) != 1))
		mexErrMsgTxt( "invalid argument: t (double vector)" );
	if (nr < 5 || !mxIsDouble( NIN ) || !mxIsScalar( NIN ))
		mexErrMsgTxt( "invalid argument: n (double scalar)" );

		/* prepare buffers */
	long int ider = mxGetScalar( NIN );
	long int m = mxGetScalar( MIN );
	long int n = fmax( mxGetM( XIN ), mxGetN( XIN ) );
	double * t = mxGetPr( TIN );
	double * x = mxGetPr( XIN );
	double * c = mxGetPr( CIN );
	long int l = 0; /* TODO*/
	double * q = (double *) mxCalloc( 2*m, sizeof( double ) );

		/* call Woltring's SPLDER */
	YOUT = mxCreateDoubleMatrix( mxGetM( TIN ), mxGetN( TIN ), 0 );
	double * y = mxGetPr( YOUT );

	for (long int ti = 0; ti < fmax( mxGetM( TIN ), mxGetN( TIN ) ); ++ti)
		y[ti] = splder_( &ider, &m, &n, &t[ti], x, c, &l, q );

		/* release buffers */
	mxFree( q );

}

