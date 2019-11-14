/* --------------------------------------------------------------------------
 * MEX wrapper for Woltring's GCVSPL function
 *
 * c = gcvsplmex( x, y, m )
 *
 * INPUT
 * x : independent variables [sorted] (double vector)
 * y : data to be smoothed (double vector)
 * m : spline half order [1-4: linear, cubic, quintic, heptic] (double scalar)
 *
 * OUTPUT
 * c : spline coefficients (double vector)
 *
 * REMARKS
 * Compile with 'mex gcvsplmex.c gcvspl.c'.
 * ------------------------------------------------------------------------ */

#include <math.h>
#include "mex.h"
#include "gcvspl.h"

#define XIN pr[0]
#define YIN pr[1]
#define MIN pr[2]

#define COUT pl[0]

void
mexFunction( int nl, mxArray ** pl, int nr, const mxArray ** pr )
{

		/* safeguard */
	if (nr < 1 || !mxIsDouble( XIN ) || (mxGetM( XIN ) != 1 && mxGetN( XIN ) != 1))
		mexErrMsgTxt( "invalid argument: x (double vector)" );
	if (nr < 2 || !mxIsDouble( YIN ) || (mxGetM( YIN ) != 1 && mxGetN( YIN ) != 1))
		mexErrMsgTxt( "invalid argument: y (double vector)" );
	if (nr < 3 || !mxIsDouble( MIN ) || !mxIsScalar( MIN ))
		mexErrMsgTxt( "invalid argument: m (double scalar)" );

		/* prepare buffers */
	double * x = mxGetPr( XIN );
	double * y = mxGetPr( YIN );
	long int n = fmax( mxGetM( XIN ), mxGetN( XIN ) );
	
	double * wx = (double *) mxCalloc( n, sizeof( double ) ); /* uniform weights */
	double * wy = (double *) mxCalloc( n, sizeof( double ) );
	for (long int wi = 0; wi < n; ++wi )
	{
		wx[wi] = 1;
		wy[wi] = 1;
	}

	long int m = mxGetScalar( MIN );
	long int k = 1; /* single data set */
	long int md = 2; /* general cross-validation */
	double val = 0; /* N/A */

	COUT = mxCreateDoubleMatrix( mxGetM( XIN ), mxGetN( XIN ), 0 );
	double * c = mxGetPr( COUT );

	double * wk = (double *) mxCalloc( 6*(n*m+1)+n, sizeof( double ) );
	long int ier;

		/* call Woltring's GCVSPL */
	gcvspl_( x, y, &n, wx, wy, &m, &n, &k, &md, &val, c, &n, wk, &ier );

		/* release buffers */
	mxFree( wx );
	mxFree( wy );
	mxFree( wk );

}

