/* --------------------------------------------------------------------------
 * MEX wrapper for Woltring's GCVSPL function
 *
 * [c, wk] = gcvsplmex( x, y, m, v, w )
 *
 * INPUT
 * x : independent variables [sorted] (double column [N, 1])
 * y : data to be smoothed (double matrix [N, K])
 * m : spline half order [1-4: linear, cubic, quintic, heptic] (double scalar)
 * v : prior given variance [negative for cross-validation] (double scalar)
 * w : weight factors (double column [N, 1])
 *
 * OUTPUT
 * c : spline coefficients (double matrix [N, K])
 * wk : GCVSPL work vector (numeric row [1, 6])
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
#define VIN pr[3]
#define WIN pr[4]

#define COUT pl[0]
#define WKOUT pl[1]

void
mexFunction( int nl, mxArray ** pl, int nr, const mxArray ** pr )
{

		/* safeguard */
	if (nr < 1 || !mxIsDouble( XIN ) || mxGetN( XIN ) != 1)
		mexErrMsgTxt( "invalid argument: x (double row [N, 1])" );
	long int n = mxGetM( XIN );
	double * x = mxGetPr( XIN );

	if (nr < 2 || !mxIsDouble( YIN ) || mxGetM( YIN ) != n)
		mexErrMsgTxt( "invalid argument: y (double matrix [N, K])" );
	long int k = mxGetN( YIN );
	double * y = mxGetPr( YIN );

	if (nr < 3 || !mxIsDouble( MIN ) || !mxIsScalar( MIN ))
		mexErrMsgTxt( "invalid argument: m (double scalar)" );
	long int m = mxGetScalar( MIN );

	if (nr < 4 || !mxIsDouble( VIN ) || !mxIsScalar( VIN ))
		mexErrMsgTxt( "invalid argument: v (double scalar)" );
	double val = mxGetScalar( VIN );

	if (nr < 5 || !mxIsDouble( WIN ) || mxGetN( WIN ) != 1 || mxGetM( WIN ) != n)
		mexErrMsgTxt( "invalid argument: w (double row [N, 1])" );
	double * wx = mxGetPr( WIN );

		/* prepare buffers */
	double * wy = (double *) mxCalloc( k, sizeof( double ) ); /* uniform weights across dimensions */
	for (long int wi = 0; wi < k; ++wi )
		wy[wi] = 1;

	long int md = 3; /* default mode: predicted variance, alternative: generalized cross-validation */
	if (val < 0)
		md = 2;

	COUT = mxCreateDoubleMatrix( n, k, 0 );
	double * c = mxGetPr( COUT );

	double * wk = (double *) mxCalloc( 6*(n*m+1)+n, sizeof( double ) );
	long int ier;

		/* call Woltring's GCVSPL */
	gcvspl_( x, y, &n, wx, wy, &m, &n, &k, &md, &val, c, &n, wk, &ier );

	WKOUT = mxCreateDoubleMatrix( 1, 6, 0 );
	double * wkp = mxGetPr( WKOUT );
	for (int wi = 0; wi < 6; wi++)
		wkp[wi] = wk[wi];

		/* release buffers */
	mxFree( wy );
	mxFree( wk );

}

