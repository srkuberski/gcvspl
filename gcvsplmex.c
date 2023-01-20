/* --------------------------------------------------------------------------
 * MEX wrapper for Woltring's GCVSPL function
 *
 * [c, wk] = gcvsplmex( x, y, m, v, w )
 *
 * INPUT
 * x : independent variables (double column [N, 1])
 * y : data to be smoothed (double matrix [N, K])
 * m : spline half order (double scalar)
 * v : prior given variance (double scalar)
 * w : weight factors (double column [N, 1])
 *
 * OUTPUT
 * c : spline coefficients (double matrix [N, K])
 * wk : internal work vector (double row [1, 6])
 *
 * REMARKS
 * - compile with 'mex -v gcvsplmex.c gcvspl.c'.
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
	long int N = mxGetM( XIN );
	double * x = mxGetPr( XIN );

	if (nr < 2 || !mxIsDouble( YIN ) || mxGetM( YIN ) != N)
		mexErrMsgTxt( "invalid argument: y (double matrix [N, K])" );
	long int K = mxGetN( YIN );
	double * y = mxGetPr( YIN );

	if (nr < 3 || !mxIsDouble( MIN ) || !mxIsScalar( MIN ))
		mexErrMsgTxt( "invalid argument: m (double scalar)" );
	long int m = mxGetScalar( MIN );

	if (nr < 4 || !mxIsDouble( VIN ) || !mxIsScalar( VIN ))
		mexErrMsgTxt( "invalid argument: v (double scalar)" );
	double val = mxGetScalar( VIN );

	if (nr < 5 || !mxIsDouble( WIN ) || mxGetN( WIN ) != 1 || mxGetM( WIN ) != N)
		mexErrMsgTxt( "invalid argument: w (double row [N, 1])" );
	double * wx = mxGetPr( WIN );

		/* prepare buffers */
	double * wy = (double *) mxCalloc( K, sizeof( double ) ); /* uniform weights across dimensions */
	for (long int wi = 0; wi < K; ++wi )
		wy[wi] = 1;

	long int md = 3; /* default mode: predicted variance, alternative: generalized cross-validation */
	if (val < 0)
		md = 2;

	COUT = mxCreateDoubleMatrix( N, K, 0 );
	double * c = mxGetPr( COUT );

	double * wk = (double *) mxCalloc( 6*(N*m+1)+N, sizeof( double ) );
	long int ier;

		/* call Woltring's GCVSPL */
	gcvspl_( x, y, &N, wx, wy, &m, &N, &K, &md, &val, c, &N, wk, &ier );

	WKOUT = mxCreateDoubleMatrix( 1, 6, 0 );
	double * wkp = mxGetPr( WKOUT );
	for (int wi = 0; wi < 6; wi++)
		wkp[wi] = wk[wi];

		/* release buffers */
	mxFree( wy );
	mxFree( wk );

}

