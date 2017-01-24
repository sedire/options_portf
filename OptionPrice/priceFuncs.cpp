#include "priceFuncs.h"

P_PRES priceEuroPut( P_PRES curT, P_PRES expT, P_PRES stockPr,
						P_PRES strikePr, P_PRES r, P_PRES sigma )
{
	P_PRES d1 = ( log( stockPr / strikePr ) + ( r + 0.5 * sigma * sigma ) * ( expT - curT ) ) / sigma / sqrt( expT - curT );
	P_PRES d2 = ( log( stockPr / strikePr ) + ( r - 0.5 * sigma * sigma ) * ( expT - curT ) ) / sigma / sqrt( expT - curT );

	return strikePr * exp( -r * ( expT - curT ) ) * calcCumDistr( -d2 ) - stockPr * calcCumDistr( -d1 );
}

HPD<P_PRES, 1> priceEuroPut( P_PRES curT, P_PRES expT, P_PRES stockPr,
						HPD<P_PRES, 1> strikePr, P_PRES r, P_PRES sigma )
{
	HPD<P_PRES, 1> d1 = ( log( stockPr / strikePr ) + ( r + 0.5 * sigma * sigma ) * ( expT - curT ) ) / sigma / sqrt( expT - curT );
	HPD<P_PRES, 1> d2 = ( log( stockPr / strikePr ) + ( r - 0.5 * sigma * sigma ) * ( expT - curT ) ) / sigma / sqrt( expT - curT );

	return strikePr * exp( -r * ( expT - curT ) ) * calcCumDistr( -d2 ) - stockPr * calcCumDistr( -d1 );
}

P_PRES dPdE( P_PRES curT, P_PRES expT, P_PRES stockPr,
						P_PRES strikePr, P_PRES r, P_PRES sigma )
{
	P_PRES d1 = ( log( stockPr / strikePr ) + ( r + 0.5 * sigma * sigma ) * ( expT - curT ) ) / sigma / sqrt( expT - curT );
	P_PRES d2 = ( log( stockPr / strikePr ) + ( r - 0.5 * sigma * sigma ) * ( expT - curT ) ) / sigma / sqrt( expT - curT );

	P_PRES dE = -1.0 / ( sigma * sqrt( expT - curT ) ) / strikePr;

	return exp( -r * ( expT - curT ) ) * calcCumDistr( -d2 ) + strikePr * exp( -r * ( expT - curT ) ) * ( -1.0 / sqrt( 2.0 * M_PI ) )
			* exp( -0.5 * d2 * d2 ) * dE + stockPr * ( -1.0 / sqrt( 2.0 * M_PI ) ) * exp( -0.5 * d1 * d1 ) * dE;
}

P_PRES calcCumDistr( P_PRES x )
{
	static const P_PRES RT2PI = sqrt( 4.0 * acos( 0.0 ) );

	static const P_PRES SPLIT = 7.07106781186547;

	static const P_PRES N0 = 220.206867912376;
	static const P_PRES N1 = 221.213596169931;
	static const P_PRES N2 = 112.079291497871;
	static const P_PRES N3 = 33.912866078383;
	static const P_PRES N4 = 6.37396220353165;
	static const P_PRES N5 = 0.700383064443688;
	static const P_PRES N6 = 3.52624965998911e-02;
	static const P_PRES M0 = 440.413735824752;
	static const P_PRES M1 = 793.826512519948;
	static const P_PRES M2 = 637.333633378831;
	static const P_PRES M3 = 296.564248779674;
	static const P_PRES M4 = 86.7807322029461;
	static const P_PRES M5 = 16.064177579207;
	static const P_PRES M6 = 1.75566716318264;
	static const P_PRES M7 = 8.83883476483184e-02;

	const P_PRES z = fabs( x );
	P_PRES c = 0.0;

	if( z <= 37.0 )
	{
		const P_PRES e = exp( -z * z / 2.0 );
		if( z < SPLIT )
		{
			const P_PRES n = ( ( ( ( ( N6 * z + N5 ) * z + N4 ) * z + N3 ) * z + N2 ) * z + N1 ) * z + N0;
			const P_PRES d = ( ( ( ( ( ( M7 * z + M6 ) * z + M5 ) * z + M4 ) * z + M3 ) * z + M2 ) * z + M1 ) * z + M0;
			c = e * n / d;
		}
		else
		{
			const P_PRES f = z + 1.0 / ( z + 2.0 / ( z + 3.0 / ( z + 4.0 / ( z + 13.0 / 20.0 ) ) ) );
			c = e / ( RT2PI * f );
		}
	}
	return x <= 0.0 ? c : 1 - c;
}

HPD<P_PRES, 1> calcCumDistr( HPD<P_PRES, 1> x )
{
	static const P_PRES RT2PI = sqrt( 4.0 * acos( 0.0 ) );

	static const P_PRES SPLIT = 7.07106781186547;

	static const P_PRES N0 = 220.206867912376;
	static const P_PRES N1 = 221.213596169931;
	static const P_PRES N2 = 112.079291497871;
	static const P_PRES N3 = 33.912866078383;
	static const P_PRES N4 = 6.37396220353165;
	static const P_PRES N5 = 0.700383064443688;
	static const P_PRES N6 = 3.52624965998911e-02;
	static const P_PRES M0 = 440.413735824752;
	static const P_PRES M1 = 793.826512519948;
	static const P_PRES M2 = 637.333633378831;
	static const P_PRES M3 = 296.564248779674;
	static const P_PRES M4 = 86.7807322029461;
	static const P_PRES M5 = 16.064177579207;
	static const P_PRES M6 = 1.75566716318264;
	static const P_PRES M7 = 8.83883476483184e-02;

	const HPD<P_PRES, 1> z = fabs( x );
	HPD<P_PRES, 1> c = 0.0;

	if( z <= 37.0 )
	{
		const HPD<P_PRES, 1> e = exp( -z * z / 2.0 );
		if( z < SPLIT )
		{
			const HPD<P_PRES, 1> n = ( ( ( ( ( N6 * z + N5 ) * z + N4 ) * z + N3 ) * z + N2 ) * z + N1 ) * z + N0;
			const HPD<P_PRES, 1> d = ( ( ( ( ( ( M7 * z + M6 ) * z + M5 ) * z + M4 ) * z + M3 ) * z + M2 ) * z + M1 ) * z + M0;
			c = e * n / d;
		}
		else
		{
			const HPD<P_PRES, 1> f = z + 1.0 / ( z + 2.0 / ( z + 3.0 / ( z + 4.0 / ( z + 13.0 / 20.0 ) ) ) );
			c = e / ( RT2PI * f );
		}
	}
	return x <= 0.0 ? c : 1.0 - c;
}

P_PRES optionPriceBench( P_PRES expT, P_PRES Spot,
						P_PRES K, P_PRES r, P_PRES sigma, P_PRES div,
						int putCall, int optionType )
{
	using namespace ublas;

	static const int Put = 0;
	static const int Call = 1;
	static const int American = 1;

	//PDE features
	static const P_PRES N = 500;					// Number of time steps
	static const P_PRES M = 500;					// Number of stock price steps
	static const P_PRES dt = expT / N;					// Time increment
	static const P_PRES mu = r - div - ( sigma * sigma ) / 2;	// Drift for stock process
	static const P_PRES dx = sigma * sqrt( 3 * dt );		// Increment for stock price

	//Probabilities
	static const P_PRES pu = -0.25 * dt * ( (sigma * sigma) / ( dx * dx ) + mu / dx );      //Up probability
	static const P_PRES pm = 1.0 + dt * ( sigma * sigma ) / 2.0 / ( dx * dx ) + r * dt / 2.0;      //Middle probability
	static const P_PRES pd = -0.25 * dt * ( ( sigma * sigma ) / ( dx * dx ) - mu / dx );      //Down probability

	//Initialize stock price and option values.
	ublas::vector<P_PRES> S( 2 * M + 1, 0 );
	matrix<P_PRES> V( 2 * M + 1, N + 1, 0 );

	// Indices for stock price step
	ublas::vector<P_PRES> J( 2 * M + 1 );
	for ( int i = M; i >= -M; --i )
	{
		J( M - i ) = i;
	}

	// Stock price at maturity
	for ( unsigned i = 0; i < S.size(); ++i )
	{
		S(i) = Spot * exp( J( i ) * dx );
	}

	// Option price at maturity
	for ( unsigned i = 0; i < V.size1(); i++ )
	{
		V( i, V.size2() - 1 ) = ( ( putCall == Put ) ? std::max( K - S( i ), 0.0 )
													: std::max( S( i ) - K, 0.0 ) ); //populate last column of V
	}	

	//Initialize required matrices
	matrix<P_PRES> pmp( 2 * M + 1, N + 1, 0 );
	matrix<P_PRES> pp( 2 * M + 1, N + 1, 0 );
	matrix<P_PRES> C( 2 * M + 1, N + 1, 0 );
	
	//Create the P' matrix
	for( unsigned j = N; j > 0; --j )
	{
		pmp( 2 * M - 1, j ) = pd + pm;
		for( unsigned i = 2 * M - 2; i > 0; --i )
		{
			pmp( i, j ) = pm - pd / ( pmp( i + 1, j ) ) * pu;
		}
	}


	//Upper and lower bounds for the American put
	P_PRES lambda_L = ( putCall == Put ) ? std::max( K - S( 2 * M ), 0.0 ) : 0.0;
	P_PRES lambda_U = ( putCall == Call ) ? std::max( S( 0 ) - K, 0.0 ) : 0.0;

	// Work backwards and obtain matrix of values V
	for( unsigned j = N; j > 0; --j )
	{
		pp( 2 * M - 1, j ) = -pu * V( 2 * M - 2, j ) - ( pm - 2.0 ) * V( 2 * M - 1, j ) - pd * V( 2 * M, j ) + pd * lambda_L;

		for ( unsigned i = 2 * M - 2; i > 0; --i )
		{
			pp(i,j) = -pu * V(i - 1, j) - (pm - 2)*V(i, j) - pd*V(i + 1, j) - pd / pmp(i + 1, j) * pp(i + 1, j);
		}

		for (unsigned i = 0; i < 2 * M + 1; ++i)
		{
			if (i == 0)
				C(i, j-1) = (pp(i + 1, j) + pmp(i + 1, j) * lambda_U) / (pmp(i + 1, j) + pu);
			else if (i < 2 * M)
				C(i, j - 1) = (pp(i, j) - pu*C(i - 1, j - 1)) / pmp(i, j);
			else
				C(i, j - 1) = C(i - 1, j - 1) - lambda_L;

			if( optionType == American )
			{
				V(i, j - 1) = std::max( ( putCall == Put ) ? K - S(i) : S(i) - K, C(i, j - 1));
			}
			else
			{
				V( i, j - 1 ) = C( i, j - 1 );
			}
		}


	}

	return V( M, 0 );
}