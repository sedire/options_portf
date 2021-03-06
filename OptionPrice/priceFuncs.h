#ifndef _PRICE_FUNCS_
#define _PRICE_FUNCS_ 1

#define _USE_MATH_DEFINES

#include "hyperDual.h"
#include "optionPriceTypes.h"
#include "math.h"
#include <fstream>
#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

using namespace System;
using std::ofstream;
using std::max;

P_PRES calcCumDistr( P_PRES x );
HPD<P_PRES, 1> calcCumDistr( HPD<P_PRES, 1> x );

P_PRES priceEuroPut( P_PRES curT, P_PRES expT, P_PRES stockPr,
						P_PRES strikePr, P_PRES r, P_PRES sigma ); 
P_PRES dPdE( P_PRES curT, P_PRES expT, P_PRES stockPr,
						P_PRES strikePr, P_PRES r, P_PRES sigma ); 

HPD<P_PRES, 1> priceEuroPut( P_PRES curT, P_PRES expT, P_PRES stockPr,
						HPD<P_PRES, 1> strikePr, P_PRES r, P_PRES sigma ); 

P_PRES optionPriceBench( P_PRES expT, P_PRES S,
						P_PRES K, P_PRES r, P_PRES sigma, P_PRES div,
						int putCall, int optionType );		//putCall: 0 - Put, 1 - Call
															//optionType: 0 - Euro, 1 - American

template <class PR_NUM>
PR_NUM payoff( PR_NUM x, P_PRES tau, P_PRES k )
{
	return exp( 0.25 * ( k + 1.0 ) * ( k + 1.0 ) * tau ) * max( exp( 0.5 * ( k - 1.0 ) * x ) - exp( 0.5 * ( k + 1.0 ) * x ), 0.0 );	
} 

template <class PR_NUM>
int priceAmerPut2( P_PRES expT, P_PRES Smin, P_PRES Smax
							, P_PRES E, P_PRES r, P_PRES sigma
							, vector<PR_NUM>* P, vector<P_PRES>* Sarr
							/*, vector<PR_NUM>* u, vector<PR_NUM>* g, vector<PR_NUM>* b*/ )
{
	P_PRES k = r / 0.5 / sigma / sigma;			//a param for the transformed equation

	P_PRES finalTau = 0.5 * sigma * sigma * expT;
	P_PRES dt = finalTau / NUM_OF_TIME_STEPS;

	//P_PRES LminRef = fabs( floor( log( LOWEST_PRICE / E ) ) );	//that's what I had initially
	P_PRES LminRef = fabs( floor( log( Smin / LOW_PRICE_FACTOR / E ) ) );
	//P_PRES LmaxRef = log( HIGH_PRICE_FACTOR );	//log( E * HIGH_PRICE_FACTOR / E );		//that's what I had initially
	P_PRES LmaxRef = log( Smax * HIGH_PRICE_FACTOR / E );

	P_PRES dx = ( LmaxRef + LminRef ) / ( NUM_OF_SPACE_STEPS - 1 );	//an alternative

	int Nmin = (int)ceil( ( LminRef - dx / 2.0 ) / dx ) + 1;
	int Nmax = (int)ceil( ( LmaxRef - dx / 2.0 ) / dx ) + 1;
	int N = Nmin + Nmax;
	P_PRES Lmin = ( Nmin - 1 ) * dx + dx / 2.0;
	P_PRES Lmax = ( Nmax - 1 ) * dx + dx / 2.0;
	
	//cout << " == " << LminRef << " " << LmaxRef << " " << dx << endl;
	//cout << " ++ " << Lmin << " " << Lmax << endl;
	//cout << " ==$ " << Nmin << " " << Nmax << " " << N << endl;

	P_PRES alpha = dt / dx / dx;
	//cout << " N " << N << endl;
	//cout << " dt: " << dt << " ; dx: " << dx << " ; alpha: " << alpha << endl;

	vector<PR_NUM> u( N, 0.0 );
	
	if( P != 0 && Sarr != 0/* && u != 0 && g != 0 && b != 0*/ )
	{
		P->resize( N, 0.0 );
		Sarr->resize( N, 0.0 );
		//u->resize( N, 0.0 );
		//g->resize( N - 2, 0.0 );
		//b->resize( N - 2, 0.0 );
	}
	else
	{
		cout << "ERROR: input arrays are zero pointers!\n";
		return 1;
	}

	PR_NUM curX = 0.0;
	for( int i = 0; i < N; ++i )
	{
		curX = -Lmin + i * dx;
		//getElem( curX, 1 ) = 1.0;
		setElem( curX, 1, 1.0 );

		(*Sarr)[i] = E * exp( getRealPart( curX ) );
		u[i] = payoff<PR_NUM>( curX, 0.0, k );
	}
	//cout << " __ " << (*Sarr)[0] << " " << (*Sarr)[1] << endl; 

	//I think I can allocate these outside of this function and re-use them
	vector<PR_NUM> g( N - 2, 0.0 );
	vector<PR_NUM> b( N - 2, 0.0 );

	P_PRES omega = 1.0;
	P_PRES dOmega = 0.05;

	int oldLoops = 10000;

	P_PRES curTau = dt;
	while( curTau <= finalTau )
	{
		//cout << " ---------------\n";
		//cout << " time " << curTau << endl;
		for( int i = 0; i < N - 2; ++i )
		{
			curX = -Lmin + ( i + 1 ) * dx;
			setElem( curX, 1, 1.0 );
			g[i] = payoff<PR_NUM>( curX, curTau, k );
			b[i] = u[i + 1] + alpha / 2.0 * ( u[i + 2] - 2.0 * u[i + 1] + u[i] );
		}

		//fix the boundary conditions
		curX = -Lmin;
		setElem( curX, 1, 1.0 );
		u[0] = payoff<PR_NUM>( curX, curTau, k );
		curX = Lmax;
		setElem( curX, 1, 1.0 );
		u[N - 1] = payoff<PR_NUM>( curX, curTau, k );

		//initial guess for the iterations
		for( int i = 1; i < N - 1; ++i )
		{
			curX = -Lmin + i * dx;
			setElem( curX, 1, 1.0 );

			u[i] = max( u[i], payoff<PR_NUM>( curX, curTau, k ) );
		}

		PR_NUM y = 0.0;

		P_PRES error = 0.0;
		int loops = 0;
		do
		{
			error = 0.0;
			for( int i = 1; i < N - 1; ++i )
			{
				y = 1.0 / ( 1.0 + alpha ) * ( b[i - 1] + 0.5 * alpha * ( u[i - 1] + u[i + 1] ) );
				y = max( u[i] + omega * ( y - u[i] ), g[i - 1] );
				//y = u[i] + omega * ( y - u[i] );
				error += ( getRealPart( u[i] ) - getRealPart( y ) ) * ( getRealPart( u[i] ) - getRealPart( y ) );
				u[i] = y;
			}
			//cout << " loop " << loops << " error " << error << endl;
			++loops;
		} while( error > PSOR_ERROR_EPS * PSOR_ERROR_EPS );

		if( loops > oldLoops )
		{
			dOmega *= -1.0;
		}
		if( omega + dOmega < 2.0 && omega + dOmega > 0.0 )
		{
			omega += dOmega;
		}
		oldLoops = loops;

		//cout << " loops " << loops << " ; omega " << omega << endl;

		curTau += dt;
	}

	for( int i = 0; i < N; ++i )
	{
		curX = -Lmin + i * dx;
		setElem( curX, 1, 1.0 );

		(*Sarr)[i] = E * exp( getRealPart( curX ) );
		(*P)[i] = E * exp( -0.5 * ( k - 1.0 ) * curX - 0.25 * ( k + 1.0 ) * ( k + 1.0 ) * finalTau ) * u[i];
		setElem( (*P)[i], 1, exp( -0.5 * ( k - 1.0 ) * getRealPart( curX ) - 0.25 * ( k + 1.0 ) * ( k + 1.0 ) * finalTau )
							* ( ( 1 + 0.5 * ( k - 1.0 ) ) * getRealPart( u[i] ) - getElem( u[i], 1 ) ) );
	}
	
	//cout << " damping results...\n";
	//ofstream off( "test2.txt" );
	//for( int i = 0; i < N; i+=10 )
	//{
	//	off << (*Sarr)[i] << " " << (*P)[i].real();// << " " << optionPriceBench( expT, (*Sarr)[i], E, r, sigma, 0.0, 0, 1 ); //<< " " << P[i].elems[1] / Sarr[i];
	//	//if( i != 0 && i != N - 1 )
	//	//{
	//	//	off << " " << ( P[i + 1].real() - P[i - 1].real() ) / ( Sarr[i + 1] - Sarr[i - 1] );
	//	//}
	//	off << endl;
	//}
	//off.close();
	//cout << "done!\n";

	return 0;
}

//template <class PR_NUM>
//int priceAmerPut( P_PRES expT, P_PRES Sinit, P_PRES E, P_PRES r, P_PRES sigma,
//					vector<PR_NUM>* P, vector<P_PRES>* Sarr )
//{
//	P_PRES k = r / 0.5 / sigma / sigma;			//a param for the transformed equation
//	//P_PRES dt = sigma * sigma / 2.0 * DT_WANTED;	//tranform the timestep -- had this originally. Probably don't need it
//
//	P_PRES finalTau = 0.5 * sigma * sigma * expT;
//
//	P_PRES dt = finalTau / NUM_OF_TIME_STEPS;
//	//cout << " orig dt is " << dt * 2.0 / sigma / sigma << endl;
//
//	//cout << " num of time steps " << finalTau / dt << endl;
//
//	//P_PRES LminRef = fabs( floor( log( LOWEST_PRICE / E ) ) );	//that's what I had initially
//	P_PRES LminRef = fabs( floor( log( Sinit / LOW_PRICE_FACTOR / E ) ) );
//	//P_PRES LmaxRef = log( HIGH_PRICE_FACTOR );	//log( E * HIGH_PRICE_FACTOR / E );		//that's what I had initially
//	P_PRES LmaxRef = log( Sinit * HIGH_PRICE_FACTOR / E );
//
//	P_PRES dx = ( LmaxRef + LminRef ) / NUM_OF_SPACE_STEPS;	//an alternative
//	//P_PRES dx = log( DS_AT_0 / E + 1 );	//was this
//
//	//cout << " num of spatial nodes " << ( LmaxRef + LminRef ) / dx << endl;
//
//	//I do this to avoid a mesh node at x == 0. There may be a singularity at x == 0 and I don't want to deal with it.
//	int Nmin = (int)ceil( ( LminRef - dx / 2.0 ) / dx ) + 1;
//	int Nmax = (int)ceil( ( LmaxRef - dx / 2.0 ) / dx ) + 1;
//	int N = Nmin + Nmax;
//	P_PRES Lmin = ( Nmin - 1 ) * dx + dx / 2.0;
//	P_PRES Lmax = ( Nmax - 1 ) * dx + dx / 2.0;
//
//	//cout << " ---- " << E * exp( LmaxRef ) << " " << E * exp( -LminRef ) << endl;
//
//	//cout << " == " << LminRef << " " << LmaxRef << " " << dx << endl;
//	//cout << " ++ " << Lmin << " " << Lmax << endl;
//	cout << " ==$ " << Nmin << " " << Nmax << " " << N << endl;
//
//	P_PRES alpha = dt / dx / dx;
//	//cout << " N " << N << endl;
//	//cout << " dt: " << dt << " ; dx: " << dx << " ; alpha: " << alpha << endl;
//
//	vector<PR_NUM> u( N, 0.0 );
//	
//	if( P != 0 && Sarr != 0 )
//	{
//		P->resize( N, 0 );
//		Sarr->resize( N, 0 );
//	}
//	else
//	{
//		cout << "ERROR: input arrays are zero pointers!\n";
//		return 1;
//	}
//
//	PR_NUM curX = 0.0;
//	for( int i = 0; i < N; ++i )
//	{
//		curX = -Lmin + i * dx;
//		curX.elems[1] = 1.0;
//
//		(*Sarr)[i] = E * exp( curX.real() );
//		u[i] = payoff<PR_NUM>( curX, 0.0, k );
//	}
//
//	vector<PR_NUM> g( N - 2, 0.0 );
//	vector<PR_NUM> b( N - 2, 0.0 );
//
//	P_PRES omega = 1.0;
//	P_PRES dOmega = 0.05;
//
//	int oldLoops = 10000;
//
//	P_PRES curTau = dt;
//	while( curTau <= finalTau )
//	{
//		//cout << " ---------------\n";
//		//cout << " time " << curTau << endl;
//		for( int i = 0; i < N - 2; ++i )
//		{
//			curX = -Lmin + ( i + 1 ) * dx;
//			curX.elems[1] = 1.0;
//			g[i] = payoff<PR_NUM>( curX, curTau, k );
//			b[i] = u[i + 1] + alpha / 2.0 * ( u[i + 2] - 2.0 * u[i + 1] + u[i] );
//		}
//
//		//fix the boundary conditions
//		curX = -Lmin;
//		curX.elems[1] = 1.0;
//		u[0] = payoff<PR_NUM>( curX, curTau, k );
//		curX = Lmax;
//		curX.elems[1] = 1.0;
//		u[N - 1] = payoff<PR_NUM>( curX, curTau, k );
//
//		//initial guess for the iterations
//		for( int i = 1; i < N - 1; ++i )
//		{
//			curX = -Lmin + i * dx;
//			curX.elems[1] = 1.0;
//
//			u[i] = max( u[i], payoff<PR_NUM>( curX, curTau, k ) );
//		}
//
//		PR_NUM y = 0.0;
//
//		P_PRES error = 0.0;
//		int loops = 0;
//		do
//		{
//			error = 0.0;
//			for( int i = 1; i < N - 1; ++i )
//			{
//				y = 1.0 / ( 1.0 + alpha ) * ( b[i - 1] + 0.5 * alpha * ( u[i - 1] + u[i + 1] ) );
//				y = max( u[i] + omega * ( y - u[i] ), g[i - 1] );
//				//y = u[i] + omega * ( y - u[i] );
//				error += ( u[i].real() - y.real() ) * ( u[i].real() - y.real() );
//				u[i] = y;
//			}
//			//cout << " loop " << loops << " error " << error << endl;
//			++loops;
//		} while( error > PSOR_ERROR_EPS * PSOR_ERROR_EPS );
//
//		if( loops > oldLoops )
//		{
//			dOmega *= -1.0;
//		}
//		if( omega + dOmega < 2.0 && omega + dOmega > 0.0 )
//		{
//			omega += dOmega;
//		}
//		oldLoops = loops;
//
//		//cout << " loops " << loops << " ; omega " << omega << endl;
//
//		curTau += dt;
//	}
//
//	for( int i = 0; i < N; ++i )
//	{
//		curX = -Lmin + i * dx;
//		curX.elems[1] = 1.0;
//
//		(*Sarr)[i] = E * exp( curX.real() );
//		(*P)[i] = E * exp( -0.5 * ( k - 1.0 ) * curX - 0.25 * ( k + 1.0 ) * ( k + 1.0 ) * finalTau ) * u[i];
//		(*P)[i].elems[1] = exp( -0.5 * ( k - 1.0 ) * curX.real() - 0.25 * ( k + 1.0 ) * ( k + 1.0 ) * finalTau )
//							* ( ( 1 + 0.5 * ( k - 1.0 ) ) * u[i].real() - u[i].elems[1] );
//	}
//	
//	//cout << " damping results...\n";
//	//ofstream off( "test2.txt" );
//	//for( int i = 0; i < N; i+=10 )
//	//{
//	//	off << (*Sarr)[i] << " " << (*P)[i].real();// << " " << optionPriceBench( expT, (*Sarr)[i], E, r, sigma, 0.0, 0, 1 ); //<< " " << P[i].elems[1] / Sarr[i];
//	//	//if( i != 0 && i != N - 1 )
//	//	//{
//	//	//	off << " " << ( P[i + 1].real() - P[i - 1].real() ) / ( Sarr[i + 1] - Sarr[i - 1] );
//	//	//}
//	//	off << endl;
//	//}
//	//cout << "done!\n";
//
//	return 0;
//}

template <class PR_NUM>
PR_NUM optPriceAt( P_PRES spot, const vector<PR_NUM>& P, const vector<P_PRES>& Sarr )
{
	if( Sarr.size() < 3 || Sarr.size() != P.size() )
	{
		cout << " something is wrong in optPriceAt\n";
		return -1;
	}

	int prevInd = 0;
	int nextInd = 1;
	PR_NUM ret = 0.0;

	for( int i = 0; i < Sarr.size() - 1; ++i )
	{
		if( nextInd >= Sarr.size() )
		{
			cout << " something is wrong in optPriceAt 2\n";
			return -1;
		}
		if( Sarr[nextInd] < spot )
		{
			++prevInd;
			++nextInd;
			continue;
		}
		else if( Sarr[nextInd] >= spot && Sarr[prevInd] <= spot )
		{
			ret = ( spot - Sarr[prevInd] ) * ( P[nextInd] - P[prevInd] ) / ( Sarr[nextInd] - Sarr[prevInd] ) + P[prevInd];
			break;
		}
		else
		{
			cout << " something is wrong in optPriceAt 3\n";
			cout << prevInd << " " << nextInd << " " << Sarr.size() << " " << Sarr[nextInd] << " " << spot << endl;
			return -1;
		}
	}

	return ret;
}

template <class PR_NUM>
int reinterpolate( P_PRES Smax,	const vector<PR_NUM>& P, const vector<P_PRES>& Sarr,
					vector<PR_NUM>* resArr )
{
	if( Smax >= Sarr[Sarr.size() - 1] )
	{
		cout << "ERROR in interpolation: Smax is too big\n";
		return 1;
	}
	P_PRES Smin = 0;
	size_t resArrSize = resArr->size();
	P_PRES dS = Smax / ( resArrSize - 1 );

	(*resArr)[0] = ( Sarr[0] * ( P[0] - P[1] ) / ( Sarr[1] - Sarr[0] ) ) + P[0];

	for( int i = 1; i < resArrSize; ++i )
	{
		P_PRES curS = i * dS;
		int leftInd = Sarr.size() - 2;

		while( leftInd > 0 && !( Sarr[leftInd] < curS && Sarr[leftInd + 1] >= curS ) )
		{
			--leftInd;
		}
		P_PRES dist = Sarr[leftInd + 1] - Sarr[leftInd];
		(*resArr)[i] = P[leftInd] * ( Sarr[leftInd + 1] - curS ) / dist + P[leftInd + 1] * ( curS - Sarr[leftInd] ) / dist;
	}

	ofstream off( "test3.txt" );
	cout << " damping results...\n";
	for( int i = 0; i < resArrSize; ++i )
	{
		P_PRES curS = i * dS;
		off << curS << " " << (*resArr)[i].real() << " " << (*resArr)[i].elems[1] << endl;
	}
	off.close();

	return 0;
}
#endif