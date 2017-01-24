#ifndef _MAIN_FUNC_
#define _MAIN_FUNC_ 1

#pragma comment( lib, "C:/Users/dchernikov/Documents/nlopt-2.4.2-dll64/libnlopt-0.lib" )

#include <iostream>
#include <iomanip>
#include "optionPriceTypes.h"
#include "priceFuncs.h"
#include "optimFuncs.h"
#include "InputReader.h"
#include <algorithm>
#include <fstream>

#include <nlopt.hpp>

using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;

int main()
{
	//StockDataPack dataPack;

	//vector<string> stockNames;
	//vector<P_PRES> stockPrices;
	//vector<P_PRES> stockVols;
	//vector<vector<P_PRES> > stockReturns;

	time_t beginT1 = time( 0 );

	InputReader inpReader;
	inpReader.readStockData( "C:\\Users\\dchernikov\\Documents\\OptionPaper\\newfile.csv"
								, &dataPack.stockNames, &dataPack.stockPrices, &dataPack.stockVols, &dataPack.stockReturns );
	dataPack.rebalanceTime = 5.0 / 252.0;	//rebalance every week
	dataPack.expTime = 1.5;				//options will mature in about a year and a half
	dataPack.riskFreeRate = 0.0072;		//12-month goverment bond yield as a rist-free rate

	vector<P_PRES> xx( 2 * dataPack.stockNames.size(), 1.0 );
	vector<P_PRES> gg( 2 * dataPack.stockNames.size(), 0.0 );
	for( int i = dataPack.stockNames.size(); i < xx.size(); ++i )
	{
		xx[i] = dataPack.stockPrices[i - dataPack.stockNames.size()];	//initial strikes are current prices of stocks
	}

	ifstream ifs( "MMAoptPoint.txt" );
	for( int i = 0; i < dataPack.stockNames.size() * 2; ++i )
	{
		ifs >> xx[i];
	}
	ifs.close();

	//cout << " Start running options algorithm:\n";
	//time_t begin = time( 0 );
	//optionPortfObjCUDA( xx, gg, ( void* )( &dataPack ) );
	//optionPortfObj2( xx, gg, ( void* )( &dataPack ) );
	//budgetConstr( xx, gg, ( void* )( &dataPack ) );
	//cout << " :: done in " << time( 0 ) - begin << endl;
	//cout << " --------------------\n";


	///////////////////////////////////
	// ACTUAL OPTIMIZAION
	///////////////////////////////////

	int dim = dataPack.stockNames.size() * 2;

	nlopt::opt opt( nlopt::LD_MMA, dim );
	//nlopt::opt opt( nlopt::LD_SLSQP, dim );
	//nlopt::opt opt( nlopt::GN_ISRES, dim );
	//nlopt::opt opt( nlopt::LD_LBFGS, dim );

	std::vector<P_PRES> lb( dim );
	for( int i = 0; i < dataPack.stockNames.size(); ++i )
	{
		lb[i] = 0.0;
	}
	for( int i = dataPack.stockNames.size(); i < dim; ++i )
	{
		lb[i] = dataPack.stockPrices[i - dataPack.stockNames.size()] / 2.0;
	}
	opt.set_lower_bounds( lb );

	std::vector<P_PRES> ub( dim );
	for( int i = 0; i < dataPack.stockNames.size(); ++i )
	{
		ub[i] = 1000.0;
	}
	for( int i = dataPack.stockNames.size(); i < dim; ++i )
	{
		ub[i] = dataPack.stockPrices[i - dataPack.stockNames.size()] * 2.0;
	}
	opt.set_upper_bounds( ub );

	opt.add_inequality_constraint( budgetConstr, (void*)&dataPack, 1e-3 );

	opt.set_max_objective( optionPortfObj2, (void*)&dataPack );

	//opt.set_xtol_rel( 1e-10 );
	//opt.set_ftol_rel( 1e-10 );
	//opt.set_stopval( 10451.8 );
	opt.set_maxtime( 7200 * 2 * 2 * 2 );

	P_PRES minf;
	nlopt::result result;
	try
	{
		result = opt.optimize( xx, minf );
	}
	catch( int e )
	{
		cout << " ERROR " << e << endl;
	}

	cout << " opt point: \n";

	for( int i = 0; i < dim; ++i )
	{
		cout << "\t" << xx[i] << endl;
	}

	cout << " opt val: \n";
	cout << "\t" << minf << endl;

	cout << " done in " << time( 0 ) - beginT1 << endl;

	ofstream off( "optimResults.txt" );
	off << result << endl;
	off << minf << endl;

	vector<P_PRES> fakeG;
	off << budgetConstr( xx, fakeG, ( void* )( &dataPack ) ) << endl;

	optionPortfObj2( xx, gg, ( void* )( &dataPack ) );
	P_PRES gradNorm = 0.0;
	for( int i = 0; i < dim; ++i )
	{
		gradNorm += gg[i] * gg[i];
	}
	gradNorm = sqrt( gradNorm );
	off << " grad l2 norm: \n";
	off << gradNorm << endl;

	off << " opt point: \n";
	for( int i = 0; i < dim; ++i )
	{
		off << xx[i] << endl;
	}
	off.close();

	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

	///////////////////////////////////
	// END OF ACTUAL OPTIMIZATION
	///////////////////////////////////

	///////////////////////////////////
	// TESTING
	///////////////////////////////////

	//int dim = 3;

	//nlopt::opt opt( nlopt::LD_MMA, dim );
	////nlopt::opt opt( nlopt::LD_SLSQP, dim );
	////nlopt::opt opt( nlopt::GN_ISRES, dim );
	////nlopt::opt opt( nlopt::LD_LBFGS, dim );

	//std::vector<P_PRES> lb( dim );
	//for( int i = 0; i < dim; ++i )
	//{
	//	lb[i] = -HUGE_VAL;
	//}
	//opt.set_lower_bounds( lb );

	//opt.set_min_objective( rosenbrockFunc, NULL );

	//opt.set_xtol_rel( 1e-14 );

	//std::vector<P_PRES> x( dim );
	//x[0] = -3.0;
	//for( int i = 1; i < dim; ++i )
	//{
	//	x[i] = -4.0;
	//}
	//P_PRES minf;
	//nlopt::result result;
	//try
	//{
	//	result = opt.optimize( x, minf );
	//}
	//catch( int e )
	//{
	//	cout << " ERROR " << e << endl;
	//}

	//cout << " opt point: \n";

	//for( int i = 0; i < dim; ++i )
	//{
	//	cout << "\t" << x[i] << endl;
	//}

	//cout << " opt val: \n";
	//cout << "\t" << minf << endl;
	//std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

	///////////////////////////////////
	// END OF TESTING
	///////////////////////////////////

	
	//ifs.open( "stockPricesForBackTest.txt" ); 
	//ofstream ofs( "resultingOptionPrices.txt" );
	//for( int i = 0; i < dataPack.stockNames.size(); ++i )
	//{
	//	cout << i << endl;

	//	P_PRES curSinit;
	//	P_PRES curSfin;
	//	P_PRES curSigma;
	//	P_PRES curE;

	//	ifs >> curSinit >> curSfin >> curSigma >> curE;	

	//	//Market option prices on CBOE and NASDAQUE can be different. 
	//	P_PRES curT = 0;//dataPack.rebalanceTime;
	//	P_PRES expT = dataPack.expTime;// - dataPack.rebalanceTime;		//It seems that the calculated price agrees more with the market price if the maturity date is father in the future
	//	//P_PRES S = 72.43;		//Spot price
	//	//P_PRES E = 57.5341;			//Strike price
	//	//P_PRES r = 0.0047;	//risk-free rate of return (3-month)
	//	//P_PRES r = 0.0058;		//risk-free rate of return (6-month)
	//	P_PRES r = dataPack.riskFreeRate;		//risk-free rate of return (12-month)
	//	//P_PRES sigma = 0.240198796117542;	//volatility of the stock

	//	//set up the derivative information
	//	//HPD<P_PRES, 1> Eh = E;
	//	//Eh.elems[1] = 1.0;

	//	vector<HPD<P_PRES, 1> > optionPrice;
	//	vector<P_PRES> Sarr;

	//	//begin = time( 0 );
	//	priceAmerPut2<HPD<P_PRES, 1> >( expT, curSinit, curSinit, curE, r, curSigma, &optionPrice, &Sarr );
	//	//cout << " option price is " << optPriceAt<HPD<P_PRES, 1> >( S, optionPrice, Sarr ) << endl;
	//	HPD<P_PRES, 1> optPriceInit = optPriceAt<HPD<P_PRES, 1> >( curSinit, optionPrice, Sarr );

	//	priceAmerPut2<HPD<P_PRES, 1> >( expT - dataPack.rebalanceTime, curSfin, curSfin, curE, r, curSigma, &optionPrice, &Sarr );
	//	HPD<P_PRES, 1> optPriceFin = optPriceAt<HPD<P_PRES, 1> >( curSfin, optionPrice, Sarr );

	//	ofs << optPriceInit.real() << "\t" << optPriceFin.real() << endl;
	//	//P_PRES dE = 5e-11;
	//	//priceAmerPut<HPD<P_PRES, 1> >( expT, S, E - dE, r, sigma, &optionPrice, &Sarr );
	//	//HPD<P_PRES, 1> priceLeft = optPriceAt<HPD<P_PRES, 1> >( S, optionPrice, Sarr );
	//	//priceAmerPut<HPD<P_PRES, 1> >( expT, S, E + dE, r, sigma, &optionPrice, &Sarr );
	//	//HPD<P_PRES, 1> priceRight = optPriceAt<HPD<P_PRES, 1> >( S, optionPrice, Sarr );
	//	//cout << " finite difference deriv is " << ( priceRight.real() - priceLeft.real() ) / 2.0 / dE << endl;

	//	//vector<HPD<P_PRES, 1> > resArr( 1000, 0.0 );
	//	//reinterpolate( std::max( 2.0 * E, 2.0 * S ), optionPrice, Sarr, &resArr );	//Reinterpolate to even intervals of S
	//
	//	//cout << "done in " << time( 0 ) - begin << endl;

	//	//cout << optionPriceBench( expT, 10.025, E, r, sigma, 0.0, 0, 1 ) << endl;

	//	//HPD<P_PRES, 1> res1 = priceEuroPut( curT, expT, S, Eh, r, sigma );
	//	//cout << ( priceEuroPut( curT, expT, S, E + dE, r, sigma ) - priceEuroPut( curT, expT, S, E - dE, r, sigma ) ) / 2.0 / dE << " " << res1.elems[1] << endl;
	//	//cout << dPdE( curT, expT, S, E, r, sigma ) << endl;
	//}
	//ifs.close();
	//ofs.close();

	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}

#endif