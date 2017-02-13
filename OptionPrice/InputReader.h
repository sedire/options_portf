#ifndef _INPUT_READER_
#define _INPUT_READER_ 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip> 
#include <algorithm>
#include "optionPriceTypes.h"
#include "priceFuncs.h"

using std::sort;
using std::vector;
using std::getline;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::stringstream;
using std::setprecision;

template<class PR_NUM>
class StockDataPack
{
public:
	StockDataPack() {};
	void doPreCalc();

	P_PRES rebalanceTime;
	P_PRES expTime;
	P_PRES riskFreeRate;
	vector<string> stockNames;
	vector<P_PRES> stockPrices;
	vector<P_PRES> stockVols;
	vector<vector<P_PRES> > stockReturns;	//first ind = stock, second ind = scen

	//int voidPack( void** pack );
	//members for option price pre-calculation
	vector<vector<P_PRES> > stockPricesScenSorted;	//first ind is stock, second is scenario
	vector<vector<int> > stockPricesScenOrder;
	vector<P_PRES> Lmins;
	vector<P_PRES> dxs;
	vector<vector<PR_NUM> > us;		//first index = stock number, second ind = position on the computational grid
									//(this is a numerical solution of the auxiliary PDE)
	//same for the initial option prices
	vector<P_PRES> Lmins0;
	vector<P_PRES> dxs0;
	vector<vector<PR_NUM> > us0;
private:
	void bubbleSortWithRank( vector<P_PRES>& v, vector<int>& r );	//sorts in ascending order and gives the rank also
};

template<class PR_NUM>
void StockDataPack<PR_NUM>::doPreCalc()
{
	//Calc all possible scenario stock prices based on returns,
	//sort them and keep the original order stored separately.
	stockPricesScenSorted.resize( stockNames.size() );
	stockPricesScenOrder.resize( stockNames.size() );
	for( int stock = 0; stock < stockNames.size(); ++stock )
	{
		for( int scen = 0; scen < stockReturns[0].size(); ++scen )
		{
			stockPricesScenSorted[stock].push_back( stockReturns[stock][scen] * stockPrices[stock] );
			//cout << " " << stockPricesScenSorted[stock][scen];
		}
		//cout << endl;

		bubbleSortWithRank( stockPricesScenSorted[stock], stockPricesScenOrder[stock] );

		//for( int scen = 0; scen < stockReturns[0].size(); ++scen )
		//{
		//	//sum += ( pricesCopy[scen] - stockPricesScenSorted[stock][scen] ) * ( pricesCopy[scen] - stockPricesScenSorted[stock][scen] );
		//	cout << " " << stockPricesScenOrder[stock][scen];
		//}
		//cout << endl;
		//cout << " --------------------------\n";
	}

	//do the pre-calc for the SCENARIO prices
	Lmins.resize( stockNames.size() );
	dxs.resize( stockNames.size() );
	us.resize( stockNames.size() );
	int ret = 0;
	for( int stock = 0; stock < stockNames.size(); ++stock )
	{
		ret += prePriceAmerPut( expTime - rebalanceTime
							, stockPricesScenSorted[stock][0]
							, stockPricesScenSorted[stock][stockPricesScenSorted[stock].size() - 1]
							, stockPrices[stock] / 2.0, stockPrices[stock] * 2.0, riskFreeRate, stockVols[stock]
							, &( us[stock] ), &( Lmins[stock] ), &( dxs[stock] ) );
	}

	//do the pre-calc for the INITIAL prices
	Lmins0.resize( stockNames.size() );
	dxs0.resize( stockNames.size() );
	us0.resize( stockNames.size() );
	for( int stock = 0; stock < stockNames.size(); ++stock )
	{
		ret += prePriceAmerPut( expTime
							, stockPricesScenSorted[stock][0]
							, stockPricesScenSorted[stock][stockPricesScenSorted[stock].size() - 1]
							, stockPrices[stock] / 2.0, stockPrices[stock] * 2.0, riskFreeRate, stockVols[stock]
							, &( us0[stock] ), &( Lmins0[stock] ), &( dxs0[stock] ) );
	}
}

template<class PR_NUM>
void StockDataPack<PR_NUM>::bubbleSortWithRank( vector<P_PRES>& v, vector<int>& r )
{
	r.resize( v.size() );
	//create the initial order
	for( int i = 0; i < r.size(); ++i )
	{
		r[i] = i;
	}

	int swapped = 0;
	do
	{
		swapped = 0;
		P_PRES tmp;
		int tmpR;
		for( int i = 0; i < v.size() - 1; ++i )
		{
			if( v[i] > v[i + 1] )
			{
				tmp = v[i];
				v[i] = v[i + 1];
				v[i + 1] = tmp;

				tmpR = r[i];
				r[i] = r[i + 1];
				r[i + 1] = tmpR;

				swapped = 1;
			}
		}
	} while( swapped != 0 );
}

class InputReader
{
public:
	void readStockData( const string& fname, vector<string>* stockNames, vector<P_PRES>* lastStockPrices, vector<P_PRES>* stockVols, vector<vector<P_PRES>>* stockHistReturns );  
};

#endif