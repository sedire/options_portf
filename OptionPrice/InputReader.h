#ifndef _INPUT_READER_
#define _INPUT_READER_ 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip> 
#include "optionPriceTypes.h"

using std::vector;
using std::getline;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::stringstream;
using std::setprecision;

class StockDataPack
{
public:
	StockDataPack() {};

	P_PRES rebalanceTime;
	P_PRES expTime;
	P_PRES riskFreeRate;
	vector<string> stockNames;
	vector<P_PRES> stockPrices;
	vector<P_PRES> stockVols;
	vector<vector<P_PRES> > stockReturns;	//first ind = stock, second ind = scen

	//int voidPack( void** pack );
};

class InputReader
{
public:
	void readStockData( const string& fname, vector<string>* stockNames, vector<P_PRES>* lastStockPrices, vector<P_PRES>* stockVols, vector<vector<P_PRES>>* stockHistReturns );  
};

#endif