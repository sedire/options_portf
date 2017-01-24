#include "InputReader.h"

void InputReader::readStockData( const string& fname, 
	vector<string>* stockNames, vector<P_PRES>* lastStockPrices, vector<P_PRES>* stockVols, vector<vector<P_PRES>>* stockHistReturns )
{
	ifstream ifs( fname );
	string line;
	if( ifs.is_open() )
	{
		//WARNING: I don't do any error checks here! I also assume that the formatting of the file is proper.
		//read the names of the stocks
		getline( ifs, line );

		stringstream parseString( line );
		string stockName;

		getline( parseString, stockName, ',' );			//the first one is empty
		while( getline( parseString, stockName, ',' ) )
		{
			stockNames->push_back( string( stockName, 0, stockName.length() ) );
		}

		//read the last Friday prices of the selected stocks
		getline( ifs, line );

		parseString = stringstream( line );
		string stockPrice;

		getline( parseString, stockPrice, ',' );			//the first one is an auxiliary number
		while( getline( parseString, stockPrice, ',' ) )
		{
			lastStockPrices->push_back( stod( stockPrice ) );
		}

		//read the volatilities of the stocks
		getline( ifs, line );

		parseString = stringstream( line );
		string stockVol;

		getline( parseString, stockVol, ',' );			//the first one is an auxiliary number
		while( getline( parseString, stockVol, ',' ) )
		{
			stockVols->push_back( stod( stockVol ) );
		}

		//now read the historical returns
		stockHistReturns->resize( stockVols->size() );
		while( getline( ifs, line ) )
		{
			parseString = stringstream( line );
			string stockRet;

			getline( parseString, stockRet, ',' );			//the first one is an auxiliary number
			int index = 0;
			while( getline( parseString, stockRet, ',' ) )
			{
				( *stockHistReturns )[index].push_back( stod( stockRet ) );
				++index;
			}
		}
		ifs.close();
	}
	else
	{
		cout << "couldn't open the file\n";
	}

	return;
}