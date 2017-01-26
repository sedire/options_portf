#include "optimFuncs.h"

P_PRES rosenbrockFunc( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data )
{
	P_PRES ret = 0;
	for( int i = 0; i < x.size() - 1; ++i )
	{
		ret += 100.0 * ( x[i + 1] - x[i] * x[i] ) * ( x[i + 1] - x[i] * x[i] ) + ( 1 - x[i] ) * ( 1 - x[i] );
	}

	grad[0] = -400.0 * x[0] * ( x[1] - x[0] * x[0] ) - 2.0 * ( 1.0 - x[0] );
	for( int i = 1; i < x.size() - 1; ++i )
	{
		grad[i] = -400.0 * x[i] * ( x[i + 1] - x[i] * x[i] ) - 2.0 * ( 1.0 - x[i] ) + 200.0 * ( x[i] - x[i - 1] * x[i - 1] );
	}
	grad[x.size() - 1] = 200.0 * ( x[x.size() - 1] - x[x.size() - 2] * x[x.size() - 2] );

	return ret;
}

P_PRES optionPortfObj( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data )	//in x put quantities first, then strikes
{
	time_t beginT = time( 0 );
	cout << " ::obj called at  " << x[0] << " ... " << x[x.size() / 2 - 1] << " " << x[x.size() / 2] << " ... " << x[x.size() - 1] << endl;
	if( f_data == 0 )
	{
		cout << "ERROR in optionPortfObj: passed data is null\n";
		return 0;
	}
	StockDataPack* dataPack = static_cast<StockDataPack*>( f_data );
	int stockNum = dataPack->stockNames.size();
	int scenNum = dataPack->stockReturns[0].size();
	P_PRES scenProb = 1.0 / scenNum;

	P_PRES maturityTime = dataPack->expTime - dataPack->rebalanceTime;

	//vector<HPD<P_PRES, 1> > optionScenPrices( stockNum, 0.0 );
	vector<HPD<P_PRES, 1> > optionPriceArr;
	vector<P_PRES> Sarr;
	vector<P_PRES> scenStockPrices( scenNum, 0.0 );

	//vector<HPD<P_PRES, 1> > priceU( NUM_OF_SPACE_STEPS + 1 );
	//vector<HPD<P_PRES, 1> > priceG( NUM_OF_SPACE_STEPS - 1 );
	//vector<HPD<P_PRES, 1> > priceB( NUM_OF_SPACE_STEPS - 1 );

	//reset the gradient
	for( int i = 0; i < grad.size(); ++i )
	{
		grad[i] = 0.0;
	}

	//ofstream of( "testVRTX3.txt" );

	P_PRES obj = 0.0;
	for( int stock = 0; stock < stockNum; ++stock )
	{
		for( int scen = 0; scen < scenNum; ++scen )
		{
			scenStockPrices[scen] = dataPack->stockPrices[stock] * dataPack->stockReturns[stock][scen];
		}
		P_PRES stockPriceMin = *min_element( scenStockPrices.begin(), scenStockPrices.end() );
		P_PRES stockPriceMax = *max_element( scenStockPrices.begin(), scenStockPrices.end() );

		//time_t begin2 = time( 0 );
		priceAmerPut2<HPD<P_PRES, 1> >( maturityTime
										, stockPriceMin, stockPriceMax
										, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
										, &optionPriceArr, &Sarr
										/*, &priceU, &priceG, &priceB*/ );
		//cout << "done in " << time( 0 ) - begin2 << endl;

		for( int scen = 0; scen < scenNum; ++scen )
		{
			//optionScenPrices[stock] = optPriceAt<HPD<P_PRES, 1> >( scenStockPrices[scen], optionPriceArr, Sarr );
			HPD<P_PRES, 1> optionScenPrice = optPriceAt<HPD<P_PRES, 1> >( scenStockPrices[scen], optionPriceArr, Sarr );
			//of << optionScenPrices[stock].real() << " " << optionScenPrices[stock].elems[1] << endl;

			if( grad.size() > 0 )
			{
				//now calculate the corresponding terms of the gradient
				//weight derivatives
				grad[stock] += scenProb * optionScenPrice.real();
				//strike derivatives
				grad[stockNum + stock] += scenProb * x[stock] * optionScenPrice.elems[1];
			}
			//now calculate the part of the objective function corresponding to the scenario
			obj += scenProb * x[stock] * optionScenPrice.real();
		}
	}

	//of.close();
	cout << " the obj is " << obj << endl;
	cout << " done in " << time( 0 ) - beginT << endl;

	return obj;
}

P_PRES optionPortfObjNoGrad( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data )	//in x put quantities first, then strikes
{
	time_t beginT = time( 0 );
	cout << " ::obj (no grad) called at  " << x[0] << " ... " << x[x.size() / 2 - 1] << " " << x[x.size() / 2] << " ... " << x[x.size() - 1] << endl;
	if( f_data == 0 )
	{
		cout << "ERROR in optionPortfObj: passed data is null\n";
		return 0;
	}
	StockDataPack* dataPack = static_cast<StockDataPack*>( f_data );
	int stockNum = dataPack->stockNames.size();
	int scenNum = dataPack->stockReturns[0].size();
	P_PRES scenProb = 1.0 / scenNum;

	P_PRES maturityTime = dataPack->expTime - dataPack->rebalanceTime;

	vector<P_PRES> optionPriceArr;
	vector<P_PRES> Sarr;
	vector<P_PRES> scenStockPrices( scenNum, 0.0 );

	//reset the gradient
	for( int i = 0; i < grad.size(); ++i )
	{
		grad[i] = 0.0;
	}

	//ofstream of( "testVRTX3.txt" );

	P_PRES obj = 0.0;
	for( int stock = 0; stock < stockNum; ++stock )
	{
		for( int scen = 0; scen < scenNum; ++scen )
		{
			scenStockPrices[scen] = dataPack->stockPrices[stock] * dataPack->stockReturns[stock][scen];
		}
		P_PRES stockPriceMin = *min_element( scenStockPrices.begin(), scenStockPrices.end() );
		P_PRES stockPriceMax = *max_element( scenStockPrices.begin(), scenStockPrices.end() );

		//cout << " _ " << stockPriceMin << " " << stockPriceMax <<  " " << dataPack->stockPrices[stock] << endl;

		//time_t begin2 = time( 0 );
		priceAmerPut2<P_PRES>( maturityTime
										, stockPriceMin, stockPriceMax
										, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
										, &optionPriceArr, &Sarr
										/*, &priceU, &priceG, &priceB*/ );
		//cout << "done in " << time( 0 ) - begin2 << endl;

		for( int scen = 0; scen < scenNum; ++scen )
		{
			//optionScenPrices[stock] = optPriceAt<HPD<P_PRES, 1> >( scenStockPrices[scen], optionPriceArr, Sarr );
			P_PRES optionScenPrice = optPriceAt<P_PRES>( scenStockPrices[scen], optionPriceArr, Sarr );
			//of << optionScenPrices[stock].real() << " " << optionScenPrices[stock].elems[1] << endl;

			if( grad.size() > 0 )
			{
				cout << "ERROR: grad should be empty in grad-free methods\n";
			}
			//now calculate the part of the objective function corresponding to the scenario
			obj += scenProb * x[stock] * optionScenPrice;
		}
	}

	//of.close();
	cout << " the obj is " << obj << endl;
	cout << " done in " << time( 0 ) - beginT << endl;

	return obj;
}

//Implementing optimization for the case when no gradient is needed
P_PRES optionPortfObjWeights( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data )		//in x put weights first, then strikes
{
	time_t beginT = time( 0 );
	cout << " ::obj called at  " << x[0] << " ... " << x[x.size() / 2 - 1] << " " << x[x.size() / 2] << " ... " << x[x.size() - 1] << endl;
	if( f_data == 0 )
	{
		cout << "ERROR in optionPortfObj: passed data is null\n";
		return 0;
	}
	StockDataPack* dataPack = static_cast<StockDataPack*>( f_data );
	int stockNum = dataPack->stockNames.size();
	int scenNum = dataPack->stockReturns[0].size();
	P_PRES scenProb = 1.0 / scenNum;

	P_PRES maturityTime0 = dataPack->expTime;
	P_PRES maturityTime = dataPack->expTime - dataPack->rebalanceTime;

	vector<P_PRES> Sarr0;
	vector<P_PRES> Sarr;
	vector<P_PRES> scenStockPrices( scenNum, 0.0 );

	vector<HPD<P_PRES, 1> > optionPriceArr0;
	vector<HPD<P_PRES, 1> > optionPriceArr;

	//for the case when no grad is requested
	vector<P_PRES> optionPriceArr0NoGrad;
	vector<P_PRES> optionPriceArrNoGrad;

	//reset the gradient
	for( int i = 0; i < grad.size(); ++i )
	{
		grad[i] = 0.0;
	}

	//ofstream of( "testVRTX3.txt" );

	P_PRES obj = 0.0;
	for( int stock = 0; stock < stockNum; ++stock )
	{
		for( int scen = 0; scen < scenNum; ++scen )
		{
			scenStockPrices[scen] = dataPack->stockPrices[stock] * dataPack->stockReturns[stock][scen];
		}
		P_PRES stockPriceMin = *min_element( scenStockPrices.begin(), scenStockPrices.end() );
		P_PRES stockPriceMax = *max_element( scenStockPrices.begin(), scenStockPrices.end() );

		HPD<P_PRES, 1> optionScenPrice0;
		P_PRES optionScenPrice0NoGrad;

		if( grad.size() > 0 )
		{
			//time_t begin2 = time( 0 );
			priceAmerPut2<HPD<P_PRES, 1> >( maturityTime0
											, dataPack->stockPrices[stock], dataPack->stockPrices[stock]
											, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
											, &optionPriceArr0, &Sarr0 );	//initial option price curve
			priceAmerPut2<HPD<P_PRES, 1> >( maturityTime
											, stockPriceMin, stockPriceMax
											, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
											, &optionPriceArr, &Sarr );		//option price curve at rebalance time
			//cout << "done in " << time( 0 ) - begin2 << endl;

			//get the initial option price from the curve
			optionScenPrice0 = optPriceAt<HPD<P_PRES, 1> >( dataPack->stockPrices[stock], optionPriceArr0, Sarr0 );	
		}
		else	//if no grad is needed
		{
			//time_t begin2 = time( 0 );
			priceAmerPut2<P_PRES>( maturityTime0
											, dataPack->stockPrices[stock], dataPack->stockPrices[stock]
											, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
											, &optionPriceArr0NoGrad, &Sarr0 );	//initial option price curve
			priceAmerPut2<P_PRES>( maturityTime
											, stockPriceMin, stockPriceMax
											, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
											, &optionPriceArrNoGrad, &Sarr );		//option price curve at rebalance time
			//cout << "done in " << time( 0 ) - begin2 << endl;

			//get the initial option price from the curve
			optionScenPrice0NoGrad = optPriceAt<P_PRES>( dataPack->stockPrices[stock], optionPriceArr0NoGrad, Sarr0 );	
		}

		for( int scen = 0; scen < scenNum; ++scen )
		{
			if( grad.size() > 0 )
			{
				HPD<P_PRES, 1> optionScenPrice = optPriceAt<HPD<P_PRES, 1> >( scenStockPrices[scen], optionPriceArr, Sarr );
				//now calculate the corresponding terms of the gradient
				//weight derivatives
				grad[stock] += scenProb * optionScenPrice.real() / optionScenPrice0.real();
				//strike derivatives
				grad[stockNum + stock] += scenProb * x[stock]
										* ( optionScenPrice.elems[1] * optionScenPrice0.real()
										- optionScenPrice0.elems[1] * optionScenPrice.real() ) 
										/ optionScenPrice0.real() / optionScenPrice0.real();
				//now calculate the part of the objective function corresponding to the scenario
				obj += scenProb * x[stock] * optionScenPrice.real() / optionScenPrice0.real();
			}
			else	//no grad case
			{
				P_PRES optionScenPrice = optPriceAt<P_PRES>( scenStockPrices[scen], optionPriceArrNoGrad, Sarr );

				//now calculate the part of the objective function corresponding to the scenario
				obj += scenProb * x[stock] * optionScenPrice / optionScenPrice0NoGrad;
			}
		}
	}

	//of.close();
	cout << " the obj is " << obj << endl;
	cout << " done in " << time( 0 ) - beginT << endl;

	return obj;
}

P_PRES budgetConstr( const vector<P_PRES> &x, vector<P_PRES> &grad, void* data )
{
	cout << " constr called at  " << x[0] << " ... " << x[x.size() / 2 - 1] << " " << x[x.size() / 2] << " ... " << x[x.size() - 1] << endl;
	if( data == 0 )
	{
		cout << "ERROR in optionPortfObj: passed data is null\n";
		return 0;
	}
	StockDataPack* dataPack = static_cast<StockDataPack*>( data );
	int stockNum = dataPack->stockNames.size();
	P_PRES maturityTime = dataPack->expTime;

	P_PRES ret = 0.0;
	//vector<HPD<P_PRES, 1> > optionPrices( stockNum, 0.0 );
	vector<HPD<P_PRES, 1> > optionPriceArr;
	vector<P_PRES> Sarr;

	//vector<HPD<P_PRES, 1> > priceU;
	//vector<HPD<P_PRES, 1> > priceG;
	//vector<HPD<P_PRES, 1> > priceB;


	//reset the gradient
	for( int i = 0; i < grad.size(); ++i )
	{
		grad[i] = 0.0;
	}

	for( int stock = 0; stock < stockNum; ++stock )
	{
		P_PRES stockPrice = dataPack->stockPrices[stock];
		//cout << "\tOption " << dataPack->stockNames[stock] << "\t";

		//time_t begin2 = time( 0 );
		priceAmerPut2<HPD<P_PRES, 1> >( maturityTime
										, stockPrice, stockPrice
										, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
										, &optionPriceArr, &Sarr 
										/*, &priceU, &priceG, &priceB*/ );
		//cout << "done in " << time( 0 ) - begin2 << endl;
		HPD<P_PRES, 1> optionPrice = 0.0;
		optionPrice = optPriceAt<HPD<P_PRES, 1> >( stockPrice, optionPriceArr, Sarr );

		ret += x[stock] * optionPrice.real();
		if( grad.size() > 0 )
		{
			grad[stock] = optionPrice.real();
			grad[stockNum + stock] = x[stock] * optionPrice.elems[1];
		}
	}
	cout << " constr value " << ret - TOTAL_BUDGET << endl;
	cout << " ====================\n";
	return ret - TOTAL_BUDGET;
}

P_PRES budgetConstrNoGrad( const vector<P_PRES> &x, vector<P_PRES> &grad, void* data )
{
	cout << " constr (no grad) called at  " << x[0] << " ... " << x[x.size() / 2 - 1] << " " << x[x.size() / 2] << " ... " << x[x.size() - 1] << endl;
	if( data == 0 )
	{
		cout << "ERROR in optionPortfObj: passed data is null\n";
		return 0;
	}
	StockDataPack* dataPack = static_cast<StockDataPack*>( data );
	int stockNum = dataPack->stockNames.size();
	P_PRES maturityTime = dataPack->expTime;

	P_PRES ret = 0.0;

	vector<P_PRES> optionPriceArr;
	vector<P_PRES> Sarr;

	//reset the gradient
	for( int i = 0; i < grad.size(); ++i )
	{
		grad[i] = 0.0;
	}

	for( int stock = 0; stock < stockNum; ++stock )
	{
		P_PRES stockPrice = dataPack->stockPrices[stock];
		//cout << "\tOption " << dataPack->stockNames[stock] << "\t";

		//time_t begin2 = time( 0 );
		priceAmerPut2<P_PRES>( maturityTime
										, stockPrice, stockPrice
										, x[stockNum + stock], dataPack->riskFreeRate, dataPack->stockVols[stock]
										, &optionPriceArr, &Sarr 
										/*, &priceU, &priceG, &priceB*/ );
		//cout << "done in " << time( 0 ) - begin2 << endl;
		P_PRES optionPrice = 0.0;
		optionPrice = optPriceAt<P_PRES>( stockPrice, optionPriceArr, Sarr );

		ret += x[stock] * optionPrice;
		if( grad.size() > 0 )
		{
			cout << "ERROR: grad should be empty in grad-free method\n";
		}
	}
	cout << " constr value " << ret - TOTAL_BUDGET << endl;
	cout << " ====================\n";
	return ret - TOTAL_BUDGET;
}

P_PRES budgetConstrWeights( const vector<P_PRES> &x, vector<P_PRES> &grad, void* data )
{
	cout << " constr called at  " << x[0] << " ... " << x[x.size() / 2 - 1] << " " << x[x.size() / 2] << " ... " << x[x.size() - 1] << endl;
	if( data == 0 )
	{
		cout << "ERROR in optionPortfObj: passed data is null\n";
		return 0;
	}
	StockDataPack* dataPack = static_cast<StockDataPack*>( data );
	int stockNum = dataPack->stockNames.size();

	if( grad.size() > 0 )
	{
		//reset the gradient
		for( int i = 0; i < stockNum; ++i )
		{
			grad[i] = 1.0;
		}
		for( int i = stockNum; i < grad.size(); ++i )
		{
			grad[i] = 0.0;
		}
	}

	P_PRES ret = 0.0;
	for( int i = 0; i < stockNum; ++i )
	{
		ret += x[i];
	}

	cout << " constr value " << ret - 1.0 << endl;
	cout << " ====================\n";
	return ret - 1.0;
}