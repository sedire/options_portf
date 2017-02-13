#ifndef _OPTIM_FUNCS_
#define _OPTIM_FUNCS_ 1

#include "optionPriceTypes.h"
#include "hyperDual.h"
#include "priceFuncs.h"
#include "InputReader.h"
#include <iostream>
#include <string>
#include <vector>
#include <time.h>  
#include <fstream>
#include <algorithm>

using std::min_element;
using std::max_element;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;

P_PRES rosenbrockFunc( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data );

//extern P_PRES cudaOptionPortfObjWrap( const vector<P_PRES> &x, vector<P_PRES> &grad, StockDataPack* dataPack );

P_PRES optionPortfObj( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data );
P_PRES preOptionPortfObj( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data );
P_PRES optionPortfObjNoGrad( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data );
P_PRES optionPortfObjWeights( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data );
P_PRES optionPortfObjCUDA( const vector<P_PRES> &x, vector<P_PRES> &grad, void* f_data );

P_PRES budgetConstr( const vector<P_PRES> &x, vector<P_PRES> &grad, void* data );
P_PRES budgetConstrNoGrad( const vector<P_PRES> &x, vector<P_PRES> &grad, void* data );
P_PRES budgetConstrWeights( const vector<P_PRES> &x, vector<P_PRES> &grad, void* data );

#endif