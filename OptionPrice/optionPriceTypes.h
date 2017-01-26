#ifndef _OPTION_PRICE_TYPES_
#define _OPTION_PRICE_TYPES_ 1

#define P_PRES double

#define LOWEST_PRICE 1e-5
#define HIGH_PRICE_FACTOR 10.0	//we assume that the highest value the price of a stock can go is this factor by the strike price
#define LOW_PRICE_FACTOR 20.0
//#define DS_AT_0 1	//was 5e-2
//#define DT_WANTED 1e-2	//was 1e-3
#define NUM_OF_TIME_STEPS 300
#define NUM_OF_SPACE_STEPS 3000

#define TOTAL_BUDGET 10000.0	//no longer need it when solving for weights, not quantities

#define PSOR_ERROR_EPS 1e-6	//was 1e-6

#endif