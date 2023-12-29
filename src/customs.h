#pragma once
#ifndef customs_h
#define customs_h

#include <vector>
#include <tuple>
#include <chrono>
#include "Xoshiro.hpp"

template<class T>
using vec = std::vector<T>;

typedef std::tuple<double, vec<double>> mypair;
typedef std::chrono::steady_clock::time_point clock_point; //now defined in simplex.h
typedef XoshiroCpp::Xoshiro256PlusPlus xpp; //now defined in simplex.h

//max number of iterations and tolerance for amoeba
const int MAX_ITER = 5000;
const double TOL = 1e-15;

const double DEX = 0.25; //the number of decades to vary the scaling regime over.

//minimum number of scaling events required to bootstrap
//set to 35, as that is when the relative error caused by sample size
//is ~20% in the ideal case.

//ideal minimum number of events is 100, since that is when the relative error is around 10%.
const size_t MIN_SIZE = 35;

//number of runs for  bootstrap
const int NUM_RUNS = 10000;

//the max number of times to try a single bootstrap run
const int CTR_MAX = 20;


//mathematical constants
const double e = std::exp(1.0);
const double pi = std::atan(1) * 4;

//can define static functions here if they are useful in every function.

//**HELPER FUNCTIONS DEFINED IN HEADERS MUST BE STATIC!**
static void hello();

static void hello()
{
	printf("hello!");
}



#endif // !customs

