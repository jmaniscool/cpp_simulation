#pragma once

#ifndef helpers_h
#define helpers_h

#include "Xoshiro.hpp"
#include "customs.h"
#include <vector>
#include <fstream>
#include <chrono>
#include <numeric> //std::iota std::accumulate std::reduce
#include <algorithm> //std::sort, std::stable_sort
//#include "gsl/gsl_interp.h" //required for log vector functions, commented out for simulation helper implementation
#include <boost/interprocess/file_mapping.hpp> //for file mapping (speed up loading and saving)
#include <boost/interprocess/mapped_region.hpp> //for mapped region object

//define helpers namespace for helper functions (i.e. generating random vectors)

namespace helpers
{

	//use this to determine if a file exists
	static int file_exists(const char* filename)
	{
		FILE* fp;
		fopen_s(&fp, filename, "r");
		int is_exist = 0;
		if (fp != NULL)
		{
			is_exist = 1;
			fclose(fp); // close the file
		}
		return is_exist;
	}

	//read in the simulation state, assumed to be in format from write_txt. Uses memory mapping to make reads more efficient.
	static vec<double> read_simulation_state(const char* fname)
	{
		boost::interprocess::file_mapping fileMapping(fname, boost::interprocess::read_only); //open file in memory map for read only

		boost::interprocess::mapped_region mappedRegion(fileMapping, boost::interprocess::read_only); //allocate file to mapped region

		std::string charData = static_cast<std::string>(static_cast<char*>(mappedRegion.get_address())); //convert the mappedRegion to a string

		std::istringstream stream(charData); //turn into a stringstream and treat it as in memory

		int count = std::count(charData.begin(), charData.end(), '\n'); //count the number of newline characters
		vec<double> out(count, -1); //initialize vector to hold that much data

		std::string hold_line;
		int i = 0;
		while (std::getline(stream, hold_line))
		{
			double val = std::stof(hold_line);
			out[i] = val;
			i++;
		}


		return out;
	}


	//do the np.linspace. From https://stackoverflow.com/questions/27028226/python-linspace-in-c
	template<typename T>
	static std::vector<double> linspace(T start_in, T end_in, int num_in);

	//write to a file fname.
	//template<class T>
	//static void write_txt(std::string& fname, std::vector<T>& indat);


	//get a random vector of type class T. Initial testing suggests this is about 10x faster than ranmarin() and may be statistically valid enough for our applications
	template<class T>
	static vec<T> ran_vec(double st, double en, int len)
	{
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		XoshiroCpp::Xoshiro256PlusPlus rng(seed);
		vec<T> out(len);

		double range = en - st;

		for (int i = 0; i < len; i++)
			out[i] = static_cast<T>(range * XoshiroCpp::DoubleFromBits(rng()) + st);

		return out;
	}

	//return a single random number
	template<class T>
	static T ran_num(double st, double en)
	{
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		XoshiroCpp::Xoshiro256PlusPlus rng(seed);
		double range = en - st;
		return static_cast<T>(range * XoshiroCpp::DoubleFromBits(rng()) + st);
	}


	static vec<double> logvec(vec<double> x)
	{
		using namespace std;
		int len = x.size();
		vec<double> out(len, -1);
		for (int i = 0; i < len; i++)
			out[i] = log(x[i]);
		return out;
	}


	//do a linear fit on vectors of points (x,y).
	//algorithm explicitly calculates the minimized sum of least squares slope
	//algorithm from https://web.archive.org/web/20150715022401/http://faculty.cs.niu.edu/~hutchins/csci230/best-fit.htm
	static double linfit(vec<double>& x, vec<double>& y)
	{
		using namespace std;
		double sumx = 0;
		double sumy = 0;
		double sumxy = 0;
		double sumx2 = 0;

		double xlen = (double)x.size(); //cast the size to a double

		//calculate the mean and sum of products of X*Y
		for (int i = 0; i < xlen; i++)
		{
			sumx += x[i];
			sumx2 += x[i] * x[i];
			sumy += y[i];
			sumxy += x[i] * y[i];
		}

		//calculate the slope
		double slope = (sumxy - sumx * sumy / xlen) / (sumx2 - sumx * sumx / xlen);
		//double yint = sumy / xlen - slope*sumx / xlen;
		return slope;
	}


	/*
	//get a linear interpolation between x and y
	//obtained from https://lists.gnu.org/archive/html/help-gsl/2007-06/msg00019.html
	static double interp(vec<double>& x, vec<double>& y, double interp_val)
	{
		//sample data
		int n = x.size(); //length of the vectors
		double* xarr = &x[0]; //this should give a pointer to the array contained in the vectors
		double* yarr = &y[0];

		//inialise and allocate the gsl objects
		gsl_interp* interpolation = gsl_interp_alloc(gsl_interp_linear, n);
		gsl_interp_init(interpolation, xarr, yarr, n);
		gsl_interp_accel* accelerator = gsl_interp_accel_alloc();

		//get interpolation for x = 1981
		double value = gsl_interp_eval(interpolation, xarr, yarr, interp_val, accelerator);

		return value;
	}

	*/

	static vec<double> log10vec(vec<double>& x)
	{
		using namespace std;
		int len = x.size();
		double log10 = log(10);
		vec<double> out(len, -1);
		for (int i = 0; i < len; i++)
			out[i] = log(x[i]) / log10;
		return out;
	}


	static double logloginterp(vec<double>& x, vec<double>& y, double interp_point)
	{
		vec<double> logx = logvec(x);
		vec<double> logy = logvec(y);
		//double val = interp(logx, logy, std::log(interp_point));
		//return std::exp(val);
		return -1;
	}

	//perform an argsort, which returns the indices that would sort vector x.
	// implemented from https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
	//template <typename T>
	//static vec<size_t> argsort(vec<T>& x);

	template <typename T>
	static vec<size_t> argsort(vec<T>& x)
	{
		// initialize original index locations
		using namespace std;
		vec<size_t> idx(x.size());
		iota(idx.begin(), idx.end(), 0);

		// sort indexes based on comparing values in v
		// using std::stable_sort instead of std::sort
		// to avoid unnecessary index re-orderings
		// when v contains elements of equal values
		stable_sort(idx.begin(), idx.end(), [&x](size_t i1, size_t i2) {return x[i1] < x[i2]; });

		return idx;
	}

	//get the closest number to a given index
	//from https://stackoverflow.com/questions/8647635/elegant-way-to-find-closest-value-in-a-vector-from-above
	//template <typename T>
	//static int closest(const std::vector<T>& vec, T value);

	template <typename T>
	static int closest(const std::vector<T>& vec, T value) {
		T tmp = value + 1e-12; //add a small amount so lower_bound() returns the first idx >= val.
		auto const it = std::lower_bound(vec.begin(), vec.end(), tmp);
		if (it == vec.end()) { return vec.size() - 1; } //if past the end, just return the end index.
		if (it == vec.begin()) { return 0; }

		return std::distance(vec.begin(), it) - 1;
	}


	//implement logbinning. From sets of scatter points x and y
	static std::tuple<vec<double>, vec<double>> logbinning(vec<double>& x, vec<double>& y, int num_bins)
	{
		using namespace std;
		int len = x.size();
		//get the indices that sort x
		vec<double> srtx(len, -1);
		vec<double> srty(len, -1);
		vec<size_t> idxs(len, 1);

		//get the sorted indices
		idxs = argsort(x);

		//sort according to x
		for (int i = 0; i < x.size(); i++)
		{
			//cout << idxs[i] << endl;
			srtx[i] = x[idxs[i]];
			srty[i] = y[idxs[i]];
		}

		double logmax = log(srtx[len - 1]);
		double logmin = log(srtx[1]); //START AT INDEX 1 BECAUSE THAT IS WHAT IS DONE IN PYTHON

		//get the linspace of log
		vec<double> binEdges = linspace(logmin, logmax, num_bins + 1);
		vec<int> binEdgesIdx(num_bins + 1, -1);

		//get the bin edges in real space
		for (int i = 0; i < num_bins + 1; i++)
		{
			binEdges[i] = exp(binEdges[i]);
			//printf("%.3f\n", binEdges[i]);
		}

		//find the values of srtx that are closes to each bin edge value.

		for (int i = 0; i < num_bins + 1; i++)
		{
			binEdgesIdx[i] = closest(srtx, binEdges[i]) - 1; //-1 MINUS ONE ACCOUNTS FOR OFF-BY-ONE ERROR IN PYTHON PROGRAM (?)
			//printf("%d\n", binEdgesIdx[i]);
		}


		//Also find the bin centers. The halfway point between bin edges in log-space.
		vec<double> binCenters(num_bins, -1);
		for (size_t i = 0; i < num_bins; i++)
		{
			binCenters[i] = exp((log(binEdges[i]) + log(binEdges[i + 1])) / 2.0);
		}


		//take the mean of the data that resides in each bin

		vec<double> avgs(num_bins, -1);
		for (size_t i = 0; i < num_bins; i++)
		{
			size_t idx_start = binEdgesIdx[i];
			size_t idx_end = binEdgesIdx[i + 1];
			//cout << "i = " << i << " starts " << idx_start << " to " << idx_end << endl;
			if (idx_start == idx_end)
				idx_end = idx_end + 1; //so there is at least one element to add
			//get the avgs[i] = sum(binCenters[idx_start:idx_en])/((double) (idx_en-idx_st))
			avgs[i] = reduce(srty.begin() + idx_start, srty.begin() + idx_end) / ((double)(idx_end - idx_start)); //get the averages
			//cout << "i = " << i << " average " << avgs[i] << endl;
		}

		auto out = make_tuple(binCenters, avgs);

		return out;






	}

	template<typename T>
	static std::vector<double> linspace(T start_in, T end_in, int num_in)
	{

		std::vector<double> linspaced;

		double start = static_cast<double>(start_in);
		double end = static_cast<double>(end_in);
		double num = static_cast<double>(num_in);

		if (num == 0) { return linspaced; }
		if (num == 1)
		{
			linspaced.push_back(start);
			return linspaced;
		}

		double delta = (end - start) / (num - 1);

		for (int i = 0; i < num - 1; ++i)
		{
			linspaced.push_back(start + delta * i);
		}
		linspaced.push_back(end); // I want to ensure that start and end
		// are exactly the same as the input
		return linspaced;
	}

	//Return if the current bootstrap limits are valid. 
	// In python, these conditions are:
	// the length of sc, dc, and vc are greater than 3
	// 
	//

	//get the confidence intervals
	//static vec<double> confidence_intervals(vec<double>& boot, double ci);

	static vec<double> confidence_intervals(vec<double>& boot, double ci)
	{
		//ci = percentage of measured values that are within (lo, hi)
		vec<double> out(2, -1); //outputs as lo, hi (no median)
		double mysize = (double)boot.size();
		vec<double> holder(mysize, -1);//holds sorted list
		for (int i = 0; i < mysize; i++)
			holder[i] = boot[i];

		std::sort(holder.begin(), holder.end()); //sort the list
		//out[0] = holder[(int)(mysize / 2.0)]; //50th percentile
		out[0] = holder[(int)(mysize * (ci / 2))]; //lower = 50 - ci/2
		out[1] = holder[(int)(mysize * (1 - ci / 2))]; //upper = 50 + ci/2

		//in ths way, we guarantee the confidence interval contains ci% of the total data
		return out;
	}

	//get 1 sigma and 2 sigma CIs more quickly
	//static vec<double> get_cis(vec<double>& boot);


	//
	static vec<double> get_cis(vec<double>& boot)
	{
		vec<double> out(4, -1);
		double mysize = (double)boot.size();
		vec<double> holder(mysize, -1);//holds sorted list
		for (int i = 0; i < mysize; i++)
			holder[i] = boot[i];

		std::sort(holder.begin(), holder.end()); //sort the list

		out[0] = holder[(int)(mysize * (0.32 / 2))];
		out[1] = holder[(int)(mysize * (0.05 / 2))];
		out[2] = holder[(int)(mysize * (1 - 0.32 / 2))];
		out[3] = holder[(int)(mysize * (1 - 0.05 / 2))];


		//returns out = [lo, lo95, hi, hi95]
		return out;
	}


	//static vec<double> read_csv(const char* fname);


	//function for reading from pandas DataFrame in C++.
	static vec<double> read_pandas(const char* fname)
	{
		std::ifstream reader(fname);
		std::string line; //holds the whole line until \n is reached.
		std::string elem; //holds the part of the line until , is reached.
		double s; //current size
		vec<double> out; //the vector of doubles.

		//ignore the first 10,000 characters, or until delimiter \n is reached. Repeat for as many times as there are headers.
		//sizes.csv has one header line, so this line is used once before the while loop.
		reader.ignore(10000, '\n');

		//while file is open
		while (std::getline(reader, line))
		{
			//initialize ss as substring object
			std::istringstream ss(line);
			ss.ignore(10000, ','); //ignore the first row (just dumb pandas stuff). 10k characters chosen to be long enought that it will always reach the comma first
			std::getline(ss, elem, ','); //get the second row and pipe it into elem until a comma is reached
			s = std::stod(elem);
			out.push_back(s);
		}

		//return the out vector
		return out;
	}

	//write to a .txt
	template<class T>
	static void write_txt(std::string& fname, std::vector<T>& indat)
	{
		//C implementation (0.32 us per write)
		auto start = std::chrono::high_resolution_clock::now();

		FILE* file; //get the file

		const char* fname_c = fname.c_str(); //get the c string

		if (fopen_s(&file, fname_c, "w") == 0) {

			for (int i = 0; i < indat.size(); i++) {
				//fprintf(file, "%.*f\n", DBL_DECIMAL_DIG-1, indat[i]); //minus 1 for the decimal point if you want full precision
				fprintf(file, "%f\n", indat[i]); //precision given to 5 or 6 decimal places (fine for pretty much anything)
			}
			fclose(file);
		}
		else {
			// Handle error if unable to open the file
			printf("Unable to open the file: %s\n", fname);
		}


		//C++ implementation (~4.8 us per write)
		/*
		std::ofstream outFile(fname);



		//if it's open, then do the cout 			
		if (outFile.is_open())
		{
			outFile << std::setprecision(std::numeric_limits<double>::max_digits10); //set the precision. Probably a bit too high.
			for (int i = 0; i < indat.size(); i++)
				 outFile << indat[i] << std::endl;
			std::cout << "Written successfully!" << std::endl;
		}
		else
		{
			std::cout << "Cannot open file! " << std::endl;
		}

		*/

		auto elapsed = std::chrono::high_resolution_clock::now() - start;

		auto time_per_write = ((double)std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count()) / indat.size();

		std::cout << "Write finished in " << std::chrono::duration_cast<std::chrono::microseconds>(
			elapsed).count() << " us." << std::endl;

		std::cout << "Write took " << time_per_write << " us per element." << std::endl;

		return;

	}


}



#endif