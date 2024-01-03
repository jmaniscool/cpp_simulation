/*
Alan Long, 6/20/2019
Last edited: Alan Long, 3/2/2021
Updated to C++ by Jordan Sickle 6/5/2023

This code simulates our mean field model exactly. 

It takes six inputs: time_max is an int that is the length of time you wish to simulate. Area is an int that is the
size of the system in cells. consv is a double that is the conservation parameter c. rate is a double is the strain
rate, set to 0 for an adiabatic system. w is a double that is the spread of arrest stresses with failure stress 
normalized to 1. weakening is a double that is the weakening parameter epsilon.

It write two files as an output. simstress.txt is the force on the system and consists of comma spaced doubles.
simstrain is the strain on the system and also consists of comma spaced doubles.
*/


//C includes
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

//C++ includes (commenting out Python linking libraries)
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h> //for std::vector imports
//#include <pybind11/numpy.h> //so we can export vectors

//local includes
#include "customs.h"
#include "Xoshiro.hpp"
#include "helpers.hpp"
//#include "progressbar.hpp"

//namespace py = pybind11;

//custom types
typedef std::tuple<vec<double>,vec<double>,vec<double>> return_tuple; //return tuple has [force,time,strs]

double wbl(double k, double lambda)
{
	double rand = helpers::ran_num<double>(0, 1);
	return lambda * std::pow(-1 * std::log(rand), 1 / k);
}

//use w = 0.1 to give a source of quenched disorder. Vary k to get different levels of system disorder.
void sim_setup(double area, double w, vec<double>& strs, vec<double>& systemstress, vec<double>& arr_strs, vec<double>& fail_strs)
{
	vec<double> randm = helpers::ran_vec<double>(0, 1, static_cast<int>(area));
	for (int i = 0; i < strs.size(); i++) {
		arr_strs[i] = w * randm[i] - w / 2;//define arrest stresses randomly
		fail_strs[i] = 1;
		strs[i] = (std::pow(i / area, .4) * 1.56585 - .56585) * (fail_strs[i] - arr_strs[i]) + arr_strs[i];//this distribution is an approximation of the steady state distribution. It will take some time to get to the steady state. This can be improved
		systemstress[0] += strs[i]; //initial system stress is sum of all strs
	}
	return;
}

bool between_slips(int& time, double& time_max, double& rate, double& area, double& modulus, double* strs, double* arr_strs, double* fail_strs, double* systemstress, double* strain)
{
	//re-iniitalize fail stress
	for (int i = 0; i < area; i++)
		fail_strs[i] = 1;

	//find the next fail
	double next_fail = 100;
	for (int i = 0; i < area; i++) {
		if (next_fail > (fail_strs[i] - strs[i]))
			next_fail = fail_strs[i] - strs[i];//chose minimum stress to next failure
	}
	if (rate > 0.0) {
		int deltat = static_cast<int> (next_fail / (rate * modulus));
		if (deltat > 0) {
			for (int i = 0; i < deltat; i++) {
				if (time == time_max) {
					return true; //flag raised that simulation is done
				}
				systemstress[time] = next_fail * area / (static_cast<double> (deltat)) + systemstress[time - 1];
				strain[time] = strain[time - 1] + rate;
				time++;
			}
		}
	}
	else {
		strain[time] = next_fail / modulus; //modulus*strain = stress
	}

	//load all cells to fail
	for (int i = 0; i < area; i++) {
		strs[i] += next_fail;
	}
	return false;
}

bool during_slip(int& time, double& rate, double& area, double& a, double& b, double& consvarea, double& invconsv, double& weakening,double& modulus, double* strs, double* arr_strs, double* fail_strs, double* systemstress, double* strain, double* failed_this_timestep)
{
	//default is false. That is, the avalanche will stop.
	bool will_fail = false;
	double fail_amount = 0.0;
	double x; //holds the amount of stress failed on a single cell taht fails, makes code easier to read
	double oneplusconsvarea = (1 + consvarea);
	double tsum = 0;
	double strainsum = 0; //sum up the strain additions, assuming that each slipping cell contributes 1/N^2 additional strain
	double invareasq = 1 / (area*area);

	for (int i = 0; i < area; i++)
	{
		tsum += strs[i];
		failed_this_timestep[i] = (strs[i] >= fail_strs[i]); //set failed_this_timestep to 1 if that cell failed
		if (strs[i] >= fail_strs[i])
		{

			strainsum += invareasq; //strain += 1/N^2 per failed cell
			will_fail = true;
			x = (fail_strs[i] - arr_strs[i]) * wbl(a, b); //wbl(a,b) > 1 some large fraction of the time
			fail_amount += x;
			strs[i] -= oneplusconsvarea * x;
			if (fail_strs[i] == 1)
				fail_strs[i] = 1 - weakening * (1 - arr_strs[i]);
		}
	}

	//redistribute all stress back into the system
	for (int i = 0; i < area; i++)
	{
		//adiabatic (works with tau = 3/2 found)
		//strs[i] += (consvarea)*fail_amount;
		strs[i]+= (consvarea)*fail_amount + (invconsv - 1.0) * rate; //alan version, for adiabatic set to rate = 0.

	}

	//update system variables
	strain[time] = strain[time - 1] + strainsum + rate;
	systemstress[time] = tsum;

	return will_fail;
}

//For most cases, set modulus = 1e-3.
return_tuple cpp_sim(double area, double time_max, double weakening, double consv, double rate, double w, double a, double modulus, int to_write_stress, int to_write_failure_times)
{
	double b = 1 / tgamma(1 + 1 / a);
	vec<double> strs(area, -1);
	vec<double> systemstress(time_max, 0); //initialize all systemstress to be 0.
	vec<double> strain(time_max, 0);
	vec<double> arr_strs(area, -1);
	vec<double> fail_strs(area, -1);
	vec<double> failed_this_timestep(area, 0); //'failed this timestep' will have 1 on all indices that have failed

	FILE* failure_file; //get the file to write failure time output to, in case the option is specified
	std::ostringstream oss_failure_file;
	oss_failure_file << "failure_time_output" << std::scientific << std::setprecision(2) << weakening << "_k_" << std::setprecision(1) << std::scientific << a << ".txt";
	const char* failure_file_name = oss_failure_file.str().c_str();//.c_str(); //convert to str() then to c_str

	if (to_write_failure_times > 0)
	{
		fopen_s(&failure_file, failure_file_name, "w"); //open file to write
	}

	//calculate once to save time
	double consvarea = consv / (area - 1);
	double invconsv = 1.0 / consv;

	double next_fail;
	int will_fail;
	double fail_amount;
	vec<double> randm(area, -1);

	int time = 0;
	int deltat = 0; //time to fail
	double x = -1; //holds amount of stress failed on the current cell

	bool simulation_done = false; //simulation is not done until this flag is pinged.

	//set up the simulation
	sim_setup(area, w, strs, systemstress, arr_strs, fail_strs);
	time = 1;

	//main loop
	while (time < time_max)
	{
		//loads all cells by the minimum amount required to the next fail. Time iterates to be when the next slip happens.
		simulation_done = between_slips(time, time_max, rate, area,modulus, strs.data(), arr_strs.data(), fail_strs.data(), systemstress.data(),strain.data());
		//if the simulation is done, break the loop.
		if (simulation_done)
			break;
		//between_slips guarantees there will be a slip.
		will_fail = 1;

		while ((time < time_max) && will_fail)
		{
			//(int& time, double& rate, double& area, double& a, double& b, double& consvarea, double& invconsv, double& weakening,double& modulus, double* strs, double* arr_strs, double* fail_strs, double* systemstress, double* strain)
			will_fail = during_slip(time, rate, area, a, b, consvarea, invconsv, weakening, modulus, strs.data(), arr_strs.data(), fail_strs.data(), systemstress.data(),strain.data(), failed_this_timestep.data());


			//Write the cells that failed this time step to a file, if that option is given.
			if (to_write_failure_times > 0)
			{
				for (int i = 0; i < area; i++)
				{
					if (failed_this_timestep[i] == 1)
					{
						fprintf(failure_file, "%d %d\n", time, i);
					}
				}

			}

			++time;
		}
	}

	//if to_write > 0, write the system stress to a file.
	if (to_write_stress > 0)
	{
		std::ostringstream oss;
		oss << "scms_avg_eps_" << std::scientific << std::setprecision(2) << weakening << "_k_" << std::setprecision(1) << std::scientific << a << ".txt";
		std::string name = oss.str();
		helpers::write_txt(name, systemstress);
	}

	if (to_write_failure_times > 0)
	{
		fclose(failure_file);
	}

	return_tuple out = std::make_tuple(systemstress, strain, strs); //stress is vector of stresses at the end of the simulation. force is the total amount of stress in the system. strain is proportional to the amount of time
	return out;
}


int main() {

	double area = 1e6; //size of the simulation (number of cells)
	double time_max = 1e2; //amount of time to run simulation for
	double consv = 1 - 1 / sqrt(area); //set to 1-1/sqrt(N), the critical point
	double weakening = 0.0025; //weakening changes the simulation from "low mode" (mostly aperiodic) to "high mode" (mostly semiperiodic) avalanches.
	double rate = 0; //set to 0 makes adiabatic simulation
	double modulus = 1e-3; //does not do anything when rate = 0
	double w = 0.1; //amount of quenched disorder
	double a = 18; //Weibull distribution shape parameter. Higher value = system is more ordered.
	int to_write = 1; //tells the simulation to write the simulation stress output.
	int to_write_failures = 1; //tells the simulation to write to a separate file the timestep at which each cell failed. This can likely be compressed, but the current implementation leads to extremely large files.

	return_tuple q = cpp_sim(area, time_max, weakening, consv, rate, w, a, modulus, to_write, to_write_failures);

	return 0;
}