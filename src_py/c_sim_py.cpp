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
#include <string.h>
#include <omp.h> //for parallel

//C++ includes
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //for std::vector imports
#include <pybind11/numpy.h> //so we can export vectors


//local includes
#include "customs.h"
#include "Xoshiro.hpp"
#include "helpers.hpp"

namespace py = pybind11;

//custom types
typedef std::tuple<vec<double>,vec<double>,vec<double>> return_tuple; //return tuple has [force,time,strs]
//double *ranmarin(int ijkl, int N);
//double *ranwbl(int ijkl, int N, double k, double lambda);
double wbl(double k, double lambda);

//void get_u(double* u, int ijkl);

//////helpers for simulation so I can avoid copying

//Function to handle during a slip. Returns True if failure continues, otherwise returns false.
//bool during_slip(int& time, double& rate, double& area, double& a, double& b, double& consvarea, double& invconsv, double& weakening, double& modulus, double* strs, double* arr_strs, double* fail_strs, double* systemstress, double* strain, double* failed_this_timestep);



double wbl(int ind, double** us, double* cs, int* i97s, int* j97s, double cm, double cd, double k, double lambda)
{
	double uni;
	double wbl;
	int i;
	uni = us[ind][i97s[ind]] - us[ind][j97s[ind]];
	if (uni < 0.0) uni += 1.0;
	us[ind][i97s[ind]] = uni;
	if (--i97s[ind] < 0) i97s[ind] = 96;
	if (--j97s[ind] < 0) j97s[ind] = 96;
	cs[ind] -= cd;
	if (cs[ind] < 0.0) cs[ind] += cm;
	uni -= cs[ind];
	if (uni < 0.0) uni += 1.0;
	//return(uni);

	  //for(i=0;i<N;i++){
	  //printf("log uni value is %.3f\n",pow(-1.0*log(uni),1/k));
	  //printf("lambda value is %.3f\n",lambda);
	if (uni != 0.0) {
		wbl = lambda * pow(-1.0 * log(uni), 1 / k);
	}
	else {
		wbl = 10.0;
	}
	//}
	return(wbl);

}

double *ranmarin(int ijkl, int N)
{
	double c,cd,cm,u[97];
	int i97,j97,y ;
	double *uni;
	uni=(double *)malloc(N*sizeof(double));
	double output;
	int i,ii,j,jj,k,l,m ;
	double s,t ;
	int BIGPRIME;
	BIGPRIME=899999963;
	ijkl=ijkl%BIGPRIME;
	int ij=ijkl/30082;
	int kl=ijkl-30082*ij;

	i=((ij/177)%177)+2 ;
	j=(ij%177)+2 ;
	k=((kl/169)%178)+1 ;
	l=kl%169 ;
	for (ii=0;ii<97;ii++) 
	{
		s=0.0 ;
		t=0.5 ;
		for (jj=0;jj<24;jj++)
		{
			m=(((i*j)%179)*k)%179 ;
			i=j;
			j=k;
			k=m;
			l=(53*l+1)%169;
			if (((l*m)%64)>=32) s+=t;
			t*=0.5;
		}
		u[ii]=s;
	}
	c=362436.0/16777216.0;
	cd=7654321.0/16777216.0;
	cm=16777213.0/16777216.0;
	i97=96;
	j97=32;
	for (y=0;y<N;y++)
	{
		uni[y]=u[i97]-u[j97];
		if (uni[y]<0.0) uni[y]+=1.0;
		u[i97]=uni[y];
		if (--i97<0) i97=96;
		if (--j97<0) j97=96;
		c-=cd;
		if (c<0.0) c+=cm;
		uni[y]-=c;
		if (uni[y]<0.0) uni[y]+=1.0;
	}
	return uni;
}
	/*
This is a makes a Weibull-ly distributed random number. It uses the uniform rng above. 

ijkl is a seed and N is the length of the resultant random array. k and lambda are the shape parameter and mean respectively. 

It outputs uni which is an array of random doubles.

It is based on the code cited below.

*/
double *ranwbl(int ijkl, int N, double k, double lambda)
{
	double *uni;
	double *wbl;
	wbl=(double *)malloc(N*sizeof(double));
	int i;	
	uni=ranmarin(ijkl,N);
	for(i=0;i<N;i++){
		if(uni[i]!=0.0){
		wbl[i]=lambda*pow(-1.0*log(uni[i]),1/k);
		}
		else{
			wbl[i]=10.0;
		}
	}
	free(uni);
	return(wbl);

}

//make a Weibully-distributed random number using rand_num.
double wbl(double k, double lambda)
{
	double rand = helpers::ran_num<double>(0, 1);
	return lambda * std::pow(-1 * std::log(rand), 1 / k);
}

void get_u(double* u, int ijkl) {
	double output;
	int i, ii, j, jj, k, l, m;
	double s, t;
	int BIGPRIME;
	BIGPRIME = 899999963;
	ijkl %= BIGPRIME;
	int ij = ijkl / 30082;
	int kl = ijkl - 30082 * ij;

	i = ((ij / 177) % 177) + 2;
	j = (ij % 177) + 2;
	k = ((kl / 169) % 178) + 1;
	l = kl % 169;
	for (ii = 0; ii < 97; ii++) {
		s = 0.0;
		t = 0.5;
		for (jj = 0; jj < 24; jj++) {
			m = (((i * j) % 179) * k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53 * l + 1) % 169;
			if (((l * m) % 64) >= 32) s += t;
			t *= 0.5;
		}
		u[ii] = s;
	}
}

//use w = 0.1 to give a source of quenched disorder. Vary k to get different levels of system disorder.
//if an initial state is given, pull random stress values from it as an initial state.
void sim_setup(double area, double w, vec<double>& strs, vec<double>& systemstress, vec<double>& arr_strs, vec<double>& fail_strs, const char* initial_state_filename) {
	//if initial_state_filename is null, then use Alan's initial state. 
	vec<double> randm = helpers::ran_vec<double>(0, 1, static_cast<int>(area));
	vec<int>idxs = {};
	vec<double> stress_pool = {};
	int is_exists = helpers::file_exists(initial_state_filename); //check if file exists. If it does, then read it in.
	if (is_exists > 0) {
		stress_pool = helpers::read_simulation_state(initial_state_filename);
		idxs = helpers::ran_vec<int>(0, stress_pool.size(), static_cast<int>(area));
	}

	for (int i = 0; i < area; i++) {
		arr_strs[i] = w * randm[i] - w / 2;//define arrest stresses randomly
		fail_strs[i] = 1;
		if (is_exists > 0) {
			strs[i] = stress_pool[idxs[i]]; //get the stresses from the input file.
		}
		else {
			strs[i] = (std::pow(i / area, .4) * 1.56585 - .56585) * (fail_strs[i] - arr_strs[i]) + arr_strs[i];//this distribution is an approximation of the steady state distribution. It will take some time to get to the steady state. This can be improved
		}

		systemstress[0] += strs[i]; //initial system stress is sum of all strs
	}
	return;
		
}


bool between_slips(int& time, double& time_max, double& rate, double& area, double& modulus, double* strs, double* arr_strs, double* fail_strs, double* systemstress, double* strain, int is_forcerate)
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
		int deltat = static_cast<int> (next_fail / (rate * (1-is_forcerate)*modulus)); //if forcerate, then modulus is not factored in.
		if (deltat > 0) {
			for (int i = 0; i < deltat; i++) {
				if (time == time_max) {
					return true; //flag raised that simulation is done
				}
				systemstress[time] = next_fail * area / (static_cast<double> (deltat)) + systemstress[time - 1];
				strain[time] = strain[time - 1] + is_forcerate*rate; //strain only increments during avalanche during a force rate controlled avalanche.
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

bool during_slip(int& time, double& rate, double& area, double& a, double& b, double& consvarea, double& invconsv, double& weakening, double& modulus, double* strs, double* arr_strs, double* fail_strs, double* systemstress, double* strain, double* failed_this_timestep)
{
	//default is false. That is, the avalanche will stop.
	bool will_fail = false;
	double fail_amount = 0.0;
	double x; //holds the amount of stress failed on a single cell taht fails, makes code easier to read
	double oneplusconsvarea = (1 + consvarea);
	double tsum = 0;
	double strainsum = 0; //sum up the strain additions, assuming that each slipping cell contributes 1/N^2 additional strain
	double invareasq = 1 / (area*area);

	int myfail = 0; //set myfail to 0, and if it's ever greater than 1, then set will_fail to true.

	//automatically parallelize if the area is greater than 2,000,000. Speedup is minimal but necessary for extremely large simulations
	#pragma omp parallel for reduction(+:tsum, strainsum, fail_amount, myfail) if (area > 2e6)
	{
		for (int i = 0; i < area; i++)
		{
			tsum += strs[i];
			if (strs[i] >= fail_strs[i])
			{
				strainsum += invareasq; //strain += 1/N^2 per failed cell
				//will_fail = true;
				myfail++;
				x = (fail_strs[i] - arr_strs[i]) * wbl(a, b); //wbl(a,b) > 1 some large fraction of the time
				fail_amount += x;
				strs[i] -= oneplusconsvarea * x;
				if (fail_strs[i] == 1)
					fail_strs[i] = 1 - weakening * (1 - arr_strs[i]);
			}
		}
	}

	//redistribute all stress back into the system
	//todo: add in rate here so the stress increases by the system loading rate.
	for (int i = 0; i < area; i++)
	{
		strs[i]+= (consvarea)*fail_amount + (invconsv - 1.0) * rate; //alan version
		//strs[i] += (consvarea)*fail_amount + modulus * rate; //possible Jordan update -- units appear to work but gives unexpected behavior
	}

	//update system variables
	strain[time] = strain[time - 1] + strainsum + rate;
	systemstress[time] = tsum;

	return myfail > 0;
}

//For most cases, set modulus = 1e-3.
//when force rate controlled, set strain = 0 outside of an avalanche, then strain += 1/N^2 inside of an avalanche
return_tuple cpp_sim(double area, double time_max, double weakening, double consv, double rate, double w, double a, double modulus, int to_write_stress, int to_write_failure_times, int is_forcerate, const char* initial_state_filename)
{
	double b = 1 / tgamma(1 + 1 / a);
	vec<double> strs(area, -1);
	vec<double> systemstress(time_max, 0); //initialize all systemstress to be 0.
	vec<double> strain(time_max, 0);
	vec<double> arr_strs(area, -1);
	vec<double> fail_strs(area, -1);
	vec<double> failed_this_timestep(area, 0); //'failed this timestep' will have 1 on all indices that have failed during the given timestep

	FILE* failure_file; //get the file to write failure time output to, in case the option is specified
	std::ostringstream oss_failure_file;
	oss_failure_file << "failure_time_output" << std::scientific << std::setprecision(2) << weakening << "_k_" << std::setprecision(1) << std::scientific << a << ".txt";
	const char* failure_file_name = oss_failure_file.str().c_str(); //convert to str() then to c_str so we can fprintf to the file.

	if (to_write_failure_times > 0)
	{
		fopen_s(&failure_file, failure_file_name, "w"); //open file to write
	}

	//calculate once to save time.
	double consvarea = consv / (area - 1);
	double invconsv = 1.0 / consv;

	double next_fail = 100;
	int will_fail = -1;
	vec<double> randm(area, -1);

	//if is_forcerate is greater than 1, set it to 1 for the calculations.
	if (is_forcerate > 1)
		is_forcerate = 1;

	int time = 0;
	int deltat = 0; //time to fail
	double x = -1; //holds amount of stress failed on the current cell

	bool simulation_done = false; //simulation is not done until this flag is pinged.

	//set up the simulation
	sim_setup(area, w, strs, systemstress, arr_strs, fail_strs, initial_state_filename);
	time = 1;

	//main loop
	while (time < time_max)
	{
		//loads all cells by the minimum amount required to the next fail. Time iterates to be when the next slip happens.
		simulation_done = between_slips(time, time_max, rate, area, modulus, strs.data(), arr_strs.data(), fail_strs.data(), systemstress.data(), strain.data(), is_forcerate);
		//if the simulation is done, break the loop.
		if (simulation_done)
			break;
		//between_slips guarantees there will be a slip.
		will_fail = 1;

		while ((time < time_max) && will_fail)
		{
			will_fail = during_slip(time, rate, area, a, b, consvarea, invconsv, weakening, modulus, strs.data(), arr_strs.data(), fail_strs.data(), systemstress.data(), strain.data(), failed_this_timestep.data());
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


//Alan's pure C simulation.
return_tuple c_sim(double area, double time_max, double weakening, double consv, double rate, double w, double a, double modulus, int to_write) {
    //define loop variables outside of pragma omp
    int idxa;
    int idxweak;
    int thread_id = 0;
    int nloops = 0;
    int intarea = (int)area;
    int inttime_max = (int)time_max;
    //#pragma omp parallel private(thread_id, nloops)
    
    ////////////////////////////////////////////////////////////////////////////////
    int count = 0;
    double b = 1 / tgamma(1 + 1 / a);
    //thread_id = omp_get_thread_num();
    //int time_max /*time period simulated*/;
    //int area /*system size*/;
    //double consv /*conservation parameter*/;
    //double rate /*strain rate*/;
    //double modulus /*shear modulus*/;
    //double w /* arrest stress spread*/;
    //double weakening /*weakening parameter*/;
    //double a;
    //double b;//this needs to be 1/Gamma(1+1/a), it's easier to just input it manually
    int numweak;
    time_t t;
    srand((unsigned)time(&t));
    //critical c = 0.9968377223398316
    //time_max = 10000; area = 10000; consv = 0.9968377223398316; rate = 0.00; w = 0.05; modulus = 0.001;//weakening=0.01;modulus=0.001;a=2.5;b=1.12838;//define parameters
    //double* strs = (double*) malloc(intarea * sizeof(double));//stress of each cell
    std::vector<double> strs = std::vector<double>(intarea, 0);
    //double* force = (double*) malloc(inttime_max * sizeof(double));//force on system
    std::vector<double> force = std::vector<double>(inttime_max, 0);
    //double* strain = (double*)malloc(inttime_max * sizeof(double));//strain on system
    std::vector<double> strain = std::vector<double>(inttime_max, 0);
    int time;//current time step
    //double* arr_strs = (double*) malloc(intarea * sizeof(double));//arrest stresses for each cell
    std::vector<double> arr_strs = std::vector<double>(intarea);
    std::vector<double> fail_strs= std::vector<double>(intarea);
    //double* fail_strs = (double*) malloc(intarea * sizeof(double));//fail stresses for each cell
    int i;//iteration index
    double next_fail;//amount of stress to next failure
    int will_fail;//will a cell fail in the next time step
    double fail_amount;//force released by failures in a time step
    double* randm = (double*) malloc(intarea * sizeof(double));//random number array
    FILE* fp;//file to write to
    double rand = -1;


    double cd = 7654321.0 / 16777216.0;
    double cm = 16777213.0 / 16777216.0;

    double** us = (double**) malloc(intarea * sizeof(double*));
    int* i97s = (int*) malloc(intarea * sizeof(int));
    int* j97s = (int*) malloc(intarea * sizeof(int));
    double* cs = (double*) malloc(intarea * sizeof(double));
    for (int ind = 0; ind < area; ind++) us[ind] = (double*) malloc(97 * sizeof(double));




    for (int ind = 0; ind < area; ind++) {
        get_u(us[ind], ind + 47);
        cs[ind] = 362436.0 / 16777216.0;
        i97s[ind] = 96;
        j97s[ind] = 32;
    }
    //printf("%.3f\n",us[100][100]);





    randm = ranmarin(47, area);
    for (i = 0; i < area; i++) {
        //double rand=wbl(i,us,cs,i97s,j97s,cm,cd,a,b);
        arr_strs[i] = w * randm[i] - w/2;//define arrest stresses randomly
    }
    free(randm);//need to free memory because it is allocated to heap
    for (i = 0; i < area; i++) {
        fail_strs[i] = 1;//normalized
        strs[i] = (std::pow( i/ area, .4) * 1.56585 - .56585) * (fail_strs[i] - arr_strs[i]) + arr_strs[i];//this distribution is an approximation of the steady state distribution. It will take some time to get to the steady state. This can be improved
    }
    strain[0] = 0.0;
    for (i = 0; i < area; i++) {
        force[0] += strs[i];
    }//force is the sum of stresses
    time = 1;
    while (time < time_max) {
        next_fail = 99.0;//this is just a large number
        for (i = 0; i < area; i++) {
            if (next_fail > fail_strs[i] - strs[i])
                next_fail = fail_strs[i] - strs[i];//chose minimum stress to next failure
        }
        if (rate > 0.0) {
            int deltat = (int)(next_fail / (modulus * rate));
            if (deltat > 0) {
                for (i = 0; i < deltat; i++) {
                    if (time == (time_max)) {
                        break;
                    }
                    force[time] = next_fail * area / deltat + force[time - 1];
                    strain[time] = rate;//+strain[time-1];
                    time++;
                }
            }
        }
        else {
            strain[time] = next_fail * modulus + strain[time - 1];
        }
        for (i = 0; i < area; i++) {
            strs[i] += next_fail;//load to failure
        }
        will_fail = 1;//says a cell will fail
        while ((time < time_max) && will_fail) {					//this loop is the failing mechanism

            //	printf("\b\b\b\b\b\b\b\b%d%% done",(time/(time_max/100)));//not necessary but tells you how far along the simulation is
            force[time] = 0.0;
            for (i = 0; i < area; i++) {

                force[time] += strs[i];
            }//force is the sum of stresses
                        
            will_fail = 0;//default to no fail
            fail_amount = 0.0;

            //			randm=ranmarin(time,area);//use this part to use a uniform distribution
            //			for (i=0;i<area;i++){
            //				randm[i]=1+w*randm[i]-w/2;
            //			}

                        //randm=ranwbl(time,area,a,b);//use this part to use a weibull distribution
			//fail loop
            strain[time] = 0;//strain[time-1];//strain is cumulative
            for (i = 0; i < area; i++) {
                if (strs[i] >= fail_strs[i]) {
                    strain[time] += 1 / (area * area);//strain increases for each failure
                    will_fail = 1;//we now have a failure
                    rand = wbl(i, us, cs, i97s, j97s, cm, cd, a, b);
                    fail_amount += (fail_strs[i] - arr_strs[i]) * rand;//randomly distribute stress after failure based on arrest stress, add to stress lost
                    strs[i] -= ((fail_strs[i] - arr_strs[i]) * rand) * (1 + consv / (area - 1));//subtract off the amount of stress lost
                    if (fail_strs[i] == 1)	fail_strs[i] = 1 - weakening * (1 - arr_strs[i]);
                }

            }
            //	free(randm);
            for (i = 0; i < area; i++) {
                strs[i] += fail_amount * consv / (area - 1) + (1.0 / consv - 1.0) * rate;//distribute lost stress across system
            }
            time++;
        } //only fail when a cell will fail and during the time period

        for (i = 0; i < area; i++)
            fail_strs[i] = 1;//reheal
    }
	if( to_write > 0)
	{
		char name[60];
		sprintf_s(name, "stress/scms_avg_eps_%f_k_%f.txt", weakening, a);
		fopen_s(&fp, name, "w");//open the stress file, change this filename to whatever you want
		for (time = 0; time < time_max; time++) { fprintf(fp, "%f\n", force[time]); }//write the force array
		fclose(fp);
	}
    count++;
    //printf("%d\n",idxa);
    //printf("sim %d done, thread_id = %d\n", count, thread_id);

    for (i = 0; i < area; i++) free(us[i]);
    free(us);
    free(cs);
    free(i97s);
    free(j97s);
    
    return_tuple out = std::make_tuple(force, strain, strs);
    return out;
}

//do the simulation, but allow the weakening to vary over the course of the simulation.
//
// Input here now also includes a list of weakenings. The weakenings will, by default, change at even intervals. 
// That is, if there are 4 weakenings, there will be a change at 25%, 50%, and 75%.
//
return_tuple simcpp_multiple_weakenings(double area, double time_max, std::vector<double> weakenings, double consv, double rate, double w, double a, double modulus, int to_write_stress, int to_write_failure_times, int is_forcerate, const char* initial_state_filename)
{

	double b = 1 / tgamma(1 + 1 / a);
	vec<double> strs(area, -1);
	vec<double> systemstress(time_max, 0); //initialize all systemstress to be 0.
	vec<double> strain(time_max, 0);
	vec<double> arr_strs(area, -1);
	vec<double> fail_strs(area, -1);
	vec<double> failed_this_timestep(area, 0); //'failed this timestep' will have 1 on all indices that have failed during the given timestep


	//calculate once to save time.
	double consvarea = consv / (area - 1);
	double invconsv = 1.0 / consv;

	double next_fail;
	int will_fail;
	vec<double> randm(area, -1);

	//if is_forcerate is greater than 1, set it to 1 for the calculations.
	if (is_forcerate > 1)
		is_forcerate = 1;

	//if weakenings is empty, set it to default values.
	if (weakenings.empty())
	{
		std::vector<double> tmp{ 0,0.005,0.007,0.1 };
		std::cout << "Using default values of weakening = [";
		for (int i = 0; i < tmp.size(); i++)
			std::cout << tmp[i] << ", ";
		std::cout << "]. " << std::endl;
		weakenings = tmp;
	}

	int time = 0;
	int deltat = 0; //time to fail
	double x = -1; //holds amount of stress failed on the current cell

	bool simulation_done = false; //simulation is not done until this flag is pinged.

	//initialize weakening to be the first weakening in the vector
	int weakidx = 0;
	double weakening = weakenings[weakidx];
	double dtimemax = static_cast<double>(time_max);
	int to_change = static_cast<int>(dtimemax / weakenings.size());

	//get mean weakening
	double meanweak = 0;
	for (int i = 0; i < weakenings.size(); i++)
		meanweak += weakenings[i] / weakenings.size();

	FILE* failure_file; //get the file to write failure time output to, in case the option is specified
	std::ostringstream oss_failure_file;
	oss_failure_file << "failure_time_output" << std::scientific << std::setprecision(2) << meanweak << "_k_" << std::setprecision(1) << std::scientific << a << ".txt";
	const char* failure_file_name = oss_failure_file.str().c_str(); //convert to str() then to c_str so we can fprintf to the file.

	if (to_write_failure_times > 0)
	{
		fopen_s(&failure_file, failure_file_name, "w"); //open file to write
	}

	//time points where the weakening is changed
	//ex) for time_max = 1000 and size of weakenings = 3,
	//changepoints = [0,333,666,999] (add time_max to the last one to ensure that the final weakening will be in effect for the rest of the sim)
	vec<int> changepoints;
	for (int i = 1; i < weakenings.size(); i++)
		changepoints.push_back(i * to_change);

	changepoints.push_back(2 * time_max); //add time_max at the end to ensure that the simulation finishes with the final time_max.

	std::cout << "Iniital weakening w = " << weakening << "." << std::endl;

	//set up the simulation
	sim_setup(area, w, strs, systemstress, arr_strs, fail_strs, initial_state_filename);
	time = 1;

	//main loop
	while (time < time_max)
	{
		if (time >= changepoints[weakidx])
		{
			weakidx++;
			weakening = weakenings[weakidx];
			std::cout << "Weakening changed to w = " << weakening << " at time step " << time << "." << std::endl;
		}
		//loads all cells by the minimum amount required to the next fail. Time iterates to be when the next slip happens.
		simulation_done = between_slips(time, time_max, rate, area, modulus, strs.data(), arr_strs.data(), fail_strs.data(), systemstress.data(), strain.data(), is_forcerate);
		//if the simulation is done, break the loop.
		if (simulation_done)
			break;
		//between_slips guarantees there will be a slip.
		will_fail = 1;

		while ((time < time_max) && will_fail)
		{
			//(int& time, double& rate, double& area, double& a, double& b, double& consvarea, double& invconsv, double& weakening,double& modulus, double* strs, double* arr_strs, double* fail_strs, double* systemstress, double* strain)
			will_fail = during_slip(time, rate, area, a, b, consvarea, invconsv, weakening, modulus, strs.data(), arr_strs.data(), fail_strs.data(), systemstress.data(), strain.data(), failed_this_timestep.data());
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

	if (to_write_stress > 0)
	{
		std::ostringstream oss;
		oss << "scms_avg_eps_" << std::scientific << std::setprecision(2) << meanweak << "_k_" << std::setprecision(1) << std::scientific << a << ".txt";
		std::string name = oss.str();
		helpers::write_txt(name, systemstress);
	}

	return_tuple out = std::make_tuple(systemstress, strain, strs); //stress is vector of stresses at the end of the simulation. force is the total amount of stress in the system. strain is proportional to the amount of time
	return out;
}

PYBIND11_MODULE(simlib, m) {
    m.doc() = R"pbdoc(
        Mean-field simulation backend. Original C version by Alan Long 2019-2022. Ported to C++ with QoL updates and Python backend by Jordan Sickle 2022-2024.
        -----------------------
        .. currentmodule:: sim
        .. autosummary::
           :toctree: _generate
           sim_debug
    )pbdoc";

	m.def("cpp_sim", &cpp_sim, "C++ Simulation written by Jordan Sickle.", py::arg("area") = 1e4, py::arg("time_max") = 1e5, py::arg("weakening") = 0, py::arg("consv") = 0.99, py::arg("rate") = 0, py::arg("w") = 0.1, py::arg("a") = 18, py::arg("modulus") = 1e-4, py::arg("to_write_stress") = 0, py::arg("to_write_failure_times") = 0, py::arg("is_forcerate") = 0, py::arg("initial_state_filename") = "null");

	m.def("c_sim", &c_sim, "C Simulation written by Alan Long. Included for testing purposes.", py::arg("area") = 1e4, py::arg("time_max") = 1e5, py::arg("weakening") = 0, py::arg("consv") = 0.99, py::arg("rate") = 0, py::arg("w") = 0.05, py::arg("a") = 18, py::arg("modulus") = 1e-4, py::arg("to_write") = 0);
	m.def("simcpp_multiple_weakenings", &simcpp_multiple_weakenings, "C++ simulation, but allows for weakening to change over the course of the simulation.", py::arg("area") = 1e4, py::arg("time_max") = 1e5, py::arg("weakenings") = std::vector<double>(), py::arg("consv") = 0.99, py::arg("rate") = 0, py::arg("w") = 0.05, py::arg("a") = 18, py::arg("modulus") = 1e-4, py::arg("to_write_stress") = 0, py::arg("to_write_failure_times") = 0, py::arg("is_forcerate") = 0, py::arg("initial_state_filename") = "null");
}