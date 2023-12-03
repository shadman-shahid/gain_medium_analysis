/*
    Copyright 2012 Lumerical Solutions, Inc. All rights reserved.
*/
#include "sixlevelexample.h"
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;
/*!
    \class SixLevelAzzamPlugin

    \brief This model implements Azzam's six level atomic system model. It additionally allows the option of
    forcing the population of the levels to be a constant. Also, it allows the initial level populations to be set
    which also allows for a different number of electrons to be set.
*/

const double SixLevelPlugin::hbar = 1.05457148e-34;
const double SixLevelPlugin::pi = 3.1415926535897931;
const double SixLevelPlugin::eps0 = 8.854187817e-12;
const double SixLevelPlugin::c = 2.99792458e8;

const char* SixLevelPlugin::names[27] = {   "w 50",          //0
                                            "gamma 50",      //1
                                            "w 40",          //2
                                            "gamma 40",      //3
                                            "w 32",          //4
                                            "gamma 32",      //5
                                            "w 31",          //6
                                            "gamma 31",      //7
                                            "t10",          //8
                                            "t20",          //9
                                            "t30",          //10
                                            "t40",          //11
                                            "t50",          //12
                                            "t31",          //13
                                            "t32",          //14
                                            "t43",          //15
                                            "t53",          //16
                                            "N density",    //17
                                            "set initial populations", //18
                                            "N0(0)",    //19
                                            "N1(0)",    //20
                                            "N2(0)",    //21
                                            "N3(0)",    //22
                                            "N4(0)",    //23
                                            "N5(0)",    //24
                                            "do not enforce electron conservation", //25
                                            0};

void SixLevelPlugin::initialize(const double** parameters, double dt)
{
    //loop over 3 axes
    for(int i=0; i<3; i++){

        //basic parameters
        double w50 =  parameters[0][i];
        double g50 =  parameters[1][i];
        double w40 =  parameters[2][i];
        double g40 =  parameters[3][i];
        double w32 =  parameters[4][i];
        double g32 =  parameters[5][i];
        double w31 =  parameters[6][i];
        double g31 =  parameters[7][i];
        t10[i] = parameters[8][i];
        t20[i] = parameters[9][i];
        t30[i] = parameters[10][i];
        t40[i] = parameters[11][i];
        t50[i] = parameters[12][i];
        t31[i] = parameters[13][i];
        t32[i] = parameters[14][i];
        t43[i] = parameters[15][i];
        t53[i] = parameters[16][i];
        double Ndensity = parameters[17][i];
        bool setInitialPopulations = bool(parameters[18][i]);
        if(setInitialPopulations) {
            N00[i] = parameters[19][i];
            N10[i] = parameters[20][i];
            N20[i] = parameters[21][i];
            N30[i] = parameters[22][i];
            N40[i] = parameters[23][i];
            N50[i] = parameters[24][i];
        } else {
            N00[i] = 1.;
            N10[i] = 0.;
            N20[i] = 0.;
            N30[i] = 0.;
            N40[i] = 0.;
            N50[i] = 0.;
        }
        totalElectrons[i] = N00[i]+N10[i]+N20[i]+N30[i]+N40[i]+N50[i];
        forceElectronConservation[i] = !bool(parameters[25][i]);

        //calculate update quantities
        ca1[i] = (2.-w50*w50*dt*dt)/(1.+0.5*g50*dt);
        ca2[i] = -(1.-0.5*g50*dt)/(1.+0.5*g50*dt);
        ca4[i] = (dt*dt*Ndensity/eps0)/(1.+0.5*g50*dt) *
                 6.*pi*eps0*c*c*c/(w50*w50*t50[i]);

        cb1[i] = (2.-w40*w40*dt*dt)/(1.+0.5*g40*dt);
        cb2[i] = -(1.-0.5*g40*dt)/(1.+0.5*g40*dt);
        cb4[i] = (dt*dt*Ndensity/eps0)/(1.+0.5*g40*dt) *
                 6.*pi*eps0*c*c*c/(w40*w40*t40[i]);

        cc1[i] = (2.-w32*w32*dt*dt)/(1.+0.5*g32*dt);
        cc2[i] = -(1.-0.5*g32*dt)/(1.+0.5*g32*dt);
        cc4[i] = (dt*dt*Ndensity/eps0)/(1.+0.5*g32*dt) *
                 6.*pi*eps0*c*c*c/(w32*w32*t32[i]);

        cd1[i] = (2.-w31*w31*dt*dt)/(1.+0.5*g31*dt);
        cd2[i] = -(1.-0.5*g31*dt)/(1.+0.5*g31*dt);
        cd4[i] = (dt*dt*Ndensity/eps0)/(1.+0.5*g31*dt) *
                 6.*pi*eps0*c*c*c/(w31*w31*t31[i]);

        //normalize the time constants to dt
        t10[i] /= dt;
        t20[i] /= dt;
        t30[i] /= dt;
        t40[i] /= dt;
        t50[i] /= dt;
        t31[i] /= dt;
        t32[i] /= dt;
        t43[i] /= dt;
        t53[i] /= dt;

        //calculate 1/w50 and 1/w40 normalized for best use in updates
        w50 *= (hbar*Ndensity)/eps0;
        w40 *= (hbar*Ndensity)/eps0;
        w32 *= (hbar*Ndensity)/eps0;
        w31 *= (hbar*Ndensity)/eps0;
        oneOverNormalizedW50[i] = 1./w50;
        oneOverNormalizedW40[i] = 1./w40;
        oneOverNormalizedW32[i] = 1./w32;
        oneOverNormalizedW31[i] = 1./w31;

    }
}

void SixLevelPlugin::initializeStorageE(int axis, float* storage)
{
    //initalize the electron level populations, which are the only non-zero storage
    storage[8] = float(N00[axis]);
    storage[9] = float(N10[axis]);
    storage[10] = float(N20[axis]);
    storage[11] = float(N30[axis]);
    storage[12] = float(N40[axis]);
    storage[13] = float(N50[axis]);
    
}

float SixLevelPlugin::calculate(int i, float U, float V, float Ef,float *storage)
{
    double E = double(Ef);
    double P50n = storage[0];
    double P50nm1 = storage[1];
    double P40n = storage[2];
    double P40nm1 = storage[3];
    double P32n = storage[4];
    double P32nm1 = storage[5];
    double P31n = storage[6];
    double P31nm1 = storage[7];
    double N0 = storage[8];
    double N1 = storage[9];
    double N2 = storage[10];
    double N3 = storage[11];
    double N4 = storage[12];
    double N5 = storage[13];

    //Step 1
    double P50np1 = ca1[i]*P50n + ca2[i]*P50nm1 + ca4[i]*(N0-N5)*E; //N0, N5 are offset in time by 1/2 time step, ignored here
    double P40np1 = cb1[i]*P40n + cb2[i]*P40nm1 + cb4[i]*(N0-N4)*E; //N0, N4 are offset in time by 1/2 time step, ignored here    
    double P32np1 = cc1[i]*P32n + cc2[i]*P32nm1 + cc4[i]*(N2-N3)*E; //N2, N3 are offset in time by 1/2 time step, ignored here
    double P31np1 = cd1[i]*P31n + cd2[i]*P31nm1 + cd4[i]*(N1-N3)*E; //N1, N3 are offset in time by 1/2 time step, ignored here


    //Step 2
    double Enp1 = (double(V) - P50np1 - P40np1 - P32np1 - P31np1)/double(U);

    //Step 3
    N5 = t53[i]*t50[i]/(2.*t50[i]*t53[i] + t50[i] + t53[i]) * 
                   (N5*(2. - 1./t50[i] - 1./t53[i]) + E*(P50np1-P50nm1)*oneOverNormalizedW50[i]);

    
    N4 = t43[i]*t40[i]/(2.*t40[i]*t43[i] + t40[i] + t43[i]) * 
                (N4*(2. - 1./t40[i] - 1./t43[i]) + E*(P40np1-P40nm1)*oneOverNormalizedW40[i]);

    //Step 4
	N2 = t20[i]/(2.*t20[i] + 1.) *
		(N2 * (2 - 1. / t20[i]) + 2.*N3 / t32[i] - E*(P32np1 - P32nm1)*oneOverNormalizedW32[i]);
	
	N1 = t10[i]/(2.*t10[i] + 1.) *
                   (N1 * (2 - 1./t10[i]) + 2.*N3/t31[i]- E*(P31np1-P31nm1)*oneOverNormalizedW31[i]);

    //Step 5
    N3 = t30[i]*t31[i]*t32[i]/(2.*t30[i]*t31[i]*t32[i] + t31[i]*t30[i] + t32[i]*t30[i] + t31[i]*t32[i]) *
                   (N3*(-1./t30[i]-1./t31[i]-1./t32[i] + 2.) + 2.*N4/t43[i] + 2.*N5/t53[i] + (Enp1+E)*(P31np1-P31n)*oneOverNormalizedW31[i] + (Enp1+E)*(P32np1-P32n)*oneOverNormalizedW32[i]);

    //Step 6, calculate N0
    if(!forceElectronConservation[i]) {
        N0 = 1./(2.) * ( 2. * N0 + 2.*N1/t10[i] + 2.*N2/t20[i] + 2.*N3/t30[i] + 2.*N4/t40[i] + 2.*N5/t50[i] - (Enp1+E)*(P40np1-P40n)*oneOverNormalizedW40[i] - (Enp1+E)*(P50np1-P50n)*oneOverNormalizedW50[i]);
    } else {
        //clamp electron probability to 0,1 for each level
        //note that we ignore here the 1/2 dt offset between levels
        N1 = max(min(1.,N1),0.);
        N2 = max(min(1.,N2),0.);
        N3 = max(min(1.,N3),0.);
        N4 = max(min(1.,N4),0.);
        N5 = max(min(1.,N5),0.);
        N0 = max(min(1.,totalElectrons[i]-N1-N2-N3-N4-N5),0.);
        N0 = max(min(1.,N0),0.);

        //redistribute any remaining electrons between levels 3 then 2 then 1
        double delta = totalElectrons[i] - N5 - N4 - N3 - N2 - N1 - N0;
        if(abs(delta) > 1e-12) {
            N5 += delta;
            N5 = max(min(1.,N5),0.);
            delta = totalElectrons[i] - N5 - N4 - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N4 += 2.*delta;
            N4 = max(min(1.,N4),0.);
            delta = totalElectrons[i] - N5 - N4 - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N3 += 2.*delta;
            N3 = max(min(1.,N3),0.);
            delta = totalElectrons[i] - N5 - N4 - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N2 += 2.*delta;
            N2 = max(min(1.,N2),0.);
            delta = totalElectrons[i] - N5 - N4 - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N1 += delta;
            N1 = max(min(1.,N1),0.);
        }
    }

    //update storage
    storage[0] = float(P50np1);
    storage[1] = float(P50n);
    storage[2] = float(P40np1);
    storage[3] = float(P40n);
    storage[4] = float(P32np1);
    storage[5] = float(P32n);
    storage[6] = float(P31np1);
    storage[7] = float(P31n);
    storage[8] = float(N0);
    storage[9] = float(N1);
    storage[10] = float(N2);
    storage[11] = float(N3);
    storage[12] = float(N4);
    storage[13] = float(N5);

    return float(Enp1);
}

float SixLevelPlugin::calculateEx( float U, float V, float Ex, float* storage )
{
    return calculate(0, U, V, Ex, storage);
}

float SixLevelPlugin::calculateEy( float U, float V, float Ey, float* storage )
{
    return calculate(1, U, V, Ey, storage);
}

float SixLevelPlugin::calculateEz( float U, float V, float Ez, float* storage )
{
    return calculate(2, U, V, Ez, storage);
}

MATERIAL_PLUGIN(SixLevelPlugin);



