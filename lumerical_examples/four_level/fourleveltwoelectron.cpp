/*
    Copyright 2012 Lumerical Solutions, Inc. All rights reserved.
*/
#include "fourleveltwoelectron.h"
#include <fstream>
#include <cmath>

using namespace std;
/*!
    \class FourLevelTwoElectrongTaflovePlugin

    \brief This model implements Taflove's four level, 2 electron model. It additionally allows the option of
    forcing the population of the levels to be a constant. Also, it allows the initial level populations to be set
    which also allows for a different number of electrons to be set.
*/

const double FourLevelTwoElectronPlugin::hbar = 1.05457148e-34;
const double FourLevelTwoElectronPlugin::pi = 3.1415926535897931;
const double FourLevelTwoElectronPlugin::eps0 = 8.854187817e-12;
const double FourLevelTwoElectronPlugin::c = 2.99792458e8;

const char* FourLevelTwoElectronPlugin::names[16] = {"w a",          //0
                                                     "gamma a",      //1
                                                     "w b",          //2
                                                     "gamma b",      //3
                                                     "t30",          //4
                                                     "t32",          //5
                                                     "t21",          //6
                                                     "t10",          //7
                                                     "N density",    //8
                                                     "set initial populations", //9
                                                     "N0(0)",    //10
                                                     "N1(0)",    //11
                                                     "N2(0)",    //12
                                                     "N3(0)",    //13
                                                     "do not enforce electron conservation", //14
                                                     0};

void FourLevelTwoElectronPlugin::initialize(const double** parameters, double dt)
{
    //loop over 3 axes
    for(int i=0; i<3; i++){

        //basic parameters
        double wa =  parameters[0][i];
        double ga =  parameters[1][i];
        double wb =  parameters[2][i];
        double gb =  parameters[3][i];
        t30[i] = parameters[4][i];
        t32[i] = parameters[5][i];
        t21[i] = parameters[6][i];
        t10[i] = parameters[7][i];
        double Ndensity = parameters[8][i];
        bool setInitialPopulations = bool(parameters[9][i]);
        if(setInitialPopulations) {
            N00[i] = parameters[10][i];
            N10[i] = parameters[11][i];
            N20[i] = parameters[12][i];
            N30[i] = parameters[13][i];
        } else {
            N00[i] = 1.;
            N10[i] = 1.;
            N20[i] = 0.;
            N30[i] = 0.;
        }
        totalElectrons[i] = N00[i]+N10[i]+N20[i]+N30[i];
        forceElectronConservation[i] = !bool(parameters[14][i]);

        //calculate update quantities
        ca1[i] = (2.-wa*wa*dt*dt)/(1.+0.5*ga*dt);
        ca2[i] = -(1.-0.5*ga*dt)/(1.+0.5*ga*dt);
        ca4[i] = (dt*dt*Ndensity/eps0)/(1.+0.5*ga*dt) *
                 6.*pi*eps0*c*c*c/(wa*wa*t21[i]);

        cb1[i] = (2.-wb*wb*dt*dt)/(1.+0.5*gb*dt);
        cb2[i] = -(1.-0.5*gb*dt)/(1.+0.5*gb*dt);
        cb4[i] = (dt*dt*Ndensity/eps0)/(1.+0.5*gb*dt) *
                 6.*pi*eps0*c*c*c/(wb*wb*t30[i]);

        //normalize the time constants to dt
        t30[i] /= dt;
        t32[i] /= dt;
        t21[i] /= dt;
        t10[i] /= dt;

        //calculate 1/wa and 1/wb normalized for best use in updates
        wa *= (hbar*Ndensity)/eps0;
        wb *= (hbar*Ndensity)/eps0;
        oneOverNormalizedWa[i] = 1./wa;
        oneOverNormalizedWb[i] = 1./wb;

    }
}

void FourLevelTwoElectronPlugin::initializeStorageE(int axis, float* storage)
{
    //initalize the electron level populations, which are the only non-zero storage
    storage[4] = float(N00[axis]);
    storage[5] = float(N10[axis]);
    storage[6] = float(N20[axis]);
    storage[7] = float(N30[axis]);
}

float FourLevelTwoElectronPlugin::calculate(int i, float U, float V, float Ef,float *storage)
{
    double E = double(Ef);
    double Pan = storage[0];
    double Panm1 = storage[1];
    double Pbn = storage[2];
    double Pbnm1 = storage[3];
    double N0 = storage[4];
    double N1 = storage[5];
    double N2 = storage[6];
    double N3 = storage[7];

    //Step 1
    double Panp1 = ca1[i]*Pan + ca2[i]*Panm1 + ca4[i]*(N1-N2)*E; //N2, N1 are offset in time by 1/2 time step, ignored here
    double Pbnp1 = cb1[i]*Pbn + cb2[i]*Pbnm1 + cb4[i]*(N0-N3)*E; //N3, N0 are offset in time by 1/2 time step, ignored here

    //Step 2
    double Enp1 = (double(V) - Panp1 - Pbnp1)/double(U);

    //Step 3
    N3 = t30[i]*t32[i]/(2.*t30[i]*t32[i]+t30[i]*(1.-N2)+t32[i]*(1.-N0)) *
                   (N3*((N0-1.)/t30[i]+(N2-1.)/t32[i] + 2.) + E*(Pbnp1-Pbnm1)*oneOverNormalizedWb[i]);

    //Step 4
    N1 = t10[i]*t21[i]/(2.*t10[i]*t21[i]+t10[i]*N2+t21[i]*(1.-N0)) *
                   (N1*((N0-1.)/t10[i]-N2/t21[i] + 2.) + 2.*N2/t21[i]- E*(Panp1-Panm1)*oneOverNormalizedWa[i]);

    //Step 5
    N2 = t21[i]*t32[i]/(2.*t21[i]*t32[i]+t21[i]*N3+t32[i]*(1.-N1)) *
                   (N2*((N1-1.)/t21[i]-N3/t32[i] + 2.) + 2.*N3/t32[i]+ (Enp1+E)*(Panp1-Pan)*oneOverNormalizedWa[i]);

    //Step 6, calculate N0
    if(!forceElectronConservation[i]) {
        N0 = t10[i]*t30[i]/(2.*t10[i]*t30[i]+t10[i]*N3+t30[i]*N1) *
                   (N0*(-N1/t10[i]-N3/t30[i] + 2.) + 2.*N1/t10[i]+ 2.*N3/t30[i]-(Enp1+E)*(Pbnp1-Pbn)*oneOverNormalizedWb[i]);
    } else {
        //clamp electron probability to 0,1 for each level
        //note that we ignore here the 1/2 dt offset between levels
        N1 = max(min(1.,N1),0.);
        N2 = max(min(1.,N2),0.);
        N3 = max(min(1.,N3),0.);
        N0 = max(min(1.,totalElectrons[i]-N1-N2-N3),0.);
        N0 = max(min(1.,N0),0.);

        //redistribute any remaining electrons between levels 3 then 2 then 1
        double delta = totalElectrons[i] - N3 - N2 - N1 - N0;
        if(abs(delta) > 1e-12) {
            N3 += delta;
            N3 = max(min(1.,N3),0.);
            delta = totalElectrons[i] - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N2 += 2.*delta;
            N2 = max(min(1.,N2),0.);
            delta = totalElectrons[i] - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N1 += delta;
            N1 = max(min(1.,N1),0.);
        }
    }

    //update storage
    storage[0] = float(Panp1);
    storage[1] = float(Pan);
    storage[2] = float(Pbnp1);
    storage[3] = float(Pbn);
    storage[4] = float(N0);
    storage[5] = float(N1);
    storage[6] = float(N2);
    storage[7] = float(N3);

    return float(Enp1);
}

float FourLevelTwoElectronPlugin::calculateEx( float U, float V, float Ex, float* storage )
{
    return calculate(0, U, V, Ex, storage);
}

float FourLevelTwoElectronPlugin::calculateEy( float U, float V, float Ey, float* storage )
{
    return calculate(1, U, V, Ey, storage);
}

float FourLevelTwoElectronPlugin::calculateEz( float U, float V, float Ez, float* storage )
{
    return calculate(2, U, V, Ez, storage);
}

MATERIAL_PLUGIN(FourLevelTwoElectronPlugin);



