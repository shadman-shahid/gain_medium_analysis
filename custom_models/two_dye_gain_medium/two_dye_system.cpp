/*
    Copyright 2012 Lumerical Solutions, Inc. All rights reserved.
*/
#include "two_dye_system.h"
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;
/*!
    \class TwoDyeSystemPlugin

    \brief This model implements Taflove's four level, 2 electron model for two dye molecules in a medium. It additionally allows the option of forcing the population of the levels to be a constant. Also, it allows the initial level populations to be set
    which also allows for a different number of electrons to be set.
*/

const double TwoDyeSystemPlugin::hbar = 1.05457148e-34;
const double TwoDyeSystemPlugin::pi = 3.1415926535897931;
const double TwoDyeSystemPlugin::eps0 = 8.854187817e-12;
const double TwoDyeSystemPlugin::c = 2.99792458e8;

const char* TwoDyeSystemPlugin::names[29] = {   "w a",          //0
                                                "gamma aa",      //1
                                                "gamma ba",      //2
                                                "w b",          //3
                                                "gamma b",      //4
                                                "w c",          //5
                                                "gamma c",      //6
                                                "at30",          //7
                                                "at32",          //8
                                                "at21",          //9
                                                "at10",          //10
                                                "bt30",          //11
                                                "bt32",          //12
                                                "bt21",          //13
                                                "bt10",          //14
                                                "N density a",    //15
                                                "N density b",    //16
                                                "set initial populations", //17
                                                "N0(0)",            //18
                                                "N1(0)",            //19
                                                "N2(0)",            //20
                                                "N3(0)",            //21
                                                "M0(0)",            //22
                                                "M1(0)",            //23
                                                "M2(0)",            //24
                                                "M3(0)",            //25
                                                "do not enforce electron conservation a", //26
                                                "do not enforce electron conservation b", //27
                                                0};

void TwoDyeSystemPlugin::initialize(const double** parameters, double dt)
{
    //loop over 3 axes
    for(int i=0; i<3; i++){

        //basic parameters
        double wa =  parameters[0][i];
        double gaa =  parameters[1][i];
        double gba =  parameters[2][i];
        double wb =  parameters[3][i];
        double gb =  parameters[4][i];
        double wc =  parameters[5][i];
        double gc =  parameters[6][i];
        at30[i] = parameters[7][i];
        at32[i] = parameters[8][i];
        at21[i] = parameters[9][i];
        at10[i] = parameters[10][i];
        bt30[i] = parameters[11][i];
        bt32[i] = parameters[12][i];
        bt21[i] = parameters[13][i];
        bt10[i] = parameters[14][i];
        double Ndensity_a = parameters[15][i];
        double Ndensity_b = parameters[16][i];
        bool setInitialPopulations = bool(parameters[17][i]);
        if(setInitialPopulations) {
            N00[i] = parameters[18][i];
            N10[i] = parameters[19][i];
            N20[i] = parameters[20][i];
            N30[i] = parameters[21][i];
            M00[i] = parameters[22][i];
            M10[i] = parameters[23][i];
            M20[i] = parameters[24][i];
            M30[i] = parameters[25][i];
        } else {
            N00[i] = 1.;
            N10[i] = 1.;
            N20[i] = 0.;
            N30[i] = 0.;
            M00[i] = 1.;
            M10[i] = 1.;
            M20[i] = 0.;
            M30[i] = 0.;
        }
        totalElectrons_a[i] = N00[i]+N10[i]+N20[i]+N30[i];
        totalElectrons_b[i] = M00[i]+M10[i]+M20[i]+M30[i];
        forceElectronConservation_1[i] = !bool(parameters[26][i]);
        forceElectronConservation_2[i] = !bool(parameters[27][i]);

        //calculate update quantities
        aa1[i] = (2.-wa*wa*dt*dt)/(1.+0.5*gaa*dt);
        aa2[i] = -(1.-0.5*gaa*dt)/(1.+0.5*gaa*dt);
        aa4[i] = (dt*dt*Ndensity_a/eps0)/(1.+0.5*gaa*dt) *
                 6.*pi*eps0*c*c*c/(wa*wa*at21[i]);

        ab1[i] = (2.-wb*wb*dt*dt)/(1.+0.5*gb*dt);
        ab2[i] = -(1.-0.5*gb*dt)/(1.+0.5*gb*dt);
        ab4[i] = (dt*dt*Ndensity_a/eps0)/(1.+0.5*gb*dt) *
                 6.*pi*eps0*c*c*c/(wb*wb*at30[i]);

        ba1[i] = (2.-wa*wa*dt*dt)/(1.+0.5*gba*dt);
        ba2[i] = -(1.-0.5*gba*dt)/(1.+0.5*gba*dt);
        ba4[i] = (dt*dt*Ndensity_b/eps0)/(1.+0.5*gba*dt) *
                 6.*pi*eps0*c*c*c/(wa*wa*bt30[i]);

        bc1[i] = (2.-wc*wc*dt*dt)/(1.+0.5*gc*dt);
        bc2[i] = -(1.-0.5*gc*dt)/(1.+0.5*gc*dt);
        bc4[i] = (dt*dt*Ndensity_b/eps0)/(1.+0.5*gc*dt) *
                 6.*pi*eps0*c*c*c/(wc*wc*bt21[i]);

        

        //normalize the time constants to dt
        at30[i] /= dt;
        at32[i] /= dt;
        at21[i] /= dt;
        at10[i] /= dt;

        //normalize the time constants to dt
        bt30[i] /= dt;
        bt32[i] /= dt;
        bt21[i] /= dt;
        bt10[i] /= dt;

        //calculate 1/wa and 1/wb normalized for best use in updates
        double waa = wa * (hbar*Ndensity_a)/eps0;
        double wba = wa * (hbar*Ndensity_b)/eps0;
        wb *= (hbar*Ndensity_a)/eps0;
        wc *= (hbar*Ndensity_b)/eps0;    
        oneOverNormalizedWaa[i] = 1./waa;
        oneOverNormalizedWba[i] = 1./wba;    
        oneOverNormalizedWb[i] = 1./wb;
        oneOverNormalizedWc[i] = 1./wc;
        
    }
}

void TwoDyeSystemPlugin::initializeStorageE(int axis, float* storage)
{
    //initalize the electron level populations, which are the only non-zero storage
	storage[8] = float(N00[axis]);
	storage[9] = float(N10[axis]);
	storage[10] = float(N20[axis]);
	storage[11] = float(N30[axis]);
	storage[12] = float(M00[axis]);
	storage[13] = float(M10[axis]);
	storage[14] = float(M20[axis]);
	storage[15] = float(M30[axis]);
}

float TwoDyeSystemPlugin::calculate(int i, float U, float V, float Ef,float *storage)
{
    double E = double(Ef);
    double Paan = storage[0];
    double Paanm1 = storage[1];
    double Pbn = storage[2];
    double Pbnm1 = storage[3];
    
    double Pbanp1 = storage[4];
    double Pban = storage[5];
    double Pcnp1 = storage[6];
    double Pcn = storage[7];
    
    double N0 = storage[8];
    double N1 = storage[9];
    double N2 = storage[10];
    double N3 = storage[11];
    double M0 = storage[12];
    double M1 = storage[13];
    double M2 = storage[14];
    double M3 = storage[15];

    //Step 1 , 1st Dye
    double Paanp1 = aa1[i]*Paan + aa2[i]*Paanm1 + aa4[i]*(N1-N2)*E; //N2, N1 are offset in time by 1/2 time step, ignored here
    double Pbnp1 = ab1[i]*Pbn + ab2[i]*Pbnm1 + ab4[i]*(N0-N3)*E; //N3, N0 are offset in time by 1/2 time step, ignored here

    //Step 2, 1st Dye
    double Enp1 = (double(V) - Paanp1 - Pbnp1)/double(U);

    //Step 3, 2nd Dye
    double Pbanp2 = ba1[i]*Pbanp1 + ba2[i]*Pban + ba4[i]*(M0-M3)*Enp1; //N2, N1 are offset in time by 1/2 time step, ignored here
    double Pcnp2 = bc1[i]*Pcnp1 + bc2[i]*Pcn + bc4[i]*(M1-M2)*Enp1; //N2, N1 are offset in time by 1/2 time step, ignored here

    //Step 4, 2nd Dye
    double Enp2 = (double(V) - Paanp1 - Pbnp1 - Pbanp2 - Pcnp2)/double(U);
    
    //Step 5
    N3 = at30[i]*at32[i]/(2.*at30[i]*at32[i]+at30[i]*(1.-N2)+at32[i]*(1.-N0)) *
                   (N3*((N0-1.)/at30[i]+(N2-1.)/at32[i] + 2.) + E*(Pbnp1-Pbnm1)*oneOverNormalizedWb[i]);

    //Step 4
    N1 = at10[i]*at21[i]/(2.*at10[i]*at21[i]+at10[i]*N2+at21[i]*(1.-N0)) *
                   (N1*((N0-1.)/at10[i]-N2/at21[i] + 2.) + 2.*N2/at21[i]- E*(Paanp1-Paanm1)*oneOverNormalizedWaa[i]);

    //Step 5
    N2 = at21[i]*at32[i]/(2.*at21[i]*at32[i]+at21[i]*N3+at32[i]*(1.-N1)) *
                   (N2*((N1-1.)/at21[i]-N3/at32[i] + 2.) + 2.*N3/at32[i] + (Enp1+E)*(Paanp1-Paan)*oneOverNormalizedWaa[i]);

    //Step 6, calculate N0
    if(!forceElectronConservation_1[i]) {
        N0 = at10[i]*at30[i]/(2.*at10[i]*at30[i]+at10[i]*N3+at30[i]*N1) *
                   (N0*(-N1/at10[i]-N3/at30[i] + 2.) + 2.*N1/at10[i]+ 2.*N3/at30[i]-(Enp1+E)*(Pbnp1-Pbn)*oneOverNormalizedWb[i]);
    } else {
        //clamp electron probability to 0,1 for each level
        //note that we ignore here the 1/2 dt offset between levels
        N1 = max(min(1.,N1),0.);
        N2 = max(min(1.,N2),0.);
        N3 = max(min(1.,N3),0.);
        N0 = max(min(1.,totalElectrons_a[i]-N1-N2-N3),0.);
        N0 = max(min(1.,N0),0.);

        //redistribute any remaining electrons between levels 3 then 2 then 1
        double delta = totalElectrons_a[i] - N3 - N2 - N1 - N0;
        if(abs(delta) > 1e-12) {
            N3 += delta;
            N3 = max(min(1.,N3),0.);
            delta = totalElectrons_a[i] - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N2 += 2.*delta;
            N2 = max(min(1.,N2),0.);
            delta = totalElectrons_a[i] - N3 - N2 - N1 - N0;
        }
        if(abs(delta) > 1e-12) {
            N1 += delta;
            N1 = max(min(1.,N1),0.);
        }
    }

    //Step 7 Second Dye
    M3 = bt30[i]*bt32[i]/(2.*bt30[i]*bt32[i]+bt30[i]*(1.-M2)+bt32[i]*(1.-M0)) *
                   (M3*((M0-1.)/bt30[i]+(M2-1.)/bt32[i] + 2.) + Enp1*(Pbanp2-Pban)*oneOverNormalizedWba[i]);

    //Step 8 Second Dye
    M1 = bt21[i]*bt21[i]/(2.*bt21[i]*bt21[i]+bt21[i]*M2+bt21[i]*(1.-M0)) *
                   (M1*((M0-1.)/bt21[i]-M2/bt21[i] + 2.) + 2.*M2/bt21[i]- Enp1*(Pcnp2-Pcn)*oneOverNormalizedWc[i]);

    //Step 9 Second Dye
    M2 = bt21[i]*bt32[i]/(2.*bt21[i]*bt32[i]+bt21[i]*M3+bt32[i]*(1.-M1)) *
                   (M2*((M1-1.)/bt21[i]-M3/bt32[i] + 2.) + 2.*M3/bt32[i]+ (Enp2+Enp1)*(Pcnp2-Pcnp1)*oneOverNormalizedWc[i]);

    //Step 10 Second Dye, calculate M0
    if(!forceElectronConservation_2[i]) {
        M0 = bt21[i]*bt30[i]/(2.*bt21[i]*bt30[i]+bt21[i]*M3+bt30[i]*M1) *
                   (M0*(-M1/bt21[i]-M3/bt30[i] + 2.) + 2.*M1/bt21[i]+ 2.*M3/bt30[i]- (Enp2+Enp1)*(Pbanp2-Pbanp1)*oneOverNormalizedWba[i]);
    } else {
        //clamp electron probability to 0,1 for each level
        //note that we ignore here the 1/2 dt offset between levels
        M1 = max(min(1.,M1),0.);
        M2 = max(min(1.,M2),0.);
        M3 = max(min(1.,M3),0.);
        M0 = max(min(1.,totalElectrons_b[i]-M1-M2-M3),0.);
        M0 = max(min(1.,M0),0.);

        //redistribute any remaining electrons between levels 3 then 2 then 1
        double delta = totalElectrons_b[i] - M3 - M2 - M1 - M0;
        if(abs(delta) > 1e-12) {
            M3 += delta;
            M3 = max(min(1.,M3),0.);
            delta = totalElectrons_b[i] - M3 - M2 - M1 - M0;
        }
        if(abs(delta) > 1e-12) {
            M2 += 2.*delta;
            M2 = max(min(1.,M2),0.);
            delta = totalElectrons_b[i] - M3 - M2 - M1 - M0;
        }
        if(abs(delta) > 1e-12) {
            M1 += delta;
            M1 = max(min(1.,M1),0.);
        }
    }

    //update storage
    storage[0] = float(Paanp1);
    storage[1] = float(Paan);
    storage[2] = float(Pbnp1);
    storage[3] = float(Pbn);
    storage[4] = float(Pbanp2);
    storage[5] = float(Pbanp1);
    storage[6] = float(Pcnp2);
    storage[7] = float(Pcnp1);

    storage[8] = float(N0);
    storage[9] = float(N1);
    storage[10] = float(N2);
    storage[11] = float(N3);
    storage[12] = float(M0);
    storage[13] = float(M1);
    storage[14] = float(M2);
    storage[15] = float(M3);

    return float(Enp2);
}

float TwoDyeSystemPlugin::calculateEx( float U, float V, float Ex, float* storage )
{
    return calculate(0, U, V, Ex, storage);
}

float TwoDyeSystemPlugin::calculateEy( float U, float V, float Ey, float* storage )
{
    return calculate(1, U, V, Ey, storage);
}

float TwoDyeSystemPlugin::calculateEz( float U, float V, float Ez, float* storage )
{
    return calculate(2, U, V, Ez, storage);
}

MATERIAL_PLUGIN(TwoDyeSystemPlugin);



