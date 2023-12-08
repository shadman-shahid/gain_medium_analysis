/*
    Copyright 2014 Lumerical Solutions, Inc. All rights reserved.
*/
#include "twolevelexample.h"
#include <cmath>
#include <algorithm>

/*!
    \class TwoLevelExample

    \brief This is a simple example to illustrate how multi-level, multi-electron
    models can be created with the material plugins.
*/


const double TwoLevelExample::C = 2.99792458e8;
const double TwoLevelExample::PI = 3.1415926535897931;
const double TwoLevelExample::HBAR = 1.05457148e-34;
const double TwoLevelExample::EPS0 = 8.854187817e-12;

const char* TwoLevelExample::names[5] = {"w0", "gamma", "gamma10","Nd", 0};

void TwoLevelExample::initialize(const double** parameters, double dt)
{
    for(int i=0; i<3; i++){
        double w0 = parameters[0][i];
        double delta = 0.5*parameters[1][i];
        double gamma10 = parameters[2][i];
        double Nd = parameters[3][i];
        double zeta = 6.*PI*C*C*C*gamma10/(w0*w0);
        a1[i] = float( (2.-dt*dt*w0*w0)/(delta*dt+1.) );
        a2[i] = float( (delta*dt-1.)/(delta*dt+1.) );
        a3[i] = float( zeta*Nd*dt*dt/(delta*dt+1.) );
        b1[i] = float( -2.*dt*gamma10 );
        b2[i] = float( EPS0/(Nd*HBAR*w0) );
    }
}

float TwoLevelExample::calculate(int axis, float U, float V, float E, float* storage)
{
    float Pn = storage[0];
    float Pnm1 = storage[1];
    float N1n = storage[2];
    float N1nm1 = storage[3];

    float Pnp1 = a1[axis]*Pn + a2[axis]*Pnm1 + a3[axis]*(1.f-2.f*N1n)*E;

    float N1np1 = b1[axis]*N1n + N1nm1 +b2[axis]*E*(Pnp1-Pnm1);
    N1np1 = std::min(1.f,std::max(0.f,N1np1));

    storage[0] = Pnp1;
    storage[1] = Pn;
    storage[2] = N1np1;
    storage[3] = N1n;

    storage[4] = a1[axis];
    storage[5] = a2[axis];
    storage[6] = a3[axis];
    storage[7] = b1[axis];
    storage[8] = b2[axis];

    return (V-Pnp1)/U;

}

float TwoLevelExample::calculateEx( float U, float V, float Ex, float* storage )
{
    return calculate(0, U, V, Ex, storage);
}

float TwoLevelExample::calculateEy( float U, float V, float Ey, float* storage )
{
    return calculate(1, U, V, Ey, storage);
}

float TwoLevelExample::calculateEz( float U, float V, float Ez, float* storage )
{
    return calculate(2, U, V, Ez, storage);
}

MATERIAL_PLUGIN(TwoLevelExample);

