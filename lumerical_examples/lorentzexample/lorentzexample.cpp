/*
    Copyright 2014 Lumerical Solutions, Inc. All rights reserved.
*/
#include "lorentzexample.h"
#include <cmath>

/*!
    \class LorentzExample

    \brief This is a simple linear Lorentz material purely for the purpose of
    demonstrating how Lumerical Material Plugins work.
*/


const char* LorentzExamplePlugin::names[4] = {"w0", "delta", "delta_eps", 0};

void LorentzExamplePlugin::initialize(const double** parameters, double dt)
{
    for(int i=0; i<3; i++){
        double w0 = parameters[0][i];
        double delta = parameters[1][i];
        double delta_eps = parameters[2][i];
        a1[i] = float( (2.-dt*dt*w0*w0)/(delta*dt+1.) );
        a2[i] = float( (delta*dt-1.)/(delta*dt+1.) );
        a3[i] = float( delta_eps*w0*w0*dt*dt/(delta*dt+1.) );
    }
}

float LorentzExamplePlugin::calculate(int axis, float U, float V, float E, float* storage)
{
    float Pn = storage[0];
    float Pnm1 = storage[1];

    float Pnp1 = a1[axis]*Pn + a2[axis]*Pnm1 + a3[axis]*E;

    storage[0] = Pnp1;
    storage[1] = Pn;

    storage[2] = a1[axis];

    return (V-Pnp1)/U;

}

float LorentzExamplePlugin::calculateEx( float U, float V, float Ex, float* storage )
{
    return calculate(0, U, V, Ex, storage);
}

float LorentzExamplePlugin::calculateEy( float U, float V, float Ey, float* storage )
{
    return calculate(1, U, V, Ey, storage);
}

float LorentzExamplePlugin::calculateEz( float U, float V, float Ez, float* storage )
{
    return calculate(2, U, V, Ez, storage);
}

MATERIAL_PLUGIN(LorentzExamplePlugin);

