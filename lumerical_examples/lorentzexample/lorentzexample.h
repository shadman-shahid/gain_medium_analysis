/*
    Copyright 2014 Lumerical Solutions, Inc. All rights reserved.
*/
#ifndef _LORENTZEXAMPLE_H
#define _LORENTZEXAMPLE_H

#include "imaterialplugin.h"

class LorentzExamplePlugin : public IMaterialPlugin
{
public:
    LorentzExamplePlugin(){};
    virtual ~LorentzExamplePlugin(){};

    const char* name() const {return "Lorentz Example";};
    const char* uniqueId() const {return "{42068DEA-F012-42b6-BAF8-B9FABAA83DD3}";};
    const char** parameterNames() const {return names;};
    float calculateEx( float U, float V, float Ex, float* storage);
    float calculateEy( float U, float V, float Ey, float* storage);
    float calculateEz( float U, float V, float Ez, float* storage);
    void initialize(const double** parameters, double dt);
    void initializeStorageEx(float* storage){};
    void initializeStorageEy(float* storage){};
    void initializeStorageEz(float* storage){};
    size_t storageSizeE() const {return 3;};

private:
    float calculate(int axis, float U, float V, float E, float* storage);

    float a1[3];
    float a2[3];
    float a3[3];

    static const char* names[4];
};

#endif


