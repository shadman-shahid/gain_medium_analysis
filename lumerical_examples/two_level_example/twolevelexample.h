/*
    Copyright 2014 Lumerical Solutions, Inc. All rights reserved.
*/
#ifndef _TWO_LEVEL_EXAMPLE_H
#define _TWO_LEVEL_EXAMPLE_H

#include "imaterialplugin.h"

class TwoLevelExample : public IMaterialPlugin
{
public:
    TwoLevelExample(){};
    virtual ~TwoLevelExample(){};

    const char* name() const {return "Two Level Example";};
    const char* uniqueId() const {return "{3E686172-9A28-4903-B7C4-8B374646904F}";};
    const char** parameterNames() const {return names;};
    float calculateEx( float U, float V, float Ex, float* storage);
    float calculateEy( float U, float V, float Ey, float* storage);
    float calculateEz( float U, float V, float Ez, float* storage);
    void initialize(const double** parameters, double dt);
    void initializeStorageEx(float* storage){};
    void initializeStorageEy(float* storage){};
    void initializeStorageEz(float* storage){};
    size_t storageSizeE() const {return 9;};

private:
    float calculate(int axis, float U, float V, float E, float* storage);

    float a1[3];
    float a2[3];
    float a3[3];
    float b1[3];
    float b2[3];

    static const char* names[5];

    static const double PI;
    static const double C;
    static const double EPS0;
    static const double HBAR;

};

#endif


