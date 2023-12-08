/*
    Copyright 2012 Lumerical Solutions, Inc. All rights reserved.
*/
#ifndef _NEW_FOUR_LEVEL_H
#define _NEW_FOUR_LEVEL_H

#include "imaterialplugin.h"
#include <vector>

class NewFourLevelPlugin : public IMaterialPlugin
{
public:
    NewFourLevelPlugin(){};
    virtual ~NewFourLevelPlugin(){};

    const char* name() const {return "New Four level model";};
    const char* uniqueId() const {return "{b611cdc7-83c4-419b-a872-b5282a75e526}";};
    const char** parameterNames() const { return names;};
    float calculateEx( float U, float V, float Ex, float* storage);
    float calculateEy( float U, float V, float Ey, float* storage);
    float calculateEz( float U, float V, float Ez, float* storage);
    void initialize(const double** parameters, double dt);
    void initializeStorageEx(float* storage) { initializeStorageE(0,storage); }
    void initializeStorageEy(float* storage) { initializeStorageE(1,storage); }
    void initializeStorageEz(float* storage) { initializeStorageE(2,storage); }
    size_t storageSizeE() const {return 14;};

private:
    float calculate(int i, float U, float V, float E, float *storage);
    void initializeStorageE(int axis, float*storage);

    double oneOverNormalizedWa[3];
    double oneOverNormalizedWb[3];

    double ca1[3];
    double ca2[3];
    double ca4[3];
    
    double cb1[3];
    double cb2[3];
    double cb4[3];

    double cn3[3];
    double cn3m1[3];
    double cn3e[3];

    double cn2[3];
    double cn2m1[3];
    double cn2e[3];

    double cn1[3];
    double cn1m1[3];
    double cn1e[3];    

    double t10[3];
    double t21[3];
    double t32[3];
    double t30[3];

    double N00[3];
    double N10[3];
    double N20[3];
    double N30[3];
    double totalElectrons[3];
    bool forceElectronConservation[3];

    static const char* names[16];

    static const double hbar;
    static const double pi;
    static const double eps0;
    static const double c;



};

#endif


