/*
    Copyright 2012 Lumerical Solutions, Inc. All rights reserved.
*/
#ifndef _SIXLEVELTWOELECTRON_H
#define _SIXLEVELTWOELECTRON_H

#include "imaterialplugin.h"
#include <vector>

class SixLevelPlugin : public IMaterialPlugin
{
public:
    SixLevelPlugin(){};
    virtual ~SixLevelPlugin(){};

    const char* name() const {return "six_level_model";};
    const char* uniqueId() const {return "{BEAEEA56-DB4F-4FAA-ACEE-36FA6041BD60}";};
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

    double oneOverNormalizedW50[3];
    double oneOverNormalizedW40[3];
    double oneOverNormalizedW32[3];
    double oneOverNormalizedW31[3];

    double ca1[3];
    double ca2[3];
    double ca4[3];
    
    double cb1[3];
    double cb2[3];
    double cb4[3];

    double cc1[3];
    double cc2[3];
    double cc4[3];

    double cd1[3];
    double cd2[3];
    double cd4[3];

    double t10[3];
    double t20[3];
    double t30[3];
    double t40[3];
    double t50[3];
    double t31[3];
    double t32[3];
    double t43[3];
    double t53[3];

    double N00[3];
    double N10[3];
    double N20[3];
    double N30[3];
    double N40[3];
    double N50[3];
    double totalElectrons[3];
    bool forceElectronConservation[3];

    static const char* names[27];

    static const double hbar;
    static const double pi;
    static const double eps0;
    static const double c;



};

#endif


