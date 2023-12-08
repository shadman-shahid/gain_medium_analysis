/*
    Copyright 2012 Lumerical Solutions, Inc. All rights reserved.
*/
#ifndef _TWODYESYSTEM_H
#define _TWODYESYSTEM_H

#include "imaterialplugin.h"
#include <vector>

class TwoDyeSystemPlugin : public IMaterialPlugin
{
public:
    TwoDyeSystemPlugin(){};
    virtual ~TwoDyeSystemPlugin(){};

    const char* name() const {return "Two Dye System";};
    const char* uniqueId() const {return "{B30704EC-9476-40AE-BB26-073291C445BD}";};
    const char** parameterNames() const { return names;};
    float calculateEx( float U, float V, float Ex, float* storage);
    float calculateEy( float U, float V, float Ey, float* storage);
    float calculateEz( float U, float V, float Ez, float* storage);
    void initialize(const double** parameters, double dt);
    void initializeStorageEx(float* storage) { initializeStorageE(0,storage); }
    void initializeStorageEy(float* storage) { initializeStorageE(1,storage); }
    void initializeStorageEz(float* storage) { initializeStorageE(2,storage); }
    size_t storageSizeE() const {return 16;};

private:
    float calculate(int i, float U, float V, float E, float *storage);
    void initializeStorageE(int axis, float*storage);

    double oneOverNormalizedWaa[3];
    double oneOverNormalizedWba[3];
    double oneOverNormalizedWb[3];
    double oneOverNormalizedWc[3];

    double aa1[3];
    double aa2[3];
    double aa4[3];
    
    double ab1[3];
    double ab2[3];
    double ab4[3];

    double ba1[3];
    double ba2[3];
    double ba4[3];
    
    double bc1[3];
    double bc2[3];
    double bc4[3];
    
    double at10[3];
    double at21[3];
    double at32[3];
    double at30[3];

    double bt10[3];
    double bt21[3];
    double bt32[3];
    double bt30[3];

    double N00[3];
    double N10[3];
    double N20[3];
    double N30[3];

    double M00[3];
    double M10[3];
    double M20[3];
    double M30[3];
    double totalElectrons_a[3];
    double totalElectrons_b[3];
    bool forceElectronConservation_1[3];
    bool forceElectronConservation_2[3];

    static const char* names[29];

    static const double hbar;
    static const double pi;
    static const double eps0;
    static const double c;



};

#endif


