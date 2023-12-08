/*
    Copyright 2012 Lumerical Solutions, Inc. All rights reserved.
*/
#ifndef _IMATERIALPLUGIN_H
#define _IMATERIALPLUGIN_H

#include <stddef.h>

/*!
    \brief The interface class definition for a material plugin for FDTD Solutions

    This pure abstract class defines the methods that must be implemented to create a
    plugin material in FDTD Solutions. A class that inherits from this interface
    class and implements all the methods can be compiled into a dynamic library that
    can be loaded into FDTD Solutions to define new materials.
*/

class IMaterialPlugin
{
public:
    virtual ~IMaterialPlugin(){};
    virtual const char* name() const = 0;
    virtual const char* uniqueId() const = 0;
    virtual const char** parameterNames() const = 0;
    virtual float calculateEx( float U, float V, float Ex, float* storage) = 0;
    virtual float calculateEy( float U, float V, float Ey, float* storage) = 0;
    virtual float calculateEz( float U, float V, float Ez, float* storage) = 0;
    virtual void initialize(const double** parameters, double dt) = 0;
    virtual void initializeStorageEx(float* storage) = 0;
    virtual void initializeStorageEy(float* storage) = 0;
    virtual void initializeStorageEz(float* storage) = 0;
    virtual size_t storageSizeE() const = 0;
};

/*!
    \brief The interface class for a magnetic material plugin

    This extends the material plugin defined above with a few more methods that need to be
    defined for a magnetic material
*/
class IMagneticMaterialPlugin : public IMaterialPlugin
{
public:
    virtual float calculateHx( float U, float V, float Ex, float* storage) = 0;
    virtual float calculateHy( float U, float V, float Ey, float* storage) = 0;
    virtual float calculateHz( float U, float V, float Ez, float* storage) = 0;
    virtual void initializeStorageHx(float* storage) = 0;
    virtual void initializeStorageHy(float* storage) = 0;
    virtual void initializeStorageHz(float* storage) = 0;
    virtual size_t storageSizeH() const = 0;
};

/*!
    \brief The interface for a factory class that creates and destroys material plugins
*/
class IMaterialPluginFactory{
public:
    virtual IMaterialPlugin* createInstance()=0;
    virtual void destroyInstance(IMaterialPlugin* i)=0;
    virtual IMagneticMaterialPlugin* toMagneticMaterialPlugin(IMaterialPlugin* p)=0;
};

/*!
    \brief A templated implementation of the IMaterialPluginFactory class

    Plugin authors do not need to write a factory class, they can just use this class. It is written
    as a template so that it can be compiled into the plugin easily. This is done in the plugin code
    usign the MATERIAL_PLUGIN(T) macro.
*/
template<class T>
class MaterialPluginFactory : public IMaterialPluginFactory
{
    IMaterialPlugin* createInstance(){return new T();}
    void destroyInstance(IMaterialPlugin* i){delete i;}
    IMagneticMaterialPlugin* toMagneticMaterialPlugin(IMaterialPlugin* p){return dynamic_cast<IMagneticMaterialPlugin*>(p);}
};

#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

//A macro to add the factory function to the plugin, instantiating the MaterialPluginFactory in the process
//All plugins should include this macro once in a source file. The argument T is the name of the user's plugin class
#define MATERIAL_PLUGIN(T) \
    extern "C" DLLEXPORT IMaterialPluginFactory* createFactoryV1(){ return new MaterialPluginFactory<T>();} \
    extern "C" DLLEXPORT void destroyFactoryV1(IMaterialPluginFactory* f){ delete f;}

#endif

