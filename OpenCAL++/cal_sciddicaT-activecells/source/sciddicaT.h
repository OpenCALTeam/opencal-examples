/*
 * Copyright (c) 2016 OpenCALTeam (https://github.com/OpenCALTeam),
 * University of Calabria, Italy.
 *
 * This file is part of an OpenCAL example.
 *
 * OpenCAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * OpenCAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OpenCAL. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef sciddicaTModel_h
#define sciddicaTModel_h

#include <OpenCAL++/calModel.h>
#include <OpenCAL++/calUnsafe.h>
#include <OpenCAL++/calVonNeumannNeighborhood.h>
#include <OpenCAL++/calRun.h>
#include <OpenCAL++/calRealConverterIO.h>
#include <OpenCAL++/calModelFunctor.h>

#include <stdlib.h>
#include <time.h>



#define P_R 0.5
#define P_EPSILON 0.001
#define STEPS 4000
#define DEM_PATH "./data/dem.txt"
#define SOURCE_PATH "./data/source.txt"
#define OUTPUT_PATH "./data/width_final.txt"

#define NUMBER_OF_OUTFLOWS 4
struct SciddicaTSubstates {
    CALSubstate<double> *z;
    CALSubstate<double> *h;
    CALSubstate<double> *f[NUMBER_OF_OUTFLOWS];
};
struct SciddicaTParameters {
    double epsilon;
    double r;
};

class SciddicaTModel
{
private:
//#define VERBOSE
    CALModel* sciddicaT;
    CALRun* sciddicaT_simulation;
    CALUnsafe* sciddicaT_unsafe;
    CALConverterIO* sciddicaConverterInputOutput;

    struct SciddicaTSubstates* Q;
    struct SciddicaTParameters* P;


    void sciddicaTLoadConfig();
public:

    SciddicaTModel (int* coordinates, size_t dimension);
    ~SciddicaTModel ();
    void sciddicaTSaveConfig();
    void sciddicaTRun();
    void sciddicaTExit();




};
#endif
