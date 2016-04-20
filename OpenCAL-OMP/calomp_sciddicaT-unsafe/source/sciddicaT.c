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

// The SciddicaT further optimized CCA debris flows model

#include <OpenCAL-OMP/cal2D.h>
#include <OpenCAL-OMP/cal2DIO.h>
#include <OpenCAL-OMP/cal2DRun.h>
#include <OpenCAL-OMP/cal2DUnsafe.h>
#include <stdlib.h>
#include <time.h>

// Some definitions...
#define ROWS 610
#define COLS 496
#define P_R 0.5
#define P_EPSILON 0.001
#define STEPS 4000
#define DEM_PATH "./data/dem.txt"
#define SOURCE_PATH "./data/source.txt"
#define OUTPUT_PATH "./data/width_final.txt"
#define NUMBER_OF_OUTFLOWS 4

// declare CCA model (sciddicaT), substates (Q), parameters (P),
// and simulation object (sciddicaT_simulation)
struct CALModel2D* sciddicaT;

struct sciddicaTSubstates {
	struct CALSubstate2Dr *z;
	struct CALSubstate2Dr *h;
} Q;

struct sciddicaTParameters {
	CALParameterr epsilon;
	CALParameterr r;
} P;

struct CALRun2D* sciddicaT_simulation;


// The sciddicaT transition function
void sciddicaTFlowsComputation(struct CALModel2D* sciddicaT, int i, int j)
{
	CALbyte eliminated_cells[5]={CAL_FALSE,CAL_FALSE,CAL_FALSE,CAL_FALSE,CAL_FALSE};
	CALbyte again;
	CALint cells_count;
	CALreal average;
	CALreal m;
	CALreal u[5];
	CALint n;
	CALreal z, h;
	CALreal f;


	m = calGet2Dr(sciddicaT, Q.h, i, j) - P.epsilon;
	u[0] = calGet2Dr(sciddicaT, Q.z, i, j) + P.epsilon;
	for (n=1; n<sciddicaT->sizeof_X; n++)
	{
		z = calGetX2Dr(sciddicaT, Q.z, i, j, n);
		h = calGetX2Dr(sciddicaT, Q.h, i, j, n);
		u[n] = z + h;
	}

	//computes outflows and updates debris thickness
	do{
		again = CAL_FALSE;
		average = m;
		cells_count = 0;

		for (n=0; n<sciddicaT->sizeof_X; n++)
			if (!eliminated_cells[n]){
				average += u[n];
				cells_count++;
			}

			if (cells_count != 0)
				average /= cells_count;

			for (n=0; n<sciddicaT->sizeof_X; n++)
				if( (average<=u[n]) && (!eliminated_cells[n]) ){
					eliminated_cells[n]=CAL_TRUE;
					again=CAL_TRUE;
				}

	}while (again);

	for (n=1; n<sciddicaT->sizeof_X; n++)
		if (!eliminated_cells[n])
		{
			f = (average-u[n])*P.r;
      calAddNext2Dr(sciddicaT,Q.h,i,j,-f);
			calAddNextX2Dr(sciddicaT,Q.h,i,j,n,f);

			//adds the cell (i, j, n) to the set of active ones
			calAddActiveCellX2D(sciddicaT, i, j, n);
		}
}


void sciddicaTRemoveInactiveCells(struct CALModel2D* sciddicaT, int i, int j)
{
	if (calGet2Dr(sciddicaT, Q.h, i, j) <= P.epsilon)
		calRemoveActiveCell2D(sciddicaT,i,j);
}


void sciddicaTSimulationInit(struct CALModel2D* sciddicaT)
{
	CALreal z, h;
	CALint i, j;

	//sciddicaT parameters setting
	P.r = P_R;
	P.epsilon = P_EPSILON;

	//sciddicaT source initialization
	for (i=0; i<sciddicaT->rows; i++)
		for (j=0; j<sciddicaT->columns; j++)
		{
			h = calGet2Dr(sciddicaT, Q.h, i, j);

			if ( h > 0.0 ) {
				z = calGet2Dr(sciddicaT, Q.z, i, j);
				calSetCurrent2Dr(sciddicaT, Q.z, i, j, z-h);

				//adds the cell (i, j) to the set of active ones
        calAddActiveCell2D(sciddicaT, i, j);
			}
		}
}


int main()
{
	time_t start_time, end_time;

	// define of the sciddicaT CA and sciddicaT_simulation simulation objects
	sciddicaT = calCADef2D (ROWS, COLS, CAL_VON_NEUMANN_NEIGHBORHOOD_2D, CAL_SPACE_TOROIDAL, CAL_OPT_ACTIVE_CELLS);
	sciddicaT_simulation = calRunDef2D(sciddicaT, 1, STEPS, CAL_UPDATE_IMPLICIT);

	//put OpenCAL - OMP in unsafe state execution(to allow unsafe operation to be used)
	calSetUnsafe2D(sciddicaT);


	// add transition function's sigma_1 and sigma_2 elementary processes
	calAddElementaryProcess2D(sciddicaT, sciddicaTFlowsComputation);
	calAddElementaryProcess2D(sciddicaT, sciddicaTRemoveInactiveCells);

	// add substates
	Q.z = calAddSingleLayerSubstate2Dr(sciddicaT);
	Q.h = calAddSubstate2Dr(sciddicaT);

	// load configuration
	calLoadSubstate2Dr(sciddicaT, Q.z, DEM_PATH);
	calLoadSubstate2Dr(sciddicaT, Q.h, SOURCE_PATH);

	// simulation run
	calRunAddInitFunc2D(sciddicaT_simulation, sciddicaTSimulationInit);
	printf ("Starting simulation...\n");
	start_time = time(NULL);
	calRun2D(sciddicaT_simulation);
	end_time = time(NULL);
	printf ("Simulation terminated.\nElapsed time: %lds\n", end_time-start_time);

	// saving configuration
	calSaveSubstate2Dr(sciddicaT, Q.h, OUTPUT_PATH);

	// finalizations
	calRunFinalize2D(sciddicaT_simulation);
	calFinalize2D(sciddicaT);

	return 0;
}
