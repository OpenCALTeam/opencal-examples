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

#include "mod2CA3D.h"
#include <stdio.h>
#include <stdlib.h>

//-----------------------------------------------------------------------
//	   THE mod2(oy) cellular automaton definition section
//-----------------------------------------------------------------------

//cadef and rundef
struct CALModel3D* mod2;							//the cellular automaton
struct life3DSubstates Q;							//the substate
struct CALRun3D* mod2_simulation;					//the simulartion run


//------------------------------------------------------------------------------
//					mod2 transition function
//------------------------------------------------------------------------------

//first elementary process
void mod2TransitionFunction(struct CALModel3D* ca, int i, int j, int k)
{
	int sum = 0, n;

	for (n=0; n<ca->sizeof_X; n++)
		sum += calGetX3Db(ca, Q.life, i, j, k, n);

	calSet3Db(ca, Q.life, i, j, k, sum%2);
}

//------------------------------------------------------------------------------
//					mod2 simulation functions
//------------------------------------------------------------------------------

void mod2SimulationInit(struct CALModel3D* ca)
{
	int i, j, k, state;

	//initializing substate to 0
	calInitSubstate3Db(ca, Q.life, 0);

	//initializing a specific cell
	calSet3Db(ca, Q.life, 12, 12, 12, 1);
}

CALbyte mod2SimulationStopCondition(struct CALModel3D* mod2)
{
	if (mod2_simulation->step >= STEPS)
		return CAL_TRUE;
	return CAL_FALSE;
}

//------------------------------------------------------------------------------
//					mod2 CADef and runDef
//------------------------------------------------------------------------------

void life3DCADef()
{
	//cadef and rundef
	mod2 = calCADef3D (ROWS, COLS, LAYERS, CAL_MOORE_NEIGHBORHOOD_3D, CAL_SPACE_TOROIDAL, CAL_NO_OPT);
	mod2_simulation = calRunDef3D(mod2, 1, CAL_RUN_LOOP, CAL_UPDATE_IMPLICIT);

	//add transition function's elementary processes
	calAddElementaryProcess3D(mod2, mod2TransitionFunction);

	//add substates
	Q.life = calAddSubstate3Db(mod2);

	//simulation run setup
	calRunAddInitFunc3D(mod2_simulation, mod2SimulationInit); calRunInitSimulation3D(mod2_simulation);
	calRunAddStopConditionFunc3D(mod2_simulation, mod2SimulationStopCondition);
}

//------------------------------------------------------------------------------
//					mod2 finalization function
//------------------------------------------------------------------------------

void life3DExit()
{
	//finalizations
	calRunFinalize3D(mod2_simulation);
	calFinalize3D(mod2);
}
