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

// mod2 3D Cellular Automaton

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DRun.h>

#define ROWS 5
#define COLS 7
#define LAYERS 3

// declare CA, substate and simulation objects
struct CALModel3D* mod2;
struct CALSubstate3Db* Q;
struct CALRun3D* mod2_simulation;

// The cell's transition function
void mod2TransitionFunction(struct CALModel3D* ca, int i, int j, int k)
{
	int sum = 0, n;

	for (n=0; n<ca->sizeof_X; n++)
		sum += calGetX3Db(ca, Q, i, j, k, n);

	calSet3Db(ca, Q, i, j, k, sum%2);
}

int main()
{
	// define of the mod2 CA and mod2_simulation simulation objects
	mod2 = calCADef3D(ROWS, COLS, LAYERS, CAL_MOORE_NEIGHBORHOOD_3D, CAL_SPACE_TOROIDAL, CAL_NO_OPT);
	mod2_simulation = calRunDef3D(mod2, 1, 1, CAL_UPDATE_IMPLICIT);

	// add the Q substate to the mod2 CA
	Q = calAddSubstate3Db(mod2);

	// add transition function's elementary process
	calAddElementaryProcess3D(mod2, mod2TransitionFunction);

	// set the whole substate to 0
	calInitSubstate3Db(mod2, Q, 0);

	// set a seed at position (2, 3, 1)
	calInit3Db(mod2, Q, 2, 3, 1, 1);

	// save the Q substate to file
	calSaveSubstate3Db(mod2, Q, "./mod2_0000.txt");

	// simulation run
	calRun3D(mod2_simulation);

	// save the Q substate to file
	calSaveSubstate3Db(mod2, Q, "./mod2_LAST.txt");

	// finalize simulation and CA objects
	calRunFinalize3D(mod2_simulation);
	calFinalize3D(mod2);

	return 0;
}
