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

#ifndef life3D_h
#define life3D_h

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>


#define ROWS 65
#define COLS 65
#define LAYERS 65
#define STEPS 1000


#define VERBOSE


//cadef and rundef
extern struct CALModel3D* mod2;
extern struct CALRun3D* mod2_simulation;

struct life3DSubstates {
	struct CALSubstate3Db *life;
};

extern struct life3DSubstates Q;


void life3DCADef();
void life3DSimulationInit(struct CALModel3D* ca);
void life3DExit();



#endif
