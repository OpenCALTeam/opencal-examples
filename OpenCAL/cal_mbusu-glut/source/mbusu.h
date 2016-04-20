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

#ifndef mbusu_h
#define mbusu_h

#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>

#define YOUT 29
#define YIN 0
#define XE 159
#define XW 0
#define ZSUP 129
#define ZFONDO 0

#define COLS YOUT+1
#define ROWS XE+1
#define LAYERS ZSUP+1

#define REFRESH 100
#define VERBOSE


//cadef and rundef
extern struct CALModel3D* mbusu;
extern struct CALRun3D* mbusuSimulation;

struct mbusuSubstates {
	struct CALSubstate3Dr *teta;
	struct CALSubstate3Dr *moist_cont;
	struct CALSubstate3Dr *psi;
	struct CALSubstate3Dr *k;
	struct CALSubstate3Dr *h;
	struct CALSubstate3Dr *dqdh;
	struct CALSubstate3Dr *convergence;
	struct CALSubstate3Dr *moist_diff;
};

extern struct mbusuSubstates Q;

extern float prm_vis;


void mbusuCADef();
void mbusuSimulationInit(struct CALModel3D* ca);
void mbusuExit();



#endif
