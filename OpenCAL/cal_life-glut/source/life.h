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

#ifndef life_h
#define life_h

#include <OpenCAL/cal2D.h>
#include <OpenCAL/cal2DRun.h>

#define ROWS 1024
#define COLS 1024
#define STATE_DEAD 0
#define STATE_ALIVE 1

//cadef and rundef
struct CellularAutomaton {
	struct CALModel2D* model;		//the cellular automaton
	struct CALSubstate2Di* Q;			//the set of call's states over the whole cellular space
	struct CALRun2D* run;		//the simulartion run
};

extern struct CellularAutomaton life;

void CADef(struct CellularAutomaton* ca);
void Init(struct CellularAutomaton* ca);
void isoExit(struct CellularAutomaton* ca);

#endif
