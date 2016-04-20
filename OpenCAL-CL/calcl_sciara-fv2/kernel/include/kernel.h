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

#ifndef kernel_h
#define kernel_h

#include <OpenCAL-CL/calcl2DActive.h>
#define F(i) (i+3)
#define SZ 0
#define SLT 1
#define ST 2
#define NUMBER_OF_OUTFLOWS 8
#define MOORE_NEIGHBORS 9
#define VON_NEUMANN_NEIGHBORS 5

typedef struct{
	int x;
	int y;
} Vent;

typedef struct{

	// Parameters
	double Pclock;	//AC clock [s]
	double Pc;		//cell side
	double Pac;		//area of the cell
	double PTsol;	//temperature of solidification
	double PTvent;	//temperature of lava at vent
	// new Paramenters
	double Pr_Tsol;
	double Pr_Tvent;
	double a;		// parametro per calcolo Pr
	double b;		// parametro per calcolo Pr
	double Phc_Tsol;
	double Phc_Tvent;
	double c;		// parametro per calcolo hc
	double d;		// parametro per calcolo hc
	double Pcool;
	double Prho;	//density
	double Pepsilon;	//emissivity
	double Psigma;	//Stephen-Boltzamnn constant
	double Pcv;		//Specific heat
	double rad2;
	unsigned int emission_time;
	double effusion_duration;

} Parameters;


#endif
