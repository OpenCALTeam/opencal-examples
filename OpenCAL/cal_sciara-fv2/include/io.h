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

//---------------------------------------------------------------------------
#ifndef io_h
#define io_h
//---------------------------------------------------------------------------
#include "GISInfo.h"
#include "Sciara.h"
#include "configurationPathLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//---------------------------------------------------------------------------
// Autosave state variables
extern bool storing;		  //se � true avviene il salvataggio automatico
extern int storing_step;   //Ogni storing_step passi salva la configurazione
extern char storing_path[]; //percorso in cui viene salvata la configurazione
extern struct TGISInfo gis_info_Sz;
extern struct TGISInfo gis_info_generic;
extern struct TGISInfo gis_info_nodata0;
//---------------------------------------------------------------------------

int loadParameters(char const* path, Sciara* sciara);
int saveParameters(char *path, Sciara* sciara);
void printParameters(Sciara* sciara);

int loadMorphology(char* path, Sciara* sciara);
int loadVents(char* path, Sciara* sciara);
int loadEmissionRate(char *path, Sciara* sciara);

int loadAlreadyAllocatedMap(char *path, int* S, int* nS, int lx, int ly);
int loadAlreadyAllocatedMap(char *path, double* S, double* nS, int lx, int ly);

int loadConfiguration(char const *path, Sciara* sciara);
int saveConfiguration(char const *path, Sciara* sciara);


//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
