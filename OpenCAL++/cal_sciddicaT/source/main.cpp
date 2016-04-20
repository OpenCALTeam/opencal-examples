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

#include "sciddicaT.h"

#define ROWS 610
#define COLS 496
time_t start_time, end_time;

int main()
{
    time_t start_time, end_time;

    //define coordinates of cellular space
    int * coordinates = new int[2] {ROWS, COLS};
    // dimension of cellular space
    int dimension = 2;
    SciddicaTModel* sciddicaTModel = new SciddicaTModel (coordinates, dimension);

    printf ("Starting simulation...\n");
    start_time = time(NULL);
    sciddicaTModel->sciddicaTRun();
    end_time = time(NULL);
    printf ("Simulation terminated.\nElapsed time: %ld\n", end_time-start_time);

    //saving configuration
    sciddicaTModel->sciddicaTSaveConfig();
    //deallocates memory
    delete sciddicaTModel;
    return 0;
}
