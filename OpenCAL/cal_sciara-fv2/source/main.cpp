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

#include "io.h"
#include "Sciara.h"
#include <time.h>

#define CONFIG_PATH "./data/2006/2006"
#define SAVE_PATH "./data/2006_SAVE/2006"
#define DEM_PATH CONFIG_PATH"_000000000000_Morphology.stt"

Sciara *sciara;
int active;
time_t start_time, end_time;

int main(int argc, char** argv) {
	start_time = time(NULL);

//	int steps = atoi(argv[1]);
//	char path[1024];
//	strcpy(path, argv[2]);
//	char * demPath = (char*)malloc(sizeof(char)*(strlen(path)+strlen("_000000000000_Morphology.stt")+1));
//	strcpy(demPath, path);
//	strcat(demPath, "_000000000000_Morphology.stt\0");
//	char * outputPath = argv[3];
//	active = atoi(argv[4]);

	int steps = 1000;
	char const* outputPath = SAVE_PATH;
	active = 0;

	char path[1024] = CONFIG_PATH;
	initSciara(DEM_PATH, steps);
	int err = loadConfiguration(path, sciara);
	if (err == 0) {
		perror("cannot load configuration\n");
		exit(EXIT_FAILURE);
	}

	runSciara();

	err = saveConfiguration(outputPath, sciara);

	if (err == 0) {
		perror("cannot save configuration\n");
		exit(EXIT_FAILURE);
	}
	end_time = time(NULL);
	printf("%ld", end_time - start_time);

	return 0;

}
