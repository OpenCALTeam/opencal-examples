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

#ifndef configurationPathLib_h
#define configurationPathLib_h
//---------------------------------------------------------------------------
void ConfigurationIdPath(char config_file_path[], char config_dir[]);
void ConfigurationFilePath(char config_file_path[], char const *name, char const *suffix, char file_path[]);
int  GetStepFromConfigurationFile(char config_file_path[]);
bool ConfigurationFileSavingPath(char config_file_path[], int step, char const * name, char const *suffix, char file_path[]);

#endif
