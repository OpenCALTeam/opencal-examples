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

#ifndef GISInfo_h
#define GISInfo_h

#include <stdio.h>

#define  GIS_FILE_OK                        0
#define  GIS_FILE_GENERIC_ERROR             1
#define  GIS_FILE_TASSELATION_ERROR         2
#define  GIS_FILE_DIMENSION_ERROR           3
#define  GIS_FILE_POSITION_ERROR            4
#define  GIS_FILE_APOTHEM_ERROR             5

struct TGISInfo {
  int ncols;
  int nrows;
  double xllcorner;
  double yllcorner;
  double cell_size;
  double NODATA_value;
};

int LeggiGISInfo(TGISInfo &gis_info, FILE* f);
int VerificaGISInfo(TGISInfo gis_info, TGISInfo gis_info_morfologia);
int SalvaGISInfo(const TGISInfo &gis_info, FILE* f);
void initGISInfoNODATA0(const TGISInfo &gis_info_source, TGISInfo &gis_info_dest);

#endif
