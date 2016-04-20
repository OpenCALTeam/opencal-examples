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

// The SciddicaT debris flows XCA transition function kernels

#include <kernel.h>

__kernel void flowsComputation(__CALCL_MODEL_2D, __global CALParameterr * Pepsilon, __global CALParameterr * Pr)
{
	calclThreadCheck2D();

	int i = calclGlobalRow();
	int j = calclGlobalColumns();

	CALbyte eliminated_cells[5] = { CAL_FALSE, CAL_FALSE, CAL_FALSE, CAL_FALSE, CAL_FALSE };
	CALbyte again;
	CALint cells_count;
	CALreal average;
	CALreal m;
	CALreal u[5];
	CALint n;
	CALreal z, h;

	CALint sizeOfX_ = calclGetNeighborhoodSize();
	CALParameterr eps = *Pepsilon;

	if (calclGet2Dr(MODEL_2D, H, i, j) <= eps)
		return;

	m = calclGet2Dr(MODEL_2D, H, i, j) - eps;
	u[0] = calclGet2Dr(MODEL_2D, Z , i, j) + eps;
	for (n = 1; n < sizeOfX_; n++) {
		z = calclGetX2Dr(MODEL_2D,Z, i, j, n);
		h = calclGetX2Dr(MODEL_2D,H, i, j, n);
		u[n] = z + h;
	}

	do {
		again = CAL_FALSE;
		average = m;
		cells_count = 0;

		for (n = 0; n < sizeOfX_; n++)
			if (!eliminated_cells[n]) {
				average += u[n];
				cells_count++;
			}

		if (cells_count != 0)
			average /= cells_count;

		for (n = 0; n < sizeOfX_; n++)
			if ((average <= u[n]) && (!eliminated_cells[n])) {
				eliminated_cells[n] = CAL_TRUE;
				again = CAL_TRUE;
			}

	} while (again);

	for (n = 1; n < sizeOfX_; n++) {
		if (eliminated_cells[n])
			calclSet2Dr(MODEL_2D, n-1, i, j, 0.0);
		else
			calclSet2Dr(MODEL_2D, n-1, i, j,(average - u[n]) * (*Pr));
	}
}

__kernel void widthUpdate(__CALCL_MODEL_2D)
{
	calclThreadCheck2D();

	int i = calclGlobalRow();
	int j = calclGlobalColumns();

	CALreal h_next;
	CALint n;

	h_next = calclGet2Dr(MODEL_2D, H, i, j);

	for (n = 1; n < calclGetNeighborhoodSize(); n++) {
		h_next += ( calclGetX2Dr(MODEL_2D, NUMBER_OF_OUTFLOWS-n, i, j, n) - calclGet2Dr(MODEL_2D, n-1, i, j) );
	}
	calclSet2Dr(MODEL_2D, H, i, j, h_next);

}

__kernel void steering(__CALCL_MODEL_2D)
{
	calclThreadCheck2D();

	CALint cols_ = calclGetColumns();
	CALint rows_ = calclGetRows();

	int i = calclGlobalRow();
	int j = calclGlobalColumns();
	int s;
	for (s = 0; s < NUMBER_OF_OUTFLOWS; ++s)
		calclInitSubstate2Dr(MODEL_2D, s, i, j, 0);
}
