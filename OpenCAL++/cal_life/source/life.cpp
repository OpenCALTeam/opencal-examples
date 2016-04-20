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

// Conway's game of Life Cellular Automaton

#include <OpenCAL++/calModel.h>
#include <OpenCAL++/calIntConverterIO.h>
#include <OpenCAL++/calMooreNeighborhood.h>
#include <OpenCAL++/calRun.h>
#include <OpenCAL++/calSubstate.h>

#include <stdlib.h>
#include <iostream>

// declare CA, substate and simulation objects

CALSubstate<int>* Q;


// The cell's transition function





class Life_transition_function :public CALElementaryProcessFunctor{

    void run(CALModel* calModel, int* indexes)
    {

        int sum = 0, n;
        int neighborhoodSize = calModel->getNeighborhoodSize();
        for (n=1; n<neighborhoodSize; n++)
        {
            sum+= Q->getX(indexes,calModel->getCoordinates(), calModel->getDimension(), n);
        }



        if ((sum == 3) || (sum == 2 && Q->getElement(indexes, calModel->getCoordinates(), calModel->getDimension()) == 1))
        {
            Q-> setElement(indexes, calModel->getCoordinates(), calModel->getDimension(), 1);
        }
        else
           Q-> setElement(indexes,calModel->getCoordinates(), calModel->getDimension(), 0);
    }
};

int main()
{

    int* coordinates = new int [2] {8,16};
    int dimension = 2;
    CALNeighborhood* neighborhood = new CALMooreNeighborhood();

    CALModel* life = new CALModel(coordinates, dimension,neighborhood, calCommon::CAL_SPACE_TOROIDAL,calCommon::CAL_NO_OPT);

    // define of the life CA and life_simulation simulation objects
    CALRun* life_simulation = new CALRun (life, 1, 1, calCommon:: CAL_UPDATE_IMPLICIT);

    // add the Q substate to the life CA
    Q = life->addSubstate<int>();

//	// add transition function's elementary process
    life->addElementaryProcess(new Life_transition_function());

//	// set the whole substate to 0
    life->initSubstate(Q, 0);

    int indexes[5][2] = {{0,2},
                        {1,0},
                        {1,2},
                        {2,1},
                        {2,2}};

    // set a glider
    for (int i = 0; i< 5; i++)
    {
        life->init(Q, indexes[i], 1);
    }

    // save the Q substate to file
    Q->saveSubstate(coordinates,dimension, new CALIntConverterIO(), (char*) "./life_000.txt");


//	// simulation run
    life_simulation->run();



//	// save the Q substate to file
    Q->saveSubstate(coordinates,dimension, new CALIntConverterIO(), (char*) "./life_LAST.txt");

//	// finalize simulation and CA objects
    delete life_simulation;
    delete life;

	return 0;
}
