#ifndef START_H_
#define START_H_

#include<model.h>
#include<build_field.h>
#include<initialise_chunk.h>
#include<generate.h>
#include<field_summary.h>
#include<ideal_gas.h>
void start(FILE* g_out){

    /**--------------------------decompose 2D grid ----------------------------**/
    //if (ops_is_root()) {
    fprintf(g_out," Setting up initial geometry\n");
    fprintf(g_out,"\n");
    //}

    clover_time  = 0.0;
    step  = 0;
    dtold = dtinit;
    dt    = dtinit;

    //cadef viene fatto dentro build_field
    build_field();

    initialise_chunk();

    /**---------------------------Generating Chunks----------------------------**/

    fprintf(g_out,"\n");
    fprintf(g_out," Generating chunks\n");
    fprintf(g_out,"\n");

    generate();

    advect_x = true;

    /**------------------------------ideal_gas---------------------------------**/
    ideal_gas(false);
    calUpdate2D(clover_model);

    //Prime all halo data for the first step
    fields[FIELD::DENSITY0]  = 1;
    fields[FIELD::ENERGY0]   = 1;
    fields[FIELD::PRESSURE]  = 1;
    fields[FIELD::VISCOSITY] = 1;
    fields[FIELD::DENSITY1]  = 1;
    fields[FIELD::ENERGY1]   = 1;
    fields[FIELD::XVEL0]     = 1;
    fields[FIELD::YVEL0]     = 1;
    fields[FIELD::XVEL1]     = 1;
    fields[FIELD::YVEL1]     = 1;

    //update_halo qui.


    fprintf(g_out,"\n");
    fprintf(g_out," Problem initialised and generated\n");
    fprintf(g_out,"\n");

    /**----------------------------field_summary-------------------------------**/

    field_summary();


}


#endif //START_H_
