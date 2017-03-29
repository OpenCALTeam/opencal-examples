#ifndef READ_INPUT_H_
#define READ_INPUT_H_
#include<model.h>
#include<cstring>
#include <stdlib.h>
#include <stdio.h>

void set_default_values(CLOVER_MODEL* clover_model){

    test_problem = 0;
    profiler_on  = 0;
    state_max = 0;
    number_of_states = 0;


    grid.xmin = 0;
    grid.ymin = 0;
    grid.xmax = 100;
    grid.ymax = 100;

    grid.x_cells = 10;
    grid.y_cells = 10;

    end_time = 10.0;
    end_step = g_ibig;
    complete = false;

    visit_frequency=10;
    summary_frequency=10;

    dtinit = 0.1;
    dtmax = 1.0;
    dtmin = 0.0000001;
    dtrise = 1.5;
    dtc_safe = 0.7;
    dtu_safe = 0.5;
    dtv_safe = 0.5;
    dtdiv_safe = 0.7;

    use_vector_loops = true;
}


void parseToken(CLOVER_MODEL* clover_model, char* token, FILE* g_out){
    printf("token: %s \n", token);
    if(strcmp(token,"initial_timestep") == 0) {
        token = strtok(NULL, " =");
        dtinit = atof(token);
        fprintf(g_out," %20s: %e\n", "initial_timestep",dtinit);
    }
    else if(strcmp(token,"max_timestep") == 0) {
        token = strtok(NULL, " =");
        dtmax = atof(token);
        fprintf(g_out," %20s: %e\n", "max_timestep",dtmax);
    }
    else if(strcmp(token,"timestep_rise") == 0) {
        token = strtok(NULL, " =");
        dtrise = atof(token);
        fprintf(g_out," %20s: %e\n", "timestep_rise", dtrise);
    }
    else if(strcmp(token,"end_time") == 0) {
        token = strtok(NULL, " =");
        end_time = atof(token);
        fprintf(g_out," %20s: %e\n", "end_time",end_time);
    }
    else if(strcmp(token,"end_step") == 0) {
        token = strtok(NULL, " =");
        end_step = atof(token);
        fprintf(g_out," %20s: %d\n", "end_step",end_step);
    }
    else if(strcmp(token,"xmin") == 0) {
        token = strtok(NULL, " =");
        grid.xmin = atof(token);
        fprintf(g_out," %20s: %e\n", "xmin",grid.xmin);
    }
    else if(strcmp(token,"xmax") == 0) {
        token = strtok(NULL, " =");
        grid.xmax = atof(token);
        fprintf(g_out," %20s: %e\n", "xmax",grid.xmax);
    }
    else if(strcmp(token,"ymin") == 0) {
        token = strtok(NULL, " =");
        grid.ymin = atof(token);
        fprintf(g_out," %20s: %e\n", "ymin",grid.ymin);
    }
    else if(strcmp(token,"ymax") == 0) {
        token = strtok(NULL, " =");
        grid.ymax = atof(token);
        fprintf(g_out," %20s: %e\n", "ymax",grid.ymax);
    }
    else if(strcmp(token,"x_cells") == 0) {
        token = strtok(NULL, " =");
        grid.x_cells = atof(token);
        fprintf(g_out," %20s: %d\n", "x_cells",grid.x_cells);
    }
    else if(strcmp(token,"y_cells") == 0) {
        token = strtok(NULL, " =");
        grid.y_cells = atof(token);
        fprintf(g_out," %20s: %d\n", "y_cells",grid.y_cells);
    }
    else if(strcmp(token,"visit_frequency") == 0) {
        token = strtok(NULL, " =");
        visit_frequency = atoi(token);
        fprintf(g_out," %20s: %d\n", "visit_frequency",visit_frequency);
    }
    else if(strcmp(token,"summary_frequency") == 0) {
        token = strtok(NULL, " =");
        summary_frequency = atoi(token);
        fprintf(g_out," %20s: %d\n", "summary_frequency",summary_frequency);
    }
    else if(strcmp(token,"test_problem") == 0) {
        token = strtok(NULL, " =");
        test_problem = atoi(token);
        fprintf(g_out," %20s: %d\n", "test_problem",test_problem);
    }
    else if(strcmp(token,"profiler_on") == 0) {
        token = strtok(NULL, " =");
        profiler_on = atoi(token);
        fprintf(g_out," %20s: %d\n", "profiler_on",profiler_on);
    }
    else if(strcmp(token,"state") == 0) {

        fprintf(g_out,"\n");
        fprintf(g_out," Reading specification for state %d\n",number_of_states+1);
        fprintf(g_out,"\n");

        token = strtok(NULL, " =");
        states =  (state_type *) realloc(states, sizeof(state_type) * (number_of_states+1));
        states[number_of_states].xvel = 0.0;
        states[number_of_states].yvel = 0.0;


        token = strtok(NULL, " =");
        while(token) {
            if(strcmp(token,"xvel") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].xvel = atof(token);
                fprintf(g_out,"xvel: %e\n", states[number_of_states].xvel);
            }
            if(strcmp(token,"yvel") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].yvel = atof(token);
                fprintf(g_out,"yvel: %e\n", states[number_of_states].yvel);
            }

            if(strcmp(token,"xmin") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].xmin = atof(token);
                fprintf(g_out," %20s: %e\n","state xmin",states[number_of_states].xmin);
            }
            if(strcmp(token,"xmax") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].xmax = atof(token);
                fprintf(g_out," %20s: %e\n","state xmax",states[number_of_states].xmax);
            }
            if(strcmp(token,"ymin") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].ymin = atof(token);
                fprintf(g_out," %20s: %e\n","state ymin",states[number_of_states].ymin);
            }
            if(strcmp(token,"ymax") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].ymax = atof(token);
                fprintf(g_out," %20s: %e\n","state ymax",states[number_of_states].ymax);
            }
            if(strcmp(token,"density") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].density = atof(token);
                fprintf(g_out," %20s: %e\n", "state density",states[number_of_states].density);
            }
            if(strcmp(token,"energy") == 0) {
                token = strtok(NULL, " =");
                states[number_of_states].energy = atof(token);
                fprintf(g_out," %20s: %e\n", "state energy",states[number_of_states].energy);
            }
            if(strcmp(token,"geometry") == 0) {
                token = strtok(NULL, " =");
                if(strcmp(token,"rectangle") == 0) {
                    states[number_of_states].geometry = g_rect;
                    fprintf(g_out," %20s: %s\n","state geometry","rectangular");
                }
                else if(strcmp(token,"circle") == 0) {
                    states[number_of_states].geometry = g_circ;
                    fprintf(g_out," %20s: %s\n","state geometry","circular");
                }
                else if(strcmp(token,"point") == 0) {
                    states[number_of_states].geometry = g_point;
                    fprintf(g_out," %20s: %s\n","state geometry","point");
                }
            }

            token = strtok(NULL, " =");

        }

        number_of_states++;
        fprintf(g_out,"\n");
    }
}

/**
 * @brief read input file and prepare state data type plus other simulation variables
 * @param g_in :it is assumes is already opened by the caller
 * @param g_out :it is assumes is already opened by the caller
 */
void read_input(CLOVER_MODEL* clover_model, FILE* const g_in,  FILE* const g_out, const std::string& in_file, const std::string& out_file){
    constexpr const static unsigned short LINESZ = 1024;
    fprintf(g_out," Reading input file\n");

    char buff[LINESZ];
    //FILE *fin = fopen (in_file.c_str(), "r");
    if (g_in != NULL) {
        while (fgets (buff, LINESZ, g_in)) {
            char* token = strtok(buff, " =");
            while (token) {
                if(strcmp(token,"*clover\n") != 0 && strcmp(token,"*endclover\n") != 0 ) {
                    parseToken(clover_model, token, g_out);
                }
                token = strtok(NULL, " =");
            }//token
        } //fgets
    }//fopen

    if(number_of_states == 0) {
      printf("read_input, No states defined.\n");
      exit(-1);
    }

    fprintf(g_out,"\n");
    fprintf(g_out," Input read finished\n");
    fprintf(g_out,"\n");

    //field = (field_type ) xmalloc(sizeof(field_type_core));
    field.x_min = 0 +2; //+2 to account for the boundary
    field.y_min = 0 +2; //+2 to account for the boundary
    field.x_max = grid.x_cells +2; //+2 to account for the boundary
    field.y_max = grid.y_cells +2; //+2 to account for the boundary
    field.left = 0;
    field.bottom = 0;

    float dx= (grid.xmax-grid.xmin)/(float)(grid.x_cells);
    float dy= (grid.ymax-grid.ymin)/(float)(grid.y_cells);

    for(int i = 0; i < number_of_states; i++)
    {
      states[i].xmin = states[i].xmin + (dx/100.00);
      states[i].ymin = states[i].ymin + (dy/100.00);
      states[i].xmax = states[i].xmax - (dx/100.00);
      states[i].ymax = states[i].ymax - (dy/100.00);
    }


}


#endif// READ_INPUT
