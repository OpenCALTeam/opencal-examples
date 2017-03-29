#ifndef INITIALISE_H_
#define INITIALISE_H_
#include<string>
#include<model.h>
#include<read_input.h>
#include<utils.h>

#include<unistd.h> //access()
using std::string;



void writeInputFileInternalTest(FILE* g_in, const string& in_file){
    g_in = open_FILE_orDie(in_file.c_str(),"w+");
    fprintf(g_in,"*clover\n");
    fprintf(g_in," state 1 density=0.2 energy=1.0\n");
    fprintf(g_in," state 2 density=1.0 energy=2.5 geometry=rectangle xmin=0.0 xmax=5.0 ymin=0.0 ymax=2.0\n");
    fprintf(g_in," x_cells=10\n");
    fprintf(g_in," y_cells=2\n");
    fprintf(g_in," xmin=0.0\n");
    fprintf(g_in," ymin=0.0\n");
    fprintf(g_in," xmax=10.0\n");
    fprintf(g_in," ymax=2.0\n");
    fprintf(g_in," initial_timestep=0.04\n");
    fprintf(g_in," timestep_rise=1.5\n");
    fprintf(g_in," max_timestep=0.04\n");
    fprintf(g_in," end_time=3.0\n");
    fprintf(g_in," test_problem 1\n");
    fprintf(g_in,"*endclover\n");
    fclose(g_in);
}

void initialise(CLOVER_MODEL* clover_model, FILE* &g_in, FILE* &g_out, const string& in_file, const string& out_file){

    g_out = open_FILE_orDie(out_file,"w");


    fprintf(g_out,"\n");
    fprintf(g_out,"Clover version %f\n", g_version);
    printf("Output file %s opened. All output will go there\n", out_file.c_str());
    fprintf(g_out,"\n");
    fprintf(g_out," Clover will run from the following input:-\n");
    fprintf(g_out,"\n");

    //check if file in_file EXISTS
    if( access(in_file.c_str(),F_OK) == -1)
        writeInputFileInternalTest(g_in,in_file);

     g_in = open_FILE_orDie(in_file.c_str(),"r");

     //read file line by line and write as it is in the out file
     char line[80];
     while(fgets(line, 80, g_in) != NULL)
     {
       fprintf(g_out,"%s", line);
     }

     fclose(g_in);

     fprintf(g_out,"\n");
     fprintf(g_out," Initialising and generating\n");
     fprintf(g_out,"\n");

    g_in = open_FILE_orDie(in_file.c_str(),"r");
     read_input(clover_model,g_in,g_out, in_file, out_file);
     fprintf(g_out," Starting the calculation\n");

     fclose(g_in);



}

#endif // INITIALISE_H_
