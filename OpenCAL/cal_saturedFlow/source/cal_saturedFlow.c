
#include <OpenCAL-OMP/cal3D.h>
#include <OpenCAL-OMP/cal3DIO.h>
#include <OpenCAL-OMP/cal3DRun.h>
#include <OpenCAL-OMP/cal2DBuffer.h>
#include <OpenCAL/cal3DReduction.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define ROWS 40
#define COLS 80
#define LAYERS 1
#define CELL_SIZE_X 5
#define CELL_SIZE_Y 5
#define CELL_SIZE_Z 20

#define SSinitial 0.00001
#define Kinitial  0.00001

#define headFirstCol 30
#define headLastCol  25

struct CALModel3D* satured;
struct CALSubstate3Dr* head; // carico
struct CALSubstate3Dr* K; // permeabilità
struct CALSubstate3Dr* SS; // immagazinamento
struct CALSubstate3Dr* Sorgente; // sorgente
struct CALSubstate3Dr* convergence; // convergence

struct CALSubstate3Dr* Mod; // immagazinamento

struct CALRun3D* satured_simulation;

int ascii_output_time_step = 180000;				//[s] in seconds
double delta_t = 1.0;
double delta_t_cum = 0.0;
double delta_t_cum_prec = 0.0;

//double convergence = ((CELL_SIZE_X*CELL_SIZE_Y)*SSinitial)/(Kinitial*4);

void saturedTransitionFunction(struct CALModel3D* ca, int i, int j, int k)
{
	// int sum = 0, n;

	// for (n=0; n<ca->sizeof_X; n++)
	// 	sum += calGetX3Db(ca, Q, i, j, k, n);

	// calSet3Db(ca, Q, i, j, k, sum%2);

    if(j==0 || j == COLS-1){
        return;
    }
    double diffHead=0.0;
   
    double sumFlows=0.0;
    double sumVelocity=0.0;
    double tmpK = 0.0;

    for (int n=1; n<ca->sizeof_X-2; n++){
            if( ( i==0 && n==1 ) || ( i == ROWS-1 && n==4 )){
                diffHead = 0;
            }
            else
            {
                diffHead = (calGetX3Dr(ca, head, i, j, k, n)- calGet3Dr(ca, head, i, j, k));
            }
            tmpK = calGetX3Dr(ca, K, i, j, k, n);
            sumFlows += (tmpK * CELL_SIZE_Z * diffHead)-calGet3Dr(ca, Sorgente, i, j, k);
            sumVelocity += calGet3Dr(ca, K, i, j, k)*abs(diffHead);


    }

    double vol = CELL_SIZE_X*CELL_SIZE_Y*CELL_SIZE_Z;
    double ht1 = (sumFlows*delta_t)/(calGet3Dr(ca, SS, i, j, k)*vol);
    calSet3Dr(ca, head, i, j, k, ht1+calGet3Dr(ca, head, i, j, k));

    // double conv = pow(CELL_SIZE_X*CELL_SIZE_Y*0.5,2)/sumVelocity;//(CELL_SIZE_X)/(CELL_SIZE_Z*calGet3Dr(ca, K, i, j, k)*maxdiffhead);
    // //  if(i == 20 && j ==50){
    // //     printf("conv = %f\n", conv);
    // // }
    // calSet3Dr(ca, convergence, i, j, k, conv);
    
}

void saturedSimulationSteering(struct CALModel3D* satured) 
{
    // double min;
    // min = calReductionComputeMin3Dr(satured, convergence);
	// printf("min = %f\n", min);
    // if (min > 105.0)
	// 	min = 105.0;

	delta_t = 0.5;//0.95*min;
	delta_t_cum_prec = delta_t_cum;
	delta_t_cum += delta_t;

    
     printf("delta_t_cum_prec = %f\n", delta_t_cum_prec);
}

CALbyte saturedSimulationStopCondition(struct CALModel3D* satured)
{

        // printf("ascii_output_time_step = %d\n", ascii_output_time_step);
        if (delta_t_cum >= ascii_output_time_step && delta_t_cum_prec <= ascii_output_time_step)
	    {
                //printf("ascii_output_time_step ===================== %d\n", ascii_output_time_step);
                return CAL_TRUE;
        }

        return CAL_FALSE;

        
}




void saturedInit(struct CALModel3D* ca)
{
    for(int i = 0; i< ROWS;i++)
        for(int j = 0; j< COLS;j++)
            for(int k = 0; k< LAYERS;k++)
            {        
              calSet3Dr(satured, head, i, j, k, headFirstCol);
            }
	
    for(int i = 0; i< ROWS;i++)
        for(int j = 0; j< COLS;j++)
            for(int k = 0; k< LAYERS;k++)
            {
                calSet3Dr(satured, K, i, j, k, Kinitial);
                calSet3Dr(satured, SS, i, j, k, SSinitial);
            }

    for(int i = 0; i< ROWS; i++)
        for(int k = 0; k< LAYERS;k++)
        {    
            calSet3Dr(satured, head, i, COLS-1, k, headLastCol);
        }

    for(int i = 0; i< ROWS;i++)
        for(int j = 0; j< COLS;j++)
            for(int k = 0; k< LAYERS;k++)
            {    
                calSet3Dr(satured, Sorgente, i, j, k, 0);
            }
	calSet3Dr(satured, Sorgente, 20, 50, 0, 0.001);
}


// CALreal* lastModFlow;
// #define STRLEN 256
// void calfLoadMatrix2Dr(CALreal* M, int rows, int columns, FILE* f)
// {
//   char str[STRLEN];
//   int i, j;

//   for (i=0; i<rows; i++)
//     for (j=0; j<columns; j++){
//       fscanf(f, "%s", str);
//       calSetMatrixElement(M, columns, i, j, atof(str));
//     }
// }

// CALbyte calLoadMatrix2Dr(CALreal* M, int rows, int columns, char* path)
// {
//   FILE *f = NULL;
//   f = fopen(path, "r");

//   if ( !f )
//     return CAL_FALSE;

//   calfLoadMatrix2Dr(M, rows, columns, f);

//   fclose(f);

//   return CAL_TRUE;
// }

int main(){


    // define of the satured CA and satured_simulation simulation objects
	satured = calCADef3D(ROWS, COLS, LAYERS, CAL_VON_NEUMANN_NEIGHBORHOOD_3D, CAL_SPACE_FLAT, CAL_NO_OPT);
	satured_simulation = calRunDef3D(satured, 1, 180000, CAL_UPDATE_IMPLICIT);

	// add the Q substate to the satured CA
	head = calAddSubstate3Dr(satured);
    K = calAddSubstate3Dr(satured);
    SS = calAddSubstate3Dr(satured);
    Sorgente = calAddSubstate3Dr(satured);
    convergence = calAddSubstate3Dr(satured);

	// add transition function's elementary process
	calAddElementaryProcess3D(satured, saturedTransitionFunction);

    // int rowsModFlow = 320;
    // int colsModFlow = 10;
    // lastModFlow = (CALreal*) malloc(sizeof(CALreal)*(rowsModFlow*colsModFlow));

    // calLoadMatrix2Dr(lastModFlow, rowsModFlow, colsModFlow, "./LastModFlow.txt");


    // CALreal* lastModFlowCorretta = (CALreal*) malloc(sizeof(CALreal)*(ROWS*COLS));

    // int cont = 0;
    // for(int i = 0; i< rowsModFlow;i++){
    //    for(int j = 0; j< colsModFlow;j++){
    //         lastModFlowCorretta[ cont ] =  lastModFlow[ i*colsModFlow+j ];
    //         cont++;
    //         //printf("%d\t%d\t%f\n",i+1,j+1,lastModFlow[i*colsModFlow+j]);
    //    }
            
    //     //printf(" \n");
    // }

    // for(int i = 0; i< ROWS;i++){
    //    for(int j = 0; j< COLS;j++)
    //         printf("%d\t%d\t%f\n",j+1,i+1,lastModFlowCorretta[i*COLS+j]);
    //     //printf(" \n");
    // }

    // for(int i = 0; i< ROWS;i++){
    //    for(int j = 0; j< COLS;j++){
    //        printf("%f ",lastModFlowCorretta[i*COLS+j]);
    //    }
    //    printf(" \n");
           
    // }


   
    // for(int i = 0; i< ROWS;i++){
    //    for(int j = 0; j< COLS;j++){
    //        printf("%f ",head->current[i*COLS+j]);
    //    }
    //    printf(" \n");
    // }
    
    
    calRunAddInitFunc3D(satured_simulation, saturedInit);
    calRunInitSimulation3D(satured_simulation);
    //calRunAddSteeringFunc3D(satured_simulation, saturedSimulationSteering);
	calRunAddStopConditionFunc3D(satured_simulation, saturedSimulationStopCondition);

	// save the Q substate to file
	// calSaveSubstate3Dr(satured, head, "./satured_head_0000.txt");
    // calSaveSubstate3Dr(satured, K, "./satured_K_0000.txt");
    // calSaveSubstate3Dr(satured, SS, "./satured_SS_0000.txt");
    
    // printf(" convergence = %f\n", convergence);
	// // simulation run
	calRun3D(satured_simulation);

    calSaveSubstate3Dr(satured, head, "./satured_head_LAST.txt");
    // calSaveSubstate3Dr(satured, K, "./satured_K_LAST.txt");
    // calSaveSubstate3Dr(satured, SS, "./satured_SS_LAST.txt");
    calSaveSubstate3Dr(satured, convergence, "./satured_convergence_LAST.txt");

    for(int i = 0; i< ROWS;i++){
       for(int j = 0; j< COLS;j++)
            printf("%d\t%d\t%f\n",j+1,i+1,head->current[i*COLS+j]);
        //printf(" \n");
    }

    // save the Q substate to file
	//calSaveSubstate3Db(satured, Q, "./satured_LAST.txt");

	// finalize simulation and CA objects
	calRunFinalize3D(satured_simulation);
	calFinalize3D(satured);

}