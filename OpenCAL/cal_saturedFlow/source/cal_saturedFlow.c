
#include <OpenCAL/cal2D.h>
#include <OpenCAL/cal2DIO.h>
#include <OpenCAL/cal2DRun.h>
#include <OpenCAL/cal2DBuffer.h>
#include <OpenCAL/cal2DReduction.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define ROWS 100//40
#define COLS 100//80
#define LAYERS 1
#define CELL_SIZE_X 10 //5
#define CELL_SIZE_Y 10 //5
#define SPESSORE 50 //20

#define Syinitial 0.1 //0.0002
#define Kinitial  0.0000125 //0.00001

#define headFixed 50 //30
#define headCalculated 50 //30

#define SogliaMinima 25

struct CALModel2D* satured;
struct CALSubstate2Dr* head; // carico
struct CALSubstate2Dr* K; // permeabilità
struct CALSubstate2Dr* Sy; // immagazinamento
struct CALSubstate2Dr* Sorgente; // sorgente
struct CALSubstate2Dr* convergence; // convergence

struct CALSubstate2Dr* Mod; // immagazinamento

struct CALRun2D* satured_simulation;

float delta_t_ = 4000;//0.5;

//double convergence = ((CELL_SIZE_X*CELL_SIZE_Y)*Syinitial)/(Kinitial*4);

void saturedTransitionFunction(struct CALModel2D* ca, int i, int j)
{
    if(j==0 || j ==COLS-1 || i ==ROWS-1  || i ==0){
        return;
    }
    double diffHead=0.0;
   
    double sumFlows=0.0;
    double tmpT = 0.0;

    for (int n=1; n<ca->sizeof_X; n++){

            diffHead = (calGetX2Dr(ca, head, i, j, n)- calGet2Dr(ca, head, i, j));

            tmpT = calGetX2Dr(ca, K, i, j, n)*SPESSORE;

            sumFlows += (tmpT * diffHead);
    }
    sumFlows+=-calGet2Dr(ca, Sorgente, i, j);

    double area = CELL_SIZE_X*CELL_SIZE_Y;
    double ht1 = (sumFlows*delta_t_)/(calGet2Dr(ca, Sy, i, j)*area);
    calSet2Dr(ca, head, i, j, ht1+calGet2Dr(ca, head, i, j));
    
}



void saturedInit(struct CALModel2D* ca)
{
for(int i = 0; i< ROWS;i++)
        for(int j = 0; j< COLS;j++)
            {        
              if( i<25 || i> 75 || j < 25 || j > 75) 
                calSet2Dr(satured, head, i, j, headFixed);
              else
                calSet2Dr(satured, head, i, j, headCalculated);
            }
	
    for(int i = 0; i< ROWS;i++)
        for(int j = 0; j< COLS;j++)
            {
                calSet2Dr(satured, K, i, j, Kinitial);
                calSet2Dr(satured, Sy, i, j, Syinitial);
            }


    for(int i = 0; i< ROWS;i++)
        for(int j = 0; j< COLS;j++)
            {    
                calSet2Dr(satured, Sorgente, i, j, 0);
            }
	calSet2Dr(satured, Sorgente, 49, 49, 0.001);
}


CALreal* lastModFlow;
#define STRLEN 256
void calfLoadMatrix2Dr(CALreal* M, int rows, int columns, FILE* f)
{
  char str[STRLEN];
  int i, j;

  for (i=0; i<rows; i++)
    for (j=0; j<columns; j++){
      fscanf(f, "%s", str);
      calSetMatrixElement(M, columns, i, j, atof(str));
    }
}

CALbyte calLoadMatrix2Dr(CALreal* M, int rows, int columns, char* path)
{
  FILE *f = NULL;
  f = fopen(path, "r");

  if ( !f )
    return CAL_FALSE;

  calfLoadMatrix2Dr(M, rows, columns, f);

  fclose(f);

  return CAL_TRUE;
}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

int main(){


    // define of the satured CA and satured_simulation simulation objects
	satured = calCADef2D(ROWS, COLS, CAL_VON_NEUMANN_NEIGHBORHOOD_2D, CAL_SPACE_FLAT, CAL_NO_OPT);
	satured_simulation = calRunDef2D(satured, 1, 260, CAL_UPDATE_IMPLICIT);

	// add the Q substate to the satured CA
	head = calAddSubstate2Dr(satured);
    K = calAddSubstate2Dr(satured);
    Sy = calAddSubstate2Dr(satured);
    Sorgente = calAddSubstate2Dr(satured);
    convergence = calAddSubstate2Dr(satured);

	// add transition function's elementary process
	calAddElementaryProcess2D(satured, saturedTransitionFunction);

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
            
    // //     //printf(" \n");
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
    
    
    calRunAddInitFunc2D(satured_simulation, saturedInit);
    calRunInitSimulation2D(satured_simulation);
    //calRunAddSteeringFunc2D(satured_simulation, saturedSimulationSteering);
	//calRunAddStopConditionFunc2D(satured_simulation, saturedSimulationStopCondition);

	// save the Q substate to file
	// calSaveSubstate2Dr(satured, head, "./satured_head_0000.txt");
    // calSaveSubstate2Dr(satured, K, "./satured_K_0000.txt");
    // calSaveSubstate2Dr(satured, SS, "./satured_SS_0000.txt");
    
    // printf(" convergence = %f\n", convergence);
	// // simulation run
	calRun2D(satured_simulation);

    calSaveSubstate2Dr(satured, head, "./satured_head_LAST.txt");
    // calSaveSubstate2Dr(satured, K, "./satured_K_LAST.txt");
    // calSaveSubstate2Dr(satured, SS, "./satured_SS_LAST.txt");
    //calSaveSubstate2Dr(satured, convergence, "./satured_convergence_LAST.txt");

    // for(int i = 0; i< ROWS;i++){
    //    for(int j = 0; j< COLS;j++)
    //         printf("%d\t%d\t%f\n",j+1,i+1,head->current[i*COLS+j]);
    //     //printf(" \n");
    // }

    // save the Q substate to file
	//calSaveSubstate2Db(satured, Q, "./satured_LAST.txt");

	// finalize simulation and CA objects
	calRunFinalize2D(satured_simulation);
	calFinalize2D(satured);

}
