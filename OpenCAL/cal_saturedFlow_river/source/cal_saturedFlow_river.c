#include <OpenCAL/cal2D.h>
#include <OpenCAL/cal2DIO.h>
#include <OpenCAL/cal2DRun.h>
#include <OpenCAL/cal2DBuffer.h>
#include <OpenCAL/cal2DReduction.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>


#define ROWS 100
#define COLS 100
#define LAYERS 1
#define CELL_SIZE_X 10
#define CELL_SIZE_Y 10
#define SPESSORE 50

#define Syinitial 0.1 
#define Kinitial  0.0000125

#define headFixed 50
#define headCalculated 50

struct CALModel2D* satured;
struct CALSubstate2Dr* head; // carico
struct CALSubstate2Dr* K; // permeabilità
struct CALSubstate2Dr* Sy; // immagazinamento
struct CALSubstate2Dr* Sorgente; // sorgente



struct CALSubstate2Dr* convergence; // convergence

struct CALSubstate2Dr* Mod; // immagazinamento

struct CALRun2D* satured_simulation;

float delta_t_ = 4000;//0.5;
double Ksb = 0.00001;//0.0000001;
double M   = 0.5;
double Hw  = 48;
double W  = 5;
int time_10 = 860000;
double riverBottom = 46.5;


// int rowMonitoringWell = 50;
// int colMonitoringWell = 65;

//FILE * monitorWellFile;
const int numberofWell = 5;
FILE * monitorWellFiles[5]; //= malloc(sizeof(FILE)*numberofWell);


const char* name = "monitorWell";
const char* extension = ".txt";

//double convergence = ((CELL_SIZE_X*CELL_SIZE_Y)*Syinitial)/(Kinitial*4);

int rowMonitoringWell[5] = {49,49,49, 49,49};
int colMonitoringWell[5] = {10,35,45, 55,65}; 

//double convergence = ((CELL_SIZE_X*CELL_SIZE_Y)*Syinitial)/(Kinitial*4);

// FILE * monitorWellFile;
FILE * HWFile;
int numberofdays =1;
int riverPosition = 24;
void saturedTransitionFunctionHw(struct CALModel2D* ca)
{
    double time_ =delta_t_*satured_simulation->step;
    if(864000 <= time_ && time_ < (15*86400))
    {
        Hw+=0.0185;
    }
    if((86400*20) > time_ && time_ >= 15*86400)
    {
        Hw-=0.0185;
    }

    if((86400*20) < time_)
        Hw = 48;
    
    int days = time_/86400;
    
    //printf("days = %d, seconds = %f, Hw = %f\n",days,time_, Hw  );
  //  printf("%d, %d, Hw = %f\n",time_10,time_10+4000, Hw );
    //time_10+=4000;


    if((numberofdays)*86400 <= time_)
    {
        for(int i = 0; i < numberofWell; i++){
           fprintf(monitorWellFiles[i], "%d %f \n", days, calGet2Dr(ca, head, rowMonitoringWell[i], colMonitoringWell[i]));
           
        }
        fprintf(HWFile, "%d %f \n", days, Hw);
        numberofdays++;
    }
   
}

void saturedTransitionFunction(struct CALModel2D* ca, int i, int j)
{
    if(j==0 || j ==COLS-1){
        return;
    }
    double diffHead=0.0;
    double Q=0.0;
    double tmpT = 0.0;  
    for (int n=1; n<satured->sizeof_X; n++){
     
      if ( i == 0 && n == 1){
        diffHead = 0;
      }else 
        if (i == ROWS -1 && n == 4){
        diffHead = 0;
      }else{
            diffHead = (calGetX2Dr(ca, head, i, j, n)- calGet2Dr(ca, head, i, j));
      }
            tmpT = calGetX2Dr(ca, K, i, j, n)*SPESSORE;
            Q += (tmpT * diffHead);
    }
    
    if(j == riverPosition){
       //  printf("Q = %f  ", Q);
        double firsttermQ = (Ksb*(CELL_SIZE_X)*W)/M;
        double Qriver = 0.0;
        if(calGet2Dr(ca, head, i, j) > riverBottom)
        {
            Qriver = firsttermQ*(Hw-calGet2Dr(ca, head, i, j));
        }else
        if(calGet2Dr(ca, head, i, j) <= riverBottom)
        {
            Qriver = firsttermQ*(Hw+M);
        }
        Q+=Qriver;
    }

    double area = CELL_SIZE_X*CELL_SIZE_Y;
    double ht1 = (Q*delta_t_)/(calGet2Dr(ca, Sy, i, j)*area);
    calSet2Dr(ca, head, i, j, ht1+calGet2Dr(ca, head, i, j));
}




void saturedInit(struct CALModel2D* ca)
{
for(int i = 0; i< ROWS;i++)
        for(int j = 0; j< COLS;j++)
            {        
              if( j == COLS-1 ) 
                calSet2Dr(satured, head, i, j, headCalculated);
              else
                calSet2Dr(satured, head, i, j, headFixed);
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
}

#define STRLEN 256
void calfLoadMatrix2Di(int* M, int rows, int columns, FILE* f)
{
  char str[STRLEN];
  int i, j;

  for (i=0; i<rows; i++)
    for (j=0; j<columns; j++){
      fscanf(f, "%s", str);
      calSetMatrixElement(M, columns, i, j, atoi(str));
    }
}

bool calLoadMatrix2Di(int* M, int rows, int columns, char* path)
{
  FILE *f = NULL;
  f = fopen(path, "r");

  if ( !f )
    return CAL_FALSE;

  calfLoadMatrix2Di(M, rows, columns, f);

  fclose(f);

  return CAL_TRUE;
}


int main(){

    for(int i = 0; i < numberofWell; i++){
        char* name_with_extension;
        
        name_with_extension = malloc(strlen(name)); /* make space for the new string (should check the return value ...) */
        
        strcpy(name_with_extension, name); /* copy name into the new var */
        
        char char_arr [100];
        printf("%i\n", colMonitoringWell[i]);
        sprintf(char_arr, "%d", colMonitoringWell[i]);
        strcat(name_with_extension, char_arr);
        strcat(name_with_extension, extension); /* add the extension */ 
        printf("%s\n", name_with_extension);

        monitorWellFiles[i] =  fopen(name_with_extension, "w");
        fclose(monitorWellFiles[i]);
        monitorWellFiles[i] =  fopen(name_with_extension, "a");    
    }
    
    HWFile =  fopen("HW.txt", "w");
    fclose(HWFile);
    HWFile =  fopen("HW.txt", "a");   
   
    // define of the satured CA and satured_simulation simulation objects
	satured = calCADef2D(ROWS, COLS, CAL_VON_NEUMANN_NEIGHBORHOOD_2D, CAL_SPACE_FLAT, CAL_NO_OPT);
    int numberOfSteps = (31*86400)/delta_t_;
    printf("number of steps %d \n", numberOfSteps);

    // return 0;

	satured_simulation = calRunDef2D(satured, 1, numberOfSteps, CAL_UPDATE_IMPLICIT);

	// add the Q substate to the satured CA
	head = calAddSubstate2Dr(satured);
    K = calAddSubstate2Dr(satured);
    Sy = calAddSubstate2Dr(satured);
    Sorgente = calAddSubstate2Dr(satured);
    convergence = calAddSubstate2Dr(satured);

	// add transition function's elementary process
    //calAddElementaryProcess2D(satured, saturedTransitionFunctionHw);
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
    calRunAddSteeringFunc2D(satured_simulation, saturedTransitionFunctionHw);
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

    for(int i = 0; i < numberofWell; i++){
       fclose(monitorWellFiles[i]); 
    }
    fclose(HWFile);

}