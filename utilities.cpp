#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/times.h>
#include <cmath>

#include "utilities.h"

using namespace std;

int N_cube; // number of particles

int fdiff, mgrid;
double hgrid, rdr, twoh, rtwoh, hsq, rfhsq, r8hcb;
double* fval = NULL; //Pointer to double
double sweight1, sweight2, sweight3, sfweight1, sfweight2, sfweight3;
double epsq, epsilon, divergence;
int sID;
int fid;
double* cpdvx  = NULL;
double* cpdvy  = NULL;
double* cpdvz  = NULL;
double** Binv = NULL;

//************************************************************************//

long getTickCount()
{
  tms tm;
  return times(&tm);
}

//************************************************************************//

// initialize a vector --- assumes size is 64

void init_vec(double vec[])
{
  for (int i=0; i<64; i++)
    vec[i]=0.0;
}

//*****************************************************************************//

void init_vec3d(double vec[][64])
{
  for (int i=0; i<64; i++)
   {
    vec[0][i]=0.0;
    vec[1][i]=0.0;
    vec[2][i]=0.0;
    vec[3][i]=0.0;
    vec[4][i]=0.0;
    vec[5][i]=0.0;
   }
}

//************************************************************************//

// read particle coordinates and weights from files

void read_particle_data(int N_cube,
			struct xyz &particles,
			double *lambda[3])
{
  FILE * fp;
  char S_data_file[64] = {0};
  sprintf(S_data_file, "./rand_%d.txt", N_cube);
  fp = fopen(S_data_file, "r");
  
  double x1, x2, x3;
  int count = -1;
  if (fp == NULL)
    {
      cout << "Cannot open random points file" << endl;
      exit(-1);
    }
  
  while (true)
    {
      count++;
      fscanf(fp,"%lf%lf%lf", &x1, &x2, &x3);
      if (feof(fp))
	break;
      if (count >= N_cube)
	{
	  cout << "Out of range" << endl;
	  exit(-1);
	}
      
      particles.x[count] = x1;
      particles.y[count] = x2;
      particles.z[count] = x3;
      particles.index[count] = -1;
      particles.old_index[count] = count;
    }
  
  
  // ********* read particle weights from a file *********
  
  char lambda_Str_data_file[64] = {0};
  sprintf(lambda_Str_data_file, "./lambda_%d.txt", N_cube);
  fp = fopen(lambda_Str_data_file, "r");
  
  count = -1;
  if (fp == NULL)
    {
      cout << "Cannot open lambda file" << endl;
      exit(-1);
    }
  
  while (true)
   {
      count++;
      fscanf(fp,"%lf%lf%lf", &x1, &x2, &x3);
      if (feof(fp))
        break;
      if (count >= N_cube)
        {
          cout << "Out of range" << endl;
          exit(-1);
        }
      
      lambda[0][count] = x1;
      lambda[1][count] = x2;
      lambda[2][count] = x3;            
    }
}

