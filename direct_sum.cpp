
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <sys/times.h>

using namespace std;

int N_cube; // N_cube points total
int N;

double epsilon;
double epsq;
double  selfint;

static const double r8pi =  0.039788735772974; 
const bool UseSleep = false; // for testing memory usage purpose

//**********************************************************//

struct vec_3d
{
    double val[3];
};

//**********************************************************//
struct xyz // particle coordinates (physical)
{
	double* x;
	double* y;
	double* z;
	size_t* index;
	size_t* old_index;
	size_t size;
	xyz(size_t N_cube_in)
	{
		size = N_cube_in;
		x = new double[size];
		y = new double[size];
		z = new double[size];
		index = new size_t[size];
		old_index = new size_t[size];
	}
	~xyz()
	{
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] index;
		delete[] old_index;
	}
};

//*****************************************************************************//
long getTickCount()
{
    tms tm;
    return times(&tm);
}

//*****************************************************************************//

int main()
{
  std::ifstream ifile("input_direct.txt");
  if (ifile.is_open())
    {
      std::string line;
      std::getline(ifile,line);
      std:stringstream stream;
      stream << line << endl;
      stream >> N_cube >> epsilon;
    }
    ifile.close();

    epsq = epsilon * epsilon;
    selfint = 2.0*r8pi/epsilon;

    string kernelName;
    kernelName="RegStokeslet";

    struct xyz particles(N_cube);
    cout << "N_cube = " << N_cube << endl;
	

    // ********* read particle coordinates from a file *********

    FILE * fp;
    char S_data_file[64] = {0};
    sprintf(S_data_file, "./rand_%d.txt", N_cube);
    fp = fopen(S_data_file, "r");
    
    double x1, x2, x3;
    int count = -1;
    if (fp == NULL)
      {
	cout << "Cannot open random points file" << endl;
	getchar();
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

    
    // ********* read particle strengths from a file *********

    double *lambda[3];
    lambda[0] = new double[N_cube];
    lambda[1] = new double[N_cube];
    lambda[2] = new double[N_cube];
    
    char lambda_Str_data_file[64] = {0};
    sprintf(lambda_Str_data_file, "./lambda_%d.txt", N_cube);
    fp = fopen(lambda_Str_data_file, "r");
    
    count = -1;
    if (fp == NULL)
      {
	cout << "Cannot open lambda file" << endl;
	getchar();
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
    
    
    //***************** Direct summation *****************
    
    vec_3d *v_true = new vec_3d[N_cube];

    for (int i = 0; i < N_cube; i++)
      {

        v_true[i].val[0] = selfint*lambda[0][i]; 
        v_true[i].val[1] = selfint*lambda[1][i];
        v_true[i].val[2] = selfint*lambda[2][i];
      }

    long Start_ds, End_ds;
    Start_ds = getTickCount(); // Get currenct CPU time
    
    for (int i = 0; i < N_cube-1; i++)
      {
        double temp_x = particles.x[i];
        double temp_y = particles.y[i];
        double temp_z = particles.z[i];
        double s0 = 0.0;
        double s1 = 0.0;
        double s2 = 0.0;
        double fi1  = lambda[0][i];
	double fi2  = lambda[1][i];
        double fi3  = lambda[2][i];

        for (int j = i+1; j < N_cube; j++)
	  {

            double fj1  = lambda[0][j];
            double fj2  = lambda[1][j];
            double fj3  = lambda[2][j];

	    double xx = temp_x - particles.x[j];
	    double yy = temp_y - particles.y[j];
	    double zz = temp_z - particles.z[j];
	
            double x2 = xx*xx;
            double xy = xx*yy;
            double xz = xx*zz;
            double y2 = yy*yy;
            double yz = yy*zz;
            double z2 = zz*zz;
	
            double Resq  = x2 + y2 + z2 + epsq;
            double grj   = Resq + epsq;
            double grjx2 = grj+x2;
            double grjy2 = grj+y2;
            double grjz2 = grj+z2;
            double H2j   =  r8pi/(Resq*sqrt(Resq));
 
	    s0 += H2j*(fj1*grjx2 + fj2*xy     + fj3*xz);
            s1 += H2j*(fj1*xy    + fj2*grjy2  + fj3*yz);
            s2 += H2j*(fj1*xz    + fj2*yz     + fj3*grjz2);
 
	    v_true[j].val[0] += H2j*(fi1*grjx2 + fi2*xy     + fi3*xz);
            v_true[j].val[1] += H2j*(fi1*xy    + fi2*grjy2  + fi3*yz);
            v_true[j].val[2] += H2j*(fi1*xz    + fi2*yz     + fi3*grjz2);
            
	  }
	
        v_true[i].val[0] += s0;
        v_true[i].val[1] += s1;
        v_true[i].val[2] += s2;

      }

    End_ds = getTickCount(); // Get current CPU time
    long ds_cpu_time;
    ds_cpu_time = End_ds - Start_ds;
    
    cout << "ds time in seconds*100 = " << ds_cpu_time << endl;
    

    // ********* save results to a file *********

    char file_name[256];
    sprintf(file_name, "exact_sum_%s_N%d", kernelName.c_str(), N_cube);
    ofstream output_file(file_name);

    for (int i = 0; i < N_cube; i++)
      {
	 output_file << i << "   " << setprecision(16) << v_true[i].val[0]
                                   <<"  "<< v_true[i].val[1] 
                                   <<"  "<< v_true[i].val[2] << endl;

      }
    
    output_file << setprecision(16) << ds_cpu_time << endl;
    output_file.close();	   


    //*********************************************************
		
    delete [] lambda[0];
    delete [] lambda[1];
    delete [] lambda[2];
    delete [] v_true;
    
    return 0;
}
