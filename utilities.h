#include <cstddef>
#include <cmath>

static const double PI = 3.14159265358979323846;
static const double rootPI = sqrt(PI);
static const double kappa = 1.0;
static const double r8pi =  0.039788735772974;
static const int P = 3; // fixed for tricubic approximation
static const int Pflat = (P + 1)*(P + 1)*(P + 1);
const bool UseSleep = false; // for testing memory usage purpose

extern int N_cube; // number of particles

extern int fdiff, mgrid;
extern double hgrid, rdr, twoh, rtwoh, hsq, rfhsq, r8hcb;
extern double* fval; //Pointer to double

extern double epsilon, epsq; // epsilon and epsilon square 
                             // regularization parameter for 
                             // regularized Stokeslet

extern double sweight,sweight1,sweight2,sweight3; //weight of source particle in cluster-particle
extern double sfweight,sfweight1,sfweight2,sfweight3; //scaled weight of source particle in cluster-particle

extern int sID; //index of source particle in cluster-particle
extern double *cpdvx; // cp treecode approximation of the x-component of the field
extern double *cpdvy; // cp treecode approximation of the y-component of the field
extern double *cpdvz; // cp treecode approximation of the z-component of the field

extern int fid; //internal parameter for regularized Stokeslet
                //Indicates which of 6 functions to be used
                //for computations at the cluster evaluation points

extern double **Binv; //B-inverse matrix

extern double divergence; // the divergence of regularized Stokeslet

//**********************************************************//

struct vec_3d
{
    double val[3];
};

//**********************************************************//

struct vec_4d
{
    double val[4];
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

//**************//

long getTickCount();

void init_vec(double vec[]);

void init_vec3d(double vec[][64]);

void read_particle_data(int N_cube,
			struct xyz &particles,
			double *lambda[3]);

