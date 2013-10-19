/**********************************************************************
 * 		liquid.h
 * the code is a direct modification of the code by Leveque[1].
 *
 * References:
 * [1] R.J. Leveque, High-resolution conservative algorithms for
 *     advection in incompressible flow.
 *
 **********************************************************************/

#include <FronTier.h>
#include <vector>
#include <petscksp.h>
#include <assert.h>

class ADVECTION {
public:

	~ADVECTION();
	void initRectGrid(double*,double*,double*,int*);

	void advection(void);
	void advection(double **q);


//private:
	enum METHOD_TYPE {
	    METHOD_DonorCell = 1, 
	    METHOD_CornerTransportUpwind = 2,
	    METHOD_SecondOrderCorrectionTermsAdded =3,
	    METHOD_TransversePropagation2 = 4,	
	    // transverse propagation of 2nd order corrections
	    METHOD_ThirdOrderCorrectionTermsAdded = 5,
	    METHOD_TransversePropagation3 = 6,	
	    // transverse propagation of 3rd order corrections
	    METHOD_LaxWendroff = 7 
	};
	METHOD_TYPE m_method;
	enum LIMITER_TYPE {
	    LIMITER_NONE = 0,
	    LIMITER_MINMOD = 1, 
	    LIMITER_SUPERBEE = 2, 
	    LIMITER_VANLEER = 3, 
	    LIMITER_MC = 4
	};
	LIMITER_TYPE m_limiter;

	double getExactSolution(double coords[2]);
	void getVelocity(double coords[2]);

	void step(double **q, double t);
	void fluxes(double *q, double *u, double *v, double *deladd, 
			double *fadd, double **gadd, int i1, int i2, 
			double dx, double dt);
	void limiter(double *dq, double *u, int i1, int i2);
	double philim(double a, double b, LIMITER_TYPE type);

	void bc(double **q, int ibc, double tn);
	void init(double **q, double **ux, double **vx, double **uy, 
			double **vy, double dx, int mitot, int mjtot);

	double gt(double t);

	// member variables
	double m_T;
	int m_Nx;
	int m_Ny;
	int m_mitot;
	int m_mjtot;
	int m_mtot;
	double m_dh[MAXD];
	double m_dt;
	double *m_x;
	double *m_y;
	double **m_q;
	double **m_ux;
	double **m_vx;
	double **m_uy;
	double **m_vy;

	double *m_q1d;
	double *m_u1d;
	double *m_v1d;
	double **m_f;
	double **m_g;
	double **m_delta;
	double *m_deladd;
	double *m_fadd;
	double **m_gadd;

	double getMin(double a, double b);
	double getMin(double a, double b, double c);
	int getMax(int a, int b);
	double getMax(double a, double b);
	double getMax(double a, double b, double c);
	double getSqr(double a);
	void printError(char *filename, double **q, int mitot, int mjtot);
	void printMaxError(double **q, int mitot, int mjtot);		
};

typedef class CARTESIAN CARTESIAN_EB;

/**************************************************************************
 *		vector/matrix operation functions 
 **************************************************************************/

void VectorIntPrint(int my_rank, char *name, int n, int *vector);
void VectorPrint(int my_rank, char *name, int n, double *vector);
void VectorZero(int n, double *vector);
void VectorZero(int n, int *vector);
void VectorCopy(int n, double *vector1, double *vector2);
void MatrixPrint(int my_rank, char *name, int m, int n, double **matrix);
void MatrixZero(int m, int n, double **matrix);
void MatrixIdentity(int n, double **matrix);
void MatrixCopy(int m, int n, double **U2, double **U1);
void ArrayZero(int m, int n, int l, double ***matrix);
void Vector2Matrix(int m, int n, double **matrix, double *vector);
void ArrayPrint(int my_rank, char*name, int m, int n, int l, double ***array);
void MatrixMultiply(int n, double **C, double **A, double **B);	// C = AB
void MatrixMultiply(int n, double **C, double **A);		// C = A'A 
void MatrixMultiply(double C[3][3], double A[3][3], double B[3][3]);  // C = AB
void MatrixVectorMultiply(int n, double *C, double **A, double *b);   // c = Ab
void SymmetricMatrixInverse(double B[3][3], double A[3][3]); // B = inverse(A); 
double MatrixMax(int m, int n, double **matrix);
double MatrixMin(int m, int n, double **matrix);

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class CARTESIAN;

//typedef class CARTESIAN CARTESIAN_EB;
//enum VISITED_TYPE {UNVISITED, VISITED, PARTIAL_VISITED, FULL_VISITED};

//-------------------------------------------------
//		STATES
// NOTE:
//      INC_STATE/INC_STATE_RECT_EDGEshould be put into 
// lcartsn.h. However, there are some trouble
// to compile in that way.
//-------------------------------------------------
// states inside a cell

class INC_STATE{
public:
	double m_U[MAXD];		// velocity vector
	double m_dU[2][2];	// used for computing velocity jump
	double m_P;			// pressure

	double m_mu;		// smoothed 
	double m_rho;		// smoothed

	INC_STATE();
	void setZero(void);
};
// states on edge

//------------------------------------------------------
//		MESH
//------------------------------------------------------
// note that the following VERTEX2D/RECT_EDGE are different 
// from those defined in MESH2D.h

enum RECT_EDGE_DIR {RECT_EDGE_U=0, RECT_EDGE_V=1, RECT_EDGE_W=2, RECT_EDGE_NULL=-1};			// not used
class RECT_EDGE {
	//static CARTESIAN *m_pCartesian;
public:
	RECT_EDGE();
	int m_comp;			 
	int m_rectangle[2];		// neighbors
	RECT_EDGE_DIR m_dir;
	INC_STATE m_state;
	void setRectangle(int i, int j);
	int getOtherNeighbor(int indexRectangle);
};

// the order of the vertices is {00, 10, 11, 01},
// then the intersection points on the edges S, E, N, W;
// the order of the edges is {S,E,N,W};
class RECTANGLE {
public:
	int m_index;			// rectangle index
	int m_vertex[4];		// 
	int m_edge[4];			
	int m_comp;			 
	INC_STATE m_state;
	double m_area;
	bool m_visited;
	double m_coords[MAXD];	
	int icoords[MAXD];

	RECTANGLE();

	void setCoords(double*,int);
};


class CARTESIAN{
	Front *front;
public:
	// outside package
	ADVECTION m_advection;
	
	CARTESIAN(Front &front);

	// member data: RECT_GRID
	int dim;

	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	double *top_L,*top_U,*top_h;
	int *top_gmax;
	COMPONENT *top_comp;

	int *lbuf,*ubuf,*gmax;
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D

	// Sweeping limites
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	enum BC_TYPE { 									// used by ADVECTION
		BC_PERIODIC = 1,
		BC_Extrapolation0 = 2,
		BC_Extrapolation1 = 3,
		BC_InflowOutflow = 4};	
	BC_TYPE m_bc[4];								// down, right, up, left 		

	// member data: mesh storage
	std::vector<RECT_EDGE>	m_edges;
	std::vector<RECTANGLE> 	m_cell_center;

	// member data: 
	int m_comp[2];
	double m_mu[2];
	double m_rho[2];		// two component at most
	double m_smoothing_radius;	// used by getSmoothingFunction()

	double m_t;                     // time
	double m_dt;			// time increment

	// constructor
	~CARTESIAN();

	// mesh
	void initRectGrid(Front*);
					// init the cartesian grid

	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setComponent(void);	// init components	
	void setSmoothedProperties(void); // smooth discontinuous properties
	double getSmoothingFunction(double r);	// Heaviside function
	void makeGridIntfc();
	void deleteGridIntfc();
	void readMeshStates(char *state_name);
	void printMeshStates(char *state_name);

	/*TMP*/
	void checkStates();

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex();

	// physics calculation
	void setInitialCondition(void);

	void setIndexMap(void);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.

	// -------------------------------------------------------
	// 		incompressible solver functions
	// -------------------------------------------------------
	void interpolateEdgeVelocity(void);
	void computeAdvectionSimple(void);
	
	double computeVelocityDivergence(int indexRectangle);
	void computePressureGradient(int indexRectangle, double *gradP);
	void compPressGradWithSmoothProperty(int indexRectangle, 
			double *gradP);

	void computeDiffusion(void);	// Crank Nicolson method, wrong
	void compDiffWithSmoothProperty(void);
					// Crank Nicolson method

	void compProjWithSmoothProperty(void);		

	void compNewVeloWithSmoothProperty(void);
					// time n+1 velocity

	void computeSourceTerm(double *coords, INC_STATE &state); 
	void computeSourceTerm(double *coords, double t, INC_STATE &state); 

	void computeExactSolution(double *coords, INC_STATE &state); 

	void getPressureJumpParameter(double *coords0, double *coords1, 
			double &theta, double &jumpPressure, 
			double &jumpDerivative);

	// main step function
	void solve(double dt);		

	// velocity field query
	void getVelocity(double *p, double *U);
	void getVelocityGradient(double *p, double *gradU, double *gradV);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	int getEdgeIndex(double i, double j);
	void getEdgeIndex(int indexEdge, double &i, double &j);

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int getRectangleComponent(int index);	// the center component
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	
	int getInteger(double i);
	bool isInteger(double i);

	double getDistance(double *coords0, double *coords1);
	
			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	void save(char *filename);
};

class SOLVER
{
public:
	SOLVER(){};
	SOLVER(int ilower, int iupper, int d_nz, int o_nz){};
	virtual ~SOLVER(){};
	virtual void Create(int ilower, int iupper, int d_nz, int o_nz){};
	
	virtual void Set_A(int i, int j, double val){};		// A[i][j]=val;
	virtual void Add_A(int i, int j, double val){};		// A[i][j]=A[i][j]+val;
	virtual void Set_x(int i, double val){};		// x[i]=val;
	virtual void Set_x(double *p){};			// x[i]=p[i];
	virtual void Add_x(int i, double val){};		// x[i]=x[i]+val;
	virtual void Get_x(double *p){};			// get the x from ij_x to p.		
	virtual void Get_x(double *p, int n, int *global_index){};
	virtual void Set_b(int i, double val){};		// b[i]=val;
	virtual void Set_b(double *b){};	
	virtual void Add_b(int i, double val){};		// b[i]=b[i]+val;

	virtual void SetMaxIter(int val){};	
	virtual void GetFinalRelativeResidualNorm(double *rel_resid_norm){};
	virtual void GetNumIterations(int *num_iterations){};

	virtual void Solve(void){};	
	virtual void Solve_withPureNeumann(void){};	
	virtual void Read_A(char *filename){};
	virtual void Print_A(char *filename){};
	virtual void Read_b(char *filename){};
	virtual void Print_b(char *filename){};
	virtual void Read_x(char *filename){};
	virtual void Print_x(char *filename){};
	virtual void test(void){};
};

class PETSc: public SOLVER
{
public:	
	MPI_Comm  comm;			// set to be MPI_COMM_WORLD.
	int iLower;
	int iUpper;			// global row range
	
	Vec x;      			/* approx solution, RHS*/
	Vec b;
  	Mat A;          		/* linear system matrix */
  	
  	KSP   ksp;          		/* Krylov subspace method context */
	PC    pc;
	MatNullSpace	nullsp;
		
	int ierr;
	int its;			// numer of iterations;

public:
	PETSc();
	PETSc(int ilower, int iupper, int d_nz, int o_nz);		
		// global row range of A, x, b on this processor
	~PETSc();
	void Create(int ilower, int iupper, int d_nz, int o_nz);	
		// same as Hypre(int, int)
	void Create(MPI_Comm Comm, int ilower, int iupper, int d_nz, int o_nz);
		// same as Hypre(int, int)
	
	void Reset_A();				// Set A[i][j]=0.0;
	void Set_A(int i, int j, double val);	// A[i][j]=val;
	void Add_A(int i, int j, double val);	// A[i][j]=A[i][j]+val;
	void Get_row_of_A(int i, int *ncol, int **cols, double **row);
	void Set_x(int i, double val);		// x[i]=val;
	void Add_x(int i, double val);		// x[i]=x[i]+val;
	void Set_b(int i, double val);		// b[i]=val;
	void Add_b(int i, double val);		// b[i]=b[i]+val;
	void Get_x(double *p);		// get the x from ij_x to p.	
	void Get_b(double *p);		// get the b from ij_x to p.
	void Get_x(double *p, int n, int *global_index);
	
	void SetMaxIter(int val); 	// Set maximum number of iterations 
	void SetTol(double val);	// Set the convergence tolerance 
	void SetKDim(int k_dim);	
			// Set the maximum size of the Krylov space 
	void GetNumIterations(int *num_iterations);	
			// Return the number of iterations taken 
	void GetFinalRelativeResidualNorm(double *rel_resid_norm);
	void Solve(void);
	void Solve_withPureNeumann(void);
};

