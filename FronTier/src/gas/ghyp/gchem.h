//#include <gdecs/gstate.h>

/* Declaration of chemical structure */
#define NUM_OF_ELEMENTS 4
#define NUM_OF_SPECIES 10
#define NUM_OF_REACTIONS 37
#define BOLTZMANN .000000000000000000000013806503  //units are in J/K
#define RGAS 83.1434	

typedef struct Vector {
	int len;        /* length of the vector */
	double *vec;         /* pointer to the elements */
} Vector;
typedef Vector *VectorPtr;

typedef struct IntVector {
	int len;
	int *vec;
} IntVector;
typedef IntVector *IntVectorPtr;

typedef struct _Species {
	int index;
	int track;
	double MolarMass;
	double MassFraction;
	double ProductionRate;
	VectorPtr NASACoeff;             /* for >1000K only */
	int NofUsedReactions;
	IntVectorPtr UsedReactions;
	IntVectorPtr SSt;		/* list of stoichiometric coefficients in the same order */
}Species;

Species species[NUM_OF_SPECIES];

struct TBcomp {
    	int index;
	double coeff;
};	/* Third body components */

typedef struct _Reaction {
	int index;
	int flag_fm;     /* 0 for basic form, 1 for TROE form */
	int flag_tb;			/* 1 if containing third body */
	int NofReactants;
	IntVectorPtr ReactantIdx;		/* list of reactant indices */
	IntVectorPtr ReactantSt;		/* list of reactant stoichiometric coefficients */
	VectorPtr ReactantMM;		/* list of reactant molarmass, ReactionMM->vec is an array */
	int NofProducts;
	IntVectorPtr ProductIdx;
	IntVectorPtr ProductSt;
	VectorPtr ProductMM;

	double HeatRelease;
	double ReactionRate;
	double RateCoeff;
	double BackwardRateCoeff;

	/* rate constant parameters */
	double PreExp;		/* preexponential constant A */
	double PreExp_i;	/* additional parameter in TROE form */
	double TempExp;		/* temperature exponent n_k */
	double TempExp_i;
	double ActEnergy;	/* activation energy */
	double ActEnergy_i;
	double fca;		/* TROE parameters */
	double fcta;
	double fcb;
	double fctb;
	double fcc;
	double fctc;

	/* third body parameters */
	int tbIndex;		/* third body index, 0 for pure composition */
	int nofcomp;
	struct TBcomp ThirdBody[6];	/* hardcoded for 6 temporarily (including 6 components) */

}Reaction;

Reaction reaction[NUM_OF_REACTIONS];

/*extern int yylex (void); //line 700 of lex.yy.c
extern FILE *yyin; //line 173 of lex.yy.c
void initDataStructures(void);
void readMolarMassData(void);
void readNASAData(void);
void putMMinReaction(void);
void freeDataStructures(void);
void ReadThermoData(void);
void ComputeChemSourceTerm(Vec_Gas *,Vec_Src *, float, float, int, int);
void ComputeHeatRelease(int, double);*/
