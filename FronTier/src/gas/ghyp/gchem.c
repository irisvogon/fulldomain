#include <unistd.h>
#include "lex.yy.c" //TMP_XY
#include <math.h>
#include <vmalloc.h>

void ComputeRateCoefficients( double *, double, double, double, double );
void ComputeBackwardRateCoefficients( double, int);
void ComputeReactionRate( double, double *, int, double);
void ComputeProductionRate( int );
double myLog10( double );

int flag = 0;

void ReadThermoData()
{	
    	
	/*In the if-statement below we check to see if the reaction array is NULL.  
		If so, we have not yet read the thermodynamic data.  If not, we 
		have already read the data and there is no need to read it again.
		In general, this type of programming is bad practice and should be
		changed in the future.*/
	if(flag == 0) {
	    	flag = 1;
		int i, j;
		int label;
		yyin = fopen("./ghyp/ThermoData/Stanfordv12-mech-no-Argon-FM.mech","r");
		initDataStructures();
		readMolarMassData();
		yylex();
		readNASAData();
		putMMinReaction();

		for ( i = 0; i < NUM_OF_REACTIONS; i++ ) {
	    		for ( j = 0; j < reaction[i].NofReactants; j++ ) {
				label = reaction[i].ReactantIdx->vec[j];
				species[label].SSt->vec[species[label].track] = reaction[i].ReactantSt->vec[j];
				species[label].track++;
	    		}		
	    		for ( j = 0; j < reaction[i].NofProducts; j++ ) {
				label = reaction[i].ProductIdx->vec[j];
				species[label].SSt->vec[species[label].track] = reaction[i].ProductSt->vec[j];
				species[label].track++;
	    		}
		}
		//freeDataStructures();
	}else{
		return;
	}
}

/*void ComputeChemSourceTerm(Vec_Gas *vst, Vec_Src *src, float temp, float dens, int cell, int num_comps)
{
	int i;
	double Y[num_comps];

	for (i=0 ;i < num_comps; i++){
		Y[i] = vst->rho0[i][cell]; 	//'cell' refers to the mesh block.
	}
	
	for (i = 0; i < NUM_OF_REACTIONS; i++) {
	    ComputeRateCoefficients( &reaction[i].RateCoeff, reaction[i].PreExp, reaction[i].TempExp, reaction[i].ActEnergy / RGAS, temp);
	    ComputeBackwardRateCoefficients( temp, i);
	    ComputeReactionRate( dens, Y, i );
	    //TMP_XY
	    //printf("For reaction %d:\n", i);
	    //printf("Temperature = %.25g\n", temp);
	    //printf("Rate coefficient = %.25g\t Backward rate coefficient = %.25g\t Reaction rate = %.25g\n", reaction[i].RateCoeff, reaction[i].BackwardRateCoeff, reaction[i].ReactionRate);
	}

	for ( i = 0; i < num_comps; i++ ) {
	    species[i].ProductionRate = 0;
	    ComputeProductionRate( i );
	    //TMP_XY
	    //printf("For species %d, the production rate = %.25g\n", i, species[i].ProductionRate);
	    //printf("Coordinates: (%lf, %lf, %lf)\n", vst->coords[cell][0], vst->coords[cell][1], vst->coords[cell][2]);

	    src->rho0[i][cell] = species[i].ProductionRate;
	}
	//printf("Leaving ComputeChemSourceTerm.\n");
}*/

void ComputeRateCoefficients( double *k, double A, double n, double eOverR, double temp )
{
    if (n) {
	if (fabs( eOverR ) > 1.0e-3 ) {
	    *k = A * exp( n * log( temp ) - eOverR / temp );
	}
	else {
	    *k = A * pow( temp, n );
	}
    }
    else {
	if ( fabs( eOverR ) > 1.0e-3 ) {
	    *k = A * exp( -eOverR / temp );
	}
	else {
	    *k = A;
	}
    }
    //TMP_XY
    //printf("In ComputeRateCoefficients, k = %lf, A = %lf, n = %lf, eOverR = %lf, temp = %lf\n", *k, A, n, eOverR, temp);
}


void ComputeBackwardRateCoefficients( double temp, int i)
{
    //MOD_XY
    //double p0 = 0.10133;	/* unit of pressure: g*cm^(-1)*ms^(-2) */
    double p0 = 1.0133;		/* unit of pressuer: g*cm^(-1)*ms^(-2) */
    double lnROverP0 = log ( RGAS/p0 );
    double RT = RGAS * temp;
    double lnT = log ( temp );
    int j;
    int n = reaction[i].NofReactants + reaction[i].NofProducts;
    double *nu;
    uni_array(&nu,n,sizeof(double));
    double sumNu = 0.0;
    double sumNuMu = 0.0;
    double *a;
    int index;
    double lnKC;
    double k_f = reaction[i].RateCoeff;

    if ( reaction[i].RateCoeff <= 1.e-200 )
    {
	reaction[i].BackwardRateCoeff = 0.0;
	free(nu);
	return;
    }
    for ( j = 0; j < reaction[i].NofReactants; j++ ) {
	nu[j] = reaction[i].ReactantSt->vec[j];
    }
    
    for ( ; j < n; j++ ) {
	nu[j] = reaction[i].ProductSt->vec[j-reaction[i].NofReactants];
    }

    for ( j = 0; j < n; j++ )
    {
	//MOD_XY
	//index = reaction[i].ReactantIdx->vec[j];
	index = ( j < reaction[i].NofReactants )? reaction[i].ReactantIdx->vec[j] : reaction[i].ProductIdx->vec[j-reaction[i].NofReactants];
	double mu = 0.0;
	a = species[index].NASACoeff->vec;
	if (temp >= 1000)
	{
	    mu = a[0] * (1.0 - lnT) + a[5] / temp - a[6];
	    mu -= 0.5 * temp * (a[1] + temp * ( a[2] / 3.0 + temp * ( a[3] / 6.0 + 0.1 * temp * a[4] ) ) );
	}
	else
	{
	    mu = a[7] * (1.0 - lnT) + a[12] / temp - a[13];
	    mu -= 0.5 * temp * (a[8] + temp * ( a[9] / 3.0 + temp * ( a[10] / 6.0 + 0.1 * temp * a[11] ) ) );
	}
	mu *= RT;
	sumNu += nu[j];
	sumNuMu += nu[j] * mu;
    }
    lnKC = -sumNu * ( lnT + lnROverP0 ) -sumNuMu / RT;
    reaction[i].BackwardRateCoeff = k_f / exp( lnKC );
    free(nu);
    //TMP_XY
    //printf("Leaving CompouteBackwardRateCoefficients.\n");
}

void ComputeReactionRate( double density, double *Y, int i, double TBConc)
{
    //TMP_XY
    //printf("Entering ComputeReactionRate.\n");
    int j;
    double r = reaction[i].RateCoeff;
    double c = reaction[i].BackwardRateCoeff;
    double r1 = 1.0;
    double c1 = 1.0;
    int nr = reaction[i].NofReactants;
    int np = reaction[i].NofProducts;
    int *nofrs = reaction[i].ReactantSt->vec;
    int *nofps = reaction[i].ProductSt->vec;

    for ( j = 0; j < nr; j++ ) {
	if( Y[reaction[i].ReactantIdx->vec[j]] <= 0 )
	{
	    r1 = 0;
	    break;
	}
	if( fabs(nofrs[j]) > 1 ) {
	    //TMP_XY
	    //printf("Testing nofrs[j]: Absolute value of %lf is: %lf\n", nofrs[j], fabs(nofrs[j]));
	    //printf("\nReactant index: %d\tReactant molar mass: %f\tReactant stochiometry: %d\n", reaction[i].ReactantIdx->vec[j], reaction[i].ReactantMM->vec[j], nofrs[j]);
	    r1 *= pow( density * Y[reaction[i].ReactantIdx->vec[j]] / reaction[i].ReactantMM->vec[j], fabs(nofrs[j]) );
	    //TMP
	    //printf("r1 = %e\t", r1);
	    //printf("Test power function: %e\n", pow( density * Y[reaction[i].ReactantIdx->vec[j]] / reaction[i].ReactantMM->vec[j], fabs(nofrs[j])));
	}
	else {
	    //TMP
	    //printf("\nReactant index: %d\tReactant molar mass: %f\tReactant stochiometry: %d\n", reaction[i].ReactantIdx->vec[j], reaction[i].ReactantMM->vec[j], nofrs[j]);
	    r1 *= density * Y[reaction[i].ReactantIdx->vec[j]] / reaction[i].ReactantMM->vec[j];
	    //TMP
	    //printf("r1 = %e\t", r1);
	}
    }
    for (j = 0; j < np; j++ ) {
	    //TMP_XY
	    /*if(i == 5)
	    {
		printf("Density: %lf\t Y: %lf\t Molarmass: %lf\t Product stoichiometric: %lf\n", density, Y[reaction[i].ProductIdx->vec[j]], reaction[i].ProductMM->vec[j], fabs(nofps[j]));
		printf("original c = %lf\t new c = %lf\t j = %d\t np = %d\n", reaction[i].BackwardRateCoeff, density*Y[reaction[i].ProductIdx->vec[j]] / reaction[i].ProductMM->vec[j], j, np);
	    } */
	if( Y[reaction[i].ProductIdx->vec[j]] <= 0 )
	{
	    c1 = 0;
	    break;
	}
	
	if( fabs(nofps[j]) > 1 ) {
	    //TMP
	    //printf("\nProduct index: %d\tProduct molar mass: %f\tProduct stochiometry: %d\n", reaction[i].ProductIdx->vec[j], reaction[i].ProductMM->vec[j], nofrs[j]);
	    c1 *= pow( density * Y[reaction[i].ProductIdx->vec[j]] / reaction[i].ProductMM->vec[j], fabs(nofps[j]) );
	    //TMP
	    //printf("\nc1 = %e\t", c1);
	}
	else {
	    //TMP
	    //printf("\nProduct index: %d\tProduct molar mass: %f\tProduct stochiometry: %d\n", reaction[i].ProductIdx->vec[j], reaction[i].ProductMM->vec[j], nofrs[j]);
	    c1 *= density * Y[reaction[i].ProductIdx->vec[j]] / reaction[i].ProductMM->vec[j];
	    //TMP
	    //printf("c1 = %e\t", c1);
	}
    }
    //TMP_XY
    //if (i == 5) printf("c = %.15lf\n", c);
    r *= r1;
    c *= c1;
    reaction[i].ReactionRate = r - c;
    reaction[i].ReactionRate *= TBConc;
    //TMP_XY
    //if ((i >= 4 && i <= 7) || i == 14 ) reaction[i].ReactionRate *= density * Y[8] / (14.007 * 2);
    //TMP_XY
    //printf("Leaving ComputeReactionRate.\n");
}

void ComputeProductionRate(int i)
{
    //TMP_XY
    //printf("Entering ComputeProductionRate.\n");
    int j;
    //TMP_1219
/*    if (i == 9)
    printf("For species OHS:\n");*/
    for ( j = 0; j < species[i].NofUsedReactions; j++ ) {
	species[i].ProductionRate += reaction[species[i].UsedReactions->vec[j]].ReactionRate * species[i].SSt->vec[j];
	//TMP_1219
/*	if ((i == 9) && (reaction[36].ReactionRate != 0))
	{
	    printf("reaction rate %d: %e\tSSt: %d\n",i, reaction[species[i].UsedReactions->vec[j]].ReactionRate, species[i].SSt->vec[j]);
	}*/
    }
    species[i].ProductionRate *= species[i].MolarMass;
/*    if (i == 9)
    {
	printf("molar mass %d: %f\tproduction rate: %e\n",i, species[i].MolarMass, species[i].ProductionRate);
    }*/
    //TMP_XY
    //printf("Leaving ComputeProductionRate.\n");
}


void ComputeHeatRelease( int i, double temp )
{
    //TMP_XY
    //printf("Entering ComputeHeatRelease.\n");
    int k, j, n = reaction[i].NofReactants + reaction[i].NofProducts;
    double *a, *h;
    //MOD_XY
    //double hr = 1.0;
    double hr = 0.0;
    //double hr1 = 0.0;	
    uni_array(&h,n,sizeof(double));
    for ( j = 0; j < reaction[i].NofReactants; j++ ) {
        k = reaction[i].ReactantIdx->vec[j];
        a = species[k].NASACoeff->vec;
	if ( temp >= 1000 )
	{
	    //TMP_XY
	    //h[j] = ( a[0] + a[1] * temp / 2 + a[2] * pow( temp, 2 ) / 3 + a[3] * pow( temp, 3 ) / 4 + a[4] * pow( temp, 4 ) / 5 + a[5] / temp ) * RGAS * temp;
	    //h[j] = ( a[0] + a[1] * temp / 2 + a[2] * pow( temp, 2 ) / 3 + a[3] * pow( temp, 3 ) / 4 + a[4] * pow( temp, 4 ) / 5 + a[5] / temp ) * RGAS * temp / species[k].MolarMass;
	    h[j] = ( a[0] + ( a[1] / 2 + ( a[2] / 3 + ( a[3] / 4 + a[4] / 5 * temp ) * temp ) * temp ) * temp + a[5] / temp ) * RGAS * temp / species[k].MolarMass;
	}
	else
	{
	    //TMP_XY
	    //h[j] = ( a[7] + a[8] * temp / 2 + a[9] * pow( temp, 2 ) / 3 + a[10] * pow( temp, 3 ) / 4 + a[11] * pow( temp, 4 ) / 5 + a[12] / temp ) * RGAS * temp;
	    //h[j] = ( a[7] + a[8] * temp / 2 + a[9] * pow( temp, 2 ) / 3 + a[10] * pow( temp, 3 ) / 4 + a[11] * pow( temp, 4 ) / 5 + a[12] / temp ) * RGAS * temp / species[k].MolarMass;
	    h[j] = ( a[7] + ( a[8] / 2 + ( a[9] / 3 + ( a[10] / 4 + a[11] / 5 * temp ) * temp ) * temp ) * temp + a[12] / temp ) * RGAS * temp / species[k].MolarMass;
	}
        //hr += h[j] * species[k].MolarMass * reaction[i].ReactantSt->vec[j];
	hr += h[j] * reaction[i].ReactionRate * reaction[i].ReactantSt->vec[j] * species[k].MolarMass;
    }
    for ( j = 0; j < reaction[i].NofProducts; j++ ) {
        k = reaction[i].ProductIdx->vec[j];
        a = species[k].NASACoeff->vec;
	if ( temp >= 1000 )
	{
	    //TMP_XY
	    //h[j] = ( a[0] + a[1] * temp / 2 + a[2] * pow( temp, 2 ) / 3 + a[3] * pow( temp, 3 ) / 4 + a[4] * pow( temp, 4 ) / 5 + a[5] / temp ) * RGAS * temp;
	    //h[j] = ( a[0] + a[1] * temp / 2 + a[2] * pow( temp, 2 ) / 3 + a[3] * pow( temp, 3 ) / 4 + a[4] * pow( temp, 4 ) / 5 + a[5] / temp ) * RGAS * temp / species[k].MolarMass;
	    h[j] = ( a[0] + ( a[1] / 2 + ( a[2] / 3 + ( a[3] / 4 + a[4] / 5 * temp ) * temp ) * temp ) * temp + a[5] / temp ) * RGAS * temp / species[k].MolarMass;
	}
	else
	{
	    //TMP_XY
	    //h[j] = ( a[7] + a[8] * temp / 2 + a[9] * pow( temp, 2 ) / 3 + a[10] * pow( temp, 3 ) / 4 + a[11] * pow( temp, 4 ) / 5 + a[12] / temp ) * RGAS * temp;
	    //h[j] = ( a[7] + a[8] * temp / 2 + a[9] * pow( temp, 2 ) / 3 + a[10] * pow( temp, 3 ) / 4 + a[11] * pow( temp, 4 ) / 5 + a[12] / temp ) * RGAS * temp / species[k].MolarMass;
	    h[j] = ( a[7] + ( a[8] / 2 + ( a[9] / 3 + ( a[10] / 4 + a[11] / 5 * temp ) * temp ) * temp ) * temp + a[12] / temp ) * RGAS * temp / species[k].MolarMass;
	}
        //hr += h[j] * species[k].MolarMass * reaction[i].ProductSt->vec[j];
	hr += h[j] * reaction[i].ReactionRate * reaction[i].ProductSt->vec[j] * species[k].MolarMass;
    }
    //MOD_XY
    //hr1 *= reaction[i].ReactionRate;
    //reaction[i].HeatRelease = hr1;
    //TMP_XY
    //reaction[i].HeatRelease = -hr;
    reaction[i].HeatRelease = -hr;
    free(h);
    //TMP_XY
    //printf("Leaving ComputeHeatRelease.\n");
}

double myLog10( double arg)
{
    static double small = 1.e-200;
    static double value = -200.0;

    return ( arg > small ) ? log10( arg ) : value;
}
