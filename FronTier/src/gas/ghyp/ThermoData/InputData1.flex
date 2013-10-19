%option	noyywrap
%top{
	#include "st.h"
	#include <math.h>
	#include <string.h>
}
	//This structure is meant to hold the data within the file containing the elemental molar mass.
	struct elementMM{
		char symbol;
		double molMass;
	}mmData[NUM_OF_ELEMENTS];

	int isNewSpecies(void);
	void computeMolarMass(void);
	void foundSpecies(void);
	void readNASAData(void);

	
	int reactionIndex =0;  				//Keeps track of reaction number.  In the end, should be equal to the total number of reactions.
	int productFlag=0;				//productFlag=0 if the program is working on the reactant side of the equation.
	char *specTok;					//Points to the species name after the array of characters has been tokenized.
	int specTokLen;
	int specIndex =0;				//Keeps track of species number.  In the end, should be equal to the total number of species.
	int reactantCount =0;				//Counts the number of reactants in each reaction.
	int productCount =0;				//Counts the number of products in each reaction.
	char *specList[NUM_OF_SPECIES];			//List of the species number and their symbols.
	int prevReactionFound[NUM_OF_SPECIES];		//Contains the previous reaction number in which a particular species was found.
	double A =0;


	
%s scPreExp
%s scPowTen
%s scTemExp
%s scActEn
%s scReactSide
%s scMoreSpec

%%
"->"	{
		productFlag=1;
	}

"[OTHER]"	;

" "[HCO0-9]+" +"	{
				foundSpecies();
				BEGIN(scMoreSpec);
			}

<scMoreSpec>" "[HCO0-9]+" "		{
						foundSpecies();
						BEGIN(0);
					}

"a = "	{
		BEGIN(scPreExp);
	}

<scPreExp>[0-9]+"."[0-9]+	A=atof(yytext); BEGIN(0);

"E+"	BEGIN(scPowTen);

<scPowTen>[0-9]+	 {	
				reaction[reactionIndex].PreExp = A*pow(10.0,atof(yytext));
				BEGIN(0);
			 }

"n = "			BEGIN(scTemExp);

<scTemExp>-?[0-9]+	{
				reaction[reactionIndex].TempExp = atof(yytext);
				BEGIN(0);
			}

<scTemExp>-?[0-9]+"."[0-9]+	{
					reaction[reactionIndex].TempExp = atof(yytext);
					BEGIN(0);
				}

"E ="			BEGIN(scActEn);

<scActEn>-?[0-9]+	{	
				reaction[reactionIndex].ActEnergy = atof(yytext);
				BEGIN(0);
			}

<scActEn>-?[0-9]+"."[0-9]+	{
					reaction[reactionIndex].ActEnergy = atof(yytext);
					BEGIN(0);
				}
"}"	{
		++reactionIndex;
		reactantCount=0;
		productCount=0;
		productFlag=0;
	}

.	;
%%
//Allocates memory for data structures.
void initDataStructures(){
	int i,j,k;
	for(i=0;i<NUM_OF_SPECIES;i++){
		species[i].MolarMass = 0.0;
		species[i].MassFraction = 0.0;
		species[i].ProductionRate = 0.0;
		species[i].NofUsedReactions = 0;
		prevReactionFound[i] = -1;	//this number needs to initialize to an impossible index.
		
		species[i].SSt = (Vector *)malloc(sizeof(Vector));
		species[i].SSt->vec = (double *)malloc(NUM_OF_REACTIONS*sizeof(double));
		species[i].SSt->len=0;
		
		species[i].NASACoeff = (Vector *)malloc(sizeof(Vector));
		species[i].NASACoeff->vec = (double *)malloc(7*sizeof(double));
		species[i].NASACoeff->len=0;	
	
		species[i].UsedReactions = (IntVector *)malloc(sizeof(IntVector));
		species[i].UsedReactions->vec = (int *)malloc(NUM_OF_REACTIONS*sizeof(int));
		species[i].UsedReactions->len=0;
	}
	for(j=0;j<NUM_OF_REACTIONS;j++){
		reaction[j].NofReactants = 0;
		reaction[j].NofProducts = 0;
		reaction[j].ReactionRate = 0.0;
		reaction[j].RateCoeff = 0.0;
		reaction[j].PreExp = 0.0;
		reaction[j].TempExp = 0.0;
		reaction[j].ActEnergy = 0.0;

		reaction[j].ReactantIdx = (IntVector *)malloc(sizeof(IntVector));
		reaction[j].ReactantIdx->vec = (int *)malloc(NUM_OF_SPECIES*sizeof(int));
		reaction[j].ReactantIdx->len = 0;
		
		reaction[j].ReactantSt = (Vector *)malloc(sizeof(Vector));
		reaction[j].ReactantSt->vec = (double *)malloc(NUM_OF_SPECIES*sizeof(double));
		reaction[j].ReactantSt->len = 0;
		
		reaction[j].ProductIdx = (IntVector *)malloc(sizeof(IntVector));
		reaction[j].ProductIdx->vec = (int *)malloc(NUM_OF_SPECIES*sizeof(int));
		reaction[j].ProductIdx->len = 0;
		
		reaction[j].ProductSt = (Vector *)malloc(sizeof(Vector));
		reaction[j].ProductSt->vec = (double *)malloc(NUM_OF_SPECIES*sizeof(double));
		reaction[j].ProductSt->len = 0;
	
		for(k=0;k<NUM_OF_SPECIES;k++){
			reaction[j].ReactantSt->vec[k]=0.0;
			reaction[j].ProductSt->vec[k]=0.0;
		}
	}
	
}

void freeDataStructures(){
	int i,j;
	for(i=0;i<NUM_OF_SPECIES;i++){
		free(species[i].NASACoeff->vec);
		free(species[i].NASACoeff);
		free(species[i].UsedReactions->vec);
		free(species[i].UsedReactions);
		free(species[i].SSt->vec);
		free(species[i].SSt);
	}
	for(j=0;j<NUM_OF_REACTIONS;j++){
		free(reaction[j].ReactantIdx->vec);
		free(reaction[j].ReactantIdx);
		free(reaction[j].ReactantSt->vec);
		free(reaction[j].ReactantSt);
		
		free(reaction[j].ProductIdx->vec);
		free(reaction[j].ProductIdx);
		free(reaction[j].ProductSt->vec);
		free(reaction[j].ProductSt);
		
		free(reaction[j].ReactantMM->vec);
		free(reaction[j].ReactantMM);
		free(reaction[j].ProductMM->vec);
		free(reaction[j].ProductMM);
	}
}
/*Everytime a species is found in the .mech file this function processes it by putting the necessary information in the appropriate data structures.*/
void foundSpecies(){
	specTok = strtok(yytext,"+ ");
	specTokLen = strlen(specTok);
	int curSpecIndex = isNewSpecies();
	int l;
	if((prevReactionFound[curSpecIndex] != reactionIndex)){  //Check if this species occurs multiple times in the same reaction.
		++species[curSpecIndex].NofUsedReactions;
		int l = species[curSpecIndex].UsedReactions->len;
		species[curSpecIndex].UsedReactions->vec[l]=reactionIndex;
		++species[curSpecIndex].UsedReactions->len;
		prevReactionFound[curSpecIndex]=reactionIndex;
	
		if(productFlag==0){
			++reaction[reactionIndex].NofReactants;
			l = reaction[reactionIndex].ReactantIdx->len;
			reaction[reactionIndex].ReactantIdx->vec[l]=curSpecIndex;			
			++reaction[reactionIndex].ReactantSt->vec[l];		
	
			++reaction[reactionIndex].ReactantIdx->len;
			++reaction[reactionIndex].ReactantSt->len;
		}else{
			++reaction[reactionIndex].NofProducts;
			l = reaction[reactionIndex].ProductIdx->len;
			reaction[reactionIndex].ProductIdx->vec[l]=curSpecIndex;
			++reaction[reactionIndex].ProductSt->vec[l];
			
			++reaction[reactionIndex].ProductIdx->len;
			++reaction[reactionIndex].ProductSt->len;
		}
	}else{
		if(productFlag==0){
			l = reaction[reactionIndex].ReactantSt->len;
			++reaction[reactionIndex].ReactantSt->vec[l-1];
		}else{
			l = reaction[reactionIndex].ProductSt->len;
			++reaction[reactionIndex].ProductSt->vec[l-1];
		}
	}
	
}
/*This function returns the index of the species that was found in the .mech file.*/
int isNewSpecies(){	
	int  i, result=1;
	for(i =0;i<specIndex;i++){
		result = strcmp(specTok,specList[i]);
		if(result==0){
			return i;
		}
	}	
	if(result != 0){			//New species has been found. Add it to list and find molar mass.
		specList[specIndex]=specTok;
		computeMolarMass();

		//printf("the species %s was added to position %d\n",specList[specIndex],specIndex);
		++specIndex;
		return specIndex-1;
	}
	printf("ERROR in isNewSpecies()\n");
}

/*Computes the molar mass of the species found based on the data in the elemental molar mass file and the symbol name found in the .mech file.*/
void computeMolarMass(){
	double mm = 0.0;
	double temp = 0.0;
	int startIndex = (specTok[0]<='9')&&(specTok[0]>='1');
	
	int i,j;
	
	for(j=startIndex;j<NUM_OF_ELEMENTS;j++){	
		if(specTok[startIndex]==mmData[j].symbol){
			temp = mmData[j].molMass;
		}	
	}
	for(i=startIndex+1;i<specTokLen;i++){
		if( (specTok[i]>='A') && (specTok[i]<='Z') ){
			mm = mm + temp;
			for(j=0;j<NUM_OF_ELEMENTS;j++){
				if(specTok[i]==mmData[j].symbol){
					temp = mmData[j].molMass;
				}
			}
		}else{
			temp *= (specTok[i]-'0');
		}
	}
	species[specIndex].MolarMass = mm+temp;
	return;
}

/*Read the elemental molar mass data file and store the data in the struct mmData.*/
void readMolarMassData()
{
	FILE *fp;
	fp = fopen("MolarMass","r");
	char buffer[256];
	char *result = NULL;
	char delims[] = " ";
	int index=0;

	fgets(buffer,256,fp);
	result = strtok(buffer,delims);
	while( result != NULL){
		if( (result[0]>='A') && (result[0]<='Z') ){
			mmData[index].symbol = result[0];
		}else{
			mmData[index].molMass = atof(result);
			++index;
		}
		result = strtok(NULL,delims);
	}		
	fclose(fp);
	return;	
}

/*Continue to put the approprate molar mass data in the data structures.*/
void putMMinReaction(){
	int i,j;

	for(i=0;i<NUM_OF_REACTIONS;i++){
		
		int len;

		//Molar Mass of reactants.
		len = reaction[i].ReactantIdx->len;
		reaction[i].ReactantMM = (Vector *)malloc(sizeof(Vector));
		reaction[i].ReactantMM->vec = (double *)malloc(len*sizeof(double));
		reaction[i].ReactantMM->len = len;
		for(j=0;j<len;j++){
			reaction[i].ReactantMM->vec[j] = species[reaction[i].ReactantIdx->vec[j]].MolarMass;
		}
	
		//Molar Mass of products.
		len = reaction[i].ProductIdx->len;
		reaction[i].ProductMM = (Vector *)malloc(sizeof(Vector));
		reaction[i].ProductMM->vec = (double *)malloc(len*sizeof(double));
		reaction[i].ProductMM->len = len;
		for(j=0;j<len;j++){
			reaction[i].ProductMM->vec[j] = species[reaction[i].ProductIdx->vec[j]].MolarMass;
		}
	}
}

/*Read the data from the NASA tables.  This routine is dependent on the specific format of the tables provided by NASA.*/
void readNASAData(){
	FILE *fp;
	int maxline = 82;
	fp = fopen("H2Thermo","r");
	char buffer[maxline];
	char *result = NULL;
	char delims[] = " \n";
	int i,j,l,offset;
	char base[11];
	char exp[3];
	double coeff;
	int foundIndex;
	
	fgets(buffer,maxline,fp);
	result = strtok(buffer,delims);
	
	if(strcmp(result,"THERMO") !=0){
		printf("ERROR: reading NASA data.\n");
		return;
	}else{
		while(fgets(buffer,maxline,fp) != NULL){
			if(buffer[79]=='1'){
				result = strtok(buffer,delims);
				for(i=0;i<NUM_OF_SPECIES;i++){
					if(strcmp(result,specList[i])==0){
						foundIndex=i;
						break;
					}
				}
			}else if(buffer[79]=='2'){
				offset=0;
				while(offset <= 60){
					for(j=offset;j<offset+11;j++){
						base[j-offset] = buffer[j];
					}
					for(j=offset;j<offset+3;j++){
						exp[j-offset]=buffer[12+j];
					}
						
					l = species[foundIndex].NASACoeff->len;
					species[foundIndex].NASACoeff->vec[l] = atof(base)*pow(10,atof(exp));
					++species[foundIndex].NASACoeff->len;
					
					offset += 15;
				}
			}else if(buffer[79]=='3'){
				offset=0;
				while(offset<=15){
					for(j=offset;j<offset+11;j++){
						base[j-offset] = buffer[j];
					}
					for(j=offset;j<offset+3;j++){
						exp[j-offset]=buffer[12+j];
					}

					l = species[foundIndex].NASACoeff->len;
					species[foundIndex].NASACoeff->vec[l] = atof(base)*pow(10,atof(exp));
					++species[foundIndex].NASACoeff->len;
					
					offset += 15;
				}
			}else{
				continue;
			}
		}
	}
}

void main(int argc, char *argv[])
{	
	if (argc >0){
		yyin = fopen(argv[1],"r");
		initDataStructures();
		printf("initDataStructures SUCCESSFUL!\n");
		readMolarMassData();
		printf("readMolarMassData SUCCESSFUL!\n");
		yylex();
		printf("yylex SUCCESSFUL!\n");
		readNASAData();
		printf("readNASAData SUCCESSFUL!\n");
		putMMinReaction();
		printf("putMMinReaction SUCCESSFUL!\n");
		
		int i;
		for(i=0;i<NUM_OF_SPECIES;i++){
			printf("Species Molar Mass: %f\n",species[i].MolarMass);
		}
		freeDataStructures();
		printf("freeDataStructures SUCCESSFUL!\n");
	}else{
		printf("ERROR: Cannot read molar mass data.\n");
		yyin = stdin;
	}
}
