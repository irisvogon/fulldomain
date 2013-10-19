/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*			navigate.c
*
*	Copywrite 1999 by The University at Stony Brook, All rights reserved.
*
*	
*	Author:  	John D. Pinezich,
*			Department of Applied Mathematics and Statistics,
*			The University at Stony Brook
*/

#if defined(NAVIGATOR)

#include <navdecs.h>

#define   MAX_NAV_NAMES  500 /* Max Number of debugging names         */
#define   MAX_NAV_CHARS   60 /* Max Number of characters used per name*/

struct _NAVIGATE_PARAMS {
        FILE	*_navigate_input;            
	int	_num_of_navigate_names;      
  	int	_psuppress[MAX_NAV_NAMES];                 
  	int	_psuppress_cnt[MAX_NAV_NAMES];                 
  	int	_pbriefly[MAX_NAV_NAMES];                    
  	int	_pbriefly_cnt[MAX_NAV_NAMES];                    
  	int	_psummary[MAX_NAV_NAMES];                    
  	int	_psummary_cnt[MAX_NAV_NAMES];                    
  	int	_pdetails[MAX_NAV_NAMES];                    
  	int	_pdetails_cnt[MAX_NAV_NAMES];                    
  	int	_pforce[MAX_NAV_NAMES];                    
  	int	_pforce_cnt[MAX_NAV_NAMES];                    
  	int	_pexplain[MAX_NAV_NAMES];                    
  	int	_pexplain_cnt[MAX_NAV_NAMES];                    
  	int	_pcall_count[MAX_NAV_NAMES]; 
  	float	_pfloat0[MAX_NAV_NAMES][3];
  	float	_pfloat1[MAX_NAV_NAMES][3];
  	float	_pfloat2[MAX_NAV_NAMES][3];
  	POINTER	_ppointer[MAX_NAV_NAMES][3];
  	int	_pinteger0[MAX_NAV_NAMES][3];
  	int	_pinteger1[MAX_NAV_NAMES][3];
  	int	_pint0[MAX_NAV_NAMES];
  	int	_pint1[MAX_NAV_NAMES];
  	int	_pdepth_limiting[MAX_NAV_NAMES];
  	int	_plimit_value[MAX_NAV_NAMES];
  	int	_plimit_cnt[MAX_NAV_NAMES];
  	int	_pnavigate[MAX_NAV_NAMES];
  	int	_pnavigate_cnt[MAX_NAV_NAMES];
  	char	_pstring_bug[MAX_NAV_NAMES][MAX_NAV_CHARS];
  	char	_padd_to_debug_string[MAX_NAV_NAMES][MAX_NAV_CHARS];
  	int	_padd_to_debug[MAX_NAV_NAMES];
  	int	_padd_to_debug_cnt[MAX_NAV_NAMES];
  	char	_premove_from_debug_string[MAX_NAV_NAMES][MAX_NAV_CHARS];
	int	_premove_from_debug[MAX_NAV_NAMES];
  	int	_premove_from_debug_cnt[MAX_NAV_NAMES];
	char 	_navigate_names[MAX_NAV_NAMES][MAX_NAV_CHARS+1];

};		/* Stores Navigator data */

typedef struct _NAVIGATE_PARAMS NAVIGATE_PARAMS;

	/* Local variables */

LOCAL   NAVIGATE_PARAMS  NavParams;
LOCAL   int nav_loaded = NO;

#define navigate_input	        NavParams._navigate_input
#define	num_of_navigate_names	NavParams._num_of_navigate_names
#define	navigate_names		NavParams._navigate_names
#define psuppress		NavParams._psuppress
#define psuppress_cnt		NavParams._psuppress_cnt
#define pbriefly		NavParams._pbriefly
#define pbriefly_cnt		NavParams._pbriefly_cnt
#define psummary		NavParams._psummary
#define psummary_cnt		NavParams._psummary_cnt
#define pdetails		NavParams._pdetails
#define pdetails_cnt		NavParams._pdetails_cnt
#define pforce			NavParams._pforce
#define pforce_cnt		NavParams._pforce_cnt
#define pexplain		NavParams._pexplain
#define pexplain_cnt		NavParams._pexplain_cnt
#define pcall_count		NavParams._pcall_count
#define pfloat0			NavParams._pfloat0
#define pfloat1			NavParams._pfloat1
#define pfloat2			NavParams._pfloat2
#define pinteger0		NavParams._pinteger0
#define pinteger1		NavParams._pinteger1
#define pint0			NavParams._pint0
#define pint1			NavParams._pint1
#define ppointer		NavParams._ppointer
#define plimit_value		NavParams._plimit_value
#define plimit_cnt		NavParams._plimit_cnt
#define pnavigate		NavParams._pnavigate
#define pnavigate_cnt		NavParams._pnavigate_cnt
#define pdepth_limiting		NavParams._pdepth_limiting
#define pstring_bug		NavParams._pstring_bug
#define padd_to_debug_string		NavParams._padd_to_debug_string
#define premove_from_debug_string 	NavParams._premove_from_debug_string
#define padd_to_debug			NavParams._padd_to_debug
#define premove_from_debug 	NavParams._premove_from_debug
#define padd_to_debug_cnt	NavParams._padd_to_debug_cnt
#define premove_from_debug_cnt 	NavParams._premove_from_debug_cnt

	/* LOCAL Function Prototypes */

LOCAL   void  white_space(int);
LOCAL	int   navigating(char*);

	/* NAVIGATOR Functions */

EXPORT	void 	read_navigator_input(
	const char *fname)
{
	int 	i,j,k;
	char 	s[Gets_BUF_SIZE];
	char 	line[MAX_NAV_CHARS+1];
	char 	*c;
	int 	integer;
	float 	flote;
	POINTER pntr;

	if (nav_loaded)
	{
	    (void) printf("\nNavigator data already loaded.\n");
	    return;
        }

	num_of_navigate_names = 0;

	if (fname != NULL)
        {
	    if ((navigate_input = fopen(fname,"r")) == NULL)
	    {
	        (void) printf("WARNING in read_navigator_input(), "
			      "navigate file %s could not be opened\n",fname);
		clean_up(ERROR);
	    }
	}
	else navigate_input = stdin;

	(void) printf("List Functions to navigate on--end  ends list.\n"
		      "Optional Keywords, one to each line following function "
		      "name and to\n"
		      "which keyword is bound -> BRIEFLY, SUMMARY, DETAILS,\n"
		      "EXPLAIN, SUPPRESS, INT0, INT1, INTEGER0, INTEGER1, "
		      "STRING_BUG,\n"
		      "NAVIGATE, LIMIT_VALUE, FLOAT0, FLOAT1, FLOAT2\n\n");

	/* read data */

	fgets(line,MAX_NAV_CHARS+1,navigate_input);
	(void) sscanf(line,"%s",s);

	for (i = 0; ; i++)
        {    
	    if (strcmp(s,"end") == 0)
            {
		(void) printf("\t: end\n\n");
		break;
            }

	    (void) printf("\t: %s\n",s); 
	    if (navigating(s))
            {
		screen("ERROR in read_navigator_input(), "
		       "already navigating on %s.\n",s);
		clean_up(ERROR);
            }

	    strncpy(navigate_names[num_of_navigate_names],s,MAX_NAV_CHARS);

	    psuppress[num_of_navigate_names] = 0;
	    psuppress_cnt[num_of_navigate_names] = -1;
	    pbriefly[num_of_navigate_names] = 0;
	    pbriefly_cnt[num_of_navigate_names] = -1;
	    psummary[num_of_navigate_names] = 0;
	    psummary_cnt[num_of_navigate_names] = -1;
	    pdetails[num_of_navigate_names] = 0;
	    pdetails_cnt[num_of_navigate_names] = -1;
	    pexplain[num_of_navigate_names] = 0;
	    pexplain_cnt[num_of_navigate_names] = -1;
	    pforce[num_of_navigate_names] = 0;
	    pforce_cnt[num_of_navigate_names] = -1;
	    pcall_count[num_of_navigate_names] = 0;
	    pint0[num_of_navigate_names] = -1;
	    pint1[num_of_navigate_names] = -1;
	    for (k = 0; k < 3; k++)
	    {
	        pinteger0[num_of_navigate_names][k] = -1;
		pinteger1[num_of_navigate_names][k] = -1;
		pfloat0[num_of_navigate_names][k] = -1;
		pfloat1[num_of_navigate_names][k] = -1;
		pfloat2[num_of_navigate_names][k] = -1;
		ppointer[num_of_navigate_names][k] = NULL;
	    }
	    pdepth_limiting[num_of_navigate_names] = 0;
	    plimit_value[num_of_navigate_names] = 1000000;
	    plimit_cnt[num_of_navigate_names] = -1;
	    pnavigate[num_of_navigate_names] = 0;
	    pnavigate_cnt[num_of_navigate_names] = -1;
	    strcpy(pstring_bug[num_of_navigate_names],"");
	    strcpy(padd_to_debug_string[num_of_navigate_names],"");
	    strcpy(premove_from_debug_string[num_of_navigate_names],"");
	    padd_to_debug[num_of_navigate_names] = 0;
	    premove_from_debug[num_of_navigate_names] = 0;
	    padd_to_debug_cnt[num_of_navigate_names] = -1;
	    premove_from_debug_cnt[num_of_navigate_names] = -1;

	    for (j = 0; ; j++)
	    {
		fgets(line,MAX_NAV_CHARS+1,navigate_input);
		sscanf(line,"%s",s);

		if (strcmp(s,"SUPPRESS") == 0)
	        {
		    (void) printf("\t\t: SUPPRESS ");
		    psuppress[num_of_navigate_names] = 1;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			psuppress_cnt[num_of_navigate_names] = integer; 
			(void) printf("%d\n",integer);
	            }
	            else
		        (void) printf("\n");
	        }

		else if (strcmp(s,"BRIEFLY") == 0)
	        {
		    (void) printf("\t\t: BRIEFLY ");
		    pbriefly[num_of_navigate_names] = 1;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pbriefly_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
		        (void) printf("\n");
	        }

		else if (strcmp(s,"SUMMARY") == 0)
	        {
		    (void) printf("\t\t: SUMMARY ");
		    psummary[num_of_navigate_names] = 1;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			psummary_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
		        (void) printf("\n");
	        }

		else if (strcmp(s,"DETAILS") == 0)
	        {
		    (void) printf("\t\t: DETAILS ");
		    pdetails[num_of_navigate_names] = 1;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pdetails_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
		        (void) printf("\n");
	        }

		else if (strcmp(s,"EXPLAIN") == 0)
	        {
		    (void) printf("\t\t: EXPLAIN ");
		    pexplain[num_of_navigate_names] = 1;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pexplain_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
		        (void) printf("\n");
	        }

		else if (strcmp(s,"FORCE") == 0)
	        {
		    (void) printf("\t\t: FORCE ");
		    pforce[num_of_navigate_names] = 1;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
		        pforce_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
		    }
		    else
		        (void) printf("\n");
		}

		else if (strcmp(s,"FLOAT0") == 0)
	        {
		    (void) printf("\t\t: FLOAT0 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
		    {
			(void) sscan_float(c,&flote);
			(void) printf("%g ",flote);
			pfloat0[num_of_navigate_names][k++] = flote;
		    }
		    (void) printf("\n");
		}

		else if (strcmp(s,"FLOAT1") == 0)
	        {
		    (void) printf("\t\t: FLOAT1 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			(void) sscan_float(c,&flote);
			(void) printf("%g ",flote);
			pfloat1[num_of_navigate_names][k++] = flote;
	            }
		    (void) printf("\n");
	        }

		else if (strcmp(s,"FLOAT2") == 0)
	        {
		    (void) printf("\t\t: FLOAT2 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			(void) sscan_float(c,&flote);
			(void) printf("%g ",flote);
			pfloat2[num_of_navigate_names][k++] = flote;
	            }
		    (void) printf("\n");
	        }

		else if (strcmp(s,"INT0") == 0)
	        {
		    (void) printf("\t\t: INT0 ");
		    if (sscanf(line,"%*s %d",&integer) == 1)
		      pint0[num_of_navigate_names] = integer;
		    (void) printf("%d\n",pint0[num_of_navigate_names]);
	        }

		else if (strcmp(s,"INT1") == 0)
	        {
		    (void) printf("\t\t: INT1 ");
		    if (sscanf(line,"%*s %d",&integer) == 1)
		      pint1[num_of_navigate_names] = integer;
		    (void) printf("%d\n",pint1[num_of_navigate_names]);
	        }

		else if (strcmp(s,"INTEGER0") == 0)
	        {
		    (void) printf("\t\t: INTEGER0 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			sscanf(c,"%d",&integer);
			(void) printf("%d ",integer);
			pinteger0[num_of_navigate_names][k++] = integer;
	            }
		    (void) printf("\n");
	        }

		else if (strcmp(s,"INTEGER1") == 0)
	        {
		    (void) printf("\t\t: INTEGER1 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			sscanf(c,"%d",&integer);
			(void) printf("%d ",integer);
			pinteger1[num_of_navigate_names][k++] = integer;
	            }
		    (void) printf("\n");
	        }

		else if (strcmp(s,"POINTER") == 0)
	        {
		    (void) printf("\t\t: POINTER ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			sscanf(c,"%p",&pntr);
			(void) printf("%p ",pntr);
			ppointer[num_of_navigate_names][k++] = pntr;
	            }
		    (void) printf("\n");
	        }

		else if (strcmp(s,"LIMIT_VALUE") == 0)
	        {
		    sscanf(line,"%*s %d ",&integer);
		    (void) printf("\t\t: LIMIT_VALUE %d ",integer);
		    pdepth_limiting[num_of_navigate_names] = 1;
		    plimit_value[num_of_navigate_names] = integer;
		    if (sscanf(line,"%*s %*d %d",&integer) == 1)
	            {
			plimit_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
		        (void) printf("\n");
	        }
		
		else if (strcmp(s,"NAVIGATE") == 0)
		{
		    (void) printf("\t\t: NAVIGATE ");
		    pnavigate[num_of_navigate_names] = 1;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pnavigate_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else (void) printf("\n");
	        }

		else if (strcmp(s,"STRING_BUG") == 0)
	        {
		    sscanf(line,"%*s %s",s);
		    (void) printf("\t\t: STRING_BUG %s\n",s);
		    strcpy(pstring_bug[num_of_navigate_names],s);
	        }

		else if (strcmp(s,"ADD_TO_DEBUG") == 0)
	        {
		    sscanf(line,"%*s %s",s);
		    (void) printf("\t\t: ADD_TO_DEBUG %s ",s);
		    padd_to_debug[num_of_navigate_names] = 1;
		    strcpy(padd_to_debug_string[num_of_navigate_names],s);
		    if (sscanf(line,"%*s %*s %d",&integer) == 1)
	            {
	        	padd_to_debug_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
		        (void) printf("\n");
	        }

		else if (strcmp(s,"REMOVE_FROM_DEBUG") == 0)
	        {
		    sscanf(line,"%*s %s",s);
		    (void) printf("\t\t: REMOVE_FROM_DEBUG %s ",s);
		    premove_from_debug[num_of_navigate_names] = 1;
		    strcpy(premove_from_debug_string[num_of_navigate_names],s);
		    if (sscanf(line,"%*s %*s %d",&integer) == 1)
	            {
			premove_from_debug_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
		        (void) printf("\n");
	        }

		else
		   break;

	    }

	    if (((int) strlen(navigate_names[num_of_navigate_names])) >= 
		MAX_NAV_CHARS)
		navigate_names[num_of_navigate_names][MAX_NAV_CHARS] = '\0';

	    num_of_navigate_names++;
	    (void) printf("\n");
	}

	if (navigate_input != stdin)
	    fclose(navigate_input);
	fflush(stdout);
	nav_loaded = YES;
}		/*end read_navigator_input*/

LOCAL	int navigating(
	char *funcname)
{
	int i;

	for (i = 0; i < num_of_navigate_names; i++) 
	    if (strncmp(funcname,navigate_names[i],MAX_NAV_CHARS) == 0)
	        return 1;
	return 0;
}		/*end navigating*/

EXPORT 	void 	print_nav_data(
	const char		*fname)
{
  	int        i, j;
	static const char *yn[2] = {"N","Y"};

	(void) printf("\n");
	(void) printf("Printout of navigate data structure from %s()\n",fname);

	for (i = 0; i < num_of_navigate_names; i++)
	{
	    for (j = 0; j < 110; j++)
	        (void) printf("=");
	    (void) printf("\n\n");
	    (void) printf("\n%s()\n",navigate_names[i]);
	    (void) printf("\n  %-8s %1s | %-8s %1s | %-8s %1s | %-8s %1s | "
			  "%-8s %1s | %-8s %1s | %-8s %1s | %-8s %1d\n",
			  "NAVIGATE",yn[pnavigate[i]],
			  "SUPPRESS",yn[psuppress[i]],
			  "BRIEFLY",yn[pbriefly[i]],
			  "SUMMARY",yn[psummary[i]],
			  "DETAILS",yn[pdetails[i]],
			  "EXPLAIN",yn[pexplain[i]],
			  "FORCE",yn[pforce[i]],
			  "LIMIT_V",plimit_value[i]);
	    (void) printf("  ");

	    for (j = 0; j < 7; j++)
	        (void) printf("-----------|-");
	    (void) printf("-----------\n");
	    (void) printf("  %-10d | %-10d | %-10d | %-10d | %-10d | %-10d | "
			  "%-10d | %-10d\n\n",
			  pnavigate_cnt[i],psuppress_cnt[i],pbriefly_cnt[i],
			  psummary_cnt[i],pdetails_cnt[i],
			  pexplain_cnt[i],pforce_cnt[i],plimit_cnt[i]); 
	    (void) printf("\n  index = %d  depth_limiting = %1s  C_COUNT = %d  "
			  "STRING_BUG = %-s\n\n",i,
			  yn[pdepth_limiting[i]],pcall_count[i],pstring_bug[i]);
	    (void) printf("  INTEGER0 (%10d %10d %10d)  INT0 = %d\n",
			  pinteger0[i][0],pinteger0[i][1],pinteger0[i][2],
			  pint0[i]); 
	    (void) printf("  INTEGER1 (%10d %10d %10d)  INT1 = %d\n",
			  pinteger1[i][0],pinteger1[i][1],pinteger1[i][2],
			  pint1[i]); 
	    (void) printf("  POINTER  (%10p %10p %10p)\n",
			  ppointer[i][0],ppointer[i][1],ppointer[i][2]);
	    (void) printf("  FLOAT0   (%10g %10g %10g)\n",
			  pfloat0[i][0],pfloat0[i][1],pfloat0[i][2]);
	    (void) printf("  FLOAT1   (%10g %10g %10g)\n",
			  pfloat1[i][0],pfloat1[i][1],pfloat1[i][2]);
	    (void) printf("  FLOAT2   (%10g %10g %10g)\n\n",
			  pfloat2[i][0],pfloat2[i][1],pfloat2[i][2]);
	    (void) printf("  %1s ADD_TO_DEBUG      %10d  %-s\n",
			  yn[padd_to_debug[i]],padd_to_debug_cnt[i],
			  padd_to_debug_string[i]);
	    (void) printf("  %1s REMOVE_FROM_DEBUG %10d  %-s\n",
			  yn[premove_from_debug[i]],premove_from_debug_cnt[i],
			  premove_from_debug_string[i]);
	}
	(void) printf("\n");
	for (j = 0; j < 110; j++)
	    (void) printf("=");
	(void) printf("\n\n\n");

	(void) nav_update(pPRINT_STACK,fname,0,0);
	(void) printf("End Printout of navigate data structure from %s()\n\n",
		      fname);
}                /*end print_nav_data*/

LOCAL 	void 	white_space(
	int 		num_spaces)
{
	int i;
  	for (i = 0; i < num_spaces; i++)
	    (void) printf(" ");
}    		/*end white_space*/

EXPORT 	int 	nav_switches(
	int		command,
	int		index) 
{
  	if (index == -1)
	    return 0;
	switch(command)
	{
	case pBRIEFLY:
	    return pbriefly[index];
	case pSUMMARY:
	    return psummary[index];
	case pDETAILS:
	    return pdetails[index];
	case pEXPLAIN:
	    return pexplain[index];
	case pNAVIGATE:
	    return pnavigate[index];
	case pDEPTH_LIMITING:
	    return pdepth_limiting[index];
	case pFORCE:
	    return pforce[index];
	case pADD_TO_DEBUG:
	    return padd_to_debug[index];
	case pREMOVE_FROM_DEBUG:
	    return premove_from_debug[index];
	case pSUPPRESS:
	    return psuppress[index];
	default:
	    (void) printf("WARNING in nav_switches(), unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	}
	return 0;
} 		/*end nav_switches*/

EXPORT 	int	nav_integers(
	int		command,
	int		index)
{
  	if (index == -1)
	    return -1;
	switch (command)
	{
	case pBRIEFLY_CNT:
	    return pbriefly_cnt[index];
	case pSUMMARY_CNT:
	    return psummary_cnt[index];
	case pDETAILS_CNT:
	    return pdetails_cnt[index];
	case pEXPLAIN_CNT:
	    return pexplain_cnt[index];
	case pFORCE_CNT:
	    return pforce_cnt[index];
	case pNAVIGATE_CNT:
	    return pnavigate_cnt[index];
	case pADD_TO_DEBUG_CNT:
	    return padd_to_debug_cnt[index];
	case pREMOVE_FROM_DEBUG_CNT:
	    return premove_from_debug_cnt[index];
	case pC_COUNT:
	    return pcall_count[index];
	case pLIMIT_VALUE:
	    return plimit_value[index];
	case pLIMIT_CNT :
	    return plimit_cnt[index];
	case pSUPPRESS_CNT:
	    return psuppress_cnt[index];
	case pINT0:
	    return pint0[index];
	case pINT1:
	    return pint1[index];
	default:
	    (void) printf("WARNING in nav_integers(), unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	}
	return 0;
}		/*end nav_integers*/

EXPORT 	int 	*nav_integer_pointers(
	int 		command,
	int		index)
{
  	switch (command)
        {
	case pINTEGER0:
	    return pinteger0[index];
	case pINTEGER1:
	    return pinteger1[index];
	default:
	    (void) printf("WARNING in nav_integer_pointers(), "
			  "unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	}
	return NULL;
}		/*end nav_integer_pointers*/

EXPORT 	float	*nav_float_pointers(
	int		command,
	int		index)
{
  	switch (command)
	{
	case pFLOAT0:
	    return pfloat0[index];
	case pFLOAT1:
	    return pfloat1[index];
	case pFLOAT2:
	    return pfloat2[index];
	default:
	    (void) printf("WARNING in nav_float_pointers(), "
			  "unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	}
	return NULL;
}		/*end nav_float_pointers*/

EXPORT 	POINTER	*nav_pointers(
	int		command,
	int		index)
{
  	switch (command)
	{
	case pPOINTER:
	    return ppointer[index];
	default:
	    (void) printf("WARNING in nav_pointers(), "
			  "unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	  }
	return NULL;
}		/*end nav_pointers*/

EXPORT	char	*nav_strings(
	int		command,
	int		index)
{
  	switch (command)
	{
	case pADD_TO_DEBUG_STRING:
	    return padd_to_debug_string[index];
	case pREMOVE_FROM_DEBUG_STRING:
	    return premove_from_debug_string[index];
	case pSTRING_BUG:
	    return pstring_bug[index];
	default:
	    (void) printf("WARNING in nav_strings(), "
			  "unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	}
	return NULL;
}		/*end nav_strings*/

EXPORT  int  depth_monitor(
	int	command)
{
	static int monitor = NO;
	switch (command)
	{
	case pSTART:
	case pSTOP:
	    monitor = command;
	case pLOOK:
	    return monitor;
	default:
	    (void) printf("WARNING in depth_monitor(), "
			  "unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	}
	return 0;
}		/*end depth_monitor*/

EXPORT  int  function_monitor(
	int 	command)
{
  	static int monitor = NO;
	switch (command)
	{
	case pSTART:
	case pSTOP:
	    monitor = command;
	case pLOOK:
	    return monitor;
	default:
	    (void) printf("WARNING in function_monitor(), "
			  "unknown command = %d\n",
			  command);
	    clean_up(ERROR);
	}
	return 0;
}		/*end function_monitor*/

EXPORT  int  nav_update(
	int	   command,
	const char *function_name,
	int	   limit_integer,
	int	   fname_index)
{
  	static int 	depth = 0;
	static int 	limit = 1000000;
	static int 	function_suppression = NO;
	static int 	reenable_depth = -1;
	static int 	reenable_index = -2;
	int 		old_limit,i;
	
	static char stack[MAX_NAV_NAMES][MAX_NAV_CHARS+1];

	switch (command)
	{
	case pINIT:
	    depth = limit_integer;
	    return 0;

	case pLIMIT:
	    old_limit = limit-depth;
	    limit = limit_integer+depth;
	    if (debugging("LIMIT_VALUE"))
	        (void) printf("LIMIT TRACE in %s( %d ), set LIMIT = %d "
			      "from lim val = %d, depth = %d; old lim = %d\n",
			      function_name,pcall_count[fname_index],limit,
			      limit_integer,depth,old_limit);
	    return old_limit;

	case pENTER:
	    if (fname_index >= 0)
	        pcall_count[fname_index]++;

	    if (!NAVIGATING)
	        return 0;

	    depth++;

	    if (depth >= 0)
	        strncpy(stack[depth],function_name,MAX_NAV_CHARS);

	    if (fname_index >= 0 && SUPPRESS && !function_suppression) 
	    {
	        function_monitor(pSTOP);
		reenable_depth = depth;
		reenable_index = fname_index;
		function_suppression = YES;
	    }

	    if (depth == limit+1)
	        depth_monitor(pSTOP);
	
	    if (LOOK)
	    {
		white_space(depth%20);
		(void) printf("\\%d %s( %d ) \n",
			      depth,function_name,pcall_count[fname_index]);
	    }

	    return 0;

	case pLEAVE:
	    if (!NAVIGATING)
	        return 0;

	    if (depth >= 0)
	    {
		if (strcmp(stack[depth],function_name) != 0) 
		    (void) printf("warning: stack misaligned, function_name "
				  "= %s(), look in stack[%d] = %s()\n",
				  function_name,depth,stack[depth]);
	    }

	    if (LOOK)
	    {
	        white_space(depth%20);
		(void) printf("/%d %s( %d )\n",
			      depth,function_name,pcall_count[fname_index]);
	    }

	    if (depth == limit+1)
	        depth_monitor(pSTART);

	    if (function_suppression && (reenable_depth == depth) &&
		(reenable_index == fname_index))
	    {
	        function_monitor(pSTART);
		reenable_depth = -1;
		reenable_index = -2;
		function_suppression = NO;
	    }

	    depth--;
	    return 0;

	case pSTATIC:
	    white_space(depth%20);
	    (void) printf(" %d>\n",depth);
	    return 0;

	case pSPACE:
	    white_space(depth%20);
	    (void) printf(" %d> ",depth);
	    return 0;

	case pPRINT_STACK:
	    if (depth < 1) 
	    { 
		(void) printf("printout of stack from %s(), depth = %d,  "
			      "stack is empty\n",function_name,depth);
		return 0; 
	    }
	    (void) printf("printout of stack from %s(), depth = %d  "
			  "limit = %d\n\n",function_name,depth,limit);
	    for (i = 1; i <= depth; i++)
	        (void) printf("\tstack[%3d] = %s()\n",i,stack[i]);
	    (void) printf("\nend of stack\n");
	    return 0;

	default:
	    (void) printf("Warning in nav_update(), command = %d unknown\n",
			  command);
	    clean_up(ERROR);
	}
	return 0;
} 		/*end nav_update*/

/*ARGSUSED*/
EXPORT	int	nav_trace(
	int		command,
	const char	*function_name,
	int		index)
{
  	static int 	trace = 0;
	static char 	start_function[MAX_NAV_CHARS+1] = "";
	static char 	stop_function[MAX_NAV_CHARS+1] = "";

	switch (command)
	{
	case pSTOP:
	case pSTART:
	    if (trace == command)
	    {
	        (void) printf("WARNING in nav_trace() command = %d "
			      "already set.\n",command);
		(void) printf("\tstart_function = %s()\n",start_function);
		(void) printf("\tstop_function  = %s()\n",stop_function);
	    }

	    if (command == pSTART)
	    {
		pindent;
		(void) printf("NAVIGATOR ON in %s\n",function_name);
		depth_monitor(pSTART);
		function_monitor(pSTART);
		strcpy(start_function,function_name);
	    }

	    if (command == pSTOP)
	    {
	        (void) printf("NAVIGATOR OFF in %s\n",function_name);
		depth_monitor(pSTOP);
		function_monitor(pSTOP);
		strcpy(stop_function,function_name);
	    }

	    trace = command;
	    return 0;

	case pLOOK:
	    return trace;

	default:
	    (void) printf("Unknown command = %d in nav_trace()\n",command);
	    clean_up(ERROR);
	}
	return 0;
}		/*end nav_trace*/

EXPORT	int	nav_index(
	const char *function_name)
{
  	int i;
	if (nav_loaded == NO)
	    return -2;
	for (i = 0; i < num_of_navigate_names; i++)
	    if (strncmp(function_name,navigate_names[i],MAX_NAV_CHARS) == 0)
	        return i;
	return -1;
}		/*end nav_index*/

EXPORT	void	nprint_long_string(
	const char *s)
{
	char	*c;
	int	len;
	int	maxlen = 79;
	int	newlen;
	int	tablen = 0;
	bool start_of_line;
	static	char	*line = NULL;
	static	size_t	allocated_length = 0;

	if (strlen(s) == 0)
	    return;

	if (strlen(s) >= allocated_length)
	{
	    if (line != NULL)
	    	free(line);
	    allocated_length = strlen(s) + 1024;
	    line = (char*) malloc(allocated_length*sizeof(char));
	}
	strcpy(line,s);

	for (len = 0, c = strtok(line," "); c != NULL; c = strtok(NULL," "))
        {
	    newlen = len+(int)strlen(c);
	    if ((newlen+1) >  maxlen)
            {
		screen("\n");
		start_of_line = YES;
		len = tablen;
		newlen = len+(int)strlen(c);
            }
	    else if (start_of_line == NO)
            {
		screen(" ");
		newlen++;
		len++;
            }

	    if (start_of_line) pindent;
	    screen("%s",c);
	    start_of_line = NO;

	    if (c[strlen(c)-1] == '\n') 
            {
		len = tablen;
		start_of_line = YES;
            }

	    else if (newlen == maxlen)
            {
		screen("\n");
		start_of_line = YES;
		len = tablen;
            }
	    else
            {
		if (c[strlen(c)-1] == '.')
		{
		    screen(" ");
		    len++;
	        }
		len += (int)strlen(c);
	    }
	}
}		/*end nprint_long_string*/


#if DONT_COMPILE /*SUGGESTED REVISION TO NAVIGATOR 19990423*/

enum {MAX_N_NAMES = 500}; /* Max Number of debugging names          */
enum {MAX_N_CHARS = 60};  /* Max Number of characters used per name */

/* commands for navigator_update() */

enum _N_UPDATE {
	pINIT 	        =  0,
	pENTER	        =  1,
	pLEAVE	        =  2,
	pSTATIC	        =  3,
	pSPACE	        =  4,
	pLIMIT          =  5,
	pPRINT_STACK    =  6,
	pEXTRA          =  7
};
typedef enum _N_UPDATE N_UPDATE;

/* commands for navigator_trace(), function_monitor(), depth_monitor() */

enum _N_MONITOR {
	pLOOK  = -1,
	pSTOP  =  0,
	pSTART =  1
};
typedef enum _N_MONITOR N_MONITOR;

/* commands for navigator_switches() */

enum _N_SWITCH {
	pBRIEFLY	   =  0,
	pSUMMARY	   =  1,
	pDETAILS	   =  2,
	pEXPLAIN	   =  3,
	pNAVIGATE	   =  4,
	pDEPTH_LIMITING    =  5,
	pFORCE		   =  6,
	pADD_TO_DEBUG	   =  7,
	pREMOVE_FROM_DEBUG =  8,
	pSUPPRESS	   =  9
};
typedef enum _N_SWITCH N_SWITCH;

/* commands for navigator_integers() */

enum _N_CNT {
	pBRIEFLY_CNT	       =   0,
	pSUMMARY_CNT	       =   1,
	pDETAILS_CNT	       =   2,
	pEXPLAIN_CNT	       =   3,
	pNAVIGATE_CNT	       =   4,
	pFORCE_CNT	       =   5,
	pADD_TO_DEBUG_CNT      =   6,
	pREMOVE_FROM_DEBUG_CNT =   7,
	pC_COUNT	       =   8,
	pLIMIT_VALUE	       =   9,
	pLIMIT_CNT	       =  10,
	pSUPPRESS_CNT	       =  11,
	pINT0		       =  12,
	pINT1		       =  13
};
typedef enum _N_CNT N_CNT;

/* commands for navigator_integer_pointers() */

enum _N_PINTEGER {
	pINTEGER0 = 0,
	pINTEGER1 = 11
};
typedef enum _N_PINTEGER N_PINTEGER;

/* commands for navigator_float_pointers() */

enum _N_PFLOAT {
	pFLOAT0 = 0,
	pFLOAT1 = 1,
	pFLOAT2 = 2
};
typedef enum _N_PFLOAT N_PFLOAT;

/* commands for navigator_pointers() */

enum _N_PPOINTER {
	pPOINTER = 0
};
typedef enum _N_PPOINTER N_PPOINTER;

/* commands for navigator_strings() */

enum _N_STRING {
	pADD_TO_DEBUG_STRING      = 0,
	pREMOVE_FROM_DEBUG_STRING = 1,
	pSTRING_BUG               = 2
};
typedef enum _N_STRING N_STRING;

enum _N_ENTER_LEAVE {
	ENTER_NAVIGATOR = 0,
	LEAVE_NAVIGATOR = 1
};
typedef enum _N_ENTER_LEAVE  N_ENTER_LEAVE;

struct _NAVIGATE_PARAMS {	/* Stores Navigator data */
        FILE		*_navigate_input;            
	int		_num_of_navigate_names;      

  	bool _pbriefly[MAX_N_NAMES];                    
  	bool _psummary[MAX_N_NAMES];                    
  	bool _pdetails[MAX_N_NAMES];                    
  	bool _pexplain[MAX_N_NAMES];                    
  	bool _pnavigate[MAX_N_NAMES];
  	bool _pdepth_limiting[MAX_N_NAMES];
  	bool _pforce[MAX_N_NAMES];                    
  	bool _padd_to_debug[MAX_N_NAMES];
	bool _premove_from_debug[MAX_N_NAMES];
  	bool _psuppress[MAX_N_NAMES];                 

  	int _pbriefly_cnt[MAX_N_NAMES];                    
  	int _psummary_cnt[MAX_N_NAMES];                    
  	int _pdetails_cnt[MAX_N_NAMES];                    
  	int _pexplain_cnt[MAX_N_NAMES];                    
  	int _pnavigate_cnt[MAX_N_NAMES];
  	int _pforce_cnt[MAX_N_NAMES];                    
  	int _padd_to_debug_cnt[MAX_N_NAMES];
  	int _premove_from_debug_cnt[MAX_N_NAMES];
  	int _pcall_count[MAX_N_NAMES]; 
  	int _plimit_value[MAX_N_NAMES];
  	int _plimit_cnt[MAX_N_NAMES];
  	int _psuppress_cnt[MAX_N_NAMES];                 
  	int _pint0[MAX_N_NAMES];
  	int _pint1[MAX_N_NAMES];

  	int _pinteger0[MAX_N_NAMES][3];
  	int _pinteger1[MAX_N_NAMES][3];

  	float _pfloat0[MAX_N_NAMES][3];
  	float _pfloat1[MAX_N_NAMES][3];
  	float _pfloat2[MAX_N_NAMES][3];

  	POINTER	 _ppointer[MAX_N_NAMES][3];

  	char _padd_to_debug_string[MAX_N_NAMES][MAX_N_CHARS];
  	char _premove_from_debug_string[MAX_N_NAMES][MAX_N_CHARS];
  	char _pstring_bug[MAX_N_NAMES][MAX_N_CHARS];

	char _navigate_names[MAX_N_NAMES][MAX_N_CHARS + 1];
};
typedef struct _NAVIGATE_PARAMS NAVIGATE_PARAMS;

	/* Local function prototypes */
LOCAL	N_MONITOR       depth_monitor(const N_MONITOR);
LOCAL	N_MONITOR       function_monitor(const N_MONITOR);
LOCAL	N_MONITOR       navigator_trace(const N_MONITOR,const char*);
LOCAL   NAVIGATE_PARAMS NavParams;
LOCAL	POINTER         *navigator_pointers(const N_PPOINTER,const int);
LOCAL	bool         ADD_TO_DEBUG(const char*);
LOCAL	bool         DEPTH_LIMITING(const char*);
LOCAL	bool         LOOK(void);
LOCAL	bool         NAVIGATING(void);
LOCAL	bool         REMOVE_FROM_DEBUG(const char*);
LOCAL	bool         SUPPRESS(const char*);
LOCAL	bool         navigator_switches(const N_SWITCH,const int);
LOCAL	char	        *ADD_TO_DEBUG_STRING(const char*);
LOCAL	char            *navigator_strings(const N_STRING,const int);
LOCAL	char	        *REMOVE_FROM_DEBUG_STRING(const char*);
LOCAL	const char	*YN(bool);
LOCAL	float           *navigator_float_pointers(const N_PFLOAT,const int);
LOCAL	int             *navigator_integer_pointers(const N_PINTEGER,const int);
LOCAL	int             navigator_index(const char*);
LOCAL	int             navigator_integers(const N_CNT,const int);
LOCAL	int             navigator_update(const N_UPDATE,const char*,
				         const int,const int);
LOCAL	void	        CHECK_ADD_TO_DEBUG(const char*);
LOCAL	void	        CHECK_REMOVE_FROM_DEBUG(const char*);
LOCAL	void	        NAVIGATE_ENTER_LEAVE(const char*,const N_ENTER_LEAVE);


	/* Local variables*/
LOCAL   int navigator_loaded = NO;

#define navigate_input	        NavParams._navigate_input
#define	num_of_navigate_names	NavParams._num_of_navigate_names
#define	navigate_names		NavParams._navigate_names
#define psuppress		NavParams._psuppress   	
#define psuppress_cnt		NavParams._psuppress_cnt
#define pbriefly		NavParams._pbriefly		
#define pbriefly_cnt		NavParams._pbriefly_cnt		
#define psummary		NavParams._psummary		
#define psummary_cnt		NavParams._psummary_cnt		
#define pdetails		NavParams._pdetails		
#define pdetails_cnt		NavParams._pdetails_cnt		
#define pforce			NavParams._pforce	
#define pforce_cnt		NavParams._pforce_cnt		
#define pexplain		NavParams._pexplain		
#define pexplain_cnt		NavParams._pexplain_cnt		
#define pcall_count		NavParams._pcall_count
#define pfloat0			NavParams._pfloat0
#define pfloat1			NavParams._pfloat1
#define pfloat2			NavParams._pfloat2
#define pinteger0		NavParams._pinteger0		
#define pinteger1		NavParams._pinteger1	
#define pint0			NavParams._pint0		
#define pint1			NavParams._pint1	
#define ppointer		NavParams._ppointer	
#define plimit_value		NavParams._plimit_value		
#define plimit_cnt		NavParams._plimit_cnt	
#define pnavigate		NavParams._pnavigate		
#define pnavigate_cnt		NavParams._pnavigate_cnt	
#define pdepth_limiting		NavParams._pdepth_limiting	
#define pstring_bug		NavParams._pstring_bug		
#define padd_to_debug_string		NavParams._padd_to_debug_string
#define premove_from_debug_string 	NavParams._premove_from_debug_string
#define padd_to_debug			NavParams._padd_to_debug
#define premove_from_debug 	NavParams._premove_from_debug
#define padd_to_debug_cnt	NavParams._padd_to_debug_cnt
#define premove_from_debug_cnt 	NavParams._premove_from_debug_cnt

	/* LOCAL Function Prototypes */

LOCAL	bool	look_and_nswitch_on(const char*,N_SWITCH,N_CNT);
LOCAL	bool navigating(const char*);
LOCAL	bool	nswitch_on(const char*,N_SWITCH,N_CNT);
LOCAL   void  	white_space(const int);

EXPORT	void read_navigator_input(
	const char *fname)
{
	int i, j, k;
	char s[Gets_BUF_SIZE];
	char filename[256], line[MAX_N_CHARS + 1];
	char *c;
	int integer;
	float flote;
	POINTER pntr;

	if (navigator_loaded)
	{
	    (void) printf("\nnavigator data already loaded\n");
	    return;
        }

	num_of_navigate_names = 0;
	
	if (fname == NULL)
        {
	    /* Look for a Filename */
	    
	    navigate_input = stdin;
	    (void) strcpy(filename,"stdin");
	    (void) printf("\nEnter optional navigate input filename: ");
	    (void) fgets(s,Gets_BUF_SIZE-1,navigate_input);
	    (void) printf("%s\n",s);
	    if (sscanf(s,"%s",filename) == 1)
	    {  
		char scfilename[1024];
		stripcomm(scfilename,filename);
		if ((navigate_input = fopen(scfilename,"r")) == NULL)
		{
		    (void) printf("WARNING navigate file %s "
				  "could not be opened\n",filename);
		    navigate_input = stdin;
		    (void) strcpy(filename,"stdin");
		    clean_up(ERROR);
	        }
	    }
	    
	    (void) printf("List Functions to navigate on--end  ends list.\n");
	    (void) printf("Optional Keywords, one to each line following "
			  "function name and to\n");
	    (void) printf("which keyword is bound -> BRIEFLY, SUMMARY, "
			  "DETAILS, EXPLAIN, SUPPRESS,\n");
	    (void) printf("INT0, INT1, STRING_BUG, NAVIGATE, LIMIT_VALUE\n\n");
	}
	else 
        {
	    (void) strcpy(filename,fname);
	    
	    if ((navigate_input = fopen(filename,"r")) == NULL)
	    {
		(void) printf("ERROR navigate file %s could not be opened\n",
			      filename);
		clean_up(ERROR);
            }
        }

	/* read data */

	(void) fgets(line,MAX_N_CHARS + 1,navigate_input);
	(void) sscanf(line,"%s",s);
	
	for (i = 0 ; ; i++)
        {    
	    if (strcmp(s,"end")==0) 
            {
		printf("\t: end\n\n");
		break;
            }
	    
	    screen("\t: ");
	    screen("%s \n",s); 
	    if (navigating(s) == YES) 
            { 
		(void) printf("ERROR in read_navigator_input() --> "
			      "already navigating on %s\n",s);
		clean_up(ERROR);
            }
	    
	    (void) strncpy(navigate_names[num_of_navigate_names],s,MAX_N_CHARS);
	    pbriefly[num_of_navigate_names] = NO;
	    pbriefly_cnt[num_of_navigate_names] = -1;
	    psummary[num_of_navigate_names] = NO;
	    psummary_cnt[num_of_navigate_names] = -1;
	    pdetails[num_of_navigate_names] = NO;
	    pdetails_cnt[num_of_navigate_names] = -1;
	    pexplain[num_of_navigate_names] = NO;
	    pexplain_cnt[num_of_navigate_names] = -1;
	    pnavigate[num_of_navigate_names] = NO;
	    pnavigate_cnt[num_of_navigate_names] = -1;
	    pdepth_limiting[num_of_navigate_names] = NO;
	    plimit_value[num_of_navigate_names] = 1000000;
	    plimit_cnt[num_of_navigate_names] = -1;
	    pforce[num_of_navigate_names] = NO;
	    pforce_cnt[num_of_navigate_names] = -1;
	    padd_to_debug[num_of_navigate_names] = NO;
	    padd_to_debug_cnt[num_of_navigate_names] = -1;
	    premove_from_debug[num_of_navigate_names] = NO;
	    premove_from_debug_cnt[num_of_navigate_names] = -1;
	    psuppress[num_of_navigate_names] = NO;
	    psuppress_cnt[num_of_navigate_names] = -1;

	    pcall_count[num_of_navigate_names] = 0;
	    pint0[num_of_navigate_names] = -1;
	    pint1[num_of_navigate_names] = -1;
	    for (k = 0 ; k < 3 ; k++)
	    {
		pinteger0[num_of_navigate_names][k] = -1;
		pinteger1[num_of_navigate_names][k] = -1;
		pfloat0[num_of_navigate_names][k] = -1;
		pfloat1[num_of_navigate_names][k] = -1;
		pfloat2[num_of_navigate_names][k] = -1;
		ppointer[num_of_navigate_names][k] = NULL;
	    }
	    (void) strcpy(pstring_bug[num_of_navigate_names],"");
	    (void) strcpy(padd_to_debug_string[num_of_navigate_names],"");
	    (void) strcpy(premove_from_debug_string[num_of_navigate_names],"");

	    for (j = 0 ; ; j++)
	    {
		(void) fgets(line,MAX_N_CHARS + 1,navigate_input);
		(void) sscanf(line,"%s",s);

		if (strncasecmp(s,"sup",3)==0)
	        {
		    (void) printf("\t\t: SUPPRESS ");
		    psuppress[num_of_navigate_names] = YES;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			psuppress_cnt[num_of_navigate_names] = integer; 
			(void) printf("%d\n",integer);
	            }
	            else printf("\n");
	        }

		else if (strncasecmp(s,"bri",3)==0)
	        {
		    printf("\t\t: BRIEFLY ");
		    pbriefly[num_of_navigate_names] = YES;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pbriefly_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
			(void) printf("\n");
	        }

		else if (strncasecmp(s,"sum",3)==0)
	        {
		    printf("\t\t: SUMMARY ");
		    psummary[num_of_navigate_names] = YES;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			psummary_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
			(void) printf("\n");
	        }

		else if (strncasecmp(s,"det",3)==0)
	        {
		    printf("\t\t: DETAILS ");
		    pdetails[num_of_navigate_names] = YES;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pdetails_cnt[num_of_navigate_names] = integer;
			printf("%d\n",integer);
	            }
		    else printf("\n");
	        }

		else if (strncasecmp(s,"exp",3)==0)
	        {
		    (void) printf("\t\t: EXPLAIN ");
		    pexplain[num_of_navigate_names] = YES;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pexplain_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
			(void) printf("\n");
	        }
		
		else if (strncasecmp(s,"for",3)==0)
	        {
		    (void) printf("\t\t: FORCE ");
		    pforce[num_of_navigate_names] = YES;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
		        pforce_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
		    }
		    else (void) printf("\n");
		}

		else if (strncasecmp(s,"float0",6)==0)
	        {
		    (void) printf("\t\t: FLOAT0 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
		    {
			(void) sscan_float(c,&flote);
			(void) printf("%g ",flote);
			pfloat0[num_of_navigate_names][k++] = flote;
		    }
		    (void) printf("\n");
		}

		else if (strncasecmp(s,"float1",6)==0)
	        {
		    (void) printf("\t\t: FLOAT1 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			sscan_float(c,&flote);
			(void) printf("%g ",flote);
			pfloat1[num_of_navigate_names][k++] = flote;
	            }
		    (void) printf("\n");
	        }

		else if (strncasecmp(s,"float2",6)==0)
	        {
		    (void) printf("\t\t: FLOAT2 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			sscan_float(c,&flote);
			(void) printf("%g ",flote);
			pfloat2[num_of_navigate_names][k++] = flote;
	            }
		    (void) printf("\n");
	        }

		else if (strncasecmp(s,"int0",4)==0)
	        {
		    (void) printf("\t\t: INT0 ");
		    if (sscanf(line,"%*s %d",&integer) == 1)
		      pint0[num_of_navigate_names] = integer;
		    (void) printf("%d\n",pint0[num_of_navigate_names]);
	        }

		else if (strncasecmp(s,"int1",4)==0)
	        {
		    (void) printf("\t\t: INT1 ");
		    if (sscanf(line,"%*s %d",&integer) == 1)
		      pint1[num_of_navigate_names] = integer;
		    (void) printf("%d\n",pint1[num_of_navigate_names]);
	        }
		    
		else if (strncasecmp(s,"integer0",8)==0)
	        {
		    (void) printf("\t\t: INTEGER0 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			(void) sscanf(c,"%d",&integer);
			(void) printf("%d ",integer);
			pinteger0[num_of_navigate_names][k++] = integer;
	            }
		    (void) printf("\n");
	        }

		else if (strncasecmp(s,"integer1",8)==0)
	        {
		    (void) printf("\t\t: INTEGER1 ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			(void) sscanf(c,"%d",&integer);
			(void) printf("%d ",integer);
			pinteger1[num_of_navigate_names][k++] = integer;
	            }
		    (void) printf("\n");
	        }

		else if (strncasecmp(s,"poin",4)==0)
	        {
		    (void) printf("\t\t: POINTER ");
		    k = 0;
		    c = strtok(line," ");
		    while((c = strtok(NULL," ")) != NULL && k < 3)
	            {
			(void) sscanf(c,"%p",&pntr);
			(void) printf("%p ",pntr);
			ppointer[num_of_navigate_names][k++] = pntr;
	            }
		    (void) printf("\n");
	        }

		else if (strncasecmp(s,"limi",4)==0)
	        {
		    (void) sscanf(line,"%*s %d ",&integer);	
		    (void) printf("\t\t: LIMIT_VALUE %d ",integer);
		    pdepth_limiting[num_of_navigate_names] = YES;
		    plimit_value[num_of_navigate_names] = integer;
		    if (sscanf(line,"%*s %*d %d",&integer) == 1)
	            {
			plimit_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
			(void) printf("\n");
	        }
		
		else if (strncasecmp(s,"nav",3)==0)
		{
		    (void) printf("\t\t: NAVIGATE ");
		    pnavigate[num_of_navigate_names] = YES;
		    if (sscanf(line,"%*s %d",&integer) == 1)
	            {
			pnavigate_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
			(void) printf("\n");
	        }

		else if (strncasecmp(s,"str",3)==0)
	        {
		    (void) sscanf(line,"%*s %s",s);	
		    (void) printf("\t\t: STRING_BUG %s\n",s);
		    (void) strcpy(pstring_bug[num_of_navigate_names],s);
	        }

		else if (strncasecmp(s,"ADD",3)==0)
	        {
		    (void) sscanf(line,"%*s %s",s);	
		    (void) printf("\t\t: ADD_TO_DEBUG %s ",s);
		    padd_to_debug[num_of_navigate_names] = YES;
		    (void) strcpy(padd_to_debug_string[num_of_navigate_names],
				  s);
		    if (sscanf(line,"%*s %*s %d",&integer) == 1)
	            {
	        	padd_to_debug_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
			(void) printf("\n");
	        }

		else if (strncasecmp(s,"rem",3)==0)
	        {
		    (void) sscanf(line,"%*s %s",s);	
		    (void) printf("\t\t: REMOVE_FROM_DEBUG %s ",s);
		    premove_from_debug[num_of_navigate_names] = YES;
		    (void) strcpy(
			premove_from_debug_string[num_of_navigate_names],s);
		    if (sscanf(line,"%*s %*s %d",&integer) == 1)
	            {
			premove_from_debug_cnt[num_of_navigate_names] = integer;
			(void) printf("%d\n",integer);
	            }
		    else
			(void) printf("\n");
	        }

		else 
		   break;

	    }
	    
	    if (((int) strlen(navigate_names[num_of_navigate_names])) >=
							   MAX_N_CHARS)
		navigate_names[num_of_navigate_names][MAX_N_CHARS]='\0';

	    num_of_navigate_names++;
	    (void) printf("\n");
	}
	
	if (navigate_input != stdin) fclose(navigate_input);

	fflush(stdout);
		
	navigator_loaded = YES;

	return;
}		/*end read_navigator_input*/


EXPORT void print_navigator_data(
	const char *fname)
{
  	int i, j;

	(void) printf("\n");
	(void) printf("Printout of navigate data structure from %s()\n",fname);

	for (i = 0 ; i < num_of_navigate_names ; i++)
	{
	    for (j = 0 ; j < 110 ; j++)
		(void) printf("=");
	    (void) printf("\n\n");
	    (void) printf("\n%s()\n",navigate_names[i]);
	    (void) printf("\n  %-8s %1s | %-8s %1s | %-8s %1s | "
			  "%-8s %1s | %-8s %1s | %-8s %1s | %-8s "
			  "%1s | %-8s %1d\n",
		          "NAVIGATE",YN(pnavigate[i]),"SUPPRESS",
			  YN(psuppress[i]),"BRIEFLY",YN(pbriefly[i]),
			  "SUMMARY",YN(psummary[i]),"DETAILS",
			  YN(pdetails[i]),"EXPLAIN",YN(pexplain[i]),
		          "FORCE",YN(pforce[i]),"LIMIT_V",
			  plimit_value[i]);
	    printf("  ");
	    for (j = 0 ; j < 7 ; j++)
		(void) printf("-----------|-");
	    (void) printf("-----------\n");
	    (void) printf("  %-10d | %-10d | %-10d | %-10d | "
			  "%-10d | %-10d | %-10d | %-10d\n\n",
		          pnavigate_cnt[i],psuppress_cnt[i],
			  pbriefly_cnt[i],psummary_cnt[i],pdetails_cnt[i],
		          pexplain_cnt[i],pforce_cnt[i],plimit_cnt[i]); 
	    (void) printf("\n  index = %d  depth_limiting = %1s  "
			  "C_COUNT = %d  STRING_BUG = %-s\n\n",
		          i,YN(pdepth_limiting[i]),pcall_count[i],
			  pstring_bug[i]);
	    (void) printf("  INTEGER0 (%10d %10d %10d)  INT0 = %d\n",
			  pinteger0[i][0],pinteger0[i][1],
			  pinteger0[i][2],pint0[i]); 
	    (void) printf("  INTEGER1 (%10d %10d %10d)  INT1 = %d\n",
			  pinteger1[i][0],pinteger1[i][1],
			  pinteger1[i][2],pint1[i]); 
	    (void) printf("  POINTER  (%10p %10p %10p)\n",
			  ppointer[i][0],ppointer[i][1],ppointer[i][2]);
	    (void) printf("  FLOAT0   (%10g %10g %10g)\n",
			  pfloat0[i][0],pfloat0[i][1],pfloat0[i][2]);
	    (void) printf("  FLOAT1   (%10g %10g %10g)\n",
			  pfloat1[i][0],pfloat1[i][1],pfloat1[i][2]);
	    (void) printf("  FLOAT2   (%10g %10g %10g)\n\n",
			  pfloat2[i][0],pfloat2[i][1],pfloat2[i][2]);
	    (void) printf("  %1s ADD_TO_DEBUG      %10d  %-s\n",
		          YN(padd_to_debug[i]),
			  padd_to_debug_cnt[i],
			  padd_to_debug_string[i]);
	    (void) printf("  %1s REMOVE_FROM_DEBUG %10d  %-s\n",
		          YN(premove_from_debug[i]),
			  premove_from_debug_cnt[i],
			  premove_from_debug_string[i]);
	    
	}
	(void) printf("\n");
	for (j = 0 ; j < 110 ; j++)
	    (void) printf("=");
	(void) printf("\n\n");
	
	(void) printf("\n");
	(void) navigator_update(pPRINT_STACK,fname,0,0); 
	(void) printf("End of Printout of navigate data "
		      "structure from %s()\n",fname);
	(void) printf("\n");
}                /* end print_navigator_data() */

LOCAL bool navigator_switches(
	const N_SWITCH command,
	const int      index) 
{
  	if (index == -1)
	    return NO;
	switch(command)
	{
	case pBRIEFLY :
	    return pbriefly[index];
	case pSUMMARY :
	    return psummary[index];
	case pDETAILS :
	    return pdetails[index ];
	case pEXPLAIN :
	    return pexplain[index];
	case pNAVIGATE :
	    return pnavigate[index];
	case pDEPTH_LIMITING :
	    return pdepth_limiting[index];
	case pFORCE:
	    return pforce[index];
	case pADD_TO_DEBUG:
	    return padd_to_debug[index];
	case pREMOVE_FROM_DEBUG:
	    return premove_from_debug[index];
	case pSUPPRESS:
	    return psuppress[index];
	default : 
	    screen("ERROR in  navigator_switches(), "
		    "unknown command = %d\n",command);
	    clean_up(ERROR);
	} 
	return NO;
} 	/* end navigator_switches */

LOCAL int  navigator_integers(
	const N_CNT  command,
	const int    index)
{
  	if (index < 0)
	    return -1;
	switch(command)
	{
	case pBRIEFLY_CNT:
	    return pbriefly_cnt[index];
	case pSUMMARY_CNT:
	    return psummary_cnt[index];
	case pDETAILS_CNT:
	    return pdetails_cnt[index];
	case pEXPLAIN_CNT:
	    return pexplain_cnt[index];
	case pFORCE_CNT:
	    return pforce_cnt[index];
	case pNAVIGATE_CNT:
	    return pnavigate_cnt[index];
	case pADD_TO_DEBUG_CNT:
	    return padd_to_debug_cnt[index];
	case pREMOVE_FROM_DEBUG_CNT:
	    return premove_from_debug_cnt[index];
	case pC_COUNT:
	    return pcall_count[index];
	case pLIMIT_VALUE :
	    return plimit_value[index ];
	case pLIMIT_CNT :
	    return plimit_cnt[index ];
	case pSUPPRESS_CNT:
	    return psuppress_cnt[index];
	case pINT0:
	    return pint0[index];
	case pINT1:
	    return pint1[index];
	}
	return -1;
}	/* end navigator_integers */

LOCAL int *navigator_integer_pointers(
	const N_PINTEGER command,
	const int        index)
{
  	if (index < 0)
	    return NULL;
  	switch(command)
        {
	case pINTEGER0:
	    return pinteger0[index];
	case pINTEGER1:
	    return pinteger1[index];
	default:
	    screen("ERROR in navigator_integer_pointers(), "
		   "no such command %d\n",command);
	    clean_up(ERROR);
	    break;
	}
	return NULL;
}	/* end navigator_integer_pointers */

LOCAL float *navigator_float_pointers(
	const N_PFLOAT command,
	const int      index)
{
  	if (index < 0)
	    return NULL;
  	switch(command)
	{
	case pFLOAT0:
	    return pfloat0[index];
	case pFLOAT1:
	    return pfloat1[index];
	case pFLOAT2:
	    return pfloat2[index];
	default:
	    screen("ERROR in navigator_float_pointers(), "
		   "no such command %d\n",command);
	    clean_up(ERROR);
	    break;
	}
	return NULL;
}	/* end navigator_float_pointers */

LOCAL POINTER *navigator_pointers(
	const N_PPOINTER command,
	const int index)
{
  	if (index < 0)
	    return NULL;
	switch(command)
	{
	case pPOINTER:
	    return ppointer[index];
	default:
	    screen("ERROR in navigator_pointers(), "
		   "no such command %d\n",command);
	    clean_up(ERROR);
	    break;
	}
	return NULL;
}	/* end navigator_pointers */

LOCAL char *navigator_strings(
	const N_STRING command,
	const int index)
{
  	if (index < 0)
	    return NULL;
  	switch(command)
	{
	case pADD_TO_DEBUG_STRING:
	    return padd_to_debug_string[index];
	case pREMOVE_FROM_DEBUG_STRING:
	    return premove_from_debug_string[index];
	case pSTRING_BUG:
	    return pstring_bug[index];
	default:
	    screen("ERROR in navigator_strings(), "
		   "no such command %d\n",command);
	    clean_up(ERROR);
	    break;
	}
	return NULL;
}	/* end navigator_strings */

LOCAL  int  navigator_update(
	const N_UPDATE  command,
	const char *fname,
	const int  limit_integer,
	const int  fname_index)
{
  	static int depth = 0;
	static int limit = 1000000;
	static int function_suppression = NO;
	static int reenable_depth = -1;
	static int reenable_index = -2;
	int old_limit, i;
	
	static char stack[MAX_N_NAMES][MAX_N_CHARS + 1];

	switch(command) 
	{
	case pINIT:
	    depth = limit_integer;
	    return 0;

	case pLIMIT:
	    old_limit = limit - depth;
	    limit = limit_integer + depth;
	    nprintf("LIMIT_TRACE in %s(%d),set LIMIT = %d  "
		    "from lim val = %d,depth = %d; old lim = %d\n",
		    fname,pcall_count[fname_index],limit,
		    limit_integer,depth,old_limit);
	    return old_limit;
      
	case pENTER:
	    if (fname_index >= 0)
		pcall_count[fname_index]++;

	    if (NAVIGATING() == NO)
		return 0;

	    depth++;

	    if (depth >= 0)
		(void) strncpy(stack[depth],fname,MAX_N_CHARS);

	    if ((fname_index >= 0) && (SUPPRESS(fname)==YES) &&
		(!function_suppression))
	    {
	        (void) function_monitor(pSTOP);
		reenable_depth = depth;
		reenable_index = fname_index;
		function_suppression = YES;
	    }

	    if (depth == limit + 1)
		(void) depth_monitor(pSTOP);
	
	    if (LOOK() == YES)
	    {
		white_space(depth % 20);
		(void) printf("\\%d %s(%d) \n",depth,fname,
						 pcall_count[fname_index]);
	    }
	    return 0;

	case pLEAVE:
	    if (NAVIGATING() == NO)
		return 0;

	    if (depth >= 0)
	    {
		if (strcmp(stack[depth],fname) != 0) 
		    (void) printf("warning: stack misaligned,"
				  "function_name = %s()  "
				  "look in stack[%d] = %s()\n",
			          fname,depth,stack[depth]);
	    }

	    if (LOOK() == YES)
	    {
	        white_space(depth % 20);
		(void) printf("/%d %s(%d) \n",depth,fname,
			      pcall_count[fname_index]);
	    }

	    if (depth == limit + 1)
		(void) depth_monitor(pSTART);
	    
	    if (function_suppression &&
		(reenable_depth==depth) &&
		(reenable_index==fname_index)) 
	    {
	        (void) function_monitor(pSTART);
		reenable_depth = -1;
		reenable_index = -2;
		function_suppression = NO;
	    }
	    depth--;
	    return 0;
      
	case pSTATIC:
	    white_space(depth % 20);
	    (void) printf(" %d>\n",depth);
	    return 0;
      
	case pSPACE:
	    white_space(depth % 20);
	    (void) printf(" %d> ",depth);
	    return 0;      
	    
	case pPRINT_STACK: 
	    if (depth < 1) 
	    { 
		(void) printf("printout of stack from %s(),depth = %d, "
			      "stack is empty\n",fname,depth);
		return 0; 
	    }
	    (void) printf("printout of stack from %s(),depth = %d  "
			 "limit = %d\n\n",fname,depth,limit);
	    for (i = 1 ; i <= depth ; i++)
		(void) printf("\tstack[%3d] = %s()\n",i,stack[i]);
	    (void) printf("\n\tend of stack\n");
	    return 0;
	}

	(void) printf("Warning in navigator_update(),"
		      "command = %d unknown\n",command);
	return 0;
} 	/* end navigator_update */

LOCAL int navigator_index(
	const char *fname)
{
	static char last_fname[1024];
  	static int i = -2;

	if (navigator_loaded == NO)
	    return -2;

	if (strcmp(fname,last_fname) == 0)
	    return i;

	(void) strcpy(last_fname,fname);
	for (i = 0 ; i < num_of_navigate_names ; i++) 
	{
	    if (strncmp(fname,navigate_names[i],MAX_N_CHARS)==0)
		return i;
	}
	i = -1;
	return i;
}	/* end navigator_index */

EXPORT	void	nprint_long_string(
	char	*s )
{
	char	*c;
	int	len;
	int	maxlen = 79;
	int	newlen;
	int	tablen = 0;
	bool start_of_line;
	static	char	*line = NULL;
	static	size_t	allocated_length = 0;

	if ((NAVIGATING() == NO) || (strlen(s) == 0))
	    return;

	if (strlen(s) >= allocated_length)
	{
	    if (line != NULL)
	    	free(line);
	    allocated_length = strlen(s) + 1024;
	    line = (char*) malloc(allocated_length*sizeof(char));
	}
	strcpy(line,s);

	for (len = 0, c = strtok(line," "); c != NULL; c = strtok(NULL," "))
        {
	    newlen = len+(int)strlen(c);
	    if ((newlen+1) >  maxlen)
            {
		screen("\n");
		start_of_line = YES;
		len = tablen;
		newlen = len+(int)strlen(c);
            }
	    else if (start_of_line == NO)
            {
		screen(" ");
		newlen++;
		len++;
            }

	    if( start_of_line )
	        (void) navigator_update(pSPACE,"",0,0);

	    screen("%s",c);
	    start_of_line = NO;

	    if( c[ strlen( c ) - 1 ] == '\n' ) 
            {
		len = tablen;
		start_of_line = YES;
            }

	    else if( newlen == maxlen )
            {
		screen("\n");
		start_of_line = YES;
		len = tablen;
            }
	    else
            {
		if (c[strlen(c)-1] == '.')
		{
		    screen(" ");
		    len++;
	        }
		len += (int)strlen(c);
	    }
	}
}		/*end nprint_long_string*/

LOCAL	bool	ADD_TO_DEBUG(
	const char *fname)
{
	return nswitch_on(fname,pADD_TO_DEBUG,pADD_TO_DEBUG_CNT);
}		/*end ADD_TO_DEBUG*/

LOCAL	bool	DEPTH_LIMITING(
	const char *fname)
{
	return nswitch_on(fname,pDEPTH_LIMITING,pLIMIT_CNT);
}		/*end DEPTH_LIMITING*/

EXPORT	bool	NAVIGATOR_DETAILS(
	const char *fname)
{
	return look_and_nswitch_on(fname,pDETAILS,pDETAILS_CNT);
}		/*end NAVIGATOR_DETAILS*/

EXPORT	bool	NAVIGATOR_EXPLAIN(
	const char *fname)
{
	return look_and_nswitch_on(fname,pEXPLAIN,pEXPLAIN_CNT);
}		/*end NAVIGATOR_EXPLAIN*/

LOCAL	bool	NAVIGATING(void)
{
	return (navigator_trace(pLOOK,"") != pSTOP) ? YES : NO;
}		/*end NAVIGATING*/

LOCAL	bool	REMOVE_FROM_DEBUG(
	const char *fname)
{
	return nswitch_on(fname,pREMOVE_FROM_DEBUG,pREMOVE_FROM_DEBUG_CNT);
}		/*end REMOVE_FROM_DEBUG*/

EXPORT	bool	NAVIGATOR_BRIEFLY(
	const char *fname)
{
	return ((look_and_nswitch_on(fname,pBRIEFLY,pBRIEFLY_CNT)==YES) ||
		(nswitch_on(fname,pFORCE,pFORCE_CNT)==YES)) ? YES : NO;
}		/*end NAVIGATOR_BRIEFLY*/

EXPORT	bool	NAVIGATOR_SUMMARY(
	const char *fname)
{
	return ((look_and_nswitch_on(fname,pSUMMARY,pSUMMARY_CNT)==YES) ||
		(nswitch_on(fname,pFORCE,pFORCE_CNT)==YES)) ? YES : NO;
}		/*end NAVIGATOR_SUMMARY*/

LOCAL	bool	SUPPRESS(
	const char *fname)
{
	return nswitch_on(fname,pSUPPRESS,pSUPPRESS_CNT);
}		/*end SUPPRESS*/

LOCAL	char	*ADD_TO_DEBUG_STRING(
	const char *fname)
{
	int fname_index;
	fname_index = navigator_index(fname);

	return navigator_strings(pADD_TO_DEBUG_STRING,fname_index);
}		/*end ADD_TO_DEBUG_STRING*/

LOCAL	char	*REMOVE_FROM_DEBUG_STRING(
	const char *fname)
{
	int fname_index;
	fname_index = navigator_index(fname);

	return navigator_strings(pREMOVE_FROM_DEBUG_STRING,fname_index);
}		/*end REMOVE_FROM_DEBUG_STRING*/

LOCAL	void	CHECK_ADD_TO_DEBUG(
	const char *fname)
{
	if (ADD_TO_DEBUG(fname) != YES)
	    return;
	add_to_debug(ADD_TO_DEBUG_STRING(fname));
}		/*end CHECK_ADD_TO_DEBUG*/

LOCAL	void	CHECK_REMOVE_FROM_DEBUG(
	const char *fname)
{
	if (REMOVE_FROM_DEBUG(fname) != YES)
	    return;
	remove_from_debug(REMOVE_FROM_DEBUG_STRING(fname));
}		/*end CHECK_REMOVE_FROM_DEBUG*/

EXPORT	void	NAVIGATOR_DEBUG_ENTER(
	const char          *fname)
{
	NAVIGATE_ENTER_LEAVE(fname,ENTER_NAVIGATOR);
}		/*end NAVIGATOR_DEBUG_ENTER*/

EXPORT	void	NAVIGATOR_DEBUG_LEAVE(
	const char          *fname)
{
	NAVIGATE_ENTER_LEAVE(fname,LEAVE_NAVIGATOR);
}		/*end NAVIGATOR_DEBUG_LEAVE*/

EXPORT	int	NAVIGATOR_INT0(
	const char	*fname)
{
	return navigator_integers(pINT0,navigator_index(fname));
}		/*end NAVIGATOR_INT0*/

EXPORT	int	NAVIGATOR_INT1(
	const char	*fname)
{
	return navigator_integers(pINT1,navigator_index(fname));
}		/*end NAVIGATOR_INT1*/

EXPORT int *NAVIGATOR_INTEGER0(
	const char	*fname)
{
	return navigator_integer_pointers(pINTEGER0,navigator_index(fname));
}		/*end NAVIGATOR_INTEGER0*/

EXPORT int *NAVIGATOR_INTEGER1(
	const char	*fname)
{
	return navigator_integer_pointers(pINTEGER1,navigator_index(fname));
}		/*end NAVIGATOR_INTEGER1*/

EXPORT float *NAVIGATOR_FLOAT0(
	const char	*fname)
{
	return navigator_float_pointers(pFLOAT0,navigator_index(fname));
}		/*end NAVIGATOR_FLOAT0*/

EXPORT float *NAVIGATOR_FLOAT1(
	const char	*fname)
{
	return navigator_float_pointers(pFLOAT1,navigator_index(fname));
}		/*end NAVIGATOR_FLOAT1*/

EXPORT float *NAVIGATOR_FLOAT2(
	const char	*fname)
{
	return navigator_float_pointers(pFLOAT2,navigator_index(fname));
}		/*end NAVIGATOR_FLOAT2*/

EXPORT POINTER *NAVIGATOR_POINTR(
	const char	*fname)
{
	return navigator_pointers(pPOINTER,navigator_index(fname));
}		/*end NAVIGATOR_POINTR*/

EXPORT	char	*NAVIGATOR_STRING_BUG(
	const char	*fname)
{
	return navigator_strings(pSTRING_BUG,navigator_index(fname));
}		/*end NAVIGATOR_STRING_BUG*/

LOCAL	void	NAVIGATE_ENTER_LEAVE(
	const char          *fname,
	const N_ENTER_LEAVE el)
{
	int fname_index;
	static int reset_limit;

	if (nswitch_on(fname,pNAVIGATE,pNAVIGATE_CNT) != YES)
	    return;

	fname_index = navigator_index(fname);
	if (el == ENTER_NAVIGATOR)
	{
	    (void) navigator_update(pENTER,fname,0,fname_index);
	    (void) navigator_trace(pSTART,fname);
	    CHECK_ADD_TO_DEBUG(fname);
	    if ((LOOK() == YES) && (DEPTH_LIMITING(fname) == YES))
	    {
		int limit_value = navigator_integers(pLIMIT_VALUE,fname_index);
		reset_limit = navigator_update(pLIMIT,fname,limit_value,
					       fname_index);
	    }
	}
	else if (el == LEAVE_NAVIGATOR)
	{
	    if ((LOOK() == YES) && (DEPTH_LIMITING(fname) == YES))
		(void) navigator_update(pLIMIT,fname,reset_limit,fname_index);
	    CHECK_REMOVE_FROM_DEBUG(fname);
	    (void) navigator_trace(pSTOP,fname);
	    (void) navigator_update(pLEAVE,fname,0,fname_index);
	}
}		/*end NAVIGATE_ENTER_LEAVE*/

LOCAL  N_MONITOR depth_monitor(
	const N_MONITOR command)
{
	static N_MONITOR monitor = pSTOP;
	switch(command)
	{
	case pSTART:
	case pSTOP:
	    monitor = command;
	    break;
	case pLOOK:
	    break;
	}
	return monitor;
}	/* end depth_monitor */

LOCAL  N_MONITOR function_monitor(
	const N_MONITOR command)
{
  	static N_MONITOR monitor = pSTOP;
	switch(command)
	{
	case pSTART:
	case pSTOP:
	    monitor = command;
	    break;
	case pLOOK:
	    break;
	}
	return monitor;
}	/* end function_monitor */

LOCAL	bool LOOK(void)
{
	return ((depth_monitor(pLOOK) != pSTOP) &&
		(function_monitor(pLOOK) != pSTOP)) ? YES : NO;
}		/*end LOOK*/

LOCAL	bool	look_and_nswitch_on(
	const char *fname,
	N_SWITCH   n_switch,
	N_CNT      n_cnt)
{
	if (LOOK() != YES)
	    return NO;
	return nswitch_on(fname,n_switch,n_cnt);
}		/*end nswitch_on*/

LOCAL	N_MONITOR navigator_trace(
	const N_MONITOR command,
	const char              *function_name)
{
  	static N_MONITOR trace = pSTOP;
	static char start_function[MAX_N_CHARS + 1] = "";
	static char stop_function[MAX_N_CHARS + 1] = "";
	
	switch(command)
	{
	case pSTOP:
	case pSTART:
	    if (trace == command)
	    {
	        (void) printf("Warning: navigator_trace() already = "
			      "command %d\n"
		              "\tstart_function = %s()\n"
		              "\tstop_function  = %s()\n",
			      command,start_function,stop_function);
	    }
	    if (command == pSTART) 
	    {
		if (NAVIGATING() == YES)
		    navigator_update(pSPACE,"",0,0);
		(void) printf("NAVIGATOR ON in %s\n",function_name);
		(void) depth_monitor(pSTART);
		(void) function_monitor(pSTART);
		(void) strcpy(start_function,function_name);
	    }
	    if (command == pSTOP) 
	    {
	        nprintf("NAVIGATOR OFF in %s\n",function_name);
		(void) depth_monitor(pSTOP);
		(void) function_monitor(pSTOP);
		(void) strcpy(stop_function,function_name);
	    }
	    trace = command;
	    return pSTOP;
     
	case pLOOK:
	    return trace;
      
	default: 
	    screen("ERROR in navigator_trace(), Unknown command = %d\n",
		   command);
	    clean_up(ERROR);
	}
	return trace;
}	/* end navigator_trace */


LOCAL	bool	nswitch_on(
	const char *fname,
	N_SWITCH   n_switch,
	N_CNT      n_cnt)
{
	int fname_index;
	int c_count;
	int cnt;
	bool nav_switches;

	fname_index = navigator_index(fname);
	nav_switches = navigator_switches(n_switch,fname_index);
	if (nav_switches != YES)
	    return NO;
	cnt = navigator_integers(n_cnt,fname_index);
	c_count = navigator_integers(pC_COUNT,fname_index);
	return ((cnt == -1) || (cnt == c_count)) ? YES : NO;
}		/*end nswitch_on*/

LOCAL	const char	*YN(
	bool	y_or_n)
{
	switch (y_or_n)
	{
	case YES:
	    return "Y";
	case NO:
	    return "N";
	default:
	    screen("ERROR in YN(), no such value %d\n",y_or_n);
	    clean_up(ERROR);
	}
	return NULL;
}		/*end YN*/


LOCAL void white_space(
	const int num_spaces)
{
	int i;
  	for (i = 0 ; i < num_spaces ; i++)
	    (void) printf(" ");
}    /* end white_space() */

LOCAL	bool navigating(
	const char *funcname)
{
	int i;

	for (i = 0 ; i < num_of_navigate_names ; i++) 
	    if (strncmp(funcname,navigate_names[i],MAX_N_CHARS)==0)
	        return YES;
	return NO;
}		/*end navigating*/

/* end of Navigator functions */

#include <stdarg.h>

/* VARARGS */
EXPORT	void	nprintf(
	const char	*fmt,
	...)
{
	va_list ap;

	if (NAVIGATING() == NO)
	    return;

	(void) navigator_update(pSPACE,"",0,0);

	va_start(ap, fmt);
	(void) vprintf(fmt,ap);
	va_end(ap);
}	/*end nprintf*/
#endif /*DONT_COMPILE*/ /*SUGGESTED REVISION TO NAVIGATOR 19990423*/

#endif /* defined(NAVIGATOR) */
