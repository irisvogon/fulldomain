/*
*
*                               gimkbub.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       This file contains functions for the construction of non-trivial
*       geometric curves needed for the initialization of an interface.
*
*/

#if defined(TWOD) || defined(THREED)

#include <ginit/ginit.h>

//#bjet2
enum _INSERT_FLAG {
	INSERT_OFF = -1,   //other comp or near a surf
        PHYSICAL_OFF = 0,  //physical no bubbles
        INSERT_ON = 1,     //can insert bubble
        NEW_INSERTED = 2,  //prepare to insert 
};
typedef enum _INSERT_FLAG INSERT_FLAG;

typedef struct {
        int     size;          // bubble diameter in mesh size
	int     spacing_size;  // dist between bubbles.
	float   r;
	float   Rbdry;
	int     compin;
	int     compout;
}  BUBBLES;

typedef struct {
        int    icoords[3];  // icrds of checked cells 
        float  press;       // pressure in the cell 
} INSERT_SITE;

LOCAL   bool    tris_near_cell(float*, RECT_GRID*, INTERFACE*);
LOCAL   bool    g_bubble_insert_cond(int*, Wave*);
LOCAL   bool    test_check_insert_bubble_position(float*,int*,int*,int*,int***,
			int*,int*,Wave*,BUBBLES*);
LOCAL   bool    check_insert_bubble_position(float*,int*,int*,int*,int***,
			int*,int*,Wave*,BUBBLES*);
LOCAL   void    average_states_inside_bubble(Locstate,int*,int*,int***,int*,int*,Wave*);
LOCAL 	void    show_pos_flag(char*,int***,int*,int*,RECT_GRID*);
LOCAL   void    set_state_on_phase_bdry(Locstate, Locstate);
LOCAL   void    get_gas_ref_state(Locstate,int,INTERFACE*);
LOCAL   void    set_bubble_params(BUBBLES*, RECT_GRID*,int,int);
LOCAL   void    set_bubble_diameter(BUBBLES*,RECT_GRID*,int);
LOCAL   int     press_ascend(const void*,const void*);

LOCAL   SURFACE*  create_bubbles3d(SURFACE*,Locstate,Locstate,Wave*,Front*,int***,
			float*,int*,int*,int*,int*,BUBBLES*);
LOCAL   void    reset_front_and_states(Wave*,Front*,int);
LOCAL   int	make_bubbles(Wave*, Front*);
LOCAL   int	make_bubbles3d(Wave*, Front*);

//outside func
bool  is_tri_outside_box(TRI *, float **);
LOCAL  int make_one_bubbles3d(Wave*,Front*);

EXPORT  void    g_prompt_for_dynamic_bubble_insertion3d(
        CHART  *chart)
{
        char      s[121];

        chart->make_bubbles = NULL;

        screen("Type 'y' to dynamically insert vapor bubbles in the"
               " liquid (y, n(dflt)): ");

        (void) Gets(s);
        if (strncmp(s,"y",1) == 0)
        {
            chart->make_bubbles = make_bubbles;
        }
}

//IMPORT  int   make_bubbles(Wave*,Front*);
//WARNING:  see  init_fuel_injection_jet3d  for predifine
LOCAL   int make_bubbles(
        Wave                 *wave,
        Front                *front)
{
        int      num_bubbles;

	if(wave == NULL || front == NULL)
	    return YES;

	//return YES;

	switch(wave->rect_grid->dim)
        {
            case 2:
	        return YES;
            case 3:
                num_bubbles = make_bubbles3d(wave,front);
                //num_bubbles = make_one_bubbles3d(wave,front);
                break;
            default:
                printf("ERROR: make_bubbles(), implement for dim %d\n",
                       wave->rect_grid->dim);
                clean_up(ERROR);
        }
       
	reset_front_and_states(wave,front,num_bubbles);

        pp_global_isum(&num_bubbles,1);
       
	//check_print_front("After bubble insertion", "bub_ins", 's', 
	//          front->interf, front->step, front->step, NO);
 
	return num_bubbles;
}

LOCAL bool  tris_near_cell(
	float	    *coords,
	RECT_GRID   *gr,
	INTERFACE   *intfc)
{
RECT_GRID      *top_grid = &topological_grid(intfc);
struct Table   *T = table_of_interface(intfc);
float          *h = gr->h;
static float   **fbox = NULL;
int	       ip0[3], ip1[3];
int	       i, j, k, n;

	if(fbox == NULL)
	    bi_array(&fbox, 2, 3, FLOAT);

	//The cell with coords as the center.
	for(i=0; i<3; i++)
	{
	    fbox[0][i] = coords[i] - 0.5*h[i];
	    fbox[1][i] = coords[i] + 0.5*h[i];
	}
	
	//near bdry 
	if((!rect_in_which(fbox[0], ip0, top_grid)) || 
	   (!rect_in_which(fbox[1], ip1, top_grid)) )
	    return YES;
	
	for(i = ip0[0]; i <= ip1[0]; i++)
	    for(j = ip0[1]; j <= ip1[1]; j++)
	        for(k = ip0[2]; k <= ip1[2]; k++)
		{
		    if(T->compon3d[k][j][i] != ONFRONT)
		        continue;
		    for(n = 0; n < T->num_of_tris[k][j][i]; n++)
		        if(!is_tri_outside_box(T->tris[k][j][i][n], fbox))
			    return YES;
		}

	return NO;
}

//physical condition for bubble insertion
LOCAL  bool  g_bubble_insert_cond(
	int    *icoords,
	Wave   *wave)
{
Locstate    state;
float	    P_sat = -10.0;

	state = Rect_state(icoords,wave);
	if(pressure(state) < P_sat)
	{
	    //printf("#press %24.16e\n", pressure(state));
	    return  YES;
	}
	return NO;

}

enum {
G_HMIN_ID     = FIRST_PHYSICS_MESSAGE_ID+10
};


LOCAL void  set_one_bubble_states(Locstate, Locstate);
void    get_default_fluid_comp(int*, int*, INTERFACE*);

LOCAL  int make_one_bubbles3d(
        Wave                 *wave,
        Front                *front)
{
	int	     compin, compout;
	BUBBLES      bubbles;
	INTERFACE    *intfc = front->interf;
	RECT_GRID    *gr = wave->rect_grid;
	SURFACE      *surf, **s, *tmp_surf;
	size_t       sizest = wave->sizest;
	static Locstate     gst=NULL, ref_st=NULL;
	float	     center[3] = {0.0, 0.0, 0.0}, ang;
	int	     i, j;

	if(gst == NULL)
	{
	    scalar(&gst, sizest);
	    scalar(&ref_st, sizest);
	}

	set_wall_flag_for_surface(intfc);
	
	surf = NULL;
	for(s=intfc->surfaces; s && *s; s++)
	    if(!is_wall_surface(*s))
	    {
	        //ASSUME one fluid surface and positive_component is the gas
	        surf = *s;
		compin = positive_component(*s);
	        compout = negative_component(*s);
		break;
	    }

	//WARNING: can not creat bubble with no fluid-gas surf, this should be fixed.
	if(surf == NULL)
	{
	    //WARNING: ASSUME first param is the fluid comp
	    get_default_fluid_comp(&compout,&compin,intfc);
	    printf("#surf==NULL when finding physical surf, "
		   "compin=%d compout=%d\n", compin, compout);
	    //debug_print("bubbles_mpi","Left make_bubbles3d()\n");
	    //printf("#left make_bubbles3d surf==NULL\n\n");
	    //clean_up(0);
	    //return 0;
	}

	set_bubble_params(&bubbles, gr, compin, compout);
	
	//set states
	get_gas_ref_state(gst, compin, intfc);
	get_gas_ref_state(ref_st, compout, intfc);
	set_one_bubble_states(gst, ref_st);

	//19
	for(i=0; i<1; i++)
	{
	    //center[2] = -0.9 + 0.1*i;
	    center[2] = 0.0;
	    center[0] = 0.0;
	    center[1] = 0.0;

	    tmp_surf = create_bubbles3d(surf,gst,ref_st,wave,front,NULL,
	        center,NULL,NULL,NULL,NULL,&bubbles);
	    if(tmp_surf != NULL)
		surf = tmp_surf;
	}

	return 1;
	
	ang = 2.0*PI/8.0;
	for(j=8; j<9; j++)
	{
	center[0] = 0.36*cos(ang*j);
	center[1] = 0.36*sin(ang*j);
	if(j == 8)
	{
	    center[0] = 0.0;
	    center[1] = 0.0;
	}
	for(i=0; i<15; i++)
	{
	    //center[0] = 0.0;
	    //center[1] = 0.0;
	    center[2] = (i-7.0)*0.23;

	    tmp_surf = create_bubbles3d(surf,gst,ref_st,wave,front,NULL,
	        center,NULL,NULL,NULL,NULL,&bubbles);
	    if(tmp_surf != NULL)
		surf = tmp_surf;
	}
	}
	
	return 1;
}

LOCAL  int make_bubbles3d(
        Wave                 *wave,
        Front                *front)
{
INTERFACE    *intfc = front->interf;
RECT_GRID    *gr = wave->rect_grid;
SURFACE      *surf, **s;
int	     ***pos_flag;
int          *gmax = gr->gmax;
int          i,j,k,n, ix,iy,iz;
int	     imin[3], imax[3], smax[3], icoords[3], inmin[3], inmax[3];
int	     compin, compout, nsite, alloc_site, num_bubble;
size_t       sizest = wave->sizest;
float	     *coords, center[3];
Locstate     state;
static Locstate     gst=NULL, ref_st=NULL;
static int   nbub = 1;   //different diameters of bubbles
static int   dbub[3] = {3, 3, 3};
BUBBLES      bubbles;
INSERT_SITE  *insert_site;
float	     spress, spmin;
int	     sendnode;
int	     picrds[3];

        debug_print("bubbles_mpi","Enter make_bubbles3d()\n");

	if(gst == NULL)
	{
	    alloc_state(intfc, &gst, sizest);
	    alloc_state(intfc, &ref_st, sizest);
	}
	//in case the intfc is modified, prepare for tris_near_cell
	if (intfc->modified || intfc->table->new_grid)
	    make_interface_topology_lists(intfc);

	printf("\n\n#make_bubbles gr=\n");
	//print_rectangular_grid(gr);
	//print_rectangular_grid(computational_grid(intfc));
	//print_rectangular_grid(&topological_grid(intfc));

        for (i = 0; i < 3; i++)
	{
            imin[i] = -gr->lbuf[i];        
	    imax[i] = gmax[i] + gr->ubuf[i];
	    
	    //avoid the subdomain bdry, otherwise, scatter_front will fail
	    if(buffered_boundary_type(rect_boundary_type(intfc,i,0)) == OPEN_BOUNDARY)
	        imin[i] += 8;
 	    if(buffered_boundary_type(rect_boundary_type(intfc,i,1)) == OPEN_BOUNDARY)
	        imax[i] -= 8;
	    smax[i] = imax[i]-imin[i];
	}

	find_Cartesian_coordinates(pp_mynode(), front->pp_grid, picrds);
	
	if(picrds[2] == 0)
	    imin[2] += 8;
	if(picrds[2] == front->pp_grid->gmax[2]-1)
	    imax[2] -= 8;

	set_wall_flag_for_surface(intfc);
	surf = NULL;
	for(s=intfc->surfaces; s && *s; s++)
	    if(!is_wall_surface(*s))
	    {
	        //ASSUME one fluid surface and positive_component is the gas
	        surf = *s;
		compin = positive_component(*s);
	        compout = negative_component(*s);
		break;
	    }
	//WARNING: can not creat bubble with no fluid-gas surf, this should be fixed.
	if(surf == NULL)
	{
	    //WARNING: ASSUME first param is the fluid comp
	    get_default_fluid_comp(&compout,&compin,intfc);
	    printf("#surf==NULL when finding physical surf, "
		   "compin=%d compout=%d\n", compin, compout);
	    //debug_print("bubbles_mpi","Left make_bubbles3d()\n");
	    //printf("#left make_bubbles3d surf==NULL\n\n");
	    //clean_up(0);
	    //return 0;
        }
	
	alloc_site = max3(imax[0],imax[1],imax[2]);
	alloc_site *= alloc_site*alloc_site;

        uni_array(&insert_site, alloc_site, sizeof(INSERT_SITE));
	tri_array(&pos_flag,smax[0],smax[1],smax[2],sizeof(int));

	printf("#chk surf %d \n", surf);
	print_int_vector("imin=", imin, 3, "\n");
	print_int_vector("imax=", imax, 3, "\n");
	//bubble insertion parameters.
	set_bubble_params(&bubbles, gr, compin, compout);

	nsite = 0;
	for (iz = imin[2]; iz < imax[2]; iz++)
        {
            icoords[2] = iz;
	    k = iz-imin[2];
            for (iy = imin[1]; iy < imax[1]; iy++)
            {
                icoords[1] = iy;
		j = iy - imin[1];
                for (ix = imin[0]; ix < imax[0]; ix++)
                {
                    icoords[0] = ix;
		    i = ix - imin[0];
		    
		    coords = Rect_coords(icoords,wave);
		    //not fluid comp
		    //no surface near block
		    if (Rect_comp(icoords,wave) != compout  || 
		        tris_near_cell(coords, gr, intfc) ) 
		    {
		        pos_flag[i][j][k] = INSERT_OFF;
		        continue;
		    }
		    
		    //physical condition satisfied
		    if(!g_bubble_insert_cond(icoords, wave))
		    {
		        pos_flag[i][j][k] = PHYSICAL_OFF;
		        continue;
		    }

		    state = Rect_state(icoords,wave);
		    
		    //assign info at the point
		    insert_site[nsite].press = pressure(state);
		    ft_assign(insert_site[nsite].icoords, icoords, 3*INT);
		    nsite++;
		    if(nsite == alloc_site)
                    {
                        INSERT_SITE *tmp_site;
                        
			alloc_site += max3(imax[0],imax[1],imax[2]);
                        uni_array(&tmp_site, alloc_site, sizeof(INSERT_SITE));
                        for(n = 0; n < nsite; n++)
                            ft_assign(&tmp_site[n], &insert_site[n], sizeof(INSERT_SITE));
                        free(insert_site);
                        insert_site = tmp_site;
                    }
	    
		    pos_flag[i][j][k] = INSERT_ON;
	        }
	    }
	}

	qsort((POINTER)insert_site, nsite, sizeof(INSERT_SITE),press_ascend);

	printf("#insert_site  press= %22.16e  %22.16e\n\n", 
	    insert_site[0].press, insert_site[nsite-1].press);
	//show_pos_flag("pos_flag", pos_flag, imin, imax, gr);
        //add_to_debug("one_bubble");
	
	set_bubble_diameter(&bubbles, gr, dbub[0]);
	for(i = 0; i < nsite; i++)
	{
            if(!test_check_insert_bubble_position(center,inmin,inmax,
	        insert_site[i].icoords, pos_flag, imin,imax, wave, &bubbles))
                    continue;
        
	    //show_pos_flag("#check_insert af", pos_flag, imin, imax, gr);
     	    print_int_vector("#bub inmin=", inmin, 3, "\n"); 
     	    print_int_vector("#bub inmax=", inmax, 3, "\n"); 
	
	    break;
        }

	if(i == nsite)
	    spress = 10.0;
	else
	    spress = insert_site[i].press;

	spmin = spress;
	pp_global_min(&spmin, 1L);
	if(spmin > 0.0)
	{
	    printf("#spmin = %24.15e press is to high\n", spmin);
	    free_these(2,pos_flag,insert_site);
	    return 0;
	}
		
	if(spress < spmin + 0.001)
	    sendnode = pp_mynode();
	else
	    sendnode = pp_numnodes();
	pp_global_imin(&sendnode, 1L);
	if(pp_mynode() == sendnode)
	    pp_send_all(G_HMIN_ID,center,3*FLOAT);
	else
	    pp_recv(G_HMIN_ID,sendnode,center,3*FLOAT);

	print_general_vector("#mpi insert pos=", center, 3, "\n");

	//make_one_bubble.
	
	//set states
	get_gas_ref_state(gst, compin, intfc);
	get_gas_ref_state(ref_st, compout, intfc);
	set_one_bubble_states(gst, ref_st);

	create_bubbles3d(surf,gst,ref_st,wave,front,NULL,
	        center,NULL,NULL,NULL,NULL,&bubbles);
	
	check_print_intfc("After bubble insertion", "bubble_af", 's', 
	       front->interf, -1, 0, NO);
   	
	printf("#left make_bubbles3d bubble made\n\n");

	free_these(2,pos_flag,insert_site);
	debug_print("bubbles_mpi","Left make_bubbles3d()\n");

	return 1;
	
	if(!rect_in_which(center,icoords,gr))
	{
	    printf("#out of grid\n");
	    free_these(2,pos_flag,insert_site);
	    return 0;
	}
	// Create bubbles 
	//inside bubble states	
	get_gas_ref_state(gst, compin, intfc);
	num_bubble = 0;
 	
	set_bubble_diameter(&bubbles, gr, dbub[0]);	
	//set center, inmin, immax.
	if(!check_insert_bubble_position(center,inmin,inmax,
	        icoords, pos_flag, imin,imax, wave, &bubbles))
	{
	    printf("#block is on bdry\n");
	    free_these(2,pos_flag,insert_site);
	    return 0;
	}
       
	num_bubble = 1;
	//show_pos_flag("#check_insert af", pos_flag, imin, imax, gr);

	print_int_vector("#bub inmin=", inmin, 3, "\n"); 
	print_int_vector("#bub inmax=", inmax, 3, "\n"); 
            
	average_states_inside_bubble(ref_st,inmin,inmax,pos_flag,imin,imax,wave);
	
	//gst was initialized in get_gas_ref_state.
	set_state_on_phase_bdry(gst, ref_st);
	    
	printf("#bubble create bf k=%d  num_bubble=%d\n", 0, num_bubble);
	verbose_print_state("ref_st", ref_st);
	verbose_print_state("gst", gst);
	
	//creat a bubble and append it to surf. Also set the pos_flag for 
	//bubble occupied region.
	surf = create_bubbles3d(surf,gst,ref_st,wave,front,pos_flag,
	        center,inmin,inmax,imin,imax,&bubbles);
                       
	//show_pos_flag("bubble create af", pos_flag, imin, imax, gr);
	
	check_print_intfc("After bubble insertion", "bubble_af", 's', 
	       front->interf, -1, 0, NO);
   	printf("#left make_bubbles3d bubble made\n\n");

	free_these(2,pos_flag,insert_site);
        
	debug_print("bubbles_mpi","Left make_bubbles3d()\n");
	return num_bubble;
}

LOCAL bool  test_check_insert_bubble_position(
      float	   *center,
      int          *inmin,	//rectangular region containing bubble
      int          *inmax,
      int          *icoords,
      int	   ***pos_flag,
      int          *dimin,	//rectangular region for insertion
      int          *dimax,
      Wave         *wave,
      BUBBLES      *bubbles)
{
int          i,j,k, ix,iy,iz;
int	     spacing_size, size;
int          ic[MAXD], Rmin[3], Rmax[3];
bool	     is_size_odd;
float	     *coords, r, Rbdry;

      size = bubbles->size;
      spacing_size = bubbles->spacing_size;
      Rbdry = bubbles->Rbdry; //Rbdry should be consistent with spacing_size
      r = bubbles->r;         //r should be consistent with size
      
      // Find the cells the bubble will occupy.
      // If the bubble diameter is odd number of mesh sizes,
      // bubble center is at the center of cell "insert_site.icoords". Otherwise,
      // It is centered at the upper corner of cell "insert_site.icoords".
      is_size_odd = (size%2 == 0 ? NO : YES);
      if(is_size_odd)
      {
          for(i = 0; i < 3; i++)
          {
              inmin[i] = icoords[i] - (size-1)/2;
              inmax[i] = icoords[i] + (size-1)/2;
              center[i] = Rect_coords(icoords, wave)[i];
          }
      }
      else
      {
          for(i = 0; i < 3; i++)
          {
              inmin[i] = icoords[i] - (size/2 - 1);
              inmax[i] = icoords[i] + size/2;
              center[i] = Rect_coords(icoords, wave)[i] +
                              0.5*wave->rect_grid->h[i];
          }
      }

      //The region contains the bubble + spacing
      for(i=0; i<3; i++)
      {
          Rmin[i] = inmin[i] - spacing_size;
          Rmax[i] = inmax[i] + spacing_size;
      }

      for(iz = Rmin[2]; iz <= Rmax[2]; iz++)
      {
          for(iy = Rmin[1]; iy <= Rmax[1]; iy++)
          {
              for(ix = Rmin[0]; ix<= Rmax[0]; ix++)
              {
                  ic[0] = ix; ic[1] = iy; ic[2] = iz;
		  if(ic[0] < dimin[0] || ic[0] >= dimax[0] ||
                     ic[1] < dimin[1] || ic[1] >= dimax[1] ||
                     ic[2] < dimin[2] || ic[2] >= dimax[2])
                      return NO;
                 
		  coords = Rect_coords(ic, wave);
		  if(distance_between_positions(coords,center,3) > Rbdry)
		      continue;
		  i = ix - dimin[0];
		  j = iy - dimin[1];
		  k = iz - dimin[2];
		 
		  //INSERT_OFF means the block is near a surface or other fluid 
		  if(pos_flag[i][j][k] == INSERT_OFF)
		      return NO;

		  if(distance_between_positions(coords,center,3) <= r)
		      if(pos_flag[i][j][k] != INSERT_ON)
		          return NO;
              }
          }
      }

      return YES;
}


LOCAL bool  check_insert_bubble_position(
      float	   *center,
      int          *inmin,	//rectangular region containing bubble
      int          *inmax,
      int          *icoords,
      int	   ***pos_flag,
      int          *dimin,	//rectangular region for insertion
      int          *dimax,
      Wave         *wave,
      BUBBLES      *bubbles)
{
int          i,j,k, ix,iy,iz;
int	     spacing_size, size;
int          ic[MAXD], Rmin[3], Rmax[3];
bool	     is_size_odd;
float	     *coords, r, Rbdry;

      size = bubbles->size;
      spacing_size = bubbles->spacing_size;
      Rbdry = bubbles->Rbdry; //Rbdry should be consistent with spacing_size
      r = bubbles->r;         //r should be consistent with size
      
      // Find the cells the bubble will occupy.
      // If the bubble diameter is odd number of mesh sizes,
      // bubble center is at the center of cell "insert_site.icoords". Otherwise,
      // It is centered at the upper corner of cell "insert_site.icoords".
      is_size_odd = (size%2 == 0 ? NO : YES);
      if(is_size_odd)
      {
          for(i = 0; i < 3; i++)
          {
              inmin[i] = icoords[i] - (size-1)/2;
              inmax[i] = icoords[i] + (size-1)/2;
              center[i] = Rect_coords(icoords, wave)[i];
          }
      }
      else
      {
          for(i = 0; i < 3; i++)
          {
              inmin[i] = icoords[i] - (size/2 - 1);
              inmax[i] = icoords[i] + size/2;
              center[i] = Rect_coords(icoords, wave)[i] +
                              0.5*wave->rect_grid->h[i];
          }
      }

      //The region contains the bubble + spacing
      //only check inside bubble states.
      for(i=0; i<3; i++)
      {
          Rmin[i] = inmin[i] - spacing_size;
          Rmax[i] = inmax[i] + spacing_size;
      }

      //see set_bubble_diameter for Rbdry.
      for(iz = Rmin[2]; iz <= Rmax[2]; iz++)
      {
          for(iy = Rmin[1]; iy <= Rmax[1]; iy++)
          {
              for(ix = Rmin[0]; ix<= Rmax[0]; ix++)
              {
                  ic[0] = ix; ic[1] = iy; ic[2] = iz;
		  if(ic[0] < dimin[0] || ic[0] >= dimax[0] ||
                     ic[1] < dimin[1] || ic[1] >= dimax[1] ||
                     ic[2] < dimin[2] || ic[2] >= dimax[2])
                      return NO;
                 
		  coords = Rect_coords(ic, wave);
		  if(distance_between_positions(coords,center,3) > Rbdry)
		      continue;
		  i = ix - dimin[0];
		  j = iy - dimin[1];
		  k = iz - dimin[2];
		 
		  //INSERT_OFF means the block is near a surface or other fluid 
		  if(pos_flag[i][j][k] == INSERT_OFF)
		      return NO;

		  if(distance_between_positions(coords,center,3) <= r)
		      if(pos_flag[i][j][k] != INSERT_ON)
		          return NO;
              }
          }
      }

      //set pos_flag for inserting position.
      for(iz = inmin[2]; iz <= inmax[2]; iz++)
      {
          for(iy = inmin[1]; iy <= inmax[1]; iy++)
          {
              for(ix = inmin[0]; ix<= inmax[0]; ix++)
              {
                  ic[0] = ix; ic[1] = iy; ic[2] = iz;
		  if(ic[0] < dimin[0] || ic[0] >= dimax[0] ||
                     ic[1] < dimin[1] || ic[1] >= dimax[1] ||
                     ic[2] < dimin[2] || ic[2] >= dimax[2])
                      return NO;
                  
		  i = ix - dimin[0];
		  j = iy - dimin[1];
		  k = iz - dimin[2];
	
		  coords = Rect_coords(ic, wave);
		  if(distance_between_positions(coords,center,3) <= r)
		      pos_flag[i][j][k] = NEW_INSERTED;
              }
	  }
      }
 
      return YES;
}

LOCAL void average_states_inside_bubble(
      Locstate   ref_st,
      int        *imin,   //box contains bubble
      int        *imax,
      int	 ***pos_flag,
      int        *dimin,  //region
      int        *dimax,
      Wave       *wave)
{
      int       i,j,k, iy,ix,iz, ic[MAXD];
      size_t    sizest = wave->sizest;
      Locstate  state;
      int       num_cells;
     
      num_cells = 0;
      for(iz = imin[2]; iz <= imax[2]; iz++)
      {
          for(iy = imin[1]; iy <= imax[1]; iy++)
          {
              for(ix = imin[0]; ix<= imax[0]; ix++)
              {
                  ic[0] = ix; ic[1] = iy; ic[2] = iz;
	          i = ix - dimin[0];
	          j = iy - dimin[1];
	          k = iz - dimin[2];

		  if(pos_flag[i][j][k] != NEW_INSERTED)
		      continue;

                  state = Rect_state(ic,wave);
                  if(num_cells == 0)
		  {
		      ft_assign(ref_st, state, sizest);
		  }
		  else
		  {
                      Dens(ref_st) += Dens(state);
                      Mom(ref_st)[0] += Mom(state)[0];
                      Mom(ref_st)[1] += Mom(state)[1];
                      Mom(ref_st)[2] += Mom(state)[2];
                      Energy(ref_st) += Energy(state);
                  }

		  num_cells++;
              }
          }
      }

      Dens(ref_st) /= num_cells;
      Mom(ref_st)[0] /= num_cells;
      Mom(ref_st)[1] /= num_cells;
      Mom(ref_st)[2] /= num_cells;
      Energy(ref_st) /= num_cells;
}


LOCAL 	void  show_pos_flag(
	char	    *msg,
	int	    ***pos_flag, 
	int	    *imin, 
	int	    *imax, 
	RECT_GRID   *gr)
{
int	ix,iy,iz, i,j,k;
int	icoords[3];

	printf("#show_pos_flag %s\n", msg);
	printf("Grid  gr:\n");
	print_rectangular_grid(gr);
	print_int_vector("imin=", imin, 3, "\n");
	print_int_vector("imax=", imax, 3, "\n");

	printf("\npos_flag=\n");
	for (iz = imin[2]; iz < imax[2]; iz++)
        {
            icoords[2] = iz;
	    k = iz-imin[2];
            printf("iz = %d\n", iz);
	    for (iy = imax[1]-1; iy >= imin[1]; iy--)
            {
                icoords[1] = iy;
		j = iy - imin[1];
                for (ix = imin[0]; ix < imax[0]; ix++)
                {
                    icoords[0] = ix;
		    i = ix - imin[0];
		    printf("%2d", pos_flag[i][j][k]);
		}
		printf("\n");
	    }
	    printf("\n\n");
	}

	printf("#ed show_pos_flag\n");

}

//gst is the gas state, ref_st is the fluid state
LOCAL void  set_state_on_phase_bdry(
	Locstate    gst, 
	Locstate    ref_st)
{
int     stype = state_type(gst);

	set_type_of_state(gst, TGAS_STATE);

	//set bubble states
	Press(gst) = 1.0;
	Dens(gst) = 0.0013;
	Vel(gst)[0] = vel(0, ref_st);
	Vel(gst)[1] = vel(1, ref_st);
	Vel(gst)[2] = vel(2, ref_st);

	set_state(gst, stype, gst);

	//set fluid states
	set_state(ref_st, TGAS_STATE, ref_st);
	Press(ref_st) = 1.0;
	set_state(ref_st, stype, ref_st);
}

LOCAL void  set_one_bubble_states(
	Locstate    gst, 
	Locstate    ref_st)
{
int     stype = state_type(gst);

	//set bubble states
	set_type_of_state(gst, TGAS_STATE);
	Press(gst) = 1.0;
	Dens(gst) = 0.0013;
	Vel(gst)[0] = 0.0;
	Vel(gst)[1] = 0.0;
	Vel(gst)[2] = 0.0;
	set_state(gst, stype, gst);

	//set fluid states
	set_type_of_state(ref_st, TGAS_STATE);
	Press(ref_st) = 1.0;
	Dens(ref_st) = 12.85;
	Vel(ref_st)[0] = 0.0;
	Vel(ref_st)[1] = 0.0;
	Vel(ref_st)[2] = 0.0;
	set_state(ref_st, stype, ref_st);
}


//NEED gas_params_for_comp new version. only Param(st) and state_type(st) 
//will be used later
LOCAL void  get_gas_ref_state(
	Locstate    st, 
	int	    comp,
	INTERFACE   *intfc)
{
Gas_param    *param;

	param = gas_params_for_comp(comp, intfc);
	if(param == NULL)
	{
	    printf("ERROR get_gas_ref_state, the param for comp %d is NULL\n", comp);
	    clean_up(ERROR);
	}
	Init_params(st, param);
	set_type_of_state(st, GAS_STATE);
}

LOCAL int press_ascend(
        const void      *c1,
        const void      *c2)
{
        if ( ((INSERT_SITE*)c1)->press >
             ((INSERT_SITE*)c2)->press ) return 1;
        else return -1;
}               /*end press_ascend*/

//normal point into the center of the ball.
EXPORT	float bubble_func(
	POINTER func_params,
	float *coords)
{
	ELLIP_PARAMS *params;
	const float *cen,*rad;
	float arg;

	params = (ELLIP_PARAMS *)func_params;
	cen = params->cen;
	rad = params->rad;

	arg = 1.0 -
                sqr(coords[0] - cen[0])/sqr(rad[0]) -
                sqr(coords[1] - cen[1])/sqr(rad[1]) -
                sqr(coords[2] - cen[2])/sqr(rad[2]);

	return arg;
}	/* end ellipsoid_func */

EXPORT HYPER_SURF       *g_make_bubble(
	ELLIPSOID       *ellip,
	COMPONENT       compin,
	COMPONENT       compout,
	Front           *front)
{
	SURFACE         *s;
        RECT_GRID       *rgr = front->rect_grid;
	ELLIP_PARAMS	ep;

	ep.cen = ellip->cen;	ep.rad = ellip->rad;
	make_level_surface(rgr,front->interf,compout,compin,bubble_func,
			(POINTER)&ep,&s);
        ellip->hs = Hyper_surf(s);

        //perturb_ellipsoid(s,ellip);
        reset_normal_on_intfc(front->interf);
        //wave_type(s) = ellip->wv_type;
        //surface_tension(s) = ellip->surf_tension;
        //layer_index(ellip->hs) = ellip->layer_index;

        //install_subdomain_bdry_curves(front->interf);
        reset_intfc_num_points(front->interf);
	//untracked_hyper_surf(Hyper_surf(s)) = ellip->untracked;
        return Hyper_surf(s);
}		/* end g_make_ellipsoid */

LOCAL SURFACE* create_bubbles3d(
      SURFACE      *surf,    //surface to append bubbles
      Locstate     gst,      //inside bubble state
      Locstate     ref_st,   //outside bubble state
      Wave         *wave,
      Front        *front,
      int	   ***pos_flag,
      float        *center,
      int          *imin,    //box contains bubbles
      int          *imax,
      int          *dimin,   //region
      int          *dimax,
      BUBBLES      *bubbles)
{
      int          ic[3], ix, iy, iz, i,j,k;
      INTERFACE    *intfc = front->interf;
      INTERFACE    *sv_intfc = current_interface();
      ELLIPSOID    ellip;
      HYPER_SURF   *hs;
      SURFACE      *ellip_surf;
      TRI          *tri;
      POINT        *pt;
      RECT_GRID	   *gr = wave->rect_grid;
      size_t       sizest = wave->sizest;
      Locstate     state;
      Gas_param    *gparam;
      int	   compin=bubbles->compin, compout=bubbles->compout;
      //this tol should be less than tol in set_bubble_diameter
      float	   tol = 0.001, shf;
      
      debug_print("bubbles_mpi","Entered create_bubbles3d()\n");

      set_current_interface(intfc);

      ellip.compin = compin;
      ellip.compout = compout;
      ellip.fpoly = NULL;
      for(k = 0; k < 3; k++)
      {
          ellip.cen[k] = center[k];
          ellip.rad[k] = bubbles->r;
      }
      ellip._make_ellipsoid = g_make_bubble;

      //make a surface
      hs = make_ellipsoid(&ellip, compin, compout, front);
      ellip_surf = Surface_of_hs(hs);

      if(ellip_surf->num_tri == 0)
      {
	    printf("#cre bub no tris.\n");
	    return NULL;
      }
      printf("#creat bubble surf %d %d\n", negative_component(ellip_surf), 
          positive_component(ellip_surf));

      //check component
      if( surf != NULL && 
         (negative_component(ellip_surf) !=  negative_component(surf) ||
          positive_component(ellip_surf) !=  positive_component(surf)) )
      {
          printf("ERROR create_bubbles3d, bubble component is inconsistent.\n");
	  clean_up(ERROR);
      }

      //shift all new points so that they are not grid based.
      shf = tol*min3(gr->h[0], gr->h[1], gr->h[2]);
      for (tri=first_tri(ellip_surf); !at_end_of_tri_list(tri,ellip_surf); tri=tri->next)
      {
          for(i = 0; i < 3; i++)
          {
              pt = Point_of_tri(tri)[i];
	      sorted(pt) = NO;
          }
      }
      for (tri=first_tri(ellip_surf); !at_end_of_tri_list(tri,ellip_surf); tri=tri->next)
      {
          for(i = 0; i < 3; i++)
          {
              pt = Point_of_tri(tri)[i];
	      if(sorted(pt))
	          continue;
	      
	      Coords(pt)[0] += shf;
	      Coords(pt)[1] += 0.842*shf;
	      Coords(pt)[2] += 0.773*shf;
	      sorted(pt) = YES;
          }
      }
      
      //init states on surf, make sure the state_type of gst and ref_st.
      for (tri=first_tri(ellip_surf); !at_end_of_tri_list(tri,ellip_surf); tri=tri->next)
      {
          for(i = 0; i < 3; i++)
          {
              pt = Point_of_tri(tri)[i];
	      ft_assign(left_state(pt), ref_st, sizest);
	      ft_assign(right_state(pt), gst, sizest);
          }
      }
      
      if(surf != NULL)
      {
          //adjoin ellip_surf tri list to surf tri list 
          last_tri(surf)->next = first_tri(ellip_surf);
          first_tri(ellip_surf)->prev = last_tri(surf);
          link_tri_list_to_surface(first_tri(surf),last_tri(ellip_surf),surf);
          //??
          ellip_surf->pos_curves = ellip_surf->neg_curves = NULL;
          delete_surface(ellip_surf);
      }
      else
      {
          surf = ellip_surf;

	  wave_type(surf) = CONTACT;
	  //ASSUME no surface tension
	  surface_tension(surf) = 0.0002;
	  //set_is_not_wall_surface(surf);
      }
      
      install_subdomain_bdry_curves(intfc);
      
      //init interior states for bubble

      gparam = gas_params_for_comp(compin,intfc);
      if(pos_flag == NULL)
      {
      float	*coords, sum;
      int	kmin[3], kmax[3];

      for(i=0; i<3; i++)
      {
          kmin[i] = -gr->lbuf[i];
          kmax[i] = gr->gmax[i] + gr->ubuf[i];
      }
      
      for (iz = kmin[2]; iz < kmax[2]; iz++)
      {
          ic[2] = iz;
          for ( iy = kmin[1]; iy < kmax[1]; iy++)
          {
              ic[1] = iy;
              for ( ix = kmin[0]; ix < kmax[0]; ix++)
              {
                  ic[0] = ix;

		  sum = 0.0;
		  coords = Rect_coords(ic,wave);
		  for(i=0; i<3; i++)
		  {
		      sum += sqr((coords[i]-center[i])/ellip.rad[i]);
		  }
		  if(sum >= 1.0)
		      continue;

                  state = Rect_state(ic,wave);
		  Rect_comp(ic,wave) = compin;
                  set_state(state, state_type(state), gst);
		  
		  //check params
		  if(Params(state) != gparam)
		  {
		      printf("ERROR create_bubbles3d, param inside bubble is inconsistent.\n");
	              clean_up(ERROR);
		  }
              }
          }
      }
      return surf;
      }

      
      for (iz = imin[2]; iz <= imax[2]; iz++)
      {
          ic[2] = iz;
          for ( iy = imin[1]; iy <= imax[1]; iy++)
          {
              ic[1] = iy;
              for ( ix = imin[0]; ix <= imax[0]; ix++)
              {
                  ic[0] = ix;
		  i = ix - dimin[0];
		  j = iy - dimin[1];
		  k = iz - dimin[2];

                  if(pos_flag[i][j][k] != NEW_INSERTED)
		      continue;
                  state = Rect_state(ic,wave);
                  
		  if (debugging("bubbles_mpi"))
                  {
                      (void) printf("Original, comp = %d, bub_comp %d ",
                          Rect_comp(ic,wave),compin);
                      verbose_print_state("state",state);
                  }
                  
		  Rect_comp(ic,wave) = compin;
                  set_state(state, state_type(state), gst);
		  
		  //check params
		  if(Params(state) != gparam)
		  {
		      printf("ERROR create_bubbles3d, param inside bubble is inconsistent.\n");
	              clean_up(ERROR);
		  }
		  //set flag, do not insert bubble again.
		  pos_flag[i][j][k] = INSERT_OFF;
              }
          }
      }
      
      set_current_interface(sv_intfc);
      debug_print("bubbles_mpi","Left create_bubbles3d()\n");
      
      return surf;
}

/*
 * |s|s|b|b|c|b|b|s|s|
 */
//set default bubble insertion params
LOCAL  void set_bubble_params(
	BUBBLES    *bubbles,
	RECT_GRID  *gr,
	int	   compin,
	int	   compout)
{
	bubbles->compin = compin;
	bubbles->compout = compout;
	bubbles->spacing_size = 2; //cells outside a bubble
	set_bubble_diameter(bubbles, gr, 9);	

	printf("#bubble_params %d %d | %d %d | %24.16e %24.16e\n", compin, compout,
	    bubbles->size, bubbles->spacing_size, bubbles->r, bubbles->Rbdry);
}

LOCAL  void set_bubble_diameter(
	BUBBLES    *bubbles,
	RECT_GRID  *gr,
	int	   d)
{
float   hmin, tol, *h = gr->h;

	tol = 0.01;
	hmin = min3(h[0], h[1], h[2]);

	bubbles->size = d;         //diameter
	bubbles->r = bubbles->size*hmin*0.5 + tol*hmin;
	bubbles->Rbdry = bubbles->r + bubbles->spacing_size*hmin;
	
	printf("#bubble diameter set  %d  %24.16e\n", bubbles->size, bubbles->r);
}


LOCAL void reset_front_and_states(
        Wave                      *wave,
        Front                     *front,
        int                       num_bubble)
{
        static int                num_bubbles = 0;
        int                       *iperm;
        INTERFACE                 *intfc = front->interf;
        int                       ix;
        
	// Local statistics 
        if ( num_bubble > 0 )
        {
            intfc->modified = YES;
            printf("%d bubbles are added at time step %d\n",
                   num_bubble, front->step);
        }

        pp_global_isum(&num_bubble,1);

	printf("num_bubble = %d\n", num_bubble);

        // Communication 
        if ( num_bubble > 0 )
	{
	    printf("#bubble inserted\n");
	    //set_stop(YES);
	    //clean_up(0);
	    
	    check_print_intfc("Before scat of bubble insertion", "bub_bf", 'g', 
	           front->interf, 1, -1, NO);

	    printf("#intfc recon %d\n", interface_reconstructed(front->interf));
	    interface_reconstructed(front->interf) = NO;

	    fflush(NULL);
           
	    add_to_debug("skip_comm_comp");
	    if (!scatter_front(front))
            {
                printf("ERROR in make_bubbles, scatter_front() failed: "
                       "number of bubbles = %d\n", num_bubble);
                clean_up(ERROR);
            }
	    remove_from_debug("skip_comm_comp");
	   
	    //clean_up(ERROR);
	    fflush(NULL);
	    
	    check_print_intfc("After scat of bubble insertion", "bubb_ins", 'g', 
	           front->interf, 0, 0, NO);

	    fflush(NULL);

            reinit_hyp_solution_function(wave,front);

            iperm = set_iperm(front->step,front->rect_grid->dim);
            for ( ix = 0; ix < front->rect_grid->dim; ix++)
            {
                if (!scatter_states(wave,front,iperm,ix))
                {
                    screen("ERROR in make_bubbles, scatter_states() "
                           "failed in %d direction, number of bubbles = "
                           "%d\n",ix,num_bubble);
                    clean_up(ERROR);
                }
            }

            num_bubbles += num_bubble;
        }

        pp_gsync();
}


#endif /* if defined(TWOD) || defined(THREED) */

