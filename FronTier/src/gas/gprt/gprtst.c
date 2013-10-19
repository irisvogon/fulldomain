//outside function
void  fprint_bond( FILE *, BOND *);

//void  print_wall_crx(char*,int*,int,int,CRXING*);
//void  print_wall_curve_crx(char*,int*,int,int,CRXING*);
//void  print_wall_curve_crx0(char*,POINT *, int,CRXING*);
//void  check_intfc_curve_connect(INTERFACE *);
//void  check_surface_curve(SURFACE *);
//void  check_curve_connect(CURVE*, SURFACE*);
//void  tecplot_tris_in_box(char*, int, int, int, INTERFACE *);
//bool  is_tri_outside_box(TRI *t, float **fbox);
LOCAL   bool  match_curve_type(CURVE *, int *);
LOCAL   int   g_which_kind_of_curve(CURVE *);
IMPORT  void  tecplot_interior_states(char *, Wave *);
IMPORT	void  test_print_Stencil(Stencil*);
IMPORT  void  print_curve_pointers(SURFACE *surf);
IMPORT  void  fprint_curve_states(FILE *fp, CURVE **c);
IMPORT  void  test_print_curve_states(INTERFACE  *intfc);
IMPORT  void  test_print_closed_curve(INTERFACE  *);
IMPORT  void  fprint_surface_states(FILE *fp, SURFACE **s);
IMPORT  void  fprint_intfc_states(FILE *, INTERFACE *);
IMPORT  void  print_intfc_in_box(INTERFACE *intfc);
IMPORT  void  out_surf(char *fname, INTERFACE *intfc);
IMPORT  void  tecplot_surface_states(const char	*, FILE*,SURFACE*);
IMPORT  void  tecplot_interface_states(const char*, INTERFACE*);
IMPORT  void  tecplot_interior_states(char *, Wave*);
IMPORT  void  check_print_nodes(INTERFACE *);

#define   MAX_BOND_TRI    10
LOCAL  bool  match_curve_type(CURVE *c, int *t)
{
BOND      *b;
BOND_TRI  **btris;
int	  surf_type, t_cnt, i, num;
int	  t_flag[MAX_BOND_TRI];

	for(i=0; i<MAX_BOND_TRI; i++)
	    t_flag[i] = 0;

	b=c->first;
	
	num = 0;
	for(btris = Btris(b); btris && *btris; btris++)
	    num++;
	
	if(num >MAX_BOND_TRI)
	{
	    printf("ERROR: in match_curve_type, too many bond_tri.\n");
	    printf("num = %d, curve = %llu\n", num, curve_number(c));
	    clean_up(ERROR);
	}

	t_cnt = 0;
	for(btris = Btris(b); btris && *btris; btris++)
	{
	    surf_type = wave_type((*btris)->surface);
	    for(i=0; i<num; i++)
	        if(surf_type == t[i] && t_flag[i] == 0)
		{
		    t_flag[i] = 1;
		    t_cnt++;
		    break;
		}
	}

	if(t_cnt == num)
	    return YES;
	return NO;

}


EXPORT  void  check_print_nodes(
	INTERFACE	*intfc)
{
CURVE	**c;
SURFACE **s;

	printf("#check_print_nodes\n");
	for(c=intfc->curves; c && *c; c++)
	{
	    fprint_hsbdry_type(stdout," ",
			       hsbdry_type(*c)," ",intfc);

	    printf("%d nse %d %d  ns %d %d  ne %d %d\n", *c, 
	        (*c)->start, (*c)->end, 
	    	(*c)->start->posn, (*c)->first->start, 
		(*c)->end->posn, (*c)->last->end);
	    if(the_point((*c)->end->posn))
	        printf("#pt found\n");
	}

	//printf("#check in check_print_nodes\n");
	//for(s=intfc->surfaces; s && *s; s++)
	//    check_surface_curve(*s);
}

#define   N_CURVE_TYPE   3
LOCAL   int   g_which_kind_of_curve(CURVE *c)
{
static  int   curve_type[N_CURVE_TYPE][3] = 
	{
	{NEUMANN_BOUNDARY, DIRICHLET_BOUNDARY, PASSIVE_BOUNDARY}, 
	{NEUMANN_BOUNDARY, NEUMANN_BOUNDARY, PASSIVE_BOUNDARY}, 
	{NEUMANN_BOUNDARY, NEUMANN_BOUNDARY, CONTACT} 
	};
int	i;

	for(i=0; i<N_CURVE_TYPE; i++)
	    if(match_curve_type(c, curve_type[i]))
	        return i;

	return -1;
}

EXPORT	void test_print_Stencil(
	Stencil		*sten)
{
	int	   i,j,imax, dim;
	Locstate   st;

	if (sten == NULL)
	{
	    (void) printf("Stencil NULL\n");
	    return;
	}
	(void) printf("reg_stencil %s, prev_reg_stencil %s\n",
		      y_or_n(sten->reg_stencil),
		      y_or_n(sten->prev_reg_stencil));
	
	(void) printf("npts = %d   newcomp = %d\n", sten->npts, sten->newcomp);
	dim = sten->fr->rect_grid->dim;
	imax = (sten->npts % 2) ? sten->npts/2 + 1 : sten->npts/2;
	for (i = -sten->npts/2; i < imax; ++i)
	{
	    (void) printf("sten->hs[%d] = %llu, ",i,
	                  hypersurface_number(sten->hs[i]));/*TODO REMOVE*/
	    (void) printf("sten->crx[%d] = %p, ",i,(POINTER)sten->crx[i]);
	    (void) printf("sten->nc[%d] = %d,\n",i,sten->nc[i]);
	    
	    (void) printf("\tsten->p[%d] = %llu",i,point_number(sten->p[i]));
	    if (sten->p[i] != NULL)
	    {
	    	(void) printf(", ");
	    	print_general_vector(NULL,Coords(sten->p[i]),dim,"");
	    }
	    (void) printf("\n");
	}
	
	printf("states\n");
	for (i = -sten->npts/2; i < imax; ++i)
	{
	    st = sten->st[i];
	    printf("%d  %d  %14.8e  %14.8e  %14.8e  %d  %d  \n", i, st, 
	        Dens(st), Press(st), Mom(st)[0], sten->comp[i], Params(st));
	}

	(void) printf("fr %p newfr %p wave %p newwave %p\n",
		      (POINTER)sten->fr,(POINTER)sten->newfr,
		      (POINTER)sten->wave,(POINTER)sten->newwave);
	(void) printf("\n");
}			/*end print_Stencil*/



EXPORT  void  print_curve_pointers(SURFACE *surf)
{
    CURVE     **c;
    SURFACE   **s;
    int 	i;

	printf("#pos %llu \n", surface_number(surf));
	for(i=0, c=surf->pos_curves; c && *c ; c++)
    	{
	    printf("curve %d ", i);
	    i++;

	    if ((s = (*c)->pos_surfaces) == NULL) 
    	        (void) printf("pos=NULL ");
	    else 
    	        while (*s)
                    (void) printf("#sp %llu ",surface_number(*s++));
	
	    if ((s = (*c)->neg_surfaces) == NULL) 
	        (void) printf("neg=NULL ");
	    else 
	        while (*s)
	            (void) printf("#sn %llu ",surface_number(*s++));
	    printf("\n");	
	}
	   
	printf("#neg\n");
	for(i=0, c=surf->neg_curves; c && *c ; c++)
    	{
	    printf("curve %d ", i);
	    i++;

	    if ((s = (*c)->pos_surfaces) == NULL) 
    	        (void) printf("pos=NULL ");
	    else 
    	        while (*s)
                (void) printf("#sp %llu ",surface_number(*s++));
	
	    if ((s = (*c)->neg_surfaces) == NULL) 
	        (void) printf("neg=NULL ");
	    else 
	        while (*s)
	            (void) printf("#sn %llu ",surface_number(*s++));
	    printf("\n");	
	}
	
}

//POINTER   Get_params(Locstate st)
//{
//	return (POINTER)Params(st);
//}

//no call
EXPORT  void	fprint_curve_states(FILE *fp, CURVE **c)
{
TRI	  *tri;
Locstate  sl, sr;
BOND	  *bond;
POINT	  *p;
BOND_TRI  **btris;
SURFACE   *s;
int	  k=0;

	fprintf(fp, "curve num %llu  curve type %d\n", 
	    curve_number(*c), g_which_kind_of_curve(*c) );
	for(bond = (*c)->first; bond; bond = bond->next)
        {		    
	    k++;
	    fprintf(fp, "##bond=%d\n", k);
	    fprint_bond(fp, bond);
	    if ((btris = Btris(bond)) != NULL)
            {
                for (btris = Btris(bond); btris && *btris; ++btris)
	        {
		    s = (*btris)->surface;
		    fprintf(fp, "#surf = %llu comp %d %d\n", surface_number(s), 
			    negative_component(s), positive_component(s));
			
	    	    sl = left_start_btri_state(*btris);
	    	    sr = right_start_btri_state(*btris);
		    fprintf(fp, "st %d  %d\n", Params(sl), Params(sr));
		    verbose_fprint_state(fp, "s sl", sl);
		    verbose_fprint_state(fp, "s sr", sr);
	        }
		for (btris = Btris(bond);  btris && *btris; ++btris)
                {
		    s = (*btris)->surface;
		    fprintf(fp, "#surf = %llu comp %d %d\n", surface_number(s), 
			    negative_component(s), positive_component(s));
			
	    	        sl = left_end_btri_state(*btris);
	    	        sr = right_end_btri_state(*btris);
			fprintf(fp, "ed %d  %d\n", Params(sl), Params(sr));
	
			verbose_fprint_state(fp, "e sl", sl);
			verbose_fprint_state(fp, "e sr", sr);

                }
		printf("#end btri  %d\n", k);
	    }
	    else
	    {
		printf("#btri NULL\n");
	    }
	    fprintf(fp, "##end\n\n");
	}
	fprintf(fp, "num_bonds=%d\n", k);

}

//no call
EXPORT  void   test_print_curve_states(INTERFACE  *intfc)
{
CURVE  **cc;
int    k=0;
	
	printf("#curve test\n");
	for(cc=intfc->curves; cc && *cc; cc++)
	{
	    if (hsbdry_type(*cc) == SUBDOMAIN_HSBDRY)
	        continue;
	    fprintf(stdout, "#cnumber=%d\n", k);
	    fprint_curve_states(stdout, cc);
	    k++;
	}
}

//no call
EXPORT  void   test_print_closed_curve(INTERFACE  *intfc)
{
CURVE       **cc;
BOND_TRI    **btri1, **btri2;
int    	    k=0;
	
	printf("#closed curve test\n");
	
	for(cc=intfc->curves; cc && *cc; cc++)
	{
	    fprintf(stdout, "#cnumber=%d\n", k);
	    if (!is_closed_curve(*cc))
	    {
	        if((*cc)->start->posn == (*cc)->end->posn)
		{
		    printf("ERROR test_print_closed_curve, same points, different nodes.\n");
		    clean_up(ERROR);
		}
		if((*cc)->first->start == (*cc)->last->end)
		{
		    printf("ERROR test_print_closed_curve, same bond points, different nodes.\n");
		    clean_up(ERROR);
		}
	        continue;
	    }

	    //Now, we have closed curve.
	    for(btri1=Btris((*cc)->first), btri2=Btris((*cc)->last); 
	        btri1 && *btri1 && btri2 && *btri2;  
		btri1++, btri2++)
	    {
	        printf("#left st   %d  %d\n", left_start_btri_state(*btri1), left_end_btri_state(*btri2));
	        printf("#right st  %d  %d\n", right_start_btri_state(*btri1), right_end_btri_state(*btri2));
	        
		if(left_start_btri_state(*btri1) != left_end_btri_state(*btri2))
		{
		    printf("ERROR test_print_closed_curve, different left btri states.\n");
		    clean_up(ERROR);
		}
	        if(right_start_btri_state(*btri1) != right_end_btri_state(*btri2))
		{
		    printf("ERROR test_print_closed_curve, different right btri states.\n");
		    clean_up(ERROR);
		}
	    }
	    k++;
	}
}

//no call
EXPORT  void  fprint_surface_states(FILE *fp, SURFACE **s)
{
TRI	  *tri;
Locstate  sl, sr;
POINT     *p;
int       i, nside, pside;
float	  pl, pr;

	fprintf(fp, "#surface num%d  Hyper_surf_index = %d\n", 
	    surface_number(*s), Hyper_surf_index((*s)));
	fprint_wave_type(fp,"surface->wave_type = ",
		         wave_type(*s),"\n",(*s)->interface);

	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    if(!the_tri(tri))
	        continue;
	    print_tri(tri, (*s)->interface);
	    print_tri_states(tri, Hyper_surf(*s));
	    continue;

	    for (i = 0; i < 3; ++i)
	    {
	    	p = Point_of_tri(tri)[i];
		nside = i;
		pside = Prev_m3(i);
		    
		{
	            fprintf(fp, "%22.16f %22.16f %22.16f\n",Coords(p)[0],Coords(p)[1],Coords(p)[2]);

		    slsr(p,Hyper_surf_element(tri),Hyper_surf(*s),&sl,&sr);
		    pl = is_obstacle_state(sl)? 0.0 : pressure(sl);
		    pr = is_obstacle_state(sr)? 0.0 : pressure(sr);
		    fprintf(fp, "%d  %25.16e  %25.16e\n", Boundary_point(p), pl, pr);
		}
	    }
	    fprintf(fp, "\n");
	}

}

//no call
EXPORT  void	fprint_intfc_states(FILE *fp, INTERFACE *intfc)
{
	SURFACE   **s;
	CURVE     **c;	
	
	for(s=intfc->surfaces; s && *s; s++)
	{
	    fprint_surface_states(fp, s);
	}	

}

EXPORT  void  print_intfc_in_box(INTERFACE *intfc)
{
	SURFACE   **s;
	TRI	  *tri;
	static float  **fbox = NULL;

	if(fbox == NULL)
	    bi_array(&fbox, 2, 3, sizeof(float));

	fbox[0][0] = -21.496;
	fbox[1][0] = -21.484;
	fbox[0][1] = -0.0626;
	fbox[1][1] = -0.0574;
	fbox[0][2] = 85.896;
	fbox[1][2] = 85.903;

	for(s=intfc->surfaces; s && *s; s++)
	{
	    for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	        if(!is_tri_outside_box(tri, fbox))
		{
		    print_tri_coords(tri);
		    printf("\n");
		}
	}	

}

//#bjet2
EXPORT  void out_surf(char *fname, INTERFACE *intfc)
{
	SURFACE    **s;

    	for(s=intfc->surfaces; s && *s; s++)
            if(wave_type(*s) == FIRST_SCALAR_PHYSICS_WAVE_TYPE)
	    	tecplot_surface(fname, NULL, *s);
}


EXPORT  void    tecplot_surface_states(
	const char	*bname,
	FILE		*file,
	SURFACE		*s)
{
	TRI	*tri;
	POINT	*p;
	int	i,npts,ntri,num_tris;
	Locstate  sl, sr;
	float	dist, lenk, max_cos;

	if (bname != NULL)//direct call
	{
	    if ((file = fopen(bname,"w")) == NULL)
	    {
		screen("WARNING in tecplot_surface(), "
		       "can't open %s\n",bname);
		return;
	    }
	    (void) fprintf(file,"TITLE = \"tecplot surface\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
			    "\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\", "
			    "\"PDL0\", \"PDL1\", \"PDR0\", \"PDR1\" \n");
	}
	
	//called from tecplot_interface
	if (file == NULL)
	{
	    screen("ERROR, in tecplot_surface, file is NULL\n");
	    clean_up(ERROR);
	}
	if (!(first_tri(s)))
	{
	    screen("WARNING, first bond of the curve is NULL\n");
	    return;
	}

	//count number of points(npts) and number of tris(ntri)
	for (tri=first_tri(s),ntri=0; !at_end_of_tri_list(tri,s); tri=tri->next,ntri++)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		}
	    }
	}
	//end counting
	
	fprint_wave_type(file, "ZONE T=\"", wave_type(s), "\"", s->interface);
    	fprintf(file, " N=%d E=%d\nF=FEPOINT, ET=TRIANGLE\n",npts,ntri);

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		    
		    slsr(p,Hyper_surf_element(tri),Hyper_surf(s),&sl,&sr);
		    sl = is_obstacle_state(sl) ? sr : sl;
		    sr = is_obstacle_state(sr) ? sl : sr;
	     	   
		    {
		    TRI		*tri2, *tri1, **ptris;
		    int		n, j, k;
		    POINT	*p1, *p2, *p3;
		    float	avep[3], hv1[3], hv2[3], ang;
			
		    for(k=0; k<3; k++)
		        avep[k] = 0.0;
		    lenk = 0.0;
		    max_cos = -1.0;

		    n = set_tri_list_around_point(p, tri, &ptris, s->interface);

		    for(j=0; j<n; j++)
		    {
		        tri1 = ptris[j];
		        p1 = Point_of_tri(tri1)[Next_m3(Vertex_of_point(tri1,p))];
		        p2 = Point_of_tri(tri1)[Prev_m3(Vertex_of_point(tri1,p))];
		        
			for(k=0; k<3; k++)
			    avep[k] += Coords(p1)[k];
			
			lenk += distance_between_positions(Coords(p1), Coords(p2), 3);

			//tst edge
			tri2 = ptris[(j-1+n)%n];
		        p3 = Point_of_tri(tri2)[Next_m3(Vertex_of_point(tri2,p))];
			triangle_height_vec(hv1, Coords(p), Coords(p2), Coords(p1));
			triangle_height_vec(hv2, Coords(p), Coords(p3), Coords(p1));
			ang = Dot3d(hv1,hv2)/(Mag3d(hv1)*Mag3d(hv2));
			if(ang > max_cos)
			    max_cos = ang;
		    }
		    
		    for(k=0; k<3; k++)
		        avep[k] /= n;
		    
		    dist = distance_between_positions(avep, Coords(p), 3);
		    }

		    if(Boundary_point(p))
		    {
		        dist = 0.0; 
			max_cos = -1.0;
		    }

		    if(wave_type(s) != PASSIVE_BOUNDARY)
		        fprintf(file,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e "
				     "%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e "
				     "%15.8e %15.8e %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    	   Coords(p)[1],Coords(p)[2], pressure(sl), pressure(sr), density(sl), density(sr), 
			   vel(0,sl), vel(1,sl), vel(2,sl), vel(0,sr), vel(1,sr), vel(2,sr),
			   pdens(sl)[0], pdens(sl)[1], pdens(sr)[0], pdens(sr)[1]);
		    else
		        fprintf(file,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e "
				     "%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e "
				     "%15.8e %15.8e %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    	 Coords(p)[1],Coords(p)[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}
	    }
	}
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		fprintf(file,"%d ",Index_of_point(Point_of_tri(tri)[i]));
	    }
	    fprintf(file,"\n");
	}

	if (ntri != s->num_tri)
	{
	    printf("WARNING, num of tri in surface is wrong\n"); 
	}
	if (bname != NULL)
	    fclose(file);
}	/* end tecplot_surface_states */

//no call
EXPORT  void    tecplot_interface_states(const char*, INTERFACE	*);

EXPORT  void    tecplot_interface_states(
	const char	*bname,
	INTERFACE	*intfc)
{
	SURFACE	**s;
	CURVE	**cc;
	char    bname1[200];
	FILE	*file;
	int	*i;

	sprintf(bname1, "%s.plt", bname);
	if ((file = fopen(bname1,"w")) == NULL)
	{
	    screen("WARNING in tecplot_interface_states(), "
	           "can't open %s\n",bname1);
	    return;
	}
	(void) fprintf(file,"TITLE = \"tecplot interface\"\n"
		   	    "VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
			    "\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\", "
			    "\"PDL0\", \"PDL1\", \"PDR0\", \"PDR1\" \n");

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    //if(wave_type(*s) != CONTACT)
	    //    continue;
	    //printf("#show contact here\n");
	    tecplot_surface_states(NULL,file,*s);
	}
	fclose(file);
}	/* end tecplot_interface */

void    tecplot_interface_in_ball(const char*, INTERFACE*);

char *  get_directory();

EXPORT  void  tecplot_interior_states(char*, Wave*);

EXPORT  void  tecplot_interior_states(
	char *bname,
	Wave *wave)
{
	const RECT_GRID 	*rgrid = wave->rect_grid;
	Locstate		state;
	int		        xmax = rgrid->gmax[0];
	int			ymax = rgrid->gmax[1];
	int 			zmax = rgrid->gmax[2];
	int 			ix, iy, iz;
	int			icoords[MAXD];
	float			*coords;
	int 			comp;
	FILE			*fp;
	char			s[200];

	sprintf(s,"%s-%d.plt", bname,pp_mynode());
	printf("tecplot_interior_states  file name %s \n",s);

	fp = fopen(s, "w");
	if(fp == NULL)
	{
	    printf("WARNING tecplot_interior_states, can not open file %s\n", s);
	    return; 
	}

	fprintf(fp, "TITLE = \"inner states\" ");
	fprintf(fp, "VARIABLES = \"x\", \"y\", \"z\", \"comp\", \"flag\", ");
	fprintf(fp, "\"dens\", \"pres\", \"entropy\", \"temp\", ");
	fprintf(fp, "\"u\", \"v\", \"w\", " );
//	fprintf(fp, "\"pd0\", \"pd1\", \"pd2\", \"pd3\", \"pd4\", \"pd5\", \"pd6\", \"pd7\", \"pd8\",\"heat flux\", \"xstress\", \"ystress\", \"zstress\",\"dens\",\"temperature\",\"C_P\", \"uvel \", \"vvel\", \"wvel\" \n " );
	fprintf(fp, "\"pd0\", \"pd1\", \"pd2\", \"pd3\", \"pd4\", \"pd5\", \"pd6\", \"pd7\", \"pd8\",\"heat flux\", \"xstress\", \"ystress\", \"zstress\" \n" );
//	fprintf(fp, "\"t0\", \"t1\", \"t2\", \"t3\", \"t4\", \"t5\", \"t6\",\"t7\", \"t8\", " );
//	fprintf(fp, "\"qh0\", \"qh1\", \"qh2\" " );
//	fprintf(fp, "\"qc00\", \"qc01\", \"qc02\", \"qc10\", \"qc11\", \"qc12\" ");
//	fprintf(fp, "\"C\", \"CI\", \"Prt\", ");
//	fprintf(fp,"\"Sct0\", \"Sct1\", \"Sct2\", \"Sct3\", \"Sct4\", \"Sct5\", \"Sct6\", \"Sct7\", \"Sct8\", \"Sct9\" \n");
	fprintf(fp, "ZONE i=%d, j=%d, k=%d\n", xmax, ymax, zmax);

	/* Calculate Kolmogorov Scale */
	/* Calculate Scalar Dissipation Rate */

        Locstate	state_xl,state_xr,state_yl,state_yr,state_zl,state_zr;
	int		icoords_xl[MAXD],icoords_yl[MAXD],icoords_zl[MAXD];
	int		icoords_xr[MAXD],icoords_yr[MAXD],icoords_zr[MAXD];
        float  		dh[3];
        dh[0] = rgrid->h[0];
        dh[1] = rgrid->h[1];
        dh[2] = rgrid->h[2];
	for(iz=0; iz < zmax; iz++)
	    for(iy=0; iy < ymax; iy++)
		  for(ix=0; ix < xmax; ix++)
		  {
		  	icoords[0] = ix;
		  	icoords[1] = iy;
		  	icoords[2] = iz;
			coords = Rect_coords(icoords, wave);
			comp = Rect_comp(icoords, wave);
			state = Rect_state(icoords, wave);

			fprintf(fp, "%f ", coords[0]);
			fprintf(fp, "%f ", coords[1]);
			fprintf(fp, "%f ", coords[2]);
			fprintf(fp, "%d ", comp);
			
			if(is_obstacle_state(state))
			{
			    fprintf(fp, "%d  %15.8e %15.8e %15.8e %15.8e", -1, -1.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e", 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%d  %15.8e %15.8e %15.8e %15.8e", 1, Dens(state),pressure(state), entropy(state),temperature(state));
			    fprintf(fp, "%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e", vel(0, state), vel(1, state), vel(2, state), pdens(state)[0], pdens(state)[1], pdens(state)[2], pdens(state)[3],pdens(state)[4],pdens(state)[5],pdens(state)[6],pdens(state)[7],pdens(state)[8],Boundary_flux(state)[0],Boundary_flux(state)[1], Boundary_flux(state)[2], Boundary_flux(state)[3]);
			}
			fprintf(fp, "\n");
		  }
	
	fclose(fp);
}

