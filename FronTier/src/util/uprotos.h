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
*			uprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Function prototypes for util libarary
*/

#if !defined(_UPROTOS_H)
#define _UPROTOS_H

#include <cdecs.h>
#include <stdint.h>

/* cleanup.c */
IMPORT	int	is_fatal_signal(int);
IMPORT	void	clean_upp(int);
IMPORT	void	init_clean_up(void(*)(void),void(*)(int));
IMPORT	void	print_call_stack(const char*);
IMPORT	void	print_storage(const char*,const char*);
IMPORT	void	set_dump_core(bool);

/* debug.c */
IMPORT	DEBUG_PARAMS	*init_debug(const DEBUG_MODE);
IMPORT  DEBUG_PARAMS    *default_debug(const DEBUG_MODE);
IMPORT	bool	debugging(const char*);
IMPORT	char	**debugging_names(int*);
IMPORT	void	add_to_debug(const char*);
IMPORT	void	debug_print(const char*,const char*,...);
IMPORT	void	debug_trace(void);
IMPORT	void	remove_from_debug(const char*);
IMPORT	void	set_debug_output(FILE*);

/* error.c */
IMPORT	void	print_errors(void);
IMPORT	void	set_error_immediate(FILE*);
IMPORT	void	log_error(const char*,int,int,const char*);

/* fgetstrin.c */
IMPORT	bool	   fgetstring(FILE*,const char*);
IMPORT	char	   *copy_until(char*,FILE*,FILE*,int);
IMPORT	const char *sgetstring(const char*,const char*);

/* fsort.c */
IMPORT	FILE	   *UncompressAndOpenFile(const char*,const char*);
IMPORT	char	   ***free_infile_list(char***);
IMPORT  const char ***sort_in_files(int,const char**);
IMPORT	void	   CloseAndCleanUpTmpFiles(FILE*);
IMPORT	void	   CloseAndCleanUpAllFiles(void);

/* machine.c */
IMPORT	FT_ENDIAN    ft_endian_type(void);
IMPORT	const char   *ft_endian_name(FT_ENDIAN);
IMPORT	double       d1_mach(int);
IMPORT	truefloat    r1_mach(int);
IMPORT	void	     reverse_string(char*,size_t);
#if defined(_HPUX_SOURCE) || defined(cray)
IMPORT	double	rint(double);
IMPORT	double	copysign(double,double);
IMPORT	double	log1p(double);
IMPORT	double	expm1(double);
#endif /* defined(_HPUX_SOURCE) || defined(cray) */
#if !defined(sun) || (defined(__SUNPRO_C) || defined(__SUNPRO_CC)) || defined(_HPUX_SOURCE) || defined(cray) || (defined(__GNUC__) && !defined(linux))
IMPORT	int irint(double);
#endif /* !defined(sun) || (defined(__SUNPRO_C) || defined(__SUNPRO_CC)) || defined(_HPUX_SOURCE)  || defined(cray) || (defined(__GNUC__) && !defined(linux)) */
IMPORT	char	*get_basename(char*);
IMPORT	char	*get_dirname(char*);

/* bi_array.c*/
IMPORT	void	rotate_matrix(float**,float**,float**,int);
IMPORT	void	rotate_vector(float*,float**,float*,int);

/* other.c */
IMPORT  const char  *ordinal_suffix(int);
IMPORT	const char  *right_flush(int,int);
IMPORT 	const char  *y_or_n(bool);
IMPORT	void	base_and_dir_name(const char*,char**,char**);
IMPORT	void	fprint_line_of_floats(FILE*,int,...);
IMPORT	void	print_bool(const char*,bool,const char*);
IMPORT	void	print_line_of_floats(int,...);

/* output.c */
IMPORT	const IO_TYPE	*open_close_file(const char* const*,int,int,int);
IMPORT	OUTPUT	*save_read_file_variables(FILE*,OUTPUT*);
IMPORT	bool	append_output(FILE*);
IMPORT	bool	check_output(FILE*);
IMPORT	bool	create_directory(const char*,bool);
IMPORT	bool	erase_last_foutput(FILE*);
IMPORT	bool	foutput(FILE*);
IMPORT	bool	is_end_output(FILE*);
IMPORT	bool	is_start_output(FILE*);
IMPORT	bool	next_output(FILE*);
IMPORT	bool	prev_output(FILE*);
IMPORT	bool	output(void);
IMPORT	bool	rewind_read_file(FILE*,OUTPUT*);
IMPORT	const char *next_output_line_containing_string(FILE*,const char*);
IMPORT	int	Fclose(FILE*);
IMPORT	uint64_t size_t2uint64_t(size_t);
IMPORT	void	determine_io_type(FILE*,IO_TYPE*);
IMPORT	void	fprint_io_type(FILE*,const char*,const IO_TYPE*);
IMPORT	void	reset_read_file_variables(OUTPUT*);
IMPORT	void	set_read_endian(FT_ENDIAN);
IMPORT	void	set_read_float_size(size_t);
IMPORT	void	set_reverse_endian(bool);
IMPORT	void	trace_foutput(FILE*);
IMPORT	void	trace_output(void);

/* ppsub.c */
IMPORT	bool	pp_global_status(bool);
IMPORT	bool	pp_iprobe(int,int);
IMPORT	bool	pp_iprobe_any(int);
IMPORT	bool	pp_max_status(bool);
IMPORT	bool	pp_min_status(bool);
IMPORT  int	pp_comm_split(int);
IMPORT	int	pp_finalize(void);
IMPORT	int	pp_init(int*,char***);
IMPORT	int	pp_mynode(void);
IMPORT	int	pp_numnodes(void);
IMPORT	void	EnsureSufficientMessageBufferSize(size_t);
IMPORT	void	pp_abort(int,bool);
IMPORT	void	pp_global_imax(int*,long);
IMPORT	void	pp_global_imin(int*,long);
IMPORT	void	pp_global_ior(int*,long);
IMPORT	void	pp_global_isum(int*,long);
IMPORT	void	pp_global_lmax(long*,long);
IMPORT	void	pp_global_lmin(long*,long);
IMPORT	void	pp_global_lor(long*,long);
IMPORT	void	pp_global_lsum(long*,long);
IMPORT	void	pp_global_max(float*,long);
IMPORT	void	pp_global_min(float*,long);
IMPORT	void	pp_global_sum(float*,long);
IMPORT	void	pp_gsync(void);
IMPORT	void	set_MSG_BUF_SIZE(size_t);
IMPORT	void	set_pp_recv_num_retries(unsigned);
IMPORT	void	set_pp_recv_timeout(unsigned);
IMPORT	void	set_pp_recv_wait_interval(unsigned);
IMPORT	void	u_pp_all_gather(POINTER,int,POINTER,int,const char*,int);
IMPORT  void    pp_f_allgatherv(float*,long,float*,long*);
IMPORT  void    pp_l_allgather(long*,long,long*,long);
IMPORT	void	u_pp_recv(int,int,POINTER,size_t,const char*,int);
IMPORT	void	u_pp_recv_any(int,POINTER,size_t,const char*,int);
IMPORT	void	u_pp_send(int,POINTER,size_t,int,const char*,int);
IMPORT	void	u_pp_send_all(int,POINTER,size_t,const char*,int);
IMPORT  void    u_pp_bcast(int,POINTER,size_t,const char*,int);
#if defined(__MPI__)
IMPORT	bool pp_test(MPI_Request*);
IMPORT	void	pp_wait(MPI_Request*);
IMPORT	void	u_pp_irecv(int,int,POINTER,size_t,MPI_Request*,const char*,int);
IMPORT	void	u_pp_isend(int,POINTER,size_t,int,MPI_Request*,const char*,int);
#endif /* defined(__MPI__) */

/* quad.c*/
IMPORT	float	dqng(float(*)(float,POINTER),POINTER,double,double,double,
	             double,float*,int*,QUADRATURE_STATUS*);
IMPORT	float	SimpRule(float(*)(float,POINTER),POINTER,double,double,double,
	                 double,float*,int*,QUADRATURE_STATUS*);

/* roots.c*/
IMPORT	bool	bisection_find_root(bool(*)(float,float*,POINTER),POINTER,
				    float,float*,float,float,float,float);
IMPORT	bool	find_root(bool(*)(float,float*,POINTER),POINTER,float,
			  float*,float,float,float,float);
IMPORT	bool	find_separation_point(bool(*)(float,float*,POINTER),
				      POINTER,float,float*,float,float,
				      float*,float*,float);
IMPORT	bool	search_harder_for_root(bool(*)(float,float*,POINTER),
				       POINTER,float,float*,float,float,
				       float*,float*,float,float,int,
				       float,float);
IMPORT	void	print_function_values(bool(*)(float,float*,POINTER),
				      POINTER,float,float,float,int,
				      const char*,FILE*);

/* runga.c*/
IMPORT	bool	runga_kutta(double,double*,double,double*,double*,int,
			    bool(*)(double,double*,double*,int,POINTER),
			    double,POINTER);

/* screen.c */
IMPORT	char	*Gets(char*);
IMPORT	int	Scanf(const char*,...);
IMPORT	void	screen(const char*,...);
IMPORT	void	screen_print_long_string(const char*);

/* simpleio.c */
IMPORT	bool	fread_bool(FILE*);
IMPORT	bool	is_binary_output(void);
IMPORT	float	fread_float(const char*,const IO_TYPE*);
IMPORT	float	read_print_float(const char*,float,const IO_TYPE*);
IMPORT	int	fscan_float(FILE*,float*);
IMPORT	int	sscan_float(const char*,float*);
IMPORT	size_t	read_binary_int_array(int*,size_t,const IO_TYPE*);
IMPORT	size_t	read_binary_real_array(float*,size_t,const IO_TYPE*);
IMPORT	size_t	read_binary_uint64_t_array(uint64_t*,size_t,const IO_TYPE*);
IMPORT	size_t	read_binary_uint_array(unsigned int*,size_t,const IO_TYPE*);
IMPORT	uint64_t u_ptr2ull(void*);
IMPORT	void	fprint_double_matrix(FILE*,const char*,int,int,
				     double**,const char*);
IMPORT	void	fprint_double_vector(FILE*,const char*,int,double*,const char*);
IMPORT	void	fprint_double_vector_as_matrix(FILE*,const char*,int,int,
					       double*,const char*);
IMPORT	void	fprint_float(FILE*,const char*,float,const char*);
IMPORT	void	fprint_float_matrix(FILE*,const char*,int,int,float**,
				    const char*);
IMPORT	void	fprint_float_vector(FILE*,const char*,int,float*,const char*);
IMPORT	void	fprint_float_vector_as_matrix(FILE*,const char*,int,int,
					      float*,const char*);
IMPORT	void	fprint_matrix(FILE*,const char*,int,int,float**,const char*);
IMPORT	void	fprint_int_matrix(FILE*,const char*,int,int,int**,const char*);
IMPORT	void	fprint_int_vector_as_matrix(FILE*,const char*,int,int,int*,
					    const char*);
IMPORT	void	fprint_vector(FILE*,const char*,int,float*,const char*);
IMPORT	void	fprint_vector_as_matrix(FILE*,const char*,int,int,float*,
					const char*);
IMPORT	void	fprint_vector_of_floats(FILE*,int,float*);
IMPORT	void	fwrite_float(FILE*,const char*,float,bool,const char*,
			     const char*);
IMPORT	void	print_machine_parameters(FILE*);
IMPORT	void	print_title(FILE*,const char*);
IMPORT	void	set_binary_output(bool);
IMPORT	void	stripcomm(char*,const char*);

/*sphhar.c*/
IMPORT	double	*SphericalHarmonic(double*,int,int,double,double);
IMPORT	double	SphericalHarmonic_i(int,int,double,double);
IMPORT	double	SphericalHarmonic_r(int,int,double,double);
IMPORT	double	SphericalHarmonic_s(int,int,double,double,double);
IMPORT	double	NALegendre(int,int,double);

/* times.c */
IMPORT	void	print_execution_times(void);
IMPORT	void	start_clock(const char*);
IMPORT	void	stop_clock(const char*);
IMPORT	void	cpu_time(const char*);
IMPORT	double	cpu_seconds(void);
IMPORT	float	real_time(void);
IMPORT	char	*date_string(void);
IMPORT  void    add_time_start(int);
IMPORT  void    add_time_end(int);
IMPORT  void    add_time_clear(int);
IMPORT  float   add_time_get(int);

/* uinit.c */
IMPORT	void	init_prompting_and_debugging(INIT_DATA*);
IMPORT  void    init_default_debugging(INIT_DATA*);

/* uni_arraymalloc.c */
IMPORT	POINTER	array_T(const char*,POINTER*,int,...);
IMPORT	int	free_from_T(POINTER);
IMPORT	int	get_vmalloc_storage_use(void);
IMPORT	void	alloc_view(FILE*);
IMPORT	void	f_ree(POINTER,const char*);
IMPORT	void	free_these(int,...);
IMPORT	void	long_alloc_view(FILE*);

/* vtk.c */
#if defined(float)
#undef float
#define not_float
#endif
IMPORT  float 	endian_float_swap(float);
#if defined(not_float)
#define float double
#undef not_float
#endif
IMPORT  double 	endian_double_swap(double);
IMPORT  int   	endian_int_swap(int);
IMPORT  bool	hardware_is_little_endian(void);
IMPORT  int	count_digits(int);

/* fft.c */
IMPORT 	bool 	fft2d(COMPLEX**,int,int,int);
IMPORT 	bool 	fft(int,int,float*,float*);
IMPORT 	bool 	Powerof2(int,int*,int*);

/* odepack prototypes */
#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
IMPORT  int    ode_solver(LSODE_FUNC*,float*,float,float,int,
                          float,float,float,int*,LSODE_JAC*);
IMPORT  int    ode_solver_VM(LSODE_FUNC*,float*,float,float,int,
                          float,float,float,int*,LSODE_JAC*,float*);
#   if defined(cray)
#	if defined(float)
#	    define lsode LSODE
#	else /* defined(float) */
#	    define lsode SLSODE
#	endif /* defined(float) */
#   else /* defined(cray) */
#	if !defined(float)
#	    define lsode slsode
#	endif /* !defined(float) */
#   endif /* defined(cray) */
    FORTRAN void    FORTRAN_NAME(lsode)(LSODE_FUNC*,int*,float*,float*,
	                                float*,int*,float*,float*,
					int*,int*,int*,float*,int*,int*,
					int*,LSODE_JAC*,int*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

/* linpak prototypes */
#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
#   if defined(cray)
#	define	dgtsl	DGTSL
#	define	sgtsl	SGTSL
#   endif /* defined(cray) */
#   if defined(float)
#	define	gtsl	FORTRAN_NAME(dgtsl)
#   else /* defined(float) */
#	define	gtsl	FORTRAN_NAME(sgtsl)
#   endif /* defined(float) */
    FORTRAN	void gtsl(int*,float*,float*,float*,float*,int*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

/* dierckx prototypes */
#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
#   if defined(cray)
#	define curfit CURFIT
#	define splev SPLEV
#   endif /* defined(cray) */
FORTRAN  void FORTRAN_NAME(curfit)(int*,int*,float*,float*,float*,
	                           float*,float*,int*,float*,
				   int*,int*,float*,float*,float*,
				   float*,int*,int*,int*);
FORTRAN  void FORTRAN_NAME(splev)(float*,int*,float*,int*,float*,float*,
	                          int*,int*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

#endif /* !defined(_UPROTOS_H) */
