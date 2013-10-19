c
c
c			fortdecs.h
c
c	Copyright 1999 by The University at Stony Brook, All rights reserved.
c
c
#if !defined(_FORTDECS_H)
#define _FORTDECS_H

#if defined(cray)

#   define amax0	max
#   define amax1	max
#   define amin0	min
#   define amin1	min
#   define alog10	log10
#   define alog	log

#endif /* defined(cray) */

#if defined(float)

#   if defined(cray)

#   	define REAL    real*16
#   	define INTEGER integer*16
#   	define COMPLEX complex*32
#   	define NWPC    5

#   else /* defined(cray) */

#   	define REAL real*8
#   	define INTEGER integer*8
#   	if defined(_AIX)
#   	    define COMPLEX double complex
#   	else /* defined(_AIX) */
#   	    define COMPLEX complex*16
#   	endif /* defined(_AIX) */
#   	define NWPC    10

#   endif /* defined(cray) */

#   define Cmplx dcmplx

#   if !defined(cray)

#   	define aint		dint
#   	define anint		dnint
#   	define nini		idnint
#   	define abs		dabs
#   	define cabs		cdabs
#   	define amod		dmod
#   	define sign		dsign
#   	define dim		ddim
#   	define amax0		dmax0
#   	define amax1		dmax1
#   	define amin0		dmin0
#   	define amin1		dmin1
#   	define Real		dble
#   	define Aimag		dimag
#   	define conjg		dconjg
#   	define sqrt		dsqrt
#   	define csqrt		cdsqrt
#   	define exp		dexp
#   	define cexp		cdexp
#   	define alog		dlog
#   	define clog		cdlog
#   	define alog10		dlog10
#   	define sin		dsin
#   	define csin		cdsin
#   	define cos		dcos
#   	define ccos		cdcos
#   	define tan		dtan
#   	define ctan		cdtan
#   	define asin		dasin
#   	define acos		dacos
#   	define atan		datan
#   	define atan2		datan2
#   	define sinh		dsinh
#   	define cosh		dcosh
#   	define tanh		dtanh

#   endif /* !defined(cray) */

#else /* defined(float) */

#   define Cmplx   cmplx
#   define Real    real
#   define Aimag   aimag

#   if defined(cray)
#   	define NWPC    10
#   else /* defined(cray) */
#   	define NWPC    20
#   endif /* defined(cray) */

#endif /* defined(float) */

#endif /* !defined(_FORTDECS_H) */
