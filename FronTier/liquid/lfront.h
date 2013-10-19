/*****************************************************************************
 *		lfront.h
 * A interface class for FronTier;	
 *****************************************************************************/
#ifndef LFRONT_H
#define LFRONT_H

typedef struct _Front Front;

class FRONTIER_EB{
	Front *m_pFront;
public:
	FRONTIER_EB(Front &front);
	// geometric parts
	void getCrossing(double *a,double *b,double *crossing,double *normal);	
			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
		
	int  getComponent(double *p);	
	void save(char *filename);
};
#endif
