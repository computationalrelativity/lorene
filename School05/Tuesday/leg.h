#ifndef __LEG_H_
#define __LEG_H_

#include "tbl.h"

using namespace Lorene ;

// Legendre n at position x
double leg (int n, double x) ;

// Norme :
double norme_leg (int) ;

// Colocation points and weights
void coloc_poids_leg(int, Tbl& coloc, Tbl& weigths) ;

// Normalisation :
Tbl gamma_leg (int) ;

// Coefficients :
Tbl coef_leg (const Tbl&) ;

// Computational routines...
void legendre (int, double&, double&, double&, double&, double&, double&, double) ;
void gauss_lobato_legendre (int, double*, double*, double, int) ;

#endif

