#ifndef __CHEBY_H_
#define __CHEBY_H_

#include "tbl.h"

using namespace Lorene ;

// Returns T_n(x)
double cheby (int n, double x) ;

// Norme of T_n
double norme_cheb (int n) ;

// Returns n colocation points :
Tbl coloc_cheb(int n) ;

// Returns the n weight functions
Tbl weight_cheb(int) ;

// Normalisation of T_n
Tbl gamma_cheb (int) ;

// Computes the coefficients where so is the value at colocation points...
Tbl coef_cheb (const Tbl& so) ;

#endif

