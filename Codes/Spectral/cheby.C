/*
 *  Solving a 1-D differential equation with spectral methods
 *  (Chebyshev polynomials)
 *
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char cheby_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/09/24 08:38:11  e_gourgoulhon
 *
 * Simple code for illustrating various Chebyshev spectral methods:
 * Galerkin, tau and collocation.
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <iostream.h>
#include <fstream.h>

// C headers
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "map.h"
#include "matrice.h"
#include "graphique.h"

int main() {

	int nn = 5 ;	// Number of Chebyshev coefficients
	
	// Grid of collocation points
	// --------------------------
	int nbr[1] ;
	int nbt[1] ;
	int nbp[1] ;
	nbr[0] = nn ;	// Number of degrees of freedom in r		
	nbt[0] = 1 ;    // Number of degrees of freedom in theta	
	nbp[0] = 1 ;    // Number of degrees of freedom in phi	
	
	int typr[1] ;
	typr[0] = FIN ;		// Type of sampling in r :
					    //   FIN = ordinary Chebychev sampling in [-1,1]
					
	int typt = SYM ;    // Type of sampling in theta
	int typp = SYM ;    // Type of sampling in theta
	
	Mg3d grid(1, nbr, typr, nbt, typt, nbp, typp) ;
	
	cout << "Grid of collocation points: " << grid << endl ;
		
    // Trivial mapping
    // ---------------

    double bornes[2] ;
    bornes[0] = -1. ;
    bornes[1] = 1. ;

    Map_af map(grid, bornes) ;

//    cout << "Mapping : " << map << endl ;

    Mtbl x = map.r ;

    cout << "Collocation points : " << x << endl ;

    // Source of the equation
    // ----------------------

    Valeur ss(grid) ;

    double cc = - exp(1.) / (1. + exp(1.)*exp(1.)) ;

    ss = exp(x) + cc ;

    cout << "Values of the source s at the collation points : "
    	 << endl << ss << endl ;

    // Chebyshev expansion of the source
    // ---------------------------------
    ss.set_base_r(0, R_CHEB) ;
    ss.set_base_t(T_COS_P) ;
    ss.set_base_p(P_COSSIN_P) ;
    ss.coef() ;
    cout << "Coef of the source : " << endl ;
    ss.affiche_seuil(cout, 0, 4, 1.e-15) ;

    // Chebyshev expansion of the source with the large number of coef
    // ---------------------------------------------------------------
	nbr[0] = 33 ;
	Mg3d grid_big(1, nbr, typr, nbt, typt, nbp, typp) ;

    Map_af map_big(grid_big, bornes) ;

    Mtbl x_big = map_big.r ;

    Valeur ss_big(grid_big) ;
    ss_big = exp(x_big) + cc ;
	ss_big.set_base_r(0, R_CHEB) ;
    ss_big.set_base_t(T_COS_P) ;
    ss_big.set_base_p(P_COSSIN_P) ;
    ss_big.coef() ;
    cout << "Coef of the source : " << endl ;
    ss_big.affiche_seuil(cout, 0, 4, 1.e-15) ;

//    des_coef_xi(ss_big, 0, 0, 0, 1.e-14, 0x0, "Coef of s") ;


    // Orthogonal projection: P_N s:
    Mtbl_cf cf_proj(grid, ss.base) ;
    cf_proj.annule_hard() ;
    for (int i=0; i<nn; i++) {
    	cf_proj.set(0,0,0,i) = (*(ss_big.c_cf))(0,0,0,i) ;
    }

    Mtbl_cf cf_alias = *(ss.c_cf) - cf_proj ;
    cout << "Aliasing error : " << cf_alias << endl ;


    ofstream file("ss.d") ;
    file << "#    x            s(x)         Is(x)      Is(x)-s(x)    Is(x)-Ps(x)"  << endl ;
    file.precision(8) ;
    int ndes = 200 ;
    double h = double(2) / double(ndes-1) ;
    for (int i=0; i<ndes; i++) {
    	double xd = -1. + h * i ;
    	double yfd = exp(xd) + cc ;
    	double yid = ss.val_point(0, xd, 0., 0.) ;
    	double ypd = cf_proj.val_point(0, xd, 0., 0.) ;
    	file << xd << "  " << yfd << "  "
    		<< yid << "  " << yid - yfd << "  " << yid - ypd << endl ;
    }
    file.close() ;
	
    file.open("colloc.d") ;
    file.precision(8) ;
    for (int i=0; i<nn; i++) {
    	file << x(0,0,0,i) << "  0." << endl ;    	
    }
    file.close() ;

    // First derivative matrix
    // -----------------------

    Matrice mat_dx(nn,nn) ;
    mat_dx.set_etat_qcq() ;

    Mtbl_cf cf_cheb(grid, ss.base) ;

    for (int i=0; i<nn; i++) {
    	cf_cheb.annule_hard() ; 	// fills with zeros
    	cf_cheb.set(0,0,0,i) = 1. ;
    	
    	cf_cheb.dsdx() ;
    	for (int j=0; j<nn; j++) {
    		mat_dx.set(j,i) = cf_cheb(0,0,0,j) ;
    	}
    }

    cout << "d/dx matrix : " << mat_dx << endl ;

    // Second derivative matrix
    // ------------------------

    Matrice mat_dx2(nn,nn) ;
    mat_dx2.set_etat_qcq() ;

    for (int i=0; i<nn; i++) {
    	cf_cheb.annule_hard() ; 	// fills with zeros
    	cf_cheb.set(0,0,0,i) = 1. ;
    	
    	cf_cheb.d2sdx2() ;
    	for (int j=0; j<nn; j++) {
    		mat_dx2.set(j,i) = cf_cheb(0,0,0,j) ;
    	}
    }

    cout << "d^2/dx^2 matrix : " << mat_dx2 << endl ;

    // Identity matrix
    // ---------------

    Matrice mat_id(nn,nn) ;
    mat_id.set_etat_qcq() ;
    for (int i=0; i<nn; i++) {
    	for (int j=0; j<nn; j++) {
			mat_id.set(i,j) = 0. ;
		}
		mat_id.set(i,i) = 1. ;
	}

    cout << "Identity matrix : " << mat_id << endl ;

    // Full operator matrix
    // --------------------

    Matrice mat_op = mat_dx2 - 4 * mat_dx + 4 * mat_id ;
    cout << "Matrix of the operator d^2/dx^2 - 4 d/dx + 4 Id : "
         << mat_op << endl ;

    //------------------
    // Galerkin method
    //------------------

    // Matrix Chebyshev basis --> Galerkin basis
    Matrice mat_gal(nn,nn-2) ;
    mat_gal.set_etat_qcq() ;
    for (int j=0; j<nn-2; j++) {
    	
    	for (int i=0; i<nn; i++) {
    		mat_gal.set(i,j) = 0 ;
    	}
    	
    	if (j%2 == 0) {
    		mat_gal.set(0,j) = -1. ; 	// - T_0(x)
    	}
    	else {
    		mat_gal.set(1,j) = -1. ;    // - T_1(x)
    	}
		
    	mat_gal.set(j+2,j) = 1. ; 		// T_{j+2}(x)    	
    	
    }

    cout << "Matrix Chebyshev basis --> Galerkin basis : "
    	<< mat_gal << endl ;
    	
    // Transpose
    Matrice mat_tg = mat_gal.transpose() ;

    // Normalization of T_0(x) :
    for (int i=0; i<nn-2; i++) {
		mat_tg.set(i,0) = 2.* mat_tg(i,0) ;
	}
    		
    cout << "Transpose Matrix  : " << mat_tg << endl ;

	// Matrix of the Galerkin linear system
	Matrice mat_lin = mat_tg * mat_op * mat_gal ;
	
    cout <<
    "Matrix of the linear system to solve in the Galerkin method: "
    << endl << mat_lin << endl ;
	
    // Right hand side
    Matrice mat_cfs(nn, 1) ;
    mat_cfs.set_etat_qcq() ;
    for (int i=0; i<nn; i++) {
    	mat_cfs.set(i,0) = (*(ss.c_cf))(0,0,0,i) ;
    }

    Matrice mat_rhs = mat_tg * mat_cfs ;

    cout << "Right-hand side : " << mat_rhs << endl ;

    Tbl rhs_gal(nn-2) ;
    rhs_gal.set_etat_qcq() ;
    for (int i=0; i<nn-2; i++) {
    	rhs_gal.set(i) = mat_rhs(i,0) ;
    }

    Tbl rhs_gal2(mat_rhs) ;

    cout << mat_rhs.get_array() << endl ;

    cout << "Right-hand side (Tbl version) : " << rhs_gal << endl ;
    cout << "Right-hand side (Tbl version 2) : " << rhs_gal2 << endl ;

    // Resolution of the linear system
    mat_lin.set_band(nn-3, nn-3) ;
    mat_lin.set_lu() ;

    Tbl resu = mat_lin.inverse(rhs_gal) ;

    cout << "resu : " << resu << endl ;

	return EXIT_SUCCESS ;

}


















