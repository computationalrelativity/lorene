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
 * Revision 1.2  2002/09/24 13:15:31  e_gourgoulhon
 *
 * First operational version.
 *
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
#include "utilitaires.h"

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

    Mtbl x = map.r ;

    cout << "Collocation points : " << x << endl ;

    // Source of the equation
    // ----------------------

    Valeur ss(grid) ;

    double cc = - 4. * exp(1.) / (1. + exp(1.)*exp(1.)) ;

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
    cout << "Coef of the source with a large number of Chebyshev polynomials: " << endl ;
    ss_big.affiche_seuil(cout, 0, 4, 1.e-15) ;

//    des_coef_xi(ss_big, 0, 0, 0, 1.e-14, 0x0, "Coef of s") ;


    // Orthogonal projection: P_N s:
    Mtbl_cf cf_proj(grid, ss.base) ;
    cf_proj.annule_hard() ;
    for (int i=0; i<nn; i++) {
    	cf_proj.set(0,0,0,i) = (*(ss_big.c_cf))(0,0,0,i) ;
    }

    Mtbl_cf cf_alias = *(ss.c_cf) - cf_proj ;
    cout << "Aliasing error : " << endl ;
    cout << "--------------" << endl ;
    cf_alias.affiche_seuil(cout) ;

    ofstream file("ss.d") ;
    file << "#    x            s(x)         Is(x)      Is(x)-s(x)    Ps(x) - s(x)     Is(x)-Ps(x)  "  << endl ;
    file.precision(8) ;
    int ndes = 200 ;
    double h = double(2) / double(ndes-1) ;
    for (int i=0; i<ndes; i++) {
    	double xd = -1. + h * i ;
    	double yfd = exp(xd) + cc ;
    	double yid = ss.val_point(0, xd, 0., 0.) ;
    	double ypd = cf_proj.val_point(0, xd, 0., 0.) ;
    	file << xd << "  " << yfd << "  "
    		<< yid << "  " << yid - yfd << "  " << ypd - yfd << " " << yid - ypd << endl ;
    }
    file.close() ;

    file.open("colloc.d") ;
    file.precision(8) ;
    for (int i=0; i<nn; i++) {
    	file << x(0,0,0,i) << "  0." << endl ;
    }
    file.close() ;


    // Exact solution of the equation
    // ------------------------------

    Valeur uu_exact(grid) ;
    uu_exact = exp(x) - sinh(1.) / sinh(2.) * exp(2*x) + cc / 4. ;
    uu_exact.set_base(ss.base) ;
    uu_exact.coef() ;



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

    arrete() ;

    //------------------
    // Galerkin method
    //------------------

    cout << "========================================================================"
         << endl << "                      Galerkin method" << endl
	 << "========================================================================"
	 << endl ;

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
    Matrice mat_ss_cf(nn, 1) ;
    mat_ss_cf.set_etat_qcq() ;
    for (int i=0; i<nn; i++) {
    	mat_ss_cf.set(i,0) = (*(ss.c_cf))(0,0,0,i) ;
    }

    Matrice mat_rhs = mat_tg * mat_ss_cf ;

    cout << "Right-hand side : " << mat_rhs << endl ;

    Tbl rhs_gal(mat_rhs) ;

    cout << "Right-hand side (Tbl version) : " << rhs_gal << endl ;

    // Resolution of the linear system
    mat_lin.set_band(nn-3, nn-3) ;
    mat_lin.set_lu() ;

    Tbl tuu_cf_gal = mat_lin.inverse(rhs_gal) ;

    cout << "Coefficients of the solution onto the Galerkin basis : "
         << endl << tuu_cf_gal << endl ;

    // Verif
    Matrice verif_gal = mat_lin * Matrice(tuu_cf_gal) - mat_rhs ;
    cout << "Check of the resolution of the linear system : " << verif_gal << endl ;

    // Chebyshev coefficients

    Matrice mat_uu_cf_gal(tuu_cf_gal) ;
    Matrice mat_uu_cf_cheb = mat_gal *  mat_uu_cf_gal ;

    Mtbl_cf uu_cf_gal(grid, ss.base) ;
    uu_cf_gal.annule_hard() ;
    for (int i=0; i<nn; i++) {
        uu_cf_gal.set(0,0,0,i) = mat_uu_cf_cheb(i,0) ;
    }

    cout << "Coefficients of the solution onto the Chebyshev basis : " << endl ;
    uu_cf_gal.affiche_seuil(cout, 4, 1.e-15) ;

    // Solution
    Valeur uu_gal(grid) ;
    uu_gal = uu_cf_gal ;
    uu_gal.coef_i() ;

    // Comparison with the exact solution
    cout << "Comparison with the exact solution : " << endl ;
    cout << "----------------------------------" << endl ;
//    cout << "uu_gal : " << uu_gal << endl ;
//    cout << "uu_exact : " << uu_exact << endl ;
    Valeur diff = uu_gal - uu_exact ;
    diff.coef() ;
    diff.affiche_seuil(cout,2) ;

    file.open("sol_gal.d") ;
    file << "#    x            uu_exact(x)        uu(x)      uu(x) - uu_exact(x)  "  << endl ;
    file.precision(8) ;
    for (int i=0; i<ndes; i++) {
    	double xd = -1. + h * i ;
    	double yuu_exact = exp(xd) - sinh(1.) / sinh(2.) * exp(2*xd) + cc / 4.  ;
    	double yuu = uu_gal.val_point(0, xd, 0., 0.) ;
    	file << xd << "  " << yuu_exact << "  " << yuu << "  "
	     << yuu - yuu_exact << endl ;
    }
    file.close() ;

    file.open("uu_gal_colloc.d") ;
    file.precision(8) ;
    for (int i=0; i<nn; i++) {
    	file << x(0,0,0,i) << "  " << uu_gal(0,0,0,i) << endl ;
    }
    file.close() ;

    arrete() ;

    //---------------
    // Tau method
    //---------------

    cout << "========================================================================"
         << endl << "                       Tau method" << endl
	 << "========================================================================"
	 << endl ;

    // Matrix of the 1-form which gives the values at x=1
    Matrice mat_valp1(1,nn) ;
    mat_valp1 = 1. ;

    // Matrix of the 1-form which gives the values at x=-1
    Matrice mat_valm1(1,nn) ;
    mat_valm1.set_etat_qcq() ;

    for (int j=0; j<nn; j++) {
        mat_valm1.set(0,j) = (j%2 == 0) ? 1. : -1. ;
    }

    cout << "Matrix of the 1-form which gives the values at x=1 : " << mat_valp1 << endl ;
    cout << "Matrix of the 1-form which gives the values at x=-1 : " << mat_valm1 << endl ;

    // Insertion in the last two rows of the operator matrix
    Matrice mat_op_tau = mat_op ;
    for (int j=0; j<nn; j++) {
        mat_op_tau.set(nn-2,j) = mat_valp1(0,j) ;
        mat_op_tau.set(nn-1,j) = mat_valm1(0,j) ;
    }

    cout << "Matrix of the operator for the tau method : " << mat_op_tau << endl ;

    // The last two coefficients of the right-hand side are set to zero
    Tbl tss_cf(mat_ss_cf) ;
    tss_cf.set(nn-2) = 0 ;
    tss_cf.set(nn-1) = 0 ;

     // Resolution of the linear system
    mat_op_tau.set_band(nn-1, nn-1) ;
    mat_op_tau.set_lu() ;

    Tbl tuu_cf_tau = mat_op_tau.inverse( tss_cf ) ;

    cout << "Coefficients of the solution by the tau method : "
         << endl << tuu_cf_tau << endl ;

    // Verif
    Matrice verif_tau = mat_op_tau * Matrice(tuu_cf_tau) - Matrice(tss_cf) ;
    cout << "Check of the resolution of the linear system : " << verif_tau << endl ;

    // Coefficients of the solution stored as a Mtbl_cf
    Mtbl_cf uu_cf_tau(grid, ss.base) ;
    uu_cf_tau.annule_hard() ;
    for (int i=0; i<nn; i++) {
        uu_cf_tau.set(0,0,0,i) = tuu_cf_tau(i) ;
    }

    cout << "Coefficients of the solution by the tau method : " << endl ;
    uu_cf_tau.affiche_seuil(cout, 4, 1.e-15) ;

    // Solution
    Valeur uu_tau(grid) ;
    uu_tau = uu_cf_tau ;
    uu_tau.coef_i() ;

    // Comparison with the exact solution
    cout << "Comparison with the exact solution : " << endl ;
    cout << "----------------------------------" << endl ;
//    cout << "uu_tau : " << uu_tau << endl ;
//    cout << "uu_exact : " << uu_exact << endl ;
    diff = uu_tau - uu_exact ;
    diff.coef() ;
    diff.affiche_seuil(cout,2) ;

    file.open("sol_tau.d") ;
    file << "#    x            uu_exact(x)        uu(x)      uu(x) - uu_exact(x)    uu(x) - uu_gal(x)"  << endl ;
    file.precision(8) ;
    for (int i=0; i<ndes; i++) {
    	double xd = -1. + h * i ;
    	double yuu_exact = exp(xd) - sinh(1.) / sinh(2.) * exp(2*xd) + cc / 4.  ;
    	double yuu = uu_tau.val_point(0, xd, 0., 0.) ;
    	double yuu_gal = uu_gal.val_point(0, xd, 0., 0.) ;
    	file << xd << "  " << yuu_exact << "  " << yuu << "  "
	     << yuu - yuu_exact << "  " << yuu - yuu_gal << endl ;
    }
    file.close() ;

    file.open("uu_tau_colloc.d") ;
    file.precision(8) ;
    for (int i=0; i<nn; i++) {
    	file << x(0,0,0,i) << "  " << uu_tau(0,0,0,i) << endl ;
    }
    file.close() ;

    arrete() ;

    //------------------------
    // Pseudo-spectral method
    //------------------------

    cout << "========================================================================"
         << endl << "                Pseudo-spectral method" << endl
	 << "========================================================================"
	 << endl ;

     // Matrix T_j(x_i)

     Matrice mat_cheb_col(nn,nn) ;
     mat_cheb_col.set_etat_qcq() ;

     for (int j=0; j<nn; j++) {
    	cf_cheb.annule_hard() ; 	// fills with zeros
    	cf_cheb.set(0,0,0,j) = 1. ;

        Valeur cheb(grid) ;
        cheb = cf_cheb ;
        cheb.coef_i() ;

        for (int i=0; i<nn; i++) {
    		mat_cheb_col.set(i,j) = cheb(0,0,0,i) ;
    	}
    }

    cout << "Matrix T_{ij} = T_j(x_i) : " << mat_cheb_col << endl ;

    // Matrix of the pseudo-spectral linear system
    Matrice mat_op_psp = mat_cheb_col * mat_op ;

    // Insertion of the boundary conditions at the first and last
    // rows of the matrix

    for (int j=0; j<nn; j++) {
        mat_op_psp.set(0,j) = mat_valm1(0,j) ;
        mat_op_psp.set(nn-1,j) = mat_valp1(0,j) ;
    }

    cout << "Matrix of the operator for the pseudo-spectral method : " << mat_op_psp << endl ;

    // Right-hand side

    Tbl rhs_psp(nn) ;
    rhs_psp.set_etat_qcq() ;

    for (int i=0; i<nn; i++) {
        rhs_psp.set(i) = ss(0,0,0,i) ;
    }

    rhs_psp.set(0) = 0. ;         // boundary conditions
    rhs_psp.set(nn-1) = 0. ;

     // Resolution of the linear system
    mat_op_psp.set_band(nn-1, nn-1) ;
    mat_op_psp.set_lu() ;

    Tbl tuu_cf_psp = mat_op_psp.inverse( rhs_psp ) ;

    cout << "Coefficients of the solution by the pseudo-spectral method : "
         << endl << tuu_cf_psp << endl ;

    // Verif
    Matrice verif_psp = mat_op_psp * Matrice(tuu_cf_psp) - Matrice(rhs_psp) ;
    cout << "Check of the resolution of the linear system : " << verif_psp << endl ;

    // Coefficients of the solution stored as a Mtbl_cf
    Mtbl_cf uu_cf_psp(grid, ss.base) ;
    uu_cf_psp.annule_hard() ;
    for (int i=0; i<nn; i++) {
        uu_cf_psp.set(0,0,0,i) = tuu_cf_psp(i) ;
    }

    cout << "Coefficients of the solution by the psp method : " << endl ;
    uu_cf_psp.affiche_seuil(cout, 4, 1.e-15) ;

    // Solution
    Valeur uu_psp(grid) ;
    uu_psp = uu_cf_psp ;
    uu_psp.coef_i() ;

    // Comparison with the exact solution
    cout << "Comparison with the exact solution : " << endl ;
    cout << "----------------------------------" << endl ;
//    cout << "uu_psp : " << uu_psp << endl ;
//    cout << "uu_exact : " << uu_exact << endl ;
    diff = uu_psp - uu_exact ;
    diff.coef() ;
    diff.affiche_seuil(cout, 2) ;

    file.open("sol_psp.d") ;
    file <<
    "#    x            uu_exact(x)        uu(x)      uu(x)-uu_exact(x)    uu(x)-uu_gal(x)   uu(x)-uu_tau(x)"  << endl ;
    file.precision(8) ;
    for (int i=0; i<ndes; i++) {
    	double xd = -1. + h * i ;
    	double yuu_exact = exp(xd) - sinh(1.) / sinh(2.) * exp(2*xd) + cc / 4.  ;
    	double yuu = uu_psp.val_point(0, xd, 0., 0.) ;
    	double yuu_gal = uu_gal.val_point(0, xd, 0., 0.) ;
    	double yuu_tau = uu_tau.val_point(0, xd, 0., 0.) ;
    	file << xd << "  " << yuu_exact << "  " << yuu << "  "
	     << yuu - yuu_exact << "  " << yuu - yuu_gal << "  " << yuu - yuu_tau << endl ;
    }
    file.close() ;

    file.open("uu_psp_colloc.d") ;
    file.precision(8) ;
    for (int i=0; i<nn; i++) {
    	file << x(0,0,0,i) << "  " << uu_psp(0,0,0,i) << endl ;
    }
    file.close() ;




    return EXIT_SUCCESS ;

}


















