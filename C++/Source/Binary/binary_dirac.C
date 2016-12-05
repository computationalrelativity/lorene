/*
 * Methods of Bin_star::dirac_gauge
 *
 * (see file star.h for documentation)
 *
 */

/*
 *   Copyright (c) 2005 Francois Limousin
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:17:47  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:52:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2006/04/11 14:25:15  f_limousin
 * New version of the code : improvement of the computation of some
 * critical sources, estimation of the dirac gauge, helical symmetry...
 *
 * Revision 1.1  2005/11/08 20:17:01  f_limousin
 * Function used to impose Dirac gauge during an iteration.
 *
 *
 * $Header$ *
 */


// Headers Lorene
#include "tenseur.h"
#include "binary.h"
#include "star.h"
#include "graphique.h"
#include "utilitaires.h"
#include "param.h"


namespace Lorene {
void Binary::dirac_gauge() {

    int nz = star1.mp.get_mg()->get_nzone() ;
    int nr = star1.mp.get_mg()->get_nr(0);
    int nt = star1.mp.get_mg()->get_nt(0);
    int np = star1.mp.get_mg()->get_np(0);

    // Importations 
    // ------------

    // Star 1

    star1.hij_comp.set_triad(star1.mp.get_bvect_cart()) ;
    Sym_tensor comp_hij1(star2.hij_auto) ;
    comp_hij1.change_triad(star2.mp.get_bvect_cart()) ;
    comp_hij1.change_triad(star1.mp.get_bvect_cart()) ;

    assert ( *(star1.hij_comp.get_triad()) == *(comp_hij1.get_triad())) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    	    star1.hij_comp.set(i,j).set_etat_qcq() ;
	    star1.hij_comp.set(i,j).import( (comp_hij1)(i,j) ) ;
	}
    star1.hij_comp.std_spectral_base() ;//set the bases for spectral expansions

     for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++)    
	    star1.hij.set(i,j) = star1.hij_auto(i,j) + star1.hij_comp(i,j) ;

    // Star 2

    star2.hij_comp.set_triad(star2.mp.get_bvect_cart()) ;
    Sym_tensor comp_hij2(star1.hij_auto) ;
    comp_hij2.change_triad(star1.mp.get_bvect_cart()) ;
    comp_hij2.change_triad(star2.mp.get_bvect_cart()) ;

    assert ( *(star2.hij_comp.get_triad()) == *(comp_hij2.get_triad())) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    	    star2.hij_comp.set(i,j).set_etat_qcq() ;
	    star2.hij_comp.set(i,j).import( (comp_hij2)(i,j) ) ;
	}
    star2.hij_comp.std_spectral_base() ;//set the bases for spectral expansions

     for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++)    
	    star2.hij.set(i,j) = star2.hij_auto(i,j) + star2.hij_comp(i,j) ;
     
     // -----------------------------------------
     // Resolution of the Poisson equation for xi
     // -----------------------------------------

     cout << "Function Binary::dirac_gauge()" << endl ;

     // Star 1
     // ----------

     int mermax = 50 ;
     double precis = 1e-5 ;
     double precis_poisson = 1e-14 ;
     double relax_poisson = 1.5 ;
     int mer_poisson = 4 ;

     Scalar rr1 (star1.mp) ;     
     rr1 = star1.mp.r ;
     Scalar rr2 (star2.mp) ;     
     rr2 = star2.mp.r ;

     Vector xi1(star1.mp, CON, star1.mp.get_bvect_cart()) ;
     xi1.set(1) = 0. ;
     xi1.set(2) = 0. ;
     xi1.set(3) = 0. ;
     xi1.std_spectral_base() ;
     Vector xi1_old(xi1) ;
     
     Scalar ssjm1_xi11 (xi1(1)) ;
     Scalar ssjm1_xi12 (xi1(2)) ;
     Scalar ssjm1_xi13 (xi1(3)) ;


     for(int mer=0; mer<mermax; mer++){

       xi1_old = xi1 ;

       // Function exp(-(r-r_0)^2/sigma^2)
       // --------------------------------
       
       double r0_1 = star1.mp.val_r(nz-2, 1, 0, 0) ;
       double sigma = 3.*r0_1 ;
         
       Scalar ff1 (star1.mp) ;
       ff1 = exp( -(rr1 - r0_1)*(rr1 - r0_1)/sigma/sigma ) ;
       for (int ii=0; ii<nz-1; ii++)
	   ff1.set_domain(ii) = 1. ;
       ff1.set_outer_boundary(nz-1, 0) ;
       ff1.std_spectral_base() ;
       
       // Source 
       
       Vector source_xi1 (star1.hij.divergence(star1.flat)) ;
       source_xi1.inc_dzpuis() ;    // dzpuis = 3
     
       double lambda = 0. ;
       Vector source_reg1 = - (1./3. - lambda) * xi1.divergence(star1.flat)
	 .derive_con(star1.flat)  ;
       source_xi1 += source_reg1 ; 
   
       // Resolution of the Poisson equations
  
       Cmp ssjm1xi11 (ssjm1_xi11) ;
       Cmp ssjm1xi12 (ssjm1_xi12) ;
       Cmp ssjm1xi13 (ssjm1_xi13) ;
       ssjm1xi11.set_etat_qcq() ;
       ssjm1xi12.set_etat_qcq() ;
       ssjm1xi13.set_etat_qcq() ;

       Param par_xi11 ;
       int niter ;
       par_xi11.add_int(mer_poisson,  0) ;  // maximum number of iterations
       par_xi11.add_double(relax_poisson,  0) ; // relaxation parameter
       par_xi11.add_double(precis_poisson, 1) ; // required precision
       par_xi11.add_int_mod(niter, 0) ; // number of iterations actually used 
       par_xi11.add_cmp_mod(ssjm1xi11) ; 
 
       Param par_xi12 ;
       par_xi12.add_int(mer_poisson,  0) ;  // maximum number of iterations
       par_xi12.add_double(relax_poisson,  0) ; // relaxation parameter
       par_xi12.add_double(precis_poisson, 1) ; // required precision
       par_xi12.add_int_mod(niter, 0) ; // number of iterations actually used 
       par_xi12.add_cmp_mod(ssjm1xi12) ; 

       Param par_xi13 ;
       par_xi13.add_int(mer_poisson,  0) ;  // maximum number of iterations
       par_xi13.add_double(relax_poisson,  0) ; // relaxation parameter
       par_xi13.add_double(precis_poisson, 1) ; // required precision
       par_xi13.add_int_mod(niter, 0) ; // number of iterations actually used 
       par_xi13.add_cmp_mod(ssjm1xi13) ; 
  
       source_xi1(1).poisson(par_xi11, xi1.set(1)) ;
       source_xi1(2).poisson(par_xi12, xi1.set(2)) ;
       source_xi1(3).poisson(par_xi13, xi1.set(3)) ;

       ssjm1_xi11 = ssjm1xi11 ;
       ssjm1_xi12 = ssjm1xi12 ;
       ssjm1_xi13 = ssjm1xi13 ;

       // Check: has the equation for xi been correctly solved ?
       // --------------------------------------------------------------

       Vector lap_xi1 = (xi1.derive_con(star1.flat)).divergence(star1.flat) 
	 + lambda* xi1.divergence(star1.flat).derive_con(star1.flat) ;

       Tbl tdiff_xi1_x = diffrel(lap_xi1(1), source_xi1(1)) ; 
       Tbl tdiff_xi1_y = diffrel(lap_xi1(2), source_xi1(2)) ; 
       Tbl tdiff_xi1_z = diffrel(lap_xi1(3), source_xi1(3)) ; 
       
       cout << 
	 "Relative error in the resolution of the equation for xi1 : "
	    << endl ; 
       cout << "x component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi1_x(l) << "  " ; 
	}
       cout << endl ;
       cout << "y component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi1_y(l) << "  " ; 
       }
       cout << endl ;
       cout << "z component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi1_z(l) << "  " ; 
       }
       cout << endl ;
       
       
       double erreur = 0 ;
       Tbl diff (diffrelmax (xi1_old(1), xi1(1))) ;
       for (int i=1 ; i<nz ; i++)
	 if (diff(i) > erreur)
	   erreur = diff(i) ;
       
       cout << "Step : " << mer << " Difference : " << erreur << endl ;
       cout << "-------------------------------------" << endl ;
       if (erreur < precis)
	 mer = mermax ;

     }

     // Star 2
     // ----------

     Vector xi2(star2.mp, CON, star2.mp.get_bvect_cart()) ;
     xi2.set(1) = 0. ;
     xi2.set(2) = 0. ;
     xi2.set(3) = 0. ;
     xi2.std_spectral_base() ;
     Vector xi2_old(xi2) ;
     
     Scalar ssjm1_xi21 (xi2(1)) ;
     Scalar ssjm1_xi22 (xi2(2)) ;
     Scalar ssjm1_xi23 (xi2(3)) ;


     for(int mer=0; mer<mermax; mer++){

       xi2_old = xi2 ;

       // Function exp(-(r-r_0)^2/sigma^2)
       // --------------------------------
       
       double r0_2 = star2.mp.val_r(nz-2, 1, 0, 0) ;
       double sigma = 3.*r0_2 ;
         
       Scalar ff2 (star2.mp) ;
       ff2 = exp( -(rr2 - r0_2)*(rr2 - r0_2)/sigma/sigma ) ;
       for (int ii=0; ii<nz-1; ii++)
	   ff2.set_domain(ii) = 1. ;
       ff2.set_outer_boundary(nz-1, 0) ;
       ff2.std_spectral_base() ;
    
       // Source

       Vector source_xi2 (star2.hij.divergence(star2.flat)) ;
       source_xi2.inc_dzpuis() ;    // dzpuis = 3
     
       double lambda = 0. ;
       Vector source_reg2 = - (1./3. - lambda) * xi2.divergence(star2.flat)
	 .derive_con(star2.flat)  ;
       source_xi2 += source_reg2 ;

       // Resolution of the Poisson equations

       Cmp ssjm1xi21 (ssjm1_xi21) ;
       Cmp ssjm1xi22 (ssjm1_xi22) ;
       Cmp ssjm1xi23 (ssjm1_xi23) ;
       ssjm1xi21.set_etat_qcq() ;
       ssjm1xi22.set_etat_qcq() ;
       ssjm1xi23.set_etat_qcq() ;

       Param par_xi21 ;
       int niter ;
       par_xi21.add_int(mer_poisson,  0) ;  // maximum number of iterations
       par_xi21.add_double(relax_poisson,  0) ; // relaxation parameter
       par_xi21.add_double(precis_poisson, 1) ; // required precision
       par_xi21.add_int_mod(niter, 0) ; // number of iterations actually used 
       par_xi21.add_cmp_mod(ssjm1xi21) ; 
 
       Param par_xi22 ;
       par_xi22.add_int(mer_poisson,  0) ;  // maximum number of iterations
       par_xi22.add_double(relax_poisson,  0) ; // relaxation parameter
       par_xi22.add_double(precis_poisson, 1) ; // required precision
       par_xi22.add_int_mod(niter, 0) ; // number of iterations actually used 
       par_xi22.add_cmp_mod(ssjm1xi22) ; 

       Param par_xi23 ;
       par_xi23.add_int(mer_poisson,  0) ;  // maximum number of iterations
       par_xi23.add_double(relax_poisson,  0) ; // relaxation parameter
       par_xi23.add_double(precis_poisson, 1) ; // required precision
       par_xi23.add_int_mod(niter, 0) ; // number of iterations actually used 
       par_xi23.add_cmp_mod(ssjm1xi23) ; 
  
       source_xi2(1).poisson(par_xi21, xi2.set(1)) ;
       source_xi2(2).poisson(par_xi22, xi2.set(2)) ;
       source_xi2(3).poisson(par_xi23, xi2.set(3)) ;

       ssjm1_xi21 = ssjm1xi21 ;
       ssjm1_xi22 = ssjm1xi22 ;
       ssjm1_xi23 = ssjm1xi23 ;

       // Check: has the equation for xi been correctly solved ?
       // --------------------------------------------------------------

       Vector lap_xi2 = (xi2.derive_con(star2.flat)).divergence(star2.flat) 
	 + lambda* xi2.divergence(star2.flat).derive_con(star2.flat) ;
      
       Tbl tdiff_xi2_x = diffrel(lap_xi2(1), source_xi2(1)) ; 
       Tbl tdiff_xi2_y = diffrel(lap_xi2(2), source_xi2(2)) ; 
       Tbl tdiff_xi2_z = diffrel(lap_xi2(3), source_xi2(3)) ; 
       
       cout << 
	 "Relative error in the resolution of the equation for xi2 : "
	    << endl ; 
       cout << "x component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi2_x(l) << "  " ; 
	}
       cout << endl ;
       cout << "y component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi2_y(l) << "  " ; 
       }
       cout << endl ;
       cout << "z component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi2_z(l) << "  " ; 
       }
       cout << endl ;
 

       double erreur = 0 ;
       Tbl diff (diffrelmax (xi2_old(1), xi2(1))) ;
       for (int i=1 ; i<nz ; i++)
	 if (diff(i) > erreur)
	   erreur = diff(i) ;
       
       cout << "Step : " << mer << " Difference : " << erreur << endl ;
       cout << "-------------------------------------" << endl ;
       if (erreur < precis)
	 mer = mermax ;

     }

     // -----------------------------
     // Computation of the new metric
     // -----------------------------

     // Star 1
     // -------

     Sym_tensor guu_dirac1 (star1.mp, CON, star1.mp.get_bvect_cart()) ;
     guu_dirac1 = star1.gamma.con().derive_lie(xi1) ;
     guu_dirac1.dec_dzpuis(2) ;
     guu_dirac1 = guu_dirac1 + star1.gamma.con() ;
     star1.gamma = guu_dirac1 ;
     
     Sym_tensor gtilde_con1(star1.mp, CON, star1.mp.get_bvect_cart()) ;
     Sym_tensor hij_dirac1(star1.mp, CON, star1.mp.get_bvect_cart()) ;

     gtilde_con1 = pow(star1.gamma.determinant(), 1./3.) * guu_dirac1 ;
     gtilde_con1.std_spectral_base() ;
     for(int i=1; i<=3; i++) 
       for(int j=i; j<=3; j++)
	   hij_dirac1.set(i,j) = gtilde_con1(i,j) - star1.flat.con()(i,j) ;

     
     star1.gtilde = gtilde_con1 ;
     star1.psi4 = pow(star1.gamma.determinant(), 1./3.) ;
     star1.psi4.std_spectral_base() ;
     
     cout << "norme de h_uu avant :" << endl ;
     for (int i=1; i<=3; i++)
	 for (int j=1; j<=i; j++) {
	     cout << "  Comp. " << i << " " << j << " :  " ;
	     for (int l=0; l<nz; l++){
		 cout << norme(star1.hij(i,j)/(nr*nt*np))(l) << " " ;
	     }
	     cout << endl ;
	 }
     cout << endl ;

     cout << "norme de h_uu en jauge de dirac :" << endl ;
     for (int i=1; i<=3; i++)
	 for (int j=1; j<=i; j++) {
	     cout << "  Comp. " << i << " " << j << " :  " ;
	     for (int l=0; l<nz; l++){
		 cout << norme(hij_dirac1(i,j)/(nr*nt*np))(l) << " " ;
	     }
	     cout << endl ;
	 }
     cout << endl ;
     

     // Check of the Dirac gauge
     // ------------------------
     
     Vector hh_dirac (star1.hij.divergence(star1.flat)) ;
     cout << "For comparaison H^i before computation = " << endl 
	  << norme(hh_dirac(1))/(nr*nt*np) 
	  << endl 
	  << norme(hh_dirac(2))/(nr*nt*np) 
	  << endl 
	  << norme(hh_dirac(3))/(nr*nt*np) 
	  << endl ; 
     
     Vector hh_dirac_new (hij_dirac1.divergence(star1.flat)) ;
     cout << "Vector H^i after the computation" << endl ;
     for (int i=1; i<=3; i++){
       cout << "  Comp. " << i << " : " << norme(hh_dirac_new(i)
					     /(nr*nt*np)) << endl ;
     }

     star1.hij_auto = star1.hij_auto + (hij_dirac1 - star1.hij) * 
       star1.decouple ;
     star1.hij_comp = star1.hij_comp + (hij_dirac1 - star1.hij) * 
       (1 - star1.decouple) ;
     star1.hij = hij_dirac1 ;


     // Star 2
     // -------

     Sym_tensor guu_dirac2 (star2.mp, CON, star2.mp.get_bvect_cart()) ;
     guu_dirac2 = star2.gamma.con().derive_lie(xi2) ;
     guu_dirac2.dec_dzpuis(2) ;
     guu_dirac2 = guu_dirac2 + star2.gamma.con() ;
     star2.gamma = guu_dirac2 ;
     
     Sym_tensor gtilde_con2(star2.mp, CON, star2.mp.get_bvect_cart()) ;
     Sym_tensor hij_dirac2(star2.mp, CON, star2.mp.get_bvect_cart()) ;

     gtilde_con2 = pow(star2.gamma.determinant(), 1./3.) * guu_dirac2 ;
     gtilde_con2.std_spectral_base() ;
     for(int i=1; i<=3; i++) 
       for(int j=i; j<=3; j++)
	   hij_dirac2.set(i,j) = gtilde_con2(i,j) - star2.flat.con()(i,j) ;

     
     star2.gtilde = gtilde_con2 ;
     star2.psi4 = pow(star2.gamma.determinant(), 1./3.) ;
     star2.psi4.std_spectral_base() ;
     

     star2.hij_auto = star2.hij_auto + (hij_dirac2 - star2.hij) * 
       star2.decouple ;
     star2.hij_comp = star2.hij_comp + (hij_dirac2 - star2.hij) * 
       (1 - star2.decouple) ;
     star2.hij = hij_dirac2 ;

     //arrete() ;
}
}
