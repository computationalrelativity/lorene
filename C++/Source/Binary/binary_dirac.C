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

char binary_dirac_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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

     int mermax = 20 ;

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
     
     for(int mer=0; mer<mermax; mer++){

       xi1_old = xi1 ;

       Vector source_xi1 (star1.hij.divergence(star1.flat)) ;
       source_xi1.inc_dzpuis(2) ;
     
       double lambda = 0. ;
       Vector source_reg1 = - (1./3. - lambda) * xi1.divergence(star1.flat)
	 .derive_con(star1.flat)  ;
       source_reg1.inc_dzpuis() ;
       source_xi1 += source_reg1 ;
       Scalar bidon1 (1e-15/rr1/rr1/rr1/rr1) ;
       bidon1.set_domain(0) = 1.e-15 ;
       bidon1.std_spectral_base() ;
       bidon1.inc_dzpuis(4) ;
       for (int i=1; i<=3 ;i++)
	 source_xi1.set(i) = source_xi1(i) + bidon1 ;

       Tenseur source1(star1.mp, 1, CON, star1.mp.get_bvect_cart() ) ;
       source1.set_etat_qcq() ;
       for (int i=0; i<3; i++) 
	 source1.set(i) = Cmp(source_xi1(i+1)) ;
       source1.set_std_base() ;

       Tenseur vect_auxi1 (star1.mp, 1, CON, star1.mp.get_bvect_cart()) ;
       vect_auxi1.set_etat_qcq() ;
       for (int i=0; i<3 ;i++){
	 vect_auxi1.set(i) = 0. ;
       }
       vect_auxi1.set_std_base() ;
       Tenseur scal_auxi1 (star1.mp) ;
       scal_auxi1.set_etat_qcq() ;
       scal_auxi1.set().annule_hard() ;
       scal_auxi1.set_std_base() ;
     
       Param par_xi1 ;
       int niter ;
       par_xi1.add_int(10,  0) ;  // maximum number of iterations
       par_xi1.add_double(1.5,  0) ; // relaxation parameter
       par_xi1.add_double(1.e-16, 1) ; // required precision
       par_xi1.add_int_mod(niter, 0) ; // number of iterations actually used 
  
       Cmp ssjm1khi1 (star1.mp) ;
       ssjm1khi1 = 0 ;
       ssjm1khi1.std_base_scal() ;
       Tenseur ssjm1wbeta1(star1.mp, 1, CON, star1.mp.get_bvect_cart()) ;
       ssjm1wbeta1.set_etat_qcq() ;
       for (int i=0; i<3; i++) 
	 ssjm1wbeta1.set(i) = 0 ;
       ssjm1wbeta1.set_std_base() ;

       par_xi1.add_cmp_mod(ssjm1khi1) ; 
       par_xi1.add_tenseur_mod(ssjm1wbeta1) ; 

       Tenseur resu1(star1.mp, 1, CON, star1.mp.get_bvect_cart() ) ;
       resu1.set_etat_qcq() ;
       resu1.set_std_base() ;
       source1.poisson_vect(lambda, par_xi1, resu1, vect_auxi1, scal_auxi1) ;

       for (int i=1; i<=3; i++) 
	 xi1.set(i) = resu1(i-1) ;

       // Check: has the equation for xi been correctly solved ?
       // --------------------------------------------------------------

       Vector lap_xi1 = (xi1.derive_con(star1.flat)).divergence(star1.flat) 
	 + lambda* xi1.divergence(star1.flat).derive_con(star1.flat) ;
       lap_xi1.inc_dzpuis() ;

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
       double precis = 1e-5 ;
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
     
     for(int mer=0; mer<mermax; mer++){

       xi2_old = xi2 ;

       Vector source_xi2 (star2.hij.divergence(star2.flat)) ;
       source_xi2.inc_dzpuis(2) ;
     
       double lambda = 0. ;
       Vector source_reg2 = - (1./3. - lambda) * xi2.divergence(star2.flat)
	 .derive_con(star2.flat) ;
       source_reg2.inc_dzpuis() ;
       source_xi2 += source_reg2 ;
       Scalar bidon2 (1e-15/rr2/rr2/rr2/rr2) ;
       bidon2.set_domain(0) = 1.e-15 ;
       bidon2.std_spectral_base() ;
       bidon2.inc_dzpuis(4) ;
       for (int i=1; i<=3 ;i++)
	 source_xi2.set(i) = source_xi2(i) + bidon2 ;

       Tenseur source2(star2.mp, 1, CON, star2.mp.get_bvect_cart() ) ;
       source2.set_etat_qcq() ;
       for (int i=0; i<3; i++) 
	 source2.set(i) = Cmp(source_xi2(i+1)) ;
       source2.set_std_base() ;

       Tenseur vect_auxi2 (star2.mp, 1, CON, star2.mp.get_bvect_cart()) ;
       vect_auxi2.set_etat_qcq() ;
       for (int i=0; i<3 ;i++){
	 vect_auxi2.set(i) = 0. ;
       }
       vect_auxi2.set_std_base() ;
       Tenseur scal_auxi2 (star2.mp) ;
       scal_auxi2.set_etat_qcq() ;
       scal_auxi2.set().annule_hard() ;
       scal_auxi2.set_std_base() ;
     
       Param par_xi2 ;
       int niter ;
       par_xi2.add_int(10,  0) ;  // maximum number of iterations
       par_xi2.add_double(1.5,  0) ; // relaxation parameter
       par_xi2.add_double(1.e-16, 1) ; // required precision
       par_xi2.add_int_mod(niter, 0) ; // number of iterations actually used 
  
       Cmp ssjm1khi2 (star2.mp) ;
       ssjm1khi2 = 0 ;
       ssjm1khi2.std_base_scal() ;
       Tenseur ssjm1wbeta2(star2.mp, 1, CON, star2.mp.get_bvect_cart()) ;
       ssjm1wbeta2.set_etat_qcq() ;
       for (int i=0; i<3; i++) 
	 ssjm1wbeta2.set(i) = 0 ;
       ssjm1wbeta2.set_std_base() ;

       par_xi2.add_cmp_mod(ssjm1khi2) ; 
       par_xi2.add_tenseur_mod(ssjm1wbeta2) ; 

       Tenseur resu2(star2.mp, 1, CON, star2.mp.get_bvect_cart() ) ;
       resu2.set_etat_qcq() ;
       source2.poisson_vect(lambda, par_xi2, resu2, vect_auxi2, scal_auxi2) ;

       for (int i=1; i<=3; i++) 
	 xi2.set(i) = resu2(i-1) ;

       // Check: has the equation for xi been correctly solved ?
       // --------------------------------------------------------------

       Vector lap_xi2 = (xi2.derive_con(star2.flat)).divergence(star2.flat) 
	 + lambda* xi2.divergence(star2.flat).derive_con(star2.flat) ;
       lap_xi2.inc_dzpuis() ;
      
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
       double precis = 1e-5 ;
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

     star2.hij_auto = star2.hij_auto + (hij_dirac2 - star2.hij) * 
       star2.decouple ;
     star2.hij_comp = star2.hij_comp + (hij_dirac2 - star2.hij) * 
       (1 - star2.decouple) ;
     star2.hij = hij_dirac2 ;

     //arrete() ;
}
