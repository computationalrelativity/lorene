/* Computes the Kerr metric in Dirac gauge.
 *
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char kerr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/06 15:09:45  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2006/01/24 16:01:33  f_limousin
 * Version ok for lambda = 0 and lambda = 1/3
 *
 * Revision 1.3  2005/11/10 16:35:30  f_limousin
 * Improve the convergence and the accuracy of the solution
 *
 * Revision 1.2  2005/11/08 14:18:57  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// headers C++
#include "headcpp.h"

// headers C
#include <cstdlib>
#include <cmath>

// headers Lorene
#include "tenseur.h" 
#include "tensor.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "metric.h"
#include "proto.h"
#include "utilitaires.h"

int main(){  
    
    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    int nt, np, nz, nr1, nrp1 ;
    double aa, hh, mm ;

    ifstream fpar("parkerr.d") ;
    fpar.ignore(1000, '\n') ;
    fpar.ignore(1000, '\n') ;
    fpar >> nz; fpar.ignore(1000, '\n');
    fpar >> nt; fpar.ignore(1000, '\n');
    fpar >> np; fpar.ignore(1000, '\n');
    fpar >> nr1; fpar.ignore(1000, '\n');
    fpar >> nrp1; fpar.ignore(1000, '\n');

     // Type of sampling in theta and phi :
    int type_t = SYM ; 
    int type_p = NONSYM ; 

    cout << "total number of domains :   nz = " << nz << endl ;
    cout << "number of points in phi :   np = " << np << endl ;
    cout << "number of points in theta : nt = " << nt << endl ;

    double relax, seuil ;
    int niter ;
    fpar >> relax; fpar.ignore(1000, '\n');
    fpar >> seuil; fpar.ignore(1000, '\n');
    fpar >> niter; fpar.ignore(1000, '\n');

    int* nr_tab = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];

    for (int l=0; l<nz; l++) {
      if (l==1) nr_tab[1] = nr1 ;
      else nr_tab[l] = nrp1 ;
      np_tab[l] = np ; 
      nt_tab[l] = nt ; 
      bornes[l] = pow(2., l-1) ;
    }
    bornes[0] = 0. ;
    bornes[nz] = __infinity ; 
    
    fpar >> hh ; fpar.ignore(1000, '\n') ;  
    fpar >> mm ; fpar.ignore(1000, '\n') ;  
    fpar >> seuil ; fpar.ignore(1000, '\n') ;  

    for (int l=0; l<nz; l++) 
	bornes[l] = bornes[l] * hh * 0.5 ;
	
    cout << "h = " << hh << endl ;
    cout << "M = " << mm << endl ;
    cout << "seuil = " << seuil << endl ;

    fpar.close();
    
    
    // Type of r sampling : 
    int* type_r = new int[nz];
    type_r[0] = RARE ; 
    for (int l=1; l<nz-1; l++) {
	type_r[l] = FIN ; 
    }
    type_r[nz-1] = UNSURR ; 
    
 
    //-----------------------------------------------------------------------
    //		Construction of multi-grid and mapping 
    //-----------------------------------------------------------------------
   
    Mg3d mg(nz, nr_tab, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af mp(mg, bornes) ;

    //-------------------------------
    // Initialisation of h_uu and xi
    //-------------------------------

    Scalar a2(mp) ;
    Scalar b2(mp) ;
 
    // We parametrize Kerr with M and h. 
    aa = sqrt (mm*mm - hh*hh) ;
   
    const Coord& rr = mp.r ;
    const Coord& theta = mp.tet ;

    a2 = 1. + 2.*mm/rr + (3.*mm*mm + aa*aa*cos(2.*theta))/(2.*rr*rr)
	+ (hh*hh*mm)/(2.*pow(rr, 3.)) + pow(hh,4.)/(16.*pow(rr,4.)) ;

    b2 = ( pow(rr,8.) + 4.*mm*pow(rr,7.) + (7.*mm*mm + 
	   aa*aa*cos(theta)*cos(theta))*pow(rr,6.) + mm*(7.*mm*mm+aa*aa)
	   *pow(rr,5.) + (4.*pow(mm,4.) + hh*hh*(3.*hh*hh/4.+aa*aa*sin(theta)
	   *sin(theta))/2.)*pow(rr,4.) + hh*hh*mm*(2.*mm*mm-hh*hh/4.)
	   *pow(rr,3.) + pow(hh,4.)/16.*(7.*mm*mm + aa*aa*cos(theta)
	   *cos(theta))*rr*rr + pow(hh,6.)*mm/16.*rr + pow(hh,8.)/256. ) 
	   / ( pow(rr,8.) + 2.*mm*pow(rr,7.) + (3.*mm*mm + aa*aa
           *cos(2.*theta))/2.*pow(rr,6.) + hh*hh*mm/2.*pow(rr,5.) 
	   + pow(hh,4.)/16.*pow(rr,4.)) ;

    b2.set_outer_boundary(nz-1 ,1.) ;
    b2.std_spectral_base() ;
    b2.set_domain(0) = 1. ;

    Sym_tensor h_uu(mp, CON, mp.get_bvect_cart()) ;
    h_uu.set(1,1) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(2,2) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(3,3) = pow(a2/b2, 2./3.) - 1 ;
    h_uu.set(1,2) = 0. ;
    h_uu.set(1,3) = 0. ;
    h_uu.set(2,3) = 0. ;
    h_uu.annule_domain(0) ;
    h_uu.std_spectral_base() ;
     
    Sym_tensor g_uu(mp, CON, mp.get_bvect_cart()) ;
    g_uu.set(1,1) = 1 / a2 ;
    g_uu.set(2,2) = 1 / a2 ;
    g_uu.set(3,3) = 1 / b2 ;
    g_uu.set(1,2).annule_hard() ;
    g_uu.set(1,3).annule_hard() ;
    g_uu.set(2,3).annule_hard() ;

    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++)
	  g_uu.set(i,j).set_domain(0) = 1. ;
    g_uu.std_spectral_base() ;

    Metric_flat flat(mp.flat_met_cart()) ;
    Metric gamma (g_uu) ;
    Sym_tensor gt_uu(mp, CON, mp.get_bvect_cart()) ;
    for(int i=1; i<=3; i++)
	for(int j=1; j<=i; j++){
	    gt_uu.set(i,j) = flat.con()(i,j) + h_uu(i,j) ;
	}
    gt_uu.std_spectral_base() ;

    // Save initial metric
    Sym_tensor huu_save (h_uu) ;
    Sym_tensor guu_save (g_uu) ;

    Vector xi(mp, CON, mp.get_bvect_cart()) ;
    xi.set(1) = 0. ;
    xi.set(2) = 0. ;
    xi.set(3) = 0. ;
    xi.std_spectral_base() ;

    //----------------------------------------------------
    // Resolution of the vectorial Poisson equation for xi
    //----------------------------------------------------

    Vector source(mp, CON, mp.get_bvect_cart()) ;
   
    Vector xi_jm1(mp, CON, mp.get_bvect_cart()) ;
    for(int i=1; i<=3; i++) 
	xi_jm1.set(i) = 0. ;
    double diff_xi = 1. ;



    for(int mer=0; mer<niter; mer++){

	cout << 
   "========================================================================"
	     << endl ;
	cout << "step = " << mer << "       diff_xi =  "<< diff_xi << endl ; 
	cout << 
   "========================================================================" 
	     << endl ;

	// Function exp(-(r-r_0)^2/sigma^2)
	// --------------------------------
	
	double r0 = mp.val_r(nz-2, 1, 0, 0) ;
	double sigma = 1.*r0 ;
	
	Scalar rrr (mp) ;
	rrr = mp.r ;
	
	Scalar ff (mp) ;
	ff = exp( -(rrr - r0)*(rrr - r0)/sigma/sigma ) ;
	for (int ii=0; ii<nz-1; ii++)
	    ff.set_domain(ii) = 1. ;
	ff.set_outer_boundary(nz-1, 0) ;
	ff.std_spectral_base() ;


	// Source 
	source = h_uu.divergence(flat) ;
	source = source * ff ;
	source.inc_dzpuis(2) ;
 
	source.annule_domain(0) ;
		
	// Printing
	//-----------

	cout << "moyenne de la source pour xi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
/*   
	des_meridian(source(1), 0, 20, "source", 0) ;
	des_meridian(source(2), 0, 20, "source", 1) ;
	des_meridian(source(3), 0, 20, "source", 2) ;
	arrete() ;
*/
	// Resolution of the Poisson equation 
	// ----------------------------------

	double lambda = 1./3. ;
	Vector source_reg = - (1./3. - lambda) * 
	    xi.divergence(flat).derive_con(flat) ;
	source_reg = source_reg * ff ;
	source_reg.inc_dzpuis() ;
	source += source_reg ;

/*
	des_profile(ff, 0, 20, 0, 0) ;
	des_meridian(source(1), 0, 20, "source", 3) ;
	des_meridian(source(2), 0, 20, "source", 4) ;
	des_meridian(source(3), 0, 20, "source", 5) ;
	arrete() ;
*/

	for (int i=1; i<=3; i++) 
	  source.set(i).set_domain(nz-1) = 1.e-15 ; // Because of the dzpuis

	double precision = 1.e-12 ;
	int num_front = 0 ; // index of the intern boundary
	int itermax = 20 ;

	Vector limite (mp, CON, mp.get_bvect_cart()) ;
	limite.std_spectral_base() ;
	Valeur lim_x (mg.get_angu()) ;
	Valeur lim_y (mg.get_angu()) ;
	Valeur lim_z (mg.get_angu()) ;
	lim_x = 0. ;
	lim_y = 0. ;
	lim_z = 0. ;
	lim_x.set_base(limite(1).get_spectral_va().get_base()) ;
	lim_y.set_base(limite(2).get_spectral_va().get_base()) ;
	lim_z.set_base(limite(3).get_spectral_va().get_base()) ;

	
	cout << "Resolution of equation for xi : poisson_vect_frontiere "  
	     << endl << "----------------------------------------------------"
	     << endl ;
	
	
	Tenseur source_p(mp, 1, CON, mp.get_bvect_cart() ) ;
	source_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) 
	  source_p.set(i) = Cmp(source(i+1)) ;
	
	Tenseur resu_p(mp, 1, CON, mp.get_bvect_cart() ) ;
	resu_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) 
	    resu_p.set(i).annule_hard() ;
	resu_p.set_std_base() ;

	poisson_vect_frontiere(lambda, source_p, resu_p, lim_x,
	
		       lim_y, lim_z, num_front, precision, itermax) ;

	for (int i=1; i<=3; i++) 
	  xi.set(i) = resu_p(i-1) ;
	
	// Test
	// ---------

	source.dec_dzpuis() ;
	Vector lap_xi = (xi.derive_con(flat)).divergence(flat) 
	    + lambda* xi.divergence(flat).derive_con(flat) ;
	
	maxabs(lap_xi - source,
	       "Absolute error in the resolution of the equation for beta") ;  
	cout << endl ;
  
       Tbl tdiff_xi_x = diffrel(lap_xi(1), source(1)) ; 
       Tbl tdiff_xi_y = diffrel(lap_xi(2), source(2)) ; 
       Tbl tdiff_xi_z = diffrel(lap_xi(3), source(3)) ; 

       cout << 
	 "Relative error in the resolution of the equation for xi1 : "
	    << endl ; 
       cout << "x component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi_x(l) << "  " ; 
	}
       cout << endl ;
       cout << "y component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi_y(l) << "  " ; 
       }
       cout << endl ;
       cout << "z component : " ;
       for (int l=0; l<nz; l++) {
	 cout << tdiff_xi_z(l) << "  " ; 
       }
       cout << endl ;

       // Xi
       cout << "xi :" << endl ;
       for (int i=1; i<=3; i++){
	 cout << "  Comp. " << i << " :  " ;
	 for (int l=0; l<nz; l++){
	   cout << norme(xi(i)/(nr_tab[0]*nt*np))(l) << " " ;
	 }
	 cout << endl ;
       }
       cout << endl ;
       
       
	// Construction of the new metric quantities :
	//--------------------------------------------
	
	Sym_tensor temp_guu (g_uu.derive_lie(xi)) ;
	temp_guu.dec_dzpuis(2) ;
	g_uu = g_uu + temp_guu ; 
	
	for (int i=1; i<=3; i++)
	    for (int j=1; j<=i; j++)
		g_uu.set(i,j).set_domain(0) = 1. ;
	
	gamma = g_uu ;
	gt_uu = pow(gamma.determinant(), 1./3.) * g_uu ;
	gt_uu.std_spectral_base() ;
	
	for(int i=1; i<=3; i++)
	    for(int j=1; j<=i; j++)
		h_uu.set(i,j) = gt_uu(i,j) - flat.con()(i,j) ;
    
	
	// End of the iteration ?
	
	diff_xi = 0. ;
	Tbl diff (max(abs(xi(1)))) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > diff_xi)
		diff_xi = diff(i) ;
	

	if(diff_xi < seuil) {
	    cout << 
   "========================================================================"
	     << endl ;
	cout << "step = " << mer << "       diff_xi =  "<< diff_xi << endl ; 
	cout << 
   "========================================================================" 
	     << endl ; 
	    mer = niter ;
	}
	
	xi = relax * xi + (1-relax) * xi_jm1 ;
	for (int i=1; i<=3; i++) {
	    xi_jm1.set(i) = xi(i) ;
	}
	
    } // End of the main iteration

    cout << 
    "========================================================================"
	 << endl ;

 
    // Compute h_uu in dirac gauge :
    //-----------------------------

    Sym_tensor guu_dirac (mp, CON, mp.get_bvect_cart()) ;

    guu_dirac = g_uu.derive_lie(xi) ;
    guu_dirac.dec_dzpuis(2) ;
    guu_dirac = guu_dirac + g_uu ;


    cout << "Relative difference between g_uu in dirac gauge and in isotropic gauge "
	 << endl ;
	cout << " Comp 1 1 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(guu_dirac(1,1), guu_save(1,1))(l) << "  " ; 
	}
	cout << endl ;
	cout << " Comp 2 2 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(guu_dirac(2,2), guu_save(2,2))(l) << "  " ; 
	}
	cout << endl ;
	cout << " Comp 3 3 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(guu_dirac(3,3), guu_save(3,3))(l) << "  " ; 
	}
	cout << endl << endl ;


       
 	cout << "norme de g_uu en jauge de dirac :" << endl ;
	for (int i=1; i<=3; i++)
	    for (int j=1; j<=i; j++) {
		cout << "  Comp. " << i << " " << j << " :  " ;
		for (int l=0; l<nz; l++){
		    cout << norme(guu_dirac(i,j)/(nr_tab[0]*nt*np))(l) << " " ;
		}
		cout << endl ;
	    }
	cout << endl ;
	
	// Check of the Dirac gauge
	// ------------------------

	for (int i=1; i<=3; i++)
	    for (int j=1; j<=i; j++)
	      guu_dirac.set(i,j).set_domain(0) = 1. ;

	Metric g_dirac (guu_dirac) ;
	Sym_tensor gtuu_dirac (mp, CON, mp.get_bvect_cart()) ;
	gtuu_dirac = pow(g_dirac.determinant(), 1./3.) * guu_dirac ;
	gtuu_dirac.std_spectral_base() ;

 	cout << "norme de gt_uu en jauge de dirac :" << endl ;
	for (int i=1; i<=3; i++)
	    for (int j=1; j<=i; j++) {
		cout << "  Comp. " << i << " " << j << " :  " ;
		for (int l=0; l<nz; l++){
		    cout << norme(gtuu_dirac(i,j)/(nr_tab[0]*nt*np))(l) <<" " ;
		}
		cout << endl ;
	    }
	cout << endl ;


	Vector d_gtuu_dirac (gtuu_dirac.divergence(flat)) ;
	d_gtuu_dirac.dec_dzpuis(2) ;

	cout << "Is Dirac gauge really satisfied ?" << endl ;
	cout << "Vector H^i" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " : " << norme(d_gtuu_dirac(i)
					     /(nr_tab[0]*nt*np)) << endl ;
	}	

	cout << "For comparaison value norme(h^11_dirac)/dist = " << endl 
	     << norme(gtuu_dirac(1,1)-1)/(nr_tab[0]*nt*np) * 2 / hh ; 

	Vector hh_dirac (huu_save.divergence(flat)) ;
	cout << "For comparaison H^i before computation = " << endl 
	     << norme(hh_dirac(1))/(nr_tab[0]*nt*np) 
	     << endl 
	     << norme(hh_dirac(2))/(nr_tab[0]*nt*np) 
	     << endl 
	     << norme(hh_dirac(3))/(nr_tab[0]*nt*np) 
	     << endl ; 
	
	// Another computation (test)
	// --------------------------
		
	Sym_tensor gtuu_dirac2 (mp, CON, mp.get_bvect_cart()) ;
	Scalar psi4 = pow(a2, 2./3.) * pow(b2, 1./3.) ;
	psi4.std_spectral_base() ;
	gtuu_dirac2 = g_uu.derive_lie(xi) * psi4 + 0.666666666666666 * 
	  xi.divergence(gamma) * gt_uu ; // xi = 0 donc... pas si top !
	gtuu_dirac2.dec_dzpuis(2) ;
	gtuu_dirac2 += gt_uu ;
	
	d_gtuu_dirac = gtuu_dirac2.divergence(flat) ;
	d_gtuu_dirac.dec_dzpuis(2) ;
	cout << "second computation" << endl ;
	cout << "Vector H^i" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " : " << norme(d_gtuu_dirac(i)
					     /(nr_tab[0]*nt*np)) << endl ;
	}	

	cout << "Relative difference between the two computation of h_uu dirac"
	 << endl ;
	cout << " Comp 1 1 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(gtuu_dirac(1,1)-1, gtuu_dirac2(1,1)-1)(l) <<"  " ; 
	}
	cout << endl ;
	cout << " Comp 2 2 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(gtuu_dirac(2,2)-1, gtuu_dirac2(2,2)-1)(l) <<"  " ; 
	}
	cout << endl ;
	cout << " Comp 3 3 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(gtuu_dirac(3,3)-1, gtuu_dirac2(3,3)-1)(l) <<"  " ; 
	}
	cout << endl ;
	
	    
    delete [] nr_tab ; 
    delete [] nt_tab ; 
    delete [] np_tab ;
    delete [] type_r ; 
    delete [] bornes ; 
    
    return EXIT_SUCCESS ; 
        
}  
  
    
