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


char Kerr_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// headers C++
#include "headcpp.h"

// headers C
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "cmp.h"
#include "tenseur.h" 
#include "tensor.h"
#include "scalar.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "metric.h"
#include "proto.h"

int main(){  
    
    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    int nt, np, nz ;
    double aa, hh, mm, seuil ;

    ifstream fpar("parkerr.d") ;
    fpar.ignore(1000, '\n') ;
    fpar.ignore(1000, '\n') ;
    fpar >> nz; fpar.ignore(1000, '\n');
    fpar >> nt; fpar.ignore(1000, '\n');
    fpar >> np; fpar.ignore(1000, '\n');

    cout << "total number of domains :   nz = " << nz << endl ;
    cout << "number of points in phi :   np = " << np << endl ;
    cout << "number of points in theta : nt = " << nt << endl ;

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];

    fpar.ignore(1000, '\n');
    for (int l=0; l<nz; l++) {
	fpar >> nr[l]; 
	fpar >> bornes[l];     fpar.ignore(1000, '\n');
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    cout << "number of points in r in domain 0 :  nr = " << nr[0] << endl ;
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
    
    // Type of sampling in theta and phi :
    int type_t = SYM ; 
    int type_p = SYM ; 
 
    //-----------------------------------------------------------------------
    //		Construction of multi-grid and mapping 
    //-----------------------------------------------------------------------
   
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af mp(mg, bornes) ;

    //-------------------------------
    // Initialisation of h_uu and xsi
    //-------------------------------

    Scalar a2(mp) ;
    Scalar b2(mp) ;
 
    // We parametrize Kerr with M and h. 
    aa = sqrt (mm*mm - hh*hh) ;
   
    const Coord& rr = mp.r ;
    const Coord& theta = mp.tet ;

    a2 = 1. + 2.*mm/rr + (3.*mm*mm + aa*aa*cos(2.*theta))/(2.*rr*rr)
	+ (hh*hh*mm)/(2.*pow(rr, 3.)) + pow(hh,4.)/(16.*pow(rr,4.)) ;

    b2 = ( pow(rr,8.) + 4*mm*pow(rr,7.) + (7.*mm*mm + 
	   aa*aa*cos(theta)*cos(theta))*pow(rr,6.) + mm*(7.*mm*mm+aa*aa)
	   *pow(rr,5.) + (4.*pow(mm,4.) + hh*hh*(3.*hh*hh/4.+aa*aa*sin(theta)
	   *sin(theta))/2.)*pow(rr,4.) + hh*hh*mm*(2.*mm*mm-hh*hh/4.)
	   *pow(rr,3.) + pow(hh,4.)/16.*(7*mm*mm + aa*aa*cos(theta)
	   *cos(theta))*rr*rr + pow(hh,6.)*mm/16.*rr + pow(hh,8.)/256 ) 
	   / ( pow(rr,8.) + 2*mm*pow(rr,7.) + (3*mm*mm + aa*aa
           *cos(2*theta))/2.*pow(rr,6.) + hh*hh*mm/2.*pow(rr,5.) 
	   + pow(hh,4)/16.*pow(rr,4)) ;


    Sym_tensor h_uu(mp, CON, mp.get_bvect_spher()) ;
    
    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++){
	    if(i != j){
		h_uu.set(i,j) = 0 ;
	    }   
	}

    h_uu.set(1,1) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(2,2) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(3,3) = pow(a2/b2, 2./3.) - 1 ;
    h_uu.annule_domain(0) ;

    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++){
	    h_uu.set(i,j).set_outer_boundary(nz-1 ,0.) ;
	}

    Metric_flat flat(mp, mp.get_bvect_spher()) ;

 
    Sym_tensor g_uu(mp, CON, mp.get_bvect_spher()) ;
    
    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++){
	    if(i != j){
		g_uu.set(i,j) = 0 ;
	    }   
	}

    g_uu.set(1,1) = 1 / a2 ;
    g_uu.set(2,2) = 1 / a2 ;
    g_uu.set(3,3) = 1 / b2 ;
    g_uu.annule_domain(0) ;

    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++){
	    if(i == j){
		g_uu.set(i,j).set_outer_boundary(nz-1 , 1.) ;
	    }   
	}
    g_uu.std_spectral_base() ;

    Sym_tensor gt_uu(mp, CON, mp.get_bvect_spher()) ;
    for(int i=1; i<=3; i++)
	for(int j=1; j<=i; j++){
	    gt_uu.set(i,j) = flat.con()(i,j) + h_uu(i,j) ;
	}


    Vector xsi(mp, CON, mp.get_bvect_spher()) ;
    xsi.set(1) = 0*0.001/(rr) ;
    xsi.set(2) = 0*0.001/(rr) ;
    xsi.set(3) = 0*0.001/(rr) ;
    xsi.std_spectral_base() ;
    h_uu.std_spectral_base() ;

    //----------------------------------------------
    // Vector Poisson equation for xsi
    //----------------------------------------------

    // Source
    //--------


    Vector source(mp, CON, mp.get_bvect_spher()) ;
    Vector source1(mp, CON, mp.get_bvect_spher()) ;
    Vector source2(mp, CON, mp.get_bvect_spher()) ;
    Vector source3(mp, CON, mp.get_bvect_spher()) ;
    Vector source4(mp, CON, mp.get_bvect_spher()) ;
    Vector source5(mp, CON, mp.get_bvect_spher()) ;
    Vector source6(mp, CON, mp.get_bvect_spher()) ;
    Vector source7(mp, CON, mp.get_bvect_spher()) ;
    Vector source8(mp, CON, mp.get_bvect_spher()) ;
    Vector source9(mp, CON, mp.get_bvect_spher()) ;

    const Tensor& dcov_huu = h_uu.derive_cov(flat) ;
    Tensor dcovdcov_huu = dcov_huu.derive_cov(flat) ;
    dcovdcov_huu.inc_dzpuis() ;

    int mer_max = 200 ;
    Vector xsi_jm1(mp, CON, mp.get_bvect_spher()) ;
    for(int i=1; i<=3; i++) xsi_jm1.set(i) = 0. ;
    double diff_xsi ;
    diff_xsi = 1 ;

    for(int mer=0; mer<mer_max; mer++){

	cout << 
   "========================================================================"
	     << endl ;
	cout << "step = " << mer << "       diff_xsi =  "<< diff_xsi << endl ; 
	cout << 
   "========================================================================" 
	     << endl ;

	xsi.annule_domain(0) ;


	const Tensor& dcov_xsi = xsi.derive_cov(flat) ;
	Tensor dcovdcov_xsi = dcov_xsi.derive_cov(flat) ;
	dcovdcov_xsi.inc_dzpuis() ;

	source1 = - contract(dcov_huu, 1, 2) ;
	source1.inc_dzpuis(2) ;
    
	source2 = contract(contract(dcov_xsi, 0, dcov_huu, 2), 0, 2) ;
    
	source3 = contract(xsi * contract(dcovdcov_huu, 1, 3), 0, 2) ;

	source4 = - contract(contract(dcov_huu, 1, 2)*dcov_xsi, 0, 2) ;
    
	source5 = - contract(h_uu, 0, 1, dcovdcov_xsi, 1, 2) ;
	source5.annule_domain(nz-1) ;

	source6 = - contract(contract(dcov_huu, 1, dcov_xsi, 1), 1, 2) ;
    
	source7 = - contract(h_uu *contract(dcovdcov_xsi, 0, 2), 1, 2) ;
	source7.annule_domain(nz-1) ;

	source8 = - 0.333333333333333 * xsi.divergence(flat)
	    .derive_con(flat) ;
	source8.inc_dzpuis() ;
/*	
	des_meridian(source8_xsi(1), 0., 4., "source8_xsi", 10) ; 
	des_meridian(source8_xsi(2), 0., 4., "source8_xsi", 10) ; 
	arrete() ;
*/

	source9 = 0.666666666666666 * contract(dcov_huu, 1, 2) *
	    xsi.divergence(flat) ;

	for (int i=1; i<=3; i++) {
	    source.set(i) = source1(i) + source2(i) 
		+ source3(i) + source4(i) + source5(i)  
		+ source6(i) + source7(i) + source8(i) + source9(i) ;
	} 

	source.annule_domain(0) ;

	
	// Printing
	//-----------

	cout << "moyenne de la source 1 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source1(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 2 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source2(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 3 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source3(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 4 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source4(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 5 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source5(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 6 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source6(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 7 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source7(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 8 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source8(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 9 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source9(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}

    
 
	// Resolution of the Poisson equation 
	// ----------------------------------

	double lambda = 0. ;
	double precision = 1.e-10 ;
	int num_front = 0 ; // index of the intern boundary
	int itermax = 20 ;

	Tenseur source_vect (mp, 1, CON, mp.get_bvect_cart()) ;
	source.change_triad(mp.get_bvect_cart()) ;
	Cmp source_1 (source(1)) ;
	Cmp source_2 (source(2)) ;
	Cmp source_3 (source(3)) ;
	source_vect.set_etat_qcq() ;
	source_vect.set(0) = source_1 ;
	source_vect.set(1) = source_2 ;
	source_vect.set(2) = source_3 ;
	source.change_triad(mp.get_bvect_spher()) ;

	Tenseur xsi_vect (mp, 1, CON, mp.get_bvect_cart()) ;
	xsi.change_triad(mp.get_bvect_cart()) ;
	Cmp xsi_vect1 (xsi(1)) ;
	Cmp xsi_vect2 (xsi(2)) ;
	Cmp xsi_vect3 (xsi(3)) ;
	xsi_vect.set_etat_qcq() ;
	xsi_vect.set(0) = xsi_vect1 ;
	xsi_vect.set(1) = xsi_vect2 ;
	xsi_vect.set(2) = xsi_vect3 ;
	xsi.change_triad(mp.get_bvect_spher()) ;

	Valeur lim_x (mg) ;
	Valeur lim_y (mg) ;
	Valeur lim_z (mg) ;
	lim_x = aa*aa ;
	lim_y = aa*aa ;
	lim_z = aa*aa ;
	lim_x.std_base_scal() ;
	lim_y.std_base_scal() ;
	lim_z.std_base_scal() ;
	
	cout << "Resolution of equation for xsi : poisson_vect_frontiere "  
	     << endl << "----------------------------------------------------"
	     << endl ;


	if (mer>=200) {
	    des_profile(source_vect(1), 0, 20, 0, 0) ;
	    des_coef_xi(source_vect(1).va, 1, 0, 0) ;
	    des_coef_xi(source_vect(1).va, 2, 0, 0) ;
	}

	poisson_vect_frontiere(lambda, source_vect, xsi_vect, lim_x,
			       lim_y, lim_z, num_front, precision, itermax) ;

	
	if (mer>=200) {
	    des_profile(xsi_vect(1), 0, 20, 0, 0) ;
	    des_coef_xi(xsi_vect(1).va, 1, 0, 0) ;
	    des_coef_xi(xsi_vect(1).va, 2, 0, 0) ;
	}

	source_vect.change_triad(mp.get_bvect_spher()) ;
	xsi_vect.change_triad(mp.get_bvect_spher()) ;
	
	xsi.set(1) = xsi_vect(0) ;
	xsi.set(2) = xsi_vect(1) ;
	xsi.set(3) = xsi_vect(2) ;

	cout << "xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(xsi(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << endl ;

	// Check: has the equation for xsi been correctly solved ?
	// --------------------------------------------------------------

	Vector lap_xsi = (xsi.derive_con(flat)).divergence(flat) 
	    + lambda * xsi.divergence(flat).derive_con(flat) ;
	lap_xsi.inc_dzpuis() ;

	cout << "moyenne de la source pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << endl ;

  	cout << "lap_xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(lap_xsi(i)/(nr[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << endl ;

	Tbl tdiff_xsi_r = diffrel(lap_xsi(1), source(1)) ; 
	Tbl tdiff_xsi_t = diffrel(lap_xsi(2), source(2)) ; 
	Tbl tdiff_xsi_p = diffrel(lap_xsi(3), source(3)) ; 

	cout << 
	    "Relative error in the resolution of the equation for xsi : "
	     << endl ; 
	cout << "r component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_xsi_r(l) << "  " ; 
	}
	cout << endl ;
	cout << "t component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_xsi_t(l) << "  " ; 
      }
	cout << endl ;
	cout << "p component : " ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_xsi_p(l) << "  " ; 
	}
	cout << endl ;
	

	int cont ;
	cont = 0 ;
	
	for(int nn=2; nn<=3; nn++)
	    for(int l=1; l<=nz-1; l++)
		for(int k=0; k<=np-1; k++)
		    for(int j=0; j<=nt-1; j++)
			for(int i=0; i<=nr[l]-1; i++){
			    if (fabs(xsi(nn).val_grid_point(l,k,j,i)) > 1.e-10){
			    if (fabs(xsi(nn).val_grid_point(l,k,j,i)) 
				> fabs(xsi_jm1(nn).val_grid_point(l,k,j,i))) {
			
				diff_xsi+=fabs((xsi(nn).val_grid_point(l,k,j,i)
				     -xsi_jm1(nn).val_grid_point(l,k,j,i)) 
				     / (xsi(nn).val_grid_point(l,k,j,i))) ;
			    }
			    else {
			    diff_xsi+=fabs((xsi(nn).val_grid_point(l,k,j,i)
				   -xsi_jm1(nn).val_grid_point(l,k,j,i)) 
				   / (xsi_jm1(nn).val_grid_point(l,k,j,i))) ;
			    }	   
			    cont++ ;
			    }
			}

	diff_xsi = diff_xsi / cont ;

	if(diff_xsi < seuil) {
	    cout << 
   "========================================================================"
	     << endl ;
	cout << "step = " << mer << "       diff_xsi =  "<< diff_xsi << endl ; 
	cout << 
   "========================================================================" 
	     << endl ; 
	    mer = mer_max ;
	}

	for (int i=1; i<=3; i++) {
	    xsi_jm1.set(i) = xsi(i) ;
	}

    }


    // Compute h_uu in dirac gauge :
    //-----------------------------

    Sym_tensor guu_dirac (mp, CON, mp.get_bvect_spher()) ;

    guu_dirac = g_uu.derive_lie(xsi) ;
    guu_dirac.dec_dzpuis(2) ;
    guu_dirac = guu_dirac + g_uu ;


    cout << "Relative difference between g_uu in dirac gauge and in isotropic gauge "
	 << endl ;
	cout << " Comp 1 1 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(guu_dirac(1,1), g_uu(1,1))(l) << "  " ; 
	}
	cout << endl ;
	cout << " Comp 2 2 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(guu_dirac(2,2), g_uu(2,2))(l) << "  " ; 
	}
	cout << endl ;
	cout << " Comp 3 3 :  " ; 
	for (int l=0; l<nz; l++) {
	    cout << diffrel(guu_dirac(3,3), g_uu(3,3))(l) << "  " ; 
	}
	cout << endl ;


       
 	cout << "norme de g_uu en jauge de dirac :" << endl ;
	for (int i=1; i<=3; i++)
	    for (int j=1; j<=i; j++) {
		cout << "  Comp. " << i << " " << j << " :  " ;
		for (int l=0; l<nz; l++){
		    cout << norme(guu_dirac(i,j)/(nr[0]*nt*np))(l) << " " ;
		}
		cout << endl ;
	    }
	cout << endl ;
	
	// Check of the Dirac gauge
	// ------------------------

	for (int i=1; i<=3; i++)
	    for (int j=1; j<=i; j++){
		guu_dirac.set(i,j) = guu_dirac(i,j) - 1. ;
		guu_dirac.set(i,j).annule_domain(0) ;
		guu_dirac.set(i,j) = guu_dirac(i,j) + 1. ;
	    }

	Metric g_dirac (guu_dirac) ;
	Sym_tensor gtuu_dirac (mp, CON, mp.get_bvect_spher()) ;
	gtuu_dirac = pow(g_dirac.determinant(), 1./3.) * guu_dirac ;
	gtuu_dirac.std_spectral_base() ;

	Vector& d_gtuu_dirac (gtuu_dirac.divergence(flat)) ;

	cout << "Is Dirac gauge really satisfied ??!" << endl ;
	cout << "Vector H^i" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " : " << norme(d_gtuu_dirac(i)
					     /(nr[0]*nt*np)) << endl ;
	}	

	cout << "For comparaison value norme(h^11_dirac)/dist = " << endl 
	     << norme(gtuu_dirac(1,1)-1)/(nr[0]*nt*np) * 2 / hh ; 
	

	    
    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ;
    delete [] type_r ; 
    delete [] bornes ; 
    
    return EXIT_SUCCESS ; 
        
}  
  
    
