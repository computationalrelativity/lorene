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
#include <stdlib.h>
#include <math.h>

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
	{ bornes[l] = bornes[l] * hh * 0.5 ;
	}
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
    h_uu.std_spectral_base() ;

    Metric_flat flat(mp.flat_met_spher()) ;
 
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
    g_uu.set(1,2) = 1e-16 ;
    g_uu.set(1,3) = 1e-16 ;
    g_uu.set(2,3) = 1e-16 ;
    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++)
	  g_uu.set(i,j).set_domain(0) = 1. ;
    g_uu.std_spectral_base() ;

    Metric gamma (g_uu) ;

    Sym_tensor gt_uu(mp, CON, mp.get_bvect_spher()) ;
    for(int i=1; i<=3; i++)
	for(int j=1; j<=i; j++){
	    gt_uu.set(i,j) = flat.con()(i,j) + h_uu(i,j) ;
	}
    gt_uu.std_spectral_base() ;


    Vector xsi(mp, CON, mp.get_bvect_spher()) ;
    xsi.set(1) = 0*0.1/(rr) ;
    xsi.set(2) = 0*0.001/(rr) ;
    xsi.set(3) = 0*0.001/(rr) ;
    xsi.std_spectral_base() ;

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
  
    const Tensor& dcov_huu = h_uu.derive_cov(flat) ;
    Tensor dcovdcov_huu = dcov_huu.derive_cov(flat) ;
    dcovdcov_huu.inc_dzpuis() ;

    Vector xsi_jm1(mp, CON, mp.get_bvect_spher()) ;
    for(int i=1; i<=3; i++) xsi_jm1.set(i) = 0. ;
    double diff_xsi ;
    diff_xsi = 1 ;

    for(int mer=0; mer<niter; mer++){

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

	source1 = contract(dcov_huu, 1, 2) ;
	source1.inc_dzpuis(2) ;
    
	source2 = contract(contract(dcov_xsi, 0, dcov_huu, 2), 0, 2) ;
    
	source3 = contract(xsi * contract(dcovdcov_huu, 1, 3), 0, 2) ;

	source4 = - contract(contract(dcov_huu, 1, 2)*dcov_xsi, 0, 2) ;
    
	source5 = - contract(h_uu, 0, 1, dcovdcov_xsi, 1, 2) ;

	source6 = - contract(contract(dcov_huu, 1, dcov_xsi, 1), 1, 2) ;
    
	source7 = - contract(h_uu *contract(dcovdcov_xsi, 0, 1), 1, 2) / 3. ;

	source8 = 0.666666666666666 * contract(dcov_huu, 1, 2) *
	    xsi.divergence(flat) ;

	for (int i=1; i<=3; i++) {
	    source.set(i) = source1(i) + 0*source2(i) 
		+ 0*source3(i) + 0*source4(i) + 0*source5(i)  
		+ 0*source6(i) + 0*source7(i) + 0*source8(i) ;
	} 

	source.annule_domain(0) ;
	for (int i=1; i<=3; i++) 
	  source.set(i).set_domain(nz-1) = 1.e-15 ;
	
	
	// Printing
	//-----------

	cout << "moyenne de la source 1 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source1(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	/*
	cout << "moyenne de la source 2 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source2(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 3 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source3(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 4 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source4(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 5 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source5(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 6 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source6(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 7 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source7(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source 8 pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source8(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << "moyenne de la source pour xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(source(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	*/
   
	// Resolution of the Poisson equation 
	// ----------------------------------

	double lambda = 0. ;
	Vector source_reg = - (1./3. - lambda) * xsi.divergence(flat)
	    .derive_con(flat) ;
	source_reg.inc_dzpuis() ;
	source += source_reg ;

	double precision = 1.e-8 ;
	int num_front = 0 ; // index of the intern boundary
	int itermax = 20 ;

	Vector limite (mp, CON, mp.get_bvect_cart()) ;
	limite.std_spectral_base() ;
	Valeur lim_x (mg) ;
	Valeur lim_y (mg) ;
	Valeur lim_z (mg) ;
	lim_x = 0. ;
	lim_y = 0. ;
	lim_z = 0. ;
	lim_x.set_base(limite(1).get_spectral_va().get_base()) ;
	lim_y.set_base(limite(2).get_spectral_va().get_base()) ;
	lim_z.set_base(limite(3).get_spectral_va().get_base()) ;


	
	cout << "Resolution of equation for xsi : poisson_vect_frontiere "  
	     << endl << "----------------------------------------------------"
	     << endl ;
	
	poisson_vect_boundary(lambda, source, xsi, lim_x,
			       lim_y, lim_z, num_front, precision, itermax) ;

	
	if (mer>=200) {
	    des_profile(xsi(1), 0, 20, 0, 0) ;
	}

	// Test
	source.dec_dzpuis() ;
	maxabs(xsi.derive_con(flat).divergence(flat) 
	       + lambda * xsi.divergence(flat)
	       .derive_con(flat) - source ,
	       "Absolute error in the resolution of the equation for beta") ;  
	cout << endl ;

	cout << "xsi :" << endl ;
	for (int i=1; i<=3; i++){
	    cout << "  Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << norme(xsi(i)/(nr_tab[0]*nt*np))(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << endl ;

	// End of the computation ?
	int cont ;
	cont = 0 ;
	
	for(int nn=2; nn<=3; nn++)
	    for(int l=1; l<=nz-1; l++)
		for(int k=0; k<=np-1; k++)
		    for(int j=0; j<=nt-1; j++)
			for(int i=0; i<=nr_tab[l]-1; i++){
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
	    mer = niter ;
	}
	
	xsi = relax * xsi + (1-relax) * xsi_jm1 ;
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
	Sym_tensor gtuu_dirac (mp, CON, mp.get_bvect_spher()) ;
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

	Vector hh_dirac (contract(dcov_huu, 1, 2)) ;
	cout << "For comparaison H^i before computation = " << endl 
	     << norme(hh_dirac(1))/(nr_tab[0]*nt*np) 
	     << endl 
	     << norme(hh_dirac(2))/(nr_tab[0]*nt*np) 
	     << endl 
	     << norme(hh_dirac(3))/(nr_tab[0]*nt*np) 
	     << endl ; 
	
	// Another calculation
		
	Sym_tensor gtuu_dirac2 (mp, CON, mp.get_bvect_spher()) ;
	Scalar psi4 = pow(a2, 2./3.) * pow(b2, 1./3.) ;
	psi4.std_spectral_base() ;
	gtuu_dirac2 = g_uu.derive_lie(xsi) * psi4 + 0.666666666666666 * 
	  xsi.divergence(gamma) * gt_uu ;
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
  
    
