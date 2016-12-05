/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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
 * Revision 1.13  2016/12/05 16:17:45  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.12  2014/10/13 08:52:40  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.11  2014/10/06 15:12:58  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.10  2007/04/24 20:14:04  f_limousin
 * Implementation of Dirichlet and Neumann BC for the lapse
 *
 * Revision 1.9  2007/02/03 07:46:30  p_grandclement
 * Addition of term kss for psi BC
 *
 * Revision 1.8  2006/04/27 09:12:31  p_grandclement
 * First try at irrotational black holes
 *
 * Revision 1.7  2006/04/25 07:21:57  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.6  2005/08/29 15:10:13  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.5  2004/03/25 10:28:57  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2003/11/25 07:11:09  k_taniguchi
 * Change some arguments from the class Etolie_bin to Et_bin_nsbh.
 *
 * Revision 1.3  2003/11/13 13:43:53  p_grandclement
 * Addition of things needed for Bhole::update_metric (const Etoile_bin&, double, double)
 *
 * Revision 1.2  2003/10/24 13:05:49  p_grandclement
 * correction of the equations for Bin_ns_bh...
 *
 * Revision 1.1  2003/02/13 16:40:25  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>
#include <cmath>

// Lorene
#include "tenseur.h"
#include "bhole.h"
#include "proto.h"
#include "utilitaires.h"
#include "et_bin_nsbh.h"
#include "graphique.h"
#include "scalar.h"

//Resolution pour le lapse pour 1 seul trou
namespace Lorene {
void Bhole::solve_lapse_with_ns (double relax, int bound_nn, double lim_nn) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "Resolution LAPSE" << endl ;
        
    // Pour la relaxation ...
    Cmp lapse_old (n_auto()) ;
    Tenseur auxi (flat_scalar_prod(tkij_tot, tkij_auto)) ;
    Tenseur kk (mp) ;
    kk = 0 ;
    Tenseur work(mp) ;
    work.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++) {
	work.set() = auxi(i, i) ;
	kk = kk + work ;
	}
  
    // La source
    Cmp psiq (pow(psi_tot(), 4.)) ;
    psiq.std_base_scal() ;
    Cmp source 
    (-2*flat_scalar_prod(psi_auto.gradient(), grad_n_tot)()/psi_tot()
	+psiq*n_tot()*kk()) ;
    source.std_base_scal() ;    
      
    Cmp soluce(n_auto()) ;
    
    if (bound_nn == 0){
      // Dirichlet
      Valeur limite (mp.get_mg()->get_angu()) ;
      limite = -0.5 + lim_nn ;
      int np = mp.get_mg()->get_np(1) ;
      int nt = mp.get_mg()->get_nt(1) ;
      for (int k=0 ; k<np ; k++) 
	for (int j=0 ; j<nt ; j++)
	limite.set(0,k,j,0) -= n_comp() (1, k, j, 0) ;
      limite.std_base_scal() ;

      soluce = source.poisson_dirichlet(limite, 0) ;
    }
    else {
      assert(bound_nn == 1);
      // Neumann
      Valeur limite (mp.get_mg()->get_angu()) ;
      limite.annule_hard() ;
      int np = mp.get_mg()->get_np(1) ;
      int nt = mp.get_mg()->get_nt(1) ;
      for (int k=0 ; k<np ; k++) 
	for (int j=0 ; j<nt ; j++)
	limite.set(0,k,j,0) -= n_tot()(1, k, j, 0)/psi_tot()(1,k,j,0)*
	  psi_tot().dsdr()(1,k,j,0) ;
      limite.std_base_scal() ;
      
      soluce = source.poisson_neumann(limite, 0) ;
    }

    soluce = soluce + 0.5 ;

    n_auto.set() = relax*soluce + (1-relax)*lapse_old ; 
    n_auto.set().raccord(3) ;
}

// Resolution sur Psi :
void Bhole::solve_psi_with_ns (double relax) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "Resolution PSI" << endl ;
    
    Cmp psi_old (psi_auto()) ;
    Tenseur auxi (flat_scalar_prod(tkij_auto, tkij_tot)) ;
    Tenseur kk (mp) ;
    kk = 0 ;
    Tenseur work(mp) ;
    work.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++) {
	work.set() = auxi(i, i) ;
	kk = kk + work ;
	}
    Cmp psic (pow(psi_tot(), 5.)) ;
    psic.std_base_scal() ;

    // La source :
    Cmp source (-psic*kk()/8.) ;
    source.std_base_scal() ;
    
    // Condition limite :
    Valeur limite (mp.get_mg()->get_angu()) ;
    limite = 1 ;


    int np = mp.get_mg()->get_np(1) ;
    int nt = mp.get_mg()->get_nt(1) ;

    double* vec_s = new double[3] ;
    Mtbl tet_mtbl (mp.get_mg()) ;
    tet_mtbl = mp.tet ;
    Mtbl phi_mtbl (mp.get_mg()) ;
    phi_mtbl = mp.phi ;

    for (int k=0 ; k<np ; k++) 
      for (int j=0 ; j<nt ; j++) {
	
	double tet = tet_mtbl(1,k,j,0) ;
        double phi = phi_mtbl(1,k,j,0) ;
        vec_s[0] = cos(phi)*sin(tet) ;
	vec_s[1] = sin(phi)*sin(tet) ;
	vec_s[2] = cos(tet) ;
	double part_ss = 0 ;
	if (tkij_tot.get_etat()==ETATQCQ) 
	for (int m=0 ; m<3 ; m++)
		for (int n=0 ; n<3 ; n++)
			part_ss += vec_s[m]*vec_s[n]*tkij_tot(m,n)(1,k,j,0) ;
	part_ss *= pow(psi_tot()(1,k,j,0),3.)/4. ;


	limite.set(0, k, j, 0) = -0.5/rayon*psi_tot()(1, k, j, 0) -
	  psi_comp().dsdr()(1, k, j, 0) - part_ss ;
}


    limite.std_base_scal() ;
    
    Cmp soluce (source.poisson_neumann(limite, 0)) ;
    soluce = soluce + 1./2. ;
    
    psi_auto.set() = relax*soluce + (1-relax)*psi_old ;    
    psi_auto.set().raccord(3) ;

}

// Le shift. Processus iteratif pour cause de CL.
void Bhole::solve_shift_with_ns (const Et_bin_nsbh& ns, 
				 double precision, double relax,
				 int bound_nn, double lim_nn) {
    
    assert (precision > 0) ;
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "resolution SHIFT" << endl ;
    
    Tenseur shift_old (shift_auto) ;
    
    Tenseur source (-6*flat_scalar_prod(taij_tot, psi_auto.gradient())/psi_tot
		    + 2*flat_scalar_prod(tkij_tot, n_auto.gradient())) ;
    source.set_std_base() ;
    
    // On verifie si les 3 composantes ne sont pas nulles :
    if (source.get_etat() == ETATQCQ) {
	int indic = 0 ;
	for (int i=0 ; i<3 ; i++)
	    if (source(i).get_etat() == ETATQCQ)
		indic = 1 ;
	if (indic ==0)
	  for (int i=0 ; i<3 ; i++)
	    source.set_etat_zero() ;
    }
 

    // On filtre les hautes frequences pour raison de stabilite :
    if (source.get_etat() == ETATQCQ)
	for (int i=0 ; i<3 ; i++)
    	    source.set(i).filtre(4) ;
    
    
    // On determine les conditions limites en fonction de omega et de NS :
    int np = mp.get_mg()->get_np(1) ;
    int nt = mp.get_mg()->get_nt(1) ;
    
    Mtbl x_mtbl (mp.get_mg()) ;
    x_mtbl.set_etat_qcq() ;
    Mtbl y_mtbl (mp.get_mg()) ;
    y_mtbl.set_etat_qcq() ;
    x_mtbl = mp.x ;
    y_mtbl = mp.y ;

    double air, theta, phi, xabs, yabs, zabs ;
    Mtbl Xabs (mp.get_mg()) ;
    Xabs = mp.xa ;
    Mtbl Yabs (mp.get_mg()) ;
    Yabs = mp.ya ;
    Mtbl Zabs (mp.get_mg()) ;
    Zabs = mp.za ;
    
    Mtbl tet_mtbl (mp.get_mg()) ;
    tet_mtbl = mp.tet ;
    Mtbl phi_mtbl (mp.get_mg()) ;
    phi_mtbl = mp.phi ;

    // Les bases pour les conditions limites :
    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;
    
    Valeur lim_x (mp.get_mg()->get_angu()) ;
    lim_x = 1 ;
    Valeur lim_y (mp.get_mg()->get_angu()) ;
    lim_y = 1 ;
    Valeur lim_z (mp.get_mg()->get_angu()) ;
    lim_z = 1 ;

    for (int k=0 ; k<np ; k++)
      for (int j=0 ; j<nt ; j++) {

	double tet = tet_mtbl(1,k,j,0) ;
        double phy = phi_mtbl(1,k,j,0) ;

	xabs = Xabs (1, k, j, 0) ;
	yabs = Yabs (1, k, j, 0) ;
	zabs = Zabs (1, k, j, 0) ;
	
	ns.get_mp().convert_absolute (xabs, yabs, zabs, air, theta, phi) ;
	
	lim_x.set(0, k, j, 0) = omega*Yabs(0, 0, 0, 0) + 
	  omega_local*y_mtbl(1,k,j,0) - 
	  ns.get_shift_auto()(0).val_point(air, theta, phi) +
	  n_tot()(1,k,j,0)/psi_tot()(1,k,j,0)/psi_tot()(1,k,j,0)*
	  cos(phy)*sin(tet) ;
	lim_x.base = *bases[0] ;
	
	
	lim_y.set(0, k, j, 0) = -omega*Xabs(0, 0, 0, 0) - 
	  omega_local*x_mtbl(1,k,j,0) - 
	  ns.get_shift_auto()(1).val_point(air, theta, phi) +
	  n_tot()(1,k,j,0)/psi_tot()(1,k,j,0)/psi_tot()(1,k,j,0)*
	  sin(phy)*sin(tet) ;
	
	lim_z.set(0, k, j, 0) = - 
	  ns.get_shift_auto()(2).val_point(air, theta, phi) +
	  n_tot()(1,k,j,0)/psi_tot()(1,k,j,0)/psi_tot()(1,k,j,0)*cos(tet) ;
      }

    lim_x.base = *bases[0] ;
    lim_y.base = *bases[1] ;
    lim_z.base = *bases[2] ;
    
    // On n'en a plus besoin
    for (int i=0 ; i<3 ; i++)
      delete bases[i] ;
    delete [] bases ;
    
    // On resout :
    poisson_vect_frontiere(1./3., source, shift_auto, lim_x, lim_y, 
			   lim_z, 0, precision, 20) ;
   
    shift_auto = relax*shift_auto + (1-relax)*shift_old ;

    for (int i=0; i<3; i++)
      shift_auto.set(i).raccord(3) ;


    // Regularisation of the shift if necessary
    // -----------------------------------------    
    if (bound_nn == 0 && lim_nn == 0)
      regul = regle (shift_auto, ns.get_shift_auto(), omega, omega_local) ;  
    else 
      regul = 0. ;
    
}


void Bhole::equilibrium (const Et_bin_nsbh& comp, 
			 double precision, double relax,
			 int bound_nn, double lim_nn) {

  // Solve for the lapse :
  solve_lapse_with_ns (relax, bound_nn, lim_nn) ;

  // Solve for the conformal factor :
  solve_psi_with_ns (relax) ;

  if (omega != 0) 
  // Solve for the shift vector :
  	solve_shift_with_ns (comp, precision, relax, bound_nn, lim_nn) ;
  
}

  
void Bhole::update_metric (const Et_bin_nsbh& comp) {

	fait_n_comp(comp) ;
	fait_psi_comp(comp) ;
	/*	
	Scalar lapse_auto (n_auto()) ;
	Scalar lapse_tot (n_tot()) ;
	Scalar lapse_comp (n_comp()) ;
	des_meridian(lapse_auto, 0, 7, "n_auto", 0) ;
	des_meridian(lapse_comp, 0, 7, "n_comp", 11) ;
	des_meridian(lapse_tot, 0, 7, "n_tot", 1) ;

	Scalar psiauto (psi_auto()) ;
	Scalar psitot (psi_tot()) ;
	des_meridian(psiauto, 0, 7, "psi_auto", 2) ;
	des_meridian(psitot, 0, 7, "psi_tot", 3) ;
	*/

}

  
}
