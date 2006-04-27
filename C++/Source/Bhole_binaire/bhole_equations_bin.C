/*
 *   Copyright (c) 2001 Philippe Grandclement
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


char bhole_equations_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2006/04/27 09:12:32  p_grandclement
 * First try at irrotational black holes
 *
 * Revision 1.2  2002/10/16 14:36:33  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2001/05/07  09:11:56  phil
 * *** empty log message ***
 *
 * Revision 2.5  2001/04/30  09:30:50  phil
 * filtre pour poisson vectoriel
 *
 * Revision 2.4  2001/04/26  12:04:34  phil
 * *** empty log message ***
 *
 * Revision 2.3  2001/04/25  15:08:48  phil
 * vire fait_tkij
 *
 * Revision 2.2  2001/04/02  12:16:11  phil
 * *** empty log message ***
 *
 * Revision 2.1  2001/03/22  10:41:36  phil
 * modification prototypage
 *
 * Revision 2.0  2001/03/01  08:18:04  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tenseur.h"
#include "bhole.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

// Resolution pour le lapse
void Bhole_binaire::solve_lapse (double precision, double relax) {
    
    assert ((relax >0) && (relax<=1)) ;
    
    cout << "-----------------------------------------------" << endl ;
    cout << "Resolution LAPSE" << endl ;
    
    Tenseur lapse_un_old (hole1.n_auto) ;
    Tenseur lapse_deux_old (hole2.n_auto) ;
    
    Tenseur auxi_un (flat_scalar_prod(hole1.tkij_tot, hole1.tkij_auto)) ;
    Cmp kk_un (auxi_un(0, 0)+auxi_un(1, 1)+auxi_un(2, 2)) ;
    
    Tenseur auxi_deux (flat_scalar_prod(hole2.tkij_tot, hole2.tkij_auto)) ;
    Cmp kk_deux (auxi_deux(0, 0)+auxi_deux(1, 1)+auxi_deux(2, 2)) ;
    
    // Les sources
    
    Cmp source_un 
(-2*flat_scalar_prod(hole1.grad_n_tot, hole1.psi_auto.gradient())()/hole1.psi_tot()
	+hole1.n_tot()*pow(hole1.psi_tot(), 4.)*kk_un) ;
    source_un.std_base_scal() ;

    Cmp source_deux   
(-2*flat_scalar_prod(hole2.grad_n_tot, hole2.psi_auto.gradient())()/hole2.psi_tot()
	+hole2.n_tot()*pow(hole2.psi_tot(), 4.)*kk_deux) ;
    source_deux.std_base_scal() ;
  	   
    //On resout
    hole1.n_auto.set() = hole1.n_auto() - 1./2. ;
    hole2.n_auto.set() = hole2.n_auto() - 1./2. ; 
    
    dirichlet_binaire (source_un, source_deux, -1., -1., hole1.n_auto.set(),
						    hole2.n_auto.set(), 0, precision) ;
    
    hole1.n_auto.set() = hole1.n_auto() + 1./2. ;
    hole2.n_auto.set() = hole2.n_auto() + 1./2. ;
    
    hole1.n_auto.set().raccord(1) ;
    hole2.n_auto.set().raccord(1) ;
    
    // La relaxation :
    hole1.n_auto.set() = relax*hole1.n_auto() + (1-relax)*lapse_un_old() ;
    hole2.n_auto.set() = relax*hole2.n_auto() + (1-relax)*lapse_deux_old() ;
 
    hole1.fait_n_comp (hole2) ;
    hole2.fait_n_comp (hole1) ;
}

//Resolution sur Psi
void Bhole_binaire::solve_psi (double precision, double relax) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "-----------------------------------------------" << endl ;
    cout << "Resolution PSI" << endl ;
    
    Tenseur psi_un_old (hole1.psi_auto) ;
    Tenseur psi_deux_old (hole2.psi_auto) ;
    
    Tenseur auxi_un (flat_scalar_prod(hole1.tkij_tot, hole1.tkij_auto)) ;
    Cmp kk_un (auxi_un(0, 0)+auxi_un(1, 1)+auxi_un(2, 2)) ;
    
    Tenseur auxi_deux (flat_scalar_prod(hole2.tkij_tot, hole2.tkij_auto)) ;
    Cmp kk_deux (auxi_deux(0, 0)+auxi_deux(1, 1)+auxi_deux(2, 2)) ;
    
    // Les sources
    Cmp source_un (-pow(hole1.psi_tot(), 5.)*kk_un/8.) ;
    source_un.std_base_scal() ;
    
    Cmp source_deux (-pow(hole2.psi_tot(), 5.)*kk_deux/8.)  ;
    source_deux.std_base_scal() ;
       
    // Les valeurs limites :
    int np_un = hole1.mp.get_mg()->get_np(1) ;
    int nt_un = hole1.mp.get_mg()->get_nt(1) ;
    Valeur lim_un (hole1.mp.get_mg()->get_angu()) ;
    lim_un = 1 ;
    for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_un.set(0, k, j, 0) = -0.5/hole1.rayon*hole1.psi_tot()(1, k, j, 0) ;
    lim_un.std_base_scal() ;
    
    int np_deux = hole2.mp.get_mg()->get_np(1) ;
    int nt_deux = hole2.mp.get_mg()->get_nt(1) ;
    Valeur lim_deux (hole2.mp.get_mg()->get_angu()) ;
    lim_deux = 1 ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	    lim_deux.set(0, k, j, 0) = -0.5/hole2.rayon*hole2.psi_tot()(1, k, j, 0) ;
    lim_deux.std_base_scal() ;
 
    //On resout
    neumann_binaire (source_un, source_deux, lim_un, lim_deux, 
	hole1.psi_auto.set(), hole2.psi_auto.set(), 0, precision) ;
    
    hole1.psi_auto.set() = hole1.psi_auto() + 1./2. ;
    hole2.psi_auto.set() = hole2.psi_auto() + 1./2. ;
    
    hole1.psi_auto.set().raccord(1) ;
    hole2.psi_auto.set().raccord(1) ;
    
    // La relaxation :
    hole1.psi_auto.set() = relax*hole1.psi_auto() + (1-relax)*psi_un_old() ;
    hole2.psi_auto.set() = relax*hole2.psi_auto() + (1-relax)*psi_deux_old() ;
    
    hole1.fait_psi_comp (hole2) ;
    hole2.fait_psi_comp (hole1) ;
}


// Resolution sur shift a omega bloque.
void Bhole_binaire::solve_shift (double precision, double relax) {
    
    cout << "------------------------------------------------" << endl ;
    cout << "Resolution shift : Omega = " << omega << endl ;
    
    // On determine les sources 
    Tenseur source_un (
	-6*flat_scalar_prod (hole1.taij_tot, hole1.psi_auto.gradient())/hole1.psi_tot
	 +2*flat_scalar_prod (hole1.tkij_tot, hole1.n_auto.gradient())) ; 
    if (source_un.get_etat() == ETATZERO) {
	source_un.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    source_un.set(i).set_etat_zero() ;
    }
    source_un.set_std_base() ;
    
    Tenseur source_deux (
	-6*flat_scalar_prod (hole2.taij_tot, hole2.psi_auto.gradient())/hole2.psi_tot
	 +2*flat_scalar_prod (hole2.tkij_tot, hole2.n_auto.gradient())) ;
      if (source_deux.get_etat() == ETATZERO) {
	source_deux.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    source_deux.set(i).set_etat_zero() ;
    }
    source_deux.set_std_base() ;
    
    // On filtre les hautes frequences.
    for (int i=0 ; i<3 ; i++) {
	if (source_un(i).get_etat() != ETATZERO)
	    source_un.set(i).filtre(4) ;
	if (source_deux(i).get_etat() != ETATZERO)
	    source_deux.set(i).filtre(4) ;
    }
    
    // Les alignemenents pour le signe des CL.
    double orientation_un = hole1.mp.get_rot_phi() ;
    assert ((orientation_un==0) || (orientation_un == M_PI)) ;
    
    double orientation_deux = hole2.mp.get_rot_phi() ;
    assert ((orientation_deux==0) || (orientation_deux == M_PI)) ;
    
    int aligne_un = (orientation_un == 0) ? 1 : -1 ;
    int aligne_deux = (orientation_deux == 0) ? 1 : -1 ;
    
    // On regarde si toutes les composantes sont nulles :
    int ind1 = 0 ;
    int ind2 = 0 ;
    for (int i=0 ; i<3 ; i++) {
	if (source_un(i).get_etat() == ETATQCQ)
	    ind1 = 1 ;
	if (source_deux(i).get_etat() == ETATQCQ)
	    ind2 = 1 ;
    }
    
    if (ind1==0)
	source_un.set_etat_zero() ;
    if (ind2==0)
	source_deux.set_etat_zero() ;
	
    // On determine les Cl en fonction de omega :
    int np_un = hole1.mp.get_mg()->get_np (1) ;
    int nt_un = hole1.mp.get_mg()->get_nt (1) ;
    
    Mtbl xa_mtbl_un (source_un.get_mp()->get_mg()) ;
    xa_mtbl_un.set_etat_qcq() ;
    Mtbl ya_mtbl_un (source_un.get_mp()->get_mg()) ;
    ya_mtbl_un.set_etat_qcq() ;
    
    xa_mtbl_un = source_un.get_mp()->xa ;
    ya_mtbl_un = source_un.get_mp()->ya ;
      
    Mtbl x_mtbl_un (source_un.get_mp()->get_mg()) ;
    x_mtbl_un.set_etat_qcq() ;
    Mtbl y_mtbl_un (source_un.get_mp()->get_mg()) ;
    y_mtbl_un.set_etat_qcq() ;
    
    x_mtbl_un = source_un.get_mp()->x ;
    y_mtbl_un = source_un.get_mp()->y ;
    
    // Les bases
    Base_val** bases_un = hole1.mp.get_mg()->std_base_vect_cart() ;
    Base_val** bases_deux = hole2.mp.get_mg()->std_base_vect_cart() ;
    
    Valeur lim_x_un (*hole1.mp.get_mg()->get_angu()) ;
    lim_x_un = 1 ; // Juste pour affecter dans espace des configs ;
    lim_x_un.set_etat_c_qcq() ;
    for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_x_un.set(0, k, j, 0) = aligne_un*omega*ya_mtbl_un(0, 0, 0, 0) + aligne_un*hole1.omega_local*y_mtbl_un(1,k,j,0) ;
    lim_x_un.base = *bases_un[0] ;
    
    Valeur lim_y_un (*hole1.mp.get_mg()->get_angu()) ;
    lim_y_un = 1 ; // Juste pour affecter dans espace des configs ;
    lim_y_un.set_etat_c_qcq() ;
    for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_y_un.set(0, k, j, 0) = -aligne_un*omega*xa_mtbl_un(0, 0, 0, 0) - aligne_un*hole1.omega_local*x_mtbl_un(1,k,j,0) ;;
    lim_y_un.base = *bases_un[1] ;
    
    Valeur lim_z_un (*hole1.mp.get_mg()->get_angu()) ;
    lim_z_un = 1 ;
     for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_z_un.set(0, k, j, 0) = 0 ;
    lim_z_un.base = *bases_un[2] ;
    
    // On determine les Cl en fonction de omega :
    int np_deux = hole2.mp.get_mg()->get_np (1) ;
    int nt_deux = hole2.mp.get_mg()->get_nt (1) ;
    
    Mtbl xa_mtbl_deux (source_deux.get_mp()->get_mg()) ;
    xa_mtbl_deux.set_etat_qcq() ;
    Mtbl ya_mtbl_deux (source_deux.get_mp()->get_mg()) ;
    ya_mtbl_deux.set_etat_qcq() ;
    
    xa_mtbl_deux = source_deux.get_mp()->xa ;
    ya_mtbl_deux = source_deux.get_mp()->ya ;

    Mtbl x_mtbl_deux (source_deux.get_mp()->get_mg()) ;
    x_mtbl_deux.set_etat_qcq() ;
    Mtbl y_mtbl_deux (source_deux.get_mp()->get_mg()) ;
    y_mtbl_deux.set_etat_qcq() ;
    
    x_mtbl_deux = source_deux.get_mp()->x ;
    y_mtbl_deux = source_deux.get_mp()->y ;
    
    Valeur lim_x_deux (*hole2.mp.get_mg()->get_angu()) ;
    lim_x_deux = 1 ; // Juste pour affecter dans espace des configs ;
    lim_x_deux.set_etat_c_qcq() ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	    lim_x_deux.set(0, k, j, 0) = aligne_deux*omega*ya_mtbl_deux(0, 0, 0, 0) + aligne_deux*hole2.omega_local*y_mtbl_deux(1,k,j,0) ;
    lim_x_deux.base = *bases_deux[0] ;
    
    Valeur lim_y_deux (*hole2.mp.get_mg()->get_angu()) ;
    lim_y_deux = 1 ; // Juste pour affecter dans espace des configs ;
    lim_y_deux.set_etat_c_qcq() ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	   lim_y_deux.set(0, k, j, 0) = -aligne_deux*omega*xa_mtbl_deux(0, 0, 0, 0) - aligne_deux*hole2.omega_local*x_mtbl_deux(1,k,j,0) ;
    lim_y_deux.base = *bases_deux[1] ;
    
    Valeur lim_z_deux (*hole2.mp.get_mg()->get_angu()) ;
    lim_z_deux = 1 ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	    lim_z_deux.set(0, k, j, 0) = 0 ;
    lim_z_deux.base = *bases_deux[2] ;
    
    for (int i=0 ; i<3 ; i++) {
	delete bases_un[i] ;
	delete bases_deux[i] ;
	}
    delete [] bases_un ;
    delete [] bases_deux ;
    
    // On resout le truc :
    Tenseur shift_un_old (hole1.shift_auto) ;
    Tenseur shift_deux_old (hole2.shift_auto) ;
    
    poisson_vect_binaire (1./3., source_un, source_deux, 
	lim_x_un, lim_y_un, lim_z_un, 
	lim_x_deux, lim_y_deux, lim_z_deux, 
	hole1.shift_auto, hole2.shift_auto, 0, precision) ;
    
    for (int i=0 ; i<3 ; i++) {
	hole1.shift_auto.set(i).raccord(1) ;
	hole2.shift_auto.set(i).raccord(1) ;
    }
    
    // On regularise les shift.
    hole1.shift_auto = relax*hole1.shift_auto + 
	(1-relax)*shift_un_old ;
    hole2.shift_auto = relax*hole2.shift_auto +
	(1-relax)*shift_deux_old ;
    
    double diff_un = regle (hole1.shift_auto, hole2.shift_auto, omega, hole1.omega_local) ;
    double diff_deux = regle (hole2.shift_auto, hole1.shift_auto, omega, hole2.omega_local) ;
    hole1.regul = diff_un ;
    hole2.regul = diff_deux ;
    
    cout << "Difference relatives : " << diff_un << " " << diff_deux << endl ;
}
