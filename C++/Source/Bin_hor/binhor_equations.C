/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose-Luis Jaramillo
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


char binhor_equations_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/02/07 10:46:28  f_limousin
 * Many changes !! The sources are written differently to minimize the
 * numerical error, the boundary conditions are changed...
 *
 * Revision 1.2  2004/12/31 15:41:26  f_limousin
 * Correction of an error
 *
 * Revision 1.1  2004/12/29 16:11:34  f_limousin
 * First version
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
#include "tensor.h"
#include "tenseur.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

// Resolution for the lapse
void Bin_hor::solve_lapse (double precision, double relax) {
    
    assert ((relax >0) && (relax<=1)) ;
    
    cout << "-----------------------------------------------" << endl ;
    cout << "Resolution LAPSE" << endl ;
    
    Scalar lapse_un_old (hole1.n_auto()) ;
    Scalar lapse_deux_old (hole2.n_auto()) ;

    Sym_tensor taa_un = hole1.aa().up_down(hole1.met_gamt) ;         
    Scalar aa_quad_un = contract(taa_un, 0, 1, hole1.aa_auto(), 0, 1) ; 

    Sym_tensor taa_deux = hole2.aa().up_down(hole2.met_gamt) ;         
    Scalar aa_quad_deux = contract(taa_deux, 0, 1, hole2.aa_auto(), 0, 1) ; 
       
    // Source 1
 
    Scalar source_un (hole1.mp) ;

    source_un = hole1.psi4()*( hole1.nn()*( aa_quad_un + 0.3333333333333333 * 
					    hole1.trK*hole1.trK*hole1.decouple)
			       - hole1.trK_point*hole1.decouple )     

       -2.*contract(hole1.dnn(), 0, hole1.psi_auto()
		    .derive_cov(hole1.ff), 0)/hole1.psi()
      - contract(hole1.hdirac(), 0, hole1.n_auto().derive_cov(hole1.ff), 0) ; 

    Scalar tmp_un (hole1.mp) ;
    
    tmp_un = hole1.psi4()* contract(hole1.beta_auto(), 0, hole1.trK.
				    derive_cov(hole1.ff), 0) 
	- contract( hole1.hh(), 0, 1, hole1.n_auto().derive_cov(hole1.ff)
		    .derive_cov(hole1.ff), 0, 1 ) ;
        
    tmp_un.inc_dzpuis() ; // dzpuis: 3 -> 4
  	   
    source_un += tmp_un ;
    source_un.std_spectral_base() ;

    // Source 2

    Scalar source_deux (hole2.mp) ;

    source_deux = hole2.psi4()*( hole2.nn()*( aa_quad_deux + 0.3333333333333333
					  * hole2.trK*hole2.trK*hole2.decouple)
			       - hole2.trK_point*hole2.decouple ) 
	-2.*contract(hole2.dnn(), 0, hole2.psi_auto()
		     .derive_cov(hole2.ff), 0)/hole2.psi()

  - contract(hole2.hdirac(), 0, hole2.n_auto().derive_cov(hole2.ff), 0) ; 

    Scalar tmp_deux (hole2.mp) ;
    
    tmp_deux = hole2.psi4()* contract(hole2.beta_auto(), 0, hole2.trK.
				    derive_cov(hole2.ff), 0) 
	- contract( hole2.hh(), 0, 1, hole2.n_auto().derive_cov(hole2.ff)
		    .derive_cov(hole2.ff), 0, 1 ) ;
        
    tmp_deux.inc_dzpuis() ; // dzpuis: 3 -> 4
  	   
    source_deux += tmp_deux ;
    source_deux.std_spectral_base() ;

    cout << "source lapse" << endl << norme(source_un) << endl ;

    // Boundary conditions :

    Valeur lim_un (hole1.boundary_nn_Dir(0.2)) ;
    Valeur lim_deux (hole2.boundary_nn_Dir(0.2)) ;

    // Scalar --> Cmp 
    Scalar n_un_temp (hole1.mp) ;
    Scalar n_deux_temp (hole2.mp) ;
    n_un_temp = hole1.n_auto() - 1./2. ;
    n_deux_temp = hole2.n_auto() - 1./2. ;
    
    Cmp sour_un (hole1.mp) ;
    Cmp sour_deux (hole2.mp) ;
    Cmp n_un_cmp (hole1.mp) ;
    Cmp n_deux_cmp (hole2.mp) ;
    
    sour_un = source_un ;
    sour_deux = source_deux ;
    n_un_cmp = n_un_temp ;
    n_deux_cmp = n_deux_temp ;


    // We solve
    //---------

    dirichlet_binaire (sour_un, sour_deux, lim_un, lim_deux, n_un_cmp,
		       n_deux_cmp, 0, precision) ;
    
    n_un_temp = n_un_cmp ;
    n_deux_temp = n_deux_cmp ;

    n_un_temp = n_un_temp + 1./2. ;
    n_deux_temp = n_deux_temp + 1./2. ;
   
    n_un_temp.raccord(1) ;
    n_deux_temp.raccord(1) ;
    
    // Relaxation :
    n_un_temp = relax*n_un_temp + (1-relax)*lapse_un_old ;
    n_deux_temp = relax*n_deux_temp + (1-relax)*lapse_deux_old ;
 
    cout << "lapse" << endl << norme (n_un_temp) << endl ;

    double ttime = hole1.the_time[hole1.jtime] ;
    hole1.n_auto_evol.update(n_un_temp, hole1.jtime, ttime) ;
    hole2.n_auto_evol.update(n_deux_temp, hole2.jtime, ttime) ;

    hole1.n_comp (hole2) ;
    hole2.n_comp (hole1) ;
}


//Resolution for Psi
void Bin_hor::solve_psi (double precision, double relax) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "-----------------------------------------------" << endl ;
    cout << "Resolution PSI" << endl ;
    
    Scalar psi_un_old (hole1.psi_auto()) ;
    Scalar psi_deux_old (hole2.psi_auto()) ;
    
    Sym_tensor taa_un = hole1.aa().up_down(hole1.met_gamt) ;         
    Scalar aa_quad_un = contract(taa_un, 0, 1, hole1.aa_auto(), 0, 1) ; 

    Sym_tensor taa_deux = hole2.aa().up_down(hole2.met_gamt) ;         
    Scalar aa_quad_deux = contract(taa_deux, 0, 1, hole2.aa_auto(), 0, 1) ; 
    
    // Source 1
       
    Scalar source_un (hole1.mp) ;
    Scalar tmp_un (hole1.mp) ;

    tmp_un = 0.125* hole1.psi_auto() * (hole1.met_gamt).ricci_scal() 
      - contract(hole1.hh(), 0, 1, hole1.psi_auto().derive_cov(hole1.ff)
		 .derive_cov(hole1.ff), 0, 1 ) ;
    tmp_un.inc_dzpuis() ; // dzpuis : 3 -> 4
    
    tmp_un -= contract(hole1.hdirac(), 0, hole1.psi_auto()
		    .derive_cov(hole1.ff), 0) ;  

    source_un = tmp_un - hole1.psi()*hole1.psi4()* ( 0.125* aa_quad_un 
       	   - 8.33333333333333e-2* hole1.trK*hole1.trK*hole1.decouple ) ;
    source_un.std_spectral_base() ;

    // Source 2
       
    Scalar source_deux (hole2.mp) ;
    Scalar tmp_deux (hole2.mp) ;

    tmp_deux = 0.125* hole2.psi_auto() * (hole2.met_gamt).ricci_scal() 
      - contract(hole2.hh(), 0, 1, hole2.psi_auto().derive_cov(hole2.ff)
		 .derive_cov(hole2.ff), 0, 1 ) ;
    tmp_deux.inc_dzpuis() ; // dzpuis : 3 -> 4
    
    tmp_deux -= contract(hole2.hdirac(), 0, hole2.psi_auto()
		    .derive_cov(hole2.ff), 0) ;  

    source_deux = tmp_deux - hole2.psi()*hole2.psi4()* ( 0.125* aa_quad_deux 
       	   - 8.33333333333333e-2* hole2.trK*hole2.trK*hole2.decouple ) ;
    source_deux.std_spectral_base() ;


    cout << "source psi" << endl << norme(source_un) << endl ;

    // Boundary conditions :

    Valeur lim_un (hole1.boundary_psi_Neu_spat()) ;
    Valeur lim_deux (hole2.boundary_psi_Neu_spat()) ;
 
/*
    int np_un = hole1.mp.get_mg()->get_np(1) ;
    int nt_un = hole1.mp.get_mg()->get_nt(1) ;
    Valeur lim_un (hole1.mp.get_mg()->get_angu()) ;
    lim_un = 1 ;
    for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_un.set(0, k, j, 0) = -0.5/hole1.radius*(hole1.psi_auto().val_grid_point(1, k, j, 0)+hole1.psi_comp().val_grid_point(1, k, j, 0)) ;
    lim_un.std_base_scal() ;
    
    int np_deux = hole2.mp.get_mg()->get_np(1) ;
    int nt_deux = hole2.mp.get_mg()->get_nt(1) ;
    Valeur lim_deux (hole2.mp.get_mg()->get_angu()) ;
    lim_deux = 1 ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	    lim_deux.set(0, k, j, 0) = -0.5/hole2.radius*(hole2.psi_auto().val_grid_point(1, k, j, 0)+hole2.psi_comp().val_grid_point(1, k, j, 0)) ;
    lim_deux.std_base_scal() ;
*/


    // Scalar --> Cmp 

    Scalar psi_un_temp (hole1.psi_auto()) ;
    Scalar psi_deux_temp (hole2.psi_auto()) ;
    
    Cmp sour_un (hole1.mp) ;
    Cmp sour_deux (hole2.mp) ;
    Cmp psi_un_cmp (hole1.mp) ;
    Cmp psi_deux_cmp (hole2.mp) ;
    
    sour_un = source_un ;
    sour_deux = source_deux ;
    psi_un_cmp = psi_un_temp ;
    psi_deux_cmp = psi_deux_temp ;

    // We solve
    // --------

    neumann_binaire (sour_un, sour_deux, lim_un, lim_deux, 
	psi_un_cmp, psi_deux_cmp, 0, precision) ;
    
    psi_un_temp = psi_un_cmp ;
    psi_deux_temp = psi_deux_cmp ;
    
    psi_un_temp = psi_un_temp + 1./2. ;
    psi_deux_temp = psi_deux_temp + 1./2. ;
     
    psi_un_temp.raccord(1) ;
    psi_deux_temp.raccord(1) ;
    
    // Relaxation :
    psi_un_temp = relax*psi_un_temp + (1-relax)*psi_un_old ;
    psi_deux_temp = relax*psi_deux_temp + (1-relax)*psi_deux_old ;
    
    cout << "psi" << endl << norme (psi_un_temp) << endl ;
 
    double ttime = hole1.the_time[hole1.jtime] ;
    hole1.psi_auto_evol.update(psi_un_temp, hole1.jtime, ttime) ;
    hole2.psi_auto_evol.update(psi_deux_temp, hole2.jtime, ttime) ;

    hole1.psi_comp (hole2) ;
    hole2.psi_comp (hole1) ;
}


// Resolution for shift with omega fixed.
void Bin_hor::solve_shift (double precision, double relax) {
    
    cout << "------------------------------------------------" << endl ;
    cout << "Resolution shift : Omega = " << omega << endl ;
    
    Sym_tensor taa_un = hole1.aa().up_down(hole1.met_gamt) ;         
    Scalar aa_quad_un = contract(taa_un, 0, 1, hole1.aa_auto(), 0, 1) ; 

    Sym_tensor taa_deux = hole2.aa().up_down(hole2.met_gamt) ;         
    Scalar aa_quad_deux = contract(taa_deux, 0, 1, hole2.aa_auto(), 0, 1) ; 

    // Source 1

    Vector source_un (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    Vector tmp_vect_un (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    
    source_un = 2.* contract(hole1.aa(), 1, 
			     hole1.n_auto().derive_cov(hole1.ff), 0 )
	  -6.*contract(hole1.aa_nn[hole1.jtime], 1, 
		       hole1.psi_auto().derive_cov(hole1.ff)/hole1.psi(), 0) ;
//			    - 6.*hole1.n_auto()*hole1.dpsi()/hole1.psi(), 0) ;
     
    tmp_vect_un = 0.66666666666666666* hole1.trK.derive_con(hole1.met_gamt)
	* hole1.decouple ;
    tmp_vect_un.inc_dzpuis() ;
   
    source_un += 2.* hole1.nn() * ( tmp_vect_un
			- contract(hole1.met_gamt.connect().get_delta(), 1, 2, 
				     hole1.aa()*hole1.decouple, 0, 1) ) ;

    Vector vtmp_un = contract(hole1.hh(), 0, 1, 
                           hole1.beta_auto().derive_cov(hole1.ff)
			   .derive_cov(hole1.ff), 1, 2)
      + 0.3333333333333333*contract(hole1.hh(), 1, hole1.beta_auto()
			       .divergence(hole1.ff).derive_cov(hole1.ff), 0) 
      - hole1.hdirac().derive_lie(hole1.beta_auto()) 
      + hole1.gamt_point.divergence(hole1.ff) ;      
    vtmp_un.inc_dzpuis() ; // dzpuis: 3 -> 4

    source_un -= vtmp_un ; 
        
    source_un += 0.66666666666666666* hole1.beta_auto().divergence(hole1.ff) 
	* hole1.hdirac() ;

    source_un.std_spectral_base() ;

    // Source 2

    Vector source_deux (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
    Vector tmp_vect_deux (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
    
    source_deux = 2.* contract(hole2.aa(), 1, 
			       hole2.n_auto().derive_cov(hole2.ff), 0) 
	-6.*contract(hole2.aa_nn[hole2.jtime], 1, hole2.psi_auto().derive_cov(hole2.ff)/hole2.psi(), 0) ;


//			    - 6.*hole2.n_auto()*hole2.dpsi()/hole2.psi(), 0) ;
     
    tmp_vect_deux = 0.66666666666666666* hole2.trK.derive_con(hole2.met_gamt)
	* hole2.decouple ;
    tmp_vect_deux.inc_dzpuis() ;
   
    source_deux += 2.* hole2.nn() * ( tmp_vect_deux
		       - contract(hole2.met_gamt.connect().get_delta(), 1, 2, 
				     hole2.aa()*hole2.decouple, 0, 1) ) ;

    Vector vtmp_deux = contract(hole2.hh(), 0, 1, 
                           hole2.beta_auto().derive_cov(hole2.ff)
			   .derive_cov(hole2.ff), 1, 2)
      + 0.3333333333333333*contract(hole2.hh(), 1, hole2.beta_auto()
			       .divergence(hole2.ff).derive_cov(hole2.ff), 0) 
      - hole2.hdirac().derive_lie(hole2.beta_auto()) 
      + hole2.gamt_point.divergence(hole2.ff) ;      
    vtmp_deux.inc_dzpuis() ; // dzpuis: 3 -> 4

    source_deux -= vtmp_deux ; 
        
    source_deux += 0.66666666666666666* hole2.beta_auto().divergence(hole2.ff) 
	* hole2.hdirac() ;

    source_deux.std_spectral_base() ;
    
    Vector source_1 (source_un) ;
    Vector source_2 (source_deux) ;
    source_1.change_triad(hole1.mp.get_bvect_cart()) ;
    source_2.change_triad(hole2.mp.get_bvect_cart()) ;

    cout << "source shift_x" << endl << norme(source_1(1)) << endl ;
    cout << "source shift_y" << endl << norme(source_1(2)) << endl ;
    cout << "source shift_z" << endl << norme(source_1(3)) << endl ;

    // Filter for high frequencies.
    for (int i=1 ; i<=3 ; i++) {
	source_un.set(i).filtre(4) ;
	source_deux.set(i).filtre(4) ;
    }

    // Les alignemenents pour le signe des CL.
    double orientation_un = hole1.mp.get_rot_phi() ;
    assert ((orientation_un==0) || (orientation_un == M_PI)) ;
    
    double orientation_deux = hole2.mp.get_rot_phi() ;
    assert ((orientation_deux==0) || (orientation_deux == M_PI)) ;
    
    int aligne_un = (orientation_un == 0) ? 1 : -1 ;
    int aligne_deux = (orientation_deux == 0) ? 1 : -1 ;

    // On determine les Cl en fonction de omega :
    int np_un = hole1.mp.get_mg()->get_np (1) ;
    int nt_un = hole1.mp.get_mg()->get_nt (1) ;
    
    tmp_vect_un = hole1.nn() * hole1.radial_vect_hor() ;
    tmp_vect_un.change_triad(hole1.mp.get_bvect_cart() ) ;


    Mtbl xa_mtbl_un (source_un.get_mp().get_mg()) ;
    xa_mtbl_un.set_etat_qcq() ;
    Mtbl ya_mtbl_un (source_un.get_mp().get_mg()) ;
    ya_mtbl_un.set_etat_qcq() ;
    
    xa_mtbl_un = source_un.get_mp().xa ;
    ya_mtbl_un = source_un.get_mp().ya ;
     
    // Les bases
    Base_val** bases_un = hole1.mp.get_mg()->std_base_vect_cart() ;
    Base_val** bases_deux = hole2.mp.get_mg()->std_base_vect_cart() ;
    
    Valeur lim_x_un (*hole1.mp.get_mg()->get_angu()) ;
    lim_x_un = 1 ; // Juste pour affecter dans espace des configs ;
    lim_x_un.set_etat_c_qcq() ;
    for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_x_un.set(0, k, j, 0) = aligne_un*omega*ya_mtbl_un(1, k, j, 0) 
		+ 0*tmp_vect_un(1).val_grid_point(1, k, j, 0) ;
    lim_x_un.base = *bases_un[0] ;
     
    Valeur lim_y_un (*hole1.mp.get_mg()->get_angu()) ;
    lim_y_un = 1 ; // Juste pour affecter dans espace des configs ;
    lim_y_un.set_etat_c_qcq() ;
    for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_y_un.set(0, k, j, 0) = -aligne_un*omega*xa_mtbl_un(1, k, j, 0) 
		+ 0*tmp_vect_un(2).val_grid_point(1, k, j, 0) ;
    lim_y_un.base = *bases_un[1] ;
    
    Valeur lim_z_un (*hole1.mp.get_mg()->get_angu()) ;
    lim_z_un = 1 ;
     for (int k=0 ; k<np_un ; k++)
	for (int j=0 ; j<nt_un ; j++)
	    lim_z_un.set(0, k, j, 0)=0*tmp_vect_un(3).val_grid_point(1, k, j, 0);
    lim_z_un.base = *bases_un[2] ;
    
    // On determine les Cl en fonction de omega :
    int np_deux = hole2.mp.get_mg()->get_np (1) ;
    int nt_deux = hole2.mp.get_mg()->get_nt (1) ;
    
    tmp_vect_deux = hole2.nn() * hole2.radial_vect_hor() ;
    tmp_vect_deux.change_triad(hole2.mp.get_bvect_cart() ) ;


    Mtbl xa_mtbl_deux (source_deux.get_mp().get_mg()) ;
    xa_mtbl_deux.set_etat_qcq() ;
    Mtbl ya_mtbl_deux (source_deux.get_mp().get_mg()) ;
    ya_mtbl_deux.set_etat_qcq() ;
    
    xa_mtbl_deux = source_deux.get_mp().xa ;
    ya_mtbl_deux = source_deux.get_mp().ya ;


    Valeur lim_x_deux (*hole2.mp.get_mg()->get_angu()) ;
    lim_x_deux = 1 ; // Juste pour affecter dans espace des configs ;
    lim_x_deux.set_etat_c_qcq() ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	    lim_x_deux.set(0, k, j, 0) = aligne_deux*omega*ya_mtbl_deux(1, k, j, 0) + 0*tmp_vect_deux(1).val_grid_point(1, k, j, 0) ;
    lim_x_deux.base = *bases_deux[0] ;
    
    Valeur lim_y_deux (*hole2.mp.get_mg()->get_angu()) ;
    lim_y_deux = 1 ; // Juste pour affecter dans espace des configs ;
    lim_y_deux.set_etat_c_qcq() ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	   lim_y_deux.set(0, k, j, 0) = -aligne_deux*omega*xa_mtbl_deux(1, k, j, 0) + 0*tmp_vect_deux(2).val_grid_point(1, k, j, 0) ;
    lim_y_deux.base = *bases_deux[1] ;
    
    Valeur lim_z_deux (*hole2.mp.get_mg()->get_angu()) ;
    lim_z_deux = 1 ;
    for (int k=0 ; k<np_deux ; k++)
	for (int j=0 ; j<nt_deux ; j++)
	    lim_z_deux.set(0, k, j, 0) = 0*tmp_vect_deux(3).val_grid_point(1, k, j, 0) ;
    lim_z_deux.base = *bases_deux[2] ;
    
    for (int i=0 ; i<3 ; i++) {
	delete bases_un[i] ;
	delete bases_deux[i] ;
	}
    delete [] bases_un ;
    delete [] bases_deux ;
    

    // Vector --> Tenseur :
    Vector beta_un_old (hole1.beta_auto()) ;
    Vector beta_deux_old (hole2.beta_auto()) ;
    
    Tenseur source1 (hole1.mp, 1, CON, hole1.mp.get_bvect_spher()) ;
    source1.set_etat_qcq() ;
    Cmp source1_r (hole1.mp) ;
    Cmp source1_t (hole1.mp) ;
    Cmp source1_p (hole1.mp) ;
    source1_r = source_un(1) ;
    source1_t = source_un(2) ;
    source1_p = source_un(3) ;
    source1.set(0) = source1_r ;
    source1.set(1) = source1_t ;
    source1.set(2) = source1_p ;
    source1.change_triad(hole1.mp.get_bvect_cart()) ;

    Tenseur source2 (hole2.mp, 1, CON, hole2.mp.get_bvect_spher()) ;
    source2.set_etat_qcq() ;
    Cmp source2_r (hole2.mp) ;
    Cmp source2_t (hole2.mp) ;
    Cmp source2_p (hole2.mp) ;
    source2_r = source_deux(1) ;
    source2_t = source_deux(2) ;
    source2_p = source_deux(3) ;
    source2.set(0) = source2_r ;
    source2.set(1) = source2_t ;
    source2.set(2) = source2_p ;
    source2.change_triad(hole2.mp.get_bvect_cart()) ;

    Tenseur beta1 (hole1.mp, 1, CON, hole1.mp.get_bvect_spher()) ;
    beta1.set_etat_qcq() ;
    Cmp beta1_r (hole1.mp) ;
    Cmp beta1_t (hole1.mp) ;
    Cmp beta1_p (hole1.mp) ;
    beta1_r = hole1.beta_auto()(1) ;
    beta1_t = hole1.beta_auto()(2) ;
    beta1_p = hole1.beta_auto()(3) ;
    beta1.set(0) = beta1_r ;
    beta1.set(1) = beta1_t ;
    beta1.set(2) = beta1_p ;
    beta1.change_triad(hole1.mp.get_bvect_cart()) ;

    Tenseur beta2 (hole2.mp, 1, CON, hole2.mp.get_bvect_spher()) ;
    beta2.set_etat_qcq() ;
    Cmp beta2_r (hole2.mp) ;
    Cmp beta2_t (hole2.mp) ;
    Cmp beta2_p (hole2.mp) ;
    beta2_r = hole2.beta_auto()(1) ;
    beta2_t = hole2.beta_auto()(2) ;
    beta2_p = hole2.beta_auto()(3) ;
    beta2.set(0) = beta2_r ;
    beta2.set(1) = beta2_t ;
    beta2.set(2) = beta2_p ;
    beta2.change_triad(hole2.mp.get_bvect_cart()) ;

    // We solve :
    poisson_vect_binaire (1./3., source1, source2, 
	lim_x_un, lim_y_un, lim_z_un, 
	lim_x_deux, lim_y_deux, lim_z_deux, 
	beta1, beta2, 0, precision) ;
    
    for (int i=0 ; i<3 ; i++) {
	beta1.set(i).raccord(1) ;
	beta2.set(i).raccord(1) ;
    }
    
    beta1.change_triad(hole1.mp.get_bvect_spher()) ;
    beta2.change_triad(hole2.mp.get_bvect_spher()) ;
    
    Vector shift1 (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    Scalar shift1_r (beta1(0)) ;
    Scalar shift1_t (beta1(1)) ;
    Scalar shift1_p (beta1(2)) ;
    shift1.set(1) = shift1_r ;
    shift1.set(2) = shift1_t ;
    shift1.set(3) = shift1_p ;

    Vector shift2 (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
    Scalar shift2_r (beta2(0)) ;
    Scalar shift2_t (beta2(1)) ;
    Scalar shift2_p (beta2(2)) ;
    shift2.set(1) = shift2_r ;
    shift2.set(2) = shift2_t ;
    shift2.set(3) = shift2_p ;


    // Regularisation
    Vector shift1_new (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    Vector shift2_new (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;

    shift1_new = relax*shift1 + (1-relax)*beta_un_old ;
    shift2_new = relax*shift2 + (1-relax)*beta_deux_old ;


    double ttime = hole1.the_time[hole1.jtime] ;
    hole1.beta_auto_evol.update(shift1_new, hole1.jtime, ttime) ;
    hole2.beta_auto_evol.update(shift2_new, hole2.jtime, ttime) ;

    shift1_new.change_triad(hole1.mp.get_bvect_cart()) ;
    shift2_new.change_triad(hole2.mp.get_bvect_cart()) ;
    cout << "shift_x" << endl << norme(shift1_new(1)) << endl ;
    cout << "shift_y" << endl << norme(shift1_new(2)) << endl ;
    cout << "shift_z" << endl << norme(shift1_new(3)) << endl ;

    // Regularisation of the shifts if necessary
    // -----------------------------------------

    int nnt = hole1.mp.get_mg()->get_nt(1) ;
    int nnp = hole1.mp.get_mg()->get_np(1) ;
    
    int check ;
    check = 0 ;
    for (int k=0; k<nnp; k++)
	for (int j=0; j<nnt; j++){
	    if ((hole1.n_auto()+hole1.n_comp()).val_grid_point(1, k, j , 0) < 1e-4){
		check = 1 ;
		break ;
	    }
	}

    if (check == 1){
	double diff_un = hole1.regularisation (hole1.beta_auto(), 
					 hole2.beta_auto(), omega) ;
	double diff_deux = hole2.regularisation (hole2.beta_auto(), 
					   hole1.beta_auto(), omega) ;
	hole1.regul = diff_un ;
	hole2.regul = diff_deux ;
    }
    
    else {
	hole1.regul = 0. ;
	hole2.regul = 0. ;
    }
    
}
