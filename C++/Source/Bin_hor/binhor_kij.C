/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo                       
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


char binhor_kij_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2006/05/24 16:56:37  f_limousin
 * Many small modifs.
 *
 * Revision 1.7  2005/09/13 18:33:15  f_limousin
 * New function vv_bound_cart_bin(double) for computing binaries with
 * berlin condition for the shift vector.
 * Suppress all the symy and asymy in the importations.
 *
 * Revision 1.6  2005/04/29 14:02:44  f_limousin
 * Important changes : manage the dependances between quantities (for
 * instance psi and psi4). New function write_global(ost).
 *
 * Revision 1.5  2005/03/10 17:21:52  f_limousin
 * Add the Berlin boundary condition for the shift.
 * Some changes to avoid warnings.
 *
 * Revision 1.4  2005/03/03 13:49:35  f_limousin
 * Add the spectral bases for both Scalars decouple.
 *
 * Revision 1.3  2005/02/07 10:48:00  f_limousin
 * The extrinsic curvature can now be computed in the case N=0 on the
 * horizon.
 *
 * Revision 1.2  2004/12/31 15:41:54  f_limousin
 * Correction of an error
 *
 * Revision 1.1  2004/12/29 16:12:03  f_limousin
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
#include "tenseur.h"
#include "tensor.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
//#include "graphique.h"

void Bin_hor::extrinsic_curvature () {
    

    double ttime = hole1.the_time[hole1.jtime] ;
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
    
    Sym_tensor aa_auto_un (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    Sym_tensor aa_auto_deux (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
       
    if (check == 0){
      
      // Computation of A^{ij}_auto
      
      aa_auto_un = ( hole1.beta_auto().ope_killing_conf(hole1.met_gamt) + 
		      hole1.gamt_point*hole1.decouple ) / (2.* hole1.nn()) ; 
      
      aa_auto_deux = ( hole2.beta_auto().ope_killing_conf(hole2.met_gamt) + 
		       hole2.gamt_point*hole2.decouple ) / (2.* hole2.nn()) ;      
      
      
      aa_auto_un.change_triad(hole1.mp.get_bvect_cart()) ;
      aa_auto_deux.change_triad(hole2.mp.get_bvect_cart()) ;
      
      for (int i=1 ; i<=3 ; i++)
	for (int j=i ; j<=3 ; j++) {
	  if (aa_auto_un(i,j).get_etat() != ETATZERO)
	    aa_auto_un.set(i, j).raccord(3) ;
	  if (aa_auto_deux(i,j).get_etat() != ETATZERO)
	    aa_auto_deux.set(i, j).raccord(3) ;
	}
      
      aa_auto_un.change_triad(hole1.mp.get_bvect_spher()) ;
      aa_auto_deux.change_triad(hole2.mp.get_bvect_spher()) ;
      
      hole1.aa_auto_evol.update(aa_auto_un, hole1.jtime, ttime) ;
      hole2.aa_auto_evol.update(aa_auto_deux, hole2.jtime, ttime) ;
      
      
      
      // Computation of A^{ij}_comp
      
      aa_auto_un.dec_dzpuis(2) ;
      aa_auto_deux.dec_dzpuis(2) ;
      
      Sym_tensor aa_comp_un (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
      aa_comp_un.set_etat_qcq() ;
      Sym_tensor aa_comp_deux (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
      aa_comp_deux.set_etat_qcq() ;
      
      aa_auto_deux.change_triad(hole2.mp.get_bvect_cart()) ;
      aa_auto_deux.change_triad(hole1.mp.get_bvect_cart()) ;
      assert(*(aa_auto_deux.get_triad()) == *(aa_comp_un.get_triad())) ;
      
      // importations :
      for (int i=1 ; i<=3 ; i++)
	for (int j=i ; j<=3 ; j++) {
	  aa_comp_un.set(i, j).import(aa_auto_deux(i, j)) ;
	  aa_comp_un.set(i, j).set_spectral_va().set_base(aa_auto_deux(i, j).
	     					  get_spectral_va().get_base()) ;
	}
      
      aa_comp_un.inc_dzpuis(2) ;
      aa_comp_un.change_triad(hole1.mp.get_bvect_spher()) ;
      
      aa_auto_un.change_triad(hole1.mp.get_bvect_cart()) ;
       aa_auto_un.change_triad(hole2.mp.get_bvect_cart()) ;
       assert(*(aa_auto_un.get_triad()) == *(aa_comp_deux.get_triad())) ;
       // importations :
       for (int i=1 ; i<=3 ; i++)
         for (int j=i ; j<=3 ; j++) {
           aa_comp_deux.set(i, j).import(aa_auto_un(i, j)) ;
           aa_comp_deux.set(i, j).set_spectral_va().set_base(aa_auto_un(i, j).
				   get_spectral_va().get_base()) ;
         }

       aa_comp_deux.inc_dzpuis(2) ;
       aa_comp_deux.change_triad(hole2.mp.get_bvect_spher()) ;
       
       /*              
       // Computation of A^{ij}_comp in the last domains
       // -----------------------------------------------
       
       int nz = hole1.mp.get_mg()->get_nzone() ;

       Sym_tensor aa_comp_un_zec (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
       aa_comp_un_zec.set_etat_qcq() ;
       Sym_tensor aa_comp_deux_zec (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
       aa_comp_deux_zec.set_etat_qcq() ;

       aa_comp_un_zec = ( hole1.beta_comp().ope_killing_conf(hole1.met_gamt) +
		  hole1.gamt_point*(1.-hole1.decouple) ) / (2.* hole1.nn()) ;

       aa_comp_deux_zec =( hole2.beta_comp().ope_killing_conf(hole2.met_gamt) +
		hole2.gamt_point*(1.-hole2.decouple) ) / (2.* hole2.nn()) ;

       for (int i=1 ; i<=3 ; i++)
         for (int j=i ; j<=3 ; j++)          
	   for (int l=nz-1 ; l<=nz-1 ; l++) {
	     if (aa_comp_un.set(i,j).get_etat() == ETATQCQ)
	       aa_comp_un.set(i,j).set_domain(l) = 
		                   aa_comp_un_zec(i,j).domain(l) ;
	     if (aa_comp_deux.set(i,j).get_etat() == ETATQCQ)
	       aa_comp_deux.set(i,j).set_domain(l)=
	                           aa_comp_deux_zec(i,j).domain(l) ;
	   }
       */      

       hole1.aa_comp_evol.update(aa_comp_un, hole1.jtime, ttime) ;
       hole2.aa_comp_evol.update(aa_comp_deux, hole2.jtime, ttime) ;
       
       
       // Computation of A^{ij}_ total
       hole1.aa_evol.update(hole1.aa_auto() + hole1.aa_comp(), 
			    hole1.jtime, ttime) ;
       hole2.aa_evol.update(hole2.aa_auto() + hole2.aa_comp(), 
			    hole2.jtime, ttime) ;
       hole1.aa_nn.update((hole1.aa_auto() + hole1.aa_comp())*2.*hole1.nn(), 
			    hole1.jtime, ttime) ;
       hole2.aa_nn.update((hole2.aa_auto() + hole2.aa_comp())*2.*hole2.nn(), 
			    hole2.jtime, ttime) ;
       
    }
   else {

       // Computation of A^{ij}_auto

       aa_auto_un = ( hole1.beta_auto().ope_killing_conf(hole1.met_gamt) + 
		      hole1.gamt_point*hole1.decouple ) ;            
       aa_auto_deux = ( hole2.beta_auto().ope_killing_conf(hole2.met_gamt) + 
			hole2.gamt_point*hole2.decouple ) ;          
       
       aa_auto_un.change_triad(hole1.mp.get_bvect_cart()) ;
       aa_auto_deux.change_triad(hole2.mp.get_bvect_cart()) ;

       for (int i=1 ; i<=3 ; i++)
	   for (int j=1 ; j<=3 ; j++) {
	       if (aa_auto_un(i,j).get_etat() != ETATZERO)
		   aa_auto_un.set(i, j).raccord(3) ;
	       if (aa_auto_deux(i,j).get_etat() != ETATZERO)
		   aa_auto_deux.set(i, j).raccord(3) ;
	   }

       // Computation of A^{ij}_comp
       
       Sym_tensor aa_auto_1 (aa_auto_un) ;
       Sym_tensor aa_auto_2 (aa_auto_deux) ;
       
       aa_auto_1.dec_dzpuis(2) ;
       aa_auto_2.dec_dzpuis(2) ;
       
       Sym_tensor aa_comp_un (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
       aa_comp_un.set_etat_qcq() ;
       Sym_tensor aa_comp_deux (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
       aa_comp_deux.set_etat_qcq() ;
       
       aa_auto_2.change_triad(hole1.mp.get_bvect_cart()) ;
       assert(*(aa_auto_2.get_triad()) == *(aa_comp_un.get_triad())) ;
       // importations :
       aa_comp_un.set(1, 1).import(aa_auto_2(1, 1)) ;
       aa_comp_un.set(1, 2).import(aa_auto_2(1, 2)) ;
       aa_comp_un.set(1, 3).import(aa_auto_2(1, 3)) ;
       aa_comp_un.set(2, 2).import(aa_auto_2(2, 2)) ;
       aa_comp_un.set(2, 3).import(aa_auto_2(2, 3)) ;
       aa_comp_un.set(3, 3).import(aa_auto_2(3, 3)) ;
       
       aa_comp_un.std_spectral_base() ;
       aa_comp_un.inc_dzpuis(2) ;
          
       aa_auto_1.change_triad(hole2.mp.get_bvect_cart()) ;
       assert(*(aa_auto_1.get_triad()) == *(aa_comp_deux.get_triad())) ;
       // importations :
       aa_comp_deux.set(1, 1).import(aa_auto_1(1, 1)) ;
       aa_comp_deux.set(1, 2).import(aa_auto_1(1, 2)) ;
       aa_comp_deux.set(1, 3).import(aa_auto_1(1, 3)) ;
       aa_comp_deux.set(2, 2).import(aa_auto_1(2, 2)) ;
       aa_comp_deux.set(2, 3).import(aa_auto_1(2, 3)) ;
       aa_comp_deux.set(3, 3).import(aa_auto_1(3, 3)) ;
       
       aa_comp_deux.std_spectral_base() ;
       aa_comp_deux.inc_dzpuis(2) ;
         
       // Computation of A^{ij}_ total
       Sym_tensor aa_tot_un (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
       Sym_tensor aa_tot_deux (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
       aa_tot_un = aa_auto_un + aa_comp_un ;
       aa_tot_deux = aa_auto_deux + aa_comp_deux ;

       Sym_tensor temp_aa_tot1 (aa_tot_un) ;
       Sym_tensor temp_aa_tot2 (aa_tot_deux) ;
       
       temp_aa_tot1.change_triad(hole1.mp.get_bvect_spher()) ;
       temp_aa_tot2.change_triad(hole2.mp.get_bvect_spher()) ;
       hole1.aa_nn.update(temp_aa_tot1, hole1.jtime, ttime) ;
       hole2.aa_nn.update(temp_aa_tot2, hole2.jtime, ttime) ;


       // Regularisation
       // --------------

       Sym_tensor aa_un (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
       Sym_tensor aa_deux (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;     

      int nz_un = hole1.mp.get_mg()->get_nzone() ;
      int nz_deux = hole2.mp.get_mg()->get_nzone() ;
      
      Scalar ntot_un (hole1.n_auto()+hole1.n_comp()) ;
      ntot_un = division_xpun (ntot_un, 0) ;
      ntot_un.raccord(1) ;
      
      Scalar ntot_deux (hole2.n_auto()+hole2.n_comp()) ;
      ntot_deux = division_xpun (ntot_deux, 0) ;
      ntot_deux.raccord(1) ;
      
      // THE TWO Aij are aligned of not !
      double orientation_un = aa_auto_un.get_mp().get_rot_phi() ;
      assert ((orientation_un==0) || (orientation_un==M_PI)) ;
      double orientation_deux = aa_auto_deux.get_mp().get_rot_phi() ;
      assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
      int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;


      // Loop on the composants :
      for (int lig = 1 ; lig<=3 ; lig++)
	  for (int col = lig ; col<=3 ; col++) {
	      
	      // The orientation
	      int ind = 1 ;
	      if (lig !=3)
		  ind *= -1 ;
	      if (col != 3)
		  ind *= -1 ;
	      if (same_orient == 1)
		  ind = 1 ;
	    
	    // Close to hole one :
	    Scalar auxi_un (aa_tot_un(lig, col)/2.) ;
	    auxi_un.dec_dzpuis(2) ;
	    auxi_un = division_xpun (auxi_un, 0) ;
	    auxi_un = auxi_un / ntot_un ;
	    if (auxi_un.get_etat() != ETATZERO)
		auxi_un.raccord(1) ;

	    // Close to hole two :
	    Scalar auxi_deux (aa_tot_deux(lig, col)/2.) ;
	    auxi_deux.dec_dzpuis(2) ;
	    auxi_deux = division_xpun (auxi_deux, 0) ;
	    auxi_deux = auxi_deux / ntot_deux ;
	    if (auxi_deux.get_etat() != ETATZERO)
		auxi_deux.raccord(1) ;
	    
	    // copy :
	    Scalar copie_un (aa_tot_un(lig, col)) ;
	    copie_un.dec_dzpuis(2) ;
	    
	    Scalar copie_deux (aa_tot_deux(lig, col)) ;
	    copie_deux.dec_dzpuis(2) ;
	    
	    double lim_un = hole1.mp.get_alpha()[1] + hole1.mp.get_beta()[1] ;
	    double lim_deux = hole2.mp.get_alpha()[1] + hole2.mp.get_beta()[1] ;
	    
	    Mtbl xabs_un (hole1.mp.xa) ;
	    Mtbl yabs_un (hole1.mp.ya) ;
	    Mtbl zabs_un (hole1.mp.za) ;
	    
	    Mtbl xabs_deux (hole2.mp.xa) ;
	    Mtbl yabs_deux (hole2.mp.ya) ;
	    Mtbl zabs_deux (hole2.mp.za) ;
	    
	    double xabs, yabs, zabs, air, theta, phi ;

	    if (auxi_un.get_etat() != ETATZERO){	    
	    // Loop on the other zones :
	    for (int l=2 ; l<nz_un ; l++) {

		int nr = hole1.mp.get_mg()->get_nr (l) ;
		
		if (l==nz_un-1)
		    nr -- ;
		
		int np = hole1.mp.get_mg()->get_np (l) ;
		int nt = hole1.mp.get_mg()->get_nt (l) ;
		
		for (int k=0 ; k<np ; k++)
		    for (int j=0 ; j<nt ; j++)
			for (int i=0 ; i<nr ; i++) {
			    
			    xabs = xabs_un (l, k, j, i) ;
			    yabs = yabs_un (l, k, j, i) ;
			    zabs = zabs_un (l, k, j, i) ;

			    // coordinates of the point in 2 :
			    hole2.mp.convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
			
			    if (air >= lim_deux)
				// Far from hole two : no pb :
				auxi_un.set_grid_point(l, k, j, i) = 
				    copie_un.val_grid_point(l, k, j, i) / 
				    ntot_un.val_grid_point(l, k, j, i)/2. ;
			    else 
				// close to hole two :
				auxi_un.set_grid_point(l, k, j, i) = 
			    ind * auxi_deux.val_point (air, theta, phi) ;
				
			}
			    
		// Case infinity :
		if (l==nz_un-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    auxi_un.set_grid_point(nz_un-1, k, j, nr) = 0 ;
	    }
	    }

	    if (auxi_deux.get_etat() != ETATZERO){	    
	    // The second hole :
	    for (int l=2 ; l<nz_deux ; l++) {
		
		int nr = hole2.mp.get_mg()->get_nr (l) ;
		
		if (l==nz_deux-1)
		    nr -- ;
		
		int np = hole2.mp.get_mg()->get_np (l) ;
		int nt = hole2.mp.get_mg()->get_nt (l) ;
		
		for (int k=0 ; k<np ; k++)
		    for (int j=0 ; j<nt ; j++)
			for (int i=0 ; i<nr ; i++) {
			    
			    xabs = xabs_deux (l, k, j, i) ;
			    yabs = yabs_deux (l, k, j, i) ;
			    zabs = zabs_deux (l, k, j, i) ;
			    
			    // coordinates of the point in 2 :
			    hole1.mp.convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
			
			    if (air >= lim_un)
				// Far from hole one : no pb :
				auxi_deux.set_grid_point(l, k, j, i) = 
				    copie_deux.val_grid_point(l, k, j, i) / 
				    ntot_deux.val_grid_point(l, k, j, i) /2.;
			    else 
			    // close to hole one :
				auxi_deux.set_grid_point(l, k, j, i) = 
			      ind * auxi_un.val_point (air, theta, phi) ;
			    }
		// Case infinity :
		if (l==nz_deux-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    auxi_un.set_grid_point(nz_deux-1, k, j, nr) = 0 ;
		}
	    }

	    auxi_un.inc_dzpuis(2) ;
	    auxi_deux.inc_dzpuis(2) ;
	    
	    aa_un.set(lig, col) = auxi_un ;
	    aa_deux.set(lig, col) = auxi_deux ;
	  }

      aa_un.change_triad(hole1.mp.get_bvect_spher()) ;
      aa_deux.change_triad(hole2.mp.get_bvect_spher()) ;

      hole1.aa_evol.update(aa_un, hole1.jtime, ttime) ;
      hole2.aa_evol.update(aa_deux, hole2.jtime, ttime) ;

      aa_auto_un.change_triad(hole1.mp.get_bvect_spher()) ;
      aa_auto_deux.change_triad(hole2.mp.get_bvect_spher()) ;

      for (int lig=1 ; lig<=3 ; lig++)
	  for (int col=lig ; col<=3 ; col++) {
	      aa_auto_un.set(lig, col) = aa_un(lig, col)*hole1.decouple ;
	      aa_auto_deux.set(lig, col) = aa_deux(lig, col)*hole2.decouple ;
	  }

      hole1.aa_auto_evol.update(aa_auto_un, hole1.jtime, ttime) ;
      hole2.aa_auto_evol.update(aa_auto_deux, hole2.jtime, ttime) ;

   }

    Sym_tensor temp (hole1.aa()) ;
    temp.change_triad(hole1.mp.get_bvect_cart()) ;

}   


void Bin_hor::decouple () {
    
    int nz_un = hole1.mp.get_mg()->get_nzone() ;
    int nz_deux = hole2.mp.get_mg()->get_nzone() ;
    
    // We determine R_limite :
    double distance = hole1.mp.get_ori_x() - hole2.mp.get_ori_x() ;
    double lim_un = distance/2. ;
    double lim_deux = distance/2. ;
    double int_un = distance/6. ;
    double int_deux = distance/6. ;
    
    // The functions used.
    Scalar fonction_f_un (hole1.mp) ;
    fonction_f_un = 0.5*pow(
	cos((hole1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.)+0.5 ;
    fonction_f_un.std_spectral_base();
    
    Scalar fonction_g_un (hole1.mp) ;
    fonction_g_un = 0.5*pow
	(sin((hole1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.) ;
    fonction_g_un.std_spectral_base();
    
    Scalar fonction_f_deux (hole2.mp) ;
    fonction_f_deux = 0.5*pow(
	cos((hole2.mp.r-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.)+0.5 ;
    fonction_f_deux.std_spectral_base();
    
    Scalar fonction_g_deux (hole2.mp) ;
    fonction_g_deux = 0.5*pow(
	sin((hole2.mp.r-int_deux)*M_PI/2./(lim_un-int_deux)), 2.) ;
    fonction_g_deux.std_spectral_base();
    
    // The functions total :
    Scalar decouple_un (hole1.mp) ;
    decouple_un.allocate_all() ;
    Scalar decouple_deux (hole2.mp) ;
    decouple_deux.allocate_all() ;
    
    Mtbl xabs_un (hole1.mp.xa) ;
    Mtbl yabs_un (hole1.mp.ya) ;
    Mtbl zabs_un (hole1.mp.za) ;
	    
    Mtbl xabs_deux (hole2.mp.xa) ;
    Mtbl yabs_deux (hole2.mp.ya) ;
    Mtbl zabs_deux (hole2.mp.za) ;
	    
    double xabs, yabs, zabs, air_un, air_deux, theta, phi ;
	    
    for (int l=0 ; l<nz_un ; l++) {
	int nr = hole1.mp.get_mg()->get_nr (l) ;
		
	if (l==nz_un-1)
	    nr -- ;
		
	int np = hole1.mp.get_mg()->get_np (l) ;
	int nt = hole1.mp.get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_un (l, k, j, i) ;
		    yabs = yabs_un (l, k, j, i) ;
		    zabs = zabs_un (l, k, j, i) ;
			    
		    // Coordinates of the point
		    hole1.mp.convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    hole2.mp.convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		
		    if (air_un <= lim_un)
			if (air_un < int_un)
			    decouple_un.set_grid_point(l, k, j, i) = 1 ;
			else
			// Close to hole 1 :
			decouple_un.set_grid_point(l, k, j, i) = 
			    fonction_f_un.val_grid_point(l, k, j, i) ;
		    else 
			if (air_deux <= lim_deux)
			    if (air_deux < int_deux)
				decouple_un.set_grid_point(l, k, j, i) = 0 ;
			    else
			// Close to hole 2 :
				decouple_un.set_grid_point(l, k, j, i) = 
		fonction_g_deux.val_point (air_deux, theta, phi) ;
		
			else
			    // Far from each holes :
			    decouple_un.set_grid_point(l, k, j, i) = 0.5 ;
		}
			    
		// Case infinity :
		if (l==nz_un-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    decouple_un.set_grid_point(nz_un-1, k, j, nr)=0.5 ;
	    }
    
    for (int l=0 ; l<nz_deux ; l++) {
	int nr = hole2.mp.get_mg()->get_nr (l) ;
		
	if (l==nz_deux-1)
	    nr -- ;
		
	int np = hole2.mp.get_mg()->get_np (l) ;
	int nt = hole2.mp.get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_deux (l, k, j, i) ;
		    yabs = yabs_deux (l, k, j, i) ;
		    zabs = zabs_deux (l, k, j, i) ;
			    
		    // coordinates of the point  :
		    hole1.mp.convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    hole2.mp.convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		    
		    if (air_deux <= lim_deux)
			if (air_deux < int_deux)
			    decouple_deux.set_grid_point(l, k, j, i) = 1 ;
			else
			// close to hole two :
			decouple_deux.set_grid_point(l, k, j, i) = 
			    fonction_f_deux.val_grid_point(l, k, j, i) ;
		    else 
			if (air_un <= lim_un)
			    if (air_un < int_un)
				decouple_deux.set_grid_point(l, k, j, i) = 0 ;
			    else
			// close to hole one :
				decouple_deux.set_grid_point(l, k, j, i) = 
			 fonction_g_un.val_point (air_un, theta, phi) ;
		
			else
			    // Far from each hole :
			    decouple_deux.set_grid_point(l, k, j, i) = 0.5 ;
		}
			    
		// Case infinity :
		if (l==nz_deux-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			 decouple_deux.set_grid_point(nz_un-1, k, j, nr)=0.5 ;
   }
   
    decouple_un.std_spectral_base() ;
    decouple_deux.std_spectral_base() ;
    
    hole1.decouple = decouple_un ;
    hole2.decouple = decouple_deux ;
    
}
