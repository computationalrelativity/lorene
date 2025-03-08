/*
 *  Methods for computing the extrinsic curvature tensor for a Bin_ns_bh
 *
 */

/*
 *   Copyright (c) 2002  Philippe Grandclement, Keisuke Taniguchi,
 *              Eric Gourgoulhon
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.12  2016/12/05 16:17:46  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.11  2014/10/13 08:52:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:13:01  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2008/09/26 08:44:04  p_grandclement
 * Mixted binaries with non vanishing spin
 *
 * Revision 1.8  2008/08/19 06:41:59  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.7  2007/04/26 14:14:59  f_limousin
 * The function fait_tkij now have default values for bound_nn and lim_nn
 *
 * Revision 1.6  2007/04/26 13:16:23  f_limousin
 * Correction of an error in the computation of grad_n_tot and grad_psi_tot
 *
 * Revision 1.5  2007/04/24 20:13:53  f_limousin
 * Implementation of Dirichlet and Neumann BC for the lapse
 *
 * Revision 1.4  2005/12/01 12:59:10  p_grandclement
 * Files for bin_ns_bh project
 *
 * Revision 1.3  2005/08/29 15:10:15  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.2  2004/05/27 12:41:00  p_grandclement
 * correction of some shadowed variables
 *
 * Revision 1.1  2003/02/13 16:40:25  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cmath>

// Lorene headers
#include "bin_ns_bh.h"
#include "proto.h"
#include "utilitaires.h"

#include "graphique.h"

namespace Lorene {

/**
 * Computes (LB)^{ij} auto.
 * To be used only when computing the total extrinsic curvature tensor
 * in the case of a Bin_ns_bh
 **/

void Et_bin_nsbh::fait_taij_auto() {

  Tenseur copie_shift (shift_auto) ;
  copie_shift.change_triad(mp.get_bvect_cart()) ;

  if (shift_auto.get_etat() == ETATZERO)
      taij_auto.set_etat_zero() ;
  else {
    Tenseur grad (copie_shift.gradient()) ;
    Tenseur trace (contract (grad,0, 1)) ;
    
    taij_auto.set_triad(mp.get_bvect_cart()) ;
    taij_auto.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++)
      for (int j=i ; j<3 ; j++)
	taij_auto.set(i, j) = grad(i, j)+grad(j, i) ;
      for (int i=0 ; i<3 ; i++)
	  taij_auto.set(i, i) -= 2./3.*trace() ;
    }

  taij_auto.change_triad(ref_triad) ;
}

void Bin_ns_bh::fait_decouple () {

  int nz_bh = hole.mp.get_mg()->get_nzone() ;
  int nz_ns = star.mp.get_mg()->get_nzone() ;
  
  // On determine R_limite (pour le moment en tout cas...) :
  double distance = star.mp.get_ori_x()-hole.mp.get_ori_x() ;
  double lim_ns = distance/2. ;
  double lim_bh = distance/2. ;
  double int_ns = lim_ns/3. ;
  double int_bh = lim_bh/3. ;
 
  // Les fonctions totales :
  Cmp decouple_ns (star.mp) ;
  decouple_ns.allocate_all() ;
  Cmp decouple_bh (hole.mp) ;
  decouple_bh.allocate_all() ;
  
  Mtbl xabs_ns (star.mp.xa) ;
  Mtbl yabs_ns (star.mp.ya) ;
  Mtbl zabs_ns (star.mp.za) ;
  
  Mtbl xabs_bh (hole.mp.xa) ;
  Mtbl yabs_bh (hole.mp.ya) ;
  Mtbl zabs_bh (hole.mp.za) ;
  
  double xabs, yabs, zabs, air_ns, air_bh, theta, phi ;
  
  // On boucle sur les autres zones :
  for (int l=0 ; l<nz_ns ; l++) {
    int nr = star.mp.get_mg()->get_nr (l) ;
    
    if (l==nz_ns-1)
      nr -- ;
    
    int np = star.mp.get_mg()->get_np (l) ;
    int nt = star.mp.get_mg()->get_nt (l) ;
    
    for (int k=0 ; k<np ; k++)
      for (int j=0 ; j<nt ; j++)
	for (int i=0 ; i<nr ; i++) {
	  
	  xabs = xabs_ns (l, k, j, i) ;
	  yabs = yabs_ns (l, k, j, i) ;
	  zabs = zabs_ns (l, k, j, i) ;
	  
	  // les coordonnees du point :
	  star.mp.convert_absolute 
	    (xabs, yabs, zabs, air_ns, theta, phi) ;
	  hole.mp.convert_absolute 
	    (xabs, yabs, zabs, air_bh, theta, phi) ;
	  
	  if (air_ns <= lim_ns)
	    if (air_ns < int_ns)
	      decouple_ns.set(l, k, j, i) = 1 ;
	    else
	      decouple_ns.set(l, k, j, i) = 
		0.5*pow(cos((air_ns-int_ns)*M_PI/2./(lim_ns-int_ns)), 2.)+0.5
		 ;
	  else 
	    if (air_bh <= lim_bh)
	      if (air_bh < int_bh)
		decouple_ns.set(l, k, j, i) = 0 ;
	      else 
		decouple_ns.set(l, k, j, i) = 0.5*
		  pow(sin((air_bh-int_bh)*M_PI/2./(lim_bh-int_bh)), 2.)
		   ;
	    else
	      // On est loin des deux trous :
	      decouple_ns.set(l, k, j, i) = 0.5 ;
	}
    
    // Cas infini :
    if (l==nz_ns-1)
      for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	  decouple_ns.set(nz_ns-1, k, j, nr) = 0.5 ;
  }
  
  for (int l=0 ; l<nz_bh ; l++) {
    int nr = hole.mp.get_mg()->get_nr (l) ;
    
    if (l==nz_bh-1)
      nr -- ;
    
    int np = hole.mp.get_mg()->get_np (l) ;
    int nt = hole.mp.get_mg()->get_nt (l) ;
    
    for (int k=0 ; k<np ; k++)
      for (int j=0 ; j<nt ; j++)
	for (int i=0 ; i<nr ; i++) {
	  
	  xabs = xabs_bh (l, k, j, i) ;
	  yabs = yabs_bh (l, k, j, i) ;
	  zabs = zabs_bh (l, k, j, i) ;
	  
	  // les coordonnees du point  :
	  star.mp.convert_absolute 
	    (xabs, yabs, zabs, air_ns, theta, phi) ;
	  hole.mp.convert_absolute 
	    (xabs, yabs, zabs, air_bh, theta, phi) ;
	  
	  if (air_bh <= lim_bh)
	    if (air_bh < int_bh)
	      decouple_bh.set(l, k, j, i) = 1 ;
	    else
	      decouple_bh.set(l, k, j, i) = 0.5*
		pow(cos((air_bh-int_bh)*M_PI/2./(lim_bh-int_bh)), 2.)+0.5 ;
	  else 
	    if (air_ns <= lim_ns)
	      if (air_ns < int_ns)
		decouple_bh.set(l, k, j, i) = 0 ;
	      else
		decouple_bh.set(l, k, j, i) = 0.5*
		  pow(sin((air_ns-int_ns)*M_PI/2./(lim_ns-int_ns)), 2.) ;
	  
	    else
	      // On est loin des deux trous :
	      decouple_bh.set(l, k, j, i) = 0.5 ;
	}
    
    // Cas infini :
    if (l==nz_bh-1)
      for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	  decouple_bh.set(nz_bh-1, k, j, nr) = 0.5 ;
  }
  
  decouple_ns.std_base_scal() ;
  decouple_bh.std_base_scal() ;
  
  star.decouple = decouple_ns ;
  hole.decouple = decouple_bh ;
}

//********************************************************
//calcul de kij total. (la regularisation ayant ete faite)
//********************************************************
void Bin_ns_bh::fait_tkij (int bound_nn, double lim_nn) {
    
  fait_decouple() ;
  
  double norme_hole = 0 ;
  double norme_star = 0 ;
  
  for (int i=0 ; i<3 ; i++) {
      norme_hole += max(norme(hole.get_shift_auto()(i))) ;
      norme_star += max(norme(star.get_shift_auto()(i))) ;
  }
  
#ifndef NDEBUG
  bool zero_shift_hole = (norme_hole <1e-14) ? true : false ; 
#endif
  bool zero_shift_star = (norme_star <1e-14) ? true : false ; 
   
  assert (zero_shift_hole == zero_shift_star) ;

  if (zero_shift_star == true) {
    hole.tkij_tot.set_etat_zero() ;
    hole.tkij_auto.set_etat_zero() ;
    
    hole.taij_tot.set_etat_zero() ;
    hole.taij_auto.set_etat_zero() ;
    hole.taij_comp.set_etat_zero() ;

    star.tkij_auto.set_etat_zero() ;
    star.tkij_comp.set_etat_zero() ;
    star.akcar_comp.set_etat_zero() ;
    star.akcar_auto.set_etat_zero() ;
  }
  else {

    if (bound_nn < 0){
      int nnt = hole.mp.get_mg()->get_nt(1) ;
      int nnp = hole.mp.get_mg()->get_np(1) ;
      
      for (int k=0; k<nnp; k++)
	for (int j=0; j<nnt; j++){
	  if (fabs((hole.n_auto+hole.n_comp)()(1, k, j , 0)) < 1e-2){
	    bound_nn = 0 ;
	    lim_nn = 0. ;
	    break ;
	  }
	}
    }  
    
    if (bound_nn != 0 || lim_nn != 0){
      
      hole.fait_taij_auto () ;  
      star.fait_taij_auto() ;

      // On trouve les trucs du compagnon
      hole.taij_comp.set_etat_qcq() ;
      // Pas de membre pour NS
      Tenseur_sym ns_taij_comp (star.mp, 2, CON, ref_triad) ;
      ns_taij_comp.set_etat_qcq() ;
      
      Tenseur_sym copie_ns (star.taij_auto) ;
      copie_ns.dec2_dzpuis() ;
      Tenseur_sym copie_bh (hole.taij_auto) ;
      copie_bh.dec2_dzpuis() ;
      
      // Les importations :
      if (hole.taij_auto.get_etat() == ETATZERO)
	ns_taij_comp.set_etat_zero() ;
      else {
	ns_taij_comp.set(0, 0).import(copie_bh(0, 0)) ;
	ns_taij_comp.set(0, 1).import(copie_bh(0, 1)) ;
	ns_taij_comp.set(0, 2).import(copie_bh(0, 2)) ;
	ns_taij_comp.set(1, 1).import(copie_bh(1, 1)) ;
	ns_taij_comp.set(1, 2).import(copie_bh(1, 2)) ;
	ns_taij_comp.set(2, 2).import(copie_bh(2, 2)) ;
	ns_taij_comp.set_triad(*copie_bh.get_triad()) ;
	ns_taij_comp.change_triad(star.ref_triad) ;
      }
      
      if (star.taij_auto.get_etat() == ETATZERO)
	hole.taij_comp.set_etat_zero() ;
      else {
	hole.taij_comp.set(0, 0).import(copie_ns(0, 0)) ;
	hole.taij_comp.set(0, 1).import(copie_ns(0, 1)) ;
	hole.taij_comp.set(0, 2).import(copie_ns(0, 2)) ;
	hole.taij_comp.set(1, 1).import(copie_ns(1, 1)) ;
	hole.taij_comp.set(1, 2).import(copie_ns(1, 2)) ;
	hole.taij_comp.set(2, 2).import(copie_ns(2, 2)) ;
	hole.taij_comp.set_triad(*copie_ns.get_triad()) ;
	hole.taij_comp.change_triad (hole.mp.get_bvect_cart()) ;
      }
      
      hole.taij_comp.set_std_base() ;
      ns_taij_comp.set_std_base() ;
      hole.taij_comp.inc2_dzpuis() ;
      ns_taij_comp.inc2_dzpuis() ;
      
      // Et enfin les trucs totaux...
      hole.taij_tot = hole.taij_auto + hole.taij_comp ;
      Tenseur_sym ns_taij_tot (star.taij_auto + ns_taij_comp) ;
      star.taij_tot = ns_taij_tot ;
      
      // Computation of tkij
      star.tkij_tot.set_etat_qcq() ;
      star.tkij_auto.set_etat_qcq() ;
      star.tkij_comp.set_etat_qcq() ;
      hole.tkij_tot.set_etat_qcq() ;
      hole.tkij_auto.set_etat_qcq() ;
      
      for (int i = 0 ; i<3 ; i++)
	for (int j = i ; j<3 ; j++) {
	  star.tkij_tot.set(i,j) = 0.5*star.taij_tot(i,j)/star.nnn() ;
	  //star.tkij_auto.set(i,j) = 0.5*star.taij_tot(i,j)/star.nnn() ;
	  //star.tkij_comp.set(i,j) = 0.5*ns_taij_comp(i,j)/star.nnn() ;
	  hole.tkij_tot.set(i,j) = 0.5*hole.taij_tot(i,j)/hole.n_tot() ;
	  //hole.tkij_auto.set(i,j) = 0.5*hole.taij_auto(i,j)/hole.n_tot() ;
	}
    
    for (int lig=0 ; lig<3 ; lig++)
      for (int col=lig ; col<3 ; col++) {
	star.tkij_auto.set(lig, col) = star.tkij_tot(lig, col)*
	  star.decouple ;
	star.tkij_comp.set(lig, col) = star.tkij_tot(lig, col)*
	  (1-star.decouple) ;
	hole.tkij_auto.set(lig, col) = hole.tkij_tot(lig, col)*
	  hole.decouple ;
      }
    star.tkij_auto.set_std_base() ;
    star.tkij_comp.set_std_base() ;
    hole.tkij_auto.set_std_base() ;
    }
    else {
      
      hole.tkij_tot.set_etat_qcq() ;
      star.tkij_tot.set_etat_qcq() ;
      
      // On construit a_ij a partir du shift ...
      // taij tot doit etre nul sur l'horizon.
      hole.fait_taij_auto () ;  
      star.fait_taij_auto() ;
      
      // On trouve les trucs du compagnon
      hole.taij_comp.set_etat_qcq() ;
      // Pas de membre pour NS
      Tenseur_sym ns_taij_comp (star.mp, 2, CON, ref_triad) ;
      ns_taij_comp.set_etat_qcq() ;
      
      Tenseur_sym copie_ns (star.taij_auto) ;
    copie_ns.dec2_dzpuis() ;
    Tenseur_sym copie_bh (hole.taij_auto) ;
    copie_bh.dec2_dzpuis() ;
    
    // Les importations :
    if (hole.taij_auto.get_etat() == ETATZERO)
      ns_taij_comp.set_etat_zero() ;
    else {
      ns_taij_comp.set(0, 0).import_asymy(copie_bh(0, 0)) ;
      ns_taij_comp.set(0, 1).import_symy(copie_bh(0, 1)) ;
      ns_taij_comp.set(0, 2).import_asymy(copie_bh(0, 2)) ;
      ns_taij_comp.set(1, 1).import_asymy(copie_bh(1, 1)) ;
      ns_taij_comp.set(1, 2).import_symy(copie_bh(1, 2)) ;
      ns_taij_comp.set(2, 2).import_asymy(copie_bh(2, 2)) ;
      ns_taij_comp.set_triad(*copie_bh.get_triad()) ;
      ns_taij_comp.change_triad(star.ref_triad) ;
    }
    
    if (star.taij_auto.get_etat() == ETATZERO)
      hole.taij_comp.set_etat_zero() ;
    else {
      hole.taij_comp.set(0, 0).import_asymy(copie_ns(0, 0)) ;
      hole.taij_comp.set(0, 1).import_symy(copie_ns(0, 1)) ;
      hole.taij_comp.set(0, 2).import_asymy(copie_ns(0, 2)) ;
      hole.taij_comp.set(1, 1).import_asymy(copie_ns(1, 1)) ;
      hole.taij_comp.set(1, 2).import_symy(copie_ns(1, 2)) ;
      hole.taij_comp.set(2, 2).import_asymy(copie_ns(2, 2)) ;
      hole.taij_comp.set_triad(*copie_ns.get_triad()) ;
      hole.taij_comp.change_triad (hole.mp.get_bvect_cart()) ;
    }
    
    hole.taij_comp.set_std_base() ;
    ns_taij_comp.set_std_base() ;
    hole.taij_comp.inc2_dzpuis() ;
    ns_taij_comp.inc2_dzpuis() ;
       
    // Et enfin les trucs totaux...
    hole.taij_tot = hole.taij_auto + hole.taij_comp ;
    Tenseur_sym ns_taij_tot (star.taij_auto + ns_taij_comp) ;
    star.taij_tot = ns_taij_tot ;
    
    int nz_ns = star.mp.get_mg()->get_nzone() ;
    Cmp ntot_ns (star.get_nnn()()) ;
    
    Cmp ntot_bh (hole.n_tot()) ;
    //des_coupe_z (ntot_bh, 0, -5, 0, -2.5, 2.5, "N") ;
    ntot_bh = division_xpun (ntot_bh, 0) ;
    ntot_bh.raccord(1) ;

    //des_coupe_z (ntot_bh, 0, -5, 0, -2.5, 2.5, "N/ (x+1)") ;
    
    // Boucle sur les composantes :
    for (int lig = 0 ; lig<3 ; lig++)
      for (int col = lig ; col<3 ; col++) {
	
	// Dans la grille du BH (pas de pb sauf pres horizon :
	Cmp auxi_bh (hole.taij_tot(lig, col)/2.) ;
	auxi_bh.dec2_dzpuis() ;
	auxi_bh = division_xpun (auxi_bh, 0) ;
	
	//des_coupe_z (auxi_bh, 0, -10, 20, -7, 7, "Aij/ (x+1)") ;
	auxi_bh = auxi_bh / ntot_bh ;
	auxi_bh.raccord(1) ;

	// Pour la NS on doit faire attention a pas etre pres du trou
	Cmp auxi_ns (star.mp) ;
	auxi_ns.allocate_all() ;
	
	// copie :
	Cmp copie_ns_bis (ns_taij_tot(lig, col)) ;
	copie_ns_bis.dec2_dzpuis() ;
	
	// Double le rayon limite :
	double lim_bh = hole.mp.get_alpha()[1] + hole.mp.get_beta()[1] ;
	
	Mtbl xabs_ns (star.mp.xa) ;
	Mtbl yabs_ns (star.mp.ya) ;
	Mtbl zabs_ns (star.mp.za) ;
	
	double xabs, yabs, zabs, air, theta, phi ;
	
	// On boucle sur les zones
	for (int l=0 ; l<nz_ns ; l++) {
	  
	  int nr = star.mp.get_mg()->get_nr (l) ;
	  
	  if (l==nz_ns-1)
	    nr -- ;
	  
	  int np = star.mp.get_mg()->get_np (l) ;
	  int nt = star.mp.get_mg()->get_nt (l) ;
	  
	  for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
	      for (int i=0 ; i<nr ; i++) {
		
		xabs = xabs_ns (l, k, j, i) ;
		yabs = yabs_ns (l, k, j, i) ;
		zabs = zabs_ns (l, k, j, i) ;
		
		// les coordonnees du point vis a vis BH :
		hole.mp.convert_absolute 
		  (xabs, yabs, zabs, air, theta, phi) ;
		
		if (air >= lim_bh)
		  // On est loin du trou 2 : pas de pb :
		  auxi_ns.set(l, k, j, i) = 
		    copie_ns_bis(l, k, j, i) / ntot_ns (l, k, j, i)/2. ;
		else 
		  // On est pres du trou (faut pas tomber dedans !) :
		  auxi_ns.set(l, k, j, i) =  auxi_bh.val_point (air, theta, phi) ;
	      }
	  
	  // Cas infini :
	  if (l==nz_ns-1)
	    for (int k=0 ; k<np ; k++)
	      for (int j=0 ; j<nt ; j++)
		auxi_ns.set(nz_ns-1, k, j, nr) = 0 ;
	}
	
	
	star.tkij_tot.set(lig, col) = auxi_ns ;
	hole.tkij_tot.set(lig, col) = auxi_bh ;
      }

    star.tkij_tot.set_std_base() ;
    hole.tkij_tot.set_std_base() ;
    star.tkij_tot.inc2_dzpuis() ;

    hole.tkij_tot.inc2_dzpuis() ;
     
    //Cmp dessin_un (hole.get_tkij_tot()(0,0)) ;
    //des_coupe_z (dessin_un, 0, -10, 20, -7, 7, "Kij tot BH") ;
    //Cmp dessin_deux (ns_tkij_tot(0,0)) ;
    //des_coupe_z (dessin_deux, 0, -10, 20, -7, 7, "Kij tot NS") ;
    
    star.tkij_auto.set_etat_qcq() ;
    star.tkij_comp.set_etat_qcq() ;
    hole.tkij_auto.set_etat_qcq() ;
    
    for (int lig=0 ; lig<3 ; lig++)
      for (int col=lig ; col<3 ; col++) {
	star.tkij_auto.set(lig, col) = star.tkij_tot(lig, col)*
	  star.decouple ;
	star.tkij_comp.set(lig, col) = star.tkij_tot(lig, col)*
	  (1-star.decouple) ;
	hole.tkij_auto.set(lig, col) = hole.tkij_tot(lig, col)*
	  hole.decouple ;
      }
    star.tkij_auto.set_std_base() ;
    star.tkij_comp.set_std_base() ;
    hole.tkij_auto.set_std_base() ;

    }

    // On doit mettre a jour les champs akcar de NS :
    star.akcar_auto.set_etat_qcq() ; 
    star.akcar_auto.set() = 0 ; 
    
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	star.akcar_auto.set() += star.tkij_auto(i, j) % star.tkij_auto(i, j) ;
    
    star.akcar_auto.set_std_base() ;
    star.akcar_auto = star.a_car % star.akcar_auto ; 
    
    star.akcar_comp.set_etat_qcq() ;
    star.akcar_comp.set() = 0 ;
    
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	star.akcar_comp.set() += star.tkij_auto(i, j) % star.tkij_comp(i, j) ;
    
    star.akcar_comp.set_std_base() ;
    star.akcar_comp = star.a_car % star.akcar_comp ;
  }
}
}
