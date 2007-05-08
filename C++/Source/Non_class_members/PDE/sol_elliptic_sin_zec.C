/*
 *   Copyright (c) 2004 Philippe Grandclement
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
 * Revision 1.6  2007/05/08 07:08:30  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.5  2007/05/06 10:48:12  p_grandclement
 * Modification of a few operators for the vorton project
 *
 * Revision 1.4  2005/11/30 11:09:08  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.3  2005/08/26 14:02:41  p_grandclement
 * Modification of the elliptic solver that matches with an oscillatory exterior solution
 * small correction in Poisson tau also...
 *
 * Revision 1.2  2004/08/24 09:14:44  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.1  2004/02/11 09:52:52  p_grandclement
 * Forgot one new file ...
 *

 * $Header$
 */

char sol_elliptic_sin_zec_C[] = "$Header$" ;

// Header C : 
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

// Headers Lorene :
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
#include "param_elliptic.h"
          

	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------

Mtbl_cf elliptic_solver_sin_zec  (const Param_elliptic& ope_var, 
				  const Mtbl_cf& source, double* amplis, double* phases) {

  // Verifications d'usage sur les zones
  int nz = source.get_mg()->get_nzone() ;
  assert (nz>1) ;
  assert (source.get_mg()->get_type_r(0) == RARE) ;
  assert (source.get_mg()->get_type_r(nz-1) == UNSURR) ;
  for (int l=1 ; l<nz-1 ; l++)
    assert(source.get_mg()->get_type_r(l) == FIN) ;
   
  // donnees sur la zone
  int nr, nt, np ;
  
  //Rangement des valeurs intermediaires 
  Tbl *so ;
  Tbl *sol_hom ;
  Tbl *sol_part ;
   
  
  // Rangement des solutions, avant raccordement
  Mtbl_cf solution_part(source.get_mg(), source.base) ;
  Mtbl_cf solution_hom_un(source.get_mg(), source.base) ;
  Mtbl_cf solution_hom_deux(source.get_mg(), source.base) ;
  Mtbl_cf resultat(source.get_mg(), source.base) ;

  solution_part.annule_hard() ;
  solution_hom_un.annule_hard() ;
  solution_hom_deux.annule_hard() ;
  resultat.annule_hard() ;
 
  // Computation of the SP and SH's in every domain but the ZEC ...
  int conte = 0 ;
  for (int zone=0 ; zone<nz-1 ; zone++) {
    nr = source.get_mg()->get_nr(zone) ;
    nt = source.get_mg()->get_nt(zone) ;
    np = source.get_mg()->get_np(zone) ;
     
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) {
	if (ope_var.operateurs[conte] != 0x0) {
	  // Calcul de la SH
	  sol_hom = new Tbl(ope_var.operateurs[conte]->get_solh()) ;

	  //Calcul de la SP
	  so = new Tbl(nr) ;
	  so->set_etat_qcq() ;
	  for (int i=0 ; i<nr ; i++)
	    so->set(i) = source(zone, k, j, i) ;
	  
	  sol_part = new Tbl(ope_var.operateurs[conte]->get_solp(*so)) ;

	  // Rangement dans les tableaux globaux ;
	  for (int i=0 ; i<nr ; i++) {
	    solution_part.set(zone, k, j, i) = (*sol_part)(i) ;
	    if (sol_hom->get_ndim()==1)
	      solution_hom_un.set(zone, k, j, i) = (*sol_hom)(i) ;
	    else
	      {
		solution_hom_un.set(zone, k, j, i) = (*sol_hom)(0,i) ;
		solution_hom_deux.set(zone, k, j, i) = (*sol_hom)(1,i) ;
	      }
	  }
	  
	  delete so ;
	  delete sol_hom ;
	  delete sol_part ;
	  
	}
	conte ++ ;
      }
  }
  
  //-------------------------------------------------
  // ON EST PARTI POUR LE RACCORD (Be carefull ....)
  //-------------------------------------------------
  
  // C'est pas simple toute cette sombre affaire...
  // Que le cas meme nombre de points dans chaque domaines...

  int start = 0 ;
  for (int k=0 ; k< source.get_mg()->get_np(0)+1; k++)
    for (int j=0 ; j<source.get_mg()->get_nt(0) ; j++) {
      if (ope_var.operateurs[start] != 0x0) {
	
	int taille = 2*nz - 2 ;
	Matrice systeme (taille, taille) ;
	systeme.set_etat_qcq() ;
	for (int i=0 ; i<taille ; i++)
	  for (int j2=0 ; j2<taille ; j2++)
	    systeme.set(i,j2) = 0 ;
	Tbl sec_membre (taille) ;
	sec_membre.set_etat_qcq() ;
	for (int i=0 ; i<taille ; i++)
	  sec_membre.set(i) = 0 ;
	//---------
	//  Noyau :
	//---------
	conte = start ;
      
	systeme.set(0,0) = ope_var.G_plus(0) * 
	  ope_var.operateurs[conte]->val_sh_one_plus() ;
	systeme.set(1,0) = 
	  ope_var.dG_plus(0) * ope_var.operateurs[conte]->val_sh_one_plus() +  
	  ope_var.G_plus(0) * ope_var.operateurs[conte]->der_sh_one_plus() ;
	
	sec_membre.set(0) -= ope_var.F_plus(0,k,j) + 
	  ope_var.G_plus(0) * ope_var.operateurs[conte]->val_sp_plus() ;
	sec_membre.set(1) -= ope_var.dF_plus(0,k,j) + 
	  ope_var.dG_plus(0) * ope_var.operateurs[conte]->val_sp_plus() + 
	  ope_var.G_plus(0) * ope_var.operateurs[conte]->der_sp_plus() ;

	//----------
	// SHELLS :
	//----------
	
	for (int l=1 ; l<nz-1 ; l++) {
	
	  // On se met au bon endroit :
	  int np_prec = source.get_mg()->get_np(l-1) ;
	  int nt_prec = source.get_mg()->get_nt(l-1) ;
	  conte += (np_prec+1)*nt_prec ;
	  
	  systeme.set(2*l-2, 2*l-1) = -ope_var.G_minus(l) * 
	    ope_var.operateurs[conte]->val_sh_one_minus() ;
	  systeme.set(2*l-2, 2*l) = - ope_var.G_minus(l) * 
	    ope_var.operateurs[conte]->val_sh_two_minus() ;
	  systeme.set(2*l-1, 2*l-1) = 
	    -ope_var.dG_minus(l)*ope_var.operateurs[conte]->val_sh_one_minus()-  
	    ope_var.G_minus(l)*ope_var.operateurs[conte]->der_sh_one_minus() ;
	  systeme.set(2*l-1, 2*l) =
	    -ope_var.dG_minus(l)*ope_var.operateurs[conte]->val_sh_two_minus()-  
	    ope_var.G_minus(l)*ope_var.operateurs[conte]->der_sh_two_minus() ;
	  
	  sec_membre.set(2*l-2) += ope_var.F_minus(l,k,j) + 
	    ope_var.G_minus(l) * ope_var.operateurs[conte]->val_sp_minus() ;
	  sec_membre.set(2*l-1) += ope_var.dF_minus(l,k,j) + 
	    ope_var.dG_minus(l) * ope_var.operateurs[conte]->val_sp_minus() + 
	    ope_var.G_minus(l) * ope_var.operateurs[conte]->der_sp_minus() ;
	  
	  // Valeurs en +1 :
	  systeme.set(2*l, 2*l-1) = ope_var.G_plus(l) * 
	    ope_var.operateurs[conte]->val_sh_one_plus() ;
	  systeme.set(2*l, 2*l) = ope_var.G_plus(l) * 
	    ope_var.operateurs[conte]->val_sh_two_plus() ;

	  systeme.set(2*l+1, 2*l-1) = 
	    ope_var.dG_plus(l)*ope_var.operateurs[conte]->val_sh_one_plus()+  
	    ope_var.G_plus(l)*ope_var.operateurs[conte]->der_sh_one_plus() ;
	  systeme.set(2*l+1, 2*l) =
	    ope_var.dG_plus(l)*ope_var.operateurs[conte]->val_sh_two_plus()+
	    ope_var.G_plus(l)*ope_var.operateurs[conte]->der_sh_two_plus() ;

	  sec_membre.set(2*l) -=  ope_var.F_plus(l,k,j) + 
	    ope_var.G_plus(l) * ope_var.operateurs[conte]->val_sp_plus();

	  sec_membre.set(2*l+1) -=  ope_var.dF_plus(l,k,j) + 
	         ope_var.dG_plus(l) * ope_var.operateurs[conte]->val_sp_plus() + 
	         ope_var.G_plus(l) * ope_var.operateurs[conte]->der_sp_plus() ;	

	}

	// On recupere la valeur de la sh : 
	double val_sh = cos(phases[start])*ope_var.operateurs[conte]->val_sh_one_plus()
			+ sin(phases[start])*ope_var.operateurs[conte]->val_sh_two_plus() ;
	double der_sh = cos(phases[start])*ope_var.operateurs[conte]->der_sh_one_plus()
			+ sin(phases[start])*ope_var.operateurs[conte]->der_sh_two_plus() ;

	systeme.set(taille-2, taille-1) = -ope_var.G_minus(nz-1) * val_sh  ;
	systeme.set(taille-1, taille-1) =  -ope_var.dG_minus(nz-1)*val_sh- ope_var.G_minus(nz-1)*der_sh ;

	sec_membre.set(taille-2) += ope_var.F_minus(nz-1,k,j) ;
	sec_membre.set(taille-1) += ope_var.dF_minus(nz-1,k,j) ;

	// On resout le systeme ...
	if (taille > 2)
	  systeme.set_band(2,2) ;
	else
	  systeme.set_band(1,1) ;
	
	systeme.set_lu() ;
	Tbl facteur (systeme.inverse(sec_membre)) ;

	amplis[start] = facteur(taille-1) ;

	// On range tout ca :
	// Noyau 
	nr = source.get_mg()->get_nr(0) ;
	for (int i=0 ; i<nr ; i++)
	  resultat.set(0,k,j,i) = solution_part(0,k,j,i) 
	    +facteur(0)*solution_hom_un(0,k,j,i) ;
  
	// Shells
	for (int l=1 ; l<nz-1 ; l++) {
	  nr = source.get_mg()->get_nr(l) ;
	  for (int i=0 ; i<nr ; i++)
	    resultat.set(l,k,j,i) = solution_part(l,k,j,i) + 
	      facteur(2*l-1)*solution_hom_un(l,k,j,i) +
	      facteur(2*l)*solution_hom_deux(l,k,j,i) ;
	}
      }
	start ++ ;
    }
  return resultat;
}


