/*
 *   Copyright (c) 2003 Philippe Grandclement
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

char sol_elliptic_only_zec_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * 
 * $Header $
 *
 */

// Header C : 
#include <stdlib.h>
#include <math.h>

// Headers Lorene :
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
#include "param_elliptic.h"
          
 
	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------



Mtbl_cf elliptic_solver_only_zec  (const Param_elliptic& ope_var, const Mtbl_cf& source, double val) {
  // Verifications d'usage sur les zones
  int nz = source.get_mg()->get_nzone() ;
  assert (nz>1) ;
  assert (source.get_mg()->get_type_r(nz-1) == UNSURR) ;
  
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

  // Computation of the SP and SH's in the ZEC ...
  int conte_start = 0 ;
  for (int l=0 ; l<nz-1 ; l++) {
    nt = source.get_mg()->get_nt(l) ;
    np = source.get_mg()->get_np(l) ;
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++)
	conte_start ++ ;
  }
  int conte = conte_start ;

  int zone = nz-1 ;

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

 
  //---------------------------------------------------------
  // ON impose la bonne CL ... CASE ONLY SPHERICAL RIGHT NOW
  //---------------------------------------------------------
  
  // C'est pas simple toute cette sombre affaire...
  // Que le cas meme nombre de points dans chaque domaines...

  int start = conte_start ;
  for (int k=0 ; k<1 ; k++)
    for (int j=0 ; j<1 ; j++) {
      if (ope_var.operateurs[start] != 0x0) {
	

	// Valeurs en -1 :
	double rlim = ope_var.get_mp().val_r_jk (zone, -1, j,k);
	
	double facteur = ((val-ope_var.variables[start]->val_F(rlim))/
			  ope_var.variables[start]->val_G(rlim) - 
			  ope_var.operateurs[start]->val_sp_minus()) /
	  ope_var.operateurs[start]->val_sh_one_minus()  ; 

	// Zec
	nr = source.get_mg()->get_nr(nz-1) ;
	for (int i=0 ; i<nr ; i++)
	  resultat.set(nz-1,k,j,i) = solution_part(nz-1,k,j,i) + 
	    facteur*solution_hom_un(nz-1,k,j,i) ;
      }
	start ++ ;
    }
    
  return resultat;
}

