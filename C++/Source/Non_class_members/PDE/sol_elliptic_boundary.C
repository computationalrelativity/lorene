/*
 *   Copyright (c) 2003 Philippe Grandclement
 *                      Jose Luis Jaramillo
 *                      Francois Limousin
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
 * Revision 1.5  2016/12/05 16:18:10  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:30  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:16:10  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/08/20 15:03:55  n_vasset
 * Correction on how the boundary condition is imposed
 *
 * Revision 1.1  2005/06/09 08:01:05  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * 
 * $Header$
 *
 */

// Header C : 
#include <cstdlib>
#include <cmath>

// Headers Lorene :
#include "param_elliptic.h"
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
          
 
	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------



namespace Lorene {
Mtbl_cf elliptic_solver_boundary  (const Param_elliptic& ope_var, const Mtbl_cf& source, 
				   const Mtbl_cf& bound, double fact_dir, double fact_neu )  {
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

  // Computation of the SP and SH's in every domain ...
  int conte = 0 ;
  for (int zone=0 ; zone<nz ; zone++) {

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
  // ON EST PARTI POUR LE RACCORD (Be careful  ....)
  //-------------------------------------------------
  
  // C'est pas simple toute cette sombre affaire...
  // Que de cas meme nombre de points dans chaque domaines...

  int start = 0 ;
  for (int k=0 ; k<source.get_mg()->get_np(0)+1 ; k++)
    for (int j=0 ; j<source.get_mg()->get_nt(0) ; j++) {
      if (ope_var.operateurs[start] != 0x0) {
	
	int taille = 2*nz - 3 ;
	Matrice systeme (taille, taille) ;
	systeme.set_etat_qcq() ;
	for (int i=0 ; i<taille ; i++)
	  for (int j2=0 ; j2<taille ; j2++)
	    systeme.set(i,j2) = 0 ;
	Tbl sec_membre (taille) ;
	sec_membre.set_etat_qcq() ;
	for (int i=0 ; i<taille ; i++)
	  sec_membre.set(i) = 0 ;
	
	//----------
	//  BOUNDARY :
	//----------
	conte = start ;

	// Boundary value at an angular point
	
	// Setting the right hand side term

	sec_membre.set(0) -= bound.val_in_bound_jk(1, j, k)/sqrt(2.) ; //Relevant value is in the first shell and only in the first coefficient


	//--------------
	// FIRST SHELL :
	//--------------
	
	int l_1=1 ;
	
	// On se met au bon endroit :
	int np_prec_1 = source.get_mg()->get_np(l_1-1) ;
	int nt_prec_1 = source.get_mg()->get_nt(l_1-1) ;
	conte += (np_prec_1+1)*nt_prec_1 ;

 	systeme.set(0, 0) = fact_dir * (-ope_var.G_minus(l_1) * 
 					ope_var.operateurs[conte]->val_sh_one_minus() ) 
	    + fact_neu * 
	 ( -ope_var.dG_minus(l_1)*ope_var.operateurs[conte]->val_sh_one_minus()-  
	  ope_var.G_minus(l_1)*ope_var.operateurs[conte]->der_sh_one_minus() );
	  
	systeme.set(0, 1) = fact_dir * (- ope_var.G_minus(l_1) * 
	  ope_var.operateurs[conte]->val_sh_two_minus() ) + fact_neu *
	  (-ope_var.dG_minus(l_1)*ope_var.operateurs[conte]->val_sh_two_minus()-  
	  ope_var.G_minus(l_1)*ope_var.operateurs[conte]->der_sh_two_minus() ) ;

	//Completing the right hand side term
 	sec_membre.set(0) += fact_dir * (ope_var.F_minus(l_1,k,j) + 
 	  ope_var.G_minus(l_1) * ope_var.operateurs[conte]->val_sp_minus() ) +
 	  fact_neu * ( ope_var.dF_minus(l_1,k,j) + 
 	  ope_var.dG_minus(l_1) * ope_var.operateurs[conte]->val_sp_minus() + 
 	  ope_var.G_minus(l_1) * ope_var.operateurs[conte]->der_sp_minus() )   ;

	
	// Valeurs en +1 :
	systeme.set(2*l_1-1, 2*l_1-2) = ope_var.G_plus(l_1) * 
	  ope_var.operateurs[conte]->val_sh_one_plus() ;
	systeme.set(2*l_1-1, 2*l_1-1) = ope_var.G_plus(l_1) * 
	  ope_var.operateurs[conte]->val_sh_two_plus() ;
	systeme.set(2*l_1, 2*l_1-2) = 
	  ope_var.dG_plus(l_1)*ope_var.operateurs[conte]->val_sh_one_plus()+  
	  ope_var.G_plus(l_1)*ope_var.operateurs[conte]->der_sh_one_plus() ;
	systeme.set(2*l_1, 2*l_1-1) =
	    ope_var.dG_plus(l_1)*ope_var.operateurs[conte]->val_sh_two_plus()+
	  ope_var.G_plus(l_1)*ope_var.operateurs[conte]->der_sh_two_plus() ;
	
	  sec_membre.set(2*l_1-1) -=  ope_var.F_plus(l_1,k,j) + 
	    ope_var.G_plus(l_1) * ope_var.operateurs[conte]->val_sp_plus();
	  sec_membre.set(2*l_1) -=  ope_var.dF_plus(l_1,k,j) + 
	    ope_var.dG_plus(l_1) * ope_var.operateurs[conte]->val_sp_plus() + 
	    ope_var.G_plus(l_1) * ope_var.operateurs[conte]->der_sp_plus() ;
	

	
	//----------
	// SHELLS :
	//----------
//	assert (nz-1 > 2) ;    // At least two shells
	for (int l=2 ; l<nz-1 ; l++) {
	  
	  // On se met au bon endroit :
	  int np_prec = source.get_mg()->get_np(l-1) ;
	  int nt_prec = source.get_mg()->get_nt(l-1) ;
	  conte += (np_prec+1)*nt_prec ;
	  	 
	  systeme.set(2*l-3, 2*l-2) = -ope_var.G_minus(l) * 
	    ope_var.operateurs[conte]->val_sh_one_minus() ;
	  systeme.set(2*l-3, 2*l-1) = - ope_var.G_minus(l) * 
	    ope_var.operateurs[conte]->val_sh_two_minus() ;
	  systeme.set(2*l-2, 2*l-2) = 
	    -ope_var.dG_minus(l)*ope_var.operateurs[conte]->val_sh_one_minus()-  
	    ope_var.G_minus(l)*ope_var.operateurs[conte]->der_sh_one_minus() ;
	  systeme.set(2*l-2, 2*l-1) =
	    -ope_var.dG_minus(l)*ope_var.operateurs[conte]->val_sh_two_minus()-  
	    ope_var.G_minus(l)*ope_var.operateurs[conte]->der_sh_two_minus() ;
	  
	  sec_membre.set(2*l-3) += ope_var.F_minus(l,k,j) + 
	    ope_var.G_minus(l) * ope_var.operateurs[conte]->val_sp_minus() ;
	  sec_membre.set(2*l-2) += ope_var.dF_minus(l,k,j) + 
	    ope_var.dG_minus(l) * ope_var.operateurs[conte]->val_sp_minus() + 
	    ope_var.G_minus(l) * ope_var.operateurs[conte]->der_sp_minus() ;
	  
	  // Valeurs en +1 :
	  systeme.set(2*l-1, 2*l-2) = ope_var.G_plus(l) * 
	    ope_var.operateurs[conte]->val_sh_one_plus() ;
	  systeme.set(2*l-1, 2*l-1) = ope_var.G_plus(l) * 
	    ope_var.operateurs[conte]->val_sh_two_plus() ;
	  systeme.set(2*l, 2*l-2) = 
	    ope_var.dG_plus(l)*ope_var.operateurs[conte]->val_sh_one_plus()+  
	    ope_var.G_plus(l)*ope_var.operateurs[conte]->der_sh_one_plus() ;
	  systeme.set(2*l, 2*l-1) =
	    ope_var.dG_plus(l)*ope_var.operateurs[conte]->val_sh_two_plus()+
	    ope_var.G_plus(l)*ope_var.operateurs[conte]->der_sh_two_plus() ;
	  
	  sec_membre.set(2*l-1) -=  ope_var.F_plus(l,k,j) + 
	    ope_var.G_plus(l) * ope_var.operateurs[conte]->val_sp_plus();
	  sec_membre.set(2*l) -=  ope_var.dF_plus(l,k,j) + 
	    ope_var.dG_plus(l) * ope_var.operateurs[conte]->val_sp_plus() + 
	    ope_var.G_plus(l) * ope_var.operateurs[conte]->der_sp_plus() ;
	}
	
	//-------
	// ZEC :
	//-------
	int np_prec = source.get_mg()->get_np(nz-2) ;
	int nt_prec = source.get_mg()->get_nt(nz-2) ;
	conte += (np_prec+1)*nt_prec ;
	
	systeme.set(taille-2, taille-1) = -ope_var.G_minus(nz-1) * 
	  ope_var.operateurs[conte]->val_sh_one_minus() ;
	systeme.set(taille-1, taille-1) = 
	  -ope_var.dG_minus(nz-1)*ope_var.operateurs[conte]->val_sh_one_minus()-  
	  ope_var.G_minus(nz-1)*ope_var.operateurs[conte]->der_sh_one_minus()  ;
	
	sec_membre.set(taille-2) += ope_var.F_minus(nz-1,k,j) + 
	  ope_var.G_minus(nz-1)*ope_var.operateurs[conte]->val_sp_minus() ;
	sec_membre.set(taille-1) += ope_var.dF_minus(nz-1,k,j) + 
	  ope_var.dG_minus(nz-1) * ope_var.operateurs[conte]->val_sp_minus() + 
	  ope_var.G_minus(nz-1) * ope_var.operateurs[conte]->der_sp_minus() ;
	
	// On resout le systeme ...
	if (taille > 2)
	  systeme.set_band(2,2) ;
	else
	  systeme.set_band(1,1) ;
	
	systeme.set_lu() ;
	Tbl facteur (systeme.inverse(sec_membre)) ;

	// On range tout ca :
  
	// Shells
	for (int l=1 ; l<nz-1 ; l++) {
	  nr = source.get_mg()->get_nr(l) ;
	  for (int i=0 ; i<nr ; i++)
	    resultat.set(l,k,j,i) = solution_part(l,k,j,i) + 
	      facteur(2*l-2)*solution_hom_un(l,k,j,i) +
	      facteur(2*l-1)*solution_hom_deux(l,k,j,i) ;
	}
	
	// Zec
	nr = source.get_mg()->get_nr(nz-1) ;
	for (int i=0 ; i<nr ; i++)
	  resultat.set(nz-1,k,j,i) = solution_part(nz-1,k,j,i) + 
	    facteur(taille-1)*solution_hom_un(nz-1,k,j,i) ;
      }
	start ++ ;
    }
    
  return resultat;
}

}
