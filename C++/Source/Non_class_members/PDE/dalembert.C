/*
 *   Copyright (c) 2000-2001 Jerome Novak
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


char dalembert_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2004/03/01 09:57:03  j_novak
 * the wave equation is solved with Scalars. It now accepts a grid with a
 * compactified external domain, which the solver ignores and where it copies
 * the values of the field from one time-step to the next.
 *
 * Revision 1.6  2003/12/19 16:21:49  j_novak
 * Shadow hunt
 *
 * Revision 1.5  2003/07/25 08:31:20  j_novak
 * Error corrected in the case of only nucleus
 *
 * Revision 1.4  2003/06/18 08:45:27  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
 * Revision 1.3  2002/01/03 13:18:41  j_novak
 * Optimization: the members set(i,j) and operator(i,j) of class Matrice are
 * now defined inline. Matrice is a friend class of Tbl.
 *
 * Revision 1.2  2002/01/02 14:07:57  j_novak
 * Dalembert equation is now solved in the shells. However, the number of
 * points in theta and phi must be the same in each domain. The solver is not
 * completely tested (beta version!).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2000/12/04  14:24:15  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Header C : 
#include <math.h>

// Headers Lorene :
#include "param.h"
#include "matrice.h"
#include "map.h"
#include "base_val.h"
#include "proto.h"


	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------

/*
 * 
 * Solution de l'equation de d'Alembert
 * 
 * Entree : mapping :   le mapping affine
 *	    source : les coefficients de la source 
 *		    La base de decomposition doit etre Ylm
 * Sortie : renvoie les coefficients de la solution dans la meme base de 
 *	    decomposition (a savoir Ylm)
 *	    
 */



Mtbl_cf sol_dalembert(Param& par, const Map_af& mapping, const Mtbl_cf& source) 
{
    
  // Verifications d'usage sur les zones
  int nz = source.get_mg()->get_nzone() ;
  bool ced = (source.get_mg()->get_type_r(nz-1) == UNSURR ) ;
  int nz0 = (ced ? nz - 1 : nz ) ;
  assert (source.get_mg()->get_type_r(0) == RARE) ;
  for (int l=1 ; l<nz0 ; l++) {
    assert(source.get_mg()->get_type_r(l) == FIN) ;
    assert(source.get_mg()->get_nt(l) == source.get_mg()->get_nt(0)) ;
    assert(source.get_mg()->get_np(l) == source.get_mg()->get_np(0)) ;
  } // Same number of points in theta and phi in all domains...
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 1) ;
  
    
  // Bases spectrales
  const Base_val& base = source.base ;
  
  // donnees sur la zone
  int nr, nt, np ;
  int base_r, type_dal ;
  double alpha, beta ;
  int l_quant, m_quant;
  nt = source.get_mg()->get_nt(0) ;
  np = source.get_mg()->get_np(0) ;
  
  //Rangement des valeurs intermediaires 
  Tbl *so ;
  Tbl *sol_hom ;
  Tbl *sol_hom2 ;
  Tbl *sol_part ;
  
  // Rangement des solutions, avant raccordement
  Mtbl_cf solution_part(source.get_mg(), base) ;
  Mtbl_cf solution_hom_un(source.get_mg(), base) ;
  Mtbl_cf solution_hom_deux(source.get_mg(), base) ;
  Mtbl_cf resultat(source.get_mg(), base) ;
  
  solution_part.set_etat_qcq() ;
  solution_hom_un.set_etat_qcq() ;
  solution_hom_deux.set_etat_qcq() ;
  resultat.set_etat_qcq() ;

  // Tbls for the boundary condition
  double* bc1 = &par.get_double_mod(1) ;
  double* bc2 = &par.get_double_mod(2) ;
  Tbl* tbc3 = &par.get_tbl_mod(1) ;
  
  for (int l=0 ; l<nz ; l++) {
    solution_part.t[l]->set_etat_qcq() ;
    solution_hom_un.t[l]->set_etat_qcq() ;
    solution_hom_deux.t[l]->set_etat_qcq() ;
    resultat.t[l]->set_etat_qcq() ;
    for (int k=0 ; k<source.get_mg()->get_np(l)+2 ; k++)
      for (int j=0 ; j<source.get_mg()->get_nt(l) ; j++)
	for (int i=0 ; i<source.get_mg()->get_nr(l) ; i++) {
	  resultat.set(l, k, j, i) = 0 ;
	  solution_part.set(l,k,j,i) = 0 ;
	  solution_hom_un.set(l,k,j,i) = 0 ;
	  solution_hom_deux.set(l,k,j,i) = 0 ;
	}
  }
  
  int lz = 0 ;
  
  
  //---------------
  //--  NUCLEUS ---
  //---------------
  
  nr = source.get_mg()->get_nr(lz) ;
  
  alpha = mapping.get_alpha()[lz] ;
  
  for (int k=0 ; k<np+1 ; k++) {
    for (int j=0 ; j<nt ; j++) {
      if (nullite_plm(j, nt, k, np, base) == 1) 
	{
	  // quantic numbers and spectral bases
	  donne_lm(nz, lz, j, k, base, m_quant, l_quant, base_r) ;

	  //Calculation of the coefficients of the operator
	  par.get_tbl_mod().set(4,lz) = 2*par.get_tbl_mod()(2,lz) ;
	  par.get_tbl_mod().set(5,lz) = 2*par.get_tbl_mod()(3,lz) ;
	  par.get_tbl_mod().set(6,lz) = 2*par.get_tbl_mod()(1,lz) ;
	  par.get_tbl_mod().set(7,lz) = 
	    -l_quant*(l_quant+1)*par.get_tbl_mod()(3,lz) ;
	  par.get_tbl_mod().set(8,lz) = 
	    -l_quant*(l_quant+1)*par.get_tbl_mod()(2,lz) ;
	  par.get_tbl_mod().set(9,lz) = 
	    -l_quant*(l_quant+1)*par.get_tbl_mod()(1,lz) ;
	  
	  Matrice operateur(nr,nr) ;
	  
	  get_operateur_dal(par, lz, base_r, type_dal, operateur) ;

	  // Getting the particular solution
	  so = new Tbl(nr) ;
	  so->set_etat_qcq() ;
	  for (int i=0 ; i<nr ; i++)
	    so->set(i) = source(lz, k, j, i) ;
	  if ((type_dal == ORDRE1_LARGE) || (type_dal == O2DEGE_LARGE)
	      || (type_dal == O2NOND_LARGE))
	    so->set(nr-1) = 0 ;
	  sol_part = new Tbl(dal_inverse(base_r, type_dal, operateur, 
					 *so, true)) ;
	  
	  // Getting the homogeneous solution
	  sol_hom = new Tbl(dal_inverse(base_r, type_dal, operateur, 
					*so, false)) ;
	  
	  
	  // Putting to Mtbl_cf
	  
	  for (int i=0 ; i<nr ; i++) {
	    solution_part.set(lz, k, j, i) = (*sol_part)(i) ;
	    solution_hom_un.set(lz, k, j, i) = (*sol_hom)(i) ;
	    solution_hom_deux.set(lz, k, j, i) = 0. ; 
	  }
	  
	  // If only one zone, the BC is set
	  if (nz0 == 1) {

	    int base_pipo = 0 ;
	    double part, dpart, hom, dhom;
	    Tbl der_part(3,1,nr) ;
	    der_part.set_etat_qcq() ;
	    for (int i=0; i<nr; i++) 
	      der_part.set(0,0,i) = (*sol_part)(i) ;
	    Tbl der_hom(3,1,nr) ;
	    der_hom.set_etat_qcq() ;
	    for (int i=0; i<nr; i++) 
	      der_hom.set(0,0,i) = (*sol_hom)(i) ;

	    if (base_r == R_CHEBP) {
	      som_r_chebp(sol_part->t, nr, 1, 1, 1., &part) ;
	      _dsdx_r_chebp(&der_part, base_pipo) ;
	      som_r_chebi(der_part.t, nr, 1, 1, 1., &dpart) ;
	      som_r_chebp(sol_hom->t, nr, 1, 1, 1., &hom) ;
	      _dsdx_r_chebp(&der_hom, base_pipo) ;
	      som_r_chebi(der_hom.t, nr, 1, 1, 1., &dhom) ;
	    }
	    else {
	      som_r_chebi(sol_part->t, nr, 1, 1, 1., &part) ;
	      _dsdx_r_chebi(&der_part, base_pipo) ;
	      som_r_chebp(der_part.t, nr, 1, 1, 1., &dpart) ;
	      som_r_chebi(sol_hom->t, nr, 1, 1, 1., &hom) ;
	      _dsdx_r_chebi(&der_hom, base_pipo) ;
	      som_r_chebp(der_hom.t, nr, 1, 1, 1., &dhom) ;
	    }
	    
	    part = part*(*bc1) + dpart*(*bc2)/alpha ;
	    hom = hom*(*bc1) + dhom*(*bc2)/alpha ;
	    double lambda = ((*tbc3)(k,j) - part) / hom ;
	    for (int i=0 ; i<nr ; i++)
	      resultat.set(lz, k, j, i) = 
		solution_part(lz, k, j, i)
		+lambda*solution_hom_un(lz, k, j, i) ; 
	  }
	  
	  delete so ;
	  delete sol_hom ;
	  delete sol_part ;
	}
    }
  }  

  //---------------------
  //--      SHELLS     --
  //---------------------
  for (lz=1 ; lz<nz0 ; lz++) {
    nr = source.get_mg()->get_nr(lz) ;
    
    alpha = mapping.get_alpha()[lz] ;
    beta = mapping.get_beta()[lz] ;
    
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) 
	if (nullite_plm(j, nt, k, np, base) == 1)
	  {
	    // calcul des nombres quantiques :
	    donne_lm(nz, lz, j, k, base, m_quant, l_quant, base_r) ;
	    
	    //Calculation of the coefficients of the operator
	    par.get_tbl_mod().set(4,lz) = 2*par.get_tbl_mod()(2,lz) ;
	    par.get_tbl_mod().set(5,lz) = 2*par.get_tbl_mod()(3,lz) ;
	    par.get_tbl_mod().set(6,lz) = 2*par.get_tbl_mod()(1,lz) ;
	    par.get_tbl_mod().set(7,lz) = 
	      -l_quant*(l_quant+1)*par.get_tbl_mod()(3,lz) ;
	    par.get_tbl_mod().set(8,lz) = 
	      -l_quant*(l_quant+1)*par.get_tbl_mod()(2,lz) ;
	    par.get_tbl_mod().set(9,lz) = 
	      -l_quant*(l_quant+1)*par.get_tbl_mod()(1,lz) ;
	  
	    Matrice operateur(nr,nr) ;

	    get_operateur_dal(par, lz, base_r, type_dal, operateur) ;
	    
	    // Calcul DES DEUX SH
	    so = new Tbl(nr) ;
	    so->set_etat_qcq() ;
	    for (int i=0; i<nr; i++) so->set(i) = 0. ;
	    so->set(nr-2) = 1. ;
	    sol_hom = new Tbl(dal_inverse(base_r, type_dal, operateur, *so,
					  false)) ;
	    so->set(nr-2) = 0. ;
	    so->set(nr-1) = 1. ;
	    sol_hom2 = new Tbl(dal_inverse(base_r, type_dal, operateur, *so,
					  false)) ;
	    
	    // Calcul de la SP
	    double *tmp = new double[nr] ;
	    for (int i=0 ; i<nr ; i++)
	      tmp[i] = source(lz, k, j, i) ;
	    if ((type_dal == O2DEGE_SMALL) || (type_dal == O2DEGE_LARGE)) {
	      for (int i=0; i<nr; i++) so->set(i) = beta*tmp[i] ;
	      multx_1d(nr, &tmp, R_CHEB) ;
	      for (int i=0; i<nr; i++) so->set(i) += alpha*tmp[i] ;
	    }
	    else {
	      for (int i=0; i<nr; i++) so->set(i) = beta*beta*tmp[i] ;
	      multx_1d(nr, &tmp, R_CHEB) ;
	      for (int i=0; i<nr; i++) so->set(i) += 2*alpha*beta*tmp[i] ;
	      multx_1d(nr, &tmp, R_CHEB) ;
	      for (int i=0; i<nr; i++) so->set(i) += alpha*alpha*tmp[i] ;
	    }
	    so->set(nr-2) = 0. ;
	    so->set(nr-1) = 0. ;
	    
	    sol_part = new Tbl (dal_inverse(base_r, type_dal, operateur, 
					    *so, true)) ;

	    // Rangement
	    for (int i=0 ; i<nr ; i++) {
	      solution_part.set(lz, k, j, i) = (*sol_part)(i) ;
	      solution_hom_un.set(lz, k, j, i) = (*sol_hom)(i) ;
	      solution_hom_deux.set(lz, k, j, i) = (*sol_hom2)(i) ;
	    }
	    
	    delete [] tmp ;
	    delete so ;
	    delete sol_hom ;
	    delete sol_hom2 ;
	    delete sol_part ;
	  }
  }
  if (nz0 > 1) {
    //--------------------------------------------------------------------
    //
    //     Combinations of particular and homogeneous solutions
    //         to verify continuity and boundary conditions
    //
    //--------------------------------------------------------------------
    int taille = 2*nz0 - 1 ; 
    Tbl deuz(taille) ;
    deuz.set_etat_qcq() ;
    Matrice systeme(taille,taille) ;
    systeme.set_etat_qcq() ;
    int sup = 2;
    int inf = (nz0>2) ? 2 : 1 ;
    for (int k=0; k<np+1; k++) {
      for (int j=0; j<nt; j++) {
	if (nullite_plm(j, nt, k, np, base)) {
	  // To get the r basis in the nucleus
	  donne_lm (nz, 0, j, k, base, m_quant, l_quant, base_r) ;
	  assert ((base_r == R_CHEBP)||(base_r == R_CHEBI)) ;
	  int parite = (base_r == R_CHEBP) ? 0 : 1 ;
	  int l ; 
	  int c ;
	  double xx = 0.;
	  for (l=0; l<taille; l++) 
	    for (c=0; c<taille; c++) systeme.set(l,c) = xx ;
	  for (l=0; l<taille; l++) deuz.set(l) = xx ;
	  
	  //---------
	  // Nucleus
	  //---------
	  nr = source.get_mg()->get_nr(0) ;
	  alpha = mapping.get_alpha()[0] ;
	  l=0 ; c=0 ;
	  for (int i=0; i<nr; i++) 
	    systeme.set(l,c) += solution_hom_un(0, k, j, i) ;

	  for (int i=0; i<nr; i++) deuz.set(l) -= solution_part(0, k, j, i) ;

	  l++ ;
	  xx = 0. ;
	  for (int i=0; i<nr; i++)
	      xx +=(2*i+parite)*(2*i+parite)
		*solution_hom_un(0, k, j, i) ;
	  systeme.set(l,c) += xx/alpha ;
	  xx = 0. ; 
	  for (int i=0; i<nr; i++) xx -= (2*i+parite)*
				     (2*i+parite)*solution_part(0, k, j, i) ;
	  deuz.set(l) += xx/alpha ;
	  
	  //----------
	  //  Shells
	  //----------
	  for (lz=1; lz<nz0; lz++) {
	    nr = source.get_mg()->get_nr(lz) ;
	    alpha = mapping.get_alpha()[lz] ;
	    l-- ; 
	    c = l+1 ;
	    for (int i=0; i<nr; i++) 
	      if (i%2 == 0)
		systeme.set(l,c) -= solution_hom_un(lz, k, j, i) ;
	      else
		systeme.set(l,c) += solution_hom_un(lz, k, j, i) ;
	    c++ ;
	    for (int i=0; i<nr; i++) 
	      if (i%2 == 0)
		systeme.set(l,c) -= solution_hom_deux(lz, k, j, i) ;
	      else
		systeme.set(l,c) += solution_hom_deux(lz, k, j, i) ;
	    for (int i=0; i<nr; i++) 
	      if (i%2 == 0) deuz.set(l) += solution_part(lz, k, j, i) ;
	      else deuz.set(l) -= solution_part(lz, k, j, i) ;
	    
	    l++ ; c-- ;
	    xx = 0. ;
	    for (int i=0; i<nr; i++) 
	      if (i%2 == 0)
		xx += i*i*solution_hom_un(lz, k, j, i) ;
	      else
		xx -= i*i*solution_hom_un(lz, k, j, i) ;
	    systeme.set(l,c) += xx/alpha ;
	    c++ ;
	    xx = 0. ;
	    for (int i=0; i<nr; i++) 
	      if (i%2 == 0)
		xx += i*i*solution_hom_deux(lz, k, j, i) ;
	      else
		xx -= i*i*solution_hom_deux(lz, k, j, i) ;
	    systeme.set(l,c) += xx/alpha ;
	    xx = 0. ;
	    for (int i=0; i<nr; i++) 
	      if (i%2 == 0) xx -= i*i*solution_part(lz, k, j, i) ;
	      else xx += i*i*solution_part(lz, k, j, i) ;
	    deuz.set(l) += xx/alpha ;

	    l++ ; c--;
	    if (lz == nz0-1) { // Last domain, the outer BC is set
	      for (int i=0; i<nr; i++) systeme.set(l,c) +=
		((*bc1)+(*bc2)*i*i/alpha)*solution_hom_un(lz, k, j, i) ;
	      c++ ;
	      for (int i=0; i<nr; i++) systeme.set(l,c) +=
		((*bc1)+(*bc2)*i*i/alpha)*solution_hom_deux(lz, k, j, i) ;
	      for (int i=0; i<nr; i++) deuz.set(l) -=
		((*bc1)+(*bc2)*i*i/alpha)*solution_part(lz, k, j, i) ;
	      deuz.set(l) += (*tbc3)(k,j) ;
	    }
	    else { // At least one more shell
	      for (int i=0; i<nr; i++) 
		systeme.set(l,c) += solution_hom_un(lz, k, j, i) ;
	      c++ ;
	      for (int i=0; i<nr; i++) 
		systeme.set(l,c) += solution_hom_deux(lz, k, j, i) ;
	      for (int i=0; i<nr; i++) 
		deuz.set(l) -= solution_part(lz, k, j, i) ;
	      l++ ; c-- ;
	      xx = 0. ;
	      for (int i=0; i<nr; i++) xx += i*i*solution_hom_un(lz, k, j, i) ;
	      systeme.set(l,c) += xx/alpha ;
	      c++ ;
	      xx = 0. ;
	      for (int i=0; i<nr; i++) 
		xx += i*i*solution_hom_deux(lz, k, j, i) ;
	      systeme.set(l,c) += xx/alpha ;
	      xx = 0. ;
	      for (int i=0; i<nr; i++) 
		xx -= i*i*solution_part(lz, k, j, i) ;	 
	      deuz.set(l) += xx/alpha ;
	    }
	  }

	  //--------------------------------------
	  //   Solution of the linear system
	  //--------------------------------------
	  
	  systeme.set_band(sup, inf) ;
	  systeme.set_lu() ;
	  Tbl facteur(systeme.inverse(deuz)) ;

	  //Linear Combination in the nucleus
	  nr = source.get_mg()->get_nr(0) ;
	  for (int i=0; i<nr; i++) 
	    resultat.set(0, k, j, i) = solution_part(0, k, j, i) 
	      + facteur(0)*solution_hom_un(0, k, j, i) ;
	  
	  //Linear combination in the shells
	  for (lz=1; lz<nz0; lz++) {
	    nr = source.get_mg()->get_nr(lz) ;
	    for (int i=0; i<nr; i++) 
	      resultat.set(lz, k, j, i) = solution_part(lz, k, j, i) 
		+ facteur(2*lz-1)*solution_hom_un(lz, k, j, i) 
		+ facteur(2*lz)*solution_hom_deux(lz, k, j, i) ;
	  }
	}
      } //End of j/theta loop   
    } //End of k/phi loop 
  } //End of case nz0>1
  
  return resultat ;

}
