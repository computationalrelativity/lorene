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
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
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
  assert (nz == 1) ; // for the moment, only the nucleus!
  assert (source.get_mg()->get_type_r(0) == RARE) ;
  //for (int l=1 ; l<nz-1 ; l++)
  //  assert(source.get_mg()->get_type_r(l) == FIN) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 1) ;
  
    
  // Bases spectrales
  const Base_val& base = source.base ;
  
  // donnees sur la zone
  int nr, nt, np ;
  int base_r, type_dal ;
  double alpha, beta ;
  int l_quant, m_quant;
  
  //Rangement des valeurs intermediaires 
  Tbl *so ;
  Tbl *sol_hom ;
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
  
  // nbre maximum de point en theta et en phi : inutile pur l'instant
  //int np_max, nt_max ;
  
  
  //---------------
  //--  NUCLEUS ---
  //---------------
  
  nr = source.get_mg()->get_nr(0) ;
  nt = source.get_mg()->get_nt(0) ;
  np = source.get_mg()->get_np(0) ;
  
  //nt_max = nt ;
  //np_max = np ;
  
  alpha = mapping.get_alpha()[0] ;
  beta = mapping.get_beta()[0] ;
  
  for (int k=0 ; k<np+1 ; k++)
    for (int j=0 ; j<nt ; j++) 
      if (nullite_plm(j, nt, k, np, base) == 1) 
	{
	  // quantic numbers and spectral bases
	  donne_lm(nz, 0, j, k, base, m_quant, l_quant, base_r) ;

	  Matrice operateur(nr,nr) ;
	  
	  get_operateur_dal(par, l_quant, base_r, alpha, beta, 
			    type_dal, operateur) ;

	  // Getting the particular solution
	  so = new Tbl(nr) ;
	  so->set_etat_qcq() ;
	  for (int i=0 ; i<nr ; i++)
	    so->set(i) = source(0, k, j, i) ;
	  if ((type_dal == ORDRE1_LARGE) || (type_dal == O2DEGE_LARGE)
	      || (type_dal == O2NOND_LARGE))
	    so->set(nr-1) = 0 ;
	  sol_part = new Tbl(dal_inverse(base_r, type_dal, alpha, beta,
					 operateur, *so, true)) ;
	  
	  // Getting the homogeneous solution
	  sol_hom = new Tbl(dal_inverse(base_r, type_dal, alpha, beta,
					operateur, *so, false)) ;
	  
	  
	  // Putting to Mtbl_cf
	  
	  for (int i=0 ; i<nr ; i++) {
	    solution_part.set(0, k, j, i) = (*sol_part)(i) ;
	    solution_hom_un.set(0, k, j, i) = (*sol_hom)(i) ;
	    solution_hom_deux.set(0, k, j, i) = 0. ; 
	  }
	  
	  // If only one zone, the BC is set
	  if (nz == 1) {

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
	    double lambda = ((*tbc3)(k,j) - part) ;
	    if (fabs(lambda)>1.e-33) {
	      lambda/= hom ;
	      for (int i=0 ; i<nr ; i++)
		resultat.set(0, k, j, i) = 
		  solution_part(0, k, j, i)
		  +lambda*solution_hom_un(0, k, j, i) ; 
	    }
	  }
	  
	  delete so ;
	  delete sol_hom ;
	  delete sol_part ;
	}

  if (nz == 1) return resultat ;
  return resultat ;

}
