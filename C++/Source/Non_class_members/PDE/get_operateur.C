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


char get_operateur_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.3  2001/10/29  10:55:28  novak
 * Error fixed for r^2 d^2/dr^2 operator
 *
 * Revision 1.2  2000/12/18 13:33:46  novak
 * *** empty log message ***
 *
 * Revision 1.1  2000/12/04 16:36:50  novak
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
#include "base_val.h"
#include "proto.h"

/*************************************************************************
 *
 * Routine used by sol_dalembert, to compute the matrix of the operator
 * to be solved. The coefficients (Ci) are stored in par.get_tbl_mod(1->9)
 * The time-step is par.get_double(0). Other inputs are:
 * l: spherical harmonic number
 * alpha, beta: coefficients of the affine mapping (see map.h)
 * Outputs are: type_dal : type of the operator (see type_parite.h)
 *              operateur: matrix of the operator
 * 
 * The operator reads: 
 * 
 *  Indentity - 0.5*dt^2 [ (C1 + C3r^2) d^2/dr^2 + (C6/r + C5r) d/dr
 *                         ({C9-l(l+1)}/r^2 + C7) Id ] 
 * 
 *************************************************************************/



		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

void _get_operateur_dal_pas_prevu(Param& , int& , double& , 
			     double& , int& , Matrice& )
{
    cout << "get_operateur_dal pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
}


void _get_operateur_dal_r_chebp(Param& par, int& l, double& alpha, 
			     double& , int& type_dal, Matrice& operateur)
{
  int nr = operateur.get_dim(0) ;
  assert (nr == operateur.get_dim(1)) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 0) ;
  assert ((par.get_tbl_mod()).get_dim(1) == 12 ) ;
  assert ((par.get_tbl_mod()).get_ndim() ==2 ) ;
  static int nap = 0 ;

  double dt = par.get_double(0) ;
  dt *= 0.5*dt ;
  double alpha2 = alpha*alpha ;

  // Copies the global coefficients to a local Tbl and adds the -l(l+1) term
  Tbl coeff(7) ;
  coeff.set_etat_qcq() ;
  coeff.set(1) = (par.get_tbl_mod())(1,0) ;
  coeff.set(2) = (par.get_tbl_mod())(3,0) ;
  coeff.set(3) = (par.get_tbl_mod())(6,0) ;
  coeff.set(4) = (par.get_tbl_mod())(5,0) ;
  coeff.set(5) = (par.get_tbl_mod())(9,0) - l*(l+1) ;
  coeff.set(6) = (par.get_tbl_mod())(7,0) ;

  //***********************************************************************
  //                Definition of the type of operator
  // For each type and a fixed time-step, if the number of points (nr) is too
  // large, the round-off error makes the matrix ill-conditioned. So one has
  // to pass the last line of the matrix to the first place (see dal_inverse).
  // The linear combinations to put the matrix into a banded form also depend
  // on the type of operator.
  //***********************************************************************

  if (fabs(coeff(1)) + fabs(coeff(2)) + fabs(coeff(5)) < 1.e-30) {
    // First order operator
    if (dt < 0.1/(fabs(coeff(3)) + fabs(coeff(4))*nr))
      type_dal = ORDRE1_SMALL ;
    else type_dal = ORDRE1_LARGE ;
  }
  else {
    // Second order degenerate (no 1/r^2 term)
    if (fabs(coeff(5)) < 1.e-24) {
      if (dt < 1./(fabs(coeff(1)) + fabs(coeff(2)) + fabs(coeff(3))*nr
		   +fabs(coeff(4)))/nr/nr) type_dal = O2DEGE_SMALL ;
      else type_dal = O2DEGE_LARGE ;
    }
    else {
      // Second order non-degenerate (most general case)
      if (dt < 1./((fabs(coeff(1)) + fabs(coeff(2)) + fabs(coeff(3))*nr
		    + fabs(coeff(4)) + fabs(coeff(5)))*nr*nr))
	type_dal = O2NOND_SMALL ;
      else type_dal = O2NOND_LARGE ;
    }
  }

  coeff.set(1) *= dt/alpha2 ;
  coeff.set(2) *= dt ;
  coeff.set(3) *= dt/alpha2 ;
  coeff.set(4) *= dt ;
  coeff.set(5) *= dt/alpha2 ;
  coeff.set(6) *= dt ;

  static Matrice Id(nr,nr) ;
  if (nap == 0) {
    Id.set_etat_qcq() ;
    for (int i=0; i<nr; i++) {
      for (int j=0; j<nr; j++) Id.set(i,j) = 0 ;
      Id.set(i,i) = 1 ;
    }
    nap = 1 ;
  }
  operateur = Id ;

  double* vect  = new double[nr] ;
  
  if (fabs(coeff(1))+fabs(coeff(2)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (nr, &vect, R_CHEBP) ;
      
      if (fabs(coeff(1)) > 1.e-30) {
	for (int j=0 ; j<=i ; j++)
	  operateur.set(j, i) -= coeff(1)*vect[j] ; 
      }
      if (fabs(coeff(2)) > 1.e-30) {
	multx2_1d (nr, &vect, R_CHEBP) ;
	for (int j=0; j<=i; j++) 
	  operateur.set(j, i) -= coeff(2)*vect[j] ;
      }
    }
  }
  
  if (fabs(coeff(3)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (nr, &vect, R_CHEBP) ;
      for (int j=0 ; j<=i ; j++)
	operateur.set(j, i) -= coeff(3)*vect[j] ; 
    }
  }
  
  if (fabs(coeff(4)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      xdsdx_1d (nr, &vect, R_CHEBP) ;
      for (int j=0 ; j<=i ; j++)
	operateur.set(j, i) -= coeff(4)*vect[j] ; 
    }
  }
  
  if (fabs(coeff(5)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      sx2_1d (nr, &vect, R_CHEBP) ;
      for (int j=0 ; j<=i ; j++)
	operateur.set(j, i) -= coeff(5)*vect[j] ; 
    }
  }
  
  if (fabs(coeff(6)) > 1.e-30) operateur = operateur - coeff(6)*Id ;
  
  delete [] vect ;

}


void _get_operateur_dal_r_chebi(Param& par, int& l, double& alpha, 
			     double& , int& type_dal, Matrice& operateur)
{
  int nr = operateur.get_dim(0) ;
  assert (nr == operateur.get_dim(1)) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 0) ;
  assert ((par.get_tbl_mod()).get_dim(1) == 12 ) ;
  assert ((par.get_tbl_mod()).get_ndim() == 2 ) ;
  static int nap = 0 ;

  double dt = par.get_double(0) ;
  dt *= 0.5*dt ;
  double alpha2 = alpha*alpha ;

  // Copies the global coefficients to a local Tbl and adds the -l(l+1) term
  Tbl coeff(7) ;
  coeff.set_etat_qcq() ;
  coeff.set(1) = (par.get_tbl_mod())(1,0) ;
  coeff.set(2) = (par.get_tbl_mod())(3,0) ;
  coeff.set(3) = (par.get_tbl_mod())(6,0) ;
  coeff.set(4) = (par.get_tbl_mod())(5,0) ;
  coeff.set(5) = (par.get_tbl_mod())(9,0) - l*(l+1) ;
  coeff.set(6) = (par.get_tbl_mod())(7,0) ;

  //***********************************************************************
  //                Definition of the type of operator
  // For each type and a fixed time-step, if the number of points (nr) is too
  // large, the round-off error makes the matrix ill-conditioned. So one has
  // to pass the last line of the matrix to the first place (see dal_inverse).
  // The linear combinations to put the matrix into a banded form also depend
  // on the type of operator.
  //***********************************************************************

  if (fabs(coeff(1)) + fabs(coeff(2)) + fabs(coeff(3)) + 
      fabs(coeff(5)) < 1.e-30) {
    // First order operator
    if (dt < 0.1/(fabs(coeff(4))*nr))
      type_dal = ORDRE1_SMALL ;
    else type_dal = ORDRE1_LARGE ;
  }
  else {
    if (fabs(coeff(5)+coeff(3)) < 1.e-6) {
    // Second order degenerate (no 1/r^2 term)
      if (dt < 1./(fabs(coeff(1)) + fabs(coeff(2)) + fabs(coeff(3))*nr
		   +fabs(coeff(4)))/nr/nr) type_dal = O2DEGE_SMALL ;
      else type_dal = O2DEGE_LARGE ;
    }
    else {
      // Second order non-degenerate (most general case)
      if (dt < 1./((fabs(coeff(1)) + fabs(coeff(2)) + fabs(coeff(3))*nr
		    + fabs(coeff(4)) + fabs(coeff(5)))*nr*nr))
	type_dal = O2NOND_SMALL ;
      else type_dal = O2NOND_LARGE ;
    }
  }

  coeff.set(1) *= dt/alpha2 ;
  coeff.set(2) *= dt ;
  coeff.set(3) *= dt/alpha2 ;
  coeff.set(4) *= dt ;
  coeff.set(5) *= dt/alpha2 ;
  coeff.set(6) *= dt ;

  static Matrice Id(nr,nr) ;
  if (nap == 0) {
    Id.set_etat_qcq() ;
    for (int i=0; i<nr; i++) {
      for (int j=0; j<nr; j++) Id.set(i,j) = 0 ;
      Id.set(i,i) = 1 ;
    }
    nap = 1 ;
  }
  operateur = Id ;

  double* vect  = new double[nr] ;
  
  if (fabs(coeff(1))+fabs(coeff(2)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (nr, &vect, R_CHEBI) ;
      
      if (fabs(coeff(1)) > 1.e-30) {
	for (int j=0 ; j<=i ; j++)
	  operateur.set(j, i) -= coeff(1)*vect[j] ; 
      }
      if (fabs(coeff(2)) > 1.e-30) {
	multx2_1d (nr, &vect, R_CHEBI) ;
	for (int j=0; j<=i; j++) 
	  operateur.set(i,j) -= coeff(2)*vect[j] ;
      }
    }
  }
  
  if (fabs(coeff(3)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (nr, &vect, R_CHEBI) ;
      for (int j=0 ; j<=i ; j++)
	operateur.set(j, i) -= coeff(3)*vect[j] ; 
    }
  }
  
  if (fabs(coeff(4)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      xdsdx_1d (nr, &vect, R_CHEBI) ;
      for (int j=0 ; j<=i ; j++)
	operateur.set(j, i) -= coeff(4)*vect[j] ; 
    }
  }
  
  if (fabs(coeff(5)) > 1.e-30) {
    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      sx2_1d (nr, &vect, R_CHEBI) ;
      for (int j=0 ; j<=i ; j++)
	operateur.set(j, i) -= coeff(5)*vect[j] ; 
    }
  }
  
  if (fabs(coeff(6)) > 1.e-30) operateur = operateur - coeff(6)*Id ;
  
  delete [] vect ;

}


		 //--------------------------
		//- La routine a appeler  ---
	       //----------------------------
void get_operateur_dal(Param& par, int& l, int& base_r, double alpha, 
		       double beta, int& type_dal, Matrice& operateur)
{

		// Routines de derivation
  static void (*get_operateur_dal[MAX_BASE])(Param&, int&, double&, double&,
					     int&, Matrice&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) 
      get_operateur_dal[i] = _get_operateur_dal_pas_prevu ;
    
    // Les routines existantes
    //get_operateur_dal[R_CHEB >> TRA_R] = _get_operateur_dal_r_cheb ;
    get_operateur_dal[R_CHEBP >> TRA_R] = _get_operateur_dal_r_chebp ;
    get_operateur_dal[R_CHEBI >> TRA_R] = _get_operateur_dal_r_chebi ;
  }
  
  get_operateur_dal[base_r](par, l, alpha, beta, type_dal,
			    operateur) ;
}
