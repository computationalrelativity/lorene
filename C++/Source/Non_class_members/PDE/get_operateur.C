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
 * Revision 1.2  2002/01/02 14:07:57  j_novak
 * Dalembert equation is now solved in the shells. However, the number of
 * points in theta and phi must be the same in each domain. The solver is not
 * completely tested (beta version!).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
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

void _get_operateur_dal_pas_prevu(const Param& , const int&, const int&, int& , Matrice& )
{
    cout << "get_operateur_dal pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
}


void _get_operateur_dal_r_cheb(const Param& par, const int& lz, const int& l, 
int& type_dal, Matrice& operateur)
{
  int nr = operateur.get_dim(0) ;
  assert (nr == operateur.get_dim(1)) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 0) ;
  assert ((par.get_tbl_mod()).get_dim(1) == 12 ) ;
  assert ((par.get_tbl_mod()).get_ndim() ==2 ) ;

  double dt = par.get_double(0) ;
  dt *= 0.5*dt ;

  // Copies the global coefficients to a local Tbl and adds the -l(l+1) term
    Tbl coeff(10) ;
    coeff.set_etat_qcq() ;
    coeff.set(1) = (par.get_tbl_mod())(1,lz) ;
    coeff.set(2) = (par.get_tbl_mod())(2,lz) ;
    coeff.set(3) = (par.get_tbl_mod())(3,lz) ;
    coeff.set(4) = (par.get_tbl_mod())(4,lz) ;
    coeff.set(5) = (par.get_tbl_mod())(5,lz) ;
    coeff.set(6) = (par.get_tbl_mod())(6,lz) ;
    coeff.set(7) = (par.get_tbl_mod())(7,lz) ;
    coeff.set(8) = (par.get_tbl_mod())(8,lz) ;
    coeff.set(9) = (par.get_tbl_mod())(9,lz) - l*(l+1) ;
    double R1 = (par.get_tbl_mod())(10,lz) ;
    double R2 = (par.get_tbl_mod())(11,lz) ;

    double a00, a01, a02, a10, a11, a12, a13, a20, a21, a22, a23, a24 ;
    bool dege = (fabs(coeff(9)) < 1.e-10) ;
    switch (dege) {
    case true:
      a00 = R1 - dt*(coeff(7)*R1 + coeff(8)) ;
      a01 = R2 - dt*R2*coeff(7) ;
      a02 = 0 ;
      a10 = -dt*(R1*coeff(4) + R1*R1*coeff(5) + coeff(6))/R2 ;
      a11 = -dt*(coeff(4) + 2*R1*coeff(5)) ;
      a12 = -dt*R2*coeff(5) ;
      a13 = 0 ;
      a20 = -dt*R1/(R2*R2)*(coeff(1) + R1*coeff(2) + R1*R1*coeff(3)) ;
      a21 = -dt/R2*(coeff(1) + 2*R1*coeff(2) + 3*R1*R1*coeff(3)) ;
      a22 = -dt*(coeff(2) + 3*R1*coeff(3)) ;
      a23 = -dt*R2*coeff(3) ;
      a24 = 0 ;
      type_dal = ((0.1*(fabs(a20)+fabs(a21)+fabs(a22)+fabs(a23))*nr*nr*nr 
		   < 1.) ? O2DEGE_SMALL : O2DEGE_LARGE ) ;
      break ;
    case false:
      a00 = R1*R1 - dt*(coeff(7)*R1*R1 + coeff(8)*R1 + coeff(9)) ;
      a01 = 2*R1*R2 - dt*(2*R1*R2*coeff(7) + R2*coeff(8)) ;
      a02 = R2*R2*(1 - dt*coeff(7)) ;
      a10 = -dt*R1/R2*(R1*coeff(4) + R1*R1*coeff(5) + coeff(6)) ;
      a11 = -dt*(2*R1*coeff(4) + 3*R1*R1*coeff(5) + coeff(6)) ;
      a12 = -dt*(R2*coeff(4) + 3*R1*R2*coeff(5)) ;
      a13 = -dt*R2*R2*coeff(5) ;
      a20 = -dt*(R1*R1)/(R2*R2)*(coeff(1) + R1*coeff(2) + R1*R1*coeff(3)) ;
      a21 = -dt*R1/R2*(2*coeff(1) + 3*R1*coeff(2) + 4*R1*R1*coeff(3)) ;
      a22 = -dt*(coeff(1) + 3*R1*coeff(2) + 6*R1*R1*coeff(3)) ;
      a23 = -dt*(R2*coeff(2) + 4*R1*R2*coeff(3)) ;
      a24 = -dt*R2*R2*coeff(3) ;
      type_dal = ((0.1*(fabs(a20)+fabs(a21)+fabs(a22)+fabs(a23)+fabs(a24))
		   *nr*nr*nr < 1.) ? O2NOND_SMALL : O2NOND_LARGE ) ;
      break ;
    }

    double* vect = new double[nr] ;
    operateur.set_etat_qcq() ;

    for (int i=0; i<nr; i++) {
      for (int j=0; j<nr; j++) operateur.set(i,j) = 0 ;
      operateur.set(i,i) = a00 ;
    }

    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      multx_1d (nr, &vect, R_CHEB) ;
      
      for (int j=0 ; j<nr ; j++)
	  operateur.set(j, i) += a01*vect[j] ; 
      if (!dege) { 
	multx_1d (nr, &vect, R_CHEB) ;
	for (int j=0 ; j<nr ; j++)
	  operateur.set(j, i) += a02*vect[j] ; 
      }
    }

    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (nr, &vect, R_CHEB) ;
      
      for (int j=0 ; j<nr ; j++)
	  operateur.set(j, i) += a10*vect[j] ; 
      multx_1d (nr, &vect, R_CHEB) ;
      
      for (int j=0 ; j<nr ; j++)
	  operateur.set(j, i) += a11*vect[j] ; 
      multx_1d (nr, &vect, R_CHEB) ;
      
      for (int j=0 ; j<nr ; j++)
	  operateur.set(j, i) += a12*vect[j] ; 
      if (!dege) {
	multx_1d (nr, &vect, R_CHEB) ;
	for (int j=0 ; j<nr ; j++)
	  operateur.set(j, i) += a13*vect[j] ; 
      }
    }

    for (int i=0 ; i<nr ; i++) {
      for (int j=0 ; j<nr ; j++) vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (nr, &vect, R_CHEB) ;
      
      for (int j=0 ; j<nr ; j++)
	  operateur.set(j, i) += a20*vect[j] ; 
      multx_1d (nr, &vect, R_CHEB) ;
      for (int j=0; j<nr; j++) 
	  operateur.set(j, i) += a21*vect[j] ;
      multx_1d (nr, &vect, R_CHEB) ;
      for (int j=0; j<nr; j++) 
	  operateur.set(j, i) += a22*vect[j] ;
      multx_1d (nr, &vect, R_CHEB) ;
      for (int j=0; j<nr; j++) 
	  operateur.set(j, i) += a23*vect[j] ;
      if (!dege) {
	multx_1d (nr, &vect, R_CHEB) ;
	for (int j=0; j<nr; j++) 
	  operateur.set(j, i) += a24*vect[j] ;
      }
    }

  delete [] vect ;
  
}

void _get_operateur_dal_r_chebp(const Param& par, const int& lzone, const int& l, int& type_dal, 
				Matrice& operateur)
{
  assert(lzone == 0) ; // Nucleus!
  int nr = operateur.get_dim(0) ;
  assert (nr == operateur.get_dim(1)) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 0) ;
  assert ((par.get_tbl_mod()).get_dim(1) == 12 ) ;
  assert ((par.get_tbl_mod()).get_ndim() ==2 ) ;
  static int nap = 0 ;

  double dt = par.get_double(0) ;
  dt *= 0.5*dt ;

  // Copies the global coefficients to a local Tbl and adds the -l(l+1) term
  Tbl coeff(7) ;
  coeff.set_etat_qcq() ;
  coeff.set(1) = (par.get_tbl_mod())(1,lzone) ;
  coeff.set(2) = (par.get_tbl_mod())(3,lzone) ;
  coeff.set(3) = (par.get_tbl_mod())(6,lzone) ;
  coeff.set(4) = (par.get_tbl_mod())(5,lzone) ;
  coeff.set(5) = (par.get_tbl_mod())(9,lzone) - l*(l+1) ;
  coeff.set(6) = (par.get_tbl_mod())(7,lzone) ;
  double alpha2 = (par.get_tbl_mod())(11,lzone)*(par.get_tbl_mod())(11,lzone) ;

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


void _get_operateur_dal_r_chebi(const Param& par, const int& lzone, const int& l, int& type_dal, 
				Matrice& operateur)
{
  assert(lzone == 0) ; // Nucleus!
  int nr = operateur.get_dim(0) ;
  assert (nr == operateur.get_dim(1)) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 0) ;
  assert ((par.get_tbl_mod()).get_dim(1) == 12 ) ;
  assert ((par.get_tbl_mod()).get_ndim() == 2 ) ;
  static int nap = 0 ;

  double dt = par.get_double(0) ;
  dt *= 0.5*dt ;

  // Copies the global coefficients to a local Tbl and adds the -l(l+1) term
  Tbl coeff(7) ;
  coeff.set_etat_qcq() ;
  coeff.set(1) = (par.get_tbl_mod())(1,lzone) ;
  coeff.set(2) = (par.get_tbl_mod())(3,lzone) ;
  coeff.set(3) = (par.get_tbl_mod())(6,lzone) ;
  coeff.set(4) = (par.get_tbl_mod())(5,lzone) ;
  coeff.set(5) = (par.get_tbl_mod())(9,lzone) - l*(l+1) ;
  coeff.set(6) = (par.get_tbl_mod())(7,lzone) ;
  double alpha2 = (par.get_tbl_mod())(11,lzone)*(par.get_tbl_mod())(11,lzone) ;

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
void get_operateur_dal(const Param& par, const int& lzone, const int& l, 
		       const int& base_r, int& type_dal, Matrice& operateur)
{

		// Routines de derivation
  static void (*get_operateur_dal[MAX_BASE])(const Param&, const int&, 
					     const int&, int&, Matrice&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) 
      get_operateur_dal[i] = _get_operateur_dal_pas_prevu ;
    
    // Les routines existantes
    get_operateur_dal[R_CHEB >> TRA_R] = _get_operateur_dal_r_cheb ;
    get_operateur_dal[R_CHEBP >> TRA_R] = _get_operateur_dal_r_chebp ;
    get_operateur_dal[R_CHEBI >> TRA_R] = _get_operateur_dal_r_chebi ;
  }
  
  get_operateur_dal[base_r](par, lzone, l, type_dal, operateur) ;
}
