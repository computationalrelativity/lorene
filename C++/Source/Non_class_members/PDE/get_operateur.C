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


 

/*
 * $Id$
 * $Log$
 * Revision 1.10  2016/12/05 16:18:09  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:53:28  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:16:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2008/08/27 08:51:15  jl_cornou
 * Added Jacobi(0,2) polynomials
 *
 * Revision 1.6  2005/01/27 10:19:43  j_novak
 * Now using Diff operators to build the matrices.
 *
 * Revision 1.5  2003/06/18 08:45:27  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
 * Revision 1.4  2002/01/03 15:30:28  j_novak
 * Some comments modified.
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
#include <cmath>

// Headers Lorene :
#include "param.h"
#include "base_val.h"
#include "diff.h"
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
 *                         (C9/r^2 + C7) Id ] 
 * 
 *************************************************************************/



		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

namespace Lorene {
void _get_operateur_dal_pas_prevu(const Param& , const int&, int& , Matrice& )
{
    cout << "get_operateur_dal pas prevu..." << endl ;
    abort() ;
    exit(-1) ;
}


void _get_operateur_dal_r_cheb(const Param& par, const int& lz, 
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

  // Copies the global coefficients to a local Tbl 
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
    coeff.set(9) = (par.get_tbl_mod())(9,lz) ;
    double R1 = (par.get_tbl_mod())(10,lz) ;
    double R2 = (par.get_tbl_mod())(11,lz) ;

    double a00 = 0. ; double a01 = 0. ; double  a02 = 0. ; 
    double a10 = 0. ; double a11 = 0. ; double  a12 = 0. ; 
    double a13 = 0. ; double a20 = 0. ; double  a21 = 0. ; 
    double a22 = 0. ; double a23 = 0. ; double  a24 = 0. ; 

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
    if (fabs(a00)<1.e-15) a00 = 0 ;
    if (fabs(a01)<1.e-15) a01 = 0 ;
    if (fabs(a02)<1.e-15) a02 = 0 ;
    if (fabs(a10)<1.e-15) a10 = 0 ;
    if (fabs(a11)<1.e-15) a11 = 0 ;
    if (fabs(a12)<1.e-15) a12 = 0 ;
    if (fabs(a13)<1.e-15) a13 = 0 ;
    if (fabs(a20)<1.e-15) a20 = 0 ;
    if (fabs(a21)<1.e-15) a21 = 0 ;
    if (fabs(a22)<1.e-15) a22 = 0 ;
    if (fabs(a23)<1.e-15) a23 = 0 ;
    if (fabs(a24)<1.e-15) a24 = 0 ;

    

    Diff_id id(R_CHEB, nr) ;
    if (fabs(a00)>1.e-15) {
	operateur = a00*id ;
    }
    else{
	operateur.set_etat_qcq() ;
	for (int i=0; i<nr; i++) 
	    for (int j=0; j<nr; j++)
		operateur.set(i,j) = 0. ;
    }
    Diff_mx op01(R_CHEB, nr) ; const Matrice& m01 = op01.get_matrice() ;
    Diff_mx2 op02(R_CHEB, nr) ; const Matrice& m02 = op02.get_matrice() ;
    Diff_dsdx op10(R_CHEB, nr) ; const Matrice& m10 = op10.get_matrice() ;
    Diff_xdsdx op11(R_CHEB, nr) ; const Matrice& m11 = op11.get_matrice() ;
    Diff_x2dsdx op12(R_CHEB, nr) ; const Matrice& m12 = op12.get_matrice() ;
    Diff_x3dsdx op13(R_CHEB, nr) ; const Matrice& m13 = op13.get_matrice() ;
    Diff_dsdx2 op20(R_CHEB, nr) ; const Matrice& m20 = op20.get_matrice() ;
    Diff_xdsdx2 op21(R_CHEB, nr) ; const Matrice& m21 = op21.get_matrice() ;
    Diff_x2dsdx2 op22(R_CHEB, nr) ; const Matrice& m22 = op22.get_matrice() ;
    Diff_x3dsdx2 op23(R_CHEB, nr) ; const Matrice& m23 = op23.get_matrice() ;
    Diff_x4dsdx2 op24(R_CHEB, nr) ; const Matrice& m24 = op24.get_matrice() ;

    for (int i=0; i<nr; i++) {
	int jmin = (i>3 ? i-3 : 0) ; 
	int jmax = (i<nr-9 ? i+10 : nr) ;
	for (int j=jmin ; j<jmax; j++) 
	    operateur.set(i,j) += a01*m01(i,j) + a02*m02(i,j) 
		+ a10*m10(i,j) + a11*m11(i,j) + a12*m12(i,j) 
		+ a13*m13(i,j) + a20*m20(i,j) + a21*m21(i,j)
		+ a22*m22(i,j) + a23*m23(i,j) + a24*m24(i,j) ;
	
    }
}

void _get_operateur_dal_r_chebp(const Param& par, const int& lzone, 
				int& type_dal, Matrice& operateur)
{
  assert(lzone == 0) ; // Nucleus!
  int nr = operateur.get_dim(0) ;
  assert (nr == operateur.get_dim(1)) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 0) ;
  assert ((par.get_tbl_mod()).get_dim(1) == 12 ) ;
  assert ((par.get_tbl_mod()).get_ndim() ==2 ) ;

  double dt = par.get_double(0) ;
  dt *= 0.5*dt ;

  // Copies the global coefficients to a local Tbl and adds the -l(l+1) term
  Tbl coeff(7) ;
  coeff.set_etat_qcq() ;
  coeff.set(1) = (par.get_tbl_mod())(1,lzone) ;
  if (fabs(coeff(1))<1.e-15) coeff.set(1) = 0 ;
  coeff.set(2) = (par.get_tbl_mod())(3,lzone) ;
  if (fabs(coeff(2))<1.e-15) coeff.set(2) = 0 ;
  coeff.set(3) = (par.get_tbl_mod())(6,lzone) ;
  if (fabs(coeff(3))<1.e-15) coeff.set(3) = 0 ;
  coeff.set(4) = (par.get_tbl_mod())(5,lzone) ;
  if (fabs(coeff(4))<1.e-15) coeff.set(4) = 0 ;
  coeff.set(5) = (par.get_tbl_mod())(9,lzone) ;
  if (fabs(coeff(5))<1.e-15) coeff.set(5) = 0 ;
  coeff.set(6) = (par.get_tbl_mod())(7,lzone) ;
  if (fabs(coeff(6))<1.e-15) coeff.set(6) = 0 ;
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

  Diff_id id(R_CHEBP, nr) ;
  if (fabs(1-coeff(6))>1.e-15) {
      operateur = (1-coeff(6))*id ;
  }
  else{
      operateur.set_etat_qcq() ;
      for (int i=0; i<nr; i++) 
	  for (int j=0; j<nr; j++)
	      operateur.set(i,j) = 0. ;
  }
  Diff_sx2 op02(R_CHEBP, nr) ; const Matrice& m02 = op02.get_matrice() ;
  Diff_xdsdx op11(R_CHEBP, nr) ; const Matrice& m11 = op11.get_matrice() ;
  Diff_sxdsdx op12(R_CHEBP, nr) ; const Matrice& m12 = op12.get_matrice() ;
  Diff_dsdx2 op20(R_CHEBP, nr) ; const Matrice& m20 = op20.get_matrice() ;
  Diff_x2dsdx2 op22(R_CHEBP, nr) ; const Matrice& m22 = op22.get_matrice() ;

    for (int i=0; i<nr; i++) {
	int jmin = (i>3 ? i-3 : 0) ; 
	int jmax = (i<nr-9 ? i+10 : nr) ;
	for (int j=jmin ; j<jmax; j++) 
	    operateur.set(i,j) -= coeff(1)*m20(i,j) + coeff(2)*m22(i,j)
		+ coeff(3)*m12(i,j) + coeff(4)*m11(i,j) + coeff(5)*m02(i,j) ;
    }


}


void _get_operateur_dal_r_chebi(const Param& par, const int& lzone,
				int& type_dal, Matrice& operateur)
{
  assert(lzone == 0) ; // Nucleus!
  int nr = operateur.get_dim(0) ;
  assert (nr == operateur.get_dim(1)) ;
  assert (par.get_n_double() > 0) ;
  assert (par.get_n_tbl_mod() > 0) ;
  assert ((par.get_tbl_mod()).get_dim(1) == 12 ) ;
  assert ((par.get_tbl_mod()).get_ndim() == 2 ) ;

  double dt = par.get_double(0) ;
  dt *= 0.5*dt ;

  // Copies the global coefficients to a local Tbl and adds the -l(l+1) term
  Tbl coeff(7) ;
  coeff.set_etat_qcq() ;
  coeff.set(1) = (par.get_tbl_mod())(1,lzone) ;
  if (fabs(coeff(1))<1.e-15) coeff.set(1) = 0 ;
  coeff.set(2) = (par.get_tbl_mod())(3,lzone) ;
  if (fabs(coeff(2))<1.e-15) coeff.set(2) = 0 ;
  coeff.set(3) = (par.get_tbl_mod())(6,lzone) ;
  if (fabs(coeff(3))<1.e-15) coeff.set(3) = 0 ;
  coeff.set(4) = (par.get_tbl_mod())(5,lzone) ;
  if (fabs(coeff(4))<1.e-15) coeff.set(4) = 0 ;
  coeff.set(5) = (par.get_tbl_mod())(9,lzone) ;
  if (fabs(coeff(5))<1.e-15) coeff.set(5) = 0 ;
  coeff.set(6) = (par.get_tbl_mod())(7,lzone) ;
  if (fabs(coeff(6))<1.e-15) coeff.set(6) = 0 ;
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

  Diff_id id(R_CHEBP, nr) ;
  if (fabs(1-coeff(6))>1.e-15) {
      operateur = (1-coeff(6))*id ;
  }
  else{
      operateur.set_etat_qcq() ;
      for (int i=0; i<nr; i++) 
	  for (int j=0; j<nr; j++)
	      operateur.set(i,j) = 0. ;
  }
  Diff_sx2 op02(R_CHEBI, nr) ; const Matrice& m02 = op02.get_matrice() ;
  Diff_xdsdx op11(R_CHEBI, nr) ; const Matrice& m11 = op11.get_matrice() ;
  Diff_sxdsdx op12(R_CHEBI, nr) ; const Matrice& m12 = op12.get_matrice() ;
  Diff_dsdx2 op20(R_CHEBI, nr) ; const Matrice& m20 = op20.get_matrice() ;
  Diff_x2dsdx2 op22(R_CHEBI, nr) ; const Matrice& m22 = op22.get_matrice() ;

    for (int i=0; i<nr; i++) {
	int jmin = (i>3 ? i-3 : 0) ; 
	int jmax = (i<nr-9 ? i+10 : nr) ;
	for (int j=jmin ; j<jmax; j++) 
	    operateur.set(i,j) -= coeff(1)*m20(i,j) + coeff(2)*m22(i,j)
		+ coeff(3)*m12(i,j) + coeff(4)*m11(i,j) + coeff(5)*m02(i,j) ;
    }

}



void _get_operateur_dal_r_jaco02(const Param& par, const int& lz, 
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

  // Copies the global coefficients to a local Tbl 
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
    coeff.set(9) = (par.get_tbl_mod())(9,lz) ;
    double R1 = (par.get_tbl_mod())(10,lz) ;
    double R2 = (par.get_tbl_mod())(11,lz) ;

    double a00 = 0. ; double a01 = 0. ; double  a02 = 0. ; 
    double a10 = 0. ; double a11 = 0. ; double  a12 = 0. ; 
    double a13 = 0. ; double a20 = 0. ; double  a21 = 0. ; 
    double a22 = 0. ; double a23 = 0. ; double  a24 = 0. ; 

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
    if (fabs(a00)<1.e-15) a00 = 0 ;
    if (fabs(a01)<1.e-15) a01 = 0 ;
    if (fabs(a02)<1.e-15) a02 = 0 ;
    if (fabs(a10)<1.e-15) a10 = 0 ;
    if (fabs(a11)<1.e-15) a11 = 0 ;
    if (fabs(a12)<1.e-15) a12 = 0 ;
    if (fabs(a13)<1.e-15) a13 = 0 ;
    if (fabs(a20)<1.e-15) a20 = 0 ;
    if (fabs(a21)<1.e-15) a21 = 0 ;
    if (fabs(a22)<1.e-15) a22 = 0 ;
    if (fabs(a23)<1.e-15) a23 = 0 ;
    if (fabs(a24)<1.e-15) a24 = 0 ;

    

    Diff_id id(R_JACO02, nr) ;
    if (fabs(a00)>1.e-15) {
	operateur = a00*id ;
    }
    else{
	operateur.set_etat_qcq() ;
	for (int i=0; i<nr; i++) 
	    for (int j=0; j<nr; j++)
		operateur.set(i,j) = 0. ;
    }
    Diff_mx op01(R_JACO02, nr) ; const Matrice& m01 = op01.get_matrice() ;
    Diff_mx2 op02(R_JACO02, nr) ; const Matrice& m02 = op02.get_matrice() ;
    Diff_dsdx op10(R_JACO02, nr) ; const Matrice& m10 = op10.get_matrice() ;
    Diff_xdsdx op11(R_JACO02, nr) ; const Matrice& m11 = op11.get_matrice() ;
    Diff_x2dsdx op12(R_JACO02, nr) ; const Matrice& m12 = op12.get_matrice() ;
    Diff_x3dsdx op13(R_JACO02, nr) ; const Matrice& m13 = op13.get_matrice() ;
    Diff_dsdx2 op20(R_JACO02, nr) ; const Matrice& m20 = op20.get_matrice() ;
    Diff_xdsdx2 op21(R_JACO02, nr) ; const Matrice& m21 = op21.get_matrice() ;
    Diff_x2dsdx2 op22(R_JACO02, nr) ; const Matrice& m22 = op22.get_matrice() ;
    Diff_x3dsdx2 op23(R_JACO02, nr) ; const Matrice& m23 = op23.get_matrice() ;
    Diff_x4dsdx2 op24(R_JACO02, nr) ; const Matrice& m24 = op24.get_matrice() ;

    for (int i=0; i<nr; i++) {
	int jmin = (i>3 ? i-3 : 0) ; 
	int jmax = (i<nr-9 ? i+10 : nr) ;
	for (int j=jmin ; j<jmax; j++) 
	    operateur.set(i,j) += a01*m01(i,j) + a02*m02(i,j) 
		+ a10*m10(i,j) + a11*m11(i,j) + a12*m12(i,j) 
		+ a13*m13(i,j) + a20*m20(i,j) + a21*m21(i,j)
		+ a22*m22(i,j) + a23*m23(i,j) + a24*m24(i,j) ;
	
    }
}




		 //--------------------------
		//- La routine a appeler  ---
	       //----------------------------
void get_operateur_dal(const Param& par, const int& lzone,
		       const int& base_r, int& type_dal, Matrice& operateur)
{

		// Routines de derivation
  static void (*get_operateur_dal[MAX_BASE])(const Param&,
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
    get_operateur_dal[R_JACO02 >> TRA_R] = _get_operateur_dal_r_jaco02 ;
  }
  
  get_operateur_dal[base_r](par, lzone, type_dal, operateur) ;
}
}
