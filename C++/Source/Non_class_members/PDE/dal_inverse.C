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


char dal_inverse_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2008/02/18 13:53:43  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.5  2004/10/05 15:44:21  j_novak
 * Minor speed enhancements.
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
 * Revision 1.1  2000/12/04  16:37:03  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <math.h>

//Headers LORENE
#include "param.h"
#include "matrice.h"
#include "proto.h"

/***************************************************************************
 *
 *      Set of routines to invert the dalembertian operator given
 * by get_operateur. The matrix is put into a banded form, following
 * the type of the operator (type_dal) and the odd/even decomposition 
 * base (base_r).
 * type_dal is supposed to be given by get_operateur (see the various cases
 * there and in type_parite.h).
 * part is a boolean saying wether one is looking for a particular sol.
 * of the system defined by operateur (i.e. X such as operateur*X = source),
 * part = true, or a homogeneous one (part = false).
 *
 ***************************************************************************/

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
Tbl _dal_inverse_pas_prevu (const Matrice&, const Tbl&, const bool) {
    cout << " Inversion du dalembertien pas prevue ..... : "<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB    -----
	      //-------------------

Tbl _dal_inverse_r_cheb_o2d_s(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(source) ;

  // Operator is put into banded form (changing the image base)

  int dirac = 2 ; // Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>1 ? i-1 : 0) ;
    int nrmax = (i<nr-9 ? i+10 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    if (part)
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>1 ? i-1 : 0) ;
    int nrmax = (i<nr-9 ? i+10 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }

  // 6 over-diagonals and 2 under...
  barre.set_band(5,1) ;

  // LU decomposition
  barre.set_lu() ;

  // Inversion using LAPACK

  return barre.inverse(aux) ;
}

Tbl _dal_inverse_r_cheb_o2d_l(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(source) ;

  // Operator is put into banded form (changing the image base)

  int dirac = 2 ; // Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>1 ? i-1 : 0) ;
    int nrmax = (i<nr-9 ? i+10 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    if (part)
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>1 ? i-1 : 0) ;
    int nrmax = (i<nr-9 ? i+10 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }

  // In this case the time-step is too large for the number of points
  // (or the number of points too large for the time-step!)
  // the matrix is ill-conditionned: the last lines are put as first
  // ones and all the others are shifted.

  double temp1, temp2 ;
  temp1 = aux(nr-1) ;
  temp2 = aux(nr-2) ;
  for (int i=nr-3; i>=0; i--) {
    int nrmin = (i>1 ? i-1 : 0) ;
    int nrmax = (i<nr-9 ? i+10 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i+2,j) = barre(i,j) ;
    aux.set(i+2) = aux(i) ;
  }
  aux.set(0) = temp2 ;
  aux.set(1) = temp1 ;
  
  barre.set(0,0) = 0.5 ;
  barre.set(0,1) = 1. ;
  barre.set(0,2) = 1. ;
  barre.set(0,3) = 1. ;
  barre.set(0,4) = 0. ;
  
  barre.set(1,0) = 0. ;
  barre.set(1,1) = 0.5 ;
  barre.set(1,2) = 1. ;
  barre.set(1,3) = -1. ;
  barre.set(1,4) = 1. ;
  barre.set(1,5) = 0. ;

  // 3 over diagonals and 3 under ...
  barre.set_band(3,3) ;

  // LU decomposition
  barre.set_lu() ;

  // Inversion using LAPACK

  return barre.inverse(aux) ;
}

Tbl _dal_inverse_r_cheb_o2_s(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(source) ;

  // Operator is put into banded form (changing the image base)

  int dirac = 2 ; // Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>2 ? i-2 : 0) ;
    int nrmax = (i<nr-10 ? i+11 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    if (part)
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>2 ? i-2 : 0) ;
    int nrmax = (i<nr-10 ? i+11 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }

  // 6 over-diagonals and 2 under...
  barre.set_band(6,2) ;

  // LU decomposition
  barre.set_lu() ;

  // Inversion using LAPACK

  return barre.inverse(aux) ;
}

Tbl _dal_inverse_r_cheb_o2_l(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(source) ;

  // Operator is put into banded form (changing the image base)

  int dirac = 2 ; // Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>2 ? i-2 : 0) ;
    int nrmax = (i<nr-10 ? i+11 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    if (part)
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmin = (i>2 ? i-2 : 0) ;
    int nrmax = (i<nr-10 ? i+11 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }

  // In this case the time-step is too large for the number of points
  // (or the number of points too large for the time-step!)
  // the matrix is ill-conditionned: the last lines are put as first
  // ones and all the others are shifted.

  double temp1, temp2 ;
  temp1 = aux(nr-1) ;
  temp2 = aux(nr-2) ;
  for (int i=nr-3; i>=0; i--) {
    int nrmin = (i>2 ? i-2 : 0) ;
    int nrmax = (i<nr-10 ? i+11 : nr) ;
    for (int j=nrmin; j<nrmax; j++) 
      barre.set(i+2,j) = barre(i,j) ;
    aux.set(i+2) = aux(i) ;
  }
  aux.set(0) = temp2 ;
  aux.set(1) = temp1 ;
  
  barre.set(0,0) = 0.5 ;
  barre.set(0,1) = 1. ;
  barre.set(0,2) = 1. ;
  barre.set(0,3) = 1. ;
  barre.set(0,4) = 1. ;
  barre.set(0,5) = 0. ;
  
  barre.set(1,0) = 0. ;
  barre.set(1,1) = 0.5 ;
  barre.set(1,2) = -1. ;
  barre.set(1,3) = 1. ;
  barre.set(1,4) = -1. ;
  barre.set(1,5) = 1. ;
  barre.set(1,6) = 0. ;

  // 4 over diagonals and 4 under ...
  barre.set_band(4,4) ;

  // LU decomposition
  barre.set_lu() ;

  // Inversion using LAPACK

  return barre.inverse(aux) ;
}
	
		//-------------------
	       //--  R_CHEBP   -----
	      //-------------------

Tbl _dal_inverse_r_chebp_o2d_s(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(nr) ;
  if (part) {
    aux = source ;
    aux.set(nr-1) = 0. ;
  }
  else {
    aux.annule_hard() ;
    aux.set(nr-1) = 1. ;
  }

  // Operator is put into banded form (changing the image base)

  int dirac = 2 ; // Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) 
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    if (part)
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) barre.set(i,j) = barre(i+1,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+1) - aux(i) ;
  }
  if (fabs(barre(nr-5,nr-1)) >= 1.e-16) {
    if (fabs(barre(nr-5,nr-1)) > fabs(barre(nr-2,nr-1))) {
      double lambda = barre(nr-2,nr-1)/barre(nr-5,nr-1) ; 
      for (int j=nr-5; j<nr; j++) barre.set(nr-5,j) = barre(nr-5,j)*lambda
				 -barre(nr-2,j) ;
      if (part) aux.set(nr-5) = aux(nr-5)*lambda - aux(nr-2) ;
    }
    else {
      double lambda = barre(nr-5,nr-1)/barre(nr-2,nr-1) ; 
      for (int j=nr-5; j<nr; j++) barre.set(nr-5,j) -= lambda*barre(nr-2,j) ;
      if (part) aux.set(nr-5) -= lambda*aux(nr-2) ;
    }
  }

  // 3 over-diagonals and 0 under...
  barre.set_band(3,0) ;

  // LU decomposition
  barre.set_lu() ;

  // Inversion using LAPACK

  return barre.inverse(aux) ;
}



Tbl _dal_inverse_r_chebp_o2d_l(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(nr) ;
  if (part) {
    aux = source ;
    aux.set(nr-1) = 0. ;
  }
  else {
    aux.annule_hard() ;
    aux.set(0) = 1. ;
  }

  // Operator is put into banded form (changing the image base)

  int dirac = 2 ; // Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) 
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    if (part)
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) barre.set(i,j) = barre(i+1,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+1) - aux(i) ;
  }
  if (fabs(barre(nr-5,nr-1)) >= 1.e-16) {
    if (fabs(barre(nr-5,nr-1)) > fabs(barre(nr-2,nr-1))) {
      double lambda = barre(nr-2,nr-1)/barre(nr-5,nr-1) ; 
      for (int j=nr-5; j<nr; j++) barre.set(nr-5,j) = barre(nr-5,j)*lambda
				 -barre(nr-2,j) ;
      if (part) aux.set(nr-5) = aux(nr-5)*lambda - aux(nr-2) ;
    }
    else {
      double lambda = barre(nr-5,nr-1)/barre(nr-2,nr-1) ; 
      for (int j=nr-5; j<nr; j++) barre.set(nr-5,j) -= lambda*barre(nr-2,j) ;
      if (part) aux.set(nr-5) -= lambda*aux(nr-2) ;
    }
  }
  
  // In this case the time-step is too large for the number of points
  // (or the number of points too large for the time-step!)
  // the matrix is ill-conditionned: the last line is put as the first
  // one and all the others are shifted.

  for (int i=nr-2; i>=0; i--) {
    for (int j=i; j<((i+5 > nr) ? nr : i+5); j++) 
      barre.set(i+1,j) = barre(i,j) ;
    if (part) aux.set(i+1) = aux(i) ;
  }
  barre.set(0,0) = 0.5 ;
  barre.set(0,1) = 1. ;
  barre.set(0,2) = 1. ;
  barre.set(0,3) = 0. ;
  
  if (part) aux.set(0) = 0 ;
  
  // 2 over diagonals and one under ...
  barre.set_band(2,1) ;

  // LU decomposition
  barre.set_lu() ;

  // Inversion using LAPACK

  return barre.inverse(aux);
}



Tbl _dal_inverse_r_chebp_o2_s(const Matrice &op, const Tbl &source, 
			      const bool part) {
  
  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(nr-1) ;
  if (part) {
    aux.set_etat_qcq() ;
    for (int i=nr-4; i<nr-1; i++) aux.set(i) = source(i) ; 
  }
  else {
    aux.annule_hard() ;
    aux.set(nr-2) = 1. ;
  }

  // Operator is put into banded form (changing the image base ...)

  int dirac = 2 ;// Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) 
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    if (part)
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) barre.set(i,j) = barre(i+1,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+1) - aux(i) ;
  }

  // ... and changing the starting base (the last column is quit) ... 
  
  Matrice tilde(nr-1,nr-1) ;
  tilde.set_etat_qcq() ;
  for (int i=0; i<nr-1; i++) 
    for (int j=0; j<nr-1;j++) 
      tilde.set(i,j) = barre(i,j+1) + barre(i,j) ;
  
  // 3 over-diagonals and 1 under...
  tilde.set_band(3,1) ;

  // LU decomposition
  tilde.set_lu() ;

  // Inversion using LAPACK
  Tbl res0(tilde.inverse(aux)) ;
  Tbl res(nr) ;
  res.set_etat_qcq() ;

  // ... finally, one has to recover the original starting base
  res.set(0) = res0(0) ;
  res.set(nr-1) = res0(nr-2) ;
  for (int i=1; i<nr-1; i++) res.set(i) = res0(i-1) + res0(i);
  
  return res ;
  
  
}
	
Tbl _dal_inverse_r_chebp_o2_l(const Matrice &op, const Tbl &source, 
			      const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(nr-1) ;
  if (part) {
    aux.set_etat_qcq() ;
    for (int i=nr-4; i<nr-1; i++) aux.set(i) = source(i) ;
  }
  else {
    aux.annule_hard() ;
    aux.set(0) = 1 ;
  }

  // Operator is put into banded form (changing the image base ...)

  int dirac = 2 ;// Don't forget the factor 2 for T_0!!
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) {
      barre.set(i,j) = (op(i+2,j) - dirac*op(i,j))/(i+1) ;
    }
    if (part) 
      aux.set(i) = (source(i+2) - dirac*source(i))/(i+1) ;
    if (i==0) dirac = 1 ;
  }
  for (int i=0; i<nr-4; i++) {
    int nrmax = (i<nr-7 ? i+8 : nr) ;
    for (int j=i; j<nrmax; j++) barre.set(i,j) = barre(i+1,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+1) - aux(i) ;
  }
  
  // ... and changing the starting base (the last column is quit) ... 
  
  Matrice tilde(nr-1,nr-1) ;
  tilde.set_etat_qcq() ;
  for (int i=0; i<nr-1; i++) 
    for (int j=0; j<nr-1;j++) 
      tilde.set(i,j) = barre(i,j+1) + barre(i,j) ;

  // In this case the time-step is too large for the number of points
  // (or the number of points too large for the time-step!)
  // the matrix is ill-conditionned: the last line is put as the first
  // one and all the others are shifted.

  for (int i=nr-3; i>=0; i--) {
    for (int j=((i>0) ? i-1 : 0); j<((i+5 > nr-1) ? nr-1 : i+5); j++) 
      tilde.set(i+1,j) = tilde(i,j) ;
    if (part) aux.set(i+1) = aux(i) ;
  }
  tilde.set(0,0) = 0.5 ;
  tilde.set(0,1) = 1. ;
  tilde.set(0,2) = 1. ;
  tilde.set(0,3) = 0. ;
  
  if (part) aux.set(0) = 0 ;
      
  // 2 over-diagonals and 2 under...
  tilde.set_band(2,2) ;
  
  // LU decomposition
  tilde.set_lu() ;
  
  // Inversion using LAPACK
  Tbl res0(tilde.inverse(aux)) ;
  Tbl res(nr) ;
  res.set_etat_qcq() ;

  // ... finally, one has to recover the original starting base

  res.set(0) = res0(0) ;
  res.set(nr-1) = res0(nr-2) ;
  for (int i=1; i<nr-1; i++) res.set(i) = res0(i-1) + res0(i);
  
  return res ;

    
}
		//-------------------
	       //--  R_CHEBI   -----
	      //-------------------

Tbl _dal_inverse_r_chebi_o2d_s(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  int nr = op.get_dim(0) ;
  Matrice barre(nr-1,nr-1) ;
  barre.set_etat_qcq() ;
  for (int i=0; i<nr-1; i++) 
    for (int j=0; j<nr-1; j++)
      barre.set(i,j) = op(i,j) ;
  Tbl aux(nr-1) ;
  if (part) {
    aux.set_etat_qcq() ;
    for (int i=nr-4; i<nr-1; i++) aux.set(i) = source(i) ; 
  }
  else {
    aux.annule_hard() ;
    aux.set(nr-2) = 1 ;
  }

  // Operator is put into banded form (changing the image base )
  // Last column is quit.

  for (int i=0; i<nr-4; i++) {
    for (int j=i; j<nr-1; j++) {
      barre.set(i,j) = (op(i+1,j) - op(i,j))/(i+1) ;
    }
    if (part) aux.set(i) = (source(i+1) - source(i))/(i+1) ;
  }
  for (int i=0; i<nr-5; i++) {
    for (int j=i; j<nr-1; j++) barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }
  if (fabs(barre(nr-6,nr-2)) >= 1.e-16) {
    if (fabs(barre(nr-6,nr-2)) > fabs(barre(nr-3,nr-2))) {
      double lambda = barre(nr-3,nr-2)/barre(nr-6,nr-2) ; 
      for (int j=0; j<nr-1; j++) barre.set(nr-6,j) = barre(nr-6,j)*lambda
				   -barre(nr-3,j) ;
      if (part) aux.set(nr-6) = aux(nr-6)*lambda - aux(nr-3) ;
    }
    else {
      double lambda = barre(nr-6,nr-2)/barre(nr-3,nr-2) ; 
      for (int j=0; j<nr-1; j++) barre.set(nr-6,j) -= lambda*barre(nr-3,j) ;
      if (part) aux.set(nr-6) -= lambda*aux(nr-3) ;
    }
  }

  // 3 over-diagonals and 0 under...
  barre.set_band(3,0) ;

  // LU decomposition
  barre.set_lu() ;

  // Inversion using LAPACK
  Tbl res0(barre.inverse(aux)) ;
  Tbl res(nr) ;
  res.set_etat_qcq() ;
  for (int i=0; i<nr-1; i++) res.set(i) = res0(i) ;
  res.set(nr-1) = 0 ;

  return res ;
}

Tbl _dal_inverse_r_chebi_o2d_l(const Matrice &op, const Tbl &source, 
			       const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(nr-1) ;
  if (part) {
    aux.set_etat_qcq() ;
    for (int i=0; i<nr-2; i++) aux.set(i) = source(i) ;
    aux.set(nr-2) = 0 ;
  }
  else {
    aux.annule_hard() ;
    aux.set(0) = 1 ;
  }
  // Operator is put into banded form (changing the image base)
  // Last column is quit.

  for (int i=0; i<nr-4; i++) {
    for (int j=i; j<nr-1; j++) {
      barre.set(i,j) = (op(i+1,j) - op(i,j))/(i+1) ;
    }
    if (part) aux.set(i) = (source(i+1) - source(i))/(i+1) ;
  }
  for (int i=0; i<nr-5; i++) {
    for (int j=i; j<nr-1; j++) barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }
  if (fabs(barre(nr-6,nr-2)) >= 1.e-16) {
    if (fabs(barre(nr-6,nr-2)) > fabs(barre(nr-3,nr-2))) {
      double lambda = barre(nr-3,nr-2)/barre(nr-6,nr-2) ; 
      for (int j=0; j<nr-1; j++) barre.set(nr-6,j) = barre(nr-6,j)*lambda
				   -barre(nr-3,j) ;
      if (part) aux.set(nr-6) = aux(nr-6)*lambda - aux(nr-3) ;
    }
    else {
      double lambda = barre(nr-6,nr-2)/barre(nr-3,nr-2) ; 
      for (int j=0; j<nr-1; j++) barre.set(nr-6,j) -= lambda*barre(nr-3,j) ;
      aux.set(nr-6) -= lambda*aux(nr-3) ;
    }
  }

  // In this case the time-step is too large for the number of points
  // (or the number of points too large for the time-step!)
  // the matrix is ill-conditionned: the last line is put as the first
  // one and all the others are shifted.

  Matrice tilde(nr-1,nr-1) ;
  tilde.set_etat_qcq() ;
  for (int i=nr-3; i>=0; i--) {
    for (int j=0; j<nr-1; j++) 
      tilde.set(i+1,j) = barre(i,j) ;
    if (part) aux.set(i+1) = aux(i) ;
  }
  tilde.set(0,0) = 1. ;
  tilde.set(0,1) = 1. ;
  tilde.set(0,2) = 1. ;
  tilde.set(0,3) = 0. ;
  
  if (part) aux.set(0) = 0 ;
      

  // 2 over-diagonals and 1 under...
  tilde.set_band(2,1) ;

  // LU decomposition
  tilde.set_lu() ;

  // Inversion using LAPACK
  Tbl res0(tilde.inverse(aux)) ;
  Tbl res(nr) ;
  res.set_etat_qcq() ;
  for (int i=0; i<nr-1; i++) res.set(i) = res0(i) ;
  res.set(nr-1) = 0 ;

  return res ;    
}
	
Tbl _dal_inverse_r_chebi_o2_s(const Matrice &op, const Tbl &source, 
			      const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(nr-2) ;
  if (part) {
    aux.set_etat_qcq() ;
    aux.set(nr-4) = source(nr-4) ;
    aux.set(nr-3) = 0 ;
  }
  else {
    aux.annule_hard() ;
    aux.set(nr-3) = 1. ;
  }

  // Operator is put into banded form (changing the image base ...)

  for (int i=0; i<nr-4; i++) {
    for (int j=i; j<nr; j++) {
      barre.set(i,j) = (op(i+1,j) - op(i,j))/(i+1) ;
    }
    if (part)
      aux.set(i) = (source(i+1) - source(i))/(i+1) ;
  }
  for (int i=0; i<nr-5; i++) {
    for (int j=i; j<nr; j++) barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }
  
  // ... and changing the starting base (first and last columns are quit)... 
  
  Matrice tilde(nr-2,nr-2) ;
  tilde.set_etat_qcq() ;
  for (int i=0; i<nr-2; i++) 
    for (int j=0; j<nr-2;j++) 
      tilde.set(i,j) = barre(i,j+1)*(2*j+1) + barre(i,j)*(2*j+3) ;
  
  // 3 over-diagonals and 1 under...
  tilde.set_band(3,1) ;
  
  // LU decomposition
  tilde.set_lu() ;

  // Inversion using LAPACK
  Tbl res0(tilde.inverse(aux)) ;
  Tbl res(nr) ;
  res.set_etat_qcq() ;

  // ... finally, one has to recover the original starting base

  res.set(0) = 3*res0(0) ;
  for (int i=1; i<nr-2; i++) res.set(i) = res0(i-1)*(2*i-1) 
			       + res0(i)*(2*i+3) ;
  res.set(nr-2) = res0(nr-3)*(2*nr-5) ;
  res.set(nr-1) = 0 ;

  return res ;
}

Tbl _dal_inverse_r_chebi_o2_l(const Matrice &op, const Tbl &source, 
			      const bool part) {

  // Operator and source are copied and prepared
  Matrice barre(op) ;
  int nr = op.get_dim(0) ;
  Tbl aux(nr-2) ;
  if (part) {
    aux.set_etat_qcq() ;
    aux.set(nr-4) = source(nr-4) ;
    aux.set(nr-3) = 0 ;
  }
  else {
    aux.annule_hard() ;
    aux.set(0) = 1. ;
  }

  // Operator is put into banded form (changing the image base ...)

  for (int i=0; i<nr-4; i++) {
    for (int j=i; j<nr; j++) {
      barre.set(i,j) = (op(i+1,j) - op(i,j))/(i+1) ;
    }
    if (part)
      aux.set(i) = (source(i+1) - source(i))/(i+1) ;
    }
  for (int i=0; i<nr-5; i++) {
    for (int j=i; j<nr; j++) barre.set(i,j) = barre(i+2,j) - barre(i,j) ;
    if (part) aux.set(i) = aux(i+2) - aux(i) ;
  }
    
  // ... and changing the starting base (first and last columns are quit) ... 
  
  Matrice tilde(nr-2,nr-2) ;
  tilde.set_etat_qcq() ;
  for (int i=0; i<nr-2; i++) 
    for (int j=0; j<nr-2;j++) 
      tilde.set(i,j) = barre(i,j+1)*(2*j+1) + barre(i,j)*(2*j+3) ;
  
  // In this case the time-step is too large for the number of points
  // (or the number of points too large for the time-step!)
  // the matrix is ill-conditionned: the last line is put as the first
  // one and all the others are shifted.

  for (int i=nr-4; i>=0; i--) {
    for (int j=((i>0) ? i-1 : 0); j<((i+5 > nr-2) ? nr-2 : i+5); j++) 
      tilde.set(i+1,j) = tilde(i,j) ;
    if (part) aux.set(i+1) = aux(i) ;
  }
  tilde.set(0,0) = 1. ;
  tilde.set(0,1) = 1. ;
  tilde.set(0,2) = 1. ;
  tilde.set(0,3) = 0. ;
  
  if (part) aux.set(0) = 0 ;
      

  // 2 over-diagonals and 2 under...
  tilde.set_band(2,2) ;

  // LU decomposition
  tilde.set_lu() ;
  
  // Inversion using LAPACK
  Tbl res0(tilde.inverse(aux)) ;
  Tbl res(nr) ;
  res.set_etat_qcq() ;
  // ... finally, one has to recover the original starting base
  res.set(0) = 3*res0(0) ;
  for (int i=1; i<nr-2; i++) res.set(i) = res0(i-1)*(2*i-1) 
			       + res0(i)*(2*i+3) ;
  res.set(nr-2) = res0(nr-3)*(2*nr-5) ;
  res.set(nr-1) = 0 ;
  
  return res ;
}


	      	//----------------------------
	       //--  Fonction a appeler   ---
	      //----------------------------
	      
	      
Tbl dal_inverse(const int& base_r, const int& type_dal, const 
		Matrice& operateur, const Tbl& source, const bool part) {

		// Routines de derivation
    static Tbl (*dal_inverse[MAX_BASE][MAX_DAL])(const Matrice&, const Tbl&, 
						 const bool) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_DAL ; i++) {
	  for (int j=0; j<MAX_BASE; j++) 
	    dal_inverse[i][j] = _dal_inverse_pas_prevu ;
	}
		// Les routines existantes
//  	dal_inverse[R_CHEB >> TRA_R] = _dal_inverse_r_cheb ; not good!!
//  	dal_inverse[ORDRE1_SMALL][R_CHEB >> TRA_R] = 
//  	  _dal_inverse_r_cheb_o1_s ;
//  	dal_inverse[ORDRE1_LARGE][R_CHEB >> TRA_R] = 
//  	  _dal_inverse_r_cheb_o1_l ;
	dal_inverse[O2DEGE_SMALL][R_CHEB >> TRA_R] = 
	  _dal_inverse_r_cheb_o2d_s ;
  	dal_inverse[O2DEGE_LARGE][R_CHEB >> TRA_R] = 
  	  _dal_inverse_r_cheb_o2d_l ;
  	dal_inverse[O2NOND_SMALL][R_CHEB >> TRA_R] = 
  	  _dal_inverse_r_cheb_o2_s ;
  	dal_inverse[O2NOND_LARGE][R_CHEB >> TRA_R] = 
  	  _dal_inverse_r_cheb_o2_l ;
//  	dal_inverse[R_CHEB >> TRA_R] = _dal_inverse_r_cheb ; not good!!
//  	dal_inverse[ORDRE1_SMALL][R_CHEBP >> TRA_R] = 
//  	  _dal_inverse_r_chebp_o1_s ;
//  	dal_inverse[ORDRE1_LARGE][R_CHEBP >> TRA_R] = 
//  	  _dal_inverse_r_chebp_o1_l ;
	dal_inverse[O2DEGE_SMALL][R_CHEBP >> TRA_R] = 
	  _dal_inverse_r_chebp_o2d_s ;
  	dal_inverse[O2DEGE_LARGE][R_CHEBP >> TRA_R] = 
  	  _dal_inverse_r_chebp_o2d_l ;
  	dal_inverse[O2NOND_SMALL][R_CHEBP >> TRA_R] = 
  	  _dal_inverse_r_chebp_o2_s ;
  	dal_inverse[O2NOND_LARGE][R_CHEBP >> TRA_R] = 
  	  _dal_inverse_r_chebp_o2_l ;
//  	dal_inverse[ORDRE1_SMALL][R_CHEBI >> TRA_R] = 
//  	  _dal_inverse_r_chebi_o1_s ;
//  	dal_inverse[ORDRE1_LARGE][R_CHEBI >> TRA_R] = 
//  	  _dal_inverse_r_chebi_o1_l ;
  	dal_inverse[O2DEGE_SMALL][R_CHEBI >> TRA_R] = 
  	  _dal_inverse_r_chebi_o2d_s ;
  	dal_inverse[O2DEGE_LARGE][R_CHEBI >> TRA_R] = 
  	  _dal_inverse_r_chebi_o2d_l ;
  	dal_inverse[O2NOND_SMALL][R_CHEBI >> TRA_R] = 
  	  _dal_inverse_r_chebi_o2_s ;
  	dal_inverse[O2NOND_LARGE][R_CHEBI >> TRA_R] = 
  	  _dal_inverse_r_chebi_o2_l ;
    }
    
    return dal_inverse[type_dal][base_r](operateur,  source, part) ;
}
