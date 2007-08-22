/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose-Luis Jaramillo
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


char binhor_hh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2007/08/22 16:10:35  f_limousin
 * Correction of many errors in binhor_hh.C
 *
 * Revision 1.1  2007/04/18 14:27:19  f_limousin
 * First version
 *
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>

// Lorene
#include "tensor.h"
#include "cmp.h"
#include "isol_hor.h"
#include "graphique.h"
#include "utilitaires.h"



Sym_tensor Bin_hor::hh_Samaya_hole1 () {

  //========
  //  Grid 1
  //========
  int nz1 = hole1.mp.get_mg()->get_nzone() ;
  int nz2 = hole2.mp.get_mg()->get_nzone() ;
  
  // General coordinate values
  const Coord& xx_1 = hole1.mp.x ; 
  const Coord& yy_1 = hole1.mp.y ;   
  const Coord& zz_1 = hole1.mp.z ; 

  //========
  //  Grid 2
  //========
  
  // General coordinate values
  const Coord& xx_2 = hole2.mp.x ; 
  const Coord& yy_2 = hole2.mp.y ;   
  const Coord& zz_2 = hole2.mp.z ; 

  //===================================
  // Definition of the relevant vectors
  //===================================

  // Coordinate vector from hole 1 in the grid 1: nn1
  //--------------------------------------------------
  Vector rr1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  rr1.set(1) = xx_1 ;
  rr1.set(2) = yy_1 ;
  rr1.set(3) = zz_1 ;
  rr1.std_spectral_base() ;  

  // Norm r1
  Scalar r1 (hole1.mp) ;
  r1 = hole1.mp.r ;
  r1.std_spectral_base() ;
  Scalar temp1 (r1) ;
  temp1.raccord(1) ;
  r1.set_domain(0) = temp1.domain(0) ;

  // Unitary vector
  Vector nn1 (rr1);
  nn1 = nn1/r1 ;

  for (int i=0; i<hole1.mp.get_mg()->get_nr(nz1-1); i++)
    for (int j=0; j<hole1.mp.get_mg()->get_nt(nz1-1); j++)
      for (int k=0; k<hole1.mp.get_mg()->get_np(nz1-1); k++)
	for (int l=0; l<nz1; l++)
	  for (int ind=1; ind<=3; ind++){
	    if (i != 0 || l != 1)
	      nn1.set(ind).set_grid_point(l,k,j,i) = nn1(ind).val_grid_point(1,k,j,0) ;
	  }

  //cout << "nn1(1)" << endl << nn1(1) << endl ;
  //cout << "nn1(2)" << endl << nn1(2) << endl ;
  //cout << "nn1(3)" << endl << nn1(3) << endl ;
  //cout << "r1" << endl << r1 << endl ;
      

  // Coordinate vector from hole 2 in the grid 2: nn2_2
  //-----------------------------------------------------
  Vector rr2_2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  rr2_2.set(1) = xx_2 ;
  rr2_2.set(2) = yy_2 ;
  rr2_2.set(3) = zz_2 ;
  rr2_2.std_spectral_base() ;  

  // Norm r2_g2 
  Scalar r2_2 (hole2.mp) ;
  r2_2 = hole2.mp.r ;
  r2_2.std_spectral_base() ;
  Scalar temp2 (r2_2) ;
  temp2.raccord(1) ;
  r2_2.set_domain(0) = temp2.domain(0) ;
   
  // Unitary vector
  Vector nn2_2 (rr2_2);
  nn2_2 = nn2_2/r2_2 ;
  
  for (int i=0; i<hole2.mp.get_mg()->get_nr(nz2-1); i++)
    for (int j=0; j<hole2.mp.get_mg()->get_nt(nz2-1); j++)
      for (int k=0; k<hole2.mp.get_mg()->get_np(nz2-1); k++)
	for (int l=0; l<nz2; l++)
	  for (int ind=1; ind<=3; ind++){
	    if (i != 0 || l != 1)
	      nn2_2.set(ind).set_grid_point(l,k,j,i) = nn2_2(ind).val_grid_point(1,k,j,0) ;
	  }

  Scalar unsr1 (hole1.mp) ;
  unsr1 = 1./hole1.mp.r ;
  unsr1.std_spectral_base() ;
  unsr1.raccord(1) ;
  
  Scalar unsr1_2 (hole2.mp) ;
  unsr1_2.set_etat_qcq() ;
  unsr1_2.import(unsr1) ;
  unsr1_2.set_spectral_va().set_base(unsr1.get_spectral_va().get_base()) ;
  
  Scalar r2sr1_2 (hole2.mp) ;
  r2sr1_2 = r2_2*unsr1_2 ;
  r2sr1_2.set_outer_boundary(nz1-1, 1.) ;

 
  //des_meridian(unsr1, 0., 20., "unsr1", 15) ;
  //des_meridian(unsr1_2, 0., 20., "unsr1_2", 16) ;
  //des_meridian(r2_2, 0., 20., "r2_2", 19) ;
  //des_meridian(r2sr1_2, 0., 20., "r2sr1_2", 18) ;
  

  // Coordinate vector from hole 2 in the grid 1: nn2
  //-----------------------------------------------------
  Vector nn2 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  nn2_2.change_triad(hole1.mp.get_bvect_cart()) ;
  nn2.set_etat_qcq() ;
  for (int i=1 ; i<=3 ; i++){ 
    nn2.set(i).import(nn2_2(i)) ;
    nn2.set(i).set_spectral_va().set_base(nn2_2(i).get_spectral_va().get_base()) ;
  }  
  
  Scalar r2sr1 (hole1.mp) ;
  r2sr1.set_etat_qcq() ;
  r2sr1.import(r2sr1_2) ;
  r2sr1.set_spectral_va().set_base(r2sr1_2.get_spectral_va().get_base()) ;
  
  //des_meridian(r2sr1, 0., 20., "r2sr1", 13) ;
  //des_meridian(r2sr1_g2, 0., 20., "r2sr1_g2", 14) ;
  //arrete() ;

  
  // Coordinate vector from hole 1 to hole 2 in the grid 1: nn12
  //----------------------------------------------------------------
  // Warning! Valid only in the symmetric case (for the general case it would
  // necessary to construct this whole function as a Bin_hor function 
  Vector rr12 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  rr12.set(1) = hole1.mp.get_ori_x() - hole2.mp.get_ori_x() ;
  rr12.set(2) = hole1.mp.get_ori_y() - hole2.mp.get_ori_y() ;
  rr12.set(3) = hole1.mp.get_ori_z() - hole2.mp.get_ori_z() ;
  rr12.std_spectral_base() ;  

  //Norm r12
  Scalar r12 (hole1.mp) ;
  r12 = sqrt( rr12(1)*rr12(1) + rr12(2)*rr12(2) + rr12(3)*rr12(3)) ;
  r12.std_spectral_base() ;

  // Unitary vector
  Vector nn12 ( rr12 );
  nn12 = nn12/ r12 ;


  Scalar f_delta (hole1.mp) ;
  Scalar f_delta_zec (hole1.mp) ;
  Scalar f_1_1 (hole1.mp) ;
  Scalar f_1_1_zec (hole1.mp) ;
  Scalar f_1_12 (hole1.mp) ;
  Scalar f_1_12_zec (hole1.mp) ;
  Scalar f_12_12 (hole1.mp) ;
  Scalar f_1_2 (hole1.mp) ;

  f_delta.set_etat_qcq() ;
  f_delta_zec.set_etat_qcq() ;
  f_1_1.set_etat_qcq() ;
  f_1_1_zec.set_etat_qcq() ;
  f_1_12.set_etat_qcq() ;
  f_1_12_zec.set_etat_qcq() ;
  f_12_12.set_etat_qcq() ;
  f_1_2.set_etat_qcq() ;  
 
  // Function exp(-(r-r_0)^2/sigma^2)
  // --------------------------------
  
  double r0 = hole1.mp.val_r(nz1-2, 1, 0, 0) ;
  double sigma = 1.*r0 ;
  
  Scalar rr (hole1.mp) ;
  rr = hole1.mp.r ;

  Scalar fexp (hole1.mp) ;
  fexp = exp( -(rr - r0)*(rr - r0)/sigma/sigma ) ;
  for (int ii=0; ii<nz1-1; ii++)
    fexp.set_domain(ii) = 1. ;
  fexp.set_outer_boundary(nz1-1, 0) ;
  fexp.std_spectral_base() ;
 
  // Conformal metric
  //=================

  // tilde{gamma}- \delta = m_1*m_2* ( f_delta \delta_{ij} 
  //                        + f_1_1 nn1*nn1 + f_1_12 nn1*nn12
  //                        + f_12_12 nn12*nn12
  //                        + f_2_2 nn2*nn2 + f_2_12 nn2*nn12
  //                        + f_1_2 nn1*nn2
 
  // f_delta
  //--------
  f_delta = -5.*r1/(8.*r12*r12*r12) - 15./(8.*r1*r12) + 
    5.*r1/(8.*r12*r12*r12*r2sr1) + 1./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1))*
    (1 + r1/r12 + r12/r1 - 1./r2sr1 - r1/(r2sr1*r12) + r12*r12/(2*r1*r1*r2sr1)) +
    1./(r1*(1.+r2sr1+r12/r1))*(-7./r1 + 2./r12) ;  

  f_delta.annule_domain(nz1-1) ;
  
  f_delta_zec = - 15./(8.*r1*r12) + 1./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1))*
    (1 + r1/r12 + r12/r1 - 1./r2sr1 - r1/(r2sr1*r12) + r12*r12/(2*r1*r1*r2sr1)) +
    1./(r1*(1.+r2sr1+r12/r1))*(-7./r1 + 2./r12) ; 
  f_delta_zec += fexp*(-5.*r1/(8.*r12*r12*r12) + 5.*r1/(8.*r12*r12*r12*r2sr1)) ;

  f_delta_zec.set_outer_boundary(nz1-1, 0.) ;
  for (int i=0 ;i<nz1-1 ; i++){
    f_delta_zec.annule_domain(i) ;
  }
  
  f_delta = f_delta + f_delta_zec ;
  
  //des_meridian(f_delta, 0., 20., "f_delta", 13) ;

  // f_1_1
  //------
  f_1_1 = r1/(8.*r12*r12*r12) + 11./(8.*r1*r12) -
    r2sr1*r2sr1*r1/(8.*r12*r12*r12) + 7./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1)) +
    7./(r1*r1*(1.+r2sr1+r12/r1)) ;

  f_1_1.annule_domain(nz1-1) ;
  
  f_1_1_zec = 11./(8.*r1*r12) + 7./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1)) +
    7./(r1*r1*(1.+r2sr1+r12/r1)) ;
  f_1_1_zec += fexp*(r1/(8.*r12*r12*r12) - r2sr1*r2sr1*r1/(8.*r12*r12*r12)) ;

  f_1_1_zec.set_outer_boundary(nz1-1, 0.) ;
  for (int i=0 ; i<nz1-1 ; i++){
    f_1_1_zec.annule_domain(i) ;
  }
  
  f_1_1 = f_1_1 + f_1_1_zec ;

  //des_meridian(f_1_1, 0., 20., "f_1_1", 14) ;

  // f_1_12
  //------
  f_1_12 = - 7./(2*r12*r12) + 8./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1)) ; // ????????
  f_1_12.annule_domain(nz1-1) ;
  
  f_1_12_zec = 8./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1)) ;
  f_1_12_zec.set_outer_boundary(nz1-1, 0.) ;

  for (int i=0 ; i<nz1-1 ; i++){
    f_1_12_zec.annule_domain(i) ;
  }
   
  f_1_12 = f_1_12 + f_1_12_zec ;
  
  //des_meridian(f_1_12, 0., 20., "f_1_12", 15) ;
 
  // f_12_12
  //-------
  f_12_12 = 2* (-4./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1)) -  // 2  ???????
		4./(r12*r1*(1.+r2sr1+r12/r1))) ;
  f_12_12.set_outer_boundary(nz1-1, 0.) ;

  //des_meridian(f_12_12, 0., 20., "f_12_12", 16) ;

  // f_1_2
  //-------
  f_1_2 = 2 * 11./(r1*r1*(1.+r2sr1+r12/r1)*(1.+r2sr1+r12/r1));        // 2  ???????
  f_1_2.set_outer_boundary(nz1-1, 0.) ;

  //des_meridian(f_1_2, 0., 20., "f_1_2", 17) ;
  
  // First part of the correction metric (needed to be complemented by the (1 <-> 2) term
  
  Sym_tensor hh_temp(hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  
  for (int i=1 ; i<= 3 ; i++){
    for (int j=i ; j<= 3 ; j++){
      hh_temp.set(i,j) =   f_delta * hole1.ff.con()(i,j) + 0*f_1_1 * nn1(i)*nn1(j) 
	+ 0*f_1_12 * 0.5 *(nn1(i) * nn12(j) + nn1(j) * nn12(i)) 
	+ 0*f_12_12 * nn12(i)*nn12(j) 
	+ 0*f_1_2 * 0.5*(nn1(i)*nn2(j) + nn1(j)*nn2(i) ) ;
    }
  }

  //des_meridian(hh_temp, 0., 20., "hh_temp", 25) ;
  //arrete() ;

  return hh_temp ;

}


Sym_tensor Bin_hor::hh_Samaya_hole2() {

  //========
  //  Grid 1
  //========
  int nz1 = hole1.mp.get_mg()->get_nzone() ;
  int nz2 = hole2.mp.get_mg()->get_nzone() ;
    
  // General coordinate values
  const Coord& xx_1 = hole1.mp.x ; 
  const Coord& yy_1 = hole1.mp.y ;   
  const Coord& zz_1 = hole1.mp.z ; 

  //========
  //  Grid 2
  //========
  
  // General coordinate values
  const Coord& xx_2 = hole2.mp.x ; 
  const Coord& yy_2 = hole2.mp.y ;   
  const Coord& zz_2 = hole2.mp.z ; 

  //===================================
  // Definition of the relevant vectors
  //===================================

  // Coordinate vector from hole 2 in the grid 2: nn2
  //--------------------------------------------------
  Vector rr2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  rr2.set(1) = xx_2 ;
  rr2.set(2) = yy_2 ;
  rr2.set(3) = zz_2 ;
  rr2.std_spectral_base() ;  

  // Norm r2
  Scalar r2 (hole2.mp) ;
  r2 = hole2.mp.r ;
  r2.std_spectral_base() ;
  Scalar temp2 (r2) ;
  temp2.raccord(1) ;
  r2.set_domain(0) = temp2.domain(0) ;

  // Unitary vector
  Vector nn2 (rr2);
  nn2 = nn2/r2 ;

  for (int i=0; i<hole2.mp.get_mg()->get_nr(nz2-1); i++)
    for (int j=0; j<hole2.mp.get_mg()->get_nt(nz2-1); j++)
      for (int k=0; k<hole2.mp.get_mg()->get_np(nz2-1); k++)
	for (int l=0; l<nz2; l++)
	  for (int ind=1; ind<=3; ind++){
	    if (i != 0 || l != 1)
	      nn2.set(ind).set_grid_point(l,k,j,i) = nn2(ind).val_grid_point(1,k,j,0) ;
	  }

  // Coordinate vector from hole 1 in the grid 1: nn1_1
  //-----------------------------------------------------
  Vector rr1_1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  rr1_1.set(1) = xx_1 ;
  rr1_1.set(2) = yy_1 ;
  rr1_1.set(3) = zz_1 ;
  rr1_1.std_spectral_base() ;  

  // Norm r1_g1 
  Scalar r1_1 (hole1.mp) ;
  r1_1 = hole1.mp.r ;
  r1_1.std_spectral_base() ;
  Scalar temp1 (r1_1) ;
  temp1.raccord(1) ;
  r1_1.set_domain(0) = temp1.domain(0) ;

  // Unitary vector
  Vector nn1_1 (rr1_1);
  nn1_1 = nn1_1/r1_1 ;
  
  for (int i=0; i<hole1.mp.get_mg()->get_nr(nz1-1); i++)
    for (int j=0; j<hole1.mp.get_mg()->get_nt(nz1-1); j++)
      for (int k=0; k<hole1.mp.get_mg()->get_np(nz1-1); k++)
	for (int l=0; l<nz1; l++)
	  for (int ind=1; ind<=3; ind++){
	    if (i != 0 || l != 1)
	      nn1_1.set(ind).set_grid_point(l,k,j,i) = nn1_1(ind).val_grid_point(1,k,j,0) ;
	  }

  Scalar unsr2 (hole2.mp) ;
  unsr2 = 1./hole2.mp.r ;
  unsr2.std_spectral_base() ;
  unsr2.raccord(1) ;
  
  Scalar unsr2_1 (hole1.mp) ;
  unsr2_1.set_etat_qcq() ;
  unsr2_1.import(unsr2) ;
  unsr2_1.set_spectral_va().set_base(unsr2.get_spectral_va().get_base()) ;
  
  Scalar r1sr2_1 (hole1.mp) ;
  r1sr2_1 = r1_1*unsr2_1 ;
  r1sr2_1.set_outer_boundary(nz2-1, 1.) ;

  // Coordinate vector from hole 1 in the grid 2: nn1
  //-----------------------------------------------------
  Vector nn1 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  nn1.change_triad(hole2.mp.get_bvect_cart()) ;
  nn1.set_etat_qcq() ;
  for (int i=1 ; i<=3 ; i++){ 
    nn1.set(i).import(nn1_1(i)) ;
    nn1.set(i).set_spectral_va().set_base(nn1_1(i).get_spectral_va().get_base()) ;
  }  
  
  Scalar r1sr2 (hole2.mp) ;
  r1sr2.set_etat_qcq() ;
  r1sr2.import(r1sr2_1) ;
  r1sr2.set_spectral_va().set_base(r1sr2_1.get_spectral_va().get_base()) ;
  

  // Coordinate vector from hole 2 to hole 1 in the grid 2: nn21
  //----------------------------------------------------------------
  // Warning! Valid only in the symmetruc case (for the general case it would
  // necessary to construct this whole function as a Bin_hor function 
  Vector rr21 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  rr21.set(1) = hole2.mp.get_ori_x() - hole1.mp.get_ori_x() ;
  rr21.set(2) = hole2.mp.get_ori_y() - hole1.mp.get_ori_y() ;
  rr21.set(3) = hole2.mp.get_ori_z() - hole1.mp.get_ori_z() ;
  rr21.std_spectral_base() ;  

  //Norm r21
  Scalar r21 (hole2.mp) ;
  r21 = sqrt(rr21(1)*rr21(1) + rr21(2)*rr21(2) + rr21(3)*rr21(3)) ;
  r21.std_spectral_base() ;

  // Unitary vector
  Vector nn21 ( rr21 );
  nn21 = nn21/ r21 ;


  Scalar f_delta (hole2.mp) ;
  Scalar f_delta_zec (hole2.mp) ;
  Scalar f_2_2 (hole2.mp) ;
  Scalar f_2_2_zec (hole2.mp) ;
  Scalar f_2_21 (hole2.mp) ;
  Scalar f_2_21_zec (hole2.mp) ;
  Scalar f_21_21 (hole2.mp) ;
  Scalar f_2_1 (hole2.mp) ;

  f_delta.set_etat_qcq() ;
  f_delta_zec.set_etat_qcq() ;
  f_2_2.set_etat_qcq() ;
  f_2_2_zec.set_etat_qcq() ;
  f_2_21.set_etat_qcq() ;
  f_2_21_zec.set_etat_qcq() ;
  f_21_21.set_etat_qcq() ;
  f_2_1.set_etat_qcq() ;

  // Function exp(-(r-r_0)^2/sigma^2)
  // --------------------------------
  
  double r0 = hole2.mp.val_r(nz2-2, 1, 0, 0) ;
  double sigma = 1.*r0 ;
  
  Scalar rr (hole2.mp) ;
  rr = hole2.mp.r ;

  Scalar fexp (hole2.mp) ;
  fexp = exp( -(rr - r0)*(rr - r0)/sigma/sigma ) ;
  for (int ii=0; ii<nz2-1; ii++)
    fexp.set_domain(ii) = 1. ;
  fexp.set_outer_boundary(nz2-1, 0) ;
  fexp.std_spectral_base() ;

  // Conformal metric
  //=================
 
   // f_delta
  //--------
  f_delta = -5.*r2/(8.*r21*r21*r21) - 15./(8.*r2*r21) + 
    5.*r2/(8.*r21*r21*r21*r1sr2) + 1./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2))*
    (1 + r2/r21 + r21/r2 - 1./r1sr2 - r2/(r1sr2*r21) + r21*r21/(2*r2*r2*r1sr2)) +
    1./(r2*(1.+r1sr2+r21/r2))*(-7./r2 + 2./r21) ;  

  f_delta.annule_domain(nz2-1) ;
  
  f_delta_zec = - 15./(8.*r2*r21) + 1./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2))*
    (1 + r2/r21 + r21/r2 - 1./r1sr2 - r2/(r1sr2*r21) + r21*r21/(2*r2*r2*r1sr2)) +
    1./(r2*(1.+r1sr2+r21/r2))*(-7./r2 + 2./r21) ;  
  f_delta_zec += fexp*(-5.*r2/(8.*r21*r21*r21) + 5.*r2/(8.*r21*r21*r21*r1sr2)) ;

  f_delta_zec.set_outer_boundary(nz2-1, 0.) ;
  for (int i=0 ;i<nz2-1 ; i++){
    f_delta_zec.annule_domain(i) ;
  }
  
  f_delta = f_delta + f_delta_zec ;

  //des_meridian(f_delta, 0., 20., "f_delta", 23) ;
  


  // f_2_2
  //------
  f_2_2 = r2/(8.*r21*r21*r21) + 11./(8.*r2*r21) -
    r1sr2*r1sr2*r2/(8.*r21*r21*r21) + 7./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2)) +
    7./(r2*r2*(1.+r1sr2+r21/r2)) ;

  f_2_2.annule_domain(nz2-1) ;
  
  f_2_2_zec = 11./(8.*r2*r21) + 7./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2)) +
    7./(r2*r2*(1.+r1sr2+r21/r2)) ;
  f_2_2_zec += fexp*(r2/(8.*r21*r21*r21) - r1sr2*r1sr2*r2/(8.*r21*r21*r21)) ;

  f_2_2_zec.set_outer_boundary(nz2-1, 0.) ;
  for (int i=0 ; i<nz2-1 ; i++){
    f_2_2_zec.annule_domain(i) ;
  }
  
  f_2_2 = f_2_2 + f_2_2_zec ;

  //des_meridian(f_2_2, 0., 20., "f_2_2", 24) ;

  // f_2_21
  //------
  f_2_21 = - 0*7./(2*r21*r21) + 8./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2)) ; // ????????
  f_2_21.annule_domain(nz2-1) ;
  
  f_2_21_zec = 8./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2)) ;
  f_2_21_zec.set_outer_boundary(nz2-1, 0.) ;

  for (int i=0 ; i<nz2-1 ; i++){
    f_2_21_zec.annule_domain(i) ;
  }
   
  f_2_21 = f_2_21 + f_2_21_zec ;
  
  //des_meridian(f_2_21, 0., 20., "f_2_21", 25) ;
 
  // f_21_21
  //-------
  f_21_21 = 2* (-4./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2)) -  // 2  ???????
		4./(r21*r2*(1.+r1sr2+r21/r2))) ;
  f_21_21.set_outer_boundary(nz2-1, 0.) ;

  //des_meridian(f_21_21, 0., 20., "f_21_21", 26) ;

  // f_2_1
  //-------
  f_2_1 = 2 * 11./(r2*r2*(1.+r1sr2+r21/r2)*(1.+r1sr2+r21/r2));        // 2  ???????
  f_2_1.set_outer_boundary(nz2-1, 0.) ;
 

  // First part of the correction metric (needed to be complemented by the (1 <-> 2) term

  
  Sym_tensor hh_temp(hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  
  for (int i=1 ; i<= 3 ; i++){
    for (int j=1 ; j<= 3 ; j++){
      hh_temp.set(i,j) =  f_delta * hole2.ff.con()(i,j) + 0*f_2_2 * nn2(i)*nn2(j) 
	+ 0*f_2_21 * 0.5 *(nn2(i) * nn21(j) + nn2(j) * nn21(i)) 
	+ 0*f_21_21 * nn21(i)*nn21(j) 
	+ 0*f_2_1 * 0.5*(nn2(i)*nn1(j) + nn2(j)*nn1(i) );
    }
  }
  
  //des_meridian(hh_temp, 0., 20., "hh_temp", 27) ;
  //arrete() ;

  return hh_temp ;

}

void Bin_hor::set_hh_Samaya() {

  Sym_tensor hh1 ( hh_Samaya_hole1() ) ;  
  Sym_tensor hh2 ( hh_Samaya_hole2() ) ;


  // Graphical outputs
  // ------------------

  des_meridian(hh1(1,1), 0., 80., "hh1", 10) ;
  arrete() ;
  des_profile(hh1(1,1), 0., 80., M_PI/2, M_PI) ;
  des_profile(hh1(1,1), 0., 80., M_PI/2, 0) ;
  des_profile(hh1(1,1), 0., 80., 0, M_PI) ;
  des_coupe_z (hh1(1,1), 0., 2) ;
  des_coupe_z (hh1(1,1), 0., 3) ;
  des_coupe_z (hh1(1,1), 0., 4) ;
  des_coupe_z (hh1(1,1), 0., 5) ;
  des_coupe_z (hh1(1,1), 0., 6) ;

  // Definition of the surface
  // -------------------------
  
  Cmp surface_un (hole1.mp) ;
  surface_un = pow(hole1.mp.r, 2.)-pow(hole1.get_radius(), 2.) ;
  surface_un.annule(hole1.mp.get_mg()->get_nzone()-1) ;
  surface_un.std_base_scal() ;
  
  Cmp surface_deux (hole2.mp) ;
  surface_deux = pow(hole2.mp.r, 2.)-pow(hole2.get_radius(), 2.) ;
  surface_deux.annule(hole1.mp.get_mg()->get_nzone()-1) ;
  surface_deux.std_base_scal() ;

  double ta = 12 ;
  Cmp dessin_un (hh1(1,1)) ;
  dessin_un.annule(0) ;
  
  Cmp dessin_deux (hh2(1,1)) ;
  dessin_deux.annule(0) ;
  
  des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "hh(1,1)", &surface_un, &surface_deux, 
		   false, 15, 300, 300) ;
 


  // Importation 
  // ----------------

  Sym_tensor hh2_1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  hh2_1.set_etat_qcq() ;  
  Sym_tensor hh1_2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  hh1_2.set_etat_qcq() ;  

  hh2.change_triad(hole1.mp.get_bvect_cart()) ;
  for (int i=1 ; i<=3 ; i++){ 
    for (int j=i ; j<=3 ; j++){ 
    hh2_1.set(i,j).import(hh2(i,j)) ;
    hh2_1.set(i,j).set_spectral_va().set_base(hh2(i,j).get_spectral_va().get_base()) ;
    }    
  }  
  hh2.change_triad(hole2.mp.get_bvect_cart()) ;

  hh1.change_triad(hole2.mp.get_bvect_cart()) ;
  for (int i=1 ; i<=3 ; i++){ 
    for (int j=i ; j<=3 ; j++){ 
      hh1_2.set(i,j).import(hh1(i,j)) ;
      hh1_2.set(i,j).set_spectral_va().set_base(hh1(i,j).get_spectral_va().get_base()) ;
    }    
  }  
  hh1.change_triad(hole1.mp.get_bvect_cart()) ;

  double m1, m2 ;
  m1 = pow(hole1.area_hor()/(16.*M_PI) + hole1.ang_mom_hor()/hole1.radius, 0.5) ;
  m2 = pow(hole2.area_hor()/(16.*M_PI) + hole2.ang_mom_hor()/hole2.radius, 0.5) ;
  
  
  hh1 = hh1 + hh2_1 ;
  hh2 = hh2 + hh1_2 ;

  des_meridian(hh1, 0., 20., "hh1 cart", 20) ;

  hh1.change_triad(hole1.mp.get_bvect_spher()) ;
  hh2.change_triad(hole2.mp.get_bvect_spher()) ;

  hole1.hh = m1*m2* hh1 ;
  hole2.hh = m1*m2* hh2 ;

  Metric tgam_1 ( hole1.ff.con() +  hh1 ) ;
  Metric tgam_2 ( hole2.ff.con() +  hh2 ) ;

  hole1.tgam = tgam_1 ;
  hole2.tgam = tgam_2 ;

  
  des_meridian(hh1, 0., 20., "hh1", 0) ;


}
