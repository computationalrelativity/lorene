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
#include "isol_hor.h"
#include "graphique.h"



Sym_tensor Bin_hor::hh_Nissanke_hole1 () {
  /*
 // Les alignemenents pour le signe des CL.
  double orientation_1 = hole1.mp.get_rot_phi() ;
  assert ((orientation_1 == 0) || (orientation_1 == M_PI)) ;
  int aligne_1 = (orientation_1 == 0) ? 1 : -1 ;
 // Les alignemenents pour le signe des CL.
  double orientation_2 = hole2.mp.get_rot_phi() ;
  assert ((orientation_2 == 0) || (orientation_2 == M_PI)) ;
  int aligne_2 = (orientation_2 == 0) ? 1 : -1 ;
  */
  //========
  //  Grid 1
  //========
  int nn_z_1 = hole1.mp.get_mg()->get_nzone() ;
  
  // General coordinate values
  const Coord& xx_g1 = hole1.mp.x ; 
  const Coord& yy_g1 = hole1.mp.y ;   
  const Coord& zz_g1 = hole1.mp.z ; 

  const Coord& xx_a_g1 = hole1.mp.xa ; 
  const Coord& yy_a_g1 = hole1.mp.ya ;   
  const Coord& zz_a_g1 = hole1.mp.za ; 

  //========
  //  Grid 2
  //========
  
  // General coordinate values
  const Coord& xx_g2 = hole2.mp.x ; 
  const Coord& yy_g2 = hole2.mp.y ;   
  const Coord& zz_g2 = hole2.mp.z ; 

  const Coord& xx_a_g2 = hole2.mp.xa ; 
  const Coord& yy_a_g2 = hole2.mp.ya ;   
  const Coord& zz_a_g2 = hole2.mp.za ; 


  //===================================
  // Definition of the relevant vectors
  //===================================

  // Coordinate vector from hole 1 in the grid 1: nn_1_g1
  //-----------------------------------------------------
  Vector rr_1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  rr_1.set(1) = xx_g1 ;
  rr_1.set(2) = yy_g1 ;
  rr_1.set(3) = zz_g1 ;
  rr_1.std_spectral_base() ;  

  // Norm r_1
  Scalar r_1 (hole1.mp) ;
  const Coord& rcord_1 = hole1.mp.r ; 
  r_1 = rcord_1 ;
  r_1.raccord(1) ;
  r_1.std_spectral_base() ;

  // Unitary vector
  Vector nn_1 ( rr_1 );
  nn_1 = nn_1/r_1 ;


  // Coordinate vector from hole 2 in the grid 2: nn_2_g2
  //-----------------------------------------------------
  Vector rr_2_g2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  rr_2_g2.set(1) = xx_g2 ;
  rr_2_g2.set(2) = yy_g2 ;
  rr_2_g2.set(3) = zz_g2 ;
  rr_2_g2.std_spectral_base() ;  


  // Norm r_2_g2 
  Scalar r_2_g2 (hole2.mp) ;
  const Coord& rcord_g2 = hole2.mp.r ; 
  r_2_g2 = rcord_g2 ;
  r_2_g2.raccord(1) ;
  r_2_g2.std_spectral_base() ;

  // Unitary vector
  Vector nn_2_g2 ( rr_2_g2 );
  nn_2_g2 = nn_2_g2/r_2_g2 ;
  for (int i=1; i<=3; i++)
    nn_2_g2.set(i).raccord(1) ;
  

  // Coordinate vector from hole 2 in the grid 1: nn_2_g1
  //-----------------------------------------------------
  Vector nn_2_g1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  nn_2_g2.change_triad(hole1.mp.get_bvect_cart()) ;
  nn_2_g1.set_etat_qcq() ;
  for (int i=1 ; i<=3 ; i++){ 
    nn_2_g1.set(i).import(nn_2_g2(i)) ;
    nn_2_g1.set(i).set_spectral_va().set_base(nn_2_g2(i).get_spectral_va().get_base()) ;
  }  
  
  Scalar r_2_g1 (hole1.mp) ;
  r_2_g1.set_etat_qcq() ;
  r_2_g1.import(r_2_g2) ;
  r_2_g1.set_spectral_va().set_base(r_2_g2.get_spectral_va().get_base()) ;
  
 

  // Coordinate vector from hole 1 to hole 2 in the grid 1: nn_12_g1
  //----------------------------------------------------------------
  // Warning! Valid only in the symmetruc case (for the general case it would
  // necessary to construct this whole function as a Bin_hor function 
  Vector rr_12 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;

  rr_12.set(1) = xx_a_g1 - xx_a_g2 ;
  rr_12.set(2) = yy_a_g1 - yy_a_g2 ; ;
  rr_12.set(3) = zz_a_g1 - zz_a_g2 ;;
  rr_12.std_spectral_base() ;  

  //Norm r_12
  Scalar r_12 (hole1.mp) ;
  r_12 = sqrt( rr_12(1)*rr_12(1) + rr_12(2)*rr_12(2) + rr_12(3)*rr_12(3)) ;
  r_12.std_spectral_base() ;

  // Unitary vector
  Vector nn_12 ( rr_12 );
  nn_12 = nn_12/ r_12 ;


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


 
  // Conformal metric
  //=================

  // tilde{gamma}- \delta = m_1*m_2* ( f_delta \delta_{ij} 
  //                        + f_1_1 nn_1*nn_1 + f_1_12 nn_1*nn_12
  //                        + f_12_12 nn_12*nn_12
  //                        + f_2_2 nn_2*nn_2 + f_2_12 nn_2*nn_12
  //                        + f_1_2 nn_1*nn_2
 
  // f_delta
  //--------
  f_delta = -5.*r_1/(8.*r_12*r_12*r_12) - 15./(8.*r_1*r_12)+ //(Beginning 1st part)
    5.*r_1*r_1/(8.*r_12*r_12*r_12)/r_2_g1 + 1./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12)*
    (1 + r_1/r_12 + r_12/r_1 - r_1/r_2_g1 - r_1*r_1/(r_2_g1*r_12) + 
     r_12*r_12/(2*r_1*r_2_g1)) +
    1./(r_1+r_2_g1+r_12)*(-7./r_1 + 2./r_12) ;  // (End 1st part)

  f_delta.annule_domain(nn_z_1-1) ;
  //  f_delta.annule_domain(0) ;
  
  f_delta_zec = -15./(8.*r_1*r_12)+ //(Beginning 1st part)
    1./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12)*
    (1 + r_1/r_12 + r_12/r_1 - r_1/r_2_g1 - r_1*r_1/(r_2_g1*r_12) + 
     r_12*r_12/(2*r_1*r_2_g1)) +
    1./(r_1+r_2_g1+r_12)*(-7./r_1 + 2./r_12) ;
  
  for (int i=0 ;i<nn_z_1-1 ; i++){
    f_delta_zec.annule_domain(i) ;
  }
  
  f_delta = f_delta + f_delta_zec ;


  // f_1_1
  //------
  f_1_1 = r_1/(8.*r_12*r_12*r_12) + 11./(8.*r_1*r_12) -
    r_2_g1*r_2_g1/(8.*r_1*r_12*r_12*r_12) + 7./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12) +
    7./r_1/(r_1+r_2_g1+r_12) ;

  f_1_1.annule_domain(nn_z_1-1) ;
  //  f_1_1.annule_domain(0) ;
  
  f_1_1_zec = + 11./(8.*r_1*r_12) -
     + 7./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12) +
    7./r_1/(r_1+r_2_g1+r_12) ;

  for (int i=0 ; i<nn_z_1-1 ; i++){
    f_1_1_zec.annule_domain(i) ;
  }
  
  f_1_1 = f_1_1 + f_1_1_zec ;


  // f_1_12
  //------
  f_1_12 = - 7./(2*r_12*r_12) + 8./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12) ;

  f_1_12.annule_domain(nn_z_1-1) ;
  //  f_1_12.annule_domain(0) ;
  
  f_1_12_zec = 8./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12) ;

  for (int i=0 ; i<nn_z_1-1 ; i++){
    f_1_12_zec.annule_domain(i) ;
  }
   
  f_1_12 = f_1_12 + f_1_12_zec ;
  
 
  // f_12_12
  //-------
  f_12_12 = 2* (-4./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12) - 4./r_12/(r_1+r_2_g1+r_12)) ;

  // f_1_2
  //-------
  f_1_2 = 2* (11./(r_1+r_2_g1+r_12)/(r_1+r_2_g1+r_12)) ;


  // First part of the correction metric (needed to be complemented by the (1 <-> 2) term

  
  Sym_tensor hh_temp(hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  
  for (int i=1 ; i<= 3 ; i++){
    for (int j=i ; j<= 3 ; j++){
      hh_temp.set(i,j) =   f_delta * hole1.ff.con()(i,j) + f_1_1 * nn_1(i)*nn_1(j) 
	+ f_1_12 * 0.5 *(nn_1(i) * nn_12(j) + nn_1(j) * nn_12(i)) 
	+ f_12_12 * nn_12(i)*nn_12(j) 
	+ f_1_2 * 0.5*(nn_1(i)*nn_2_g1(j) + nn_1(j)*nn_2_g1(i) ) ;
    }
  }
  

  return hh_temp ;

}






Sym_tensor Bin_hor::hh_Nissanke_hole2 () {
  /*
// Les alignemenents pour le signe des CL.
  double orientation_1 = hole1.mp.get_rot_phi() ;
  assert ((orientation_1 == 0) || (orientation_1 == M_PI)) ;
  int aligne_1 = (orientation_1 == 0) ? 1 : -1 ;
 // Les alignemenents pour le signe des CL.
  double orientation_2 = hole2.mp.get_rot_phi() ;
  assert ((orientation_2 == 0) || (orientation_2 == M_PI)) ;
  int aligne_2 = (orientation_2 == 0) ? 1 : -1 ;
  */
 
  //========
  //  Grid 1
  //========
    
  // General coordinate values
  const Coord& xx_g1 = hole1.mp.x ; 
  const Coord& yy_g1 = hole1.mp.y ;   
  const Coord& zz_g1 = hole1.mp.z ; 

  const Coord& xx_a_g1 = hole1.mp.xa ; 
  const Coord& yy_a_g1 = hole1.mp.ya ;   
  const Coord& zz_a_g1 = hole1.mp.za ; 

  //========
  //  Grid 2
  //========
  int nn_z_2 = hole2.mp.get_mg()->get_nzone() ;
  
  // General coordinate values
  const Coord& xx_g2 = hole2.mp.x ; 
  const Coord& yy_g2 = hole2.mp.y ;   
  const Coord& zz_g2 = hole2.mp.z ; 

  const Coord& xx_a_g2 = hole2.mp.xa ; 
  const Coord& yy_a_g2 = hole2.mp.ya ;   
  const Coord& zz_a_g2 = hole2.mp.za ; 


  //===================================
  // Definition of the relevant vectors
  //===================================

  // Coordinate vector from hole 2 in the grid 2: nn_2
  //-----------------------------------------------------
  Vector rr_2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;

  rr_2.set(1) = xx_g2 ;
  rr_2.set(2) = yy_g2 ;
  rr_2.set(3) = zz_g2 ;
  rr_2.std_spectral_base() ;  

  // Norm r_2
  Scalar r_2 (hole2.mp) ;
  const Coord& rcord_2 = hole2.mp.r ; 
  r_2 = rcord_2 ;
  r_2.raccord(1) ;
  r_2.std_spectral_base() ;

  // Unitary vector
  Vector nn_2 ( rr_2 );
  nn_2 = nn_2/r_2 ;


  // Coordinate vector from hole 1 in the grid 1: nn_1_g1
  //-----------------------------------------------------
  Vector rr_1_g1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;

  rr_1_g1.set(1) = xx_g1 ;
  rr_1_g1.set(2) = yy_g1 ;
  rr_1_g1.set(3) = zz_g1 ;
  rr_1_g1.std_spectral_base() ;  


  // Norm r_1_g1 
  Scalar r_1_g1 (hole1.mp) ;
  const Coord& rcord_g1 = hole1.mp.r ; 
  r_1_g1 = rcord_g1 ;
  r_1_g1.raccord(1) ;
  r_1_g1.std_spectral_base() ;

  // Unitary vector
  Vector nn_1_g1 ( rr_1_g1 );
  nn_1_g1 = nn_1_g1/r_1_g1 ;
  for (int i=1; i<=3; i++)
    nn_1_g1.set(i).raccord(1) ;

  // Coordinate vector from hole 1 in the grid 2: nn_1_g2
  //-----------------------------------------------------
  Vector nn_1_g2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  nn_1_g1.change_triad(hole2.mp.get_bvect_cart()) ;
  nn_1_g2.set_etat_qcq() ;
  for (int i=1 ; i<=3 ; i++){ 
    nn_1_g2.set(i).import(nn_1_g1(i)) ;
    nn_1_g2.set(i).set_spectral_va().set_base(nn_1_g1(i).get_spectral_va().get_base()) ;
  }  
  
  Scalar r_1_g2 (hole1.mp) ;
  r_1_g2.set_etat_qcq() ;
  r_1_g2.import(r_1_g1) ;
  r_1_g2.set_spectral_va().set_base(r_1_g1.get_spectral_va().get_base()) ;
  
 

  // Coordinate vector from hole 2 to hole 1 in the grid 2: nn_21_g2
  //----------------------------------------------------------------
  // Warning! Valid only in the symmetruc case (for the general case it would
  // necessary to construct this whole function as a Bin_hor function 
  Vector rr_21 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;

  rr_21.set(1) = xx_a_g2 - xx_a_g1 ;
  rr_21.set(2) = yy_a_g2 - yy_a_g1 ; ;
  rr_21.set(3) = zz_a_g2 - zz_a_g1 ;;
  rr_21.std_spectral_base() ;  

  //Norm r_21
  Scalar r_21 (hole2.mp) ;
  r_21 = sqrt( rr_21(1)*rr_21(1) + rr_21(2)*rr_21(2) + rr_21(3)*rr_21(3)) ;
  r_21.std_spectral_base() ;

  // Unitary vector
  Vector nn_21 ( rr_21 );
  nn_21 = nn_21/ r_21 ;


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



  // Conformal metric
  //=================

 
  // f_delta
  //--------
  f_delta = -5.*r_2/(8.*r_21*r_21*r_21) - 15./(8.*r_2*r_21)+ 
    5.*r_2*r_2/(8.*r_21*r_21*r_21)/r_1_g2 + 1./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21)*
    (1 + r_2/r_21 + r_21/r_2 - r_2/r_1_g2 - r_2*r_2/(r_1_g2*r_21) + 
     r_21*r_21/(2*r_2*r_1_g2)) +
    1./(r_2+r_1_g2+r_21)*(-7./r_2 + 2./r_21) ;  // (End 1st part)

  f_delta.annule_domain(nn_z_2-1) ;
  //  f_delta.annule_domain(0) ;
  
  f_delta_zec = -15./(8.*r_2*r_21)+ //(Beginning 1st part)
    1./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21)*
    (1 + r_2/r_21 + r_21/r_2 - r_2/r_1_g2 - r_2*r_2/(r_1_g2*r_21) + 
     r_21*r_21/(2*r_2*r_1_g2)) +
    1./(r_2+r_1_g2+r_21)*(-7./r_2 + 2./r_21) ;
  
  for (int i=0 ;i<nn_z_2-1 ; i++){
    f_delta_zec.annule_domain(i) ;
  }
  
  f_delta = f_delta + f_delta_zec ;


  // f_2_2
  //------
  f_2_2 = r_2/(8.*r_21*r_21*r_21) + 11./(8.*r_2*r_21) -
    r_1_g2*r_1_g2/(8.*r_2*r_21*r_21*r_21) + 7./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21) +
    7./r_2/(r_2+r_1_g2+r_21) ;

  f_2_2.annule_domain(nn_z_2-1) ;
  //  f_2_2.annule_domain(0) ;
  
  f_2_2_zec = + 11./(8.*r_2*r_21) -
     + 7./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21) +
    7./r_2/(r_2+r_1_g2+r_21) ;

  for (int i=0 ; i<nn_z_2-1 ; i++){
    f_2_2_zec.annule_domain(i) ;
  }
  
  f_2_2 = f_2_2 + f_2_2_zec ;


  // f_2_21
  //------
  f_2_21 = - 7./(2*r_21*r_21) + 8./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21) ;

  f_2_21.annule_domain(nn_z_2-1) ;
  //  f_2_21.annule_domain(0) ;
  
  f_2_21_zec = 8./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21) ;

  for (int i=0 ; i<nn_z_2-1 ; i++){
    f_2_21_zec.annule_domain(i) ;
  }
   
  f_2_21 = f_2_21 + f_2_21_zec ;
  
 
  // f_21_21
  //-------
  f_21_21 = 2* (-4./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21) - 4./r_21/(r_2+r_1_g2+r_21)) ;

  // f_2_1
  //-------
  f_2_1 = 2* (11./(r_2+r_1_g2+r_21)/(r_2+r_1_g2+r_21)) ;


  // First part of the correction metric (needed to be complemented by the (1 <-> 2) term

  
  Sym_tensor hh_temp(hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  
  for (int i=1 ; i<= 3 ; i++){
    for (int j=1 ; j<= 3 ; j++){
      hh_temp.set(i,j) =  f_delta * hole2.ff.con()(i,j) + f_2_2 * nn_2(i)*nn_2(j) 
	+ f_2_21 * 0.5 *(nn_2(i) * nn_21(j) + nn_2(j) * nn_21(i)) 
	+ f_21_21 * nn_21(i)*nn_21(j) 
	+ f_2_1 * 0.5*(nn_2(i)*nn_1_g2(j) + nn_2(j)*nn_1_g2(i) );
    }
  }
  
  return hh_temp ;


}



void Bin_hor::set_hh_Nissanke () {

  Sym_tensor hh1 ( hh_Nissanke_hole1() ) ;  
  Sym_tensor hh2 ( hh_Nissanke_hole2() ) ;


  Sym_tensor hh2_g1 (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
  hh2_g1.set_etat_qcq() ;  
  Sym_tensor hh1_g2 (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
  hh1_g2.set_etat_qcq() ;  

  for (int i=1 ; i<=3 ; i++){ 
    for (int j=i ; j<=3 ; j++){ 
    hh2_g1.set(i,j).import(hh2(i,j)) ;
    hh2_g1.set(i,j).set_spectral_va().set_base(hh2(i,j).get_spectral_va().get_base()) ;
    }    
  }  
  for (int i=1 ; i<=3 ; i++){ 
    for (int j=i ; j<=3 ; j++){ 
      hh1_g2.set(i,j).import(hh1(i,j)) ;
      hh1_g2.set(i,j).set_spectral_va().set_base(hh1(i,j).get_spectral_va().get_base()) ;
    }    
  }  


  double m1, m2 ;
  m1 = pow(hole1.area_hor()/(16.*M_PI) + hole1.ang_mom_hor()/hole1.radius, 0.5) ;
  m2 = pow(hole2.area_hor()/(16.*M_PI) + hole2.ang_mom_hor()/hole2.radius, 0.5) ;
  
  
  hh1 =  hh1 + hh2_g1 ;
  hh2 = hh2 + hh1_g2 ;

  hole1.hh = m1*m2* hh1 ;
  hole2.hh = m1*m2* hh2 ;

  Metric tgam_1 ( hole1.ff.con() +  hh1 ) ;
  Metric tgam_2 ( hole2.ff.con() +  hh2 ) ;

  hole1.tgam = tgam_1 ;
  hole2.tgam = tgam_2 ;

}
