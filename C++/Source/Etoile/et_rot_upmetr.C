/*
 * Method Etoile_rot::update_metric
 *
 * (see file etoile.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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
 * Revision 1.4  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:52:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2004/10/11 15:09:02  j_novak
 * The radial manipulation functions take Scalar as arguments, instead of Cmp.
 * Added a conversion operator from Scalar to Cmp.
 * The Cmp radial manipulation function make conversion to Scalar, call to the
 * Map_radial version with a Scalar argument and back.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2001/01/25  13:01:28  eric
 * Appel de set_std_base() sur bbb et b_car.
 *
 * Revision 2.3  2000/11/20  21:43:08  eric
 * Ajout de bbb.set_etat_qcq() avant bbb.set().
 *
 * Revision 1.2  2000/09/18  16:15:52  eric
 * *** empty log message ***
 *
 * Revision 1.1  2000/07/20  15:33:06  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "etoile.h"


namespace Lorene {
void Etoile_rot::update_metric() {
 
    // Lapse function N
    // ----------------
    
    nnn = exp( unsurc2 * logn ) ; 

    nnn.set_std_base() ;   // set the bases for spectral expansions
    
    
    // Metric factor A^2
    // -----------------
    
    a_car = exp( 2*unsurc2*( dzeta - logn ) ) ; 

    a_car.set_std_base() ;   // set the bases for spectral expansions

    // Metric factor B
    // ---------------
    
    Cmp tmp = tggg() ; 
    tmp.div_rsint() ;	        //... Division of tG by r sin(theta)

    bbb.set_etat_qcq() ;
    bbb.set() = tmp ;
    bbb = (bbb + 1) / nnn ; 

    bbb.set_std_base() ;   // set the bases for spectral expansions
    
    // Metric factor B^2
    // -----------------
    
    b_car = bbb * bbb ; 
    
    b_car.set_std_base() ;   // set the bases for spectral expansions

    // Tensor B^{-2} K_{ij} and Scalar A^2 K_{ij} K^{ij}
    // -------------------------------------------------
    
    extrinsic_curvature() ; 
    
  
    // The derived quantities are no longer up to date : 
    // -----------------------------------------------

    del_deriv() ;  

}
}
