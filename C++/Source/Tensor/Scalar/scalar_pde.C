/*
 * Methods of the class Scalar for various partial differential equations
 *
 *  See file scalar.h for documentation. 
 */

/*
 *   Copyright (c) 2003-2005 Eric Gourgoulhon & Jerome Novak
 *   Copyright (c) 2004 Philippe Grandclement
 *
 *   Copyright (c) 1999-2001 Eric Gourgoulhon (for preceding class Cmp)
 *   Copyright (c) 1999-2001 Philippe Grandclement (for preceding class Cmp)
 *   Copyright (c) 2000-2001 Jerome Novak (for preceding class Cmp)
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


char scalar_pde_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.13  2005/04/04 21:34:44  e_gourgoulhon
 * Added argument lambda to method poisson_angu
 * to deal with the generalized angular Poisson equation:
 *     Lap_ang u + lambda u = source.
 *
 * Revision 1.12  2004/08/24 09:14:52  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.11  2004/06/22 08:50:00  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.10  2004/05/25 14:30:48  f_limousin
 * Minor modif.
 *
 * Revision 1.9  2004/03/17 15:58:50  p_grandclement
 * Slight modification of sol_elliptic_no_zec
 *
 * Revision 1.8  2004/03/01 09:57:04  j_novak
 * the wave equation is solved with Scalars. It now accepts a grid with a
 * compactified external domain, which the solver ignores and where it copies
 * the values of the field from one time-step to the next.
 *
 * Revision 1.7  2004/02/11 09:47:46  p_grandclement
 * Addition of a new elliptic solver, matching with the homogeneous solution
 * at the outer shell and not solving in the external domain (more details
 * coming soon ; check your local Lorene dealer...)
 *
 * Revision 1.6  2004/01/28 16:59:14  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.5  2004/01/28 16:46:24  p_grandclement
 * Addition of the sol_elliptic_fixe_der_zero stuff
 *
 * Revision 1.4  2004/01/14 10:11:51  f_limousin
 * Corrected bug in poisson with parameters.
 *
 * Revision 1.3  2003/12/11 14:48:51  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * Revision 1.2  2003/10/15 21:14:23  e_gourgoulhon
 * Added method poisson_angu().
 *
 * Revision 1.1  2003/09/25 08:06:56  e_gourgoulhon
 * First versions (use Cmp as intermediate quantities).
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "scalar.h"
#include "tensor.h"
#include "param.h"
#include "cmp.h"
#include "param_elliptic.h"

  
		    //-----------------------------------//
		    //      Scalar Poisson equation      //
		    //-----------------------------------//

// Version without parameters
// --------------------------

Scalar Scalar::poisson() const {
    
    Param bidon ;
    Cmp csource(*this) ; 
    Cmp cresu(mp) ;     
    cresu = 0. ;

    mp->poisson(csource, bidon, cresu) ; 

    Scalar resu(cresu) ; 
    return resu ;          
}

// Version with parameters
// -----------------------

void Scalar::poisson(Param& par, Scalar& uu) const {
    
    Cmp csource(*this) ; 
    Cmp cuu(uu) ;     

    mp->poisson(csource, par, cuu) ;     
    
    uu = cuu ; 
}


		    //-----------------------------------//
		    //      Angular Poisson equation	 //
		    //-----------------------------------//


Scalar Scalar::poisson_angu(double lambda) const {
    
    Param bidon ;

    Scalar resu(*mp) ; 
    resu = 0. ;
		
    mp->poisson_angu(*this, bidon, resu, lambda) ; 

    return resu ;          
}


		    //-----------------------------------//
		    //      Scalar d'Alembert equation	 //
		    //-----------------------------------//

Scalar Scalar::avance_dalembert(Param& par, const Scalar& fjm1, 
	const Scalar& source) const {

  Scalar fjp1(*mp) ;
  
  mp->dalembert(par, fjp1, *this, fjm1, source) ;
	
  return fjp1 ;

}


		    //-----------------------------------//
		    //      General elliptic equation	 //
		    //-----------------------------------//


Scalar Scalar::sol_elliptic(Param_elliptic& ope_var) const {

  // Right now, only applicable with affine mapping or log one
  const Map_af* map_affine = dynamic_cast <const Map_af*> (mp) ;
  const Map_log* map_log = dynamic_cast <const Map_log*> (mp) ;

  if ((map_affine == 0x0) && (map_log == 0x0))  {
    cout << "sol_elliptic only defined for affine or log mapping" << endl ;
    abort() ;
  }
  
  Scalar res (*mp) ;
  res.set_etat_qcq() ;
  
  if (map_affine != 0x0)
    map_affine->sol_elliptic (ope_var, *this, res) ;
  else
    map_log->sol_elliptic (ope_var, *this, res) ;

  return (res) ;
}

     
		    //-----------------------------------//
		    //      General elliptic equation	 //
                    //             with no ZEC           //
		    //-----------------------------------//

Scalar Scalar::sol_elliptic_no_zec(Param_elliptic& ope_var, double val) const {

  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (mp) ;
  const Map_log* map_log = dynamic_cast <const Map_log*> (mp) ;

  if ((map_affine == 0x0) && (map_log == 0x0)) {
    cout << "sol_elliptic_no_zec only defined for affine or log mapping" << endl ;
    abort() ;
  }
  
  Scalar res (*mp) ;
  res.set_etat_qcq() ;
  
  if (map_affine != 0x0)
    map_affine->sol_elliptic_no_zec (ope_var, *this, res, val) ;
  else 
    map_log->sol_elliptic_no_zec (ope_var, *this, res, val) ;

  return (res) ;
}		
 
                    //-----------------------------------//  
                    //      General elliptic equation	 //
                    //             with no ZEC           //
		    //-----------------------------------//

Scalar Scalar::sol_elliptic_only_zec(Param_elliptic& ope_var, double val) const {

  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (mp) ;

  if (map_affine == 0x0) {
    cout << "sol_elliptic_no_zec only defined for affine or log mapping" << endl ;
    abort() ;
  }
  
  Scalar res (*mp) ;
  res.set_etat_qcq() ;
  
  map_affine->sol_elliptic_only_zec (ope_var, *this, res, val) ;
  return (res) ;
}
    		    //-----------------------------------//
		    //      General elliptic equation	 //
                    //         with no ZEC and a         //
                    //        matching with sin(r)/r     //
		    //-----------------------------------//

Scalar Scalar::sol_elliptic_sin_zec(Param_elliptic& ope_var, 
				    double freq, int nbr_phase, 
				    double& ampli_min ,double& phase_min) 
  const {

  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (mp) ;
  
  if (map_affine == 0x0) {
    cout << "sol_elliptic_sin_zec only defined for affine mapping" << endl ;
    abort() ;
  }
  
  Scalar res (*mp) ;
  res.set_etat_qcq() ;
  
  map_affine->sol_elliptic_sin_zec (ope_var, *this, res, freq, 
				    nbr_phase, ampli_min, phase_min) ;

  return (res) ;
}
		    //-----------------------------------//
		    //      General elliptic equation	 //
                    //      fixing the radial derivative //
		    //-----------------------------------//

Scalar Scalar::sol_elliptic_fixe_der_zero (double valeur, 
					   Param_elliptic& ope_var) const {

  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (mp) ;
  
  if (map_affine == 0x0) {
    cout << "sol_elliptic_no_zec only defined for affine mapping" << endl ;
    abort() ;
  }
  
  Scalar res (*mp) ;
  res.set_etat_qcq() ;
  
  map_affine->sol_elliptic_fixe_der_zero (valeur, ope_var, *this, res) ;

  return (res) ;
}

		    //-----------------------------------//
		    //     Two-dimensional Poisson eq.   //
		    //-----------------------------------//

Scalar Scalar::sol_elliptic_2d (Param_elliptic& ope_var) const {

  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (mp) ;
  
  if (map_affine == 0x0) {
    cout << "Poisson 2D only defined for affine mapping" << endl ;
    abort() ;
  }
  
  Scalar res (*mp) ;
  res.set_etat_qcq() ;
  
  map_affine->sol_elliptic_2d(ope_var, *this, res) ;

  return (res) ;
}
		    //-----------------------------------//
		    //     Pseudo-1dimensional eq.   //
		    //-----------------------------------//

Scalar Scalar::sol_elliptic_pseudo_1d (Param_elliptic& ope_var) const {

  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (mp) ;
  
  if (map_affine == 0x0) {
    cout << "Pseudo_1d only defined for affine mapping" << endl ;
    abort() ;
  }
  
  Scalar res (*mp) ;
  res.set_etat_qcq() ;
  
  map_affine->sol_elliptic_pseudo_1d(ope_var, *this, res) ;

  return (res) ;
}
