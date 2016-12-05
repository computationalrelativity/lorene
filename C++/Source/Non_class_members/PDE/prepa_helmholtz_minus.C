/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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
 * Revision 1.9  2016/12/05 16:18:10  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:53:30  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:16:10  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2008/07/09 06:51:58  p_grandclement
 * some corrections to helmholtz minus in the nucleus
 *
 * Revision 1.5  2008/07/08 11:45:28  p_grandclement
 * Add helmholtz_minus in the nucleus
 *
 * Revision 1.4  2008/02/18 13:53:43  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.3  2004/08/24 09:14:44  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.2  2004/01/15 09:15:37  p_grandclement
 * Modification and addition of the Helmholtz operators
 *
 * Revision 1.1  2003/12/11 14:48:49  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

//fichiers includes
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"




		//------------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

namespace Lorene {
Matrice _prepa_helmholtz_minus_nondege_pas_prevu(const Matrice &so) {
  
    cout << "Unknown case for prepa_helmholtz_minus_nondege" << endl ;
    abort() ;
    exit(-1) ;
    return so;
}
  
	     	//-------------------
	       //--  R_CHEB   -------
	      //--------------------

Matrice _prepa_helmholtz_minus_nondege_r_cheb (const Matrice &lap) {
    
  int n = lap.get_dim(0) ;
  int non_dege = 2 ;
  
  Matrice res(n-non_dege, n-non_dege) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-non_dege ; i++)
    for (int j=0 ; j<n-non_dege ; j++)
      res.set(i, j) = lap(i, j+non_dege) ;
  
  res.set_band (4,4) ;
  res.set_lu() ;
  return res ;
} 

	     	//-------------------
	       //--  R_CHEBU  -------
	      //--------------------
Matrice _prepa_helmholtz_minus_nondege_r_chebu (const Matrice &lap) {
    
  int n = lap.get_dim(0) ;
  int non_dege = 1 ;
  
  Matrice res(n-non_dege, n-non_dege) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-non_dege ; i++)
    for (int j=0 ; j<n-non_dege ; j++)
      res.set(i, j) = lap(i, j+non_dege) ;
  
  res.set_band (5,3) ;
  res.set_lu() ;
  return res ;
} 
	     	//-------------------
	       //--  R_CHEBP  -------
	      //--------------------
Matrice _prepa_helmholtz_minus_nondege_r_chebp (const Matrice &lap) {
    
  int n = lap.get_dim(0) ;
  int non_dege = 1 ;
  
  Matrice res(n-non_dege, n-non_dege) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-non_dege ; i++)
    for (int j=0 ; j<n-non_dege ; j++)
      res.set(i, j) = lap(i, j+non_dege) ;
  
  res.set_band (4,2) ;
  res.set_lu() ;
  return res ;
} 
	     	//-------------------
	       //--  R_CHEBI  -------
	      //--------------------
Matrice _prepa_helmholtz_minus_nondege_r_chebi (const Matrice &lap) {
    
  int n = lap.get_dim(0) ;
  int non_dege = 1 ;
  
  Matrice res(n-non_dege, n-non_dege) ;
  res.set_etat_qcq() ;
  for (int i=0 ; i<n-non_dege ; i++)
    for (int j=0 ; j<n-non_dege ; j++)
      res.set(i, j) = lap(i, j+non_dege) ;
  
  res.set_band (4,2) ;
  res.set_lu() ;
  return res ;
} 

        	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
Matrice prepa_helmholtz_minus_nondege(const Matrice &ope, int base_r) {

  // Routines de derivation
  static Matrice (*prepa_helmholtz_minus_nondege[MAX_BASE])
    (const Matrice&) ;
  static int nap = 0 ;
  
  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      prepa_helmholtz_minus_nondege[i] = 
	_prepa_helmholtz_minus_nondege_pas_prevu ;
    }
    // Les routines existantes
    prepa_helmholtz_minus_nondege[R_CHEB >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_cheb ;
    prepa_helmholtz_minus_nondege[R_CHEBU >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_chebu ; 
    prepa_helmholtz_minus_nondege[R_CHEBP >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_chebp ;
    prepa_helmholtz_minus_nondege[R_CHEBI >> TRA_R] = 
      _prepa_helmholtz_minus_nondege_r_chebi ;
  }
  
  Matrice res(prepa_helmholtz_minus_nondege[base_r](ope)) ;
  return res ;
}

}
