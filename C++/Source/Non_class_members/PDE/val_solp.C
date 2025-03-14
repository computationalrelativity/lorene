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
 * Revision 1.7  2016/12/05 16:18:10  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:31  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:16:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/02/18 13:53:45  j_novak
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
 * Revision 1.2  2003/12/11 15:37:09  p_grandclement
 * sqrt(2) ----> sqrt(double(2))
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

#include "proto.h"
#include "matrice.h"
#include "type_parite.h"


		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
namespace Lorene {
Tbl _val_solp_pas_prevu (const Tbl&, double) {

    cout << " Base_r unknown in val_solp."<< endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _val_solp_r_cheb (const Tbl& sp, double alpha) {
  
  int nr = sp.get_dim(0) ;
  Tbl res(4) ;
  res.annule_hard() ;
  
  // Solution en + 1 
  for (int i=0 ; i<nr ; i++)
    res.set(0) += sp(i) ;

  // Solution en -1 :
  for (int i=0 ; i<nr ; i++)
    if (i%2 == 0)
      res.set(1) += sp(i) ;
    else
      res.set(1) -= sp(i) ;

  // Derivee en +1 :
  for (int i=0 ; i<nr ; i++)
    res.set(2) += sp(i)*i*i/alpha ;

  // Derivee en -1 :
  for (int i=0 ; i<nr ; i++)
    if (i%2 == 0)
      res.set(3) -= sp(i)*i*i/alpha ;
    else
      res.set(3) += sp(i)*i*i/alpha ;

  res /= sqrt(double(2)) ;
  return res ;
}	
	
		//-------------------
	       //--  R_CHEBP  ------
	      //-------------------

Tbl _val_solp_r_chebp (const Tbl& sp, double alpha) {
  
  int nr = sp.get_dim(0) ;
  Tbl res(4) ;
  res.annule_hard() ;
  
  // Solution en +1 :
  for (int i=0 ; i<nr ; i++)
    res.set(0) += sp(i) ;

  // Solution en 0 (a priori pas trop utilise)
  for (int i=0 ; i<nr ; i++)
    if (i%2==0)
      res.set(1) += sp(i) ;
    else
      res.set(1) -= sp(i) ;
  
  // Derivee en +1 :
  for (int i=0 ; i<nr ; i++) 
    res.set(2) += sp(i)*(2*i)*(2*i)/alpha ;

  // Derivee en 0
  res.set(3) = 0 ;

  res /= sqrt(double(2)) ;
  return res ;
}
	
	
	      	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	
Tbl _val_solp_r_chebi (const Tbl& sp, double alpha) {
     
  int nr = sp.get_dim(0) ;
  Tbl res(4) ;
  res.annule_hard() ;
  
  // Solution en +1 :
  for (int i=0 ; i<nr ; i++)
    res.set(0) += sp(i) ;

  // Solution en 0 :
  res.set(1) = 0 ;

  // Derivee en +1 :
  for (int i=0 ; i<nr ; i++) 
    res.set(2) += sp(i)*(2*i+1)*(2*i+1)/alpha ;
  
  // Derivee en 0 :
  for (int i=0 ; i<nr ; i++)
    if (i%2==0)
      res.set(3) += (2*i+1)*sp(i) ;
    else
      res.set(3) -= (2*i+1)*sp(i) ;

  res /= sqrt(double(2)) ;
  return res ;   
}
	
	
	
	       	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
	
Tbl _val_solp_r_chebu (const Tbl& sp, double alpha) {
 
  int nr = sp.get_dim(0) ;
  Tbl res(4) ;
  res.annule_hard() ;

  // Solution en + 1 
  for (int i=0 ; i<nr ; i++)
    res.set(0) += sp(i) ;

  // Solution en -1 :
  for (int i=0 ; i<nr ; i++)
    if (i%2==0)
      res.set(1) += sp(i) ;
    else
      res.set(1) -= sp(i) ;

  // Derivee en +1 c'est zero ca !

  // Derivee en -1 :
  for (int i=0 ; i<nr ; i++)
    if (i%2==0)
      res.set(3) += 4.*alpha*i*i*sp(i) ;
    else
      res.set(3) -= 4.*alpha*i*i*sp(i) ;
 
  res /= sqrt(double(2)) ;
  return res ;
}
	
	
	
	
	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl val_solp (const Tbl& sp, double alpha, int base_r) {

		// Routines de derivation
    static Tbl (*val_solp[MAX_BASE])(const Tbl&, double) ;
    static int nap = 0 ;
    
    // Premier appel
    if (nap==0) {
      nap = 1 ;
      for (int i=0 ; i<MAX_BASE ; i++) {
	val_solp[i] = _val_solp_pas_prevu ;
      }
      // Les routines existantes
      val_solp[R_CHEB >> TRA_R] = _val_solp_r_cheb ;
      val_solp[R_CHEBU >> TRA_R] = _val_solp_r_chebu ;
      val_solp[R_CHEBP >> TRA_R] = _val_solp_r_chebp ;
      val_solp[R_CHEBI >> TRA_R] = _val_solp_r_chebi ;
    }
    
    Tbl res(val_solp[base_r](sp, alpha)) ;
    return res ;
}
}
