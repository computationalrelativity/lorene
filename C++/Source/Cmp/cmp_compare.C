/*
 *  File to compare two Cmp
 *
 */

/*
 *   Copyright (c) 2003 Francois Limousin
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char cmp_compare_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/06/20 15:04:07  f_limousin
 * *** empty log message ***
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 * template files
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include "math.h"

// Lorene headers
#include "cmp.h"


void Cmp::compare(const Cmp& comp, const char* name, int ii, int jj) {

  int nzone = (*((*mp).get_mg())).get_nzone() ;
  double nr = (*((*mp).get_mg())).get_nr(0) ;
  double nt = (*((*mp).get_mg())).get_nt(0) ;
  double np = (*((*mp).get_mg())).get_np(0) ;
 

  double diff_moy[nzone], diff_max[nzone] ;
  double compt[nzone], tmp ;

  for (int l=0; l<nzone; l++) {
    diff_moy[l] = 0 ;
    diff_max[l] = 0 ;
    compt[l] = 0 ;
  }

  for (int i=0; i<nr; i++) {
    for (int j=0; j<nt; j++) {
      for (int k=0; k<np; k++) {
	for (int l=0; l<nzone; l++) {
	  
	  if (fabs((*this)(l,k,j,i))>1.e-10) {
	    tmp = fabs(((*this)(l,k,j,i) - comp(l,k,j,i))
		       /(*this)(l,k,j,i)) ;
	    diff_moy[l] += tmp  ;
	    compt[l]++ ;
	    if (diff_max[l]<tmp) {
	      diff_max[l] = tmp ;
	    }
	  }
	}
      }
    }
  }


  for (int l=0; l<nzone; l++) {
    if (compt[l] != 0) {
    diff_moy[l] /= compt[l] ;
    }
  }


  if ( ii == -1 && jj == -1 ) {
    for (int l=0; l<nzone; l++) {
      if (compt[l] != 0) {    
	cout << "moy relative difference for " << name << " in zone " << l 
	     << " = " << diff_moy[l] << endl ;
	cout << "max relative difference for " << name << " in zone " << l 
	     << " = " << diff_max[l] << endl << endl ;
      }
      else {
	cout << name << " in zone " << l << " < 1.e-10" << endl << endl ;
      }
    }
  }

  if ( ii != -1 && jj == -1 ) {
    for (int l=0; l<nzone; l++) {
      if (compt[l] != 0) {    
	cout << "moy relative difference for " << name << "(" << ii 
	     << ") in zone " << l << " = " << diff_moy[l] << endl ;
	cout << "max relative difference for " << name << "(" << ii
	     << ") in zone " << l << " = " << diff_max[l] << endl << endl ;
      }
      else {
	cout << name << "(" << ii << ") in zone " << l << " < 1.e-10" 
	     << endl << endl ;
      }
    }
  }

  if ( ii != -1 && jj != -1 ) {
    for (int l=0; l<nzone; l++) {
      if (compt[l] != 0) {    
	cout << "moy relative difference for " << name << "(" << ii << "," 
	     << jj << ") in zone " << l << " = " << diff_moy[l] << endl ;
	cout << "max relative difference for " << name << "(" << ii << "," 
	     << jj << ") in zone " << l << " = " << diff_max[l] << endl << endl ;
      }
      else {
	cout << name << "(" << ii << "," << jj << ") in zone " << l 
	     << " < 1.e-10" << endl << endl ;
      }
    }
  }

}


void Cmp::compare(FILE* fich, const char* name_i) {


  Mg3d mg(fich) ;
  Map_et mp(mg, fich) ;

  Cmp comp(mp, mg, fich) ;

  cout << "----------------------------------------" << endl ;
  cout << "Comparison of " << name_i << endl << endl ;

  compare(comp, name_i) ;

}
