/*
 * Methods for making interpolations with Godunov-type arrays.
 *
 * See the file tbl_val.h for documentation
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


char TBL_VAL_INTER_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2003/12/19 15:05:15  j_novak
 * Trying to avoid shadowed variables
 *
 * Revision 1.5  2003/10/03 16:17:17  j_novak
 * Corrected some const qualifiers
 *
 * Revision 1.4  2002/11/13 11:22:57  j_novak
 * Version "provisoire" de l'interpolation (sommation depuis la grille
 * spectrale) aux interfaces de la grille de Valence.
 *
 * Revision 1.3  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/11/23 16:03:07  j_novak
 *
 *  minor modifications on the grid check.
 *
 * Revision 1.1  2001/11/22 13:41:54  j_novak
 * Added all source files for manipulating Valencia type objects and making
 * interpolations to and from Meudon grids.
 *
 *
 * $Header$
 *
 */

// headers C
#include <math.h>

// headers Lorene
#include "headcpp.h"
#include "tbl_val.h"


Cmp Tbl_val::to_spectral(const Map& mp, const int lmax, const int lmin, 
			     int type_inter) const {

  assert(etat != ETATNONDEF) ;
  assert( gval->compatible(&mp, lmax, lmin) ) ;
  Cmp resu(&mp) ;

  if (etat == ETATZERO) {
    resu.annule(lmin, lmax-1) ;
    return resu ;
  }
  else {
    int nzin = lmax - lmin ;
    int dim_spec = 1 ;
    const Mg3d* mgrid = mp.get_mg() ;
    for (int i=lmin; i<lmax; i++) {
      if ((mgrid->get_nt(i) > 1)&&(dim_spec==1)) dim_spec = 2; 
      if (mgrid->get_np(i) > 1) dim_spec = 3;
    }
    const Coord& rr = mp.r ;
    
    int* ntet = new int[nzin] ;
    int ntetmax = 1 ;
    int* nphi = new int[nzin] ;
    int nphimax = 1 ;
    for (int i=lmin; i<lmax; i++) {
      int tmp = mgrid->get_np(i) ;
      nphi[i-lmin] = tmp ;
      nphimax = (tmp > nphimax ? tmp : nphimax) ;
      tmp = mgrid->get_nt(i) ;
      ntet[i-lmin] = tmp ;
      ntetmax = (tmp > ntetmax ? tmp : ntetmax) ;
    }
    if (dim_spec > 1) {
      for (int i=lmin; i<lmax; i++) {
	if ((nphimax % nphi[i-lmin]) != 0) {
	cout << "Tbl_val::to_spectral: The numbers of points in phi" << endl ;
	cout << "in the different domains of Meudon grid are not" << endl;
	cout << "well defined; see the documentation." << endl ;
	abort() ;
	}
	assert((ntet[i-lmin]-1) > 0) ;
	if (((ntetmax-1) % (ntet[i-lmin]-1)) != 0) {
	cout <<"Tbl_val::to_spectral: The numbers of points in theta"<< endl ;
	cout << "in the different domains of Meudon grid are not" << endl;
	cout << "well defined; see the documentation." << endl ;
	abort() ;
	}
      }
    }
    
    resu.allocate_all() ;
    if (lmin>0) resu.annule(0,lmin-1) ;
    if (lmax < mgrid->get_nzone()) resu.annule(lmax, mgrid->get_nzone()-1) ;
    
    int fant = gval->get_fantome() ;
    int flag = 1 ; 
    int nrarr = 0 ;
    for (int l=lmin; l<lmax; l++) nrarr += mgrid->get_nr(l) -1 ;
    nrarr++ ;
    switch (dim_spec) {
      
    case 1: {
      int tsize = dim->dim[0] + 2*fant ;
      Tbl fdep(tsize) ; 
      fdep.set_etat_qcq() ;
      for (int i=0; i<tsize; i++) fdep.set(i) = t[i] ;
      Tbl farr(nrarr) ;
      Tbl rarr(nrarr) ;
      rarr.set_etat_qcq() ;
      int inum = 0 ;
      for (int l=lmin; l<lmax; l++) {
	for (int i=0; i<mgrid->get_nr(l); i++) {
	  rarr.set(inum) = (+rr)(l,0,0,i) ;
	  inum++ ;
	}
	inum--;
      }
      farr = gval->interpol1(*gval->zr, rarr, fdep, flag, type_inter) ;
      inum = 0 ;
      for (int l=lmin; l<lmax; l++) {
	for (int i=0; i<mgrid->get_nr(l); i++) {
	  resu.set(l,0,0,i) = farr(inum) ;
	inum++ ;
	}
	inum--;
      }
      break ;
    }

    case 2: {
      int tsizex = dim->dim[1] + 2*fant ;
      int tsizez = dim->dim[0] + 2*fant ;
      Tbl fdep(tsizex, tsizez) ;
      fdep.set_etat_qcq() ;
      for (int j=0; j<tsizex; j++) {
	for (int i=0; i<tsizez; i++) {
	  int l = tsizez*j + i ;
	  fdep.t[l] = t[l] ;
	}
      }
      Tbl farr(ntetmax, nrarr) ;
      Tbl rarr(nrarr) ;
      rarr.set_etat_qcq() ;
      Tbl tetarr(ntetmax) ;
      tetarr.set_etat_qcq() ;
      int ltmax = 0 ;
      int inum = 0 ;
      for (int l=lmin; l<lmax; l++) {
	if (ntetmax == ntet[l-lmin]) ltmax = l ;
	for (int i=0; i<mgrid->get_nr(l); i++) {
	  rarr.set(inum) = (+rr)(l,0,0,i) ;
	  inum++ ;
	}
	inum--;
      }
      const Coord& tet = mp.tet ;
      for (int j=0; j<ntetmax; j++) 
	tetarr.set(j) = (+tet)(ltmax,0,j,0) ;
      farr = gval->interpol2(fdep, rarr, tetarr, type_inter) ;
      inum = 0 ;
      for (int l=lmin; l<lmax; l++) {
	for (int j=0; j<ntet[l-lmin]; j++) {
	  int itet = (ntetmax-1)/(ntet[l-lmin]-1)*j ;
	  for (int i=0; i<mgrid->get_nr(l); i++) {
	    resu.set(l,0,j,i) = farr(itet,inum) ;
	    inum++ ;
	  }
	  inum -= mgrid->get_nr(l) ;
	}
	inum += mgrid->get_nr(l) - 1;
      }
      break ;
    }
    
    case 3: {
      if (type_inter == 0) {
	cout << "The use of routine INSMTS is not well suited" << endl ;
	cout << "for 3D interpolation." << endl ;
	//	abort() ;
      }
      int tsizey = dim->dim[2] + 2*fant ;
      int tsizex = dim->dim[1] + 2*fant ;
      int tsizez = dim->dim[0] + 2*fant ;
      Tbl fdep(tsizey, tsizex, tsizez) ;
      fdep.set_etat_qcq() ;
      for (int k=0; k<tsizey; k++) {
	for (int j=0; j<tsizex; j++) {
	  for (int i=0; i<tsizez; i++) {
	    int l = (k*tsizex+j)*tsizez+i ;
	    fdep.t[l] = t[l];
	  }
	}
      }
      Tbl farr(nphimax, ntetmax, nrarr) ;
      Tbl rarr(nrarr) ;
      rarr.set_etat_qcq() ;
      Tbl tetarr(ntetmax) ;
      tetarr.set_etat_qcq() ;
      Tbl phiarr(nphimax) ;
      phiarr.set_etat_qcq() ;
      int lpmax = 0 ;
      int ltmax = 0 ;
      int inum = 0 ;
      for (int l=lmin; l<lmax; l++) {
	if (ntetmax == ntet[l-lmin]) ltmax = l ;
	if (nphimax == nphi[l-lmin]) lpmax = l ;
	for (int i=0; i<mgrid->get_nr(l); i++) {
	  rarr.set(inum) = (+rr)(l,0,0,i) ;
	  inum++ ;
	}
	inum-- ;
      }
      const Coord& tet = mp.tet ;
      const Coord& phi = mp.phi ;
      for (int k=0; k<nphimax; k++) {
      phiarr.set(k) = (+phi)(lpmax,k,0,0) ;
      }
      for (int j=0; j<ntetmax; j++) 
	tetarr.set(j) = (+tet)(ltmax,0,j,0) ;
      farr = gval->interpol3(fdep, rarr, tetarr, phiarr, type_inter) ;
      inum = 0 ;
      for (int l=lmin; l<lmax; l++) {
	for (int k=0; k<nphi[l-lmin]; k++) {
	  int iphi = (nphimax-1)/(nphi[l-lmin]-1)*k ;
	  for (int j=0; j<ntet[l-lmin]; j++) {
	    int itet = (ntetmax-1)/(ntet[l-lmin]-1)*j ;
	    for (int i=0; i<mgrid->get_nr(l); i++) {
	      resu.set(l,k,j,i) = farr(iphi,itet,inum) ;
	      inum++ ;
	    }
	    inum -= mgrid->get_nr(l) ;
	  }
	}
	inum += mgrid->get_nr(l) - 1 ;
      }
      break ;
    }
    
    default:
      cout << "Tbl_val::to_spectral:Strange error..." << endl ;
      abort() ;
      break ;
      
    }
    
    delete [] ntet ;
    delete [] nphi ;
    return resu ;
  }
}

void Tbl_val::from_spectral(const Cmp& meudon, int lmax, int lmin,
			    bool interfr, bool interft)
{
  assert(meudon.get_etat() != ETATNONDEF) ;
  const Map* mp = meudon.get_mp() ;
  assert( gval->contenue_dans(mp, lmax, lmin) ) ;

  if (meudon.get_etat() == ETATZERO) {
    set_etat_zero() ;
    return ;
  }
  else {
    assert(meudon.get_etat() == ETATQCQ) ;
    set_etat_qcq() ;
    int dim_spec = 1 ;
    const Mg3d* mgrid = mp->get_mg() ;
    for (int i=lmin; i<lmax; i++) {
      if ((mgrid->get_nt(i) > 1)&&(dim_spec==1)) dim_spec = 2; 
      if (mgrid->get_np(i) > 1) dim_spec = 3;
    }
    
    switch (dim_spec) {
      
    case 1: {
      delete [] t ;
      t = gval->somme_spectrale1(meudon) ;
    break ;
    }
    
    case 2: {
      delete [] t ;
      t = gval->somme_spectrale2(meudon) ;
      if (interfr) {
	delete [] tzri ;
	const Gval_spher* gvs = dynamic_cast<const Gval_spher*>(gval) ; //## A modifier
	assert (gvs != 0x0) ;
	tzri = gvs->somme_spectrale2ri(meudon) ;
      }
      if (interft) {
	delete [] txti ;
	const Gval_spher* gvs = dynamic_cast<const Gval_spher*>(gval) ; //## A modifier
	assert (gvs != 0x0) ;
	txti = gvs->somme_spectrale2ti(meudon) ;
      }
      break ;
    }
    
    case 3: {
      delete [] t ;
      t = gval->somme_spectrale3(meudon) ;
      break ;
    }
    
    default:
      cout << "Tbl_val::from_spectral:Strange error..." << endl ;
      abort() ;
      break ;
      
    }    
    return ;
  }
}















