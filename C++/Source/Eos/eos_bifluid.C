/*
 * Methods of the class Eos_bifluid.
 *
 * (see file eos_bifluid.h for documentation).
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


char eos_bifluid_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2003/10/03 15:58:46  j_novak
 * Cleaning of some headers
 *
 * Revision 1.7  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.6  2002/05/31 16:13:36  j_novak
 * better inversion for eos_bifluid
 *
 * Revision 1.5  2002/05/10 09:55:27  j_novak
 * *** empty log message ***
 *
 * Revision 1.4  2002/05/10 09:26:52  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 * Revision 1.3  2002/01/03 15:30:27  j_novak
 * Some comments modified.
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.5  2001/10/10  13:49:53  eric
 * Modif Joachim: &(Eos_bifluid::...) --> &Eos_bifluid::...
 *  pour conformite au compilateur HP.
 *
 * Revision 1.4  2001/08/31  15:48:11  novak
 * The flag tronc has been added to nbar_ent... functions
 *
 * Revision 1.3  2001/08/27 12:23:40  novak
 * The Cmp arguments delta2 put to const
 *
 * Revision 1.2  2001/08/27 09:52:49  novak
 * Use of new variable delta2
 *
 * Revision 1.1  2001/06/21 15:21:47  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Headers Lorene
#include "eos_bifluid.h"
#include "cmp.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//


// Standard constructor without name
// ---------------------------------
Eos_bifluid::Eos_bifluid(){
        
    set_name("") ; 
    
}

// Standard constructor with name
// ---------------------------------
Eos_bifluid::Eos_bifluid(const char* name_i){
        
    set_name(name_i) ; 
    
}

// Copy constructor
// ----------------
Eos_bifluid::Eos_bifluid(const Eos_bifluid& eos_i){
        
    set_name(eos_i.name) ; 
    
}

// Constructor from a binary file
// ------------------------------
Eos_bifluid::Eos_bifluid(FILE* fich){
        
    fread(name, sizeof(char), 100, fich) ;		
    
}

// Constructor from a formatted file
// ---------------------------------
Eos_bifluid::Eos_bifluid(ifstream& fich){
        
    fich.getline(name, 100) ;
    
}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_bifluid::~Eos_bifluid(){
    
    // does nothing
        
}

			//-------------------------//
			//  Manipulation of name   //
			//-------------------------//
			
			
void Eos_bifluid::set_name(const char* name_i) {

    strncpy(name, name_i,  100) ; 
    
}

const char* Eos_bifluid::get_name() const {
    
    return name ; 
    
}

			//------------//
			//  Outputs   //
			//------------//

void Eos_bifluid::sauve(FILE* fich) const {

    int ident = identify() ; 
    fwrite_be(&ident, sizeof(int), 1, fich) ;	
    	
    fwrite(name, sizeof(char), 100, fich) ;		
   
}
    



ostream& operator<<(ostream& ost, const Eos_bifluid& eqetat)  {
    ost << eqetat.get_name() << endl ; 
    eqetat >> ost ;
    return ost ;
}


			//-------------------------------//
			//    Computational routines     //
			//-------------------------------//

// Complete computational routine giving all thermo variables
//-----------------------------------------------------------

void Eos_bifluid::calcule_tout(const Cmp& ent1, const Cmp& ent2, 
			       const Cmp& delta2, Cmp& nbar1, Cmp& nbar2,  
			       Cmp& ener, Cmp& press, 
			       int nzet, int l_min) const {
  
  assert(ent1.get_etat() != ETATNONDEF) ; 
  assert(ent2.get_etat() != ETATNONDEF) ; 
  assert(delta2.get_etat() != ETATNONDEF) ;
  
  const Map* mp = ent1.get_mp() ;	// Mapping
  assert(mp == ent2.get_mp()) ;
  assert(mp == delta2.get_mp()) ;
  assert(mp == ener.get_mp()) ;
  
  if ((ent1.get_etat() == ETATZERO)&&(ent2.get_etat() == ETATZERO)) {
    nbar1.set_etat_zero() ; 
    nbar2.set_etat_zero() ; 
    ener.set_etat_zero() ; 
    press.set_etat_zero() ; 
    return ; 
  }
  nbar1.allocate_all() ;
  nbar2.allocate_all() ;
  ener.allocate_all() ;
  press.allocate_all() ;

  const Mg3d* mg = mp->get_mg() ;	// Multi-grid
  
  int nz = mg->get_nzone() ;		// total number of domains
  
  // resu is set to zero in the other domains :
  
  if (l_min > 0) {
    nbar1.annule(0, l_min-1) ; 
    nbar2.annule(0, l_min-1) ; 
    ener.annule(0, l_min-1) ; 
    press.annule(0, l_min-1) ; 
  }
  
  if (l_min + nzet < nz) {
    nbar1.annule(l_min + nzet, nz - 1) ; 
    nbar2.annule(l_min + nzet, nz - 1) ; 
    ener.annule(l_min + nzet, nz - 1) ; 
    press.annule(l_min + nzet, nz - 1) ; 
  }

  int np0 = mp->get_mg()->get_np(0) ;
  int nt0 = mp->get_mg()->get_nt(0) ;
  for (int l=l_min; l<l_min+nzet; l++) {
    assert(mp->get_mg()->get_np(l) == np0) ;
    assert(mp->get_mg()->get_nt(l) == nt0) ; 
  }

  for (int k=0; k<np0; k++) {
    for (int j=0; j<nt0; j++) {
      bool inside = true ;
      bool inside1 = true ;
      bool inside2 = true ;
      for (int l=l_min; l< l_min + nzet; l++) {
	for (int i=0; i<mp->get_mg()->get_nr(l); i++) {
	  double xx, xx1, xx2, n1, n2 ;
	  xx1 = ent1(l,k,j,i) ; 
	  xx2 = ent2(l,k,j,i) ; 
	  xx = delta2(l,k,j,i) ;
	  if (inside) {
	    inside = (!nbar_ent_p(xx1, xx2, xx, n1, n2)) ; 
	    inside1 = ((n1 > 0.)&&(xx1>0.)) ;
	    inside2 = ((n2 > 0.)&&(xx2>0.)) ;
	  }
	  if (inside) {
	    nbar1.set(l,k,j,i) = n1 ;
	    nbar2.set(l,k,j,i) = n2 ;
	    ener.set(l,k,j,i) = ener_nbar_p(n1, n2, xx) ;
	    press.set(l,k,j,i) = press_nbar_p(n1, n2, xx) ;
	  }
	  else {
	    if (inside1) {
	      n1 = nbar_ent_p1(xx1) ;
	      inside1 = (n1 > 0.) ;
	    }
	    if (!inside1) n1 = 0. ;
	    if (inside2) {
	      n2 = nbar_ent_p2(xx2) ;
	      inside2 = (n2 > 0.) ;
	    }
	    if (!inside2) n2 = 0. ;
	    nbar1.set(l,k,j,i) = n1 ;
	    nbar2.set(l,k,j,i) = n2 ;
	    
	    ener.set(l,k,j,i) = ener_nbar_p(n1, n2, 0.) ;
	    press.set(l,k,j,i) = press_nbar_p(n1, n2, 0. ) ;
	  }
	}
      }
    }
  }

}

// Baryon density from enthalpy 
//------------------------------

void Eos_bifluid::nbar_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2,
			   Cmp& nbar1, Cmp& nbar2, int nzet, int l_min) const {
  
  assert(ent1.get_etat() != ETATNONDEF) ; 
  assert(ent2.get_etat() != ETATNONDEF) ; 
  assert(delta2.get_etat() != ETATNONDEF) ;
  
  const Map* mp = ent1.get_mp() ;	// Mapping
  assert(mp == ent2.get_mp()) ;
  assert(mp == delta2.get_mp()) ;
  assert(mp == nbar1.get_mp()) ;
  
  if ((ent1.get_etat() == ETATZERO)&&(ent2.get_etat() == ETATZERO)) {
    nbar1.set_etat_zero() ; 
    nbar2.set_etat_zero() ; 
    return ; 
  }
  nbar1.allocate_all() ;
  nbar2.allocate_all() ;

  const Mg3d* mg = mp->get_mg() ;	// Multi-grid
  
  int nz = mg->get_nzone() ;		// total number of domains
  
  // resu is set to zero in the other domains :
  
  if (l_min > 0) {
    nbar1.annule(0, l_min-1) ; 
    nbar2.annule(0, l_min-1) ; 
  }
  
  if (l_min + nzet < nz) {
    nbar1.annule(l_min + nzet, nz - 1) ; 
    nbar2.annule(l_min + nzet, nz - 1) ; 
  }

  int np0 = mp->get_mg()->get_np(0) ;
  int nt0 = mp->get_mg()->get_nt(0) ;
  for (int l=l_min; l<l_min+nzet; l++) {
    assert(mp->get_mg()->get_np(l) == np0) ;
    assert(mp->get_mg()->get_nt(l) == nt0) ; 
  }

  for (int k=0; k<np0; k++) {
    for (int j=0; j<nt0; j++) {
      bool inside = true ;
      bool inside1 = true ;
      bool inside2 = true ;
      for (int l=l_min; l< l_min + nzet; l++) {
	for (int i=0; i<mp->get_mg()->get_nr(l); i++) {
	  double xx, xx1, xx2, n1, n2 ;
	  xx1 = ent1(l,k,j,i) ; 
	  xx2 = ent2(l,k,j,i) ; 
	  xx = delta2(l,k,j,i) ;
	  if (inside) {
	    inside = (!nbar_ent_p(xx1, xx2, xx, n1, n2)) ; 
	    inside1 = ((n1 > 0.)&&(xx1>0.)) ;
	    inside2 = ((n2 > 0.)&&(xx2>0.)) ;
	  }
	  if (inside) {
	    nbar1.set(l,k,j,i) = n1 ;
	    nbar2.set(l,k,j,i) = n2 ;
	  }
	  else {
	    if (inside1) {
	      n1 = nbar_ent_p1(xx1) ;
	      inside1 = (n1 > 0.) ;
	    }
	    if (!inside1) n1 = 0. ;
	    if (inside2) {
	      n2 = nbar_ent_p2(xx2) ;
	      inside2 = (n2 > 0.) ;
	    }
	    if (!inside2) n2 = 0. ;
	    nbar1.set(l,k,j,i) = n1 ;
	    nbar2.set(l,k,j,i) = n2 ;
	  }
	}
      }
    }
  }

}


// Energy density from enthalpy 
//------------------------------

Cmp Eos_bifluid::ener_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
			  int nzet, int l_min) const {
    
  Cmp ener(ent1.get_mp()) ; 
  assert(ent1.get_etat() != ETATNONDEF) ; 
  assert(ent2.get_etat() != ETATNONDEF) ; 
  assert(delta2.get_etat() != ETATNONDEF) ;
    
  const Map* mp = ent1.get_mp() ;	// Mapping
  assert(mp == ent2.get_mp()) ;
  assert(mp == delta2.get_mp()) ;
    
  if ((ent1.get_etat() == ETATZERO)&&(ent2.get_etat() == ETATZERO)) {
    ener.set_etat_zero() ; 
    return ener; 
  }
  
  ener.allocate_all() ;

  const Mg3d* mg = mp->get_mg() ;	// Multi-grid
  
  int nz = mg->get_nzone() ;		// total number of domains
  
  // resu is set to zero in the other domains :
  
  if (l_min > 0) 
    ener.annule(0, l_min-1) ; 
  
  if (l_min + nzet < nz) 
    ener.annule(l_min + nzet, nz - 1) ; 

  int np0 = mp->get_mg()->get_np(0) ;
  int nt0 = mp->get_mg()->get_nt(0) ;
  for (int l=l_min; l<l_min+nzet; l++) {
    assert(mp->get_mg()->get_np(l) == np0) ;
    assert(mp->get_mg()->get_nt(l) == nt0) ; 
  }

  for (int k=0; k<np0; k++) {
    for (int j=0; j<nt0; j++) {
      bool inside = true ;
      bool inside1 = true ;
      bool inside2 = true ;
      for (int l=l_min; l< l_min + nzet; l++) {
	for (int i=0; i<mp->get_mg()->get_nr(l); i++) {
	  double xx, xx1, xx2, n1, n2 ;
	  xx1 = ent1(l,k,j,i) ; 
	  xx2 = ent2(l,k,j,i) ; 
	  xx = delta2(l,k,j,i) ;
	  if (inside) {
	    inside = (!nbar_ent_p(xx1, xx2, xx, n1, n2)) ; 
	    inside1 = ((n1 > 0.)&&(xx1>0.)) ;
	    inside2 = ((n2 > 0.)&&(xx2>0.)) ;
	  }
	  if (inside) {
	    ener.set(l,k,j,i) = ener_nbar_p(n1, n2, xx) ;
	  }
	  else {
	    if (inside1) {
	      n1 = nbar_ent_p1(xx1) ;
	      inside1 = (n1 > 0.) ;
	    }
	    if (!inside1) n1 = 0. ;
	    if (inside2) {
	      n2 = nbar_ent_p2(xx2) ;
	      inside2 = (n2 > 0.) ;
	    }
	    if (!inside2) n2 = 0. ;
	    ener.set(l,k,j,i) = ener_nbar_p(n1, n2, 0.) ;
	  }
	}
      }
    }
  }
  return ener ;
}

// Pressure from enthalpies 
//-------------------------

Cmp Eos_bifluid::press_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
		      int nzet, int l_min ) const {
    
  Cmp press(ent1.get_mp()) ; 
  assert(ent1.get_etat() != ETATNONDEF) ; 
  assert(ent2.get_etat() != ETATNONDEF) ; 
  assert(delta2.get_etat() != ETATNONDEF) ;
    
  const Map* mp = ent1.get_mp() ;	// Mapping
  assert(mp == ent2.get_mp()) ;
    
  if ((ent1.get_etat() == ETATZERO)&&(ent2.get_etat() == ETATZERO)) {
    press.set_etat_zero() ; 
    return press; 
  }
  press.allocate_all() ;

  const Mg3d* mg = mp->get_mg() ;	// Multi-grid
  
  int nz = mg->get_nzone() ;		// total number of domains
  
  // resu is set to zero in the other domains :
  
  if (l_min > 0) 
    press.annule(0, l_min-1) ; 
  
  if (l_min + nzet < nz) 
    press.annule(l_min + nzet, nz - 1) ; 

  int np0 = mp->get_mg()->get_np(0) ;
  int nt0 = mp->get_mg()->get_nt(0) ;
  for (int l=l_min; l<l_min+nzet; l++) {
    assert(mp->get_mg()->get_np(l) == np0) ;
    assert(mp->get_mg()->get_nt(l) == nt0) ; 
  }

  for (int k=0; k<np0; k++) {
    for (int j=0; j<nt0; j++) {
      bool inside = true ;
      bool inside1 = true ;
      bool inside2 = true ;
      for (int l=l_min; l< l_min + nzet; l++) {
	for (int i=0; i<mp->get_mg()->get_nr(l); i++) {
	  double xx, xx1, xx2, n1, n2 ;
	  xx1 = ent1(l,k,j,i) ; 
	  xx2 = ent2(l,k,j,i) ; 
	  xx = delta2(l,k,j,i) ;
	  if (inside) {
	    inside = (!nbar_ent_p(xx1, xx2, xx, n1, n2)) ; 
	    inside1 = ((n1 > 0.)&&(xx1>0.)) ;
	    inside2 = ((n2 > 0.)&&(xx2>0.)) ;
	  }
	  if (inside) 
	    press.set(l,k,j,i) = press_nbar_p(n1, n2, xx) ;
	  else {
	    if (inside1) {
	      n1 = nbar_ent_p1(xx1) ;
	      inside1 = (n1 > 0.) ;
	    }
	    if (!inside1) n1 = 0. ;
	    if (inside2) {
	      n2 = nbar_ent_p2(xx2) ;
	      inside2 = (n2 > 0.) ;
	    }
	    if (!inside2) n2 = 0. ;
	    press.set(l,k,j,i) = press_nbar_p(n1, n2, 0. ) ;
	  }
	}
      }
    }
  }
  return press ;
}

// Generic computational routine for get_Kxx
//------------------------------------------

void Eos_bifluid::calcule(const Cmp& nbar1, const Cmp& nbar2, const Cmp&
		     x2, int nzet, int l_min, double
		     (Eos_bifluid::*fait)(double, double, double) const, 
			  Cmp& resu) const {
    
    assert(nbar1.get_etat() != ETATNONDEF) ; 
    assert(nbar2.get_etat() != ETATNONDEF) ; 
    assert(x2.get_etat() != ETATNONDEF) ; 
    
    const Map* mp = nbar1.get_mp() ;	// Mapping
    assert(mp == nbar2.get_mp()) ;
    
    
    if ((nbar1.get_etat() == ETATZERO)&&(nbar2.get_etat() == ETATZERO)) {
	resu.set_etat_zero() ; 
	return ; 
    }
    
    bool nb1 = nbar1.get_etat() == ETATQCQ ; 
    bool nb2 = nbar2.get_etat() == ETATQCQ ; 
    bool bx2 = x2.get_etat() == ETATQCQ ; 
    const Valeur* vnbar1 = 0x0 ;
    const Valeur* vnbar2 = 0x0 ;
    const Valeur* vx2 = 0x0 ;
    if (nb1) { vnbar1 = &nbar1.va ;
    vnbar1->coef_i() ; }
    if (nb2) { vnbar2 = &nbar2.va ;
    vnbar2->coef_i() ; }
    if (bx2) {vx2 = & x2.va ;
    vx2->coef_i() ; }
   
    const Mg3d* mg = mp->get_mg() ;	// Multi-grid
    
    int nz = mg->get_nzone() ;		// total number of domains
    
    // Preparations for a point by point computation:
    resu.set_etat_qcq() ;
    Valeur& vresu = resu.va ; 
    vresu.set_etat_c_qcq() ;
    vresu.c->set_etat_qcq() ;

    // Loop on domains where the computation has to be done :
    for (int l = l_min; l< l_min + nzet; l++) {
	
      assert(l>=0) ; 
      assert(l<nz) ; 
      
      Tbl* tnbar1 = 0x0 ;
      Tbl* tnbar2 = 0x0 ;
      Tbl* tx2 = 0x0 ;
      
      if (nb1) tnbar1 = vnbar1->c->t[l] ;
      if (nb2) tnbar2 = vnbar2->c->t[l] ;
      if (bx2) tx2 = vx2->c->t[l] ;
      Tbl* tresu = vresu.c->t[l] ; 
      
      bool nb1b = false ;
      if (nb1) nb1b = tnbar1->get_etat() == ETATQCQ ;
      bool nb2b = false ;
      if (nb2) nb2b = tnbar2->get_etat() == ETATQCQ ;
      bool bx2b = false ;
      if (bx2) bx2b = tx2->get_etat() == ETATQCQ ;
      tresu->set_etat_qcq() ;
      
      double n1, n2, xx ;
      
      for (int i=0; i<tnbar1->get_taille(); i++) {
	
	n1 = nb1b ? tnbar1->t[i] : 0 ;
	n2 = nb2b ? tnbar2->t[i] : 0 ;
	xx = bx2b ? tx2->t[i] : 0 ;
	tresu->t[i] = (this->*fait)(n1, n2, xx ) ;
      }  
      
    }  // End of the loop on domains where the computation had to be done
    
    // resu is set to zero in the other domains :
    
    if (l_min > 0) {
      resu.annule(0, l_min-1) ; 
    }
    
    if (l_min + nzet < nz) {
      resu.annule(l_min + nzet, nz - 1) ; 
    }
}

Cmp Eos_bifluid::get_Knn(const Cmp& nbar1, const Cmp& nbar2, const Cmp&
			 delta2, int nzet, int l_min) const {
    
    Cmp resu(nbar1.get_mp()) ; 
    
    calcule(nbar1, nbar2, delta2, nzet, l_min, &Eos_bifluid::get_K11, resu) ;
    
    return resu ; 
    
}

Cmp Eos_bifluid::get_Knp(const Cmp& nbar1, const Cmp& nbar2, const Cmp& 
			 delta2, int nzet, int l_min) const {
    
    Cmp resu(delta2.get_mp()) ; 
    
    calcule(nbar1, nbar2, delta2, nzet, l_min, &Eos_bifluid::get_K12, resu) ;
    
    return resu ; 
    
}

Cmp Eos_bifluid::get_Kpp(const Cmp& nbar1, const Cmp& nbar2, const Cmp& 
			 delta2, int nzet, int l_min) const {
    
    Cmp resu(nbar2.get_mp()) ; 
    
    calcule(nbar1, nbar2, delta2, nzet, l_min, &Eos_bifluid::get_K22, resu) ;
    
    return resu ; 
    
}
