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


 

/*
 * $Id$
 * $Log$
 * Revision 1.23  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.22  2015/06/12 12:38:24  j_novak
 * Implementation of the corrected formula for the quadrupole momentum.
 *
 * Revision 1.21  2015/06/11 14:41:59  a_sourie
 * Corrected minor bug
 *
 * Revision 1.20  2015/06/11 13:50:19  j_novak
 * Minor corrections
 *
 * Revision 1.19  2015/06/10 14:39:17  a_sourie
 * New class Eos_bf_tabul for tabulated 2-fluid EoSs and associated functions for the computation of rotating stars with such EoSs.
 *
 * Revision 1.18  2014/10/13 08:52:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.17  2014/04/25 10:43:51  j_novak
 * The member 'name' is of type string now. Correction of a few const-related issues.
 *
 * Revision 1.16  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.15  2006/03/10 08:38:55  j_novak
 * Use of C++-style casts.
 *
 * Revision 1.14  2004/09/01 16:12:30  r_prix
 * (hopefully) fixed seg-fault bug with my inconsistent treatment of eos-bifluid 'name'
 * (was char-array, now char*)
 *
 * Revision 1.13  2004/09/01 09:50:34  r_prix
 * adapted to change in read_variable() for strings
 *
 * Revision 1.12  2003/12/17 23:12:32  r_prix
 * replaced use of C++ <string> by standard ANSI char* to be backwards compatible
 * with broken compilers like MIPSpro Compiler 7.2 on SGI Origin200. ;-)
 *
 * Revision 1.11  2003/12/05 15:09:47  r_prix
 * adapted Eos_bifluid class and subclasses to use read_variable() for
 * (formatted) file-reading.
 *
 * Revision 1.10  2003/12/04 14:17:26  r_prix
 * new 2-fluid EOS subtype 'typeos=5': this is identical to typeos=0
 * (analytic EOS), but we perform the EOS inversion "slow-rot-style",
 * i.e. we don't switch to a 1-fluid EOS when one fluid vanishes.
 *
 * Revision 1.9  2003/11/18 18:28:38  r_prix
 * moved particle-masses m_1, m_2 of the two fluids into class eos_bifluid (from eos_bf_poly)
 *
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
 * The flag tronc has been added to ar_ent... functions
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
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "eos_bifluid.h"
#include "cmp.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//


// Standard constructor without name
// ---------------------------------
namespace Lorene {
Eos_bifluid::Eos_bifluid() :
  name("EoS bi-fluid"), m_1(1), m_2(1)
{ }

// Standard constructor with name and masses
// ---------------------------------
Eos_bifluid::Eos_bifluid(const char* name_i, double mass1, double mass2) :
  name(name_i), m_1(mass1), m_2(mass2)
{ }

// Copy constructor
// ----------------
Eos_bifluid::Eos_bifluid(const Eos_bifluid& eos_i) :
  name(eos_i.name), m_1(eos_i.m_1), m_2(eos_i.m_2)
{ }

// Constructor from a binary file
// ------------------------------
Eos_bifluid::Eos_bifluid(FILE* fich)
{
  char dummy [MAX_EOSNAME];
  if (fread(dummy, sizeof(char),MAX_EOSNAME, fich) == 0) 
    cerr << "Reading problem in " << __FILE__ << endl ;
  name = dummy ;
  fread_be(&m_1, sizeof(double), 1, fich) ;		
  fread_be(&m_2, sizeof(double), 1, fich) ;	
    
}

// Constructor from a formatted file
// ---------------------------------
Eos_bifluid::Eos_bifluid(const char *fname)
{
  int res=0;
  char* dummy = 0x0 ;
  char* char_name = const_cast<char*>("name");
  char* char_m1 = const_cast<char*>("m_1") ;
  char* char_m2 = const_cast<char*>("m_2") ;
  res += read_variable (fname, char_name, &dummy);
  res += read_variable (fname, char_m1, m_1);
  res += read_variable (fname, char_m2, m_2);

  name = dummy ;

  free(dummy) ;

  if (res)
    {
      cerr << "ERROR: Eos_bifluid(const char*): could not read either of \
the variables 'name', 'm_1, or 'm_2' from file " << fname << endl;
      exit (-1);
    }

}

// Constructor from a formatted file
// ---------------------------------
  Eos_bifluid::Eos_bifluid(ifstream& fich){
      
  string aname ;

  // EOS identificator : 
  char blabla[80] ;
  fich >> aname; 
  name = aname ;
  fich.getline(blabla, 80) ;

  fich >> m_1 ; fich.getline(blabla, 80) ;
  fich >> m_2 ; fich.getline(blabla, 80) ;
}




			//--------------//
			//  Destructor  //
			//--------------//

Eos_bifluid::~Eos_bifluid()
{ }

			//--------------//
			//  Assignment  //
			//--------------//
void Eos_bifluid::operator=(const Eos_bifluid& eosi) {
    
    name = eosi.name ; 
    m_1 = eosi.m_1;
    m_2 = eosi.m_2;
    
}


			//------------//
			//  Outputs   //
			//------------//

void Eos_bifluid::sauve(FILE* fich) const 
{
  assert(name.size() < MAX_EOSNAME) ;
  char dummy [MAX_EOSNAME];
  int ident = identify() ; 

  fwrite_be(&ident, sizeof(int), 1, fich) ;	

  strncpy (dummy, name.c_str(), MAX_EOSNAME);
  dummy[MAX_EOSNAME-1] = 0;
  if (fwrite(dummy, sizeof(char), MAX_EOSNAME, fich ) == 0)
    cerr << "Writing problem in " << __FILE__ << endl ;

  fwrite_be(&m_1, sizeof(double), 1, fich) ;	
  fwrite_be(&m_2, sizeof(double), 1, fich) ;	
   
}
    

ostream& operator<<(ostream& ost, const Eos_bifluid& eqetat)  {
    ost << eqetat.get_name() << endl ; 
    ost << "   Mean particle 1 mass : " << eqetat.get_m1() << " m_B" << endl ;
    ost << "   Mean particle 2 mass : " << eqetat.get_m2() << " m_B" << endl ;

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

  //**********************************************************************
  //RP: for comparison with slow-rotation, we might have to treat the 
  // 1-fluid region somewhat differently...
  bool slow_rot_style = false;  // off by default

  if ( identify() == 2 )  // only applies if newtonian 2-fluid polytrope
    {
      const Eos_bf_poly_newt *this_eos = dynamic_cast<const Eos_bf_poly_newt*>(this);
      if (this_eos -> get_typeos() == 5)
	{
	  slow_rot_style = true;
	  cout << "DEBUG: using EOS-inversion slow-rot-style!! " << endl;
	}
    }

  //**********************************************************************

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

	    //	    inside1 = ((n1 > 0.)&&(xx1>0.)) ;
	    inside1 = (n1 > 0.) ;
	    //	    inside2 = ((n2 > 0.)&&(xx2>0.)) ;
	    inside2 = (n2 > 0.) ;

	    // slowrot special treatment follows here.
	    if (slow_rot_style)
	      {
		inside = true;  // no 1-fluid transition!
		n1 = (n1 > 0) ? n1: 0;  // make sure only positive densities
		n2 = (n2 > 0) ? n2: 0;
	      }
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
	    if (inside2) {
	      n2 = nbar_ent_p2(xx2) ;
	      inside2 = (n2 > 0.) ;
	    }
	    if (!inside1) n1 = 0. ;
	    if (!inside2) n2 = 0. ;
	    nbar1.set(l,k,j,i) = n1 ;
	    nbar2.set(l,k,j,i) = n2 ;
	  }

	  ener.set(l,k,j,i) = ener_nbar_p(n1, n2, xx) ;
	  press.set(l,k,j,i) = press_nbar_p(n1, n2, xx) ;

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

// new functions --> ok ?

// Generic computational routine for get_Kxx_ent
//----------------------------------------------
//Rk : nbar1 and nbar2 have been replaced by ent1 and ent2
void Eos_bifluid::calcule_ent(const Cmp& ent1, const Cmp& ent2, const Cmp&
		     x2, int nzet, int l_min, double
		     (Eos_bifluid::*fait)(double, double, double) const, 
			  Cmp& resu) const {
    
    assert(ent1.get_etat() != ETATNONDEF) ; 
    assert(ent2.get_etat() != ETATNONDEF) ; 
    assert(x2.get_etat() != ETATNONDEF) ; 
    
    const Map* mp = ent1.get_mp() ;	// Mapping
    assert(mp == ent2.get_mp()) ;    
    
    
    if ((ent1.get_etat() == ETATZERO)&&(ent2.get_etat() == ETATZERO)) {
	resu.set_etat_zero() ; 
	return ; 
    }
    
    bool le1 = ent1.get_etat() == ETATQCQ ; 
    bool le2 = ent2.get_etat() == ETATQCQ ; 
    bool bx2 = x2.get_etat() == ETATQCQ ; 
    const Valeur* vent1 = 0x0 ;
    const Valeur* vent2 = 0x0 ;
    const Valeur* vx2 = 0x0 ;
    if (le1) { vent1 = &ent1.va ;
    vent1->coef_i() ; }
    if (le2) { vent2 = &ent2.va ;
    vent2->coef_i() ; }
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
      
      Tbl* tent1 = 0x0 ;
      Tbl* tent2 = 0x0 ;
      Tbl* tx2 = 0x0 ;
      
      if (le1) tent1 = vent1->c->t[l] ;
      if (le2) tent2 = vent2->c->t[l] ;
      if (bx2) tx2 = vx2->c->t[l] ;
      Tbl* tresu = vresu.c->t[l] ; 
      
      bool le1b = false ;
      if (le1) le1b = tent1->get_etat() == ETATQCQ ;
      bool le2b = false ;
      if (le2) le2b = tent2->get_etat() == ETATQCQ ;
      bool bx2b = false ;
      if (bx2) bx2b = tx2->get_etat() == ETATQCQ ;
      tresu->set_etat_qcq() ;
      
      double H1, H2, xx ;
      
      for (int i=0; i<tent1->get_taille(); i++) {
	
	H1 = le1b ? tent1->t[i] : 0 ;
	H2 = le2b ? tent2->t[i] : 0 ;
	xx = bx2b ? tx2->t[i] : 0 ;

//cout << " LA" << H1 << "   " << H2 << "   " << xx <<  endl;
	tresu->t[i] = (this->*fait)(xx, H1, H2) ;
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

// rajout de ces fonctions
Cmp Eos_bifluid::get_Knn_ent(const Cmp& ent1, const Cmp& ent2, const Cmp&
			 delta2, int nzet, int l_min) const {
    
    Cmp resu(ent1.get_mp()) ; 

    calcule_ent(ent1, ent2, delta2, nzet, l_min, &Eos_bifluid::get_K11, resu) ;

    return resu ; 
    
}

Cmp Eos_bifluid::get_Knp_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& 
			 delta2, int nzet, int l_min) const {
    
    Cmp resu(delta2.get_mp()) ; 
    
    calcule_ent(ent1, ent2, delta2, nzet, l_min, &Eos_bifluid::get_K12, resu) ;
    
    return resu ; 
    
}

Cmp Eos_bifluid::get_Kpp_ent(const Cmp& ent1, const Cmp& ent2, const Cmp& 
			 delta2, int nzet, int l_min) const {
    
    Cmp resu(ent2.get_mp()) ; 
    
    calcule_ent(ent1, ent2, delta2, nzet, l_min, &Eos_bifluid::get_K22, resu) ;
    
    return resu ; 
    
}

}

