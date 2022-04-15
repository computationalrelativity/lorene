/*
 *  Methods of class Eos_compose_fit
 *
 *  (see file eos_compose_fit.h for documentation).
 *
 */

/*
 *   Copyright (c) 2022 Jerome Novak
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
 * Revision 1.1  2022/04/15 13:39:24  j_novak
 * New class Eos_compose_fit to generate fitted EoSs from CompOSE tables.
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "eos_compose_fit.h"
#include "scalar.h"
#include "utilitaires.h"
#include "unites.h"


namespace Lorene {

    
			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
Eos_compose_fit::Eos_compose_fit(const string& par_file_name) :
  Eos("EoS fitted from CompOSE data"), mg(nullptr),
  mp(nullptr), log_p(nullptr), log_e(nullptr), log_nb(nullptr), log_cs2(nullptr)
{
  ifstream para_file(par_file_name) ;
  read_and_compute(para_file) ;

}


// Constructor from binary file
// ----------------------------
  Eos_compose_fit::Eos_compose_fit(FILE* fich) : Eos(fich), mg(nullptr),
  mp(nullptr), log_p(nullptr), log_e(nullptr), log_nb(nullptr), log_cs2(nullptr)
  {

  char tmp_string[160] ;
  fread(tmp_string, sizeof(char), 160, fich) ;
  tablename = tmp_string ;
  
  //### to be done...

}



// Constructor from a formatted file
// ---------------------------------

  Eos_compose_fit::Eos_compose_fit(ifstream& para_file) :
    Eos(para_file), mg(nullptr), mp(nullptr), log_p(nullptr), log_e(nullptr),
    log_nb(nullptr), log_cs2(nullptr) 
{
  read_and_compute(para_file) ;
}


			//--------------//
			//  Destructor  //
			//--------------//

Eos_compose_fit::~Eos_compose_fit(){
  if (mg != nullptr) delete mg ;
  if (mp != nullptr) delete mp ;
  if (log_p != nullptr) delete log_p ;
  if (log_e != nullptr) delete log_e ;
  if (log_nb != nullptr) delete log_nb ;
  if (log_cs2 != nullptr) delete log_cs2 ;
}

			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Eos_compose_fit::operator==(const Eos& eos_i) const {
    
    bool resu = true ;
    
    if ( eos_i.identify() != identify() ) {
      cout << "The second EOS is not of type Eos_compose_fit !" << endl ;
      resu = false ;
    }
    return resu ; 
  }

  bool Eos_compose_fit::operator!=(const Eos& eos_i) const {
    
    return !(operator==(eos_i)) ;
  }
  
			//------------//
			//  Outputs   //
			//------------//

void Eos_compose_fit::sauve(FILE* fich) const {

  Eos::sauve(fich) ;
  
  char tmp_string[160] ;
  strcpy(tmp_string, tablename.c_str()) ;
  fwrite(tmp_string, sizeof(char), 160, fich) ;

  mg->sauve(fich) ;

  //### to be done...
}

  ostream& Eos_compose_fit::operator>>(ostream & ost) const {

    ost << "EOS of class Eos_compose_fit." << endl ;
    ost << "Built from files in directory " << tablename << endl ;
    ost << "Number of coefficients used for the polynomial fit : "
	<< n_coefs << endl ;
    ost << "nb_min = " << nb_min << ", nb_max = " << nb_max << endl ;
    ost << *mg << endl ;
    ost << *mp << endl ;
    return ost ;

}

		//----------------------------------------------//
		//  Reading of the table and computing the fit  //
		//----------------------------------------------//
			
void Eos_compose_fit::read_and_compute(ifstream& para_file) {

  para_file >> tablename ;
  para_file >> nb_max >> nb_min ;
  para_file.ignore(1000, '\n') ;
  assert(nb_min < nb_max) ;
  para_file >> n_coefs ;
  para_file.ignore(1000, '\n') ;
  int nr ; para_file >> nr ;
  para_file.ignore(1000, '\n') ;
  
  Tbl* logh_read = nullptr ;
  Tbl* logp_read = nullptr ;
  Tbl* loge_read = nullptr ;
  Tbl* lognb_read = nullptr ;
  Tbl* dlpsdlnb = nullptr ;
  int nbp = 0 ;
  
  read_compose_data(nbp, logh_read, logp_read, loge_read, lognb_read, dlpsdlnb) ;
  
  int i_min = 0 ; int i_max = 0 ;

  Tbl nb_read = exp((*lognb_read)) ;

  for (int i=0; i<nbp; i++) {
    if (nb_read(i) < nb_min) i_min = i ;
    if (nb_read(i) < nb_max) i_max = i ;
  }

  int nz = 3 ;
  double x_0 = (*logh_read)(0) ;
  double x_lim_min = (*logh_read)(i_min) ;
  double x_lim_max = (*logh_read)(i_max) ;
  double x_end = (*logh_read)(nbp-1) ;

  // Boundaries of spectral grid
  //------------------------------
  Tbl xtab(nz+1) ; xtab.set_etat_qcq() ;
  xtab.set(0) = x_0 ;
  xtab.set(1) = x_lim_min ;
  xtab.set(2) = x_lim_max ;
  xtab.set(nz) = x_end*(1. + 1.e-12) ;

  // Grid and xx = log(h) coordinate
  mg = new Mg3d(nz, nr, 1, 1, SYM, SYM) ;
  mp = new Map_af(*mg, xtab) ;

  adiabatic_index_fit(i_min, i_max, *logh_read, *logp_read, *loge_read, *lognb_read,
    		      *dlpsdlnb) ;
  
}

			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_compose_fit::nbar_ent_p(double ent, const Param* ) const {

  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_compose_fit::nbar_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logh0 = log( ent ) ;
    return exp(log_nb->val_point(logh0, 0. ,0.)) ;
  }
  else{
    return 0 ;
  }
}

// Energy density from enthalpy
//------------------------------

double Eos_compose_fit::ener_ent_p(double ent, const Param* ) const {
  
  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_compose_fit::ener_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logh0 = log( ent ) ;
    return exp(log_e->val_point(logh0, 0., 0.)) ;
  }
  else{
    return 0 ;
  }
}

// Pressure from enthalpy
//------------------------

double Eos_compose_fit::press_ent_p(double ent, const Param* ) const {
  
    if ( ent > hmin ) {
      if (ent > hmax) {
	cout << "Eos_compose_fit::press_ent_p : ent > hmax !" << endl ;
	abort() ;
      }
      
      double logh0 = log( ent ) ;
      return exp(log_p->val_point(logh0, 0., 0.)) ;
    }
    else{
      return 0 ;
    }
}

// dln(n)/ln(H) from enthalpy 
//---------------------------

double Eos_compose_fit::der_nbar_ent_p(double ent, const Param* ) const {

  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_compose_fit::der_nbar_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    
    double logh0 = log(ent) ;
    return (log_nb->dsdr()).val_point(logh0, 0., 0.) ; 
  }
  else 
    return (log_nb->dsdr()).val_grid_point(0, 0, 0, 0) ;  // to ensure continuity
}


// dln(e)/ln(H) from enthalpy 
//---------------------------

double Eos_compose_fit::der_ener_ent_p(double ent, const Param* ) const {

  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_compose_fit::der_ener_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logh0 = log(ent) ;
    return (log_e->dsdr()).val_point(logh0, 0., 0.) ;
    }
  else
    return (log_e->dsdr()).val_grid_point(0, 0, 0, 0) ; 
}


// dln(p)/ln(H) from enthalpy 
//---------------------------

double Eos_compose_fit::der_press_ent_p(double ent, const Param* ) const {

  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_compose_fit::der_press_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logh0 = log(ent) ;
    return  (log_p->dsdr()).val_point(logh0, 0., 0.) ;
  }
  else  
    return (log_p->dsdr()).val_grid_point(0, 0, 0, 0) ;
}


// dln(p)/dln(n) from enthalpy 
//---------------------------

double Eos_compose_fit::der_press_nbar_p(double ent, const Param*) const {

  double dlpsdlh0 = der_press_ent_p(ent) ;
  double dlnsdlh0 = der_nbar_ent_p(ent) ;
  
  return dlpsdlh0 / dlnsdlh0  ;

}

double Eos_compose_fit::csound_square_ent_p(double ent, const Param*) const {
  
  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_compose_fit::csound_square_ent_p : ent>hmax !" << endl ;
      abort() ;
    }
    double logh0 = log(ent) ;    
    return exp(log_cs2->val_point(logh0, 0., 0.)) ;
  }
  else
      return exp(log_cs2->val_grid_point(0, 0, 0, 0)) ;
}
}
