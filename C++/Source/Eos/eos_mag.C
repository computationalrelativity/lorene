/*
 *  Methods of class Eos_mag
 *
 *  (see file eos_mag.h for documentation).
 *
 */

/*
 *   Copyright (c) 2011 Thomas Elghozi & Jerome Novak
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


char eos_mag_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2011/06/16 10:49:18  j_novak
 * New class Eos_mag for EOSs depending on density and magnetic field.
 *
 *
 * $Header$
 *
 */

// headers C
#include <cmath>

// Headers Lorene
#include "eos.h"
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"


void interpol_herm(const Tbl& , const Tbl&, const Tbl&, double, int&,
		   double&, double& ) ;

void interpol_linear(const Tbl&, const Tbl&, double, int&, double&) ;

			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
Eos_mag::Eos_mag(const char* name_i, const char* table,
		 const char* path) : Eos(name_i), tablename(path) {	

  tablename += '/' ;
  tablename += table ;
	
  read_table() ; 	

}

// Standard constructor with full filename
// ---------------------------------------
Eos_mag::Eos_mag(const char* name_i, const char* file_name) 
  : Eos(name_i), tablename(file_name) {	

	read_table() ; 	

}


// Constructor from binary file
// ----------------------------
Eos_mag::Eos_mag(FILE* fich) : Eos(fich) {

  char tmp_name[160] ;

  fread(tmp_name, sizeof(char), 160, fich) ;		
  tablename += tmp_name ;

  read_table() ;

}



// Constructor from a formatted file
// ---------------------------------
Eos_mag::Eos_mag(ifstream& fich, const char* table) : Eos(fich) {

  fich >> tablename ;
  tablename += '/' ;
  tablename += table ;
  
  read_table() ; 	

}

Eos_mag::Eos_mag(ifstream& fich) : Eos(fich) {

  fich >> tablename ;

  read_table() ; 	

}


			//--------------//
			//  Destructor  //
			//--------------//

Eos_mag::~Eos_mag(){
  delete lognb ;
  delete logp ;
  delete loge ;
  delete logmu ;
  delete magB ;
  delete magM ;
  delete chi ;
}

			//------------//
			//  Outputs   //
			//------------//

void Eos_mag::sauve(FILE* fich) const {
  
  Eos::sauve(fich) ;
  
  fwrite(tablename.c_str(), sizeof(char), 160, fich) ;		

}
			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_mag::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_mag !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_mag::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_mag::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_mag : tabulated EOS depending on two variables: enthalpy and magnetic field."
    	<< '\n' ;
    ost << "Table read from " << tablename << endl ;
    	
    return ost ;

}

			

			//------------------------//
			//  Reading of the table  //
			//------------------------//
			
void Eos_mag::read_table() {

  using namespace Unites_mag ;

  ifstream fich(tablename.data()) ;

  if (!fich) {
    cout << "Eos_mag::read_table(): " << endl ;
    cout << "Problem in opening the EOS file!" << endl ;
    cout << "Aborting..." << endl ;
    abort() ;
  }

  for (int i=0; i<5; i++) {		//  jump over the file
    fich.ignore(1000, '\n') ;             // header
  }                                       //

  int nbp1, nbp2 ;
  fich >> nbp1 >> nbp2 ; fich.ignore(1000, '\n') ;   // number of data
  if ((nbp1<=0) || (nbp2<=0)) {
    cout << "Eos_mag::read_table(): " << endl ;
    cout << "Wrong value for the number of lines!" << endl ;
    cout << "nbp1 = " << nbp1 << ", nbp2 = " << nbp2 << endl ;
    cout << "Aborting..." << endl ;
    abort() ;
  }

  for (int i=0; i<3; i++) {		//  jump over the table
    fich.ignore(1000, '\n') ;
  }                                      
  
  lognb = new Tbl(nbp2, nbp1) ;
  logp = new Tbl(nbp2, nbp1) ;
  loge = new Tbl(nbp2, nbp1) ;
  logmu = new Tbl(nbp2, nbp1) ;
  magB = new Tbl(nbp2, nbp1) ;
  magM = new Tbl(nbp2, nbp1) ;
  chi = new Tbl(nbp2, nbp1) ;
    	
  lognb->set_etat_qcq() ;
  logp->set_etat_qcq() ;
  loge->set_etat_qcq() ;
  logmu->set_etat_qcq() ;
  magB->set_etat_qcq() ;
  magM->set_etat_qcq() ;
  chi->set_etat_qcq() ;
    	
  double rhonuc_cgs = rhonuc_si * 1e-3 ;
  double c2_cgs = c_si * c_si * 1e4 ;
  double mag_PG = mag_unit / 100. ;
    	
  int no1, no2 ;
  double nb_fm3, rho_cgs, p_cgs, mu_MeV, magB_PG, magM_PG, chi_PGpMeV ;
    	
  double ww = 0. ;
    	
  for (int j=0; j<nbp2; j++) {
    
    for (int i=0; i<nbp1; i++) {
      fich >> no1 >> no2 >> nb_fm3 >> rho_cgs >> p_cgs >> mu_MeV 
	   >> magB_PG >> magM_PG >> chi_PGpMeV ; 
      fich.ignore(1000,'\n') ;
      if ( (nb_fm3<0) || (rho_cgs<0)) { // || (p_cgs < 0) ){
	cout << "Eos_mag::read_table(): " << endl ;
	cout << "Negative value in table!" << endl ;
	cout << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	  ", p = " << p_cgs << endl ;
	cout << "Aborting..." << endl ;
	abort() ;
      }
      double psc2_cgs = p_cgs / c2_cgs ;
      // double h_comp = log( (rho_cgs + psc2_cgs) /
      // 			   (10 * nb_fm3 * rhonuc_cgs) ) ;
      double h_read = log(mu_MeV) ;
      
      if ( (i==0) && (j==0) ) ww = h_read ;
    		
      lognb->set(j, i) = log10(nb_fm3) ;
      logp->set(j, i) = log10( psc2_cgs / rhonuc_cgs ) ; 
      loge->set(j, i) = log10( rho_cgs / rhonuc_cgs ) ;	
      logmu->set(j, i) = log10( h_read - ww + 1.e-14 ) ;
      magB->set(j, i) = magB_PG / mag_PG ;
      magM->set(j, i) = magM_PG / mag_PG ;
      chi->set(j, i) = chi_PGpMeV / (mag_PG*exp(ww)) ;
    }
  }
            
  hmin = pow( double(10), (*logmu)(0, 0) ) ;
  hmax = pow( double(10), (*logmu)(0, nbp1-1) ) ;
  Bmax = (*magB)(nbp2-1, 0) ;

  fich.close();
 
}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_mag::nbar_ent_p(double ent, const Param* ) const {

  c_est_pas_fait("Eos_mag::nbar_ent_p" ) ;

  return ent ;

}

// Energy density from enthalpy
//------------------------------

double Eos_mag::ener_ent_p(double ent, const Param* ) const {

  c_est_pas_fait("Eos_mag::ener_ent_p" ) ;

  return ent ;

}

// Pressure from enthalpy
//------------------------

double Eos_mag::press_ent_p(double ent, const Param* ) const {

  c_est_pas_fait("Eos_mag::press_ent_p" ) ;

  return ent ;

}

// dln(n)/ln(H) from enthalpy 
//---------------------------

double Eos_mag::der_nbar_ent_p(double ent, const Param* ) const {

  c_est_pas_fait("Eos_mag::der_nbar_ent_p" ) ;

  return ent ;

}


// dln(e)/ln(H) from enthalpy 
//---------------------------

double Eos_mag::der_ener_ent_p(double ent, const Param* ) const {


  c_est_pas_fait("Eos_mag::der_ener_enr_p" ) ;

  return ent ;

}


// dln(p)/ln(H) from enthalpy 
//---------------------------

double Eos_mag::der_press_ent_p(double ent, const Param* ) const {

  c_est_pas_fait("Eos_mag::" ) ;

  return ent ;

}


// dln(p)/dln(n) from enthalpy 
//---------------------------

double Eos_mag::der_press_nbar_p(double ent, const Param*) const {

  c_est_pas_fait("Eos_mag::der_press_nbar_p" ) ;

  return ent ;
          
}
