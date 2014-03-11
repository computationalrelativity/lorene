 /*  Methods of class Eos_mag
 *
 *  (see file eos_mag.h for documentation).
 *
 */

/*
 *   Copyright (c) 2011 Thomas Elghozi & Jerome Novak
 *   Copyright (c) 2013 Debarati Chatterjee
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
 * Revision 1.8  2014/03/11 14:27:26  j_novak
 * Corrected a missing 4pi term.
 *
 * Revision 1.7  2014/03/03 16:23:08  j_novak
 * Updated error message
 *
 * Revision 1.6  2013/12/12 16:07:30  j_novak
 * interpol_herm_2d outputs df/dx, used to get the magnetization.
 *
 * Revision 1.5  2013/11/25 15:00:52  j_novak
 * Correction of memory error.
 *
 * Revision 1.4  2013/11/14 16:12:55  j_novak
 * Corrected a mistake in the units.
 *
 * Revision 1.2  2011/10/04 16:05:19  j_novak
 * Update of Eos_mag class. Suppression of loge, re-definition of the derivatives
 * and use of interpol_herm_2d.
 *
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
#include "cmp.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"


void interpol_herm_2d(const Tbl&, const Tbl&, const Tbl&, const Tbl&, const Tbl&, const Tbl&, double, double, double&, double&, double&) ;


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
  delete d2lp ;
  delete dlpsdB ;
  delete dlpsdlh ;
  delete Bfield ;
  delete logh ;
  delete logp ;
}

			//------------//
			//  Outputs   //
			//------------//

void Eos_mag::sauve(FILE* fich) const {
  
  Eos::sauve(fich) ;
  
  fwrite(tablename.c_str(), sizeof(char), tablename.size(), fich) ;		

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
  
  logp = new Tbl(nbp2, nbp1) ;
  logh = new Tbl(nbp2, nbp1) ;
  Bfield = new Tbl(nbp2, nbp1) ;
  dlpsdlh = new Tbl(nbp2, nbp1) ;
  dlpsdB = new Tbl(nbp2, nbp1) ;
  d2lp = new Tbl(nbp2, nbp1) ;
    	
  logp->set_etat_qcq() ;
  logh->set_etat_qcq() ;
  Bfield->set_etat_qcq() ;
  dlpsdlh->set_etat_qcq() ;
  dlpsdB->set_etat_qcq() ;
  d2lp->set_etat_qcq() ;
  
  double c2 = c_si * c_si ;
  double mag_PG = mag_unit / 100. ; 
    	
  int no1 ;
  double nb_fm3, rho_cgs, p_cgs, mu_MeV, magB_PG, magM_PG, chi_PGpMeV ;
    	
  double ww = 0. ;

  for (int j=0; j<nbp2; j++) {
    
    for (int i=0; i<nbp1; i++) {
      fich >> no1 >> nb_fm3 >> rho_cgs >> p_cgs >> mu_MeV 
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
      
      double p_si = p_cgs*0.1 ; // in N m^-2
      double psc2 = p_si/c2 ; //  in kg m^-3
      double rho_si = rho_cgs*1000. ; // in kg m^-3

      double h_read = log(mu_MeV) ;
      if ( (i==0) && (j==0) ) ww = h_read ;
      double h_new = h_read - ww ;

      logp->set(j, i) = psc2/rhonuc_si ; 
      logh->set(j, i) = h_new ;
      Bfield->set(j, i) = magB_PG / mag_PG ; // in Lorene units
      dlpsdlh->set(j, i) = (rho_si + psc2)/rhonuc_si ; 
      dlpsdB->set(j, i) = magM_PG*mag_PG/c2/rhonuc_si/(4*M_PI) ; 
      d2lp->set(j, i) = mu_MeV*chi_PGpMeV*mag_PG/c2/rhonuc_si ; 

    }
  }

  hmin = (*logh)(0, 0) ;
  hmax = (*logh)(0, nbp1-1) ;

  Bmax = (*Bfield)(nbp2-1, 0) ;


  fich.close();
 
}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_mag::nbar_ent_p(double ent, const Param* par ) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if ( ent >= hmin) {
    if (ent > hmax) {
      cout << "Eos_tabul::nbar_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    // recuperer magB0 (input)
    double magB0 = 0. ;
    if (par != 0x0) 
      if (par->get_n_double_mod() > 0) {
	magB0 = par->get_double_mod() ;
      }
    
    double p_int, dpdb_int, dpdh_int ;
    
    interpol_herm_2d(*Bfield, *logh, *logp, *dlpsdB, *dlpsdlh, *d2lp, magB0, ent, 
		     p_int, dpdb_int, dpdh_int) ;
    
    double nbar_int = dpdh_int * exp(-ent) ;
    
    return nbar_int ;
  }
  else{
    return 0 ;
  }
}


// Energy density from enthalpy
//------------------------------

double Eos_mag::ener_ent_p(double ent, const Param* par ) const {

  
  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if ( ent >= hmin ) {
    if (ent > hmax) {
      cout << "Eos_tabul::ener_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    
    // recuperer magB0 (input)
    double magB0 = 0. ;
    if (par != 0x0) 
      if (par->get_n_double_mod() > 0) {
	magB0 = par->get_double_mod() ;
      }

    double p_int, dpdb_int, dpdh_int ;

    interpol_herm_2d(*Bfield, *logh, *logp, *dlpsdB, *dlpsdlh, *d2lp, magB0, ent, 
		     p_int, dpdb_int, dpdh_int) ;
    
    double nbar_int = dpdh_int * exp(-ent) ;

    double f_int = - p_int + exp(ent) * nbar_int;

    return f_int ;
  }
  else{
    return 0 ;
  }
}

// Pressure from enthalpy
//------------------------

double Eos_mag::press_ent_p(double ent, const Param* par ) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if ( ent >= hmin ) {                            
    if (ent > hmax) {
      cout << "Eos_mag::press_ent_p : ent > hmax !" << endl ;
      abort() ;
    }

    // recuperer magB0 (input)
    double magB0 = 0. ;
    if (par != 0x0) 
      if (par->get_n_double_mod() > 0) {
	magB0 = par->get_double_mod() ;
      }
    
    double p_int, dpdb_int, dpdh_int ;

    interpol_herm_2d(*Bfield, *logh, *logp, *dlpsdB, *dlpsdlh, *d2lp, magB0, ent, 
		     p_int, dpdb_int, dpdh_int) ; 
    
    return p_int;
  }
  else{
    return 0 ;
  }
}

// Magnetization from enthalpy
//----------------------------

double Eos_mag::mag_ent_p(double ent, const Param* par) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if ( ent >= hmin ) {                            
    if (ent > hmax) {
      cout << "Eos_mag::mag_ent_p : ent > hmax !" << endl ;
      abort() ;
    }

    // recuperer magB0 (input)
    double magB0 = 0. ;
    if (par != 0x0) 
      if (par->get_n_double_mod() > 0) {
	magB0 = par->get_double_mod() ;
      }
    
    double p_int, dpdb_int, dpdh_int ;

    interpol_herm_2d(*Bfield, *logh, *logp, *dlpsdB, *dlpsdlh, *d2lp, magB0, ent, 
		     p_int, dpdb_int, dpdh_int) ; 
    
    double magnetization = 4*M_PI * dpdb_int ;

    if (magB0 == 0.) //## Test to be improved...
      return 0 ;
    else 
      return  magnetization / magB0 ;
  }
  
  else 
    return 0. ;

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
