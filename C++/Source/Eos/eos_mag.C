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


void interpol_herm_2d(const Tbl&, const Tbl&, const Tbl&, const Tbl&, const Tbl&, 
		      const Tbl&, double, double, double&, double&) ;


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
  
  double rhonuc_cgs = rhonuc_si * 1e-3 ;
  double c2_cgs = c_si * c_si * 1e4 ;
  double MeVpfm3_cgs = 1.6022e33 ;
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
      double h_read = log(mu_MeV) ;
      
      if ( (i==0) && (j==0) ) ww = h_read ;
    		
      double h_new = h_read - ww + 1.e-14 ;
      logp->set(j, i) = log10( psc2_cgs / rhonuc_cgs ) ; 
      logh->set(j, i) = log10( h_new ) ;
      Bfield->set(j, i) = magB_PG / mag_PG ;
      dlpsdlh->set(j, i) = h_new * (rho_cgs + psc2_cgs) / psc2_cgs ; 
      dlpsdB->set(j, i) = (magM_PG * rhonuc_cgs) / (mag_PG * psc2_cgs) ;
      d2lp->set(j, i) = h_new * mu_MeV * 
	( (chi_PGpMeV / mag_PG) * (rhonuc_cgs / p_cgs)
	  - nb_fm3 * (MeVpfm3_cgs / p_cgs) * (*dlpsdB)(j, i) ) ;
    }
  }
            
  hmin = pow( double(10), (*logh)(0, 0) ) ;
  hmax = pow( double(10), (*logh)(0, nbp1-1) ) ;

  Bmax = (*Bfield)(nbp2-1, 0) ;

  fich.close();
 
}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_mag::nbar_ent_p(double ent, const Param* par ) const {

  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_tabul::nbar_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logent0 = log10( ent ) ;
    // recuperer magB0 (input)
    const Cmp& par_mag = par->get_cmp();
    int lg = par->get_int();
    int kg = par->get_int();
    int jg = par->get_int();
    int ig = par->get_int();
    double magB0 = par_mag(lg,kg,jg,ig);
    
    double p_int, dp_int ;
    interpol_herm_2d(*Bfield, *logh, *logp, *dlpsdB, *dlpsdlh, *d2lp, magB0, logent0, 
		     p_int, dp_int) ;

    double nbar_int = pow(double(10), p_int) * dp_int * exp(-ent) / ent ;

    return nbar_int ;
    }
    else{
	return 0 ;
    }
}


// Energy density from enthalpy
//------------------------------

double Eos_mag::ener_ent_p(double ent, const Param* par ) const {

  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_tabul::ener_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    
    double logent0 = log10( ent ) ;
    // recuperer magB0 (input)
    const Cmp& par_mag = par->get_cmp();
    int lg = par->get_int();
    int kg = par->get_int();
    int jg = par->get_int();
    int ig = par->get_int();
    double magB0 = par_mag(lg,kg,jg,ig);
    
    double logp_int, dlogp_int ;
    interpol_herm_2d(*Bfield, *logh, *logp, *dlpsdB, *dlpsdlh, *d2lp, magB0, logent0, 
		     logp_int, dlogp_int) ;
    
    double p_int = pow(double(10), logp_int) ;
    double nbar_int = p_int * dlogp_int * exp(-ent) / ent ;

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

  if ( ent > hmin ) {                              /////////   indices i et j ??
    if (ent > hmax) {
      cout << "Eos_mag::press_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logent0 = log10( ent ) ;
    // recuperer magB0 (input)
    const Cmp& par_mag = par->get_cmp();
    int lg = par->get_int();
    int kg = par->get_int();
    int jg = par->get_int();
    int ig = par->get_int();
    double magB0 = par_mag(lg,kg,jg,ig);
    
    double p_int, dp_int ;
    interpol_herm_2d(*Bfield, *logh, *logp, *dlpsdB, *dlpsdlh, *d2lp, magB0, logent0, 
		     p_int, dp_int) ;

    return pow(10., p_int);
  }
  else{
    return 0 ;
  }
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
