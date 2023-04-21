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
 * Revision 1.6  2023/04/21 14:51:51  j_novak
 * The analytic part for low densities has been changed to a new model with one parameter more than the polytrope, allowing for a smooth transition of all quantities (p, e, nB AND cs2). More checks need to be done...
 *
 * Revision 1.5  2023/03/17 08:36:59  j_novak
 * Output of "middle" enthalpy (boundary between domains).
 *
 * Revision 1.4  2023/01/27 16:10:35  j_novak
 * A polytrope (class Eos_poly) is used for low and high enthalpies.
 *
 * Revision 1.3  2022/07/21 12:34:18  j_novak
 * Corrected units
 *
 * Revision 1.2  2022/07/20 12:59:43  j_novak
 * Added methods for saving into files and constructor from a file
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
    Eos("EoS fitted from CompOSE data"), p_eos_high(nullptr),
    mg(nullptr), mp(nullptr), log_p(nullptr), log_e(nullptr), log_nb(nullptr),
    log_cs2(nullptr)
  {
    ifstream para_file(par_file_name) ;
    para_file.ignore(1000, '\n') ;
    para_file.ignore(1000, '\n') ;
    read_and_compute(para_file) ;
    
  }


// Constructor from binary file
// ----------------------------
  Eos_compose_fit::Eos_compose_fit(FILE* fich) :
    Eos(fich), p_eos_high(nullptr), mg(nullptr),
    mp(nullptr), log_p(nullptr), log_e(nullptr), log_nb(nullptr), log_cs2(nullptr)
  {

    char tmp_string[160] ;
    size_t res = fread(tmp_string, sizeof(char), 160, fich) ;
    if (res != 0)
      tablename = tmp_string ;
    else {
      cerr << "Eos_compose_fit: constructor from a binary file: " << endl ;
      cerr << "Error reading tablename... aborting." << endl ;
      abort() ;
    }
    fread_be(&n_coefs, sizeof(int), 1, fich) ;
    fread_be(&nb_min, sizeof(double), 1, fich) ;
    fread_be(&nb_mid, sizeof(double), 1, fich) ;
    fread_be(&nb_max, sizeof(double), 1, fich) ;
    fread_be(&hmin, sizeof(double), 1, fich) ;
    fread_be(&hmax, sizeof(double), 1, fich) ;
    double gamma_high, kappa_high, m0_high, mu0_high ;
    fread_be(&gamma_high, sizeof(double), 1, fich) ;
    fread_be(&kappa_high, sizeof(double), 1, fich) ;
    fread_be(&m0_high, sizeof(double), 1, fich) ;
    fread_be(&mu0_high, sizeof(double), 1, fich) ;
    fread_be(&k_low, sizeof(double), 1, fich) ;
    fread_be(&c_low, sizeof(double), 1, fich) ;
    fread_be(&alpha_low, sizeof(double), 1, fich) ;

    p_eos_high = new Eos_poly(gamma_high, kappa_high, m0_high, mu0_high) ;

    mg = new Mg3d(fich) ;
    mp = new Map_af(*mg, fich) ;
    
    log_p = new Scalar(*mp, *mg, fich) ;
    log_e = new Scalar(*mp, *mg, fich) ;
    log_nb = new Scalar(*mp, *mg, fich) ;
    log_cs2 = new Scalar(*mp, *mg, fich) ;
  }



// Constructor from a formatted file
// ---------------------------------

  Eos_compose_fit::Eos_compose_fit(ifstream& para_file) :
    Eos(para_file), p_eos_high(nullptr),
    mg(nullptr), mp(nullptr), log_p(nullptr), log_e(nullptr), log_nb(nullptr),
    log_cs2(nullptr) 
  {
    read_and_compute(para_file) ;
  }
  

			//--------------//
			//  Destructor  //
			//--------------//

Eos_compose_fit::~Eos_compose_fit(){
  if (p_eos_high != nullptr) delete p_eos_high ;
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
  fwrite_be(&n_coefs, sizeof(int), 1, fich) ;
  fwrite_be(&nb_min, sizeof(double), 1, fich) ;
  fwrite_be(&nb_mid, sizeof(double), 1, fich) ;
  fwrite_be(&nb_max, sizeof(double), 1, fich) ;
  fwrite_be(&hmin, sizeof(double), 1, fich) ;
  fwrite_be(&hmax, sizeof(double), 1, fich) ;
  double gamma_high = p_eos_high->get_gam() ;
  double kappa_high = p_eos_high->get_kap() ;
  double m0_high = p_eos_high->get_m_0() ;
  double mu0_high = p_eos_high->get_mu_0() ;
  fwrite_be(&gamma_high, sizeof(double), 1, fich) ;
  fwrite_be(&kappa_high, sizeof(double), 1, fich) ;
  fwrite_be(&m0_high, sizeof(double), 1, fich) ;
  fwrite_be(&mu0_high, sizeof(double), 1, fich) ;
  fwrite_be(&k_low, sizeof(double), 1, fich) ;
  fwrite_be(&c_low, sizeof(double), 1, fich) ;
  fwrite_be(&alpha_low, sizeof(double), 1, fich) ;
  
  mg->sauve(fich) ;
  mp->sauve(fich) ;

  log_p->sauve(fich) ;
  log_e->sauve(fich) ;
  log_nb->sauve(fich) ;
  log_cs2->sauve(fich) ;
  
}

  ostream& Eos_compose_fit::operator>>(ostream & ost) const {

    ost << "EOS of class Eos_compose_fit." << endl ;
    ost << "Built from files in directory " << tablename << endl ;
    ost << "Number of coefficients used for the polynomial fit : "
	<< n_coefs << endl ;
    ost << "nb_min = " << nb_min << ", nb_mid = " << nb_mid
    	<< ", nb_max = " << nb_max << " [fm^-3]" << endl ;
    ost << "hmin = " << hmin << ", hmid = " << exp(mp->val_r_jk(0, 1., 0, 0))
	<< ", hmax = " << hmax << endl ;
    ost << "EoS for low density part: " << endl ;
    ost << " of the form p = kappa(1 - exp(alpha H))^(1/(alpha C))" << endl ;
    ost << "kappa = " << k_low << ", C = " << c_low << ", alpha = " << alpha_low
	<< endl ;
    ost << "EoS for high density part: " << *p_eos_high ;
    ost << *mg << endl ;
    return ost ;

}

		//----------------------------------------------//
		//  Reading of the table and computing the fit  //
		//----------------------------------------------//
			
void Eos_compose_fit::read_and_compute(ifstream& para_file) {

  cout << para_file.good() << endl ;
  para_file >> tablename ;
  para_file >> nb_mid >> nb_min ;
  para_file.ignore(1000, '\n') ;
  assert(nb_min < nb_mid) ;
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
  
  int i_min = 0 ; int i_mid = 0 ;

  Tbl nb_read = exp((*lognb_read)) ;

  for (int i=0; i<nbp; i++) {
    if (nb_read(i) < 10.*nb_min) i_min = i ; // nb_min & nb_mid are in fm^-3
    if (nb_read(i) < 10.*nb_mid) i_mid = i ; // Lorene units are 0.1 fm^-3
  }

  int i_max = nbp - 1 ; //## to improve, depending on the sound speed values

  int nz = 2 ;
  double x_min = (*logh_read)(i_min) ;
  double x_mid = (*logh_read)(i_mid) ;
  double x_max = (*logh_read)(i_max) ;

  hmin = exp(x_min) ;
  hmax = exp(x_max) ;
  
  // Boundaries of spectral grid
  //------------------------------
  Tbl xtab(nz+1) ; xtab.set_etat_qcq() ;
  xtab.set(0) = x_min ;
  xtab.set(1) = x_mid ;
  xtab.set(nz) = x_max ;

  // Grid and xx = log(h) coordinate
  mg = new Mg3d(nz, nr, 1, 1, SYM, SYM) ;
  mp = new Map_af(*mg, xtab) ;

  adiabatic_index_fit(i_min, i_mid, *logh_read, *logp_read, *loge_read, *lognb_read,
    		      *dlpsdlnb) ;

  // Cleaning
  if (logh_read != nullptr) delete logh_read ;
  if (logp_read != nullptr) delete logp_read ;
  if (loge_read != nullptr) delete loge_read ;
  if (lognb_read != nullptr) delete lognb_read ;
  if (dlpsdlnb != nullptr) delete dlpsdlnb ;  
}

void Eos_compose_fit::write_lorene_table(const string& name_file, int nlines)
  const {

  cout << "Writing the Eos_compose_fit object into a Lorene-format file ("
       << name_file << ") ..." << flush ;
  if (nlines < 10) {
    cerr << "Eos_compose_fit::write_lorene_table" << endl ;
    cerr << "The number of lines to be outputted is too small!" << endl ;
    cerr << " nlines = " << nlines << endl ;
    cerr << "Aborting..." << endl ;
    abort() ;
  }
  double rhonuc_cgs = Unites::rhonuc_si * 1.e-3 ;
  double c2_cgs = Unites::c_si * Unites::c_si * 1.e4 ;

  ofstream lorene_file(name_file) ;
  Scalar press = exp(*log_p) * c2_cgs * rhonuc_cgs ;
  Scalar ener = exp(*log_e) * rhonuc_cgs ;
  Scalar nbar = exp(*log_nb) * 10. ;

  lorene_file << "#" << endl ;
  lorene_file << "# Built from Eos_compose_fit object" << endl ;
  lorene_file << "# From the table: " << tablename << endl ;
  lorene_file << "#\n#" << endl ;
  lorene_file << nlines << "\t Number of lines" << endl ;
  lorene_file << "#\n#\t  n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]" << endl ;
  lorene_file << "#" << endl ;

  //  const Coord& xx = mp->r ;
#ifndef NDEBUG
  int nz = mg->get_nzone() ;
  assert (nz > 1) ;
#endif
  double xmin0 = log(1.e-14) ;
  double xmax0 = log(hmin) ;
  int nlines0 = nlines / 10 ;
  double dx = (xmax0 - xmin0)/double(nlines0-1) ;
  assert(dx>0.) ;
  double logh = xmin0 ;
  for (int i=0; i<nlines0; i++) {
    double ent = exp(logh) ;
    lorene_file << setprecision(16) ;
    lorene_file << i << '\t' << nbar_ent_p(ent)/10. << '\t'
		<< ener_ent_p(ent)*rhonuc_cgs / c2_cgs << '\t'
		<< press_ent_p(ent) *c2_cgs * rhonuc_cgs << endl ;
    logh += dx ;
  }
  logh -= dx ;
  xmin0 = xmax0 ;
  xmax0 = log(10.*hmax) ; // ## to improve?
  dx = (xmax0 - xmin0) / double(nlines - nlines0-1) ;
  for (int i=nlines0; i<nlines; i++) {
    double ent = exp(logh) ;
    lorene_file << setprecision(16) ;
    lorene_file << i << '\t' << nbar_ent_p(ent)/10. << '\t'
		<< ener_ent_p(ent)*rhonuc_cgs / c2_cgs << '\t'
		<< press_ent_p(ent) *c2_cgs * rhonuc_cgs << endl ;
    logh += dx ;
  }
  lorene_file.close() ;
  cout << " done!" << endl ;
}




			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_compose_fit::nbar_ent_p(double ent, const Param* ) const {

  if ( ent > hmin ) {
    if (ent < hmax) {
      double logh0 = log( ent ) ;
      return exp(log_nb->val_point(logh0, 0. ,0.)) ;
    }
    else {
      assert(p_eos_high != nullptr) ;
      return p_eos_high->nbar_ent_p(ent) ;
    }
  }
  else{
    if (ent > 0.) {
      double expam1 = expm1(alpha_low*ent) ;
      double expmam1 = expm1(-alpha_low*ent) ;
      double ppp = k_low*pow(fabs(expam1), 1./(alpha_low*c_low)) ;
      double yyy = -c_low*expmam1 ;
      double nba = ppp / (yyy*exp(ent)) ;
      return nba ;
    }
    else
      return 0. ;
  }
}

// Energy density from enthalpy
//------------------------------

double Eos_compose_fit::ener_ent_p(double ent, const Param* ) const {
  
  if ( ent > hmin ) {
    if (ent < hmax) {
      double logh0 = log( ent ) ;
      return exp(log_e->val_point(logh0, 0., 0.)) ;
    }
    else {
      assert(p_eos_high != nullptr) ;
      return p_eos_high->ener_ent_p(ent) ;
    }
  }
  else{
    if (ent > 0.) {
      double expam1 = expm1(alpha_low*ent) ;
      double expmam1 = expm1(-alpha_low*ent) ;
      double ppp = k_low*pow(fabs(expam1), 1./(alpha_low*c_low)) ;
      double yyy = -c_low*expmam1 ;
      return ppp*(1-yyy) / yyy ;
    }
    else
      return 0. ;
  }
}

// Pressure from enthalpy
//------------------------

double Eos_compose_fit::press_ent_p(double ent, const Param* ) const {
  
  if ( ent > hmin ) {
    if (ent < hmax) {
      double logh0 = log( ent ) ;
      return exp(log_p->val_point(logh0, 0., 0.)) ;
    }
    else {
      assert(p_eos_high != nullptr) ;
      return p_eos_high->press_ent_p(ent) ;
    }
  }
  else{
    if (ent > 0.) {
      double expam1 = fabs(expm1(alpha_low*ent)) ;
      return k_low*pow(expam1, 1./(alpha_low*c_low)) ;
    }
    else
      return 0. ;
  }
}

// dln(n)/ln(H) from enthalpy 
//---------------------------

double Eos_compose_fit::der_nbar_ent_p(double ent, const Param* ) const {

  if (ent > 0.) 
    return ent / csound_square_ent_p(ent) ;
  else
    return 1./(c_low*alpha_low) - 1. ;
}


// dln(e)/ln(H) from enthalpy 
//---------------------------

double Eos_compose_fit::der_ener_ent_p(double ent, const Param* ) const {

  if (ent > 0.) {
    double ener = ener_ent_p(ent) ;
    double press = press_ent_p(ent) ;
    double cs2 = csound_square_ent_p(ent) ;
    
    return (ener + press)*ent / (ener * cs2) ;
  }
  else
    return  1./(c_low*alpha_low) - 1. ;
}


// dln(p)/ln(H) from enthalpy 
//---------------------------

double Eos_compose_fit::der_press_ent_p(double ent, const Param* ) const {
  
  if (ent > 0.) {
    double ener = ener_ent_p(ent) ;
    double press = press_ent_p(ent) ;
    
    return ent*(ener + press)/press ;
  }
  else
    return 1./(c_low*alpha_low) ;
  
}


// dln(p)/dln(n) from enthalpy 
//---------------------------

double Eos_compose_fit::der_press_nbar_p(double ent, const Param*) const {

  if (ent > 0.) {
    double dlpsdlh0 = der_press_ent_p(ent) ;
    double dlnsdlh0 = der_nbar_ent_p(ent) ;
    
    return dlpsdlh0 / dlnsdlh0  ;
  }
  else
    return 1./ (1. - alpha_low*c_low) ;

}

double Eos_compose_fit::csound_square_ent_p(double ent, const Param*) const {
  
  if ( ent > hmin ) {
    if (ent < hmax) {
      double logh0 = log( ent ) ;    
      return exp(log_cs2->val_point(logh0, 0., 0.)) ;
    }
    else {
      assert(p_eos_high != nullptr) ;
      return p_eos_high->csound_square_ent_p(ent) ;
    }
  }
  else{
    if (ent > 0.) {
      double expmam1 = expm1(-alpha_low*ent) ;
      double yyy = -c_low*expmam1 ;
      return yyy / (1 - c_low*(1 + alpha_low-1)*exp(-alpha_low*ent)) ;
    }
    else
      return 0. ;
  }
}

} // end of namespace Lorene
