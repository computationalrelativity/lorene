/*
 *  Methods of class Eos_consistent
 *
 *  (see file eos_compose.h for documentation).
 *
 */

/*
 *   Copyright (c) 2015 Jerome Novak
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
 * Revision 1.7  2022/03/10 16:38:39  j_novak
 * log(cs^2) is tabulated instead of cs^2.
 *
 * Revision 1.6  2021/05/31 11:31:23  g_servignat
 * Added csound_square_ent routine to calculate the sound speed from enthalpy to Eos_consistent and corrected error outputs
 *
 * Revision 1.5  2019/03/28 13:41:02  j_novak
 * Improved managed of saved EoS (functions sauve and constructor form FILE*)
 *
 * Revision 1.4  2017/08/11 13:42:04  j_novak
 * Suppression of spurious output
 *
 * Revision 1.3  2017/08/10 15:14:27  j_novak
 * Now Eos_consistent is also reading Lorene format tables.
 *
 * Revision 1.2  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.1  2015/08/04 14:41:29  j_novak
 * Back to previous version for Eos_CompOSE. Enthalpy-consistent EoS can be accessed using Eos_consistent class (derived from Eos_CompOSE).
 *
 *
 * $Header$
 *
 */

#include <string>

// Headers Lorene
#include "headcpp.h"
#include "eos.h"
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"
#include<string>

namespace Lorene {
  void interpol_herm(const Tbl& , const Tbl&, const Tbl&, double, int&,
		   double&, double& ) ;
  void interpol_linear(const Tbl&, const Tbl&, double, int&, double&) ;

			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
  Eos_consistent::Eos_consistent(const char* file_name)
    : Eos_CompOSE(file_name) { }
  

// Constructor from binary file
// ----------------------------
  Eos_consistent::Eos_consistent(FILE* fich) : Eos_CompOSE(fich) {
    if (format == 0) read_table() ;
    else read_compose_data() ;
  }
  

// Constructor from a formatted file
// ---------------------------------
  Eos_consistent::Eos_consistent(ifstream& fich) : Eos_CompOSE(fich) {
    read_table() ; 
  }
  

// Constructor from CompOSE data files
// ------------------------------------
  Eos_consistent::Eos_consistent(const string& files) : Eos_CompOSE(files) 
  {
    read_compose_data() ;
  }
  
  void Eos_consistent::read_table() {
    
    using namespace Unites ;

    cout << "Computing a thermodynamically - consistent table." << endl ;
    
    ifstream tfich(tablename.data()) ;
    
    for (int i=0; i<5; i++) {		//  jump over the file
      tfich.ignore(1000, '\n') ;             // header
    }                                       //
    
    int nbp ;
    tfich >> nbp ; tfich.ignore(1000, '\n') ;   // number of data
    if (nbp<=0) {
      cerr << "Eos_consistent::read_table(): " << endl ;
      cerr << "Wrong value for the number of lines!" << endl ;
      cerr << "nbp = " << nbp << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }
    
    for (int i=0; i<3; i++) {		//  jump over the table
      tfich.ignore(1000, '\n') ;             // description
    }
    
    press = new double[nbp] ;
    nb    = new double[nbp] ;
    ro    = new double[nbp] ; 
    
    double rhonuc_cgs = rhonuc_si * 1e-3 ;
    double c2_cgs = c_si * c_si * 1e4 ;    	
    
    int no ;
    double nb_fm3, rho_cgs, p_cgs ;
    
    for (int i=0; i<nbp; i++) {
      
      tfich >> no ;
      tfich >> nb_fm3 ;
      tfich >> rho_cgs ;
      tfich >> p_cgs ; tfich.ignore(1000,'\n') ;    		
      if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	cout << "Eos_consistent(ifstream&): " << endl ;
	cout << "Negative value in table!" << endl ;
	cout << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	  ", p = " << p_cgs << endl ;
	cout << "Aborting..." << endl ;
	abort() ;
      }
      
      press[i] = p_cgs / c2_cgs ;
      nb[i]    = nb_fm3 ;
      ro[i]    = rho_cgs ;
    }
    
    Tbl pp(nbp) ; pp.set_etat_qcq() ;
    Tbl dh(nbp) ; dh.set_etat_qcq() ;
    for (int i=0; i<nbp; i++) {
      pp.set(i) = log(press[i] / rhonuc_cgs) ;
      dh.set(i) = press[i] / (ro[i] + press[i]) ;
    }
    
    Tbl hh  = integ1D(pp, dh) + 1.e-14 ;
    
    for (int i=0; i<nbp; i++) {
      logh->set(i) = log10( hh(i) ) ;
      logp->set(i) = log10( press[i] / rhonuc_cgs ) ;
      dlpsdlh->set(i) = hh(i) / dh(i) ;
      lognb->set(i) = log10(nb[i]) ;
    }
    
    hmin = pow( double(10), (*logh)(0) ) ;
    hmax = pow( double(10), (*logh)(nbp-1) ) ;
    
    // Cleaning
    //---------
    delete [] press ; 
    delete [] nb ; 
    delete [] ro ; 
    
  }


  void Eos_consistent::read_compose_data() {

    using namespace Unites ;

    cout << "Computing a thermodynamically - consistent table." << endl ;

    // Files containing data and a test
    //---------------------------------
    string file_nb = tablename + ".nb" ;
    string file_thermo = tablename + ".thermo" ;
    
    ifstream in_nb(file_nb.data()) ;
    
    // obtaining the size of the tables for memory allocation
    //-------------------------------------------------------
    int index1, index2 ;
    in_nb >> index1 >> index2 ;
    int nbp = index2 - index1 + 1 ;
    assert( nbp == logh->get_taille() ) ;
    
    press = new double[nbp] ;
    nb    = new double[nbp] ;
    ro    = new double[nbp] ; 
    
    // Variables and conversion
    //-------------------------
    double nb_fm3, rho_cgs, p_cgs, p_over_nb_comp, eps_comp ;
    double dummy_x ;
    int dummy_n ;
    
    double rhonuc_cgs = rhonuc_si * 1e-3 ;
    double c2_cgs = c_si * c_si * 1e4 ;    
    double m_neutron_MeV, m_proton_MeV ;
    
    ifstream in_p_rho (file_thermo.data()) ;
    in_p_rho >> m_neutron_MeV >> m_proton_MeV ; //Neutron and proton masses
    in_p_rho.ignore(1000, '\n') ;
    
    double p_convert = mev_si * 1.e45 * 10. ; // Conversion from MeV/fm^3 to cgs
    double eps_convert = mev_si * 1.e42 / (c_si*c_si) ; //From meV/fm^3 to g/cm^3
    
    // Main loop reading the table
    //----------------------------
    for (int i=0; i<nbp; i++) {
      in_nb >> nb_fm3 ;
      in_p_rho >> dummy_n >> dummy_n >> dummy_n >> p_over_nb_comp ;
      in_p_rho >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> eps_comp ;
      in_p_rho.ignore(1000, '\n') ;
      p_cgs = p_over_nb_comp * nb_fm3 * p_convert ;
      rho_cgs = ( eps_comp + 1. ) * m_neutron_MeV * nb_fm3 * eps_convert ;
      
      press[i] = p_cgs / c2_cgs ; 
      nb[i]    = nb_fm3 ;
      ro[i]    = rho_cgs ; 
  }
    
    Tbl pp(nbp) ; pp.set_etat_qcq() ;
    Tbl dh(nbp) ; dh.set_etat_qcq() ;
    for (int i=0; i<nbp; i++) {
      pp.set(i) = log(press[i] / rhonuc_cgs) ;
      dh.set(i) = press[i] / (ro[i] + press[i]) ;
    }
    
    Tbl hh  = integ1D(pp, dh) + 1.e-14 ;
    
    for (int i=0; i<nbp; i++) {
      logh->set(i) = log10( hh(i) ) ;
      logp->set(i) = log10( press[i] / rhonuc_cgs ) ;
      dlpsdlh->set(i) = hh(i) / dh(i) ;
      lognb->set(i) = log10(nb[i]) ;
    }
    
    hmin = pow( double(10), (*logh)(0) ) ;
    hmax = pow( double(10), (*logh)(nbp-1) ) ;
    
    // Cleaning
    //---------
    delete [] press ; 
    delete [] nb ; 
    delete [] ro ; 

  }

  
			//--------------//
			//  Destructor  //
			//--------------//

Eos_consistent::~Eos_consistent(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_consistent::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_consistent !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_consistent::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_consistent::nbar_ent_p(double ent, const Param* ) const {

  static int i_near = logh->get_taille() / 2 ;
  
  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_consistent::nbar_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logh0 = log10( ent ) ;
    double logp0 ;
    double dlpsdlh0 ;
    interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
		  dlpsdlh0) ;
    
    double pp = pow(double(10), logp0) ;
    
    double resu = pp / ent * dlpsdlh0 * exp(-ent) ;
    if (i_near == 0) 
      { // Use of linear interpolation for the first interval
	double pp_near = pow(double(10), (*logp)(i_near)) ;
	double ent_near = pow(double(10), (*logh)(i_near)) ;
	resu = pp_near / ent_near * (*dlpsdlh)(i_near) * exp(-ent_near) ;
      }
    return resu ;
  }
  else{
    return 0 ;
  }
}

// Energy density from enthalpy
//------------------------------

double Eos_consistent::ener_ent_p(double ent, const Param* ) const {

  static int i_near = logh->get_taille() / 2 ;
  
  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_consistent::ener_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    double logh0 = log10( ent ) ;
    double logp0 ;
    double dlpsdlh0 ;
    interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
		  dlpsdlh0) ;
    
    double pp = pow(double(10), logp0) ;
    
    double resu = pp / ent * dlpsdlh0 - pp ;
    if (i_near == 0)
      {
	double p_near = pow(double(10), (*logp)(i_near)) ;
	double p_nearp1 = pow(double(10), (*logp)(i_near+1)) ;
	double ener_near = p_near/ pow(double(10), (*logh)(i_near))
	  * (*dlpsdlh)(i_near) - p_near ;
	double ener_nearp1 = p_nearp1/ pow(double(10), (*logh)(i_near+1))
	  * (*dlpsdlh)(i_near+1) - p_nearp1 ; 
	double delta = (*logh)(i_near+1) - (*logh)(i_near) ;
	resu = (ener_nearp1*(logh0 - (*logh)(i_near)) 
		 - ener_near*(logh0 - (*logh)(i_near+1))) / delta  ;
      }
    return resu ;
  }
  else{
    return 0 ;
  }
}

// Pressure from enthalpy
//------------------------

double Eos_consistent::press_ent_p(double ent, const Param* ) const {

  static int i_near = logh->get_taille() / 2 ;
  
  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_consistent::press_ent_p : ent > hmax !" << endl ;
      abort() ;
    }
    
    double logh0 = log10( ent ) ;
    double logp0 ;
    double dlpsdlh0 ;
    interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
		  dlpsdlh0) ;
    if (i_near == 0)
      {
	double logp_near = (*logp)(i_near) ;
	double logp_nearp1 = (*logp)(i_near+1) ;
	double delta = (*logh)(i_near+1) - (*logh)(i_near) ;
	logp0 = (logp_nearp1*(logh0 - (*logh)(i_near)) 
		 - logp_near*(logh0 - (*logh)(i_near+1))) / delta  ;
      }	   
    return pow(double(10), logp0) ;
  }
  else{
    return 0 ;
  }
}

// Square of sound speed from enthalpy
double Eos_consistent::csound_square_ent_p(double ent, const Param*) const {
  static int i_near = lognb->get_taille() / 2 ;
  
  if ( ent > hmin ) {
    if (ent > hmax) {
      cout << "Eos_consistent::csound_square_ent_p : ent>hmax !" << endl ;
      abort() ;
    }
    double log_ent0 = log10( ent ) ;
    double log_csound0 ;
        
    interpol_linear(*logh, *log_cs2, log_ent0, i_near, log_csound0) ;
      
    return pow(10., log_csound0) ;
  }
  else
    {
      return pow(10., (*log_cs2)(0)) ; 
    }
}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_consistent::operator>>(ostream & ost) const {

    ost << "EOS of class Eos_consistent." << endl ;
    ost << "Built from file " << tablename << endl ;
    ost << "Authors : " << authors << endl ;
    ost << "Number of points in file : " << logh->get_taille() << endl ;
    ost << "Table eventually slightly modified to ensure the relation" << endl ;
    ost << "dp = (e+p) dh" << endl ;
    return ost ;
}

			
}
