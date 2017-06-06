/*
 *  Methods of class Hoteos_tabul
 *
 *  (see file hoteos.h for documentation).
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
 * Revision 1.4  2017/06/06 15:36:42  j_novak
 * Reads a number as first column of tabulated EoS
 *
 * Revision 1.3  2016/12/05 16:17:52  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2015/12/08 15:42:17  j_novak
 * Low/zero entropy is set to the lowest value in the table in computational functions.
 *
 * Revision 1.1  2015/12/08 10:52:18  j_novak
 * New class Hoteos_tabul for tabulated temperature-dependent EoSs.
 *
 *
 *
 * $Header$
 *
 */

// headers C
#include <cstdlib>
#include <cstring>
#include <cmath>

// Headers Lorene
#include "hoteos.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"


namespace Lorene {
  void interpol_herm_2d(const Tbl&, const Tbl&, const Tbl&, const Tbl&, const Tbl&, 
			const Tbl&, double, double, double&, double&, double&) ;


			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

  // Standard constructor
  // --------------------	
  Hoteos_tabul::Hoteos_tabul(const string& filename): 
    Hot_eos("Tabulated hot EoS"), tablename(filename), authors("Unknown"), 
    hmin(0.), hmax(1.), sbmin(0.), sbmax(1.)
  {
    set_arrays_0x0() ;
    read_table() ;
  }
  
  
  // Constructor from binary file
  // ----------------------------
  Hoteos_tabul::Hoteos_tabul(FILE* fich) : Hot_eos(fich) {
    
    char tmp_string[160] ;
    fread(tmp_string, sizeof(char), 160, fich) ;
    tablename = tmp_string ;
    set_arrays_0x0() ;
    read_table() ;  
  }
  
  // Constructor from a formatted file
  // ---------------------------------
  Hoteos_tabul::Hoteos_tabul(ifstream& fich) : 
    Hot_eos(fich) {
    
    fich >> tablename ;
    set_arrays_0x0() ;
    read_table() ; 	 
  }

  //Sets the arrays to the null pointer
  void Hoteos_tabul::set_arrays_0x0() {
    hhh = 0x0 ;
    s_B = 0x0 ;
    ppp = 0x0 ;
    dpdh = 0x0 ;
    dpds = 0x0 ;
    d2p = 0x0 ;
  }

			//--------------//
			//  Destructor  //
			//--------------//

  Hoteos_tabul::~Hoteos_tabul(){
    if (hhh != 0x0) delete hhh ;
    if (s_B != 0x0) delete s_B ;
    if (ppp != 0x0) delete ppp ;
    if (dpdh != 0x0) delete dpdh ;
    if (dpds != 0x0) delete dpds ;
    if (d2p != 0x0) delete d2p ;
  }

			//------------//
			//  Outputs   //
			//------------//

  void Hoteos_tabul::sauve(FILE* fich) const {
    
    Hot_eos::sauve(fich) ;
    
    char tmp_string[160] ;
    strcpy(tmp_string, tablename.c_str()) ;
    fwrite(tmp_string, sizeof(char), 160, fich) ;		
  }

  ostream& Hoteos_tabul::operator>>(ostream & ost) const {
    
    ost << "Hot EOS of class Hoteos_tabul (tabulated hot beta-equilibrium EoS) : "
	<< endl ; 
    ost << "Built from file " << tablename << endl ;
    ost << "Authors : " << authors << endl ;
    ost << "Number of points in file : " << hhh->get_dim(0) 
	<< " in enthalpy, and " << hhh->get_dim(1) 
	<< " in entropy." << endl ;

    return ost ;
}

			//------------------------//
			//  Reading of the table  //
			//------------------------//
			
  void Hoteos_tabul::read_table() {

    using namespace Unites ;
    
    ifstream fich(tablename.data()) ;
    
    if (!fich) {
      cerr << "Hoteos_tabul::read_table(): " << endl ;
      cerr << "Problem in opening the EOS file!" << endl ;
      cerr << "While trying to open " << tablename << endl ; 
      cerr << "Aborting..." << endl ;
      abort() ;
    }
    
    fich.ignore(1000, '\n') ;             // Jump over the first header
    fich.ignore(1) ;
    getline(fich, authors, '\n') ;
    for (int i=0; i<3; i++) {		//  jump over the file
      fich.ignore(1000, '\n') ;             // header
    }                                       //
    
    int nbp1, nbp2 ;
    fich >> nbp1 >> nbp2 ; fich.ignore(1000, '\n') ;   // number of data
    if ( (nbp1<=0) || (nbp2<=0) ) {
      cerr << "Hoteos_tabul::read_table(): " << endl ;
      cerr << "Wrong value for the number of lines!" << endl ;
      cerr << "nbp1 = " << nbp1 << ", nbp2 = " << nbp2 << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }
    
    for (int i=0; i<3; i++) {		//  jump over the table
      fich.ignore(1000, '\n') ;             // description
    }                                      
    
    ppp = new Tbl(nbp2, nbp1) ;
    hhh = new Tbl(nbp2, nbp1) ;
    s_B = new Tbl(nbp2, nbp1) ;
    dpdh = new Tbl(nbp2, nbp1) ;
    dpds = new Tbl(nbp2, nbp1) ;
    d2p = new Tbl(nbp2, nbp1) ;
    
    ppp->set_etat_qcq() ;
    hhh->set_etat_qcq() ;
    s_B->set_etat_qcq() ;
    dpdh->set_etat_qcq() ;
    dpds->set_etat_qcq() ;
    d2p->set_etat_qcq() ;
    
    double c2 = c_si * c_si ;
    double dummy, nb_fm3, rho_cgs, p_cgs, mu_MeV, entr, temp, der2 ;	
    double ww = 0. ;
    int no;

    for (int j=0; j<nbp2; j++) {
      for (int i=0; i<nbp1; i++) {
	fich >> no >> nb_fm3>> rho_cgs >> p_cgs>> mu_MeV >> entr >> temp >> der2 
	     >> dummy;
	fich.ignore(1000,'\n') ;
	if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	  cerr << "Eos_mag::read_table(): " << endl ;
	  cerr << "Negative value in table!" << endl ;
	  cerr << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	    ", p = " << p_cgs << endl ;
	  cerr << "Aborting..." << endl ;
	  abort() ;
	}
      
	double psc2 = 0.1*p_cgs/c2 ; //  in kg m^-3
	double rho_si = rho_cgs*1000. ; // in kg m^-3
	
	double h_read = log(mu_MeV) ;
	if ( (i==0) && (j==0) ) ww = h_read ;
	double h_new = h_read - ww ;
	
	ppp->set(j, i) = psc2/rhonuc_si ; 
	hhh->set(j, i) = h_new ;
	s_B->set(j, i) = entr ; // in Lorene units (k_B)
	dpdh->set(j, i) = (rho_si + psc2)/rhonuc_si ; 
	dpds->set(j, i) = -temp*nb_fm3*mevpfm3 ; 
	d2p->set(j, i) = 0.1*der2*mu_MeV/(c2*rhonuc_si) ;	
      }
    }

    hmin = (*hhh)(0, 0) ;
    hmax = (*hhh)(0, nbp1-1) ;
    
    sbmin = (*s_B)(0, 0) ;
    sbmax = (*s_B)(nbp2-1, 0) ;
    
    cout << "hmin: " << hmin << ", hmax: " << hmax << endl ;
    cout << "sbmin: " << sbmin << ", sbmax: " << sbmax << endl ;
    
    fich.close();
}

			//-------------------------------//
			//  The corresponding cold Eos   //
			//-------------------------------//

  const Eos& Hoteos_tabul::new_cold_Eos() const {
    
    if (p_cold_eos == 0x0) {
      cerr << "Warning: Hoteos_tabul::new_cold_Eos " <<
	"The corresponding cold EoS is likely not to work" << endl ;
      cout << "read from file:   "<< tablename.c_str() << endl;
      p_cold_eos = new Eos_CompOSE(tablename.c_str()) ;
    }
    
    return *p_cold_eos ;
  }



			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Hoteos_tabul::operator==(const Hot_eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
      cout << "The second EOS is not of type Hoteos_tabul !" << endl ; 
      resu = false ; 
    }

    return resu ; 
  }


  bool Hoteos_tabul::operator!=(const Hot_eos& eos_i) const {
    return !(operator==(eos_i)) ;  
  }
  
  
			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy and entropy
//------------------------------------------

double Hoteos_tabul::nbar_Hs_p(double ent, double sb) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if (sb < sbmin) sb = sbmin ;

  if ( ent >= hmin ) {
    if (ent > hmax) {
      cout << "Hoteos_tabul::nbar_Hs_p : ent > hmax !" << endl ;
      abort() ;
    }

    if (sb > sbmax) {
      cerr << "Hoteos_tabul::nbar_Hs_p : s_B not in the tabulated interval !" 
	   << endl ;
      cerr << "s_B = " << sb << ", sbmin = " << sbmin << ", sbmax = " << sbmax 
	   << endl ;
      abort() ;
    }

    double p_int, dpds_int, dpdh_int ;
    interpol_herm_2d(*s_B, *hhh, *ppp, *dpds, *dpdh, *d2p, sb, ent, p_int,
		     dpds_int, dpdh_int) ;

    double nbar_int = dpdh_int * exp(-ent) ;
    return nbar_int ;
  }
  else{
    return 0 ;
  }
}

// Energy density from enthalpy and entropy
//-----------------------------------------

double Hoteos_tabul::ener_Hs_p(double ent, double sb) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if (sb < sbmin) sb = sbmin ;

  if ( ent >= hmin ) {
    if (ent > hmax) {
      cout << "Hoteos_tabul::ener_Hs_p : ent > hmax !" << endl ;
      abort() ;
    }

    if (sb > sbmax) {
      cerr << "Hoteos_tabul::ener_Hs_p : s_B not in the tabulated interval !" 
	   << endl ;
      cerr << "s_B = " << sb << ", sbmin = " << sbmin << ", sbmax = " << sbmax 
	   << endl ;
      abort() ;
    }

    double p_int, dpds_int, dpdh_int ;
    interpol_herm_2d(*s_B, *hhh, *ppp, *dpds, *dpdh, *d2p, sb, ent, p_int,
		     dpds_int, dpdh_int) ;

    double nbar_int = dpdh_int * exp(-ent) ;

    double f_int = - p_int + exp(ent) * nbar_int;

    return f_int ;
  }
  else{
    return 0 ;
  }
}

// Pressure from enthalpy and entropy
//-----------------------------------

double Hoteos_tabul::press_Hs_p(double ent, double sb) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if (sb < sbmin) sb = sbmin ;

  if ( ent >= hmin ) {
    if (ent > hmax) {
      cout << "Hoteos_tabul::press_Hs_p : ent > hmax !" << endl ;
      abort() ;
    }

    if (sb > sbmax) {
      cerr << "Hoteos_tabul::press_Hs_p : s_B not in the tabulated interval !" 
	   << endl ;
      cerr << "s_B = " << sb << ", sbmin = " << sbmin << ", sbmax = " << sbmax 
	   << endl ;
      abort() ;
    }

    double p_int, dpds_int, dpdh_int ;
    interpol_herm_2d(*s_B, *hhh, *ppp, *dpds, *dpdh, *d2p, sb, ent, p_int,
		     dpds_int, dpdh_int) ;

    return p_int ;
    }
    else{
      return 0 ;
    }
}

  // Temperature from enthalpy and entropy
  //--------------------------------------
  double Hoteos_tabul::temp_Hs_p(double ent, double sb) const  {
    
    if ((ent > hmin - 1.e-12) && (ent < hmin))
      ent = hmin ;

    if (sb < sbmin) sb = sbmin ;
    
    if ( ent >= hmin ) {
      if (ent > hmax) {
	cout << "Hoteos_tabul::temp_Hs_p : ent > hmax !" << endl ;
	abort() ;
      }
      
      if (sb > sbmax) {
	cerr << "Hoteos_tabul::temp_Hs_p : s_B not in the tabulated interval !" 
	     << endl ;
	cerr << "s_B = " << sb << ", sbmin = " << sbmin << ", sbmax = " << sbmax 
	     << endl ;
	abort() ;
      }
      
      double p_int, dpds_int, dpdh_int ;
      interpol_herm_2d(*s_B, *hhh, *ppp, *dpds, *dpdh, *d2p, sb, ent, p_int,
		       dpds_int, dpdh_int) ;
      
      double nbar_int = dpdh_int * exp(-ent) ;
      
      double temp_int = -dpds_int / nbar_int ;
      
      return temp_int ;
    }
    else {
      return 0 ;
    }
  }

} //End of namespace Lorene
