/*
 *  Methods of class Dyn_eos_tab
 *
 *  (see file dyneos.h for documentation).
 *
 */

/*
 *   Copyright (c) 2019 Jerome Novak
 *             (c) 2000 Eric Gourgoulhon for Eos classes
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
 * Revision 1.3  2022/03/22 13:36:00  j_novak
 * Added declaration of compute_derivative to utilitaires.h
 *
 * Revision 1.2  2020/12/17 17:00:27  j_novak
 * Output of sound speed squared, instead of sound speed.
 *
 * Revision 1.1  2019/12/06 14:30:50  j_novak
 * New classes Dyn_eos... for cold Eos's with baryon density as input.
 *
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "dyneos.h"
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"


namespace Lorene {
  void interpol_herm(const Tbl& , const Tbl&, const Tbl&, double, int&,
		     double&, double& ) ;
  
  void interpol_linear(const Tbl&, const Tbl&, double, int&, double&) ;
  

    
			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
  Dyn_eos_tab::Dyn_eos_tab(const string& name_i, const string& tablename_i,
			   bool compose) : Dyn_eos(name_i), tablename(tablename_i),
					   compose_format(compose), lognb(0x0),
					   loge(0x0), dlesdlnb(0x0), c_sound2(0x0)
  {
    if (compose_format)
      read_table_compose() ;
    else
      read_table_lorene() ;
  }

// Default constructor (protected, to be used by derived classes)
// ---------------------------------------------------------------
  Dyn_eos_tab::Dyn_eos_tab() : Dyn_eos()
  {
  }
  
// Constructor from binary file
// ----------------------------
  Dyn_eos_tab::Dyn_eos_tab(FILE* fich) : Dyn_eos(fich), lognb(0x0),
					   loge(0x0), dlesdlnb(0x0), c_sound2(0x0)
  {
    const int nc = 160 ;
    char tmp_string[nc] ;
    size_t ret = fread(tmp_string, sizeof(char), nc, fich) ;
    if (int(ret) == nc)
      tablename = tmp_string ;
    int comp ;
    fread_be(&comp, sizeof(int), 1, fich) ;
    compose_format = comp ;
    
    if (compose_format)
      read_table_compose() ;
    else
      read_table_lorene() ;
  }


// Constructor from a formatted file
// ---------------------------------
  Dyn_eos_tab::Dyn_eos_tab(ifstream& fich) : Dyn_eos(fich), lognb(0x0),
					   loge(0x0), dlesdlnb(0x0), c_sound2(0x0)
  {  
    fich.seekg(0, fich.beg) ;
    fich.ignore(1000, '\n') ;
    fich >> compose_format ;
    fich.ignore(1000, '\n') ;
    getline(fich, name, '\n') ;
    fich >> tablename ;
	
    if (compose_format)
      read_table_compose() ;
    else
      read_table_lorene() ;
  }

			//--------------//
			//  Destructor  //
			//--------------//

  Dyn_eos_tab::~Dyn_eos_tab()
  {
    if (lognb != 0x0) delete lognb ;
    if (loge != 0x0) delete loge ;
    if (dlesdlnb != 0x0) delete dlesdlnb ;
    if (c_sound2 != 0x0) delete c_sound2 ;
  }

			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Dyn_eos_tab::operator==(const Dyn_eos& eos_i) const {
    
    bool resu = true ;
    
    if ( eos_i.identify() != identify() ) {
      cout << "The second EOS is not of type Dyn_eos_tab !" << endl ;
      resu = false ;
    }
    
    return resu ;
    
  }

  bool Dyn_eos_tab::operator!=(const Dyn_eos& eos_i) const {
    
    return !(operator==(eos_i)) ;
    
  }

			//------------//
			//  Outputs   //
			//------------//

  void Dyn_eos_tab::sauve(FILE* fich) const
  {
    Dyn_eos::sauve(fich) ;
  
    char tmp_string[160] ;
    strcpy(tmp_string, tablename.c_str()) ;
    fwrite(tmp_string, sizeof(char), 160, fich) ;
    int comp = int(compose_format) ;
    fwrite_be(&comp, sizeof(int), 1, fich) ;
  }

  ostream& Dyn_eos_tab::operator>>(ostream & ost) const
  {
    ost << "EOS of class Dyn_eos_tab." << endl ;
    ost << "Built from file " << tablename << endl ;
    ost << ((compose_format == 0) ? "Old LORENE format" : "CompOSE format") << endl ; 
    ost << "Authors : " << authors << endl ;
    ost << "Number of points in file : " << lognb->get_taille() << endl ;
    return ost ;
  }
			//------------------------//
			//  Reading of the table  //
			//------------------------//

  // LORENE format
  //---------------
  void Dyn_eos_tab::read_table_lorene() {

    using namespace Unites ;
    
    ifstream fich(tablename.data()) ;
    
    if (!fich) {
      cerr << "Dyn_eos_tab::read_table_lorene(): " << endl ;
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
    
    int nbp ;
    fich >> nbp ; fich.ignore(1000, '\n') ;   // number of data
    if (nbp<=0) {
      cerr << "Dyn_eos_tab::read_table_lorene(): " << endl ;
      cerr << "Wrong value for the number of lines!" << endl ;
      cerr << "nbp = " << nbp << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }
    
    for (int i=0; i<3; i++) {		//  jump over the table
      fich.ignore(1000, '\n') ;             // description
    }                                      
    
    Tbl press(nbp) ; press.set_etat_qcq() ;
    Tbl nb(nbp) ; nb.set_etat_qcq() ;
    Tbl ro(nbp) ; ro.set_etat_qcq() ;
    
    lognb = new Tbl(nbp) ;
    loge = new Tbl(nbp) ;
    dlesdlnb = new Tbl(nbp) ;
    
    lognb->set_etat_qcq() ;
    loge->set_etat_qcq() ;
    dlesdlnb->set_etat_qcq() ;
    
    double rhonuc_cgs = rhonuc_si * 1e-3 ;
    double c2_cgs = c_si * c_si * 1e4 ;    	
    
    int no ;
    double nb_fm3, rho_cgs, p_cgs ;
    
    cout << "Dyn_eos_tab: reading Lorene format table from the file "
	 << tablename << endl ;
    for (int i=0; i<nbp; i++) {
      
      fich >> no ;
      fich >> nb_fm3 ;
      fich >> rho_cgs ;
      fich >> p_cgs ; fich.ignore(1000,'\n') ;    		
      if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	cout << "Dyn_eos_tab::read_table_lorene(): " << endl ;
	cout << "Negative value in table!" << endl ;
	cout << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	  ", p = " << p_cgs << endl ;
	cout << "Aborting..." << endl ;
	abort() ;
      }
      
      press.set(i) = p_cgs / (c2_cgs*rhonuc_cgs) ;
      nb.set(i)    = 10.*nb_fm3 ; // Units 0.1 fm^-3
      ro.set(i)    = rho_cgs / rhonuc_cgs ; 
    }
    
    *lognb = log10(nb) ;
    *loge = log10(ro) ;
    *dlesdlnb = (ro + press) / ro ;        
    Tbl tmp(nbp) ; tmp.set_etat_qcq() ;
    compute_derivative(ro, press, tmp) ;
    c_sound2 = new Tbl(tmp) ; // c_s^2 = dp/de
    
    nbmin = pow( double(10), (*lognb)(0) ) ;
    nbmax = pow( double(10), (*lognb)(nbp-1) ) ;
    
    cout << "nbmin, nbmax : " << nbmin << "  " << nbmax << endl ;
    
    fich.close();  
  }

  // CompOSE format
  //----------------
  void Dyn_eos_tab::read_table_compose()
  {
    using namespace Unites ;
    
    // Files containing data and a test
    //---------------------------------
    string file_nb = tablename + ".nb" ;
    string file_thermo = tablename + ".thermo" ;

    ifstream in_nb(file_nb.data()) ;
    if (!in_nb) {
      cerr << "Dyn_eos_tab::read_table_compose(): " << endl ;
      cerr << "Problem in opening the EOS file!" << endl ;
      cerr << "While trying to open " << file_nb << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }

    // obtaining the size of the tables for memory allocation
    //-------------------------------------------------------
    int index1, index2 ;
    in_nb >> index1 >> index2 ;
    int nbp = index2 - index1 + 1 ;
    assert(nbp > 0) ;

    Tbl press(nbp) ; press.set_etat_qcq() ;
    Tbl nb(nbp) ; nb.set_etat_qcq() ;
    Tbl ro(nbp) ; ro.set_etat_qcq() ;
    
    lognb = new Tbl(nbp) ;
    loge = new Tbl(nbp) ;
    dlesdlnb = new Tbl(nbp) ;
    
    lognb->set_etat_qcq() ;
    loge->set_etat_qcq() ;
    dlesdlnb->set_etat_qcq() ;

    // Variables and conversion
    //-------------------------
    double nb_fm3, rho_cgs, p_cgs, p_over_nb_comp, eps_comp ;
    double dummy_x ;
    int dummy_n ;
    
    double rhonuc_cgs = rhonuc_si * 1e-3 ;
    double c2_cgs = c_si * c_si * 1e4 ;
    double m_neutron_MeV, m_proton_MeV ;
    
    ifstream in_p_rho (file_thermo.data()) ;
    if (!in_p_rho) {
      cerr << "Dyn_eos_tab::read_table_compose(): " << endl ;
      cerr << "Problem in opening the EOS file!" << endl ;
      cerr << "While trying to open " << file_thermo << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }
    in_p_rho >> m_neutron_MeV >> m_proton_MeV ; //Neutron and proton masses
    in_p_rho.ignore(1000, '\n') ;
    
    double p_convert = mev_si * 1.e45 * 10. ; // Conversion from MeV/fm^3 to cgs
    double eps_convert = mev_si * 1.e42 / (c_si*c_si) ; //From meV/fm^3 to g/cm^3
    
    // Main loop reading the table
    //----------------------------
    cout << "Dyn_eos_tab: reading CompOSE format table from the file "
	 << tablename << ".thermo" << endl ;
    for (int i=0; i<nbp; i++) {
      in_nb >> nb_fm3 ;
      in_p_rho >> dummy_n >> dummy_n >> dummy_n >> p_over_nb_comp ;
      in_p_rho >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> eps_comp ;
      in_p_rho.ignore(1000, '\n') ;
      p_cgs = p_over_nb_comp * nb_fm3 * p_convert ;
      rho_cgs = ( eps_comp + 1. ) * m_neutron_MeV * nb_fm3 * eps_convert ;
      
      if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	cout << "Dyn_eos_tab::read_table_compose(): " << endl ;
	cout << "Negative value in table!" << endl ;
	cout << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	  ", p = " << p_cgs << endl ;
	cout << "Aborting..." << endl ;
	abort() ;
      }
      
      press.set(i) = (p_cgs / c2_cgs) / rhonuc_cgs ;
      nb.set(i) = 10. * nb_fm3 ; // Units : 0.1 fm^-3
      ro.set(i)    = rho_cgs / rhonuc_cgs ;
    }
    
    *lognb = log10(nb) ;
    *loge = log10(ro) ;
    *dlesdlnb = (ro + press) / ro ;    
    Tbl tmp(nbp) ; tmp.set_etat_qcq() ;
    compute_derivative(ro, press, tmp) ;
    c_sound2 = new Tbl(tmp) ; // c_s^2 = dp/de
    
    nbmin = pow( double(10), (*lognb)(0) ) ;
    nbmax = pow( double(10), (*lognb)(nbp-1) ) ;
    
    cout << "nbmin, nbmax : " << nbmin << "  " << nbmax << endl ;
  }
			//------------------------------//
			//    Computational routines    //
			//------------------------------//

  // Enthalpy from baryon density
  //------------------------------
  double Dyn_eos_tab::ent_nbar_p(double nbar, const Param* ) const
  {
    static int i_near = lognb->get_taille() / 2 ;

    if ( nbar > nbmin ) {
      if (nbar > nbmax) {
	cout << "Dyn_eos_tab::ent_nbar_p : nbar > nbmax !" << endl ;
	abort() ;
      }
      double lognb0 = log10( nbar ) ;
      double loge0 ;
      double dlesdlnb0 ;
      interpol_herm(*lognb, *loge, *dlesdlnb, lognb0, i_near, loge0,
		    dlesdlnb0) ;
      double ee = pow(double(10), loge0) ;
      double resu = dlesdlnb0*ee / nbar ;
      return log(resu) ;
    }
    else
      return log((*dlesdlnb)(0)*pow(10.,(*loge)(0)) / nbmin) ;
  }

  // Energy density from baryon density
  //------------------------------------

  double Dyn_eos_tab::ener_nbar_p(double nbar, const Param* ) const
  {
    static int i_near = lognb->get_taille() / 2 ;

    if ( nbar > nbmin ) {
      if (nbar > nbmax) {
	cout << "Dyn_eos_tab::ener_nbar_p : nbar > nbmax !" << endl ;
	abort() ;
      }
      double lognb0 = log10( nbar ) ;
      double loge0 ;
      double dlesdlnb0 ;
      interpol_herm(*lognb, *loge, *dlesdlnb, lognb0, i_near, loge0,
		    dlesdlnb0) ;
      return pow(double(10), loge0) ;
    }
    else
      return pow(10.,(*loge)(0)) ;
  }

  // Pressure from baryon density
  //------------------------------

  double Dyn_eos_tab::press_nbar_p(double nbar, const Param* ) const
  {
    static int i_near = lognb->get_taille() / 2 ;

    if ( nbar > nbmin ) {
      if (nbar > nbmax) {
	cout << "Dyn_eos_tab::press_nbar_p : nbar > nbmax !" << endl ;
	abort() ;
      }
      double lognb0 = log10( nbar ) ;
      double loge0 ;
      double dlesdlnb0 ;
      interpol_herm(*lognb, *loge, *dlesdlnb, lognb0, i_near, loge0,
		    dlesdlnb0) ;
      double ee = pow(double(10), loge0) ;
      double hnb = ee * dlesdlnb0 ;
      return hnb - ee ;
    }
    else{
      return pow(10.,(*loge)(0))*((*dlesdlnb)(0) - 1.) ;
    }
  }


  // Sound speed from baryon density 
  //---------------------------------

  double Dyn_eos_tab::csound_square_nbar_p(double nbar, const Param*) const {

    static int i_near = lognb->get_taille() / 2 ;

    if ( nbar > nbmin ) {
      if (nbar > nbmax) {
	cout << "Dyn_eos_tab::press_nbar_p : nbar > nbmax !" << endl ;
	abort() ;
      }
      double lognb0 = log10( nbar ) ;
      double csound0 ;
    
      interpol_linear(*lognb, *c_sound2, lognb0, i_near, csound0) ;
    
      return csound0 ;
    }
    else
      {
	return (*c_sound2)(0) ; 
      }
  
  
  }
}
