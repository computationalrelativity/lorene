/*
 *  Methods of class Dyn_eos_cons
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
    
  void compute_derivative(const Tbl&, const Tbl&, Tbl&) ;

                        //----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
  Dyn_eos_cons::Dyn_eos_cons(const string& name_i, const string& tablename_i,
			   bool compose) : Dyn_eos_tab()
  {
    name = name_i ;
    tablename = tablename_i ;
    compose_format = compose ;
    if (compose_format)
      read_table_compose() ;
    else
      read_table_lorene() ;
  }

// Constructor from binary file
// ----------------------------
  Dyn_eos_cons::Dyn_eos_cons(FILE* fich) : Dyn_eos_tab()
  {
    const int nc1 = 100 ;
    const int nc2 = 160 ;
    char tmp_s1[nc1] ;
    size_t ret = fread(tmp_s1, sizeof(char), nc1, fich) ;
    if (int(ret) == nc1)
      name = tmp_s1 ;
    char tmp_s2[nc2] ;
    ret = fread(tmp_s2, sizeof(char), nc2, fich) ;
    if (ret == nc2)
      tablename = tmp_s2 ;
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
  Dyn_eos_cons::Dyn_eos_cons(ifstream& fich) : Dyn_eos_tab()
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

  Dyn_eos_cons::~Dyn_eos_cons()
  {
   }

			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Dyn_eos_cons::operator==(const Dyn_eos& eos_i) const {
    
    bool resu = true ;
    
    if ( eos_i.identify() != identify() ) {
      cout << "The second EOS is not of type Dyn_eos_cons !" << endl ;
      resu = false ;
    }
    
    return resu ;
    
  }

  bool Dyn_eos_cons::operator!=(const Dyn_eos& eos_i) const {
    
    return !(operator==(eos_i)) ;
    
  }

			//------------//
			//  Outputs   //
			//------------//

  ostream& Dyn_eos_cons::operator>>(ostream & ost) const
  {
    ost << "EOS of class Dyn_eos_cons." << endl ;
    ost << "Built from file " << tablename << endl ;
    ost << ((compose_format == 0) ? "LORENE format" : "CompOSE format") << endl ; 
    ost << "Authors : " << authors << endl ;
    ost << "Number of points in file : " << lognb->get_taille() << endl ;
    ost << "Table eventually slightly modified to ensure the relation" << endl ;
    ost << "de/dn = (e+p)/n" << endl ;
    return ost ;
  }
			//------------------------//
			//  Reading of the table  //
			//------------------------//

  // LORENE format
  //---------------
  void Dyn_eos_cons::read_table_lorene() {

    using namespace Unites ;
    
    ifstream fich(tablename.data()) ;
    
    if (!fich) {
      cerr << "Dyn_eos_cons::read_table_lorene(): " << endl ;
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
      cerr << "Dyn_eos_cons::read_table_lorene(): " << endl ;
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
    
    cout << "Dyn_eos_cons: reading Lorene format table from the file "
	 << tablename << endl ;
    for (int i=0; i<nbp; i++) {
      
      fich >> no ;
      fich >> nb_fm3 ;
      fich >> rho_cgs ;
      fich >> p_cgs ; fich.ignore(1000,'\n') ;    		
      if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	cout << "Dyn_eos_cons::read_table_lorene(): " << endl ;
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

    Tbl ee(nbp) ; ee.set_etat_qcq() ;
    Tbl dn(nbp) ; dn.set_etat_qcq() ;
    ee = log(ro) ;
    dn = ro / ( ro + press ) ;

    Tbl nb2 = exp( integ1D(ee, dn)  + log(nb(0)) ) ;
    double maxdif = diffrelmax(nb2, nb) ;
    cout << "Dyn_eos_cons: max. relative difference between new (consistent)" << endl ;
    cout << "and old density data: " << maxdif << endl ;
   
    *lognb = log10(nb2) ;
    *loge = log10(ro) ;
    *dlesdlnb = (ro + press) / ro ;    
    Tbl tmp(nbp) ; tmp.set_etat_qcq() ;
    compute_derivative(ro, press, tmp) ;
    c_sound2 = new Tbl(tmp) ; // c_s^2 = dp/de
    
    nbmin = pow( double(10), (*lognb)(0) ) ;
    nbmax = pow( double(10), (*lognb)(nbp-1) ) ;
    
    cout << "nbmin, nbmax : " << 0.1*nbmin << " " << 0.1*nbmax << " fm^-3" << endl ;
    
    fich.close();  
  }

  // CompOSE format
  //----------------
  void Dyn_eos_cons::read_table_compose()
  {
    using namespace Unites ;
    
    // Files containing data and a test
    //---------------------------------
    string file_nb = tablename + ".nb" ;
    string file_thermo = tablename + ".thermo" ;

    ifstream in_nb(file_nb.data()) ;
    if (!in_nb) {
      cerr << "Dyn_eos_cons::read_table_compose(): " << endl ;
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
      cerr << "Dyn_eos_cons::read_table_compose(): " << endl ;
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
    cout << "Dyn_eos_cons: reading CompOSE format table from the file "
	 << tablename << ".thermo" << endl ;
    for (int i=0; i<nbp; i++) {
      in_nb >> nb_fm3 ;
      in_p_rho >> dummy_n >> dummy_n >> dummy_n >> p_over_nb_comp ;
      in_p_rho >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> eps_comp ;
      in_p_rho.ignore(1000, '\n') ;
      p_cgs = p_over_nb_comp * nb_fm3 * p_convert ;
      rho_cgs = ( eps_comp + 1. ) * m_neutron_MeV * nb_fm3 * eps_convert ;
      
      if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	cout << "Dyn_eos_cons::read_table_compose(): " << endl ;
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
    
    Tbl ee(nbp) ; ee.set_etat_qcq() ;
    Tbl dn(nbp) ; dn.set_etat_qcq() ;
    ee = log(ro) ;
    dn = ro / ( ro + press ) ;

    Tbl nb2 = exp( integ1D(ee, dn)  + log(nb(0)) ) ;
    double maxdif = diffrelmax(nb2, nb) ;
    cout << "Dyn_eos_cons: max. relative difference between new (consistent)" << endl ;
    cout << "and old density data: " << maxdif << endl ;
   
    *lognb = log10(nb2) ;
    *loge = log10(ro) ;
    *dlesdlnb = (ro + press) / ro ;    
    Tbl tmp(nbp) ; tmp.set_etat_qcq() ;
    compute_derivative(ro, press, tmp) ;
    c_sound2 = new Tbl(tmp) ; // c_s^2 = dp/de
    
    nbmin = pow( double(10), (*lognb)(0) ) ;
    nbmax = pow( double(10), (*lognb)(nbp-1) ) ;
    
    cout << "nbmin, nbmax : " << nbmin << "  " << nbmax << endl ;
  }

}
