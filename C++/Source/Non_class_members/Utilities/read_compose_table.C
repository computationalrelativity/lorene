/*
 * Function to read a CompOSE table (general purpose or neutron star matter).
 *
 */

/*
 *   Copyright (c) 2024 Jerome Novak
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

// Lorene headers
#include "tbl.h"
#include "unites.h"


namespace Lorene {

  bool read_compose_table(const string& tablename, Tbl*& p_ental, Tbl*& p_entro,
			  Tbl*& p_press, Tbl*& p_ener, Tbl*& p_nb, Tbl*& p_temp,
			  Tbl*& p_ye) {

    bool resu = true ;
    if (p_ental != nullptr) delete p_ental ; 
    if (p_entro != nullptr) delete p_entro ; 
    if (p_press != nullptr) delete p_press ; 
    if (p_ener != nullptr) delete p_ener ; 
    if (p_nb != nullptr) delete p_nb ; 
    if (p_temp != nullptr) delete p_temp ; 
    if (p_ye != nullptr) delete p_ye ; 
    
    // Files containing data and a test
    //---------------------------------
    string file_nb = tablename + ".nb" ;
    string file_t = tablename + ".t" ;
    string file_ye = tablename + ".yq" ;
    string file_thermo = tablename + ".thermo" ;
  
    ifstream in_nb(file_nb.data()) ;
    ifstream in_t(file_t.data()) ;
    ifstream in_ye(file_ye.data()) ;
    if ((!in_nb)||(!in_t)||(!in_ye)) {
      cerr << "Problem in opening the EOS files!" << endl ;
      cerr << "While trying to open " << tablename << endl ;
      cerr << "Aborting ... " << endl ;
      abort() ;
    }

    // reading density, temperature & Ye tables
    //------------------------------------------
    int index1_n, index2_n ;
    in_nb >> index1_n >> index2_n ;
    int nbn = index2_n - index1_n + 1 ;
    p_nb = new Tbl(nbn) ; p_nb->set_etat_qcq() ;
    for (int i=0; i<nbn; i++) 
      in_nb >> p_nb->set(i) ;
    
    int index1_t, index2_t ;
    in_t >> index1_t >> index2_t ;
    int nbt = index2_t - index1_t + 1 ;
    p_temp = new Tbl(nbt) ; p_temp->set_etat_qcq() ;
    for (int j=0; j<nbt; j++)
      in_t >> p_temp->set(j) ;
  
    int index1_y, index2_y ;
    in_ye >> index1_y >> index2_y ;
    int nby = index2_y - index1_y + 1 ;
    p_ye = new Tbl(nby) ; p_ye->set_etat_qcq() ;
    for (int k=0; k<nby; k++)
      in_ye >> p_ye->set(k) ;

    // Arrays for other thermo quantities
    //------------------------------------
    p_ener = new Tbl(nby, nbn, nbt) ; p_ener->set_etat_qcq() ; (*p_ener) = NAN ;
    p_press = new Tbl(nby, nbn, nbt) ; p_press->set_etat_qcq() ;
    p_ental = new Tbl(nby, nbn, nbt) ; p_ental->set_etat_qcq() ; // Enthalpy
    p_entro = new Tbl(nby, nbn, nbt) ; p_entro->set_etat_qcq() ; // Entropy / baryon
  
    // Variables and conversion
    //-------------------------
    double e_compose, p_compose, p_over_nb_comp, s_per_bar, mub_shifted, eps_comp ;
    double dummy_x ;
    int dummy_n ;

    // Conversion from MeV/fm^3 to cgs
    double p_convert = Unites::mev_si * 1.e45 * 10. ;
    // From meV/fm^3 to g/cm^3
    double eps_convert = Unites::mev_si * 1.e42 / (Unites::c_si*Unites::c_si) ;
    double rhonuc_cgs = Unites::rhonuc_si * 1e-3 ;
    double c2_cgs = Unites::c_si * Unites::c_si * 1e4 ;

    double m_neutron_MeV, m_proton_MeV ;

    ifstream in_thermo(file_thermo.data()) ;
    if (!in_thermo) {
      cerr << "Problem in opening the EOS file!" << endl ;
      cerr << "While trying to open " << file_thermo << endl ;
      cerr << "Aborting ... " << endl ;
      abort() ;
    }

    // neutron and proton masses
    in_thermo >> m_neutron_MeV >> m_proton_MeV >> dummy_n;

    long int ntot = nby*nbn*nbt ;

    // Main loop reading the table
    //----------------------------
#ifndef NDEBUG
    cout << "Reading " << file_thermo << " ... " << flush ;
#endif

    for (int i=0; i<ntot; i++) {
#ifndef NDEBUG
      if (!in_thermo) {
	cerr << "Problem in EOS file: seems too small!" << endl ;
	cerr << "Name: " << file_thermo << endl ;
	cerr << "Aborting ... " << endl ;
	abort() ;
      }
#endif
      int it, inb, iye ;
      in_thermo >> it >> inb >> iye >> p_over_nb_comp ; // T, nB, Ye
      in_thermo >> s_per_bar >> mub_shifted >> dummy_x >> dummy_x >> dummy_x
		>> eps_comp ;
      in_thermo.ignore(1000, '\n') ;
      double nb0 = (*p_nb)(inb - index1_n) ;
      p_compose = p_over_nb_comp * nb0 ;
      e_compose = ( eps_comp + 1. ) * m_neutron_MeV * nb0 ;
      //double mub_compose = (mub_shifted + 1.)*m_neutron_MeV ;
      if ( (e_compose<0) || (p_compose < 0) ){
	cout << "Negative value in table!" << endl ;
	cout << "e = " << e_compose << ", p = " << p_compose << endl ;
	cout << "x = " << it  << ", y = " << inb << ", z = " << iye << endl ;      
	resu = false ;
      }
      p_ener->set(iye-index1_y, inb-index1_n, it-index1_t)
	= e_compose*eps_convert / rhonuc_cgs;
      p_press->set(iye-index1_y, inb-index1_n, it-index1_t)
	= p_compose * p_convert / (c2_cgs*rhonuc_cgs) ;
      double ent = ( e_compose + p_compose) / (m_neutron_MeV*nb0) ;
      p_ental->set(iye-index1_y, inb-index1_n, it-index1_t)
	= ent ;
      p_entro->set(iye-index1_y, inb-index1_n, it-index1_t)
	= s_per_bar ;
    }
    cout << "done! " << endl ;
#ifndef NDEBUG
    for (int i=0; i<ntot; i++) {
      if (isnan(p_ener->t[i])) {
	cerr << "Missing points in EoS table for " << file_thermo << endl ;
	resu = false ;
      }
    }
#endif  

    return resu ;
  } // end of read_compose_table()

} // end of namespace Lorene
