/*
 *  Methods of class Eos_bf_tabul.
 *
 *  (see file eos_bifluid.h for documentation).
 *
 */

/*
 *   Copyright (c) 2015 Aurelien Sourie
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
 * Revision 1.5  2017/10/06 12:36:33  a_sourie
 * Cleaning of tabulated 2-fluid EoS class + superfluid rotating star model.
 *
 * Revision 1.4  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2015/06/11 14:41:59  a_sourie
 * Corrected minor bug
 *
 * Revision 1.2  2015/06/11 13:50:19  j_novak
 * Minor corrections
 *
 * Revision 1.1  2015/06/10 14:39:17  a_sourie
 * New class Eos_bf_tabul for tabulated 2-fluid EoSs and associated functions for the computation of rotating stars with such EoSs.
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
#include "eos_bifluid.h"
#include "param.h"
#include "eos.h" 
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"
#include "cmp.h"
#include "nbr_spx.h"

namespace Lorene {


void interpol_herm(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
		   double x, int& i, double& y, double& dy) ;


			//----------------------------//
			//   	Constructors	      	//
			//----------------------------//


// Standard constructor
// --------------------			
Eos_bf_tabul::Eos_bf_tabul(const char* name_i, const char* table,
			   const char* path, double mass1, double mass2) : 
Eos_bifluid(name_i, mass1, mass2) 
{
	tablename = path ;
	tablename += "/" ;
	tablename += table ;
	
	read_table() ; 	
}

// Standard constructor with full filename
// ---------------------------------------
Eos_bf_tabul::Eos_bf_tabul(const char* name_i, const char* file_name, 
			   double mass1, double mass2) :
Eos_bifluid(name_i, mass1, mass2)
{
	tablename = file_name ;
	
	read_table() ; 	
}

// Constructor from binary file
// ----------------------------
Eos_bf_tabul::Eos_bf_tabul(FILE* fich) : 
Eos_bifluid(fich) 
{
  char tmp_string[160] ;
  fread(tmp_string, sizeof(char), 160, fich) ;
  tablename = tmp_string ;

  read_table() ;
}

// Constructor from a formatted file
// ---------------------------------
Eos_bf_tabul::Eos_bf_tabul(ifstream& fich, const char* table) : 
Eos_bifluid(fich) 
{
	fich >> tablename ;
   tablename += "/" ;
   tablename += table ;

   read_table() ;    
}

Eos_bf_tabul::Eos_bf_tabul(ifstream& fich) : 
Eos_bifluid(fich) 
{
   fich.ignore(1000, '\n') ;
   fich >> tablename ; 		//  directory in which the EoSs are stored
   tablename += name ; 		// the name of the EoS is set thanks to the call to Eos_bifluid(fich) 
   tablename += "_2f.d" ; 	// table for the two-fluid EoS
   read_table() ;    

}
                    

			//--------------//
			//  Destructor  //
			//--------------//

Eos_bf_tabul::~Eos_bf_tabul(){

	delete mu1_tab ;
	delete mu2_tab ;
	delete delta_car_tab ; 
	delete press_tab ;
	delete n1_tab;
   delete n2_tab ;
	delete d2psdmu1dmu2_tab   ;
	delete dpsddelta_car_tab ;
	delete dn1sddelta_car_tab ;
	delete dn2sddelta_car_tab ;
	delete delta_car_n0 ;
	delete mu1_n0 ;
	delete mu2_n0 ;
	delete delta_car_p0 ;
	delete mu1_p0 ;
	delete mu2_p0 ;
	delete mu1_N ;
	delete n_n_N ;
	delete press_N ;
	delete mu2_P ;
	delete n_p_P ;
	delete press_P ;

}

			//--------------//
			//  Assignment  //
			//--------------//

void Eos_bf_tabul::operator=(const Eos_bf_tabul& eosi) {
    
  	// Assignment of quantities common to all the derived classes of Eos_bifluid
    Eos_bifluid::operator=(eosi) ;	    
    
    tablename 				= eosi.tablename; 
    mu1_min 				= eosi.mu1_min ;
    mu1_max 				= eosi.mu1_max ; 
    mu2_min 				= eosi.mu2_min ; 
    mu2_max 				= eosi.mu2_max ; 
    delta_car_min 		= eosi.delta_car_min ;
    delta_car_max 		= eosi.delta_car_max ;
    mu1_tab 				= eosi.mu1_tab ;  
    mu2_tab 				= eosi.mu2_tab ;
    delta_car_tab 		= eosi.delta_car_tab ; 
    press_tab 				= eosi.press_tab ;
    n1_tab					= eosi.n1_tab;
    n2_tab 					= eosi.n2_tab ;
    d2psdmu1dmu2_tab   	= eosi.d2psdmu1dmu2_tab   ;
    dpsddelta_car_tab 	= eosi.dpsddelta_car_tab ;
    dn1sddelta_car_tab 	= eosi.dn1sddelta_car_tab ;
    dn2sddelta_car_tab	= eosi.dn2sddelta_car_tab;
    delta_car_n0			= eosi.delta_car_n0 ;
    mu1_n0					= eosi.mu1_n0 ;
    mu2_n0					= eosi.mu2_n0 ;
    delta_car_p0			= eosi.delta_car_p0 ;
    mu1_p0 					= eosi.mu1_p0;
    mu2_p0 					= eosi.mu2_p0 ;
    mu1_N  					= eosi.mu1_N  ;
    n_n_N  					= eosi.n_n_N;
    press_N  				= eosi.press_N;
    mu2_P  					= eosi.mu2_P ;
    n_p_P  					= eosi.n_p_P;
    press_P  				= eosi.press_P;
    
}

         //------------//
			//  Outputs   //
			//------------//

void Eos_bf_tabul::sauve(FILE* fich) const {

  Eos_bifluid::sauve(fich) ;
  
  char tmp_string[160] ;
  strcpy(tmp_string, tablename.c_str()) ;
  fwrite(tmp_string, sizeof(char), 160, fich) ;		

}


ostream& Eos_bf_tabul::operator>>(ostream & ost) const {

  ost <<
  "EOS of class Eos_bf_tabul : tabulated EOS depending on three variables (relative velocity and chemical potentials of neutrons and protons)."
  		<< '\n' ;
  ost << "Table read from " << tablename << endl ;
    	
  return ost ;
}

			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_bf_tabul::operator==(const Eos_bifluid& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
			cout << "The second EOS is not of type Eos_bf_tabul !" << endl ;
			resu = false ;
    }
    else { 
      	const Eos_bf_tabul& eos = dynamic_cast<const Eos_bf_tabul&>( eos_i ) ; 
    
      	if (eos.tablename != tablename){
				cout << 
	  			"The two Eos_bf_tabul have different names and not refer to the same tables! " << endl ; 
				resu = false ; 
      	}
   }
   
  return resu ; 
  
}


bool Eos_bf_tabul::operator!=(const Eos_bifluid& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------------------//
			//  Reading of the tables //
			//------------------------//
			

void Eos_bf_tabul::read_table() {

    using namespace Unites ;
    double m_b_si = rhonuc_si / 1e44; 	

    //------------------------------------
    //---------- 2 fluid table -----------
    //------------------------------------

    ifstream fich(tablename.data()) ;

    if (!fich) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Problem in opening the EOS file!" << endl ;
	  		cerr << "While trying to open " << tablename << endl ; 
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }

    fich.ignore(1000, '\n') ;             // Jump over the first header
    fich.ignore(2) ;
    getline(fich, authors, '\n') ;

    for (int i=0; i<5; i++) {					//  Jump over the file header
	  		fich.ignore(1000, '\n') ;     	// 	
    }                                   	//

    int nbp ;
    fich >> nbp ; fich.ignore(1000, '\n') ;   // number of data

    int n_delta, n_mu1, n_mu2; 	// number of values in delta_car, mu_n and  mu_p
    fich >> n_delta;  fich.ignore(1000, '\n') ;
    fich >> n_mu1;  	 fich.ignore(1000, '\n') ;
    fich >> n_mu2;    fich.ignore(1000, '\n') ;

    if (nbp<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of lines!" << endl ;
	  		cerr << "nbp = " << nbp << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    if (n_delta<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of values of delta_car!" << endl ;
	  		cerr << "n_delta = " << n_delta << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    if (n_mu1<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of values of mu_n!" << endl ;
	  		cerr << "n_mu1 = " << n_mu1 << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    if (n_mu2<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of values of mu_p!" << endl ;
	  		cerr << "n_mu2 = " << n_mu2 << endl ;
	  		cerr << "Aborting..." << endl ;
	 	 	abort() ;
    }
    if (n_mu2*n_mu1*n_delta != nbp ) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of lines!" << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }

    for (int i=0; i<3; i++) {				 // jump over the table
	  		fich.ignore(1000, '\n') ;      // description
    }                                      
	
    mu1_tab 				= new Tbl(n_delta, n_mu1, n_mu2) ;
    mu2_tab 				= new Tbl(n_delta, n_mu1, n_mu2) ;
    delta_car_tab 		= new Tbl(n_delta, n_mu1, n_mu2) ; 
    press_tab 				= new Tbl(n_delta, n_mu1, n_mu2) ;
    n1_tab					= new Tbl(n_delta, n_mu1, n_mu2) ;
    n2_tab 					= new Tbl(n_delta, n_mu1, n_mu2) ;
    d2psdmu1dmu2_tab   	= new Tbl(n_delta, n_mu1, n_mu2) ; 
    dpsddelta_car_tab 	= new Tbl(n_delta, n_mu1, n_mu2) ;
    dn1sddelta_car_tab 	= new Tbl(n_delta, n_mu1, n_mu2) ;
    dn2sddelta_car_tab 	= new Tbl(n_delta, n_mu1, n_mu2) ;

    mu1_tab->set_etat_qcq() ;
    mu2_tab->set_etat_qcq() ;
    delta_car_tab->set_etat_qcq() ; 
    press_tab->set_etat_qcq() ;
    n1_tab->set_etat_qcq() ;
    n2_tab->set_etat_qcq() ;
    d2psdmu1dmu2_tab  ->set_etat_qcq() ; 
    dpsddelta_car_tab->set_etat_qcq() ;
    dn1sddelta_car_tab->set_etat_qcq() ;
    dn2sddelta_car_tab->set_etat_qcq() ;

    //------------------------------------------------------
    // We have to choose the right unites (SI, cgs , LORENE)
    //------------------------------------------------------

    // Quantities given by the tabulated EoS (be careful to the units!)
    double mu1_MeV, mu2_MeV, n1_fm3, n2_fm3; 
    double Knp_Mev_2, press_MeV_fm3;
    double d2press_dmu1dmu2_MeV_fm3, dn1_ddelta_car_fm3, dn2_ddelta_car_fm3;

    // Stored quantities (in Lorene units)
    double delta_car_adim, mu1_adim, mu2_adim, n1_01fm3, n2_01fm3, Knp_adim  ;
    double press_adim,  dpress_ddelta_car_adim, dn1_ddelta_car_adim, dn2_ddelta_car_adim ;
    double d2press_dmu1dmu2_adim;
    
    double hbarc = 197.33 ; // Mev*fm
    double hbarc3 = hbarc * hbarc * hbarc ;

    // reading of the table.
    for (int j=0 ; j < n_delta ; j++) {
	  		for (int k=0 ; k < n_mu1 ; k++) {
					for (int l=0 ; l < n_mu2 ; l++) {
	      
		  				fich >> delta_car_adim ;
		  				fich >> mu1_MeV ;
		  				mu1_adim = mu1_MeV * mev_si / (m_b_si * v_unit * v_unit ) ;
		  				fich >> mu2_MeV ;
		  				mu2_adim = mu2_MeV * mev_si / (m_b_si * v_unit * v_unit ) ;
		  				fich >> n1_fm3 ;
		 				n1_01fm3 = 10. * n1_fm3 ;
		  				fich >> n2_fm3 ;
		  				n2_01fm3 = 10. * n2_fm3 ;
		  				fich >> Knp_Mev_2 ; 
		  				Knp_adim = Knp_Mev_2 / ( m_b_si * v_unit * v_unit *10. ) * (mev_si *hbarc3 )  ;
		 				dpress_ddelta_car_adim = - Knp_adim * n1_01fm3 * n2_01fm3  * pow( 1.-delta_car_adim, -1.5)  / 2. ;
		  				fich >> press_MeV_fm3 ;
		  				press_adim = press_MeV_fm3 * mevpfm3 ;	
		  				fich >>  d2press_dmu1dmu2_MeV_fm3 ;
		  				d2press_dmu1dmu2_adim = d2press_dmu1dmu2_MeV_fm3 * (10. * m_b_si * v_unit * v_unit ) / mev_si ;
		  				fich >> dn1_ddelta_car_fm3 ;
		  				dn1_ddelta_car_adim = dn1_ddelta_car_fm3 * 10. ;
		  				fich >> dn2_ddelta_car_fm3 ;
		  				dn2_ddelta_car_adim = dn2_ddelta_car_fm3 * 10. ;
	      
		  				fich.ignore(1000,'\n') ;    		
	      
		  				mu1_tab->set(j,k,l) = mu1_adim ;
		  				mu2_tab->set(j,k,l) = mu2_adim ;
		  				delta_car_tab->set(j,k,l) = delta_car_adim ;
		  				if ((n1_01fm3 <=0) && (n2_01fm3 <=0)) {press_adim = 0. ;}
		  				press_tab->set(j,k,l) = press_adim ;  	   
		  				n1_tab->set(j,k,l) = n1_01fm3 ;  
		  				n2_tab->set(j,k,l) = n2_01fm3 ; 
 			    	   if ((n1_01fm3 <= 0 ) || (n2_01fm3 <=0)) {
							d2press_dmu1dmu2_adim  	= 0. ;
							dpress_ddelta_car_adim	= 0. ;
							dn1_ddelta_car_adim 		= 0. ;
							dn2_ddelta_car_adim  	= 0. ;
		  				}
		  				d2psdmu1dmu2_tab ->set(j,k,l)  = d2press_dmu1dmu2_adim  ; 
		  				dpsddelta_car_tab->set(j,k,l)  = dpress_ddelta_car_adim  ;
		  				dn1sddelta_car_tab->set(j,k,l) = dn1_ddelta_car_adim  ;
		  				dn2sddelta_car_tab->set(j,k,l) = dn2_ddelta_car_adim  ;
					}  

					fich.ignore(1000, '\n') ; 

	  		}

	  		fich.ignore(1000, '\n') ;
	  
    }

    delta_car_min = (*delta_car_tab)(0, 0, 0)  ;
    delta_car_max = (*delta_car_tab)(n_delta-1, 0, 0) ;
    mu1_min = (*mu1_tab)(0, 0, 0)  ;  
    mu1_max = (*mu1_tab)(0, n_mu1-1, 0) ;
    mu2_min = (*mu2_tab)(0, 0, 0) ;
    mu2_max = (*mu2_tab)(0, 0, n_mu2-1);
	
    fich.close();

    //---------------------------------------------------------
    //---------- Table with a single fluid: fluid N -----------
    //---------------------------------------------------------

	 string tablename_1f_n = tablename.c_str() ;
    tablename_1f_n.replace(tablename_1f_n.end()-5,tablename_1f_n.end(),"_1f_n.d");

    ifstream fichN  ;
    fichN.open(tablename_1f_n.c_str()) ;
 
    if (!fichN) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Problem in opening the EOS file!" << endl ;
	  		cerr << "While trying to open " << tablename << endl ; 
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }

    fichN.ignore(1000, '\n') ;             // Jump over the first header	
    fichN.ignore(2) ;
    fichN.ignore(1000, '\n') ;  

    for (int i=0; i<5; i++) {	
	  		fichN.ignore(1000, '\n') ;    
    }                               

    int nbp_N ;		// number of data
    fichN >> nbp_N ; 
	 fichN.ignore(1000, '\n') ;   	
    int n_mu1_N; 		// number of different mu_n
    fichN >> n_mu1_N ;
	 fichN.ignore(1000, '\n') ;
	
    if (nbp_N<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of lines!" << endl ;
	  		cerr << "nbp = " << nbp << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    if (n_mu1_N<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of values of mu_n!" << endl ;
	  		cerr << "n_mu1 = " << n_mu1 << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    for (int i=0; i<3; i++) {					//  jump over the table
	  		fichN.ignore(1000, '\n') ;       // description
    }                                      
	
    mu1_N = new Tbl(n_mu1_N) ;
    n_n_N = new Tbl(n_mu1_N) ;
    press_N = new Tbl(n_mu1_N) ;

    mu1_N ->set_etat_qcq() ;
    n_n_N->set_etat_qcq() ;
    press_N ->set_etat_qcq() ;

    // Quantities given by the tabulated EoS (be careful to the units!)
    double mu1_MeV_N, n1_fm3_N, press_MeV_fm3_N;

    // Stored quantities (in Lorene units)
    double mu1_adimN, n1_01fm3N, press_adimN;

    for (int k=0 ; k < n_mu1_N ; k++) {
   	
   		fichN >> mu1_MeV_N ;
			mu1_adimN = mu1_MeV_N * mev_si /  (m_b_si * v_unit * v_unit ) ;
			fichN >> n1_fm3_N ; 
			n1_01fm3N = 10. * n1_fm3_N ;
    		fichN >> press_MeV_fm3_N;
			press_adimN = press_MeV_fm3_N * mevpfm3 ;
 			fichN.ignore(1000,'\n') ;  

			if ( (n1_01fm3N<0)  || (press_adimN < 0)){ 
		  		cout << "Eos_tabul::read_table(): " << endl ;
		  		cout << "Negative value in table!" << endl ;
		 		cout << "n_neutrons = " << n1_01fm3N <<  ", p = " << press_adimN << ", "<< endl ;
		 		cout << "Aborting..." << endl ;
		  		abort() ;		
			}

			mu1_N ->set(k) = mu1_adimN ;
			n_n_N->set(k) =  n1_01fm3N ; 
			press_N ->set(k) = press_adimN ;
   }

   fichN.close();

   //---------------------------------------------------------
   //---------- Table with a single fluid: fluid P -----------
   //---------------------------------------------------------
    
	 string tablename_1f_p = tablename.c_str() ;
    tablename_1f_p.replace(tablename_1f_p.end()-5,tablename_1f_p.end(),"_1f_p.d");

    ifstream fichP ;
    fichP.open(tablename_1f_p.c_str()) ;

    if (!fichP) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Problem in opening the EOS file!" << endl ;
	  		cerr << "While trying to open " << tablename << endl ; 
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }

    fichP.ignore(1000, '\n') ;             // Jump over the first header
    fichP.ignore(2) ;
    fichP.ignore(1000, '\n') ;  

    for (int i=0; i<5; i++) {			 			//  jump over the file
	  		fichP.ignore(1000, '\n') ;          // header
    }                                       	

    int nbp_P ;				// number of data
    fichP >> nbp_P ; 
	 fichP.ignore(1000, '\n') ;   	
    int n_mu2_P; 				// number of different values in mu_2
    fichP >> n_mu2_P ;
    fichP.ignore(1000, '\n') ;
    
    if (nbp_P<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of lines!" << endl ;
	  		cerr << "nbp = " << nbp << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    if (n_mu2_P<=0) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Wrong value for the number of values of mu_p!" << endl ;
	  		cerr << "n_mu2 = " << n_mu2 << endl ;
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }

    for (int i=0; i<3; i++) {				//  jump over the table
	  		fichP.ignore(1000, '\n') ;    // description
    }                                      
	
    mu2_P = new Tbl(n_mu2_P) ;
    n_p_P = new Tbl(n_mu2_P) ;
    press_P = new Tbl(n_mu2_P) ;

    mu2_P ->set_etat_qcq() ;
    n_p_P->set_etat_qcq() ;
    press_P ->set_etat_qcq() ;
	  
    // Quantities given by the tabulated EoS (be careful to the units!)
    double mu2_MeV_P, n2_fm3_P, press_MeV_fm3_P;
    	
    // Stored quantities (in Lorene units)
    double mu2_adimP,  n2_01fm3P, press_adimP;

    for (int l=0 ; l < n_mu2_P ; l++) {
  	
			fichP >> mu2_MeV_P ;
			mu2_adimP = mu2_MeV_P * mev_si /  (m_b_si * v_unit * v_unit ) ;
    		fichP >> n2_fm3_P ; 
			n2_01fm3P = 10. * n2_fm3_P ;
    		fichP >> press_MeV_fm3_P;
			press_adimP = press_MeV_fm3_P * mevpfm3 ;
    	   fichP.ignore(1000,'\n') ;    		

			if ( (n2_01fm3P<0) || (press_adimP < 0)){ 
		  		cout << "Eos_tabul::read_table(): " << endl ;
		  		cout << "Pegative value in table!" << endl ;
		  		cout <<", n_protons " << n2_01fm3P << ", p = " << press_adimP << ", "<< endl ;
		  		cout << "Aborting..." << endl ;
		  		abort() ;		
			}

			mu2_P ->set(l) = mu2_adimP ;
			n_p_P->set(l) =  n2_01fm3P ;
			press_P ->set(l) = press_adimP ;

    }

    fichP.close();

    //----------------------------------------------------
    //---------- Curve corresponding to np = 0 -----------
    //----------------------------------------------------
    // Rk: the table is sorted with increasing values of mu_n
    
	 string tablename_np_0 = tablename.c_str() ;
    tablename_np_0.replace(tablename_np_0.end()-5,tablename_np_0.end(),"_np=0.d");

 	 ifstream fich1 ;
    fich1.open(tablename_np_0.c_str()) ;

    if (!fich1) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Problem in opening the EOS file!" << endl ;
	  		cerr << "While trying to open " << tablename << endl ; 
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    
    int n_delta_n0, n_mu1_n0;
    fich1 >> n_delta_n0;fich1.ignore(1000, '\n') ;
    fich1 >> n_mu1_n0;fich1.ignore(1000, '\n') ;
    fich1.ignore(1000, '\n') ;             // Jump over the first header

    delta_car_n0 = new Tbl(n_delta_n0, n_mu1_n0) ;
    mu1_n0 = new Tbl(n_delta_n0, n_mu1_n0) ;
    mu2_n0 = new Tbl(n_delta_n0, n_mu1_n0) ;
	
    delta_car_n0 ->set_etat_qcq() ;
    mu1_n0->set_etat_qcq() ;
    mu2_n0->set_etat_qcq() ;

    double delta_car_nn0, mu1_MeV_nn0, mu2_MeV_nn0; 
    	
    for (int o = 0; o < n_delta_n0 ; o++ ) { 
  
      for (int p = 0 ; p < n_mu1_n0 ; p++) {
   	
	    		fich1 >> delta_car_nn0 ; 
	    		fich1 >> mu1_MeV_nn0 ;
	    		fich1 >> mu2_MeV_nn0 ;
	    		fich1.ignore(1000,'\n') ;    		
	
	    		delta_car_n0 ->set(o,p) = delta_car_nn0;
	    		mu1_n0 ->set(o,p) = mu1_MeV_nn0 * mev_si / (m_b_si * v_unit * v_unit ) ;
	    		mu2_n0 ->set(o,p) = mu2_MeV_nn0 * mev_si / (m_b_si * v_unit * v_unit ) ;
		  
      }
      fich1.ignore(1000,'\n') ; 
	     
    }

    fich1.close();
 
    //----------------------------------------------------
    //---------- Curve corresponding to nn = 0 -----------
    //----------------------------------------------------
    // Rk: the table is sorted with increasing values of mu_p

	 string tablename_nn_0 = tablename.c_str() ;
    tablename_nn_0.replace(tablename_nn_0.end()-5,tablename_nn_0.end(),"_nn=0.d");

 	 ifstream fich2 ;
    fich2.open(tablename_nn_0.c_str()) ;

    if (!fich2) {
	  		cerr << "Eos_bf_tabul::read_table(): " << endl ;
	  		cerr << "Problem in opening the EOS file!" << endl ;
	  		cerr << "While trying to open " << tablename << endl ; 
	  		cerr << "Aborting..." << endl ;
	  		abort() ;
    }
    
    int n_delta_p0, n_mu2_p0;
    fich2 >> n_delta_p0;fich2.ignore(1000, '\n') ;
    fich2 >> n_mu2_p0;fich2.ignore(1000, '\n') ;
    fich2.ignore(1000, '\n') ;            
    
    delta_car_p0 = new Tbl(n_delta_p0, n_mu2_p0) ;
    mu1_p0 = new Tbl(n_delta_p0, n_mu2_p0) ;
    mu2_p0 = new Tbl(n_delta_p0, n_mu2_p0) ;

    delta_car_p0 ->set_etat_qcq() ;
    mu1_p0->set_etat_qcq() ;
    mu2_p0 ->set_etat_qcq() ;
	
    double delta_car_np0, mu1_MeV_np0, mu2_MeV_np0; 
    	
    for (int o = 0; o < n_delta_p0 ; o++ ) { 
  
      for (int p = 0 ; p < n_mu2_p0 ; p++) {
   	
	    		fich2 >> delta_car_np0 ; 
	    		fich2 >> mu1_MeV_np0 ;
	    		fich2 >> mu2_MeV_np0 ;
	    		fich2.ignore(1000,'\n') ;    		

	   		delta_car_p0 ->set(o,p) = delta_car_np0;
	    		mu1_p0 ->set(o,p) = mu1_MeV_np0 * mev_si / (m_b_si * v_unit * v_unit ) ;
	    		mu2_p0 ->set(o,p) = mu2_MeV_np0 * mev_si / (m_b_si * v_unit * v_unit ) ;
		
      }
      fich2.ignore(1000,'\n') ; 
  
    }

    fich2.close();

}


		//------------------------------//
		//    Computational routines    //
		//------------------------------//


// Complete computational routine giving all thermo variables
//-----------------------------------------------------------

void Eos_bf_tabul::calcule_interpol(const Cmp& ent1, const Cmp& ent2, const Cmp& delta2, 
				    							Cmp& nbar1, Cmp& nbar2, Cmp& ener, Cmp& press, 
				    							Cmp& K_nn, Cmp& K_np, Cmp& K_pp, Cmp& alpha_eos,
				    							int nzet, int l_min) const {

    	assert(ent1.get_etat() != ETATNONDEF) ; 
    	assert(ent2.get_etat() != ETATNONDEF) ; 
   	assert(delta2.get_etat() != ETATNONDEF) ;
  
    	const Map* mp = ent1.get_mp() ;	// Mapping
    	assert(mp == ent2.get_mp()) ;
   	assert(mp == delta2.get_mp()) ;
    	assert(mp == ener.get_mp()) ;
  
    	if ((ent1.get_etat() == ETATZERO)&&(ent2.get_etat() == ETATZERO)) {
	
			nbar1.set_etat_zero() ; 
			nbar2.set_etat_zero() ; 
			ener.set_etat_zero() ; 
			press.set_etat_zero() ; 
			K_nn.set_etat_zero() ; 
			K_np.set_etat_zero() ; 
			K_pp.set_etat_zero() ; 
			alpha_eos.set_etat_zero() ; 
	
			return ; 
 
   	}
    
    	nbar1.allocate_all() ;
    	nbar2.allocate_all() ;
   	ener.allocate_all() ;
    	press.allocate_all() ;
    	K_nn.allocate_all() ;
    	K_np.allocate_all() ;
    	K_pp.allocate_all() ;
    	alpha_eos.allocate_all() ;
  
    	const Mg3d* mg = mp->get_mg() ;	// Multi-grid
  
    	int nz = mg->get_nzone() ;		// total number of domains
  
  		// resu is set to zero in the other domains :
  		if (l_min > 0) {
  
    		nbar1.annule(0, l_min-1) ; 
    		nbar2.annule(0, l_min-1) ; 
    		ener.annule(0, l_min-1) ; 
    		press.annule(0, l_min-1) ; 
   		K_nn.annule(0, l_min-1) ; 
    		K_np.annule(0, l_min-1) ; 
    		K_pp.annule(0, l_min-1) ; 
    		alpha_eos.annule(0, l_min-1) ; 
 
  		}
  
  		if (l_min + nzet < nz) {
 
    		nbar1.annule(l_min + nzet, nz - 1) ; 
    		nbar2.annule(l_min + nzet, nz - 1) ; 
    		ener.annule(l_min + nzet, nz - 1) ; 
   		press.annule(l_min + nzet, nz - 1) ; 
    		K_nn.annule(l_min + nzet, nz - 1) ; 
    		K_np.annule(l_min + nzet, nz - 1) ; 
    		K_pp.annule(l_min + nzet, nz - 1) ; 
    		alpha_eos.annule(l_min + nzet, nz - 1) ; 

  		}

  		int np0 = mp->get_mg()->get_np(0) ;
  		int nt0 = mp->get_mg()->get_nt(0) ;
  
  		for (int l=l_min; l<l_min+nzet; l++) {
      		assert(mp->get_mg()->get_np(l) == np0) ;
      		assert(mp->get_mg()->get_nt(l) == nt0) ; 
  		}

  		for (int k=0; k<np0; k++) {
    		for (int j=0; j<nt0; j++) {
				for (int l=l_min; l< l_min + nzet; l++) {
	      		for (int i=0; i<mp->get_mg()->get_nr(l); i++) {

	  					double xx, xx1, xx2; // xx1 = H1 = ln(mu1/m1)
						xx1 = ent1(l,k,j,i) ; 
						xx2 = ent2(l,k,j,i) ; 
						xx = delta2(l,k,j,i) ;

						if (xx < delta_car_min) {
		     					cout << "Eos_bf_tabul::calcule_tout : delta2 < delta_car_min !" << endl ;
		      				abort() ; 
						}
						if (xx > delta_car_max) {	
		      				cout << "Eos_bf_tabul::calcule_tout : delta2 > delta_car_max !" << endl ;
		      				abort() ; 
						}		
						if (m_1 * exp(xx1) > mu1_max) {
		      				cout << "Eos_bf_tabul::calcule_tout : ent1 > mu1_max !" << endl ;
		      				abort() ; 
						}
						if (m_2 *exp(xx2) > mu2_max) {
		      				cout << "Eos_bf_tabul::calcule_tout : ent2 > mu2_max !" << endl ;
		      				abort() ; 
						}
	
						double n1 = 0 , n2 = 0, pressure = 0 ;
						double alpha_int = 0, K11 = 0, K12 = 0, K22 = 0 ;
	  	
						double mu1 = m_1 * exp(xx1);
						double mu2 = m_2 * exp(xx2);  
	    
						if ( (exp(xx1) < 1.) && (exp(xx2) < 1.) ) {	// no fluds 
		     					n1 = 0. ;
		     					n2 = 0. ;
		      				pressure = 0.;
		      				alpha_int = 0 ;
		      				K11  = 0 ;
		      				K12 =  0 ;
		      				K22 = 0 ;
						}
	    
						else { // at least one fluid is present !
    
		      				double p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha_interpol ;
	
		     					Eos_bf_tabul::interpol_3d_bifluid(xx, mu1, mu2, p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha_interpol) ;

		      				n1 = dpsdent1_interpol ;
		      				n2 = dpsdent2_interpol ;
		      				pressure = p_interpol;
		      				alpha_int = alpha_interpol ; 	    
						}
  
	  					if (n1 < 0 ) {
		    				cout << "Eos_bf_tabul::calcule_tout : nbar1<0 !" << endl ;
		    				// abort() ; 
		    				n1 = 0 ;
						}
						if (n2 < 0 ) {
		    				cout << "Eos_bf_tabul::calcule_tout : nbar2<0 !" << endl ;
		    				// abort() ; 
		    				n2 = 0 ;		
						}
						if (pressure < 0 ) {
		    				cout << "Eos_bf_tabul::calcule_tout : pressure < 0 !" << endl ;
		    				// abort() ; 
		    				pressure = 0 ;
						}
	
						// Knn
						if (n1 >0.) { 
		      			K11 = mu1 / n1 - double(2) * alpha_int * ( 1. - xx) / ( n1 * n1 ) ; 
						}	
						// Kpp
						if (n2 >0.) { 
		     	 			K22 = mu2 / n2 - double(2) * alpha_int * ( 1. - xx) / ( n2 * n2 ) ; 
						}
						// Knp
						if ((n1 <= 0.) || (n2 <= 0.) ) { 
		      			K12 = 0. ;
		      			alpha_int = 0. ;
						}
						else {
		      			K12 = double(2) * alpha_int * pow(1.-xx, 1.5)/ ( n1 * n2 );
						}

						nbar1.set(l,k,j,i) = n1 ;
						nbar2.set(l,k,j,i) = n2 ;
						press.set(l,k,j,i) = pressure ;
						ener.set(l,k,j,i) = mu1 * n1 + mu2 * n2 - pressure ;
						K_nn.set(l,k,j,i) = K11 ;   
						K_np.set(l,k,j,i) = K12;  
						K_pp.set(l,k,j,i) = K22 ;   
						alpha_eos.set(l,k,j,i) = alpha_int ;   

	  				}
				} 
      	}
    	}
 
}


// Baryon densities from enthalpies 
//---------------------------------

// this function is not called anymore but should be implemented (virtual function)
bool Eos_bf_tabul::nbar_ent_p(const double ent1, const double ent2, 
				const double delta2, double& nbar1, 
			       double& nbar2) const { 

	bool one_fluid = false; 

	if (delta2 < delta_car_min) {
				cout << "Eos_bf_tabul::nbar_ent_p : delta2 < delta_car_min !" << endl ;
           	abort() ; 
	}
 	if (delta2 > delta_car_max) {	
				cout << "Eos_bf_tabul::nbar_ent_p : delta2 > delta_car_max !" << endl ;
           	abort() ; 
	}
 	if (m_1 * exp(ent1) > mu1_max) {
				cout << "Eos_bf_tabul::nbar_ent_p : ent1 > mu1_max !" << endl ;
           	abort() ; 
	}
 	if (m_2 *exp(ent2) > mu2_max) {
				cout << "Eos_bf_tabul::nbar_ent_p : ent2 > mu2_max !" << endl ;
           	abort() ; 
	}
	
	if ( (exp(ent1) < 1.) && (exp(ent2) < 1.) ) {	
			nbar1 = 0. ;
			nbar2 = 0. ;
	}
	else {

	  		double mu1 = m_1 * exp(ent1);
	  		double mu2 = m_2 * exp(ent2);

	  		double p_interpol ;
	  		double dpsdent1_interpol ; 
	  		double dpsdent2_interpol ;
	  		double alpha ;
	
	  		Eos_bf_tabul::interpol_3d_bifluid(delta2, mu1, mu2, 
			p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;

	  		nbar1 = dpsdent1_interpol ;
	  		nbar2 = dpsdent2_interpol ;
	
	}
  
   if (nbar1 < 0 ) {
		cout << "Eos_bf_tabul::nbar_ent_p : nbar1<0 !" << endl ;
      //	abort() ; 
		nbar1 = 0 ;
	}
  	if (nbar2 < 0 ) {
		cout << "Eos_bf_tabul::nbar_ent_p : nbar2<0 !" << endl ;
     	// abort() ; 
		nbar2 = 0 ;		
	}

    one_fluid = ((nbar1 <= 0.)||(nbar2 <= 0.)) ;
    
    return one_fluid; 

}
   
// One fluid sub-EOSs
//-------------------

// this function is not called anymore but should be implemented (virtual function)
double Eos_bf_tabul::nbar_ent_p1(const double ent1) const {

	c_est_pas_fait("Eos_bf_tabul::nbar_ent_p1" ) ;

   return ent1 ;

	/*
	double pressN_interpol ; 
	double n_n_N_interpol ; 
	double mu1 = m_1 * exp(ent1);
	int i =0;
	
	if (exp(ent1) < 1. ) {
			n_n_N_interpol = 0. ;
	}
	else {
			interpol_herm( *mu1_N, *press_N, *n_n_N,
		   					mu1,  i, pressN_interpol , n_n_N_interpol) ;
	}
	return n_n_N_interpol ;
	*/

}

// this function is not called anymore but should be implemented (virtual function)
 double Eos_bf_tabul::nbar_ent_p2(const double ent2) const {

	c_est_pas_fait("Eos_bf_tabul::nbar_ent_p2" ) ;

   return ent2 ;

	/*
	double pressP_interpol ; 
	double n_p_P_interpol ; 
	double mu2 = m_2 * exp(ent2);
	int i =0;
	if (exp(ent2) < 1. ) {
			n_p_P_interpol = 0. ;
	}
	else {
			interpol_herm( *mu2_P, *press_P, *n_p_P,
		   					mu2,  i, pressP_interpol , n_p_P_interpol) ;
	}
	return n_p_P_interpol ;
	*/
}

 // Energy density from baryonic densities 
//-----------------------------------------

// this function is not called anymore but should be implemented (virtual function)
double Eos_bf_tabul::ener_nbar_p(const double nbar1, const double nbar2, 
		   const double delta2) const{

  c_est_pas_fait("Eos_bf_tabul::ener_nbar_p" ) ;

  return nbar1 + nbar2 + delta2;

 }

// Pressure from baryonic densities 
//----------------------------------

// this function is not called anymore but should be implemented (virtual function)
double Eos_bf_tabul::press_nbar_p(const double nbar1, const double nbar2, 
		    const double delta2) const {

  c_est_pas_fait("Eos_bf_tabul::press_nbar_p" ) ;

    return nbar1 + nbar2 + delta2;

} 
     
	
// Pressure from baryonic log-enthalpies
//--------------------------------------

// this function is not called  but can be useful if necessary
double Eos_bf_tabul::press_ent_p(const double ent1, const double ent2, const double delta2) const {
	
	if (delta2 < delta_car_min) {
				cout << "Eos_bf_tabul::press_ent_p : delta2 < delta_car_min !" << endl ;
           	abort() ; 
	}
 	if (delta2 > delta_car_max) {
				cout << "Eos_bf_tabul::press_ent_p : delta2 > delta_car_max !" << endl ;
           	abort() ; 
	}
 	if (m_1 * exp(ent1) > mu1_max) {
				cout << "Eos_bf_tabul::press_ent_p : ent1 > mu1_max !" << endl ;
           	abort() ; 
	}
 	if (m_2 * exp(ent2) > mu2_max) {
				cout << "Eos_bf_tabul::press_ent_p : ent2 > mu2_max !" << endl ;
           	abort() ; 
	}

	double pressure ;

	if ( (exp(ent1) < 1.) && (exp(ent2) < 1.)) {
			//abort();
			pressure = 0. ;
	}
	else {
	
	  		double mu1 = m_1 * exp(ent1);
	  		double mu2 = m_2 * exp(ent2);

	  		double p_interpol ;
	  		double dpsdent1_interpol ; 
	  		double dpsdent2_interpol ;
	  		double alpha ;
	
	  		Eos_bf_tabul::interpol_3d_bifluid(delta2, mu1, mu2, 
				p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;

 			pressure = p_interpol;

	}
 
	if (pressure < 0 ) {
			cout << "Eos_bf_tabul::press_ent_p : pressure < 0 !" << endl ;
         // abort() ; 
			pressure = 0 ;
	}
	return pressure ;
 }


// Energy from baryonic log - enthalpies
//--------------------------------------

// this function is not called but can be useful if necessary
  double Eos_bf_tabul::ener_ent_p(const double ent1, const double ent2, 
				  const double delta2) const {
   double energy= 0.; 

   if ( (exp(ent1) < 1.) && ( exp(ent2) < 1.)) {
     	energy = 0. ;
   }
   else {
     	double mu1 = m_1 * exp(ent1);
		double mu2 = m_2 * exp(ent2);

   	double p_interpol ;
		double dpsdent1_interpol ; 
   	double dpsdent2_interpol ;
		double alpha ;
	
		Eos_bf_tabul::interpol_3d_bifluid( delta2, mu1, mu2, p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;

		energy = mu1 * dpsdent1_interpol + mu2 * dpsdent2_interpol - p_interpol ;
	}
	
	if (energy < 0 ) {
		cout << "Eos_bf_tabul::ener_ent_p : energy < 0 !" << endl ;
      //    	abort() ; 
	energy = 0 ;
	}
	return energy;

}


// Alpha  from baryonic log - enthalpies
//--------------------------------------- 

// this function is not called but can be useful if necessary
double Eos_bf_tabul::alpha_ent_p(const double ent1, const double ent2, 
				 const double delta2) const {

	if (delta2 < delta_car_min) {
		cout << "Eos_bf_tabul::alpha_ent_p : delta2 < delta_car_min !" 
		     << endl ;
      abort() ; 
	}
 	if (delta2 > delta_car_max) {
		cout << "Eos_bf_tabul::alpha_ent_p : delta2 > delta_car_max !" 
		     << endl ;
       abort() ; 
	}
	if (m_1 * exp(ent1) > mu1_max) {
		cout << "Eos_bf_tabul::alpha_ent_p : ent1 > mu1_max !" << endl ;
      abort() ; 
	}
	if (m_2 * exp(ent2) > mu2_max) {  
		cout << "Eos_bf_tabul::alpha_ent_p : ent2 > mu2_max !" << endl ;
      abort() ; 
	}

	double alpha;
	
	if ((exp(ent1) <= 1.) && (exp(ent2) <= 1.) ) {
		alpha = 0. ;
	}
	else {
	  	double mu1 = m_1 * exp(ent1);
	  	double mu2 = m_2 * exp(ent2);
	   double p_interpol ;
	  	double dpsdent1_interpol ; 
	  	double dpsdent2_interpol ;

		Eos_bf_tabul::interpol_3d_bifluid( delta2, mu1, mu2, p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;

	}
 
	return alpha;
}


// Derivatives of energy
//----------------------

// this function is not called but can be useful if necessary
double Eos_bf_tabul::get_K11(const double delta2, const double ent1, const double ent2)  const  {

	double xx = 0.; // K_nn
	double mu_1 = m_1 * exp(ent1);
	double nbar1;
	double nbar2;

	if ((exp(ent1) <= 1.) && (exp(ent2) <= 1.) ){
		xx = 0. ;
	}
	else {
		
		Eos_bf_tabul::nbar_ent_p(ent1,ent2, delta2, nbar1, nbar2) ; 
 
		double alpha = Eos_bf_tabul::alpha_ent_p(ent1,ent2,delta2) ;
		if (nbar1 >0.) {
				xx = mu_1 / nbar1 - double(2) * alpha * ( 1. - delta2) / ( nbar1 * nbar1 ) ;
		}
	}

   return xx;
}

// this function is not called but can be useful if necessary
double Eos_bf_tabul::get_K22( const double delta2, const double ent1, const double ent2)  const {
	
	double xx=0.;
	double mu_2 = m_2  * exp (ent2) ;
	double nbar1;
	double nbar2;
	
	if ((exp(ent1) <= 1.) && (exp(ent2) <= 1.) ){
		xx = 0. ;
	}
	else {

	  Eos_bf_tabul::nbar_ent_p(ent1,ent2, delta2, nbar1,nbar2) ;

	  double alpha = Eos_bf_tabul::alpha_ent_p(ent1,ent2,delta2) ;
	  if (nbar2 >0.) {
	    		xx = mu_2 / nbar2 - double(2) * alpha * ( 1. - delta2) / ( nbar2 * nbar2 ) ;
	  }
   }
  
   return xx;
}

double Eos_bf_tabul::get_K12(const double delta2, const double ent1, const double ent2)  const {
	double xx =0.;
	double nbar1;
	double nbar2;
	
	if ((exp(ent1) <= 1.) && (exp(ent2) <= 1.) ){
		xx = 0. ;
	}
	else {

	  	Eos_bf_tabul::nbar_ent_p(ent1,ent2, delta2, nbar1,nbar2) ;

	 	double alpha = Eos_bf_tabul::alpha_ent_p(ent1,ent2,delta2) ;
	  	if ((nbar1 <= 0.) || (nbar2 <= 0.) ) { 
	      xx = 0. ;
	  	}
	  	else {
		  	xx = double(2) * alpha * pow(1.-delta2, 1.5)/ ( nbar1 * nbar2 );
		}
	}

   return xx;
}

// Computes the interpolated values of the pressure, the baryon densities and alpha at the point under consideration from tabulated EoSs.
// This routine uses the following 3D interpolation scheme : Hermite interpolation in the chemical potentials
// and linear interpolation in the relative speed.
// --------------------------------------------------------------------------------------------------------------------------------------
void Eos_bf_tabul::interpol_3d_bifluid(const double delta2, const double mu1, const double mu2, double& press, double& nbar1, double& nbar2, double& alpha) const
{
	 
	assert((*mu1_tab).dim == (*delta_car_tab).dim) ; 
   assert((*mu2_tab).dim == (*delta_car_tab).dim) ;
   assert((*press_tab).dim == (*delta_car_tab).dim) ;
   assert((*n1_tab).dim == (*delta_car_tab).dim) ;
   assert((*n2_tab).dim == (*delta_car_tab).dim) ;
   assert((*d2psdmu1dmu2_tab  ).dim == (*delta_car_tab).dim) ;
  
   int nbp1, nbp2, nbp3;
   nbp1 = (*delta_car_tab).get_dim(2) ; // number of values of \Delta^2 in the table
   nbp2 = (*delta_car_tab).get_dim(1) ; // number of values of \mu_n in the table
   nbp3 = (*delta_car_tab).get_dim(0) ; // number of values of \mu_p in the table

	Tbl* null_tab = new  Tbl(nbp1,nbp2,nbp3) ; // Table whose components are all equal to zero
 	null_tab->set_etat_zero() ;
	
   int i_near = 0 ; 
   int j_near = 0 ;
   int k_near = 0 ;

   // looking for the positions of (delta2,mu1,mu2) in the tables
	while ( ( (*delta_car_tab)(i_near,0,0) <= delta2 ) && ( ( nbp1-1 ) > i_near ) ) {
		i_near++ ;
  	}
  	if (i_near != 0) { 
		i_near -- ; 
   }
	while ( ( (*mu1_tab)(i_near,j_near, 0) <= mu1 ) && ( ( nbp2-1 ) > j_near ) ) {
		j_near++ ;
   }
   if (j_near != 0) {
		j_near -- ; 
   }
	while ( ( (*mu2_tab)( i_near, j_near, k_near) <= mu2) && ( ( nbp3-1 ) > k_near ) ) {
		k_near++ ;
   }
   if (k_near != 0) {
		k_near-- ; 
   }
	int i1 = i_near + 1 ;
   int j1 = j_near + 1 ;
   int k1 = k_near + 1 ;

	// The location in the table is refined if necessary
   if ( ( (*delta_car_tab)( i_near, j_near, k_near) > delta2 ) && (i_near !=0 ) )   {
      i_near -= 1;
      i1 -= 1; 
   }
   if (  (delta2 > (*delta_car_tab)( i1, j_near, k_near) ) && (i_near != nbp1 ) )  {
      i_near += 1;
      i1 += 1; 
   }
   if (  ( (*mu1_tab)( i_near, j_near, k_near) > mu1 ) && (j_near !=0 ) )   {
      j_near -= 1;
      j1 -= 1; 
   }
   if (  ( mu1 > (*mu1_tab)( i1, j1, k_near) ) && ( j_near != nbp2) )   {
      j_near += 1;
      j1 += 1; 
   }
   if (  ( (*mu2_tab)( i_near, j_near, k_near) > mu2 ) && (k_near !=0 ) )   {
      k_near -= 1;
      k1 -= 1; 
   }
   if (  ( mu2 > (*mu2_tab)( i1, j_near, k1) ) && ( k_near != nbp3) )   {
      k_near += 1;
      k1 += 1; 
   } 
     
	// Check of the location
  	if ( ( (*delta_car_tab)( i_near, j_near, k_near) > delta2 ) || (delta2 > (*delta_car_tab)( i1, j_near, k_near) ) ) {
      cout << "bad location of delta2 in *delta_car_tab " << endl ;
      cout << (*delta_car_tab)( i_near, j_near, k_near) << "  " << delta2 <<  "  " << (*delta_car_tab)( i1, j_near, k_near) << endl;
      abort();
   }
   if ( ( (*mu1_tab)( i_near, j_near, k_near) > mu1 ) || (mu1 > (*mu1_tab)( i1, j1, k_near) ) ) {
      cout << "bad location of mu1 in *mu1_tab " << endl ;
      cout << (*mu1_tab)( i_near, j_near, k_near) << "  " << mu1 <<   "  " << (*mu1_tab)( i1, j1, k_near) << endl;
      abort();
   }
   if ( ( (*mu2_tab)( i_near, j_near, k_near) > mu2 ) || ( mu2 > (*mu2_tab)( i1, j_near, k1) ) ){
      cout << "bad location of mu2 in *mu2_tab "<< endl ;
      cout << (*mu2_tab)( i_near, j_near, k_near) << "  " << mu2 <<  "  " << (*mu2_tab)( i1, j_near, k1) << endl;
      abort();
   }
	
   // Values in the slice i 
	double press_i_near 	= 0. ;
  	double nbar1_i_near 	= 0. ;
 	double nbar2_i_near 	= 0. ;
  	double malpha_i_near = 0. ;
	// Values in the slice i+1  	
	double press_i1  = 0. ;
  	double nbar1_i1  = 0. ;
   double nbar2_i1  = 0. ;
   double malpha_i1 = 0. ; 
   // -alpha  
	double malpha = 0. ;
 
  	int n_deltaN = (*delta_car_n0).get_dim(1) ;
  	int n_mu1N 	= (*delta_car_n0).get_dim(0) ;
  	int n_deltaP = (*delta_car_p0).get_dim(1) ;
  	int n_mu2P 	= (*delta_car_p0).get_dim(0) ;


  /*******************************
 	* 2D interpolation on Slice i *
	*******************************/ 
 
  // Looking for the table to be used (concerning protons only, neutrons only or both fluids).
  // -----------------------------------------------------------------------------------------

  int Placei = 0 ; // 0 = two fluids, 1 = only neutrons (fluid 1), 2 = only protons (fluid 2), 3 = no fluids

  int i_nearN_i = 0;
  int j_nearN_i = 0;
  int i_nearP_i = 0;
  int j_nearP_i = 0;
  
  if ( mu1 > m_1 ) // both fluids are present or only the neutron fluid (fluid 1)
  {
		// to find if one or two fluid(s) is (are) present, we compare the position under consideration  
      // with the curve n_p = 0 for the EoS which is used.
		// Note that the following procedure is adapted to DDH and DDHdelta EoSs (close to beta eq. and corotation)
		// but the curve n_p = 0 could be possibly much more complicated depending on the EoS...
		while ( ( (*delta_car_n0)(i_nearN_i,0) <= (*delta_car_tab)(i_near, j_near, k_near) ) && ( ( n_deltaN-1 ) > i_nearN_i ) ) {
  	  		i_nearN_i++ ;
  		}
  		if (i_nearN_i != 0) { 
	  		i_nearN_i -- ; 
 		}
		while ( ( (*mu1_n0)(i_nearN_i,j_nearN_i) <= mu1 ) && ( ( n_mu1N-1 ) > j_nearN_i ) ) {
	  		j_nearN_i++ ;
  		}
  		if (j_nearN_i != 0) { 
	  		j_nearN_i -- ; 
  		}

 		// some checks
		if ( ( (*delta_car_n0)(i_nearN_i,0) > (*delta_car_tab)(i_near, j_near, k_near)  ) || ((*delta_car_tab)(i_near, j_near, k_near)  > (*delta_car_n0)(i_nearN_i+1,0) ) ) 
		{
 	  		cout << " bad location of delta_car_tab_i in *delta_car_n0 (courbe limite np = 0) " << endl ;
 	  		cout << (*delta_car_n0)(i_nearN_i,0) << "  " << (*delta_car_tab)(i_near, j_near, k_near) <<  "  " <<  (*delta_car_n0)(i_nearN_i+1,0) << endl;
 			abort();
 		}
	 	if ( ( (*mu1_n0)(i_nearN_i,j_nearN_i) > mu1  ) || (mu1 > (*mu1_n0)(i_nearN_i,j_nearN_i+1) ) ) 
		{
 	  		cout << " bad location of mu_n in  *mu1_n0 (limit curve np = 0) " << endl ;
 	  		cout << (*mu1_n0)(i_nearN_i,j_nearN_i) << "  " << mu1 <<  "  " <<  (*mu1_n0)(i_nearN_i,j_nearN_i+1) << endl;
 			abort();
 		}	
		
     	double aN_i, bN_i;
		aN_i = ((*mu2_n0)(i_nearN_i,j_nearN_i+1) -  (*mu2_n0)(i_nearN_i,j_nearN_i) ) / ((*mu1_n0)(i_nearN_i,j_nearN_i+1) -  (*mu1_n0)(i_nearN_i,j_nearN_i) ) ;
		bN_i = (*mu2_n0)(i_nearN_i,j_nearN_i)  - aN_i * (*mu1_n0)(i_nearN_i,j_nearN_i) ;
		double zN_i = aN_i * mu1 + bN_i ;

	  	if (zN_i <  mu2) 
		{ 
				Placei = 0; 		// two fluids
	  	}
	  	else 
		{	
				Placei = 1 ;		// fluid 1 only
	  	}	
  }
  
  else	// both fluids are present or only the proton fluid (fluid 2) or no fluids !
  {
    	if ( mu2 <= m_2) 
		{
				Placei = 3 ; 		// no fluids
		}
		else 
		{
		
			// to find if one or two fluid(s) is (are) present, we compare the position under consideration  
     		// with the curve n_n = 0 for the EoS which is used.
			// Note that the following procedure is adapted to DDH and DDHdelta EoSs (close to beta eq. and corotation)
			// but the curve n_n = 0 could be possibly much more complicated depending on the EoS...
			while ( ( (*delta_car_p0)(i_nearP_i,0) <= (*delta_car_tab)(i_near, j_near, k_near)) && ( ( n_deltaP-1 ) > i_nearP_i ) ) {
      			i_nearP_i++ ;
	  		}
	  		if (i_nearP_i != 0) { 
      			i_nearP_i -- ; 
	  		}
			while ( ( (*mu2_p0)(i_nearP_i,j_nearP_i) <= mu2 ) && ( ( n_mu2P-1 ) > j_nearP_i ) ) {
		 		j_nearP_i++ ;
	 		}
	 		if (j_nearP_i != 0) { 
					j_nearP_i -- ; 
	  		}

			// some checks
 	      if ( ( (*delta_car_p0)(i_nearP_i,0) > (*delta_car_tab)(i_near, j_near, k_near)  ) || ((*delta_car_tab)(i_near, j_near, k_near)  > (*delta_car_p0)(i_nearP_i+1,0) ) ) 
			{
 	  			cout << " bad location of delta_car_tab_i in *delta_car_p0 (courbe limite nn = 0) " << endl ;
	 	      cout << (*delta_car_p0)(i_nearP_i,0) << "  " << (*delta_car_tab)(i_near, j_near, k_near) <<  "  " <<  (*delta_car_p0)(i_nearP_i+1,0) << endl;
 				abort();
 			}
	 		if ( ( (*mu2_p0)(i_nearP_i,j_nearP_i) > mu2  ) || (mu2 > (*mu2_p0)(i_nearP_i,j_nearP_i+1) ) ) {
 	  			cout << " bad location of mu_p in  *mu2_p0 (limit curve nn = 0) " << endl ;
	 	      cout << (*mu2_p0)(i_nearP_i,j_nearP_i) << "  " << mu2 <<  "  " <<  (*mu2_p0)(i_nearP_i,j_nearP_i+1) << endl;
 				abort();
 			}	
		
			double aP_i, bP_i;
	  		aP_i = ( (*mu2_p0)(i_nearP_i,j_nearP_i+1) -  (*mu2_p0)(i_nearP_i,j_nearP_i) ) / ( (*mu1_p0)(i_nearP_i,j_nearP_i+1) -  (*mu1_p0)(i_nearP_i,j_nearP_i) ) ;
	  		bP_i = (*mu2_p0)(i_nearP_i,j_nearP_i)  - aP_i * (*mu1_p0)(i_nearP_i,j_nearP_i) ;
	  		double yP_i = (mu2- bP_i) /aP_i ;
		
	  		if (yP_i <  mu1) 
			{ 
	    			Placei = 0; 		// two fluids
	  		}
	  		else 
			{	
	    			Placei = 2 ;   	// fluid 2 only
	  		}
	  
		}

	}	// end of the search of the table to be used.

	/*	
	cout << mu2* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "   "  << mu1 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13) << endl ;
	cout << " Placei = " << Placei << endl ;
	*/

	// Interpolation in the different areas
	if (Placei == 3 ) // no fluids
	{ 			
  			press_i_near = 0. ;
			nbar1_i_near = 0. ;
			nbar2_i_near = 0. ;
			malpha_i_near = 0.;
  	}
 	else if (Placei == 1 )  // fluid 1 only
	{	
			malpha_i_near = 0.;
			nbar2_i_near = 0. ;
			int i = 0;
			interpol_herm(*mu1_N, *press_N, *n_n_N, mu1,  i, press_i_near , nbar1_i_near ) ;
 			if (press_i_near < 0.) { 
		  			cout << " INTERPOLATION FLUID N --> negative pressure " << endl ;
		  			abort();
			}
 			if (nbar1_i_near < 0.) { 
		  			cout << " INTERPOLATION FLUID N --> negative density " << endl ;
		  			abort();
			}
  	}
  	else if (Placei == 2 ) // fluid 2 only
	{	
			malpha_i_near = 0.;
			nbar1_i_near = 0. ;
			int i =0;
			interpol_herm( *mu2_P, *press_P, *n_p_P, mu2,  i, press_i_near,  nbar2_i_near) ;
 			if (press_i_near < 0.) { 
		  			cout << " INTERPOLATION FLUID P --> negative pressure " << endl ;
		  			abort();
			}
 			if (nbar2_i_near < 0.) { 
		  			cout << " INTERPOLATION FLUID P --> negative density " << endl ;
		  			abort();
			}
  	}
 	else if (Placei == 0 ) {	// two fluids
   
			// interpolation of press, nbar1 and nbar2	
		   
			Eos_bf_tabul::interpol_2d_Tbl3d(i_near, j_near, k_near, *mu1_tab, *mu2_tab,
												  *press_tab, *n1_tab, *n2_tab, *d2psdmu1dmu2_tab  , 
		      								  mu1, mu2, press_i_near, nbar1_i_near , nbar2_i_near) ; 
   
			// interpolation of malpha
			double der1 = 0., der2 = 0.	;

		  	Eos_bf_tabul::interpol_2d_Tbl3d(i_near, j_near, k_near, *mu1_tab, *mu2_tab,
												  *dpsddelta_car_tab, *dn1sddelta_car_tab, *dn2sddelta_car_tab, *null_tab,
		      								  mu1, mu2, malpha_i_near,  der1 , der2) ; 
   
 			if (press_i_near < 0.) { 
		  			//cout << press_i_near << " --> negative pressure " << endl ;
		  			press_i_near = 0 ; // interpolation for very small values of press can lead to press <0...
					// to check the smallness of press -> abort();
			}
 			if (nbar1_i_near < 0.) { 
		  			//cout <<  nbar1_i_near << " --> negative density nbar1 " << endl ;
		  			nbar1_i_near = 0 ; // interpolation for very small values of nbar1 can lead to nbar1  <0...
			}
			if (nbar2_i_near < 0.) { 
		  			//cout <<  nbar2_i_near << " --> negative density nbar2 " << endl ;
		  			nbar2_i_near = 0 ; // interpolation for very small values of nbar2 can lead to nbar2 <0...
			}

	}


  /***********************************
 	* 2D interpolation on Slice i + 1 *
	***********************************/ 
   
	// Looking for the table to be used (concerning protons only, neutrons only or both fluids).
  	// -----------------------------------------------------------------------------------------

	int Placei1 = 0 ; // 0 = two fluids, 1 = only neutrons (fluid 1), 2 = only protons (fluid 2), 3 = no fluids

  	int i_nearN_i1 = 0;
  	int j_nearN_i1 = 0;
  	int i_nearP_i1 = 0;
  	int j_nearP_i1 = 0; 
  
  	if ( mu1 > m_1 ) // both fluids are present or only the neutron fluid (fluid 1)
  	{
		// to find if one or two fluid(s) is (are) present, we compare the position under consideration  
      // with the curve n_p = 0 for the EoS which is used.
		// Note that the following procedure is adapted to DDH and DDHdelta EoSs (close to beta eq. and corotation)
		// but the curve n_p = 0 could be possibly much more complicated depending on the EoS...
		while ( ( (*delta_car_n0)(i_nearN_i1,0) <= (*delta_car_tab)(i_near+1, j_near, k_near) ) && ( ( n_deltaN-1 ) > i_nearN_i1 ) ) {
  	  		i_nearN_i1++ ;
  		}
  		if (i_nearN_i1 != 0) { 
	  		i_nearN_i1 -- ; 
		}
		while ( ( (*mu1_n0)(i_nearN_i1,j_nearN_i1) <= mu1 ) && ( ( n_mu1N-1 ) > j_nearN_i1 ) ) {
	  		j_nearN_i1++ ;
  		}
  		if (j_nearN_i1 != 0) { 
	  		j_nearN_i1 -- ; 
  		}

		// some checks
		if ( ( (*delta_car_n0)(i_nearN_i1,0) > (*delta_car_tab)(i_near+1, j_near, k_near)  ) || ((*delta_car_tab)(i_near+1, j_near, k_near)  > (*delta_car_n0)(i_nearN_i1+1,0) ) )
		{
 	  		cout << " bad location of delta_car_tab_i+1 in *delta_car_n0 (courbe limite np = 0) " << endl ;
 	  		cout << (*delta_car_n0)(i_nearN_i1,0) << "  " << (*delta_car_tab)(i_near+1, j_near, k_near) <<  "  " <<  (*delta_car_n0)(i_nearN_i1+1,0) << endl;
 			abort();
 		}
	 	if ( ( (*mu1_n0)(i_nearN_i1,j_nearN_i1) > mu1  ) || (mu1 > (*mu1_n0)(i_nearN_i1,j_nearN_i1+1) ) ) 
		{
 	  		cout << " bad location of mu_n in  *mu1_n0 (limit curve np = 0) " << endl ;
 	  		cout << (*mu1_n0)(i_nearN_i1,j_nearN_i1) << "  " << mu1 <<  "  " <<  (*mu1_n0)(i_nearN_i1,j_nearN_i1+1) << endl;
 			abort();
 		}	

		double aN_i1, bN_i1;
		aN_i1 = ((*mu2_n0)(i_nearN_i1,j_nearN_i1+1) -  (*mu2_n0)(i_nearN_i1,j_nearN_i1) ) / ((*mu1_n0)(i_nearN_i1,j_nearN_i1+1) -  (*mu1_n0)(i_nearN_i1,j_nearN_i1) ) ;
		bN_i1 = (*mu2_n0)(i_nearN_i1,j_nearN_i1)  - aN_i1 * (*mu1_n0)(i_nearN_i1,j_nearN_i1) ;
		double zN_i1 = aN_i1 * mu1 + bN_i1 ;

		if (zN_i1 <  mu2) 
		{	
	  			Placei1 = 0 ;		// two fluids
		}
		else {			
	 			Placei1 = 1 ;		// fluid 1 only
		}	

	}
	else	// both fluids are present or only the proton fluid (fluid 2) or no fluids !
  	{
    	if ( mu2 <= m_2) 
		{
				Placei = 3 ; 		// no fluids
		}
		else 
		{
		
			// to find if one or two fluid(s) is (are) present, we compare the position under consideration  
     		// with the curve n_n = 0 for the EoS which is used.
			// Note that the following procedure is adapted to DDH and DDHdelta EoSs (close to beta eq. and corotation)
			// but the curve n_n = 0 could be possibly much more complicated depending on the EoS... 	  
		  	while ( ( (*delta_car_p0)(i_nearP_i1,0) <= (*delta_car_tab)(i_near+1, j_near, k_near)) && ( ( n_deltaP-1 ) > i_nearP_i1) ) {
	    		i_nearP_i1++ ;
	  		}
	  		if (i_nearP_i1 != 0) { 
	    		i_nearP_i1 -- ; 
	  		}
			while ( ( (*mu2_p0)(i_nearP_i1,j_nearP_i1) <= mu2 ) && ( ( n_mu2P-1 ) > j_nearP_i1 ) ) {
	    		j_nearP_i1++ ;
	  		}
	  		if (j_nearP_i1 != 0) { 
	    		j_nearP_i1 -- ; 
	  		}

			// some checks
 	      if  ( ( (*delta_car_p0)(i_nearP_i1,0) > (*delta_car_tab)(i_near+1, j_near, k_near)  ) || ((*delta_car_tab)(i_near+1, j_near, k_near)  > (*delta_car_p0)(i_nearP_i1+1,0) ) )
			{
 	  			cout << " bad location of delta_car_tab_i+1 in *delta_car_p0 (courbe limite nn = 0) " << endl ;
				cout << (*delta_car_p0)(i_nearP_i1,0) << "  " << (*delta_car_tab)(i_near+1, j_near, k_near) <<  "  " <<  (*delta_car_p0)(i_nearP_i1+1,0) << endl;
 				abort();
 			}
	 		if ( ( (*mu2_p0)(i_nearP_i1,j_nearP_i1) > mu2  ) || (mu2 > (*mu2_p0)(i_nearP_i1,j_nearP_i1+1) ) ) {
 	  			cout << " bad location of mu_p in  *mu2_p0 (limit curve nn = 0) " << endl ;
	 	      cout << (*mu2_p0)(i_nearP_i1,j_nearP_i1) << "  " << mu2 <<  "  " <<  (*mu2_p0)(i_nearP_i1,j_nearP_i1+1) << endl;
 				abort();
 			}	

	  		double aP_i1, bP_i1;
	  		aP_i1 = ( (*mu2_p0)(i_nearP_i1,j_nearP_i1+1) -  (*mu2_p0)(i_nearP_i1,j_nearP_i1) ) / ( (*mu1_p0)(i_nearP_i1,j_nearP_i1+1) -  (*mu1_p0)(i_nearP_i1,j_nearP_i1) ) ;
	  		bP_i1 = (*mu2_p0)(i_nearP_i1,j_nearP_i1)  - aP_i1 * (*mu1_p0)(i_nearP_i1,j_nearP_i1) ;
	  		double yP_i1 = (mu2- bP_i1) /aP_i1 ;
		
	 		if (yP_i1 <  mu1) 
			{ 
	    			Placei = 0; 		// two fluids
	  		}
	  		else 
			{	
	    			Placei = 2 ;   	// fluid 2 only
	  		}
	  
		}

	}	// end of the search of the table to be used.
 
	/*	
	cout << mu2* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "   "  << mu1 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13) << endl ;
	cout << " Placei = " << Placei << endl ;
	*/

	// Interpolation in the different areas
	if (Placei == 3 ) // no fluids
	{
	  		press_i1 = 0. ;
			nbar1_i1 = 0. ;
			nbar2_i1 = 0. ;
			malpha_i1 = 0.;
  	}
  	else if (Placei1 == 1 ) // fluid 1 only
	{	
			malpha_i1 = 0. ;
			nbar2_i1 = 0. ;
			int i =0;
			interpol_herm(*mu1_N, *press_N, *n_n_N, mu1,  i, press_i1 , nbar1_i1 ) ;
			if (press_i1 < 0.) { 
		  		cout << " INTERPOLATION FLUID N i+1 --> negative pressure " << endl ;
		  		abort();
			}
 			if (nbar1_i1 < 0.) { 
		  		cout << " INTERPOLATION FLUID N i+1--> negative density " << endl ;
		  		abort();
			}
  }
  else if (Placei1 == 2 ) // fluid 2 only
  {	
			malpha_i1 = 0.;
			nbar1_i1 = 0. ;
			int i =0;
			interpol_herm( *mu2_P, *press_P, *n_p_P, mu2,  i, press_i1,  nbar2_i1) ;
			if (press_i1 < 0.) { 
		  		cout << " INTERPOLATION FLUID P i+1--> negative pressure " << endl ;
		  		abort();
		  	}
 			if (nbar2_i1 < 0.) { 
		  		cout << " INTERPOLATION FLUID P i+1 --> negative density " << endl ;
		  		abort();
			}
  }
  else if (Placei1 == 0 ) { // two fluids
  
		// interpolation of press, nbar1 and nbar2	

		Eos_bf_tabul::interpol_2d_Tbl3d(i1, j_near, k_near, *mu1_tab, *mu2_tab,
												  *press_tab, *n1_tab, *n2_tab, *d2psdmu1dmu2_tab  , 
		      								  mu1, mu2, press_i1, nbar1_i1 , nbar2_i1) ; 
	
		// interpolation of malpha
		double der1b = 0., der2b = 0.;

		Eos_bf_tabul::interpol_2d_Tbl3d(i1, j_near, k_near, *mu1_tab, *mu2_tab,
												  *dpsddelta_car_tab, *dn1sddelta_car_tab, *dn2sddelta_car_tab, *null_tab,
		      								  mu1, mu2, malpha_i1,  der1b , der2b ) ; 					


 		if (press_i1 < 0.) { 
		  		//cout << press_i1 << " --> negative pressure " << endl ;
		  		press_i1 = 0 ; // interpolation for very small values of press can lead to press <0...
				// to check the smallness of press -> abort();
		}
 		if (nbar1_i1 < 0.) { 
		  		//cout <<  nbar1_i1 << " --> negative density nbar1 " << endl ;
		  		nbar1_i1 = 0 ; // interpolation for very small values of nbar1 can lead to nbar1  <0...
		}
		if (nbar2_i1 < 0.) { 
		  		//cout <<  nbar2_i1 << " --> negative density nbar2 " << endl ;
		  		nbar2_i1 = 0 ; // interpolation for very small values of nbar2 can lead to nbar2 <0...
		}

	}


  /***********************************
 	* linear interpolation in Delta^2 *
	***********************************/ 
 
 	double x1  = (*delta_car_tab)(i_near, 0, 0) ;
 	double x2  = (*delta_car_tab)(i1, 0, 0) ;
 	double x12 = x1-x2 ;
   
 	//pressure 
 	double y1 = press_i_near;
 	double y2 = press_i1;
 	double a  = (y1-y2)/x12 ;
 	double b  = (x1*y2-y1*x2)/x12 ;
	press  = delta2*a+b ; 


  	//nbar1
  	double y1_y = nbar1_i_near;
  	double y2_y = nbar1_i1;
  	double a_y  = (y1_y-y2_y)/x12 ;
  	double b_y  = (x1*y2_y-y1_y*x2)/x12 ;
	nbar1  = delta2*a_y+b_y ; 

  	//nbar2
  	double y1_z = nbar2_i_near;
  	double y2_z = nbar2_i1;
  	double a_z  = (y1_z-y2_z)/x12 ;
  	double b_z  = (x1*y2_z-y1_z*x2)/x12 ;
	nbar2  = delta2*a_z+b_z ; 
  
  	// for alpha 
  	double y1_alpha = malpha_i_near;
  	double y2_alpha = malpha_i1;
  	double a_alpha  = (y1_alpha-y2_alpha)/x12 ;
  	double b_alpha  = (x1*y2_alpha-y1_alpha*x2)/x12 ;
	malpha  = delta2*a_alpha+b_alpha ; 
	alpha = -malpha ;
 
	delete null_tab; 
   return;

}

// routine used in interpol_3d_bifluid to perform a 2D thermodynamically consistent interpolation
// on each slice of constant delta_car
//----------------------------------------------------------------------------------------------
void Eos_bf_tabul::interpol_2d_Tbl3d(const int i, const int j, const int k,  const Tbl& ytab, const Tbl& ztab,
												 const Tbl& ftab, const Tbl& dfdytab, const Tbl& dfdztab, const Tbl& d2fdydztab,
		      								 const double y, const double z, double& f, double& dfdy, double& dfdz) const
{

  assert(dfdytab.dim == ftab.dim ) ;
  assert(dfdztab.dim == ftab.dim ) ;
  assert(d2fdydztab.dim == ftab.dim ) ;
    
  int j1 = j+1 ; 
  int k1 = k+1 ;

  double dy = ytab(i, j1, k) - ytab(i, j, k) ;
  double dz = ztab(i, j, k1) - ztab(i, j, k) ;

  double u = (y - ytab(i, j, k)) / dy ;
  double v = (z - ztab(i, j, k)) / dz ;

  double u2 = u*u ; double v2 = v*v ;
  double u3 = u2*u ; double v3 = v2*v ;

  double psi0_u = 2.*u3 - 3.*u2 + 1. ;
  double psi0_1mu = -2.*u3 + 3.*u2 ;
  double psi1_u = u3 - 2.*u2 + u ;
  double psi1_1mu = -u3 + u2 ;
  double psi0_v = 2.*v3 - 3.*v2 + 1. ;
  double psi0_1mv = -2.*v3 + 3.*v2 ;
  double psi1_v = v3 - 2.*v2 + v ;
  double psi1_1mv = -v3 + v2 ;

  f = ftab(i, j,  k)  * psi0_u * psi0_v
    + ftab(i, j1, k)  * psi0_1mu * psi0_v 
    + ftab(i, j,  k1) * psi0_u * psi0_1mv
    + ftab(i, j1, k1) * psi0_1mu * psi0_1mv ;

  f += (dfdytab(i, j,  k)  * psi1_u * psi0_v
		- dfdytab(i, j1, k)  * psi1_1mu * psi0_v
		+ dfdytab(i, j,  k1) * psi1_u * psi0_1mv
		- dfdytab(i, j1, k1) * psi1_1mu * psi0_1mv) * dy ;

  f += (dfdztab(i, j,  k)  * psi0_u * psi1_v
		+ dfdztab(i, j1, k)  * psi0_1mu * psi1_v
		- dfdztab(i, j,  k1) * psi0_u * psi1_1mv
		- dfdztab(i, j1, k1) * psi0_1mu * psi1_1mv) * dz ;
  
  f += (d2fdydztab(i, j,  k) * psi1_u * psi1_v
		- d2fdydztab(i, j1, k) * psi1_1mu * psi1_v
		- d2fdydztab(i, j, k1) * psi1_u * psi1_1mv 
		+ d2fdydztab(i, j1, k1) * psi1_1mu * psi1_1mv) * dy * dz ;

  double dpsi0_u = 6.*(u2 - u) ;
  double dpsi0_1mu = 6.*(u2 - u) ;
  double dpsi1_u = 3.*u2 - 4.*u + 1. ;
  double dpsi1_1mu = 3.*u2 - 2.*u ;
 
  dfdy = (ftab(i, j, k) 	* dpsi0_u * psi0_v
	  	  - ftab(i, j1, k) 	* dpsi0_1mu * psi0_v
	  	  + ftab(i, j, k1) 	* dpsi0_u * psi0_1mv
	  	  - ftab(i, j1, k1)  * dpsi0_1mu * psi0_1mv ) / dy;

  dfdy += (dfdytab(i, j, k) 	* dpsi1_u * psi0_v
			+ dfdytab(i, j1, k) 	* dpsi1_1mu * psi0_v
		   + dfdytab(i, j, k1) 	* dpsi1_u * psi0_1mv
		   + dfdytab(i, j1, k1) * dpsi1_1mu * psi0_1mv) ;

  dfdy += (dfdztab(i, j, k) 	* dpsi0_u* psi1_v
		   - dfdztab(i, j1, k) 	* dpsi0_1mu * psi1_v
		   - dfdztab(i, j, k1) 	* dpsi0_u * psi1_1mv
			+ dfdztab(i, j1, k1) * dpsi0_1mu * psi1_1mv) * dz /dy ;
  
 dfdy += (d2fdydztab(i, j, k) 	* dpsi1_u * psi1_v
			+ d2fdydztab(i, j1, k) 	* dpsi1_1mu * psi1_v
			- d2fdydztab(i, j, k1) 	* dpsi1_u * psi1_1mv
			- d2fdydztab(i, j1, k1) * dpsi1_1mu * psi1_1mv) * dz ;

  double dpsi0_v = 6.*(v2 - v) ;
  double dpsi0_1mv = 6.*(v2 - v) ;
  double dpsi1_v = 3.*v2 - 4.*v + 1. ;
  double dpsi1_1mv = 3.*v2 - 2.*v ;

  dfdz = (ftab(i, j, k) * psi0_u * dpsi0_v
    		+ ftab(i, j1, k) * psi0_1mu * dpsi0_v
    		- ftab(i, j, k1) * psi0_u * dpsi0_1mv
    		- ftab(i, j1, k1)  * psi0_1mu * dpsi0_1mv) / dz ;

  dfdz += (dfdytab(i, j, k) * psi1_u * dpsi0_v
			- dfdytab(i, j1, k) * psi1_1mu * dpsi0_v
			- dfdytab(i, j, k1) * psi1_u * dpsi0_1mv
			+ dfdytab(i, j1, k1) * psi1_1mu * dpsi0_1mv) * dy / dz ;

  dfdz += (dfdztab(i, j, k) * psi0_u * dpsi1_v
			+ dfdztab(i, j1, k) * psi0_1mu * dpsi1_v
			+ dfdztab(i, j, k1) * psi0_u * dpsi1_1mv
			+ dfdztab(i, j1, k1) * psi0_1mu * dpsi1_1mv) ;
  
  dfdz += (d2fdydztab(i, j, k) * psi1_u* dpsi1_v
			- d2fdydztab(i, j1, k) * psi1_1mu * dpsi1_v
			+ d2fdydztab(i, j, k1) * psi1_u* dpsi1_1mv
			- d2fdydztab(i, j1, k1) * psi1_1mu * dpsi1_1mv) * dy;

	return ;

}

// Conversion functions 
// ---------------------

//This function is necessary for "Et_rot", which needs an eos of type "Eos"
//But this eos is not used in the code, except for the construction of the star
Eos* Eos_bf_tabul::trans2Eos() const {

  Eos_poly* eos_simple = new Eos_poly(2.,0.016,1.008) ; // we can take whatever we want that makes sense as parameters
  
  return eos_simple ;
}

}
