/*
 *  Methods of class Eos_bf_tabul.
 *
 *  (see file eos_bifluid.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char eos_bf_tabul_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2015/06/10 14:39:17  a_sourie
 * New class Eos_bf_tabul for tabulated 2-fluid EoSs and associated functions for the computation of rotating stars with such EoSs.
 *
 * Revision 1.12  2014/03/06 15:53:35  j_novak
 * Eos_compstar is now Eos_compOSE. Eos_tabul uses strings and contains informations about authors.
 *
 * Revision 1.11  2012/11/09 10:32:44  m_bejger
 * Implementing der_ener_ent_p
 *
 *
 * Revision 2.0  2000/11/22  19:31:30  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers C // Are they all useful ?
#include <cstdlib>  
#include <cstring>
#include <cmath>
#include <stdlib.h>


// Headers Lorene // Are they all useful ??
#include "eos_bifluid.h"
#include "param.h"
#include "eos.h" 
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"
#include "cmp.h"
#include "nbr_spx.h"

namespace Lorene {

void interpol_linear(const Tbl& xtab, const Tbl& ytab, 
                     double x, int& i, double& y) ;

void interpol_herm(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
		   double x, int& i, double& y, double& dy) ;

void interpol_mixed_3d(const Tbl&, const Tbl&, const Tbl&, const Tbl&, 
		      const Tbl&, const Tbl&, const Tbl&,
		      double, double, double, double&, double&, double&) ;


void interpol_mixed_3d_new(double m1, double m2, const Tbl& xtab, const Tbl& ytab, const Tbl& ztab, 
			const Tbl& ftab, const Tbl& dfdytab, const Tbl& dfdztab, const Tbl& d2fdydztab,
			const Tbl& dlpsddelta_car, const Tbl&  d2lpsdlent1ddelta_car, const Tbl& d2lpsdlent2ddelta_car, 
			const Tbl&  mu2_P, const Tbl&  n_p_P,   const Tbl& press_P,
			const Tbl& mu1_N, const Tbl& n_n_N, const Tbl& press_N,
			const Tbl& delta_car_n0, const Tbl& mu1_n0, const Tbl& mu2_n0,
			const Tbl& delta_car_p0, const Tbl& mu1_p0, const Tbl& mu2_p0, 
		      double x, double y, double z, double& f, double& dfdy, double& dfdz, double& alpha) ;

void interpolation_brutale(double x,double y, double z, // localisation
			// sommet vitesse inferieure, mun inf, mup_inf
			double delta111, double delta211,
			double mu1_111,  double mu1_121,  double mu2_111,  double mu2_112,  double mu1_211,  double mu1_221,  double mu2_211, double  mu2_212,
			double p_111,  double p_121,  double p_112,  double p_122,  double p_211,  double p_221,  double p_212,  double p_222,  
			double n1_111,  double n1_121,  double n1_112,  double n1_122,  double n1_211,  double n1_221,  double n1_212, double  n1_222, 
			double n2_111,  double n2_121,  double n2_112,  double n2_122,  double n2_211,  double n2_221,  double n2_212,  double n2_222, 
			double cross_111,  double cross_121,  double cross_112,  double cross_122,  double cross_211, double  cross_221,  double cross_212,  double cross_222,
			double& f, double& dfdy, double& dfdz)  ;

void interpolation_brutale_mod(double x,double y, double z, // localisation
			// sommet vitesse inferieure, mun inf, mup_inf
			double delta111, double delta211,
			double mu1_111,  double mu1_121,  double mu2_111,  double mu2_112,  double mu1_211,  double mu1_221,  double mu2_211, double  mu2_212,
			double p_111,  double p_121,  double p_112,  double p_122,  double p_211,  double p_221,  double p_212,  double p_222,  
			double n1_111,  double n1_121,  double n1_112,  double n1_122,  double n1_211,  double n1_221,  double n1_212, double  n1_222, 
			double n2_111,  double n2_121,  double n2_112,  double n2_122,  double n2_211,  double n2_221,  double n2_212,  double n2_222, 
			double& f, double& dfdy, double& dfdz)  ;

			
/*void interpol_mixed_3d_test3(const Tbl&, const Tbl&, const Tbl&, const Tbl&, 
		      const Tbl&, const Tbl&, const Tbl&,
		      double, double, double, double&, double&, double&) ;

void interpol_mixed_3d_test2(const Tbl&, const Tbl&, const Tbl&, const Tbl&, 
		      const Tbl&, const Tbl&, const Tbl&, const Tbl&, const Tbl&,
		      double, double, double, double&, double&, double&) ;
*/
void interpol_mixed_3d_mod(const Tbl&, const Tbl&, const Tbl&, const Tbl&, 
		      const Tbl&, const Tbl&,
		      double, double, double, double&,  double&, double&) ;


			//----------------------------//
			//   	Constructors	      //
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
  fich >> tablename ;

  read_table() ;    
}
                    

			//--------------//
			//  Destructor  //
			//--------------//

Eos_bf_tabul::~Eos_bf_tabul(){

	delete logent1 ;
	delete logent2 ;
	delete delta_car ; 
	delete logp ;
	delete dlpsdlent1 ;
    	delete dlpsdlent2 ;
	delete d2lpsdlent1dlent2 ;
	delete d2lpsdlent1dlent1 ;
	delete d2lpsdlent2dlent2 ;
	delete dlpsddelta_car ;
	delete d2lpsdlent1ddelta_car ;
	delete d2lpsdlent2ddelta_car ;
	delete d3lpsdlent1dlent2ddelta_car ; // if necessary
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
    
    tablename = eosi.tablename; 
    ent1_min = eosi.ent1_min ;
    ent1_max = eosi.ent1_max ; 
    ent2_min = eosi.ent2_min ; 
    ent2_max = eosi.ent2_max ; 
    delta_car_min = eosi.delta_car_min ;
    delta_car_max = eosi.delta_car_max ;
    logent1 = eosi.logent1 ;  
    logent2 = eosi.logent2 ;
    delta_car = eosi.delta_car ; 
    logp = eosi.logp ;
    dlpsdlent1 = eosi.dlpsdlent1 ;
    dlpsdlent2 = eosi.dlpsdlent2 ;
    d2lpsdlent1dlent2 = eosi.d2lpsdlent1dlent2 ;
    dlpsddelta_car = eosi.dlpsddelta_car ;
    d2lpsdlent1ddelta_car = eosi.d2lpsdlent1ddelta_car ;
    d2lpsdlent2ddelta_car = eosi.d2lpsdlent2ddelta_car;
    d3lpsdlent1dlent2ddelta_car = eosi.d3lpsdlent1dlent2ddelta_car;
    delta_car_n0= eosi.delta_car_n0 ;
    mu1_n0= eosi.mu1_n0 ;
    mu2_n0= eosi.mu2_n0 ;
    delta_car_p0= eosi.delta_car_p0 ;
    mu1_p0 = eosi.mu1_p0;
    mu2_p0 = eosi.mu2_p0 ;
    mu1_N  = eosi.mu1_N  ;
    n_n_N  = eosi.n_n_N;
    press_N  = eosi.press_N;
    mu2_P  = eosi.mu2_P ;
    n_p_P  = eosi.n_p_P;
    press_P  = eosi.press_P;
    
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
    "EOS of class Eos_bf_tabul : tabulated EOS depending on three variables (relative velocity and enthalpies of neutrons and protons)."
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
  double m_b_si = 1.66E-27; 	

  //-------------------------------------------
  //---------- Table à deux fluides -----------
  //-------------------------------------------

  ifstream fich(tablename.data()) ;
	//DEBUG : cout << tablename << endl;
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
	//cout << authors << endl ;
	for (int i=0; i<5; i++) {		//  jump over the file
	  fich.ignore(1000, '\n') ;             // header
    	}                                       //

    	int nbp ;
    	fich >> nbp ; fich.ignore(1000, '\n') ;   // number of data
	//cout << nbp << endl ;
	
	int n_delta, n_mu1, n_mu2; // number of different values for : delta_car, mu_n, mu_p
	fich >> n_delta;  fich.ignore(1000, '\n') ;
	fich >> n_mu1;  fich.ignore(1000, '\n') ;
	fich >> n_mu2;  fich.ignore(1000, '\n') ;
	//cout << "n_delta : " << n_delta << endl ;
	//cout << "n_mu1 : " << n_mu1 << endl ;
	//cout << "n_mu2 : " <<n_mu2 << endl ;

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

	for (int i=0; i<3; i++) {		//  jump over the table
	  fich.ignore(1000, '\n') ;             // description
    	}                                      
	
	logent1 = new Tbl(n_delta, n_mu1, n_mu2) ;
	logent2 = new Tbl(n_delta, n_mu1, n_mu2) ;
	delta_car = new Tbl(n_delta, n_mu1, n_mu2) ; 
	logp = new Tbl(n_delta, n_mu1, n_mu2) ;
	dlpsdlent1 = new Tbl(n_delta, n_mu1, n_mu2) ;
        dlpsdlent2 = new Tbl(n_delta, n_mu1, n_mu2) ;
	d2lpsdlent1dlent2 = new Tbl(n_delta, n_mu1, n_mu2) ; 
      	dlpsddelta_car = new Tbl(n_delta, n_mu1, n_mu2) ;
	d2lpsdlent1ddelta_car = new Tbl(n_delta, n_mu1, n_mu2) ;
	d2lpsdlent2ddelta_car = new Tbl(n_delta, n_mu1, n_mu2) ;
	d3lpsdlent1dlent2ddelta_car = new Tbl(n_delta, n_mu1, n_mu2) ;
	d2lpsdlent1dlent1 = new Tbl(n_delta, n_mu1, n_mu2) ; 
	d2lpsdlent2dlent2 = new Tbl(n_delta, n_mu1, n_mu2) ; 

	logent1->set_etat_qcq() ;
	logent2->set_etat_qcq() ;
	delta_car->set_etat_qcq() ; 
	logp->set_etat_qcq() ;
	dlpsdlent1->set_etat_qcq() ;
	dlpsdlent2->set_etat_qcq() ;
	d2lpsdlent1dlent2->set_etat_qcq() ; 
      	dlpsddelta_car->set_etat_qcq() ;
	d2lpsdlent1ddelta_car->set_etat_qcq() ;
	d2lpsdlent2ddelta_car->set_etat_qcq() ;
	d3lpsdlent1dlent2ddelta_car->set_etat_qcq()  ;
 	d2lpsdlent1dlent1->set_etat_qcq() ; 
	d2lpsdlent2dlent2->set_etat_qcq() ; 

	//------------------------------------------------------
	// We have to choose the right unites (SI, cgs , LORENE)
	//------------------------------------------------------
	// sor far, all the calculations are done with Lorene units

    	int no ;

	//Grandeurs issues de la table à deux fluides (unités physiques adaptées)
	double mu1_MeV, mu2_MeV, n1_fm3, n2_fm3; 
	double Knp_Mev_2, press_MeV_fm3;
	double d2press_dmu1dmu2_MeV_fm3, dn1_ddelta_car_fm3, dn2_ddelta_car_fm3;

	//Grandeurs stockées (unité Lorene)
    	double delta_car_adim, mu1_adim, mu2_adim, n1_01fm3, n2_01fm3, Knp_adim  ;
	double press_adim,  dpress_ddelta_car_adim, dn1_ddelta_car_adim, dn2_ddelta_car_adim ;
	double d2press_dmu1dmu2_adim;
	double d3press_dmu1dmu2ddelta_car_adim; // ce terme n'est donné que pour la table analytique. Je le stocke en unité Lorene pour plus de simplicité

	double hbarc = 197.33 ; //en Mev*fm
	double hbarc3 = hbarc * hbarc * hbarc ;


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
			   // cout << "neg" << dpress_ddelta_car_adim << endl ;
    		            fich >> press_MeV_fm3 ;
			    press_adim = press_MeV_fm3 * mevpfm3 ;	
			    fich >>  d2press_dmu1dmu2_MeV_fm3 ;
			    d2press_dmu1dmu2_adim = d2press_dmu1dmu2_MeV_fm3 * (10. * m_b_si * v_unit * v_unit ) / mev_si ;
		            fich >> dn1_ddelta_car_fm3 ;
			    dn1_ddelta_car_adim = dn1_ddelta_car_fm3 * 10. ;
			    fich >> dn2_ddelta_car_fm3 ;
			    dn2_ddelta_car_adim = dn2_ddelta_car_fm3 * 10. ;
			   // fich >> d3press_dmu1dmu2ddelta_car_adim ;  // depend de la table -------------------------------------------------------------------
			   
	 	            //To check the reading is well done :
	 		    //cout << j << "   " << k << "    " << l << "   " << delta_car_adim << "  "<< mu1_MeV << "   "<< mu2_MeV << "   " <<  Knp_Mev_2 << "  " <<  dn2_ddelta_car_fm3  <<  endl;
		
			    fich.ignore(1000,'\n') ;    		
		
			    /* if ( (n1_fm3<0) || (n2_fm3<0) || (press_MeV_fm3  < 0)){ 
			      cout << "Eos_tabul::read_table(): " << endl ;
			      cout << "Negative value in table!" << endl ;
			      cout << "n_neutrons = " << n1_fm3 << ", n_protons " << n2_fm3 <<
				", p = " << press_MeV_fm3  << ", "<< endl ;
			      cout << "Aborting..." << endl ;
			      abort() ;		
			      }*/ // il peut y avoir des densités négatives selon la table à deux fluides utilisés...



			    /* redefinition de la masse --> permet de redefinir le minimum en log enthalpie !
			       double ww1 = 1.;
			       double ww2 = 1.;
			       if ((k==0) && (l==0) ) {ww1 = mu1_adim; ww2 = mu2_adim;} ;
			       cout << setprecision(16) ; 
			       cout << ww1 ;
			       mu1_adim = mu1_adim/ww1 ; //+1e-14 ;
			       mu2_adim = mu2_adim/ww2  ;
			    */

			    // on enregistre les grandeurs dans des tableaux sans passer par des logarithmes !!! --> si on garde il faudra choisir des meilleurs noms de tableaux !
			    logent1->set(j,k,l) = mu1_adim ;
			    logent2->set(j,k,l) = mu2_adim ;
			    delta_car->set(j,k,l) = delta_car_adim ;
// 			    logp->set(j,k,l) = log(press_adim) ;
 			    if (press_adim<=0) {press_adim = 0. ;}
 			    logp->set(j,k,l) = press_adim ;					// non en log	   
// 			    dlpsdlent1->set(j,k,l) = n1_01fm3 /press_adim ;
 			    dlpsdlent1->set(j,k,l) = n1_01fm3 ; 				// non en log					    
// 			    dlpsdlent2->set(j,k,l) = n2_01fm3 /press_adim ;
 			    dlpsdlent2->set(j,k,l) = n2_01fm3 ;					// non en log
 			    
			    if ((n1_01fm3 < 1e-16) && (n2_01fm3 <1e-16 )) {
			      d2press_dmu1dmu2_adim  = 0. ;
			      dpress_ddelta_car_adim = 0. ;
			       dn1_ddelta_car_adim =0. ;
			      dn2_ddelta_car_adim  = 0. ;
			    }
// 			    d2lpsdlent1dlent2->set(j,k,l) = d2press_dmu1dmu2_adim / press_adim - n1_01fm3 /press_adim * n2_01fm3 /press_adim    ; 
 			    d2lpsdlent1dlent2->set(j,k,l) = d2press_dmu1dmu2_adim  ; 		// non en log
			    dlpsddelta_car->set(j,k,l) = dpress_ddelta_car_adim  ;
			    d2lpsdlent1ddelta_car->set(j,k,l) =  dn1_ddelta_car_adim  ;
			    d2lpsdlent2ddelta_car->set(j,k,l) =  dn2_ddelta_car_adim  ;
			    d3lpsdlent1dlent2ddelta_car->set(j,k,l) =  d3press_dmu1dmu2ddelta_car_adim ;
		
    			}  

			  fich.ignore(1000, '\n') ;  // BIEN VERIFIER QU ON SAUTE UNE LIGNE DANS TOUTES LES TABLES MAIS NORMALEMENT OUI !

 	 	}

		fich.ignore(1000, '\n') ;

	  }
	//abort();
 	delta_car_min = (*delta_car)(0, 0, 0)  ;
 	delta_car_max = (*delta_car)(n_delta-1, 0, 0) ;

	// si on passe par des log
	/*ent1_min = pow( double(10), (*logent1)(0, 0, 0))  ;  
	  ent1_max = pow( double(10), (*logent1)(0, n_mu1-1, 0)) ;
	  ent2_min = pow( double(10), (*logent2)(0, 0, 0)) ;
	  ent2_max = pow( double(10), (*logent2)(0, 0, n_mu2-1));
	*/

	/*ent1_min = log((*logent1)(0, 0, 0)) ;		 
	  ent1_max = log((*logent1)(0, n_mu1-1, n_mu2-1)) ;	
	  ent2_min = log((*logent2)(0, 0, 0)) ;		
	  ent2_max = log((*logent2)(0, n_mu1-1, n_mu2-1));
	*/

	// sans passer par les log !
	ent1_min = (*logent1)(0, 0, 0)  ;  // en fait il s'agit de mu
	ent1_max = (*logent1)(0, n_mu1-1, 0) ;
	ent2_min = (*logent2)(0, 0, 0) ;
	ent2_max = (*logent2)(0, 0, n_mu2-1);
	
	//cout <<"DEBUG : ent1_min "<< ent1_min << endl ;
	//cout <<"DEBUG : ent1_max "<< ent1_max << endl ;   
	//cout <<"DEBUG : ent2_min "<< ent2_min << endl ; 
	//cout <<"DEBUG : ent2_max "<< ent2_max << endl ; 
	
	
	
	
	//cout <<"DEBUG : ent1_min en MeV "<< ent1_min / mev_si * (m_b_si * v_unit * v_unit ) << endl ;
	//cout <<"DEBUG : ent1_max en MeV "<< ent1_max / mev_si * (m_b_si * v_unit * v_unit ) << endl ;   
	//cout <<"DEBUG : ent2_min en MeV "<< ent2_min / mev_si * (m_b_si * v_unit * v_unit ) << endl ; 
	//cout <<"DEBUG : ent2_max en MeV "<< ent2_max / mev_si * (m_b_si * v_unit * v_unit )<< endl ; 
	//abort();
	fich.close();


  //---------------------------------------------------
  //---------- Table à un fluide : fluide N -----------
  //---------------------------------------------------

	ifstream fichN ;
	fichN.open("eos_bf_test_1_fluide_N.d") ;

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

    	int nbp_N ;
    	fichN >> nbp_N ; fichN.ignore(1000, '\n') ;   // number of data
	cout<< "nbp_N = " << nbp_N ; 
	int n_mu1_N; // number of different values for : delta_car, n_n, n_p // SEULS n_n sert le reste est inutile et pourra etre enlever dans une version finale
	fichN >> n_mu1_N;fichN.ignore(1000, '\n') ;
	
	//cout << n_mu1_N endl ; 

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
	for (int i=0; i<3; i++) {		//  jump over the table
	  fichN.ignore(1000, '\n') ;             // description
    	}                                      
	
	mu1_N = new Tbl(n_mu1_N) ;
	n_n_N = new Tbl(n_mu1_N) ;
	press_N = new Tbl(n_mu1_N) ;

	mu1_N ->set_etat_qcq() ;
	n_n_N->set_etat_qcq() ;
	press_N ->set_etat_qcq() ;

 	//Grandeurs issues de la table à un fluide (unités physiques adaptées)
	double mu1_MeV_N, n1_fm3_N,press_MeV_fm3_N;

	//Grandeurs stockées (unité Lorene)
   	double mu1_adimN, n1_01fm3N,press_adimN;

       for (int k=0 ; k < n_mu1_N ; k++) {
   	
       		fichN >> mu1_MeV_N ;
		mu1_adimN = mu1_MeV_N * mev_si /  (m_b_si * v_unit * v_unit ) ;
    		fichN >> n1_fm3_N ; 
		n1_01fm3N = 10. * n1_fm3_N ;
    		fichN >> press_MeV_fm3_N;
		press_adimN = press_MeV_fm3_N * mevpfm3 ;
 		fichN.ignore(1000,'\n') ;  
	
		//cout  << k << "   " << press_MeV_fm3_N << endl; 	
		
		if ( (n1_01fm3N<0)  || (press_adimN < 0)){ 
		  cout << "Eos_tabul::read_table(): " << endl ;
		  cout << "Negative value in table!" << endl ;
		  cout << "n_neutrons = " << n1_01fm3N << 
		    ", p = " << press_adimN << ", "<< endl ;
		  cout << "Aborting..." << endl ;
		  abort() ;		
		}

		mu1_N ->set(k) = mu1_adimN ;
		n_n_N->set(k) =  n1_01fm3N ; 
		press_N ->set(k) = press_adimN ;

	}
	//abort();
	//cout << setprecision(16) ;
	//cout << endl ;
	//cout << m_b_si * v_unit * v_unit / mev_si << endl  ;
	//cout << m_b_si << endl  ;
	//cout << v_unit << endl  ;
	//cout <<  mev_si << endl  ;
	//cout << rhonuc_si << endl ;
	//abort() ;
	fichN.close();



  //---------------------------------------------------
  //---------- Table à un fluide : fluide P -----------
  //---------------------------------------------------


	ifstream fichP ;
	fichP.open("eos_bf_test_1_fluide_P.d") ;

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

	for (int i=0; i<5; i++) {		//  jump over the file
	  fichP.ignore(1000, '\n') ;             // header
    	}                                       //

    	int nbp_P ;
    	fichP >> nbp_P ; fichP.ignore(1000, '\n') ;   // number of data

	int n_mu2_P; // number of different values for : delta_car, n_n, n_p // SEUL n_mu2_P a un sens !!!!!!!!!!!! ENLEVER LES AUTRES DANS UNE FUTURE VERSION
	   fichP >> n_mu2_P;fichP.ignore(1000, '\n') ;

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
	
	for (int i=0; i<3; i++) {		//  jump over the table
	  fichP.ignore(1000, '\n') ;             // description
    	}                                      
	
	mu2_P = new Tbl(n_mu2_P) ;
	n_p_P = new Tbl(n_mu2_P) ;
	press_P = new Tbl(n_mu2_P) ;

	mu2_P ->set_etat_qcq() ;
	n_p_P->set_etat_qcq() ;
	press_P ->set_etat_qcq() ;
	  
	//Grandeurs issues de la table à un fluide (unités physiques adaptées)
	double mu2_MeV_P, n2_fm3_P,press_MeV_fm3_P;
    	
	//Grandeurs stockées (unité Lorene)
	double mu2_adimP,  n2_01fm3P, press_adimP;

 	for (int l=0 ; l < n_mu2_P ; l++) {
  	
		fichP >> mu2_MeV_P ;
		mu2_adimP = mu2_MeV_P * mev_si /  (m_b_si * v_unit * v_unit ) ;
    		fichP >> n2_fm3_P ; 
		n2_01fm3P = 10. * n2_fm3_P ;
    		fichP >> press_MeV_fm3_P;
		press_adimP = press_MeV_fm3_P * mevpfm3 ;
    	
		//cout << l << "   " << press_MeV_fm3_P << endl;
    		 
		fichP.ignore(1000,'\n') ;    		
		if ( (n2_01fm3P<0) || (press_adimP < 0)){ 
		  cout << "Eos_tabul::read_table(): " << endl ;
		  cout << "Pegative value in table!" << endl ;
		  cout <<", n_protons " << n2_01fm3P <<
		    ", p = " << press_adimP << ", "<< endl ;
		  cout << "Aborting..." << endl ;
		  abort() ;		
		}

	mu2_P ->set(l) = mu2_adimP ;
	n_p_P->set(l) =  n2_01fm3P ;
	press_P ->set(l) = press_adimP ;

	}
	//abort();
	fichP.close();


  //------------------------------------------------------------------
  //---------- COURBE LIMITE TABLE A DEUX FLUIDES : np = 0 -----------
  //------------------------------------------------------------------

	ifstream fich1 ;
	fich1.open("np=0.dat") ;
	// table classee en mun croissant !
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
	//cout << n_delta_n0 << "    " << n_mu1_n0 << endl;
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
		//  cout << o << "  " << p << "  " << delta_car_nn0  << "  " << mu1_MeV_nn0  << "  " <<  mu2_MeV_nn0<< endl ;
		  delta_car_n0 ->set(o,p) = delta_car_nn0;
		  mu1_n0 ->set(o,p) = mu1_MeV_nn0 * mev_si / (m_b_si * v_unit * v_unit ) ;
		  mu2_n0 ->set(o,p) = mu2_MeV_nn0 * mev_si / (m_b_si * v_unit * v_unit ) ;
		  
	       	}
		fich1.ignore(1000,'\n') ; 
	}
	//abort();
	fich1.close();
 
  //------------------------------------------------------------------
  //---------- COURBE LIMITE TABLE A DEUX FLUIDES : nn = 0 -----------
  //------------------------------------------------------------------

	ifstream fich2 ;
	fich2.open("nn=0.dat") ;
	// table classe en mup croissant ---------------> attention

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
	fich2.ignore(1000, '\n') ;             // Jump over the first header
	//cout << n_delta_p0 << "    " << n_mu2_p0 << endl;
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
	//	  cout << o << "  " << p << "  " << delta_car_np0  << "  " << mu1_MeV_np0  << "  " <<  mu2_MeV_np0<< endl ;
		  delta_car_p0 ->set(o,p) = delta_car_np0;
		  mu1_p0 ->set(o,p) = mu1_MeV_np0 * mev_si / (m_b_si * v_unit * v_unit ) ;
		  mu2_p0 ->set(o,p) = mu2_MeV_np0 * mev_si / (m_b_si * v_unit * v_unit ) ;
		
		  //cout << o << "  " << p << "  " << mu2_MeV_np0 << endl;
		}
		fich2.ignore(1000,'\n') ; 
	}
	//abort();
	fich2.close();

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//


// Complete computational routine giving all thermo variables
//-----------------------------------------------------------

void Eos_bf_tabul::calcule_tout(const Cmp& ent1, const Cmp& ent2, 
			       const Cmp& delta2, Cmp& nbar1, Cmp& nbar2,  
			       Cmp& ener, Cmp& press, Cmp& K_nn, Cmp& K_np, Cmp& K_pp,
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
    return ; 
  }
  nbar1.allocate_all() ;
  nbar2.allocate_all() ;
  ener.allocate_all() ;
  press.allocate_all() ;
  K_nn.allocate_all() ;
  K_np.allocate_all() ;
  K_pp.allocate_all() ;
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
  }
  
  if (l_min + nzet < nz) {
    nbar1.annule(l_min + nzet, nz - 1) ; 
    nbar2.annule(l_min + nzet, nz - 1) ; 
    ener.annule(l_min + nzet, nz - 1) ; 
    press.annule(l_min + nzet, nz - 1) ; 
    K_nn.annule(l_min + nzet, nz - 1) ; 
    K_np.annule(l_min + nzet, nz - 1) ; 
    K_pp.annule(l_min + nzet, nz - 1) ; 
  }

  int np0 = mp->get_mg()->get_np(0) ;
  int nt0 = mp->get_mg()->get_nt(0) ;
  for (int l=l_min; l<l_min+nzet; l++) {
    assert(mp->get_mg()->get_np(l) == np0) ;
    assert(mp->get_mg()->get_nt(l) == nt0) ; 
  }


  for (int k=0; k<np0; k++) {
    for (int j=0; j<nt0; j++) {

      bool inside = true ;  
      bool inside1 = true ; 
      bool inside2 = true ; 
     
      for (int l=l_min; l< l_min + nzet; l++) {
	for (int i=0; i<mp->get_mg()->get_nr(l); i++) {

	  /*cout << " k =" << k <<  " np0 =" << np0 << endl ;
	    cout << " j =" << j << " nt0 =" << nt0 << endl ;
	    cout << " l =" << l << " lmax =" << l_min + nzet << endl ;
	    cout << " i =" << i << " imax =" << mp->get_mg()->get_nr(l) << endl ;
	  */
	  
	  double xx, xx1, xx2; // xx1 = H1 = ln(mu1/m1)
	  xx1 = ent1(l,k,j,i) ; 
	  xx2 = ent2(l,k,j,i) ; 
	  xx = delta2(l,k,j,i) ;

	
	  
/*
 *  AVANT MODIFICATION
 *	  
	  // ZONE a zero fluide ----------------------------------------------------------------------------------------------
	  if ((exp(xx1) <= 1. ) && (exp(xx2) <= 1. ) ) { 
	    nbar1.set(l,k,j,i) = 0. ;
	    nbar2.set(l,k,j,i) = 0. ;
	    press.set(l,k,j,i) = 0. ;
	    ener.set(l,k,j,i) = 0.;
	    K_nn.set(l,k,j,i) = 0.;
	    K_np.set(l,k,j,i) = 0.;
 	    K_pp.set(l,k,j,i) = 0.;
	  }

	  // ZONE a deux fluides ou fluide seul -----------------------------------------------------------------------------
	  else {

	    double n1 = 0. ;
	    double n2 = 0.;
	    double p = 0. ;
	    
	    inside = (!nbar_ent_p(xx1, xx2, xx, n1, n2)) ;
	    nbar1.set(l,k,j,i) = n1 ;
	    nbar2.set(l,k,j,i) = n2 ;
	    
	    	//---si n1<0 -------------------------
// 		if (n1 <0. )
// 		{	
		// soit on met a zero
// 		n1 = 0.; 
		// soit on fait une interpolation lineaire (au choix)
		//int i_tri1 = 0;
		//interpol_linear(*mu1_N,  *n_1_N, mu1_adim, i_tri1 , n1) ;
// 		}
		//---si n1<0 -------------------------	
		//---si n2<0 -------------------------
// 		if (n2 <0. )
// 		{	
		// soit on met a zero
// 		n2 = 0.; 
		// soit on fait une interpolation lineaire (au choix)
		//int i_tri2 = 0;
		//interpol_linear(*mu2_P,  *n_2_P, mu2_adim, i_tri2 , n2) ;
// 		}
		//---si n2<0 -------------------------	
		
		//---si p<0 -------------------------
// 		if (p <0. )
// 		{	
		// soit on met a zero
// 		p = 0.; 
		// soit on fait une interpolation lineaire (au choix) // pas terrible
		//int i_trip = 0;
		//interpol_linear(*mu2_P,  *press_P, mu2_adim, i_trip , p) ;
// 		}
		//---si p<0 -------------------------
	   
	    p = press_ent_p(xx1, xx2, xx) ;
	    press.set(l,k,j,i) = p ;
	    double mu1 = m_1 * exp (xx1 );
	    double mu2 = m_2 * exp( xx2) ;
	    ener.set(l,k,j,i) = mu1 * n1 + mu2 * n2 - p ;//ener_ent_p(xx1, xx2, xx) ;
 	    K_nn.set(l,k,j,i) = get_K11(xx, xx1, xx2);   
 	    K_np.set(l,k,j,i) = get_K12(xx, xx1, xx2);   
 	    K_pp.set(l,k,j,i) = get_K22(xx, xx1, xx2);   
	    }
	    
*
*  	AVANT MODIFICATION
*/
	    
	    //---------------------------------------------------------------------------------------------------------------
	    // plus d'appel à plein de fonctions, tout est fait d'un coup! 
	    //-----------------------------------------------------------------------------------------------------------------

	    bool one_fluid = false; 
	    if (xx < delta_car_min) {
		cout << "Eos_bf_tabul::calcule_tout : delta2 < delta_car_min !" << endl ;
           	abort() ; 
	    }
	    if (xx > delta_car_max) {	
		cout << "Eos_bf_tabul::calcule_tout : delta2 > delta_car_max !" << endl ;
           	abort() ; 
	    }
	    if (m_1 * exp(xx1) > ent1_max) {
		cout << "Eos_bf_tabul::calcule_tout : ent1 > ent1_max !" << endl ;
           	abort() ; 
	    }
	    if (m_2 *exp(xx2) > ent2_max) {
		cout << "Eos_bf_tabul::calcule_tout : ent2 > ent2_max !" << endl ;
           	abort() ; 
	    }
	
	    double n1 = 0 ;
	    double n2 = 0 ;
	    double pressure = 0 ;
	    double alpha = 0 ;
	    double K11 = 0 ;
	    double K12 = 0 ; 
	    double K22 = 0 ;
	  	
	    double mu1 = m_1 * exp(xx1);
	    double mu2 = m_2 * exp(xx2);  
	    
	    if ( (exp(xx1) < 1.) && (exp(xx2) < 1.) ) {	
		n1 = 0. ;
		n2 = 0. ;
		pressure = 0.;
		alpha = 0 ;
		K11  = 0 ;
		K12 =  0 ;
		K22 = 0 ;
	    }
	    else {



	    double p_interpol ;
	    double dpsdent1_interpol ; 
	    double dpsdent2_interpol ;
	    double alpha_interpol ;
	
	   // cout << "DEBUG:" << xx << "  " << mu1 << "  " << mu2 << endl;
	    
	    interpol_mixed_3d_new(m_1, m_2, *delta_car, *logent1, *logent2, 
			*logp,  *dlpsdlent1, *dlpsdlent2,  *d2lpsdlent1dlent2,
			*dlpsddelta_car,  *d2lpsdlent1ddelta_car, *d2lpsdlent2ddelta_car, 
			*mu2_P, *n_p_P,   *press_P,
			*mu1_N, *n_n_N, *press_N,
			*delta_car_n0, *mu1_n0, *mu2_n0,
			*delta_car_p0, *mu1_p0, *mu2_p0, 
			xx, mu1, mu2, 
			p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha_interpol) ;

	   // en log : 
// 	   pressure = exp(p_interpol) ; 
// 	   n1 = dpsdent1_interpol * pressure ;
// 	   n2 = dpsdent2_interpol * pressure ;			
 	    n1 = dpsdent1_interpol ;
 	    n2 = dpsdent2_interpol ;
 	    pressure = p_interpol;
 	    alpha = alpha_interpol ;
	    
	  }
  
	  if (n1 < 0 ) {
		cout << "Eos_bf_tabul::nbar_ent_p : nbar1<0 !" << endl ;
//           	abort() ; 
//  		n1 = 0 ;
	  }
	  if (n2 < 0 ) {
		cout << "Eos_bf_tabul::nbar_ent_p : nbar2<0 !" << endl ;
// //           	abort() ; 
// 		n2 = 0 ;		
	  }
    //cout << "nbar1  " << nbar1 << endl ;
    //cout << "nbar2  " << nbar2 << endl ;

	  if (pressure < 0 ) {
		cout << "Eos_bf_tabul::press_ent_p : pressure < 0 !" << endl ;
//           	abort() ; 
//  		pressure = 0 ;
	  }
	
	
	  if (n1 >0.) {

	    K11 = mu1 / n1 - double(2) * alpha * ( 1. - xx) / ( n1 * n1 ) ; //OK

// 	    K11 = mu1 *n1  - double(2) * alpha * ( 1. - xx) ; // KNN----> KNN * Nn^2
	    
	  }
	  if (n2 >0.) {
 	      
 	    K22 = mu2 / n2 - double(2) * alpha * ( 1. - xx) / ( n2 * n2 ) ; //OK
  
// 	    K22 = mu2 * n2 - double(2) * alpha * ( 1. - xx)  ;	// KNN----> KNN * Nn^2	
	  
	  }
	  if ((n1 <= 0.) || (n2 <= 0.) ) { 
	  
	    K12 = 0. ;
	  
	  }
	  else {
 	  
	    K12 = double(2) * alpha * pow(1.-xx, 1.5)/ ( n1 * n2 );	//OK
// 	    K12 = double(2) * alpha * pow(1.-xx, 1.5) ;// KNN----> KNN * Nn^2
	    
	  }

	    nbar1.set(l,k,j,i) = n1 ;
	    nbar2.set(l,k,j,i) = n2 ;
	    press.set(l,k,j,i) = pressure ;
	    ener.set(l,k,j,i) = mu1 * n1 + mu2 * n2 - pressure ;//ener_ent_p(xx1, xx2, xx) ;
	    K_nn.set(l,k,j,i) = K11 ;   
	    K_np.set(l,k,j,i) = K12;  
	    K_pp.set(l,k,j,i) = K22 ;   
	  
	//--------------------------------------------------------------------------------------------------------------------
	    
	    
	    
	  //cout<< setprecision(16) ;
	  //cout << "DEBUG-------> delta2 = " << xx << " test de la fonction pow : pow(1-delta2,1.5) = " << pow(1.-xx,1.5) << endl;
	  	  
	  //Comparaison entre interpolation et calcul direct
	 // cout<< setprecision(16) ;
	 // cout << m_1 * exp (xx1 )* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "   " <<  m_2 * exp (xx2 )* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "  :" << endl;
	 // cout << "H1 = " << xx1 << " H2 =  " << xx2 << " n1 =  "<< nbar1(l,k,j,i)  << " n2 =  " << nbar2(l,k,j,i) <<   "  P =   " << press(l,k,j,i) <<    "  E =  " << ener(l,k,j,i) <<  endl; 
	  
	  //abort();
	  
	}
      }
      //abort();
      
    }
  }
  
}


// Baryon densities from enthalpies 
//---------------------------------

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
 	if (m_1 * exp(ent1) > ent1_max) {
		cout << "Eos_bf_tabul::nbar_ent_p : ent1 > ent1_max !" << endl ;
           	abort() ; 
	}
 	if (m_2 *exp(ent2) > ent2_max) {
		cout << "Eos_bf_tabul::nbar_ent_p : ent2 > ent2_max !" << endl ;
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
	
	  interpol_mixed_3d_new(m_1, m_2, *delta_car, *logent1, *logent2, 
			*logp,  *dlpsdlent1, *dlpsdlent2,  *d2lpsdlent1dlent2,
			*dlpsddelta_car,  *d2lpsdlent1ddelta_car, *d2lpsdlent2ddelta_car, 
			*mu2_P, *n_p_P,   *press_P,
			*mu1_N, *n_n_N, *press_N,
			*delta_car_n0, *mu1_n0, *mu2_n0,
			*delta_car_p0, *mu1_p0, *mu2_p0, 
			delta2, mu1, mu2, 
			p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;

	  nbar1 = dpsdent1_interpol ;
	  nbar2 = dpsdent2_interpol ;
	
	}
  
    if (nbar1 < 0 ) {
		cout << "Eos_bf_tabul::nbar_ent_p : nbar1<0 !" << endl ;
//           	abort() ; 
		nbar1 = 0 ;
	}
     if (nbar2 < 0 ) {
		cout << "Eos_bf_tabul::nbar_ent_p : nbar2<0 !" << endl ;
//           	abort() ; 
		nbar2 = 0 ;		
	}
    //cout << "nbar1  " << nbar1 << endl ;
    //cout << "nbar2  " << nbar2 << endl ;
    one_fluid = ((nbar1 <= 0.)||(nbar2 <= 0.)) ;
    
    return one_fluid; 
}
   
// One fluid sub-EOSs
//-------------------

//PLUS APPELE ACTUELLEMENT!!!!!!!!!!!!!!!!!!!
double Eos_bf_tabul::nbar_ent_p1(const double ent1) const {
	//cout << " appel a nbar_ent_p1" << endl;
	
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

}

//PLUS APPELE ACTUELLEMENT!!!!!!!!!!!!!!!!!!!
 double Eos_bf_tabul::nbar_ent_p2(const double ent2) const {
 //cout << " appel a nbar_ent_p2" << endl;
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

}

 // Energy density from baryonic densities 
//-----------------------------------------
double Eos_bf_tabul::ener_nbar_p(const double nbar1, const double nbar2, 
		   const double delta2) const{

  c_est_pas_fait("Eos_bf_tabul::ener_nbar_p" ) ;

  return nbar1 + nbar2 + delta2;

 }

// Pressure from baryonic densities 
//----------------------------------
double Eos_bf_tabul::press_nbar_p(const double nbar1, const double nbar2, 
		    const double delta2) const {

  c_est_pas_fait("Eos_bf_tabul::press_nbar_p" ) ;

    return nbar1 + nbar2 + delta2;

} 
     
	
// Pressure from baryonic log-enthalpies
//--------------------------------------
 double Eos_bf_tabul::press_ent_p(const double ent1, const double ent2, const double delta2) const {
	
	if (delta2 < delta_car_min) {
	cout << "Eos_bf_tabul::press_ent_p : delta2 < delta_car_min !" << endl ;
           	abort() ; 
	}
 	if (delta2 > delta_car_max) {
	cout << "Eos_bf_tabul::press_ent_p : delta2 > delta_car_max !" << endl ;
           	abort() ; 
	}
 	if (m_1 * exp(ent1) > ent1_max) {
	cout << "Eos_bf_tabul::press_ent_p : ent1 > ent1_max !" << endl ;
           	abort() ; 
	}
 	if (m_2 * exp(ent2) > ent2_max) {
	cout << "Eos_bf_tabul::press_ent_p : ent2 > ent2_max !" << endl ;
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
	
	  interpol_mixed_3d_new(m_1, m_2, *delta_car, *logent1, *logent2, 
			*logp,  *dlpsdlent1, *dlpsdlent2,  *d2lpsdlent1dlent2,
			*dlpsddelta_car,  *d2lpsdlent1ddelta_car, *d2lpsdlent2ddelta_car, 
			*mu2_P, *n_p_P,   *press_P,
			*mu1_N, *n_n_N, *press_N,
			*delta_car_n0, *mu1_n0, *mu2_n0,
			*delta_car_p0, *mu1_p0, *mu2_p0, 
			delta2, mu1, mu2, 
			p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;

 	pressure = p_interpol;

	}
 
	if (pressure < 0 ) {
		cout << "Eos_bf_tabul::press_ent_p : pressure < 0 !" << endl ;
//           	abort() ; 
		pressure = 0 ;
	}
	return pressure ;
 }

// PLUS APPELE MAINTEANTN
 double Eos_bf_tabul::press_ent_p1(const double ent1) const {
 	 
	double pressure ;
	if (m_1 * exp(ent1) > ent1_max) {
	cout << "Eos_bf_tabul::press_ent_p : ent1 > ent1_max !" << endl ;
           	abort() ; 
	}

	else {
 	
	double pressN_interpol ; 
	double n_n_N_interpol ; 
	
    	double mu1 = m_1 * exp(ent1);
	int i =0;
	interpol_herm( *mu1_N, *press_N, *n_n_N,
		   mu1,  i, pressN_interpol , n_n_N_interpol) ;
 
 	pressure =  pressN_interpol ;
 
	}

	return pressure ;
 }

 //PLUS APPELE MAINTENANT
double Eos_bf_tabul::press_ent_p2(const double ent2) const {
	double pressure ;
	  
	if (m_2*exp(ent2) > ent2_max) {
	cout << "Eos_bf_tabul::press_ent_p : ent1 > ent1_max !" << endl ;
           	abort() ; 
	}
	else {

	double pressP_interpol ; 
	double n_p_P_interpol ; 
	double mu2 = m_2* exp(ent2);
	int i =0;
	interpol_herm( *mu2_P, *press_P, *n_p_P,
		   mu2,  i, pressP_interpol , n_p_P_interpol) ;
 
 	pressure =  pressP_interpol ; 
 
	}
	return pressure ;
}


// Energy from baryonic log - enthalpies
//--------------------------------------

 double Eos_bf_tabul::ener_ent_p(const double ent1, const double ent2, const double delta2) const {

   double energy= 0.; 
   double pressure = 0. ; 
   double nbar1= 0. ;
   double nbar2=0.;

if ( (exp(ent1) < 1.) && ( exp(ent2) < 1.)) {
	energy = 0. ;
}
else {
	/*Eos_bf_tabul::nbar_ent_p( ent1,  ent2, delta2, nbar1, nbar2);
	
	if ( (nbar1 > 0) && (nbar2 > 0) ) { 
	
	pressure = Eos_bf_tabul::press_ent_p(ent1,ent2,delta2);      
  	double mu_1 = m_1 * exp ( ent1 );
   	double mu_2 = m_2 * exp ( ent2 );
	energy = nbar1 * mu_1 + nbar2 * mu_2 - pressure;

	}
	else if  ( (nbar1 <= 0) && (nbar2 > 0) ) { 
	
	//nbar1 = 0.;
	//nbar2 = Eos_bf_tabul::nbar_ent_p2(ent2);
	pressure = Eos_bf_tabul::press_ent_p(ent1,ent2,delta2); 
	//pressure = Eos_bf_tabul::press_ent_p2(ent2); 
	double mu_2 = m_2 * exp ( ent2 );
	energy =   nbar2 * mu_2 - pressure;
	//cout << nbar1 << "   " << pressure <<    "    "  << energy << endl ;
	
	}
	else if  ((nbar2 <= 0) && (nbar1 > 0)) { 
	//nbar2 = 0.;
	//nbar1 = Eos_bf_tabul:: nbar_ent_p1(ent1);
	//pressure = Eos_bf_tabul::press_ent_p1(ent1); 
	pressure = Eos_bf_tabul::press_ent_p(ent1,ent2,delta2);  
 	double mu_1 = m_1 * exp (ent1);
	energy =   nbar1 * mu_1 - pressure;
	//cout << nbar1 << "   " << pressure <<    "    "  << energy << endl ;

	}
	else if  ((nbar2 <= 0) && (nbar1 <= 0)) { 
	
	energy =  0.;*/

	double mu1 = m_1 * exp(ent1);
	double mu2 = m_2 * exp(ent2);

        double p_interpol ;
	double dpsdent1_interpol ; 
   	double dpsdent2_interpol ;
	double alpha ;
	
	interpol_mixed_3d_new(m_1, m_2, *delta_car, *logent1, *logent2, 
			*logp,  *dlpsdlent1, *dlpsdlent2,  *d2lpsdlent1dlent2,
			*dlpsddelta_car,  *d2lpsdlent1ddelta_car, *d2lpsdlent2ddelta_car, 
			*mu2_P, *n_p_P,   *press_P,
			*mu1_N, *n_n_N, *press_N,
			*delta_car_n0, *mu1_n0, *mu2_n0,
			*delta_car_p0, *mu1_p0, *mu2_p0, 
	 delta2, mu1, mu2, p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;

	energy = mu1 * dpsdent1_interpol + mu2 * dpsdent2_interpol - p_interpol ;
 
	}
	
	if (energy < 0 ) {
		cout << "Eos_bf_tabul::ener_ent_p : energy < 0 !" << endl ;
//           	abort() ; 
	}
	return energy;

}


// Alpha  from baryonic log - enthalpies
//--------------------------------------- 

double Eos_bf_tabul::alpha_ent_p(const double ent1, const double ent2, const double delta2) const {

	if (delta2 < delta_car_min) {
		cout << "Eos_bf_tabul::alpha_ent_p : delta2 < delta_car_min !" << endl ;
           	abort() ; 
	}
 	if (delta2 > delta_car_max) {
		cout << "Eos_bf_tabul::alpha_ent_p : delta2 > delta_car_max !" << endl ;
           	abort() ; 
	}
	if (m_1 * exp(ent1) > ent1_max) {
		cout << "Eos_bf_tabul::alpha_ent_p : ent1 > ent1_max !" << endl ;
           	abort() ; 
	}
	if (m_2 * exp(ent2) > ent2_max) {  
		cout << "Eos_bf_tabul::alpha_ent_p : ent2 > ent2_max !" << endl ;
           	abort() ; 
	}


	bool test= false ;
	double alpha;
	
	if ((exp(ent1) <= 1.) && (exp(ent2) <= 1.) ) {
		alpha = 0. ;
		test = false ;
	}

	else {

	  double mu1 = m_1 * exp(ent1);
	  double mu2 = m_2 * exp(ent2);

	  double p_interpol ;
	  double dpsdent1_interpol ; 
	  double dpsdent2_interpol ;
//	  double alpha ;
	
	  interpol_mixed_3d_new(m_1, m_2, *delta_car, *logent1, *logent2, 
			*logp,  *dlpsdlent1, *dlpsdlent2,  *d2lpsdlent1dlent2,
			*dlpsddelta_car,  *d2lpsdlent1ddelta_car, *d2lpsdlent2ddelta_car, 
			*mu2_P, *n_p_P,   *press_P,
			*mu1_N, *n_n_N, *press_N,
			*delta_car_n0, *mu1_n0, *mu2_n0,
			*delta_car_p0, *mu1_p0, *mu2_p0, 
	  delta2, mu1, mu2, p_interpol, dpsdent1_interpol, dpsdent2_interpol, alpha) ;
//cout << dpsdent1_interpol << "  " <<  dpsdent2_interpol << "   " << alpha <<endl;
	  //cout << " alpha   " << alpha << endl ;
	}
	//alpha = -alpha ; // NON NON NON LA CONVERSION EN -alpha est déjà faite dans l'interpolation
    //cout << ent1 << "   " << ent2 << "    " << alpha << endl;
   return alpha;
	

 }



// Derivatives of energy
//----------------------

double Eos_bf_tabul::get_K11(const double delta2, const double ent1, const double ent2)  const  {

	double xx=0.;
	double mu_1 = m_1 * exp(ent1);
	double nbar1;
	double nbar2;

	if ((exp(ent1) <= 1.) && (exp(ent2) <= 1.) ){
		xx = 0. ;
	}
	else {
		Eos_bf_tabul::nbar_ent_p(ent1,ent2, delta2, nbar1, nbar2) ; 
 
		/*if (nbar1 <= 0.) { 
		xx = 0. ;
		}
		else if ( (nbar1 > 0.) && ( nbar2 > 0. )) {
		//double mu_1 = m_1 * exp ( ent1 );
		double alpha = Eos_bf_tabul::alpha_ent_p(ent1,ent2,delta2) ;
 		xx = mu_1 / nbar1 - double(2) * alpha / ( nbar1 * nbar1 ) * ( 1. - delta2) ;
		}
		else if ( (nbar1 > 0.) && ( nbar2 <= 0. )) {
		//double mu_1 = m_1 * exp ( ent1 );
		//nbar1 = Eos_bf_tabul::nbar_ent_p1(ent1) ; 
	
		if (nbar1 !=0.) {
		xx = mu_1 / nbar1 ;
		}
		else { xx = 0.; }
		}*/
		double alpha = Eos_bf_tabul::alpha_ent_p(ent1,ent2,delta2) ;
		if (nbar1 >0.) {
		  
		  xx = mu_1 / nbar1 - double(2) * alpha * ( 1. - delta2) / ( nbar1 * nbar1 ) ;
		//cout << alpha << "   " << mu_1 / nbar1  << "   " << double(2) * alpha / ( nbar1 * nbar1 ) * ( 1. - delta2) << endl;
		}
	//cout << " test : " << " nbar =  " << nbar1 << "      " << mu_1 / nbar1  << "     " << double(2) * alpha * ( 1. - delta2) / ( nbar1 * nbar1 ) << "    " << "  Knn= " <<
// 	    xx << endl;	
		//cout << "get_K11  : ent1 =" << ent1 << " ent2 = " << ent2 << " delta2 =  " << delta2 << " nbar1 =  " <<  nbar1 << " nbar2 =  " << nbar2 << " alpha =  " << alpha << " Knn = " << xx << endl;
	}
	


   return xx;

}

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

	/*if  (nbar2 <= 0.) { 
	xx = 0. ;
	}
	else if ( (nbar2 > 0.) && ( nbar1 > 0. )) {
	double alpha = Eos_bf_tabul::alpha_ent_p(ent1,ent2,delta2) ;
 	xx = mu_2 / nbar2 - double(2) * alpha / ( nbar2 * nbar2 ) * ( 1. - delta2);
	}
	else if ( (nbar2 > 0.) && ( nbar1 <= 0. )) {

	//double mu_2 = m_2 * exp ( ent2);
	//nbar2 = Eos_bf_tabul::nbar_ent_p2(ent2) ; 
	if (nbar2 !=0.) {
	xx = mu_2 / nbar2 ;
	}
	else { xx = 0. ;}
	}*/
	  double alpha = Eos_bf_tabul::alpha_ent_p(ent1,ent2,delta2) ;
	  if (nbar2 >0.) {

	    xx = mu_2 / nbar2 - double(2) * alpha * ( 1. - delta2) / ( nbar2 * nbar2 ) ;

	//cout << alpha << "   " << mu_2 / nbar2  << "   " << double(2) * alpha / ( nbar2 * nbar2 ) * ( 1. - delta2) << endl;
	  }
//   	    cout << " test : " << " nbar2 =  " << nbar2 << "      " << mu_2 / nbar2  << "     " << double(2) * alpha * ( 1. - delta2) / ( nbar2 * nbar2 ) << "    " << "  Kpp= " <<
// 	    xx << endl;  	  
	  //cout << "get_K22  : ent1 =" << ent1 << " ent2 = " << ent2 << " delta2 =  " << delta2 << " nbar1 =  " <<  nbar1 << " nbar2 =  " << nbar2 << " alpha =  " << alpha << " Kpp = " << xx << endl;
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
	//cout << alpha << "   " << double(2) * alpha / ( nbar1 * nbar2 ) * pow(1.-delta2, 1.5)<< endl;	
	  }
	  //cout << "get_K12  : ent1 =" << ent1 << " ent2 = " << ent2 << " delta2 =  " << delta2 << " nbar1 =  " <<  nbar1 << " nbar2 =  " << nbar2 << " alpha =  " << alpha << " Knp = " << xx << endl;
	}
   // cout << ent1 << "   " << ent2 << "   " << xx << endl;
   return xx;
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
