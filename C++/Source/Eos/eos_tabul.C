/*
 *  Methods of class Eos_tabul
 *
 *  (see file eos.h for documentation).
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.22  2021/05/31 11:31:23  g_servignat
 * Added csound_square_ent routine to calculate the sound speed from enthalpy to Eos_consistent and corrected error outputs
 *
 * Revision 1.21  2021/05/14 15:39:22  g_servignat
 * Added sound speed computation from enthalpy to Eos class and tabulated+polytropic derived classes
 *
 * Revision 1.20  2019/12/02 14:51:36  j_novak
 * Moved piecewise parabolic computation of dpdnb to a separate function.
 *
 * Revision 1.19  2019/03/28 13:41:02  j_novak
 * Improved managed of saved EoS (functions sauve and constructor form FILE*)
 *
 * Revision 1.18  2017/12/15 15:36:38  j_novak
 * Improvement of the MEos class. Implementation of automatic offset computation accross different EoSs/domains.
 *
 * Revision 1.17  2016/12/05 16:17:52  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.16  2015/08/04 14:41:29  j_novak
 * Back to previous version for Eos_CompOSE. Enthalpy-consistent EoS can be accessed using Eos_consistent class (derived from Eos_CompOSE).
 *
 * Revision 1.15  2015/01/27 14:22:38  j_novak
 * New methods in Eos_tabul to correct for EoS themro consistency (optional).
 *
 * Revision 1.14  2014/10/13 08:52:54  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2014/06/30 16:13:18  j_novak
 * New methods for reading directly from CompOSE files.
 *
 * Revision 1.12  2014/03/06 15:53:35  j_novak
 * Eos_compstar is now Eos_compOSE. Eos_tabul uses strings and contains informations about authors.
 *
 * Revision 1.11  2012/11/09 10:32:44  m_bejger
 * Implementing der_ener_ent_p
 *
 * Revision 1.10  2010/02/16 11:14:50  j_novak
 * More verbose opeining of the file.
 *
 * Revision 1.9  2010/02/02 13:22:16  j_novak
 * New class Eos_Compstar.
 *
 * Revision 1.8  2004/03/25 10:29:02  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.7  2003/11/25 13:42:50  m_bejger
 * read_table written in more ordered way
 *
 * Revision 1.6  2003/11/21 16:11:16  m_bejger
 * added the computation dlnp/dlnn_b, dlnn/dlnH
 *
 * Revision 1.5  2003/05/30 07:50:06  e_gourgoulhon
 *
 * Reformulate the computation of der_nbar_ent
 * Added computation of der_press_ent_p.
 *
 * Revision 1.4  2003/05/15 09:42:47  e_gourgoulhon
 *
 * Now computes d ln / dln H.
 *
 * Revision 1.3  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/04/09 14:32:15  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2001/09/13  13:35:49  eric
 * Suppression des affichages dans read_table().
 *
 * Revision 2.2  2001/02/07  09:48:05  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 *
 * Revision 2.1  2000/11/23  00:10:20  eric
 * Enthalpie minimale fixee a 1e-14.
 *
 * Revision 2.0  2000/11/22  19:31:30  eric
 * *** empty log message ***
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
#include "eos.h"
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"


namespace Lorene {
  void interpol_herm(const Tbl& , const Tbl&, const Tbl&, double, int&,
		     double&, double& ) ;
  
  void interpol_linear(const Tbl&, const Tbl&, double, int&, double&) ;
  
  void interpol_quad(const Tbl&, const Tbl&, double, int&, double&) ;
  
  void compute_derivative(const Tbl&, const Tbl&, Tbl&) ;
  

    
			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
Eos_tabul::Eos_tabul(const char* name_i, const char* table,
		     const char* path) : Eos(name_i)
{	
	tablename = path ;
	tablename += "/" ;
	tablename += table ;
	
	read_table() ; 	
}

// Standard constructor with full filename
// ---------------------------------------
  Eos_tabul::Eos_tabul(const char* name_i, const char* file_name) 
    : Eos(name_i) {	

	tablename = file_name ;
	
	read_table() ; 	

}


// Constructor from binary file
// ----------------------------
  Eos_tabul::Eos_tabul(FILE* fich) : Eos(fich) {

  char tmp_string[160] ;
  fread(tmp_string, sizeof(char), 160, fich) ;
  tablename = tmp_string ;

  read_table() ;

}



// Constructor from a formatted file
// ---------------------------------
  Eos_tabul::Eos_tabul(ifstream& fich, const char* table) : 
    Eos(fich) {

  fich >> tablename ;
  tablename += "/" ;
  tablename += table ;
	
  read_table() ; 	

}

  Eos_tabul::Eos_tabul(ifstream& fich) : Eos(fich) 
{
  fich >> tablename ;
	
  read_table() ; 	
}

// Standard constructor with a name
// ---------------------------------
  Eos_tabul::Eos_tabul(const char* name_i) : Eos(name_i), 
					     logh(0x0), logp(0x0), dlpsdlh(0x0), 
					     lognb(0x0), dlpsdlnb(0x0), c_sound2(0x0)
{}


			//--------------//
			//  Destructor  //
			//--------------//

Eos_tabul::~Eos_tabul(){
  if (logh != 0x0) delete logh ;
  if (logp != 0x0) delete logp ;
  if (dlpsdlh != 0x0) delete dlpsdlh ;
  if (lognb != 0x0) delete lognb ;
  if (dlpsdlnb != 0x0) delete dlpsdlnb ;
  if (c_sound2 != 0x0) delete c_sound2 ;
}

			//------------//
			//  Outputs   //
			//------------//

void Eos_tabul::sauve(FILE* fich) const {

  Eos::sauve(fich) ;
  
  char tmp_string[160] ;
  strcpy(tmp_string, tablename.c_str()) ;
  fwrite(tmp_string, sizeof(char), 160, fich) ;		
}

			//------------------------//
			//  Reading of the table  //
			//------------------------//
			
void Eos_tabul::read_table() {

  using namespace Unites ;

  ifstream is_meos("meos_is_being_built.d") ;
    	
  ifstream fich(tablename.data()) ;

  if (!fich) {
    cerr << "Eos_tabul::read_table(): " << endl ;
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
    cerr << "Eos_tabul::read_table(): " << endl ;
    cerr << "Wrong value for the number of lines!" << endl ;
    cerr << "nbp = " << nbp << endl ;
    cerr << "Aborting..." << endl ;
    abort() ;
  }
  
  for (int i=0; i<3; i++) {		//  jump over the table
    fich.ignore(1000, '\n') ;             // description
  }                                      
  
  press = new double[nbp] ;
  nb    = new double[nbp] ;
  ro    = new double[nbp] ; 
  
  Tbl press_tbl(nbp) ; press_tbl.set_etat_qcq() ;
  Tbl nb_tbl(nbp)    ; nb_tbl.set_etat_qcq()    ;
  Tbl ro_tbl(nbp)    ; ro_tbl.set_etat_qcq()    ;
  
  logh = new Tbl(nbp) ;
  logp = new Tbl(nbp) ;
  dlpsdlh = new Tbl(nbp) ;
  lognb = new Tbl(nbp) ;
  dlpsdlnb = new Tbl(nbp) ;
  
  logh->set_etat_qcq() ;
  logp->set_etat_qcq() ;
  dlpsdlh->set_etat_qcq() ;
  lognb->set_etat_qcq() ;
  dlpsdlnb->set_etat_qcq() ;   	
  
  double rhonuc_cgs = rhonuc_si * 1e-3 ;
  double c2_cgs = c_si * c_si * 1e4 ;    	
  
  int no ;
  double nb_fm3, rho_cgs, p_cgs ;
  
  for (int i=0; i<nbp; i++) {
    
    fich >> no ;
    fich >> nb_fm3 ;
    fich >> rho_cgs ;
    fich >> p_cgs ; fich.ignore(1000,'\n') ;    		
    if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
      cout << "Eos_tabul::read_table(): " << endl ;
      cout << "Negative value in table!" << endl ;
      cout << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	", p = " << p_cgs << endl ;
      cout << "Aborting..." << endl ;
      abort() ;
    }
    
    press[i] = p_cgs / c2_cgs ;
    nb[i]    = nb_fm3 ;
    ro[i]    = rho_cgs ; 
    
    press_tbl.set(i) = press[i] ;
    nb_tbl.set(i) = nb[i] ;
    ro_tbl.set(i) = ro[i] ;
  }

  double ww = 0. ;
  bool store_offset = false ;
  bool compute_offset = true ;
  if (is_meos) {
    string temp ;
    string ok("ok") ;
    is_meos >> temp ;
    if (temp.find(ok) != string::npos) {
      is_meos >> ww ;
      compute_offset = false ;
    }
    else {
      is_meos.close() ;
      store_offset = true ;
    }
  }

  for (int i=0; i<nbp; i++) {
    double h = log( (ro[i] + press[i]) /
		    (10 * nb[i] * rhonuc_cgs) ) ;
    
    if ((i==0) && compute_offset) { ww = h ; }
    h = h - ww + 1.e-14 ;    		
    
    logh->set(i) = log10( h ) ;
    logp->set(i) = log10( press[i] / rhonuc_cgs ) ;
    dlpsdlh->set(i) = h * (ro[i] + press[i]) / press[i] ;
    lognb->set(i) = log10(nb[i]) ;
  } 

  if (store_offset == true) {
    ofstream write_meos("meos_is_being_built.d", ios_base::app) ;
    write_meos << "ok" << endl ;
    write_meos << setprecision(16) << ww << endl ;
  }

  Tbl logn0 = *lognb * log(10.) ;
  Tbl logp0 = ((*logp) + log10(rhonuc_cgs)) * log(10.) ;
  compute_derivative(logn0, logp0, *dlpsdlnb) ;

  hmin = pow( double(10), (*logh)(0) ) ;
  hmax = pow( double(10), (*logh)(nbp-1) ) ;
	
  Tbl tmp(nbp) ; tmp.set_etat_qcq() ;
  compute_derivative(ro_tbl,press_tbl,tmp) ;
  c_sound2 = new Tbl(tmp) ;// c_s^2 = dp/de
  // // Smoothing tmp
  // Tbl y(nbp) ; y.set_etat_qcq() ;
  // y.set(0) = tmp(0) ;
  // double RC = 1e-1 ; // Space constant, 1/2\piRC is cutoff space freq
  // double dx = (hmax-hmin)/nbp   ;
  // double alpha = dx/(RC+dx) ;
  // for (int i=1 ; i<5*nbp/8 ; i++){
  //   y.set(i) = y(i-1) + alpha*(tmp(i)-y(i-1)) ;
  // }
  // for (int i = 5*nbp/8 ; i<nbp ; i++) {
  //   y.set(i) = tmp(i) ;
  // }
  // c_sound2 = new Tbl(y) ;// c_s^2 = dp/de2
  
  cout << "hmin, hmax : " << hmin << "  " << hmax << endl ;

  fich.close();
  
  delete [] press ; 
  delete [] nb ; 
  delete [] ro ; 
  
   	
}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_tabul::nbar_ent_p(double ent, const Param* ) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::nbar_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }
           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;

           double pp = pow(double(10), logp0) ;

           double resu = pp / ent * dlpsdlh0 * exp(-ent) ;
	   return resu ;
    }
    else{
      return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_tabul::ener_ent_p(double ent, const Param* ) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::ener_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }
           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;

           double pp = pow(double(10), logp0) ;

	   double resu = pp / ent * dlpsdlh0 - pp ;
	   return resu ;
    }
    else{
      return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Eos_tabul::press_ent_p(double ent, const Param* ) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::press_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;
	   return pow(double(10), logp0) ;
    }
    else{
      return 0 ;
    }
}

// dln(n)/ln(H) from enthalpy 
//---------------------------

double Eos_tabul::der_nbar_ent_p(double ent, const Param* ) const {

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_nbar_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

	   double zeta = der_press_ent_p(ent) / der_press_nbar_p(ent) ; 

	   return zeta ; 
	  
    }
    else 

    return 1./(der_press_nbar_p(hmin)-1.) ;  // to ensure continuity

}


// dln(e)/ln(H) from enthalpy 
//---------------------------

double Eos_tabul::der_ener_ent_p(double ent, const Param* ) const {

    if ( ent > hmin ) {
		
           if (ent > hmax) {
           	cout << "Eos_tabul::der_ener_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

		   return ( der_nbar_ent_p(ent) 
			     *( double(1.) + press_ent_p(ent)/ener_ent_p(ent) )) ; 

    } else return der_nbar_ent_p(hmin) ;
    
}


// dln(p)/ln(H) from enthalpy 
//---------------------------

double Eos_tabul::der_press_ent_p(double ent, const Param* ) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_press_ent_p : ent > hmax !" << endl ;
           	abort() ;
           }

           double logh0 = log10( ent ) ;
           double logp0 ;
           double dlpsdlh0 ;
           interpol_herm(*logh, *logp, *dlpsdlh, logh0, i_near, logp0,
           		 dlpsdlh0) ;

           return dlpsdlh0 ;

    }
    else{
        
        return 0 ; 
	// return der_press_ent_p(hmin) ; // to ensure continuity
    }
}


// dln(p)/dln(n) from enthalpy 
//---------------------------

double Eos_tabul::der_press_nbar_p(double ent, const Param*) const {

    static int i_near = logh->get_taille() / 2 ;

    if ( ent > hmin ) {
           if (ent > hmax) {
           	cout << "Eos_tabul::der_press_nbar_p : ent > hmax !" << endl ;
           	abort() ;
           }

           double logh0 = log10( ent ) ;
           double dlpsdlnb0 ;

           interpol_linear(*logh, *dlpsdlnb, logh0, i_near, dlpsdlnb0) ;
 
          return dlpsdlnb0 ;

    }
    else{
        
         return (*dlpsdlnb)(0) ; 
	 // return der_press_nbar_p(hmin) ; // to ensure continuity

    }
}

double Eos_tabul::csound_square_ent_p(double ent, const Param*) const {
      static int i_near = lognb->get_taille() / 2 ;

      if ( ent > hmin ) {
        if (ent > hmax) {
    cout << "Eos_tabul::csound_square_ent_p : ent>hmax !" << endl ;
    abort() ;
        }
        double logent0 = log10( ent ) ;
        double csound0 ;
        
        
        interpol_linear(*logh, *c_sound2, logent0, i_near, csound0) ;
      
        return csound0 ;
      }
      else
        {
    return (*c_sound2)(0) ; 
        }
	}
}
