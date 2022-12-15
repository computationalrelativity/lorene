/*
 *  Methods of class Ye_eos_tabul
 *
 *  (see file hoteos.h for documentation).
 *
 */

/*
 *   Copyright (c) 2022 Jerome Novak, Gael Servignat
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
  void interpol_linear_2d(const Tbl&, const Tbl&, const Tbl&, double, double, int&, int&, double&) ;
  void huntm(const Tbl& xx, double& x, int& i_low) ;
  Tbl extract_column(const Tbl&, const Tbl&, double) ;


			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

  // Standard constructor
  // --------------------	 
  Ye_eos_tabul::Ye_eos_tabul(const string& filename): 
    Hot_eos("Tabulated Y_e EoS"), tablename(filename), authors("Unknown"), 
    hmin(0.), hmax(1.), yemin(0.), yemax(1.)
  {
    set_arrays_0x0() ;
    read_table() ;
  }
  
  
  // Constructor from binary file
  // ----------------------------
  Ye_eos_tabul::Ye_eos_tabul(FILE* fich) : Hot_eos(fich) {
    
    const int nc = 160 ;
    char tmp_string[nc] ;
    size_t ret = fread(tmp_string, sizeof(char), nc, fich) ;
    if (int(ret) == nc)
      tablename = tmp_string ;
    else {
      cerr << "Ye_tabul: constructor from a binary file:" << endl ;
      cerr << "Problems in reading the table name." << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }
    set_arrays_0x0() ;
    read_table() ;  
  }
  
  // Constructor from a formatted file
  // ---------------------------------
  Ye_eos_tabul::Ye_eos_tabul(ifstream& fich) : 
    Hot_eos(fich) {
    
    fich >> tablename ;
    set_arrays_0x0() ;
    read_table() ; 	 
  }

  //Sets the arrays to the null pointer
  void Ye_eos_tabul::set_arrays_0x0() {
    hhh = 0x0 ;
    Y_e = 0x0 ;
    c_sound2 = 0x0 ;
    mu_l = 0x0 ;
    ppp = 0x0 ;
    dpdh = 0x0 ;
    dpdye = 0x0 ;
    d2p = 0x0 ;
    Sourcetbl = 0x0 ;
  }

			//--------------//
			//  Destructor  //
			//--------------//

  Ye_eos_tabul::~Ye_eos_tabul(){
    if (hhh != 0x0) delete hhh ;
    if (Y_e != 0x0) delete Y_e ;
    if (ppp != 0x0) delete ppp ;
    if (dpdh != 0x0) delete dpdh ;
    if (dpdye != 0x0) delete dpdye ;
    if (d2p != 0x0) delete d2p ;
    if (Sourcetbl != 0x0) delete Sourcetbl ;
  }

			//------------//
			//  Outputs   //
			//------------//

  void Ye_eos_tabul::sauve(FILE* fich) const {
    
    Hot_eos::sauve(fich) ;
    
    char tmp_string[160] ;
    strcpy(tmp_string, tablename.c_str()) ;
    fwrite(tmp_string, sizeof(char), 160, fich) ;		
  }

  ostream& Ye_eos_tabul::operator>>(ostream & ost) const {
    
    ost << "Non-beta equilibrium EOS of class Ye_eos_tabul (tabulated out of beta-equilibrium EoS) : "
	<< endl ; 
    ost << "Built from file " << tablename << endl ;
    ost << "Authors : " << authors << endl ;
    ost << "Number of points in file : " << hhh->get_dim(0) 
	<< " in enthalpy, and " << hhh->get_dim(1) 
	<< " in electronic fraction." << endl ;

    return ost ;
}

			//------------------------//
			//  Reading of the table  //
			//------------------------//
			
  void Ye_eos_tabul::read_table() {

    using namespace Unites ;
    
    ifstream fich(tablename.data()) ;
    
    if (!fich) {
      cerr << "Ye_eos_tabul::read_table(): " << endl ;
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
      cerr << "Ye_eos_tabul::read_table(): " << endl ;
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
    Y_e = new Tbl(nbp2, nbp1) ;
    mu_l = new Tbl(nbp2, nbp1) ;
    c_sound2 = new Tbl(nbp2, nbp1) ;
    dpdh = new Tbl(nbp2, nbp1) ;
    dpdye = new Tbl(nbp2, nbp1) ;
    d2p = new Tbl(nbp2, nbp1) ;
    Sourcetbl = new Tbl(nbp2, nbp1) ; 
    
    ppp->set_etat_qcq() ;
    hhh->set_etat_qcq() ;
    Y_e->set_etat_qcq() ;
    mu_l->set_etat_qcq() ;
    c_sound2->set_etat_qcq() ;
    dpdh->set_etat_qcq() ;
    dpdye->set_etat_qcq() ;
    d2p->set_etat_qcq() ;
    Sourcetbl->set_etat_qcq() ;
    
    double c2 = c_si * c_si ;
    double nb_fm3, rho_cgs, p_cgs, ent, ye, mul, der2, cs2, source_d ;	
    double ww = 0. ;
    int no;

    for (int j=0; j<nbp2; j++) {
      for (int i=0; i<nbp1; i++) {
	fich >> no >> nb_fm3>> rho_cgs >> p_cgs>> ent >> ye >> cs2 >> mul >> der2 >> source_d ;

	fich.ignore(1000,'\n') ;
	if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	  cerr << "Ye_eos_tabul::read_table(): " << endl ;
	  cerr << "Negative value in table!" << endl ;
	  cerr << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	    ", p = " << p_cgs << endl ;
	  cerr << "Aborting..." << endl ;
	  abort() ;
	}
      
	double psc2 = 0.1*p_cgs/c2 ; //  in kg m^-3
	double rho_si = rho_cgs*1000. ; // in kg m^-3
	
	double h_read = ent ;
  if ( (i==0) ) ww = h_read ;
	double h_new = h_read - ww + 1e-14;
	
	ppp->set(j, i) = psc2/rhonuc_si ; 
	hhh->set(j, i) = h_new ;
	Y_e->set(j, i) = ye ;
  c_sound2->set(j, i) = cs2 ;
  mu_l->set(j, i) = mul*mev_si*1e44/(rhonuc_si*c2) ;
	dpdh->set(j, i) = (rho_si + psc2)/rhonuc_si ;
	dpdye->set(j, i) = -mul*nb_fm3*mevpfm3 ; 
	d2p->set(j, i) = 0.1*der2/(c2*rhonuc_si) ;	// Bien revérifier l'expression !!
  Sourcetbl->set(j, i) = source_d*t_unit ; 
      }
    }

    hmin = (*hhh)(0, 0) ;
    hmax = (*hhh)(0, nbp1-1) ;
    
    yemin = (*Y_e)(0, 0) ;
    yemax = (*Y_e)(nbp2-1, 0) ;
    
    cout << "hmin: " << hmin << ", hmax: " << hmax << endl ;
    cout << "yemin: " << yemin << ", yemax: " << yemax << endl ;
    
    fich.close();
}

			//-------------------------------//
			//  The corresponding cold Eos   //
			//-------------------------------//

  const Eos& Ye_eos_tabul::new_cold_Eos() const {
    
    if (p_cold_eos == 0x0) {
      cerr << "Warning: Ye_eos_tabul::new_cold_Eos " <<
	"The corresponding cold EoS is likely not to work" << endl ;
      cout << "read from file:   "<< tablename.c_str() << endl;
      p_cold_eos = new Eos_CompOSE(tablename.c_str()) ;
    }
    
    return *p_cold_eos ;
  }
  
  



			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Ye_eos_tabul::operator==(const Hot_eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
      cout << "The second EOS is not of type Ye_eos_tabul !" << endl ; 
      resu = false ; 
    }

    return resu ; 
  }


  bool Ye_eos_tabul::operator!=(const Hot_eos& eos_i) const {
    return !(operator==(eos_i)) ;  
  }
  
  
			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy and electronic fraction
//------------------------------------------

double Ye_eos_tabul::nbar_Hs_p(double ent, double ye) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if (ye < yemin) ye = yemin ;

  if ( ent >= hmin ) {
    if (ent > hmax) {
      cout << "Ye_eos_tabul::nbar_Hs_p : ent > hmax !" << endl ;
      abort() ;
    }

    if (ye > yemax) {
      cerr << "Ye_eos_tabul::nbar_Hs_p : Y_e not in the tabulated interval !" 
	   << endl ;
      cerr << "Y_e = " << ye << ", yemin = " << yemin << ", yemax = " << yemax 
	   << endl ;
      abort() ;
    }

    double p_int, dpdye_int, dpdh_int ;

    interpol_herm_2d(*Y_e, *hhh, *ppp, *dpdye, *dpdh, *d2p, ye, ent, p_int,
		     dpdye_int, dpdh_int) ;

    double nbar_int = dpdh_int * exp(-ent) ;
    return nbar_int ;
  }
  else{
    return 0 ;
  }
}


// Energy density from enthalpy and electronic fraction
//-----------------------------------------

double Ye_eos_tabul::ener_Hs_p(double ent, double ye) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if (ye < yemin) ye = yemin ;

  if ( ent >= hmin ) {
    if (ent > hmax) {
      cout << "Ye_eos_tabul::ener_Hs_p : ent > hmax !" << endl ;
      abort() ;
    }

    if (ye > yemax) {
      cerr << "Ye_eos_tabul::ener_Hs_p : Y_e not in the tabulated interval !" 
	   << endl ;
      cerr << "Y_e = " << ye << ", yemin = " << yemin << ", yemax = " << yemax 
	   << endl ;
      abort() ;
    }

    double p_int, dpdye_int, dpdh_int ;
    interpol_herm_2d(*Y_e, *hhh, *ppp, *dpdye, *dpdh, *d2p, ye, ent, p_int,
		     dpdye_int, dpdh_int) ;


    double f_int = - p_int + dpdh_int;

    return f_int ;
  }
  else{
    return 0 ;
  }
}

double Ye_eos_tabul::press_Hs_p(double ent, double ye) const {

  if ((ent > hmin - 1.e-12) && (ent < hmin))
    ent = hmin ;

  if (ye < yemin) ye = yemin ;

  if ( ent >= hmin ) {
    if (ent > hmax) {
      cout << "Ye_eos_tabul::press_Hs_p : ent > hmax !" << endl ;
      abort() ;
    }

    if (ye > yemax) {
      cerr << "Ye_eos_tabul::press_Hs_p : Y_e not in the tabulated interval !" 
	   << endl ;
      cerr << "Y_e = " << ye << ", yemin = " << yemin << ", yemax = " << yemax 
	   << endl ;
      abort() ;
    }

    double p_int, dpdye_int, dpdh_int ;
    interpol_herm_2d(*Y_e, *hhh, *ppp, *dpdye, *dpdh, *d2p, ye, ent, p_int,
		     dpdye_int, dpdh_int) ;

    return p_int ;
  }
  else{
    return 0 ;
  }
}

double Ye_eos_tabul::csound_square_Hs_p(double ent, double ye) const {
      static int i_near = hhh->get_dim(0) / 2 ;
      static int j_near = Y_e->get_dim(1) / 2 ;
      
      if ((ent > hmin - 1.e-12) && (ent < hmin))
        ent = hmin ;

      if (ye < yemin) ye = yemin ;
      
      Tbl ye_1D(Y_e->get_dim(1)) ; ye_1D.set_etat_qcq() ;
        for (int i=0 ; i<Y_e->get_dim(1) ; i++){
          ye_1D.set(i) = (*Y_e)(i,0) ;
      }
      Tbl enthalpy = extract_column(*hhh,ye_1D,ye) ;
      if ( ent > hmin ) {
        if (ent > hmax) {
    cout << "Eos_tabul::csound_square_Hs_p : ent>hmax !" << endl ;
    abort() ;
        }
        if (ye > yemax) {
          cerr << "Ye_eos_tabul::csound_square_Hs_p : Y_e not in the tabulated interval !" 
        << endl ;
          cerr << "Y_e = " << ye << ", yemin = " << yemin << ", yemax = " << yemax 
        << endl ;
          abort() ;
        }
        double csound0 ;
            
        
        interpol_linear_2d(enthalpy, ye_1D, *c_sound2, ent, ye, i_near, j_near, csound0) ;

      
        return csound0 ;
      }
      else
        {
          double csound0 ;
          interpol_linear_2d(enthalpy, ye_1D, *c_sound2, hmin, ye, i_near, j_near, csound0) ;
    return csound0 ; 
        }
	}

double Ye_eos_tabul::chi2_Hs_p(double, double) const {
    
    cerr << "Warning : (H,Y_e) EoS have no contribution from chi^2 ; Ye_eos_tabul::chi2_Hs_p :function not implemented." << endl;
    cerr << "Aborting ..." << endl;
    abort() ;
    
    return 0. ;
	}

double Ye_eos_tabul::mul_Hs_p(double ent, double ye) const {
      static int i_near = hhh->get_dim(0) / 2 ;
      static int j_near = Y_e->get_dim(1) / 2 ;
      
      if ((ent > hmin - 1.e-12) && (ent < hmin))
        ent = hmin ;

      if (ye < yemin) ye = yemin ;
      
      Tbl ye_1D(Y_e->get_dim(1)) ; ye_1D.set_etat_qcq() ;
      for (int i=0 ; i<Y_e->get_dim(1) ; i++){
        ye_1D.set(i) = (*Y_e)(i,0) ;
      }
      Tbl enthalpy = extract_column(*hhh,ye_1D,ye) ;
      
      if ( ent > hmin ) {
        if (ent > hmax) {
    cout << "Eos_tabul::mul_Hs_p : ent>hmax !" << endl ;
    abort() ;
        }
        
        if (ye > yemax) {
          cerr << "Ye_eos_tabul::mul_Hs_p : Y_e not in the tabulated interval !" 
        << endl ;
          cerr << "Y_e = " << ye << ", yemin = " << yemin << ", yemax = " << yemax 
        << endl ;
          abort() ;
        }
        
        double mul0 ;

        interpol_linear_2d(enthalpy, ye_1D, *mu_l, ent, ye, i_near, j_near, mul0) ;
      
        return mul0 ;
      }
      else
        {
          double mul0 ;
          interpol_linear_2d(enthalpy, ye_1D, *mu_l, hmin, ye, i_near, j_near, mul0) ;
    return mul0 ; 
        }
	}

  double Ye_eos_tabul::sigma_Hs_p(double ent, double ye) const {
      static int i_near = hhh->get_dim(0) / 2 ;
      static int j_near = Y_e->get_dim(1) / 2 ;
      
      if ((ent > hmin - 1.e-12) && (ent < hmin))
        ent = hmin ;

      if (ye < yemin) ye = yemin ;
      
      Tbl ye_1D(Y_e->get_dim(1)) ; ye_1D.set_etat_qcq() ;
      for (int i=0 ; i<Y_e->get_dim(1) ; i++){
        ye_1D.set(i) = (*Y_e)(i,0) ;
      }
      Tbl enthalpy = extract_column(*hhh,ye_1D,ye) ;
      
      if ( ent > hmin ) {
        if (ent > hmax) {
    cout << "Eos_tabul::sigma_Hs_p : ent>hmax !" << endl ;
    abort() ;
        }
        
        if (ye > yemax) {
          cerr << "Ye_eos_tabul::sigma_Hs_p : Y_e not in the tabulated interval !" 
        << endl ;
          cerr << "Y_e = " << ye << ", yemin = " << yemin << ", yemax = " << yemax 
        << endl ;
          abort() ;
        }
        
        double sigma0 ;

        interpol_linear_2d(enthalpy, ye_1D, *Sourcetbl, ent, ye, i_near, j_near, sigma0) ;
      
        return sigma0 ;
      }
      else
        {
          double sigma0 ;
          interpol_linear_2d(enthalpy, ye_1D, *Sourcetbl, hmin, ye, i_near, j_near, sigma0) ;
    return sigma0 ; 
        }
	}
  
  double Ye_eos_tabul::temp_Hs_p(double, double) const {
    
    cerr << "Warning : (H,Y_e) EoS does not use T as a parameter ; Ye_eos_tabul::temp_Hs_p :function not implemented." << endl;
    cerr << "Aborting ..." << endl;
    abort() ;
    
    return 0. ;
  }
} //End of namespace Lorene


