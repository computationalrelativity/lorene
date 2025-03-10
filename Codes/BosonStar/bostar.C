/*
 *  Code for computing an equilibrium configuration of a rotating boson star
 *
 *    (see file boson_star.h for documentation).
 *
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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
 * Revision 1.4  2016/12/05 16:18:24  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:55  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2012/12/03 15:28:03  c_some
 * First call to Boson_star::equilibrium
 *
 * Revision 1.1  2012/11/23 15:41:44  c_some
 * First version
 *
 *
 * $Header$
 *
 */

// C++ headers

// C headers

// Lorene headers
#include "boson_star.h"
#include "nbr_spx.h"


using namespace Lorene ;

int main() {
	
    // Identification of all the subroutines called by the code : 
    
    system("ident bostar > identif.d") ; 

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    int kkk, mer_max,  mermax_poisson, graph, nz, nzadapt, nt, np, mer_triax ; 
    double mmm, rphi_c, iphi_c, precis, thres_adapt,  relax, relax_poisson, ampli_triax, 
	   precis_adapt ;  
    
    ifstream fpar("param.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file param.d ! " << endl ;
        abort() ;
    }
    
    fpar.ignore(1000,'\n') ;    // skip title
    fpar >> mmm ; fpar.ignore(1000,'\n') ;
    fpar >> kkk ; fpar.ignore(1000,'\n') ;
    fpar >> rphi_c ; fpar.ignore(1000,'\n') ;
    fpar >> iphi_c ; fpar.ignore(1000,'\n') ;
    fpar.ignore(1000,'\n') ;	// skip title
    fpar >> mer_max ; fpar.ignore(1000,'\n') ;
    fpar >> precis ; fpar.ignore(1000,'\n') ;
    fpar >> thres_adapt ; fpar.ignore(1000,'\n') ;
    fpar >> mer_triax ; fpar.ignore(1000,'\n') ;
    fpar >> ampli_triax ; fpar.ignore(1000,'\n') ;
    fpar >> relax ; fpar.ignore(1000,'\n') ;
    fpar >> mermax_poisson ; fpar.ignore(1000,'\n') ;
    fpar >> relax_poisson ; fpar.ignore(1000,'\n') ;
    fpar >> precis_adapt ; fpar.ignore(1000,'\n') ;
    fpar >> graph ; fpar.ignore(1000,'\n') ;
    fpar.ignore(1000,'\n') ; // skip title
    fpar >> nz ; fpar.ignore(1000,'\n') ;
    fpar >> nzadapt; fpar.ignore(1000,'\n') ;
    fpar >> nt; fpar.ignore(1000,'\n') ;
    fpar >> np; fpar.ignore(1000,'\n') ;

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
     
    fpar.ignore(1000,'\n'); 	// skip title
    for (int l=0; l<nz; l++) {
		fpar >> nr[l]; 
		fpar >> bornes[l]; fpar.ignore(1000,'\n') ;
		np_tab[l] = np ; 
		nt_tab[l] = nt ; 
    }
    bornes[nz] = __infinity ;

    fpar.close();

    //-----------------------------------------------------------------------
    //		Construction of the multi-grid and the mapping
    //-----------------------------------------------------------------------


    // Type of r sampling :
    int* type_r = new int[nz];
    type_r[0] = RARE ; 
    for (int l=1; l<nz-1; l++) {
		type_r[l] = FIN ; 
    }
    type_r[nz-1] = UNSURR ; 
    
    // Type of sampling in theta and phi :
    int type_t = SYM ; 
    int type_p = SYM ; 
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_et mp(mg, bornes) ;
   
    // Cleaning
    // --------

    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] type_r ; 
    delete [] bornes ; 
       
    //-----------------------------------------------------------------------
    //		Construction of the boson star
    //-----------------------------------------------------------------------

	Boson_star star(mp, mmm, kkk) ; 
	
	//-----------------------------------------------------------------------
    //		Initialization of Phi
    //-----------------------------------------------------------------------

    const Coord& r = mp.r ;
    double ray0 = mp.val_r(0, 1., 0., 0.) ;  
    Scalar rphi0(mp) ; 
    rphi0 = rphi_c * exp(- r*r / (ray0*ray0) ) ; 
 	rphi0.std_spectral_base() ;  // sets the standard bases for spectral expansions
 
    Scalar iphi0(mp) ; 
    iphi0 = iphi_c * exp(- r*r / (ray0*ray0) ) ; 
 	iphi0.std_spectral_base() ;  // sets the standard bases for spectral expansions

    star.set_rphi() = rphi0 ;  
    star.set_iphi() = iphi0 ;  
    
    // Initialization of the energy-momentum tensor
    star.update_ener_mom() ; 

    cout << endl << "Initial star : " 
	 << endl << "==========   " << endl ;
	cout <<  star << endl ; 	


    //-----------------------------------------------------------------------
    //		Computation of the rotating equilibrium
    //-----------------------------------------------------------------------

    Itbl icontrol(8) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
//##    icontrol.set(1) = mer_rot ; 
//    icontrol.set(2) = mer_change_omega ; 
//    icontrol.set(3) = mer_fix_omega ; 
//##    icontrol.set(4) = mer_mass ; 
    icontrol.set(5) = mermax_poisson ; 
    icontrol.set(6) = mer_triax ; 
//##    icontrol.set(7) = delta_mer_kep ; 
    
    Tbl control(7) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
//##    control.set(1) = omega_ini ; 
    control.set(2) = relax ; 
    control.set(3) = relax_poisson ; 
//##    control.set(4) = thres_adapt ; 
    control.set(5) = ampli_triax ; 
    control.set(6) = precis_adapt ; 

    Tbl diff(8) ;     

    Tbl phi_limit(1) ;
    phi_limit.set_etat_qcq() ;
    phi_limit.set(0) = 1e-3*rphi_c  ; 	// Phi at the stellar "surface"

	star.equilibrium(rphi_c, iphi_c, nzadapt, phi_limit, icontrol, control, diff) ;


    cout << endl << "Final star : " 
	 << endl << "==========   " << endl ;
	cout <<  star << endl ; 	

}
