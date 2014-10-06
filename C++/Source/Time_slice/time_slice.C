/*
 *  Methods of class Time_slice
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon, Jose Luis Jaramillo & Jerome Novak
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

char time_slice_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.15  2014/10/06 15:13:21  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.14  2008/12/02 15:02:21  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.13  2004/05/31 09:08:18  e_gourgoulhon
 * Method sauve and constructor from binary file are now operational.
 *
 * Revision 1.12  2004/05/27 15:25:04  e_gourgoulhon
 * Added constructors from binary file, as well as corresponding
 * functions sauve and save.
 *
 * Revision 1.11  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.10  2004/05/10 09:08:34  e_gourgoulhon
 * Added "adm_mass_evol.downdate(jtime)" in method del_deriv.
 * Added printing of ADM mass in operator>>(ostream&).
 *
 * Revision 1.9  2004/05/09 20:57:34  e_gourgoulhon
 * Added data member adm_mass_evol.
 *
 * Revision 1.8  2004/05/05 14:26:25  e_gourgoulhon
 * Minor modif. in operator>>(ostream& ).
 *
 * Revision 1.7  2004/04/07 07:58:21  e_gourgoulhon
 * Constructor as Minkowski slice: added call to std_spectral_base()
 * after setting the lapse to 1.
 *
 * Revision 1.6  2004/04/01 16:09:02  j_novak
 * Trace of K_ij is now member of Time_slice (it was member of Time_slice_conf).
 * Added new methods for checking 3+1 Einstein equations (preliminary).
 *
 * Revision 1.5  2004/03/29 11:59:23  e_gourgoulhon
 * Added operator>>.
 *
 * Revision 1.4  2004/03/28 21:29:45  e_gourgoulhon
 * Evolution_std's renamed with suffix "_evol"
 * Method gam() modified
 * Added special constructor for derived classes.
 *
 * Revision 1.3  2004/03/26 13:33:02  j_novak
 * New methods for accessing/updating members (nn(), beta(), gam_uu(), k_uu(), ...)
 *
 * Revision 1.2  2004/03/26 08:22:56  e_gourgoulhon
 * Modifications to take into account the new setting of class
 * Evolution.
 *
 * Revision 1.1  2004/03/24 14:57:17  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
#include <cassert>

// Lorene headers
#include "time_slice.h"
#include "utilitaires.h"



			    //--------------//
			    // Constructors //
			    //--------------//


// Standard constructor (Hamiltonian-like)
// ---------------------------------------
Time_slice::Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor& kk_in,
               int depth_in) 
               : depth(depth_in),
		 scheme_order(depth_in-1),
                 jtime(0),
                 the_time(0., depth_in),
                 gam_dd_evol(depth_in),
                 gam_uu_evol(depth_in),
                 k_dd_evol(depth_in),
                 k_uu_evol(depth_in),
                 n_evol(lapse_in, depth_in),
                 beta_evol(shift_in, depth_in),
		 trk_evol(depth_in),
                 adm_mass_evol() {
                                  
    set_der_0x0() ; 

    double time_init = the_time[jtime] ; 

    p_gamma = new Metric(gamma_in) ; 
                    
    if (gamma_in.get_index_type(0) == COV) {
        gam_dd_evol.update(gamma_in, jtime, time_init) ; 
    }
    else {
        gam_uu_evol.update(gamma_in, jtime, time_init) ; 
    }
                 
    if (kk_in.get_index_type(0) == COV) {
        k_dd_evol.update(kk_in, jtime, time_init) ; 
    }
    else {
        k_uu_evol.update(kk_in, jtime, time_init) ; 
    }
                 
    trk_evol.update( kk_in.trace(*p_gamma), jtime, the_time[jtime] ) ; 
    
}
                 
                 
// Standard constructor (Lagrangian-like)
// ---------------------------------------
Time_slice::Time_slice(const Scalar& lapse_in, const Vector& shift_in,
               const Evolution_std<Sym_tensor>& gamma_in) 
               : depth(gamma_in.get_size()), 
 		 scheme_order(gamma_in.get_size()-1),
		 jtime(0),
                 the_time(0., gamma_in.get_size()),
                 gam_dd_evol( gamma_in.get_size() ),
                 gam_uu_evol( gamma_in.get_size() ),
                 k_dd_evol( gamma_in.get_size() ),
                 k_uu_evol( gamma_in.get_size() ),
                 n_evol(lapse_in, gamma_in.get_size() ),
                 beta_evol(shift_in, gamma_in.get_size() ),
		 trk_evol(gamma_in.get_size() ),
                 adm_mass_evol() {

    cerr << 
    "Time_slice constuctor from evolution of gamma not implemented yet !\n" ;
    abort() ; 
                 
    set_der_0x0() ; 

}                 
                 
// Constructor as standard time slice of flat spacetime (Minkowski)
//-----------------------------------------------------------------               
Time_slice::Time_slice(const Map& mp, const Base_vect& triad, int depth_in)  
               : depth(depth_in),
		 scheme_order(depth_in-1),
                 jtime(0),
                 the_time(0., depth_in),
                 gam_dd_evol(depth_in),
                 gam_uu_evol(depth_in),
                 k_dd_evol(depth_in),
                 k_uu_evol(depth_in),
                 n_evol(depth_in),
                 beta_evol(depth_in),
		 trk_evol(depth_in),
                 adm_mass_evol() {
                 
    double time_init = the_time[jtime] ; 
    
    const Base_vect_spher* ptriad_s = 
        dynamic_cast<const Base_vect_spher*>(&triad) ;                  
    bool spher = (ptriad_s != 0x0) ; 
    
    if (spher) {
        gam_dd_evol.update( mp.flat_met_spher().cov(), jtime, time_init) ;  
    }         
    else {
        assert( dynamic_cast<const Base_vect_cart*>(&triad) != 0x0) ; 
        gam_dd_evol.update( mp.flat_met_cart().cov(), jtime, time_init) ;                           
    }


    // K_ij identically zero:
    Sym_tensor ktmp(mp, COV, triad) ;
    ktmp.set_etat_zero() ; 
    k_dd_evol.update(ktmp, jtime, time_init) ;  

    // Lapse identically one:
    Scalar tmp(mp) ; 
    tmp.set_etat_one() ; 
    tmp.std_spectral_base() ; 
    n_evol.update(tmp, jtime, time_init) ; 
    
    // shift identically zero:
    Vector btmp(mp, CON, triad) ;
    btmp.set_etat_zero() ; 
    beta_evol.update(btmp, jtime, time_init) ;  
    
    // trace(K) identically zero:
    tmp.set_etat_zero() ; 
    trk_evol.update(tmp, jtime, time_init) ; 
    
    set_der_0x0() ; 
}
    
// Constructor from binary file             
// ----------------------------

Time_slice::Time_slice(const Map& mp, const Base_vect& triad, FILE* fich, 
                       bool partial_read, int depth_in) 
               : depth(depth_in),
                 the_time(depth_in),
                 gam_dd_evol(depth_in),
                 gam_uu_evol(depth_in),
                 k_dd_evol(depth_in),
                 k_uu_evol(depth_in),
                 n_evol(depth_in),
                 beta_evol(depth_in),
		         trk_evol(depth_in),
                 adm_mass_evol() {

    // Reading various integer parameters
    // ----------------------------------
    
    int depth_file ; 
    fread_be(&depth_file, sizeof(int), 1, fich) ;
    if (depth_file != depth_in) {
        cout << 
        "Time_slice constructor from file: the depth read in file \n"
        << " is different from that given in the argument list : \n"
        << "   depth_file = " << depth_file 
        << " <-> depth_in " << depth_in << " !" << endl ; 
        abort() ;  
    }	
    fread_be(&scheme_order, sizeof(int), 1, fich) ;	
    fread_be(&jtime, sizeof(int), 1, fich) ;	
    
    // Reading the_time
    // ----------------
    int jmin = jtime - depth + 1 ; 
    int indicator ; 
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;	
        if (indicator == 1) {	
            double xx ; 
            fread_be(&xx, sizeof(double), 1, fich) ;	
            the_time.update(xx, j, xx) ; 
        }
    }

    // Reading of various fields
    // -------------------------
    
    // N
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;	
        if (indicator == 1) {
            Scalar nn_file(mp, *(mp.get_mg()), fich) ; 
            n_evol.update(nn_file, j, the_time[j]) ; 
        }
    }

    // beta
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;	
        if (indicator == 1) {
            Vector beta_file(mp, triad, fich) ; 
            beta_evol.update(beta_file, j, the_time[j]) ; 
        }
    }

    // Case of a full reading
    // -----------------------
    if (!partial_read) {
        cout << 
        "Time_slice constructor from file: the case of full reading\n"
        << " is not ready yet !" << endl ; 
        abort() ; 
    }

    set_der_0x0() ; 

} 


// Copy constructor
// ----------------
Time_slice::Time_slice(const Time_slice& tin) 
               : depth(tin.depth),
		 scheme_order(tin.scheme_order),
                 jtime(tin.jtime),
                 the_time(tin.the_time),
                 gam_dd_evol(tin.gam_dd_evol),
                 gam_uu_evol(tin.gam_uu_evol),
                 k_dd_evol(tin.k_dd_evol),
                 k_uu_evol(tin.k_uu_evol),
                 n_evol(tin.n_evol),
                 beta_evol(tin.beta_evol),
		 trk_evol(tin.trk_evol),
                 adm_mass_evol(tin.adm_mass_evol) {
                 
    set_der_0x0() ; 
}

// Special constructor for derived classes
//----------------------------------------               
Time_slice::Time_slice(int depth_in)  
               : depth(depth_in),
		 scheme_order(depth_in-1),
                 jtime(0),
                 the_time(0., depth_in),
                 gam_dd_evol(depth_in),
                 gam_uu_evol(depth_in),
                 k_dd_evol(depth_in),
                 k_uu_evol(depth_in),
                 n_evol(depth_in),
                 beta_evol(depth_in),
		 trk_evol(depth_in),
                 adm_mass_evol() {
                 
    set_der_0x0() ; 
}
                 



			    //--------------//
			    //  Destructor  //
			    //--------------//

Time_slice::~Time_slice(){

    Time_slice::del_deriv() ; 

}

                //---------------------//
                //  Memory management  //
                //---------------------//

void Time_slice::del_deriv() const {

    if (p_gamma != 0x0) delete p_gamma ; 
    
    set_der_0x0() ;
    
    adm_mass_evol.downdate(jtime) ; 
}


void Time_slice::set_der_0x0() const {

    p_gamma = 0x0 ; 
    
}


                //-----------------------//
                // Mutators / assignment //
                //-----------------------//

void Time_slice::operator=(const Time_slice& tin) {

    depth = tin.depth;
    scheme_order = tin.scheme_order ;
    jtime = tin.jtime;
    the_time = tin.the_time;
    gam_dd_evol = tin.gam_dd_evol;
    gam_uu_evol = tin.gam_uu_evol;
    k_dd_evol = tin.k_dd_evol;
    k_uu_evol = tin.k_uu_evol;
    n_evol = tin.n_evol;
    beta_evol = tin.beta_evol; 
    trk_evol = tin.trk_evol ;
    
    del_deriv() ; 
    
}


                //------------------//
                //      output      //
                //------------------//

ostream& Time_slice::operator>>(ostream& flux) const {

    flux << "\n------------------------------------------------------------\n" 
         << "Lorene class : " << typeid(*this).name() << '\n' ; 
    flux << "Number of stored slices : " << depth  
        << "     order of time scheme : " << scheme_order << '\n' 
        << "Time label t = " << the_time[jtime]  
        << "               index of time step j = " << jtime << '\n' << '\n' ; 
    if (adm_mass_evol.is_known(jtime)) {
        flux << "ADM mass : " << adm_mass() << endl ;
    }
         
    flux << "Max. of absolute values of the various fields in each domain: \n" ;
    if (gam_dd_evol.is_known(jtime)) {
        maxabs( gam_dd_evol[jtime], "gam_{ij}", flux) ;
    }
    if (gam_uu_evol.is_known(jtime)) {
        maxabs( gam_uu_evol[jtime], "gam^{ij}", flux) ;
    }
    if (k_dd_evol.is_known(jtime)) {
        maxabs( k_dd_evol[jtime], "K_{ij}", flux) ;
    }
    if (k_uu_evol.is_known(jtime)) {
        maxabs( k_uu_evol[jtime], "K^{ij}", flux) ;
    }
    if (n_evol.is_known(jtime)) {
        maxabs( n_evol[jtime], "N", flux) ;
    }
    if (beta_evol.is_known(jtime)) {
        maxabs( beta_evol[jtime], "beta^i", flux) ;
    }
    if (trk_evol.is_known(jtime)) {
        maxabs( trk_evol[jtime], "tr K", flux) ;
    }

    if (p_gamma != 0x0) flux << "Metric gamma is up to date" << endl ; 
    
    return flux ; 

}


ostream& operator<<(ostream& flux, const Time_slice& sigma) {

    sigma >> flux ;     
    return flux ; 

}


void Time_slice::save(const char* rootname) const {

    // Opening of file 
    // ---------------
    char* filename = new char[ strlen(rootname)+10 ] ; 
    strcpy(filename, rootname) ; 
    char nomj[7] ; 
    sprintf(nomj, "%06d", jtime) ; 
    strcat(filename, nomj) ; 
    strcat(filename, ".d") ; 
        
    FILE* fich = fopen(filename, "w") ; 
    if (fich == 0x0) {
    	cout << "Problem in opening file " << filename << " ! " << endl ; 
	    perror(" reason") ; 
	    abort() ; 
    }

    // Write grid, mapping, triad and depth
    // ------------------------------------
    const Map& map = nn().get_mp() ;
    const Mg3d& mgrid = *(map.get_mg()) ; 
    const Base_vect& triad = *(beta().get_triad()) ; 
    
    mgrid.sauve(fich) ; 
    map.sauve(fich) ; 
    triad.sauve(fich) ;  
   
    fwrite_be(&depth, sizeof(int), 1, fich) ;	

    // Write all binary data by means of virtual function sauve
    // --------------------------------------------------------
    bool partial_save = false ;
    sauve(fich, partial_save) ; 
    
    // Close the file
    // --------------
    
    fclose(fich) ; 
        
    delete [] filename ;        
        
}



void Time_slice::sauve(FILE* fich, bool partial_save) const {

    // Writing various integer parameters
    // ----------------------------------
    
    fwrite_be(&depth, sizeof(int), 1, fich) ;	
    fwrite_be(&scheme_order, sizeof(int), 1, fich) ;	
    fwrite_be(&jtime, sizeof(int), 1, fich) ;	
    
    // Writing the_time
    // ----------------
    int jmin = jtime - depth + 1 ; 
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (the_time.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) {	
            double xx = the_time[j] ; 
            fwrite_be(&xx, sizeof(double), 1, fich) ;	
        }
    }

    // Writing of various fields
    // -------------------------
    
    // N
    nn() ;     // forces the update at the current time step
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (n_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) n_evol[j].sauve(fich) ; 
    }

    // beta
    beta() ;     // forces the update at the current time step
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (beta_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) beta_evol[j].sauve(fich) ; 
    }

    // Case of a complete save
    // -----------------------
    if (!partial_save) {
    
        cout << "Time_slice::sauve: the full writing is not ready yet !" 
             << endl ; 
        abort() ; 
    }


}









                
                
                
                

