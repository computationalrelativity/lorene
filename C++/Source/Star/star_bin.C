/*
 * Methods for the class Star_bin
 *
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char star_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/02/21 17:05:13  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.4  2004/01/22 10:07:18  f_limousin
 * Add methods set_logn_comp() and set_shift_auto().
 *
 * Revision 1.3  2004/01/20 15:17:34  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h" 
#include "star.h"
#include "eos.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Star_bin::Star_bin(Map& mpi, int nzet_i, const Eos& eos_i, 
		     bool irrot, bool conf_flat0)
    : Star(mpi, nzet_i, eos_i), 
      irrotational(irrot), 
      psi0(mpi), 
      d_psi(mpi, COV, mpi.get_bvect_spher()), 
      wit_w(mpi, CON, mpi.get_bvect_spher()), 
      loggam(mpi), 
      bsn(mpi, CON, mpi.get_bvect_spher()), 
      pot_centri(mpi), 
      logn_auto(mpi),
      logn_comp(mpi), 
      dcov_logn(mpi, COV, mpi.get_bvect_spher()),
      dcon_logn(mpi, CON, mpi.get_bvect_spher()),
      shift_auto(mpi, CON, mpi.get_bvect_spher()), 
      shift_comp(mpi, CON, mpi.get_bvect_spher()), 
      qq_auto(mpi),
      qq_comp(mpi),
      psi4(mpi),
      dcov_lnpsi(mpi, COV, mpi.get_bvect_spher()),
      dcon_lnpsi(mpi, CON, mpi.get_bvect_spher()),
      flat(mpi, mpi.get_bvect_spher()),
      gtilde(flat),
      hij(mpi, CON, mpi.get_bvect_spher()),
      hij_auto(mpi, CON, mpi.get_bvect_spher()),
      hij_comp(mpi, CON, mpi.get_bvect_spher()), 
      tkij_auto(mpi, CON, mpi.get_bvect_spher()), 
      tkij_comp(mpi, CON, mpi.get_bvect_spher()), 
      kcar_auto(mpi), 
      kcar_comp(mpi), 
      ssjm1_logn(mpi),
      ssjm1_qq(mpi),
      ssjm1_h00(mpi),
      ssjm1_h10(mpi),
      ssjm1_h20(mpi),
      ssjm1_h11(mpi),
      ssjm1_h21(mpi),
      ssjm1_h22(mpi),
      decouple(mpi),
      conf_flat(conf_flat0){
    
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    // All quantities are initialized to zero : 
    psi0 = 0 ; 
    d_psi.set_etat_zero() ; 
    wit_w.set_etat_zero() ; 
    loggam = 0 ; 
    bsn.set_etat_zero() ; 
    pot_centri = 0 ;
  
    logn_auto = 0 ;
    logn_comp = 0 ; 
    dcov_logn.set_etat_zero() ;
    dcon_logn.set_etat_zero() ;
    shift_auto.set_etat_zero() ; 
    shift_comp.set_etat_zero() ; 
    qq_auto = 1 ;
    qq_comp = 0 ;
    psi4 = 1 ;
    dcov_lnpsi.set_etat_zero() ;
    dcon_lnpsi.set_etat_zero() ;
    hij.set_etat_zero() ;
    hij_auto.set_etat_zero() ;
    hij_comp.set_etat_zero() ;

    tkij_auto.set_etat_zero() ; 
    tkij_comp.set_etat_zero() ; 
    kcar_auto = 0 ;
    kcar_comp = 0 ; 
}

// Copy constructor
// ----------------
Star_bin::Star_bin(const Star_bin& star)
		       : Star(star), 
			 irrotational(star.irrotational), 
			 psi0(star.psi0), 
			 d_psi(star.d_psi), 
			 wit_w(star.wit_w), 
			 loggam(star.loggam), 
			 bsn(star.bsn), 
			 pot_centri(star.pot_centri), 
			 logn_auto(star.logn_auto),
			 logn_comp(star.logn_comp), 
			 dcov_logn(star.dcov_logn),
			 dcon_logn(star.dcon_logn),
			 shift_auto(star.shift_auto), 
			 shift_comp(star.shift_comp), 
			 qq_auto(star.qq_auto),
			 qq_comp(star.qq_comp),
			 psi4(star.psi4),
			 dcov_lnpsi(star.dcov_lnpsi),
			 dcon_lnpsi(star.dcon_lnpsi),
			 flat(star.flat),
			 gtilde(star.gtilde),
			 hij(star.hij),
			 hij_auto(star.hij_auto),
			 hij_comp(star.hij_comp),
			 tkij_auto(star.tkij_auto), 
			 tkij_comp(star.tkij_comp), 
			 kcar_auto(star.kcar_auto), 
			 kcar_comp(star.kcar_comp), 
			 ssjm1_logn(star.ssjm1_logn),
			 ssjm1_qq(star.ssjm1_qq),
			 ssjm1_h00(star.ssjm1_h00),
			 ssjm1_h10(star.ssjm1_h10),
			 ssjm1_h20(star.ssjm1_h20),
			 ssjm1_h11(star.ssjm1_h11),
			 ssjm1_h21(star.ssjm1_h21),
			 ssjm1_h22(star.ssjm1_h22),
			 decouple(star.decouple),
			 conf_flat(star.conf_flat)
{
    set_der_0x0() ;    

}    

// Constructor from a file
// -----------------------
Star_bin::Star_bin(Map& mpi, const Eos& eos_i, FILE* fich)
		       : Star(mpi, eos_i, fich), 
			 psi0(mpi), 
			 d_psi(mpi, COV, mpi.get_bvect_spher()), 
			 wit_w(mpi, CON, mpi.get_bvect_spher()), 
			 loggam(mpi), 
			 bsn(mpi, CON, mpi.get_bvect_spher()), 
			 pot_centri(mpi), 
			 logn_auto(mpi, *(mpi.get_mg()), fich),
			 logn_comp(mpi), 
			 dcov_logn(mpi, COV, mpi.get_bvect_spher()),
			 dcon_logn(mpi, CON, mpi.get_bvect_spher()),
			 shift_auto(mpi, mpi.get_bvect_spher(), fich), 
			 shift_comp(mpi, CON, mpi.get_bvect_spher()), 
			 qq_auto(mpi, *(mpi.get_mg()), fich),
			 qq_comp(mpi),
			 psi4(mpi),
			 dcov_lnpsi(mpi, COV, mpi.get_bvect_spher()),
			 dcon_lnpsi(mpi, CON, mpi.get_bvect_spher()),
			 flat(mpi, mpi.get_bvect_spher()),
			 gtilde(flat),
			 hij(mpi, CON, mpi.get_bvect_spher()),
			 hij_auto(mpi, mpi.get_bvect_spher(), fich),
			 hij_comp(mpi, CON, mpi.get_bvect_spher()),
     			 tkij_auto(mpi, CON, mpi.get_bvect_spher()), 
			 tkij_comp(mpi, CON, mpi.get_bvect_spher()), 
			 kcar_auto(mpi), 
			 kcar_comp(mpi), 
			 ssjm1_logn(mpi),
			 ssjm1_qq(mpi),
			 ssjm1_h00(mpi),
			 ssjm1_h10(mpi),
			 ssjm1_h20(mpi),
			 ssjm1_h11(mpi),
			 ssjm1_h21(mpi),
			 ssjm1_h22(mpi),
			 decouple(mpi){

    // Etoile parameters
    // -----------------

    // irrotational and conf_flat is read in the file:     
    fread(&irrotational, sizeof(bool), 1, fich) ;
    fread(&conf_flat, sizeof(bool), 1, fich) ;
    	  
   
    // Read of the saved fields:
    // ------------------------

    if (irrotational) {
	Scalar gam_euler_file(mp, *(mp.get_mg()), fich) ; 
	gam_euler = gam_euler_file ; 	

	Scalar psi0_file(mp, *(mp.get_mg()), fich) ; 
	psi0 = psi0_file ; 
    }

    // All other fields are initialized to zero : 
    // ----------------------------------------

    d_psi.set_etat_zero() ; 
    wit_w.set_etat_zero() ; 
    loggam = 0 ; 
    bsn.set_etat_zero() ; 
    pot_centri = 0 ;
    logn_comp = 0 ; 
    dcov_logn.set_etat_zero() ;
    dcon_logn.set_etat_zero() ;
    shift_comp.set_etat_zero() ; 
    qq_comp = 0 ;
    psi4 = 1 ;
    dcov_lnpsi.set_etat_zero() ;
    dcon_lnpsi.set_etat_zero() ;
    hij.set_etat_zero() ;
    hij_comp.set_etat_zero() ;

    tkij_auto.set_etat_zero() ; 
    tkij_comp.set_etat_zero() ; 
    kcar_auto = 0 ;
    kcar_comp = 0 ; 
 
    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}

			    //------------//
			    // Destructor //
			    //------------//

Star_bin::~Star_bin(){

    del_deriv() ; 

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Star_bin::del_deriv() const {

    Star::del_deriv() ; 

    if (p_xa_barycenter != 0x0) delete p_xa_barycenter ; 
    
    set_der_0x0() ; 
}			    




void Star_bin::set_der_0x0() const {

    Star::set_der_0x0() ;

    p_xa_barycenter = 0x0 ; 

}			    

void Star_bin::del_hydro_euler() {

    Star::del_hydro_euler() ; 

    del_deriv() ; 

}			    


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Star_bin
// --------------------------------
void Star_bin::operator=(const Star_bin& star) {

    // Assignment of quantities common to the derived classes of Star
    Star::operator=(star) ;	    

    // Assignement of proper quantities of class Star_bin
    irrotational = star.irrotational ; 
    psi0 = star.psi0 ; 
    d_psi = star.d_psi ;
    wit_w = star.wit_w ; 
    loggam = star.loggam ;
    bsn = star.bsn ;
    pot_centri = star.pot_centri ;
    logn_auto = star.logn_auto ;    
    logn_comp = star.logn_comp ;
    dcov_logn = star.dcov_logn ;
    dcon_logn = star.dcon_logn ;
    shift_auto = star.shift_auto ;
    shift_comp = star.shift_comp ;
    qq_auto = star.qq_auto ;
    qq_comp = star.qq_comp ;
    psi4 = star.psi4 ;
    dcov_lnpsi = star.dcov_lnpsi ;
    dcon_lnpsi = star.dcon_lnpsi ;
    flat = star.flat ;
    gtilde = star.gtilde ;
    hij = star.hij ;
    hij_auto = star.hij_auto ;
    hij_comp = star.hij_comp ; 
    tkij_auto = star.tkij_auto ;
    tkij_comp = star.tkij_comp ;
    kcar_auto = star.kcar_auto ;
    kcar_comp = star.kcar_comp ;
    ssjm1_logn = star.ssjm1_logn ;
    ssjm1_qq = star.ssjm1_qq ;
    ssjm1_h00 = star.ssjm1_h00 ;
    ssjm1_h10 = star.ssjm1_h10 ;
    ssjm1_h20 = star.ssjm1_h20 ;
    ssjm1_h11 = star.ssjm1_h11 ;
    ssjm1_h21 = star.ssjm1_h21 ;
    ssjm1_h22 = star.ssjm1_h22 ;
    decouple = star.decouple ;
    conf_flat = star.conf_flat ;
    
    del_deriv() ;  // Deletes all derived quantities

}	

Scalar& Star_bin::set_pot_centri() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return pot_centri ;
    
} 

Scalar& Star_bin::set_logn_comp() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return logn_comp ;
    
} 


Vector& Star_bin::set_shift_auto() {

    return shift_auto ;
    
} 


			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Star_bin::sauve(FILE* fich) const {
    
    Star::sauve(fich) ; 
    
    logn_auto.sauve(fich) ;
    shift_auto.sauve(fich) ;
    qq_auto.sauve(fich) ;
    hij_auto.sauve(fich) ;

    fwrite(&irrotational, sizeof(bool), 1, fich) ;		
    fwrite(&conf_flat, sizeof(bool), 1, fich) ;		
    
    if (irrotational) {
	gam_euler.sauve(fich) ; // required to construct d_psi from psi0
	psi0.sauve(fich) ; 
    }
 
}

// Printing
// --------

ostream& Star_bin::operator>>(ostream& ost) const {
    
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << qpig << msol << f_unit << mevpfm3 << endl ; 
    }    

    Star::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Star in a binary system" << endl ; 
    ost << "-----------------------" << endl ; 
    
    if (irrotational) {
	ost << "irrotational configuration" << endl ; 
    }
    else {
	ost << "corotating configuration" << endl ; 
    }
       
    ost << "Absolute abscidia of the stellar center: " <<
	mp.get_ori_x() / km << " km" << endl ; 
    
    ost << "Absolute abscidia of the barycenter of the baryon density : " <<
	xa_barycenter() / km << " km" << endl ; 
    
    double r_0 = 0.5 * ( ray_eq() + ray_eq_pi() ) ; 
    double d_ns = fabs( mp.get_ori_x() ) + ray_eq_pi() - r_0 ;
    double d_tilde = 2 * d_ns / r_0 ;  
    
    ost << "d_tilde : " << d_tilde << endl ; 

    ost << "Central value of gam_euler : " 
        << gam_euler.val_grid_point(0, 0, 0, 0)  << endl ; 

    ost << "Central u_euler (U^X, U^Y, U^Z) [c] : " 
	<< u_euler(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< u_euler(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< u_euler(3).val_grid_point(0, 0, 0, 0) << endl ; 

    if (irrotational) {
    ost << "Central d_psi (X, Y, Z) [c] :         " 
	    << d_psi(1).val_grid_point(0, 0, 0, 0) << "  " 
	    << d_psi(2).val_grid_point(0, 0, 0, 0) << "  " 
	    << d_psi(3).val_grid_point(0, 0, 0, 0) << endl ; 

	ost << "Central vel. / co-orb. (W^X, W^Y, W^Z) [c] : " 
	    << wit_w(1).val_grid_point(0, 0, 0, 0) << "  " 
	    << wit_w(2).val_grid_point(0, 0, 0, 0) << "  " 
	    << wit_w(3).val_grid_point(0, 0, 0, 0) << endl ; 

	ost << "Max vel. / co-orb. (W^X, W^Y, W^Z) [c] : " 
	    << max(max(wit_w(1))) << "  " 
	    << max(max(wit_w(2))) << "  " 
	    << max(max(wit_w(3))) << endl ; 

	ost << "Min vel. / co-orb. (W^X, W^Y, W^Z) [c] : " 
	    << min(min(wit_w(1))) << "  " 
	    << min(min(wit_w(2))) << "  " 
	    << min(min(wit_w(3))) << endl ; 

	double r_surf = mp.val_r(0,1.,M_PI/4,M_PI/4) ;

	ost << "Velocity at (r_surf,pi/4,pi/4) / co-orb. [c] : "
	    << wit_w(1).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(2).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(3).val_point(r_surf,M_PI/4,M_PI/4) << endl ;

	ost << "Central value of loggam : " 
	    << loggam.val_grid_point(0, 0, 0, 0)  << endl ; 	
    }


    ost << "Central value of log(N) auto, comp :         " 
	<< logn_auto.val_grid_point(0, 0, 0, 0) << "  " 
	<< logn_comp.val_grid_point(0, 0, 0, 0) << endl ; 

    ost << "Central value of shift (N^X, N^Y, N^Z) [c] : " 
	<< shift(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift(3).val_grid_point(0, 0, 0, 0) << endl ; 

    ost << "  ... shift_auto part of it [c] :            " 
	<< shift_auto(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift_auto(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift_auto(3).val_grid_point(0, 0, 0, 0) << endl ; 

    ost << "  ... shift_comp part of it [c] :            " 
	<< shift_comp(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift_comp(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< shift_comp(3).val_grid_point(0, 0, 0, 0) << endl ; 

    ost << endl << "Central value of (B^X, B^Y, B^Z)/N [c] : " 
	<< bsn(1).val_grid_point(0, 0, 0, 0) << "  " 
	<< bsn(2).val_grid_point(0, 0, 0, 0) << "  " 
	<< bsn(3).val_grid_point(0, 0, 0, 0) << endl ; 


    ost << endl << "Central A^2 K^{ij} [c/km] : " << endl ; 
    ost << "  A^2 K^{xx} auto, comp : " 
	<< tkij_auto(1, 1).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 1).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{xy} auto, comp : " 
	<< tkij_auto(1, 2).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 2).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{xz} auto, comp : " 
	<< tkij_auto(1, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 3).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{yy} auto, comp : " 
	<< tkij_auto(2, 2).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(2, 2).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{yz} auto, comp : " 
	<< tkij_auto(2, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(2, 3).val_grid_point(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{zz} auto, comp : " 
	<< tkij_auto(3, 3).val_grid_point(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(3, 3).val_grid_point(0, 0, 0, 0) * km << endl ; 

    ost << endl << "Central A^2 K_{ij} K^{ij} [c^2/km^2] : " << endl ; 
    ost << "   A^2 K_{ij} K^{ij}  auto, comp : " 
	<< kcar_auto.val_grid_point(0, 0, 0, 0) * km*km  << "  "
	<< kcar_comp.val_grid_point(0, 0, 0, 0) * km*km << endl ; 

    
    return ost ; 
}

			    //-------------------------//
			    //	Computational routines //
			    //-------------------------//


Tensor Star_bin::sprod(const Tensor& t1, const Tensor& t2) const {
     
    Tensor* p_tens_metr = 0x0 ;

   // Both indices should be contravariant or both covariant : 
    if (t1.get_index_type(t1.get_valence()-1) == CON) {
      assert( t2.get_index_type(0) == CON ) ;
      
      p_tens_metr = new Tensor(gamma.cov()) ;
	}
    
    if (t1.get_index_type(t1.get_valence()-1) == COV) {
      assert( t2.get_index_type(0) == COV ) ;

       p_tens_metr = new Tensor(gamma.con()) ;
       }

    // Verifs :
    assert (t1.get_mp() == t2.get_mp()) ;
    
    // Contraction possible ?
     if ( (t1.get_valence() != 0) && (t2.get_valence() != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    int val_res = t1.get_valence() + t2.get_valence() - 2;

    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;

    for (int i=0 ; i<t1.get_valence() - 1 ; i++)
	tipe.set(i) = t1.get_index_type(i) ;
    for (int i = t1.get_valence()-1 ; i<val_res ; i++)
	tipe.set(i) = t2.get_index_type(i-t1.get_valence()+2) ;
	
    Tensor res(t1.get_mp(), val_res, tipe, t1.get_triad()) ;
	
    Scalar work(t1.get_mp()) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_t1(t1.get_valence()) ;
    Itbl jeux_indice_t2(t2.get_valence()) ;
    jeux_indice_t1.set_etat_qcq() ;
    jeux_indice_t2.set_etat_qcq() ;
    
    for (int ir=0 ; ir<res.get_n_comp() ; ir++) {    // Boucle sur les composantes
					       // du resultat 

	// Indices du resultat correspondant a la position ir : 
	Itbl jeux_indice_res = res.indices(ir) ;

	// Premiers indices correspondant dans t1 : 
	for (int i=0 ; i<t1.get_valence() - 1 ; i++) {
	    jeux_indice_t1.set(i) = jeux_indice_res(i) ;
	}
	
	// Derniers indices correspondant dans t2 : 
	for (int i=1 ; i<t2.get_valence() ; i++) {
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.get_valence()+i-2) ;
	}
	
	work.set_etat_zero() ;

	// Sommation sur le dernier indice de t1 et le premier de t2 : 
	
	for (int i=1 ; i<=3 ; i++) {
	  for (int j=1 ; j<=3 ; j++) {
	    jeux_indice_t1.set(t1.get_valence() - 1) = i ;
	    jeux_indice_t2.set(0) = j ;
	    work = work + (*p_tens_metr)(i,j)*t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	    }
	}
	res.set(jeux_indice_res) = work ;
    
    }	// fin de la boucle sur les composantes du resultat
    delete p_tens_metr ; 
    return res ;

}



void Star_bin::fait_d_psi() {

    if (!irrotational) {
	d_psi.set_etat_nondef() ; 
	return ; 
    }
/*
    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Scalar hhh = exp(ent) ;  // = 1 at the Newtonian limit
 
    //  Computation of W^i = - psi^4 h Gamma_n B^i/N
    //----------------------------------------------

    Vector www = - psi4 * hhh * gam_euler * bsn ; 
    
    
    // Constant value of W^i at the center of the star
    //-------------------------------------------------
    
    Vector v_orb(mp, COV, mp.get_bvect_spher()) ; 
    
    for (int i=1; i<=3; i++) {
	v_orb.set(i) = www(i).val_grid_point(0, 0, 0, 0) ; 
    }
    
    // Gradient of psi 
    //----------------

    Vector d_psi0 = psi0.derive_cov(flat) ; 
    
    d_psi0.change_triad( mp.get_bvect_spher() ) ; 

    d_psi = d_psi0 + v_orb ; 
    

    // C^1 continuation of d_psi outside the star
    //  (to ensure a smooth enthalpy field accross the stellar surface)
    // ----------------------------------------------------------------
    
    int nzm1 = mp.get_mg()->get_nzone() - 1 ;    

    if (d_psi0.get_etat() == ETATQCQ ) {
	d_psi.annule(nzet, nzm1) ;	 
	for (int i=0; i<3; i++) {
	    d_psi.set(i).set_spectral_va().set_base( d_psi0(i).
						   get_spectral_va().base ) ; 
	    d_psi.set(i) = raccord_c1(d_psi(i), nzet) ; 
	}
    }
    else{ 
	d_psi.annule_domain(nzm1) ;	 
    }

*/ 
} 


void Star_bin::relaxation(const Star_bin& star_jm1, double relax_ent, 
			    double relax_met, int mer, int fmer_met) {
				
    double relax_ent_jm1 = 1. - relax_ent ; 
    double relax_met_jm1 = 1. - relax_met ; 

    ent = relax_ent * ent + relax_ent_jm1 * star_jm1.ent ; 

    if ( (mer != 0) && (mer % fmer_met == 0)) {

	logn_auto = relax_met * logn_auto + relax_met_jm1 * star_jm1.logn_auto ;
	qq_auto = relax_met * qq_auto + relax_met_jm1 * star_jm1.qq_auto ;
	
	shift_auto = relax_met * shift_auto 
					+ relax_met_jm1 * star_jm1.shift_auto ;
	
	hij_auto = relax_met * hij_auto + relax_met_jm1 * star_jm1.hij_auto ;
	
    }

    del_deriv() ; 
    
    equation_of_state() ; 

}

