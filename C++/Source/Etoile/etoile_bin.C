/*
 * Methods for the class Etoile_bin
 *
 * (see file etoile.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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


char etoile_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.31  2001/06/25  12:53:15  eric
 * Ajout du membre p_companion et des fonctions associees set_companion() et get_companion().
 *
 * Revision 2.30  2000/12/22  13:09:09  eric
 * Modif fait_d_psi : prolongement C^1 de dpsi en dehors de l'etoile.
 *
 * Revision 2.29  2000/09/29  11:54:40  keisuke
 * Add the relaxations for logn_auto_div and d_logn_auto_div.
 *
 * Revision 2.28  2000/09/29  09:57:26  keisuke
 * Add the relaxation for logn_auto_regu.
 *
 * Revision 2.27  2000/09/22  15:51:19  keisuke
 * d_logn_auto_div devient desormais un membre de la classe Etoile
 * et non plus de la classe derivee Etoile_bin.
 *
 * Revision 2.26  2000/09/07  14:35:31  keisuke
 * Ajout du membre d_logn_auto_regu.
 *
 * Revision 2.25  2000/08/29  11:38:03  eric
 * Ajout du membre d_logn_auto_div.
 *
 * Revision 2.24  2000/07/06  09:53:22  eric
 * Ajout de xa_barycenter dans l'affichage.
 *
 * Revision 2.23  2000/07/06  09:40:11  eric
 *  Ajout du membre derive p_xa_barycenter.
 *
 * Revision 2.22  2000/06/15  15:50:02  eric
 * Ajout du calcul de d_tilde dans l'affichage.
 *
 * Revision 2.21  2000/05/23  12:28:00  phil
 * changement apres modification de skxk
 *
 * Revision 2.20  2000/03/15  11:04:58  eric
 * Ajout des fonctions Etoile_bin::set_w_shift() et Etoile_bin::set_khi_shift()
 *
 * Revision 2.19  2000/02/24  09:12:56  keisuke
 * Add an output for the velocity field in the corotating frame.
 *
 * Revision 2.18  2000/02/21  16:31:58  eric
 * Modif affichage.
 *
 * Revision 2.17  2000/02/21  14:54:13  eric
 * fait_d_psi appel d_psi.set_etat_nondef() et return dans le cas
 * corotation.
 *
 * Revision 2.16  2000/02/21  14:31:43  eric
 * gam_euler est desormais sauve dans le fichier (cas irrotationnel seulement)
 * psi0 n'est sauve que dans le cas irrotationnel.
 *
 * Revision 2.15  2000/02/21  13:58:22  eric
 * Suppression du membre psi: remplacement par psi0.
 *
 * Revision 2.14  2000/02/17  15:30:04  eric
 * Ajout de la fonction Etoile_bin::relaxation.
 *
 * Revision 2.13  2000/02/17  14:42:09  eric
 * Modif affichage.
 *
 * Revision 2.12  2000/02/16  17:12:25  eric
 * Modif initialisation de w_shift, khi_shift et ssjm1_wshift dans le
 * constructeur standard.
 *
 * Revision 2.11  2000/02/16  16:08:10  eric
 * w_shift et ssjm1_wshift sont desormais definis sur la triade du mapping.
 *
 * Revision 2.10  2000/02/16  15:06:11  eric
 *  Ajout des membres w_shift et khi_shift.
 * (sauves dans les fichiers a la place de shift_auto).
 * Ajout de la fonction Etoile_bin::fait_shift_auto.
 *
 * Revision 2.9  2000/02/16  13:47:22  eric
 * Ajout des membres ssjm1_khi et ssjm1_wshift.
 * (sauvegardes dans les fichiers).
 *
 * Revision 2.8  2000/02/16  11:54:40  eric
 * Ajout des membres ssjm1_logn et ssjm1_beta (sauvegarde dans les fichiers).
 *
 * Revision 2.7  2000/02/10  18:56:52  eric
 * Modifs routine d'affichage.
 *
 * Revision 2.6  2000/02/10  16:12:05  eric
 * Ajout de la fonction fait_d_psi
 *
 * Revision 2.5  2000/02/09  20:24:03  eric
 * Appel de set_triad(ref_triad) sur u_euler et shift dans les
 * constructeurs.
 * ,.
 *
 * Revision 2.4  2000/02/09  19:31:33  eric
 * La triade de decomposition doit desormais figurer en argument des
 *  constructeurs de Tenseur.
 *
 * Revision 2.3  2000/02/08  19:29:50  eric
 * La fonction Etoile_bin::scal_prod est rebaptisee Etoile_bin::sprod
 *
 * Revision 2.2  2000/02/04  17:15:36  eric
 * Ajout du membre ref_triad.
 *
 * Revision 2.1  2000/02/01  16:00:00  eric
 * Ajout de la fonction Etoile_bin::scal_prod.
 *
 * Revision 2.0  2000/01/31  15:57:12  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "etoile.h"
#include "eos.h"

// Local prototype
Cmp raccord_c1(const Cmp& uu, int l1) ; 

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Etoile_bin::Etoile_bin(Map& mpi, int nzet_i, bool relat, const Eos& eos_i, 
		       bool irrot, const Base_vect& ref_triad_i)
		       : Etoile(mpi, nzet_i, relat, eos_i), 
			 irrotational(irrot), 
			 ref_triad(ref_triad_i), 
			 psi0(mpi), 
			 d_psi(mpi, 1, COV, ref_triad), 
			 wit_w(mpi, 1, CON, ref_triad), 
			 loggam(mpi), 
			 logn_comp(mpi), 
			 d_logn_auto(mpi, 1, COV, ref_triad), 
			 d_logn_auto_regu(mpi, 1, COV, ref_triad), 
			 d_logn_comp(mpi, 1, COV, ref_triad), 
			 beta_comp(mpi), 
			 d_beta_auto(mpi, 1, COV, ref_triad), 
			 d_beta_comp(mpi, 1, COV, ref_triad), 
			 shift_auto(mpi, 1, CON, ref_triad), 
			 shift_comp(mpi, 1, CON, ref_triad), 
			 w_shift(mpi, 1, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 tkij_auto(mpi, 2, CON, ref_triad), 
			 tkij_comp(mpi, 2, CON, ref_triad), 
			 akcar_auto(mpi), 
			 akcar_comp(mpi), 
			 bsn(mpi, 1, CON, ref_triad), 
			 pot_centri(mpi), 
			 ssjm1_logn(mpi), 
			 ssjm1_beta(mpi), 			 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, 1, CON, mp.get_bvect_cart())			 
{
    
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    // Pointer to the companion star is initialized to 0x0 (the function
    //  set_companion() should be called later)
    p_companion = 0x0 ; 
    
    // The reference triad is assigned to the vectors u_euler and shift :
    u_euler.set_triad(ref_triad) ; 
    shift.set_triad(ref_triad) ;

    // All quantities are initialized to zero : 
    psi0 = 0 ; 
    d_psi = 0 ; 
    wit_w = 0 ; 
    loggam = 0 ; 
    logn_comp = 0 ; 
    d_logn_auto = 0 ; 
    d_logn_auto_regu = 0 ; 
    d_logn_comp = 0 ; 
    beta_comp = 0 ; 
    d_beta_auto = 0 ; 
    d_beta_comp = 0 ; 
    shift_auto = 0 ; 
    shift_comp = 0 ; 

    w_shift.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	w_shift.set(i) = 0 ; 
    }

    khi_shift.set_etat_qcq() ; 
    khi_shift.set() = 0 ; 

    tkij_auto.set_etat_zero() ; 
    tkij_comp.set_etat_zero() ; 
    akcar_auto = 0 ;
    akcar_comp = 0 ; 
    bsn = 0 ; 
    pot_centri = 0 ;
    ssjm1_logn = 0 ; 
    ssjm1_beta = 0 ; 
    ssjm1_khi = 0 ; 
    
    ssjm1_wshift.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	ssjm1_wshift.set(i) = 0 ; 
    }
    
}

// Copy constructor
// ----------------
Etoile_bin::Etoile_bin(const Etoile_bin& et)
		       : Etoile(et), 
			 irrotational(et.irrotational), 
			 ref_triad(et.ref_triad), 
			 psi0(et.psi0), 
			 d_psi(et.d_psi), 
			 wit_w(et.wit_w), 
			 loggam(et.loggam), 
			 logn_comp(et.logn_comp), 
			 d_logn_auto(et.d_logn_auto), 
			 d_logn_auto_regu(et.d_logn_auto_regu), 
			 d_logn_comp(et.d_logn_comp), 
			 beta_comp(et.beta_comp), 
			 d_beta_auto(et.d_beta_auto), 
			 d_beta_comp(et.d_beta_comp), 
			 shift_auto(et.shift_auto), 
			 shift_comp(et.shift_comp), 
			 w_shift(et.w_shift), 
			 khi_shift(et.khi_shift), 
			 tkij_auto(et.tkij_auto), 
			 tkij_comp(et.tkij_comp), 
			 akcar_auto(et.akcar_auto), 
			 akcar_comp(et.akcar_comp), 
			 bsn(et.bsn), 
			 pot_centri(et.pot_centri), 
			 ssjm1_logn(et.ssjm1_logn), 
			 ssjm1_beta(et.ssjm1_beta), 
			 ssjm1_khi(et.ssjm1_khi), 
			 ssjm1_wshift(et.ssjm1_wshift)
{
    set_der_0x0() ;    

    // Pointer to the companion star is initialized to 0x0 (the function
    //  set_companion() should be called later)
    p_companion = 0x0 ; 
    
}    

// Constructor from a file
// -----------------------
Etoile_bin::Etoile_bin(Map& mpi, const Eos& eos_i, const Base_vect& ref_triad_i,
		       FILE* fich)
		       : Etoile(mpi, eos_i, fich), 
			 ref_triad(ref_triad_i), 
			 psi0(mpi), 
			 d_psi(mpi, 1, COV, ref_triad), 
			 wit_w(mpi, 1, CON, ref_triad), 
			 loggam(mpi), 
			 logn_comp(mpi), 
			 d_logn_auto(mpi, 1, COV, ref_triad), 
			 d_logn_auto_regu(mpi, 1, COV, ref_triad), 
			 d_logn_comp(mpi, 1, COV, ref_triad), 
			 beta_comp(mpi), 
			 d_beta_auto(mpi, 1, COV, ref_triad), 
			 d_beta_comp(mpi, 1, COV, ref_triad), 
			 shift_auto(mpi, 1, CON, ref_triad), 
			 shift_comp(mpi, 1, CON, ref_triad), 
			 w_shift(mpi, 1, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 tkij_auto(mpi, 2, CON, ref_triad), 
			 tkij_comp(mpi, 2, CON, ref_triad), 
			 akcar_auto(mpi), 
			 akcar_comp(mpi), 
			 bsn(mpi, 1, CON, ref_triad), 
			 pot_centri(mpi), 
			 ssjm1_logn(mpi), 
			 ssjm1_beta(mpi), 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, 1, CON, mp.get_bvect_cart())
{

    // The reference triad is assigned to the vectors u_euler and shift :
    u_euler.set_triad(ref_triad) ; 
    shift.set_triad(ref_triad) ;

    // Etoile parameters
    // -----------------

    // irrotational is read in the file:     
    fread(&irrotational, sizeof(bool), 1, fich) ;		
    	  
   
    // Read of the saved fields:
    // ------------------------

    if (irrotational) {
	Tenseur gam_euler_file(mp, fich) ; 
	gam_euler = gam_euler_file ; 	

	Tenseur psi0_file(mp, fich) ; 
	psi0 = psi0_file ; 
    }

    
    Tenseur w_shift_file(mp, mp.get_bvect_cart(), fich) ; 
    w_shift = w_shift_file ;
    
    Tenseur khi_shift_file(mp, fich) ; 
    khi_shift = khi_shift_file ;
    
    fait_shift_auto() ;	    // constructs shift_auto from w_shift and khi_shift
    
    Cmp ssjm1_logn_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_logn = ssjm1_logn_file ; 

    Cmp ssjm1_beta_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_beta = ssjm1_beta_file ; 

    Cmp ssjm1_khi_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_khi = ssjm1_khi_file ; 

    Tenseur ssjm1_wshift_file(mp, mp.get_bvect_cart(), fich) ; 
    ssjm1_wshift = ssjm1_wshift_file ; 

    // All other fields are initialized to zero : 
    // ----------------------------------------
    d_psi = 0 ; 
    wit_w = 0 ; 
    loggam = 0 ; 
    logn_comp = 0 ; 
    d_logn_auto = 0 ; 
    d_logn_auto_regu = 0 ; 
    d_logn_comp = 0 ; 
    beta_comp = 0 ; 
    d_beta_auto = 0 ; 
    d_beta_comp = 0 ; 
    shift_comp = 0 ; 
    tkij_auto.set_etat_zero() ; 
    tkij_comp.set_etat_zero() ; 
    akcar_auto = 0 ;
    akcar_comp = 0 ; 
    bsn = 0 ; 
    pot_centri = 0 ;

    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
    // Pointer to the companion star is initialized to 0x0 (the function
    //  set_companion() should be called later)
    p_companion = 0x0 ; 
    
}

			    //------------//
			    // Destructor //
			    //------------//

Etoile_bin::~Etoile_bin(){

    del_deriv() ; 

}

			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Etoile_bin::del_deriv() const {

    Etoile::del_deriv() ; 

    if (p_xa_barycenter != 0x0) delete p_xa_barycenter ; 
    
    set_der_0x0() ; 
}			    




void Etoile_bin::set_der_0x0() const {

    Etoile::set_der_0x0() ;

    p_xa_barycenter = 0x0 ; 

}			    

void Etoile_bin::del_hydro_euler() {

    Etoile::del_hydro_euler() ; 

    del_deriv() ; 

}			    


			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Etoile_bin
// --------------------------------
void Etoile_bin::operator=(const Etoile_bin& et) {

    // Assignment of quantities common to the derived classes of Etoile
    Etoile::operator=(et) ;	    

    // Assignement of proper quantities of class Etoile_bin
    irrotational = et.irrotational ; 

    assert(et.ref_triad == ref_triad) ; 

    psi0 = et.psi0 ; 
    d_psi = et.d_psi ;
    wit_w = et.wit_w ; 
    loggam = et.loggam ;
    logn_comp = et.logn_comp ;
    d_logn_auto = et.d_logn_auto ;
    d_logn_auto_regu = et.d_logn_auto_regu ;
    d_logn_comp = et.d_logn_comp ;
    beta_comp = et.beta_comp ;
    d_beta_auto = et.d_beta_auto ;
    d_beta_comp = et.d_beta_comp ;
    shift_auto = et.shift_auto ;
    shift_comp = et.shift_comp ;
    w_shift = et.w_shift ;
    khi_shift = et.khi_shift ;
    tkij_auto = et.tkij_auto ;
    tkij_comp = et.tkij_comp ;
    akcar_auto = et.akcar_auto ;
    akcar_comp = et.akcar_comp ;
    bsn = et.bsn ;
    pot_centri = et.pot_centri ;
    ssjm1_logn = et.ssjm1_logn ;
    ssjm1_beta = et.ssjm1_beta ; 
    ssjm1_khi = et.ssjm1_khi ;
    ssjm1_wshift = et.ssjm1_wshift ; 
    
    del_deriv() ;  // Deletes all derived quantities

}	

Tenseur& Etoile_bin::set_logn_comp() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return logn_comp ;
    
} 

Tenseur& Etoile_bin::set_pot_centri() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return pot_centri ;
    
} 

Tenseur& Etoile_bin::set_w_shift() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return w_shift ;
    
} 

Tenseur& Etoile_bin::set_khi_shift() {

    del_deriv() ;	// sets to 0x0 all the derived quantities
    return khi_shift ;
    
} 

void Etoile_bin::set_companion(const Etoile_bin& comp) {
    
    p_companion = &comp ; 
    
}


			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Etoile_bin::sauve(FILE* fich) const {
    
    Etoile::sauve(fich) ; 
    
    fwrite(&irrotational, sizeof(bool), 1, fich) ;		
    
    if (irrotational) {
	gam_euler.sauve(fich) ; // required to construct d_psi from psi0
	psi0.sauve(fich) ; 
    }
    
    w_shift.sauve(fich) ; 
    khi_shift.sauve(fich) ; 
    
    ssjm1_logn.sauve(fich) ; 
    ssjm1_beta.sauve(fich) ; 
    ssjm1_khi.sauve(fich) ; 
    ssjm1_wshift.sauve(fich) ;
        
}

// Printing
// --------

ostream& Etoile_bin::operator>>(ostream& ost) const {
    
    #include "unites.h"	    
    // To avoid some compilation warnings
    if (&ost == 0x0) {
	cout << qpig << msol << f_unit << mevpfm3 << endl ; 
    }    

    Etoile::operator>>(ost) ; 
    
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

    ost << "Orientation with respect to the absolute frame : " <<
	mp.get_rot_phi() << " rad" << endl ; 

    ost << "Central value of gam_euler : " 
        << gam_euler()(0, 0, 0, 0)  << endl ; 

    ost << "Central u_euler (U^X, U^Y, U^Z) [c] : " 
	<< u_euler(0)(0, 0, 0, 0) << "  " 
	<< u_euler(1)(0, 0, 0, 0) << "  " 
	<< u_euler(2)(0, 0, 0, 0) << endl ; 

    if (irrotational) {
    ost << "Central d_psi (X, Y, Z) [c] :         " 
	    << d_psi(0)(0, 0, 0, 0) << "  " 
	    << d_psi(1)(0, 0, 0, 0) << "  " 
	    << d_psi(2)(0, 0, 0, 0) << endl ; 

	ost << "Central vel. / co-orb. (W^X, W^Y, W^Z) [c] : " 
	    << wit_w(0)(0, 0, 0, 0) << "  " 
	    << wit_w(1)(0, 0, 0, 0) << "  " 
	    << wit_w(2)(0, 0, 0, 0) << endl ; 

	ost << "Max vel. / co-orb. (W^X, W^Y, W^Z) [c] : " 
	    << max(max(wit_w(0))) << "  " 
	    << max(max(wit_w(1))) << "  " 
	    << max(max(wit_w(2))) << endl ; 

	ost << "Min vel. / co-orb. (W^X, W^Y, W^Z) [c] : " 
	    << min(min(wit_w(0))) << "  " 
	    << min(min(wit_w(1))) << "  " 
	    << min(min(wit_w(2))) << endl ; 

	double r_surf = mp.val_r(0,1.,M_PI/4,M_PI/4) ;

	ost << "Velocity at (r_surf,pi/4,pi/4) / co-orb. [c] : "
	    << wit_w(0).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(1).val_point(r_surf,M_PI/4,M_PI/4) << "  "
	    << wit_w(2).val_point(r_surf,M_PI/4,M_PI/4) << endl ;

	ost << "Central value of loggam : " 
	    << loggam()(0, 0, 0, 0)  << endl ; 	
    }


    ost << "Central value of log(N) auto, comp :         " 
	<< logn_auto()(0, 0, 0, 0) << "  " 
	<< logn_comp()(0, 0, 0, 0) << endl ; 

    ost << "Central value of beta=log(AN) auto, comp :   " 
	<< beta_auto()(0, 0, 0, 0) << "  " 
	<< beta_comp()(0, 0, 0, 0) << endl ; 

    ost << "Central value of shift (N^X, N^Y, N^Z) [c] : " 
	<< shift(0)(0, 0, 0, 0) << "  " 
	<< shift(1)(0, 0, 0, 0) << "  " 
	<< shift(2)(0, 0, 0, 0) << endl ; 

    ost << "  ... shift_auto part of it [c] :            " 
	<< shift_auto(0)(0, 0, 0, 0) << "  " 
	<< shift_auto(1)(0, 0, 0, 0) << "  " 
	<< shift_auto(2)(0, 0, 0, 0) << endl ; 

    ost << "  ... shift_comp part of it [c] :            " 
	<< shift_comp(0)(0, 0, 0, 0) << "  " 
	<< shift_comp(1)(0, 0, 0, 0) << "  " 
	<< shift_comp(2)(0, 0, 0, 0) << endl ; 

    ost << "  ... w_shift (NB: components in the star Cartesian frame) [c] :  "
	<< endl  
	<< w_shift(0)(0, 0, 0, 0) << "  " 
	<< w_shift(1)(0, 0, 0, 0) << "  " 
	<< w_shift(2)(0, 0, 0, 0) << endl ; 

    ost << "Central value of khi_shift [km c] : " 
        << khi_shift()(0, 0, 0, 0) / km << endl ; 

    ost << endl << "Central value of (B^X, B^Y, B^Z)/N [c] : " 
	<< bsn(0)(0, 0, 0, 0) << "  " 
	<< bsn(1)(0, 0, 0, 0) << "  " 
	<< bsn(2)(0, 0, 0, 0) << endl ; 

    ost << endl << 
	"Central (d/dX,d/dY,d/dZ)(logn_auto) [km^{-1}] : " 
	<< d_logn_auto(0)(0, 0, 0, 0) * km << "  " 
	<< d_logn_auto(1)(0, 0, 0, 0) * km  << "  " 
	<< d_logn_auto(2)(0, 0, 0, 0) * km  << endl ; 

    ost << "Central (d/dX,d/dY,d/dZ)(logn_comp) [km^{-1}] : " 
	<< d_logn_comp(0)(0, 0, 0, 0) * km  << "  " 
	<< d_logn_comp(1)(0, 0, 0, 0) * km  << "  " 
	<< d_logn_comp(2)(0, 0, 0, 0) * km  << endl ; 

    ost << endl << 
	"Central (d/dX,d/dY,d/dZ)(beta_auto) [km^{-1}] : " 
	<< d_beta_auto(0)(0, 0, 0, 0) * km << "  " 
	<< d_beta_auto(1)(0, 0, 0, 0) * km  << "  " 
	<< d_beta_auto(2)(0, 0, 0, 0) * km  << endl ; 

    ost << "Central (d/dX,d/dY,d/dZ)(beta_comp) [km^{-1}] : " 
	<< d_beta_comp(0)(0, 0, 0, 0) * km  << "  " 
	<< d_beta_comp(1)(0, 0, 0, 0) * km  << "  " 
	<< d_beta_comp(2)(0, 0, 0, 0) * km  << endl ; 


    ost << endl << "Central A^2 K^{ij} [c/km] : " << endl ; 
    ost << "  A^2 K^{xx} auto, comp : " 
	<< tkij_auto(0, 0)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(0, 0)(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{xy} auto, comp : " 
	<< tkij_auto(0, 1)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(0, 1)(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{xz} auto, comp : " 
	<< tkij_auto(0, 2)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(0, 2)(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{yy} auto, comp : " 
	<< tkij_auto(1, 1)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 1)(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{yz} auto, comp : " 
	<< tkij_auto(1, 2)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 2)(0, 0, 0, 0) * km << endl ; 
    ost << "  A^2 K^{zz} auto, comp : " 
	<< tkij_auto(2, 2)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(2, 2)(0, 0, 0, 0) * km << endl ; 

    ost << endl << "Central A^2 K_{ij} K^{ij} [c^2/km^2] : " << endl ; 
    ost << "   A^2 K_{ij} K^{ij}  auto, comp : " 
	<< akcar_auto()(0, 0, 0, 0) * km*km  << "  "
	<< akcar_comp()(0, 0, 0, 0) * km*km << endl ; 

    
    return ost ; 
}

			    //-------------------------//
			    //	Computational routines //
			    //-------------------------//
			    
Tenseur Etoile_bin::sprod(const Tenseur& t1, const Tenseur& t2) const {

    int val1 = t1.get_valence() ; 
    
    // Both indices should be contravariant : 
    assert( t1.get_type_indice(val1-1) == CON ) ; 
    assert( t2.get_type_indice(0) == CON ) ; 

    return a_car * flat_scalar_prod(t1, t2) ; 
    
} 

void Etoile_bin::fait_d_psi() {

    if (!irrotational) {
	d_psi.set_etat_nondef() ; 
	return ; 
    }

    // Specific relativistic enthalpy		    ---> hhh
    //----------------------------------
    
    Tenseur hhh = exp(unsurc2 * ent) ;  // = 1 at the Newtonian limit
 
    //  Computation of W^i = - A^2 h Gamma_n B^i/N
    //----------------------------------------------

    Tenseur www = - a_car * hhh * gam_euler * bsn ; 
    
    
    // Constant value of W^i at the center of the star
    //-------------------------------------------------
    
    Tenseur v_orb(mp, 1, COV, ref_triad) ; 
    
    v_orb.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	v_orb.set(i) = www(i)(0, 0, 0, 0) ; 
    }
    
    v_orb.set_triad( *(www.get_triad()) ) ;     

    // Gradient of psi 
    //----------------

    Tenseur d_psi0 = psi0.gradient() ; 
    
    d_psi0.change_triad( ref_triad ) ; 

    d_psi = d_psi0 + v_orb ; 
    

    // C^1 continuation of d_psi outside the star
    //  (to ensure a smooth enthalpy field accross the stellar surface)
    // ----------------------------------------------------------------
    
    int nzm1 = mp.get_mg()->get_nzone() - 1 ;    

    if (d_psi0.get_etat() == ETATQCQ ) {
	d_psi.annule(nzet, nzm1) ;	 
	for (int i=0; i<3; i++) {
	    d_psi.set(i).va.set_base( d_psi0(i).va.base ) ; 
	    d_psi.set(i) = raccord_c1(d_psi(i), nzet) ; 
	}
    }
    else{ 
	d_psi.annule(nzm1) ;	 
    }   
} 


void Etoile_bin::fait_shift_auto() {

    Tenseur d_khi = khi_shift.gradient() ; 
    
    if (d_khi.get_etat() == ETATQCQ) { 
	d_khi.dec2_dzpuis() ;   // divide by r^2 in the external compactified
				// domain
    }
 
    // x_k dW^k/dx_i
    
    Tenseur x_d_w = skxk( w_shift.gradient() ) ;
    x_d_w.dec_dzpuis() ;
    
    double lambda = double(1) / double(3) ; 

    // The final computation is done component by component because
    // d_khi and x_d_w are covariant comp. whereas w_shift is
    // contravariant
    
    shift_auto.set_etat_qcq() ; 
    
    for (int i=0; i<3; i++) {
	shift_auto.set(i) = (lambda+2)/2./(lambda+1) * w_shift(i)
		- (lambda/2./(lambda+1))  * (d_khi(i) + x_d_w(i)) ;      
    }
    
    shift_auto.set_triad( *(w_shift.get_triad()) ) ; 
    
    // The final components of shift_auto should be those with respect
    // to the absolute frame (X,Y,Z) :
    
    shift_auto.change_triad( ref_triad ) ;  
    
} 


void Etoile_bin::relaxation(const Etoile_bin& star_jm1, double relax_ent, 
			    double relax_met, int mer, int fmer_met) {
				
    double relax_ent_jm1 = 1. - relax_ent ; 
    double relax_met_jm1 = 1. - relax_met ; 

    ent = relax_ent * ent + relax_ent_jm1 * star_jm1.ent ; 

    if ( (mer != 0) && (mer % fmer_met == 0)) {

	logn_auto = relax_met * logn_auto + relax_met_jm1 * star_jm1.logn_auto ;

	logn_auto_regu = relax_met * logn_auto_regu
	  + relax_met_jm1 * star_jm1.logn_auto_regu ;

	logn_auto_div = relax_met * logn_auto_div
	  + relax_met_jm1 * star_jm1.logn_auto_div ;

	d_logn_auto_div = relax_met * d_logn_auto_div
	  + relax_met_jm1 * star_jm1.d_logn_auto_div ;

	beta_auto = relax_met * beta_auto + relax_met_jm1 * star_jm1.beta_auto ;
	
	shift_auto = relax_met * shift_auto 
					+ relax_met_jm1 * star_jm1.shift_auto ;
	
    }

    del_deriv() ; 
    
    equation_of_state() ; 

}
