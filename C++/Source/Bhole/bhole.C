/*
 *  Methods of class Bhole
 *
 *   (see file bhole.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char bhole_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/02/13 16:40:25  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.4  2003/01/31 16:57:12  p_grandclement
 * addition of the member Cmp decouple used to compute the K_ij auto, once
 * the K_ij total is known
 *
 * Revision 1.3  2002/10/16 14:36:31  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/12/04 21:27:52  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.23  2001/06/28  15:08:27  eric
 * Ajout de la fonction area().
 *
 * Revision 2.22  2001/05/16  14:22:14  phil
 * modif init_bhole_seul
 *
 * Revision 2.21  2001/05/16  11:33:08  phil
 * ajout de init_bhole_seul
 *
 * Revision 2.20  2001/05/07  12:45:29  phil
 * *** empty log message ***
 *
 * Revision 2.19  2001/05/07  09:30:34  phil
 * *** empty log message ***
 *
 * Revision 2.18  2001/05/07  09:11:36  phil
 * *** empty log message ***
 *
 * Revision 2.17  2001/04/26  12:04:27  phil
 * *** empty log message ***
 *
 * Revision 2.16  2001/04/02  12:17:21  phil
 * *** empty log message ***
 *
 * Revision 2.15  2001/02/28  13:22:43  phil
 * vire kk_auto
 *
 * Revision 2.14  2001/01/29  14:30:37  phil
 * ajout type_rotation
 *
 * Revision 2.13  2001/01/10  09:31:31  phil
 * vire fait_kk_auto
 *
 * Revision 2.12  2001/01/09  13:38:26  phil
 * ajout memebre kk_auto
 *
 * Revision 2.11  2000/12/21  10:48:59  phil
 * retour a version 2.9
 *
 * Revision 2.9  2000/12/18  16:40:45  phil
 * *** empty log message ***
 *
 * Revision 2.8  2000/12/14  10:44:48  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.7  2000/12/04  14:29:47  phil
 * ajout de grad_n_tot
 *
 * Revision 2.6  2000/12/01  16:39:25  phil
 * correction impoetation grad_phi_comp
 *
 * Revision 2.5  2000/12/01  16:21:25  phil
 * modification init_bhole
 *
 * Revision 2.4  2000/12/01  16:17:42  phil
 * correction base import de grad_phi_comp
 *
 * Revision 2.3  2000/12/01  16:14:33  phil
 * *** empty log message ***
 *
 * Revision 2.2  2000/12/01  16:12:55  phil
 * *** empty log message ***
 *
 * Revision 2.1  2000/10/23  09:22:07  phil
 * rearrangement dans constructeurs.
 *
 * Revision 2.0  2000/10/20  09:18:47  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "tenseur.h"
#include "bhole.h"
#include "proto.h"
#include "utilitaires.h"
#include "etoile.h"
#include "graphique.h"

// Constructeur standard
Bhole::Bhole (Map_af& mpi) : mp(mpi),
			rayon ((mp.get_alpha())[0]), 
			omega(0), boost (new double[3]), regul(0), 
			n_auto(mpi), n_comp(mpi), n_tot(mpi), 
			psi_auto(mpi), psi_comp(mpi), psi_tot(mpi),  
			grad_n_tot (mpi, 1, COV, mpi.get_bvect_cart()),
			grad_psi_tot (mpi, 1, COV, mpi.get_bvect_cart()),
			shift_auto(mpi, 1, CON, mpi.get_bvect_cart()), 
			taij_auto(mpi, 2, CON,  mpi.get_bvect_cart()), 
			taij_comp(mpi, 2, CON,  mpi.get_bvect_cart()), 
			taij_tot (mpi, 2, CON,  mpi.get_bvect_cart()),
			tkij_auto (mpi, 2, CON,  mpi.get_bvect_cart()),  
			tkij_tot (mpi, 2, CON,  mpi.get_bvect_cart()), 
			decouple (mpi) {
     for (int i=0 ; i<3 ; i++)
	boost[i] = 0 ;
}


// Constructeur par recopie
Bhole::Bhole (const Bhole& source) :
    mp(source.mp), rayon(source.rayon), omega (source.omega), 
    boost (new double [3]), regul(source.regul), n_auto(source.n_auto),
    n_comp(source.n_comp), n_tot(source.n_tot), psi_auto(source.psi_auto), 
    psi_comp(source.psi_comp), psi_tot(source.psi_tot),
    grad_n_tot (source.grad_n_tot), grad_psi_tot (source.grad_psi_tot), 
    shift_auto(source.shift_auto), 
    taij_auto (source.taij_auto), 
    taij_comp(source.taij_comp), taij_tot(source.taij_tot), 
    tkij_auto(source.tkij_auto), tkij_tot(source.tkij_tot), decouple(source.decouple) {
    
    for (int i=0 ; i<3 ; i++)
	boost[i] = source.boost[i] ;
}

// Destructeur
Bhole::~Bhole() {
    delete [] boost ;
}

// Affectation
void Bhole::operator= (const Bhole& source) {
    
    assert (&mp == &source.mp) ;
    
    rayon = source.rayon ;
    omega = source.omega ;
    for (int i=0 ; i<3 ; i++)
	boost[i] = source.boost[i] ;
    regul = regul ;
    
    n_auto = source.n_auto ;
    n_comp = source.n_comp ;
    n_tot = source.n_tot ;
    
    psi_auto = source.psi_auto ;
    psi_comp = source.psi_comp ;
    psi_tot = source.psi_tot ;
    
    grad_n_tot = source.grad_n_tot ;
    grad_psi_tot = source.grad_psi_tot ;
    
    shift_auto = source.shift_auto ;
    
    taij_auto = source.taij_auto ;
    taij_comp = source.taij_comp ;
    taij_tot = source.taij_tot ;
    
    tkij_auto = source.tkij_auto ;
    tkij_tot = source.tkij_tot ;
    decouple = source.decouple ;
}

// Importe le lapse du compagnon (Bhole case)

void Bhole::fait_n_comp (const Bhole& comp) {
     // Alignes ou non ?
    double orientation = mp.get_rot_phi() ;
    assert ((orientation==0) || (orientation==M_PI)) ;
    double orientation_comp = comp.mp.get_rot_phi() ;
    assert ((orientation_comp==0) || (orientation_comp==M_PI)) ;
    int same_orient = (orientation == orientation_comp) ? 1 : -1 ;
    
    n_comp.set_etat_qcq() ;
    n_comp.set().import_symy(comp.n_auto()) ;
    n_comp.set_std_base() ;
    
    n_tot = n_comp + n_auto ;
    n_tot.set_std_base() ;
    
    Tenseur grad_comp (mp, 1, COV, mp.get_bvect_cart()) ;
    Tenseur auxi (comp.n_auto.gradient()) ;
    auxi.dec2_dzpuis() ;
    grad_comp.set_etat_qcq() ;
    grad_comp.set(0).import_symy(same_orient*auxi(0)) ;
    grad_comp.set(1).import_asymy(same_orient*auxi(1)) ;
    grad_comp.set(2).import_symy(auxi(2)) ;
    grad_comp.set_std_base() ;
    grad_comp.inc2_dzpuis() ;
    
    grad_n_tot = n_auto.gradient() + grad_comp ;
}

// Importe le facteur conforme du compagnon (Bhole case)

void Bhole::fait_psi_comp (const Bhole& comp) {
  
    // Alignes ou non ?
    double orientation = mp.get_rot_phi() ;
    assert ((orientation==0) || (orientation==M_PI)) ;
    double orientation_comp = comp.mp.get_rot_phi() ;
    assert ((orientation_comp==0) || (orientation_comp==M_PI)) ;
    int same_orient = (orientation == orientation_comp) ? 1 : -1 ;
    
    psi_comp.set_etat_qcq() ;
    psi_comp.set().import_symy(comp.psi_auto()) ;
    psi_comp.set_std_base() ;
    
    psi_tot = psi_comp + psi_auto ;
    psi_tot.set_std_base() ;
    
    Tenseur grad_comp (mp, 1, COV, mp.get_bvect_cart()) ;
    Tenseur auxi (comp.psi_auto.gradient()) ;
    auxi.dec2_dzpuis() ;
    grad_comp.set_etat_qcq() ;
    grad_comp.set(0).import_symy(same_orient*auxi(0)) ;
    grad_comp.set(1).import_asymy(same_orient*auxi(1)) ;
    grad_comp.set(2).import_symy(auxi(2)) ;
    grad_comp.set_std_base() ;
    grad_comp.inc2_dzpuis() ;
    
    grad_psi_tot = psi_auto.gradient() + grad_comp ;
}


// Importe le lapse du compagnon (NS case)
void Bhole::fait_n_comp (const Etoile_bin& comp) {
     // Alignes ou non ?
    double orientation = mp.get_rot_phi() ;
    assert ((orientation==0) || (orientation==M_PI)) ;
    double orientation_comp = comp.get_mp().get_rot_phi() ;
    assert ((orientation_comp==0) || (orientation_comp==M_PI)) ;
    int same_orient = (orientation == orientation_comp) ? 1 : -1 ;
    
    n_comp.set_etat_qcq() ;
    n_comp.set().import_symy(comp.get_logn_auto()()) ;
    n_comp.set_std_base() ;
    
    n_tot = n_comp + n_auto ;
    n_tot.set_std_base() ;
    
    Tenseur grad_comp (mp, 1, COV, mp.get_bvect_cart()) ;
    Tenseur auxi (comp.get_d_logn_auto()) ;
    auxi.dec2_dzpuis() ;
    grad_comp.set_etat_qcq() ;
    grad_comp.set(0).import_symy(same_orient*auxi(0)) ;
    grad_comp.set(1).import_asymy(same_orient*auxi(1)) ;
    grad_comp.set(2).import_symy(auxi(2)) ;
    grad_comp.set_std_base() ;
    grad_comp.inc2_dzpuis() ;
    
    grad_n_tot = n_auto.gradient() + grad_comp ;
}

// Importe le facteur conforme du compagnon (Bhole case)

void Bhole::fait_psi_comp (const Etoile_bin& comp) {
  
    // Alignes ou non ?
    double orientation = mp.get_rot_phi() ;
    assert ((orientation==0) || (orientation==M_PI)) ;
    double orientation_comp = comp.get_mp().get_rot_phi() ;
    assert ((orientation_comp==0) || (orientation_comp==M_PI)) ;
    int same_orient = (orientation == orientation_comp) ? 1 : -1 ;
    
    psi_comp.set_etat_qcq() ;
    psi_comp.set().import_symy(comp.get_beta_auto()()) ;
    psi_comp.set_std_base() ;
    
    psi_tot = psi_comp + psi_auto ;
    psi_tot.set_std_base() ;
    
    Tenseur grad_comp (mp, 1, COV, mp.get_bvect_cart()) ;
    Tenseur auxi (comp.get_d_beta_auto()) ;
    auxi.dec2_dzpuis() ;
    grad_comp.set_etat_qcq() ;
    grad_comp.set(0).import_symy(same_orient*auxi(0)) ;
    grad_comp.set(1).import_asymy(same_orient*auxi(1)) ;
    grad_comp.set(2).import_symy(auxi(2)) ;
    grad_comp.set_std_base() ;
    grad_comp.inc2_dzpuis() ;
    
    grad_psi_tot = psi_auto.gradient() + grad_comp ;
}

// Calcul Taij auto (nul sur H que si la regularisation a ete faite)
void Bhole::fait_taij_auto () {
    
    if (shift_auto.get_etat() == ETATZERO)
	taij_auto.set_etat_zero() ;
    else {
	Tenseur grad (shift_auto.gradient()) ;
	Tenseur trace (mp) ;
	trace = grad(0, 0)+grad(1, 1)+grad(2, 2) ;
    
	taij_auto.set_etat_qcq() ;
	taij_auto.set_std_base() ;
	for (int i=0 ; i<3 ; i++) {
	    for (int j=i+1 ; j<3 ; j++)
		taij_auto.set(i, j) = grad(i, j)+grad(j, i) ;
	    taij_auto.set(i, i) = 2*grad(i, i) -2./3.*trace() ;
	}
    
	for (int i=0 ; i<3 ; i++)
	  for (int j=0 ; j<3 ; j++)
		taij_auto.set(i, j).raccord(1) ;
    }
}

// Regularise le shift, untilisant le compagnon.
void Bhole::regularise_shift (const Tenseur& shift_comp) {
    regul = regle (shift_auto, shift_comp, omega) ;
}

//Initialise a Schwartz
// Compagnon a 0
void Bhole::init_bhole () {
    
    Cmp auxi(mp) ;
    
    auxi = 1./2.-rayon/mp.r ;
    auxi.annule(0);
    auxi.set_dzpuis(0) ;
    n_auto = auxi;
    n_auto.set_std_base() ;
    n_auto.set().raccord(1) ;
    n_comp =0.5;
    n_comp.set_std_base() ;
    n_tot = n_comp+n_auto ;
  
    auxi = 0.5+rayon/mp.r ;
    auxi.annule(0);
    auxi.set_dzpuis(0) ;
    psi_auto = auxi;
    psi_auto.set_std_base() ;
    psi_auto.set().raccord(1) ;
    psi_comp = 0.5;
    psi_comp.set_std_base() ;
    psi_tot = psi_comp+psi_auto ;
    
    grad_n_tot = n_tot.gradient() ;
    grad_psi_tot = psi_tot.gradient() ;
    
    shift_auto.set_etat_zero() ;
    taij_auto.set_etat_zero();
    taij_comp.set_etat_zero();
    taij_tot.set_etat_zero() ;
    tkij_auto.set_etat_zero() ;
    tkij_tot.set_etat_zero();
    decouple.set_etat_zero() ;
}

void Bhole::init_bhole_seul () {
    
    Cmp auxi(mp) ;
    
    auxi = (1-rayon/mp.r)/(1+rayon/mp.r) ;
    auxi.annule(0);
    auxi.set_val_inf(1) ;
    auxi.set_dzpuis(0) ;
    n_auto = auxi;
    n_auto.set_std_base() ;
    n_auto.set().raccord(1) ;
    n_comp.set_etat_zero();
    n_tot.set_etat_zero() ;
    
    auxi = 1+rayon/mp.r ;
    auxi.annule(0);
    auxi.set_val_inf(1) ;
    auxi.set_dzpuis(0) ;
    psi_auto = auxi;
    psi_auto.set_std_base() ;
    psi_auto.set().raccord(1) ;
    psi_comp.set_etat_zero() ;
    psi_tot.set_etat_zero() ;
    
    grad_n_tot = n_tot.gradient() ;
    grad_psi_tot = psi_tot.gradient() ;
    
    shift_auto.set_etat_zero() ;
    taij_auto.set_etat_zero();
    taij_comp.set_etat_zero();
    taij_tot.set_etat_zero() ;
    tkij_auto.set_etat_zero() ;
    tkij_tot.set_etat_zero();
    decouple.set_etat_zero() ;
}

// Sauvegarde dans fichier
void Bhole::sauve (FILE* fich) const {
    
    fwrite_be (&omega, sizeof(double), 1, fich) ;
    fwrite_be (boost, sizeof(double), 3, fich) ;
    fwrite_be (&regul, sizeof(double), 1, fich) ;
    
    n_auto.sauve(fich) ;
    psi_auto.sauve(fich) ;
    shift_auto.sauve(fich) ;
}

//Constructeur par fichier
Bhole::Bhole(Map_af& mpi, FILE* fich) :
	    mp(mpi), rayon ((mp.get_alpha())[0]), 
	    boost (new double[3]), n_auto(mpi), n_comp(mpi), n_tot(mpi), 
			psi_auto(mpi), psi_comp(mpi), psi_tot(mpi), 
			grad_n_tot (mpi, 1, COV, mpi.get_bvect_cart()), 
			grad_psi_tot (mpi, 1, COV, mpi.get_bvect_cart()), 
			shift_auto(mpi, 1, CON, mpi.get_bvect_cart()), 
			taij_auto(mpi, 2, CON,  mpi.get_bvect_cart()), 
			taij_comp(mpi, 2, CON,  mpi.get_bvect_cart()), 
			taij_tot (mpi, 2, CON,  mpi.get_bvect_cart()), 
			tkij_auto (mpi, 2, CON,  mpi.get_bvect_cart()) , 
			tkij_tot (mpi, 2, CON,  mpi.get_bvect_cart()) , 
			decouple (mpi) {

	fread_be(&omega, sizeof(double), 1, fich) ;
	fread_be(boost, sizeof(double), 3, fich) ;
	fread_be(&regul, sizeof(double), 1, fich) ;
	
	Tenseur n_auto_file (mp, fich) ;
	n_auto = n_auto_file ;
	
	Tenseur psi_auto_file (mp, fich) ;
	psi_auto = psi_auto_file ;
	
	Tenseur shift_auto_file (mp, mp.get_bvect_cart(), fich) ;
	shift_auto = shift_auto_file ;
	
	grad_n_tot.set_etat_qcq() ;
	grad_psi_tot.set_etat_zero() ;
	
	taij_auto.set_etat_zero();
	taij_comp.set_etat_zero();
	taij_tot.set_etat_zero() ;
	
	tkij_auto.set_etat_zero() ;
	tkij_tot.set_etat_zero() ;
	decouple.set_etat_zero() ;
}


		    //---------------------------------//
		    //	Computation of throat area     //
		    //---------------------------------//

double Bhole::area() const {
    
    Cmp integrant( pow(psi_tot(), 4) ) ;	// Psi^4
    
    integrant.std_base_scal() ;
    integrant.raccord(1) ;
    
    return mp.integrale_surface(integrant, rayon) ;
    
}

		    
