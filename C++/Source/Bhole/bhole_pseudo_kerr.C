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


char bhole_pseudo_kerr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/03/25 10:28:57  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.3  2003/10/03 15:58:43  j_novak
 * Cleaning of some headers
 *
 * Revision 1.2  2002/10/16 14:36:32  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2001/05/16  14:51:36  phil
 * correction calcul de kij
 *
 * Revision 2.8  2001/05/07  12:24:30  phil
 * correction calcul de J
 *
 * Revision 2.7  2001/05/07  09:28:33  phil
 * *** empty log message ***
 *
 * Revision 2.6  2001/02/12  15:36:58  phil
 * ajout calcul de J a l infini
 *
 * Revision 2.5  2000/12/14  14:11:54  phil
 * correction cl sur psi
 *
 * Revision 2.4  2000/12/14  12:41:54  phil
 * on met les bases dans les sources
 *
 * Revision 2.3  2000/12/14  10:57:47  phil
 * corections diverses et sans importances
 *
 * Revision 2.2  2000/12/14  10:45:00  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.1  2000/11/03  12:56:33  phil
 * ajout de const
 *
 * Revision 2.0  2000/10/20  09:19:09  phil
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

//Resolution pour le lapse pour 1 seul trou
void Bhole::solve_lapse_seul (double relax) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "Resolution LAPSE" << endl ;
    
    // Pour la relaxation ...
    Cmp lapse_old (n_auto()) ;
    Tenseur auxi (flat_scalar_prod(tkij_auto, tkij_auto)) ;
    Tenseur kk (mp) ;
    kk = 0 ;
    Tenseur work(mp) ;
    work.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++) {
	work.set() = auxi(i, i) ;
	kk = kk + work ;
	}
    
    // La source
    Cmp source 
    (-2*flat_scalar_prod(psi_auto.gradient(), n_auto.gradient())()/psi_auto()
	+pow(psi_auto(), 4.)*n_auto()*kk()) ;
    source.std_base_scal() ;
     
    // On resout pour N-1 :
    Valeur limite (mp.get_mg()->get_angu()) ;
    limite = -1 ;
    limite.std_base_scal() ;
    
    Cmp soluce (source.poisson_dirichlet(limite, 0)) ;
    soluce = soluce + 1 ;   // Permet de trouver N
    soluce.raccord(1) ;
    
    n_auto.set() = relax*soluce + (1-relax)*lapse_old ;
}


// Resolution sur Psi :
void Bhole::solve_psi_seul (double relax) {
    
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "Resolution PSI" << endl ;
    
    Cmp psi_old (psi_auto()) ;
    Tenseur auxi (flat_scalar_prod(tkij_auto, tkij_auto)) ;
    Tenseur kk (mp) ;
    kk = 0 ;
    Tenseur work(mp) ;
    work.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++) {
	work.set() = auxi(i, i) ;
	kk = kk + work ;
	}
    
    // La source :
    Cmp source (-pow(psi_auto(), 5.)*kk()/8.) ;
    source.std_base_scal() ;
    
    // La condition limite de type neumann :
    int np = mp.get_mg()->get_np(1) ;
    int nt = mp.get_mg()->get_nt(1) ;
    Valeur limite (mp.get_mg()->get_angu()) ;
    limite = 1 ;
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    limite.set(0, k, j, 0) = -0.5/rayon*psi_auto()(1, k, j, 0) ;
    limite.std_base_scal() ;
    
    Cmp soluce (source.poisson_neumann(limite, 0)) ;
    soluce = soluce + 1 ;
    soluce.raccord(1) ;
    
    psi_auto.set() = relax*soluce + (1-relax)*psi_old ;
    
}


// Le shift. Processus iteratif pour cause de CL.
void Bhole::solve_shift_seul (double precision, double relax) {
    
    assert (precision > 0) ;
    assert ((relax>0) && (relax<=1)) ;
    
    cout << "resolution SHIFT" << endl ;
    
    Tenseur shift_old (shift_auto) ;
    
    Tenseur source (-6*flat_scalar_prod(taij_auto, psi_auto.gradient())/psi_auto
		    + 2*flat_scalar_prod(tkij_auto, n_auto.gradient())) ;
    source.set_std_base() ;
    
    // On verifie si les 3 composantes ne sont pas nulles :
    if (source.get_etat() == ETATQCQ) {
	int indic = 0 ;
	for (int i=0 ; i<3 ; i++)
	    if (source(i).get_etat() == ETATQCQ)
		indic = 1 ;
	if (indic ==0)
	    source.set_etat_zero() ;
    }
    
    // On filtre les hautes frequences pour raison de stabilite :
    if (source.get_etat() == ETATQCQ)
	for (int i=0 ; i<3 ; i++)
    	    source.set(i).filtre(4) ;
    
    
    // On determine les conditions limites en fonction du omega et du boost :
    int np = mp.get_mg()->get_np(1) ;
    int nt = mp.get_mg()->get_nt(1) ;
    
    Mtbl x_mtbl (mp.get_mg()) ;
    x_mtbl.set_etat_qcq() ;
    Mtbl y_mtbl (mp.get_mg()) ;
    y_mtbl.set_etat_qcq() ;
    x_mtbl = mp.x ;
    y_mtbl = mp.y ;
    
    // Les bases pour les conditions limites :
    Base_val** bases = mp.get_mg()->std_base_vect_cart() ;
    
    Valeur lim_x (mp.get_mg()->get_angu()) ;
    lim_x = 1 ;
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    lim_x.set(0, k, j, 0) = omega*y_mtbl(1, k, j, 0)-boost[0] ;
    lim_x.base = *bases[0] ;
    
    Valeur lim_y (mp.get_mg()->get_angu()) ;
    lim_y = 1 ;
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    lim_y.set(0, k, j, 0) = - omega*x_mtbl(1, k, j, 0)-boost[1] ;
    lim_y.base = *bases[1] ;
   
    Valeur lim_z (mp.get_mg()->get_angu()) ;
    lim_z = 1 ;
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    lim_z.set(0, k, j, 0) = -boost[2] ;
    lim_z.base = *bases[2] ;
     
     // On n'en a plus besoin
    for (int i=0 ; i<3 ; i++)
	delete bases[i] ;
    delete [] bases ;
    
    // On resout :
    poisson_vect_frontiere(1./3., source, shift_auto, lim_x, lim_y, 
				    lim_z, 0, precision, 20) ;
      
    shift_auto = relax*shift_auto + (1-relax)*shift_old ;
}


// La regularisation si un seul trou noir :
void Bhole::regularise_seul () {
    
    // Le vecteur B (non tournant et non boostant)
    Tenseur tbi (shift_auto) ;
    for (int i=0 ; i<3 ; i++) {
	tbi.set(i).va.coef_i() ;
	tbi.set(i).va.set_etat_c_qcq() ;
	}
	
    for (int i=0 ; i<3 ; i++)
	shift_auto(i).va.coef_i() ;
	
    tbi.set(0) = *shift_auto(0).va.c - omega*shift_auto.get_mp()->y + boost[0];
    tbi.set(1) = *shift_auto(1).va.c + omega*shift_auto.get_mp()->x + boost[1];
    tbi.set(2) = *shift_auto(2).va.c + boost[2];
    tbi.set_std_base() ;
    
    // On evite soucis a l'infini (on a besoin que de la valeur sur horizon
    tbi.set(0).annule(mp.get_mg()->get_nzone()-1) ;
    tbi.set(1).annule(mp.get_mg()->get_nzone()-1) ;
      
    Tenseur derive_r (mp, 1, CON, mp.get_bvect_cart()) ;
    derive_r.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++)
	derive_r.set(i) = tbi(i).dsdr() ;
    
    Valeur val_hor (mp.get_mg()) ;
    Valeur fonction_radiale (mp.get_mg()) ;
    Cmp enleve (mp) ;
    
    double erreur = 0 ;
    int nz = mp.get_mg()->get_nzone() ;
    int np = mp.get_mg()->get_np(1) ;
    int nt = mp.get_mg()->get_nt(1) ;
    int nr = mp.get_mg()->get_nr(1) ;
    
    double r_0 = mp.val_r(1, -1, 0, 0) ;
    double r_1 = mp.val_r(1, 1, 0, 0) ;
    
    for (int comp=0 ; comp<3 ; comp++) {
	val_hor.annule_hard() ;
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++)
		    val_hor.set(1, k, j, i) = derive_r(comp)(1, k, j, 0) ;
	
	fonction_radiale = pow(r_1-mp.r, 3.)* (mp.r-r_0)/pow(r_1-r_0, 3.) ;
	fonction_radiale.annule(0) ;
	fonction_radiale.annule(2, nz-1) ;
    
	enleve = fonction_radiale*val_hor ;
	enleve.va.base = shift_auto(comp).va.base ;
	
	// Ca devrai annuler la derivee de B sur H et donc rendre K regulier
	Cmp copie (shift_auto(comp)) ;
	shift_auto.set(comp) = shift_auto(comp)-enleve ;
	
	assert (shift_auto(comp).check_dzpuis(0)) ;
	
	// On regarde l'intensite de la correction si non nul !
	    Tbl norm (norme(shift_auto(comp))) ;
	    if (norm(1) > 1e-5) {
		Tbl diff (diffrelmax (copie, shift_auto(comp))) ;
		if (erreur<diff(1)) ;
		erreur = diff(1) ;
	    }
	}
    regul = erreur ;
}

// On calcul Kij sachant que la regulatisation sur le shift doit 
// etre faite avant.
void Bhole::fait_tkij () {
    
    fait_taij_auto() ;
    
    Cmp lapse_non_sing (division_xpun(n_auto(), 0)) ;
    Cmp auxi (mp) ;
    
    tkij_auto.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++)
	for (int j=i ; j<3 ; j++) {
	    auxi = taij_auto(i, j) ;
	    auxi = division_xpun (auxi, 0) ;
	    tkij_auto.set(i, j) = auxi/lapse_non_sing/2. ;
	    tkij_auto.set(i, j).raccord(1) ;
	}
}

// Calcul masse ADM via valeur asymptotiques
double Bhole::masse_adm_seul () const {
    
    Cmp integrant (psi_auto().dsdr()) ;
    double masse = mp.integrale_surface_infini (integrant) ;
    masse /= -2*M_PI ;
    return masse ;
}

double Bhole::masse_komar_seul() const {
    
    Cmp integrant (n_auto().dsdr()) ;
    double masse = mp.integrale_surface_infini (integrant) ;
    masse /= 4*M_PI ;
    return masse ;
}

// Calcul du moment angulaire via integrale a l infini
// Non valable si le boost != 0 ;

double Bhole::moment_seul_inf() const {

    // On verifie si le boost est bien nul :
    double indic = 0 ;
    for (int i=0 ; i<3 ; i++)
	if (boost[i] != 0)
	    indic = 1 ;
    if (indic == 1) {
	cout << "Calcul du moment non valable pour un boost != 0" << endl ;
	abort() ;
    }
    
    if (omega == 0)
	return 0 ;
    else {
	Tenseur vecteur_un (mp, 1, CON, mp.get_bvect_cart()) ;
	vecteur_un.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_un.set(i) = tkij_auto(0, i) ;
	vecteur_un.change_triad (mp.get_bvect_spher()) ;
	Cmp integrant_un (vecteur_un(0)) ;
	integrant_un.mult_r_zec() ;
	
	Tenseur vecteur_deux (mp, 1, CON, mp.get_bvect_cart()) ;
	vecteur_deux.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_deux.set(i) = tkij_auto(1, i) ;
	vecteur_deux.change_triad (mp.get_bvect_spher()) ;
	Cmp integrant_deux (vecteur_deux(0)) ;
	integrant_deux.mult_r_zec() ;
	
	// On multiplie par y et x :
	integrant_un.va = integrant_un.va.mult_st() ;
	integrant_un.va = integrant_un.va.mult_sp() ;
	
	integrant_deux.va = integrant_deux.va.mult_st() ;
	integrant_deux.va = integrant_deux.va.mult_cp() ;
	
	double moment = mp.integrale_surface_infini (-integrant_un+integrant_deux) ;
	moment /= 8*M_PI ;
	return moment ;
	}
}

// Calcul du moment angulaire via integrale sur l horizon
// Non valable si le boost != 0 ;

double Bhole::moment_seul_hor() const {

    // On verifie si le boost est bien nul :
    double indic = 0 ;
    for (int i=0 ; i<3 ; i++)
	if (boost[i] != 0)
	    indic = 1 ;
    if (indic == 1) {
	cout << "Calcul du moment non valable pour un boost != 0" << endl ;
	abort() ;
    }
    
    if (omega == 0)
	return 0 ;
    else {
	// Integrale sur l horizon :
	Cmp xa (mp) ;
	xa = mp.xa ;
	xa.std_base_scal() ;
	
	Cmp ya (mp) ;
	ya = mp.ya ;
	ya.std_base_scal() ;
	
	Tenseur vecteur (mp, 1, CON, mp.get_bvect_cart()) ;
	vecteur.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur.set(i) = (-ya*tkij_auto(0, i)+xa * tkij_auto(1, i)) ;
	vecteur.set_std_base() ;
	vecteur.annule(mp.get_mg()->get_nzone()-1) ;
	vecteur.change_triad (mp.get_bvect_spher()) ;
	
	Cmp integrant (pow(psi_auto(), 6)*vecteur(0)) ;
	integrant.std_base_scal() ;
	double moment = mp.integrale_surface (integrant, rayon)/8/M_PI ;
	return moment ;
	}
}

