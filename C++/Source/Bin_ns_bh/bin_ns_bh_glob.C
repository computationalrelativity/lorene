/*
 *   Copyright (c) 2005 Philippe Grandclement
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


char bin_ns_bh_glob_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2006/06/23 07:09:24  p_grandclement
 * Addition of spinning black hole
 *
 * Revision 1.4  2006/06/01 12:47:52  p_grandclement
 * update of the Bin_ns_bh project
 *
 * Revision 1.3  2005/12/01 12:59:10  p_grandclement
 * Files for bin_ns_bh project
 *
 * Revision 1.2  2005/11/30 11:09:06  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.1  2005/08/29 15:10:15  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 *
 * $Header$
 *
 */



//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tenseur.h"
#include "bhole.h"
#include "bin_ns_bh.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"
#include "unites.h"
#include "metrique.h"

double Bin_ns_bh::adm_systeme() const {
    Cmp der_un (hole.psi_auto().dsdr()) ;
    
    Map_af auxi_mp (star.get_mp()) ;
    Cmp der_deux (star.confpsi_auto().dsdr()) ;
    
    double masse = hole.mp.integrale_surface_infini(der_un) + 
	auxi_mp.integrale_surface_infini(der_deux) ;
    
    masse /= -2*M_PI ;
    return masse ;
}

double Bin_ns_bh::adm_systeme_volume() const {
 
    using namespace Unites ;

    Tenseur auxi_bh (flat_scalar_prod(hole.tkij_tot, hole.tkij_auto)) ;
    Tenseur kk_bh (hole.mp) ;
    kk_bh = 0 ;
    Tenseur work_bh(hole.mp) ;
    work_bh.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++) {
	work_bh.set() = auxi_bh(i, i) ;
	kk_bh = kk_bh + work_bh ;
	}
    Cmp integ_bh (pow(hole.psi_tot(), 5.)*kk_bh()) ;
    integ_bh.annule(0) ;
    integ_bh.std_base_scal() ;
    
    Cmp integ_hor (hole.psi_tot()) ;
    
    Tenseur auxi_ns (flat_scalar_prod(star.tkij_tot, star.tkij_auto)) ;
    Tenseur kk_ns (star.mp) ;
    kk_ns = 0 ;
    Tenseur work_ns(star.mp) ;
    work_ns.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++) {
	work_ns.set() = auxi_ns(i, i) ;
	kk_ns = kk_ns + work_ns ;
	}
    Cmp integ_ns (pow(star.confpsi_comp() + star.confpsi_auto(), 5.)*kk_ns()) ;
    integ_ns.std_base_scal() ;
    
    Cmp integ_matter (pow(star.confpsi_comp()+star.confpsi_auto(), 5.)*star.ener_euler()) ;
    integ_matter.std_base_scal() ;
    
    double masse = (integ_bh.integrale()+integ_ns.integrale())/16/M_PI  + 
    	hole.mp.integrale_surface(integ_hor, hole.rayon)/hole.rayon/4/M_PI 
	+ integ_matter.integrale()*ggrav ;
	
    return masse ;
}

double Bin_ns_bh::komar_systeme() const {
    Cmp der_un (hole.n_auto().dsdr()) ;
    
    Map_af auxi_mp (star.get_mp()) ;
    Cmp der_deux (star.n_auto().dsdr()) ;

    double masse = hole.mp.integrale_surface_infini(der_un) + 
	auxi_mp.integrale_surface_infini(der_deux) ;

    masse /= 4*M_PI ;
    
    return masse ;
}

double Bin_ns_bh::viriel() const {
    double adm = adm_systeme() ;
    double komar = komar_systeme() ;
    
    return (adm-komar)/adm ;
}

double Bin_ns_bh::moment_systeme_inf() const {
    
    if (omega == 0)
	return 0 ;
    else {
	
	// On construit une grille et un mapping auxiliaire :
	int nzones = hole.mp.get_mg()->get_nzone() ;
	double* bornes = new double [nzones+1] ;
	double courant = fabs(hole.mp.get_ori_x()-star.mp.get_ori_x())+1 ;
	for (int i=nzones-1 ; i>0 ; i--) {
	    bornes[i] = courant ;
	    courant /= 2. ;
	    }
	bornes[0] = 0 ;
	bornes[nzones] = __infinity ;
	
	Map_af mapping (*hole.mp.get_mg(), bornes) ;
	
	delete [] bornes ; 
	
	// On construit k_total dessus :
	Tenseur_sym k_total (mapping, 2, CON, mapping.get_bvect_cart()) ;
	k_total.set_etat_qcq() ;
	
	Tenseur shift_un (mapping, 1, CON, mapping.get_bvect_cart()) ;
	shift_un.set_etat_qcq() ;
	
	Tenseur shift_deux (mapping, 1, CON, mapping.get_bvect_cart()) ;
	shift_deux.set_etat_qcq() ;
	
	shift_un.set_triad (*hole.shift_auto.get_triad()) ;
	shift_un.set(0).import_asymy(hole.shift_auto(0)) ;
	shift_un.set(1).import_symy(hole.shift_auto(1)) ;
	shift_un.set(2).import_asymy(hole.shift_auto(2)) ;
	shift_un.change_triad (mapping.get_bvect_cart()) ;
	
	shift_deux.set_triad (*star.shift_auto.get_triad()) ;
	shift_deux.set(0).import_asymy(star.shift_auto(0)) ;
	shift_deux.set(1).import_symy(star.shift_auto(1)) ;
	shift_deux.set(2).import_asymy(star.shift_auto(2)) ;
	shift_deux.change_triad(mapping.get_bvect_cart()) ;
	
	Tenseur shift_tot (shift_un+shift_deux) ;
	shift_tot.set_std_base() ;
	shift_tot.annule(0, nzones-2) ;
	// On enleve les residus
	shift_tot.inc2_dzpuis() ;
	shift_tot.dec2_dzpuis() ;
	
	Tenseur grad (shift_tot.gradient()) ;
	Tenseur trace (grad(0, 0)+grad(1, 1)+grad(2, 2)) ;
	for (int i=0 ; i<3 ; i++) {
	    k_total.set(i, i) = grad(i, i)-trace()/3. ;
	    for (int j=i+1 ; j<3 ; j++)
		k_total.set(i, j) = (grad(i, j)+grad(j, i))/2. ;
	    }	
	
	for (int lig=0 ; lig<3 ; lig++)
	   for (int col=lig ; col<3 ; col++)
		k_total.set(lig, col).mult_r_zec() ;
	
	Tenseur vecteur_un (mapping, 1, CON, mapping.get_bvect_cart()) ;
	vecteur_un.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_un.set(i) = k_total(0, i) ;
	vecteur_un.change_triad (mapping.get_bvect_spher()) ;
	Cmp integrant_un (vecteur_un(0)) ;
	
	Tenseur vecteur_deux (mapping, 1, CON, mapping.get_bvect_cart()) ;
	vecteur_deux.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_deux.set(i) = k_total(1, i) ;
	vecteur_deux.change_triad (mapping.get_bvect_spher()) ;
	Cmp integrant_deux (vecteur_deux(0)) ;
	
	// On multiplie par y et x :
	integrant_un.va = integrant_un.va.mult_st() ;
	integrant_un.va = integrant_un.va.mult_sp() ;
	
	integrant_deux.va = integrant_deux.va.mult_st() ;
	integrant_deux.va = integrant_deux.va.mult_cp() ;
	
	double moment = mapping.integrale_surface_infini (-integrant_un+integrant_deux) ;
	
	moment /= 8*M_PI ;
	return moment ;
	}
}

double Bin_ns_bh::moment_systeme_hor() const {
    
    using namespace Unites ;

    if (omega == 0)
	return 0 ;
    else {	
	//Contribution du trou noir :
	Cmp xa_bh (hole.mp) ;
	xa_bh = hole.mp.xa ;
	xa_bh.std_base_scal() ;
	
	Cmp ya_bh (hole.mp) ;
	ya_bh = hole.mp.ya ;
	ya_bh.std_base_scal() ;
	
	Tenseur vecteur_bh (hole.mp, 1, CON, hole.mp.get_bvect_cart()) ;
	vecteur_bh.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_bh.set(i) = (-ya_bh*hole.tkij_tot(0, i)+
		xa_bh * hole.tkij_tot(1, i)) ;
	vecteur_bh.set_std_base() ;
	vecteur_bh.annule(hole.mp.get_mg()->get_nzone()-1) ;
	vecteur_bh.change_triad (hole.mp.get_bvect_spher()) ;
	
	Cmp integrant_bh (pow(hole.psi_tot(), 6)*vecteur_bh(0)) ;
	integrant_bh.std_base_scal() ;
	double moment_bh = hole.mp.integrale_surface
	    (integrant_bh, hole.rayon)/8/M_PI ;
	
	// Contribution de l'étoile :
        Cmp xa_ns (star.mp) ;
	xa_ns = star.mp.xa ;
	xa_ns.std_base_scal() ;
	
	Cmp ya_ns (star.mp) ;
	ya_ns = star.mp.ya ;
	ya_ns.std_base_scal() ;
	
	Cmp integrant_ns (pow(star.confpsi_auto()+star.confpsi_comp(), 10)*(star.ener_euler()+star.press())*
	              (xa_ns*star.u_euler(1) - ya_ns*star.u_euler(0))) ;
	integrant_ns.std_base_scal() ;
	
	double moment_ns = integrant_ns.integrale() * ggrav ;
	return moment_ns + moment_bh ;
	}
}

Tbl Bin_ns_bh::linear_momentum_systeme_inf() const {
    
	Tbl res (3) ;
	res.set_etat_qcq() ;
	
	// On construit une grille et un mapping auxiliaire :
	int nzones = hole.mp.get_mg()->get_nzone() ;
	double* bornes = new double [nzones+1] ;
	double courant = fabs(hole.mp.get_ori_x()-star.mp.get_ori_x())+1 ;
	for (int i=nzones-1 ; i>0 ; i--) {
	    bornes[i] = courant ;
	    courant /= 2. ;
	    }
	bornes[0] = 0 ;
	bornes[nzones] = __infinity ;
	
	Map_af mapping (*hole.mp.get_mg(), bornes) ;
	
	delete [] bornes ; 
	
	// On construit k_total dessus :
	Tenseur_sym k_total (mapping, 2, CON, mapping.get_bvect_cart()) ;
	k_total.set_etat_qcq() ;
	
	Tenseur shift_un (mapping, 1, CON, mapping.get_bvect_cart()) ;
	shift_un.set_etat_qcq() ;
	
	Tenseur shift_deux (mapping, 1, CON, mapping.get_bvect_cart()) ;
	shift_deux.set_etat_qcq() ;
	
	shift_un.set_triad (*hole.shift_auto.get_triad()) ;
	shift_un.set(0).import_asymy(hole.shift_auto(0)) ;
	shift_un.set(1).import_symy(hole.shift_auto(1)) ;
	shift_un.set(2).import_asymy(hole.shift_auto(2)) ;
	shift_un.change_triad (mapping.get_bvect_cart()) ;
	
	shift_deux.set_triad (*star.shift_auto.get_triad()) ;
	shift_deux.set(0).import_asymy(star.shift_auto(0)) ;
	shift_deux.set(1).import_symy(star.shift_auto(1)) ;
	shift_deux.set(2).import_asymy(star.shift_auto(2)) ;
	shift_deux.change_triad(mapping.get_bvect_cart()) ;
	
	shift_un.set_std_base() ;
	shift_deux.set_std_base() ;
	
	Tenseur shift_tot (shift_un+shift_deux) ;
	shift_tot.set_std_base() ;
	shift_tot.annule(0, nzones-2) ;
	
	Cmp compy (shift_tot(1)) ;
	compy.mult_r_zec() ;
	
	int nr = mapping.get_mg()->get_nr(nzones-1) ;
	int nt = mapping.get_mg()->get_nt(nzones-1) ;
	int np = mapping.get_mg()->get_np(nzones-1) ;
	Tbl val_inf (nt*np) ;
	val_inf.set_etat_qcq() ;
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
	        val_inf.set(k*nt + j) = fabs(compy (nzones-1, k, j, nr-1)) ;
	
	Tenseur grad (shift_tot.gradient()) ;
	Tenseur trace (grad(0, 0)+grad(1, 1)+grad(2, 2)) ;
	for (int i=0 ; i<3 ; i++) {
	    k_total.set(i, i) = grad(i, i)-trace()/3. ;
	    for (int j=i+1 ; j<3 ; j++)
		k_total.set(i, j) = (grad(i, j)+grad(j, i))/2. ;
	    }
	
	for (int comp=0 ; comp<3 ; comp++) {
		Tenseur vecteur (mapping, 1, CON, mapping.get_bvect_cart()) ;
		vecteur.set_etat_qcq() ;
		for (int i=0 ; i<3 ; i++)
	    		vecteur.set(i) = k_total(i, comp) ;
		vecteur.change_triad (mapping.get_bvect_spher()) ;
		Cmp integrant (vecteur(0)) ;
	
		res.set(comp) = mapping.integrale_surface_infini (integrant)/8/M_PI ;
	}
	return res ;
}

double Bin_ns_bh::distance_propre_axe_bh (const int nr) const {
   
    double x_bh = hole.mp.get_ori_x() + hole.rayon ;
    
    // Les coefficients du changement de variable :
    double pente = -2./x_bh ;
    double constante = - 1. ;
    
    double ksi ; // variable d'integration.
    double xabs ; // x reel.
    double air_un, theta_un, phi_un ; // coordonnee spheriques 1
    double air_deux, theta_deux, phi_deux ; // coordonnee spheriques 2
    
    double* coloc = new double[nr] ;
    double* coef = new double[nr] ;
    int* deg = new int[3] ;
    deg[0] = 1 ; deg[1] = 1 ; deg[2] = nr ;
    
    for (int i=0 ; i<nr ; i++) {
	ksi = -cos (M_PI*i/(nr-1)) ;
	xabs = (ksi+constante)/pente ;
	
	hole.mp.convert_absolute (xabs, 0, 0, air_un, theta_un, phi_un) ;
	star.mp.convert_absolute (xabs, 0, 0, air_deux, theta_deux, phi_deux) ;
	
	coloc[i] = pow (hole.psi_auto().val_point (air_un, theta_un, phi_un) +
		   star.confpsi_auto().val_point (air_deux, theta_deux, phi_deux), 2.) ;
    }
    
    // On prend les coefficients de la fonction
    cfrcheb(deg, deg, coloc, deg, coef) ;
    
    // On integre
    double* som = new double[nr] ;
    som[0] = 2 ;
    for (int i=2 ; i<nr ; i+=2)
	som[i] = 1./(i+1)-1./(i-1) ;
    for (int i=1 ; i<nr ; i+=2)
	som[i] = 0 ;
    
    double res = 0 ;
    for (int i=0 ; i<nr ; i++)
	res += som[i]*coef[i] ;
	
    res /= pente ;
    
    delete [] deg ;
    delete [] coef ;
    delete [] coloc ;
    delete [] som ;
    
    return res ;
}

double Bin_ns_bh::distance_propre_axe_ns (const int nr) const {
   
    double x_ns = star.mp.get_ori_x() - star.mp.val_r (star.nzet, -1, M_PI/2, M_PI) ;
    
    // Les coefficients du changement de variable :
    double pente = 2./x_ns ;
    double constante = - 1. ;
    
    double ksi ; // variable d'integration.
    double xabs ; // x reel.
    double air_un, theta_un, phi_un ; // coordonnee spheriques 1
    double air_deux, theta_deux, phi_deux ; // coordonnee spheriques 2
    
    double* coloc = new double[nr] ;
    double* coef = new double[nr] ;
    int* deg = new int[3] ;
    deg[0] = 1 ; deg[1] = 1 ; deg[2] = nr ;
    
    for (int i=0 ; i<nr ; i++) {
	ksi = -cos (M_PI*i/(nr-1)) ;
	xabs = (ksi-constante)/pente ;
	
	hole.mp.convert_absolute (xabs, 0, 0, air_un, theta_un, phi_un) ;
	star.mp.convert_absolute (xabs, 0, 0, air_deux, theta_deux, phi_deux) ;
	
	coloc[i] = pow (hole.psi_auto().val_point (air_un, theta_un, phi_un) +
		   star.confpsi_auto().val_point (air_deux, theta_deux, phi_deux), 2.) ;
    }
    
    // On prend les coefficients de la fonction
    cfrcheb(deg, deg, coloc, deg, coef) ;
    
    // On integre
    double* som = new double[nr] ;
    som[0] = 2 ;
    for (int i=2 ; i<nr ; i+=2)
	som[i] = 1./(i+1)-1./(i-1) ;
    for (int i=1 ; i<nr ; i+=2)
	som[i] = 0 ;
    
    double res = 0 ;
    for (int i=0 ; i<nr ; i++)
	res += som[i]*coef[i] ;
	
    res /= pente ;
    
    delete [] deg ;
    delete [] coef ;
    delete [] coloc ;
    delete [] som ;
    
    return res ;
}

double Bin_ns_bh::smarr() const {

    using namespace Unites ;

    // The tests
    Tenseur psiq_t (pow(star.get_confpsi()(), 4.)) ;
    psiq_t.set_std_base() ;
    
    Tenseur_sym furmet (star.mp, 2, CON, star.mp.get_bvect_cart()) ;
    furmet.set_etat_qcq() ;
    for (int i=0 ; i< 3 ; i++) {
        furmet.set(i,i) = 1/psiq_t() ;
	for(int j=i+1 ; j<3 ; j++)
	   furmet.set(i,j) = 0 ;
    }
    Metrique met (furmet, false) ;
    
    Tenseur_sym kij (star.get_tkij_tot()/psiq_t) ;
    kij.change_triad(star.mp.get_bvect_cart()) ;
    kij.dec2_dzpuis() ;
    Tenseur_sym kij_cov (manipule (kij, met)) ;
    Tenseur shift (star.get_shift()) ;
    shift.change_triad(star.mp.get_bvect_cart()) ;
    
    Tenseur aime (star.mp, 1, CON, star.mp.get_bvect_cart()) ;
    aime.set_etat_qcq() ;
    aime.set(0) = -omega*star.mp.ya ; 
    aime.set(1) = omega*star.mp.xa ;
    aime.set(2) = 0 ;
    aime.annule(star.mp.get_mg()->get_nzone()-1) ;
    aime.set_std_base() ;
    shift = shift - aime ;
 
    // La matière :
    Tenseur u_euler (star.get_u_euler()) ;
    u_euler.change_triad(star.mp.get_bvect_cart()) ;
    Tenseur u_i_bas (manipule(u_euler, met)) ;
    Tenseur mat (qpig*(star.get_nnn()*(star.get_ener_euler() + star.get_s_euler()) - 2*(star.get_ener_euler()+star.get_press())*contract(u_i_bas, 0, shift, 0))) ;
    
    // La partie avec la matière :
    Cmp psiq (pow(star.get_confpsi()(), 4.)) ;
    
    Cmp integ_matter (star.get_nnn()()*(star.get_ener_euler()() + star.get_s_euler()()) 
    		- 2*(star.get_ener_euler()()+star.get_press()())*psiq*flat_scalar_prod(u_euler, shift)()) ;
    integ_matter = integ_matter * pow(star.get_confpsi()(),6.) ;
    integ_matter.std_base_scal() ;
    integ_matter.annule(star.get_nzet(), star.get_mp().get_mg()->get_nzone()-1) ;
    double matter_term = integ_matter.integrale()*qpig/4/M_PI ;
    
    // Integrale sur horizon :
    Cmp integ_hor (hole.get_n_tot()().dsdr()*
	pow(hole.get_psi_tot()(), 2)) ;
    integ_hor.std_base_scal() ;
    integ_hor.raccord(1) ;
    double hor_term = hole.mp.integrale_surface(integ_hor, hole.get_rayon()) ;	
    hor_term /= 4*M_PI ;
   
    double m_test = hor_term + matter_term + 2*omega*moment_systeme_inf() + 
      2*(hole.omega_local-omega)*hole.local_momentum() ;
    
    return m_test ;
    }
