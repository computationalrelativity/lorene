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


char regularise_shift_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2006/04/27 09:12:32  p_grandclement
 * First try at irrotational black holes
 *
 * Revision 1.4  2005/08/29 15:10:14  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.3  2003/10/03 15:58:44  j_novak
 * Cleaning of some headers
 *
 * Revision 1.2  2003/02/13 16:40:25  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2001/05/07  09:11:48  phil
 * *** empty log message ***
 *
 * Revision 2.5  2001/01/29  14:30:48  phil
 * ajout type rotation
 *
 * Revision 2.4  2000/11/02  10:18:05  phil
 * modification du degre de l;a fonction\
 * de correction. On passe de 2 a 3
 *
 * Revision 2.3  2000/10/31  10:04:28  phil
 * correction importation
 *
 * Revision 2.2  2000/10/26  12:34:17  phil
 * *** empty log message ***
 *
 * Revision 2.1  2000/10/26  12:31:55  phil
 * correction orientation pour calcul du shift total
 *
 * Revision 2.0  2000/10/19  10:08:03  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


//Standard
#include <stdlib.h>
#include <math.h>

//Lorene
#include "nbr_spx.h"
#include "tenseur.h"

double regle (Tenseur& shift_auto, const Tenseur& shift_comp, double omega, double omega_local) {
    
    Tenseur shift_old (shift_auto) ;
    
    double orientation = shift_auto.get_mp()->get_rot_phi() ;
    assert ((orientation==0) || (orientation == M_PI)) ;
    double orientation_autre = shift_comp.get_mp()->get_rot_phi() ;
    assert ((orientation_autre==0) || (orientation_autre == M_PI)) ;
    
    int alignes = (orientation == orientation_autre) ? 1 : -1 ;
    
    // Cas triades identiques
    if (*shift_comp.get_triad() == *shift_auto.get_triad())
         alignes = 1 ;
    
    int nz = shift_auto.get_mp()->get_mg()->get_nzone() ;
    int np = shift_auto.get_mp()->get_mg()->get_np(1) ;
    int nt = shift_auto.get_mp()->get_mg()->get_nt(1) ;
    int nr = shift_auto.get_mp()->get_mg()->get_nr(1) ;
    
    // On minimise la valeur de la derivee de B sur R :
    Tenseur shift_tot (*shift_auto.get_mp(), 1, CON, *shift_auto.get_triad()) ;
    shift_tot.set_etat_qcq() ;
    shift_tot.set(0).import_asymy (alignes*shift_comp(0)) ;
    shift_tot.set(1).import_symy (alignes*shift_comp(1)) ;
    shift_tot.set(2).import_asymy (shift_comp(2)) ;

    shift_tot = shift_tot + shift_auto ;
 
    double indic = (orientation == 0) ? 1 : -1 ;
    
    Mtbl xa_mtbl (shift_tot.get_mp()->get_mg()) ;
    xa_mtbl = shift_tot.get_mp()->xa ;
    Mtbl ya_mtbl (shift_tot.get_mp()->get_mg()) ;
    ya_mtbl = shift_tot.get_mp()->ya ;
    
    Tenseur tbi (shift_tot) ;
    if (omega != 0) {
	for (int i=0 ; i<3 ; i++) {
	    tbi.set(i).va.coef_i() ;
	    tbi.set(i).va.set_etat_c_qcq() ;
	    }
	    
	tbi.set(0) = *shift_tot(0).va.c - indic *omega * ya_mtbl(0,0,0,0) - indic*omega_local* shift_tot.get_mp()->y ;
	tbi.set(1) = *shift_tot(1).va.c + indic *omega * xa_mtbl(0,0,0,0) + indic*omega_local* shift_tot.get_mp()->x ;
	tbi.set_std_base() ;
	tbi.set(0).annule(nz-1) ;
	tbi.set(1).annule(nz-1) ;
    }
     
    Tenseur derive_r (*shift_auto.get_mp(), 1, CON, *shift_auto.get_triad()) ;
    derive_r.set_etat_qcq() ;
    for (int i=0 ; i<3 ; i++) {
	if (tbi(i).get_etat() != ETATZERO)
	    derive_r.set(i) = tbi(i).dsdr() ;
	else
	    derive_r.set(i).set_etat_zero() ;
	}
    
    // On enleve un fonction pour rendre Kij regulier !
    
    Valeur val_hor (shift_auto.get_mp()->get_mg()) ;
    Valeur fonction_radiale (shift_auto.get_mp()->get_mg()) ;
    Cmp enleve (shift_auto.get_mp()) ;
  
    double erreur = 0 ;
    for (int comp=0 ; comp<3 ; comp++)
	{
	    val_hor.annule_hard() ; // Pour initialiser les tableaux 
	    for (int k=0 ; k<np ; k++)
		for (int j=0 ; j<nt ; j++)
		    for (int i=0 ; i<nr ; i++)
		    val_hor.set(1, k, j, i) = derive_r(comp) (1, k, j, 0) ;
			     
	    double r_0 = shift_auto.get_mp()->val_r (1, -1, 0, 0) ;
	    double r_1 = shift_auto.get_mp()->val_r (1, 1, 0, 0) ;
	    
	    fonction_radiale = pow(r_1-shift_auto.get_mp()->r, 3.)*
		    (shift_auto.get_mp()->r-r_0)/pow(r_1-r_0, 3.) ;
	    fonction_radiale.annule(0) ;
	    fonction_radiale.annule(2, nz-1) ;
	      
	    enleve = fonction_radiale * val_hor ;
	    enleve.va.base = shift_auto(comp).va.base ;
	    
	    if (norme(enleve)(1) != 0)
		shift_auto.set(comp) = shift_auto(comp) - enleve ;
	    if (norme(shift_auto(comp))(1) > 1e-5) {
		Tbl diff (diffrelmax (shift_auto(comp), shift_old(comp))) ;
		if (erreur < diff(1))
		    erreur = diff(1) ;
	    }
    }
    return erreur ;
}
