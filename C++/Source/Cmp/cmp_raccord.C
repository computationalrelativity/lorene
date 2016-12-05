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


 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:17:49  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/10/03 15:58:45  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/09/07  13:19:58  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/06/06  12:18:27  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>
#include <cmath>

// LORENE
#include "matrice.h"
#include "cmp.h"
#include "proto.h"


namespace Lorene {
Matrice matrice_raccord_pair (int cont, double alpha_kernel) {
    
    Matrice systeme (cont, cont) ;
    systeme.set_etat_qcq() ;
    for (int i=0 ; i<cont ; i++)
	for (int j=0 ; j<cont ; j++)
	    systeme.set(i, j) = 0 ;
    
    double somme ;
    for (int i=0 ; i<cont ; i++)
	for (int k=0 ; k<cont ; k++)
	    if (k<= 2*i) {
	    somme = 1 ;
	    for (int boucle=0 ; boucle<k ; boucle++)
		somme *= (4*i*i-boucle*boucle)/(2.*boucle+1.)/alpha_kernel ;
	    systeme.set(k, i) = somme ;
	}
    int inf = (cont%2 == 1) ? (cont-1)/2 : (cont-2)/2 ;
    systeme.set_band (cont-1, inf) ;
    systeme.set_lu() ;
    return systeme ;
}

Matrice matrice_raccord_impair (int cont, double alpha_kernel) {
    
    Matrice systeme (cont, cont) ;
    systeme.set_etat_qcq() ;
    for (int i=0 ; i<cont ; i++)
	for (int j=0 ; j<cont ; j++)
	    systeme.set(i, j) = 0 ;
    
    double somme ;
    for (int i=0 ; i<cont ; i++)
	for (int k=0 ; k<cont ; k++)
	    if (k<= 2*i+1) {
	    somme = 1 ;
	    for (int boucle=0 ; boucle<k ; boucle++)
		somme *= (pow(2*i+1, 2.)-boucle*boucle)/(2.*boucle+1.)/alpha_kernel ;
	    systeme.set(k, i) = somme ;
	}
    int inf = (cont%2 == 0) ? cont/2 : (cont-1)/2 ;
    systeme.set_band (cont-1, inf) ;
    systeme.set_lu() ;
    return systeme ;
}


Tbl sec_membre_raccord (Tbl coef, int cont, double alpha_shell) {
    
    assert (coef.get_etat() != ETATNONDEF) ;
    int nr = coef.get_dim(0) ;
    
    Tbl sec_membre(cont) ;
    sec_membre.set_etat_qcq() ;
    for (int i=0 ; i<cont ; i++)
	sec_membre.set(i) = 0 ;
    
    double somme ;
    for (int i=0 ; i<nr ; i++)
	for (int k=0 ; k<cont ; k++)
	    if (k<= i) {
	    somme = 1 ;
	    for (int boucle=0 ; boucle<k ; boucle++)
		somme *= (i*i-boucle*boucle)/(2.*boucle+1.)/alpha_shell ;
	    if ((i+k)%2 == 0)
		sec_membre.set(k) += coef(i)*somme ;
	    else
		sec_membre.set(k) -= coef(i)*somme ;
	}
    
    return sec_membre ;
}


Tbl regularise (Tbl coef, int nr, int base_r) {
    
    assert ((base_r == R_CHEBI) || (base_r == R_CHEBP)) ;
    assert (coef.get_etat() != ETATNONDEF) ;
    int cont = coef.get_dim(0) ;
    assert (nr >= cont) ;
    
    Tbl resu (nr) ;
    resu.set_etat_qcq() ;
    
    double* x4coef = new double[nr] ;
    for (int i=0 ; i<cont ; i++)
	x4coef[i] = coef(i) ;
    for (int i=cont ; i<nr ; i++)
	x4coef[i] = 0 ;
    double* x6coef = new double[nr] ;
    
    multx2_1d (nr, &x4coef, base_r) ;
    multx2_1d (nr, &x4coef, base_r) ;
    for (int i=0 ; i<nr ; i++)
	x6coef[i] = x4coef[i] ;
    multx2_1d (nr, &x6coef, base_r) ;
    
    for (int i=0 ; i<nr ; i++)
	resu.set(i) = 3*x4coef[i]-2*x6coef[i] ;
    
    delete [] x4coef ;
    delete [] x6coef ;
    
    return resu ;
}



void Cmp::raccord (int aux) {
    assert (etat != ETATNONDEF) ;
    
    assert (aux >=0) ;
    int cont = aux+1 ;
    
    const Map_af* mapping = dynamic_cast<const Map_af*>(get_mp() ) ; 

    if (mapping == 0x0) {
	cout << 
	"raccord : The mapping does not belong to the class Map_af !"
	    << endl ; 
	abort() ;
    }
    
    assert (mapping->get_mg()->get_type_r(1) == FIN) ;
    assert (mapping->get_mg()->get_type_r(0) == RARE) ;
    
    // On passe en Ylm et vire tout dans la zone interne...
    va.coef() ;
    va.ylm() ;
    va.set_etat_cf_qcq() ;
    va.c_cf->t[0]->annule_hard() ;
    
    // Confort :
    int nz = mapping->get_mg()->get_nzone() ;
    int nbrer_kernel = mapping->get_mg()->get_nr(0) ;
    int nbrer_shell  = mapping->get_mg()->get_nr(1) ;
    
    int nbret_kernel = mapping->get_mg()->get_nt(0) ;
    int nbret_shell  = mapping->get_mg()->get_nt(1) ;
    
    int nbrep_kernel = mapping->get_mg()->get_np(0) ;
    int nbrep_shell  = mapping->get_mg()->get_np(1) ;
    
    double alpha_kernel = mapping->get_alpha()[0] ;
    double alpha_shell  = mapping->get_alpha()[1] ;
    
    int base_r, m_quant, l_quant ;
    
    for (int k=0 ; k<nbrep_kernel+1 ; k++)
	for (int j=0 ; j<nbret_kernel ; j++)
	    if (nullite_plm(j, nbret_kernel, k,nbrep_kernel, va.base) == 1)
		 if (nullite_plm(j, nbret_shell, k, nbrep_shell, va.base) == 1)
	{
		// calcul des nombres quantiques :
	    donne_lm(nz, 0, j, k, va.base, m_quant, l_quant, base_r) ;
	    assert ((base_r == R_CHEBP) || (base_r == R_CHEBI)) ;
	    
	    Matrice systeme(cont, cont) ;
	    
	    Tbl facteur (nbrer_kernel) ;
	    facteur.annule_hard() ;
	    for (int i=0 ; i<nbrer_shell ; i++)
		if (i<nbrer_kernel)
		    facteur.set(i) = (*va.c_cf)(1, k, j, i) ;
	    
	    Tbl sec_membre (sec_membre_raccord (facteur, cont, alpha_shell)) ;
	   
	    if (base_r == R_CHEBP)
		systeme = matrice_raccord_pair (cont, alpha_kernel) ;	    
	    else
		systeme = matrice_raccord_impair (cont, alpha_kernel) ;
	    
	    Tbl soluce (systeme.inverse(sec_membre)) ;
	    Tbl regulier (nbrer_kernel) ;
	    
	    if (l_quant == 0)
		for (int i=0 ; i<cont ; i++)
		    va.c_cf->set(0, k, j, i) = soluce(i) ;
	    else {
		if (l_quant %2 == 0)
		    regulier = regularise (soluce, nbrer_kernel, R_CHEBP) ;
		else
		    regulier = regularise (soluce, nbrer_kernel, R_CHEBI) ;
		
		for (int i=0 ; i<nbrer_kernel ; i++)
		    va.c_cf->set(0, k, j, i) = regulier(i) ;
		}
	    }
    va.ylm_i() ;
}
}
