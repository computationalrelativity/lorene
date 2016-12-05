/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 2000-2001 Philippe Grandclement (for preceding Cmp version)
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
 * Revision 1.5  2016/12/05 16:18:19  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:16:16  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/10/01 13:04:44  e_gourgoulhon
 * The method Tensor::get_mp() returns now a reference (and not
 * a pointer) onto a mapping.
 *
 * Revision 1.1  2003/09/25 08:58:10  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>

// LORENE
#include "tensor.h"
#include "proto.h"
#include "matrice.h"

namespace Lorene {
Matrice matrice_raccord_pair (int cont, double alpha_kernel) ;
Matrice matrice_raccord_impair (int cont, double alpha_kernel) ;
Tbl sec_membre_raccord (Tbl coef, int cont, double alpha_shell) ;
Tbl regularise (Tbl coef, int nr, int base_r) ;

void Scalar::raccord (int aux) {

    assert (etat != ETATNONDEF) ;
    
    assert (aux >=0) ;
    int cont = aux+1 ;
    
    const Map_af* mapping = dynamic_cast<const Map_af*>( mp ) ; 

    if (mapping == 0x0) {
	cout << 
	"Scalar::raccord : The mapping does not belong to the class Map_af !"
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
