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
 * Revision 1.5  2018/11/16 14:34:36  j_novak
 * Changed minor points to avoid some compilation warnings.
 *
 * Revision 1.4  2016/12/05 16:18:10  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:29  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:16:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/05/15  15:46:43  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/04/27  15:19:52  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Header C : 
#include <cstdlib>
#include <cmath>

// Headers Lorene :
#include "matrice.h"
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
#include "base_val.h"
#include "proto.h"
#include "type_parite.h"
#include "utilitaires.h"




	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------

namespace Lorene {
Mtbl_cf sol_poisson_frontiere_double (const Map_af& mapping, 
    const Mtbl_cf& source, const Mtbl_cf& lim_func, const Mtbl_cf& lim_der, 
					int num_zone)

{
    
    // Verifications d'usage sur les zones
    int nz = source.get_mg()->get_nzone() ;
    assert (nz>1) ;
    assert ((num_zone>0) && (num_zone<nz-1)) ;
    assert(source.get_mg()->get_type_r(num_zone) == FIN) ;
    
    assert (lim_func.get_mg() == source.get_mg()->get_angu()) ;
    assert (lim_der.get_mg() == source.get_mg()->get_angu()) ;
    assert (source.get_etat() != ETATNONDEF) ;
    assert (lim_func.get_etat() != ETATNONDEF) ;
    assert (lim_der.get_etat() != ETATNONDEF) ;
     
    // Bases spectrales
    const Base_val& base = source.base ;
    
    // donnees sur la zone
    int nr = source.get_mg()->get_nr(num_zone) ;
    int nt = source.get_mg()->get_nt(num_zone) ;
    int np = source.get_mg()->get_np(num_zone) ;;
    int base_r ;
    int l_quant, m_quant;
    
    double alpha = mapping.get_alpha()[num_zone] ;
    double beta = mapping.get_beta()[num_zone] ;
    double echelle = beta/alpha ;
    double facteur ;
    
    //Rangement des valeurs intermediaires 
    Tbl *so ;
    Tbl *sol_hom ;
    Tbl *sol_part ;
    Matrice *operateur ;
    Matrice *nondege ;
    
    
    Mtbl_cf resultat(source.get_mg(), base) ;
    resultat.annule_hard() ;
    
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1)
	    {
		// calcul des nombres quantiques :
		donne_lm(nz, num_zone, j, k, base, m_quant, l_quant, base_r) ;
	    
		// Construction de l'operateur
		operateur = new Matrice(laplacien_mat
				    (nr, l_quant, echelle, 0, base_r)) ;
		
		(*operateur) = combinaison(*operateur, l_quant, echelle, 0, 
									 base_r) ;
		
		 // Operateur inversible
		nondege = new Matrice(prepa_nondege(*operateur, l_quant, 
							echelle, 0, base_r)) ;		
		
		// Calcul DES DEUX SH
		sol_hom = new Tbl(solh(nr, l_quant, echelle, base_r)) ;
		
		// Calcul de la SP
		so = new Tbl(nr) ;
		so->set_etat_qcq() ;
		for (int i=0 ; i<nr ; i++)
		    so->set(i) = source(num_zone, k, j, i) ;
		
		sol_part = new Tbl (solp(*operateur, *nondege, alpha,
					 beta, *so, 0, base_r)) ;
		
		 //-------------------------------------------
		// On est parti pour imposer la boundary
		//-------------------------------------------    
		// Conditions de raccord type Dirichlet :
		// Pour la sp :
		double somme = 0 ;
		for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0)
		    somme += (*sol_part)(i) ;
		  else
		    somme -= (*sol_part)(i) ;
		
		facteur = (lim_func(num_zone-1, k, j, 0)-somme)
		  * pow(echelle-1, l_quant+1) ;
		
		for (int i=0 ; i<nr ; i++)
		  sol_part->set(i) +=
		    facteur*(*sol_hom)(1, i) ;
		
		// pour l'autre solution homogene :
		facteur = - pow(echelle-1, 2*l_quant+1) ;
		for (int i=0 ; i<nr ; i++)
		  sol_hom->set(0, i) +=
		    facteur*(*sol_hom)(1, i) ;
		
		// Condition de raccord de type Neumann :
		double val_der_solp = 0 ;
		for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0)
		    val_der_solp -= i*i*(*sol_part)(i)/alpha ;
		  else
		    val_der_solp += i*i*(*sol_part)(i)/alpha ;
		
		double val_der_solh = 0 ;
		for (int i=0 ; i<nr ; i++)
		  if (i%2 == 0)
		    val_der_solh -= i*i*(*sol_hom)(0, i)/alpha ;
		  else
		    val_der_solh += i*i*(*sol_hom)(0, i)/alpha ;
		
		assert (val_der_solh != 0) ;
		
		facteur = (lim_der(num_zone-1, k, j, 0)-val_der_solp) /
		  val_der_solh ;
		
		for (int i=0 ; i<nr ; i++)
		  sol_part->set(i) +=
		    facteur*(*sol_hom)(0, i) ;
		
		// solp contient le bon truc (normalement ...)
		for (int i=0 ; i<nr ; i++)
		  resultat.set(num_zone, k, j, i) = (*sol_part)(i) ;
		
		delete operateur ;
		delete nondege ;
		delete so ;
		delete sol_hom ;
		delete sol_part ;
	    }
    return resultat ;
}
}
