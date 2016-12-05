/*
 *  Method of Non_class_members for solving Poisson equations
 *   with a multipole falloff condition at the outer boundary
 *
 */

/*
 *   Copyright (c) 2004 Joshua A. Faber
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:10  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:30  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:16:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/29 16:32:02  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cmath>

// Lorene headers
#include "matrice.h"
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
#include "base_val.h"
#include "proto.h"
#include "type_parite.h"



	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------

/*
 * 
 * Solution de l'equation de poisson
 * 
 * Entree : mapping :   le mapping affine
 *	    source : les coefficients de la source qui a ete multipliee par
 *		    r^4 ou r^2 dans la ZEC.
 *		    La base de decomposition doit etre Ylm
 *	    k_ylm: exponent in radial dependence of field: phi \propto r^{-k}
 * Sortie : renvoie les coefficients de la solution dans la meme base de 
 *	    decomposition (a savoir Ylm)
 *	    
 */



namespace Lorene {
Mtbl_cf sol_poisson_ylm(const Map_af& mapping, const Mtbl_cf& source, const int nylm, const double* intvec)
{
 
    // Verifications d'usage sur les zones
    int nz = source.get_mg()->get_nzone() ;
    assert (nz>1) ;
    assert (source.get_mg()->get_type_r(0) == RARE) ;
    //    assert (source.get_mg()->get_type_r(nz-1) == UNSURR) ;
    for (int l=1 ; l<nz ; l++)
	assert(source.get_mg()->get_type_r(l) == FIN) ;

    //     assert ((dzpuis==4) || (dzpuis==2) || (dzpuis==3)) ;
    
    
    // Bases spectrales
    const Base_val& base = source.base ;
    
    // donnees sur la zone
    int nr, nt, np ;
    int base_r ;
    double alpha, beta, echelle ;
    int l_quant, m_quant;
    
    //Rangement des valeurs intermediaires 
    Tbl *so ;
    Tbl *sol_hom ;
    Tbl *sol_part ;
    Matrice *operateur ;
    Matrice *nondege ;
    
    
    // Rangement des solutions, avant raccordement
    Mtbl_cf solution_part(source.get_mg(), base) ;
    Mtbl_cf solution_hom_un(source.get_mg(), base) ;
    Mtbl_cf solution_hom_deux(source.get_mg(), base) ;
    Mtbl_cf resultat(source.get_mg(), base) ;
    
    solution_part.set_etat_qcq() ;
    solution_hom_un.set_etat_qcq() ;
    solution_hom_deux.set_etat_qcq() ;
    resultat.set_etat_qcq() ;
    
    for (int l=0 ; l<nz ; l++) {
	solution_part.t[l]->set_etat_qcq() ;
	solution_hom_un.t[l]->set_etat_qcq() ;
	solution_hom_deux.t[l]->set_etat_qcq() ;
	resultat.t[l]->set_etat_qcq() ;
	for (int k=0 ; k<source.get_mg()->get_np(l)+2 ; k++)
	    for (int j=0 ; j<source.get_mg()->get_nt(l) ; j++)
		for (int i=0 ; i<source.get_mg()->get_nr(l) ; i++)
		    resultat.set(l, k, j, i) = 0 ;
    }
    
    // nbre maximum de point en theta et en phi :
    int np_max, nt_max ;
    
		   //---------------
		  //--  NOYAU -----
		 //---------------
		 
    nr = source.get_mg()->get_nr(0) ;
    nt = source.get_mg()->get_nt(0) ;
    np = source.get_mg()->get_np(0) ;
    
    nt_max = nt ;
    np_max = np ;
    
    alpha = mapping.get_alpha()[0] ;
    beta = mapping.get_beta()[0] ;
    
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1) 
	    {
	    // calcul des nombres quantiques :
	    donne_lm(nz, 0, j, k, base, m_quant, l_quant, base_r) ;
	    
	    // Construction operateur
	    operateur = new Matrice(laplacien_mat(nr, l_quant, 0., 0, base_r)) ;
	    (*operateur) = combinaison(*operateur, l_quant, 0., 0, base_r) ;
	    
	    //Operateur inversible
	    nondege = new Matrice(prepa_nondege(*operateur, l_quant, 0., 0, base_r)) ;
	    
	    // Calcul de la SH
	    sol_hom = new Tbl(solh(nr, l_quant, 0., base_r)) ;
	 
	    //Calcul de la SP
	    so = new Tbl(nr) ;
	    so->set_etat_qcq() ;
	    for (int i=0 ; i<nr ; i++)
	      so->set(i) = source(0, k, j, i) ;
		
	    sol_part = new Tbl(solp(*operateur, *nondege, alpha, beta, 
				    *so, 0, base_r)) ;
	    
	    // Rangement dans les tableaux globaux ;
	    
	    for (int i=0 ; i<nr ; i++) {
		solution_part.set(0, k, j, i) = (*sol_part)(i) ;
		solution_hom_un.set(0, k, j, i) = (*sol_hom)(i) ;
		solution_hom_deux.set(0, k, j, i) = 0. ; 
	    }	    
	    
	    
	    
	    delete operateur ;
	    delete nondege ;
	    delete so ;
	    delete sol_hom ;
	    delete sol_part ;
	}


    
		   //---------------
		  //- COQUILLES ---
		 //---------------

    for (int zone=1 ; zone<nz ; zone++) {
	nr = source.get_mg()->get_nr(zone) ;
	nt = source.get_mg()->get_nt(zone) ;
	np = source.get_mg()->get_np(zone) ;
	
	if (np > np_max) np_max = np ;
	if (nt > nt_max) nt_max = nt ;
	
	alpha = mapping.get_alpha()[zone] ;
	beta = mapping.get_beta()[zone] ;

	for (int k=0 ; k<np+1 ; k++)
	    for (int j=0 ; j<nt ; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1)
	    {
		// calcul des nombres quantiques :
		donne_lm(nz, zone, j, k, base, m_quant, l_quant, base_r) ;
	    
		// Construction de l'operateur
		operateur = new Matrice(laplacien_mat
				    (nr, l_quant, beta/alpha, 0, base_r)) ;
		
		(*operateur) = combinaison(*operateur, l_quant, beta/alpha, 0, 
									 base_r) ;
		
		 // Operateur inversible
		nondege = new Matrice(prepa_nondege(*operateur, l_quant, 
							beta/alpha, 0, base_r)) ;		
		
		// Calcul DES DEUX SH
		sol_hom = new Tbl(solh(nr, l_quant, beta/alpha, base_r)) ;
		
		// Calcul de la SP
		so = new Tbl(nr) ;
		so->set_etat_qcq() ;
		for (int i=0 ; i<nr ; i++)
		  so->set(i) = source(zone, k, j, i) ;
		
		sol_part = new Tbl (solp(*operateur, *nondege, alpha,
					 beta, *so, 0, base_r)) ;
		
		
		// Rangement
		for (int i=0 ; i<nr ; i++) {
		    solution_part.set(zone, k, j, i) = (*sol_part)(i) ;
		    solution_hom_un.set(zone, k, j, i) = (*sol_hom)(0, i) ;
		    solution_hom_deux.set(zone, k, j, i) = (*sol_hom)(1, i) ;
		}
			
		
		delete operateur ;
		delete nondege ;
		delete so ;
		delete sol_hom ;
		delete sol_part ;
	    }
    }

	     //-------------------------------------------
	    // On est parti pour le raccord des solutions
	   //-------------------------------------------
    // Tableau de 0 et 1 sur les zones, indiquant si le Plm considere
    // intervient dans le developpement de cette zone.
    int * indic = new int [nz] ;
    int taille ;
    Tbl *sec_membre ; // termes constants du systeme
    Matrice *systeme ; //le systeme a resoudre

    // des compteurs pour le remplissage du systeme
    int ligne ;
    int colonne ;

    // compteurs pour les diagonales du systeme :
    int sup ;
    int inf ;
    int sup_new, inf_new ;

    // on boucle sur le plus grand nombre possible de Plm intervenant...
    for (int k=0 ; k<np_max+1 ; k++)
	for (int j=0 ; j<nt_max ; j++)
	    if (nullite_plm(j, nt_max, k, np_max, base) == 1) {
		
		ligne = 0 ;
		colonne = 0 ;
		sup = 0 ;
		inf = 0 ;
		
		for (int zone=0 ; zone<nz ; zone++)
			indic[zone] = nullite_plm(j, source.get_mg()->get_nt(zone), 
					 k,  source.get_mg()->get_np(zone), base);
		
		// taille du systeme a resoudre pour ce Plm 
		taille = indic[0] ;
		for (int zone=1 ; zone<nz ; zone++)
		    taille+=2*indic[zone] ;
		
		// on verifie que la taille est non-nulle.
		// cas pouvant a priori se produire...
		
		if (taille != 0) {
		    
		    sec_membre = new Tbl(taille) ;
		    sec_membre->set_etat_qcq() ;
		    for (int l=0 ; l<taille ; l++)
			sec_membre->set(l) = 0 ;
			
		    systeme = new Matrice(taille, taille) ;
		    systeme->set_etat_qcq() ;
		    for (int l=0 ; l<taille ; l++)
			for (int c=0 ; c<taille ; c++)
			    systeme->set(l, c) = 0 ;
		    
		    //Calcul des nombres quantiques
		    //base_r est donne dans le noyau, sa valeur dans les autres
		    //zones etant inutile.
		    
		    donne_lm (nz, 0, j, k, base, m_quant, l_quant, base_r) ;
			
			
		     //--------------------------
		    //		NOYAU 
		   //---------------------------
		    
		    if (indic[0] == 1) {
			nr = source.get_mg()->get_nr(0) ;
			
			alpha = mapping.get_alpha()[0] ;
			// valeur de x^l en 1 :
			systeme->set(ligne, colonne) = 1. ; /* ligne=0, colonne=0 */
			// coefficient du Plm dans la solp
			for (int i=0 ; i<nr ; i++)
			    sec_membre->set(ligne) -= solution_part(0, k, j, i) ; /* ligne=0 */
			    
			ligne++ ;  /* ligne=1, colonne=0 */
			// on prend les derivees que si Plm existe 
			//dans la zone suivante
			
			if (indic[1] == 1) {
			    // derivee de x^l en 1 :
			    systeme->set(ligne, colonne) = 1./alpha*l_quant ; /* ligne=1, colonne=0 */
			
			    // coefficient de la derivee du Plm dans la solp
			    if (base_r == R_CHEBP)
				// cas en R_CHEBP 
				for (int i=0 ; i<nr ; i++)
				    sec_membre->set(ligne) -= 
					    4*i*i/alpha
						*solution_part(0, k, j, i) ; /* ligne=1 */
			    else
				// cas en R_CHEBI
				for (int i=0 ; i<nr ; i++)
				    sec_membre->set(ligne) -=
					(2*i+1)*(2*i+1)/alpha
					    *solution_part(0, k, j, i) ; /* ligne=1 */
			    		    
			    // on a au moins un diag inferieure dans ce cas ...
			    inf = 1 ;
			}
		    colonne++ ; /* ligne=1, colonne=1 */
		    }
		    
		    //-----------------------------
		   //        COQUILLES
		  //------------------------------
		  
		  for (int zone=1 ; zone<nz ; zone++)
		    if (indic[zone] == 1) {
			
			nr = source.get_mg()->get_nr(zone) ;
			alpha = mapping.get_alpha()[zone] ;
			echelle = mapping.get_beta()[zone]/alpha ;
			
			//Frontiere avec la zone precedente :
			if (indic[zone-1] == 1) ligne -- ; /* ligne=0, colonne=1 */
			
			//valeur de (x+echelle)^l en -1 :
			systeme->set(ligne, colonne) = 
			    -pow(echelle-1, double(l_quant)) ; /* ligne=0, colonne=1 */
			
			// valeur de 1/(x+echelle) ^(l+1) en -1 :
			systeme->set(ligne, colonne+1) = 
			    -1/pow(echelle-1, double(l_quant+1)) ; /* ligne=0, colonne=1, colonne+1=2 */
			
			// la solution particuliere :
			for (int i=0 ; i<nr ; i++)
			    if (i%2 == 0)
				sec_membre->set(ligne) += solution_part(zone, k, j, i) ;
			    else sec_membre->set(ligne) -= solution_part(zone, k, j, i) ; /* ligne=0 */
			
			// les diagonales :
			sup_new = colonne+1-ligne ; /* ligne=0, colonne=1, colonne+1-ligne=2 */
			if (sup_new > sup) sup = sup_new ;
			
			ligne++ ; /* ligne=1 */
			
			// on prend les derivees si Plm existe dans la zone 
			// precedente :
			
			if (indic[zone-1] == 1) {
			// derivee de (x+echelle)^l en -1 :
			    systeme->set(ligne, colonne) = 
			      -l_quant*pow(echelle-1, double(l_quant-1))/alpha ; /* ligne=1, colonne=1 */
			// derivee de 1/(x+echelle)^(l+1) en -1 :
			    systeme->set(ligne, colonne+1) = 
				(l_quant+1)/pow(echelle-1, double(l_quant+2))/alpha ; /* ligne=1, colonne=1, colonne+1=2 */
			
			
			
			// la solution particuliere :
			    for (int i=0 ; i<nr ; i++)
				if (i%2 == 0)
				    sec_membre->set(ligne) -= 
				      i*i/alpha*solution_part(zone, k, j, i) ; /* ligne=1 */
				else
				    sec_membre->set(ligne) +=
					i*i/alpha*solution_part(zone, k, j, i) ; /* ligne=1 */
			
			// les diagonales :
			sup_new = colonne+1-ligne ; /* ligne=1, colonne=1, colonne+1-ligne=1 */
			if (sup_new > sup) sup = sup_new ;
			
			ligne++ ; /* ligne=2 */
			}
			

			if(zone < nz-1) {

			// Frontiere avec la zone suivante :
			//valeur de (x+echelle)^l en 1 :
			systeme->set(ligne, colonne) = 
			  pow(echelle+1, double(l_quant)) ; /* ligne=2, colonne=1 */ 
			
			// valeur de 1/(x+echelle)^(l+1) en 1 :
			systeme->set(ligne, colonne+1) =
			    1/pow(echelle+1, double(l_quant+1)) ; /* ligne=2, colonne=1, colonne+1=2 */
			    
			// la solution particuliere :
			for (int i=0 ; i<nr ; i++)
			  sec_membre->set(ligne) -= solution_part(zone, k, j, i) ; /* ligne=2 */
			
			// les diagonales inf :
			inf_new = ligne-colonne ; /* ligne=2, colonne=1, ligne-colonne=1 */
			if (inf_new > inf) inf = inf_new ;  
			
			ligne ++ ; /* ligne=3 */
			    
			// Utilisation des derivees ssi Plm existe dans la
			//zone suivante :
			if (indic[zone+1] == 1) {
			    
			    //derivee de (x+echelle)^l en 1 :
			    systeme->set(ligne, colonne) = 
			      l_quant*pow(echelle+1, double(l_quant-1))/alpha ; /* ligne=3, colonne=1 */
				
			    //derivee de 1/(echelle+x)^(l+1) en 1 :
			    systeme->set(ligne, colonne+1) =
				-(l_quant+1)/pow(echelle+1, double(l_quant+2))/alpha ; /* ligne=3, colonne=1, colonne+1=2 */
			    
			    // La solution particuliere :
			    for (int i=0 ; i<nr ; i++)
				sec_membre->set(ligne) -= 
				  i*i/alpha*solution_part(zone, k, j, i) ; /* ligne=3 */
				    
			// les diagonales inf :
			inf_new = ligne-colonne ; /* ligne=3, colonne=1, ligne-colonne=2 */
			if (inf_new > inf) inf = inf_new ;  
			    
			     }
			colonne += 2 ; /* ligne=colonne=3 -> ligne=colonne=1 next block of two */
			} else {



   		    //-------------------------
		   //  Ylm outer boundary		    
		  //---------------------------

			  /* ligne=2, colonne=1 */
			  
			  // The coefficient for A_n is (k+l)(echelle+1)^l
			 systeme->set(ligne, colonne) = 
			   pow(echelle+1, double(l_quant)) ;  /* ligne=2, colonne=1 */

			 // The coefficient for B_n is (k-(l+1))(echelle+1)^{-(l+1)}
			systeme->set(ligne, colonne+1) =
			    pow(echelle+1, double(-1-l_quant)) ; /* ligne=2, colonne=1, colonne+1=2 */

			// Here we have -f - multipole integral 
			//(pos source->neg field!!)

			int indylm;
			double intterm=0.;
			if(l_quant*l_quant < nylm && m_quant<=l_quant) {
			  indylm=l_quant*l_quant+2*m_quant;
			  if(k%2==0 && k>=2)indylm-=1;
			  intterm=intvec[indylm];
			  cout <<"l,m:"<<l_quant<<" "<<m_quant<<" "<<j<<" "<<k<<" "<<indylm<<" "<<intterm<<endl;
			}
			// This is just for testing!!!
			//if(l_quant==1 && m_quant==1&& k%2==0) {
			//  intterm=1.0 ;
			//} else {
			//  intterm=0.;
			//}
			
			sec_membre->set(ligne) -= intterm;
			
			// La solution particuliere :
			for (int i=0 ; i<nr ; i++)
			  sec_membre->set(ligne) -= 
			    solution_part(zone, k, j, i) ; /* ligne=2 */

			// les diagonales inf :
			inf_new = ligne-colonne ; /* ligne=2, colonne=1, ligne-colonne=1 */
			if (inf_new > inf) inf = inf_new ;  
			    
			}
		    }
		    

		    //-------------------------
		   //  resolution du systeme 
		  //--------------------------
		    
		    systeme->set_band(sup, inf) ;
		    systeme->set_lu() ;

		    Tbl facteurs(systeme->inverse(*sec_membre)) ;
		    int conte = 0 ;
		    
		 
		// rangement dans le noyau :
			
		    if (indic[0] == 1) {
			nr = source.get_mg()->get_nr(0) ;
			for (int i=0 ; i<nr ; i++)
			    resultat.set(0, k, j, i) = solution_part(0, k, j, i)
				+facteurs(conte)*solution_hom_un(0, k, j, i) ;
			conte++ ;
			}
			
			// rangement dans les coquilles :
		    for (int zone=1 ; zone<nz ; zone++)
			if (indic[zone] == 1) {
			    nr = source.get_mg()->get_nr(zone) ;
			    for (int i=0 ; i<nr ; i++)
			    resultat.set(zone, k, j, i) = 
				solution_part(zone, k, j, i)
				+facteurs(conte)*solution_hom_un(zone, k, j, i) 
				+facteurs(conte+1)*solution_hom_deux(zone, k, j, i) ;
			    conte+=2 ;
			}
			
		    delete sec_membre ;
		    delete systeme ;
		}
		
	    }
    
    delete [] indic ;
    
    return resultat;
}
}
