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


char poisson_frontiere_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:29  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:16:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2004/11/23 12:50:44  f_limousin
 * Intoduce function poisson_dir_neu(...) to solve a scalar poisson
 * equation with a mixed boundary condition (Dirichlet + Neumann).
 *
 * Revision 1.2  2004/09/08 15:12:16  f_limousin
 * Delete some assert.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/05/22  16:03:32  phil
 * ajout du cas dzpuis = 3
 *
 * Revision 2.3  2000/05/15  15:46:35  phil
 * *** empty log message ***
 *
 * Revision 2.2  2000/04/27  12:28:56  phil
 * correction pour le raccord des differents domaines
 *
 * Revision 2.1  2000/03/20  13:06:32  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/03/17  17:24:49  phil
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
#include "valeur.h"



	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------

/*
 * 
 * Solution de l'equation de poisson avec Boundary condition a 
 * l'interieur d'une coquille.
 * 
 * Entree : mapping :   le mapping affine
 *	    source : les coefficients de la source qui a ete multipliee par
 *		    r^4 ou r^2 dans la ZEC.
 *		    La base de decomposition doit etre Ylm
 *	    limite : la CL (fonction angulaire) sur une frontiere spherique
 *	    type_raccord : 1 pour Dirichlet et 2 pour Neumann
 *	    num_front : indique la frontiere ou on impose la CL : 1 pour la 
 * frontiere situee entre le domain 1 et 2.
 *	    dzpuis : exposant de r dans le factor multiplicatif dans la ZEC
 * Sortie : renvoie les coefficients de la solution dans la meme base de 
 *	    decomposition (a savoir Ylm)
 *	    
 */



namespace Lorene {
Mtbl_cf sol_poisson_frontiere(const Map_af& mapping, const Mtbl_cf& source,
		const Mtbl_cf& limite, int type_raccord, int num_front, int dzpuis, double fact_dir, double fact_neu)

{
    
    // Verifications d'usage sur les zones
    int nz = source.get_mg()->get_nzone() ;
    assert (nz>1) ;
    assert ((num_front>=0) && (num_front<nz-2)) ;
    assert (source.get_mg()->get_type_r(nz-1) == UNSURR) ;
    for (int l=num_front+1 ; l<nz-1 ; l++)
	assert(source.get_mg()->get_type_r(l) == FIN) ;
    
    assert (source.get_etat() != ETATNONDEF) ;
    assert (limite.get_etat() != ETATNONDEF) ;
     
    assert ((dzpuis==4) || (dzpuis==2) || (dzpuis==3)) ;
    assert ((type_raccord == 1) || (type_raccord == 2)|| (type_raccord == 3)) ;
    
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
    int np_max = 0 ;
    int nt_max = 0 ;
    

		   //---------------
		  //--  ZEC   -----
		 //---------------
		 
    nr = source.get_mg()->get_nr(nz-1) ;
    nt = source.get_mg()->get_nt(nz-1) ;
    np = source.get_mg()->get_np(nz-1) ;
    
    if (np > np_max) np_max = np ;
    if (nt > nt_max) nt_max = nt ;
   
    alpha = mapping.get_alpha()[nz-1] ;
    beta = mapping.get_beta()[nz-1] ;
    
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1)
	    {
	    
	    // calcul des nombres quantiques :
	    donne_lm(nz, nz-1, j, k, base, m_quant, l_quant, base_r) ;
	    
	    //Construction operateur
	    operateur = new Matrice(laplacien_mat(nr, l_quant, 0., dzpuis,
								 base_r)) ;
	    (*operateur) = combinaison(*operateur, l_quant, 0., dzpuis, base_r) ;
	    
	    // Operateur inversible
	    nondege = new Matrice(prepa_nondege(*operateur, l_quant, 0.,
							     dzpuis, base_r)) ;
	   
	    // Calcul de la SH
	    sol_hom = new Tbl(solh(nr, l_quant, 0., base_r)) ;
	   
	    // Calcul de la SP
	    so = new Tbl(nr) ;
	    so->set_etat_qcq() ;
	    for (int i=0 ; i<nr ; i++)
		so->set(i) = source(nz-1, k, j, i) ;
	    sol_part = new Tbl(solp(*operateur, *nondege, alpha, beta, 
				*so, dzpuis, base_r)) ;
	
	    // Rangement
	    for (int i=0 ; i<nr ; i++) {
		solution_part.set(nz-1, k, j, i) = (*sol_part)(i) ;
		solution_hom_un.set(nz-1, k, j, i) = (*sol_hom)(i) ;
		solution_hom_deux.set(nz-1, k, j, i) = 0. ;
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

    for (int zone=num_front+1 ; zone<nz-1 ; zone++) {
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
	    // On est parti pour imposer la boundary
	   //-------------------------------------------
	   
    nr = source.get_mg()->get_nr(num_front+1) ;
    nt = source.get_mg()->get_nt(num_front+1) ;
    np = source.get_mg()->get_np(num_front+1) ;
    double facteur ;
    double somme ;
    
    alpha = mapping.get_alpha()[num_front+1] ;
    beta = mapping.get_beta()[num_front+1] ;
    echelle = beta/alpha ;
    
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1)
	    {
		// calcul des nombres quantiques :
		donne_lm(nz, num_front+1, j, k, base, m_quant, l_quant, base_r) ;
		
		switch (type_raccord) {
		case 1 : 
		    // Conditions de raccord type Dirichlet :
		    // Pour la sp :
		    somme = 0 ;
		    for (int i=0 ; i<nr ; i++)
			if (i%2 == 0)
			    somme += solution_part(num_front+1, k, j, i) ;
			else
			    somme -= solution_part(num_front+1, k, j, i) ;
		    facteur = (limite(num_front, k, j, 0)-somme)
				* pow(echelle-1, l_quant+1) ;
		    
		    for (int i=0 ; i<nr ; i++)
			solution_part.set(num_front+1, k, j, i) +=
			    facteur*solution_hom_deux(num_front+1, k, j, i) ;
		    
		    // pour la solution homogene :
		    facteur = - pow(echelle-1, 2*l_quant+1) ;
		    for (int i=0 ; i<nr ; i++)
			solution_hom_un.set(num_front+1, k, j, i) +=
			    facteur*solution_hom_deux(num_front+1, k, j, i) ;
		    break ;
		 
		 
		 case 2 :
		    //Conditions de raccord de type Neumann :
		     // Pour la sp :
		    somme = 0 ;
		    for (int i=0 ; i<nr ; i++)
			if (i%2 == 0)
			    somme -= i*i/alpha*solution_part(num_front+1, k, j, i) ;
			else
			    somme += i*i/alpha*solution_part(num_front+1, k, j, i) ;
		    facteur = (somme-limite(num_front, k, j, 0))
				* alpha*pow(echelle-1, l_quant+2)/(l_quant+1) ;
		    for (int i=0 ; i<nr ; i++)
			solution_part.set(num_front+1, k, j, i) +=
			    facteur*solution_hom_deux(num_front+1, k, j, i) ;
		    
		    // pour la solution homogene :
		    facteur = pow(echelle-1, 2*l_quant+1)*l_quant/(l_quant+1) ;
		    for (int i=0 ; i<nr ; i++)
			solution_hom_un.set(num_front+1, k, j, i) +=
			    facteur*solution_hom_deux(num_front+1, k, j, i) ;
		    break ;
		    
		case 3 : 
		    // Conditions de raccord type Dirichlet-Neumann :
		    somme = 0 ;
		    for (int i=0 ; i<nr ; i++)
			if (i%2 == 0)
			    somme += solution_part(num_front+1, k, j, i) *
				fact_dir - fact_neu *
			        i*i/alpha*solution_part(num_front+1, k, j, i) ;
			else
			    somme += - solution_part(num_front+1, k, j, i) *
				fact_dir + fact_neu *
				i*i/alpha*solution_part(num_front+1, k, j, i) ;

		    double somme2 ;
		    somme2 = fact_dir * pow(echelle-1, -l_quant-1) - 
			fact_neu/alpha*pow(echelle-1, -l_quant-2)*(l_quant+1) ;

		    facteur = (limite(num_front, k, j, 0)- somme) / somme2 ;

		    for (int i=0 ; i<nr ; i++)
			solution_part.set(num_front+1, k, j, i) +=
			    facteur*solution_hom_deux(num_front+1, k, j, i) ;
		    
		    // pour la solution homogene :
		    double somme1 ;
		    somme1 = fact_dir * pow(echelle-1, l_quant) + 
			fact_neu / alpha * pow(echelle-1, l_quant-1) * 
			l_quant ;
		    facteur = - somme1 / somme2 ;
		    for (int i=0 ; i<nr ; i++)
			solution_hom_un.set(num_front+1, k, j, i) +=
			    facteur*solution_hom_deux(num_front+1, k, j, i) ;

		    break ;

		 default :
		    cout << "Diantres nous ne devrions pas etre ici ! " << endl ;
		    abort() ;
		    break ;
		}
		
		// Securite.
		for (int i=0 ; i<nr ; i++)
		    solution_hom_deux.set(num_front+1, k, j, i) = 0 ;
	    }
	    
	    
	     //-------------------------------------------
	    // On est parti pour le raccord des solutions
	   //-------------------------------------------
	   
    // On 
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
		
		for (int zone=num_front+1 ; zone<nz ; zone++)
			indic[zone] = nullite_plm(j, source.get_mg()->get_nt(zone), 
					 k,  source.get_mg()->get_np(zone), base);
		
		// taille du systeme a resoudre pour ce Plm 
		taille = indic[nz-1]+indic[num_front+1] ;
		for (int zone=num_front+2 ; zone<nz-1 ; zone++)
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
			
			
		     //----------------------------------
		    //		COQUILLE LA + INTERNE 
		   //------------------------------------
		    
		    if (indic[num_front+1] == 1) {
			nr = source.get_mg()->get_nr(num_front+1) ;
			
			alpha = mapping.get_alpha()[num_front+1] ;
			beta = mapping.get_beta()[num_front+1] ;
			echelle = beta/alpha ;
			
			// valeur de la solhomogene en 1 :
			somme = 0 ;
			for (int i=0 ; i<nr ; i++)
			    somme += solution_hom_un (num_front+1, k, j, i) ;
			systeme->set(ligne, colonne) = somme ;
			
			// coefficient du Plm dans la solp
			for (int i=0 ; i<nr ; i++)
			    sec_membre->set(ligne) -= solution_part(num_front+1, k, j, i) ;
			    
			ligne++ ;
			// on prend les derivees que si Plm existe 
			//dans la zone suivante
			
			if (indic[num_front+1] == 1) {

			    // derivee de la solhomogene en 1 :
			    somme = 0 ;
			    for (int i=0 ; i<nr ; i++)
				somme += i*i/alpha*
				    solution_hom_un(num_front+1, k, j, i) ;
			    systeme->set(ligne, colonne) = somme ; 
			
			    // coefficient de la derivee du Plm dans la solp
			    for (int i=0 ; i<nr ; i++)
				sec_membre->set(ligne) -= i*i/alpha
						*solution_part(num_front+1, k, j, i) ;
							    
			    // on a au moins un diag inferieure dans ce cas ...
			    inf = 1 ;
			}
		    colonne++ ; 
		    }
		    
		    //-----------------------------
		   //        COQUILLES "normales"
		  //------------------------------
		  
		  for (int zone=num_front+2 ; zone<nz-1 ; zone++)
		    if (indic[zone] == 1) {
			
			nr = source.get_mg()->get_nr(zone) ;
			alpha = mapping.get_alpha()[zone] ;
			echelle = mapping.get_beta()[zone]/alpha ;
			
			//Frontiere avec la zone precedente :
			if (indic[zone-1] == 1) ligne -- ;
			
			//valeur de (x+echelle)^l en -1 :
			systeme->set(ligne, colonne) = 
			    -pow(echelle-1, double(l_quant)) ;
			
			// valeur de 1/(x+echelle) ^(l+1) en -1 :
			systeme->set(ligne, colonne+1) = 
			    -1/pow(echelle-1, double(l_quant+1)) ;
			
			// la solution particuliere :
			for (int i=0 ; i<nr ; i++)
			    if (i%2 == 0)
				sec_membre->set(ligne) += solution_part(zone, k, j, i) ;
			    else sec_membre->set(ligne) -= solution_part(zone, k, j, i) ;
			
			// les diagonales :
			sup_new = colonne+1-ligne ;
			if (sup_new > sup) sup = sup_new ;
			
			ligne++ ;
			
			// on prend les derivees si Plm existe dans la zone 
			// precedente :
			
			if (indic[zone-1] == 1) {
			// derivee de (x+echelle)^l en -1 :
			    systeme->set(ligne, colonne) = 
				-l_quant*pow(echelle-1, double(l_quant-1))/alpha ;
			// derivee de 1/(x+echelle)^(l+1) en -1 :
			    systeme->set(ligne, colonne+1) = 
				(l_quant+1)/pow(echelle-1, double(l_quant+2))/alpha ;
			
			
			
			// la solution particuliere :
			    for (int i=0 ; i<nr ; i++)
				if (i%2 == 0)
				    sec_membre->set(ligne) -= 
					i*i/alpha*solution_part(zone, k, j, i) ;
				else
				    sec_membre->set(ligne) +=
					i*i/alpha*solution_part(zone, k, j, i) ;
			
			// les diagonales :
			sup_new = colonne+1-ligne ;
			if (sup_new > sup) sup = sup_new ;
			
			ligne++ ;
			}
			
			// Frontiere avec la zone suivante :
			//valeur de (x+echelle)^l en 1 :
			systeme->set(ligne, colonne) = 
			    pow(echelle+1, double(l_quant)) ;
			
			// valeur de 1/(x+echelle)^(l+1) en 1 :
			systeme->set(ligne, colonne+1) =
			    1/pow(echelle+1, double(l_quant+1)) ;
			    
			// la solution particuliere :
			for (int i=0 ; i<nr ; i++)
			    sec_membre->set(ligne) -= solution_part(zone, k, j, i) ;
			
			// les diagonales inf :
			inf_new = ligne-colonne ;
			if (inf_new > inf) inf = inf_new ;  
			
			ligne ++ ;
			    
			// Utilisation des derivees ssi Plm existe dans la
			//zone suivante :
			if (indic[zone+1] == 1) {
			    
			    //derivee de (x+echelle)^l en 1 :
			    systeme->set(ligne, colonne) = 
				l_quant*pow(echelle+1, double(l_quant-1))/alpha ;
				
			    //derivee de 1/(echelle+x)^(l+1) en 1 :
			    systeme->set(ligne, colonne+1) =
				-(l_quant+1)/pow(echelle+1, double(l_quant+2))/alpha ;
			    
			    // La solution particuliere :
			    for (int i=0 ; i<nr ; i++)
				sec_membre->set(ligne) -= 
				    i*i/alpha*solution_part(zone, k, j, i) ;
				    
			// les diagonales inf :
			inf_new = ligne-colonne ;
			if (inf_new > inf) inf = inf_new ;  
			    
			     }
			 colonne += 2 ;
			 }
		    
		    
		    //--------------------------------
		   //             ZEC
		  //---------------------------------
		  
		  if (indic[nz-1] == 1) {
		      
		      nr = source.get_mg()->get_nr(nz-1) ;
		      
		      
		      alpha = mapping.get_alpha()[nz-1] ;
		      
		      if (indic[nz-2] == 1) ligne -- ;
		      
		      //valeur de (x-1)^(l+1) en -1 :
		      systeme->set(ligne, colonne) = -pow(-2, double(l_quant+1)) ;
		      //solution particuliere :
		      for (int i=0 ; i<nr ; i++)
			if (i%2 == 0)
			    sec_membre->set(ligne) += solution_part(nz-1, k, j, i) ;
			else sec_membre->set(ligne) -= solution_part(nz-1, k, j, i) ;
			
		    //on prend les derivees ssi Plm existe dans la zone precedente :
		    if (indic[nz-2] == 1) {
			
			//derivee de (x-1)^(l+1) en -1 :
			systeme->set(ligne+1, colonne) = 
			    alpha*(l_quant+1)*pow(-2., double(l_quant+2)) ;
			
			// Solution particuliere :
			for (int i=0 ; i<nr ; i++)
			    if (i%2 == 0)
				sec_membre->set(ligne+1) -=
				    -4*alpha*i*i*solution_part(nz-1, k, j, i) ;
			    else
				sec_membre->set(ligne+1) +=
				    -4*alpha*i*i*solution_part(nz-1, k, j, i) ;
			
			//les diags :
			if (sup == 0) sup = 1 ;
			}
		    }
		    
		    //-------------------------
		   //  resolution du systeme 
		  //--------------------------
		    
		    systeme->set_band(sup, inf) ;
		    systeme->set_lu() ;

		    Tbl facteurs(systeme->inverse(*sec_membre)) ;
		    int conte = 0 ;
		    
		 
		// rangement dans la coquille interne :
			
		    if (indic[num_front+1] == 1) {
			nr = source.get_mg()->get_nr(num_front+1) ;
			for (int i=0 ; i<nr ; i++)
			    resultat.set(num_front+1, k, j, i) = 
				 solution_part(num_front+1, k, j, i)
			+facteurs(conte)*solution_hom_un(num_front+1, k, j, i) ;
			conte++ ;
			}
			
			// rangement dans les coquilles :
		    for (int zone=num_front+2 ; zone<nz-1 ; zone++)
			if (indic[zone] == 1) {
			    nr = source.get_mg()->get_nr(zone) ;
			    for (int i=0 ; i<nr ; i++)
			    resultat.set(zone, k, j, i) = 
				solution_part(zone, k, j, i)
				+facteurs(conte)*solution_hom_un(zone, k, j, i) 
				+facteurs(conte+1)*solution_hom_deux(zone, k, j, i) ;
			    conte+=2 ;
			}
			
			//rangement dans la ZEC :
		    if (indic[nz-1] == 1) {
			 nr = source.get_mg()->get_nr(nz-1) ;
			 for (int i=0 ; i<nr ; i++)
			 resultat.set(nz-1, k, j, i) = 
			   solution_part(nz-1, k, j, i)
			    +facteurs(conte)*solution_hom_un(nz-1, k, j, i) ;
			}
		    
		    delete sec_membre ;
		    delete systeme ;
	    }
    
    }
    
    delete [] indic ;
    
    // Les trucs les plus internes sont mis a zero ...
    for (int l=0 ; l<num_front+1 ; l++)
	resultat.t[l]->set_etat_zero() ;
    return resultat;
}
}
