/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/02/06 10:53:54  j_novak
 * New dzpuis = 5 -> dzpuis = 3 case (not ready yet).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.24  2000/07/18  10:23:09  eric
 * Suppression d'une erreur sans consequence : dans le noyau, remplacement de
 *   beta = mapping.get_alpha()[0]    par
 *   beta = mapping.get_beta()[0]
 * Cette erreur etait sans consequence car beta n'intervient pas dans la
 * suite des calculs pour le noyau.
 *
 * Revision 2.23  2000/05/22  13:43:46  phil
 * ajout du cas dzpuis =3
 *
 * Revision 2.22  2000/02/09  14:40:49  eric
 * Ajout de l'argument dzpuis a sol_poisson.
 * (dzpuis n'est plus lu sur le Mtbl_cf).
 *
 * Revision 2.21  1999/12/02  14:29:05  eric
 * Suppression de la base en argument de la version Mtbl_cf.
 * La version Map se trouve desormais dans le fichier map_af_poisson.C
 * La version Cmp se trouve desormais dans le fichier cmp_pde.C
 *
 * Revision 2.20  1999/11/30  12:56:45  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.19  1999/11/24  14:34:25  eric
 * Accession des membres alpha et beta de Map_af (desormais prives) par
 *   Map_af::get_alpha() et Map_af::get_beta().
 *
 * Revision 2.18  1999/10/27  15:47:19  eric
 * Suppression du membre Cmp::c.
 *
 * Revision 2.17  1999/10/14  14:28:07  eric
 * Methode const.
 *
 * Revision 2.16  1999/10/13  15:52:22  eric
 * Ajout de la base dans l'appel au constructeur de Mtbl_cf.
 *
 * Revision 2.15  1999/10/11  16:27:59  phil
 * suppression du sur echantillonnage
 *
 * Revision 2.14  1999/10/11  14:28:50  phil
 * & -> &&
 *
 * Revision 2.13  1999/09/06  16:25:42  phil
 * ajout de la version Cmp
 *
 * Revision 2.12  1999/07/02  15:05:01  phil
 * *** empty log message ***
 *
 * Revision 2.11  1999/06/23  12:35:18  phil
 * ajout de dzpuis = 2
 *
 * Revision 2.10  1999/04/27  13:12:28  phil
 * *** empty log message ***
 *
 * Revision 2.9  1999/04/14  13:56:03  phil
 * *** empty log message ***
 *
 * Revision 2.8  1999/04/14  10:22:28  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/04/14  10:20:09  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/04/13  17:18:50  phil
 * Ajout de la version Valeur
 *
 * Revision 2.5  1999/04/12  15:21:37  phil
 * correction de la construction du systeme
 *
 * Revision 2.4  1999/04/12  14:55:16  phil
 * correction matrices
 *
 * Revision 2.3  1999/04/12  14:28:46  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/07  14:56:33  phil
 * Changement prototypage
 *
 * Revision 2.1  1999/04/07  14:40:55  phil
 * mise au point des includes.
 *
 * Revision 2.0  1999/04/07  14:10:55  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header C : 
#include <stdlib.h>
#include <math.h>

// Headers Lorene :
#include "matrice.h"
#include "map.h"
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
 *	    dzpuis : exposant de r dans le factor multiplicatif dans la ZEC
 * Sortie : renvoie les coefficients de la solution dans la meme base de 
 *	    decomposition (a savoir Ylm)
 *	    
 */



Mtbl_cf sol_poisson(const Map_af& mapping, const Mtbl_cf& source, int dzpuis)
{
    
    // Verifications d'usage sur les zones
    int nz = source.get_mg()->get_nzone() ;
    assert (nz>1) ;
    assert (source.get_mg()->get_type_r(0) == RARE) ;
    assert (source.get_mg()->get_type_r(nz-1) == UNSURR) ;
    for (int l=1 ; l<nz-1 ; l++)
	assert(source.get_mg()->get_type_r(l) == FIN) ;

     assert ((dzpuis==4) || (dzpuis==2) || (dzpuis==3) || (dzpuis==5)) ;
    
    
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
	    if (dzpuis !=5) 
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
		if (dzpuis !=5) {
		  solution_hom_un.set(nz-1, k, j, i) = (*sol_hom)(i) ;
		  solution_hom_deux.set(nz-1, k, j, i) = 0. ;
		}
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

    for (int zone=1 ; zone<nz-1 ; zone++) {
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


    if (dzpuis !=5) {

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
		taille = indic[nz-1]+indic[0] ;
		for (int zone=1 ; zone<nz-1 ; zone++)
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
			systeme->set(ligne, colonne) = 1. ;
			// coefficient du Plm dans la solp
			for (int i=0 ; i<nr ; i++)
			    sec_membre->set(ligne) -= solution_part(0, k, j, i) ;
			    
			ligne++ ;
			// on prend les derivees que si Plm existe 
			//dans la zone suivante
			
			if (indic[1] == 1) {
			    // derivee de x^l en 1 :
			    systeme->set(ligne, colonne) = 1./alpha*l_quant ; 
			
			    // coefficient de la derivee du Plm dans la solp
			    if (base_r == R_CHEBP)
				// cas en R_CHEBP 
				for (int i=0 ; i<nr ; i++)
				    sec_membre->set(ligne) -= 
					    4*i*i/alpha
						*solution_part(0, k, j, i) ;
			    else
				// cas en R_CHEBI
				for (int i=0 ; i<nr ; i++)
				    sec_membre->set(ligne) -=
					(2*i+1)*(2*i+1)/alpha
					    *solution_part(0, k, j, i) ;
			    		    
			    // on a au moins un diag inferieure dans ce cas ...
			    inf = 1 ;
			}
		    colonne++ ; 
		    }
		    
		    //-----------------------------
		   //        COQUILLES
		  //------------------------------
		  
		  for (int zone=1 ; zone<nz-1 ; zone++)
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
		    
		 
		// rangement dans le noyau :
			
		    if (indic[0] == 1) {
			nr = source.get_mg()->get_nr(0) ;
			for (int i=0 ; i<nr ; i++)
			    resultat.set(0, k, j, i) = solution_part(0, k, j, i)
				+facteurs(conte)*solution_hom_un(0, k, j, i) ;
			conte++ ;
			}
			
			// rangement dans les coquilles :
		    for (int zone=1 ; zone<nz-1 ; zone++)
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

    } 
    else { //Si dzpuis =5, seule la sol. part est renvoyee

      for (int zone = 0; zone < nz; zone++) 
	for (int k=0 ; k<k<source.get_mg()->get_np(zone)+1 ; k++)
	  for (int j=0 ; j<source.get_mg()->get_nt(zone) ; j++)
	    if (nullite_plm(j,source.get_mg()->get_nt(zone) , 
			    k, source.get_mg()->get_np(zone), base) == 1) 
	      for (int i=0; i<source.get_mg()->get_nr(zone); i++)
		resultat.set(zone, k, j, i) = solution_part(zone, k, j, i) ;

    }

      


    return resultat;
}
