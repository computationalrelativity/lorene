/*
 *  Test code for the classe Metconf and the associated covariant 
 *  derivatives, using Schwarzschild metric.
 *
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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

char metconf_derive_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:37:18  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/09/10 14:16:14  j_novak
 * Modif commentaires
 *
 *
 * $Header$
 *
 */

//standard
#include <math.h>

// Headers Lorene :
#include "utilitaires.h"
#include "param.h"
#include "metconf.h" 
#include "graphique.h"
#include "proto.h"
#include "nbr_spx.h"

int main() {

    int nz = 2 ;
    double R = 3. ;

    // Nombre points en r 
    int nr = 21 ;
    
    // echantillonnage en phi :
    int* np = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	np[l] = 4 ;
    int type_p = SYM ;
    
    // echantillonnage en theta :
    int* nt = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	nt[l] = 5 ;
    int type_t = SYM ;

    int* type_r = new int[nz] ;
    type_r[0] = FIN ;
    for (int l=1 ; l<nz ; l++)
	type_r[l] = FIN ;
    type_r[nz-1] = UNSURR ;
    
   // Construction de la grille 

    // echantillonage en r :
    int* nr0 = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
      nr0[l] = nr ;
	
    Mg3d grille (nz, nr0, type_r, nt, type_t, np, type_p) ;
    cout << grille ;

    //Construction du mapping :
    double* bornes = new double[nz+1] ;
    for (int i=0 ; i<nz ; i++)
	bornes[i] = R*double(i+nz)/double(nz) ;
    bornes[nz] = __infinity ;
    Map_af mapping(grille, bornes) ;

    Coord& r = mapping.r ;
    Cmp source(mapping) ;
    source = 1. - 2/r ; //Schwarzschild avec M = 1

    Tenseur_sym gspher(mapping, 2, COV, mapping.get_bvect_spher() ) ;
    gspher.set_etat_qcq() ;
    for (int i=0; i<3; i++) 
      for (int j=0; j<i; j++) 
	gspher.set(i,j) = 0 ;
    gspher.set(0,0) = 1./source ;
    gspher.set(1,1) = 1. ;
    gspher.set(2,2) = 1. ;
    gspher.set_std_base() ;

    Tenseur_sym gcart(gspher) ;
    gcart.change_triad(mapping.get_bvect_cart() ) ;
    gcart.set_std_base() ;

    Metrique gij(gcart) ;
    cout << "----------------------------------------------------" << endl ;
    cout << "Norme du scalaire de courbure de la 3-metrique gij: " << endl ;
    cout << "----------------------------------------------------" << endl ;
    Cmp courb(gij.ricci_scal()()) ;
    courb.dec2_dzpuis() ;
    courb.dec_dzpuis() ;
    cout << norme(courb) ;
    arrete() ;  

    Tenseur_sym plat(mapping, 2, COV, mapping.get_bvect_cart() ) ;
    plat.set_etat_qcq() ;
    for (int i=0; i<3; i++) {
      for (int j=0; j<i; j++) 
        plat.set(i,j) = 0 ;
      plat.set(i,i) = 1 ;
    }
    plat.set_std_base() ;

    Metrique fij(plat, true) ;
    Metconf gtilde(pow(gij.determinant(),(-1./3.))*gcart, fij) ;

   cout << "--------------------" << endl ;
   cout << "Norme de D_k g~_ij: " << endl ;
   cout << "--------------------" << endl ;
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) 
	for (int k=j; k<3; k++) {
	  cout << "k, i, j : " << i << ", " << j << ", " << k << endl ;
	  cout << norme(gtilde.cov().derive_cov(gij)(i,j,k)) ;
	}
    arrete() ;

   cout << "---------------------" << endl ;
   cout << "Norme de D~_k g~_ij: " << endl ;
   cout << "---------------------" << endl ;
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) 
	for (int k=j; k<3; k++) {
	  cout << "k, i, j : " << i << ", " << j << ", " << k << endl ;
	  cout << norme(gtilde.cov().derive_cov(gtilde)(i,j,k)) ;
	}
    arrete() ;

    cout << "------------" << endl ;
    cout << "Delta^k_ik: " << endl ;
    cout << "------------" << endl ;
    Tenseur beta(contract(gtilde.delta(),0,1)) ;
    for (int i=0; i<3; i++) {
      cout << "i: " << i <<  endl ;
      cout << norme(beta(i)) ;
    }
    arrete() ;

    cout << "---------------------------------------------------------" 
	 << endl ;
    cout << "Erreur relative sur  C~ - gamma + gamma~ [formule (169)]: " 
	 << endl ;
    cout << "---------------------------------------------------------" 
	 << endl ;
    Tenseur T1(contract(gij.con(),1,gij.cov().derive_cov(gtilde),1)) ;
    Tenseur T2(contract(gij.con(),1,gij.cov().derive_cov(gtilde),0)) ;
    Tenseur_sym ctilde(gij.gamma()) ;
    for (int k=0; k<3; k++) 
      for (int i=0; i<3; i++) 
	for (int j=0; j<3; j++) {
	  ctilde.set(k,i,j) = 0.5*(T1(k,i,j) + T1(k,j,i) - T2(k,i,j)) ;
	}

    Tenseur diff(ctilde - gij.gamma() + gtilde.gamma()) ;
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) 
	for (int k=j; k<3; k++) {
	  cout << "k, i, j : " << i << ", " << j << ", " << k << endl ;
	  cout << norme(diff(i,j,k))/norme(ctilde(i,j,k)) ;
	}
    arrete() ;

    cout << "----------------------------------------------------" << endl ;
    cout << "Verification de la relation entre Rij et R~ij (178):" << endl ;
    cout << "----------------------------------------------------" << endl ;
    Tenseur lnh(log(gij.determinant())) ;
    lnh.set_std_base() ;
    Tenseur gradh(lnh.gradient()) ;
    Tenseur gradh2(gradh) ;
    gradh2.dec2_dzpuis() ;
    Tenseur_sym ggh(gradh2.derive_cov(gtilde)) ;
    ggh.inc_dzpuis() ;

    Tenseur_sym R1(fij.cov()) ;
    for(int i=0; i<3; i++) 
      for (int j=i; j<3; j++) 
	R1.set(i,j) = gradh(i)*gradh(j)/6. ;
    R1.dec_dzpuis() ;
    Tenseur Raux(contract(gtilde.con(),0,ggh+R1,0)) ;
    Tenseur_sym R2(contract(Raux,0,1)*gtilde.cov()) ;
    Tenseur_sym Rij(gij.ricci() - gtilde.ricci() + (ggh - R1 + R2)/6.) ;
    for(int i=0; i<3; i++) 
      for (int j=i; j<3; j++) {
	cout << "Erreur relative/composante: " << i << ", " << j << endl ;
	cout << norme(Rij(i,j))/norme(gij.ricci()(i,j)) ;}
    arrete() ;

    cout << "------------------------------------------------" << endl ;
    cout << "Verification de la relation entre R et R~ (182):" << endl ;
    cout << "------------------------------------------------" << endl ;
    Raux = contract(gij.con(), 0,(ggh - R1 + R2)/6., 0) ;
    Tenseur R3(contract(Raux,0,1)) ;
    Tenseur R4(gtilde.ricci_scal()) ;
    if (gtilde.dirac_gauge()) R4.dec_dzpuis() ;
    Tenseur RR(gij.ricci_scal() - pow(gij.determinant(), -1./3.)*R4 + R3) ;
    cout << norme(RR())/norme(gtilde.ricci_scal()()) ;
    arrete() ;
    
    cout << "--------------------------------------------" << endl ;
    cout << "Verification de la formule (199) pour R~_ij:" << endl ;
    cout << "--------------------------------------------" << endl ;
    Tenseur_sym L1(contract(contract(gtilde.delta(),0, gtilde.delta(),2),
			    1,2) );
    L1.dec_dzpuis() ;
    Tenseur_sym del2(gtilde.delta()) ;
    del2.dec2_dzpuis() ;
    Tenseur_sym L2(contract(del2.derive_cov(fij), 0, 1)) ;
    L2.inc_dzpuis() ;

    Tenseur_sym fin(gtilde.ricci() - L2 + L1) ;
    for (int i=0; i<3; i++) 
      for (int j=i; j<3; j++) {
	cout << "i, j: " << i << ", " << j << endl ;
	cout << norme(fin(i,j))/norme(gtilde.ricci()(i,j)) ;
      }

    delete [] np ;
    delete [] nt ;
    delete [] nr0 ;
    delete [] type_r ;
    delete [] bornes ;
    return EXIT_SUCCESS ; 
}    






