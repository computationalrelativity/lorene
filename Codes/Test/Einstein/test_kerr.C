/*
 *  Test code for LORENE: given the Kerr metric in quasi-isotropic
 *  gauge, it tests it within the 3+1 formalism, in cartesian components.
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

char test_kerr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:37:17  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/09/11 08:31:26  j_novak
 * *** empty log message ***
 *
 * Revision 1.1  2002/09/10 13:23:19  j_novak
 * A test code for LORENE on Kerr metric is now available.
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
#include "tenseur.h"
#include "metconf.h" 
#include "graphique.h"
#include "proto.h"
#include "nbr_spx.h"

int main() {

  double M = 1. ;
  double a = 0.99 ;
  
  int nz = 2 ;
  double R = 2*M ;
  
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
    nt[l] = 9 ;
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
  bornes[0] = R ;
  for (int i=1 ; i<nz ; i++)
    bornes[i] = R*double(i+nz)/double(nz) ;
  bornes[nz] = __infinity ;
  Map_af mapping(grille, bornes) ;
  
  Coord& r = mapping.r ;
  Coord& sint = mapping.sint ;
  Coord& cost = mapping.cost ;
  Mtbl erre(1.+M/r + (M*M - a*a)/(4*r*r)) ;
  Mtbl rs(r + M + (M*M - a*a)/(4*r)) ;
  Mtbl sigmasr(rs + (a*a*cost*cost)/rs) ;
  Cmp adeux(mapping) ; Cmp bdeux(mapping) ;
  adeux = erre*erre + (a*a*cost*cost)/(r*r) ; //Kerr dans la jauge quasi-
  bdeux = erre*erre + a*a/(r*r) + 2*a*a*sint*sint*M/(sigmasr*r*r) ; //isotrope
  
  Tenseur_sym gspher(mapping, 2, COV, mapping.get_bvect_spher() ) ;
  gspher.set_etat_qcq() ;
  for (int i=0; i<3; i++) 
    for (int j=0; j<i; j++) 
      gspher.set(i,j) = 0 ;
  gspher.set(0,0) = adeux ;
  gspher.set(1,1) = adeux ;
  gspher.set(2,2) = bdeux ;
  gspher.set_std_base() ;
  
  Tenseur shift(mapping, 1, CON, mapping.get_bvect_spher() ) ;
  shift.set_etat_qcq() ;
  shift.set(0) = 0 ;
  shift.set(1) = 0 ; 
  Cmp bphi(mapping) ;
  bphi = 2*a*M*sint/(sigmasr*(rs*erre + a*a/r) + 2*a*a*M*sint*sint/r) ;
  shift.set(2) = bphi ;
  shift.set_std_base() ;

  shift.change_triad(mapping.get_bvect_cart()) ;
  shift.set_std_base() ;

  Tenseur lapse(mapping) ;
  lapse.set_etat_qcq() ;
  Cmp enne2(mapping) ;
  enne2 = 1 - 2*M/sigmasr + 
    1./(sigmasr*sigmasr*(rs*rs + a*a)/(4*a*a*M*M*sint*sint) 
	+ sigmasr/(2*M)) ;
  lapse.set() = sqrt(enne2 ) ;
  lapse.set_std_base() ;
  
  Tenseur_sym gcart(gspher) ;
  gcart.change_triad(mapping.get_bvect_cart() ) ;
  gcart.set_std_base() ;
  
  Metrique gij(gcart) ;
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

  cout << "----------------------" << endl ;
  cout << "D_k(metrique conforme)" << endl ;
  cout << "----------------------" << endl ;  
  Tenseur_sym dcov(gtilde.cov().derive_cov(gij)) ;
  Tenseur_sym dcon(gtilde.con().derive_cov(gij)) ;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) 
       for (int k=j; k<3; k++) {
	 cout << "i, j, k: " << i << ", " << j << ", " << k << endl ; 
	 cout << "Rep. covariante: " << max(abs(dcov(i,j,k))) ;
	 cout << "Rep. contravariante: " << max(abs(dcon(i,j,k))) ;
 	cout << "-----------------------------------------" << endl ;
       }
  
  cout << "---------------------" << endl ;
  cout << "Norme de la trace K  " << endl ;
  cout << "---------------------" << endl ;
  Tenseur_sym Kij(0.5*lie_derive(gij.cov(), shift)/lapse) ;
  Tenseur tra(contract(gij.con(), 0, Kij, 0)) ;
  Cmp trace(mapping) ; trace.set_etat_zero() ;
  for (int i=0; i<3; i++) 
    trace += tra(i,i) ;
  cout << norme(trace)/norme(Kij(0,0)) ;
  arrete() ;
  
  Tenseur trK(trace) ;
  trK.dec_dzpuis() ;

  cout << "-------------------------------" << endl ;
  cout << "Verification du calcul de Kij  " << endl ;
  cout << "-------------------------------" << endl ;
  Tenseur betad(manipule(shift, gij, 0)) ;
  Tenseur dbeta(betad.derive_cov(gij)) ;
  Tenseur ver(0.5*dbeta/lapse) ;
  for (int i=0; i<3; i++) 
    for (int j=i; j<3; j++) {
      cout << "i,j : " << i << ", " << j << endl ;
      cout << norme(Kij(i,j) - ver(i,j) - ver(j,i)) ;
    }
  arrete() ;

  cout << "**********************************************" << endl ;
  cout << "*** Verifications des equations d'Einstein ***" << endl ;
  cout << "***    en 3+1 ecrites avec la derivee D~   ***" << endl ;
  cout << "**********************************************" << endl ;

  cout << "-------------------------------" << endl ;
  cout << "Membre de droite des equations " << endl ;
  cout << "       d'evolution de Kij      " << endl ;
  cout << "-------------------------------" << endl ;
  Tenseur h13(pow(gij.determinant(), -1./3.)) ;
  Tenseur_sym Aij(h13*Kij) ;
  Tenseur_sym Aij_copie(Aij) ;
  Aij_copie.dec2_dzpuis() ;
  Tenseur_sym L1(lie_derive(Aij_copie, shift)) ;
  L1.inc_dzpuis() ;
  Tenseur_sym KK(contract(contract(gtilde.con(),0, Aij,1),0, Aij, 1)) ;
  KK.dec_dzpuis() ;
  Tenseur dN(lapse.derive_cov(gij)) ;
  dN.dec2_dzpuis() ;
  Tenseur_sym ddN(dN.derive_cov(gij)) ;
  ddN.inc_dzpuis() ;
  Tenseur_sym L2(lapse*(h13*sans_trace(gij.ricci(),gij) - 2*KK + trK*Aij)
		 + L1) ;
  Tenseur_sym dKdt(L2 - h13*sans_trace(ddN,gij)) ;

  Tenseur_sym Kij_copie(Kij) ;
  Kij_copie.dec2_dzpuis() ;
  Tenseur_sym L10(lie_derive(Kij_copie, shift)) ;
  L10.inc_dzpuis() ;
  Tenseur_sym KK0(contract(contract(gij.con(), 0, Kij, 1), 0, Kij, 0)) ;
  KK0.dec_dzpuis() ;
  Tenseur_sym L20(lapse*(gij.ricci() - 2*KK0 + trK*Kij) + L10) ;
  Tenseur dN0(lapse.derive_cov(gij)) ;
  dN0.dec2_dzpuis() ;
  Tenseur_sym ddN0(dN0.derive_cov(gij)) ;
  ddN0.inc_dzpuis() ;
  ddN0 = ddN0 - L20 ;

  for (int i=0; i<3; i++) 
    for (int j=i; j<3; j++) {
      cout << "i,j: " << i << ", " << j << endl ;
      cout << "Formulation avec ~Aij:" << max(abs(dKdt(i,j)))
	/max(abs(L2(i,j))) ;
      cout << "Formulation avec Kij:" <<max(abs(ddN0(i,j)))
	/max(abs(L20(i,j))) ;
    }
  arrete() ;

  cout << "--------------------------------" << endl ;
  cout << "Eq. de contrainte hamiltonienne " << endl ;
  cout << "--------------------------------" << endl ;
  Tenseur AA(Aij.carre_scal(gtilde)) ;
  AA.dec_dzpuis() ;
  Tenseur Ham(gij.ricci_scal() - AA) ;
  Ham.dec_dzpuis() ; Ham.dec2_dzpuis() ;
  AA.dec_dzpuis() ; AA.dec2_dzpuis() ;
  cout << max(abs(Ham()))/max(abs(AA())) ;
  arrete() ;

  cout << "---------------------------------" << endl ;
  cout << "Eq. de contrainte impulsionnelle " << endl ;
  cout << "---------------------------------" << endl ;
  Tenseur_sym Aup(contract(contract(Aij,0,gtilde.con(),0), 0, gtilde.con(), 0)) ;
  Aup.dec2_dzpuis() ;
  Tenseur_sym dA(Aup.derive_cov(gij)) ;
  Tenseur divA(contract(dA,0,2)) ; divA.dec2_dzpuis() ;
  for (int i=0; i<3; i++) {
    cout << "i:" << i << max(divA(i)) ;
  }

  delete [] np ;
  delete [] nt ;
  delete [] nr0 ;
  delete [] type_r ;
  delete [] bornes ;
  return EXIT_SUCCESS ; 
}    






