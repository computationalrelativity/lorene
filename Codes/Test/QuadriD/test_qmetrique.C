/*
 * Main code for testing the classes Qmetrique and Qtenseur.
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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

char test_qmetrique_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:12:55  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/01/09 11:07:59  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1  2002/09/19 11:41:29  j_novak
 * Added main program for testing the classes Qmetrique and Qtenseur
 *
 *
 * $Header$
 *
 */

//standard
#include <cmath>

// Headers Lorene :
#include "qmetrique.h"
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
  double dt = 0.1 ;
  Param para ;
  para.add_double(dt) ;
  Qmetrique kerr(lapse, shift, gij) ;
  cout << "---------------------------------" << endl ;
  cout << "4-Ricci tensor of the Kerr metric" << endl ;
  cout << "---------------------------------" << endl ;
  kerr.avance_temps(lapse, shift, gij, dt) ;
  kerr.avance_temps(lapse, shift, gij, dt) ;
  kerr.avance_temps(lapse, shift, gij, dt) ;
  for (int i=0; i<4; i++)
    for (int j=i; j<4; j++) {
      cout << "i, j " << i << ", " << j << endl ; 
      cout << max(abs(kerr.ricci(dt)(i,j))) ; 
    }

  delete [] np ;
  delete [] nt ;
  delete [] nr0 ;
  delete [] type_r ;
  delete [] bornes ;
  return EXIT_SUCCESS ; 
}    






