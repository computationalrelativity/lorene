/*
 * Test code for LORENE class Valeur and PGPLOT
 */
 
/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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
 * Revision 1.2  2002/10/16 14:37:19  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "valeur.h"
#include "map.h"
#include "graphique.h"

int main(){

  // Nombre de points
  // ----------------
  int np ; 
  cout << "Nombre de points ? " << endl ; 
  cin >> np ; 


  // Construction de la grille
  // -------------------------
  int nz = 1 ; 
  int nbr[] = {1} ; 
  int type_r[] = {RARE} ; 
  int nbt[] = {1} ; 
  int type_t = SYM ; 
  int nbp[] = {np} ; 
  int type_p = NONSYM ;
  
  Mg3d mg(nz, nbr, type_r, nbt, type_t, nbp, type_p) ;
  
  // Construction du mapping associe
  // -------------------------------

  double r_limits[] = {0,1} ; 
  Map_af mp(mg,r_limits) ; 

  // Phi
  // ---

  const Coord& phi = mp.phi ; 

  // Valeur de la fonction
  // ---------------------

  Valeur ff(mg) ; 
  ff = pow(cos(phi),5) ; 

  cout << "f : " << endl ; 
  cout << ff << endl ; 

  // Transformation de Fourier
  // -------------------------

  ff.std_base_scal() ; // definit la base spectrale a utiliser (Fourier)

  ff.coef() ; // effectue la transformation de Fourier

  cout << "Coefficients de la transformation de Fourier de f : " << endl ; 

  ff.affiche_seuil(cout) ; 

  // Dessin des coefficients de Fourier
  // ----------------------------------
  des_coef_phi(ff, 0, 0, 0, 1.e-14, "log|C_k|","Fourier coefficients") ; 


  // Fin du programme
  // ----------------

  return EXIT_SUCCESS ; 
 
}
