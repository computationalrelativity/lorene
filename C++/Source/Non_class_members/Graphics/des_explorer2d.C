/*
 * Prepares a file for an Explorer visualisation of a "slice" (at $\varphi$ =
 * const) of a Cmp in a given domain. The result has to be read by the
 * Lit_scal2D module of Explorer (see graphique.h).
 *
 */

/*
 *   Copyright (c) 2000-2001 Jerome Novak
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


char des_explorer2d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:57  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/12/04  14:15:44  novak
 * Phi = NONSYM case added
 *
 * Revision 1.1  2000/11/27 15:16:33  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

// LORENE
#include "cmp.h" 

void des_explorer2D(const Cmp& uu, int nz, const int k_phi,
		     const char* filename, const double scale = 1.)
{
  const Map* mp = uu.get_mp() ;
  const Mg3d* grille = mp->get_mg() ;
  int np = grille->get_np(nz) ;
  
  assert (nz < grille->get_nzone()) ;
	  
  assert (grille->get_type_r(nz) != UNSURR) ;
  assert (grille->get_type_t() == SYM) ;
  assert (k_phi < np) ;
  if (np%2 != 0) cout << "des_explorer2D: Warning: np is not even!" ;
  
  // Definition of the grid
  Mtbl erre (mp->r) ;
  Mtbl theta (mp->tet) ;
  Mtbl xx = erre*sin(theta) ;
  Mtbl yy = erre*cos(theta) ;
  
  int nr = grille->get_nr(nz) ;
  int nt = grille->get_nt(nz) ;  
  
  int nt4 = 4*(nt-1) + 1 ;
  
  int ntot = nr * nt4 ;
  
  float* x_exp = new float[ntot] ;
  float* y_exp = new float[ntot] ;
  float* z_exp = new float[ntot] ;
  
  float* fuu = new float[ntot] ;
  
  float* px = x_exp ;
  float* py = y_exp ;
  float* pz = z_exp ;
  float* puu = fuu ;
  
  for (int j=0; j<nt; j++) {
    for (int i=0; i<nr; i++) {
      *px = xx(nz, k_phi, j, i) ;
      *py = yy(nz, k_phi, j, i) ;
      *pz = scale*uu(nz, k_phi, j, i) ;
      
      *puu = uu(nz, k_phi, j, i) ;
      
      px++ ;
      py++ ;
      pz++ ;
      puu++ ;
    }
  }
  
  for (int j=nt-2; j>=0; j--) {
    for (int i=0; i<nr; i++) {
      *px =  xx(nz, k_phi, j, i) ;
      *py = -yy(nz, k_phi, j, i) ;
      *pz = scale*uu(nz, k_phi, j, i) ;
      
      *puu = uu(nz, k_phi, j, i) ;
      
      px++ ;
      py++ ;
      pz++ ;
      puu++ ;
    }
  }
  
  if (grille->get_type_p() == SYM) {
    for (int j=1; j<nt; j++) {
      for (int i=0; i<nr; i++) {
	*px = - xx(nz, k_phi, j, i) ;
	*py = - yy(nz, k_phi, j, i) ;
	*pz = scale*uu(nz, k_phi, j, i) ;
	
	*puu = uu(nz, k_phi, j, i) ;
	
	px++ ;
	py++ ;
	pz++ ;
	puu++ ;
      }
    }
    
    for (int j=nt-2; j>=0; j--) {
      for (int i=0; i<nr; i++) {
	*px = - xx(nz, k_phi, j, i) ;
	*py = yy(nz, k_phi, j, i) ;
	*pz = scale*uu(nz, k_phi, j, i) ;
	
	*puu = uu(nz, k_phi, j, i) ;
	
	px++ ;
	py++ ;
	pz++ ;
	puu++ ;
      }
    }
  }

  else {
    int k_phi2 = k_phi+np/2 ;
    for (int j=1; j<nt; j++) {
      for (int i=0; i<nr; i++) {
	*px = -xx(nz, k_phi2, j, i) ;
	*py = -yy(nz, k_phi2, j, i) ;
	*pz = scale*uu(nz, k_phi2, j, i) ;
	
	*puu = uu(nz, k_phi2, j, i) ;
	
	px++ ;
	py++ ;
	pz++ ;
	puu++ ;
      }
    }
    
    for (int j=nt-2; j>=0; j--) {
      for (int i=0; i<nr; i++) {
	*px = -xx(nz, k_phi2, j, i) ;
	*py = yy(nz, k_phi2, j, i) ;
	*pz = scale*uu(nz, k_phi2, j, i) ;
	
	*puu = uu(nz, k_phi2, j, i) ;
	
	px++ ;
	py++ ;
	pz++ ;
	puu++ ;
      }
    }
  }
    
  
  FILE* fich = fopen(filename, "w" ) ;
  fwrite (&nr, sizeof(int), 1, fich);
  fwrite (&nt4, sizeof(int), 1, fich);
  for (int i=0 ; i<ntot ; i++ ) {
    fwrite (&(x_exp[i]), sizeof(float), 1, fich);
    fwrite (&(y_exp[i]), sizeof(float), 1, fich);
    fwrite (&(z_exp[i]), sizeof(float), 1, fich);
  }
  fwrite (fuu, sizeof(float), ntot, fich);
  fclose(fich) ;
  
  delete[] x_exp ;
  delete[] y_exp ;
  delete[] z_exp ;
  delete[] fuu ;
      
  
}
