/*
 * Methods for the class Grille_val, and its derivative classes.
 *
 * See the file grille_val.h for documentation
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


char Grille_val_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/12/19 15:05:14  j_novak
 * Trying to avoid shadowed variables
 *
 * Revision 1.3  2003/10/03 16:17:17  j_novak
 * Corrected some const qualifiers
 *
 * Revision 1.2  2001/12/04 21:27:54  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1  2001/11/22 13:41:54  j_novak
 * Added all source files for manipulating Valencia type objects and making
 * interpolations to and from Meudon grids.
 *
 *
 * $Header$
 *
 */

// Fichier includes
#include "grille_val.h"
#include "utilitaires.h"

			//---------------//
			// Constructeurs //
			//---------------//

// Fonction auxilliaire
Tbl* Grille_val::fait_grille1D(const double rmin, const double rmax, const
			       int n) 
{
  assert(rmin<rmax) ;
  Tbl* resu = new Tbl(n) ;
  double step = (rmax - rmin)/double(n-1) ;
  resu->set_etat_qcq() ;
  for (int i=0; i<n; i++) resu->set(i) = rmin + i*step ;
  return resu ;
}

// Constructeur 1D
Grille_val::Grille_val(const double izrmin, const double izrmax, const
		       int n1, const int fantome): 
  dim(n1), nfantome(fantome), type_t(SYM), type_p(SYM)
{
  assert (n1 > 0) ;
  zrmin = new double(izrmin) ;
  zrmax = new double(izrmax) ;
  double amin = ((n1 + nfantome - 0.5)*izrmin -
		 (nfantome-0.5)*izrmax) / n1 ;
  double amax = ((n1 + nfantome - 0.5)*izrmax -
		 (nfantome-0.5)*izrmin) / n1 ;
  zr = fait_grille1D(amin, amax, n1+2*fantome) ;
  amin = ((n1 + nfantome)*izrmin - nfantome*izrmax) / n1 ;
  amax = ((n1 + nfantome)*izrmax - nfantome*izrmin) / n1 ;
  zri = fait_grille1D(amin, amax, n1+2*fantome+1) ;
}

// Constructeur 2D
Grille_val::Grille_val(const double izrmin, const double izrmax, const
		       int n2, const int n1, const int itype_t, 
		       const int fantome): 
  dim(n2,n1), nfantome(fantome), type_t(itype_t), type_p(SYM)

{
  zrmin = new double(izrmin) ;
  zrmax = new double(izrmax) ;
  double amin = ((n1 + nfantome - 0.5)*izrmin -
		 (nfantome-0.5)*izrmax) / n1 ;
  double amax = ((n1 + nfantome - 0.5)*izrmax -
		 (nfantome-0.5)*izrmin) / n1 ;
  zr = fait_grille1D(amin, amax, n1+2*fantome) ;
  amin = ((n1 + nfantome)*izrmin - nfantome*izrmax) / n1 ;
  amax = ((n1 + nfantome)*izrmax - nfantome*izrmin) / n1 ;
  zri = fait_grille1D(amin, amax, n1+2*fantome+1) ;
  
}

// Constructeur 3D
Grille_val::Grille_val(const double izrmin, const double izrmax, 
		       const int n3, const int n2, const int n1, 
		       const int itype_t, const int itype_p, 
		       const int fantome ): dim(n3,n2,n1), nfantome(fantome), 
  type_t(itype_t), type_p(itype_p) 
{
  zrmin = new double(izrmin) ;
  zrmax = new double(izrmax) ;
  double amin = ((n1 + nfantome - 0.5)*izrmin -
		 (nfantome-0.5)*izrmax) / n1 ;
  double amax = ((n1 + nfantome - 0.5)*izrmax -
		 (nfantome-0.5)*izrmin) / n1 ;
  zr = fait_grille1D(amin, amax, n1+2*fantome) ;
  amin = ((n1 + nfantome)*izrmin - nfantome*izrmax) / n1 ;
  amax = ((n1 + nfantome)*izrmax - nfantome*izrmin) / n1 ;
  zri = fait_grille1D(amin, amax, n1+2*fantome+1) ;
}


// Constructeur par recopie
Grille_val::Grille_val(const Grille_val & titi): dim(titi.dim), 
  nfantome(titi.nfantome), type_t(titi.type_t), type_p(titi.type_p)
{
  assert(titi.zr != 0x0) ;
  assert(titi.zri != 0x0) ;
  assert(titi.zrmin != 0x0) ;
  assert(titi.zrmax != 0x0) ;
  zr = new Tbl(*titi.zr) ;
  zri = new Tbl(*titi.zri) ;
  zrmin = new double(*titi.zrmin) ;
  zrmax = new double(*titi.zrmax) ;

}
	
// Depuis un fichier
Grille_val::Grille_val(FILE* fd):dim(fd) {

  fread_be(&nfantome, sizeof(int), 1, fd) ;		
  fread_be(&type_t, sizeof(int), 1, fd) ;		
  fread_be(&type_p, sizeof(int), 1, fd) ;		
  
  double amin, amax ;
  fread_be(&amin, sizeof(double), 1, fd) ;		
  fread_be(&amax, sizeof(double), 1, fd) ;		
  zrmin = new double(amin) ;
  zrmax = new double(amax) ;
  zr = new Tbl(fd) ;
  zri = new Tbl(fd) ;
  
}

			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Grille_val::~Grille_val() {

  assert(zr != 0x0) ;
  assert(zri != 0x0) ;
  assert(zrmin != 0x0) ; 
  assert(zrmax != 0x0) ; 
  delete zr ;
  delete zri ;
  delete zrmin ;
  delete zrmax ; 

}

			//-------------//
			// Affectation //
			//-------------//

// Depuis une autre Grille_val
void Grille_val::operator=(const Grille_val & titi) {

  dim = titi.dim ;
  nfantome = titi.nfantome ;
  type_t = titi.type_t ;
  type_p = titi.type_p ;

  assert(titi.zr != 0x0) ;
  assert(titi.zri != 0x0) ;
  assert(titi.zrmin != 0x0) ;
  assert(titi.zrmax != 0x0) ;
  
  for (int i=0; i<dim.dim[0]+2*nfantome; i++) 
    zr->t[i] = titi.zr->t[i] ;
  for (int i=0; i<dim.dim[0] + 2*nfantome + 1; i++) 
    zri->t[i] = titi.zri->t[i] ;
  *zrmin = *titi.zrmin ;
  *zrmax = *titi.zrmax ;
}
    
			//------------//
			// Sauvegarde //
			//------------//

// Sauve dans un fchier
void Grille_val::sauve(FILE* fd) const {

  dim.sauve(fd) ;
  fwrite_be(&nfantome, sizeof(int), 1, fd) ;		
  fwrite_be(&type_t, sizeof(int), 1, fd) ;		
  fwrite_be(&type_p, sizeof(int), 1, fd) ;		
  
  fwrite_be(zrmin, sizeof(double), 1, fd) ;		
  fwrite_be(zrmax, sizeof(double), 1, fd) ;		
  
  zr->sauve(fd) ; zri->sauve(fd) ;

}
    
			//------------//
			// Impression //
			//------------//

// Operateurs <<
ostream& operator<<(ostream& o, const Grille_val & titi) {
  titi >> o ;
  return o ;
}

ostream& Grille_val::operator>>(ostream& o) const {
  int ndim = dim.ndim ;
  int nfant = nfantome ;
  o.precision(4);
  o.setf(ios::showpoint);
  o << "*** Grille_val " << ndim << "D" << "   size: " ; 
  for (int i = 0; i<ndim-1; i++) {
    o << dim.dim[ndim-1-i] ;
    if (ndim-i == 3) o << " x " ;
    if (ndim-i == 2) o << " x " ;
  } 
  o << dim.dim[0] << iendl ;
  o << nfant << " hidden cells on each side " << iendl ;
  return o ;
}
    
		    //------------------------------------//
		    //		class Gval_cart		  //
		    //------------------------------------//

/*********************************************************************
 *
 *         Cartesian grid for Godunov-type integration schemes
 *
 *********************************************************************/

			//---------------//
			// Constructeurs //
			//---------------//

// Constructeur 1D
Gval_cart::Gval_cart(const double izmin, const double izmax, const int nz,
		     const int fantome)
  :Grille_val(izmin, izmax, nz, fantome), 
  xmin(0x0), xmax(0x0),
  ymin(0x0), ymax(0x0),
  x(0x0), xi(0x0),
  y(0x0), yi(0x0){
}

// Constructeur 2D
Gval_cart::Gval_cart(const double ixmin, const double ixmax, const 
		     double izmin, const double izmax, const int nx, 
		     const int nz, const int itype_t, const int fantome)
  :Grille_val(izmin, izmax, nx, nz, itype_t, fantome),
   ymin(0x0), ymax(0x0),
   y(0x0), yi(0x0) 
{
  assert ( (type_t!=SYM) || (izmin >= double(0)) ) ;
  
  xmin = new double(ixmin) ;
  xmax = new double(ixmax) ;
  double amin = ((nx + nfantome - 0.5)*ixmin -
	  (nfantome-0.5)*ixmax) / nx ;
  double amax = ((nx + nfantome - 0.5)*ixmax -
	  (nfantome-0.5)*ixmin) / nx ;
  x = fait_grille1D(amin, amax, nx+2*fantome) ;
  amin = ((nx + nfantome)*ixmin - nfantome*ixmax) / nx ;
  amax = ((nx + nfantome)*ixmax - nfantome*ixmin) / nx ;
  xi = fait_grille1D(amin, amax, nx+2*fantome+1) ;
  
}
  
// Constructeur 3D
Gval_cart::Gval_cart(const double iymin, const double iymax, 
		     const double ixmin, const double ixmax, const 
		     double izmin, const double izmax, const int ny,
		     const int nx, const int nz, const int itype_t, 
		     const int itype_p, const int fantome)
  :Grille_val(izmin, izmax, ny, nx, nz, itype_t, itype_p, fantome)
{
  assert ( (type_t!=SYM) || (izmin >= double(0)) ) ;
  assert ( (type_p!=SYM) || (iymin >= double(0)) ) ; 
  
  xmin = new double(ixmin) ;
  xmax = new double(ixmax) ;
  double amin = ((nx + nfantome - 0.5)*ixmin -
	  (nfantome-0.5)*ixmax) / nx ;
  double amax = ((nx + nfantome - 0.5)*ixmax -
	  (nfantome-0.5)*ixmin) / nx ;
  x = fait_grille1D(amin, amax, nx+2*fantome) ;
  amin = ((nx + nfantome)*ixmin - nfantome*ixmax) / nx ;
  amax = ((nx + nfantome)*ixmax - nfantome*ixmin) / nx ;
  xi = fait_grille1D(amin, amax, nx+2*fantome+1) ;
  
  ymin = new double(iymin) ;
  ymax = new double(iymax) ;
  amin = ((ny + nfantome - 0.5)*iymin -
	  (nfantome-0.5)*iymax) / ny ;
  amax = ((ny + nfantome - 0.5)*iymax -
	  (nfantome-0.5)*iymin) / ny ;
  y = fait_grille1D(amin, amax, ny+2*fantome) ;
  amin = ((ny + nfantome)*iymin - nfantome*iymax) / ny ;
  amax = ((ny + nfantome)*iymax - nfantome*iymin) / ny ;
  yi = fait_grille1D(amin, amax, ny+2*fantome+1) ;
 
}

// Constructeur par recopie
Gval_cart::Gval_cart(const Gval_cart& titi)
  :Grille_val(titi) 
{
  if (titi.x != 0x0) x = new Tbl(*titi.x) ;
  if (titi.xi != 0x0) xi = new Tbl(*titi.xi) ;
  if (titi.xmin != 0x0) xmin = new double(*titi.xmin) ;
  if (titi.xmax != 0x0) xmax = new double(*titi.xmax) ;
  if (titi.y != 0x0) y = new Tbl(*titi.y) ;
  if (titi.yi != 0x0) yi = new Tbl(*titi.yi) ;
  if (titi.ymin != 0x0) ymin = new double(*titi.ymin) ;
  if (titi.ymax != 0x0) ymax = new double(*titi.ymax) ;

}
   
// Depuis un fichier
Gval_cart::Gval_cart(FILE* fd)
  : Grille_val(fd)
{
  double amin, amax ;
  if (dim.ndim >= 2) {
    fread_be(&amin, sizeof(double), 1, fd) ;		
    fread_be(&amax, sizeof(double), 1, fd) ;		
    xmin = new double(amin) ;
    xmax = new double(amax) ;
    x = new Tbl(fd) ;
    xi = new Tbl(fd) ;
  }
  if (dim.ndim >= 3) {
    fread_be(&amin, sizeof(double), 1, fd) ;		
    fread_be(&amax, sizeof(double), 1, fd) ;		
    ymin = new double(amin) ;
    ymax = new double(amax) ;
    y = new Tbl(fd) ;
    yi = new Tbl(fd) ;
  }
}
  

			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Gval_cart::~Gval_cart() {

  if (x != 0x0) delete x ;
  if (xi != 0x0) delete xi ;
  if (xmin != 0x0) delete xmin ;
  if (xmax != 0x0) delete xmax ; 
  if (y != 0x0) delete y ;
  if (yi != 0x0) delete yi ;
  if (ymin != 0x0) delete ymin ;
  if (ymax != 0x0) delete ymax ; 
  
}

			//-------------//
			// Affectation //
			//-------------//

// Depuis une autre Grille_val
void Gval_cart::operator=(const Gval_cart& titi) {

  Grille_val::operator=(titi) ;

  if (titi.x != 0x0) *x = *titi.x ;
  if (titi.xi != 0x0) *xi = *titi.xi ;
  if (titi.xmin != 0x0) *xmin = *titi.xmin ;
  if (titi.xmax != 0x0) *xmax = *titi.xmax ;
  if (titi.y != 0x0) *y = *titi.y ;
  if (titi.yi != 0x0) *yi = *titi.yi ;
  if (titi.ymin != 0x0) *ymin = *titi.ymin ;
  if (titi.ymax != 0x0) *ymax = *titi.ymax ;
}

			//------------//
			// Sauvegarde //
			//------------//

// save onto a file
void Gval_cart::sauve(FILE* fd) const {

  Grille_val::sauve(fd) ;

  if (dim.ndim >= 2) {
    fwrite_be(xmin, sizeof(double), 1, fd) ;		
    fwrite_be(xmax, sizeof(double), 1, fd) ;		
    x->sauve(fd) ; xi->sauve(fd) ;
  }
  if (dim.ndim >= 3) {
    fwrite_be(ymin, sizeof(double), 1, fd) ;		
    fwrite_be(ymax, sizeof(double), 1, fd) ;		
    y->sauve(fd) ; yi->sauve(fd) ;
  }

}
    

			//------------//
			// Impression //
			//------------//

// Operateurs <<
ostream& Gval_cart::operator>>(ostream& o) const {

  int ndim = dim.ndim;
  Grille_val::operator>>(o) ;

  o << "*** Cartesian grid ***" << endl ;

  switch (ndim) {

  case 1 : {
    o << "Z nodes: " << iendl ;
    for (int i=0; i<dim.dim[0]; i++) {
      o << zr->set(i+nfantome) << " " ;
    }
    o << endl ;
    break ;
  }
  
  
  case 2 : {
    o << "X nodes: " << iendl ;
    for (int j=0 ; j<dim.dim[1] ; j++) {
      o << " " << x->set(j+nfantome) ;
    }
    o << endl ;
    
    o << "Z nodes: " << iendl ;
    
    for (int i=0 ; i<dim.dim[0] ; i++) {
      o << " " << zr->set(i+nfantome) ;
    }
    o << endl ;
    break ;    
  }
  
  case 3 : {
    o << "Y nodes: " << iendl ;
    for (int k=0 ; k<dim.dim[2] ; k++) {
      o << " " << y->set(k+nfantome) ;
    }
    o << endl ;

    o << "X nodes: " << iendl ;
    for (int j=0 ; j<dim.dim[1] ; j++) {
      o << " " << x->set(j+nfantome) ;
    }
    o << endl ;
    
    o << "Z nodes: " << iendl ;
    
    for (int i=0 ; i<dim.dim[0] ; i++) {
      o << " " << zr->set(i+nfantome) ;
    }
    o << endl ;
    break ;
  }
  
  default : {
    cout << "operator>> Gval_cart : unexpected dimension !" << endl ;
    cout << " ndim = " << ndim << endl ;        
    abort() ;
    break ;
  }
  }
  
  return o ;
}
		    //------------------------------------//
		    //		class Gval_spher	  //
		    //------------------------------------//

/*********************************************************************
 *
 *         Spherical grids for Godunov-type integration schemes
 *
 *********************************************************************/

			//---------------//
			// Constructeurs //
			//---------------//

// Constructeur 1D
Gval_spher::Gval_spher(const double irmin, const double irmax, const int nr,
		     const int fantome)
  :Grille_val(irmin, irmax, nr, fantome), 
  tet(0x0), teti(0x0),
  phi(0x0), phii(0x0){
  assert(irmin>=double(0)) ;
}

// Constructeur 2D
Gval_spher::Gval_spher(const double irmin, const double irmax, const int nt, 
		     const int nr, const int itype_t, const int fantome)
  :Grille_val(irmin, irmax, nt, nr, itype_t, fantome),
   phi(0x0), phii(0x0) 
{
  assert (irmin >= double(0)) ;
  
  double tetmin = 0. ;
  double tetmax = (type_t == SYM ? M_PI_2 : M_PI) ;
  double amin = ((nt + nfantome - 0.5)*tetmin -
	  (nfantome-0.5)*tetmax) / nt ;
  double amax = ((nt + nfantome - 0.5)*tetmax -
	  (nfantome-0.5)*tetmin) / nt ;
  tet = fait_grille1D(amin, amax, nt+2*fantome) ;
  amin = ((nt + nfantome)*tetmin - nfantome*tetmax) / nt ;
  amax = ((nt + nfantome)*tetmax - nfantome*tetmin) / nt ;
  teti = fait_grille1D(amin, amax, nt+2*fantome+1) ;
  
}
  
// Constructeur 3D
Gval_spher::Gval_spher(const double irmin, const double irmax, const int np,
		     const int nt, const int nr, const int itype_t, 
		     const int itype_p, const int fantome)
  :Grille_val(irmin, irmax, np, nt, nr, itype_t, itype_p, fantome)
{
  assert (irmin >= double(0))  ;
  

  double tetmin = 0. ;
  double tetmax = (type_t == SYM ? M_PI_2 : M_PI) ;
  double amin = ((nt + nfantome - 0.5)*tetmin -
	  (nfantome-0.5)*tetmax) / nt ;
  double amax = ((nt + nfantome - 0.5)*tetmax -
	  (nfantome-0.5)*tetmin) / nt ;
  tet = fait_grille1D(amin, amax, nt+2*fantome) ;
  amin = ((nt + nfantome)*tetmin - nfantome*tetmax) / nt ;
  amax = ((nt + nfantome)*tetmax - nfantome*tetmin) / nt ;
  teti = fait_grille1D(amin, amax, nt+2*fantome+1) ;
  
  double phimin = 0. ;
  double phimax = ( type_p == SYM ? M_PI : 2.*M_PI) ; //??? a verifier!
  amin = ((np + nfantome - 0.5)*phimin -
	  (nfantome-0.5)*phimax) / np ;
  amax = ((np + nfantome - 0.5)*phimax -
	  (nfantome-0.5)*phimin) / np ;
  phi = fait_grille1D(amin, amax, np+2*fantome) ;
  amin = ((np + nfantome)*phimin - nfantome*phimax) / np ;
  amax = ((np + nfantome)*phimax - nfantome*phimin) / np ;
  phii = fait_grille1D(amin, amax, np+2*fantome+1) ;
 
}

// Constructeur par recopie
Gval_spher::Gval_spher(const Gval_spher& titi)
  :Grille_val(titi) 
{
  if (titi.tet != 0x0) tet = new Tbl(*titi.tet) ;
  if (titi.teti != 0x0) teti = new Tbl(*titi.teti) ;
  if (titi.phi != 0x0) phi = new Tbl(*titi.phi) ;
  if (titi.phii != 0x0) phii = new Tbl(*titi.phii) ;

}
   
// Depuis un fichier
Gval_spher::Gval_spher(FILE* fd)
  : Grille_val(fd)
{
  if (dim.ndim >= 2) {
    tet = new Tbl(fd) ;
    teti = new Tbl(fd) ;
  }
  if (dim.ndim >= 3) {
    phi = new Tbl(fd) ;
    phii = new Tbl(fd) ;
  }
}
  
			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Gval_spher::~Gval_spher() {

  if (tet != 0x0) delete tet ;
  if (teti != 0x0) delete teti ;
  if (phi != 0x0) delete phi ;
  if (phii != 0x0) delete phii ;
  
}

			//-------------//
			// Affectation //
			//-------------//

// Depuis une autre Grille_val
void Gval_spher::operator=(const Gval_spher& titi) {

  Grille_val::operator=(titi) ;

  if (titi.tet != 0x0) *tet = *titi.tet ;
  if (titi.teti != 0x0) *teti = *titi.teti ;
  if (titi.phi != 0x0) *phi = *titi.phi ;
  if (titi.phii != 0x0) *phii = *titi.phii ;
}

			//------------//
			// Sauvegarde //
			//------------//

// save onto a file
void Gval_spher::sauve(FILE* fd) const {

  Grille_val::sauve(fd) ;

  if (dim.ndim >= 2) {
    tet->sauve(fd) ; 
    teti->sauve(fd) ;
  }
  if (dim.ndim >= 3) {
    phi->sauve(fd) ; 
    phii->sauve(fd) ;
  }

}
    
			//------------//
			// Impression //
			//------------//

// Operateurs <<
ostream& Gval_spher::operator>>(ostream& o) const {

  int ndim = dim.ndim;
  Grille_val::operator>>(o) ;

  o << "*** Spherical grid ***" << endl ;

  switch (ndim) {

  case 1 : {
    o << "R nodes: " << iendl ;
    for (int i=0; i<dim.dim[0]; i++) {
      o << zr->set(i+nfantome) << " " ;
    }
    o << endl ;
    break ;
  }
  
  
  case 2 : {
    o << "THETA nodes: " << iendl ;
    for (int j=0 ; j<dim.dim[1] ; j++) {
      o << " " << tet->set(j+nfantome) ;
    }
    o << endl ;
    
    o << "R nodes: " << iendl ;
    
    for (int i=0 ; i<dim.dim[0] ; i++) {
      o << " " << zr->set(i+nfantome) ;
    }
    o << endl ;
    break ;    
  }
  
  case 3 : {
    o << "PHI nodes: " << iendl ;
    for (int k=0 ; k<dim.dim[2] ; k++) {
      o << " " << phi->set(k+nfantome) ;
    }
    o << endl ;

    o << "THETA nodes: " << iendl ;
    for (int j=0 ; j<dim.dim[1] ; j++) {
      o << " " << tet->set(j+nfantome) ;
    }
    o << endl ;
    
    o << "R nodes: " << iendl ;
    
    for (int i=0 ; i<dim.dim[0] ; i++) {
      o << " " << zr->set(i+nfantome) ;
    }
    o << endl ;
    break ;
  }
  
  default : {
    cout << "operator>> Gval_spher : unexpected dimension !" << endl ;
    cout << " ndim = " << ndim << endl ;        
    abort() ;
    break ;
  }
  }
  
  return o ;
}
