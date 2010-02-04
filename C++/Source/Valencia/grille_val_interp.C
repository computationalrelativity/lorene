/*
 * Methods for interpolating with class Grille_val, and its derivative classes.
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


char grille_val_interp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.12  2010/02/04 16:44:35  j_novak
 * Reformulation of the parabolic interpolation, to have better accuracy
 *
 * Revision 1.11  2005/06/23 13:40:08  j_novak
 * The tests on the number of dimensions have been changed to handle better the
 * axisymmetric case.
 *
 * Revision 1.10  2005/06/22 09:11:17  lm_lin
 *
 * Grid wedding: convert from the old C++ object "Cmp" to "Scalar".
 *
 * Revision 1.9  2004/05/07 12:32:13  j_novak
 * New summation from spectral to FD grid. Much faster!
 *
 * Revision 1.8  2004/03/25 14:52:33  j_novak
 * Suppressed some documentation/
 *
 * Revision 1.7  2003/12/19 15:05:14  j_novak
 * Trying to avoid shadowed variables
 *
 * Revision 1.6  2003/12/05 14:51:54  j_novak
 * problem with new SGI compiler
 *
 * Revision 1.5  2003/10/03 16:17:17  j_novak
 * Corrected some const qualifiers
 *
 * Revision 1.4  2002/11/13 11:22:57  j_novak
 * Version "provisoire" de l'interpolation (sommation depuis la grille
 * spectrale) aux interfaces de la grille de Valence.
 *
 * Revision 1.3  2002/09/09 13:00:40  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.2  2001/11/23 16:03:07  j_novak
 *
 *  minor modifications on the grid check.
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
#include "proto_f77.h"

                        //------------------
                        // Compatibilite
                        //------------------

//Compatibilite entre une grille valencienne cartesienne et une meudonaise
bool Gval_cart::compatible(const Map* mp, const int lmax, const int lmin) 
  const {

  //Seulement avec des mappings du genre affine
  assert( dynamic_cast<const Map_af*>(mp) != 0x0) ; 

  const Mg3d* mgrid = mp->get_mg() ;
  assert(lmin >= 0 && lmax <= mgrid->get_nzone()) ;
  int dim_spec = 1 ;
  for (int i=lmin; i<lmax; i++) {
    if ((mgrid->get_nt(i) > 1)&&(dim_spec==1)) dim_spec = 2; 
    if (mgrid->get_np(i) > 1) dim_spec = 3;
  }
  if (dim_spec != dim.ndim) {
    cout << "Grille_val::compatibilite: the number of dimensions" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }
  if (type_t != mgrid->get_type_t()) {
    cout << "Grille_val::compatibilite: the symmetries in theta" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }
  if (type_p != mgrid->get_type_p()) {
    cout << "Grille_val::compatibilite: the symmetries in phi" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }

  bool dimension = true ;
  const Coord& rr = mp->r ;

  double rout = (+rr)(lmax-1, 0, 0, mgrid->get_nr(lmax-1) - 1) ;

  dimension &= (rout <= *zrmax) ;
  switch (dim_spec) {
  case 1:{
    dimension &= ((+rr)(lmin,0,0,0) >= *zrmin) ;
    break ;
  }
  case 2: {
    if (mgrid->get_type_t() == SYM) 
      {dimension &= (*zrmin <= 0.) ;}
    else {
      dimension &= (*zrmin <= -rout ) ;}
    dimension &= (*xmin <= 0.) ;
    dimension &= (*xmax >= rout ) ;
    break ;
  }
  case 3: {
    if (mgrid->get_type_t() == SYM) 
      {dimension &= (*zrmin <= 0.) ;}
    else {
      dimension &= (*zrmin <= -rout) ;}
    if (mgrid->get_type_p() == SYM) {
      dimension &= (*ymin <= 0.) ;
      dimension &= (*xmin <= -rout) ;
    }
    else {
      dimension &= (*xmin <= -rout ) ;
      dimension &= (*ymin <= -rout ) ;
    }
    dimension &= (*xmax >= rout) ;
    dimension &= (*ymax >= rout) ;
    break ;
  }
  }
  return dimension ;

}
//Compatibilite entre une grille valencienne spherique  et une meudonaise
bool Gval_spher::compatible(const Map* mp, const int lmax, const int lmin) 
  const {

  //Seulement avec des mappings du genre affine.
  assert( dynamic_cast<const Map_af*>(mp) != 0x0) ;
 
  int dim_spec = 1 ;

  const Mg3d* mgrid = mp->get_mg() ;
  for (int i=lmin; i<lmax; i++) {
    if ((mgrid->get_nt(i) > 1)&&(dim_spec==1)) dim_spec = 2; 
    if (mgrid->get_np(i) > 1) dim_spec = 3;
  }
  if (dim_spec > dim.ndim) {
    cout << "Grille_val::compatibilite: the number of dimensions" << endl ;
    cout << "of both grids do not coincide!" << endl;
    cout << "Spectral : " << dim_spec << "D,   FD: " << dim.ndim << "D" << endl ;
    abort() ;
  }
  if (type_t != mgrid->get_type_t()) {
    cout << "Grille_val::compatibilite: the symmetries in theta" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }
  if (type_p != mgrid->get_type_p()) {
    cout << "Grille_val::compatibilite: the symmetries in phi" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }

  const Coord& rr = mp->r ;

  int i_b = mgrid->get_nr(lmax-1) - 1 ;

  double rmax = (+rr)(lmax-1, 0, 0, i_b) ;
  double rmin = (+rr)(lmin, 0, 0, 0) ;
  double valmax = get_zr(dim.dim[0]+nfantome - 1) ;
  double valmin = get_zr(-nfantome) ;

  bool dimension = ((rmax <= (valmax)) && (rmin>= (valmin))) ;

  return dimension ;
}

// Teste si la grille valencienne cartesienne est contenue dans le mapping
// de Meudon (pour le passage Meudon->Valence )
bool Gval_cart::contenue_dans(const Map& mp, const int lmax, const int lmin)
 const {
  //Seulement avec des mappings du genre affine
  assert( dynamic_cast<const Map_af*>(&mp) != 0x0) ; 

  const Mg3d* mgrid = mp.get_mg() ;
  assert(lmin >= 0 && lmax <= mgrid->get_nzone()) ;
  int dim_spec = 1 ;
  for (int i=lmin; i<lmax; i++) {
    if ((mgrid->get_nt(i) > 1)&&(dim_spec==1)) dim_spec = 2; 
    if (mgrid->get_np(i) > 1) dim_spec = 3;
  }
  if (dim_spec != dim.ndim) {
    cout << "Grille_val::contenue_dans: the number of dimensions" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }
  if (type_t != mgrid->get_type_t()) {
    cout << "Grille_val::contenue_dans: the symmetries in theta" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }
  if (type_p != mgrid->get_type_p()) {
    cout << "Grille_val::contenue_dans: the symmetries in phi" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }

  bool dimension = true ;
  const Coord& rr = mp.r ;

  //For an affine mapping:
  double radius = (+rr)(lmax-1,0,0,mgrid->get_nr(lmax-1)-1) ;
  double radius2 = radius*radius ;

  if (dim_spec ==1) {
    dimension &= ((+rr)(lmin,0,0,0) <= *zrmin) ;
    dimension &= (radius >= *zrmax) ;
  }
  if (dim_spec ==2) { //## a transformer en switch...
    dimension &= ((+rr)(lmin,0,0,0)/radius < 1.e-6) ;
    dimension &= (*xmin >= 0.) ;
    if (mgrid->get_type_t() == SYM) dimension &= (*zrmin >= 0.) ;
    double x1 = *xmax ;
    double z1 = (fabs(*zrmax)>fabs(*zrmin)? *zrmax : *zrmin) ;
    dimension &= (x1*x1+z1*z1 <= radius2) ;
  }
  if (dim_spec == 3) {
    if (mgrid->get_type_t() == SYM) dimension &= (*zrmin >= 0.) ;
    if (mgrid->get_type_p() == SYM) dimension &= (*ymin >= 0.) ;
    double x1 = (fabs(*xmax)>fabs(*xmin)? *xmax : *xmin) ;
    double y1 = (fabs(*ymax)>fabs(*ymin)? *ymax : *ymin) ;
    double z1 = (fabs(*zrmax)>fabs(*zrmin)? *zrmax : *zrmin) ;
    dimension &= (x1*x1+y1*y1+z1*z1 <= radius2) ;
  }
  return dimension ;
}

// Teste si la grille valencienne spherique est contenue dans le mapping
// de Meudon  (pour le passage Meudon->Valence )
bool Gval_spher::contenue_dans(const Map& mp, const int lmax, const int lmin)
  const {

  //Seulement avec des mappings du genre affine.
  assert( dynamic_cast<const Map_af*>(&mp) != 0x0) ;
 
  const Mg3d* mgrid = mp.get_mg() ;
  
  if (type_t != mgrid->get_type_t()) {
    cout << "Grille_val::contenue_dans: the symmetries in theta" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }
  if (type_p != mgrid->get_type_p()) {
    cout << "Grille_val::contenue_dans: the symmetries in phi" << endl ;
    cout << "of both grids do not coincide!" << endl;
    abort() ;
  }

  const Coord& rr = mp.r ;

  int i_b = mgrid->get_nr(lmax-1) - 1 ;

  double rmax = (+rr)(lmax-1, 0, 0, i_b) ;
  double rmin = (+rr)(lmin, 0, 0, 0) ;
  double valmin = get_zr(0) ;
  double valmax = get_zr(dim.dim[0] - 1) ;

  bool dimension = ((rmax >= valmax) && (rmin<= valmin)) ;

  return dimension ;
}

                        //------------------
                        // Interpolation 1D
                        //------------------

// Interpolation pour la classe de base
Tbl Grille_val::interpol1(const Tbl& rdep, const Tbl& rarr, const Tbl& fdep, 
			  int flag, const int type_inter) const {
  assert(rdep.get_ndim() == 1) ;
  assert(rarr.get_ndim() == 1) ;
  assert(rdep.dim == fdep.dim) ;

  Tbl farr(rarr.dim) ;
  farr.set_etat_qcq() ;
  
  int ndep = rdep.get_dim(0) ;
  int narr = rarr.get_dim(0) ;

  switch (type_inter) {
  case 0: {
    int ndeg[4] ;
    ndeg[0] = ndep ;
    ndeg[1] = narr ;
    double* err0 = new double[ndep+narr] ;
    double* err1 = new double[ndep+narr] ;
    double* den0 = new double[ndep+narr] ;
    double* den1 = new double[ndep+narr] ;
    for (int i=0; i<ndep; i++) {
      err0[i] = rdep(i) ;
      den0[i] = fdep(i) ; 
    }
    for (int i=0; i<narr; i++) err1[i] = rarr(i) ;
    F77_insmts(ndeg, &flag, err0, err1, den0, den1) ;
    for (int i=0; i<narr; i++) farr.set(i) = den1[i] ;
    delete[] err0 ;
    delete[] den0 ;
    delete[] err1 ;
    delete[] den1 ;
    break ;
  }
   case 1: {
     int ip = 0 ;
     int is = 1 ;
     assert(ndep > 1);
     for (int i=0; i<narr; i++) {
       while(rdep(is) < rarr(i)) is++ ;
       assert(is<ndep) ;
       ip = is - 1 ;
       farr.t[i] = (fdep(is)*(rdep(ip)-rarr(i)) + 
 				    fdep(ip)*(rarr(i)-rdep(is))) /
 	     (rdep(ip)-rdep(is)) ;
     }
     break ;
   }
    
   case 2:
     int i1, i2, i3 ;
     double xr, x1, x2, x3, y1, y2, y3 ;
     i2 = 0 ;
     i3 = 1 ;
    assert(ndep > 2) ;
    for (int i=0; i<narr; i++) {
      xr = rarr(i) ;
      while(rdep.t[i3] < xr) i3++ ;
      assert(i3<ndep) ;
      if (i3 == 1) {
	  i1 = 0 ;
	  i2 = 1 ;
	  i3 = 2 ;
      }
      else {
	  i2 = i3 - 1 ;
	  i1 = i2 - 1 ;
      }
      x1 = rdep(i1) ;
      x2 = rdep(i2) ;
      x3 = rdep(i3) ;
      y1 = fdep(i1) ;
      y2 = fdep(i2) ;
      y3 = fdep(i3) ;
      double c = y1 ;
      double b = (y2 - y1) / (x2 - x1) ;
      double a = ( (y3 - y2)/(x3 - x2) - (y2 - y1)/(x2 - x1) ) / (x3 - x1) ;
      farr.t[i] = c + b*(xr - x1) + a*(xr - x1)*(xr - x2) ;
    }
    break ;
    
  case 3:
    cout << "Spline interpolation not implemented yet!" << endl ;
    abort() ;
    break ;
    
   default:
     cout << "Unknown type of interpolation!" << endl ;
     abort() ;
     break ;
   }
   return farr ;
}
  
                        //------------------
                        // Interpolation 2D
                        //------------------

// Interpolation pour les classes derivees
Tbl Gval_spher::interpol2(const Tbl& fdep, const Tbl& rarr, 
		       const Tbl& tarr, const int type_inter) const 
{
  assert(dim.ndim >= 2) ;
  assert(fdep.get_ndim() == 2) ;
  assert(rarr.get_ndim() == 1) ;
  assert(tarr.get_ndim() == 1) ;

  int ntv = tet->get_dim(0) ;
  int nrv = zr->get_dim(0) ;
  int ntm = tarr.get_dim(0) ;
  int nrm = rarr.get_dim(0) ;

  Tbl *fdept = new Tbl(nrv) ;
  fdept->set_etat_qcq() ;
  Tbl intermediaire(ntv, nrm) ;
  intermediaire.set_etat_qcq() ;

  Tbl farr(ntm, nrm) ;
  farr.set_etat_qcq() ;

  int job = 1 ;
  for (int i=0; i<ntv; i++) {
    for (int j=0; j<nrv; j++) fdept->t[j] = fdep.t[i*nrv+j] ;
    Tbl fr(interpol1(*zr, rarr, *fdept, job, type_inter)) ;
    job = 0 ;
    for (int j=0; j<nrm; j++) intermediaire.t[i*nrm+j] = fr.t[j] ;
  }
  delete fdept ;

  fdept = new Tbl(ntv) ;
  fdept->set_etat_qcq() ;
  job = 1 ;
  for (int i=0; i<nrm; i++) {
    for (int j=0; j<ntv; j++) fdept->t[j] = intermediaire.t[j*nrm+i] ;
    Tbl fr(interpol1(*tet, tarr, *fdept, job, type_inter)) ;
    job = 0 ;
    for (int j=0; j<ntm; j++) farr.set(j,i) = fr(j) ;
  }
  delete fdept ;
  return farr ;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

struct Point {
  double x ;
  int l ;
  int k ;
};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */ 

int copar(const void* a, const void* b) {
  double x = (reinterpret_cast<const Point*>(a))->x ;
  double y = (reinterpret_cast<const Point*>(b))->x ;
  return x > y ? 1 : -1 ;
} 

Tbl Gval_cart::interpol2(const Tbl& fdep, const Tbl& rarr, 
		       const Tbl& tetarr, const int type_inter) const 
{
  return interpol2c(*zr, *x, fdep, rarr, tetarr, type_inter) ;
}

Tbl Gval_cart::interpol2c(const Tbl& zdep, const Tbl& xdep, const Tbl& fdep, 
	      const Tbl& rarr, const Tbl& tarr, const int inter_type) const {
  
  assert(fdep.get_ndim() == 2) ;
  assert(zdep.get_ndim() == 1) ;
  assert(xdep.get_ndim() == 1) ;
  assert(rarr.get_ndim() == 1) ;
  assert(tarr.get_ndim() == 1) ;
  
  int nz = zdep.get_dim(0) ;
  int nx = xdep.get_dim(0) ;
  int nr = rarr.get_dim(0) ;
  int nt = tarr.get_dim(0) ;
 
  Tbl farr(nt, nr) ;
  farr.set_etat_qcq() ;
  
  int narr = nt*nr ;
  Point* zlk = new Point[narr] ;
  int inum = 0 ;
  int ir, it ;
  for (it=0; it < nt; it++) {
    for (ir=0; ir < nr; ir++) {
      zlk[inum].x = rarr(ir)*cos(tarr(it)) ; 
      zlk[inum].l = ir ;
      zlk[inum].k = it ;
      inum++ ;
    }
  }
  
  void* base = reinterpret_cast<void*>(zlk) ;
  size_t nel = size_t(narr) ;
  size_t width = sizeof(Point) ;
  qsort (base, nel, width, copar) ; 
  
  Tbl effdep(nz) ; effdep.set_etat_qcq() ;
 
  double x12 = 1e-6*(zdep(nz-1) - zdep(0)) ; 
  // Attention!! x12 doit etre compatible avec son equivalent dans insmts
  int ndistz = 0;
  inum = 0 ;
  do  {
    inum++ ;
    if (inum < narr) {
      if ( (zlk[inum].x - zlk[inum-1].x) > x12 ) {ndistz++ ; }
    }
  } while (inum < narr) ;
  ndistz++ ;
  Tbl errarr(ndistz) ; 
  errarr.set_etat_qcq() ;
  Tbl effarr(ndistz) ; 
  ndistz = 0 ;
  inum = 0 ;
  do  {
    errarr.set(ndistz) = zlk[inum].x ;
    inum ++ ;
    if (inum < narr) {
      if ( (zlk[inum].x - zlk[inum-1].x) > x12 ) {ndistz++ ; }
    }
  } while (inum < narr) ;
  ndistz++ ;

  int ijob = 1 ;

  Tbl tablo(nx, ndistz) ;
  tablo.set_etat_qcq() ;
  for (int j=0; j<nx; j++) {
      for (int i=0; i<nz; i++) effdep.set(i) = fdep(j,i) ;
      effarr = interpol1(zdep, errarr, effdep, ijob, inter_type) ;
      ijob = 0 ;
      for (int i=0; i<ndistz; i++) tablo.set(j,i) = effarr(i) ;
  }

  inum = 0 ;
  int indz = 0 ;
  Tbl effdep2(nx) ;
  effdep2.set_etat_qcq() ;
  while (inum < narr) {
    Point* xlk = new Point[3*nr] ;
    int nxline = 0 ;
    int inum1 ;
    do { 
      ir = zlk[inum].l ;
      it = zlk[inum].k ;
      xlk[nxline].x = rarr(ir)*sin(tarr(it)) ;
      xlk[nxline].l = ir ;
      xlk[nxline].k = it ;
      nxline ++ ; inum ++ ;
      inum1 = (inum < narr ? inum : 0) ;
    } while ( ( (zlk[inum1].x - zlk[inum-1].x) < x12 ) && (inum < narr)) ;
    void* bas2 = reinterpret_cast<void*>(xlk) ;
    size_t ne2 = size_t(nxline) ;
    qsort (bas2, ne2, width, copar) ;
    
    int inum2 = 0 ;
    int ndistx = 0 ;
    do  {
      inum2 ++ ;
      if (inum2 < nxline) {
	if ( (xlk[inum2].x - xlk[inum2-1].x) > x12 ) {ndistx++ ; }
      }
    } while (inum2 < nxline) ;
    ndistx++ ;

    Tbl errarr2(ndistx) ;
    errarr2.set_etat_qcq() ;
    Tbl effarr2(ndistx) ;
    inum2 = 0 ; 
    ndistx = 0 ;
    do  {
      errarr2.set(ndistx) = xlk[inum2].x ;
      inum2 ++ ;
      if (inum2 < nxline) {
	if ( (xlk[inum2].x - xlk[inum2-1].x) > x12 ) {ndistx++ ; }
      }
    } while (inum2 < nxline) ;
    ndistx++ ;

    for (int j=0; j<nx; j++) {
      effdep2.set(j) = tablo(j,indz) ;
    } 
    indz++ ;
    ijob = 1 ;
    effarr2 = interpol1(xdep, errarr2, effdep2, ijob, inter_type) ;
    int iresu = 0 ;
    if (ijob == -1) {
      for (int i=0; i<nxline; i++) {
	while (fabs(xlk[i].x - xdep(iresu)) > x12 ) {
	  iresu++ ;
	}
	ir = xlk[i].l ;
	it = xlk[i].k ;
	farr.set(it,ir) = effdep2(iresu) ;
      }
    }
    else {
      double resu ;
      for (int i=0; i<nxline; i++) {
	resu = effarr2(iresu) ;
	if (i<nxline-1) {
	  if ((xlk[i+1].x-xlk[i].x) > x12) {
	    iresu++ ;
	  }
	}
	ir = xlk[i].l ;
	it = xlk[i].k ;
	farr.set(it,ir) = resu ;
      }
    }
    delete [] xlk ;
  }

  delete [] zlk ;
  return farr ;
}


                        //------------------
                        // Interpolation 3D
                        //------------------

// Interpolation pour les classes derivees
Tbl Gval_spher::interpol3(const Tbl& fdep, const Tbl& rarr, const Tbl& tarr, 
			  const Tbl& parr, const int type_inter) const {
  assert(dim.ndim == 3) ;
  assert(fdep.get_ndim() == 3) ;
  assert(rarr.get_ndim() == 1) ;
  assert(tarr.get_ndim() == 1) ;
  assert(parr.get_ndim() == 1) ;

  int npv = phi->get_dim(0) ;
  int ntv = tet->get_dim(0) ;
  int nrv = zr->get_dim(0) ;
  int npm = parr.get_dim(0) ;
  int ntm = tarr.get_dim(0) ;
  int nrm = rarr.get_dim(0) ;

  Tbl *fdept = new Tbl(ntv, nrv) ;
  fdept->set_etat_qcq() ;
  Tbl intermediaire(npv, ntm, nrm) ;
  intermediaire.set_etat_qcq() ;

  Tbl farr(npm, ntm, nrm) ;
  farr.set_etat_qcq() ;

  for (int i=0; i<npv; i++) {
    for (int j=0; j<ntv; j++) 
      for (int k=0; k<nrv; k++) fdept->t[j*nrv+k] = fdep.t[(i*ntv+j)*nrv+k] ;
    Tbl fr(interpol2(*fdept, rarr, tarr, type_inter)) ;
    for (int j=0; j<ntm; j++)
      for (int k=0; k<nrm; k++) intermediaire.set(i,j,k) = fr(j,k) ;
  }
  delete fdept ;

  int job = 1 ;
  fdept = new Tbl(npv) ;
  fdept->set_etat_qcq() ;
  for (int i=0; i<ntm; i++) {
    for (int j=0; j<nrm; j++) {
      for (int k=0; k<npv; k++) fdept->set(k) = intermediaire(k,i,j) ;
      Tbl fr(interpol1(*phi, parr, *fdept, job, type_inter)) ;
      job = 0 ;
      for (int k=0; k<npm; k++) farr.set(k,i,j) = fr(k) ;
    }
  }
  delete fdept ;
  return farr ;
}

Tbl Gval_cart::interpol3(const Tbl& fdep, const Tbl& rarr, 
			  const Tbl& tarr, const Tbl& parr, const 
			  int inter_type) const {

  assert(fdep.get_ndim() == 3) ;
  assert(rarr.get_ndim() == 1) ;
  assert(tarr.get_ndim() == 1) ;
  assert(parr.get_ndim() == 1) ;
  
  int nz = zr->get_dim(0) ;
  int nx = x->get_dim(0) ;
  int ny = y->get_dim(0) ;
  int nr = rarr.get_dim(0) ;
  int nt = tarr.get_dim(0) ;
  int np = parr.get_dim(0) ;
  Tbl farr(np, nt, nr) ;
  farr.set_etat_qcq() ;

  bool coq = (rarr(0)/rarr(nr-1) > 1.e-6) ;
  Tbl* rarr2(0x0);
  if (coq) {     // If the spectral grid is only made of shells
    rarr2 = new Tbl(2*nr) ;
    rarr2->set_etat_qcq() ;
    double dr = rarr(0)/nr ;
    for (int i=0; i<nr; i++) rarr2->set(i) = i*dr ;
    for (int i=nr; i<2*nr; i++) rarr2->set(i) = rarr(i-nr) ;
  }

  int nr2 = coq ? 2*nr : nr ;

  Tbl cylindre(nz, np, nr2) ;
  cylindre.set_etat_qcq() ;
  for(int iz=0; iz<nz; iz++) {
    Tbl carre(ny,nx) ;
    carre.set_etat_qcq() ;
    Tbl cercle(np, nr2) ;
    for (int iy=0; iy<ny; iy++) 
      for (int ix=0; ix<nx; ix++) 
	carre.set(iy,ix) = fdep(iy,ix,iz) ; // This should be optimized...
    cercle = interpol2c(*x, *y, carre, coq ? *rarr2 : rarr, parr, inter_type) ;
      
    for (int ip=0; ip<np; ip++) 
      for (int ir=0; ir<nr2; ir++) 
	cylindre.set(iz,ip,ir) = cercle(ip,ir) ;
  }

 for (int ip=0; ip<np; ip++) {
    Tbl carre(nr2, nz) ;
    carre.set_etat_qcq() ;
    Tbl cercle(nt, nr) ;
    for (int ir=0; ir<nr2; ir++) 
      for (int iz=0; iz<nz; iz++) 
	carre.set(ir,iz) = cylindre(iz,ip,ir) ;
    cercle = interpol2c(*zr, coq ? *rarr2 : rarr , carre, rarr, tarr, 
			inter_type) ;
    for (int it=0; it<nt; it++) 
      for (int ir=0; ir<nr; ir++) 
	farr.set(ip,it,ir) = cercle(it,ir) ;
  }

 if (coq) delete rarr2 ;
 return farr ;

}


