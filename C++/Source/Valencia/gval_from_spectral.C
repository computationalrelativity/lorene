/*
 *  Functions for spectral summation to a Valencia-type grid (see grille_val.h)
 *
 */

/*
 *   Copyright (c) 2001 and 2004 Jerome Novak
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
 * Revision 1.15  2016/12/05 16:18:20  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.14  2014/10/13 08:53:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.12  2009/10/28 13:40:23  j_novak
 * General case for the theta symmetry (now should work).
 *
 * Revision 1.11  2009/10/21 13:19:04  j_novak
 * Going back (temporary) to previous version.
 *
 * Revision 1.9  2007/12/21 10:46:29  j_novak
 * In "from_spectral..." functions: better treatment of ETATZERO case.
 *
 * Revision 1.8  2007/11/02 16:49:12  j_novak
 * Suppression of intermediate array for spectral summation.
 *
 * Revision 1.7  2006/10/02 07:41:03  j_novak
 * Corrected an error in the case r=0, when exporting to a cartesian grid.
 *
 * Revision 1.6  2005/06/23 13:44:18  j_novak
 * Removed some old comments.
 *
 * Revision 1.5  2005/06/23 13:40:08  j_novak
 * The tests on the number of dimensions have been changed to handle better the
 * axisymmetric case.
 *
 * Revision 1.4  2005/06/22 09:11:17  lm_lin
 *
 * Grid wedding: convert from the old C++ object "Cmp" to "Scalar".
 *
 * Revision 1.3  2004/12/17 13:35:04  m_forot
 * Add the case T_LEG
 *
 * Revision 1.2  2004/05/07 13:19:24  j_novak
 * Prevention of warnings
 *
 * Revision 1.1  2004/05/07 12:32:13  j_novak
 * New summation from spectral to FD grid. Much faster!
 *
 *
 * $Header$
 *
 */

#include <cmath>

// Lorene headers
#include "grille_val.h"
#include "proto_f77.h"


                 //--------------------------------------
                 // Sommation depuis une grille spectrale
                 //--------------------------------------

namespace Lorene {
void Grille_val::somme_spectrale1(const Scalar& meudon, double* resu, int taille_in) const {

  int taille = dim.dim[0]+2*nfantome ;
  if (taille != taille_in) {
      cout << "Gval_spher::somme_spectral2():\n" ;
      cout << "grid size incompatible with array size... exiting!" << endl ;
      abort() ;
  }
  int nrv = dim.dim[0]+nfantome ;
  const Map& mp = meudon.get_mp() ;
  int l ;
  double xi ;
  for (int i=0; i<nfantome; i++) resu[i] = 0 ;
  for (int i=nfantome; i<nrv; i++) {
    mp.val_lx(zr->t[i],0.,0.,l,xi) ;
    resu[i] = meudon.get_spectral_va().val_point_jk(l, xi, 0, 0) ;
  }
  for (int i=nrv; i<taille; i++) resu[i] = 0 ;
}
 
void Gval_cart::somme_spectrale2(const Scalar& meudon, double* resu, int taille_in) const {
  int nzv = dim.dim[0] + nfantome ;
  int nxv = dim.dim[1] + nfantome ;
  int nzv2 = dim.dim[0] + 2*nfantome ;
  int nxv2 = dim.dim[1] + 2*nfantome ;
  int taille = nxv2*nzv2 ;
  if (taille != taille_in) {
      cout << "Gval_spher::somme_spectral2():\n" ;
      cout << "grid size incompatible with array size... exiting!" << endl ;
      abort() ;
  }
  const Map& mp = meudon.get_mp() ;
  int l ;
  double xi0, rr, theta ;
  double phi = 0 ;
  int inum = 0 ;
  for (int ix=0; ix<nfantome; ix++) {
    for (int iz=0; iz<nzv2; iz++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int ix=nfantome; ix<nxv; ix++) {
    for (int iz=0; iz<nfantome; iz++) {
      resu[inum] = 0. ;
      inum++ ;
    }
    double xx2 = (x->t[ix])*(x->t[ix]) ;
    for (int iz=nfantome; iz<nzv; iz++) {
      rr = sqrt((zr->t[iz])*(zr->t[iz]) + xx2) ;
      theta = (rr != 0. ? acos((zr->t[iz])/rr) : 0) ;
      mp.val_lx(rr, theta, phi, l, xi0) ;
      resu[inum] = meudon.get_spectral_va().val_point(l, xi0, theta, phi) ;
      inum++ ;
    }
    for (int iz=nzv; iz<nzv2; iz++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int ix=nxv; ix<nxv2; ix++) {
    for (int iz=0; iz<nzv2; iz++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }  
}

void Gval_cart::somme_spectrale3(const Scalar& meudon, double* resu, int taille_in) const{
  int nzv = dim.dim[0] + nfantome ;
  int nxv = dim.dim[1] + nfantome ;
  int nyv = dim.dim[2] + nfantome ;
  int nzv2 = dim.dim[0] + 2*nfantome ;
  int nxv2 = dim.dim[1] + 2*nfantome ;
  int nyv2 = dim.dim[2] + 2*nfantome ;
  int taille = nyv2*nxv2*nzv2 ;
  if (taille != taille_in) {
      cout << "Gval_spher::somme_spectral2():\n" ;
      cout << "grid size incompatible with array size... exiting!" << endl ;
      abort() ;
  }
  const Map& mp = meudon.get_mp() ;
  int l ;
  double xi0, rr, theta, phi ;
  int inum = 0 ;
  for (int iy=0; iy<nfantome; iy++) {
    for (int ix=0; ix<nxv2; ix++) {
      for (int iz=0; iz<nzv2; iz++){
	resu[inum] = 0. ;
	inum++ ;
      }
    }
  }
  for (int iy=nfantome; iy<nyv; iy++) { 
    double yy = x->t[iy] ;
    double yy2 = yy*yy ;
    for (int ix=0; ix<nfantome; ix++) {
      for (int iz=0; iz<nzv2; iz++) {
	resu[inum] = 0. ;
	inum++ ;
      }
    }
    for (int ix=nfantome; ix<nxv; ix++) {
      for (int iz=0; iz<nfantome; iz++) {
	resu[inum] = 0. ;
	inum++ ;
      }
      double xx = x->t[ix] ;
      double xx2 = xx*xx ;
      for (int iz=nfantome; iz<nzv; iz++) {
	rr = sqrt((zr->t[iz])*(zr->t[iz]) + xx2 + yy2) ;
	theta = (rr != 0. ? acos((zr->t[iz])/rr) : 0. );
	phi = (rr != 0. ? atan2(yy, xx) : 0. ) ; // return value in [-M_PI,M_PI], should work
	mp.val_lx(rr, theta, phi, l, xi0) ;
	resu[inum] = meudon.get_spectral_va().val_point(l, xi0, theta, phi) ;
	inum++ ;
      }
      for (int iz=nzv; iz<nzv2; iz++) {
	resu[inum] = 0. ;
	inum++ ;
      }
    }
    for (int ix=nxv; ix<nxv2; ix++) {
      for (int iz=0; iz<nzv2; iz++) {
	resu[inum] = 0. ;
	inum++ ;
      }
    }  
  }
  for (int iy=nyv; iy<nyv2; iy++) {
    for (int ix=0; ix<nxv2; ix++) {
      for (int iz=0; iz<nzv2; iz++){
	resu[inum] = 0. ;
	inum++ ;
      }
    }
  }
}

void Gval_spher::somme_spectrale2(const Scalar& meudon, double* resu, int taille_in) const {

    assert (dim.ndim >=2) ;
  int nrv = dim.dim[0] + nfantome ;
  int ntv = dim.dim[1] + nfantome ;
  int nrv2 = dim.dim[0] + 2*nfantome ;
  int ntv2 = dim.dim[1] + 2*nfantome ;
  int taille = ntv2*nrv2 ;
  if (taille != taille_in) {
      cout << "Gval_spher::somme_spectral2():\n" ;
      cout << "grid size incompatible with array size... exiting!" << endl ;
      abort() ;
  }
  const Map& mp = meudon.get_mp() ;
  int l ;
  double xi, rr, theta ;
  double phi0 = 0 ;
  int inum = 0 ;
  for (int it=0; it<nfantome; it++) {
    for (int ir=0; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int it=nfantome; it<ntv; it++) {
    for (int ir=0; ir<nfantome; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
    theta = tet->t[it] ;
    for (int ir=nfantome; ir<nrv; ir++) {
      rr = zr->t[ir] ;
      mp.val_lx(rr, theta, phi0, l, xi) ;
      resu[inum] = meudon.get_spectral_va().val_point(l, xi, theta, phi0) ;
      inum++ ;
    }
    for (int ir=nrv; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int it=ntv; it<ntv2; it++) {
    for (int ir=0; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }  
}

double* Gval_spher::somme_spectrale2ri(const Scalar& meudon) const {
  int nrv = dim.dim[0] + 1 + nfantome ;
  int ntv = dim.dim[1] + nfantome ;
  int nrv2 = dim.dim[0] + 1 + 2*nfantome ;
  int ntv2 = dim.dim[1] + 2*nfantome ;
  int taille = ntv2*nrv2 ;
  const Map& mp = meudon.get_mp() ;
  double* resu = new double[taille] ;
  int l ;
  double xi, rr, theta ;
  double phi0 = 0 ;
  int inum = 0 ;
  for (int it=0; it<nfantome; it++) {
    for (int ir=0; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int it=nfantome; it<ntv; it++) {
    for (int ir=0; ir<nfantome; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
    theta = tet->t[it] ;
    for (int ir=nfantome; ir<nrv; ir++) {
      rr = zri->t[ir] ;
      mp.val_lx(rr, theta, phi0, l, xi) ;
      resu[inum] = meudon.get_spectral_va().val_point(l, xi, theta, phi0) ;
      inum++ ;
    }
    for (int ir=nrv; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int it=ntv; it<ntv2; it++) {
    for (int ir=0; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }  
  return resu ;
}

double* Gval_spher::somme_spectrale2ti(const Scalar& meudon) const {
  int nrv = dim.dim[0] + nfantome ;
  int ntv = dim.dim[1] + 1 + nfantome ;
  int nrv2 = dim.dim[0] + 2*nfantome ;
  int ntv2 = dim.dim[1] + 1 + 2*nfantome ;
  int taille = ntv2*nrv2 ;
  const Map& mp = meudon.get_mp() ;
  double* resu = new double[taille] ;
  int l ;
  double xi, rr, theta ;
  double phi0 = 0 ;
  int inum = 0 ;
  for (int it=0; it<nfantome; it++) {
    for (int ir=0; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int it=nfantome; it<ntv; it++) {
    for (int ir=0; ir<nfantome; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
    theta = teti->t[it] ;
    for (int ir=nfantome; ir<nrv; ir++) {
      rr = zr->t[ir] ;
      mp.val_lx(rr, theta, phi0, l, xi) ;
      resu[inum] = meudon.get_spectral_va().val_point(l, xi, theta, phi0) ;
      inum++ ;
    }
    for (int ir=nrv; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }
  for (int it=ntv; it<ntv2; it++) {
    for (int ir=0; ir<nrv2; ir++) {
      resu[inum] = 0. ;
      inum++ ;
    }
  }  
  return resu ;
}

void Gval_spher::somme_spectrale3(const Scalar& meudon, double* resu, int taille_in) const{

  assert(meudon.get_etat() == ETATQCQ) ;
  meudon.get_spectral_va().coef() ;

  //Sizes of both grids
  //-------------------
  int nrv0 = dim.dim[0] ;
  int ntv0 = dim.dim[1] ;
  int nrv = dim.dim[0] + nfantome ;
  int ntv = dim.dim[1] + nfantome ;
  int npv = dim.dim[2] + nfantome ;
  int nrv2 = dim.dim[0] + 2*nfantome ;
  int ntv2 = dim.dim[1] + 2*nfantome ;
  int npv2 = dim.dim[2] + 2*nfantome ;
  int taille = npv2*ntv2*nrv2 ;
  if (taille != taille_in) {
      cout << "Gval_spher::somme_spectral3():\n" ;
      cout << "grid size incompatible with array size... exiting!" << endl ;
      abort() ;
  }
  const Map& mp = meudon.get_mp() ;
#ifndef NDEBUG
  const Map_af* mpaff = dynamic_cast<const Map_af*>(&mp) ;
  assert(mpaff != 0x0) ;
#endif
  const Mg3d* mg = mp.get_mg() ;
  int ntm = mg->get_nt(0) ;
  int npm = mg->get_np(0) ;
  int nz = mg->get_nzone() ;
#ifndef NDEBUG
  for (int lz=1; lz<nz; lz++) { 
    assert (ntm == mg->get_nt(lz)) ; //Same angular grids in all domains...
    assert (npm == mg->get_np(lz)) ;
  }
#endif

  //Intermediate quantities
  //-----------------------
  double* alpha = new double[nrv0*(npm+2)*ntm] ;
  double* p_coef = alpha ;
  double* chebnri = 0x0 ; //size ~ nrv0 * (npm+2) * nr ...
  int* idom = 0x0 ;
  initialize_spectral_r(mp, meudon.get_spectral_va().get_base(), idom, chebnri) ;
  double* p_func = chebnri ;
  Mtbl_cf& mtbcf = *meudon.get_spectral_va().c_cf ;
  double** coefm = new double*[nz] ;
  for (int lz=0; lz<nz; lz++) {
      assert((mtbcf.t[lz])->get_etat() != ETATNONDEF) ;
      coefm[lz] = (mtbcf.t[lz])->t ;
      if (coefm[lz] == 0x0) {
	  int sizem = mg->get_nr(lz)*ntm*(npm+2) ;
	  coefm[lz] = new double[sizem] ;
	  double* pcf = coefm[lz] ;
	  for (int i=0; i<sizem; i++)
	      pcf[i] = 0. ;
      }
  }
 
  //First partial summation
  //-----------------------
  for (int irv=0; irv<nrv0; irv++) {
    int lz = idom[irv] ;
    double* tbcf = coefm[lz] ;
    int nrm = mg->get_nr(lz) ;
    for (int mpm=0; mpm<npm+2; mpm++) {
      for (int ltm=0; ltm<ntm; ltm++) {
	*p_coef = 0 ;
	for (int irm=0; irm<nrm; irm++) {
	  *p_coef += (*tbcf)*(*p_func) ;
	  tbcf++ ;
	  p_func++ ;
//	  cout << *p_func << ", " << *tbcf << ", " << *p_coef << endl ;
	}
	p_coef++ ;
      }
    }
  }

  for (int lz=0; lz<nz; lz++) {
      if ((mtbcf.t[lz])->t == 0x0) delete [] coefm[lz] ;
  }
  delete [] coefm ;
  delete [] chebnri ;
  delete [] idom ;

  double* beta = new double[ntv0*nrv0*(npm+2)] ;
  p_coef = beta ;
  double* tetlj = 0x0 ;
  initialize_spectral_theta(mp, meudon.get_spectral_va().get_base(), tetlj) ;
  p_func = tetlj ;
  double* p_interm = alpha ;

  //Second partial summation
  //------------------------
  for (int jtv=0; jtv<ntv0; jtv++) {
    for (int irv=0; irv<nrv0; irv++) {
      for (int mpm=0; mpm<npm+2; mpm++) {
	*p_coef = 0 ;	
	for (int ltm=0; ltm<ntm; ltm++) {
	  *p_coef += (*p_interm) * (*p_func) ;
	  p_interm++ ;
	  p_func++ ;
	}
	p_coef++ ;
      } // Loop on m      
      p_func -= (npm+2)*ntm ;
    } //Loop on irv
    p_interm = alpha ;
    p_func += (npm+2)*ntm ;
  } //Loop on jtv

  delete [] alpha ;
  delete [] tetlj ;



  // Final summation
  //----------------
  p_interm = beta ;
  double* expmk = 0x0 ;
  initialize_spectral_phi(mp, meudon.get_spectral_va().get_base(), expmk) ;
  p_func = expmk ;
  p_coef = resu ;
  for (int ip=0; ip<nfantome; ip++) {
    for (int it=0; it<ntv2; it++) {
      for (int ir=0; ir<nrv2; ir++){
	*p_coef = 0. ;
	p_coef++ ;
      }
    }
  }
  for (int ip=nfantome; ip<npv; ip++) { 
    for (int it=0; it<nfantome; it++) {
      for (int ir=0; ir<nrv2; ir++) {
	*p_coef = 0. ;
	p_coef++ ;
      }
    }
    for (int it=nfantome; it<ntv; it++) {
      for (int ir=0; ir<nfantome; ir++) {
	*p_coef = 0. ;
	p_coef++ ;
      }
      for (int ir=nfantome; ir<nrv; ir++) {
	*p_coef = 0. ;
	for (int mpm=0; mpm<npm+2; mpm++) {
	  *p_coef += (*p_interm) * (*p_func) ;
	  p_interm++ ;
	  p_func++ ;
	}
	p_coef++ ;
	p_func -= (npm+2) ;
      }
      for (int ir=nrv; ir<nrv2; ir++) {
	*p_coef = 0. ;
	p_coef++ ;
      }
    }
    for (int it=ntv; it<ntv2; it++) {
      for (int ir=0; ir<nrv2; ir++) {
	*p_coef = 0. ;
	p_coef++ ;
      }
    }
    p_func += npm+2 ; //Next point in phi
    p_interm = beta ;
  }
  for (int ip=npv; ip<npv2; ip++) {
    for (int it=0; it<ntv2; it++) {
      for (int ir=0; ir<nrv2; ir++){
	*p_coef = 0. ;
	p_coef++ ;
      }
    }
  }
  delete [] expmk ;
  delete [] beta ;
}


void Gval_spher::initialize_spectral_r(const Map& mp, const Base_val& base,
				       int*& idom, double*& chebnri) const {
    
    int nrv0 = dim.dim[0] ;
    const Mg3d* mg = mp.get_mg() ;
    int npm = mg->get_np(0) ;
    int ntm = mg->get_nt(0) ;
    
    assert (idom == 0x0) ;
    idom = new int[nrv0] ;
    double* xi = new double[nrv0] ;
    int nrmax = 0 ;

    for (int i=0; i<nrv0; i++) {
	mp.val_lx(zr->t[i+nfantome], 0., 0., idom[i], xi[i]) ;
	nrmax += mg->get_nr(idom[i]) ;
    }
    
    assert (chebnri == 0x0) ;
    chebnri = new double[(npm+2)*ntm*nrmax] ;
    double* p_out = chebnri ;
    for (int irv=0; irv<nrv0; irv++) {
	bool nucleus = (mg->get_type_r(idom[irv]) == RARE) ;
	int nmax = (nucleus ? 2*mg->get_nr(idom[irv]) + 1 
		    : mg->get_nr(idom[irv])) ;
	double* cheb = new double[nmax] ;
	cheb[0] = 1. ;
	cheb[1] = xi[irv] ;
	for (int ir=2; ir<nmax; ir++) {
	    cheb[ir] = 2*xi[irv]*cheb[ir-1] - cheb[ir-2] ;
	}

	int base_r = base.get_base_r(idom[irv]) ;
    
	for (int ip=0; ip<npm+2; ip++) {
	    for (int it=0; it<ntm; it++) {
		int fact = 1 ;
		int par = 0 ;
		if (nucleus) {
		    fact = 2 ;
		    switch (base_r) {
			
			case R_CHEBP : {
			    break ;
			}
			    
			case R_CHEBI : {
			    par = 1 ;
			    break ;
			}
			    
			case R_CHEBPI_P : {
			    par = it % 2 ;
			    break ;
			}
			    
			case R_CHEBPI_I : {
			    par = 1 - (it % 2) ;
			    break ;
			}
			case R_CHEBPIM_P : {
			    par = (ip/2) % 2 ;
			    break ;
			}
			    
			case R_CHEBPIM_I : {
			    par = 1 - ((ip/2) % 2) ;
			    break ;
			}
			
			default : {
			    cout << "Gval_spher::initialize_spectral_r : " << '\n' 
				 << "Unexpected radial base !" << '\n' 
				 << "Base : " << base_r << endl ;
			    abort() ;
			break ;
			}
		}
		}
		for (int ir=0; ir<mg->get_nr(idom[irv]); ir++) {
		    *p_out = cheb[fact*ir+par] ;
		    p_out++ ;
		}
	    
	    } // Loop on it
	} // Loop on ip
	delete [] cheb ;
	
    }// Loop on irv
    
    delete [] xi ;

}

void Gval_spher::initialize_spectral_theta(const Map& mp, const Base_val& base,
				       double*& tetlj) const {
 
  int ntv0 = dim.dim[1] ;
  const Mg3d* mg = mp.get_mg() ;
  int npm = mg->get_np(0) ;
  int ntm = mg->get_nt(0) ;
  int base_t = base.get_base_t(0) ;

  assert (tetlj == 0x0) ;
  tetlj = new double[(npm+2)*ntv0*ntm] ;
  double* p_out = tetlj ;

  for (int jtv=0; jtv<ntv0; jtv++) {
    double teta = tet->t[jtv+nfantome] ;
    for (int mpm=0; mpm < npm+2; mpm++) {
      for (int ltm=0; ltm<ntm; ltm++) {
	switch (base_t)  { //## One should use array of functions...
	case T_COS : {
	  *p_out = cos(ltm*teta) ;
	  break ;
	}
	case T_SIN : {
	  *p_out  = sin(ltm*teta) ;
	  break ;
	}
	case T_COS_P : {
	  *p_out  = cos(2*ltm*teta) ;
	  break ;
	}
	case T_COS_I : {
	  *p_out = cos((2*ltm+1)*teta) ;
	  break ;
	}
	case T_SIN_P : {
	  *p_out  = sin(2*ltm*teta) ;
	  break ;
	}
	case T_SIN_I : {
	  *p_out = sin((2*ltm+1)*teta) ;
	  break ;
	}
	case T_COSSIN_CP : {
	  *p_out = ( ((mpm/2) % 2  == 0) ? cos(2*ltm*teta) 
		     : sin((2*ltm+1)*teta)) ;
	  break ;
	}
	case T_COSSIN_CI : {
	  *p_out = ( ((mpm/2) % 2  == 0) ? cos((2*ltm+1)*teta) 
		     : sin(2*ltm*teta)) ;
	  break ;
	}
	case T_COSSIN_SP : {
	  *p_out = ( ((mpm/2) % 2  == 0) ? sin(2*ltm*teta) 
		     : cos((2*ltm+1)*teta)) ;
	  break ;
	}
	case T_COSSIN_SI : {
	  *p_out = ( ((mpm/2) % 2  == 0) ? sin((2*ltm+1)*teta) 
		     : cos(2*ltm*teta)) ;
	  break ;
	}
	case T_COSSIN_C : {
	  *p_out = ( ((mpm/2) % 2  == 0) ? cos(ltm*teta) 
		     : sin(ltm*teta)) ;
	  break ;
	}
	case T_COSSIN_S : {
	  *p_out = ( ((mpm/2) % 2  == 0) ? sin(ltm*teta) 
		     : cos(ltm*teta)) ;
	  break ;  
	}
	default : {
	  cout << "Gval_spher::initialize_spectral_theta : " << '\n' 
	       << "Unexpected theta base !" << '\n' 
	       << "Base : " << base_t << endl ;
	  abort() ;
	  break ;
	}
	}
	p_out++ ;
      }
      if ( (base_t == T_COS_I) || (base_t == T_SIN_P) || (base_t == T_SIN_I) )
	{
	  p_out-- ;
	  *p_out = 0. ;
	  p_out++ ;
	}
    } //Loop on mpm
  } //Loop on jtv

}


void Gval_spher::initialize_spectral_phi(const Map& mp, const Base_val& base,
				       double*& expmk) const {
 
  int npv0 = dim.dim[2] ;
  const Mg3d* mg = mp.get_mg() ;
  int npm = mg->get_np(0) ;
  int base_p = base.get_base_p(0) ;

  assert (expmk == 0x0) ;
  expmk = new double[(npm+2)*npv0] ;
  double* p_out = expmk ;

  for (int kpv=0; kpv<npv0; kpv++) {
    double fi = phi->t[kpv+nfantome] ;
    for (int mpm=0; mpm < npm+2; mpm++) {
	switch (base_p)  { //## One should use array of functions...
	case P_COSSIN : {
	  int m = mpm / 2 ;
	  *p_out  = ( (mpm%2 == 0) ? cos(m*fi) : sin(m*fi) ) ;
	  break ;
	}
	case P_COSSIN_P : {
	  int m = mpm / 2 ;
	  *p_out  = ( (mpm%2 == 0) ? cos(2*m*fi) : sin(2*m*fi) ) ;
	  break ;
	}
	case P_COSSIN_I : {
	  int m = mpm / 2 ;
	  *p_out  = ( (mpm%2 == 0) ? cos((2*m+1)*fi) : sin((2*m+1)*fi) ) ;
	  break ;
	}
	default : {
	  cout << "Gval_spher::initialize_spectral_phi : " << '\n' 
	       << "Unexpected phi base !" << '\n' 
	       << "Base : " << base_p << endl ;
	  abort() ;
	  break ;
	}
	}
	p_out++ ;
    } //Loop on mpm
  } //Loop on kpv

}
}
