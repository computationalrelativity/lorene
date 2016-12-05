/*
 *  Small utility for determining the surface for 2-fluid stars.
 *
 */

/*
 *   Copyright (c) 2002  Jerome Novak
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

 

#include "cmp.h"

namespace Lorene {
Cmp prolonge_c1(const Cmp& uu, const int nzet) {

  const Map_radial* mpi = dynamic_cast<const Map_radial*>( uu.get_mp() ) ; 
  
  if (mpi == 0x0) {
    cout << 
      "prolonge_c1 : The mapping does not belong to the class Map_radial !"
	 << endl ; 
    abort() ;
  }
  
  const Map_radial& mp = *mpi ; 

  const Mg3d& mg = *(mp.get_mg()) ; 
    
  int nz = mg.get_nzone() ; 
  assert((nzet > 0)&&(nzet<nz)) ; 

  Cmp resu(mp) ;
  resu.allocate_all() ;
  double nbc = uu(0, 0, 0, 0) ;
  double resu_ext =  - 0.2 * nbc ;
  const Coord& rr = mp.r ;
  double rout = (+rr)(nzet,0,0,mg.get_nr(nzet)-1) ;
  double rin = 0 ; double resu1 = 0 ; double dresu1 = 0 ;
  double a = 0 ; double b = 0 ;
  for (int k=0; k<mg.get_np(0); k++) 
    for (int j=0; j<mg.get_nt(0); j++) { 
      double ir = 0 ;
      double nm2 = nbc ;
      double nm1 = nbc ;
      double nb0 = nbc ;
      double np1 = nbc ;
      double rm2 = 0 ;
      double rm1 = 0 ;
      double r0 = 0 ;
      double rp1 = 0 ;
      bool dedans = true ;
      for (int lz=0; lz <= nzet; lz++) {
	for (int i=1; i<mg.get_nr(lz); i++) {
	  ir++ ;
	  rm2 = rm1 ;
	  rm1 = r0 ;
	  r0 = rp1 ;
	  rp1 = (+rr)(lz,k,j,i) ;
	  nm2 = nm1 ;
	  nm1 = nb0 ;
	  nb0 = np1 ;
	  np1 = uu(lz,k,j,i) ;
	  if ((np1<=0.) && dedans) {
	    if (ir<2) {
	      cout << "Problem prolonge_c1!" << endl ;
	      abort() ;
	    }
	    resu1 = nm1 * (r0-rp1)*(rm2-rp1) 
	      / ((r0-rm1)*(rm2-rm1))
	      + nm2 * (r0-rp1)*(rm1-rp1) / ((r0-rm2)*(rm1-rm2))
	      + nb0 * (rm1-rp1)*(rm2-rp1) / ((rm1-r0)*(rm2-r0)) ;
	    resu.set(lz,k,j,i) = resu1 ; 
	    dresu1 = nm1 * ((rp1-r0) + (rp1-rm2)) 
	      / ((r0-rm1)*(rm2-rm1))
	      + nm2 * ((rp1-r0) + (rp1-rm1)) / ((r0-rm2)*(rm1-rm2))
	      + nb0 * ((rp1-rm1) + (rp1-rm2)) / ((rm1-r0)*(rm2-r0)) ;
	    a = (dresu1 - 2*(resu_ext - resu1)/(rout - rp1))/
	      ((rout-rp1)*(rout-rp1)) ;
	    b = 0.5*(-dresu1/(rout-rp1) - a*(rout+rp1)) ;
	    rin = rp1 ;
	    dedans = false ;
	  }
	  else {
	    resu.set(lz,k,j,i) = (dedans ? np1 : 
		  resu1*(rout-rp1)/(rout-rin) + resu_ext*(rin-rp1)/(rin-rout)
				  +(rp1-rin)*(rp1-rout)*(a*rp1+b) );
	  }
	}
	resu.set(lz,k,j,0) = (lz==0 ? nbc : 
			      resu(lz-1, k, j, mg.get_nr(lz-1)-1 ) ) ;
      }
    }
  Cmp resu2(mp) ;
  resu2 = resu_ext ;
  resu2.annule(0,nzet) ;
  resu.annule(nzet+1, nz-1) ;
  resu += resu2 ;
  resu.std_base_scal() ;
  return resu ;
}
}
