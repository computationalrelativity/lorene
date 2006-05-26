/*
 *  Definition of methods for the class Spheroid and its subclass App_hor
 *
 */

/*
 *   Copyright (c) 2006  Nicolas Vasset, Jerome Novak & Jose-Luis Jaramillo
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

char spheroid_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2006/05/26 13:49:10  j_novak
 * Method for computing the Ricci tensor on the Spheroid.
 *
 * Revision 1.1  2006/05/26 13:20:43  j_novak
 * New class spheroid.
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "spheroid.h"
#include "utilitaires.h"
#include "param.h"

                     //---------------//
                     //  Constructors //
                     //--------------// 
 
Spheroid::Spheroid(const Map_af& map, double radius):
    h_surf(map), 
    proj (map, COV, map.get_bvect_spher()),
    qab(map.flat_met_spher()), 
    hh(map, COV, map.get_bvect_spher()),
    trk(map),
    ll(map, COV, map.get_bvect_spher()), 
    jj(map, COV, map.get_bvect_spher()),    
    fff(map),
    ggg(map) 

{  
    assert(radius > 1.e-15) ;
    assert(map.get_mg()->get_nzone() == 1) ; // one domain
    assert(map.get_mg()->get_nr(0) == 1) ; // angular grid
    assert(map.get_mg()->get_type_r(0) == FIN) ; //considered as a shell

    h_surf = radius ;
    proj.set_etat_zero();
    hh.set_etat_zero() ;
    for (int i=1; i<=3; i++)
	hh.set(i,i) = 2./radius ;
    trk.set_etat_zero() ;
    ll.set_etat_zero() ;
    jj.set_etat_zero() ;
    fff.set_etat_zero();
    ggg.set_etat_zero();
    set_der_0x0() ;
}

Spheroid::Spheroid(const Scalar& h_in, const Metric& gamij, const Sym_tensor& Kij):
    h_surf(h_in),
    proj(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()),
    qab(h_in.get_mp().flat_met_spher()),
     hh(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()),
    trk(h_in.get_mp()), 
    ll(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()), 
    jj(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()),
    fff(h_in.get_mp()),
    ggg(h_in.get_mp())

{
    const Map& map = h_in.get_mp() ;
    const Mg3d& gri2d = *map.get_mg() ;

    const Map_radial* mp_rad = dynamic_cast<const Map_radial*>(&Kij.get_mp()) ;
    assert(mp_rad != 0x0) ;
    
    assert(&gri2d == Kij.get_mp().get_mg()->get_angu_1dom()) ;

    int np = gri2d.get_np(0) ;
    int nt = gri2d.get_nt(0) ;
    Param pipo ;
    int lz = 0 ;
    double xi = 0. ; 
  
  // Scalar field which annulation characterizes the 2-surface
    Scalar carac = h_surf-map.r ; 
    // Computation of the normal vector (covariant form) 
    Vector ss (h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()) ;
    Vector ss3d= carac.derive_con(gamij) ;
    ss.allocate_all() ; 
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    Scalar temp= ss(l);
	    temp.set_grid_point(0, k, j, 0) = 
		ss3d(l).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	    ss.set(l)=temp;
	}
  
	Sym_tensor sxss = ss*ss ; 
	/* Computation of the metric linked to the 2-surface. 
           It will be designed as a block-diagonal 3-metric, with 1 for 
           the first coordinate and the concerned 2-d  metric as a 
	   second block */

	Sym_tensor qq3d (h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher());
        Sym_tensor qab2 (h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher());
	qq3d.set_etat_zero();
  
        qq3d.set(1,1) = 1;
        qq3d.set(2,2)= gamij.cov()(1,1) * (h_surf.srdsdt())* (h_surf.srdsdt()) + 2*  gamij.cov()(1,2)* (h_surf.srdsdt()) +  gamij.cov()(2,2);
        qq3d.set(3,3)= gamij.cov()(1,1)* (h_surf.srstdsdp())*(h_surf.srstdsdp())+2*gamij.cov()(1,3)* (h_surf.srstdsdp()) +gamij.cov()(3,3);
        qq3d.set(2,3)= gamij.cov()(1,1)* (h_surf.srdsdt())* (h_surf.srstdsdp())+
	    gamij.cov()(1,2)* (h_surf.srstdsdp())+gamij.cov()(1,3)* (h_surf.srdsdt()) + gamij.cov()(2,3) ; 
        qq3d.set(3,2)= gamij.cov()(1,1)* (h_surf.srdsdt())* (h_surf.srstdsdp())+
	    gamij.cov()(1,2)* (h_surf.srstdsdp())+gamij.cov()(1,3)* (h_surf.srdsdt()) + gamij.cov()(2,3) ; 
  
         qab2.allocate_all() ; 
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    Scalar temp3= qab2(l,m);
	    temp3.set_grid_point(0, k, j, 0) = 
		qq3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	    qab2.set(l,m)=temp3;
	
            qab= qq3d;
	}
	    // Computation of the trivial projector on the two last components
   
	    Sym_tensor proj3d(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher());
 
        proj3d.set_etat_zero() ;
       	proj3d.set(3,3)= 1. ;
        proj3d.set(2,2)= 1. ;    
        
	proj.allocate_all() ; 
         for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    Scalar temp4= proj(l,m);
	    temp4.set_grid_point(0, k, j, 0) = 
		proj3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	    proj.set(l,m)=temp4;
       
	}   


    // Computation of the trace of the extrinsic curvature of 3-slice
    Scalar trk3d = Kij.trace(gamij) ;
    trk.allocate_all() ; 
    for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    trk.set_grid_point(0, k, j, 0) = 
		trk3d.get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}
   
    // Computation of the normalization factor of the outgoing null vector.

    Scalar fff3d (map);
    fff3d = 1. ; 
    fff.allocate_all() ; 
  
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    fff.set_grid_point(0, k, j, 0) = 
		fff3d.get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}


    // Computation of the normalization factor of the ingoing null vector.

    Scalar ggg3d (map);
    ggg3d = 1. ; 
    ggg.allocate_all() ; 
  
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    ggg.set_grid_point(0, k, j, 0) = 
		ggg3d.get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}

	/* Computation of the jacobian matrix linked to the mapping from the
            spheroid to a coordinate sphere. All quantities will then be calculated
             as from a real coordinate sphere 
	*/

	Tensor jac (h_in.get_mp(), 2,COV, h_in.get_mp().get_bvect_spher());
 
	jac.set_etat_zero(); 

	jac.set(1,1) = 1. ;
        jac.set(2,2)= 1. ;
        jac.set(3,3) = 1. ; 
        jac.set(1,2) = h_surf.srdsdt() ; 
        jac.set(1,3) = h_surf.srstdsdp() ;

        jac.std_spectral_base() ;

    /* Computation of the tangent part of the extrinsic curvature of
     * the 2 surface embedded in the 3 slice. The reference vector
     used is the vector field s */
   
    Tensor temp = ss * qab.cov() ; 
    Vector ll3d = contract( Kij, 1, temp, 1);
    Vector ll2 (h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher());
    ll2.set_etat_zero();
            ll.allocate_all() ; 
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    Scalar temp1= ll2(l);
	    temp1.set_grid_point(0, k, j, 0) = 
		ll3d(l).get_spectral_va().val_point_jk(lz, xi, j, k) ;
            ll2.set(l)=temp1;

	    ll = (contract(jac.up(0, gamij),0, ll,0));
	} 
    
	/* computation of the tangent components of the extrinsic curvature 
         *of the 3-slice 
	 *(extracted from the curvature of the timeslice)
         * Note: this is not the actual 2d_ extrinsic curvature, but the 
         *tangent part of the time-slice extrinsic curvature (nuance)*/
 
    Tensor jj3d= Kij + ss*ll + ll*ss + contract(Kij, 0,1,sxss, 0, 1);
    Tensor jj2  (h_in.get_mp(), 2,COV, h_in.get_mp().get_bvect_spher());
    jj2.set_etat_zero();
                     jj.allocate_all() ; 
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    Scalar temp2= jj2(l,m);
	    temp2.set_grid_point(0, k, j, 0) = 
		jj3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
            jj2.set(l,m)=temp2;
 
	    jj = contract(jac.up(0, gamij),0,(contract(jac.up(0, gamij),0, jj,1)),0);
	} 

	// COmputation of 2d extrinsic curvature in the 3-slice
    
	Tensor hh3d = qab.cov().up(1, gamij)*ss.derive_con(gamij) ;
        Tensor hh2  (h_in.get_mp(), 2,COV, h_in.get_mp().get_bvect_spher());
	hh2.set_etat_zero();
                     hh.allocate_all() ; 
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    Scalar temp3= hh2(l,m);
	    temp3.set_grid_point(0, k, j, 0) = 
		hh3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	    hh2.set(l,m)=temp3;

	    hh= contract(jac.up(0, gamij),0,(contract(jac.up(0, gamij),0, hh,1)),0); 
	}
       

          
    set_der_0x0() ;
}





//Copy constructor//

Spheroid::Spheroid (const Spheroid &sph_in) :h_surf(sph_in.h_surf),
                                             proj(sph_in.proj),
					    qab(sph_in.qab),
					    hh(sph_in.hh),
					    trk(sph_in.trk),
					    ll(sph_in.ll),
					    jj(sph_in.jj),
					     fff(sph_in.fff),
                                             ggg(sph_in.ggg)
              
{
    set_der_0x0() ; 

}
                    //------------//
                    //Destructor //
                    //-----------//

Spheroid::~Spheroid()
{
    del_deriv() ;
}

                    // -----------------//
                    // Memory management//
                    //------------------//
void Spheroid::del_deriv() const {
    if (p_sqrt_q != 0x0) delete p_sqrt_q ;
    if (p_area != 0x0) delete p_area ;
    if (p_angu_mom != 0x0) delete p_angu_mom ;
    if (p_theta_plus != 0x0) delete p_theta_plus ;
    if (p_theta_minus != 0x0) delete p_theta_minus ;
    if (p_shear != 0x0) delete p_shear ;
    if (p_ricci != 0x0) delete p_ricci ;

    set_der_0x0() ;
}

void Spheroid::set_der_0x0() const {
    p_sqrt_q = 0x0 ;
    p_area = 0x0 ;
    p_angu_mom = 0x0 ;
    p_theta_plus = 0x0 ;
    p_theta_minus = 0x0 ;
    p_shear = 0x0 ;
    p_ricci = 0x0 ;
} 


 
                    //---------//
                    //Accessors//
                    //---------//

// Computation of the 2-dimensional Jacobian amplitude for the surface
  const  Scalar& Spheroid::sqrt_q() const { 
      if (p_sqrt_q == 0x0) {
        p_sqrt_q = new Scalar(sqrt(qab.determinant())) ;
	p_sqrt_q->std_spectral_base() ;
      }
      return *p_sqrt_q ; 
  }

// Computation of the 2-dimensional area of the surface

   double  Spheroid::area() const {
      if (p_area == 0x0) {
   
      p_area = new double ((sqrt_q()*h_surf*h_surf).integrale()) ;
  } 
     return *p_area;
} 

// Outgoing null expansion of 2-surface

const Scalar &Spheroid::theta_plus() const {

    if (p_theta_plus == 0x0) {
        p_theta_plus = new Scalar(fff*(hh.trace() - jj.trace())) ;
	p_theta_plus->std_spectral_base() ;
      }

    return *p_theta_plus; 
				}
// ingoing null expansion of 2-surface

const Scalar& Spheroid::theta_minus() const {

    if (p_theta_minus == 0x0) {
        p_theta_minus = new Scalar(ggg*(-hh.trace() - jj.trace())) ;
	p_theta_minus->std_spectral_base() ;
      }

    return *p_theta_minus; 
 
}

// null-oriented shear of 2-surface

const Sym_tensor& Spheroid::shear() const { 

    if (p_shear == 0x0) {
        p_shear = new Sym_tensor(fff*(jj + hh - (qab.cov()/2) *(hh.trace() - jj.trace()))) ;
	p_shear->std_spectral_base() ;
      }

    return *p_shear; 
 
}

 
void Spheroid::sauve(FILE* ) const {

    cout << "c'est pas fait!" << endl ;
    return ; 

}
