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
 * Revision 1.4  2006/06/01 11:47:50  j_novak
 * Memory error hunt.
 *
 * Revision 1.3  2006/06/01 08:18:16  n_vasset
 * Further implementation of Spheroid class definitions
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "metric.h"
#include "spheroid.h"
#include "utilitaires.h"
#include "param.h"
#include "itbl.h"
#include "map.h"
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
    const Map& map3 = Kij.get_mp();
    const Mg3d& gri2d = *map.get_mg() ;

    const Map_radial* mp_rad = dynamic_cast<const Map_radial*>(&Kij.get_mp()) ;
    assert(mp_rad != 0x0) ; 
    
    const Mg3d& gri3d = *map3.get_mg();
    assert(&gri2d == Kij.get_mp().get_mg()->get_angu_1dom()) ;


    int np = gri2d.get_np(0) ;
    int nt = gri2d.get_nt(0) ;
    //  int nr = gri2d.get_nr(0) ;

    int nr3 = gri3d.get_nr(0) ;
    int nt3 = gri3d.get_nt(0) ;    
    int np3 = gri3d.get_np(0) ;
    int nz = gri3d.get_nzone() ;
    assert( nt == nt3 ) ; 
    assert ( np == np3 ); 
    assert ( nr3 != 1 );  
 

    Param pipo ;
    double xi = 0. ;
    int lz = 0 ;
    // Copy of the 2-dimensional h_surf to a 3_d h_surf (calculus commodity, in order to match grids)
    Scalar h_surf3 (map3); 
 
    h_surf3.allocate_all();
    h_surf3.std_spectral_base();
    for (int f= 0; f< nz; f++)
        for (int k=0; k<np3; k++)
	for (int j=0; j<nt3; j++) {
	    for (int l=0; l<nr3; l++) {
 
		h_surf3.set_grid_point(f,k,j,l) = h_surf.val_grid_point(0,k,j,0);
	}
	}

	/* Computation of the jacobian matrix linked to the mapping from the
            spheroid to a coordinate sphere. All quantities will then be calculated
             as from a real coordinate sphere 
	*/

	Tensor jac (h_in.get_mp(), 2,COV, h_in.get_mp().get_bvect_spher());
 
	jac.set_etat_zero(); 
        jac.std_spectral_base();
	jac.set(1,1) = 1. ;
        jac.set(2,2)= 1. ;
        jac.set(3,3) = 1. ; 
        jac.set(1,2) = h_surf.srdsdt() ; 
        jac.set(1,3) = h_surf.srstdsdp() ;

        jac.std_spectral_base() ;

	/* Computation of the jacobian projector, on the 3-D grid*/


	Tensor tempjac (Kij.get_mp(), 2, COV, Kij.get_mp().get_bvect_spher());
	tempjac.std_spectral_base();
        tempjac.allocate_all(); 
	for (int f=0; f<nz; f++)          
       for (int k=0; k<np3; k++)
	for (int j=0; j<nt3; j++) {
	    for (int l=0; l<nr3; l++) {
                for( int m=1; m<4; m++)
		    for( int n=1; n<4; n++){
		      tempjac.set (m,n).set_grid_point(f,k,j,l) =jac(m,n).val_grid_point(0,k,j,0);
		    }
	    }
	}



		Tensor jac3d = tempjac.up(0, gamij);

		/* Computation of the jacobian projector on the 2-D grid*/
		Itbl type (2);
                type.set(0)= CON;
		type.set(1)= COV;
             
		Tensor jacproj (h_in.get_mp(), 2, type, h_in.get_mp().get_bvect_spher()); 
                jacproj.allocate_all();
                jacproj.std_spectral_base();
                jacproj.set_etat_qcq();
    for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    jacproj.set(l,m).set_grid_point(0, k, j, 0) = 
		jac3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
       
	}   

		    // Computation of the trivial projector on the two last components
   
	    Tensor proj3d(Kij.get_mp(),2,  type, Kij.get_mp().get_bvect_spher());
	proj3d.allocate_all();
        proj3d.set_etat_zero() ;
       	proj3d.set(3,3)= 1. ;
        proj3d.set(2,2)= 1. ;    
        
        proj3d.std_spectral_base();
	proj.allocate_all() ; 

         for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    proj.set(l,m).set_grid_point(0, k, j, 0) = 
		proj3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
       
      
	}   


  // Scalar field which annulation characterizes the 2-surface
    Scalar carac = h_surf3-Kij.get_mp().r ; 
    // Computation of the normal vector (covariant form) 
    Vector ss (h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()) ;
    Vector ss3d= carac.derive_cov(gamij) ;

    ss3d = contract( proj3d, 0, contract(jac3d, 0, ss3d, 0), 0); 
    ss.allocate_all() ; 
    ss.std_spectral_base();

        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    ss.set(l).set_grid_point(0, k, j, 0) = 
		ss3d(l).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}
	Sym_tensor sxss = ss*ss ; 


	/* Computation of the metric linked to the 2-surface. 
           It will be designed as a block-diagonal 3-metric, with 1 for 
           the first coordinate and the concerned 2-d  metric as a 
	   second block */

	Sym_tensor qq3d (Kij.get_mp(), COV, Kij.get_mp().get_bvect_spher());
        Sym_tensor qab2 (h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher());
	qq3d.set_etat_zero();
        qq3d.set(1,1) = 1;
        qq3d.set(2,2)= gamij.cov()(1,1) * (h_surf3.srdsdt())* (h_surf3.srdsdt()) + 2*  gamij.cov()(1,2)* (h_surf3.srdsdt()) +  gamij.cov()(2,2);
        qq3d.set(3,3)= gamij.cov()(1,1)* (h_surf3.srstdsdp())*(h_surf3.srstdsdp())+2*gamij.cov()(1,3)* (h_surf3.srstdsdp()) +gamij.cov()(3,3);
        qq3d.set(2,3)= gamij.cov()(1,1)* (h_surf3.srdsdt())* (h_surf3.srstdsdp())+
	    gamij.cov()(1,2)* (h_surf3.srstdsdp())+gamij.cov()(1,3)* (h_surf3.srdsdt()) + gamij.cov()(2,3) ; 
        qq3d.set(3,2)= gamij.cov()(1,1)* (h_surf3.srdsdt())* (h_surf3.srstdsdp())+
	    gamij.cov()(1,2)* (h_surf3.srstdsdp())+gamij.cov()(1,3)* (h_surf3.srdsdt()) + gamij.cov()(2,3) ; 
         qq3d.std_spectral_base();
         qab2.allocate_all() ;
         qab2.std_spectral_base(); 
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    qab2.set(l,m).set_grid_point(0, k, j, 0) = 
		qq3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}
            qab= qab2;

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

    Scalar fff3d (map3);
    fff3d = 1. ; 
    fff.allocate_all() ; 
    fff3d.std_spectral_base();
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    fff.set_grid_point(0, k, j, 0) = 
		fff3d.get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}


    // Computation of the normalization factor of the ingoing null vector.

    Scalar ggg3d (map3);
    ggg3d = 1. ; 
    ggg.allocate_all() ; 
    ggg3d.std_spectral_base();
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    ggg.set_grid_point(0, k, j, 0) = 
		ggg3d.get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}


	
     
    /* Computation of the tangent part of the extrinsic curvature of
     * the 2 surface embedded in the 3 slice. The reference vector
     used is the vector field s */
   
    Tensor tempq= qab.cov() ;
    
  
    // now we need to convert this tensor in the 3d grid in order to 
    //perform following contraction. (new variable temp2q)

    Tensor temp2q(Kij.get_mp(),2, COV, Kij.get_mp().get_bvect_spher());
    temp2q.allocate_all();
    for( int m=1; m<4; m++)
	for( int n=1; n<4; n++){
	    temp2q.set(m,n) = -1.111 ;
	    for (int f=0; f< nz; f++)
		for (int k=0; k<np3; k++)
		    for (int j=0; j<nt3; j++) {
			for (int l=0; l<nr3; l++) {
			    temp2q.set(m,n).set_grid_point(f,k,j,l) =
				tempq(m,n).val_grid_point(0,k,j,0);
			}
		    }
	}
    temp2q.std_spectral_base();
    
    Tensor sxss3 = ss3d * ss3d;
    sxss3.std_spectral_base();
 

	Tensor temp = ss3d.up(0, gamij) * (temp2q.up(0, gamij)) ;
	Vector ll3d = contract(contract(Kij, 1, temp ,1), 0, 1); 

	/* Computation of the 3-d l vector on the "round" coordinates
	 * using jacobian
         * NB: only for intermediate calculus purpose */
      Vector ll3 = contract(proj3d,0, (contract(jac3d,0, ll3d,0)),0);
            ll.allocate_all() ;
            ll.std_spectral_base(); 
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    ll.set(l).set_grid_point(0, k, j, 0) = 
		ll3(l).get_spectral_va().val_point_jk(lz, xi, j, k) ;

	} 
 
	// note: check whether projection method is correct here or not
    
	/* computation of the tangent components of the extrinsic curvature 
         *of the 3-slice 
	 *(extracted from the curvature of the timeslice)
         * Note: this is not the actual 2d_ extrinsic curvature, but the 
         *tangent part of the time-slice extrinsic curvature (nuance)*/
 
    Tensor jj3d= Kij + ss3d*ll3d + ll3d*ss3d + contract((Kij.up(0, gamij)).up(1,gamij), 0,1,sxss3, 0, 1)* sxss3;
    Tensor jj3 =contract(proj3d,0,( contract(proj3d,0,contract(jac3d,0,(contract(jac3d,0, jj3d,1)),0),1)),0);
                     jj.allocate_all() ; 
		     jj.std_spectral_base();
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    jj.set(l,m).set_grid_point(0, k, j, 0) = 
		jj3(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
 
	} 

	// COmputation of 2d extrinsic curvature in the 3-slice
    
	Tensor hh3d = contract(qq3d.up(0, gamij), 0,ss3d.derive_cov(gamij),1) ;

        Tensor hh3 =contract(proj3d,0,( contract(proj3d,0,contract(jac3d,0,(contract(jac3d,0, hh3d,1)),0),1)),0);

                     hh.allocate_all() ; 
                     hh.std_spectral_base();
        for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
        for (int l=1; l<4; l++)
        for (int m=1; m<4; m++) {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0), j, k, pipo,
				   lz, xi) ;
	    hh.set(l,m).set_grid_point(0, k, j, 0) = 
		hh3(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;

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

    set_der_0x0() ;
}

void Spheroid::set_der_0x0() const {
    p_sqrt_q = 0x0 ;
    p_area = 0x0 ;
    p_angu_mom = 0x0 ;
    p_theta_plus = 0x0 ;
    p_theta_minus = 0x0 ;
    p_shear = 0x0 ;
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
   const Map_af& mp_ang = dynamic_cast<const Map_af&>(h_surf.get_mp()) ;
      p_area = new double (mp_ang.integrale_surface((sqrt_q()) * h_surf *h_surf, 1)) ;
  } 
     return *p_area;
} 

// Computation of the angular momentum of the surface (G is set to be identically one)

   double Spheroid::angu_mom(const Vector& phi) const {
    if (p_angu_mom == 0x0) { 
   
	Scalar temp = contract(ll.up(0, qab),0, phi,0);
	p_angu_mom = new double ((sqrt_q()*h_surf*h_surf*temp).integrale()) ;
    }
    return *p_angu_mom;
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

// Ricci tensor of the 2-surface 
 
const Sym_tensor& Spheroid::ricci() const {
    if (p_ricci == 0x0) {
   	Tensor riccia(h_surf.get_mp(), 2, COV, h_surf.get_mp().get_bvect_spher()) ;
        riccia.set_etat_zero();
	for (int i=1; i<=3; i++)
	    riccia.set(1,i) = 0 ;
	Vector ei(h_surf.get_mp(), CON, h_surf.get_mp().get_bvect_spher()) ;
	ei.set_etat_zero() ;
	ei.set(2) = 1 ;
	ei.set(2).std_spectral_base() ;
	ei.set(2).mult_cost() ;
	ei.set(2).mult_sint() ;

	Tensor dei = ei.derive_cov(qab) ;
	dei = contract(proj, 1, contract(proj, 0, dei, 1), 1) ;
	Tensor d2ei = dei.derive_cov(qab) ;
	d2ei = contract(proj, 1, contract(proj, 0, contract(proj, 0, d2ei, 2), 2), 2) ;
	Vector r1 = d2ei.trace(0,2) - d2ei.trace(0,1) ;
	for (int i=1; i<=3; i++) {
	    r1.set(i).div_sint() ;
	    r1.set(i).div_cost() ;
	}
	for (int j=1; j<=3; j++){
	    riccia.set(j,2)= r1(j);
	}
	
	ei.set_etat_zero() ;
	ei.set(3) = 1 ;
	ei.set(3).std_spectral_base() ;
	ei.set(3).mult_sint() ;

	dei = ei.derive_cov(qab) ;
	dei = contract(proj, 1, contract(proj, 0, dei, 1), 1) ;
	d2ei = dei.derive_cov(qab) ;
	d2ei = contract(proj, 1, contract(proj, 0, contract(proj, 0, d2ei, 2), 2), 2) ;
	Vector r2 = d2ei.trace(0,2) - d2ei.trace(0,1) ;
	for (int i=1; i<=3; i++) {
	    r2.set(i).div_sint() ;
	}
	for (int j=1; j<=3; j++) {
	    riccia.set(j,3) = r2(j) ;
	}
	p_ricci = new Sym_tensor(riccia);
}
    return * p_ricci;
}

void Spheroid::sauve(FILE* ) const {

    cout << "c'est pas fait!" << endl ;
    return ; 

}
