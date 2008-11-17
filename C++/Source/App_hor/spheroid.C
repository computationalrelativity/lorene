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
 * Revision 1.14  2008/11/17 08:30:01  jl_jaramillo
 * Nicolas's changes on multipoles. Sign correction in outgoing shear expression
 *
 * Revision 1.13  2008/11/12 15:17:47  n_vasset
 * Bug-hunting, and new definition for computation of the Ricci scalar
 * (instead of Ricci tensor previously)
 *
 * Revision 1.12  2008/07/09 08:47:33  n_vasset
 * new version for multipole calculation. Function zeta implemented.
 * Revision 1.11  2008/06/04 12:31:23  n_vasset
 * New functions multipole_mass and multipole_angu. first version.
 *
 * Revision 1.9  2006/09/07 08:39:45  j_novak
 * Minor changes.
 *
 * Revision 1.8  2006/07/03 10:13:48  n_vasset
 *  More efficient method for calculation of ricci tensor. Adding of flag issphere
 *
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
#include <assert.h>
#include "nbr_spx.h"
#include "math.h"
#include "param_elliptic.h"
#include "tensor.h"
#include "sym_tensor.h"
#include "diff.h"
#include "proto.h"
#include "param.h"


 


//---------------//
//  Constructors //
//--------------// 
 
Spheroid::Spheroid(const Map_af& map, double radius):
  h_surf(map), 
  jac2d(map, 2, COV, map.get_bvect_spher()),
  proj(map, 2, COV, map.get_bvect_spher()), 
  qq(map, COV, map.get_bvect_spher()),
  ss ( map, COV, map.get_bvect_spher ()),
  qab(map.flat_met_spher()), 
  ricci(map),
  hh(map, COV, map.get_bvect_spher()),
  trk(map),
  ll(map, COV, map.get_bvect_spher()), 
  jj(map, COV, map.get_bvect_spher()),    
  fff(map),
  ggg(map), 
  zeta(map),
  issphere(true)
{    

  //    Itbl type (2) ; 
  //     type.set(1) = CON ; 
  //     type.set(2) = COV ;
  //    Tensor proj(map, 2, type, map.get_bvect_spher());


  assert(radius > 1.e-15) ;
  assert(map.get_mg()->get_nzone() == 1) ; // one domain
  assert(map.get_mg()->get_nr(0) == 1) ; // angular grid
  assert(map.get_mg()->get_type_r(0) == FIN) ; //considered as a shell

  // Setting of real index types for jacobian and projector (first contravariant, second covariant)
  jac2d.set_index_type(0) = CON ;
  proj.set_index_type(0) = CON ; 
 
  jac2d.set_etat_zero() ; 

 
  h_surf = radius ;
  proj.set_etat_zero();
  hh.set_etat_zero() ;
  for (int i=1; i<=3; i++)
    hh.set(i,i) = 2./radius ; // WTF?
  trk.set_etat_zero() ;
  ll.set_etat_zero() ;
  jj.set_etat_zero() ;
  fff.set_etat_zero();
  ggg.set_etat_zero();
  zeta.set_etat_zero();
  set_der_0x0() ;
}

Spheroid::Spheroid(const Scalar& h_in, const Metric& gamij, const Sym_tensor& Kij):
  h_surf(h_in),
  jac2d(h_in.get_mp(),2, COV, h_in.get_mp().get_bvect_spher()),
  proj(h_in.get_mp(),2, COV, h_in.get_mp().get_bvect_spher()),
  qq(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()),
  ss (h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()),
  qab(h_in.get_mp().flat_met_spher()),
  ricci(h_in.get_mp()),
  hh(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()),
  trk(h_in.get_mp()), 
  ll(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()), 
  jj(h_in.get_mp(), COV, h_in.get_mp().get_bvect_spher()),
  fff(h_in.get_mp()),
  ggg(h_in.get_mp()),
  zeta(h_in.get_mp()),
  issphere(true)

{
  const Map& map = h_in.get_mp() ; // The 2-d 1-domain map
  const Map& map3 = Kij.get_mp(); // The 3-d map
  const Mg3d& gri2d = *map.get_mg() ;

  const Map_radial* mp_rad = dynamic_cast<const Map_radial*>(&Kij.get_mp()) ;
  assert(mp_rad != 0x0) ; 
    
  const Mg3d& gri3d = *map3.get_mg();
  assert(&gri2d == Kij.get_mp().get_mg()->get_angu_1dom()) ;

 

  int np = gri2d.get_np(0) ;
  int nt = gri2d.get_nt(0) ;

  int nr3 = gri3d.get_nr(0) ;
  int nt3 = gri3d.get_nt(0) ;    
  int np3 = gri3d.get_np(0) ;
  int nz = gri3d.get_nzone() ;
  assert( nt == nt3 ) ; 
  assert ( np == np3 ); 
  assert ( nr3 != 1 );  
 

  Param pipo ;
  double xi = 0. ;
  int lz = 1 ;


  // Setting of real index types forjacobian and projector(first contravariant, other covariant)
  proj.set_index_type(0) = CON; 
  jac2d.set_index_type(0) = CON; 

  // Copy of the 2-dimensional h_surf to a 3_d h_surf (calculus commodity, in order to match grids)
  Scalar h_surf3 (map3); 
 
  h_surf3.allocate_all();
  h_surf3.std_spectral_base();
  for (int f= 0; f<nz; f++)
    for (int k=0; k<np3; k++)
      for (int j=0; j<nt3; j++) {
	for (int l=0; l<nr3; l++) {
		
	  h_surf3.set_grid_point(f,k,j,l) = h_surf.val_grid_point(0,k,j,0);
	
	}
      }

   h_surf3.annule_domain(0);
   h_surf3.annule_domain(nz - 1);

   //  h_surf3.std_spectral_base();

 


  /* Computation of the jacobian projector linked to the mapping from the
     spheroid to a coordinate sphere. All quantities will then be calculated
     as from a real coordinate sphere 
  */

  Itbl type (2); 
  type.set(0) = CON ; 
  type.set(1) = COV ;


  Tensor jac (Kij.get_mp(),2,type, Kij.get_mp().get_bvect_spher());
    
  jac.set_etat_zero(); 
  jac.std_spectral_base();
  jac.set(1,1) = 1. ;
  jac.set(2,2)= 1. ;
  jac.set(3,3) = 1. ; 
  jac.set(1,2) = -h_surf3.srdsdt() ; 
  jac.set(1,3) = -h_surf3.srstdsdp() ;
 
  jac.std_spectral_base() ; 

  // Copy on the 2-d grid

  jac2d.allocate_all() ; 
  jac2d.std_spectral_base();
  for (int l=1; l<4; l++)
    for (int m=1; m<4; m++)     
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  {
	    mp_rad->val_lx_jk((h_surf.val_grid_point(0, k, j, 0))*1.00000000001, j, k, pipo,
			      lz, xi) ;
	    jac2d.set(l,m).set_grid_point(0, k, j, 0) = 
	      jac(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;

	  }




  // Inverse jacobian (on 3-d grid)
  Tensor jac_inv = jac ;
  jac_inv.set(1,2) = - jac_inv(1,2);  
  jac_inv.set(1,3) = - jac_inv(1,3) ;
    
 


  // Scalar field which annulation characterizes the 2-surface
  Scalar carac = Kij.get_mp().r - h_surf3; 
  carac.std_spectral_base();
 // Computation of the normal vector (covariant form) on both grids

  Vector ss3d= carac.derive_cov(gamij) ;
  Vector ss3dcon=  carac.derive_con(gamij) ;
  Scalar ssnorm =  contract (ss3d.up(0, gamij), 0, ss3d, 0); 
  ssnorm.std_spectral_base() ;
  ss3d =  ss3d / sqrt (ssnorm) ; 
  ss3dcon =  ss3dcon / sqrt (ssnorm) ; 
  ss3d.annule_domain(0);
  ss3dcon.annule_domain(0);
  ss3d.annule_domain(nz-1);
  ss3dcon.annule_domain(nz -1);
  ss3d.std_spectral_base();
  ss3dcon.std_spectral_base();


  // Provisory handling of dzpuis problems 
  h_surf3.annule_domain(nz-1);
  for (int ii=1; ii <=3; ii++){
    ss3d.set(ii).dec_dzpuis(ss3d(ii).get_dzpuis());
    ss3dcon.set(ii).dec_dzpuis(ss3dcon(ii).get_dzpuis());
  }

  for (int ii=1; ii <=3; ii++)
    for (int jjy = 1; jjy <=3; jjy ++){
      jac_inv.set(ii, jjy).dec_dzpuis(jac_inv(ii, jjy).get_dzpuis());
      jac.set(ii, jjy).dec_dzpuis(jac(ii, jjy).get_dzpuis());
    }

  // End of dzpuis handling.


  Sym_tensor sxss3d = ss3d * ss3d ;
  Sym_tensor sxss3dcon = ss3dcon * ss3dcon ; 
  Vector ss3 (Kij.get_mp(), COV, Kij.get_mp().get_bvect_spher()) ;
  Vector ss3con(Kij.get_mp(), CON, Kij.get_mp().get_bvect_spher()) ;


  ss3 = contract(jac_inv, 0, ss3d, 0); 
  ss.allocate_all() ; 
  ss.std_spectral_base();

  for (int l=1; l<4; l++)
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++)
	{
	  mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.00000000001, j, k, pipo,
			    lz, xi) ;
	  ss.set(l).set_grid_point(0, k, j, 0) = 
	    ss3(l).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	}



  // Computation of the 3-d projector on the 2-sphere
   
   
  Tensor proj3d(Kij.get_mp(),2,  type, Kij.get_mp().get_bvect_spher());
  Tensor proj_prov = gamij.con().down(1, gamij) - ss3dcon*ss3d;
  proj.allocate_all();
  proj.std_spectral_base();
  proj3d.allocate_all();
  proj3d = contract(jac, 1, contract( jac_inv, 0, proj_prov , 1), 1 ); 
       
  proj3d.std_spectral_base();
  
  for (int l=1; l<4; l++)
    for (int m=1; m<4; m++)
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
			      lz, xi) ;
	    proj.set(l,m).set_grid_point(0, k, j, 0) = 
	      proj3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
       
      
	  }   



  /* Computation of the metric linked to the 2-surface (linked to covariant form 
     of the degenerated 2-metric). 
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
  for (int l=1; l<4; l++)
    for (int m=1; m<4; m++) 
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
			      lz, xi) ;
	    qab2.set(l,m).set_grid_point(0, k, j, 0) = 
	      qq3d(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	  }
  qab= qab2;



  // Computation of the degenerated 3d degenerated covariant metric on the 2-surface 


  Sym_tensor qqq = contract(jac_inv, 0, contract( jac_inv, 0, (gamij.cov() - ss3d * ss3d) , 0), 1) ; 
 

  qq.allocate_all() ; 
  qq.std_spectral_base();
  for (int l=1; l<4; l++)
    for (int m=1; m<4; m++)     
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
			      lz, xi) ;
	    qq.set(l,m).set_grid_point(0, k, j, 0) = 
	      qqq(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;

	  }
       
  





  // Computation of the trace of the extrinsic curvature of 3-slice
  Scalar trk3d = Kij.trace(gamij) ;
  trk.allocate_all() ; 
  for (int k=0; k<np; k++)
    for (int j=0; j<nt; j++) {
      mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
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
      mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
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
      mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
			lz, xi) ;
      ggg.set_grid_point(0, k, j, 0) = 
	ggg3d.get_spectral_va().val_point_jk(lz, xi, j, k) ;
    }


  //Computation of function zeta, changing spheroid coordinates to get a round measure.

  Scalar ftilde = sqrt_q();
    double rayon = sqrt(area()/(4.*M_PI));
    ftilde = -ftilde/(rayon*rayon);
  ftilde.set_spectral_va().ylm();


  Base_val base = ftilde.get_spectral_base() ;

  Mtbl_cf *coefftilde = ftilde.get_spectral_va().c_cf;  int nombre = 2*nt; // Doubled in SYM base!!

  double *a_tilde = new double[nombre];
  

  lz = 0; // Now we work with 2d map associated with sqrt(q)

  
 
 
         for (int k=0; k<np+1; k++)
             for (int j=0; j<nt; j++) {
                int l_q, m_q, base_r ;
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
                 if (nullite_plm(j, nt, k, np, base) == 1) {
		   donne_lm(1, lz, j, k, base, m_q, l_q,base_r) ;
			   if (m_q ==0) {
		    a_tilde[l_q] = (*coefftilde)(0, k, j, 0);
	
			   }}}

	 Scalar zeta2(map); zeta2= 3.; zeta2.std_spectral_base(); zeta2.mult_cost();
	zeta2.set_spectral_va().ylm();

		Base_val base2 = zeta2.get_spectral_base() ; 
		Mtbl_cf *dzeta = zeta2.set_spectral_va().c_cf;

		for (int k=0; k<np+1; k++)
		  for (int j=0; j<nt; j++) {
		    int l_q2, m_q2, base_r2 ;
		base.give_quant_numbers(lz, k, j, m_q2, l_q2, base_r2) ;
                 if (nullite_plm(j, nt, k, np, base2) == 1) {
                donne_lm(1, lz, j, k, base2, m_q2, l_q2,base_r2) ;
		if (m_q2 ==0){
 		  if(l_q2 ==0){
 		    (*dzeta).set(0,k,j,0) = a_tilde[0] + (a_tilde[1]/3.)*sqrt(2.*1. + 1.);
		    } 
		  if (l_q2 >0){ 
		    (*dzeta).set(0,k,j,0) =  (a_tilde[l_q2 +1]/(2.*l_q2 + 3.))*sqrt((2.*(l_q2 +1.)+1.)/(2.*l_q2 + 1.)) -  (a_tilde[l_q2 -1]/(2.*l_q2 - 1.))*sqrt((2.*(l_q2 - 1.)+1.)/(2.*l_q2 + 1.));
		       }
                }
		 }
		  }
           
	
		zeta2.set_spectral_va().coef();
		zeta = zeta2;

     
  /* Computation of the tangent part of the extrinsic curvature of
   * the 2 surface embedded in the 3 slice. The reference vector
   used is the vector field s */

  Vector ll3d = contract( proj_prov, 0, contract(Kij, 1, ss3dcon, 0), 0) ; 

  Vector ll3 = contract( jac_inv, 0, ll3d, 0) ; 

  ll.allocate_all() ;
  ll.std_spectral_base(); 
  for (int l=1; l<4; l++)      
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++)
	{
	  mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
			    lz, xi) ;
	  ll.set(l).set_grid_point(0, k, j, 0) = 
	    ll3(l).get_spectral_va().val_point_jk(lz, xi, j, k) ;

	} 
 






    
  /* computation of the tangent components of the extrinsic curvature 
   *of the 3-slice 
   *(extracted from the curvature of the timeslice)
   * Note: this is not the actual 2d_ extrinsic curvature, but the 
   *tangent part of the time-slice extrinsic curvature */
 
  Tensor jj3d = Kij - ss3d*ll3d - ll3d*ss3d - contract(Kij, 0 , 1, sxss3dcon , 0, 1)* sxss3d ;
  Tensor jj3 =contract(jac_inv, 0 , contract(jac_inv,0 , jj3d,1),1);
  jj.allocate_all() ; 
  jj.std_spectral_base();
  for (int l=1; l<4; l++)
    for (int m=1; m<4; m++)     
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
			      lz, xi) ;
	    jj.set(l,m).set_grid_point(0, k, j, 0) = 
	      jj3(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;
	    
	  } 
  
  
  
  
  
  // Computation of 2d extrinsic curvature in the 3-slice
    

  Tensor hh3d = contract(proj_prov, 0, contract(proj_prov, 0,ss3d.derive_cov(gamij),1), 1 ) ;

  Tensor hh3 =contract(jac_inv, 0 , contract(jac_inv,0 , hh3d,1),1);

  hh.allocate_all() ; 
  hh.std_spectral_base();
  for (int l=1; l<4; l++)
    for (int m=1; m<4; m++)     
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  {
	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
			      lz, xi) ;
	    hh.set(l,m).set_grid_point(0, k, j, 0) = 
	      hh3(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;

	  }


 
 Tensor hh3dupdown = hh3d.up(0, gamij);
 Scalar ricciscal3 = gamij.ricci().trace(gamij); ricciscal3.annule_domain(nz-1); ricciscal3.std_spectral_base(); // Ricci scalar on the 3-surface.
 Tensor hh3dupup = hh3dupdown.up(1,gamij); hh3dupup.annule_domain(nz-1); hh3dupup.std_spectral_base();
 Scalar ricci22 = ricciscal3 - 2.*contract(contract(gamij.ricci(), 1, ss3dcon, 0),0, ss3dcon, 0);
 ricci22.annule_domain(nz-1); ricci22.std_spectral_base();
 ricci22 += (hh3d.trace(gamij)*hh3d.trace(gamij)) - contract(contract(hh3dupup,0, hh3d,0),0,1);
 ricci22.annule_domain(nz-1); ricci22.std_spectral_base();
 
  ricci.allocate_all();
  ricci.std_spectral_base();
      for (int k=0; k<np; k++)
        for (int j=0; j<nt; j++)
 	 {
 	   mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
 			     lz, xi) ;
 	   ricci.set_grid_point(0, k, j, 0) = 
 	     ricci22.get_spectral_va().val_point_jk(lz, xi, j, k) ;

	 }
 
 
 set_der_0x0() ;

}





//Copy constructor//

Spheroid::Spheroid (const Spheroid &sph_in) :h_surf(sph_in.h_surf),
					     jac2d(sph_in.jac2d),
					     proj(sph_in.proj),                          
					     qq(sph_in.qq),
					     ss (sph_in.ss),
					     qab(sph_in.qab),
					     ricci(sph_in.ricci),
					     hh(sph_in.hh),
					     trk(sph_in.trk),
					     ll(sph_in.ll),
					     jj(sph_in.jj),
					     fff(sph_in.fff),
					     ggg(sph_in.ggg),
					     zeta(sph_in.zeta),
					     issphere(sph_in.issphere)
  
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
  if (p_mass != 0x0) delete p_mass ;
  if (p_multipole_mass != 0x0) delete p_multipole_mass ;
  if (p_multipole_angu != 0x0) delete p_multipole_angu ;
  if (p_theta_plus != 0x0) delete p_theta_plus ;
  if (p_theta_minus != 0x0) delete p_theta_minus ;
  if (p_shear != 0x0) delete p_shear ;
  if (p_delta != 0x0) delete p_delta ;
  set_der_0x0() ;
}

void Spheroid::set_der_0x0() const {
  p_sqrt_q = 0x0 ;
  p_area = 0x0 ;
  p_angu_mom = 0x0 ;
  p_mass = 0x0 ;
  p_multipole_mass = 0x0;
  p_multipole_angu = 0x0;
  p_theta_plus = 0x0 ;
  p_theta_minus = 0x0 ;
  p_shear = 0x0 ;
  p_delta = 0x0;

} 


 
//---------//
//Accessors//
//---------//





// Computation of the 2-dimensional Jacobian amplitude for the surface
const  Scalar& Spheroid::sqrt_q() const { 
  if (p_sqrt_q == 0x0) {
    p_sqrt_q = new Scalar(sqrt((get_qq()(2,2)*get_qq()(3,3))- (get_qq()(2,3)*get_qq()(2,3)))) ;
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
    const Map_af& mp_ang = dynamic_cast<const Map_af&>(h_surf.get_mp()) ;
    Scalar tmp = contract(ll,0, contract (jac2d, 1,phi,0), 0 );
    p_angu_mom = new double (mp_ang.integrale_surface((sqrt_q()*h_surf*h_surf*tmp),1)) ;
    *p_angu_mom = *p_angu_mom /(8. * M_PI) ; 
  }

  return *p_angu_mom;

}


double Spheroid::mass(const Vector& phi) const {
  if (p_mass == 0x0) {
    assert( issphere == true ) ;
    double rayon_coord= h_surf.val_grid_point (0,0,0,0); // This is only the coordinate radius...
    double rayon = sqrt(area()/(4.*M_PI));

    p_mass = new double ((1/(2.*rayon))*sqrt(rayon*rayon*rayon*rayon + 4.*angu_mom(phi)*angu_mom(phi)));

  }
  return *p_mass;

}



double Spheroid::multipole_mass(const int order, const Vector& phi) const{
 
    assert( issphere == true ) ;
    const Map_af& mp_ang = dynamic_cast<const Map_af&>(h_surf.get_mp()) ;
    double rayon = sqrt(area()/(4.*M_PI));
    // Multiplicative factor before integral.
    double factor = mass(phi)/(8. * M_PI); // To check later
    if (order >0)
      { for (int compte=0; compte <=order -1; compte++)
	factor = factor*rayon;
      }
    // Calculus of legendre polynomial of order n, as function of cos(theta)
    Scalar Pn(mp_ang); Pn=1; Pn.std_spectral_base(); Pn.set_spectral_va().ylm();
    if (order >0)
      { Pn = Pn*zeta;     
      }
    if (order >1)
      { Scalar Pnold(mp_ang); Pnold = 1; Pnold.std_spectral_base(); Pnold.set_spectral_va().ylm();

      for (int nn=1; nn<order; nn++){
  
	Scalar Pnnew = (2.*nn +1.)*Pn;
	Pnnew = Pnnew*zeta;
	Pnnew = Pnnew - nn*Pnold;
	Pnnew = Pnnew/(double(nn) + 1.);
 
	Pnold = Pn;
	Pn = Pnnew;

      }
      }
 
    // Calculus of Ricci Scalar over the surface
    Scalar ricciscal(mp_ang);
    ricciscal = get_ricci(); // As it is 2-dimensional, ricci tensor is only proportional to ricci scalar.
    ricciscal.std_spectral_base();
    ricciscal.set_spectral_va().ylm();
  
    Scalar rayyon = h_surf;
    rayyon.std_spectral_base();
    rayyon.set_spectral_va().ylm();
    
    

    Scalar sqq = sqrt_q();
    Scalar integrande = sqq * rayyon *rayyon*ricciscal*Pn; integrande.std_spectral_base();

    
    p_multipole_mass = new double (factor*mp_ang.integrale_surface(integrande, 1)) ;
    
  
 

  return *p_multipole_mass;
}



double Spheroid::multipole_angu(const int order, const Vector& phi) const{

    assert( issphere == true ) ;
    assert (order >=1) ;
    const Map_af& mp_ang = dynamic_cast<const Map_af&>(h_surf.get_mp()) ;

  double rayon = sqrt(area()/(4.*M_PI));

    double factor = 1./(8. * M_PI);
    if (order >1)
      { for (int compte=0; compte <=order -2; compte++)
	factor = factor*rayon;
      }

    // Calculus of legendre polynomial of order n, as function of cos(theta)
    Scalar Pn(mp_ang); Pn=1; Pn.std_spectral_base(); Pn.set_spectral_va().ylm();
    Scalar dPn = Pn;

    Pn = Pn*zeta;     
  
    if (order >1)
   
      { Scalar Pnold(mp_ang); Pnold = 1; Pnold.std_spectral_base(); Pnold.set_spectral_va().ylm();

      for (int nn=1; nn<order; nn++){
  
	Scalar Pnnew = (2.*nn +1.)*Pn;
	Pnnew = Pnnew*zeta;
	Pnnew = Pnnew - nn*Pnold;
	Pnnew = Pnnew/(double(nn) + 1.);
 
	Pnold = Pn;
	Pn = Pnnew; // Pn is now P(n+1)

      }

      //  Calculus of functional derivative of order N legendre polynomial.
    
      dPn = Pn* zeta; dPn = dPn - Pnold; dPn = double(order)*dPn;
    
      Scalar quotient(mp_ang); quotient = 1.; quotient.std_spectral_base();
      quotient = quotient*zeta*zeta; quotient = quotient -1.;
    
      dPn = dPn/quotient; 

      }

    // Computation of the multipole;
    Scalar tmp = contract(ll,0, contract (jac2d, 1,phi,0), 0 ); tmp.std_spectral_base();
    Scalar tmp2 = (sqrt_q()) * h_surf *h_surf*tmp*dPn; tmp2.std_spectral_base();

     

    p_multipole_angu = new double (factor*mp_ang.integrale_surface(tmp2, 1)) ;  
   

  return *p_multipole_angu;
  
}




// Outgoing null expansion of 2-surface

const Scalar &Spheroid::theta_plus() const {

  if (p_theta_plus == 0x0) {
    p_theta_plus = new Scalar(fff*(hh.trace(qab) - jj.trace(qab))) ;
    p_theta_plus->std_spectral_base() ;
  }

  return *p_theta_plus; 
}








// ingoing null expansion of 2-surface

const Scalar& Spheroid::theta_minus() const {

  if (p_theta_minus == 0x0) {
    p_theta_minus = new Scalar(ggg*(-hh.trace(qab) - jj.trace(qab))) ;
    p_theta_minus->std_spectral_base() ;
  }

  return *p_theta_minus; 
 
}




//outer null-oriented shear of 2-surface

const Sym_tensor& Spheroid::shear() const { 
  
  if (p_shear == 0x0) {
     p_shear = new Sym_tensor( fff*(hh - jj) - (qab.cov()/2) *(hh.trace(qab) - jj.trace(qab))) ;  
// This is associated with the null vector: "l = n + s";
// For a null vector, "l = f (n + s)", multiply by the global factor "f" (e.g. the lapse "N" for "l=N(n+s)".  
                                                                                                
 
    p_shear->std_spectral_base() ;
  }

  return *p_shear; 
 
}











//-------------------------------------------//
// Covariant flat derivative, returning a pointer.//
//-------------------------------------------//

Tensor Spheroid::derive_cov2dflat(const Tensor& uu) const{

  // Notations: suffix 0 in name <=> input tensor
  //            suffix 1 in name <=> output tensor

  int valence0 = uu.get_valence() ; 
  int valence1 = valence0 + 1 ; 
  int valence1m1 = valence1 - 1 ; // same as valence0, but introduced for 
  // the sake of clarity


  // Protections
  // -----------
  if (valence0 >= 1) {

  }

  // Creation of the result (pointer)
  // --------------------------------
  Tensor *resu ;

  // If uu is a Scalar, the result is a vector
  if (valence0 == 0) {
    resu = new Vector(uu.get_mp(), COV, uu.get_mp().get_bvect_spher()) ;
  }
  else {

    // Type of indices of the result :
    Itbl tipe(valence1) ; 
    const Itbl& tipeuu = uu.get_index_type() ;  
    for (int id = 0; id<valence0; id++) {
      tipe.set(id) = tipeuu(id) ;   // First indices = same as uu
    }
    tipe.set(valence1m1) = COV ;  // last index is the derivation index

    // if uu is a Tensor_sym, the result is also a Tensor_sym:
    const Tensor* puu = &uu ; 
    const Tensor_sym* puus = dynamic_cast<const Tensor_sym*>(puu) ; 
    if ( puus != 0x0 ) {    // the input tensor is symmetric
      resu = new Tensor_sym(uu.get_mp(), valence1, tipe, *uu.get_triad(),
			    puus->sym_index1(), puus->sym_index2()) ;
    }
    else {  
      resu = new Tensor(uu.get_mp(), valence1, tipe, *uu.get_triad()) ;  // no symmetry  
    }

  }

  int ncomp1 = resu->get_n_comp() ; 
	
  Itbl ind1(valence1) ; // working Itbl to store the indices of resu
  Itbl ind0(valence0) ; // working Itbl to store the indices of uu
  Itbl ind(valence0) ;  // working Itbl to store the indices of uu
	
  Scalar tmp(uu.get_mp()) ;	// working scalar

  // Determination of the dzpuis parameter of the result  --> dz_resu
  // --------------------------------------------------
        
  int dz_resu = 0;

  // (We only work here on a non-compactified shell) // 




  // Loop on all the components of the output tensor
  // -----------------------------------------------
  /* Note: we have here preserved all the non-useful terms in this case(typically christoffel symbols) for the sake of understandng what's going on... 
   */
 
 
  for (int ic=0; ic<ncomp1; ic++) {
    
    // indices corresponding to the component no. ic in the output tensor
    ind1 = resu->indices(ic) ; 
    
    // Component no. ic:
    Scalar& cresu = resu->set(ind1) ; 
		
    // Indices of the input tensor
    for (int id = 0; id < valence0; id++) {
      ind0.set(id) = ind1(id) ; 
    }
         
    // Value of last index (derivation index)
    int k = ind1(valence1m1) ; 
        
    switch (k) {
        
    case 1 : {  // Derivation index = r
      //---------------------
	
      cresu = 0; //(uu(ind0)).dsdr() ; 	// d/dr
		
      // all the connection coefficients Gamma^i_{jk} are zero for k=1
      break ; 
    }   

    case 2 : {  // Derivation index = theta
      //-------------------------
			
      cresu = (uu(ind0)).srdsdt() ;  // 1/r d/dtheta 
		
      // Loop on all the indices of uu
      for (int id=0; id<valence0; id++) {
		
	switch ( ind0(id) ) {
				
	case 1 : {	// Gamma^r_{l theta} V^l 
	  // or -Gamma^l_{r theta} V_l 
	  /*   ind = ind0 ; 
	       ind.set(id) = 2 ;   // l = theta

	       // Division by r :
	       tmp = uu(ind) ; 
	       tmp.div_r_dzpuis(dz_resu) ;

	       cresu -= tmp ; */
	  break ; 
	}
		    		
	case 2 : {	// Gamma^theta_{l theta} V^l 
	  // or -Gamma^l_{theta theta} V_l
	  /*  ind = ind0 ; 
	      ind.set(id) = 1 ;   // l = r
	      tmp = uu(ind) ; 
	      tmp.div_r_dzpuis(dz_resu) ;

	      cresu += tmp ; */
	  break ; 
	}
				
	case 3 : {	// Gamma^phi_{l theta} V^l 
	  // or -Gamma^l_{phi theta} V_l
	  break ; 
	}
				
	default : {
	  cerr << "Connection_fspher::p_derive_cov : index problem ! "
	       << endl ; 
	  abort() ;  
	}
	}

      }
      break ; 
    }


    case 3 : {  // Derivation index = phi
      //-----------------------
					
      cresu = (uu(ind0)).srstdsdp() ;  // 1/(r sin(theta)) d/dphi 	
		
      // Loop on all the indices of uu
      for (int id=0; id<valence0; id++) {
		
	switch ( ind0(id) ) {
				
	case 1 : {	// Gamma^r_{l phi} V^l 
	  // or -Gamma^l_{r phi} V_l 
	  /* ind = ind0 ; 
	     ind.set(id) = 3 ;   // l = phi
	     tmp = uu(ind) ; 
	     tmp.div_r_dzpuis(dz_resu) ;

	     cresu -= tmp ; */
	  break ; 
	}
		    	
	case 2 : {	// Gamma^theta_{l phi} V^l 
	  // or -Gamma^l_{theta phi} V_l
	  ind = ind0 ; 
	  ind.set(id) = 3 ;   // l = phi
	  tmp = uu(ind) ; 
	  tmp.div_r_dzpuis(dz_resu) ;

	  tmp.div_tant() ; 	// division by tan(theta)
					
	  cresu -= tmp ; 
	  break ; 
	}
				
	case 3 : {	// Gamma^phi_{l phi} V^l 
	  // or -Gamma^l_{phi phi} V_l
						
	  ind = ind0 ; 
	  // 	            ind.set(id) = 1 ;   // l = r
	  // 	            tmp = uu(ind) ; 
	  // 		    tmp.div_r_dzpuis(dz_resu) ;

	  // 	            cresu += tmp ; 

	  ind.set(id) = 2 ;   // l = theta
	  tmp = uu(ind) ; 
	  tmp.div_r_dzpuis(dz_resu) ;

	  tmp.div_tant() ; 	// division by tan(theta)

	  cresu += tmp ; 
	  break ; 
	}
				
	default : {
	  cerr << "Connection_fspher::p_derive_cov : index problem ! \n"
	       << endl ; 
	  abort() ;  
	}
	}

      }
            
      break ; 
    }

    default : {
      cerr << "Connection_fspher::p_derive_cov : index problem ! \n" ;
      abort() ;  
    }

    } // End of switch on the derivation index


  } // End of loop on all the components of the output tensor

    // C'est fini !
    // -----------
  return *resu ; 

}

void Spheroid::sauve(FILE* ) const {

  cout << "c'est pas fait!" << endl ;
  return ; 

}








// Computation of the delta coefficients
 
 
const Tensor& Spheroid::delta() const {

  // Tensor tg 
  if (p_delta == 0x0) {
   
    Tensor christoflat(qab.get_mp(), 3, COV, qab.get_mp().get_bvect_spher()); 
    christoflat.set_index_type(0) = CON; 
    christoflat.set_etat_zero() ;

    // assert(flat_met != 0x0) ; 
    Tensor dgam = derive_cov2dflat(qab.cov()) ; 
 
    for (int k=1; k<=3; k++) {
      for (int i=1; i<=3; i++) {
	for (int j=1; j<=i; j++) {
	  Scalar& cc= christoflat.set(k,i,j); 
	  for (int l=1; l<=3; l++) {
	    cc += qab.con()(k,l) * ( 
				    dgam(l,j,i) + dgam(i,l,j) - dgam(i,j,l) ) ; 
                        
	  }
	  cc = 0.5 * cc ; 
	}
      }
    }

    p_delta = new Tensor (christoflat) ;
    
  }
  return *p_delta; 
}






// Computation of global derivative on 2-sphere 
Tensor Spheroid::derive_cov2d(const Tensor& uu) const {
  
  if(uu.get_valence()>=1){
    int nbboucle =  uu.get_valence(); 
    Tensor resu = derive_cov2dflat(uu);
    for (int y=1; y<=nbboucle; y++){
           
      int df = uu.get_index_type(y-1); 
      if (df == COV) {
	resu -= contract(delta(),0, uu, y-1);
      }
      else {resu += contract(delta(),1,  uu, y-1);}
   
      return resu;
    }
  }
  else return derive_cov2dflat(uu);  

  return derive_cov2dflat(uu); // to avoid warnings...
}
   








// // COmputation of the ricci tensor  

// const Sym_tensor& Spheroid::ricci() const {

//   if (p_ricci == 0x0) {  // a new computation is necessary
    
//     assert( issphere == true ) ;
//     Sym_tensor riccia(h_surf.get_mp(), CON, h_surf.get_mp().get_bvect_spher()) ;
//     riccia.set_etat_zero(); 
        
//     const Tensor& d_delta = derive_cov2dflat(delta()) ; 
                
//     for (int i=1; i<=3; i++) {
        
//       int jmax = 3 ; 
            
//       for (int j=1; j<=jmax; j++) {

// 	Scalar tmp1(h_surf.get_mp()) ;
// 	tmp1.set_etat_zero() ; 
// 	for (int k=1; k<=3; k++) {
// 	  tmp1 += d_delta(k,i,j,k) ; 
// 	} 
                
// 	Scalar tmp2(h_surf.get_mp()) ;
// 	tmp2.set_etat_zero() ; 
// 	for (int k=1; k<=3; k++) {
// 	  tmp2 += d_delta(k,i,k,j) ; 
// 	} 
                
// 	Scalar tmp3(h_surf.get_mp()) ;
// 	tmp3.set_etat_zero() ; 
// 	for (int k=1; k<=3; k++) {
// 	  for (int m=1; m<=3; m++) {
// 	    tmp3 += delta()(k,k,m) * delta()(m,i,j) ; 
// 	  }
// 	} 
// 	tmp3.dec_dzpuis() ;  // dzpuis 4 -> 3
                
// 	Scalar tmp4(h_surf.get_mp()) ;
// 	tmp4.set_etat_zero() ; 
// 	for (int k=1; k<=3; k++) {
// 	  for (int m=1; m<=3; m++) {
// 	    tmp4 += delta()(k,j,m) * delta()(m,i,k) ; 
// 	  }
// 	} 
// 	tmp4.dec_dzpuis() ;  // dzpuis 4 -> 3
              

// 	riccia.set(i,j) = tmp1 - tmp2 + tmp3 - tmp4 ; 
                  
        
//       }
//     }
//     /* Note: Here we must take into account the fact that a round metric on a spheroid doesn't give zero as "flat" ricci part. Then a diagonal scalar term must be added. 
//        WARNING: this only works with "round" horizons!! */ 
 
//     double rayon = sqrt(area()/(4.*M_PI));
//     Scalar rayon2  = h_surf;
//     rayon2 = rayon;
//     rayon2.std_spectral_base();
    
//     for (int hi=1; hi<=3; hi++){
 
//       riccia.set(hi,hi) += 2/(rayon2 * rayon2) ; // Plutot 1/hsurf^2, non? 
//     }
//     p_ricci = new Sym_tensor(riccia);
//   }
	
//   return *p_ricci ; 
	
// }




// COmputation of the ricci tensor  

// const Sym_tensor& Spheroid::ricci() const {

//   if (p_ricci == 0x0) {  // a new computation is necessary
//     Sym_tensor riccia(h_surf.get_mp(), CON, h_surf.get_mp().get_bvect_spher()) ;
//     Sym_tensor ricci3 = gamij.ricci();
    
//     Sym_tensor ricci3-2d(h_surf.get_mp(), COV, h_surf.get_mp().get_bvect_spher());
//       ricci3-2d.allocate_all() ; 
//   ricci3-2d.std_spectral_base();
//   for (int l=1; l<4; l++)
//     for (int m=1; m<4; m++)     
//       for (int k=0; k<np; k++)
// 	for (int j=0; j<nt; j++)
// 	  {
// 	    mp_rad->val_lx_jk(h_surf.val_grid_point(0, k, j, 0)*1.000000000001, j, k, pipo,
// 			      lz, xi) ;
// 	    ricci3-2d.set(l,m).set_grid_point(0, k, j, 0) = 
// 	      ricci3(l,m).get_spectral_va().val_point_jk(lz, xi, j, k) ;

// 	  }
       


//   }
//   return *p_ricci ; 
  
// }





