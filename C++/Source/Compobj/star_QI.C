/*
 *  Methods of the class Star_QI
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
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

char star_QI_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:52:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2012/12/03 15:27:30  c_some
 * Small changes
 *
 * Revision 1.3  2012/11/22 16:04:51  c_some
 * Minor modifications
 *
 * Revision 1.2  2012/11/21 14:55:27  c_some
 * Added methods fait_shift and fait_nphi
 *
 * Revision 1.1  2012/11/20 16:29:50  c_some
 * New class Star_QI
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>

// Lorene headers
#include "compobj.h"
#include "unites.h"

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
namespace Lorene {
Star_QI::Star_QI(Map& mpi) :
			 Compobj_QI(mpi) ,
			 logn(mpi), 
			 tnphi(mpi), 
			 nuf(mpi), 
			 nuq(mpi), 
			 dzeta(mpi), 
			 tggg(mpi), 
			 w_shift(mpi, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 ssjm1_nuf(mpi), 
			 ssjm1_nuq(mpi), 
			 ssjm1_dzeta(mpi), 
			 ssjm1_tggg(mpi), 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, CON, mp.get_bvect_cart())	
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

    // Initialization to a flat metric : 
    logn = 0 ;   
    tnphi = 0 ;   
    nuf = 0 ; 
    nuq = 0 ;
    dzeta = 0 ; 
    tggg = 0 ;   

    w_shift.set_etat_zero() ; 
    khi_shift =  0 ; 

    beta.set_etat_zero() ; 

    ssjm1_nuf = 0 ; 
    ssjm1_nuq = 0 ; 
    ssjm1_dzeta = 0 ; 
    ssjm1_tggg = 0 ; 
    ssjm1_khi = 0 ; 

    ssjm1_wshift.set_etat_zero() ; 

}

// Copy constructor
// --------------------
Star_QI::Star_QI(const Star_QI& st) :
			 Compobj_QI(st), 
			 logn(st.logn), 
			 tnphi(st.tnphi), 
			 nuf(st.nuf), 
			 nuq(st.nuq),
			 dzeta(st.dzeta), 
			 tggg(st.tggg), 
			 w_shift(st.w_shift), 
			 khi_shift(st.khi_shift), 
			 ssjm1_nuf(st.ssjm1_nuf), 
			 ssjm1_nuq(st.ssjm1_nuq), 
			 ssjm1_dzeta(st.ssjm1_dzeta), 
			 ssjm1_tggg(st.ssjm1_tggg), 
			 ssjm1_khi(st.ssjm1_khi), 
			 ssjm1_wshift(st.ssjm1_wshift)			 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}


// Constructor from a file
// -----------------------
Star_QI::Star_QI(Map& mpi, FILE* fich) :
			 Compobj_QI(mpi) , 
			 logn(mpi), 
			 tnphi(mpi), 
			 nuf(mpi), 
			 nuq(mpi), 
			 dzeta(mpi), 
			 tggg(mpi), 
			 w_shift(mpi, CON, mp.get_bvect_cart()), 
			 khi_shift(mpi), 
			 ssjm1_nuf(mpi), 
			 ssjm1_nuq(mpi), 
			 ssjm1_dzeta(mpi), 
			 ssjm1_tggg(mpi), 
			 ssjm1_khi(mpi), 
			 ssjm1_wshift(mpi, CON, mp.get_bvect_cart())			 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    // Read of the saved fields:
    // ------------------------

    Scalar nuf_file(mp, *(mp.get_mg()), fich) ; 
    nuf = nuf_file ; 
        
    Scalar nuq_file(mp, *(mp.get_mg()), fich) ; 
    nuq = nuq_file ; 
        
    logn = nuf + nuq ; //## to be checked !
    
    Scalar dzeta_file(mp, *(mp.get_mg()), fich) ; 
    dzeta = dzeta_file ; 
        
    Scalar tggg_file(mp, *(mp.get_mg()), fich) ; 
    tggg = tggg_file ; 
        
    Vector w_shift_file(mp, mp.get_bvect_cart(), fich) ; 
    w_shift = w_shift_file ;
    
    Scalar khi_shift_file(mp, *(mp.get_mg()), fich) ; 
    khi_shift = khi_shift_file ;
    
    fait_shift() ;	    // constructs shift from w_shift and khi_shift
    fait_nphi() ;       // constructs N^phi from (N^x,N^y,N^z)
    
    Scalar ssjm1_nuf_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_nuf = ssjm1_nuf_file ; 

    Scalar ssjm1_nuq_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_nuq = ssjm1_nuq_file ; 

    Scalar ssjm1_dzeta_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_dzeta = ssjm1_dzeta_file ; 

    Scalar ssjm1_tggg_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_tggg = ssjm1_tggg_file ; 

    Scalar ssjm1_khi_file(mp, *(mp.get_mg()), fich) ; 
    ssjm1_khi = ssjm1_khi_file ; 

    Vector ssjm1_wshift_file(mp, mp.get_bvect_cart(), fich) ; 
    ssjm1_wshift = ssjm1_wshift_file ; 

     
    // Initialization of N, A^2, B^2, gamma_ij, tkij and ak_car
    update_metric() ; 
    
}

			    //------------//
			    // Destructor //
			    //------------//

Star_QI::~Star_QI(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Star_QI::del_deriv() const {

    Compobj_QI::del_deriv() ; 

   	if (p_grv2 != 0x0) delete p_grv2 ; 
    if (p_grv3 != 0x0) delete p_grv3 ;
    if (p_mom_quad != 0x0) delete p_mom_quad ;
    if (p_mass_g != 0x0) delete p_mass_g ;

    Star_QI::set_der_0x0() ; 
}			    


void Star_QI::set_der_0x0() const {

    p_grv2 = 0x0 ; 
    p_grv3 = 0x0 ;
    p_mom_quad = 0x0 ;
    p_mass_g = 0x0 ;
 	 
}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Star_QI
// --------------------------------
void Star_QI::operator=(const Star_QI& st) {

    // Assignment of quantities common to all the derived classes of Compobj_QI
    Compobj_QI::operator=(st) ;	    
    
    logn = st.logn ; 
    tnphi = st.tnphi ; 
    nuf = st.nuf ; 
    nuq = st.nuq ; 
    dzeta = st.dzeta ; 
    tggg = st.tggg ; 
    w_shift = st.w_shift ;
    khi_shift = st.khi_shift ;
    ssjm1_nuf = st.ssjm1_nuf ;
    ssjm1_nuq = st.ssjm1_nuq ;
    ssjm1_dzeta = st.ssjm1_dzeta ; 
    ssjm1_tggg = st.ssjm1_tggg ;
    ssjm1_khi = st.ssjm1_khi ;
    ssjm1_wshift = st.ssjm1_wshift ; 

    del_deriv() ;  // Deletes all derived quantities
}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Star_QI::sauve(FILE* fich) const {

    nuf.sauve(fich) ; 
    nuq.sauve(fich) ; 
    dzeta.sauve(fich) ; 
    tggg.sauve(fich) ; 
    w_shift.sauve(fich) ; 
    khi_shift.sauve(fich) ; 
    
    ssjm1_nuf.sauve(fich) ; 
    ssjm1_nuq.sauve(fich) ; 
    ssjm1_dzeta.sauve(fich) ; 
    ssjm1_tggg.sauve(fich) ; 
    ssjm1_khi.sauve(fich) ; 
    ssjm1_wshift.sauve(fich) ; 
    
}

// Printing
// --------

ostream& Star_QI::operator>>(ostream& ost) const {

 	using namespace Unites ;
	
	Compobj_QI::operator>>(ost) ; 
    
    ost << endl << "Axisymmetric stationary compact star in quasi-isotropic coordinates (class Star_QI) " << endl ; 

    ost << "Central values of various fields : " << endl ; 
    ost << "-------------------------------- " << endl ; 
    ost << "   ln(N) : " << logn.val_grid_point(0,0,0,0) << endl ; 
    ost << "   nuf : " << nuf.val_grid_point(0,0,0,0) << endl ; 
    ost << "   nuq : " << nuq.val_grid_point(0,0,0,0) << endl ; 
    ost << "   zeta = ln(AN): " << dzeta.val_grid_point(0,0,0,0) << endl << endl ;
    
    ost << "Error on the virial identity GRV2 : " <<  grv2() << endl ; 
    ost << "Error on the virial identity GRV3 : " <<  grv3(&ost) << endl ; 

    double mom_quad_38si = mom_quad() * rho_unit * (pow(r_unit, double(5.)) 
 							/ double(1.e38) ) ;
    ost << "Quadrupole moment Q : " << mom_quad_38si << " 10^38 kg m^2"
         << endl ; 
    ost << "c^4 Q / (G^2 M^3) :      " 
	 << mom_quad() / ( pow(qpig/(4*M_PI), 2.) * pow(mass_g(), 3.) ) 
	 << endl ; 
    
    ost << "Total angular momentum J :      " 
	 << angu_mom()/( qpig / (4* M_PI) * msol*msol) << " G M_sol^2 / c"
	 << endl ; 
    ost << "c J / (G M^2) :           " 
	 << angu_mom()/( qpig / (4* M_PI) * pow(mass_g(), 2.) ) << endl ;    	 	
    	
    return ost ; 
      
}

			    //-------------------------//
			    //	Computational methods  //
			    //-------------------------//
			    
void Star_QI::fait_shift() {

    Vector d_khi = khi_shift.derive_con( mp.flat_met_cart() ) ; 
    
    d_khi.dec_dzpuis(2) ;   // divide by r^2 in the external compactified domain  
 
    // x_k dW^k/dx_i
      Scalar xx(mp) ; 
      Scalar yy(mp) ; 
      Scalar zz(mp) ; 
      Scalar sintcosp(mp) ;
      Scalar sintsinp(mp) ;
      Scalar cost(mp) ;
      xx = mp.x ; 
      yy = mp.y ; 
      zz = mp.z ; 
      sintcosp = mp.sint * mp.cosp ;
      sintsinp = mp.sint * mp.sinp ;
      cost = mp.cost ;

      int nz = mp.get_mg()->get_nzone() ;
      Vector xk(mp, COV, mp.get_bvect_cart()) ;
      xk.set(1) = xx ;
      xk.set(1).set_domain(nz-1) = sintcosp.domain(nz-1) ;
      xk.set(1).set_dzpuis(-1) ;
      xk.set(2) = yy ;
      xk.set(2).set_domain(nz-1) = sintsinp.domain(nz-1) ;
      xk.set(2).set_dzpuis(-1) ;
      xk.set(3) = zz ;
      xk.set(3).set_domain(nz-1) = cost.domain(nz-1) ;
      xk.set(3).set_dzpuis(-1) ;
      xk.std_spectral_base() ;
      
      Tensor d_w = w_shift.derive_con( mp.flat_met_cart() ) ;
      
      Vector x_d_w = contract(xk, 0, d_w, 0) ;
      x_d_w.dec_dzpuis() ;

      double lambda = double(1) / double(3) ; 

      beta = - (lambda+2)/2./(lambda+1) * w_shift
		+ (lambda/2./(lambda+1))  * (d_khi + x_d_w) ;      
    
} 


void Star_QI::fait_nphi() {

    if ( (beta(1).get_etat() == ETATZERO) && (beta(2).get_etat() == ETATZERO) ) {
		tnphi = 0 ; 
		nphi = 0 ;
	return ;  
    }

    // Computation of tnphi
    // --------------------
    
    mp.comp_p_from_cartesian( -beta(1), -beta(2), tnphi ) ; 

    // Computation of nphi
    // -------------------
    
    nphi = tnphi ; 
    nphi.div_rsint() ; 
    
}


// Updates the 3-metric and the shift

void Star_QI::update_metric() {

    // Lapse function N
    // ----------------
    
    nn = exp( logn ) ; 

    nn.std_spectral_base() ;   // set the bases for spectral expansions
    
    
    // Metric factor A^2
    // -----------------
    
    a_car = exp( 2*( dzeta - logn ) ) ; 

    a_car.std_spectral_base() ;   // set the bases for spectral expansions

    // Metric factor B
    // ---------------
    
    Scalar tmp = tggg ; 
    tmp.div_rsint() ;	        //... Division of tG by r sin(theta)

    bbb = (1 + tmp) / nn ; 

    bbb.std_spectral_base() ;   // set the bases for spectral expansions
        
    b_car = bbb * bbb ; 
	

	Compobj_QI::update_metric() ;  // updates gamma_{ij} 
	
	
    // Tensor B^{-2} K_{ij} and Scalar A^2 K_{ij} K^{ij}
    // -------------------------------------------------
    
    //## extrinsic_curvature() ; // should be done by Compobj_QI::update_metric()
    
  
    // The derived quantities are no longer up to date : 
    // -----------------------------------------------

    del_deriv() ;  
	
}



}
