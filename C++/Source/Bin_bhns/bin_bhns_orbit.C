/*
 *  Methods of class Bin_bhns to compute the orbital angular velocity
 *
 *    (see file bin_bhns.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005-2007 Keisuke Taniguchi
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

char bin_bhns_orbit_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:52:41  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:00  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/05/15 19:01:28  k_taniguchi
 * Change of some parameters.
 *
 * Revision 1.1  2007/06/22 01:10:20  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C++ headers
//#include <>

// C headers
#include <cmath>

// Lorene headers
#include "bin_bhns.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
double func_binbhns_orbit_ks(double , const Param& ) ;
double func_binbhns_orbit_is(double , const Param& ) ;

//**********************************************************************

void Bin_bhns::orbit_omega(double fact_omeg_min, double fact_omeg_max) {

    using namespace Unites ;

    if (hole.is_kerrschild()) {

      //--------------------------------------------------------------------
      // Evaluation of various quantities at the center of the neutron star
      //--------------------------------------------------------------------

      double dnulg, p6sl2, dp6sl2 ;
      double shiftx, shifty, dshiftx, dshifty, shift2, dshift2 ;
      double x_orb, y_orb, y_separ, xbh_orb, mhsr ;

      const Map& mp = star.get_mp() ;

      const Scalar& lapconf = star.get_lapconf_tot() ;
      const Scalar& lapconf_auto = star.get_lapconf_auto() ;
      const Scalar& confo = star.get_confo_tot() ;
      const Scalar& confo_auto = star.get_confo_auto() ;
      const Scalar& loggam = star.get_loggam() ;
      const Vector& shift = star.get_shift_tot() ;
      const Vector& shift_auto = star.get_shift_auto() ;

      const Vector& dlapconf_comp = star.get_d_lapconf_comp() ;
      const Vector& dconfo_comp = star.get_d_confo_comp() ;
      const Tensor& dshift_comp = star.get_d_shift_comp() ;

      const double& massbh = hole.get_mass_bh() ;
      double mass = ggrav * massbh ;

      //----------------------------------------------------------
      // Calculation of d/dX( ln(lapconf) - ln(psi) + ln(Gamma) )
      //  at the center of NS
      //----------------------------------------------------------

      // Factor to translate x --> X
      double factx ;
      if (fabs(mp.get_rot_phi()) < 1.e-14) {
        factx = 1. ;
      }
      else {
        if (fabs(mp.get_rot_phi() - M_PI) < 1.e-14) {
	  factx = - 1. ;
	}
	else {
	  cout << "Bin_bhns::orbit_omega : unknown value of rot_phi !"
	       << endl ;
	  abort() ;
	}
      }

      Scalar tmp1(mp) ;
      tmp1 = loggam ;
      tmp1.std_spectral_base() ;

      // d/dX tmp1
      dnulg = factx * ( ((lapconf_auto.dsdx()).val_grid_point(0,0,0,0)
			 + dlapconf_comp(1).val_grid_point(0,0,0,0))
			/ lapconf.val_grid_point(0,0,0,0)
			- ((confo_auto.dsdx()).val_grid_point(0,0,0,0)
			   + dconfo_comp(1).val_grid_point(0,0,0,0))
			/ confo.val_grid_point(0,0,0,0)
			+ (tmp1.dsdx()).val_grid_point(0,0,0,0) ) ;


      //----------------------------------------------------
      // Calculation of psi^6/lapconf^2 at the center of NS
      //----------------------------------------------------

      double lapconf_c = lapconf.val_grid_point(0,0,0,0) ;
      double confo_c = confo.val_grid_point(0,0,0,0) ;

      p6sl2 = pow(confo_c,6.) / lapconf_c / lapconf_c ;


      //----------------------------------------------------------
      // Calculation of d/dX(psi^6/lapconf^2) at the center of NS
      //----------------------------------------------------------

      double dlapconf_c = factx *
	( (lapconf_auto.dsdx()).val_grid_point(0,0,0,0)
	  + dlapconf_comp(1).val_grid_point(0,0,0,0) ) ;

      double dpsi6_c = 6. * factx * pow(confo_c,5.)
	* ((confo_auto.dsdx()).val_grid_point(0,0,0,0)
	   + dconfo_comp(1).val_grid_point(0,0,0,0)) ;

      dp6sl2 = (dpsi6_c - 2.*pow(confo_c,6.)*dlapconf_c/lapconf_c)
	/ lapconf_c / lapconf_c ;


      //--------------------------------------------------------
      // Calculation of shift^X and shift^Y at the center of NS
      //--------------------------------------------------------

      shiftx = shift(1).val_grid_point(0,0,0,0) ;
      shifty = shift(2).val_grid_point(0,0,0,0) ;


      //------------------------------------------------------------------
      // Calculation of d shift^X/dX and d shift^Y/dX at the center of NS
      //------------------------------------------------------------------

      dshiftx = factx * ((shift_auto(1).dsdx()).val_grid_point(0,0,0,0)
			 + dshift_comp(1,1).val_grid_point(0,0,0,0)) ;

      dshifty = factx * ((shift_auto(2).dsdx()).val_grid_point(0,0,0,0)
			 + dshift_comp(1,2).val_grid_point(0,0,0,0)) ;


      //--------------------------------------------------------
      // Calculation of (shift^X)^2 + (shift^Y)^2 + (shift^Z)^2
      //  at the center of NS
      //--------------------------------------------------------

      Scalar tmp2(mp) ;
      tmp2 = shift(1)%shift(1) + shift(2)%shift(2) + shift(3)%shift(3) ;
      tmp2.std_spectral_base() ;

      shift2 = tmp2.val_grid_point(0,0,0,0) ;


      //----------------------------------------------------------------
      // Calculation of d/dX( (shift^X)^2 + (shift^Y)^2 + (shift^Z)^2 )
      //  at the center of NS
      //----------------------------------------------------------------

      dshift2 = 2.*factx*((shift(1).val_grid_point(0,0,0,0))
			  * ((shift_auto(1).dsdx()).val_grid_point(0,0,0,0)
			     + dshift_comp(1,1).val_grid_point(0,0,0,0))
			  +(shift(2).val_grid_point(0,0,0,0))
			  * ((shift_auto(2).dsdx()).val_grid_point(0,0,0,0)
			     + dshift_comp(1,2).val_grid_point(0,0,0,0))
			  +(shift(3).val_grid_point(0,0,0,0))
			  * ((shift_auto(3).dsdx()).val_grid_point(0,0,0,0)
			     + dshift_comp(1,3).val_grid_point(0,0,0,0)) ) ;


      //-------------------------
      // Information of position
      //-------------------------

      x_orb = (star.get_mp()).get_ori_x() - x_rot ;
      y_orb = (star.get_mp()).get_ori_y() - y_rot ;
      y_separ = (star.get_mp()).get_ori_y() - (hole.get_mp()).get_ori_y() ;

      xbh_orb = (hole.get_mp()).get_ori_x() - x_rot ;

       //------------------------------
      // Calculation of H_BH / r_bh^4
      //------------------------------

      mhsr = mass / pow( separ*separ+y_separ*y_separ, 2.5 ) ;

      // Output
      // ------

      cout << "Bin_bhns::orbit_omega: central d(log(lap)+log(Gam))/dX : "
	   << dnulg << endl ;
      cout << "Bin_bhns::orbit_omega: central psi^6/lapconf^2 :         "
	   << p6sl2 << endl ;
      cout << "Bin_bhns::orbit_omega: central d(psi^6/lapconf^2)/dX :   "
	   << dp6sl2 << endl ;
      cout << "Bin_bhns::orbit_omega: central shift^X :                 "
	   << shiftx << endl ;
      cout << "Bin_bhns::orbit_omega: central shift^Y :                 "
	   << shifty << endl ;
      cout << "Bin_bhns::orbit_omega: central d(shift^X)/dX :           "
	   << dshiftx << endl ;
      cout << "Bin_bhns::orbit_omega: central d(shift^Y)/dX :           "
	   << dshifty << endl ;
      cout << "Bin_bhns::orbit_omega: central shift^i shift_i :         "
	   << shift2 << endl ;
      cout << "Bin_bhns::orbit_omega: central d(shift^i shift_i)/dX :   "
	   << dshift2 << endl ;


      //---------------------------------------------------------------//
      //          Calculation of the orbital angular velocity          //
      //---------------------------------------------------------------//

      Param parorb ;
      parorb.add_double(dnulg, 0) ;
      parorb.add_double(p6sl2, 1) ;
      parorb.add_double(dp6sl2, 2) ;
      parorb.add_double(shiftx, 3) ;
      parorb.add_double(shifty, 4) ;
      parorb.add_double(dshiftx, 5) ;
      parorb.add_double(dshifty, 6) ;
      parorb.add_double(shift2, 7) ;
      parorb.add_double(dshift2, 8) ;
      parorb.add_double(x_orb, 9) ;
      parorb.add_double(y_orb, 10) ;
      parorb.add_double(separ, 11) ;
      parorb.add_double(y_separ, 12) ;
      parorb.add_double(xbh_orb, 13) ;
      parorb.add_double(mhsr, 14) ;

      double omega1 = fact_omeg_min * omega ;
      double omega2 = fact_omeg_max * omega ;

      cout << "Bin_bhns::orbit_omega: omega1,  omega2 [rad/s] : "
	   << omega1 * f_unit << "  " << omega2 * f_unit << endl ;

      // Search for the various zeros in the interval [omega1,omega2]
      // ------------------------------------------------------------

      int nsub = 50 ;  // total number of subdivisions of the interval
      Tbl* azer = 0x0 ;
      Tbl* bzer = 0x0 ;
      zero_list(func_binbhns_orbit_ks, parorb, omega1, omega2, nsub,
		azer, bzer) ;

      // Search for the zero closest to the previous value of omega
      // ----------------------------------------------------------

      double omeg_min, omeg_max ;
      int nzer = azer->get_taille() ; // number of zeros found by zero_list

      cout << "Bin_bhns::orbit_omega: " << nzer
	   << " zero(s) found in the interval [omega1,  omega2]." << endl ;
      cout << "omega, omega1, omega2 : " << omega << "  " << omega1
	   << "  " << omega2 << endl ;
      cout << "azer : " << *azer << endl ;
      cout << "bzer : " << *bzer << endl ;

      if (nzer == 0) {
        cout <<
	  "Bin_bhns::orbit_omega: WARNING : no zero detected in the interval"
	     << endl << "   [" << omega1 * f_unit << ", "
	     << omega2 * f_unit << "]  rad/s  !" << endl ;
	omeg_min = omega1 ;
	omeg_max = omega2 ;
      }
      else {
        double dist_min = fabs(omega2 - omega1) ;
	int i_dist_min = -1 ;
	for (int i=0; i<nzer; i++) {
	  // Distance of previous value of omega from the center of the
	  //  interval [azer(i), bzer(i)]

	  double dist = fabs( omega - 0.5 * ( (*azer)(i) + (*bzer)(i) ) ) ;

	  if (dist < dist_min) {
	    dist_min = dist ;
	    i_dist_min = i ;
	  }
	}
	omeg_min = (*azer)(i_dist_min) ;
	omeg_max = (*bzer)(i_dist_min) ;
      }

      delete azer ; // Tbl allocated by zero_list
      delete bzer ;

      cout << "Bin_bhns::orbit_omega: interval selected for the search of the zero : "
	   << endl << "  [" << omeg_min << ", " << omeg_max << "]  =  ["
	   << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s " << endl ;

      // Computation of the zero in the selected interval by the secant method
      // ---------------------------------------------------------------------

      int nitermax = 200 ;
      int niter ;
      double precis = 1.e-13 ;
      omega = zerosec_b(func_binbhns_orbit_ks, parorb, omeg_min, omeg_max,
			precis, nitermax, niter) ;

      cout << "Bin_bhns::orbit_omega: Number of iterations in zerosec for omega : "
	   << niter << endl ;

      cout << "Bin_bhns::orbit_omega: omega [rad/s] : "
	   << omega * f_unit << endl ;


    }
    else { // Isotropic coordinates with the maximal slicing

      //--------------------------------------------------------------------
      // Evaluation of various quantities at the center of the neutron star
      //--------------------------------------------------------------------

      double dnulg, p6sl2, dp6sl2 ;
      double shiftx, shifty, dshiftx, dshifty, shift2, dshift2 ;
      double x_orb, y_orb ;

      const Map& mp = star.get_mp() ;

      const Scalar& lapconf = star.get_lapconf_tot() ;
      const Scalar& lapconf_auto = star.get_lapconf_auto() ;
      const Scalar& confo = star.get_confo_tot() ;
      const Scalar& confo_auto = star.get_confo_auto() ;
      const Scalar& loggam = star.get_loggam() ;
      const Vector& shift = star.get_shift_tot() ;
      const Vector& shift_auto = star.get_shift_auto() ;

      const Vector& dlapconf_comp = star.get_d_lapconf_comp() ;
      const Vector& dconfo_comp = star.get_d_confo_comp() ;
      const Tensor& dshift_comp = star.get_d_shift_comp() ;

      //----------------------------------------------------------
      // Calculation of d/dX( ln(lapconf) - ln(psi) + ln(Gamma) )
      //  at the center of NS
      //----------------------------------------------------------

      // Factor to translate x --> X
      double factx ;
      if (fabs(mp.get_rot_phi()) < 1.e-14) {
        factx = 1. ;
      }
      else {
        if (fabs(mp.get_rot_phi() - M_PI) < 1.e-14) {
	  factx = - 1. ;
	}
	else {
	  cout << "Bin_bhns::orbit_omega: unknown value of rot_phi !"
	       << endl ;
	  abort() ;
	}
      }

      Scalar tmp1(mp) ;
      tmp1 = loggam ;
      tmp1.std_spectral_base() ;

      // d/dX tmp1
      dnulg = factx * ( ((lapconf_auto.dsdx()).val_grid_point(0,0,0,0)
			 + dlapconf_comp(1).val_grid_point(0,0,0,0))
			/ lapconf.val_grid_point(0,0,0,0)
			- ((confo_auto.dsdx()).val_grid_point(0,0,0,0)
			 + dconfo_comp(1).val_grid_point(0,0,0,0))
			/ confo.val_grid_point(0,0,0,0)
			+ (tmp1.dsdx()).val_grid_point(0,0,0,0) ) ;


      //----------------------------------------------------
      // Calculation of psi^6/lapconf^2 at the center of NS
      //----------------------------------------------------

      double lapconf_c = lapconf.val_grid_point(0,0,0,0) ;
      double confo_c = confo.val_grid_point(0,0,0,0) ;

      p6sl2 = pow(confo_c,6.) / lapconf_c / lapconf_c ;


      //----------------------------------------------------------
      // Calculation of d/dX(psi^6/lapconf^2) at the center of NS
      //----------------------------------------------------------

      double dlapconf_c = factx *
	( (lapconf_auto.dsdx()).val_grid_point(0,0,0,0)
	  + dlapconf_comp(1).val_grid_point(0,0,0,0) ) ;

      double dpsi6_c = 6. * factx * pow(confo_c,5.)
	* ((confo_auto.dsdx()).val_grid_point(0,0,0,0)
	   + dconfo_comp(1).val_grid_point(0,0,0,0)) ;

      dp6sl2 = (dpsi6_c - 2.*pow(confo_c,6.)*dlapconf_c/lapconf_c)
	/ lapconf_c / lapconf_c ;


      //--------------------------------------------------------
      // Calculation of shift^X and shift^Y at the center of NS
      //--------------------------------------------------------

      shiftx = shift(1).val_grid_point(0,0,0,0) ;
      shifty = shift(2).val_grid_point(0,0,0,0) ;


      //------------------------------------------------------------------
      // Calculation of d shift^X/dX and d shift^Y/dX at the center of NS
      //------------------------------------------------------------------

      dshiftx = factx * ( (shift_auto(1).dsdx()).val_grid_point(0,0,0,0)
			  + dshift_comp(1,1).val_grid_point(0,0,0,0) ) ;

      dshifty = factx * ( (shift_auto(2).dsdx()).val_grid_point(0,0,0,0)
			  + dshift_comp(1,2).val_grid_point(0,0,0,0) ) ;


      //--------------------------------------------------------
      // Calculation of (shift^X)^2 + (shift^Y)^2 + (shift^Z)^2
      //  at the center of NS
      //--------------------------------------------------------

      Scalar tmp2(mp) ;
      tmp2 = shift(1)%shift(1) + shift(2)%shift(2) + shift(3)%shift(3) ;
      tmp2.std_spectral_base() ;

      shift2 = tmp2.val_grid_point(0,0,0,0) ;


      //----------------------------------------------------------------
      // Calculation of d/dX( (shift^X)^2 + (shift^Y)^2 + (shift^Z)^2 )
      //  at the center of NS
      //----------------------------------------------------------------

      dshift2 = 2.*factx*( (shift(1).val_grid_point(0,0,0,0))
			   * ((shift_auto(1).dsdx()).val_grid_point(0,0,0,0)
			      + dshift_comp(1,1).val_grid_point(0,0,0,0))
			   +(shift(2).val_grid_point(0,0,0,0))
			   * ((shift_auto(2).dsdx()).val_grid_point(0,0,0,0)
			      + dshift_comp(1,2).val_grid_point(0,0,0,0))
			   +(shift(3).val_grid_point(0,0,0,0))
			   * ((shift_auto(3).dsdx()).val_grid_point(0,0,0,0)
			      + dshift_comp(1,3).val_grid_point(0,0,0,0)) ) ;


      //-------------------------
      // Information of position
      //-------------------------

      x_orb = (star.get_mp()).get_ori_x() - x_rot ;
      y_orb = (star.get_mp()).get_ori_y() - y_rot ;


      // Output
      // ------

      cout << "Bin_bhns::orbit_omega: central d(log(lap)+log(Gam))/dX : "
	   << dnulg << endl ;
      cout << "Bin_bhns::orbit_omega: central psi^6/lapconf^2 :         "
	   << p6sl2 << endl ;
      cout << "Bin_bhns::orbit_omega: central d(psi^6/lapconf^2)/dX :   "
	   << dp6sl2 << endl ;
      cout << "Bin_bhns::orbit_omega: central shift^X :                 "
	   << shiftx << endl ;
      cout << "Bin_bhns::orbit_omega: central shift^Y :                 "
	   << shifty << endl ;
      cout << "Bin_bhns::orbit_omega: central d(shift^X)/dX :           "
	   << dshiftx << endl ;
      cout << "Bin_bhns::orbit_omega: central d(shift^Y)/dX :           "
	   << dshifty << endl ;
      cout << "Bin_bhns::orbit_omega: central shift^i shift_i :         "
	   << shift2 << endl ;
      cout << "Bin_bhns::orbit_omega: central d(shift^i shift_i)/dX :   "
	   << dshift2 << endl ;


      //---------------------------------------------------------------//
      //          Calculation of the orbital angular velocity          //
      //---------------------------------------------------------------//

      Param parorb ;
      parorb.add_double(dnulg, 0) ;
      parorb.add_double(p6sl2, 1) ;
      parorb.add_double(dp6sl2, 2) ;
      parorb.add_double(shiftx, 3) ;
      parorb.add_double(shifty, 4) ;
      parorb.add_double(dshiftx, 5) ;
      parorb.add_double(dshifty, 6) ;
      parorb.add_double(shift2, 7) ;
      parorb.add_double(dshift2, 8) ;
      parorb.add_double(x_orb, 9) ;
      parorb.add_double(y_orb, 10) ;

      double omega1 = fact_omeg_min * omega ;
      double omega2 = fact_omeg_max * omega ;

      cout << "Bin_bhns::orbit_omega: omega1,  omega2 [rad/s] : "
	   << omega1 * f_unit << "  " << omega2 * f_unit << endl ;

      // Search for the various zeros in the interval [omega1,omega2]
      // ------------------------------------------------------------

      int nsub = 50 ;  // total number of subdivisions of the interval
      Tbl* azer = 0x0 ;
      Tbl* bzer = 0x0 ;
      zero_list(func_binbhns_orbit_is, parorb, omega1, omega2, nsub,
		azer, bzer) ;

      // Search for the zero closest to the previous value of omega
      // ----------------------------------------------------------

      double omeg_min, omeg_max ;
      int nzer = azer->get_taille() ; // number of zeros found by zero_list

      cout << "Bin_bhns::orbit_omega: " << nzer
	   << " zero(s) found in the interval [omega1,  omega2]." << endl ;
      cout << "omega, omega1, omega2 : " << omega << "  " << omega1
	   << "  " << omega2 << endl ;
      cout << "azer : " << *azer << endl ;
      cout << "bzer : " << *bzer << endl ;

      if (nzer == 0) {
        cout <<
	  "Bin_bhns::orbit_omega: WARNING : no zero detected in the interval"
	     << endl << "   [" << omega1 * f_unit << ", "
	     << omega2 * f_unit << "]  rad/s  !" << endl ;
	omeg_min = omega1 ;
	omeg_max = omega2 ;
      }
      else {
        double dist_min = fabs(omega2 - omega1) ;
	int i_dist_min = -1 ;
	for (int i=0; i<nzer; i++) {
	  // Distance of previous value of omega from the center of the
	  //  interval [azer(i), bzer(i)]

	  double dist = fabs( omega - 0.5 * ( (*azer)(i) + (*bzer)(i) ) ) ;

	  if (dist < dist_min) {
	    dist_min = dist ;
	    i_dist_min = i ;
	  }
	}
	omeg_min = (*azer)(i_dist_min) ;
	omeg_max = (*bzer)(i_dist_min) ;
      }

      delete azer ; // Tbl allocated by zero_list
      delete bzer ;

      cout << "Bin_bhns::orbit_omega: interval selected for the search of the zero : "
	   << endl << "  [" << omeg_min << ", " << omeg_max << "]  =  ["
	   << omeg_min * f_unit << ", " << omeg_max * f_unit << "] rad/s " << endl ;

      // Computation of the zero in the selected interval by the secant method
      // ---------------------------------------------------------------------

      int nitermax = 200 ;
      int niter ;
      double precis = 1.e-13 ;
      omega = zerosec_b(func_binbhns_orbit_is, parorb, omeg_min, omeg_max,
			precis, nitermax, niter) ;

      cout << "Bin_bhns::orbit_omega: Number of iterations in zerosec for omega : "
	   << niter << endl ;

      cout << "Bin_bhns::orbit_omega: omega [rad/s] : "
	   << omega * f_unit << endl ;

    }
}

//******************************
// Function for searching omega
//******************************

double func_binbhns_orbit_ks(double om, const Param& parorb) {

    double dnulg = parorb.get_double(0) ;
    double p6sl2 = parorb.get_double(1) ;
    double dp6sl2 = parorb.get_double(2) ;
    double shiftx = parorb.get_double(3) ;
    double shifty = parorb.get_double(4) ;
    double dshiftx = parorb.get_double(5) ;
    double dshifty = parorb.get_double(6) ;
    double shift2 = parorb.get_double(7) ;
    double dshift2 = parorb.get_double(8) ;
    double x_orb = parorb.get_double(9) ;
    double y_orb = parorb.get_double(10) ;
    double x_separ = parorb.get_double(11) ;
    double y_separ = parorb.get_double(12) ;
    double xbh_orb = parorb.get_double(13) ;
    double mhsr = parorb.get_double(14) ;

    double om2 = om * om ;
    /*
    double bpb = om2 * (x_orb * x_orb + y_orb * y_orb)
      + 2. * om * (shifty * x_orb - shiftx * y_orb) + shift2
      + 2. * mhsr * (x_separ*x_separ+y_separ*y_separ)
      * (x_separ*shiftx + y_separ*shifty + om * y_separ * xbh_orb)
      * (x_separ*shiftx + y_separ*shifty + om * y_separ * xbh_orb) ;
      */

    double bpb = om2 * (x_orb * x_orb + y_orb * y_orb
			+ 2.*mhsr*(x_separ*x_separ+y_separ*y_separ)
			* y_separ * y_separ * xbh_orb * xbh_orb)
      + 2.*om*(shifty * x_orb - shiftx * y_orb
	       + 2.*mhsr*(x_separ*x_separ+y_separ*y_separ)
	       * (shiftx*x_separ + shifty*y_separ) * y_separ * xbh_orb)
      + shift2 + 2.*mhsr*(x_separ*x_separ+y_separ*y_separ)
      * (shiftx*x_separ + shifty*y_separ)
      * (shiftx*x_separ + shifty*y_separ) ;

    double dlngam0 =
      ( 0.5 * dp6sl2 * bpb
	+ p6sl2 * (0.5*dshift2
		   + om * (shifty - dshiftx*y_orb + dshifty*x_orb)
		   + om2 * x_orb
		   - mhsr * x_separ * (x_separ*shiftx+y_separ*shifty)
		   * (x_separ*shiftx+y_separ*shifty)
		   + 2. * mhsr * (x_separ*shiftx+y_separ*shifty)
		   * (y_separ*y_separ*shiftx - x_separ*y_separ*shifty
		      +(x_separ*x_separ+y_separ*y_separ)*(x_separ*dshiftx
							  +y_separ*dshifty))
		   + 2. * mhsr * om * y_separ * xbh_orb
		   * ( (y_separ*y_separ-2.*x_separ*x_separ)*shiftx
		       - 3. * x_separ * y_separ * shifty
		       +(x_separ*x_separ+y_separ*y_separ)*(x_separ*dshiftx
							   +y_separ*dshifty) )
		   - 3. * mhsr * om2 * x_separ * y_separ * y_separ * xbh_orb
		   * xbh_orb)
	) / (1 - p6sl2 * bpb) ;

    return dnulg - dlngam0 ;

}

double func_binbhns_orbit_is(double om, const Param& parorb) {

    double dnulg = parorb.get_double(0) ;
    double p6sl2 = parorb.get_double(1) ;
    double dp6sl2 = parorb.get_double(2) ;
    double shiftx = parorb.get_double(3) ;
    double shifty = parorb.get_double(4) ;
    double dshiftx = parorb.get_double(5) ;
    double dshifty = parorb.get_double(6) ;
    double shift2 = parorb.get_double(7) ;
    double dshift2 = parorb.get_double(8) ;
    double x_orb = parorb.get_double(9) ;
    double y_orb = parorb.get_double(10) ;

    double om2 = om * om ;

    double bpb = om2 * (x_orb * x_orb + y_orb * y_orb)
      + 2. * om * (shifty * x_orb - shiftx * y_orb) + shift2 ;

    double dlngam0 = ( 0.5 * dp6sl2 * bpb
		       + p6sl2 * (0.5*dshift2
				  + om *
				  (shifty - dshiftx*y_orb + dshifty*x_orb)
				  + om2 * x_orb)
		       ) / (1 - p6sl2 * bpb) ;

    return dnulg - dlngam0 ;

}
}
