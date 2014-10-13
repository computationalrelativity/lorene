/*
 *  Definition of Lorene class App_hor
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
 *
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

#ifndef __APP_HOR_H_ 
#define __APP_HOR_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/13 08:52:31  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2012/01/02 13:52:57  j_novak
 * New parameter 'verbose' to get less output if needed.
 *
 * Revision 1.4  2005/12/09 09:35:59  lm_lin
 *
 * Minor fix in the documentation.
 *
 * Revision 1.3  2005/12/07 11:11:30  lm_lin
 *
 * Add option to turn off screen output during iterations.
 *
 * Revision 1.2  2005/11/17 14:19:49  lm_lin
 *
 * Check the expansion function evaluated on the apparent horizon after the
 * iteration of the 2-surface converges.
 *
 * Revision 1.1  2005/10/13 08:51:14  j_novak
 * New stuff for apparent horizon finder. For the moment, there is only an
 * external function. A class should come soon...
 *
 * $Header$
 *
 */


// Headers Lorene
#include "metric.h" 

namespace Lorene {
/**
 * Class for apparent horizon (under heavy development)
 * 
 */

    // Function (Apparent horizon finder) 
    //-----------------------------------------------
	
         /**
	 * Apparent horizon finder.  \ingroup (star)
	 *
	 * Find the apparent horizon (AH) on a spatial slice for given
	 * 3-metric \f$\gamma_{ij}\f$ and extrinsic curvature \f$K_{ij}\f$.
	 * Method: We solve the apparent-horizon equation 
	 * \f$\Theta=0\f$ (where \f$\Theta\f$ is the expansion function
	 * of the outgoing null geodesic) as a generalized angular Poisson
	 * equation (\f$\Delta_{\theta\phi}u+\lambda u=\sigma\f$) in the form:
	 * 
	 * \f[ 
	 *   \Delta_{\theta\phi}h - 2h = \Psi^4 |DF| h^2\Theta + 
	 *   \Delta_{\theta\phi}h - 2h
	 * \f]
	 *
	 * where \f$h=h(\theta,\phi)\f$ (\c h) is the 2-surface of the AH, 
	 * \f$\Psi\f$ is the conformal factor, and 
	 * \f$|DF|:= (\gamma^{ij}D_iFD_jF)^{1/2} \f$ (where \f$F\f$ is the 
	 * level set function \f$ F:=r-h(\theta,\phi) \f$ and \f$D_i\f$ is the 
	 * covariant derivative w.r.t. \f$\gamma_{ij}\f$).	 
       	 * We solve the above equation iteratively with the right hand side 
         * of the equation as a source term. 
	 * 
	 * @param gamma : [input] the 3-metric \f$\gamma_{ij}\f$ w.r.t. which 
	 * the AH is to be found. 
	 * @param k_dd_in  : [input] the extrinsic curvature \f$K_{ij}\f$.
	 * @param Valeur h  : [output] the 2-surface of the apparent horizon
	 * @param Scalar ex_fcn : [output] the expansion function defined from 
	 *                        the level set function \f$F:= r - h(\theta,\phi)\f$
	 * @param a_axis : [input] the initial guess for \f$h\f$ is a triaxial 
	 * ellipsoidal surface with a_axis the axis along the x-axis of  
	 * the Cartesian grid
	 * @param b_axis : [input] axis along the y-axis (cf a_axis)
	 * @param c_axis : [input] axis along the z-axis (cf a_axis)
	 * @param bool print : [input] screen printout during iterations? (default: false)
	 * @param precis : [input] threshold in the relative difference between 
	 * the 2-surface \f$h\f$ of two consecutive steps to stop
	 * the iterative procedure (default value: 1.e-8)
	 * @param precis_exp : [input] maximum error of the expansion function evaluated on 
	 * the 2-surface \f$h\f$ (should be zero by definition) after the iteration
	 * is stopped (default value: 1.e-6)
	 * @param it_max : [input] maximum number of steps (default value: 200)
	 * @param it_relax : [input] step at which relaxation is switched on 
	 *                   (default value: 200)
	 * @param relax_fac : [input] relaxation factor (default value: 
	 *                    1 for no relaxation)
	 * @return bool ah_flag : a flag to indicate whether an apparent horizon is found
	 */
bool ah_finder(const Metric& gamma, const Sym_tensor& k_dd_in, Valeur& h, Scalar& ex_fcn,
	       double a_axis, double b_axis, double c_axis, bool verbose = true, 
	       bool print = false, double precis = 1.e-8, double precis_exp = 1.e-6,  
	       int it_max = 200, int it_relax = 200, double relax_fac = 1.) ;
		
}
#endif
