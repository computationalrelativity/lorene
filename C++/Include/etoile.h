/*
 *  Definition of Lorene classes Etoile
 *				 Etoile_bin
 *				 Etoile_rot
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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


#ifndef __ETOILE_H_ 
#define __ETOILE_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.27  2005/10/05 15:14:47  j_novak
 * Added a Param* as parameter of Etoile_rot::equilibrium
 *
 * Revision 1.26  2005/08/29 15:10:12  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.25  2005/01/18 22:35:51  k_taniguchi
 * Delete a pointer for Etoile::ray_eq(int kk).
 *
 * Revision 1.24  2005/01/18 20:34:14  k_taniguchi
 * Addition of Etoile::ray_eq(int kk).
 *
 * Revision 1.23  2005/01/17 20:39:32  k_taniguchi
 * Addition of Etoile::ray_eq_3pis2().
 *
 * Revision 1.22  2004/11/30 20:40:06  k_taniguchi
 * Addition of the method for calculating a spherical star with falloff
 * condition at the outer boundary.
 *
 * Revision 1.21  2004/05/10 10:18:33  f_limousin
 * Change to avoid warning in the compilation of Lorene
 *
 * Revision 1.20  2004/05/07 12:37:12  f_limousin
 * Add new member ssjm1_psi
 *
 * Revision 1.19  2004/04/08 16:42:31  f_limousin
 * Add a function velocity_potential with argument ssjm1_psi for the
 * case of strange stars.
 *
 * Revision 1.18  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.17  2003/10/24 12:24:41  k_taniguchi
 * Suppress the methods of update metric for NS-BH
 *
 * Revision 1.16  2003/10/21 11:44:43  k_taniguchi
 * Delete various things for the Bin_ns_bh project.
 * They are moved to et_bin_nsbh.h.
 *
 * Revision 1.15  2003/10/20 14:50:04  k_taniguchi
 * Addition of various things for the Bin_ns_bh project
 * which are related with the part of the neutron star.
 *
 * Revision 1.14  2003/10/20 13:11:03  k_taniguchi
 * Back to version 1.12
 *
 * Revision 1.13  2003/10/20 12:15:55  k_taniguchi
 * Addition of various things for the Bin_ns_bh project
 * which are related with the part of the neutron star.
 *
 * Revision 1.12  2003/06/20 14:13:16  f_limousin
 * Change to virtual the functions equilibrium_spher() and kinematics().
 *
 * Revision 1.11  2003/02/13 16:40:24  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.10  2003/02/04 16:20:35  f_limousin
 * Change to virtual the routine extrinsic_curvature
 *
 * Revision 1.9  2003/01/31 16:57:12  p_grandclement
 * addition of the member Cmp decouple used to compute the K_ij auto, once
 * the K_ij total is known
 *
 * Revision 1.8  2002/12/19 14:48:00  e_gourgoulhon
 *
 * Class Etoile_bin: added the new functions:
 * 	void update_metric(const Bhole& comp)
 *  	void update_metric_der_comp(const Bhole& comp)
 * to treat the case where the companion is a black hole
 *
 * Revision 1.7  2002/12/17 21:17:08  e_gourgoulhon
 * Class Etoile_bin: suppression of the member p_companion
 *                   as well as the associated functions set_companion
 *   		  and get_companion.
 *
 * Revision 1.6  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
 * Revision 1.5  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.4  2002/04/05 09:09:36  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.3  2002/01/11 14:09:34  j_novak
 * Added newtonian version for 2-fluid stars
 *
 * Revision 1.2  2001/12/06 15:11:43  jl_zdunik
 * Introduction of the new function f_eq() in the class Etoile_rot
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.61  2001/10/25  08:32:57  eric
 * Etoile_rot::display_poly passe de protected a public.
 *
 * Revision 2.60  2001/10/24  15:35:55  eric
 * Classe Etoile_rot: ajout de la fonction display_poly.
 *
 * Revision 2.59  2001/10/16  14:48:00  eric
 * La fonction Etoile_rot::get_omega() devient
 *   virtual Etoile_rot::get_omega_c()
 *  (retourne la valeur centrale de Omega).
 *
 * Revision 2.58  2001/10/11  09:24:00  eric
 * La fonction Etoile_rot::equilibrium est desormais virtuelle.
 *
 * Revision 2.57  2001/08/06  15:39:04  keisuke
 * Addition of a new argument to Etoile_bin::equilibrium and equil_regular.
 *
 * Revision 2.56  2001/06/25  12:52:33  eric
 * Classe Etoile_bin : ajout du membre p_companion et des fonctions
 *  associees set_companion() et get_companion().
 *
 * Revision 2.55  2001/06/13  14:11:55  eric
 * Modif commentaires (mise en conformite Doc++ 3.4.7)
 *
 * Revision 2.54  2001/03/26  09:29:59  jlz
 * Classe Etoile_rot: new members p_espec_isco and p_lspec_isco.
 *
 * Revision 2.53  2001/02/08  15:37:09  eric
 * *** empty log message ***
 *
 * Revision 2.52  2001/02/08  15:12:42  eric
 * Ajout de la fonction Etoile_rot::f_eccentric.
 *
 * Revision 2.51  2000/11/23  15:43:24  eric
 * Ajout de l'argument ent_limit a Etoile_rot::equilibrium.
 *
 * Revision 2.50  2000/11/19  18:51:13  eric
 * Etoile_rot: ajout de la fonction (static) lambda_grv2
 *
 * Revision 2.49  2000/11/18  23:17:32  eric
 * Ajout de l'argument ost a Etoile_rot::r_isco.
 *
 * Revision 2.48  2000/11/18  21:08:33  eric
 * Classe Etoile_rot: ajout des fonctions r_isco() et f_isco()
 *   ainsi que des membres associes p_r_isco et p_f_isco.
 *
 * Revision 2.47  2000/11/10  15:15:38  eric
 * Modif arguments Etoile_rot::equilibrium.
 *
 * Revision 2.46  2000/10/20  13:10:23  eric
 * Etoile_rot::equilibrium: ajout de l'argument nzadapt.
 *
 * Revision 2.45  2000/10/17  15:59:14  eric
 * Modif commentaires Etoile_rot::equilibrium
 *
 * Revision 2.44  2000/10/12  15:22:29  eric
 * Modif commentaires Etoile_rot.
 *
 * Revision 2.43  2000/10/12  10:20:12  eric
 * Ajout de la fonction Etoile_rot::fait_nphi().
 *
 * Revision 2.42  2000/09/22  15:50:13  keisuke
 * d_logn_auto_div devient desormais un membre de la classe Etoile
 * et non plus de la classe derivee Etoile_bin.
 *
 * Revision 2.41  2000/09/18  16:14:37  eric
 * Classe Etoile_rot: ajout du membre tkij et de la fonction
 *   extrinsic curvature().
 *
 * Revision 2.40  2000/09/07  14:31:09  keisuke
 * Ajout des membres logn_auto_regu et get_logn_auto_regu() a la classe Etoile.
 * Ajout des membres d_logn_auto_regu et get_d_logn_auto_regu()
 *  a la classe Etoile_bin.
 *
 * Revision 2.39  2000/09/07  10:25:50  keisuke
 * Ajout du membre get_logn_auto_div() a la classe Etoile.
 * Ajout du membre get_d_logn_auto_div() a la classe Etoile_bin.
 *
 * Revision 2.38  2000/08/31  11:25:24  eric
 * Classe Etoile_rot: ajout des membres tnphi et ak_car.
 *
 * Revision 2.37  2000/08/29  11:37:06  eric
 * Ajout des membres k_div et logn_auto_div a la classe Etoile.
 * Ajout du membre d_logn_auto_div a la classe Etoile_bin.
 *
 * Revision 2.36  2000/08/25  10:25:57  keisuke
 * Ajout de Etoile_bin::equil_regular.
 *
 * Revision 2.35  2000/08/18  14:01:07  eric
 * Modif commentaires.
 *
 * Revision 2.34  2000/08/17  12:38:30  eric
 * Modifs classe Etoile_rot : ajout des membres nuf, nuq, ssjm1_nuf et ssjm1_nuq
 * Modif arguments Etoile_rot::equilibrium.
 * .\
 *
 * Revision 2.33  2000/08/07  12:11:13  keisuke
 * Ajout de Etoile::equil_spher_regular.
 *
 * Revision 2.32  2000/07/21  12:02:55  eric
 * Suppression de Etoile_rot::relaxation.
 *
 * Revision 2.31  2000/07/20  15:32:28  eric
 * *** empty log message ***
 *
 * Revision 2.30  2000/07/06  09:39:12  eric
 * Ajout du membre p_xa_barycenter a Etoile_bin, ainsi que de la
 * fonction associee xa_barycenter().
 *
 * Revision 2.29  2000/05/25  13:47:38  eric
 * Modif Etoile_bin::equilibrium: ajout de l'argument thres_adapt
 *
 * Revision 2.28  2000/03/22  16:41:45  eric
 * Ajout de la sortie de nouvelles erreurs dans Etoile_bin::equilibrium.
 *
 * Revision 2.27  2000/03/22  12:54:44  eric
 * Nouveau prototype d'Etoile_bin::velocity_potential : l'erreur est
 * retournee en double.
 * Nouveau prototype d'Etoile_bin::equilibrium : diff_ent est remplace
 *  par le Tbl diff.
 *
 * Revision 2.26  2000/03/15  11:04:15  eric
 * Ajout des fonctions Etoile_bin::set_w_shift() et Etoile_bin::set_khi_shift()
 * Amelioration des commentaires.
 *
 * Revision 2.25  2000/03/08  12:12:49  eric
 * Ajout de la fonction Etoile_bin::is_irrotational().
 *
 * Revision 2.24  2000/03/07  14:48:01  eric
 * Ajout de la fonction Etoile_bin::extrinsic_curvature().
 *
 * Revision 2.23  2000/02/21  13:57:57  eric
 * classe Etoile_bin: suppression du membre psi: remplacement par psi0.
 *
 * Revision 2.22  2000/02/17  16:51:22  eric
 * Ajout de l'argument diff_ent dans Etoile_bin::equilibrium.
 *
 * Revision 2.21  2000/02/17  15:29:38  eric
 * Ajout de la fonction Etoile_bin::relaxation.
 *
 * Revision 2.20  2000/02/17  13:54:21  eric
 * Ajout de la fonction Etoile_bin::velocity_potential.
 *
 * Revision 2.19  2000/02/16  15:05:14  eric
 * Ajout des membres w_shift et khi_shift.
 * (sauves dans les fichiers a la place de shift_auto).
 * Ajout de la fonction Etoile_bin::fait_shift_auto.
 *
 * Revision 2.18  2000/02/16  13:47:02  eric
 * Classe Etoile_bin: ajout des membres ssjm1_khi et ssjm1_wshift.
 *
 * Revision 2.17  2000/02/16  11:54:13  eric
 * Classe Etoile_bin : ajout des membres ssjm1_logn et ssjm1_beta.
 *
 * Revision 2.16  2000/02/15  15:40:07  eric
 * Ajout de Etoile_bin::equilibrium.
 *
 * Revision 2.15  2000/02/12  18:40:15  eric
 * Modif commentaires.
 *
 * Revision 2.14  2000/02/12  14:44:26  eric
 * Ajout des fonctions Etoile_bin::set_logn_comp et set_pot_centri.
 *
 * Revision 2.13  2000/02/10  20:22:25  eric
 * Modif commentaires.
 *
 * Revision 2.12  2000/02/10  16:11:24  eric
 * Classe Etoile_bin : ajout des accesseurs get_psi, etc...
 *                     ajout de la fonction fait_d_psi
 *
 * Revision 2.11  2000/02/08  19:28:29  eric
 * La fonction Etoile_bin::scal_prod est rebaptisee Etoile_bin::sprod
 *
 * Revision 2.10  2000/02/04  17:15:15  eric
 * Classe Etoile_bin: ajout du membre ref_triad.
 *
 * Revision 2.9  2000/02/04  16:36:48  eric
 * Ajout des fonctions update_metric* et kinematics.
 *
 * Revision 2.8  2000/02/02  10:12:37  eric
 * Ajout des fonctions de lecture/ecriture mp, nzet, eos, etc...
 *
 * Revision 2.7  2000/02/01  15:59:43  eric
 * Ajout de la fonction Etoile_bin::scal_prod.
 *
 * Revision 2.6  2000/01/31  15:56:45  eric
 * Introduction de la classe derivee Etoile_bin.
 *
 * Revision 2.5  2000/01/28  17:17:45  eric
 * Ajout des fonctions de calcul des quantites globales.
 *
 * Revision 2.4  2000/01/27  16:46:59  eric
 * Ajout des fonctions get_ent(), etc....
 *
 * Revision 2.3  2000/01/24  17:19:48  eric
 * Modif commentaires.
 *
 * Revision 2.2  2000/01/24  17:13:04  eric
 * Le mapping mp n'est plus constant.
 * Ajout de la fonction equilibrium_spher.
 *
 * Revision 2.1  2000/01/24  13:37:19  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/01/20  17:04:33  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tenseur.h"
class Eos ;
class Bhole ;


			//---------------------------//
			//    base class Etoile      //
			//---------------------------//

/**
 * Base class for stars. \ingroup (star)
 * 
 * An \c Etoile  is constructed upon (i) a mapping 
 * (derived class of \c Map ), the center of which defines the center of the 
 * star, and (ii) an equation of state (derived class of \c Eos ).  
 * It contains tensor fields (class \c Tenseur ) which describle the
 * hydrodynamical quantities as well as the gravitational field (spacetime
 * metric). 
 * 
 * According to the 3+1 formalism, the spacetime metric is written
 * \f[ \label{eetoilemetrique}
 *   ds^2 = - (N^2 - N_i N^i) dt^2 - 2 N_i \,  dt\,  dx^i
 *	    + A^2  \,  {\tilde h}_{ij} \,  dx^i dx^j
 * \f]
 * where \f${\tilde h}_{ij}\f$ is a 3-metric, the exact form of which is specified
 * in the derived classes of \c Etoile . The base class \c Etoile  by
 * itself provides only storage for the lapse function \e N  (member \c nnn ), 
 * the shift vector \f$N^i\f$ (member \c shift ) and the conformal factor
 * \f$A^2\f$ (member \c a_car ). 
 * 
 * The 3+1 formalism introduces two kinds of priviledged observers: the
 * fluid comoving observer and the Eulerian observer, whose 4-velocity
 * is the unit future directed normal to the \e t = const hypersurfaces. 
 * The hydrodynamical quantities measured by the fluid observer correspond
 * to the members \c ent , \c nbar , \c ener , and \c press . 
 * The hydrodynamical quantities measured by the Eulerian observer correspond
 * to the members \c ener_euler , \c s_euler , \c gam_euler , and 
 * \c u_euler .  
 * 
 * A star of class \c Etoile  can be either relativistic or Newtonian, 
 * depending on the boolean indicator \c relativistic . For a Newtonian
 * star, the metric coefficients \e N  and \e A  are set to 1,  and \f$N^i\f$ is
 * set to zero; the only relevant gravitational quantity in this case is
 * \c logn_auto  which represents the (regular part of) the
 * Newtonian gravitational potential
 * generated by the star. 
 *
 */
class Etoile {

    // Data : 
    // -----
    protected:
	Map& mp ;	    ///< Mapping associated with the star

	/** Number of domains of \c *mp  occupied by the star
	 * 
	 */
	int nzet ;
	
	/** Indicator of relativity: \c true  for a relativistic star,
	 *	\c false  for a Newtonian one. 
	 */
	bool relativistic ;
	
	/** \f$1/c^2\f$ : \c unsurc2=1 for a relativistic star,
	 *  0 for a Newtonian one. 
	 */
	double unsurc2 ; 	     
    
	/** Index of regularity of the gravitational potential \c logn_auto .
	 *  If \c k_div=0 , \c logn_auto  contains the total potential
	 *  generated principaly by the star, otherwise it should be
	 *  supplemented by \c logn_auto_div . 
	 */
	int k_div ; 
    
	const Eos& eos ;   ///< Equation of state of the stellar matter
    
	// Fluid quantities with respect to the fluid frame
	// ------------------------------------------------
	/// Log-enthalpy (relativistic case) or specific enthalpy (Newtonian case)
	Tenseur ent ;	  
	
	Tenseur nbar ; 	   ///< Baryon density in the fluid frame
	Tenseur ener ;	   ///< Total energy density in the fluid frame
	Tenseur press ;	   ///< Fluid pressure
    
	// Fluid quantities with respect to the Eulerian frame
	// ---------------------------------------------------
	Tenseur ener_euler ; ///< Total energy density in the Eulerian frame

	/// Trace of the stress tensor in the Eulerian frame
	Tenseur s_euler ;   
	
	/// Lorentz factor between the fluid and Eulerian observers 
	Tenseur gam_euler ; 
	
	/// Fluid 3-velocity with respect to the Eulerian observer
	Tenseur u_euler ; 
	
	// Metric potentials
	// -----------------
	
	/** Total of the logarithm of the part of the lapse \e N  
	 *   generated principaly by the star. In the Newtonian case, 
	 *   this is the Newtonian gravitational potential
	 *   (in units of \f$c^2\f$). 
	 */
	Tenseur logn_auto ;

	/** Regular part of the logarithm of the part of the lapse \e N  
	 *   generated principaly by the star. In the Newtonian case, 
	 *   this is the Newtonian gravitational potential
	 *   (in units of \f$c^2\f$). 
	 */
	Tenseur logn_auto_regu ;
	
	/** Divergent part (if \c k_div!=0 ) 
	 *  of the logarithm of the part of the lapse \e N  
	 *   generated principaly by the star. 
	 */
	Tenseur logn_auto_div ; 
	
	/** Gradient of \c logn_auto_div  (if \c k_div!=0 )
	 */
	Tenseur d_logn_auto_div ; 

	/** Logarithm of the part of the product \e AN  generated principaly by
	 *   by the star
	 */
	Tenseur beta_auto ; 
	
	/// Total lapse function 
	Tenseur nnn ; 
	
	/// Total shift vector
	Tenseur shift ;
	
	/// Total conformal factor \f$A^2\f$
	Tenseur a_car ; 
	
    // Derived data : 
    // ------------
    protected:
	/// Coordinate radius at \f$\phi=0\f$, \f$\theta=\pi/2\f$. 
	mutable double* p_ray_eq ; 

	/// Coordinate radius at \f$\phi=\pi/2\f$, \f$\theta=\pi/2\f$. 
	mutable double* p_ray_eq_pis2 ;

	/// Coordinate radius at \f$\phi=\pi\f$, \f$\theta=\pi/2\f$. 
	mutable double* p_ray_eq_pi ;
	
	/// Coordinate radius at \f$\phi=3\pi/2\f$, \f$\theta=\pi/2\f$. 
	mutable double* p_ray_eq_3pis2 ;

	/// Coordinate radius at \f$\theta=0\f$. 
	mutable double* p_ray_pole ;

	/** Description of the stellar surface: 2-D \c Itbl  containing the 
	 *	values of the domain index \e l  on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$
	 */
	mutable Itbl* p_l_surf ; 
	
	/** Description of the stellar surface: 2-D \c Tbl  containing the 
	 *	values of the radial coordinate \f$\xi\f$ on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$
	 */
	mutable Tbl* p_xi_surf ; 
	
	mutable double* p_mass_b ;	///< Baryon mass
	mutable double* p_mass_g ;	///< Gravitational mass


    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be \c true  for a relativistic
	 *			star,  \c false  for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 * 
	 */
	Etoile(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i) ;			
	
	
	Etoile(const Etoile& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	\c sauve )
	 */
	Etoile(Map& mp_i, const Eos& eos_i, FILE* fich) ;    		

	virtual ~Etoile() ;			///< Destructor

	
    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0  all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to \c ETATNONDEF  (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Etoile 
	void operator=(const Etoile&) ;	
	
	/// Read/write of the mapping
	Map& set_mp() {return mp; } ; 

	/// Assignment of the enthalpy field.
	void set_enthalpy(const Cmp& ) ; 
	
	/** Computes the proper baryon and energy density, as well as
	 *  pressure from the enthalpy.
	 */
	void equation_of_state() ; 
	
	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame (\c nbar , \c ener 
	 *  and \c press ).
	 */
	virtual void hydro_euler() ; 
	
	/** Computes a spherical static configuration. 
	 * 
	 *  @param ent_c [input] central value of the enthalpy
	 *  @param precis [input] threshold in the relative difference between 
	 *	the enthalpy fields of two consecutive steps
	 *	to stop the iterative procedure (default value: 1.e-14)
	 */
	virtual void equilibrium_spher(double ent_c, double precis = 1.e-14) ; 

	/** Computes a spherical static configuration. 
	 *  The sources for Poisson equations are regularized
	 *  by extracting analytical diverging parts.
	 * 
	 *  @param ent_c [input] central value of the enthalpy
	 *  @param precis [input] threshold in the relative difference between 
	 *	the enthalpy fields of two consecutive steps
	 *	to stop the iterative procedure (default value: 1.e-14)
	 */
	void equil_spher_regular(double ent_c, double precis = 1.e-14) ; 

	/** Computes a spherical static configuration with the outer
	 *  boundary condition at a finite radius
	 * 
	 *  @param ent_c [input] central value of the enthalpy
	 *  @param precis [input] threshold in the relative difference between 
	 *	the enthalpy fields of two consecutive steps
	 *	to stop the iterative procedure (default value: 1.e-14)
	 */
	virtual void equil_spher_falloff(double ent_c,
					 double precis = 1.e-14) ;

    // Accessors
    // ---------
    public:
	/// Returns the mapping
	const Map& get_mp() const {return mp; } ; 

	/// Returns the number of domains occupied by the star
	int get_nzet() const {return nzet; } ; 

	/** Returns \c true  for a relativistic star, \c false  for 
	 *  a Newtonian one
	 */
	bool is_relativistic() const {return relativistic; } ; 

	/// Returns the equation of state
	const Eos& get_eos() const {return eos; } ; 

	/// Returns the enthalpy field 
	const Tenseur& get_ent() const {return ent;} ;

	/// Returns the proper baryon density
	const Tenseur& get_nbar() const {return nbar;} ;

	/// Returns the proper total energy density
	const Tenseur& get_ener() const {return ener;} ;

	/// Returns the fluid pressure
	const Tenseur& get_press() const {return press;} ;

	/// Returns the total energy density with respect to the Eulerian observer
	const Tenseur& get_ener_euler() const {return ener_euler;} ;

	/// Returns the trace of the stress tensor in the Eulerian frame
	const Tenseur& get_s_euler() const {return s_euler;} ;

	/// Returns the Lorentz factor between the fluid and Eulerian observers
	const Tenseur& get_gam_euler() const {return gam_euler;} ;

	/// Returns the fluid 3-velocity with respect to the Eulerian observer
	const Tenseur& get_u_euler() const {return u_euler;} ;

	/** Returns the logarithm of the part of the lapse \e N  generated 
	 *   principaly by the star.
	 *   In the Newtonian case, this is the Newtonian
	 *   gravitational potential (in units of \f$c^2\f$). 
	 */
	const Tenseur& get_logn_auto() const {return logn_auto;} ;

	/** Returns the regular part of the logarithm of the part of
	 *   the lapse \e N  generated principaly by the star.
	 *   In the Newtonian case, this is the Newtonian
	 *   gravitational potential (in units of \f$c^2\f$). 
	 */
	const Tenseur& get_logn_auto_regu() const {return logn_auto_regu;} ;

	/** Returns the divergent part of the logarithm of the part of
	 *   the lapse \e N  generated principaly by the star.
	 *   In the Newtonian case, this is the diverging part of
	 *   the Newtonian gravitational potential (in units of \f$c^2\f$). 
	 */
	const Tenseur& get_logn_auto_div() const {return logn_auto_div;} ;

	/** Returns the gradient of \c logn_auto_div 
	 */
	const Tenseur& get_d_logn_auto_div() const {return d_logn_auto_div;} ;

	/** Returns the logarithm of the part of the product \e AN  generated 
	 *  principaly by the star.
	 */
	const Tenseur& get_beta_auto() const {return beta_auto;} ;

	/// Returns the total lapse function \e N 
	const Tenseur& get_nnn() const {return nnn;} ;

	/// Returns the total shift vector \f$N^i\f$
	const Tenseur& get_shift() const {return shift;} ;

	/// Returns the total conformal factor \f$A^2\f$
	const Tenseur& get_a_car() const {return a_car;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Etoile& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	/// Coordinate radius at \f$\phi=0\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq() const ; 
	
	/// Coordinate radius at \f$\phi=\pi/2\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq_pis2() const ; 
	
	/// Coordinate radius at \f$\phi=\pi\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq_pi() const ; 
	
	/// Coordinate radius at \f$\phi=3\pi/2\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq_3pis2() const ;

	/// Coordinate radius at \f$\theta=0\f$ [r_unit]. 
	double ray_pole() const ; 
    
	/// Coordinate radius at \f$\phi=2k\pi/np\f$, \f$\theta=\pi/2\f$ [r_unit].
	double ray_eq(int kk) const ;

	/** Description of the stellar surface: returns a 2-D \c Itbl  
	 *	containing the 
	 *	values of the domain index \e l  on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$.
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member \c ent ) vanishes.
	 */
	virtual const Itbl& l_surf() const ; 
	
	/** Description of the stellar surface: returns a 2-D \c Tbl  
	 *	containing the values of the radial coordinate \f$\xi\f$ 
	 *	on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$. 
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member \c ent ) vanishes.
	 */
	const Tbl& xi_surf() const ; 

	/// Baryon mass
    	virtual double mass_b() const ;
	
	/// Gravitational mass
    	virtual double mass_g() const ;
	
};
ostream& operator<<(ostream& , const Etoile& ) ;	


			//---------------------------//
			//    class Etoile_bin       //
			//---------------------------//

/**
 * Class for stars in binary system. \ingroup (star)
 *
 * This class implements the formalism for corotating or irrotational
 * systems presented in Bonazzola, Gourgoulhon \& Marck \a Phys. \a Rev. \a Lett.
 * \b 82 , 892 (1999). In particular, the conformal 3-metric 
 * \f$\tilde h_{ij}\f$ introduced in Eq.~(\ref{eetoilemetrique}) is flat. 
 *
 * An \c Etoile_bin  can be construted in two states, represented by
 * the \c bool  member \c irrotational : (i) irrotational
 * (i.e. the fluid motion is irrotational) or (ii) rigidly corotating 
 * with respect to the orbital motion (synchronized binary). 
 *
 */
class Etoile_bin : public Etoile {

    // Data : 
    // -----
    protected:
	/** \c true  for an irrotational star, \c false  for a
	 *  corotating one
	 */
	bool irrotational ; 
	
	/** Reference triad ("absolute frame"), with
	 *  respect to which the components of all the member \c Tenseur 's
	 *  are defined, except for \c w_shift  and \c ssjm1_wshift .
	 */
	const Base_vect& ref_triad ; 
	
	/** Scalar potential \f$\Psi_0\f$ of the non-translational part of
	 *  fluid 4-velocity (in the irrotational case)
	 */
	Tenseur psi0 ; 

	/** Gradient of \f$\Psi\f$ (in the irrotational case)
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur d_psi ; 
	
	/** Spatial projection of the fluid 3-velocity with respect to  
	 *  the co-orbiting observer. 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur wit_w ; 
	
	/** Logarithm of the Lorentz factor between the fluid and 
	 *  the co-orbiting observer.
	 */
	Tenseur loggam ; 

	/** Part of the lapse logarithm (gravitational potential at the
	 *  Newtonian limit) generated principaly by the companion star. 
	 */
	Tenseur logn_comp ; 
	
	/** Gradient of \c logn_auto 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur d_logn_auto ; 

	/** Gradient of \c logn_auto_regu 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur d_logn_auto_regu ; 

	/** Gradient of \c logn_comp 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur d_logn_comp ; 

	/** Part of the logarithm of \e AN  generated principaly by the 
	 *  companion star. 
	 */
	Tenseur beta_comp ; 
	
	/** Gradient of \c beta_auto 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur d_beta_auto ; 

	/** Gradient of \c beta_comp 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur d_beta_comp ; 

	/** Part of the shift vector \f$N^i\f$ generated principaly by the star.
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur shift_auto ; 
	
	/** Part of the shift vector \f$N^i\f$ generated principaly by the 
	 *  companion star. 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur shift_comp ; 
	
	/** Vector \f$W^i\f$ used in the decomposition of \c shift_auto ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101, 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	Tenseur w_shift ; 
	
	/** Scalar \f$\chi\f$ used in the decomposition of \c shift_auto ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 */
	Tenseur khi_shift ; 

	/** Part of the extrinsic curvature tensor 
	 *  \f$\tilde K^{ij} = A^2 K^{ij}\f$
	 *  generated by \c shift_auto . 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur_sym tkij_auto ;
	
	/** Part of the extrinsic curvature tensor 
	 *  \f$\tilde K^{ij}  = A^2 K^{ij}\f$
	 *  generated by \c shift_comp . 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur_sym tkij_comp ;
	
	/** Part of the scalar \f$A^2 K_{ij} K^{ij}\f$
	 *  generated by \c shift_auto , i.e. 
	 *    \f$A^2 K_{ij}^{\rm auto} K^{ij}_{\rm auto}\f$ 
	 */
	Tenseur akcar_auto ;
	
	/** Part of the scalar \f$A^2 K_{ij} K^{ij}\f$
	 *  generated by \c shift_auto  and \c shift_comp , i.e. 
	 *    \f$A^2 K_{ij}^{\rm auto} K^{ij}_{\rm comp}\f$ 
	 */
	Tenseur akcar_comp ;
	
	/** 3-vector shift, divided by \e N , of the rotating coordinates,
	 *  \f$B^i/N\f$. 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	Tenseur bsn ; 
	
	/// Centrifugal potential
	Tenseur pot_centri ; 	
	
	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c logn_auto  by means of
	 *  \c Map_et::poisson .
	 */
	Cmp ssjm1_logn ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c beta_auto  by means of
	 *  \c Map_et::poisson .
	 */
	Cmp ssjm1_beta ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for the scalar \f$\chi\f$ by means of
	 *  \c Map_et::poisson . 
	 *  \f$\chi\f$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector \f$N^i\f$
	 */
	 Cmp ssjm1_khi ; 
	 
	/** Effective source at the previous step for the resolution of 
	 *  the vector Poisson equation for \f$W^i\f$ by means of
	 *  \c Map_et::poisson . 
	 *  \f$W^i\f$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector \f$N^i\f$
	 *  (Components with respect to the Cartesian triad associated with 
	 *   the mapping \c mp )
	 */
	 Tenseur ssjm1_wshift ; 
	 
	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for the scalar \f$\psi\f$ by means of
	 *  \c Map_et::poisson_interne . 
	 */
	 Cmp ssjm1_psi ; 

	 /**
	  * Function used to construct the part of \f$K^{ij}\f$ generated by 
	  * the star  from the total \f$K^{ij}\f$. Only used for a binary system
	  * where the other member is a black hole.
	  * 
	  * Mainly this \c Cmp  is 1 around the hole and 0 around the companion
	  * and the sum of \c decouple for the hole and his companion is 1 
	  * everywhere.
	  */
	 Cmp decouple ;
	
    // Derived data : 
    // ------------
    protected:
	/// Absolute coordinate X of the barycenter of the baryon density
	mutable double* p_xa_barycenter ; 
	 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be \c true  for a relativistic
	 *			star,  \c false  for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 * @param irrot should be \c true  for an irrotational star, 
	 *		    \c false  for a corotating one
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    \c Tenseur 's are defined, except for \c w_shift 
	 *	    and \c ssjm1_wshift  whose components are defined
	 *	    with respect to the mapping \c mp  Cartesian triad. 
	 */
	Etoile_bin(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		   bool irrot, const Base_vect& ref_triad_i) ;			
	
	
	Etoile_bin(const Etoile_bin& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    \c Tenseur 's are defined, except for \c w_shift 
	 *	    and \c ssjm1_wshift  whose components are defined
	 *	    with respect to the mapping \c mp  Cartesian triad. 
	 * @param fich	input file (must have been created by the function
	 *	\c sauve )
	 */
	Etoile_bin(Map& mp_i, const Eos& eos_i, const Base_vect& ref_triad_i, 
		   FILE* fich) ;    		

	virtual ~Etoile_bin() ;			///< Destructor


    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0  all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to \c ETATNONDEF  (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Etoile_bin 
	void operator=(const Etoile_bin& ) ;	
	
	/** Read/write the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the companion star. 
	 */
	Tenseur& set_logn_comp() ;

	/// Read/write the centrifugal potential
	Tenseur& set_pot_centri() ;
	
	/// Read/write of \c w_shift 
	Tenseur& set_w_shift() ;
	
	/// Read/write of \c khi_shift 
	Tenseur& set_khi_shift() ;
		
    // Accessors
    // ---------
    public:
	/** Returns \c true  for an irrotational motion, \c false  for 
	 *  a corotating one. 
	 */
	bool is_irrotational() const {return irrotational; } ; 

	/// Returns the non-translational part of the velocity potential
	const Tenseur& get_psi0() const {return psi0;} ;

	/** Returns the gradient of the velocity potential 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_d_psi() const {return d_psi;} ;

	/** Returns the spatial projection of the fluid 3-velocity with 
	 *  respect to the co-orbiting observer. 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_wit_w() const {return wit_w;} ;

	/** Returns the logarithm of the Lorentz factor between the fluid and 
	 *  the co-orbiting observer.
	 */
	const Tenseur& get_loggam() const {return loggam;} ;

	/** Returns the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the companion star. 
	 */
	const Tenseur& get_logn_comp() const {return logn_comp;} ;

	/** Returns the gradient of \c logn_auto 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_d_logn_auto() const {return d_logn_auto;} ;

	/** Returns the gradient of \c logn_auto_regu 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_d_logn_auto_regu() const {return d_logn_auto_regu;} ;

	/** Returns the gradient of \c logn_comp 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_d_logn_comp() const {return d_logn_comp;} ;

	/** Returns the part of the logarithm of \e AN  generated principaly 
	 *  by the companion star. 
	 */
	const Tenseur& get_beta_comp() const {return beta_comp;} ;

	/** Returns the gradient of \c beta_auto 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_d_beta_auto() const {return d_beta_auto;} ;

	/** Returns the gradient of \c beta_comp  
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_d_beta_comp() const {return d_beta_comp;} ;

	/** Returns the part of the shift vector \f$N^i\f$ generated principaly 
	 * by the star.
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_shift_auto() const {return shift_auto;} ;

	/** Returns the part of the shift vector \f$N^i\f$ generated principaly 
	 *   by the companion star. 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_shift_comp() const {return shift_comp;} ;

	/** Returns the vector \f$W^i\f$ used in the decomposition of 
	 *  \c shift_auto ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	const Tenseur& get_w_shift() const {return w_shift;} ; 
	
	/** Returns the scalar \f$\chi\f$ used in the decomposition of 
	 *  \c shift_auto  
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	const Tenseur& get_khi_shift() const {return khi_shift;} ; 

	/** Returns the part of the extrinsic curvature tensor 
	 *  \f$\tilde K^{ij} = A^2 K^{ij}\f$ generated by \c shift_auto . 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur_sym& get_tkij_auto() const {return tkij_auto;} ;

	/** Returns the part of the extrinsic curvature tensor 
	 *  \f$\tilde K^{ij} = A^2 K^{ij}\f$ generated by \c shift_comp . 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur_sym& get_tkij_comp() const {return tkij_comp;} ;

	/** Returns the part of the scalar \f$A^2 K_{ij} K^{ij}\f$
	 *  generated by \c shift_auto , i.e. 
	 *    \f$A^2 K_{ij}^{\rm auto} K^{ij}_{\rm auto}\f$ 
	 */
	const Tenseur& get_akcar_auto() const {return akcar_auto;} ;

	/** Returns the part of the scalar \f$A^2 K_{ij} K^{ij}\f$
	 *  generated by \c shift_auto  and \c shift_comp , i.e. 
	 *    \f$A^2 K_{ij}^{\rm auto} K^{ij}_{\rm comp} \f$ 
	 */
	const Tenseur& get_akcar_comp() const {return akcar_comp;} ;

	/** Returns the shift vector, divided by \e N , of the rotating 
	 *   coordinates, \f$B^i/N\f$. 
	 *  (Cartesian components with respect to \c ref_triad )
	 */
	const Tenseur& get_bsn() const {return bsn;} ;

	/// Returns the centrifugal potential
	const Tenseur& get_pot_centri() const {return pot_centri;} ;
	/**
	 * Returns the function used to construct \c tkij_auto  
	 from \c tkij_tot .
	*/
	const Cmp get_decouple() const {return decouple ;}

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    ///< Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	/// Baryon mass
    	virtual double mass_b() const ;
	
	/// Gravitational mass
    	virtual double mass_g() const ;
	
	/** Absolute coordinate X of the barycenter of the baryon density, 
	 *  defined according to the formula
	 *  \f[
	 *    X_G := \int A^3 \Gamma_{\rm n} \,  n \,  X \, d^3x \ ,  
	 *  \f]
	 *  where \f$\Gamma_{\rm n}\f$ is the Lorentz factor between the fluid 
	 *  and Eulerian observers.
	 */
    	virtual double xa_barycenter() const ;
	

    // Computational routines
    // ----------------------
    public: 
	/** Performs the scalar product of two tensors by contracting
	 *  the last index of \c t1  with the first index of \c t2 .
	 *  Both indices are supposed to be contravariant, so that a 
	 *  multiplication by \f$A^2\f$ is performed to lower one index. 
	 *  For instance, for two vectors \f$V^i\f$ and \f$W^i\f$, this function
	 *  returns the scalar \f$h_{ij} V^i W^j = A^2 f_{ij} V^i W^j\f$.  
	 */
	virtual Tenseur sprod(const Tenseur& t1, const Tenseur& t2) const ; 

	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame, as well as 
	 *  \c wit_w  and \c loggam .  
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c ent , \c ener , \c press , \c a_car  and \c bsn ,  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  \c gam_euler , \c u_euler , \c ener_euler , \c s_euler , 
	 *  \c wit_w  and \c loggam . 
	 * 
	 */
	virtual void hydro_euler() ; 
	
	/** Computes metric coefficients from known potentials,
	 * when the companion is another star.
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c logn_auto ,  \c beta_auto , \c shift_auto ,
	 *  \c comp.logn_auto ,  \c comp.beta_auto ,
	 *  \c comp.shift_auto 
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  \c logn_comp , \c beta_comp , \c shift_comp ,
	 *  \c nnn ,  \c a_car ,  \c shift ,
	 *  \c d_logn_auto , \c d_beta_auto , \c tkij_auto ,
	 *  \c akcar_auto .
	 *
	 *  @param comp companion star.
	 *
	 */
	void update_metric(const Etoile_bin& comp) ;

	/** Same as \c update_metric(const Etoile_bin\& \c ) but with
	 *  relaxation.
	 *
	 *  @param comp companion star.
	 *  @param star_prev previous value of the star. 
	 *  @param relax relaxation parameter. 
	 * 
	 */
	void update_metric(const Etoile_bin& comp, const Etoile_bin& star_prev, 
			   double relax) ; 
	
	/** Computes the derivative of metric functions related to the
	 *  companion star.
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c comp.d_logn_auto ,  \c comp.d_beta_auto ,
	 *  \c comp.tkij_auto 
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  \c d_logn_comp , \c d_beta_comp , \c tkij_comp ,
	 *  \c akcar_comp .
	 *
	 *  @param comp companion star.
	 *
	 */
	 void update_metric_der_comp(const Etoile_bin& comp) ;

	/** Computes the quantities \c bsn  and \c pot_centri .
	 * 
	 *  The calculation is performed starting from the quantities
	 *  \c nnn , \c shift ,  \c a_car ,  
	 *  which are supposed to be up to date.  
	 * 
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 *  @param x_axe  absolute X coordinate of the rotation axis
	 */
	virtual void kinematics(double omega, double x_axe) ; 
	
	/** Computes the gradient of the total velocity potential \f$\psi\f$. 
	 * 
	 */
	void fait_d_psi() ; 
	
	/** Computes \c shift_auto  from \c w_shift  and \c khi_shift 
	 *  according to Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 */
	void fait_shift_auto() ; 
	
	/** Computes \c tkij_auto  and \c akcar_auto  from 
	 *  \c shift_auto , \c nnn  and \c a_car .
	 */
	virtual void extrinsic_curvature() ; 
		
	
	/** Computes an equilibrium configuration.
	 * 
	 *  The values of \c logn_comp , \c beta_comp , \c pot_centri 
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param fact [input]    1-D \c Tbl  for the input
	 *                          of some factors : 
	 *          \li \c fact(0)  : A resizing factor for the first shell
	 *  @param diff [output]   1-D \c Tbl  for the storage of some
	 *			    error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error returned by the routine
	 *				\c Etoile_bin::velocity_potential   
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c logn_auto    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c beta_auto    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (x comp.)   
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (y comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (z comp.)   
	 */
	void equilibrium(double ent_c, int mermax, int mermax_poisson, 
			 double relax_poisson, int mermax_potvit, 
			 double relax_potvit, double thres_adapt, 
			 const Tbl& fact, Tbl& diff) ;

	/** Computes an equilibrium configuration by regularizing
	 *  the diverging source.
	 * 
	 *  The values of \c logn_comp , \c beta_comp , \c pot_centri 
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map_radial::poisson_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map_radial::poisson_compact
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param fact [input]    1-D \c Tbl  for the input
	 *                          of some factors : 
	 *          \li \c fact(0)  : A resizing factor for the first shell
	 *  @param diff [output]   1-D \c Tbl  for the storage of some
	 *			    error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error returned by the routine
	 *				\c Etoile_bin::velocity_potential   
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c logn_auto    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c beta_auto    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (x comp.)   
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (y comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift_auto  (z comp.)   
	 */
	void equil_regular(double ent_c, int mermax, int mermax_poisson, 
			   double relax_poisson, int mermax_potvit, 
			   double relax_potvit, double thres_adapt, 
			   const Tbl& fact, Tbl& diff) ;

	/** Computes the non-translational part of the velocity scalar potential
	 *  \f$\psi0\f$ by solving the continuity equation.
	 *  
	 *  @param mermax  [input] Maximum number of steps in the iteration
	 *  @param precis  [input] Required precision: the iteration will
	 *			   be stopped when the relative difference
	 *			   on \f$\psi0\f$ between two successive steps
	 *			   is lower than \c precis .
	 *  @param relax   [input] Relaxation factor.  
	 *
	 *  @return Relative error of the resolution obtained by comparing
	 *	    the operator acting on the solution with the source.
	 */
	double velocity_potential(int mermax, double precis, double relax) ;
	
	/** Performs a relaxation on \c ent , \c logn_auto ,
	 *  \c beta_auto  and \c shift_auto . 
	 * 
	 *  @param star_prev   [input] star at the previous step.
	 *  @param relax_ent   [input] Relaxation factor for \c ent  
	 *  @param relax_met   [input] Relaxation factor for \c logn_auto ,
	 *			       \c beta_auto , \c shift_auto , 
	 *			       only if \c (mer \% fmer_met == 0) .
	 *  @param mer	       [input] Step number
	 *  @param fmer_met    [input] Step interval between metric updates
	 */
	void relaxation(const Etoile_bin& star_prev, double relax_ent, 
			double relax_met, int mer, int fmer_met) ;
	
	friend class Bin_ns_bh ; ///< Friend class Bin_ns_bh
};



			//---------------------------//
			//    class Etoile_rot       //
			//---------------------------//

/**
 * Class for isolated rotating stars. \ingroup (star)
 * 
 * The metric is
 * \f[
 *   ds^2 = - N^2 dt^2 + A^2 (dr^2 + r^2 d\theta^2)
 *		       + B^2 r^2 \sin^2\theta (d\varphi - N^\varphi dt)^2
 * \f]
 *
 *
 */
class Etoile_rot : public Etoile {

    // Data : 
    // -----
    protected:
	double omega ;	    ///< Rotation angular velocity (\c [f_unit] ) 

	/// Metric factor \e B 
	Tenseur bbb ; 

	/// Square of the metric factor \e B 
	Tenseur b_car ; 

	/// Metric coefficient \f$N^\varphi\f$
	Tenseur nphi ; 

	/** Component \f$\tilde N^\varphi = N^\varphi r\sin\theta\f$ of the
	 *  shift vector
	 */
	Tenseur tnphi ; 

	/// Norm of \c u_euler 
	Tenseur uuu ;		
	
	/// Metric potential \f$\nu = \ln N\f$ = \c logn_auto 
	Tenseur& logn ;	

	/** Part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the matter terms
	 */
	Tenseur nuf ;	

	/** Part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the quadratic terms
	 */
	Tenseur nuq ;	

	/// Metric potential \f$\zeta = \ln(AN)\f$ = \c beta_auto 
	Tenseur& dzeta ;	

	/// Metric potential \f$\tilde G = (NB-1) r\sin\theta\f$
	Tenseur tggg ; 

	/** Vector \f$W^i\f$ used in the decomposition of \c shift ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	Tenseur w_shift ; 
	
	/** Scalar \f$\chi\f$ used in the decomposition of \c shift ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 */
	Tenseur khi_shift ; 

	/** Tensor \f${\tilde K_{ij}}\f$ related to the extrinsic curvature
	 *  tensor by \f${\tilde K_{ij}} = B^{-2} K_{ij}\f$.
	 *  \c tkij  contains the Cartesian components of
	 *  \f${\tilde K_{ij}}\f$. 
	 */
	Tenseur_sym tkij ; 

	/** Scalar \f$A^2 K_{ij} K^{ij}\f$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of \f$N^\varphi\f$ by
	 * \f[
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \f]
	 * In particular it is related to the quantities \f$k_1\f$ and \f$k_2\f$
	 * introduced by Eqs.~(3.7) and (3.8) of 
	 * Bonazzola et al. \a Astron. \a Astrophys. \b 278 , 421 (1993)
	 * by 
	 * \f[
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \f]
	 */
	 Tenseur ak_car ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c nuf  by means of
	 *  \c Map_et::poisson .
	 */
	Cmp ssjm1_nuf ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c nuq  by means of
	 *  \c Map_et::poisson .
	 */
	Cmp ssjm1_nuq ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c dzeta .
	 */
	Cmp ssjm1_dzeta ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for \c tggg .
	 */
	Cmp ssjm1_tggg ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for the scalar \f$\chi\f$ by means of
	 *  \c Map_et::poisson . 
	 *  \f$\chi\f$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector \f$N^i\f$
	 */
	 Cmp ssjm1_khi ; 
	 
	/** Effective source at the previous step for the resolution of 
	 *  the vector Poisson equation for \f$W^i\f$.
	 *  \f$W^i\f$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector \f$N^i\f$
	 *  (Components with respect to the Cartesian triad associated with 
	 *   the mapping \c mp )
	 */
	 Tenseur ssjm1_wshift ; 
	 
    // Derived data : 
    // ------------
    protected:
	
	mutable double* p_angu_mom ;	///< Angular momentum 
	mutable double* p_tsw ;		///< Ratio T/W
	mutable double* p_grv2 ;	///< Error on the virial identity GRV2
	mutable double* p_grv3 ;	///< Error on the virial identity GRV3
	mutable double* p_r_circ ;	///< Circumferential radius
	mutable double* p_aplat ;	///< Flatening r_pole/r_eq
	mutable double* p_z_eqf ;	///< Forward redshift factor at equator
	mutable double* p_z_eqb ;	///< Backward redshift factor at equator
	mutable double* p_z_pole ;	///< Redshift factor at North pole
	mutable double* p_mom_quad ;	///< Quadrupole moment
	mutable double* p_r_isco ;	///< Circumferential radius of the ISCO
	mutable double* p_f_isco ;	///< Orbital frequency of the ISCO
	/// Specific energy of a particle on the ISCO 
	mutable double* p_espec_isco ;	
	/// Specific angular momentum of a particle on the ISCO
	mutable double* p_lspec_isco ;	
        mutable double* p_f_eq ;        ///< Orbital frequency at the equator
	
	 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be \c true  for a relativistic
	 *			star,  \c false  for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 */
	Etoile_rot(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i) ;			
	
	
	Etoile_rot(const Etoile_rot& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	\c sauve )
	 */
	Etoile_rot(Map& mp_i, const Eos& eos_i, FILE* fich) ;    		

	virtual ~Etoile_rot() ;			///< Destructor


    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0  all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to \c ETATNONDEF  (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Etoile_rot 
	void operator=(const Etoile_rot& ) ;	
	
    // Accessors
    // ---------
    public:
	/** Returns the central value of the rotation angular velocity 
	 *  (\c [f_unit] )
	 */ 
	virtual double get_omega_c() const ;	    

	/// Returns the metric factor \e B 
	const Tenseur& get_bbb() const {return bbb;} ; 

	/// Returns the square of the metric factor \e B 
	const Tenseur& get_b_car() const {return b_car;} ; 

	/// Returns the metric coefficient \f$N^\varphi\f$
	const Tenseur& get_nphi() const {return nphi;} ; 

	/** Returns the component \f$\tilde N^\varphi = N^\varphi r\sin\theta\f$ 
	 *  of the shift vector
	 */
	const Tenseur& get_tnphi() const {return tnphi;} ; 

	/// Returns the norm of \c u_euler 
	const Tenseur& get_uuu() const {return uuu;} ;		
	
	/// Returns the metric potential \f$\nu = \ln N\f$ = \c logn_auto 
	const Tenseur& get_logn() const {return logn;} ;	

	/** Returns the part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the matter terms
	 */
	const Tenseur& get_nuf() const {return nuf;} ;	

	/** Returns the Part of the Metric potential \f$\nu = \ln N\f$ = \c logn 
	 *  generated by the quadratic terms
	 */
	const Tenseur& get_nuq() const {return nuq;} ;	

	/// Returns the Metric potential \f$\zeta = \ln(AN)\f$ = \c beta_auto 
	const Tenseur& get_dzeta() const {return dzeta;} ;	

	/// Returns the Metric potential \f$\tilde G = (NB-1) r\sin\theta\f$
	const Tenseur& get_tggg() const {return tggg;} ; 

	/** Returns the vector \f$W^i\f$ used in the decomposition of 
	 *  \c shift ,
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	const Tenseur& get_w_shift() const {return w_shift;} ; 
	
	/** Returns the scalar \f$\chi\f$ used in the decomposition of 
	 *  \c shift  
	 *  following Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 * NB: \c w_shift  contains the components of \f$W^i\f$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping \c mp . 
	 */
	const Tenseur& get_khi_shift() const {return khi_shift;} ; 

	/** Returns the tensor \f${\tilde K_{ij}}\f$ related to the extrinsic 
	 *  curvature tensor by \f${\tilde K_{ij}} = B^{-2} K_{ij}\f$.
	 *  \c tkij  contains the Cartesian components of
	 *  \f${\tilde K_{ij}}\f$. 
	 */
	const Tenseur_sym& get_tkij() const {return tkij;} ; 

	/** Returns the scalar \f$A^2 K_{ij} K^{ij}\f$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of \f$N^\varphi\f$ by
	 * \f[
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \f]
	 * In particular it is related to the quantities \f$k_1\f$ and \f$k_2\f$
	 * introduced by Eqs.~(3.7) and (3.8) of 
	 * Bonazzola et al. \a Astron. \a Astrophys. \b 278 , 421 (1993)
	 * by 
	 * \f[
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \f]
	 */
	 const Tenseur& get_ak_car() const {return ak_car;} ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    ///< Save in a file
    
	/// Display in polytropic units
	virtual void display_poly(ostream& ) const ; 

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

	/// Printing of some informations, excluding all global quantities
	virtual void partial_display(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	
	/** Description of the stellar surface: returns a 2-D \c Itbl 
	 *	containing the 
	 *	values of the domain index \e l  on the surface at the 
	 *	collocation points in \f$(\theta', \phi')\f$.
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member \c ent ) vanishes.
	 */
	virtual const Itbl& l_surf() const ; 
	
	virtual double mass_b() const ;	    ///< Baryon mass
	virtual double mass_g() const ;	    ///< Gravitational mass
	virtual double angu_mom() const ;	///< Angular momentum 
	virtual double tsw() const ;		///< Ratio T/W

	/** Error on the virial identity GRV2.
	 *  This indicator is only valid for relativistic computations.
	 */
	virtual double grv2() const ;	

	/** Error on the virial identity GRV3.
	 *  The error is computed as the integral defined
	 *  by Eq. (43) of [Gourgoulhon and Bonazzola, 
	 *  \a Class. \a Quantum \a Grav. \b 11, 443 (1994)] divided by
	 *  the integral of the matter terms.
	 * 
	 *  @param ost output stream to give details of the computation;
	 *		if set to 0x0 [default value], no details will be
	 *		given.
	 *   
	 */
	virtual double grv3(ostream* ost = 0x0) const ;	

	virtual double r_circ() const ;	///< Circumferential radius
	virtual double aplat() const ;	///< Flatening r_pole/r_eq
	virtual double z_eqf() const ;	///< Forward redshift factor at equator
	virtual double z_eqb() const ;	///< Backward redshift factor at equator
	virtual double z_pole() const ;	///< Redshift factor at North pole
    
	/** Quadrupole moment.
	 *  The quadrupole moment \e Q is defined according to Eq. (7) of
	 *  [Salgado, Bonazzola, Gourgoulhon and Haensel, \a Astron. \a Astrophys.
	 *   \b 291 , 155 (1994)]. At the Newtonian limit it is related to
	 *  the component \f${\bar I}_{zz}\f$ of the MTW (1973) reduced quadrupole 
	 *  moment \f${\bar I}_{ij}\f$ by: \f$Q = -3/2 {\bar I}_{zz}\f$. 
	 *  Note that \e Q is the negative of the quadrupole moment defined 
	 *  by Laarakkers and Poisson, \a Astrophys. \a J. \b 512 , 282 (1999).
	 */
	virtual double mom_quad() const ;	

	/** Circumferential radius of the innermost stable circular orbit (ISCO).	
	 *
	 *  @param ost output stream to give details of the computation;
	 *		if set to 0x0 [default value], no details will be
	 *		given.
	 */
	virtual double r_isco(ostream* ost = 0x0) const ;	
 	
 	/// Orbital frequency at the innermost stable circular orbit (ISCO).	
 	virtual double f_isco() const ;	

	/// Energy of a particle on the ISCO 
 	virtual double espec_isco() const ;	
	
	/// Angular momentum of a particle on the ISCO
 	virtual double lspec_isco() const ;	


	/** Computation of frequency of eccentric orbits.
	 * 
	 *  @param ecc eccentricity of the orbit
	 *  @param periasrt periastron of the orbit
	 *  @param ost output stream to give details of the computation;
	 *		if set to 0x0 [default value], no details will be
	 *		given.
	 * 
	 *  @return orbital frequency
	 */
	virtual double f_eccentric(double ecc, double periast, 
				   ostream* ost = 0x0) const ; 

        /// Orbital frequency at the equator.
	virtual double f_eq() const ;
	

    // Computational routines
    // ----------------------
    public: 
	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame.
	 *
	 *  The calculation is performed starting from the quantities
	 *  \c ent , \c ener , \c press , and \c a_car ,  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  \c gam_euler , \c u_euler , \c ener_euler , \c s_euler . 
	 * 
	 */
	virtual void hydro_euler() ; 
	
	/** Computes metric coefficients from known potentials. 
	 * 
	 *  The calculation is performed starting from the quantities
	 *  \c logn ,  \c dzeta , \c tggg  and \c shift , 
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  \c nnn , \c a_car ,  \c bbb  and \c b_car . 
	 * 
	 */
	void update_metric() ; 
		
	/** Computes \c shift  from \c w_shift  and \c khi_shift 
	 *  according to Shibata's prescription 
	 *  [\a Prog. \a Theor. \a Phys. \b 101 , 1199 (1999)] :
	 * \f[
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \f]
	 */
	void fait_shift() ; 
	
	/** Computes \c tnphi  and \c nphi  from the Cartesian 
	 *   components of the shift, stored in \c shift .
	 */
	void fait_nphi() ; 
		
	/** Computes \c tkij  and \c ak_car  from 
	 *  \c shift , \c nnn  and \c b_car .
	 */
	void extrinsic_curvature() ;
	
	/** Computes the coefficient \f$\lambda\f$ which ensures that the
	 *	GRV2 virial identity is satisfied.
	 *  \f$\lambda\f$ is the coefficient by which one must multiply
	 *  the quadratic source term \f$\sigma_q\f$ of the 2-D Poisson equation
	 *	\f[
	 *		\Delta_2 u = \sigma_m + \sigma_q
	 *	\f]
	 *  in order that the total source does not contain any monopolar term,
	 *  i.e. in order that
	 *  \f[
	 *		\int_0^{2\pi} \int_0^{+\infty} \sigma(r, \theta)
	 *				\, r \, dr \, d\theta = 0	    \ ,
	 *  \f]
	 *  where \f$\sigma = \sigma_m + \sigma_q\f$.
	 *	\f$\lambda\f$ is computed according to the formula
	 *  \f[
	 *		\lambda = - { \int_0^{2\pi} \int_0^{+\infty} \sigma_m(r, \theta)
	 *				\, r \, dr \, d\theta	    \over
	 * 			\int_0^{2\pi} \int_0^{+\infty} \sigma_q(r, \theta)
	 *				\, r \, dr \, d\theta } \ .
	 *  \f]
	 *  Then, by construction, the new source
	 *	\f$\sigma' = \sigma_m + \lambda \sigma_q\f$ has a vanishing monopolar
	 *  term.
	 *
	 *	@param sou_m [input] matter source term \f$\sigma_m\f$
	 *	@param sou_q [input] quadratic source term \f$\sigma_q\f$
	 *  @return	value of \f$\lambda\f$
	 */
	static double lambda_grv2(const Cmp& sou_m, const Cmp& sou_q) ;
		
	/** Computes an equilibrium configuration.
	 *  
	 *  @param ent_c  [input] Central enthalpy 
	 *  @param omega0  [input] Requested angular velocity 
	 *			     (if \c fact_omega=1. )
	 *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
	 *			      1. = otherwise.
	 *  @param nzadapt  [input] Number of (inner) domains where the mapping 
	 *			    adaptation to an iso-enthalpy surface
	 *			    should be performed
	 *  @param ent_limit [input] 1-D \c Tbl  of dimension \c nzet  which
	 *				defines the enthalpy at the outer boundary
	 *				of each domain
	 *  @param icontrol [input] Set of integer parameters (stored as a
	 *			    1-D \c Itbl  of size 8) to control the 
	 *			    iteration: 
	 *	\li \c icontrol(0) = mer_max  : maximum number of steps 
	 *	\li \c icontrol(1) = mer_rot  : step at which the rotation is 
	 *				      switched on 
	 *	\li \c icontrol(2) = mer_change_omega  : step at which the rotation
	 *			  velocity is changed to reach the final one  
	 *	\li \c icontrol(3) = mer_fix_omega  :  step at which the final
	 *			    rotation velocity must have been reached  
	 *	\li \c icontrol(4) = mer_mass  : the absolute value of 
	 *			    \c mer_mass  is the step from which the 
	 *			    baryon mass is forced to converge, 
	 *			    by varying the central enthalpy 
	 *			    (\c mer_mass>0 ) or the angular 
	 *			    velocity (\c mer_mass<0 ) 
	 *	\li \c icontrol(5) = mermax_poisson  : maximum number of steps in 
	 *				\c Map_et::poisson  
	 *	\li \c icontrol(6) = mer_triax  : step at which the 3-D 
	 *				perturbation is switched on 
	 *	\li \c icontrol(7) = delta_mer_kep  : number of steps
	 *			    after \c mer_fix_omega  when \c omega 
	 *			    starts to be increased by \c fact_omega 
	 *			    to search for the Keplerian velocity
	 * 	 
	 *  @param control [input] Set of parameters (stored as a 
	 *			    1-D \c Tbl  of size 7) to control the 
	 *			    iteration: 
	 *	\li \c control(0) = precis  : threshold on the enthalpy relative 
	 *				change for ending the computation 
	 *	\li \c control(1) = omega_ini  : initial angular velocity, 
	 *			    switched on only if \c mer_rot<0 , 
	 *			    otherwise 0 is used  
	 *	\li \c control(2) = relax  : relaxation factor in the main 
	 *				   iteration  
	 *	\li \c control(3) = relax_poisson  : relaxation factor in 
	 *				   \c Map_et::poisson  
	 *	\li \c control(4) = thres_adapt  :  threshold on dH/dr for 
	 *			    freezing the adaptation of the mapping 
	 *	\li \c control(5) = ampli_triax  :  relative amplitude of 
	 *			    the 3-D perturbation 
	 *	\li \c control(6) = precis_adapt  : precision for 
	 *			    \c Map_et::adapt 
	 *
	 *  @param mbar_wanted [input] Requested baryon mass (effective only 
	 *				if \c mer_mass > \c mer_max )
	 *  @param aexp_mass [input] Exponent for the increase factor of the 
	 *			      central enthalpy to converge to the 
	 *			      requested baryon mass
	 *  @param diff [output]   1-D \c Tbl  of size 7 for the storage of 
	 *			    some error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuf    
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuq    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c dzeta    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    Poisson equation for \c tggg    
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift  (x comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift  (y comp.)   
	 */
	virtual void equilibrium(double ent_c, double omega0, double fact_omega, 
			 int nzadapt, const Tbl& ent_limit,
			 const Itbl& icontrol, const Tbl& control,
			 double mbar_wanted, double aexp_mass, 
			 Tbl& diff, Param* = 0x0) ;
	

};




#endif
