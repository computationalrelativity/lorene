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
 * Base class for stars.
 * 
 * An {\tt Etoile} is constructed upon (i) a mapping 
 * (derived class of {\tt Map}), the center of which defines the center of the 
 * star, and (ii) an equation of state (derived class of {\tt Eos}).  
 * It contains tensor fields (class {\tt Tenseur}) which describle the
 * hydrodynamical quantities as well as the gravitational field (spacetime
 * metric). 
 * 
 * According to the 3+1 formalism, the spacetime metric is written
 * \begin{equation} \label{eetoilemetrique}
 *   ds^2 = - (N^2 - N_i N^i) dt^2 - 2 N_i \,  dt\,  dx^i
 *	    + A^2  \,  {\tilde h}_{ij} \,  dx^i dx^j
 * \end{equation}
 * where ${\tilde h}_{ij}$ is a 3-metric, the exact form of which is specified
 * in the derived classes of {\tt Etoile}. The base class {\tt Etoile} by
 * itself provides only storage for the lapse function {\it N} (member {\tt nnn}), 
 * the shift vector $N^i$ (member {\tt shift}) and the conformal factor
 * $A^2$ (member {\tt a\_car}). 
 * 
 * The 3+1 formalism introduces two kinds of priviledged observers: the
 * fluid comoving observer and the Eulerian observer, whose 4-velocity
 * is the unit future directed normal to the {\it t} = const hypersurfaces. 
 * The hydrodynamical quantities measured by the fluid observer correspond
 * to the members {\tt ent}, {\tt nbar}, {\tt ener}, and {\tt press}. 
 * The hydrodynamical quantities measured by the Eulerian observer correspond
 * to the members {\tt ener\_euler}, {\tt s\_euler}, {\tt gam\_euler}, and 
 * {\tt u\_euler}.  
 * 
 * A star of class {\tt Etoile} can be either relativistic or Newtonian, 
 * depending on the boolean indicator {\tt relativistic}. For a Newtonian
 * star, the metric coefficients {\it N} and {\it A} are set to 1,  and $N^i$ is
 * set to zero; the only relevant gravitational quantity in this case is
 * {\tt logn\_auto} which represents the (regular part of) the
 * Newtonian gravitational potential
 * generated by the star. 
 *
 * @version #$Id$#
 */
class Etoile {

    // Data : 
    // -----
    protected:
	Map& mp ;	    /// Mapping associated with the star

	/** Number of domains of {\tt *mp} occupied by the star
	 * 
	 */
	int nzet ;
	
	/** Indicator of relativity: {\tt true} for a relativistic star,
	 *	{\tt false} for a Newtonian one. 
	 */
	bool relativistic ;
	
	/** $1/c^2$ : {\tt unsurc2 = 1} for a relativistic star,
	 *  0 for a Newtonian one. 
	 */
	double unsurc2 ; 	     
    
	/** Index of regularity of the gravitational potential {\tt logn\_auto}.
	 *  If {\tt k\_div=0}, {\tt logn\_auto} contains the total potential
	 *  generated principaly by the star, otherwise it should be
	 *  supplemented by {\tt logn\_auto\_div}. 
	 */
	int k_div ; 
    
	const Eos& eos ;   /// Equation of state of the stellar matter
    
	// Fluid quantities with respect to the fluid frame
	// ------------------------------------------------
	/// Log-enthalpy (relativistic case) or specific enthalpy (Newtonian case)
	Tenseur ent ;	  
	
	Tenseur nbar ; 	   /// Baryon density in the fluid frame
	Tenseur ener ;	   /// Total energy density in the fluid frame
	Tenseur press ;	   /// Fluid pressure
    
	// Fluid quantities with respect to the Eulerian frame
	// ---------------------------------------------------
	Tenseur ener_euler ; /// Total energy density in the Eulerian frame

	/// Trace of the stress tensor in the Eulerian frame
	Tenseur s_euler ;   
	
	/// Lorentz factor between the fluid and Eulerian observers 
	Tenseur gam_euler ; 
	
	/// Fluid 3-velocity with respect to the Eulerian observer
	Tenseur u_euler ; 
	
	// Metric potentials
	// -----------------
	
	/** Total of the logarithm of the part of the lapse {\it N} 
	 *   generated principaly by the star. In the Newtonian case, 
	 *   this is the Newtonian gravitational potential
	 *   (in units of $c^2$). 
	 */
	Tenseur logn_auto ;

	/** Regular part of the logarithm of the part of the lapse {\it N} 
	 *   generated principaly by the star. In the Newtonian case, 
	 *   this is the Newtonian gravitational potential
	 *   (in units of $c^2$). 
	 */
	Tenseur logn_auto_regu ;
	
	/** Divergent part (if {\tt k\_div != 0} ) 
	 *  of the logarithm of the part of the lapse {\it N} 
	 *   generated principaly by the star. 
	 */
	Tenseur logn_auto_div ; 
	
	/** Gradient of {\tt logn\_auto\_div} (if {\tt k\_div != 0} )
	 */
	Tenseur d_logn_auto_div ; 

	/** Logarithm of the part of the product {\it AN} generated principaly by
	 *   by the star
	 */
	Tenseur beta_auto ; 
	
	/// Total lapse function 
	Tenseur nnn ; 
	
	/// Total shift vector
	Tenseur shift ;
	
	/// Total conformal factor $A^2$
	Tenseur a_car ; 
	
    // Derived data : 
    // ------------
    protected:
	/// Coordinate radius at $\phi=0$, $\theta=\pi/2$. 
	mutable double* p_ray_eq ; 

	/// Coordinate radius at $\phi=\pi/2$, $\theta=\pi/2$. 
	mutable double* p_ray_eq_pis2 ;

	/// Coordinate radius at $\phi=\pi$, $\theta=\pi/2$. 
	mutable double* p_ray_eq_pi ;
	
	/// Coordinate radius at $\theta=0$. 
	mutable double* p_ray_pole ;
	
	/** Description of the stellar surface: 2-D {\tt Itbl} containing the 
	 *	values of the domain index {\it l} on the surface at the 
	 *	collocation points in $(\theta', \phi')$
	 */
	mutable Itbl* p_l_surf ; 
	
	/** Description of the stellar surface: 2-D {\tt Tbl} containing the 
	 *	values of the radial coordinate $\xi$ on the surface at the 
	 *	collocation points in $(\theta', \phi')$
	 */
	mutable Tbl* p_xi_surf ; 
	
	mutable double* p_mass_b ;	/// Baryon mass
	mutable double* p_mass_g ;	/// Gravitational mass


    // Constructors - Destructor
    // -------------------------
    public:
    
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be {\tt true} for a relativistic
	 *			star,  {\tt false} for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 * 
	 */
	Etoile(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i) ;			
	
	
	Etoile(const Etoile& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Etoile(Map& mp_i, const Eos& eos_i, FILE* fich) ;    		

	virtual ~Etoile() ;			/// Destructor

	
    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to {\tt ETATNONDEF} (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another {\tt Etoile}
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
	 *  observer from those in the fluid frame ({\tt nbar}, {\tt ener}
	 *  and {\tt press}).
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

    // Accessors
    // ---------
    public:
	/// Returns the mapping
	const Map& get_mp() const {return mp; } ; 

	/// Returns the number of domains occupied by the star
	int get_nzet() const {return nzet; } ; 

	/** Returns {\tt true} for a relativistic star, {\tt false} for 
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

	/** Returns the logarithm of the part of the lapse {\it N} generated 
	 *   principaly by the star.
	 *   In the Newtonian case, this is the Newtonian
	 *   gravitational potential (in units of $c^2$). 
	 */
	const Tenseur& get_logn_auto() const {return logn_auto;} ;

	/** Returns the regular part of the logarithm of the part of
	 *   the lapse {\it N} generated principaly by the star.
	 *   In the Newtonian case, this is the Newtonian
	 *   gravitational potential (in units of $c^2$). 
	 */
	const Tenseur& get_logn_auto_regu() const {return logn_auto_regu;} ;

	/** Returns the divergent part of the logarithm of the part of
	 *   the lapse {\it N} generated principaly by the star.
	 *   In the Newtonian case, this is the diverging part of
	 *   the Newtonian gravitational potential (in units of $c^2$). 
	 */
	const Tenseur& get_logn_auto_div() const {return logn_auto_div;} ;

	/** Returns the gradient of {\tt logn\_auto\_div}
	 */
	const Tenseur& get_d_logn_auto_div() const {return d_logn_auto_div;} ;

	/** Returns the logarithm of the part of the product {\it AN} generated 
	 *  principaly by the star.
	 */
	const Tenseur& get_beta_auto() const {return beta_auto;} ;

	/// Returns the total lapse function {\it N}
	const Tenseur& get_nnn() const {return nnn;} ;

	/// Returns the total shift vector $N^i$
	const Tenseur& get_shift() const {return shift;} ;

	/// Returns the total conformal factor $A^2$
	const Tenseur& get_a_car() const {return a_car;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    /// Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Etoile& ) ;	

    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
	/// Coordinate radius at $\phi=0$, $\theta=\pi/2$ [r\_unit].
	double ray_eq() const ; 
	
	/// Coordinate radius at $\phi=\pi/2$, $\theta=\pi/2$ [r\_unit].
	double ray_eq_pis2() const ; 
	
	/// Coordinate radius at $\phi=\pi$, $\theta=\pi/2$ [r\_unit].
	double ray_eq_pi() const ; 
	
	/// Coordinate radius at $\theta=0$ [r\_unit]. 
	double ray_pole() const ; 
    
	/** Description of the stellar surface: returns a 2-D {\tt Itbl} 
	 *	containing the 
	 *	values of the domain index {\it l} on the surface at the 
	 *	collocation points in $(\theta', \phi')$.
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member {\tt ent}) vanishes.
	 */
	virtual const Itbl& l_surf() const ; 
	
	/** Description of the stellar surface: returns a 2-D {\tt Tbl} 
	 *	containing the values of the radial coordinate $\xi$ 
	 *	on the surface at the 
	 *	collocation points in $(\theta', \phi')$. 
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member {\tt ent}) vanishes.
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
 * Class for stars in binary system.
 *
 * This class implements the formalism for corotating or irrotational
 * systems presented in Bonazzola, Gourgoulhon \& Marck {\sl Phys. Rev. Lett.}
 * {\bf 82}, 892 (1999). In particular, the conformal 3-metric 
 * $\tilde h_{ij}$ introduced in Eq.~(\ref{eetoilemetrique}) is flat. 
 *
 * An {\tt Etoile\_bin} can be construted in two states, represented by
 * the {\tt bool} member {\tt irrotational}: (i) irrotational
 * (i.e. the fluid motion is irrotational) or (ii) rigidly corotating 
 * with respect to the orbital motion (synchronized binary). 
 *
 * @version #$Id$#
 */
class Etoile_bin : public Etoile {

    // Data : 
    // -----
    protected:
	/** {\tt true} for an irrotational star, {\tt false} for a
	 *  corotating one
	 */
	bool irrotational ; 
	
	/** Reference triad ("absolute frame"), with
	 *  respect to which the components of all the member {\tt Tenseur}'s
	 *  are defined, except for {\tt w\_shift} and {\tt ssjm1\_wshift}.
	 */
	const Base_vect& ref_triad ; 
	
	/** Scalar potential $\Psi_0$ of the non-translational part of
	 *  fluid 4-velocity (in the irrotational case)
	 */
	Tenseur psi0 ; 

	/** Gradient of $\Psi$ (in the irrotational case)
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur d_psi ; 
	
	/** Spatial projection of the fluid 3-velocity with respect to  
	 *  the co-orbiting observer. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
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
	
	/** Gradient of {\tt logn\_auto}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur d_logn_auto ; 

	/** Gradient of {\tt logn\_auto\_regu}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur d_logn_auto_regu ; 

	/** Gradient of {\tt logn\_comp}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur d_logn_comp ; 

	/** Part of the logarithm of {\it AN} generated principaly by the 
	 *  companion star. 
	 */
	Tenseur beta_comp ; 
	
	/** Gradient of {\tt beta\_auto}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur d_beta_auto ; 

	/** Gradient of {\tt beta\_comp}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur d_beta_comp ; 

	/** Part of the shift vector $N^i$ generated principaly by the star.
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur shift_auto ; 
	
	/** Part of the shift vector $N^i$ generated principaly by the 
	 *  companion star. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur shift_comp ; 
	
	/** Vector $W^i$ used in the decomposition of {\tt shift\_auto},
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 * NB: {\tt w\_shift} contains the components of $W^i$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping {\tt mp}. 
	 */
	Tenseur w_shift ; 
	
	/** Scalar $\chi$ used in the decomposition of {\tt shift\_auto},
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 */
	Tenseur khi_shift ; 

	/** Part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$
	 *  generated by {\tt shift\_auto}. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur_sym tkij_auto ;
	
	/** Part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$
	 *  generated by {\tt shift\_comp}. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur_sym tkij_comp ;
	
	/** Part of the scalar $A^2 K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto}, i.e. 
	 *    $A^2 K_{ij}^{\rm auto} K^{ij}_{\rm auto}$ 
	 */
	Tenseur akcar_auto ;
	
	/** Part of the scalar $A^2 K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto} and {\tt shift\_comp}, i.e. 
	 *    $A^2 K_{ij}^{\rm auto} K^{ij}_{\rm comp}$ 
	 */
	Tenseur akcar_comp ;
	
	/** 3-vector shift, divided by {\it N}, of the rotating coordinates,
	 *  $B^i/N$. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	Tenseur bsn ; 
	
	/// Centrifugal potential
	Tenseur pot_centri ; 	
	
	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt logn\_auto} by means of
	 *  {\tt Map\_et::poisson}.
	 */
	Cmp ssjm1_logn ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt beta\_auto} by means of
	 *  {\tt Map\_et::poisson}.
	 */
	Cmp ssjm1_beta ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for the scalar $\chi$ by means of
	 *  {\tt Map\_et::poisson}. 
	 *  $\chi$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector $N^i$
	 */
	 Cmp ssjm1_khi ; 
	 
	/** Effective source at the previous step for the resolution of 
	 *  the vector Poisson equation for $W^i$ by means of
	 *  {\tt Map\_et::poisson}. 
	 *  $W^i$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector $N^i$
	 *  (Components with respect to the Cartesian triad associated with 
	 *   the mapping {\tt mp})
	 */
	 Tenseur ssjm1_wshift ; 
	 
	 /**
	  * Function used to construct the part of $K^{ij}$ generated by 
	  * the star  from the total $K^{ij}$. Only used for a binary system
	  * where the other member is a black hole.
	  * 
	  * Mainly this {\tt Cmp} is 1 around the hole and 0 around the companion
	  * and the sum of {tt decouple} for the hole and his companion is 1 
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
	 * @param relat should be {\tt true} for a relativistic
	 *			star,  {\tt false} for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 * @param irrot should be {\tt true} for an irrotational star, 
	 *		    {\tt false} for a corotating one
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    {\tt Tenseur}'s are defined, except for {\tt w\_shift}
	 *	    and {\tt ssjm1\_wshift} whose components are defined
	 *	    with respect to the mapping {\tt mp} Cartesian triad. 
	 */
	Etoile_bin(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		   bool irrot, const Base_vect& ref_triad_i) ;			
	
	
	Etoile_bin(const Etoile_bin& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    {\tt Tenseur}'s are defined, except for {\tt w\_shift}
	 *	    and {\tt ssjm1\_wshift} whose components are defined
	 *	    with respect to the mapping {\tt mp} Cartesian triad. 
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Etoile_bin(Map& mp_i, const Eos& eos_i, const Base_vect& ref_triad_i, 
		   FILE* fich) ;    		

	virtual ~Etoile_bin() ;			/// Destructor


    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to {\tt ETATNONDEF} (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another {\tt Etoile\_bin}
	void operator=(const Etoile_bin& ) ;	
	
	/** Read/write the part of the lapse logarithm (gravitational potential 
	 *  at the Newtonian limit) generated principaly by the companion star. 
	 */
	Tenseur& set_logn_comp() ;

	/// Read/write the centrifugal potential
	Tenseur& set_pot_centri() ;
	
	/// Read/write of {\tt w\_shift}
	Tenseur& set_w_shift() ;
	
	/// Read/write of {\tt khi\_shift}
	Tenseur& set_khi_shift() ;
		
    // Accessors
    // ---------
    public:
	/** Returns {\tt true} for an irrotational motion, {\tt false} for 
	 *  a corotating one. 
	 */
	bool is_irrotational() const {return irrotational; } ; 

	/// Returns the non-translational part of the velocity potential
	const Tenseur& get_psi0() const {return psi0;} ;

	/** Returns the gradient of the velocity potential 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_d_psi() const {return d_psi;} ;

	/** Returns the spatial projection of the fluid 3-velocity with 
	 *  respect to the co-orbiting observer. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
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

	/** Returns the gradient of {\tt logn\_auto}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_d_logn_auto() const {return d_logn_auto;} ;

	/** Returns the gradient of {\tt logn\_auto\_regu}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_d_logn_auto_regu() const {return d_logn_auto_regu;} ;

	/** Returns the gradient of {\tt logn\_comp}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_d_logn_comp() const {return d_logn_comp;} ;

	/** Returns the part of the logarithm of {\it AN} generated principaly 
	 *  by the companion star. 
	 */
	const Tenseur& get_beta_comp() const {return beta_comp;} ;

	/** Returns the gradient of {\tt beta\_auto}
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_d_beta_auto() const {return d_beta_auto;} ;

	/** Returns the gradient of {\tt beta\_comp} 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_d_beta_comp() const {return d_beta_comp;} ;

	/** Returns the part of the shift vector $N^i$ generated principaly 
	 * by the star.
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_shift_auto() const {return shift_auto;} ;

	/** Returns the part of the shift vector $N^i$ generated principaly 
	 *   by the companion star. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_shift_comp() const {return shift_comp;} ;

	/** Returns the vector $W^i$ used in the decomposition of 
	 *  {\tt shift\_auto},
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 * NB: {\tt w\_shift} contains the components of $W^i$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping {\tt mp}. 
	 */
	const Tenseur& get_w_shift() const {return w_shift;} ; 
	
	/** Returns the scalar $\chi$ used in the decomposition of 
	 *  {\tt shift\_auto} 
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 * NB: {\tt w\_shift} contains the components of $W^i$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping {\tt mp}. 
	 */
	const Tenseur& get_khi_shift() const {return khi_shift;} ; 

	/** Returns the part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$ generated by {\tt shift\_auto}. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur_sym& get_tkij_auto() const {return tkij_auto;} ;

	/** Returns the part of the extrinsic curvature tensor 
	 *  $\tilde K^{ij} = A^2 K^{ij}$ generated by {\tt shift\_comp}. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur_sym& get_tkij_comp() const {return tkij_comp;} ;

	/** Returns the part of the scalar $A^2 K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto}, i.e. 
	 *    $A^2 K_{ij}^{\rm auto} K^{ij}_{\rm auto}$ 
	 */
	const Tenseur& get_akcar_auto() const {return akcar_auto;} ;

	/** Returns the part of the scalar $A^2 K_{ij} K^{ij}$
	 *  generated by {\tt shift\_auto} and {\tt shift\_comp}, i.e. 
	 *    $A^2 K_{ij}^{\rm auto} K^{ij}_{\rm comp}$ 
	 */
	const Tenseur& get_akcar_comp() const {return akcar_comp;} ;

	/** Returns the shift vector, divided by {\it N}, of the rotating 
	 *   coordinates, $B^i/N$. 
	 *  (Cartesian components with respect to {\tt ref\_triad})
	 */
	const Tenseur& get_bsn() const {return bsn;} ;

	/// Returns the centrifugal potential
	const Tenseur& get_pot_centri() const {return pot_centri;} ;
	/**
	 * Returns the function used to construct {\tt tkij\_auto} 
	 from {\tt tkij\_tot}.
	*/
	const Cmp get_decouple() const {return decouple ;}

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    /// Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

	/// Computes the auto part of $\left(L\beta\right)^{ij} (not stored)
	Tenseur_sym fait_taij_auto() const ;

    // Global quantities
    // -----------------
    public:
	/// Baryon mass
    	virtual double mass_b() const ;
	
	/// Gravitational mass
    	virtual double mass_g() const ;
	
	/** Absolute coordinate X of the barycenter of the baryon density, 
	 *  defined according to the formula
	 *  \begin{equation}
	 *    X_G := \int A^3 \Gamma_{\rm n} \,  n \,  X \, d^3x \ ,  
	 *  \end{equation}
	 *  where $\Gamma_{\rm n}$ is the Lorentz factor between the fluid 
	 *  and Eulerian observers.
	 */
    	virtual double xa_barycenter() const ;
	

    // Computational routines
    // ----------------------
    public: 
	/** Performs the scalar product of two tensors by contracting
	 *  the last index of {\tt t1} with the first index of {\tt t2}.
	 *  Both indices are supposed to be contravariant, so that a 
	 *  multiplication by $A^2$ is performed to lower one index. 
	 *  For instance, for two vectors $V^i$ and $W^i$, this function
	 *  returns the scalar $h_{ij} V^i W^j = A^2 f_{ij} V^i W^j$.  
	 */
	virtual Tenseur sprod(const Tenseur& t1, const Tenseur& t2) const ; 

	/** Computes the hydrodynamical quantities relative to the Eulerian
	 *  observer from those in the fluid frame, as well as 
	 *  {\tt wit\_w} and {\tt loggam}.  
	 *
	 *  The calculation is performed starting from the quantities
	 *  {\tt ent}, {\tt ener}, {\tt press}, {\tt a\_car} and {\tt bsn},  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  {\tt gam\_euler}, {\tt u\_euler}, {\tt ener\_euler}, {\tt s\_euler}, 
	 *  {\tt wit\_w} and {\tt loggam}. 
	 * 
	 */
	virtual void hydro_euler() ; 
	
	/** Computes metric coefficients from known potentials,
	 * when the companion is another star.
	 *
	 *  The calculation is performed starting from the quantities
	 *  {\tt logn\_auto},  {\tt beta\_auto}, {\tt shift\_auto},
	 *  {\tt comp.logn\_auto},  {\tt comp.beta\_auto},
	 *  {\tt comp.shift\_auto}
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  {\tt logn\_comp}, {\tt beta\_comp}, {\tt shift\_comp},
	 *  {\tt nnn},  {\tt a\_car},  {\tt shift},
	 *  {\tt d\_logn\_auto}, {\tt d\_beta\_auto}, {\tt tkij\_auto},
	 *  {\tt akcar\_auto}.
	 *
	 *  @param comp companion star.
	 *
	 */
	void update_metric(const Etoile_bin& comp) ;

	/** Computes metric coefficients from known potentials,
	 *  when the companion is a black hole.
	 *
	 *  @param comp companion black hole
	 *
	 */
	void update_metric(const Bhole& comp) ;

	/** Same as {\tt update\_metric(const Etoile\_bin\& )} but with
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
	 *  {\tt comp.d\_logn\_auto},  {\tt comp.d\_beta\_auto},
	 *  {\tt comp.tkij\_auto}
	 *  which are supposed to be up to date.
	 *  From these,  the following fields are updated:
	 *  {\tt d\_logn\_comp}, {\tt d\_beta\_comp}, {\tt tkij\_comp},
	 *  {\tt akcar\_comp}.
	 *
	 *  @param comp companion star.
	 *
	 */
	 void update_metric_der_comp(const Etoile_bin& comp) ;

	/** Computes the derivative of metric functions related to the
	 *  companion black hole.
	 *
	 *  @param comp companion BH.
	 *
	 */
	void update_metric_der_comp(const Bhole& comp) ;

	/** Computes the quantities {\tt bsn} and {\tt pot\_centri}.
	 * 
	 *  The calculation is performed starting from the quantities
	 *  {\tt nnn}, {\tt shift},  {\tt a\_car},  
	 *  which are supposed to be up to date.  
	 * 
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 *  @param x_axe  absolute X coordinate of the rotation axis
	 */
	virtual void kinematics(double omega, double x_axe) ; 
	
	/** Computes the gradient of the total velocity potential $\psi$. 
	 * 
	 */
	void fait_d_psi() ; 
	
	/** Computes {\tt shift\_auto} from {\tt w\_shift} and {\tt khi\_shift}
	 *  according to Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 */
	void fait_shift_auto() ; 
	
	/** Computes {\tt tkij\_auto} and {\tt akcar\_auto} from 
	 *  {\tt shift\_auto}, {\tt nnn} and {\tt a\_car}.
	 */
	virtual void extrinsic_curvature() ; 
		
	
	/** Computes an equilibrium configuration.
	 * 
	 *  The values of {\tt logn\_comp}, {\tt beta\_comp}, {\tt pot\_centri}
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map\_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map\_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map\_radial::poisson\_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map\_radial::poisson\_compact
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param fact [input]    1-D {\tt Tbl} for the input
	 *                          of some factors : \\
	 *          {\tt fact(0)} : A resizing factor for the first shell
	 *  @param diff [output]   1-D {\tt Tbl} for the storage of some
	 *			    error indicators : \\
	 *	    {\tt diff(0)} : Relative change in the enthalpy field
	 *			      between two successive steps \\
	 *	    {\tt diff(1)} : Relative error returned by the routine
	 *				{\tt Etoile\_bin::velocity\_potential}  
	 *	    {\tt diff(2)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt logn\_auto} \\  
	 *	    {\tt diff(3)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt beta\_auto} \\  
	 *	    {\tt diff(4)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (x comp.) \\  
	 *	    {\tt diff(5)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (y comp.) \\  
	 *	    {\tt diff(6)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (z comp.)   
	 */
	void equilibrium(double ent_c, int mermax, int mermax_poisson, 
			 double relax_poisson, int mermax_potvit, 
			 double relax_potvit, double thres_adapt, 
			 const Tbl& fact, Tbl& diff) ;


	/** Computes an equilibrium configuration by regularizing
	 *  the diverging source.
	 * 
	 *  The values of {\tt logn\_comp}, {\tt beta\_comp}, {\tt pot\_centri}
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map\_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map\_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map\_radial::poisson\_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map\_radial::poisson\_compact
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param fact [input]    1-D {\tt Tbl} for the input
	 *                          of some factors : \\
	 *          {\tt fact(0)} : A resizing factor for the first shell
	 *  @param diff [output]   1-D {\tt Tbl} for the storage of some
	 *			    error indicators : \\
	 *	    {\tt diff(0)} : Relative change in the enthalpy field
	 *			      between two successive steps \\
	 *	    {\tt diff(1)} : Relative error returned by the routine
	 *				{\tt Etoile\_bin::velocity\_potential}  
	 *	    {\tt diff(2)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt logn\_auto} \\  
	 *	    {\tt diff(3)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt beta\_auto} \\  
	 *	    {\tt diff(4)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (x comp.) \\  
	 *	    {\tt diff(5)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (y comp.) \\  
	 *	    {\tt diff(6)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (z comp.)   
	 */
	void equil_regular(double ent_c, int mermax, int mermax_poisson, 
			   double relax_poisson, int mermax_potvit, 
			   double relax_potvit, double thres_adapt, 
			   const Tbl& fact, Tbl& diff) ;


	/** Computes the non-translational part of the velocity scalar potential
	 *  $\psi0$ by solving the continuity equation.
	 *  
	 *  @param mermax  [input] Maximum number of steps in the iteration
	 *  @param precis  [input] Required precision: the iteration will
	 *			   be stopped when the relative difference
	 *			   on $\psi0$ between two successive steps
	 *			   is lower than {\tt precis}.
	 *  @param relax   [input] Relaxation factor.  
	 *
	 *  @return Relative error of the resolution obtained by comparing
	 *	    the operator acting on the solution with the source.
	 */
	double velocity_potential(int mermax, double precis, double relax) ;

	/** Performs a relaxation on {\tt ent}, {\tt logn\_auto},
	 *  {\tt beta\_auto} and {\tt shift\_auto}. 
	 * 
	 *  @param star_prev   [input] star at the previous step.
	 *  @param relax_ent   [input] Relaxation factor for {\tt ent} 
	 *  @param relax_met   [input] Relaxation factor for {\tt logn\_auto},
	 *			       {\tt beta\_auto}, {\tt shift\_auto}, 
	 *			       only if {\tt (mer \% fmer\_met == 0)}.
	 *  @param mer	       [input] Step number
	 *  @param fmer_met    [input] Step interval between metric updates
	 */
	void relaxation(const Etoile_bin& star_prev, double relax_ent, 
			double relax_met, int mer, int fmer_met) ;
	
	friend class Bin_ns_bh ; /// Friend class Bin_ns_bh
};



			//---------------------------//
			//    class Etoile_rot       //
			//---------------------------//

/**
 * Class for isolated rotating stars. 
 * 
 * The metric is
 * \begin{equation}
 *   ds^2 = - N^2 dt^2 + A^2 (dr^2 + r^2 d\theta^2)
 *		       + B^2 r^2 \sin^2\theta (d\varphi - N^\varphi dt)^2
 * \end{equation}
 *
 *
 * @version #$Id$#
 */
class Etoile_rot : public Etoile {

    // Data : 
    // -----
    protected:
	double omega ;	    /// Rotation angular velocity ({\tt [f\_unit]}) 

	/// Metric factor {\it B}
	Tenseur bbb ; 

	/// Square of the metric factor {\it B}
	Tenseur b_car ; 

	/// Metric coefficient $N^\varphi$
	Tenseur nphi ; 

	/** Component $\tilde N^\varphi = N^\varphi r\sin\theta$ of the
	 *  shift vector
	 */
	Tenseur tnphi ; 

	/// Norm of {\tt u\_euler}
	Tenseur uuu ;		
	
	/// Metric potential $\nu = \ln N$ = {\tt logn\_auto}
	Tenseur& logn ;	

	/** Part of the Metric potential $\nu = \ln N$ = {\tt logn}
	 *  generated by the matter terms
	 */
	Tenseur nuf ;	

	/** Part of the Metric potential $\nu = \ln N$ = {\tt logn}
	 *  generated by the quadratic terms
	 */
	Tenseur nuq ;	

	/// Metric potential $\zeta = \ln(AN)$ = {\tt beta\_auto}
	Tenseur& dzeta ;	

	/// Metric potential $\tilde G = (NB-1) r\sin\theta$
	Tenseur tggg ; 

	/** Vector $W^i$ used in the decomposition of {\tt shift},
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 * NB: {\tt w\_shift} contains the components of $W^i$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping {\tt mp}. 
	 */
	Tenseur w_shift ; 
	
	/** Scalar $\chi$ used in the decomposition of {\tt shift},
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 */
	Tenseur khi_shift ; 

	/** Tensor ${\tilde K_{ij}}$ related to the extrinsic curvature
	 *  tensor by ${\tilde K_{ij}} = B^{-2} K_{ij}$.
	 *  {\tt tkij} contains the Cartesian components of
	 *  ${\tilde K_{ij}}$. 
	 */
	Tenseur_sym tkij ; 

	/** Scalar $A^2 K_{ij} K^{ij}$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of $N^\varphi$ by
	 * \begin{equation}
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \end{equation}
	 * In particular it is related to the quantities $k_1$ and $k_2$
	 * introduced by Eqs.~(3.7) and (3.8) of 
	 * Bonazzola et al. {\sl Astron. Astrophys.} {\bf 278}, 421 (1993)
	 * by 
	 * \begin{equation}
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \end{equation}
	 */
	 Tenseur ak_car ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt nuf} by means of
	 *  {\tt Map\_et::poisson}.
	 */
	Cmp ssjm1_nuf ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt nuq} by means of
	 *  {\tt Map\_et::poisson}.
	 */
	Cmp ssjm1_nuq ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt dzeta}.
	 */
	Cmp ssjm1_dzeta ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for {\tt tggg}.
	 */
	Cmp ssjm1_tggg ; 

	/** Effective source at the previous step for the resolution of 
	 *  the Poisson equation for the scalar $\chi$ by means of
	 *  {\tt Map\_et::poisson}. 
	 *  $\chi$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector $N^i$
	 */
	 Cmp ssjm1_khi ; 
	 
	/** Effective source at the previous step for the resolution of 
	 *  the vector Poisson equation for $W^i$.
	 *  $W^i$ is an intermediate quantity for the resolution of the
	 *  elliptic equation for the shift vector $N^i$
	 *  (Components with respect to the Cartesian triad associated with 
	 *   the mapping {\tt mp})
	 */
	 Tenseur ssjm1_wshift ; 
	 
    // Derived data : 
    // ------------
    protected:
	
	mutable double* p_angu_mom ;	/// Angular momentum 
	mutable double* p_tsw ;		/// Ratio T/W
	mutable double* p_grv2 ;	/// Error on the virial identity GRV2
	mutable double* p_grv3 ;	/// Error on the virial identity GRV3
	mutable double* p_r_circ ;	/// Circumferential radius
	mutable double* p_aplat ;	/// Flatening r\_pole/r\_eq
	mutable double* p_z_eqf ;	/// Forward redshift factor at equator
	mutable double* p_z_eqb ;	/// Backward redshift factor at equator
	mutable double* p_z_pole ;	/// Redshift factor at North pole
	mutable double* p_mom_quad ;	/// Quadrupole moment
	mutable double* p_r_isco ;	/// Circumferential radius of the ISCO
	mutable double* p_f_isco ;	/// Orbital frequency of the ISCO
	/// Specific energy of a particle on the ISCO 
	mutable double* p_espec_isco ;	
	/// Specific angular momentum of a particle on the ISCO
	mutable double* p_lspec_isco ;	
        mutable double* p_f_eq ;        /// Orbital frequency at the equator
	
	 

    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be {\tt true} for a relativistic
	 *			star,  {\tt false} for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 */
	Etoile_rot(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i) ;			
	
	
	Etoile_rot(const Etoile_rot& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Etoile_rot(Map& mp_i, const Eos& eos_i, FILE* fich) ;    		

	virtual ~Etoile_rot() ;			/// Destructor


    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to {\tt 0x0} all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

	/** Sets to {\tt ETATNONDEF} (undefined state) the hydrodynamical 
	 *  quantities relative to the Eulerian observer.
	 */
	virtual void del_hydro_euler() ; 
	

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another {\tt Etoile\_rot}
	void operator=(const Etoile_rot& ) ;	
	
    // Accessors
    // ---------
    public:
	/** Returns the central value of the rotation angular velocity 
	 *  ({\tt [f\_unit]})
	 */ 
	virtual double get_omega_c() const ;	    

	/// Returns the metric factor {\it B}
	const Tenseur& get_bbb() const {return bbb;} ; 

	/// Returns the square of the metric factor {\it B}
	const Tenseur& get_b_car() const {return b_car;} ; 

	/// Returns the metric coefficient $N^\varphi$
	const Tenseur& get_nphi() const {return nphi;} ; 

	/** Returns the component $\tilde N^\varphi = N^\varphi r\sin\theta$ 
	 *  of the shift vector
	 */
	const Tenseur& get_tnphi() const {return tnphi;} ; 

	/// Returns the norm of {\tt u\_euler}
	const Tenseur& get_uuu() const {return uuu;} ;		
	
	/// Returns the metric potential $\nu = \ln N$ = {\tt logn\_auto}
	const Tenseur& get_logn() const {return logn;} ;	

	/** Returns the part of the Metric potential $\nu = \ln N$ = {\tt logn}
	 *  generated by the matter terms
	 */
	const Tenseur& get_nuf() const {return nuf;} ;	

	/** Returns the Part of the Metric potential $\nu = \ln N$ = {\tt logn}
	 *  generated by the quadratic terms
	 */
	const Tenseur& get_nuq() const {return nuq;} ;	

	/// Returns the Metric potential $\zeta = \ln(AN)$ = {\tt beta\_auto}
	const Tenseur& get_dzeta() const {return dzeta;} ;	

	/// Returns the Metric potential $\tilde G = (NB-1) r\sin\theta$
	const Tenseur& get_tggg() const {return tggg;} ; 

	/** Returns the vector $W^i$ used in the decomposition of 
	 *  {\tt shift},
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 * NB: {\tt w\_shift} contains the components of $W^i$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping {\tt mp}. 
	 */
	const Tenseur& get_w_shift() const {return w_shift;} ; 
	
	/** Returns the scalar $\chi$ used in the decomposition of 
	 *  {\tt shift} 
	 *  following Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 * NB: {\tt w\_shift} contains the components of $W^i$
	 *      with respect to the Cartesian triad associated with the 
	 *	mapping {\tt mp}. 
	 */
	const Tenseur& get_khi_shift() const {return khi_shift;} ; 

	/** Returns the tensor ${\tilde K_{ij}}$ related to the extrinsic 
	 *  curvature tensor by ${\tilde K_{ij}} = B^{-2} K_{ij}$.
	 *  {\tt tkij} contains the Cartesian components of
	 *  ${\tilde K_{ij}}$. 
	 */
	const Tenseur_sym& get_tkij() const {return tkij;} ; 

	/** Returns the scalar $A^2 K_{ij} K^{ij}$.
	 *  For axisymmetric stars, this quantity is related to the 
	 *  derivatives of $N^\varphi$ by
	 * \begin{equation}
	 *	A^2 K_{ij} K^{ij} = {B^2 \over 2 N^2} \, r^2\sin^2\theta \,  
	 *    \left[ \left( {\partial N^\varphi \over \partial r} \right) ^2
	 *	    + {1\over r^2} \left( {\partial N^\varphi \over 
	 *		    \partial \theta} \right) ^2 \right] \ . 
	 * \end{equation}
	 * In particular it is related to the quantities $k_1$ and $k_2$
	 * introduced by Eqs.~(3.7) and (3.8) of 
	 * Bonazzola et al. {\sl Astron. Astrophys.} {\bf 278}, 421 (1993)
	 * by 
	 * \begin{equation}
	 *	A^2 K_{ij} K^{ij} = 2 A^2 (k_1^2 + k_2^2) \ . 
	 * \end{equation}
	 */
	 const Tenseur& get_ak_car() const {return ak_car;} ; 

    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    /// Save in a file
    
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
	
	/** Description of the stellar surface: returns a 2-D {\tt Itbl} 
	 *	containing the 
	 *	values of the domain index {\it l} on the surface at the 
	 *	collocation points in $(\theta', \phi')$.
	 *	The stellar surface is defined as the location where
	 *	the enthalpy (member {\tt ent}) vanishes.
	 */
	virtual const Itbl& l_surf() const ; 
	
	virtual double mass_b() const ;	    /// Baryon mass
	virtual double mass_g() const ;	    /// Gravitational mass
	virtual double angu_mom() const ;	/// Angular momentum 
	virtual double tsw() const ;		/// Ratio T/W

	/** Error on the virial identity GRV2.
	 *  This indicator is only valid for relativistic computations.
	 */
	virtual double grv2() const ;	

	/** Error on the virial identity GRV3.
	 *  The error is computed as the integral defined
	 *  by Eq. (43) of [Gourgoulhon and Bonazzola, 
	 *  Class. Quantum Grav. {\bf 11}, 443 (1994)] divided by
	 *  the integral of the matter terms.
	 * 
	 *  @param ost output stream to give details of the computation;
	 *		if set to 0x0 [default value], no details will be
	 *		given.
	 *   
	 */
	virtual double grv3(ostream* ost = 0x0) const ;	

	virtual double r_circ() const ;	/// Circumferential radius
	virtual double aplat() const ;	/// Flatening r\_pole/r\_eq
	virtual double z_eqf() const ;	/// Forward redshift factor at equator
	virtual double z_eqb() const ;	/// Backward redshift factor at equator
	virtual double z_pole() const ;	/// Redshift factor at North pole
    
	/** Quadrupole moment.
	 *  The quadrupole moment {\it Q} is defined according to Eq. (7) of
	 *  [Salgado, Bonazzola, Gourgoulhon and Haensel, Astron. Astrophys.
	 *   {\bf 291}, 155 (1994)]. At the Newtonian limit it is related to
	 *  the component ${\bar I}_{zz}$ of the MTW (1973) reduced quadrupole 
	 *  moment ${\bar I}_{ij}$ by: $Q = -3/2 {\bar I}_{zz}$. 
	 *  Note that {\it Q} is the negative of the quadrupole moment defined 
	 *  by Laarakkers and Poisson, Astrophys. J. {\bf 512}, 282 (1999).
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
	 *  {\tt ent}, {\tt ener}, {\tt press}, and {\tt a\_car},  
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  {\tt gam\_euler}, {\tt u\_euler}, {\tt ener\_euler}, {\tt s\_euler}. 
	 * 
	 */
	virtual void hydro_euler() ; 
	
	/** Computes metric coefficients from known potentials. 
	 * 
	 *  The calculation is performed starting from the quantities
	 *  {\tt logn},  {\tt dzeta}, {\tt tggg} and {\tt shift}, 
	 *  which are supposed to be up to date.  
	 *  From these,  the following fields are updated:
	 *  {\tt nnn}, {\tt a\_car},  {\tt bbb} and {\tt b\_car}. 
	 * 
	 */
	void update_metric() ; 
		
	/** Computes {\tt shift} from {\tt w\_shift} and {\tt khi\_shift}
	 *  according to Shibata's prescription 
	 *  [Prog. Theor. Phys. {\bf 101}, 1199 (1999)] :
	 * \begin{equation}
	 *  N^i = {7\over 8} W^i - {1\over 8} 
	 *			\left(\nabla^i\chi+\nabla^iW^kx_k\right)
	 * \end{equation}
	 */
	void fait_shift() ; 
	
	/** Computes {\tt tnphi} and {\tt nphi} from the Cartesian 
	 *   components of the shift, stored in {\tt shift}.
	 */
	void fait_nphi() ; 
		
	/** Computes {\tt tkij} and {\tt ak\_car} from 
	 *  {\tt shift}, {\tt nnn} and {\tt b\_car}.
	 */
	void extrinsic_curvature() ;
	
	/** Computes the coefficient $\lambda$ which ensures that the
	 *	GRV2 virial identity is satisfied.
	 *  $\lambda$ is the coefficient by which one must multiply
	 *  the quadratic source term $\sigma_q$ of the 2-D Poisson equation
	 *	\begin{equation}
	 *		\Delta_2 u = \sigma_m + \sigma_q
	 *	\end{equation}
	 *  in order that the total source does not contain any monopolar term,
	 *  i.e. in order that
	 *  \begin{equation}
	 *		\int_0^{2\pi} \int_0^{+\infty} \sigma(r, \theta)
	 *				\, r \, dr \, d\theta = 0	    \ ,
	 *  \end{equation}
	 *  where $\sigma = \sigma_m + \sigma_q$.
	 *	$\lambda$ is computed according to the formula
	 *  \begin{equation}
	 *		\lambda = - { \int_0^{2\pi} \int_0^{+\infty} \sigma_m(r, \theta)
	 *				\, r \, dr \, d\theta	    \over
	 * 			\int_0^{2\pi} \int_0^{+\infty} \sigma_q(r, \theta)
	 *				\, r \, dr \, d\theta } \ .
	 *  \end{equation}
	 *  Then, by construction, the new source
	 *	$\sigma' = \sigma_m + \lambda \sigma_q$ has a vanishing monopolar
	 *  term.
	 *
	 *	@param sou_m [input] matter source term $\sigma_m$
	 *	@param sou_q [input] quadratic source term $\sigma_q$
	 *  @return	value of $\lambda$
	 */
	static double lambda_grv2(const Cmp& sou_m, const Cmp& sou_q) ;
		
	/** Computes an equilibrium configuration.
	 *  
	 *  @param ent_c  [input] Central enthalpy 
	 *  @param omega0  [input] Requested angular velocity 
	 *			     (if {\tt fact\_omega=1.})
	 *  @param fact_omega [input] 1.01 = search for the Keplerian frequency,
	 *			      1. = otherwise.
	 *  @param nzadapt  [input] Number of (inner) domains where the mapping 
	 *			    adaptation to an iso-enthalpy surface
	 *			    should be performed
	 *  @param ent_limit [input] 1-D {\tt Tbl} of dimension {\tt nzet} which
	 *				defines the enthalpy at the outer boundary
	 *				of each domain
	 *  @param icontrol [input] Set of integer parameters (stored as a
	 *			    1-D {\tt Itbl} of size 8) to control the 
	 *			    iteration: \\
	 *	{\tt icontrol(0) = mer\_max} : maximum number of steps \\
	 *	{\tt icontrol(1) = mer\_rot} : step at which the rotation is 
	 *				      switched on \\
	 *	{\tt icontrol(2) = mer\_change\_omega} : step at which the rotation
	 *			  velocity is changed to reach the final one  \\
	 *	{\tt icontrol(3) = mer\_fix\_omega} :  step at which the final
	 *			    rotation velocity must have been reached  \\
	 *	{\tt icontrol(4) = mer\_mass} : the absolute value of 
	 *			    {\tt mer\_mass} is the step from which the 
	 *			    baryon mass is forced to converge, 
	 *			    by varying the central enthalpy 
	 *			    ({\tt mer\_mass > 0}) or the angular 
	 *			    velocity ({\tt mer\_mass < 0}) \\
	 *	{\tt icontrol(5) = mermax\_poisson} : maximum number of steps in 
	 *				{\tt Map\_et::poisson} \\
	 *	{\tt icontrol(6) = mer\_triax} : step at which the 3-D 
	 *				perturbation is switched on \\
	 *	{\tt icontrol(7) = delta\_mer\_kep} : number of steps
	 *			    after {\tt mer\_fix\_omega} when {\tt omega}
	 *			    starts to be increased by {\tt fact\_omega}
	 *			    to search for the Keplerian velocity
	 * 	 
	 *  @param control [input] Set of parameters (stored as a 
	 *			    1-D {\tt Tbl} of size 7) to control the 
	 *			    iteration: \\
	 *	{\tt control(0) = precis} : threshold on the enthalpy relative 
	 *				change for ending the computation \\
	 *	{\tt control(1) = omega\_ini} : initial angular velocity, 
	 *			    switched on only if {\tt mer\_rot < 0}, 
	 *			    otherwise 0 is used \\ 
	 *	{\tt control(2) = relax} : relaxation factor in the main 
	 *				   iteration \\ 
	 *	{\tt control(3) = relax\_poisson} : relaxation factor in 
	 *				   {\tt Map\_et::poisson}\\ 
	 *	{\tt control(4) = thres\_adapt} :  threshold on dH/dr for 
	 *			    freezing the adaptation of the mapping \\
	 *	{\tt control(5) = ampli\_triax} :  relative amplitude of 
	 *			    the 3-D perturbation \\
	 *	{\tt control(6) = precis\_adapt} : precision for 
	 *			    {\tt Map\_et::adapt}
	 *
	 *  @param mbar_wanted [input] Requested baryon mass (effective only 
	 *				if {\tt mer\_mass > mer\_max})
	 *  @param aexp_mass [input] Exponent for the increase factor of the 
	 *			      central enthalpy to converge to the 
	 *			      requested baryon mass
	 *  @param diff [output]   1-D {\tt Tbl} of size 7 for the storage of 
	 *			    some error indicators : \\
	 *	    {\tt diff(0)} : Relative change in the enthalpy field
	 *			      between two successive steps \\
	 *	    {\tt diff(1)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt nuf} \\  
	 *	    {\tt diff(2)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt nuq} \\  
	 *	    {\tt diff(3)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt dzeta} \\  
	 *	    {\tt diff(4)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt tggg} \\  
	 *	    {\tt diff(5)} : Relative error in the resolution of the
	 *			    equation for {\tt shift} (x comp.) \\  
	 *	    {\tt diff(6)} : Relative error in the resolution of the
	 *			    equation for {\tt shift} (y comp.) \\  
	 */
	virtual void equilibrium(double ent_c, double omega0, double fact_omega, 
			 int nzadapt, const Tbl& ent_limit,
			 const Itbl& icontrol, const Tbl& control,
			 double mbar_wanted, double aexp_mass, 
			 Tbl& diff) ;
	

};




#endif
