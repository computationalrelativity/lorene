/*
 *  Definition of Lorene classes Bhole
 *				 Bhole_binaire
 *
 */

/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


#ifndef __BHOLE_H_ 
#define __BHOLE_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.45  2001/06/28  15:08:16  eric
 * Ajout de la fonction area().
 *
 * Revision 2.44  2001/05/16  11:32:45  phil
 * ajout de init_bhole_seul
 *
 * Revision 2.43  2001/05/07  12:35:28  phil
 * mmodifi commentaires
 *
 * Revision 2.42  2001/05/07  11:43:43  phil
 * *** empty log message ***
 *
 * Revision 2.41  2001/05/07  11:40:50  phil
 * mise a jour des commentaires
 *
 * Revision 2.40  2001/05/07  09:30:37  phil
 * *** empty log message ***
 *
 * Revision 2.39  2001/05/07  09:28:37  phil
 * *** empty log message ***
 *
 * Revision 2.38  2001/05/07  09:11:31  phil
 * *** empty log message ***
 *
 * Revision 2.37  2001/04/27  09:25:11  phil
 * fait_decouple est devenu public
 *
 * Revision 2.36  2001/04/26  12:04:11  phil
 * *** empty log message ***
 *
 * Revision 2.35  2001/04/05  13:42:40  phil
 * *** empty log message ***
 *
 * Revision 2.34  2001/04/05  13:34:29  phil
 * ajout resolution en utilisant phi
 *
 * Revision 2.33  2001/03/22  10:49:34  phil
 * *** empty log message ***
 *
 * Revision 2.32  2001/03/22  10:41:25  phil
 * pleins de modif
 *
 * Revision 2.31  2001/02/28  13:22:31  phil
 * modif vire kk_auto
 *
 * Revision 2.30  2001/02/12  15:37:42  phil
 * ajout calcul de J a linfini
 *
 * Revision 2.29  2001/01/29  14:29:56  phil
 * ajout type rotation
 *
 * Revision 2.28  2001/01/24  10:57:08  phil
 * ajout de set_rayonm
 *
 * Revision 2.27  2001/01/22  09:29:08  phil
 * vire les procedures qui servent pas
 *
 * Revision 2.26  2001/01/10  09:31:05  phil
 * modification de fait_kk_auto (membre de binaire maintenant)
 *
 * Revision 2.25  2001/01/09  13:38:15  phil
 * ajout memebre kk_auto
 *
 * Revision 2.24  2000/12/21  10:47:07  phil
 * retour eventuel a 2.21
 *
 * Revision 2.21  2000/12/20  09:07:56  phil
 * ajout set_statiques
 *
 * Revision 2.20  2000/12/18  16:40:39  phil
 * *** empty log message ***
 *
 * Revision 2.19  2000/12/18  16:38:04  phil
 * ajout convergence vers une masse donnee
 *
 * Revision 2.18  2000/12/15  16:41:28  phil
 * ajout calcul de la separation
 *
 * Revision 2.17  2000/12/14  10:44:28  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.16  2000/12/13  15:35:18  phil
 * ajout calcul bare_masse
 *
 * Revision 2.15  2000/12/04  14:29:37  phil
 * ajout de grad_n_tot
 *
 * Revision 2.14  2000/12/01  16:12:52  phil
 * *** empty log message ***
 *
 * Revision 2.13  2000/12/01  14:16:20  phil
 * *** empty log message ***
 *
 * Revision 2.12  2000/11/24  15:14:44  phil
 * ajout de find_horizon
 *
 * Revision 2.11  2000/11/24  09:57:18  phil
 * ajout calcul masse et moment pour systeme
 *
 * Revision 2.10  2000/11/17  10:03:13  phil
 * ajout de coal
 *
 * Revision 2.9  2000/11/15  18:26:29  phil
 * simplifaction resolution du shift : on bosse a omega bloque
 *
 * Revision 2.8  2000/11/15  12:59:50  phil
 * changement solve_shift_omega
 *
 * Revision 2.7  2000/11/15  09:41:03  phil
 * *** empty log message ***
 *
 * Revision 2.6  2000/11/15  09:39:26  phil
 * ajout viriel
 *
 * Revision 2.5  2000/11/03  12:56:41  phil
 * ajout de const
 *
 * Revision 2.4  2000/10/26  08:20:45  phil
 * ajout verifie_shift
 *
 * Revision 2.3  2000/10/23  09:15:12  phil
 * modif commentaires
 *
 * Revision 2.2  2000/10/20  10:51:13  phil
 * Modif commentaires (minimale !)
 *
 * Revision 2.1  2000/10/20  09:27:28  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/10/20  09:22:04  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

class Tenseur ;
class Tenseur_sym ;
class Map_af ;
class Bhole_binaire ;

/**
 * Black hole.
 * 
 * This class represents a black hole in a comformaly flat approximation. It is 
 * defined with an affine mapping with a nucleus, not used in calculations, and 
 * must have,  at least two shells. The horizon (i.e. the surface where {\it N=0}) 
 * is located at the inner boudary of the innermost shell (i.e. the domain
 * number 1).
 * 
 * It can described either :
 * \begin{itemize}
 * \item {\bf an isolated black hole}. The fields of the companion and the total
 * fields being zero.
 * \item {\bf a black hole in a binary system}. The boost having no direct 
 * meaning in this case.
 * \end{itemize}
 * The tensor $K^{ij}$ is the the extrinsic curavature tensor with the conformal 
 * factor extracted.
 * 
 * The tensor $A^{ij}$ denotes $\nabla^i N^j + \nabla^j N^i-\frac{2}{3}
 * \nabla_k N^k \delta^{ij}$,  that is $2NK^{ij}$.
 *
 * @version #$Id$#
 */
class Bhole {

    // Data :
    protected:
	
	Map_af& mp ;  /// Affine mapping.
	double rayon ; /// Radius of the horizon in LORENE's units.
	double omega ; /// Angular velocity in LORENE's units.
	/**
	 * Cartesian components of the boost in the reference frame.
	 */
	double* boost ;
	double regul ; /// Intensity of the correction on the shift vector. 
	
	Tenseur n_auto ; /// Part of {\it N} generated by the hole.
	Tenseur n_comp ; ///Part of {\it N} generated by the companion hole.
	Tenseur n_tot ; /// Total {\it N}.
	
	Tenseur psi_auto ; /// Part of $\Psi$ generated by the hole.
	Tenseur psi_comp ; /// Part of $\Psi$ generated by the companion hole.
	Tenseur psi_tot ; /// Total $\Psi$.
	
	Tenseur grad_n_tot ; /// Total gradient of {\it N}.
	Tenseur grad_psi_tot ; /// Total gradient of $\Psi$.
	
	Tenseur shift_auto ; /// Part of $\beta^i$ generated by the hole.
	
	Tenseur_sym taij_auto ; /// Part of $A^{ij}$ generated by the hole.
	Tenseur_sym taij_comp ;/// Part of $A^{ij}$ generated by the companion hole.
	/**
	 * Total $A^{ij}$,  which must be zero on the horizon of the
	 * regularisation on the shift has been done.
	 */
	Tenseur_sym taij_tot ;
	
	Tenseur_sym tkij_auto ; /// Auto $K^{ij}$.
	Tenseur_sym tkij_tot ; /// Total $K^{ij}$.
	
	/**
	 * Function used to construct the part of $K^{ij}$ generated by the hole 
	 * from the total $K^{ij}$. Only used for a binary system.
	 * 
	 * Mainly this {\tt Cmp} is 1 around the hole and 0 around the companion
	 * and the sum of {tt decouple} for the hole and his companion is 1 
	 * everywhere.
	 */
	Cmp decouple ;
	
    //Constructors :
    public:
	/**
	 * Standard constructor. All the fields are set to zero, except
	 * {\tt rayon},  which value is deducted from the {\tt mapping}
	 */
	Bhole (Map_af& mapping) ;
	Bhole (const Bhole&) ;	/// Constructor by copy
	Bhole (Map_af&, FILE*) ; /// Constructor from a {\tt Map\_af} and a file
	~Bhole() ; /// Destructor
	
    public:
    //sauvegarde
	void sauve (FILE* fich) const ; /// Write on a file
	
    public:
	void operator= (const Bhole&) ; /// Affectation
	
	const Map_af& get_mp() const {return mp;} ; /// Returns the mapping.
	/**
	 * Returns the radius of the horizon.
	 */
	double get_rayon() const {return rayon;} ;
	
	/**
	 * Sets the radius of the horizon to {\tt ray}.
	 */
	void set_rayon(double ray) {rayon = ray ;} ;
	
	/**
	 * Returns the angular velocity.
	 */
	double get_omega() const {return omega;} ;
	/**
	 * Sets the angular velocity to {\tt ome}.
	 */
	void set_omega(double ome) {omega = ome ;} ;
	
	/**
	 * Returns the cartesian components of the boost with respect to the 
	 * reference frame.
	 */
	double* get_boost () const {return boost;} ;
	void set_boost (double vx, double vy, double vz) {
	    boost[0] = vx ; boost[1] = vy ; boost[2] = vz ;
	}
	
	/**
	 * Returns the intensity of the regularisation on $\vec{N}$.
	 */
	double get_regul() const {return regul;} ;
	
	/**
	 * Returns the part of {\it N} generated by the hole.
	 */
	const Tenseur& get_n_auto() const {return n_auto ;} ;
	/**
	 * Returns the part of {\it N} generated by the companion hole.
	 */
	const Tenseur& get_n_comp() const {return n_comp ;} ;
	/**
	 * Returns the total {\it N}.
	 */
	const Tenseur& get_n_tot() const {return n_tot ;} ;
	
	/**
	 * Returns the part of $\Psi$ generated by the hole.
	 */
	const Tenseur& get_psi_auto() const {return psi_auto ;} ;
	/**
	 * Returns the part of $\Psi$ generated by the companion hole.
	 */
	const Tenseur& get_psi_comp() const {return psi_comp ;} ;
	/**
	 * Returns the total $\Psi$.
	 */
	const Tenseur& get_psi_tot() const {return psi_tot ;} ;
	
	/**
	 * Returns the gradient of $\Psi$.
	 */
	const Tenseur& get_grad_psi_tot() const {return grad_psi_tot ;} ;
	
	/**
	 * Returns the gradient of {\it N}.
	 */
	const Tenseur& get_grad_n_tot() const {return grad_n_tot ;} ;
	
	
	/**
	 * Returns the part of $\beta^i$ generated by the hole.
	 */
	const Tenseur& get_shift_auto() const {return shift_auto ;} ;
	
	/**
	 * Returns the part of $A^{ij}$ generated by the hole.
	 */
	const Tenseur& get_taij_auto() const {return taij_auto ;} ;
	/**
	 * Returns the part of $A^{ij}$ generated by the companion hole.
	 */
	const Tenseur& get_taij_comp() const {return taij_comp ;} ;
	/**
	 * Returns the total $A^{ij}$.
	 */
	const Tenseur& get_taij_tot() const {return taij_tot ;}
	
	/**
	 * Returns the total $K^{ij}$.
	 */
	const Tenseur& get_tkij_tot() const {return tkij_tot ;}
	/**
	 * Returns the part of $K^{ij}$ generated by the hole.
	 */
	const Tenseur& get_tkij_auto() const {return tkij_auto ;}
	/**
	 * Returns the function used to construct {\tt tkij\_auto} from {\tt tkij\_tot}.
	 */
	const Cmp get_decouple() const {return decouple ;}
    public:
	
	/**
	 * Imports the part of {\it N} due to the companion hole {\tt comp}. The 
	 * total {\it N} is then calculated.
	 * 
	 * It also imports the gradient of {\it N} and construct the total $\nabla N$.
	 */
	void fait_n_comp (Bhole comp) ;
	
	/**
	 * Imports the part of $\Psi$ due to the companion hole {\tt comp}. The 
	 * total $\Psi$ is then calculated.
	 * 
	 * It also imports the gradient of $\Psi$ and construct the total $\nabla \Psi$.
	 */
	void fait_psi_comp (Bhole comp) ;
	
	/**
	 * Calculates the part of $A^{ij}$ generated by {\tt shift\_auto}.
	 */
	void fait_taij_auto () ;
	

	/**
	 * Corrects {\tt shift\_auto} in such a way that the total $A^{ij}$ is 
	 * equal to zero in the horizon,  which should ensure the regularity 
	 * of $K^{ij}$,  using the stored values of the boost and the angular 
	 * velocity.
	 * 
	 * {\bf WARNING : } this should only be used for a binary black hole.
	 * 
	 * @param comp [input] : the part of $\beta^i$ generated by the companion 
	 * hole.
	 */
	void regularise_shift (Tenseur& comp) ;
	
	/**
	 * Computes the viriel error, that is the difference between the ADM and the Komar 
	 * masses,  calculated by the asymptotic behaviours of respectively $\Psi$ and {\it N}.
	 * 
	 * {\bf WARNING} this should only be used for an isolated black hole.
	 */
	 double viriel_seul () const ;
	 
	/**
	 * Sets the values of the fields to :
	 * \begin {itemize}
	 * \item {\tt n\_auto} $= \frac{1}{2}-2\frac{a}{r}$
	 * \item {\tt n\_comp}  $= \frac{1}{2}$
	 * \item {\tt psi\_auto} $= \frac{1}{2}+\frac{a}{r}$
	 * \item {\tt psi\_comp}  $= \frac{1}{2}$
	 * \end{itemize}
	 * {\it a} being the radius of the hole, the other fields being set to zero.
	 */
	void init_bhole () ;
	
	/**
	 *  Set the inital values to those of Kerr
	 * @param masse [input] : ADM mass in LORENE's units.
	 * @param moment [input] : $\frac{J}{M}=a$, in LORENE's units.
	 * 
	 * The radius $r_0$ of {\tt *this} is supposed to be equal to 
	 * $\frac{1}{2}\sqrt{M^2-a^2}$.
	 * 
	 * The fields have then the following values~:
	 * \begin{itemize}
	 * \item $\Omega = \frac{a}{2M\sqrt{M^2-a^2}}$
	 * 
	 * \item $\Psi^4 = 
	 * 1+\frac{2M}{r}+\frac{3M^2+a^2\cos^2\theta}{2r^2}+\frac{2Mr_0^2}{r^3}+
	 * \frac{r_0^4}{r^4}$.
	 * \item $N = \left(1-\frac{2MR}{\Sigma}+\frac{4a^2M^2R^2\sin^2\theta}
	 * {\Sigma^2\left(R^2+a^2\right)+2a^2\Sigma R\sin^2\theta}\right)^{\frac{1}{2}}$
	 * \item $N^\phi = \frac{2aMR}{\Sigma\left(R^2+a^2\right)+2a^2MR\sin^2\theta}$
	 * \end {itemize}
	 * where~:
	 * \begin{itemize}
	 * \item $R = r+\frac{M^2-a^2}{4r}+M$.
	 * \item $\Sigma = R^2+a^2-2MR$
	 * \item $r_0$ is equal to {\tt rayon}.
	 * \end{itemize}
	 */
	void init_kerr (double masse, double moment) ;
	
	/**
	 * Solves the equation for {\it N}~:
	 *    \begin{equation}
	 *\Delta N = -\frac{2}{\Psi}\nabla_i \Psi \nabla^i N 
	 * + N \Psi^4 K_{ij} K^{ij}
	 * \end{equation}
	 * with the condition that {\tt N}=0 on the horizon.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * {\bf WARNING} this should only be used for an isolated black hole.
	 */
	void solve_lapse_seul (double relax) ;
	
	/**
	 * Solves the equation for $\Psi$~:
	 *    \begin{equation}
	 (\Delta \Psi = - \frac{\Psi^5} {8} K_{ij}K^{ij}
	 *  \end{equation}
	 * with the condition that $\partial_r \Psi=-\frac{1}{2 {\tt rayon}}
	 * f\left(\theta, \phi\right)$ on the horizon, where {\tt f} is the value of $\Psi$ on 
	 * the horizon at the preceeding step.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * {\bf WARNING} this should only be used for an isolated black hole.
	 */
	 
	void solve_psi_seul (double relax) ;
	
	/**
	 * Solves the equation for $\vec{\beta}$~:
	 *    \begin{equation}
	 *\Delta \beta^i +\frac{1}{3} \nabla^i \left(\nabla_j \beta^j\right)
	 * = -6A^{ij}\frac{\nabla_j \Psi}{\Psi} + 2 K^{ij}\nabla_j N
	 * \end{equation}
	 * with $\vec{\beta} = -\Omega \vec{m} - \vec{V}$ on the horizon, 
	 * $\vec{V}$ being the boost and $\vec{m} = \frac{\partial}{\partial \phi}$.
	 * The solution is solved using Oohara-scheme and an iteration.
	 * @param precis [input] : parameter for the Oohara-solver which is an iterative 
	 scheme.
	 * @param relax [input] : the relaxation parameter.
	 * 
	 * {\bf WARNING} this should only be used for an isolated black hole.
	 */
	void solve_shift_seul (double precis, double relax) ;
	
	/**
	 * Corrects the shift in the innermost shell, so that it remains $
	 * {\mathcal{C}}^2$ and that $A^{ij}$ equals zero on the horizon.
	 * 
	 * {\tt regul} is then,  the relative difference between the shift before
	 * and after the regularisation.
	 * 
	 * {\bf WARNING} this should only be used for an isolated black hole.
	 */
	void regularise_seul () ;
	
	/**
	 * Calculates the total $K^{ij}$. The regularisation of the shift must be done
	 * before to ensure regularity.
	 */
	void fait_tkij() ;
	
	/// Computes the area of the throat. 
	double area() const ; 
	
	/**
	 *  Calculates the ADM mass of the black hole using :
	 * $M = -\frac{1}{2 \pi} \oint_{\infty} \nabla^i \Psi {\mathrm d} S_i$.
	 * 
	 * {\bf WARNING} this should only be used for an isolated black hole.
	 */
	double masse_adm_seul () const ;
	
	/**
	 *  Calculates the Komar mass of the black hole using :
	 * $M = \frac{1}{4 \pi} \oint_{\infty} \nabla^i N {\mathrm d} S_i$.
	 * 
	 * {\bf WARNING} this should only be used for an isolated black hole.
	 */
	double masse_komar_seul() const ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula at infinity :
	 * $ J = -\frac{1}{8 \pi} \oint_{\infty} K_j^i m^j {\mathrm d}S_i$
	 * where $\vec{m} = \frac{\partial}{\partial \phi}$.
	 * 
	 *{\bf WARNING} It supposes that the boost is zero and should only be 
	 * used for an isolated black hole..
	 */
	double moment_seul_inf() const ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula on the horizon :
	 * $ J = -\frac{1}{8 \pi} \oint_{H} \Psi^6 K_j^i m^j {\mathrm d}S_i$
	 * where $\vec{m} = \frac{\partial}{\partial \phi}$ and {\it H} denotes the horizon.
	 * 
	 *{\bf WARNING} It supposes that the boost is zero and should only be 
	 * used for an isolated black hole..
	 */
	double moment_seul_hor() const ;
	
	/**
	 * Initiates the black hole for a resolution with $\Phi = \log \Psi$.
	 */
	void init_bhole_phi () ;
	
	/**
	 * Initiates for a single the black hole.
	 * 
	 * {\bf WARNING} It supposes that the boost is zero and should only be 
	 * used for an isolated black hole..
	 */
	void init_bhole_seul () ;
	
	friend class Bhole_binaire ; /// Binary black hole system.

} ;

/**
 * Binary black holes system.
 * 
 * This class is intended for dealing with binary black holes configurations 
 * in the conformaly flat approximation.
 */
class Bhole_binaire {
    
    // data :
    private:
	// les deux trous noirs.
	Bhole hole1 ;	/// Black hole one
	Bhole hole2 ;	/// Black hole two
	
	// Tableau sur les deux trous.
	Bhole* holes[2] ; /// Array on the black holes
	
	double omega ;	/// Angular velocity
	
    public:
	// constructeurs & destructeur
	Bhole_binaire(Map_af& mp1, Map_af& mp2) ;   /// Standard constructor
	Bhole_binaire(const Bhole_binaire& ) ;	/// Copy
	~Bhole_binaire() ;  /// Destructor
	
    public:
	// trucs pour modifier
	void operator=(const Bhole_binaire&) ; /// Affectation operator
	
	/**
	 * Read/write of a component of the system. {\tt i} must be equal to
	 * 1 or 2.
	 */
	Bhole& set(int i) 
	    { assert( (i==1) || (i==2) ); 
	      return *holes[i-1] ;} ; 
	/**
	 * Sets the orbital velocity to {\tt ome}
	 */
	void set_omega(double ome) {omega = ome ; 
			     hole1.set_omega (ome) ;
			     hole2.set_omega (ome) ;} ;
	
    public:
	// trucs pour lire :
	/**
	 * Read only of a component of the system. {\tt i} must be equal to
	 * 1 or 2.
	 */
	const Bhole& operator()(int i) const 
	    { assert( (i==1) || (i==2) ); 
	      return *holes[i-1] ;} ;
	      
	double get_omega() const {return omega; } ; /// Returns the angular velocity 
    
	/**
	 * Initialisation of the system. Each hole is set close to a Schwarzschild one 
	 * and the parts of the fields generated by 
	 * the companion are calculated.
	 * 
	 * The angular velocity is set to zero.
	 */
	void init_bhole_binaire() ;
	
	/**
	 * Computes the viriel error, that is the difference between the ADM and the Komar 
	 * masses,  calculated by the asymptotic behaviours of respectively $\Psi$ and {\it N}.
	 */
	double viriel() const ;
	/**
	 * Solves the equation for the lapse~:
	 *    \begin{equation}
	 * \Delta N_a = -\frac{2}{\Psi}\nabla_i \Psi_a \nabla^i N + N \Psi^4 K_{ij}K_a^{ij}
	 *  \end{equation}
	 * The fields are the total values excpet those with subscript $_a$, which are 
	 * the fields generated by each holes ({\it a} = 1, 2). The boundary conditions are 
	 * such that {\it N}=0 on both horizons.
	 * The companion contributions are then obtained.
	 * @param precis [input] : precision,  for the boudary conditions are
	 * obtained by iteration.
	 * @param relax [input] : relaxation parameter.
	 */
	void solve_lapse (double precis, double relax) ;
	
	/**
	 * Solves the equation for the conformal factor~:
	 *    \begin{equation}
	 *\Delta \Psi_a = -\frac {\Psi^5}{8} K_{ij}K_a^{ij}
	 * \end{equation}
	 * The fields are the total values excpet those with subscript $_a$, which are 
	 * the fields generated by each holes ({\it a} = 1, 2). The boundary conditions are such that 
	 * $\partial_r \Psi = -\frac{1}{2 r_0} f\left(\theta, \phi\right)$ on both horizons,  $r_0$ being 
	 * the radii of those horizons and {\it f} the value of $\Psi$ on the 
	 * horizons at the preceeding step. The companion contributions are then 
	 * obtained.
	 * @param precis [input] : precision,  for the boudary conditions are being 
	 * obtained by iteration.
	 * @param relax [input] : relaxation parameter.
	 */
	void solve_psi (double precis, double relax) ;
   
	/**
	 * Solves the equation for the shift, using the Oohara-Nakarmure scheme~:
	 *    \begin{equation}
	 *\Delta \beta_a^i +\frac{1}{3} \nabla^i \left(\nabla_j \beta_a^j\right) = 
	 * -6A^{ij}\frac{\nabla_i \Psi_a}{\Psi} + 2K^{ij}\nabla_j N_a
	 * \end{equation}
	 * using the current $\Omega$ for the boudary conditions~:
	 * $\vec{N} = -\Omega \vec{m}$ on both horizons. 
	 * The fields are the total values excpet those with subscript $_a$, which are 
	 * the fields generated by each holes ({\tt a} = 1, 2). The companion contributions are then 
	 * obtained.
	 * @param precis [input] : precision for the solver, the boundary 
	 * conditions and the inversion of the operator being 
	 * obtained by iteration.
	 * @param relax [input] : relaxation parameter.
	 */
	void solve_shift (double precis, double relax) ;
	
	/**
	 * Calculation af the extrinsic curvature tensor.
	 * 
	 * The regularisation of the shifts must have been done. All the 
	 * contributions of $A^{ij}$ are then calculated and the total tensor 
	 * must be zero on both horizons. The computation is done to avoid every singularity 
	 * close to the horizons (division by {\it N}=0) for it is done in the coefficient space 
	 * in the two regions surrounding the holes.
	 */
	void fait_tkij () ;
	
	/**
	 * Calculates {tt decouple} which is used to obtain {\tt tkij\_auto} by the formula : 
	 * {\tt tkij\_auto} = {\tt decouple} * {\tt tkij\_tot}.
	 * (see the membre {tt Cmp decouple} for more precisions about its value).
	 * 
	 */
	void fait_decouple () ;
    
    public:
	 /**
	  * Initialize the systeme to Misner Lindquist solution, that is solving for {\it N} and 
	  * $Psi$ in the case $\Omega = 0$.
	  * @param precis [input] : precision for the convergence (on {\it N}).
	  * @param relax [input] : relaxation parameter.
	  */
	 
	 void set_statiques (double precis, double relax) ;
	 
	 /**
	  * Solves the equation for a particular angular velocity, the systeme being 
	  * initialized to Misner-Lindquist solution.
	  * @param angu [input] : angular velocity used for the boundary condition on 
	  * $\vec{\beta}$.
	  * @param precis [input] : precision for the convergence (on $\beta$).
	  * @param relax [input] : relaxation parameter.
	  * @param nbre_ome [input] : number of intermediates velocities to go from 0 to 
	  * {\tt omega}, typically 10.
	  * @param sortie [input] : flag for the output on files (0 no output files).
	  * @returns : the virial error.
	  */
	  double coal (double angu, double precis, double relax,
	   	    double nbre_ome,const int sortie = 0) ;
			    
	/**
	 *  Calculates the ADM mass of the system using :
	 * $M = -\frac{1}{2 \pi} \oint_{\infty} \nabla^i \Psi {\mathrm d} S_i$.
	 */
	double adm_systeme() const ;
	
	/**
	 *  Calculates the Komar mass of the system using :
	 * $M = \frac{1}{4 \pi} \oint_{\infty} \nabla^i N {\mathrm d} S_i$.
	 */
	double komar_systeme() const ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula at infinity :
	 * $ J = -\frac{1}{8 \pi} \oint_{\infty} K_j^i m^j {\mathrm d}S_i$
	 * where $\vec{m} = \frac{\partial}{\partial \phi}$.
	 */
	double moment_systeme_inf() ;
	
	/**
	 * Calculates the angular momentum of the black hole using the formula on the horizon :
	 * $ J = -\sum_{a= 1, 2} \frac{1}{8 \pi} \oint_{H_a} \Psi^6 K_j^i m^j {\mathrm d}S_i$
	 * where $\vec{m} = \frac{\partial}{\partial \phi}$ and $H_a$ denotes the horizon 
	 * of the hole {\it a}.
	 * 
	 */
	double moment_systeme_hor() const ;
	
	/**
	 * Calculation of the proper distance between the two spheres of inversion, 
	 * along the x axis.
	 * @param nr [input] : number of points used for the calculation.
	 */
	 double distance_propre(const int nr = 65) const ;
	
	/**
	 * Solve the equation for the logarithm of $\Psi$.
	 */
	void solve_phi (double precision, double relax) ;
	/**
	 * Initiates the system for a resolution using the logarithm of $\Psi$.
	 */
	void init_phi() ;

} ;

#endif
