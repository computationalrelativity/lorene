/*
 *  Prototypes of non class-member functions
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
 *   Copyright (c) 2002-2003 Jerome Novak
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


#ifndef	__PROTO_H_
#define	__PROTO_H_


/*
 * $Id$
 * $Log$
 * Revision 1.11  2003/09/16 13:07:40  j_novak
 * New files for coefficient trnasformation to/from the T_LEG_II base.
 *
 * Revision 1.10  2003/09/16 08:53:05  j_novak
 * Addition of the T_LEG_II base (odd in theta, only for odd m) and the
 * transformation functions from and to T_SIN_P.
 *
 * Revision 1.9  2003/06/18 08:45:26  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
 * Revision 1.8  2003/02/13 16:40:24  p_grandclement
 * Addition of various things for the Bin_ns_bh project, non of them being
 * completely tested
 *
 * Revision 1.7  2002/11/12 17:45:19  j_novak
 * Added transformation function for T_COS basis.
 *
 * Revision 1.6  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.5  2002/05/11 12:39:08  e_gourgoulhon
 * Added declaration of som_tet_cossin_si.
 *
 * Revision 1.4  2002/05/05 16:24:48  e_gourgoulhon
 * Added som_tet_cossin_sp
 *
 * Revision 1.3  2002/01/03 15:30:27  j_novak
 * Some comments modified.
 *
 * Revision 1.2  2002/01/02 14:07:56  j_novak
 * Dalembert equation is now solved in the shells. However, the number of
 * points in theta and phi must be the same in each domain. The solver is not
 * completely tested (beta version!).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.69  2001/05/07  09:11:26  phil
 * *** empty log message ***
 *
 * Revision 2.68  2001/04/03  12:41:23  phil
 * modification de itemax dans separation
 *
 * Revision 2.67  2001/03/22  10:40:13  phil
 * modification prototypage se separation
 *
 * Revision 2.66  2001/03/02  10:21:26  phil
 * *** empty log message ***
 *
 * Revision 2.65  2001/03/02  10:18:47  phil
 * modif parametrage separation
 *
 * Revision 2.64  2001/02/28  11:23:30  phil
 * ajout separation
 *
 * Revision 2.63  2001/01/29  14:30:10  phil
 * ajout type rotation
 *
 * Revision 2.62  2000/12/13  15:42:14  phil
 * ajout des trucs relatifs a Lindquist
 *
 * Revision 2.61  2000/12/04  13:33:28  novak
 * *** empty log message ***
 *
 * Revision 2.60  2000/12/04 13:29:08  novak
 * Added prototypes for the dalembertian
 *
 * Revision 2.59  2000/10/19 10:07:49  phil
 * ajout de regle
 *
 * Revision 2.58  2000/10/19  09:35:44  phil
 * *** empty log message ***
 *
 * Revision 2.57  2000/10/04  14:40:34  eric
 * *** empty log message ***
 *
 * Revision 2.56  2000/09/29  14:02:06  eric
 * *** empty log message ***
 *
 * Revision 2.55  2000/09/28  10:01:36  eric
 * *** empty log message ***
 *
 * Revision 2.54  2000/09/08  16:26:26  eric
 * *** empty log message ***
 *
 * Revision 2.53  2000/09/08  16:07:43  eric
 * *** empty log message ***
 *
 * Revision 2.52  2000/09/07  15:07:40  eric
 * *** empty log message ***
 *
 * Revision 2.51  2000/09/06  13:59:53  eric
 * *** empty log message ***
 *
 * Revision 2.50  2000/06/06  12:42:55  phil
 * ajout de Cmp division_xpun (const Cmp&,int)
 *
 * Revision 2.49  2000/05/22  13:33:15  phil
 * ajout des trucs pour poisson avec dzpuis == 3
 *
 * Revision 2.48  2000/04/03  17:01:01  phil
 * ajout de sxpun_1d
 *
 * Revision 2.47  2000/03/16  16:28:30  phil
 * *** empty log message ***
 *
 * Revision 2.46  2000/03/09  13:52:55  phil
 * *** empty log message ***
 *
 * Revision 2.45  2000/03/09  13:42:34  phil
 * vire les trucs relatifs a Poisson compacts
 *
 * Revision 2.44  2000/03/06  10:27:07  eric
 * Ajout des protos som_*_symy et som_*_asymy.
 *
 * Revision 2.43  2000/01/20  14:07:59  phil
 * vire poisson_vect et xksk
 *
 * Revision 2.42  1999/12/15  09:41:52  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

class Param ;
class Tbl ;
class Mtbl ;
class Mtbl_cf ;
class Map_af ;
class Matrice ;
class Valeur ;
class Base_val ;
class Cmp ;
class Tenseur ;

// Routines calcul de coefficients
// -------------------------------
double* cheb_ini(const int) ;
double* chebimp_ini(const int) ;
void four1d(const int, double* ) ;
void chebyf1d(const int, double* ) ; 
void chebyr1d(const int, double* ) ;
void cfpcossin(const int* ,const int* ,  double* ) ;
void cfpcossini(const int* ,const int* ,  double* ) ;
void cftcos(const int*, const int*, double*, const int*, double*) ;
void cftcosp(const int*, const int*, double*, const int*, double*) ;
void cftcosi(const int*, const int*, double*, const int*, double*) ;
void cftsinp(const int*, const int*, double*, const int*, double*) ;
void cftsini(const int*, const int*, double*, const int*, double*) ;
void cftcossincp(const int*, const int*, double*, const int*, double*) ;
void cftcossinsi(const int*, const int*, double*, const int*, double*) ;
void cftcossinsp(const int*, const int*, double*, const int*, double*) ;
void cftcossinci(const int*, const int*, double*, const int*, double*) ;
void cftlegp(const int*, const int*, double*, const int*, double*) ;
void cftlegpp(const int*, const int*, double*, const int*, double*) ;
void cftlegi(const int*, const int*, double*, const int*, double*) ;
void cftlegip(const int*, const int*, double*, const int*, double*) ;
void cftlegpi(const int*, const int*, double*, const int*, double*) ;
void cftlegii(const int*, const int*, double*, const int*, double*) ;
void cfrcheb(const int*, const int*, double*, const int*, double*) ;
void cfrchebp(const int*, const int*, double*, const int*, double*) ;
void cfrchebi(const int*, const int*, double*, const int*, double*) ;
void cfrchebpimp(const int*, const int*, double*, const int*, double*) ;
void cfrchebpimi(const int*, const int*, double*, const int*, double*) ;
void cipcossin(const int* , const int* , const int* , double* , double* ) ;
void cipcossini(const int* , const int* , const int* , double* , double* ) ;
void citcos(const int*, const int*, double*, const int*, double*) ;
void citcosp(const int*, const int*, double*, const int*, double*) ;
void citcosi(const int*, const int*, double*, const int*, double*) ;
void citsinp(const int*, const int*, double*, const int*, double*) ;
void citsini(const int*, const int*, double*, const int*, double*) ;
void citcossincp(const int*, const int*, double*, const int*, double*) ;
void citcossinsi(const int*, const int*, double*, const int*, double*) ;
void citcossinsp(const int*, const int*, double*, const int*, double*) ;
void citcossinci(const int*, const int*, double*, const int*, double*) ;
void citlegp(const int*, const int*, double*, const int*, double*) ;
void citlegpp(const int*, const int*, double*, const int*, double*) ;
void citlegi(const int*, const int*, double*, const int*, double*) ;
void citlegip(const int*, const int*, double*, const int*, double*) ;
void citlegpi(const int*, const int*, double*, const int*, double*) ;
void citlegii(const int*, const int*, double*, const int*, double*) ;
void circheb(const int*, const int*, double*, const int*, double*) ;
void circhebp(const int*, const int*, double*, const int*, double*) ;
void circhebi(const int*, const int*, double*, const int*, double*) ;
void circhebpimp(const int*, const int*, double*, const int*, double*) ;
void circhebpimi(const int*, const int*, double*, const int*, double*) ;
double* legendre(int , int ) ;
double* legendre_norm(int , int ) ;
double* mat_cossincp_legp(int, int) ;   
double* mat_cossinci_legi(int, int) ;   
double* mat_cosp_legpp(int, int) ;   
double* mat_cosi_legip(int, int) ;   
double* mat_sini_legpi(int, int) ;   
double* mat_sinp_legii(int, int) ;   
double* mat_legp_cossincp(int,  int) ;
double* mat_legi_cossinci(int,  int) ;
double* mat_legpp_cosp(int, int) ;   
double* mat_legip_cosi(int, int) ;   
double* mat_legpi_sini(int, int) ;   
double* mat_legii_sinp(int, int) ;   
void chb_cossincp_legp(const int* , const double* , double* ) ;
void chb_legp_cossincp(const int* , const double* , double* ) ;
void chb_cosp_legpp(const int* , const double* , double* ) ;
void chb_legpp_cosp(const int* , const double* , double* ) ;
void chb_cosi_legip(const int* , const double* , double* ) ;
void chb_legip_cosi(const int* , const double* , double* ) ;
void chb_sini_legpi(const int* , const double* , double* ) ;
void chb_legpi_sini(const int* , const double* , double* ) ;
void chb_cossinci_legi(const int* , const double* , double* ) ;
void chb_legi_cossinci(const int* , const double* , double* ) ;
void chb_sinp_legii(const int* , const double* , double* ) ;
void chb_legii_sinp(const int* , const double* , double* ) ;

double int1d_chebp(int, const double* ) ;

// Routines calcul de coef inverse
void cipcossin(const int* , const int* , const int* , double* , double* ) ;
void citcosp(const int*, const int*, double*, const int*, double*) ;
void citcosi(const int*, const int*, double*, const int*, double*) ;
void citsinp(const int*, const int*, double*, const int*, double*) ;
void citsini(const int*, const int*, double*, const int*, double*) ;
void citcossincp(const int*, const int*, double*, const int*, double*) ;
void citcossinsi(const int*, const int*, double*, const int*, double*) ;
void citcossinsp(const int*, const int*, double*, const int*, double*) ;
void citcossinci(const int*, const int*, double*, const int*, double*) ;
void citlegp(const int*, const int*, double*, const int*, double*) ;
void citlegpp(const int*, const int*, double*, const int*, double*) ;
void citlegi(const int*, const int*, double*, const int*, double*) ;
void circheb(const int*, const int*, double*, const int*, double*) ;
void circhebp(const int*, const int*, double*, const int*, double*) ;
void circhebi(const int*, const int*, double*, const int*, double*) ;
void circhebpimp(const int*, const int*, double*, const int*, double*) ;
void circhebpimi(const int*, const int*, double*, const int*, double*) ;

// Routines calculant la matrice du laplacien
Matrice _laplacien_mat_pas_prevu(int, int, double, int) ;
Matrice _laplacien_mat_r_chebp(int, int, double, int) ;
Matrice _laplacien_mat_r_chebi(int, int, double, int) ;
Matrice _laplacien_mat_r_chebu(int, int, double, int) ;
Matrice _laplacien_mat_r_chebu_deux(int, int) ;
Matrice _laplacien_mat_r_chebu_trois(int, int) ;
Matrice _laplacien_mat_r_chebu_quatre(int, int) ;
Matrice _laplacien_mat_r_cheb(int, int, double, int) ;
Matrice laplacien_mat(int , int , double , int, int ) ;

//Routines de passage a bande versions Matrice et Tbl
Matrice _cl_pas_prevu (const Matrice&, int, double, int) ;
Matrice _cl_r_cheb (const Matrice&, int, double, int) ;
Matrice _cl_r_chebi (const Matrice&, int, double, int) ;
Matrice _cl_r_chebu (const Matrice&, int, double, int) ;
Matrice _cl_r_chebu_quatre (const Matrice&, int) ;
Matrice _cl_r_chebu_trois (const Matrice&, int) ;
Matrice _cl_r_chebu_deux (const Matrice&, int) ;
Matrice _cl_r_chebp (const Matrice&, int, double, int) ;
Matrice combinaison (const Matrice&, int, double, int, int) ;

Tbl _cl_pas_prevu (const Tbl&, int) ;
Tbl _cl_r_cheb (const Tbl&, int) ;
Tbl _cl_r_chebi (const Tbl&, int) ;
Tbl _cl_r_chebu (const Tbl&, int) ;
Tbl _cl_r_chebu_deux (const Tbl&) ;
Tbl _cl_r_chebu_trois (const Tbl&) ;
Tbl _cl_r_chebu_quatre (const Tbl&) ;
Tbl _cl_r_chebp (const Tbl&, int) ;
Tbl combinaison (const Tbl&, int, int) ;


// Routines de preparation du laplacien inversible
Matrice _prepa_nondege_pas_prevu(const Matrice &, int , double, int) ;
Matrice _prepa_nondege_r_cheb (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebp (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebi (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebu (const Matrice&, int, double, int) ;
Matrice _prepa_nondege_r_chebu_deux (const Matrice&, int) ;
Matrice _prepa_nondege_r_chebu_trois (const Matrice&, int) ;
Matrice _prepa_nondege_r_chebu_quatre (const Matrice&, int) ;
Matrice prepa_nondege (const Matrice&, int, double, int, int) ;

//Routines de calcul de la solution particuliere
Tbl _solp_pas_prevu(const Matrice&, const Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_cheb (const Matrice&, const  Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebp (const Matrice&, const  Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebi (const Matrice&, const Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebu (const Matrice&, const Matrice&, double, double, const Tbl&, int) ;
Tbl _solp_r_chebu_deux (const Matrice&, const Matrice&, const Tbl&) ;
Tbl _solp_r_chebu_trois (const Matrice&, const Matrice&, double, const Tbl&) ;
Tbl _solp_r_chebu_quatre (const Matrice&, const Matrice&, double, const Tbl&) ;
Tbl solp (const Matrice&, const Matrice&, double, double, const Tbl&, int, int) ;

//Routines de calcul des solutions homogenes
Tbl _solh_pas_prevu (int, int, double) ;
Tbl _solh_r_cheb (int, int, double) ;
Tbl _solh_r_chebp (int, int, double) ;
Tbl _solh_r_chebi (int, int, double) ;
Tbl _solh_r_chebu (int, int, double) ;
Tbl solh (int, int, double, int) ;

//Routines de derivations version 1d
void _d2sdx2_1d_pas_prevu(int, double*, double* ) ;
void _d2sdx2_1d_r_chebu(int, double*, double* ) ;
void _d2sdx2_1d_r_cheb(int, double*, double* ) ;
void _d2sdx2_1d_r_chebp(int, double*, double* ) ;
void _d2sdx2_1d_r_chebi(int, double*, double * ) ;
void d2sdx2_1d(int, double** , int) ;

void _multx_1d_pas_prevu(int, double*, double* ) ;
void _multx_1d_r_cheb(int, double*, double* ) ;
void multx_1d(int, double **, int) ;

void _sx2_1d_pas_prevu(int, double*, double*) ;
void _sx2_1d_identite(int, double*, double*) ;
void _sx2_1d_r_chebp(int, double*, double*) ;
void _sx2_1d_r_chebi(int, double*, double*) ;
void _sxm12_1d_r_chebu(int, double *, double*) ;
void sx2_1d(int, double**, int) ;

void _sxdsdx_1d_pas_prevu(int, double*, double*) ;
void _dsdx_1d_r_cheb(int, double*, double*) ;
void _sxdsdx_1d_r_chebi(int, double*, double*) ;
void _sxdsdx_1d_r_chebp(int, double*, double*) ;
void sxdsdx_1d(int, double** , int) ;

//Routines de derivations (pour sol_dalembert)
void _dsdx_r_chebp(Tbl *, int &) ;
void _dsdx_r_chebi(Tbl *, int &) ;

// Resolution de l'equation de Poisson
int nullite_plm_sym (int, int, int, int) ;
int nullite_plm_nonsym (int, int, int, int) ;
int nullite_plm_nonsym_anti (int, int, int, int) ;
int nullite_plm (int, int, int, int, Base_val) ;

void donne_lm_sym (int, int, int, int, int&, int&, int&) ;
void donne_lm_nonsym (int, int, int, int, int&, int&, int&) ;
void donne_lm_nonsym_anti (int, int, int, int, int&, int&, int&) ;
void donne_lm (int, int, int, int, Base_val, int&, int&, int&) ;


// Les sommations en r :
void som_r_pas_prevu
    (double*, const int, const int, const int, const double, double*) ;
void som_r_cheb
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebi
    (double*, const int, const int, const int, const double, double*) ;    
void som_r_chebp
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebu
    (double*, const int, const int, const int, const double, double*) ; 
void som_r_chebpim_p
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpim_i
    (double*, const int, const int, const int, const double, double*) ;
void som_r_cheb_symy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebu_symy
    (double*, const int, const int, const int, const double, double*) ; 
void som_r_chebpim_p_symy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpim_i_symy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_cheb_asymy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebu_asymy
    (double*, const int, const int, const int, const double, double*) ; 
void som_r_chebpim_p_asymy
    (double*, const int, const int, const int, const double, double*) ;
void som_r_chebpim_i_asymy
    (double*, const int, const int, const int, const double, double*) ;
    
// Les sommations en theta :
void som_tet_pas_prevu
    (double*, const int, const int, const double, double*) ;
void som_tet_cos
    (double*, const int, const int, const double, double* ) ;
void som_tet_cos_p
    (double*, const int, const int, const double, double* ) ;
void som_tet_cos_i
    (double*, const int, const int, const double, double* ) ;
void som_tet_sin
    (double*, const int, const int, const double, double* ) ;
void som_tet_sin_p
    (double*, const int, const int, const double, double* ) ;
void som_tet_sin_i
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_cp
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_ci
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_sp
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_si
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_cp_symy
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_ci_symy
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_cp_asymy
    (double*, const int, const int, const double, double* ) ;
void som_tet_cossin_ci_asymy
    (double*, const int, const int, const double, double* ) ;

// Les sommations en phi : 
void som_phi_pas_prevu
    (double*, const int, const double, double* ) ;
void som_phi_cossin
    (double*, const int, const double, double* ) ;
void som_phi_cossin_p
    (double*, const int, const double, double* ) ;
void som_phi_cossin_i
    (double*, const int, const double, double* ) ;
void som_phi_cossin_symy
    (double*, const int, const double, double* ) ;
void som_phi_cossin_asymy
    (double*, const int, const double, double* ) ;

// les divisions et multiplications par x-1 :
void sxm1_1d_cheb(int, double*) ;
void mult_xm1_1d_cheb(int, const double*, double*) ;
void mult2_xm1_1d_cheb(int, const double*, double*) ;

// x * dsdx ...
void _xdsdx_1d_pas_prevu (int, double*, double*) ;
void _xdsdx_1d_r_cheb (int, double*, double*) ;
void _xdsdx_1d_r_chebp (int, double*, double*) ;
void _xdsdx_1d_r_chebi (int, double*, double*) ;
void xdsdx_1d(int, double**, int) ;

// Multiplication par x^2
void multx2_1d(int, double **, int)	;
void _multx2_1d_r_cheb(int, double* , double *);
void _multx2_1d_r_chebp(int, double* , double *);
void _multx2_1d_r_chebi(int, double* , double *);
void _multx2_1d_pas_prevu(int, double* , double *);

// division pas (x+1)
void sxpun_1d(int, double **, int)	;
void _sxpun_1d_r_cheb(int, double* , double *);
void _sxpun_1d_pas_prevu(int, double* , double *);
Cmp division_xpun (const Cmp&, int) ;

// Fonctions utilisees par le laplacien :
void _lapang_pas_prevu(Mtbl_cf*, int) ;
void _lapang_t_leg_p(Mtbl_cf*, int) ;
void _lapang_t_leg_pp(Mtbl_cf*, int) ;
void _lapang_t_leg_i(Mtbl_cf*, int) ;

// Fonctions liees au dalembertien
void get_operateur_dal(const Param&, const int&, const int&, 
		       int&, Matrice& );
Tbl dal_inverse(const int&, const int&, const Matrice&, const Tbl&, 
		const bool) ;
Mtbl_cf sol_dalembert(Param&, const Map_af&, const Mtbl_cf&) ;

// Fonctions diverses : 
void c_est_pas_fait(char * ) ;

// Trucs utilises pour poisson_compact :
Matrice lap_cpt_mat(int, int, int) ;
Matrice xdsdx_mat(int, int, int) ;
Matrice combinaison_cpt (const Matrice &, int, int) ;
Tbl combinaison_cpt (const Tbl &, int) ;

// Trucs binaires :
void dirichlet_binaire (const Cmp& source_un, const Cmp& source_deux, 
			const Valeur& boundary_un, const Valeur& boundary_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) ;
void dirichlet_binaire (const Cmp& source_un, const Cmp& source_deux, 
			double bound_un, double bound_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) ;
void neumann_binaire (const Cmp& source_un, const Cmp& source_deux, 
			const Valeur& boundary_un, const Valeur& boundary_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) ;
void neumann_binaire (const Cmp& source_un, const Cmp& source_deux, 
			double bound_un, double bound_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) ;
void poisson_vect_frontiere (double lambda, const Tenseur& source, Tenseur& shift, 
	    const Valeur& lim_x, const Valeur& lim_y, const Valeur& lim_z, 
	    int num_front, double precision, int itermax) ;
void poisson_vect_binaire ( double lambda, 
		const Tenseur& source_un, const Tenseur& source_deux, 
		const Valeur& bound_x_un, const Valeur& bound_y_un, 
		const Valeur& bound_z_un, const Valeur& bound_x_deux, 
		const Valeur& bound_y_deux, const Valeur& bound_z_deux, 
		Tenseur& sol_un, Tenseur& sol_deux, int num_front, double precision) ;

// Regularisation du shift :
double regle (Tenseur& shift_auto, const Tenseur& shift_comp, double omega) ;

// Trucs pour la solution de Misner-Lindquist
double serie_lindquist_plus (double rayon, double distance, double xa, double ya, 
	double za, double precision, double itemax) ;

double serie_lindquist_moins (double rayon, double distance, double xa, double ya, 
	double za, double precision, double itemax) ;

double adm_serie (double rayon, double distance, double precision) ;

double bare_serie (double rayon, double distance, double precision) ;

void set_lindquist (Cmp& psi_un, Cmp& psi_deux, double rayon, double precision) ;

void separation (const Cmp& c1, const Cmp& c2, Cmp& res1, Cmp& res2, int decrois, 
    int puiss, int lmax, double precision, const double relax = 0.5, const int itemax = 100, const int flag = 1) ;


#endif
