/*
 *  Lorene's macros
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2003 Jerome Novak
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


#ifndef	__TYPE_PARITE_H_
#define	__TYPE_PARITE_H_

/*
 * Constantes utilisees dans les types de grilles et les parites
 */

/*
 * $Id$
 * $Log$
 * Revision 1.6  2004/11/04 15:40:14  e_gourgoulhon
 * Added definition of symbol T_LEG.
 *
 * Revision 1.5  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.4  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.3  2003/09/16 08:53:05  j_novak
 * Addition of the T_LEG_II base (odd in theta, only for odd m) and the
 * transformation functions from and to T_SIN_P.
 *
 * Revision 1.2  2002/03/07 15:41:12  n_chamel
 * New class for dealing with Cartesian grids
 * Added the sampling type UNIFORM in type_parite.h
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/12/04  13:09:55  novak
 * Added constants for the dalembertian
 *
 * Revision 2.2  2000/09/28 08:55:52  eric
 * Ajout de T_LEG_IP et T_LEG_PI.
 *
 * Revision 2.1  1999/10/01  15:22:13  eric
 * Vire les jacobiennes.
 *
 * Revision 2.0  1999/02/15  10:41:51  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */
/**
 * \defgroup def_mac LORENE's macros.
 *
 * @{
 */

/* Etat des tbl et autres */
/* ---------------------- */
#define	    ETATZERO        0
#define	    ETATUN	    1
#define	    ETATQCQ	    2
#define	    ETATNONDEF	    3

/* Uniform sampling on a Cartesian grid */
#define         UNIFORM         3

/* Echantillonage fin en r */
/* ----------------------- */
#define	    RARE    1
#define	    FIN	    0
#define	    UNSURR  2

/* Symetries en theta */
/* ------------------ */
#define	    SYM	    1
#define	    NONSYM  0

/* Les bases de developement */
/* ------------------------- */
/// Nombre max. de bases differentes 
#define	    MAX_BASE	32		

    /* Divers (masques, nulls,...) */
/// base inconnue 
#define	    NONDEF	0x00000000	
/// Extraction de l'info sur R 
#define	    MSQ_R	0x000000ff	
/// Extraction de l'info sur Theta 
#define	    MSQ_T	0x0000ff00	
/// Extraction de l'info sur Phi 
#define	    MSQ_P	0x00ff0000	
/// Translation en R 
#define	    TRA_R	0		
/// Translation en Theta 
#define	    TRA_T	8		
/// Translation en Phi 
#define	    TRA_P	16		

    /* R */
/// base de Chebychev ordinaire (fin) 
#define	    R_CHEB	0x00000001	
/// base de Cheb. paire (rare) seulement 
#define	    R_CHEBP	0x00000002	
/// base de Cheb. impaire (rare) seulement 
#define	    R_CHEBI	0x00000003	
/// Cheb. pair-impair suivant l pair pour l=0 
#define	    R_CHEBPI_P	0x00000004	
/// Cheb. pair-impair suivant l impair pour l=0 
#define	    R_CHEBPI_I	0x00000005	
/// Cheb. pair-impair suivant m, pair pour m=0 
#define	    R_CHEBPIM_P 0x00000006	
/// Cheb. pair-impair suivant m, impair pour m=0 
#define	    R_CHEBPIM_I 0x00000007	
/// base de Chebychev ordinaire (fin), dev. en 1/r 
#define	    R_CHEBU	0x00000008	
/// \f$r^l\f$.Cheb. pair, l et m pair 
#define	    R_RLCHEB_PP	0x00000009	
/// \f$r^l\f$.Cheb. pair, l pair ou impair suivant m
#define	    R_RLCHEB_P	0x0000000a	

    /* Theta */
/// dev. cos-sin alternes, cos pour m=0 
#define	    T_COSSIN_C	0x00000100	
/// dev. cos-sin alternes, sin pour m=0 
#define	    T_COSSIN_S	0x00000200	
/// dev. cos seulement 
#define	    T_COS	0x00000300	
/// dev. sin seulement 
#define	    T_SIN	0x00000400	
/// dev. cos seulement, harmoniques paires 
#define	    T_COS_P	0x00000500	
/// dev. sin seulement, harmoniques paires 
#define	    T_SIN_P	0x00000600	
/// dev. cos seulement, harmoniques impaires 
#define	    T_COS_I	0x00000700	
/// dev. sin seulement, harmoniques impaires 
#define	    T_SIN_I	0x00000800	
/// cos pair-sin impair alternes, cos pour m=0 
#define	    T_COSSIN_CP	0x00000900	
/// sin pair-cos impair alternes, sin pour m=0 
#define	    T_COSSIN_SP	0x00000a00	
/// cos impair-sin pair alternes, cos pour m=0 
#define	    T_COSSIN_CI	0x00000b00	
/// sin impair-cos pair alternes, sin pour m=0 
#define	    T_COSSIN_SI	0x00000c00	
/// fct. de Legendre associees paires 
#define	    T_LEG_P	0x00000d00	
/// fct. de Legendre associees paires avec m pair 
#define	    T_LEG_PP	0x00000e00	
/// fct. de Legendre associees impaires 
#define	    T_LEG_I	0x00000f00	
/// fct. de Legendre associees impaires avec m pair 
#define	    T_LEG_IP	0x00001000	
/// fct. de Legendre associees paires avec m impair 
#define	    T_LEG_PI	0x00001100	
/// fct. de Legendre associees impaires avec m impair 
#define	    T_LEG_II	0x00001200	
/// CL of even cosines
#define T_CL_COS_P 0x00001300
/// CL of even sines
#define T_CL_SIN_P 0x00001400
/// CL of odd cosines
#define T_CL_COS_I 0x00001500
/// CL of odd sines.
#define T_CL_SIN_I 0x00001600
/// fct. de Legendre associees 
#define T_LEG 0x00001700


    /* Phi */
/// dev. standart 
#define	    P_COSSIN	0x00010000	
/// dev. sur Phi = 2*phi, freq. paires 
#define	    P_COSSIN_P  0x00020000	
/// dev. sur Phi = 2*phi, freq. impaires 
#define	    P_COSSIN_I  0x00030000	
/// dev. cos seulement 
#define	    P_COS	0x00040000	
/// dev. sin seulement 
#define	    P_SIN	0x00050000	


/************/
/* Type EOS */
/************/
/// eos polytropique 
#define	    POLYTROPE	    0x000000001	
/// eos incompressible 
#define	    INCOMP	    0x000000002	
/// eos polytropique (cas newtonien) 
#define	    POLYTROPE_NEWT  0x000000003	
/// eos incompressible (cas newtonien) 
#define	    INCOMP_NEWT	    0x000000004	

/*******************************/
/* Type operateur dalembertien */
/* (uniquement pour le noyau   */
/*       pour l'instant)       */
/*******************************/
/// Nombre max d'operateurs (pour l'instant)
#define     MAX_DAL       6           
/// Operateur du premier ordre, \f$\delta < \delta_{crit}\f$
#define     ORDRE1_SMALL  0x000000001 
/// Operateur du premier ordre \f$\delta > \delta_{crit}\f$
#define     ORDRE1_LARGE  0x000000002 
/// Operateur du deuxieme ordre degenere \f$\delta < \delta_{crit}\f$
#define     O2DEGE_SMALL  0x000000003 
/// Operateur du deuxieme ordre degenere \f$\delta > \delta_{crit}\f$
#define     O2DEGE_LARGE  0x000000004 
/// Operateur du deuxieme ordre non degenere \f$\delta < \delta_{crit}\f$
#define     O2NOND_SMALL  0x000000005 
/// Operateur du deuxieme ordre non degenere \f$\delta > \delta_{crit}\f$
#define     O2NOND_LARGE  0x000000006 
/**@} */
#endif
