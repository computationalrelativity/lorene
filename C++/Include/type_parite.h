/*
 *  Lorene's macros
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Jerome Novak
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
#define	    MAX_BASE	32		/* Nombre max. de bases differentes */

    /* Divers (masques, nulls,...) */
#define	    NONDEF	0x00000000	/* base inconnue */
#define	    MSQ_R	0x000000ff	/* Extraction de l'info sur R */
#define	    MSQ_T	0x0000ff00	/* Extraction de l'info sur Theta */
#define	    MSQ_P	0x00ff0000	/* Extraction de l'info sur Phi */
#define	    TRA_R	0		/* Translation en R */
#define	    TRA_T	8		/* Translation en Theta */
#define	    TRA_P	16		/* Translation en Phi */

    /* R */
#define	    R_CHEB	0x00000001	/* base de Chebychev ordinaire (fin) */
#define	    R_CHEBP	0x00000002	/* base de Cheb. paire (rare) seulement */
#define	    R_CHEBI	0x00000003	/* base de Cheb. impaire (rare) seulement */
#define	    R_CHEBPI_P	0x00000004	/* Cheb. pair-impair suivant l pair pour l=0 */
#define	    R_CHEBPI_I	0x00000005	/* Cheb. pair-impair suivant l impair pour l=0 */
#define	    R_CHEBPIM_P 0x00000006	/* Cheb. pair-impair suivant m, pair pour m=0 */
#define	    R_CHEBPIM_I 0x00000007	/* Cheb. pair-impair suivant m, impair pour m=0 */
#define	    R_CHEBU	0x00000008	/* base de Chebychev ordinaire (fin), dev. en 1/r */
#define	    R_RLCHEB_PP	0x00000009	/* r^l.Cheb. pair, l et m pair */
#define	    R_RLCHEB_P	0x0000000a	/* r^l.Cheb. pair, l pair ou impair suivant m*/

    /* Theta */
#define	    T_COSSIN_C	0x00000100	/* dev. cos-sin alternes, cos pour m=0 */
#define	    T_COSSIN_S	0x00000200	/* dev. cos-sin alternes, sin pour m=0 */
#define	    T_COS	0x00000300	/* dev. cos seulement */
#define	    T_SIN	0x00000400	/* dev. sin seulement */
#define	    T_COS_P	0x00000500	/* dev. cos seulement, harmoniques paires */
#define	    T_SIN_P	0x00000600	/* dev. sin seulement, harmoniques paires */
#define	    T_COS_I	0x00000700	/* dev. cos seulement, harmoniques impaires */
#define	    T_SIN_I	0x00000800	/* dev. sin seulement, harmoniques impaires */
#define	    T_COSSIN_CP	0x00000900	/* cos pair-sin impair alternes, cos pour m=0 */
#define	    T_COSSIN_SP	0x00000a00	/* sin pair-cos impair alternes, sin pour m=0 */
#define	    T_COSSIN_CI	0x00000b00	/* cos impair-sin pair alternes, cos pour m=0 */
#define	    T_COSSIN_SI	0x00000c00	/* sin impair-cos pair alternes, sin pour m=0 */
#define	    T_LEG_P	0x00000d00	/* fct. de Legendre associees paires */
#define	    T_LEG_PP	0x00000e00	/* fct. de Legendre associees paires avec m pair */
#define	    T_LEG_I	0x00000f00	/* fct. de Legendre associees impaires */
#define	    T_LEG_IP	0x00001000	/* fct. de Legendre associees impaires avec m pair */
#define	    T_LEG_PI	0x00001100	/* fct. de Legendre associees paires avec m impair */

    /* Phi */
#define	    P_COSSIN	0x00010000	/* dev. standart */
#define	    P_COSSIN_P  0x00020000	/* dev. sur Phi = 2*phi, freq. paires */
#define	    P_COSSIN_I  0x00030000	/* dev. sur Phi = 2*phi, freq. impaires */
#define	    P_COS	0x00040000	/* dev. cos seulement */
#define	    P_SIN	0x00050000	/* dev. sin seulement */


/************/
/* Type EOS */
/************/
#define	    POLYTROPE	    0x000000001	/* eos polytropique */
#define	    INCOMP	    0x000000002	/* eos incompressible */
#define	    POLYTROPE_NEWT  0x000000003	/* eos polytropique (cas newtonien) */
#define	    INCOMP_NEWT	    0x000000004	/* eos incompressible (cas newtonien) */

/*******************************/
/* Type operateur dalembertien */
/* (uniquement pour le noyau   */
/*       pour l'instant)       */
/*******************************/
#define     MAX_DAL       6           /* Nombre max d'operateurs (pour l'instant)*/
#define     ORDRE1_SMALL  0x000000001 /* Operateur du premier ordre, $\delta < \delata_{crit}$*/
#define     ORDRE1_LARGE  0x000000002 /* Operateur du premier ordre $\delta > \delata_{crit}$*/
#define     O2DEGE_SMALL  0x000000003 /* Operateur du deuxieme ordre degenere $\delta < \delta_{crit}$*/
#define     O2DEGE_LARGE  0x000000004 /* Operateur du deuxieme ordre degenere $\delta > \delta_{crit}$*/
#define     O2NOND_SMALL  0x000000005 /* Operateur du deuxieme ordre non degenere $\delta < \delta_{crit}$*/
#define     O2NOND_LARGE  0x000000006 /* Operateur du deuxieme ordre non degenere $\delta > \delta_{crit}$*/

#endif
