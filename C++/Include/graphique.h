/*
 *  Prototypes of graphical routines
 *
 */

/*
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


#ifndef	__GRAPHIQUE_H_
#define	__GRAPHIQUE_H_

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/06/03 10:00:37  e_gourgoulhon
 * Added a new version of des_profile for Cmp with scale and nomx
 * specified in the argument list
 *
 * Revision 1.3  2003/01/17 13:48:17  f_limousin
 * Add des_explorer and des_explorer_symz for a Bin_ns_ncp
 *
 * Revision 1.2  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.24  2001/06/21  07:35:44  novak
 * Added two routines for 2-surface star drawing (des_bi_coupe_y)
 *
 * Revision 1.23  2001/05/22 13:31:54  eric
 * Ajout de des_explorer_coef
 *
 * Revision 1.22  2001/03/07  10:47:09  eric
 * Ajout de des_explorer_symz
 *
 * Revision 1.21  2000/12/04  14:16:55  novak
 * des_explorer2D added
 *
 * Revision 1.20  2000/06/22 16:09:03  eric
 * Retour a la version 1.18 (1.19 etait une erreur).
 *
 * Revision 1.18  2000/03/02  10:33:32  eric
 * Ajout des routines des_vect_bin_*
 *
 * Revision 1.17  2000/03/01  16:11:14  eric
 * Ajout des dessins de champs vectoriels.
 *
 * Revision 1.16  2000/02/12  11:17:46  eric
 * Ajout des versions de des_coupe_* avec determination automatique des
 * bornes de la fenetre graphique.
 *
 * Revision 1.15  2000/02/11  18:43:27  eric
 * Ajout de l'argument draw_bound aux routines des_coupe*.
 *
 * Revision 1.14  2000/02/11  17:47:33  eric
 * Ajout des routines des_coupe_bin_*
 *
 * Revision 1.13  2000/02/11  16:51:49  eric
 * Les routines de dessins de Cmp utilisent desormais les coordonnees
 * cartesiennes abolues (X,Y,Z) et non plus relatives (x,y,z).
 *
 * Revision 1.12  2000/02/11  09:58:12  eric
 * *** empty log message ***
 *
 * Revision 1.11  2000/02/11  09:56:14  eric
 * Ajout des sorties pour Explorer.
 *
 * Revision 1.10  1999/12/27  12:22:25  eric
 * *** empty log message ***
 *
 * Revision 1.9  1999/12/27  12:17:11  eric
 * Ajout des routines des_domaine_*.
 * Les valeurs par defaut du nombre de mailles pour le quadrillage des
 * dans des_coupe_* passent de 80x80 a 100x100.
 *
 * Revision 1.8  1999/12/24  12:59:38  eric
 * Ajout des routines des_surface_*
 *
 * Revision 1.7  1999/12/23  16:14:33  eric
 * Les routines des_coupe_* dessine desormais egalement la surface
 *  de l'objet (ajout de l'argument defsurf).
 *
 * Revision 1.6  1999/12/20  11:04:24  eric
 * Modif commentaires.
 *
 * Revision 1.5  1999/12/20  11:00:38  eric
 * *** empty log message ***
 *
 * Revision 1.4  1999/12/20  10:53:52  eric
 * Ajout des arguments device, newgraph, nxpage et nypage
 *  a des_coef_xi, des_coef_theta et des_coef_phi.
 * Ajout de la routine des_map_et.
 *
 * Revision 1.3  1999/12/15  09:42:02  eric
 * *** empty log message ***
 *
 * Revision 1.2  1999/12/10  12:14:09  eric
 * Ajout des fonctions des_coef.
 *
 * Revision 1.1  1999/12/09  16:37:57  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */
 
class Valeur ; 
class Map ; 
class Map_et ; 
class Cmp ; 
class Tenseur ; 
class Etoile ; 
class Binaire ; 
class Bin_ns_ncp ;

    /** @name Basic routines.
     */
    //@{     


/** Basic routine for drawing profiles.
 *  A profile is a function y=y(x). 
 *
 *  @param uutab [input] Array (size: {\tt nx}) of y values to be drawn
 *			 (the x sampling is supposed to be uniform).
 *  @param nx [input] Number of points
 *  @param xmin [input] lowest value of x
 *  @param xmax [input] highest value of x
 *  @param nomx [input] x legend of the figure
 *  @param nomy [input] y legend of the figure
 *  @param title [input] title of the figure
 *  @param device [input] PGPLOT device (default value = 0x0)
 *
 */
void des_profile(float* uutab, int nx, float xmin, float xmax, 
		 char* nomx, char* nomy, char* title, char* device = 0x0) ;


/** Basic routine for drawing isocontours.
 * 
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. 
 *
 *  @param uutab [input] field to be drawn;  
 *	    the value of the field a the point of coordinates \\
 *		      x\_i = xmin + i (xmax-xmin)/(nx-1)     0 <= i <= nx-1 \\  
 *		      y\_j = ymin + j (ymax-ymin)/(ny-1)     0 <= j <= ny-1 \\
 *     must be stored at the following position in the float 1-D array uu : \\
 *			index = j * nx + i 
 *  @param nx [input]  number of points in the x direction
 *  @param ny [input]  number of points in the y direction
 *  @param xmin [input] lowest value of x 
 *  @param xmax [input] highest value of x 
 *  @param ymin [input] lowest value of y 
 *  @param ymax [input] highest value of y
 *  @param ncour [input] number of isocontour lines
 *  @param nomx [input] x legend of the figure
 *  @param nomy [input] y legend of the figure
 *  @param title [input] title of the figure
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *			    
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_equipot(float* uutab, int nx, int ny, float xmin, float xmax, 
		 float ymin, float ymax, int ncour, char* nomx, char* nomy, 
		 char* title, char* device = 0x0, int newgraph = 3, 
		 int nxpage = 1, int nypage = 1) ;

/** Basic routine for plotting vector field. 
 * 
 *
 *  @param vvx [input] x-component of the vector field to be drawn ;  
 *	    the value of the field a the point of coordinates \\
 *		      x\_i = xmin + i (xmax-xmin)/(nx-1)     0 <= i <= nx-1 \\  
 *		      y\_j = ymin + j (ymax-ymin)/(ny-1)     0 <= j <= ny-1 \\
 *     must be stored at the following position in the float 1-D array uu : \\
 *			index = j * nx + i 
 *  @param vvy [input] y-component of the vector field to be drawn ;  
 *		       same storage as {\tt vvx}.
 *  @param nx [input]  number of points in the x direction
 *  @param ny [input]  number of points in the y direction
 *  @param xmin [input] lowest value of x 
 *  @param xmax [input] highest value of x 
 *  @param ymin [input] lowest value of y 
 *  @param ymax [input] highest value of y
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param nomx [input] x legend of the figure
 *  @param nomy [input] y legend of the figure
 *  @param title [input] title of the figure
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *			    
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_vect(float* vvx, float* vvy, int nx, int ny, float xmin, float xmax, 
		 float ymin, float ymax, double scale,  double sizefl, 
		 char* nomx, char* nomy, char* title, char* device = 0x0, 
		 int newgraph = 3, int nxpage = 1, int nypage = 1) ;


/** Basic routine for drawing spectral coefficients.
 *  
 *  @param cf  [input] 1-D array of the coefficients to be drawn (size: {\tt n})
 *  @param n  [input] number of coefficients
 *  @param pzero [input] positive number under which (in absolute value)
 *		        a coefficient will be considered as zero 
 *  @param nomx [input] x legend of the figure 
 *  @param nomy [input] y legend of the figure 
 *  @param title [input] title of the figure
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *			    
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_coef(const double* cf, int n, double pzero,
	      char* nomx, char* nomy, char* title, char* device = 0x0, 
	      int newgraph = 3, int nxpage = 1, int nypage = 1) ;

    //@}
    
    /** @name Plots of spectral coefficients
     */
    //@{     

/** Plots the coefficients of the spectral expansion in $\xi$ of a {\tt Valeur}.
 * 
 *  This routine performs a logarithmic plot of the coefficients of the
 *  $\xi$ expansion of the coefficient which is in front of a given
 *  $\phi'$ basis function (index {\tt k}) and a given $\theta'$
 *  basis function (index {\tt j}) in the general spectral expansion of a
 *  field. The plotted quantities are the logarithm of the absolute value
 *  of the coefficients, using solid lines (resp. dashed lines) for 
 *  positive coefficients (resp. negative coefficients). 
 * 
 *  @param uu [input] {\tt Valeur} the $\xi$ coefficients of which are to 
 *		      be plotted
 *  @param l [input] index of the domain 
 *  @param k [input] index of the considered basis function in $\phi'$
 *  @param j [input] index of the considered basis function in $\theta'$
 *  @param pzero [input] positive number under which (in absolute value)
 *		        a coefficient will be considered as zero 
 *		        (default value = 1.e-14)
 *  @param nomy [input] y legend of the figure (default value = 0x0,  
 *		        corresponds to no y legend)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *			    
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_coef_xi(const Valeur& uu, int l, int k, int j, double pzero = 1.e-14, 
		 char* nomy = 0x0, char* title = 0x0, char* device = 0x0, 
	         int newgraph = 3, int nxpage = 1, int nypage = 1) ;


/** Plots the coefficients of the spectral expansion in $\theta'$ of a 
 * {\tt Valeur}.
 * 
 *  This routine performs a logarithmic plot of the coefficients of the
 *  $\theta'$ expansion of the coefficient which is in front of a given
 *  $\phi'$ basis function (index {\tt k}) and a given $\xi$
 *  basis function (index {\tt i}) in the general spectral expansion of a
 *  field. The plotted quantities are the logarithm of the absolute value
 *  of the coefficients, using solid lines (resp. dashed lines) for 
 *  positive coefficients (resp. negative coefficients). 
 * 
 *  @param uu [input] {\tt Valeur} the $\xi$ coefficients of which are to 
 *		      be plotted
 *  @param l [input] index of the domain 
 *  @param k [input] index of the considered basis function in $\phi'$
 *  @param i [input] index of the considered basis function in $\xi$
 *  @param pzero [input] positive number under which (in absolute value)
 *		        a coefficient will be considered as zero 
 *		        (default value = 1.e-14)
 *  @param nomy [input] y legend of the figure (default value = 0x0,  
 *		        corresponds to no y legend)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *			    
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_coef_theta(const Valeur& uu, int l, int k, int i, double pzero = 1.e-14, 
		    char* nomy = 0x0, char* title = 0x0, char* device = 0x0, 
	            int newgraph = 3, int nxpage = 1, int nypage = 1) ;
		 
		 
/** Plots the coefficients of the spectral expansion in $\phi'$ of a 
 *  {\tt Valeur}.
 * 
 *  This routine performs a logarithmic plot of the coefficients of the
 *  $\phi'$ expansion of the coefficient which is in front of a given
 *  $\theta'$ basis function (index {\tt j}) and a given $\xi$
 *  basis function (index {\tt i}) in the general spectral expansion of a
 *  field. The plotted quantities are the logarithm of the absolute value
 *  of the coefficients, using solid lines (resp. dashed lines) for 
 *  positive coefficients (resp. negative coefficients). 
 * 
 *  @param uu [input] {\tt Valeur} the $\xi$ coefficients of which are to 
 *		      be plotted
 *  @param l [input] index of the domain 
 *  @param j [input] index of the considered basis function in $\theta'$
 *  @param i [input] index of the considered basis function in $\xi$
 *  @param pzero [input] positive number under which (in absolute value)
 *		        a coefficient will be considered as zero 
 *		        (default value = 1.e-14)
 *  @param nomy [input] y legend of the figure (default value = 0x0,  
 *		        corresponds to no y legend)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *			    
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_coef_phi(const Valeur& uu, int l, int j, int i, double pzero = 1.e-14, 
		  char* nomy = 0x0, char* title = 0x0, char* device = 0x0, 
	          int newgraph = 3, int nxpage = 1, int nypage = 1) ;

/** Plots the coefficients of the functions $F_l(\theta', \phi')$ and
 *  $G_l(\theta', \phi')$ of a mapping of class {\tt Map\_et}. 
 * 
 *  This routine performs a logarithmic plot of the coefficients of the
 *  $\theta'$ (resp. $\phi'$) expansion of the coefficient which is in front 
 *  of a given $\phi'$ (resp. $\theta'$) basis function (index {\tt k}
 *  (resp. index {\tt j})).
 *  The plotted quantities are the logarithm of the absolute value
 *  of the coefficients, using solid lines (resp. dashed lines) for 
 *  positive coefficients (resp. negative coefficients). 
 * 
 *  @param mp [input] Mapping of class {\tt Map\_et}
 *  @param lz [input] Index of the domain where the plot is to performed
 * 
 */
void des_map_et(const Map_et& mp, int lz) ;


    //@}



    /** @name Plot of a scalar field
     */
    //@{     


/** Draws the profile of a {\tt Cmp} along some radial axis determined by
 *  a fixed value of $(\theta, \phi)$. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param r_min [input] Minimal value of {\it r} for the drawing
 *  @param r_max [input] Maximal value of {\it r} for the drawing
 *  @param theta [input] Value of $\theta$ which defines the profile axis
 *  @param phi [input] Value of $\phi$ which defines the profile axis
 *  @param nomy [input] y legend of the figure (default value = 0x0,  
 *		        corresponds to no y legend)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 * 
 */
 
void des_profile(const Cmp& uu, double r_min, double r_max, 
		     double theta, double phi, char* nomy = 0x0,  
		     char* title = 0x0 ) ;


/** Draws the profile of a {\tt Cmp} along some radial axis determined by
 *  a fixed value of $(\theta, \phi)$. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param r_min [input] Minimal value of {\it r} for the drawing
 *  @param r_max [input] Maximal value of {\it r} for the drawing
 *  @param scale scale factor for the radius in the plot
 *  @param theta [input] Value of $\theta$ which defines the profile axis
 *  @param phi [input] Value of $\phi$ which defines the profile axis
 *  @param nomx [input] x legend of the figure (default value = 0x0,  
 *		        corresponds to no x legend)
 *  @param nomy [input] y legend of the figure (default value = 0x0,  
 *		        corresponds to no y legend)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 * 
 */
 
void des_profile(const Cmp& uu, double r_min, double r_max, double scale,
		     double theta, double phi, char* nomx = 0x0, 
		     char* nomy = 0x0, char* title= 0x0) ;


/** Basic routine for drawing a stellar surface in a plane X=constant.
 * 
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *  The surface to be drawn is defined as the location where a given scalar
 *  field vanishes. 
 *  
 *  @param defsurf [input] field defining the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes 
 *  @param x0 [input] value of the absolute coordinate X which defines the plane
 *		      of the drawing
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *  @param y_min [input] lowest value of absol. coord. Y (default value = -1)
 *  @param y_max [input] highest value of absol. coord. Y (default value = 1)
 *  @param z_min [input] lowest value of absol. coord. Z (default value = -1)
 *  @param z_max [input] highest value of absol. coord. Z (default value = 1)
 *  @param nomy [input] y legend of the figure (default value = 0x0)
 *  @param nomz [input] z legend of the figure (default value = 0x0)
 *  @param title [input] title of the figure (default value = 0x0)
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_surface_x(const Cmp& defsurf, double x0, char* device = 0x0, 
		   int newgraph = 3, double y_min = -1, double y_max = 1, 
		   double z_min = -1, double z_max = 1, 
		   char* nomy = 0x0, char* nomz = 0x0, char* title = 0x0, 
		   int nxpage = 1, int nypage = 1) ;

/** Basic routine for drawing a stellar surface in a plane Y=constant.
 * 
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *  The surface to be drawn is defined as the location where a given scalar
 *  field vanishes. 
 *  
 *  @param defsurf [input] field defining the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes 
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *  @param x_min [input] lowest value of absol. coord. X (default value = -1)
 *  @param x_max [input] highest value of absol. coord. X (default value = 1)
 *  @param z_min [input] lowest value of absol. coord. Z (default value = -1)
 *  @param z_max [input] highest value of absol. coord. Z (default value = 1)
 *  @param nomx [input] x legend of the figure (default value = 0x0)
 *  @param nomz [input] z legend of the figure (default value = 0x0)
 *  @param title [input] title of the figure (default value = 0x0)
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_surface_y(const Cmp& defsurf, double y0, char* device = 0x0, 
		   int newgraph = 3, double x_min = -1, double x_max = 1, 
		   double z_min = -1, double z_max = 1, 
		   char* nomx = 0x0, char* nomz = 0x0, char* title = 0x0, 
		   int nxpage = 1, int nypage = 1) ;


/** Basic routine for drawing a stellar surface in a plane Z=constant.
 * 
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *  The surface to be drawn is defined as the location where a given scalar
 *  field vanishes. 
 *  
 *  @param defsurf [input] field defining the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes 
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *  @param x_min [input] lowest value of absol. coord. X (default value = -1)
 *  @param x_max [input] highest value of absol. coord. X (default value = 1)
 *  @param y_min [input] lowest value of absol. coord. Y (default value = -1)
 *  @param y_max [input] highest value of absol. coord. Y (default value = 1)
 *  @param nomx [input] x legend of the figure (default value = 0x0)
 *  @param nomy [input] y legend of the figure (default value = 0x0)
 *  @param title [input] title of the figure (default value = 0x0)
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_surface_z(const Cmp& defsurf, double z0, char* device = 0x0, 
		   int newgraph = 3, double x_min = -1, double x_max = 1, 
		   double y_min = -1, double y_max = 1, 
		   char* nomx = 0x0, char* nomz = 0x0, char* title = 0x0, 
		   int nxpage = 1, int nypage = 1) ;


/** Basic routine for drawing the outer boundary of a given domain 
 *  in a plane X=constant.
 * 
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *  The domain outer boundary is defined by $\xi = 1$. 
 *  
 *  @param mp [input] Mapping defining the various domains
 *  @param l0 [input] Index of the domain, the outer boundary of which is
 *		      to be drawn
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *  @param y_min [input] lowest value of absol. coord. Y (default value = -1)
 *  @param y_max [input] highest value of absol. coord. Y (default value = 1)
 *  @param z_min [input] lowest value of absol. coord. Z (default value = -1)
 *  @param z_max [input] highest value of absol. coord. Z (default value = 1)
 *  @param nomy [input] y legend of the figure (default value = 0x0)
 *  @param nomz [input] z legend of the figure (default value = 0x0)
 *  @param title [input] title of the figure (default value = 0x0)
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_domaine_x(const Map& mp, int l0, double x0, char* device = 0x0, 
		   int newgraph = 3, double y_min = -1, double y_max = 1, 
		   double z_min = -1, double z_max = 1, 
		   char* nomy = 0x0, char* nomz = 0x0, char* title = 0x0, 
		   int nxpage = 1, int nypage = 1) ;


/** Basic routine for drawing the outer boundary of a given domain 
 *  in a plane Y=constant.
 * 
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *  The domain outer boundary is defined by $\xi = 1$. 
 *  
 *  @param mp [input] Mapping defining the various domains
 *  @param l0 [input] Index of the domain, the outer boundary of which is
 *		      to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *  @param x_min [input] lowest value of absol. coord. X (default value = -1)
 *  @param x_max [input] highest value of absol. coord. X (default value = 1)
 *  @param z_min [input] lowest value of absol. coord. Z (default value = -1)
 *  @param z_max [input] highest value of absol. coord. Z (default value = 1)
 *  @param nomx [input] x legend of the figure (default value = 0x0)
 *  @param nomz [input] z legend of the figure (default value = 0x0)
 *  @param title [input] title of the figure (default value = 0x0)
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_domaine_y(const Map& mp, int l0, double y0, char* device = 0x0, 
		   int newgraph = 3, double x_min = -1, double x_max = 1, 
		   double z_min = -1, double z_max = 1, 
		   char* nomx = 0x0, char* nomz = 0x0, char* title = 0x0, 
		   int nxpage = 1, int nypage = 1) ;


/** Basic routine for drawing the outer boundary of a given domain 
 *  in a plane Z=constant.
 * 
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *  The domain outer boundary is defined by $\xi = 1$. 
 *  
 *  @param mp [input] Mapping defining the various domains
 *  @param l0 [input] Index of the domain, the outer boundary of which is
 *		      to be drawn
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param device [input] PGPLOT device (default value = 0x0)
 *  @param newgraph [input] controls the opening/closing of the graphic device: \\
 *			0 : does nothing (the device must be already opened) \\
 *			1 : opens the device but does not close it at the end \\
 *			2 : closes the device at the end but does not open it at
 *			    the beginning \\
 *		        3 (default value) : opens and closes the device 
 *  @param x_min [input] lowest value of absol. coord. X (default value = -1)
 *  @param x_max [input] highest value of absol. coord. X (default value = 1)
 *  @param y_min [input] lowest value of absol. coord. Y (default value = -1)
 *  @param y_max [input] highest value of absol. coord. Y (default value = 1)
 *  @param nomx [input] x legend of the figure (default value = 0x0)
 *  @param nomy [input] y legend of the figure (default value = 0x0)
 *  @param title [input] title of the figure (default value = 0x0)
 *  @param nxpage [input] number of graphs in the horizontal direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 *  @param nypage [input] number of graphs in the vertical direction of the
 *			  display window (meaningfull only if 
 *			  {\tt newgraph} = 1 or 3) (default value = 1)
 */
void des_domaine_z(const Map& mp, int l0, double z0, char* device = 0x0, 
		   int newgraph = 3, double x_min = -1, double x_max = 1, 
		   double y_min = -1, double y_max = 1, 
		   char* nomx = 0x0, char* nomz = 0x0, char* title = 0x0, 
		   int nxpage = 1, int nypage = 1) ;



/** Draws isocontour lines of a {\tt Cmp} in a plane X=constant.
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the {\tt nzdes} innermost domains fit in it (for
 *			{\tt zoom} = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the {\tt nzdes} innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_x(const Cmp& uu, double x0, int nzdes, char* title = 0x0, 
		 const Cmp* defsurf = 0x0, double zoom = 1.2, 
		 bool draw_bound = true, int ncour = 15, int ny = 100, 
		 int nz = 100) ; 


/** Draws isocontour lines of a {\tt Cmp} in a plane X=constant
 *  within a specified graphic window. 
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param y_min [input] lowest value of absol. coord. Y 
 *  @param y_max [input] highest value of absol. coord. Y 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_x(const Cmp& uu, double x0, double y_min, double y_max, 
		 double z_min, double z_max, char* title = 0x0, 
		 const Cmp* defsurf = 0x0, bool draw_bound = true,
		 int ncour = 15, int ny = 100, int nz = 100) ; 


/** Draws isocontour lines of a {\tt Cmp} in a plane Y=constant.
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the {\tt nzdes} innermost domains fit in it (for
 *			{\tt zoom} = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the {\tt nzdes} innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_y(const Cmp& uu, double y0, int nzdes, char* title = 0x0, 
		 const Cmp* defsurf = 0x0, double zoom = 1.2,
		 bool draw_bound = true, int ncour = 15, int nx = 100, 
		 int nz = 100) ; 

/** Draws isocontour lines of a {\tt Cmp} in a plane Y=constant
 *  within a specified graphic window. 
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param z_min [input] lowest value of absol. coord. Z 
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_y(const Cmp& uu, double y0, double x_min, double x_max, 
		 double z_min, double z_max, char* title = 0x0, 
		 const Cmp* defsurf = 0x0, bool draw_bound = true,
		 int ncour = 15, int nx = 100, int nz = 100) ; 


/** Draws isocontour lines of a {\tt Cmp} in a plane Z=constant.
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the {\tt nzdes} innermost domains fit in it (for
 *			{\tt zoom} = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the {\tt nzdes} innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 */
void des_coupe_z(const Cmp& uu, double z0, int nzdes, char* title = 0x0, 
		 const Cmp* defsurf = 0x0, double zoom = 1.2, 
		 bool draw_bound = true, int ncour = 15, int nx = 100, 
		 int ny = 100) ;

/** Draws isocontour lines of a {\tt Cmp} in a plane Z=constant
 *  within a specified graphic window. 
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. 
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param y_min [input] lowest value of absol. coord. Y
 *  @param y_max [input] highest value of absol. coord. Y
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 */
void des_coupe_z(const Cmp& uu, double z0, double x_min, double x_max, 
		 double y_min, double y_max, char* title = 0x0, 
		 const Cmp* defsurf = 0x0, bool draw_bound = true,
		 int ncour = 15, int nx = 100, int ny = 100) ;


/** Draws isocontour lines of a {\tt Cmp} in a plane Y=constant.
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. The routine allows for the drawing of two
 *  surfaces given by two enthalpies (defsurf and defsurf2).
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the {\tt nzdes} innermost domains fit in it (for
 *			{\tt zoom} = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface for the first fluid (see et_rot_biluid.h): 
 *                         the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of 
 *			   the surface for the second fluid (analog to defsurf)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the {\tt nzdes} innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_bi_coupe_y(const Cmp& uu, double y0, int nzdes, char* title = 0x0, 
                 const Cmp* defsurf = 0x0, const Cmp* defsurf2 = 0x0, 
		 double zoom = 1.2,
                 bool draw_bound = true, int ncour = 15, int nx = 100, 
                 int nz = 100) ; 

/** Draws isocontour lines of a {\tt Cmp} in a plane Y=constant
 *  within a specified graphic window. 
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field. The routine allows for the drawing of two
 *  surfaces given by two enthalpies (defsurf and defsurf2).
 *
 *  @param uu [input] {\tt Cmp} to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param z_min [input] lowest value of absol. coord. Z 
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface for the first fluid (see et_rot_biluid.h): 
 *                         the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of  
 *			   the surface for the second fluid (analog to defsurf)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_bi_coupe_y(const Cmp& uu, double y0, double x_min, double x_max, 
		 double z_min, double z_max, char* title = 0x0, 
		 const Cmp* defsurf = 0x0, const Cmp* defsurf2 = 0x0, 
		 bool draw_bound = true,
		 int ncour = 15, int nx = 100, int nz = 100) ; 


/** Draws isocontour lines of a the sum of two {\tt Cmp}'s in a plane X=constant.
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field {\tt uu1 + uu2}. 
 *
 *  @param uu1 [input] first {\tt Cmp} to define the field {\tt uu1 + uu2} to
 *		       be drawn.  
 *  @param uu2 [input] second {\tt Cmp} to define the field {\tt uu1 + uu2} to
 *		       be drawn.  
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param y_min [input] lowest value of absol. coord. Y 
 *  @param y_max [input] highest value of absol. coord. Y 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf1 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   first surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   second surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_bin_x(const Cmp& uu1, const Cmp& uu2, double x0, double y_min, 
		     double y_max, double z_min, double z_max, char* title, 
		     const Cmp* defsurf1 = 0x0,  const Cmp* defsurf2 = 0x0, 
		     bool draw_bound = true, int ncour = 15, int ny = 100, 
		     int nz = 100) ; 


/** Draws isocontour lines of a the sum of two {\tt Cmp}'s in a plane Y=constant.
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field {\tt uu1 + uu2}. 
 *
 *  @param uu1 [input] first {\tt Cmp} to define the field {\tt uu1 + uu2} to
 *		       be drawn.  
 *  @param uu2 [input] second {\tt Cmp} to define the field {\tt uu1 + uu2} to
 *		       be drawn.  
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf1 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   first surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   second surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_bin_y(const Cmp& uu1, const Cmp& uu2, double y0, double x_min, 
		     double x_max, double z_min, double z_max, char* title, 
		     const Cmp* defsurf1 = 0x0,  const Cmp* defsurf2 = 0x0, 
		     bool draw_bound = true, int ncour = 15, int nx = 100, 
		     int nz = 100) ; 


/** Draws isocontour lines of a the sum of two {\tt Cmp}'s in a plane Z=constant.
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *  Solid (resp. dashed) lines correspond to positive (resp. negative)
 *  values of the field {\tt uu1 + uu2}. 
 *
 *  @param uu1 [input] first {\tt Cmp} to define the field {\tt uu1 + uu2} to
 *		       be drawn.  
 *  @param uu2 [input] second {\tt Cmp} to define the field {\tt uu1 + uu2} to
 *		       be drawn.  
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param y_min [input] lowest value of absol. coord. Y
 *  @param y_max [input] highest value of absol. coord. Y
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf1 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   first surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   second surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ncour [input] number of isocontour lines (default value = 15)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 */
void des_coupe_bin_z(const Cmp& uu1, const Cmp& uu2, double z0, double x_min, 
		     double x_max, double y_min, double y_max, char* title, 
		     const Cmp* defsurf1 = 0x0,  const Cmp* defsurf2 = 0x0, 
		     bool draw_bound = true, int ncour = 15, int nx = 100, 
		     int ny = 100) ; 


    //@}


    /** @name Plot of a vector field
     */
    //@{     

/** Plots a vector field in a plane X=constant.
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv [input] vector field to be drawn
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the {\tt nzdes} innermost domains fit in it (for
 *			{\tt zoom} = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the {\tt nzdes} innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ny [input]  number of points in the Y direction (default value = 20)
 *  @param nz [input]  number of points in the Z direction (default value = 20)
 */
void des_coupe_vect_x(const Tenseur& vv, double x0, double scale, double sizefl,
		      int nzdes, char* title = 0x0, const Cmp* defsurf = 0x0, 
		      double zoom = 1.2, bool draw_bound = true, 
		      int ny = 20, int nz = 20) ; 


/** Plots a vector field in a plane X=constant
 *  within a specified graphic window. 
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv [input] vector field to be drawn
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param y_min [input] lowest value of absol. coord. Y 
 *  @param y_max [input] highest value of absol. coord. Y 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ny [input]  number of points in the Y direction (default value = 20)
 *  @param nz [input]  number of points in the Z direction (default value = 20)
 */
void des_coupe_vect_x(const Tenseur& vv, double x0, double scale, double
		      sizefl, double y_min, double y_max, double z_min, 
		      double z_max, char* title = 0x0, const Cmp* defsurf = 0x0,
		      bool draw_bound = true, int ny = 20, int nz = 20) ;

/** Plots a vector field in a plane Y=constant.
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv [input] vector field to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the {\tt nzdes} innermost domains fit in it (for
 *			{\tt zoom} = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the {\tt nzdes} innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param nx [input]  number of points in the X direction (default value = 20)
 *  @param nz [input]  number of points in the Z direction (default value = 20)
 */
void des_coupe_vect_y(const Tenseur& vv, double y0, double scale, double sizefl,
		      int nzdes, char* title = 0x0, const Cmp* defsurf = 0x0, 
		      double zoom = 1.2, bool draw_bound = true, 
		      int nx = 20, int nz = 20) ; 


/** Plots a vector field in a plane Y=constant
 *  within a specified graphic window. 
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv [input] vector field to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param nx [input]  number of points in the X direction (default value = 20)
 *  @param nz [input]  number of points in the Z direction (default value = 20)
 */
void des_coupe_vect_y(const Tenseur& vv, double y0, double scale, double
		      sizefl, double x_min, double x_max, double z_min, 
		      double z_max, char* title = 0x0, const Cmp* defsurf = 0x0,
		      bool draw_bound = true, int nx = 20, int nz = 20) ;

/** Plots a vector field in a plane Z=constant.
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv [input] vector field to be drawn
 *  @param z0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the {\tt nzdes} innermost domains fit in it (for
 *			{\tt zoom} = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the {\tt nzdes} innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param nx [input]  number of points in the X direction (default value = 20)
 *  @param ny [input]  number of points in the Y direction (default value = 20)
 */
void des_coupe_vect_z(const Tenseur& vv, double z0, double scale, double sizefl,
		      int nzdes, char* title = 0x0, const Cmp* defsurf = 0x0, 
		      double zoom = 1.2, bool draw_bound = true, 
		      int nx = 20, int ny = 20) ; 


/** Plots a vector field in a plane Z=constant
 *  within a specified graphic window. 
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv [input] vector field to be drawn
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param y_min [input] lowest value of absol. coord. Y
 *  @param y_max [input] highest value of absol. coord. Y
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param nx [input]  number of points in the X direction (default value = 20)
 *  @param ny [input]  number of points in the Y direction (default value = 20)
 */
void des_coupe_vect_z(const Tenseur& vv, double z0, double scale, double
		      sizefl, double x_min, double x_max, double y_min, 
		      double y_max, char* title = 0x0, const Cmp* defsurf = 0x0,
		      bool draw_bound = true, int nx = 20, int ny = 20) ;


/** Plots the sum of two vectors in a plane X=constant.
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv1 [input] first vector to define the field {\tt vv1 + vv2} to
 *		       be drawn.  
 *  @param vv2 [input] second vector to define the field {\tt vv1 + vv2} to
 *		       be drawn.  
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param y_min [input] lowest value of absol. coord. Y 
 *  @param y_max [input] highest value of absol. coord. Y 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf1 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   first surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   second surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param ny [input]  number of points in the Y direction (default value = 20)
 *  @param nz [input]  number of points in the Z direction (default value = 20)
 */
void des_vect_bin_x(const Tenseur& vv1, const Tenseur& vv2, double x0, 
		    double scale, double sizefl, double y_min, double y_max, 
		    double z_min, double z_max, char* title, 
		    const Cmp* defsurf1 = 0x0,  const Cmp* defsurf2 = 0x0, 
		    bool draw_bound = true, int ny = 20, int nz = 20) ;


/** Plots the sum of two vectors in a plane Y=constant.
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv1 [input] first vector to define the field {\tt vv1 + vv2} to
 *		       be drawn.  
 *  @param vv2 [input] second vector to define the field {\tt vv1 + vv2} to
 *		       be drawn.  
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf1 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   first surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   second surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param nx [input]  number of points in the X direction (default value = 20)
 *  @param nz [input]  number of points in the Z direction (default value = 20)
 */
void des_vect_bin_y(const Tenseur& vv1, const Tenseur& vv2, double x0, 
		    double scale, double sizefl, double x_min, double x_max, 
		    double z_min, double z_max, char* title, 
		    const Cmp* defsurf1 = 0x0,  const Cmp* defsurf2 = 0x0, 
		    bool draw_bound = true, int nx = 20, int nz = 20) ;


/** Plots the sum of two vectors in a plane Z=constant.
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param vv1 [input] first vector to define the field {\tt vv1 + vv2} to
 *		       be drawn.  
 *  @param vv2 [input] second vector to define the field {\tt vv1 + vv2} to
 *		       be drawn.  
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param scale [input] controls the length of the drawn arrows; if {\tt scale}
 *			 is negative, the length is determined automatically by
 *			 the routine: \\
 *			{\tt scale = -1} : max. length = step of the rectangular
 *					   grid \\ 
 *			{\tt scale = -2} : max. length = 2* step of the rectangular
 *					   grid \\ 
 *			 etc...	   
 *  @param sizefl [input] size of the arrows extremities (standard value: 1)
 *  @param x_min [input] lowest value of absol. coord. X 
 *  @param x_max [input] highest value of absol. coord. X 
 *  @param y_min [input] lowest value of absol. coord. Y
 *  @param y_max [input] highest value of absol. coord. Y
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param defsurf1 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   first surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param defsurf2 [input] pointer on a {\tt Cmp} giving the definition of the 
 *			   second surface: the surface is defined as the location
 *			   where this {\tt Cmp} vanishes (default value = 0x0, 
 *			   corresponds to no surface plot)
 *  @param draw_bound [input] true for drawing the boundaries of the various
 *			      domains (default value = true)
 *  @param nx [input]  number of points in the X direction (default value = 20)
 *  @param ny [input]  number of points in the Y direction (default value = 20)
 */
void des_vect_bin_z(const Tenseur& vv1, const Tenseur& vv2, double x0, 
		    double scale, double sizefl, double x_min, double x_max, 
		    double y_min, double y_max, char* title, 
		    const Cmp* defsurf1 = 0x0,  const Cmp* defsurf2 = 0x0, 
		    bool draw_bound = true, int nx = 20, int ny = 20) ;


    //@}




    /** @name 3-D visualization via Iris Explorer
     */
    //@{     

/** Prepare a file for visualizing a {\tt Cmp} with Iris Explorer, in a 
 *  given domain. 
 * 
 *  @param uu	{\tt Cmp} to examine
 *  @param lz	index of the domain  
 *  @param filename  name of the output file, to be read by Explorer. 
 * 
 */
void des_explorer(const Cmp& uu, int lz, const char* filename) ; 

/** Prepare a file for visualizing a {\tt Cmp} with Iris Explorer, in a 
 *  given domain, including data for z<0, assuming that the {\tt Cmp}
 *  is symmetric with respect to z=0. 
 * 
 *  @param uu	{\tt Cmp} to examine
 *  @param lz	index of the domain  
 *  @param filename  name of the output file, to be read by Explorer. 
 * 
 */
void des_explorer_symz(const Cmp& uu, int lz, const char* filename) ; 

/** Prepare a file for visualizing spectral coefficients with Iris Explorer, 
 *  in a given domain. 
 * 
 *  @param uu	{\tt Cmp}, the coefficient of which are to be examined
 *  @param lz	index of the domain  
 *  @param filename  name of the output file, to be read by Explorer. 
 * 
 */
void des_explorer_coef(const Cmp& uu, int lz, const char* filename) ; 




/** Prepare a file for visualizing a star (class {\tt Etoile}) with 
 *  Iris Explorer.
 * 
 *  @param star	{\tt Etoile} to examine
 *  @param filename  beginning of the names of the output files for Explorer.
 *		     {\tt nz-1} files will be created, where {\tt nz} is the
 *		     total number of domains. Each file name is indexed by 
 *		     the domain index.  
 */
void des_explorer(const Etoile& star, const char* name) ;

/** Prepare a file for visualizing a star (class {\tt Etoile}) with 
 *  Iris Explorer, including data for z<0, assuming that all {\tt Cmp}
 *  are symmetric with respect to z=0.
 * 
 *  @param star	{\tt Etoile} to examine
 *  @param filename  beginning of the names of the output files for Explorer.
 *		     {\tt nz-1} files will be created, where {\tt nz} is the
 *		     total number of domains. Each file name is indexed by 
 *		     the domain index.  
 */
void des_explorer_symz(const Etoile& star, const char* name) ;


/** Prepare a file for visualizing a binary star (class {\tt Binaire}) with 
 *  Iris Explorer.
 * 
 *  @param bibi	    {\tt Binaire} to examine
 *  @param filename  beginning of the names of the output files for Explorer.
 *		     {\tt nz1 + nz2 -2} files will be created, where {\tt nz1} 
 *		     (resp. {\tt nz2}) is the total number of domains in the 
 *		     mapping of star 1 (resp. star 2). Each file name is 
 *		     indexed first by the star number (after the suffixe 
 *		     "st") and then by the domain index (after the suffixe
 *		     "dm").  
 */
void des_explorer(const Binaire& bibi, const char* name) ;

/** Prepare a file for visualizing a binary star (class {\tt Binaire}) with 
 *  Iris Explorer,  including data for z<0, assuming that all {\tt Cmp}
 *  are symmetric with respect to z=0.
 * 
 *  @param bibi	    {\tt Binaire} to examine
 *  @param filename  beginning of the names of the output files for Explorer.
 *		     {\tt nz1 + nz2 -2} files will be created, where {\tt nz1} 
 *		     (resp. {\tt nz2}) is the total number of domains in the 
 *		     mapping of star 1 (resp. star 2). Each file name is 
 *		     indexed first by the star number (after the suffixe 
 *		     "st") and then by the domain index (after the suffixe
 *		     "dm").  
 */
void des_explorer_symz(const Binaire& bibi, const char* name) ;


/** Prepare a file for visualizing a meridian slice of a {\tt Cmp} 
 *  in a given domain with Iris Explorer.
 *  This is to draw a function $z=f(r,\theta)$, where {\it f} is a
 *  ``slice'' in the plane $\phi$ = const of the function defined by
 *  uu (to use with the {\bf Lit\_scal2D} module of Explorer).
 *
 *  @param uu       [input]: {\tt Cmp} to examine
 *  @param nz       [input]: index of the domain
 *  @param k_phi    [input]: $\phi $ index
 *  @param filename [input]: name of the output file, to be read by Explorer
 *  @param scale    [input]: scale by which the value of the {\tt Cmp}
 *                           is to be multiplied
 */

void des_explorer(const Bin_ns_ncp& bibi, const char* name) ;

/** Prepare a file for visualizing a binary star (class {\tt Bin_ns_ncp}) with 
 *  Iris Explorer,  including data for z<0, assuming that all {\tt Cmp}
 *  are symmetric with respect to z=0.
 * 
 *  @param bibi	    {\tt Bin_ns_ncp} to examine
 *  @param filename  beginning of the names of the output files for Explorer.
 *		     {\tt nz1 + nz2 -2} files will be created, where {\tt nz1} 
 *		     (resp. {\tt nz2}) is the total number of domains in the 
 *		     mapping of star 1 (resp. star 2). Each file name is 
 *		     indexed first by the star number (after the suffixe 
 *		     "st") and then by the domain index (after the suffixe
 *		     "dm").  
 */
void des_explorer_symz(const Bin_ns_ncp& bibi, const char* name) ;


/** Prepare a file for visualizing a meridian slice of a {\tt Cmp} 
 *  in a given domain with Iris Explorer.
 *  This is to draw a function $z=f(r,\theta)$, where {\it f} is a
 *  ``slice'' in the plane $\phi$ = const of the function defined by
 *  uu (to use with the {\bf Lit\_scal2D} module of Explorer).
 *
 *  @param uu       [input]: {\tt Cmp} to examine
 *  @param nz       [input]: index of the domain
 *  @param k_phi    [input]: $\phi $ index
 *  @param filename [input]: name of the output file, to be read by Explorer
 *  @param scale    [input]: scale by which the value of the {\tt Cmp}
 *                           is to be multiplied
 */


void des_explorer2D(const Cmp& uu, int nz, const int k_phi,
		    const char* filename, const double scale = 1.) ;

    //@}
    



#endif
