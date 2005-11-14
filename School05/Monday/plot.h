/*
 * Declaration of elementary graphical functions
 */

/*
 *   Copyright (c) 2005 Eric Gourgoulhon
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

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/11/14 01:57:00  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

#ifndef __PLOT_H_ 
#define __PLOT_H_ 

/**
 * \defgroup graphbasic Basic graphical routines.
 * Different figures can created simultaneously. They are characterized
 * by the label \c nfig. 
 *
 * @{
 */

/** Drawing a point as a small circle. 
 *
 *  @param x point abscissa
 *  @param y point ordinate
 *  @param color color of the small circle: \anchor color_list
 *      \li 0 : black (background)
 *      \li 1 : white (default)
 *      \li 2 : red
 *      \li 3 : green
 *      \li 4 : blue
 *      \li 5 : cyan
 *      \li 6 : magenta
 *      \li 7 : yellow
 *      \li 8 : orange
 *      \li 9 : green + yellow
 *      \li 10 : green + cyan
 *      \li 11 : blue + cyan
 *      \li 12 : blue + magenta
 *      \li 13 : red + magenta
 *      \li 14 : dark gray 
 *      \li 15 : light gray
 *  @param nfig index of the figure (in the range [0,99])
 *  to be used for the plot: if this figure does not exist, 
 *    it will be created with the device name \c device  provided by the last
 *      argument. 
 *  @param ymin lower bound on y of the graphical window (used only if a new 
 *      figure must be created)
 *  @param ymax upper bound on y of the graphical window (used only if a new 
 *      figure must be created)
 *  @param title title of the figure (used only if a new figure must be created)
 *  @param label_y y legend of the figure (used only if a new 
 *      figure must be created)
 *  @param device type of graphical device (default value = 0x0, will result in
 *  interactive choice) (used only if a new 
 *      figure must be created)
 */
void plot_point(double x, double y, int color = 1, int nfig = 0, 
                double ymin = -1., double ymax = 1., const char* title = 0x0, 
                const char* label_y = 0x0, const char* device = 0x0) ;

/** Drawing a set of points as small circles. 
 *
 *  @param np Number of points
 *  @param xx Array (size: \c np ) of abscissas of the points
 *  @param yy Array (size: \c np ) of ordinates of the points
 *  @param color color of the small circles: see \ref color_list 
 *  @param nfig index of the figure (in the range [0,99])
 *  to be used for the plot: if this figure does not exist, 
 *    it will be created with the device name \c device  provided by the last
 *      argument. 
 *  @param ymin lower bound on y of the graphical window (used only if a new 
 *      figure must be created)
 *  @param ymax upper bound on y of the graphical window (used only if a new 
 *      figure must be created)
 *  @param title title of the figure (used only if a new figure must be created)
 *  @param label_y y legend of the figure (used only if a new 
 *      figure must be created)
 *  @param device type of graphical device (default value = 0x0, will result in
 *  interactive choice) (used only if a new 
 *      figure must be created)
 */
void plot_point_set(int np, const double* xx, const double* yy, int color = 1, 
                    int nfig = 0, double ymin = -1., double ymax = 1., 
                  const char* title = 0x0, const char* label_y = 0x0, 
                  const char* device = 0x0) ;

/** Drawing a profile with uniform x sampling.
 *  A profile is a curve y=y(x). It is drawn on the specified figure. If the 
 *  latter is not opened, it will be opened with the device name \c device  
 *  provided by the last argument. 
 *
 *  @param yy Array (size: \c nx ) of y values to be drawn
 *			 (the x sampling is supposed to be uniform in [-1,1]).
 *  @param nx Number of points 
 *  @param color color of the profile: see \ref color_list 
 *  @param style style of the profile:
 *       the possible values are \c line_style[i] = 1  
 * (full line), 2 (dashed), 3 (dot-dash-dot-dash), 4 (dotted), 5 
 * (dash-dot-dot-dot).  
 *  @param nfig index of the figure (in the range [0,99])
 *  to be used for the plot: if this figure does not exist, 
 *    it will be created with the device name \c device  provided by the last
 *      argument. 
 *  @param ymin lower bound on y of the graphical window (used only if a new 
 *      figure must be created)
 *  @param ymax upper bound on y of the graphical window (used only if a new 
 *      figure must be created)
 *  @param title title of the figure (used only if a new figure must be created)
 *  @param label_y y legend of the figure (used only if a new 
 *      figure must be created)
 *  @param device type of graphical device (default value = 0x0, will result in
 *  interactive choice) (used only if a new 
 *      figure must be created)
 */
void plot_profile(const double* yy, int nx, int color = 1, int style = 1,    
                  int nfig = 0, double ymin = -1., double ymax = 1., 
                  const char* title = 0x0, const char* label_y = 0x0, 
                  const char* device = 0x0) ;

/** Closing a figure.
  *
  *  @param nfig [input] Index of the figure (in the range [0,99])
  * 
  */
void plot_close(int nfig = 0) ;

/** Closing all opened figures.
  *
  */
void plot_close_all() ;


/** @} */

/** Initialization of the graphical displays.
 * This function is called only by the other graphical routines
 * and is required if the graphical engine is PGPLOT.
 */
void plot_init() ; 

/** Opening a figure. 
 *  To be called only by the other graphical routines
 *
 */
void plot_open(int nfig, double ymin, double ymax, const char* title,    
    const char* label_y, const char* device) ; 

#endif
