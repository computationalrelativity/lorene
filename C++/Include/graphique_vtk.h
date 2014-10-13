/*
 *  Prototypes of graphical routines
 *
 */


#ifndef	__GRAPHIQUE_VTK_H_
#define	__GRAPHIQUE_VTK_H_

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace Lorene {
class Scalar ;
class Map ; 

#endif /* DOXYGEN_SHOULD_SKIP_THIS */ 

/**
 * \defgroup grafbasic  Low level graphical routines.
 *  \ingroup (graphics)
 * @{
 */

/*
 * The following are routines used to save Lorene data in the vtk format
 * for use with VisIt
 *
 */

/** Saves the data for a \c Scalar in the plane X=constant.
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the \c nzdes  innermost domains fit in it (for
 *			\c zoom  = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the \c nzdes  innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_vtk_x(const Scalar& uu, double x0, int nzdes, const char* title = 0x0,
                 double zoom = 1.2, int ny = 100, int nz = 100) ;

/** Saves the data for a \c Scalar in the plane X=constant
 *  within a specified graphic window. 
 *
 *  X is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param x0 [input] value of the coordinate X which defines the plane
 *		      of the drawing
 *  @param y_min [input] lowest value of absol. coord. Y 
 *  @param y_max [input] highest value of absol. coord. Y 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_vtk_x(const Scalar& uu, double x0, double y_min, double y_max,
                 double z_min, double z_max, const char* title = 0x0,
                 int ny = 100, int nz = 100) ;

/** Saves the data for a \c Scalar in the plane Y=constant.
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the \c nzdes  innermost domains fit in it (for
 *			\c zoom  = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the \c nzdes  innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_vtk_y(const Scalar& uu, double y0, int nzdes, const char* title = 0x0,
                 double zoom = 1.2, int nx = 100, int nz = 100) ;

/** Saves the data for a \c Scalar in the plane Y=constant
 *  within a specified graphic window. 
 *
 *  Y is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param y0 [input] value of the coordinate Y which defines the plane
 *		      of the drawing
 *  @param x_min [input] lowest value of absol. coord. x 
 *  @param x_max [input] highest value of absol. coord. x 
 *  @param z_min [input] lowest value of absol. coord. Z
 *  @param z_max [input] highest value of absol. coord. Z
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_coupe_vtk_y(const Scalar& uu, double y0, double x_min, double x_max,
                 double z_min, double z_max, const char* title = 0x0,
                 int nx = 100, int nz = 100) ;

/** Saves the data for a \c Scalar in the plane Z=constant.
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the \c nzdes  innermost domains fit in it (for
 *			\c zoom  = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param zoom [input] Factor by which the size of the graphic window 
 *			(determined from the size of the \c nzdes  innermost 
 *			 domains) is multiplied (default value = 1.2)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 */
void des_coupe_vtk_z(const Scalar& uu, double z0, int nzdes, const char* title = 0x0,
                 double zoom = 1.2, int nx = 100, int ny = 100) ;

/** Saves the data for a \c Scalar in the plane Z=constant
 *  within a specified graphic window. 
 *
 *  Z is the Cartesian coordinate relative to the absolute frame.  
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param z0 [input] value of the coordinate Z which defines the plane
 *		      of the drawing
 *  @param x_min [input] lowest value of absol. coord. X
 *  @param x_max [input] highest value of absol. coord. X
 *  @param y_min [input] lowest value of absol. coord. Y 
 *  @param y_max [input] highest value of absol. coord. Y 
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 */
void des_coupe_vtk_z(const Scalar& uu, double z0, double x_min, double x_max, 
                 double y_min, double y_max, const char* title = 0x0,
                 int nx = 100, int ny = 100) ;

/** Saves a three dimensional \c Scalar.
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param nzdes [input] number of domains for which the plot is performed:
 *			the size of the graphic window is determined so that
 *			the \c nzdes  innermost domains fit in it (for
 *			\c zoom  = 1.)
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_vtk_xyz(const Scalar& uu, int nzdes, const char* title = 0x0,
                int nx = 100, int ny = 100, int nz = 100) ;

/** Saves the data for a \c Scalar within a specified graphic window. 
 *
 *  @param uu [input] \c Scalar  to be drawn
 *  @param x_min [input] lowest value of absol. coord. X
 *  @param x_max [input] highest value of absol. coord. X
 *  @param y_min [input] lowest value of absol. coord. Y 
 *  @param y_max [input] highest value of absol. coord. Y 
 *  @param z_min [input] lowest value of absol. coord. Z 
 *  @param z_max [input] highest value of absol. coord. Z 
 *  @param title [input] title of the figure (default value = 0x0, 
 *			corresponds to no title)
 *  @param nx [input]  number of points in the X direction (default value = 100)
 *  @param ny [input]  number of points in the Y direction (default value = 100)
 *  @param nz [input]  number of points in the Z direction (default value = 100)
 */
void des_vtk_xyz(const Scalar& uu, double x_min, double x_max,
                 double y_min, double y_max, double z_min, double z_max,
                 const char* title = 0x0, int nx = 100, int ny = 100, int nz = 100) ;

/** @} */

}
#endif
