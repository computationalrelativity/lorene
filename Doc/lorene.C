/*!\mainpage LORENE --- Reference manual
 *
 * <A HREF="http://www.lorene.obspm.fr/"> LORENE </A>
 * (Langage Objet pour la RElativit&eacute; Num&eacute;riquE)
 * is a C++ based language for numerical relativity.
 *
 * Lorene home page is 
 * <A HREF="http://www.lorene.obspm.fr/"> http://www.lorene.obspm.fr/ </A>
 *
 */

/**
 * \brief Main namespace containing all \b LORENE classes and functions.
 *
 * The namespace \c Lorene gathers all function and classes defined in the 
 * LORENE library. To use it, simply put a 
 \code{.cpp} 
 using namespace Lorene ;
 \endcode
 * after the include lines in your code.
 * \copyright E Gourgoulhon, P Grandclement, J-A Marck, J Novak, K Taniguchi
 */
namespace Lorene {}

/**
 * \defgroup util Utilities.
 *
 * The classes \c Tbl and \c Itbl implement the 1-D, 2-D and 3-D
 * array representation in Lorene, whereas the class \c Matrice
 * is devoted to matrix representation and \c Param to generic
 * parameter storage.  
 * All these classes are independent of the numerical method 
 * (spectral method).
 *
 */

/**
 * \defgroup spec Spectral representation.
 *
 * These classes correspond to the implementation of spectral
 * methods in Lorene. They describe the collocation points
 * (classes \c Grille3d  and \c Mg3d ), the values
 * of a function at these points (class \c Mtbl ), and the coefficients
 * of the spectral expansions of a function (class \c Mtbl_cf ). 
 *
 */

/**
 * \defgroup map Mapping grid -> physical space (spherical coordinates)
 *
 * These classes implement the mapping between the
 * grid coordinates \f$(\xi, \theta', \phi')\f$ (described by the \c Mg3d
 * class) and the physical coordinates
 * \f$(r, \theta, \phi)\f$ [cf. Bonazzola, Gourgoulhon \& Marck, \e Phys. 
 * \e Rev. \e D \b 58 , 104020 (1998)].
 *
 * The class \c Map and its derived classes determine the methods for
 * partial derivatives with respect to the physical coordinates, as well
 * as resolution of basic partial differential equations (e.g. Poisson
 * equations).
 *
 */
  

/**\defgroup tensor Tensorial fields
 *
 * These classes implement the tensorial calculus in Lorene. 
 * They are high level classes and therefore are independent of the
 * actual numerical method (spectral method). 
 *
 */
    

/**\defgroup evol Time evolution (***under development***)
 *
 * Classes for time evolving fields and spacetimes.
 * Two families of classes are provided: (i) the \c Evolution family, 
 * which are template classes to store and manipulate (e.g. taking time 
 * derivatives) any evolving \e Lorene structure which has some 
 * arithmetics (e.g. \c Tbl, \c Tensor, etc...); (ii) the
 * \c Time_slice family which is devoted to the evolution of a hypersurface
 * \e t = const of the 3+1 formalism of General Relativity. 
 */


/**\defgroup otens Old tensorial fields (*** Deprecated ***)
 *
 * These classes have been used up to 2003 to treat scalar 
 * (class \c Cmp ) and tensorial fields (class \c Tenseur ). 
 * They are now deprecated and
 * have been replaced by the class \c Tensor and its various derived
 * classes, among which \c Scalar .
 */

/**\defgroup eos Equations of state
 *
 */
 
/**\defgroup star Stars and black holes
 *
 */
 
/**\defgroup compactobjects Stationary compact objects (***under development***)
 *
 */
 

/**\defgroup ellip General PDE solvers (***under development***)
 *
 * These classes are needed for using the general partial differential
 * equation solver, for
 * which the variables and the operators can be different from one domain
 * to the other.
 */


/**\defgroup mdm Grid wedding
 *
 * These classes are used to make the interface with Godunov-type
 * methods, used by the Garching and Valencia groups.
 *
 */

/**
 * \defgroup unites Physical units.
 *
 */


/**
 * \defgroup graphics Graphical outputs.
 *
 * The 2-D graphical outputs of various Lorene objects 
 * are performed via the PGPLOT library:
 * <A HREF="http://astro.caltech.edu/~tjp/pgplot/">http://astro.caltech.edu/~tjp/pgplot/</A>
 *
 * To open an X11 display area, use the command 'pgdisp'.
 *
 * 3-D visualization of various Lorene objects is performed by
 * means of OpenDX:
 * <A HREF="http://www.opendx.org/">http://www.opendx.org/</A>
 *
 * OpenDX is called from the methods Scalar::visu_section, 
 * \c Scalar::visu_box, \c Vector::visu_arrows , etc...
 * The corresponding OpenDX scripts are provided in the 
 * directory \c Lorene/Visu/OpenDX. They must be copied to the
 * working directory in order to use the above methods.  
 *
 */


