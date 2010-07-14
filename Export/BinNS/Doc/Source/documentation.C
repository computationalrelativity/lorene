/**@name Lorene initial data for binary neutron stars
 *
 *  Lorene data represents quasistationary binary neutron stars configurations,
 *  obtained by
 *   \begin{itemize}
 *      \item E. Gourgoulhon, P. Grandcl\'ement, K. Taniguchi, J.-A. Marck,
 *               S. Bonazzola, Phys. Rev. D {\bf 63}, 064029 (2001)
 *       \item K. Taniguchi, E. Gourgoulhon, S. Bonazzola,
 *               Phys.  Rev. D  {\bf 64}, 064012 (2001)
 *       \item K. Taniguchi, E. Gourgoulhon, Phys. Rev. D {\bf 65}, 044027 (2002)
 *		 \item K. Taniguchi, E. Gourgoulhon, Phys. Rev. D {\bf 66}, 104019 (2002)
 *		 \item K. Taniguchi, E. Gourgoulhon, Phys. Rev. D {\bf 68}, 124025 (2003) 
 *               \item M. Bejger, D. Gondek-Rosinska, E. Gourgoulhon, P. Haensel, 
 *                     K. Taniguchi, J.L. Zdunik,  Astron. Astrophys. {\bf 431}, 297 (2005)
 *    \end{itemize}
 *
 *  The exportation of this data, computed by means of
 *  \URL[LORENE]{http://www.lorene.obspm.fr/}
 *  on a multi-domain spectral grid, onto a Cartesian grid
 *  (e.g. for CACTUS), is performed by means of the C++ class {\tt Bin\_NS}.
 *  The class {\tt Bin\_NS} comes along with
 *  \URL[LORENE distribution]{http://www.lorene.obspm.fr/}.
 *  This class is very simple, with all data members being public.
 *  A typical example of use is the following one
 *
 *  \begin{verbatim}
 *	    // Define the Cartesian grid by means of the arrays xg, yg, zg:
 *	    for (int i=0; i<nb_points; i++) {
 *           xg[i] = ...
 *           yg[i] = ...
 *           zg[i] = ...
 *	    }
 *
 *	    // Read the file containing the spectral data and evaluate
 *	    //  all the fields on the Cartesian grid :
 *
 *	    Bin_NS binary_system(nb_points, xg, yg, zg, datafile) ;
 *
 *	    // Extract what you need :
 *
 *	    double* gamma_xx = binary_system.g_xx ; // metric coefficient g_xx
 *
 *	    double* shift_x = binary_system.beta_x ; // x comp. of shift vector
 *
 *	    ...
 *
 *	    // Save everything in an ASCII file :
 *
 *	    ofstream file_ini("ini.d") ;
 *	    binary_system.save_form(file_ini) ;
 *	    file_ini.close() ;
 *
 *  \end{verbatim}
 *
 */

//@{
    	//@Include:bin_ns.h_r
//@}
