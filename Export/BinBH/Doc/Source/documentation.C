/**@name Meudon initial data for binary black holes
 *
 *  Meudon data represents quasistationary binary black configurations,
 *  obtained by P. Grandcl\'ement, E. Gourgoulhon \& S. Bonazzola,
 *   Phys. Rev. D {\bf 65}, 044021 (2002).
 *
 *  The exportation of this data, computed by means of
 *  \URL[LORENE]{http://www.lorene.obspm.fr/}
 *  on a multi-domain spectral grid, onto a Cartesian grid
 *  (e.g. for CACTUS), is performed by means of the C++ class {\tt Bin\_BH}.
 *  The class {\tt Bin\_BH} comes along with
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
 *	    Bin_BH binary_system(nb_points, xg, yg, zg, fill, datafile) ;
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
    	//@Include:bin_bh.h_r
//@}
