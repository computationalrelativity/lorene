/**@name Meudon initial data for binary black holes
 *
 *  The transfert of Meudon data, computed by means of LORENE
 *  on a multi-domain spectral grid, onto a Cartesian grid
 *  for CACTUS, is performed by means of the C++ class {\tt Bin\_BH}. 
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
 *	    Bin_BH binary_system(nb_points, xg, yg, zg, datafile) ;
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
    	//@Include:../bin_bh.h
//@}
