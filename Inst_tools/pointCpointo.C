/*
 * Code to generate *.o file names from *.C
 *
 */
#include <stdlib.h>
#include <stdio.h>

#ifdef OBSOLETE_HEADERS

#include <iostream.h>

#else

#include <iostream>
using namespace std ;
#endif 

int main() {
    char c ;
    int i ;

    cout << "OBJ =\t" ;
    
    while( (i = getchar()) != EOF) {
	c = char(i) ;
	switch (c) {
	    case '.':
	    c = getchar() ; // peut-etre le C
	    if (c == 'C') {
		cout << ".o" ;
	    }
	    else {
		cout << "." << c ;
	    }
	    break ;
	    
	    case ' ':
	    cout << " \\\n\t" ;
	    break ;
	    
	    default:
	    cout << c ;
	    break ;
	}
    }

    return 0 ;
}
