/*
 * Code to generate *.o file names from *.C
 *
 */

#include <iostream.h>

#include <stdlib.h>
#include <stdio.h>

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
