#include "readinteger.h"

// Read a positive integer from the file
// Returns INT32_MAX if EOF
INT readinteger(FILE * file)
{	
	INT ID=0;
	char carac;
	do 
	{	
		carac=getc(file);
		if (carac!=EOF && carac>47 && carac<58) ID=ID*10+carac-48;
	}	while (carac!=EOF && carac>47 && carac<58);
	if (carac==EOF) ID=MAXINT;
	return(ID); 	
}


