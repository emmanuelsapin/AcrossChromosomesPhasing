#include "readreal.h"

// Read a positive real number from the file
// Returns FLT_MAX if EOF
float readreal(FILE * file)
{	
	float ID=0;
	char carac;
	do 
	{	
		carac=getc(file);
		if (carac!=EOF && carac>47 && carac<58) ID=ID*10+carac-48;
	}	while (carac!=EOF && carac>47 && carac<58);
	if (carac==EOF) ID=FLT_MAX;
	else
	{	
		if (carac=='.')
		{	
			float dec=10;
			do 
			{	
				carac=getc(file);
				if (carac!=EOF && carac>47 && carac<58) ID=ID+(carac-48)/dec;
				dec=dec*10;
			}	while (carac!=EOF && carac>47 && carac<58);
		};
	}
	if (carac==EOF) ID=FLT_MAX;
	return(ID); 	
}


