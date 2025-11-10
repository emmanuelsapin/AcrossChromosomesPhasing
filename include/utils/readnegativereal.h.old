#include <inttypes.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>   

//
//  read an positive integer from the file file
// return INT32_MAX if eof
//	
double readnegativereal(FILE * file)
{	double ID=0;
	char carac;
	int nega=0;
	carac=getc(file);
	if (carac=='-') 
	{	nega=1; 	
		carac='0';
	};
	if (carac!=EOF) do 
	{	if (carac!=EOF && carac>47 && carac<58) ID=ID*10+carac-48;	
		carac=getc(file);
	} while (carac!=EOF && carac>47 && carac<58);
	if (carac!=EOF) if (carac=='.')
	{	float dec=10;
		do 
		{	carac=getc(file);
			if (carac!=EOF && carac>47 && carac<58) ID=ID+(carac-48)/dec;
			dec=dec*10;
		}	while (carac!=EOF && carac>47 && carac<58);
	};
	if (carac!=EOF) if (carac=='e')
	{	carac=getc(file);
		int negaexp=0;
		if (carac=='-') 
		{	negaexp=1; 	
			carac='0';
		};
		int exp=0;
		do 
		{	if (carac!=EOF && carac>47 && carac<58) exp=exp*10+carac-48;	
			carac=getc(file);
		}	while (carac!=EOF && carac>47 && carac<58);
		if (negaexp) exp=-1*exp;
		ID = ID * pow(10,exp);
	}	
	if (carac==EOF) ID=FLT_MAX;
	if (nega) ID=-1*ID;
	return(ID); 	
}