#include <inttypes.h>
#include <stdio.h>

#define INT int32_t
#define MAXINT INT32_MAX

//
//  read an positive integer from the file file
// return INT32_MAX if eof
//	
INT readintegerbase94(FILE * file)
{	INT ID=0;
	unsigned char carac;
	do 
	{	carac=getc(file);
		if (carac!=EOF && carac>32 && carac<127) ID=ID*94+carac-33;
	}	while (carac!=EOF && carac>32 && carac<127);
	if (carac==EOF) ID=MAXINT;
	return(ID); 	
}

void basetento94(int pos, char* str)
{	 char unit0=pos%94+33;
	pos=pos-pos%94;
	pos=pos/94;
	if (pos) 
	{	 char unit1=pos%94+33;
		pos=pos-pos%94;
		pos=pos/94;
		if (pos) 
		{	 char unit2=pos%94+33;
			pos=pos-pos%94;
			pos=pos/94;
			if (pos) 
			{	 char unit3=pos%94+33;
				pos=pos-pos%94;
				pos=pos/94;
				if (pos) 
				{	 char unit4=pos%94+33;
					str[0]=unit4;str[1]=unit3;str[2]=unit2;str[3]=unit1;str[4]=unit0;str[5]='\0';
				} else
				{	str[0]=unit3;str[1]=unit2;str[2]=unit1;str[3]=unit0;str[4]='\0';
				};
			} else
			{	str[0]=unit2;str[1]=unit1;str[2]=unit0;str[3]='\0';
			};
		}
		else
		{	str[0]=unit1;str[1]=unit0;str[2]='\0';
		};
	} else
	{	str[0]=unit0;str[1]='\0';
	};	
};