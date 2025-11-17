/***************************************************************************************************************************************************
*
*                                 This program outputs the common ancestor for any two people 
*(The speed of this program can be greatly decreased by writing on file a big bunch of data from one command instead of how it is done here)
*                                 
****************************************************************************************************************************************************/
#include <time.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "readinteger.h"
#include "readintegerbase94.h"
#include "readreal.h"

#define MAXINT INT32_MAX
#define MAXPOP 110000
#define MAXLASTGENERATION 51000
#define INT int32_t
#define SIZEMAXLINE 100000

INT subsets[2][MAXLASTGENERATION];

char line[SIZEMAXLINE];

//
// 
//	

INT readsubset(char path[],INT generation)
{	FILE *fp; 
    char filename[200]; 
    strcpy(filename,path);
	printf("The file to read the subset is %s", filename); 
	fp = fopen(filename, "r"); 
    if (fp == NULL) 
    {   printf("Could not open file %s", filename); 
        return (0); 
    }; 	
	INT ID=0;
	do
	{	ID=readinteger(fp);
		if (ID<MAXINT) if (ID<MAXLASTGENERATION) 
		{	subsets[generation][ID]=1;
		} 
		else 
		{	printf("ID %d higher than %d", ID,MAXLASTGENERATION);
		};
	} while (ID<MAXINT);
	fclose (fp);
	return (0);
};

//
// 
//

INT extractsample(INT chr,char arm,char pathorigin[],char pathoutput[], INT gen, INT last, INT set)
{	FILE *fp; 
    char filename[100]; 
    strcpy(filename,pathorigin);
	char number[100];
    sprintf(number, "%d", gen);
	strcat(filename,number);
	strcat(filename,".chr");
    sprintf(number, "%d", chr);
	strcat(filename,number);
	strcat(filename,".");
	filename[strlen(filename)+1]='\0';
	filename[strlen(filename)]=arm;
	strcat(filename,".ped");
	char fileoutput[100]; 
    strcpy(fileoutput,pathoutput);
	sprintf(number, "%d", gen-last);
	strcat(fileoutput,number);
	sprintf(number, "%d", gen-last+1);
	strcat(fileoutput,number);
	strcat(fileoutput,".39subchr.chr");
	sprintf(number, "%d", chr);
	strcat(fileoutput,number);
	strcat(fileoutput,".");
	fileoutput[strlen(fileoutput)+1]='\0';
	fileoutput[strlen(fileoutput)]=arm;
	strcat(fileoutput,".");
	sprintf(number, "%d", set);
	strcat(fileoutput,number);
	strcat(fileoutput,".ped");
	printf("\nThe file to output is %s\n", fileoutput); 
	FILE *fo; 
	if (last)
	{	fo = fopen(fileoutput, "a"); 
	}
	else
	{	fo = fopen(fileoutput, "w"); 
	}	
	if (fo == NULL) 
    {   printf("Could not open file %s", fileoutput);
        return 0; 
    };
    // Open the file 
    printf("The file to read is %s\n", filename); 
	fp = fopen(filename, "r"); 
    // Check if file exists 
    if (fp == NULL) 
    {   printf("Could not open file %s", filename);
        return 0; 
    };
	char c;
	char line[SIZEMAXLINE];
	INT lenghtline=0;
	INT numline=1;
	//subsets[last][1]=1;
	do 
	{	if (subsets[last][numline]==1)
		{	//fprintf(fo,	"g%d %d ", gen, numline);
			for (c = getc(fp); (c!='\n' && c!=EOF); c = getc(fp)) 
			{   if (c!=EOF)
				{	//printf("%c",c);
					line[lenghtline]=c;
					line[lenghtline+1]='\0';
					if (lenghtline+1==SIZEMAXLINE) 
					{	printf("line too long ( longer than %d caracteres",SIZEMAXLINE);
						return (0);
					};
					lenghtline++;
					
				}
			};
			printf("end line %d\n",numline);
			if (c == '\n') 
			{	line[lenghtline]=c;
				line[lenghtline+1]='\0';
				fprintf(fo,	"%s",line);
				lenghtline=0;
			};
			numline++;
		// Increment count if this character is newline 
        } 
		else 
		{	for (c = getc(fp); (c!='\n' && c!=EOF); c = getc(fp)) {};
			numline++;
		};
	} while (c != EOF);
	fclose (fp);
	fclose (fo);
	return (0);
}
//
//  main program
//	
int main(int argc, char *argv[])
{	 clock_t begin = clock();
	printf("start..\n");
	//
	printf("Reading arguments...\n");
	//
	INT set;
	if (argc<3)
    {   printf("argument set missing. To run the program the arguments --set set need to be added. Default is 1000\n");
		set=1000;
	} else 
	{	set = atoi (argv[2]);
	};
	printf("set is %d\n",set);
	char pathorigin[200];
    if (argc<5)
	{   printf("argument path origin is missing. To run the program the arguments --pathorigin pathorigin\n");
		strcpy(pathorigin,"/rc_scratch/emsa1620/GeneEvolveoutput/Output_files50000.2.save/ped.full.phased/");
	} else 
	{	strcpy(pathorigin,argv[4]);
	};
	printf("The path where are the ped file is %s\n",pathorigin);
	char pathsamplelastlast[200];
	if (argc<7)
	{   printf("The path of the run is missing. No default.\n");
		return 0;
	} else 
	{	strcpy(pathsamplelastlast,argv[6]);
	};
	printf("The path of the sample to take from one generation before the last is %s\n",pathsamplelastlast);
	char pathsamplelast[200];
	if (argc<9)
	{   printf("The path of the run is missing. No default.\n");
		return 0;
	} else 
	{	strcpy(pathsamplelast,argv[8]);
	};
	printf("The path of the sample to take from the last generation is %s\n",pathsamplelast);
	char pathoutput[200];
	if (argc<11)
	{   printf("The path of the output is missing. No default.\n");
		return 0;
	} else 
	{	strcpy(pathoutput,argv[10]);
	};
	printf("The path of the output is %s\n",pathoutput);
	INT gen;
	if (argc<13)
	{   printf("gen is missing. Default is 5.\n");
		gen=5;
	} else 
	{	gen = atoi (argv[12]);
	};
	printf("gen is %d\n",gen);
	INT chr;
	if (argc<15)
	{   printf("chr is missing. Default is 5.\n");
		chr=5;
	} else 
	{	chr = atoi (argv[14]);
	};
	printf("chr is %d\n",chr);
	char arm='q';
	if (argc<17)
	{   printf("arm is missing. Default is p.\n");
	} else 
	{	arm = argv[16][0];
	};
	printf("arm is %c\n",arm);
	INT generation;
	for (generation=0;generation<2;generation++)
	{	INT ID;
		for(ID=0;ID<MAXLASTGENERATION;ID++)
		{	subsets[generation][ID]=0;
		};
	}; 
	readsubset(pathsamplelastlast,0);
	readsubset(pathsamplelast,1);
	
	extractsample(chr,arm,pathorigin,pathoutput,gen-1,0,set);
	extractsample(chr,arm,pathorigin,pathoutput,gen,1,set);
	
	clock_t end = clock();
	float elapsed_secs = (float)(end - begin) ;
    printf("\nTotal time:%f cpu click so %f seconds",elapsed_secs,elapsed_secs/CLOCKS_PER_SEC );
	return (0);
} 