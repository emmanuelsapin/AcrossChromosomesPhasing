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
#define MAXLASTGENERATION 50050
#define INT int32_t
#define TOTMAXNUMSEG 8000000
#define MAXFILEDESTINATION 100
#define MAXNUMSEG TOTMAXNUMSEG/MAXFILEDESTINATION

INT individual[MAXPOP][MAXPOP];

typedef struct 
{	INT ID1;
	INT ID2;
	INT start; 
	INT end; 
	INT chr;
	INT lenght;	
} structsegment;

structsegment popsegment[MAXFILEDESTINATION][MAXNUMSEG];

typedef struct
{	INT lenght;	
	INT nbsegment;
} structarms;

typedef struct
{	structarms arms[2];	
} structchrs;

//
//
//

//
//  read file that resulkt from germline
//	574 33574 3076339.000000

//
//  save file sibling  
//
INT emptyinfile(INT filedestination,char path[],INT nbsegmentperfile[])
{	char filename[200]; 
	FILE * file;
	strcpy (filename,path);
	strcat(filename,".filenumber");
	char number[100];
    sprintf(number, "%d", filedestination);
	strcat(pathfile,number);
	strcat(filename,".match");
    if ((file = fopen(filename, "a")) == NULL) return 1;
    INT nbsegment;
	for(nbsegment=0;nbsegment<nbsegmentperfile[filedestination];nbsegment++)
	{	popsegment[nbsegment].ID1

	};
	fprintf(file,"%d %d %f\n",indiv1,indiv2,size);
	fclose (file);
	return (0);
}
//
//  get the number of the file destination when the two ID are integers
//	
INT getfiledestinationfrominteger(INT ID1,INT ID2)
{	INT sumID=ID1+ID2;
	do 
	{	sumID=sumID+ID1%10;
		ID1=trunc(ID1/10);
	} while (ID1>10);
	do 
	{	sumID=sumID+ID2%10;
		ID2=trunc(ID2/10);
	} while (ID2>10);
	return (sumID%MAXFILEDESTINATION);
}
//
//  read file that resulkt from germline
//	
INT readfilegermline(INT cmfin, INT set,char path[],INT chr,char arm,INT nbsegmentperfile[],structchrs chrs[])
{	FILE * file;
	char pathfile[200];
	strcpy(pathfile,path);
	char number[100];
    sprintf(number, "%d", chr);
	strcat(pathfile,number);
	strcat(pathfile,".");
	number[0]=arm;
	number[1]='\0';
	strcat(pathfile,number);
	strcat(pathfile,".");
	sprintf(number, "%d", set);
	strcat(pathfile,number);
	strcat(pathfile,".");
	sprintf(number, "%d", cmfin);
	strcat(pathfile,number);
	strcat(pathfile,".GERMLINE2.match");
	if ((file = fopen(pathfile, "r")) == NULL) 
	{	printf("file %s no found\n",pathfile);		
		return (1);
	}
    printf("file %s opened\n",pathfile);
	char carac;
  //  do {carac=getc(file);} while (carac!='\n' && carac!= EOF);
	INT ID1=0;
	INT ID2=0;
	INT print=0; 
	INT newindivi=0;
	INT nbline=0;
	do
	{	nbline++;
		ID1=readinteger(file);
		if (ID1<MAXINT)
		{	
			if (ID1<500000)
			{	ID1=ID1%100000;
			} else 
			{	ID1=ID1%100000+55000;
			};
			//printf("ID1 %d",ID1);
			if (ID1>MAXPOP)
			{	printf("ID1 is higher than MAXPOP on line %d",chrs[chr].arms[arm=='p'?1:0].nbsegment);
				return (0);
			};
			ID2=readinteger(file);
			if (ID2<500000)
			{	ID2=ID2%100000;
			} else 
			{	ID2=ID2%100000+55000;
			};
			//printf(" ID2 %d",ID2);
			if (ID2>MAXPOP)
			{	printf("ID2 is higher than MAXPOP on line %d",chrs[chr].arms[arm=='p'?1:0].nbsegment);
				return 0;
			};
			
			INT filedestination=getfiledestinationfrominteger(ID1,ID2);
		//	printf(" ID2 %d",ID2);
			popsegment[filedestination][nbsegmentperfile[filedestination]].ID1=ID1;
			popsegment[filedestination][nbsegmentperfile[filedestination]].ID2=ID2;
			popsegment[filedestination][nbsegmentperfile[filedestination]].start=readintegerbase94(file);
			popsegment[filedestination][nbsegmentperfile[filedestination]].end=readintegerbase94(file);
			popsegment[filedestination][nbsegmentperfile[filedestination]].chr=readintegerbase94(file);
			popsegment[filedestination][nbsegmentperfile[filedestination]].lenght=readintegerbase94(file);
	//		printf(" lenght %d\n",popsegment[filedestination][nbsegmentperfile[filedestination]].lenght);
			
			if (chrs[chr].arms[arm=='p'?1:0].lenght>chrs[chr].arms[arm=='p'?1:0].lenght+popsegment[filedestination][nbsegmentperfile[filedestination]].lenght)
			{	printf("error lenght can not be added anymore. Going over integer limit");
				return 0;
			};
			if (popsegment[filedestination][nbsegmentperfile[filedestination]].lenght<1000)
			{	if (individual[ID1][ID2]==0) 
				{	newindivi++;
					individual[ID1][ID2]=1;
				};
				chrs[chr].arms[arm=='p'?1:0].lenght=chrs[chr].arms[arm=='p'?1:0].lenght+popsegment[filedestination][nbsegmentperfile[filedestination]].lenght;
				chrs[chr].arms[arm=='p'?1:0].nbsegment++;
			};
			nbsegmentperfile[filedestination]++;
			if (nbsegmentperfile[filedestination]==MAXNUMSEG) 
			{	printf("Too many segments for file number %d",filedestination);
				emptyinfile(filedestination,path);
				nbsegmentperfile[filedestination]--;
			};
		};
	} while(ID1<MAXINT);
	printf("Arm %c of chromosome %d done:\n",arm,chr);
	printf("    %d segments found\n",chrs[chr].arms[arm=='p'?1:0].nbsegment);
	printf("    %d lenght of all segments\n",chrs[chr].arms[arm=='p'?1:0].lenght);
	printf("    %f average lenght of segments\n",((float) chrs[chr].arms[arm=='p'?1:0].lenght)/chrs[chr].arms[arm=='p'?1:0].nbsegment);
	printf("    %d new individuals found\n",newindivi);
	INT nbfile;
	for(nbfile=0;nbfile<MAXFILEDESTINATION;nbfile++) printf("%d segments put in file number %d\n",nbsegmentperfile[nbfile],nbfile);
	return (0);
}
	

//
//  main program
//	
int main(int argc, char *argv[])
{	 clock_t begin = clock();
	printf("start....\n");
	//
	printf("Reading arguments...\n");
	//
	INT cmfin;
	if (argc<3)
    {   printf("argument cmfin missing. To run the program the arguments --cmfin cmfin 	need to be added. Default is 5\n");
		cmfin=5;
	} else 
	{	cmfin = atoi (argv[2]);
	};
	printf("cmfin is %d\n",cmfin);
	INT set;
    if (argc<5)
	{   printf("argument set missing. To run the program the arguments --set set where set is the size of the sample. Default is 50000\n");
		set=50000;
	} else 
	{	set = atoi (argv[4]);
	};
	printf("set is %d\n",set);
	char path[200];
	if (argc<7)
	{   printf("The path of the run is missing. Default is /rc_scratch/emsa1620/GeneEvolve/Output_files300/.\n");
		strcpy(path,"/rc_scratch/emsa1620/GeneEvolve/Output_files300/");
	} else 
	{	strcpy(path,argv[6]);
	};
	printf("The path is %s\n",path);
	INT chr;
	INT nbsegmentperfile[MAXFILEDESTINATION];
	for(chr=0;chr<MAXFILEDESTINATION;chr++) nbsegmentperfile[chr]=0;
	INT chr1;
	for(chr=0;chr<MAXPOP;chr++) for(chr1=0;chr1<MAXPOP;chr1++) individual[chr][chr1]=0;
	structchrs chrs[23];
	for(chr=1;chr<23;chr++)
	{	chrs[chr].arms[0].lenght=0;
		chrs[chr].arms[1].lenght=0;
		chrs[chr].arms[0].nbsegment=0;
		chrs[chr].arms[1].nbsegment=0;
	};
	float totlenght=0;
	INT totngssegments=0;
	for(chr=1;chr<23;chr++)
	{	char arm='q';
		readfilegermline(cmfin,set,path,chr,arm,nbsegmentperfile,chrs);	 
		totlenght=totlenght+chrs[chr].arms[0].lenght;
		totngssegments=totngssegments+chrs[chr].arms[0].nbsegment;
		printf("total lenght of segment is %f\n",totlenght );
		printf("total number of segment is %d\n",totngssegments );	
		INT nbindividus=0;
		INT chr2;
		//for(chr2=0;chr2<MAXPOP;chr2++)  for(chr1=0;chr1<MAXPOP;chr1++) if (individual[chr2][chr1]) nbindividus++;
	//	printf("Average number of segments per pair is %f\n",((float) totngssegments)/(nbindividus>0?nbindividus:1));
		if (chr!=13 && chr!=14 && chr!=15 && chr!=21 && chr!=22) 
		{	arm='p';
			readfilegermline(cmfin,set,path,chr,arm,nbsegmentperfile,chrs);	
			totlenght=totlenght+chrs[chr].arms[1].lenght;
			totngssegments=totngssegments+chrs[chr].arms[1].nbsegment;
			printf("total lenght of segment is %f\n",totlenght );
			printf("total number of segment is %d\n",totngssegments );	
			nbindividus=0;
			//for(chr2=0;chr2<MAXPOP;chr2++)  for(chr1=0;chr1<MAXPOP;chr1++) if (individual[chr2][chr1]) nbindividus++;
		//	printf("Average number of segments per pair is %f (%d,%d)\n",((float) totngssegments)/(nbindividus>0?nbindividus:1),totngssegments,nbindividus);
		};
	};
	printf("Total number of segment is %d\n",totngssegments);
	printf("Total lenght of segment is %f\n",totlenght);
	printf("Average lenght of segment is %f\n",totlenght/(totngssegments>0?totngssegments:1));
	INT nbindividus=0;
	for(chr=0;chr<MAXPOP;chr++)  for(chr1=0;chr1<MAXPOP;chr1++) if (individual[chr][chr1]) 
	{	//printf("pair %d,%d\n",chr,chr1);
		nbindividus++;
	}
	printf("Average number of segments per pair is %f\n",((float) totngssegments)/(nbindividus>0?nbindividus:1));
	
	/*
	char path[100];
	if (argc<7)
	{   printf("The path of the run is missing. Default is /rc_scratch/emsa1620/GeneEvolve/Output_files300/.\n");
	} else 
	{	strcpy(path,argv[6]);
	};
	printf("The path is %s\n",path);
	INT numgen;
	if (argc<9)
	{   printf("All pairs will be considered from $generation$ to $generation-numgen+1$. Default is $numgen=1$.\n");
		numgen=1;
	} else 
	{	numgen = atoi (argv[8]);
	};*/	
	
	clock_t end = clock();
	float elapsed_secs = (float)(end - begin) ;
    printf("\nTotal time:%f cpu click so %f seconds",elapsed_secs,elapsed_secs/CLOCKS_PER_SEC );
	return (0);
} 