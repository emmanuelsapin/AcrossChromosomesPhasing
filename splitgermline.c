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
#define NMAX 700 
#define MAXPOP NMAX*NMAX
#define NBGROUPPERPERSON NMAX+1

#define MAXLASTGENERATION 51000
#define INT32 int32_t
#define INT64 int64_t
#define TOTMAXNUMSEG 8000000


INT32 groups[MAXPOP][NBGROUPPERPERSON];
INT32 nbgroups[MAXPOP];
INT32 groupnumber=0;
INT32 n;
INT32 ntop;
INT32 sqrtntop;
INT32 nbgroupperperson;

//
//
//

//
//  read file that resulkt from germline
//	574 33574 3076339.000000

//
//  save file sibling  
//

INT32 determinenbpeople(char path[], INT32 chr, INT32 set, INT32 gen)
{	char filename[200]; 
	FILE * file;
	strcpy (filename,path);
	strcat(filename,"/GE.gen");
	char number[100];
    sprintf(number, "%d", gen);
	strcat(filename,number);
	strcat(filename,".39subchr.chr");
	sprintf(number, "%d", chr+1);
	strcat(filename,number);
	strcat(filename,".f.");
	sprintf(number, "%d", set);
	strcat(filename,number);
	strcat(filename,".ped");
	file = fopen(filename, "r"); 
    // Check if file exists 
    if (file == NULL) 
    {   printf("Could not open file %s", filename); 
        return 0; 
    }; 
    // Extract characters from file and store in character c 
    char c;
	INT32 count;
	for (c = getc(file); c != EOF; c = getc(file)) 
    {   if (c == '\n') // Increment count if this character is newline 
        count = count + 1; 
	}
	// Close the file 
	fclose(file); 
	return (count-1);	
}



INT32 fillgroups()
{	INT32 n1;
	INT32 group;
	for(n1=0;n1<ntop;n1=n1+sqrtntop)
	{	INT32 n2;
		for(n2=n1;n2<n1+sqrtntop;n2++)
		{	groups[n2][nbgroups[n2]]=groupnumber;
			nbgroups[n2]++;
		/*	if (n2<100)
			{	printf("\n %d and %d",n2,groupnumber);
				printf("Creation of the blocks. There are %d blocks",groupnumber );
				printf("person 0 are in the %d groups:\n",nbgroups[0]);
				for(group=0;group<NBGROUPPERPERSON;group++)
				{	printf("%d ",groups[0][group]);
				};
				printf("\nperson 1 are in the %d groups:\n",nbgroups[1]);
				for(group=0;group<NBGROUPPERPERSON;group++)
				{	printf("%d ",groups[1][group]);
				};
				printf("\n");
			};*/
			
		};
		groupnumber++;
	};
	
	for(n1=0;n1<sqrtntop;n1++)
	{	INT32 n2;
		for(n2=0;n2<sqrtntop;n2++)
		{	groups[n1+n2*sqrtntop][nbgroups[n1+n2*sqrtntop]]=groupnumber;
			nbgroups[n1+n2*sqrtntop]++;
		};
		groupnumber++;
	};
				/*printf("Creation of the blocks. There are %d blocks",groupnumber );
				printf("person 0 are in the %d groups:\n",nbgroups[0]);
				for(group=0;group<NBGROUPPERPERSON;group++)
				{	printf("%d ",groups[0][group]);
				};
				printf("\nperson 1 are in the %d groups:\n",nbgroups[1]);
				for(group=0;group<NBGROUPPERPERSON;group++)
				{	printf("%d ",groups[1][group]);
				};
				printf("\n");*/
				
	INT32 moveup[sqrtntop];
	for(n1=0;n1<sqrtntop;n1++)
	{	moveup[n1]=0;
	};
	INT32 increment=1;
	do 
	{	for(n1=1;n1<sqrtntop;n1++)
		{	moveup[n1]=(n1*increment )%(sqrtntop);
		//	printf ("%d ",moveup[n1]);
		};
	//	printf("\n");
		for(n1=0;n1<sqrtntop;n1++)
		{	INT32 n2;
			for(n2=0;n2<sqrtntop;n2++)
			{	groups[(n1+moveup[n2])%sqrtntop+n2*sqrtntop ][nbgroups[(n1+moveup[n2])%sqrtntop+n2*sqrtntop]]=groupnumber;
				nbgroups[(n1+moveup[n2])%sqrtntop+n2*sqrtntop]++;
				//if (n1<10) printf("%d ", (n1+moveup[n2])%sqrtntop+n2*sqrtntop);
			};
		/*	if (n1<10) 
			{	printf("\n");
				for(n2=1;n2<sqrtntop;n2++)
				{	printf("%d ",moveup[n2] );
				};
				printf(" \n");
			}*/
			
			groupnumber++;
		};
		increment++;
	/*	printf("Creation of the blocks. There are %d blocks",groupnumber );
				printf("person 0 are in the %d groups:\n",nbgroups[0]);
				for(group=0;group<NBGROUPPERPERSON;group++)
				{	printf("%d ",groups[0][group]);
				};
				printf("\nperson 1 are in the %d groups:\n",nbgroups[1]);
				for(group=0;group<NBGROUPPERPERSON;group++)
				{	printf("%d ",groups[1][group]);
				};
				printf("\n");*/
				
	} while (increment<sqrtntop);
};

createfill(char path[], INT32 lowerbond, INT32 higherbond, INT32 set, INT32 gen, INT32 chr)
{	char filename[22][100]; 
	FILE * file[22];
	strcpy (filename[chr],path);
	strcat(filename[chr],"/GE.gen");
	char number[100];
	sprintf(number, "%d", gen);
	strcat(filename[chr],number);
	strcat(filename[chr],".39subchr.chr");
	sprintf(number, "%d", chr+1);
	strcat(filename[chr],number);
	strcat(filename[chr],".f.");
	sprintf(number, "%d", set);
	strcat(filename[chr],number);
	strcat(filename[chr],".ped");
	file[chr] = fopen(filename[chr], "r"); 
// Check if file exists 
	if (file[chr] == NULL) 
	{   printf("Could not open file %s", filename); 
		return 0; 
	};
	char filenamechunk[200]; 
	strcpy(filenamechunk,path);
	strcat(filenamechunk,"/chunk");
	INT32 person;
	for(person=0;person<n;person++)
	{	if (person%1000==0) printf("%d\n",person);
		char str[110000];
		fgets (str, 110000, file[chr]);
		//printf("%s",str);
		INT32 group;
	//	printf("person of line %d put in subgroup: ",person);
		for(group=0;group<nbgroupperperson;group++)
		{	//printf("%d\n",groups[person][group]);
			if (groups[person][group]<higherbond && groups[person][group]>lowerbond-1)
			{	//printf("%d, ",groups[person][group]);
				char filenamechunkn[200]; 
				strcpy (filenamechunkn,filenamechunk);
				char number[100];
				sprintf(number, "%d", groups[person][group]);
				strcat(filenamechunkn,number);
				strcat(filenamechunkn,"/chr");
				sprintf(number, "%d",chr+1);
				strcat(filenamechunkn,number);
				strcat(filenamechunkn,".ped");
				FILE * chunk;
				chunk = fopen(filenamechunkn, "a"); 
				if (chunk == NULL) 
				{   printf("Could not open file %s", filenamechunkn); 
					return 0; 
				};
				fputs (str,chunk);
				fclose(chunk);
			};
		};
	//	printf("\n");
	};
	//for(chr=0;chr<22;chr++)
	{	fclose(file[chr]);
	};
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
	{   printf("The path of the ped files is missing. Default is /rc_scratch/emsa1620/pedigree/testoutput2g50bits3Hextend1/STEP4_TGF/ped.phased/2000.\n");
		strcpy(path,"/rc_scratch/emsa1620/pedigree/testoutput2g50bits3Hextend1/STEP4_TGF/ped.phased/2000");
	} else 
	{	strcpy(path,argv[6]);
	};
	printf("The path of the ped file is %s\n",path);
	INT32 full;
    if (argc<9)
	{   printf("argument full missing. To run the program the arguments --set set where set is the size of the sample. full is 1\n");
		full=1;
	} else 
	{	full = atoi (argv[8]);
	};
	printf("full is %d\n",full);
	INT32 lowerbond;
    if (argc<11)
	{   printf("argument lowerbond missing. Default is 0\n");
		lowerbond=0;
	} else 
	{	lowerbond = atoi (argv[10]);
	};
	printf("lowerbond is %d\n",lowerbond);
	INT32 higherbond;
    if (argc<13)
	{   printf("argument higherbond missing. Default is 1\n");
		higherbond=1;
	} else 
	{	higherbond = atoi (argv[12]);
	};
	printf("higherbond is %d\n",higherbond);
	INT chr;
    if (argc<15)
	{   printf("argument chr missing. Default is 1\n");
		chr=1;
	} else 
	{	chr = atoi (argv[14])-1;
	};
	printf("chr is %d\n",chr);
	INT nfiction;
    if (argc<17)
	{   printf("argument nfiction missing. Default is 1\n");
		nfiction=1;
	} else 
	{	nfiction = atoi (argv[16]);
	};
	printf("nfiction is %d\n",nfiction);
	
	INT32 person;
	INT32 group;
	n=determinenbpeople(path,chr,set,45);
//	n=nfiction;
//	INT32 fact=0;
//	do 
//	{	fact++;
	//	n=nfiction*fact;
		printf("number of individuals is %d\n",n);
		INT32 i=ceil(sqrt(n));
		INT32 flag=0;
		do
		{	flag=0;
			INT32 j;
			for(j=2;j<ceil(i/2)+1;j++)
			{	if(i%j==0)
				{	flag=1;
					j=i-1;
				};
			};
			if (flag==1) i++;
		} while (flag==1 && i<10000); 
		printf("lowest square of a prime number higher than %d is %d\n",n,i*i);
		ntop=i*i;
		sqrtntop=i;
		nbgroupperperson=sqrtntop+1;
		printf("n = %d , ntop = %d , sqrtntop or nbpersonpergroup = %d , nbgroupperperson = %d, nbgroups = %d \n", n, ntop, sqrtntop, nbgroupperperson,ntop+sqrtntop);
		if (lowerbond>higherbond)
		{	char filename[200]; 
			FILE * file;
			strcpy (filename,path);
			strcat(filename,"/nbgroups.txt");
			file = fopen(filename, "w"); 
			// Check if file exists 
			if (file == NULL) 
			{   printf("Could not open file %s", filename); 
				return 0; 
			}; 
			fprintf(file,"%d\n",ntop+sqrtntop);
			fclose (file);
		};
		for(person=0;person<MAXPOP;person++)
		{	nbgroups[person]=0;
			for(group=0;group<NBGROUPPERPERSON;group++)
			{	groups[person][group]=0;
			};
		};
		printf("Creation of the blocks\n");
		fillgroups();
		printf("Creation of the blocks. There are %d blocks\n",groupnumber );
	/*	INT nbwrong=0;
		for(group=0;group<nbgroupperperson;group++)
		{	INT32 group1;
			for(group1=0;group1<nbgroupperperson;group1++)
			{	if (groups[2084][group]==groups[963][group1]) printf ("Individuals 23803 and 2028 interact in group %d which is the %d for 2048 and the %d for 963\n",groups[2084][group],group,group1 );
			};
		};		*/	
				
				 
	/*	for(person=0;person<n;person++)
		{	
			printf("Check if person %d interact one and only one time with the rest of the dataset\n",person);
			INT32 person2;
			for(person2=person;person2<n;person2++)
			{	INT32 nbinter=0;
				for(group=0;group<nbgroupperperson;group++)
				{	INT32 group1;
					for(group1=0;group1<nbgroupperperson;group1++)
					{	if (person!=person2 && groups[person][group]==groups[person2][group1]) nbinter++;
						if (person==person2 && group!=group1 && groups[person][group]==groups[person2][group1]) nbinter++;
							
						
					};
				};
				if (person2!=person)
				{	if (nbinter!=1) 
					{	printf("the pair %d and %d interacts %d times\n",person,person2,nbinter);
						nbwrong++;
					};
				} else 
				{	if (nbinter!=0) 
					{	printf("the pair %d and %d interacts %d times\n",person,person2,nbinter);
						nbwrong++;
					};
				};
				if (nbwrong>100) return (0);
			};
		};*/
	//} while (	nfiction*fact<1000000);
	printf("%d,%d",lowerbond,higherbond);
	createfill(path,lowerbond,higherbond,set, 45,chr);
	clock_t end = clock();
	float elapsed_secs = (float)(end - begin);
    printf("\nTotal time:%f cpu click so %f seconds",elapsed_secs,elapsed_secs/CLOCKS_PER_SEC );
	return (0);
} 