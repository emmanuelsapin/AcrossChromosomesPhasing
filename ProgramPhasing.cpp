#include <time.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <errno.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sched.h>
#include <sys/syscall.h>
#include <omp.h>
#include "readinteger.h"
#include "readreal.h"
#include "readnegativereal.h"
#define NSNPPERCHR 60000
#define MAXPOP 435188
#define MAXCLOSERELAT 6000
#define MAXCLOSERELATTEMP 450000
#define MAXNBDIVISOR 25
#define NBINDIVMAX 100000
#define NBINDIV 100000
#define NBINDIVEA 1
#define MAXGEN 1
#define INCREMENTLOOP 1
#define MAXBREAK 1000

typedef struct
{
	int start;
	int end;
	int phasing;
	int nbright;
	int nbwrong;
	int segment;
} typechrdivider;

typedef struct
{	int chi;
	unsigned char haptohap;
} resultphaing;

typedef struct
{	int ID;
	unsigned char coef;
} indivclose;

typedef struct
{	int ID;
	int PIHAT10000;
} structrelatif;

typedef struct
{
	int chr;
	int snpstart;
	int snpend;
	float cor1;
	float cor2;
	int nblastsplit;
} structphaseerror;

typedef struct
{
	int IDoffspring;
	int IDp1;
	int IDp2;
	float score;
	int nbbad;
	int nbkeep;
} structtrio;

typedef struct
{	double cor;
	int nbgroup1;
	int nbgroup2;
} pointdecision;

typedef struct
{
	int IDoffspring;
	int IDp1;
	float score;
	int nbphaseerror;
} structPO;

typedef struct
{	int ID;
	int parent1;
	int parent2;
} structparentinfo;

typedef struct
{
	int ID;
	int start;
	int end;
	int chr;
	int hapfocal;
	int haprelat;
	int averageseg;
	int nbhet;
} structseg;

class PhasingProgram
{
public:
	int bestpihatagainstallID[100];
	int IDbestpihat;
	int IDbestpihat2;
	float pihatagainstall[MAXPOP];
	int nbchrdivider[51][23];
	int IDjob;
	int placefirttoconsider;
	typechrdivider chrdivider[50][23][20];
	float seuilpihat[3];
	int result[4][6];
	float bestpihat[7];
	int printdetail;
	int nbsnpperchr[23];
	int nbsnpperchrinfile[23];
	int MAF[NSNPPERCHR][23];
	int genomeoffpss[3][NSNPPERCHR][23];
	int NbIndiv;
	unsigned char * genomes[23];
	pointdecision tappointdec[MAXNBDIVISOR*22];
	structparentinfo parentinfo[NBINDIVMAX];
	int nbparentinfo;
	int listIndivToProcess[NBINDIVMAX];
	int nbIndivToProcess;
	int useListIndiv;
	int nbbreak;
	int distrigametic[12];
	int distrigametickeep[23];
	int checkrelat[MAXPOP];
	int segstart[24];
	int nbrelatpihat;
	
	PhasingProgram() : IDjob(0), placefirttoconsider(0), printdetail(0), NbIndiv(0),
		nbparentinfo(0), nbIndivToProcess(0), useListIndiv(0), nbbreak(0), nbrelatpihat(0)
	{
		for (int i = 0; i < NSNPPERCHR; i++)
		{
			for (int j = 0; j < 23; j++)
			{
				MAF[i][j] = 0;
			}
		}
		for (int i = 0; i < 23; i++)
		{
			genomes[i] = nullptr;
		}
	}
	
	int readgenomelocal(char pathfile[], int chr, int run, int step, unsigned char* gentomodify)
	{
		FILE* endfile;
		char number[100];
		char pathfilechr[300];
		strcpy(pathfilechr, pathfile);
		sprintf(number, "%d", chr);
		strcat(pathfilechr, number);
		strcat(pathfilechr, ".hap");
		printf("%s\n", pathfilechr);
		if ((endfile = fopen(pathfilechr, "r")) == NULL)
		{
			printf("file end is not found\n");
			return 1;
		}
		char first;
		int snp = 0;
		printf("Reading file: %s corresponding to chromosome %d\n", pathfilechr, chr);
		do
		{
			first = getc(endfile);
			if (first != EOF)
			{
				first = getc(endfile);
				int i = readinteger(endfile);
				
				do { first = getc(endfile); } while (first != 32 && first != EOF);
				i = readinteger(endfile);
				first = getc(endfile);
				getc(endfile);
				getc(endfile);
				if (first != EOF)
				{
					for (int relat = 0; relat < NbIndiv; relat++)
					{
						first = getc(endfile);
						int geno = 0;
						first = getc(endfile);
						if (first == 49) geno = geno + 2;
						first = getc(endfile);
						first = getc(endfile);
						first = getc(endfile);
						first = getc(endfile);
						if (first == 49) geno = geno + 1;
						*(gentomodify + (unsigned long long)(((((unsigned long long)relat) * (nbsnpperchr[chr] / 4 + ((nbsnpperchr[chr] % 4) > 0))) + snp / 4))) =
							(*(gentomodify + (unsigned long long)(((((unsigned long long)relat) * (nbsnpperchr[chr] / 4 + ((nbsnpperchr[chr] % 4) > 0))) + snp / 4))) &
							(~(3 << ((snp % 4) * 2)))) |
							((geno) << ((snp % 4) * 2));
						getc(endfile);
						getc(endfile);
					}
					snp++;
				}
			}
		} while (first != EOF);
		nbsnpperchrinfile[chr] = snp - 2;
		printf("%d SNPs in file\n", nbsnpperchrinfile[chr]);
		fclose(endfile);
		return 0;
	}
	int readParentInfo(const char* filename)
	{
		FILE* fp = fopen(filename, "r");
		if (fp == NULL)
		{
			printf("ERROR: Could not open parent info file %s\n", filename);
			return 1;
		}
		
		nbparentinfo = 0;
		char line[256];
		int lineNum = 0;
		
		while (fgets(line, sizeof(line), fp) != NULL && nbparentinfo < NBINDIVMAX)
		{
			lineNum++;
			if (line[0] == '\n' || line[0] == '\r')
				continue;
			
			int id, p1, p2;
			if (sscanf(line, "%d %d %d", &id, &p1, &p2) == 3)
			{
				if (id < 0 || id >= NBINDIVMAX)
				{
					printf("WARNING: Individual ID %d on line %d is out of range [0, %d). Skipping.\n",
						id, lineNum, NBINDIVMAX);
					continue;
				}
				parentinfo[nbparentinfo].ID = id;
				parentinfo[nbparentinfo].parent1 = p1;
				parentinfo[nbparentinfo].parent2 = p2;
				nbparentinfo++;
			}
			else
			{
				printf("WARNING: Invalid format on line %d. Expected: ID parent1 parent2\n", lineNum);
			}
		}
		
		fclose(fp);
		printf("Successfully loaded %d parent information entries from file %s\n", nbparentinfo, filename);
		return 0;
	}

	int getParentInfo(int indivID, int* parent1, int* parent2)
	{
		for (int i = 0; i < nbparentinfo; i++)
		{
			if (parentinfo[i].ID == indivID)
			{
				*parent1 = parentinfo[i].parent1;
				*parent2 = parentinfo[i].parent2;
				return 0;
			}
		}
		*parent1 = -1;
		*parent2 = -1;
		return 1;
	}

	int readListIndiv(const char* filename)
	{
		FILE* fp = fopen(filename, "r");
		if (fp == NULL)
		{
			printf("ERROR: Could not open file %s for reading individual list\n", filename);
			return 1;
		}
		
		nbIndivToProcess = 0;
		char line[256];
		int lineNum = 0;
		
		while (fgets(line, sizeof(line), fp) != NULL && nbIndivToProcess < NBINDIVMAX)
		{
			lineNum++;
			if (line[0] == '\n' || line[0] == '\r' || line[0] == '#')
				continue;
			
			int indivID = atoi(line);
			
			if (indivID < 0 || indivID >= NbIndiv)
			{
				printf("WARNING: Individual ID %d on line %d is out of range [0, %d). Skipping.\n",
					indivID, lineNum, NbIndiv);
				continue;
			}
			
			listIndivToProcess[nbIndivToProcess] = indivID;
			nbIndivToProcess++;
		}
		
		fclose(fp);
		
		if (nbIndivToProcess == 0)
		{
			printf("ERROR: No valid individual IDs found in file %s\n", filename);
			return 1;
		}
		
		printf("Successfully loaded %d individuals from file %s\n", nbIndivToProcess, filename);
		useListIndiv = 1;
		return 0;
	}

	int findrelative(int ID,int numtrio,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[])
{	int parametercombine=0; 
	int parametercalculfromparent=0; 
	int parameterfitness=2;
	int parameterlocal=0; 
	int tabseuil[10]; 
	tabseuil[0]=1;tabseuil[1]=2;tabseuil[2]=3;tabseuil[3]=4;tabseuil[4]=5;tabseuil[5]=6;tabseuil[6]=8;tabseuil[7]=11;tabseuil[8]=15;tabseuil[9]=20;
	int seuilacrosschr=0;
	int seuilinchr=3;
	int hetsquare=0;
	int prodadd=0;
	int prodaddpersnp=1;
	int restriclenseg=125;
	int keepone=0;
	float seuilpihat=0.33;
	int LD=0;
	int takelenasweight=0;
	int howcalculfreq=2;
	int MAFintoaccount=1;
	int limitweight=1;
	int limitweightproduct=13000; 
	int removelastchrsassoc=3;
	int addiftwicelimit=0;
	int hetsnp=0;
	int limitsumweightforindiv=0;
	int relatpihatP1[MAXCLOSERELATTEMP];
	int relatpihatP2[MAXCLOSERELAT];
	int relatpihatID[MAXPOP];
	int relatpihatchr[MAXCLOSERELAT][23][2];
	int seuilscore=1;
	int powersegment=1;	
	for(int IDrun=0;IDrun<MAXPOP;IDrun++)
	{	relatpihatID[IDrun]=-1;	
	};
	
	int64_t segnum=0;
	

	uint64_t averageofaverage[23];
	char number[100];
	char pathfile[300]; 								
	int nbseg=0;
	int segstarttemp[24];
	float bestpihatagainstall[100]={-1};
	for(int relat=0;relat<MAXPOP;relat++)
	{	pihatagainstall[relat]=0;
	};
	
	printf("Search for relatives\n");
	for(int chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	
		#pragma omp parallel for	
		for(int snp=0;snp<nbsnpperchrinfile[chrtemp1];snp++)		
		{	MAF[snp][chrtemp1]=0;
		};
	};
	for(int chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	
		#pragma omp parallel for	
		for(int snp=0;snp<nbsnpperchrinfile[chrtemp1];snp++)		
		{	for(int relat=0;relat<NbIndiv;relat++)
			{	int snpvalue0=(*((genomes[chrtemp1]+(unsigned long long) relat*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
				MAF[snp][chrtemp1]=MAF[snp][chrtemp1]+(snpvalue0>>1)+(snpvalue0&1);
			};
		};
	};
		int relat=ID;
	for(int chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	for(int snp=0;snp<nbsnpperchrinfile[chrtemp1];snp++)		
		{	
			float maffloat=1.0*MAF[snp][chrtemp1]/(NbIndiv)/2;
			if (maffloat>0.5) maffloat=1-maffloat;
			if (maffloat>0.01)
			{	int snpvalue0=(*((genomes[chrtemp1]+(unsigned long long) relat*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
				int parent0indiv0=(snpvalue0>>1);
				int parent1indiv0=(snpvalue0&1);
				float pcontribu10=((parent0indiv0)-maffloat);
				float pcontribu11=((parent1indiv0)-maffloat);
				float maffloatdiviseur=2*(maffloat)*(.5-maffloat/2);		
			
				float tabpihatcontri[3];
				tabpihatcontri[0]=pcontribu10*(-maffloat)*2/maffloatdiviseur+pcontribu11*(-maffloat)*2/maffloatdiviseur;
				tabpihatcontri[1]=pcontribu10*	(1-maffloat*2)/maffloatdiviseur+pcontribu11*	(1-maffloat*2)/maffloatdiviseur;		
				tabpihatcontri[2]=pcontribu10*	((2)-maffloat*2)/maffloatdiviseur+pcontribu11*	((2)-maffloat*2)/maffloatdiviseur;
				int multiplerelat=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
				int snpdiv4=snp/4;
				int snpdecal=(((snp%4)*2)); 
				printf("%d %d %f %f %f 	%f %f %f %f %f %f 	%d %d %d 	%f %f %f \n",chrtemp1,snpvalue0,pcontribu10,pcontribu11,pcontribu10+pcontribu11,	
											maffloat,maffloatdiviseur,-maffloat*2,2-maffloat*2,(pcontribu10+pcontribu11)*(2-maffloat*2),(pcontribu10+pcontribu11)*(-maffloat*2),	
											multiplerelat,snpdiv4,snpdecal,	tabpihatcontri[0],tabpihatcontri[1],tabpihatcontri[2]);
				#pragma omp parallel for			
				for(int relat2=0;relat2<NbIndiv;relat2++) 
				{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat2*multiplerelat )+snpdiv4)>>snpdecal)&3;
					int parent0indiv1=(snpvalue1&1);
					int parent1indiv1=(snpvalue1>>1);	
					pihatagainstall[relat2]=pihatagainstall[relat2]+tabpihatcontri[parent0indiv1+parent1indiv1];	
				};	
			};
		};	
	};
	for(int relat2=0;relat2<NbIndiv;relat2++)
	{	
		pihatagainstall[relat2]=pihatagainstall[relat2]/2/(330005)-1;
		 	printf("%d %f\n",relat2 ,pihatagainstall[relat2]);
	};
	return 0;
	for(int relat2=0;relat2<NbIndiv;relat2++)
	{	if (pihatagainstall[relat2]>0.1)
		{	printf("Relative of %d found: %d with pihat of %f\n",relat,relat2,pihatagainstall[relat2]);
		}
		int comparepihat=99;
		while (pihatagainstall[relat2]>bestpihatagainstall[comparepihat] && comparepihat>-1)
		{	comparepihat--;
		};
		if (comparepihat<99)
		{	for(int movepihat=99;movepihat>comparepihat;movepihat--)
			{	bestpihatagainstall[movepihat]=bestpihatagainstall[movepihat-1];
				bestpihatagainstallID[movepihat]=bestpihatagainstallID[movepihat-1];
			};
		}
		bestpihatagainstall[comparepihat]=pihatagainstall[relat2];
		bestpihatagainstallID[comparepihat]=relat2;
	};
	int compteur=0;
	for(int comp=0;comp<7;comp++) bestpihat[comp]=0;
	bestpihat[0]=bestpihatagainstall[compteur];
	int nbtimeIDchange=0;
	do {compteur++;} while (bestpihatagainstall[compteur]>seuilpihat);
	placefirttoconsider=compteur;
	bestpihat[1]=bestpihatagainstall[compteur];
	IDbestpihat=bestpihatagainstallID[compteur];
	IDbestpihat2=bestpihatagainstallID[compteur+1];
	for(int comp=compteur;comp<compteur+2;comp++)	bestpihat[2]=bestpihat[2]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+5;comp++)	bestpihat[3]=bestpihat[3]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+10;comp++)	bestpihat[4]=bestpihat[4]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+20;comp++)	bestpihat[5]=bestpihat[5]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+50;comp++)	bestpihat[6]=bestpihat[6]+bestpihatagainstall[comp];	
	segstarttemp[0]=0;
	
	float pihatagainstallchrMP[MAXPOP][2]={0};
	for(int relat=0;relat<MAXPOP;relat++) 
	{	pihatagainstallchrMP[relat][0]=0;pihatagainstallchrMP[relat][1]=0;
	}
	return (1);
}
	int predict_phasing_without_parents(int ID,int chr1, int chr2,int numtrio,float limitnbsnp,int run,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[])
{	printf("%d\n",IDjob);
	seuilpihat[0]=0.33;
	seuilpihat[1]=0.33;
	seuilpihat[2]=0.33;
	int increaseincrementwhenneg=1;
	int takesumintoaccout=1;
	int limitweight=105;
	int howtocalculatenotimprove=0;
	int howtocalculatenotimprove2=0;
	int nbmutdependonnotimprove=0;
	int mutateallhapwhennoimprove=0;
	int changeonlyifneg=2;
	int crossclever=0;
	int numgentowait=17;
	int howtocalculaterelattochr=0; 
	int nbmutvari=1;
	int64_t tabresmut[50][50];
	int nbtimetabresmut[50][50];
	for(int mut1=0;mut1<50;mut1++)
	{	for(int mut2=0;mut2<50;mut2++)
		{	tabresmut[mut1][mut2]=0;
			nbtimetabresmut[mut1][mut2]=1;
		};
	};
	int mostrelatinscore=0;
	int howcalculmostrelaetd=0;
	int howcalculmostrelaetd2=0;
	int addoneandzero=1;
	int localsearch=1;
	int randmuthap=0;
	unsigned char relatsuperpose[MAXCLOSERELAT][MAXCLOSERELAT];
	int64_t relatpihatchr[MAXCLOSERELAT][23][2];
	int seuilscore=1;
	int powersegment=1;		
	int64_t segnum=0;
	unsigned char relatbreak[MAXBREAK][MAXCLOSERELAT];
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	for(int  break1=0;break1<nbbreak;break1++)		
		{	relatbreak[break1][relat]=0;
		};
	};
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	int highestnbhet=0;
		int highestnbhetbreak=0;
		for(int  break1=0;break1<nbbreak;break1++)		
		{	if (highestnbhet<relatbreak[break1][relat]) 
			{	highestnbhet=relatbreak[break1][relat];
				highestnbhetbreak=break1;
			};
			if (relatbreak[break1][relat]) relatbreak[break1][relat]=1;
		};
		for(int  break1=0;break1<nbbreak;break1++)		
		{	relatbreak[break1][relat]=(highestnbhetbreak==break1);
		};
	};	
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	printf("%d ",relat);
		for(int  break1=0;break1<nbbreak;break1++)		
		{	printf("%d ",relatbreak[break1][relat]);
		};
	};
	int relattochr[MAXCLOSERELAT];
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	int bestchr=1;
		for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
		{	int hapmost=(relatpihatchr[relat][chrtemp1][0]>relatpihatchr[relat][chrtemp1][1])?0:1;
			int pihatnorm=relatpihatchr[relat][chrtemp1][hapmost]*nbsnpperchrinfile[1]/nbsnpperchrinfile[chrtemp1];
			if (chrtemp1==1 || relattochr[relat]<pihatnorm)  
			{	relattochr[relat]=pihatnorm;
				bestchr=chrtemp1+22*hapmost;
			};
		};
		relattochr[relat]=bestchr;
	};
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	for(int snp=0;snp<nbsnpperchrinfile[chrtemp1];snp++)		
		{	genomeoffpss[0][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			genomeoffpss[1][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			genomeoffpss[2][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) IDp2loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			printf("%d %d %d\n",chrtemp1,snp,genomeoffpss[0][snp][chrtemp1]);
		};	
	};
	int maxpihat=0;
	int64_t nbseghap[2][23][26229];
	int64_t nbsegoverchr[23][2];
	int ratio=2;	
	printf("maxpihat %d\n",maxpihat);
	FILE * endfile;
	uint64_t averageofaverage[23];
	char number[100];
	char pathfile[300]; 								
	uint64_t sumsumnincrementplus=0;
	uint64_t sumsumnincrementmoins=0;
	sumsumnincrementplus=0;
	sumsumnincrementmoins=0;
	float plusmoins=1.0*sumsumnincrementplus/sumsumnincrementmoins;
	float moinsplus=1.0*sumsumnincrementmoins/sumsumnincrementplus;
	printf("%f %f\n",plusmoins,moinsplus);
	sumsumnincrementplus=0;
	sumsumnincrementmoins=0;
	int64_t sumallpihatpos=0;
	int64_t sumallpihatneg=0;
	int nbpos=0;
	int nbneg=0;
	int64_t minpihat=0;
	int minpihatID1=0;
	int minpihatID2=0;
	unsigned char indivEA[NBINDIVEA][MAXCLOSERELAT];
	unsigned char bestindiv[MAXCLOSERELAT];
	int64_t scoreEA[NBINDIVEA];	
	int64_t bestscore=-1000000;
	int64_t nbsumerrorbestscore=0;
	int64_t bestscoreofbestfitness=0;
	int64_t bestphasing=0;
	int same;
	int gen=gentostart;	
	int lastgen=100;	
	int nuincreasenewtop=50;
	int64_t averageparentold5=-1000000;
	double averageparent[MAXGEN];
	int print=1;
	int lastgenimprove=gen;
	int genintiindiv=gen;
	int64_t bestscoresinceinit=bestscore;
	int64_t bestscorefitness=-100000;
	int64_t bestscoreclustering=-100000;
	do
	{	int64_t scorerelattab[NBINDIVEA][MAXCLOSERELAT];
		int64_t scorerelattaballchr[NBINDIVEA][MAXCLOSERELAT][MAXBREAK];
		double scorehap[NBINDIVEA][MAXBREAK];
		{	int groumdtruthdone=0;
			int bestwellphased=0;
			int bestunwellphased=0;
			int64_t score1=0;
			int64_t score2=0;
			int64_t score1tab[MAXBREAK];
			int64_t score2tab[MAXBREAK];			
			for(int  chrtemp1=0;chrtemp1<MAXBREAK;chrtemp1++)		
			{	score1tab[chrtemp1]=0;
				score2tab[chrtemp1]=0;
			}
			int nbonebalanced=0;
		
			int64_t score1averageneg=0;
			int score1nbneg=1;
			int64_t score2averageneg=0;
			int score2nbneg=1;
			int chagemade=1;
			int64_t iter=0;
			int iterbreak=0;
			int64_t maxscorerelattab=0;
			int64_t plus[MAXCLOSERELAT]={0};
			int64_t minus[MAXCLOSERELAT]={0};
			maxscorerelattab=0;
			int64_t maxscorerelattabID=0;
			int64_t score1temp=0;
			 int64_t  score2temp=0;
			printf("%d %d %d %d\n",score1,score2,score1nbneg,score1nbneg);
			int64_t score;
			int64_t averagescore1=1;
			int64_t averagescore2=1;
			if (nbonebalanced*2>nbrelatpihat) nbonebalanced=1+nbonebalanced*2-nbrelatpihat; else nbonebalanced=1+nbrelatpihat-nbonebalanced*2;
			if (bestscoresinceinit<score) 
			{
				lastgenimprove=gen;
				bestscoresinceinit=score;
			};
			int64_t highest=score;
			{	
				bestscoreclustering=score;
				if (bestscoreofbestfitness==0) bestscoreofbestfitness=1;
				char number[100];
				char pathfile[300]; 
				int64_t scorechr1[23];
				int64_t scorechr2[23];
				int64_t scorechr3[23];
				int64_t scorechr4[23];
				int keep[23][4];
				int64_t nbonefromhap[2][23][26229];
				int64_t nbzerofromhap[2][23][26229];
				int64_t nbonequalified[2];
				int64_t nbzeroqualified[2];
				nbonequalified[0]=0;
				nbonequalified[1]=0;
				nbzeroqualified[0]=0;
				nbzeroqualified[1]=0;
				int64_t sumdiff[23];
				int64_t nbcontracditperchr[23];
				for(int nbkickout=1;nbkickout<23;nbkickout++)
				{	keep[nbkickout][ratio]=1;
				};
				{	nbsumerrorbestscore=1;
					int phase1[23][1][1];
					int guessphase1[23][1][1];
					int phase1guessphase1[23][1][1];
					int phase1guessphase2[23][1][1];
					int phase2[23][1][1];
					int guessphase2[23][1][1];
					int phase2guessphase2[23][1][1];
					int phase2guessphase1[23][1][1];
					int phaseguessright[23][1][1];
					int phaseguesswrong[23][1][1];
					int phasenotguessed[23][1][1];
					int nbsum=1;
					int wrong[23];
					int right[23];
					unsigned int x, y, width, height;
					for(int prod=0;prod<1;prod++)
					{	for(int chr=0;chr<23;chr++)
						{	for(int ratio=0;ratio<1;ratio++)
							{	phaseguessright[chr][prod][ratio]=0;
								phaseguesswrong[chr][prod][ratio]=0;
								phasenotguessed[chr][prod][ratio]=0;
								phase1[chr][prod][ratio]=0;
								guessphase1[chr][prod][ratio]=0;
								phase1guessphase1[chr][prod][ratio]=0;
								phase1guessphase2[chr][prod][ratio]=0;
								phase2[chr][prod][ratio]=0;
								guessphase2[chr][prod][ratio]=0;
								phase2guessphase2[chr][prod][ratio]=0;
								phase2guessphase1[chr][prod][ratio]=0;
							};	
						};
					};
					int lenshowchr=height*0.9;
					int widhtshowchr=5;
					int gapshowchry=63;
					int widthshowchry=9;
					int gapshowchrx=5;
					int64_t sumdifftot=0;
					int phaseparent1=-1;
					int phaseparent1tap[23][26229];
					int phaseguesedparent1tap[23][26229];
					for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
					{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
						{	phaseparent1tap[chrtemp1][snp]=0;
							phaseguesedparent1tap[chrtemp1][snp]=0;
						};
					};
					int phaseunknow=0;	
					int bestpur=0;
					double dataperchr[23][110];		
					int nbdivision=2;
					int powerpihat=0;
					int exponentcorp=0;
					double sumoftabcor=0;
					int firstexponent;
					int thirdindice;
					double pihatagainstallchrMPphaseerror[23][MAXPOP][200][2];
					int purcentage;
					int powerpihatDEGREE4=0;
					for(int size=0;size<51;size++)
					{	for(int chr=1;chr<23;chr++)
						{	for(int window=0;window<1;window++)
							{	chrdivider[size][chr][window].segment=0;
							};
						};
					};
					chrdivider[25][1][0].end=2203;chrdivider[25][1][1].start=2203;
					chrdivider[25][1][1].end=7136;chrdivider[25][1][2].start=7136;
					chrdivider[25][1][2].end=14606;chrdivider[25][1][3].start=14606;
					chrdivider[25][1][3].end=18289;chrdivider[25][1][4].start=18289;
					chrdivider[25][1][4].end=nbsnpperchrinfile[1];
					nbchrdivider[25][1]=5;
					nbchrdivider[25][0]=nbchrdivider[25][1];
					chrdivider[25][2][0].end=2097;chrdivider[25][2][1].start=2097;
					chrdivider[25][2][1].end=6555;chrdivider[25][2][2].start=6555;
					chrdivider[25][2][2].end=12181;chrdivider[25][2][3].start=12181;
					chrdivider[25][2][3].end=18340;chrdivider[25][2][4].start=18340;
					chrdivider[25][2][4].end=nbsnpperchrinfile[2];
					nbchrdivider[25][2]=5;
					if (nbchrdivider[25][2]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][2];
					chrdivider[25][3][0].end=4737;chrdivider[25][3][1].start=4737;
					chrdivider[25][3][1].end=9553;chrdivider[25][3][2].start=9553;
					chrdivider[25][3][2].end=14838;chrdivider[25][3][3].start=14838;
					chrdivider[25][3][3].end=18519;chrdivider[25][3][4].start=18519;
					chrdivider[25][3][4].end=nbsnpperchrinfile[3];
					nbchrdivider[25][3]=5;
					if (nbchrdivider[25][3]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][3];
					chrdivider[25][4][0].end=3754;chrdivider[25][4][1].start=3754;
					chrdivider[25][4][1].end=6615;chrdivider[25][4][2].start=6615;
					chrdivider[25][4][2].end=9749;chrdivider[25][4][3].start=9749;
					chrdivider[25][4][3].end=16124;chrdivider[25][4][4].start=16124;
					chrdivider[25][4][4].end=nbsnpperchrinfile[4];
					nbchrdivider[25][4]=5;
					if (nbchrdivider[25][4]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][4];
					chrdivider[25][5][0].end=6946;chrdivider[25][5][1].start=6946;
					chrdivider[25][5][1].end=11357;chrdivider[25][5][2].start=11357;
					chrdivider[25][5][2].end=15071;chrdivider[25][5][3].start=15071;
					chrdivider[25][5][3].end=nbsnpperchrinfile[5];
					nbchrdivider[25][5]=4;
					if (nbchrdivider[25][5]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][5];
					chrdivider[25][6][0].end=3976;chrdivider[25][6][1].start=3976;
					chrdivider[25][6][1].end=11907;chrdivider[25][6][2].start=11907;
					chrdivider[25][6][2].end=nbsnpperchrinfile[6];
					nbchrdivider[25][6]=3;
					if (nbchrdivider[25][6]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][6];
					chrdivider[25][7][0].end=3082;chrdivider[25][7][1].start=3082;
					chrdivider[25][7][1].end=7188;chrdivider[25][7][2].start=7188;
					chrdivider[25][7][2].end=12684;chrdivider[25][7][3].start=12684;
					chrdivider[25][7][3].end=nbsnpperchrinfile[7];
					nbchrdivider[25][7]=4;
					if (nbchrdivider[25][7]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][7];
					chrdivider[25][8][0].end=7271;chrdivider[25][8][1].start=7271;
					chrdivider[25][8][1].end=12465;chrdivider[25][8][2].start=12465;
					chrdivider[25][8][2].end=nbsnpperchrinfile[8];
					nbchrdivider[25][8]=3;
					if (nbchrdivider[25][8]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][8];
					chrdivider[25][9][0].end=3384;chrdivider[25][9][1].start=3384;
					chrdivider[25][9][1].end=7020;chrdivider[25][9][2].start=7020;
					chrdivider[25][9][2].end=9803;chrdivider[25][9][3].start=9803;
					chrdivider[25][9][3].end=nbsnpperchrinfile[9];
					nbchrdivider[25][9]=4;
					if (nbchrdivider[25][9]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][9];
					chrdivider[25][10][0].end=2407;chrdivider[25][10][1].start=2407;
					chrdivider[25][10][1].end=8424;chrdivider[25][10][2].start=8424;
					chrdivider[25][10][2].end=12288;chrdivider[25][10][3].start=12288;
					chrdivider[25][10][3].end=nbsnpperchrinfile[10];
					nbchrdivider[25][10]=4;
					if (nbchrdivider[25][10]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][10];
					chrdivider[25][11][0].end=2953;chrdivider[25][11][1].start=2953;
					chrdivider[25][11][1].end=5344;chrdivider[25][11][2].start=5344;
					chrdivider[25][11][2].end=nbsnpperchrinfile[11];
					nbchrdivider[25][11]=3;
					if (nbchrdivider[25][11]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][11];
					chrdivider[25][12][0].end=5937;chrdivider[25][12][1].start=5937;
					chrdivider[25][12][1].end=9107;chrdivider[25][12][2].start=9107;
					chrdivider[25][12][2].end=14609;chrdivider[25][12][3].start=14609;
					chrdivider[25][12][3].end=nbsnpperchrinfile[12];
					nbchrdivider[25][12]=4;
					if (nbchrdivider[25][12]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][12];
					chrdivider[25][13][0].end=2524;chrdivider[25][13][1].start=2524;
					chrdivider[25][13][1].end=6229;chrdivider[25][13][2].start=6229;
					chrdivider[25][13][2].end=nbsnpperchrinfile[13];
					nbchrdivider[25][13]=3;
					if (nbchrdivider[25][13]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][13];
					chrdivider[25][14][0].end=3228;chrdivider[25][14][1].start=3228;
					chrdivider[25][14][1].end=6437;chrdivider[25][14][2].start=6437;
					chrdivider[25][14][2].end=nbsnpperchrinfile[14];
					nbchrdivider[25][14]=3;
					if (nbchrdivider[25][14]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][14];
					chrdivider[25][15][0].end=3288;chrdivider[25][15][1].start=3288;
					chrdivider[25][15][1].end=6043;chrdivider[25][15][2].start=6043;
					chrdivider[25][15][2].end=8861;chrdivider[25][15][3].start=8861;
					chrdivider[25][15][3].end=nbsnpperchrinfile[15];
					nbchrdivider[25][15]=4;
					if (nbchrdivider[25][15]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][15];
					chrdivider[25][16][0].end=3860;chrdivider[25][16][1].start=3860;
					chrdivider[25][16][1].end=nbsnpperchrinfile[16];
					nbchrdivider[25][16]=2;
					if (nbchrdivider[25][16]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][16];
					chrdivider[25][17][0].end=2196;chrdivider[25][17][1].start=2196;
					chrdivider[25][17][1].end=4676;chrdivider[25][17][2].start=4676;
					chrdivider[25][17][2].end=8643;chrdivider[25][17][3].start=8643;
					chrdivider[25][17][3].end=nbsnpperchrinfile[17];
					nbchrdivider[25][17]=4;
					if (nbchrdivider[25][17]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][17];
					chrdivider[25][18][0].end=3734;chrdivider[25][18][1].start=3734;
					chrdivider[25][18][1].end=7591;chrdivider[25][18][2].start=7591;
					chrdivider[25][18][2].end=nbsnpperchrinfile[18];
					nbchrdivider[25][18]=3;
					if (nbchrdivider[25][18]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][18];
					chrdivider[25][19][0].end=2982;chrdivider[25][19][1].start=2982;
					chrdivider[25][19][1].end=6035;chrdivider[25][19][2].start=6035;
					chrdivider[25][19][2].end=nbsnpperchrinfile[19];
					nbchrdivider[25][19]=3;
					if (nbchrdivider[25][19]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][19];
					chrdivider[25][20][0].end=2054;chrdivider[25][20][1].start=2054;
					chrdivider[25][20][1].end=5072;chrdivider[25][20][2].start=5072;
					chrdivider[25][20][2].end=nbsnpperchrinfile[20];
					nbchrdivider[25][20]=3;
					if (nbchrdivider[25][20]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][20];
					chrdivider[25][21][0].end=3025;chrdivider[25][21][1].start=3025;
					chrdivider[25][21][1].end=nbsnpperchrinfile[21];
					nbchrdivider[25][21]=2;
					if (nbchrdivider[25][21]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][21];
					chrdivider[25][22][0].end=1615;chrdivider[25][22][1].start=1615;
					chrdivider[25][22][1].end=nbsnpperchrinfile[22];
					nbchrdivider[25][22]=2;
					if (nbchrdivider[25][22]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][22];
					int nbdivisionDEGREE4=0; 
{	char segwithav[23][NSNPPERCHR];
						double seuilpihatcorrection=0.022;
						printf("%d %f\n",placefirttoconsider,pihatagainstall[placefirttoconsider]);
						for(int  relattocompare=0;relattocompare<20;relattocompare++)	
						{	int countsnp=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchrinfile[chrtemp1];chrdividerrun++)
								{	segwithav[chrtemp1][chrdividerrun]=0; 
								};
							};
							int IDrelattocompare=bestpihatagainstallID[placefirttoconsider+relattocompare];
							if (pihatagainstall[IDrelattocompare]>seuilpihatcorrection	)
							{	
								#pragma omp parallel for 	
								for(int  chrtemp1=22;chrtemp1>0;chrtemp1--)		
								{	printf("Start correcting chr %d\n",chrtemp1);
									int nbindivmatchseg[MAXPOP][4];
									for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++)
									{	nbindivmatchseg[relat][0]=0;
										nbindivmatchseg[relat][1]=0;
										nbindivmatchseg[relat][2]=0;
										nbindivmatchseg[relat][3]=0;
									};
									int nbsegmentofthislength[NSNPPERCHR]={0};
									int typesegment=-1;
									int lasttypesegment=-1;
									int endlastsegment[4]={-1};
									double ratetoconsiderseg=0.0001;
									double ratetoconsidersegPE=ratetoconsiderseg*2; 
									ratetoconsiderseg=0.0002;   
									ratetoconsidersegPE=2*ratetoconsiderseg; 
									int phaseerrorpossible[4]={0};
									int breaknubercm=25;
									for(int snp=0;snp<nbsnpperchrinfile[chrtemp1];snp++)
									{	
										int snpvalue0=genomeoffpss[0][snp][chrtemp1];		
										int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
										int snpmod4=(((snp%4)*2));
										int snpdiv4=snp/4;
										for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++) if (relat!=ID && pihatagainstall[relat]<seuilpihat[0] )
										{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
											if ((snpvalue0&1)==(snpvalue1&1))
											{	if (nbindivmatchseg[relat][0]<10) nbindivmatchseg[relat][0]++; else 
												nbsegmentofthislength[++nbindivmatchseg[relat][0]]++;
											} else 
											{	if (nbindivmatchseg[relat][0]<10) nbindivmatchseg[relat][0]=0; else 
												{	for(int lenght=10;lenght<nbindivmatchseg[relat][0]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][0]=0;
												}
											};
											if ((snpvalue0&1)==(snpvalue1>>1))
											{	if (nbindivmatchseg[relat][1]<10) nbindivmatchseg[relat][1]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][1]]++;
											} else 
											{	if (nbindivmatchseg[relat][1]<10) nbindivmatchseg[relat][1]=0; else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][1]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][1]=0;
												};
											};
											if ((snpvalue0>>1)==(snpvalue1&1))
											{	if (nbindivmatchseg[relat][2]<10) nbindivmatchseg[relat][2]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][2]]++;
											} else 
											{	if (nbindivmatchseg[relat][2]<10) nbindivmatchseg[relat][2]=0;
												else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][2]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][2]=0;
												}
											};
											if ((snpvalue0>>1)==(snpvalue1>>1))
											{	if (nbindivmatchseg[relat][3]<10) nbindivmatchseg[relat][3]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][3]]++;
											} else 
											{	if (nbindivmatchseg[relat][3]<10) nbindivmatchseg[relat][3]=0;	
												else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][3]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][3]=0;	
												};
											};
										};
										int relat=IDrelattocompare;
										if (nbindivmatchseg[relat][0]==0) phaseerrorpossible[0]=0;
										if (nbindivmatchseg[relat][1]==0) phaseerrorpossible[1]=0;
										if (nbindivmatchseg[relat][2]==0) phaseerrorpossible[2]=0;
										if (nbindivmatchseg[relat][3]==0) phaseerrorpossible[3]=0;
										int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
										int lenght;
										if (nbindivmatchseg[relat][0]>10 && 
												nbsegmentofthislength[nbindivmatchseg[relat][0]]<(phaseerrorpossible[0]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg) && 
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][1] &&
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][2] &&
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 1 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][0], nbsegmentofthislength[nbindivmatchseg[relat][0]]);
											lenght=nbindivmatchseg[relat][0];
											typesegment=1;
											for(int spnrun=snp-nbindivmatchseg[relat][0];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
											endlastsegment[0]=snp;	
										} else if (nbindivmatchseg[relat][1]>10 && nbsegmentofthislength[nbindivmatchseg[relat][1]]<(phaseerrorpossible[1]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg) &&
													nbindivmatchseg[relat][1]>nbindivmatchseg[relat][2] &&
													nbindivmatchseg[relat][1]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 3 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][1], nbsegmentofthislength[nbindivmatchseg[relat][1]]);
											lenght=nbindivmatchseg[relat][1];
											typesegment=3;
											endlastsegment[1]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][1];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
										} else if (nbindivmatchseg[relat][2]>10 && nbsegmentofthislength[nbindivmatchseg[relat][2]]<(phaseerrorpossible[2]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg) &&
													nbindivmatchseg[relat][2]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 2 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][2], nbsegmentofthislength[nbindivmatchseg[relat][2]]);
											lenght=nbindivmatchseg[relat][2];
											typesegment=2;
											endlastsegment[2]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][2];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=2;
										} else if (nbindivmatchseg[relat][3]>10 && nbsegmentofthislength[nbindivmatchseg[relat][3]]<(phaseerrorpossible[3]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg)
												)
										{	printf("chr %d segment type 4 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][3], nbsegmentofthislength[nbindivmatchseg[relat][3]]);
											lenght=nbindivmatchseg[relat][3];
											typesegment=4;
											endlastsegment[3]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][3];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=2;
										};
										if (endlastsegment[0]==snp-1 || endlastsegment[1]==snp-1 || endlastsegment[2]==snp-1 || endlastsegment[3]==snp-1)
										{	if ((snpvalue1&1)!=(snpvalue1>>1))
											{	phaseerrorpossible[0]=1;
												phaseerrorpossible[1]=1;
												phaseerrorpossible[2]=1;
												phaseerrorpossible[3]=1;
											}
										}
										if (typesegment>-1 && lasttypesegment>-1 && (typesegment&1)!=(lasttypesegment&1))
										{	int end;
											if (typesegment==1) end=endlastsegment[0];
											if (typesegment==2) end=endlastsegment[2];
											if (typesegment==3) end=endlastsegment[1];
											if (typesegment==4) end=endlastsegment[3];
											printf("change from %d to %d\n",(end+snp-lenght)/2,nbsnpperchrinfile[chrtemp1]);
											for(int snprun=(snp-lenght);snprun<nbsnpperchrinfile[chrtemp1];snprun++)		
											{
												genomeoffpss[0][snprun][chrtemp1]=((genomeoffpss[0][snprun][chrtemp1]&1)<<1)+(genomeoffpss[0][snprun][chrtemp1]>>1);
											};
											snp=end;
											for(int snprun=end;snprun<nbsnpperchrinfile[chrtemp1];snprun++)		
											{	segwithav[chrtemp1][snprun]=0;
											};
											typesegment=-1;
											for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++) if (pihatagainstall[relat]<seuilpihat[0] )
											{	nbindivmatchseg[relat][0]=0;
												nbindivmatchseg[relat][1]=0;
												nbindivmatchseg[relat][2]=0;
												nbindivmatchseg[relat][3]=0;
											};
											for(int snprun=0;snprun<nbsnpperchr[chrtemp1];snprun++)
											{	nbsegmentofthislength[snprun]=0;
											};
										} else if (typesegment>-1)
										{	lasttypesegment=typesegment;
										};									
									};		
								};
							};
						};
						int breaknubercm=25;
						int breaknubercmloop=25;
						{	for(int prod=0;prod<1;prod++)
							{	for(int chr=0;chr<23;chr++)
								{	for(int ratio=0;ratio<1;ratio++)
									{	phaseguessright[chr][prod][ratio]=0;
										phaseguesswrong[chr][prod][ratio]=0;
										phasenotguessed[chr][prod][ratio]=0;
										phase1[chr][prod][ratio]=0;
										guessphase1[chr][prod][ratio]=0;
										phase1guessphase1[chr][prod][ratio]=0;
										phase1guessphase2[chr][prod][ratio]=0;
										phase2[chr][prod][ratio]=0;
										guessphase2[chr][prod][ratio]=0;
										phase2guessphase2[chr][prod][ratio]=0;
										phase2guessphase1[chr][prod][ratio]=0;
									};	
								};
							};
							for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
							{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
								{	phaseparent1tap[chrtemp1][snp]=0;
									phaseguesedparent1tap[chrtemp1][snp]=0;
								};
							};							
							char segwithav[23][NSNPPERCHR];
							int countsnp=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchr[chrtemp1];chrdividerrun++)
								{	segwithav[chrtemp1][chrdividerrun]=0; 
								};
							};
							int64_t purcentage=0;
							ratio=2;
							printf("%d %d %d\n",6,phaseguessright[6][0][0],phaseguesswrong[6][0][0]);
							purcentage=(int64_t) 10000*
													(1+(phaseguessright[ratio][0][0]>phaseguesswrong[ratio][0][0]?
															phaseguessright[ratio][0][0]:
															phaseguesswrong[ratio][0][0]))/ 
													(phaseguesswrong[ratio][0][0]+phaseguessright[ratio][0][0]+1.0);
							if (purcentage>bestphasing) bestphasing=purcentage;
							int nbphaseerror=0;
							double pihatthrehjold[23];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	pihatthrehjold[chrtemp1]=0.10;
							};
							int nbphaseright=0;
							int nbphasewrong=0;
							int tabnbphaseright[23][25];
							int tabnbphasewrong[23][25];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	tabnbphaseright[chrtemp1][chrdividerrun]=0;
									tabnbphasewrong[chrtemp1][chrdividerrun]=0;
								};
							};
							double pihatagainstallchrMP[MAXPOP][2]={0};
							 #pragma omp parallel for
							for(int relat=0;relat<MAXPOP;relat++)
							{	pihatagainstallchrMP[relat][0]=0;
								pihatagainstallchrMP[relat][1]=0;
							}
							int sizebin=450;
							#pragma omp parallel for					
							for(int relat=0;relat<MAXPOP;relat++)
							{	pihatagainstallchrMPphaseerror[0][relat][0][0]=0;
								pihatagainstallchrMPphaseerror[0][relat][0][1]=0;
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][0]*4;chrdividerrun++)
								{	pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][0]=0;
									pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][1]=0;
								};
							}
							int hetsnp=0;
							double exponentinterchr=1;
							for(int relat=0;relat<MAXPOP;relat++)
							{	
							}
							for(int relat=0;relat<MAXPOP;relat++)
							{	if (pihatagainstall[relat]>0)  printf("PIHAT %d %f %f\n",relat,pihatagainstall[relat],seuilpihat[0]);
							}
							float firstexponent=2;
							int powersnppihat=3;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	
								{	hetsnp=0;
									#pragma omp parallel for 				
									for(int relat=0;relat<MAXPOP;relat++)
									{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1]*4;chrdividerrun++)
										{	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][0]=0;
											pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][1]=0;
										};
									};
									for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
									{	double temp[MAXPOP][4];
										for(int64_t relat=0;relat<(int64_t) MAXPOP;relat++) 
										{	temp[relat][0]=0;
											temp[relat][1]=0;
											temp[relat][2]=0;
											temp[relat][3]=0;
										};
										for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
										{	int snpvalue0=genomeoffpss[0][snp][chrtemp1];
											if (snpvalue0!=0 && snpvalue0!=3 )
											{
												hetsnp++;
												double maffloat=1.0*MAF[snp][chrtemp1]/MAXPOP/2;
												int parent0indiv0=(snpvalue0>>1);
												int parent1indiv0=(snpvalue0>>1);
												double pcontribu10=((parent0indiv0)-1.0*maffloat);
												double pcontribu11=((parent1indiv0)-1.0*maffloat);
												double maffloatdiviseur=(maffloat)*(2-maffloat*2)/3*2;
												double pconttab[4][2];
												pconttab[0][0]=pow(((0)-maffloat)/maffloatdiviseur*pcontribu10,1.0/powersnppihat);
												pconttab[0][1]=pow(((0)-maffloat)/maffloatdiviseur*pcontribu11,1.0/powersnppihat);
												pconttab[1][0]=pow(((1)-maffloat)/maffloatdiviseur*pcontribu10,1.0/powersnppihat);
												pconttab[1][1]=pow(((1)-maffloat)/maffloatdiviseur*pcontribu11,1.0/powersnppihat); 
												int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
												int snpmod4=(((snp%4)*2));
												int snpdiv4=snp/4;
												#pragma omp parallel for
												for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++) 
												{	if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )
													{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
														if (pconttab[snpvalue1&1][0]>0) temp[relat][0]=temp[relat][0]+(pconttab[snpvalue1&1][0]);
														else temp[relat][0]=0;
														if (temp[relat][0]<0) temp[relat][0]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]<temp[relat][0]) 
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]=temp[relat][0];
														}
														if (pconttab[snpvalue1>>1][0]>0) temp[relat][1]=temp[relat][1]+(pconttab[snpvalue1>>1][0]);
														else temp[relat][1]=0;
														if (temp[relat][1]<0) temp[relat][1]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]<temp[relat][1])
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]=temp[relat][1];
														};
														if (pconttab[snpvalue1&1][1]>0) temp[relat][2]=temp[relat][2]+(pconttab[snpvalue1&1][1]);
														else temp[relat][2]=0;
														if (temp[relat][2]<0) temp[relat][2]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]<temp[relat][2]) 
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]=temp[relat][2];
														}
														if (pconttab[snpvalue1>>1][1]>0) temp[relat][3]=temp[relat][3]+(pconttab[snpvalue1>>1][1]);
														else temp[relat][3]=0;
														if (temp[relat][3]<0) temp[relat][3]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]<temp[relat][3])
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]=temp[relat][3];
														};
													} else 
													{	printf("%d %f %f\n",relat,pihatagainstall[relat],pihatagainstall2[relat]);
													};
												};
											};
										};

										#pragma omp parallel for
										for(int relat=0;relat<NbIndiv;relat++) 
										{	if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 )
											{	int value=0;
												{	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]=
																					(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]>pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][value]?
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]:
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][value]);
													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]=
																					(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]>pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][value]?
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]:
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][value]);
												};
											};
										};
									
									};
								};	
							};
							float thirdindice=3;
							double corglobmax=0;						
							int chrmax1;						
							int chrmax2;						
							int dividmax1;						
							int dividmax2;						
							double sumall[23][25][2][2];
							double sumallsquare[23][25][2][2];
							#pragma omp parallel for 
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int value=0;value<2;value++)
									{	sumall[chrtemp1][chrdividerrun][0][value]=0;
										sumall[chrtemp1][chrdividerrun][1][value]=0;
										sumallsquare[chrtemp1][chrdividerrun][0][value]=0;
										sumallsquare[chrtemp1][chrdividerrun][1][value]=0;
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	sumall[chrtemp1][chrdividerrun][0][value]=sumall[chrtemp1][chrdividerrun][0][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]),firstexponent)*
													(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]>0?1:-1);
											sumall[chrtemp1][chrdividerrun][1][value]=sumall[chrtemp1][chrdividerrun][1][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]),firstexponent)*
													(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]>0?1:-1);
											sumallsquare[chrtemp1][chrdividerrun][0][value]=sumallsquare[chrtemp1][chrdividerrun][0][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]),firstexponent*2);
											sumallsquare[chrtemp1][chrdividerrun][1][value]=sumallsquare[chrtemp1][chrdividerrun][1][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]),firstexponent*2);
										};
									};
								};
							};
							double allcor[23][MAXNBDIVISOR][23][MAXNBDIVISOR];
							double allseg[23][MAXNBDIVISOR][23][MAXNBDIVISOR];
							#pragma omp parallel for		
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
									{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) 
										{	allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=0;
										};
									};
								};
							};			
							int countseg[23][MAXNBDIVISOR];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	countseg[chrtemp1][chrdividerrun]=1;
								};
							};
							float seuil1=0.5;
							#pragma omp parallel for		
							for(int relat=0;relat<NbIndiv;relat++) if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )
							{	
								#pragma omp parallel for		
								for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
								{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
									{	
										double thresholdseg=seuil1*(chrdivider[breaknubercm][chrtemp1][chrdividerrun].end-chrdivider[breaknubercm][chrtemp1][chrdividerrun].start);				
										if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]>thresholdseg || pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]>thresholdseg)
										{	if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]<thresholdseg || pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]<thresholdseg)
											{	countseg[chrtemp1][chrdividerrun]++;
											};
										};
									};
								};
							};
							double sumallseg=0;
							float penalty=0.7;
							double coefforcor=1.0/4;
							#pragma omp parallel for		
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
									{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2)
										{	double sumproduct; 
											sumproduct=0;
											double nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
											{	nbelem++;
												sumproduct=sumproduct+
															pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
															pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
											};
											double  cor1=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
											double sumproduct1=sumproduct;
											sumproduct=0;
											nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
											};
											double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
											sumproduct=0;
											 nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent)*
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
											};
											double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
											sumproduct=0;
											 nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 ) 
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent)*	
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
											};
											double  cor4=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
											double corglob=(pow(fabs(cor1),thirdindice-1)*cor1 -
																							pow(fabs(cor2),thirdindice-1)*cor2+ 
																							pow(fabs(cor3),thirdindice-1)*cor3-
																							pow(fabs(cor4),thirdindice-1)*cor4);
											if (corglob<0 && chrtemp1==chrtemp2) corglob=corglob*penalty;
											allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=corglob;
											allcor[chrtemp2][chrdividerrun2][chrtemp1][chrdividerrun]=corglob;
											if (fabs(corglob)>fabs(corglobmax))
											{	corglobmax=corglob;
												chrmax1=chrtemp1;						
												chrmax2=chrtemp2;						
												dividmax1=chrdividerrun;						
												dividmax2=chrdividerrun2;
												printf("%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f  \n",chrtemp1,chrtemp2,chrdividerrun,chrdividerrun2,cor1,cor2,cor3,cor4,
																							sumproduct,nbelem,
																							sumallsquare[chrtemp1][chrdividerrun][0][0],sumallsquare[chrtemp2][chrdividerrun2][0][0],
																							sumall[chrtemp1][chrdividerrun][0][0],sumall[chrtemp2][chrdividerrun2][0][0],
																							sumproduct1,allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2],corglobmax);
											};											
										};
									};
								};
							};
							printf("MAX %d %d %d %d %f\n",chrmax1,dividmax1,chrmax2,dividmax2,corglobmax);
							int group[23][MAXNBDIVISOR];
							int havemerged[23][MAXNBDIVISOR];
							int nbingroup[23][MAXNBDIVISOR];
							#pragma omp parallel for 	
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	havemerged[chrtemp1][chrdividerrun]=0;
									nbingroup[chrtemp1][chrdividerrun]=1;
									group[chrtemp1][chrdividerrun]=chrtemp1+chrdividerrun*23; 
								};
							};	
							int nbmerge=0;
							float mergingexponent[2];
							mergingexponent[0]=2;
							mergingexponent[1]=5; 
							double sumcorglobal=0;
							int summingroup=0;
							int summaxgroup=0;
							int sumgroup=0;
							int nbpair=0;
							do 
							{	group[chrmax2][dividmax2]=chrmax1+dividmax1*23;
								printf("MERGE chr %d chunk %d to chr %d chnk %d :%d %d %d %d corglobmax %f \n",chrmax2,dividmax2,chrmax1,dividmax1,tabnbphaseright[chrmax1][dividmax1],
											tabnbphasewrong[chrmax1][dividmax1],tabnbphaseright[chrmax2][dividmax2],tabnbphasewrong[chrmax2][dividmax2],corglobmax);
								tappointdec[nbmerge].nbgroup1=nbingroup[chrmax1][dividmax1];
								tappointdec[nbmerge].nbgroup2=nbingroup[chrmax2][dividmax2];
								tappointdec[nbmerge].cor=fabs(corglobmax);
								nbmerge++;
								sumcorglobal=sumcorglobal+fabs(corglobmax);
								sumgroup=sumgroup+nbingroup[chrmax1][dividmax1]+nbingroup[chrmax2][dividmax2];
								summingroup=summingroup+((nbingroup[chrmax1][dividmax1]<nbingroup[chrmax2][dividmax2])?nbingroup[chrmax1][dividmax1]:nbingroup[chrmax2][dividmax2]);
								summaxgroup=summaxgroup+((nbingroup[chrmax1][dividmax1]>nbingroup[chrmax2][dividmax2])?nbingroup[chrmax1][dividmax1]:nbingroup[chrmax2][dividmax2]);
								nbingroup[chrmax1][dividmax1]=nbingroup[chrmax1][dividmax1]+nbingroup[chrmax2][dividmax2];
								if (corglobmax>0)
								{	havemerged[chrmax2][dividmax2]=1;
									for(int  chrtemp3=1;chrtemp3<23;chrtemp3++)	
									{	for(int chrdividerrun3=0;chrdividerrun3<nbchrdivider[breaknubercm][chrtemp3];chrdividerrun3++) 
										{	
										};
									};
									if (nbdivisionDEGREE4==1) 
										countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==0) countseg[chrmax1][dividmax1]=(countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==2) countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2]+
																							((countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2]);
									#pragma omp parallel for		
									for(int relat=0;relat<NbIndiv;relat++) 
										if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 )										
									{	
										for(int value=0;value<2;value++)
										{	double p1=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
											double p2=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
											double p3=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]>0?1:-1);
											double p4=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]=pow(fabs(p1+p3),1.0/mergingexponent[value])*(p1+p3>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]=pow(fabs(p2+p4),1.0/mergingexponent[value])*(p2+p4>0?1:-1);
										};
									}	
								} else
								{	havemerged[chrmax2][dividmax2]=-1;
									#pragma omp parallel for
									for(int  chrtemp3=1;chrtemp3<23;chrtemp3++)		
									{	for(int chrdividerrun3=0;chrdividerrun3<nbchrdivider[breaknubercm][chrtemp3];chrdividerrun3++)
										{
										};
									};
									if (nbdivisionDEGREE4==2) 
										countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==0) countseg[chrmax1][dividmax1]=(countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==1) countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2]+
																							((countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2]);
									#pragma omp parallel for		
									for(int relat=0;relat<NbIndiv;relat++)
										if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
									{	for(int value=0;value<2;value++)
										{
											double p1=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
											double p2=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
											double p3=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]>0?1:-1);
											double p4=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]=pow(fabs(p1+p4),1.0/mergingexponent[value])*(p1+p4>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]=pow(fabs(p2+p3),1.0/mergingexponent[value])*(p2+p3>0?1:-1);
										};
									};					
								};
								sumall[chrmax1][dividmax1][0][0]=0;
								sumall[chrmax1][dividmax1][1][0]=0;
								sumallsquare[chrmax1][dividmax1][0][0]=0;
								sumallsquare[chrmax1][dividmax1][1][0]=0;
								sumall[chrmax1][dividmax1][0][1]=0;
								sumall[chrmax1][dividmax1][1][1]=0;
								sumallsquare[chrmax1][dividmax1][0][1]=0;
								sumallsquare[chrmax1][dividmax1][1][1]=0;
								for(int relat=0;relat<NbIndiv;relat++) 
									if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 ) 
								{	for(int value=0;value<2;value++)
									{	sumall[chrmax1][dividmax1][0][value]=sumall[chrmax1][dividmax1][0][value]+
											pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*(0.3+0.7))*
												(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
										sumall[chrmax1][dividmax1][1][value]=sumall[chrmax1][dividmax1][1][value]+
											pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*(0.3+0.7))*
											(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
										sumallsquare[chrmax1][dividmax1][0][value]=sumallsquare[chrmax1][dividmax1][0][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*2);
										sumallsquare[chrmax1][dividmax1][1][value]=sumallsquare[chrmax1][dividmax1][1][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*2);
									};
								};
								int chrtemp1=chrmax1;
								int chrdividerrun=dividmax1;
								corglobmax=0;	
								double highestnew=0;		
								for(int  chrtemp2=1;chrtemp2<23;chrtemp2++)		
								{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (havemerged[chrtemp2][chrdividerrun2]==0 && (chrtemp1!=chrtemp2 || chrdividerrun!=chrdividerrun2))
									{	double sumproduct; 
										sumproduct=0;
										int nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )											
										{	nbelem++;
											sumproduct=sumproduct+
														pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
														pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
										};
										double  cor1=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
										sumproduct=0;
										sumproduct=0;
										nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
										};
										double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
										sumproduct=0;
										sumproduct=0;
										 nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
										};
										double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
										sumproduct=0;
										sumproduct=0;
										 nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent*(0.3+0.7));
										};
										double  cor4=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
										sumproduct=0;
										float exponantnumer=1;
										double corglob=(pow(fabs(cor1),thirdindice-1)*cor1 -
																						pow(fabs(cor2),thirdindice-1)*cor2+ 
																						pow(fabs(cor3),thirdindice-1)*cor3-
																						pow(fabs(cor4),thirdindice-1)*cor4)
																				*pow((nbingroup[chrtemp1][chrdividerrun]+nbingroup[chrtemp2][chrdividerrun2]),exponantnumer);
										if (highestnew<fabs(corglob)) 
										{	highestnew=fabs(corglob); 
											printf("%f",corglob);
											printf("NB1 %d MB2 %d \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2]);  
										};
										allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=corglob;										
										allcor[chrtemp2][chrdividerrun2][chrtemp1][chrdividerrun]=corglob;										
									};
								};
								nbpair=0;
								for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
								{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++) if (havemerged[chrtemp1][chrdividerrun]==0)
									{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
										{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (havemerged[chrtemp2][chrdividerrun2]==0 && (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2))
											{	nbpair++;
												double corglob=allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2];
												if (fabs(corglob)>fabs(corglobmax))
												{
													corglobmax=corglob;
													chrmax1=chrtemp1;						
													chrmax2=chrtemp2;						
													dividmax1=chrdividerrun;						
													dividmax2=chrdividerrun2;
													printf("NB1 %d MB2 %d c %f \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2],corglob);  
													if (nbingroup[chrtemp1][chrdividerrun]==1 || nbingroup[chrtemp2][chrdividerrun2]==1)
													{	int chrone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrtemp2:chrtemp1;
														int chrdividerrunone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrdividerrun2:chrdividerrun;
														int chrnotone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrtemp2:chrtemp1;
														int chrdividerrunnotone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrdividerrun2:chrdividerrun;
														printf("1 is %d %d\n",chrnotone,chrdividerrunnotone);
														int numberattached=0;
														for(int chrdividerrunrun=0;chrdividerrunrun<nbchrdivider[breaknubercm][chrone];chrdividerrunrun++) 
														{
															if (chrdividerrunrun==chrdividerrunone)
															{	
															} else
															{	if (havemerged[chrone][chrdividerrunrun]==0)
																{
																} else 
																{	int hrdividerloop=chrdividerrunrun;	
																	int chrloop=chrone;
																	int phaseknown=havemerged[chrloop][hrdividerloop];
																	while (iter<1000 && havemerged[chrloop][hrdividerloop]!=0)
																	{	
																		iter++;
																		phaseknown=phaseknown*(havemerged[chrloop][hrdividerloop]!=0?havemerged[chrloop][hrdividerloop]:1);
																		int savechr=chrloop;
																		chrloop=group[chrloop][hrdividerloop]%23;
																		hrdividerloop=group[savechr][hrdividerloop]/23;
																	};
																	if 	(chrloop==chrnotone && hrdividerloop==chrdividerrunnotone)
																	{	printf("window n%d is attached to %d %d with phaseknown %d\n", chrdividerrunrun,chrnotone,chrdividerrunnotone,phaseknown);
																		numberattached=numberattached+phaseknown;
																	} else 
																	{	printf("window n%d is not attached to %d %d\n", chrdividerrunrun,chrnotone,chrdividerrunnotone);
																	};																		
																};
															};
														};
														if (corglobmax*numberattached<0) corglobmax=corglobmax*penalty;
													};
												};	
											};
										};
									};
								};
								printf("Nb pairs: %d\n",nbpair);
							} while (nbpair>0 && nbmerge<22*19-1);
							 nbphaseright=0;
							 nbphasewrong=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		 
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	int chr=chrtemp1;
									int chrdivi=chrdividerrun;
									int phaseknown=1;
									int iter=0;
									while (iter<1000 && havemerged[chr][chrdivi]!=0)
									{	
										iter++;
										phaseknown=phaseknown*(havemerged[chr][chrdivi]!=0?havemerged[chr][chrdivi]:1);
										int savechr=chr;
										chr=group[chr][chrdivi]%23;
										chrdivi=group[savechr][chrdivi]/23;
									} 
									if (iter>900) exit(0);
									printf("pahse chr %d div %d is %d\n",chrtemp1,chrdividerrun,phaseknown);
									chrdivider[breaknubercm][chrtemp1][chrdividerrun].phasing=phaseknown;
									for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
									{	if ((phaseknown+3)/2==2) 
										{	if (genomeoffpss[0][snp][chrtemp1]==1) 
											{	genomeoffpss[0][snp][chrtemp1]=2;
												*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
													(*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
														(~(3<<((snp%4)*2))))
															|	((2)<<((snp%4)*2));
											}
											else if (genomeoffpss[0][snp][chrtemp1]==2) 
											{	genomeoffpss[0][snp][chrtemp1]=1;
												*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
													(*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
														(~(3<<((snp%4)*2))))
															|	((1)<<((snp%4)*2));
											};
										};
										if ((phaseknown+3)/2==phaseparent1tap[chrtemp1][snp]) 
										{	nbphaseright++;
											chrdivider[breaknubercm][0][0].nbright++;
											chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbright++;
										} else 
										{	nbphasewrong++;
											chrdivider[breaknubercm][0][0].nbwrong++;
											chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbwrong++;
										};
									};
									printf("chr %d div %d nbphaseright %d nbphasewrong %d\n",chrtemp1,chrdividerrun,nbphaseright,nbphasewrong);
								};
								printf("chr %d nbphaseright %d nbphasewrong %d ratio %f \n",chrtemp1,nbphaseright,nbphasewrong,
																							(nbphaseright>nbphasewrong?nbphaseright:nbphasewrong)/(1.0+nbphaseright+nbphasewrong));
							};
							purcentage=(int64_t) 10000*
											(nbphaseright>nbphasewrong?nbphaseright:nbphasewrong)/(nbphaseright+nbphasewrong>0?nbphaseright+nbphasewrong:1);						
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)	
							{	dataperchr[chrtemp1][52+nbdivisionDEGREE4-18+3*powerpihatDEGREE4+breaknubercmloop-23]=purcentage;		
							};
						}
					};
				};					
			};	
		};
	} while (gen<MAXGEN);	
	return (0);
}
	int predict_phasing_with_parents_providing_GT(int ID,int chr1, int chr2,int numtrio,float limitnbsnp,int run,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[])
{	printf("%d\n",IDjob);
	seuilpihat[0]=0.33;
	seuilpihat[1]=0.33;
	seuilpihat[2]=0.33;
	int increaseincrementwhenneg=1;
	int takesumintoaccout=1;
	int limitweight=105;
	int howtocalculatenotimprove=0;
	int howtocalculatenotimprove2=0;
	int nbmutdependonnotimprove=0;
	int mutateallhapwhennoimprove=0;
	int changeonlyifneg=2;
	int crossclever=0;
	int numgentowait=17;
	int howtocalculaterelattochr=0; 
	int nbmutvari=1;
	int64_t tabresmut[50][50];
	int nbtimetabresmut[50][50];
	for(int mut1=0;mut1<50;mut1++)
	{	for(int mut2=0;mut2<50;mut2++)
		{	tabresmut[mut1][mut2]=0;
			nbtimetabresmut[mut1][mut2]=1;
		};
	};
	int mostrelatinscore=0;
	int howcalculmostrelaetd=0;
	int howcalculmostrelaetd2=0;
	int addoneandzero=1;
	int localsearch=1;
	int randmuthap=0;
	unsigned char relatsuperpose[MAXCLOSERELAT][MAXCLOSERELAT];
	int64_t relatpihatchr[MAXCLOSERELAT][23][2];
	int seuilscore=1;
	int powersegment=1;		
	int64_t segnum=0;
	unsigned char relatbreak[MAXBREAK][MAXCLOSERELAT];
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	for(int  break1=0;break1<nbbreak;break1++)		
		{	relatbreak[break1][relat]=0;
		};
	};
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	int highestnbhet=0;
		int highestnbhetbreak=0;
		for(int  break1=0;break1<nbbreak;break1++)		
		{	if (highestnbhet<relatbreak[break1][relat]) 
			{	highestnbhet=relatbreak[break1][relat];
				highestnbhetbreak=break1;
			};
			if (relatbreak[break1][relat]) relatbreak[break1][relat]=1;
		};
		for(int  break1=0;break1<nbbreak;break1++)		
		{	relatbreak[break1][relat]=(highestnbhetbreak==break1);
		};
	};	
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	printf("%d ",relat);
		for(int  break1=0;break1<nbbreak;break1++)		
		{	printf("%d ",relatbreak[break1][relat]);
		};
	};
	int relattochr[MAXCLOSERELAT];
	for(int relat=0;relat<nbrelatpihat;relat++)
	{	int bestchr=1;
		for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
		{	int hapmost=(relatpihatchr[relat][chrtemp1][0]>relatpihatchr[relat][chrtemp1][1])?0:1;
			int pihatnorm=relatpihatchr[relat][chrtemp1][hapmost]*nbsnpperchrinfile[1]/nbsnpperchrinfile[chrtemp1];
			if (chrtemp1==1 || relattochr[relat]<pihatnorm)  
			{	relattochr[relat]=pihatnorm;
				bestchr=chrtemp1+22*hapmost;
			};
		};
		relattochr[relat]=bestchr;
	};
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	for(int snp=0;snp<nbsnpperchrinfile[chrtemp1];snp++)		
		{	genomeoffpss[0][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			genomeoffpss[1][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			genomeoffpss[2][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) IDp2loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			printf("%d %d %d\n",chrtemp1,snp,genomeoffpss[0][snp][chrtemp1]);
		};	
	};
	int maxpihat=0;
	int64_t nbseghap[2][23][26229];
	int64_t nbsegoverchr[23][2];
	int ratio=2;	
	printf("maxpihat %d\n",maxpihat);
	FILE * endfile;
	uint64_t averageofaverage[23];
	char number[100];
	char pathfile[300]; 								
	uint64_t sumsumnincrementplus=0;
	uint64_t sumsumnincrementmoins=0;
	sumsumnincrementplus=0;
	sumsumnincrementmoins=0;
	float plusmoins=1.0*sumsumnincrementplus/sumsumnincrementmoins;
	float moinsplus=1.0*sumsumnincrementmoins/sumsumnincrementplus;
	printf("%f %f\n",plusmoins,moinsplus);
	sumsumnincrementplus=0;
	sumsumnincrementmoins=0;
	int64_t sumallpihatpos=0;
	int64_t sumallpihatneg=0;
	int nbpos=0;
	int nbneg=0;
	int64_t minpihat=0;
	int minpihatID1=0;
	int minpihatID2=0;
	unsigned char indivEA[NBINDIVEA][MAXCLOSERELAT];
	unsigned char bestindiv[MAXCLOSERELAT];
	int64_t scoreEA[NBINDIVEA];	
	int64_t bestscore=-1000000;
	int64_t nbsumerrorbestscore=0;
	int64_t bestscoreofbestfitness=0;
	int64_t bestphasing=0;
	int same;
	int gen=gentostart;	
	int lastgen=100;	
	int nuincreasenewtop=50;
	int64_t averageparentold5=-1000000;
	double averageparent[MAXGEN];
	int print=1;
	int lastgenimprove=gen;
	int genintiindiv=gen;
	int64_t bestscoresinceinit=bestscore;
	int64_t bestscorefitness=-100000;
	int64_t bestscoreclustering=-100000;
	do
	{	int64_t scorerelattab[NBINDIVEA][MAXCLOSERELAT];
		int64_t scorerelattaballchr[NBINDIVEA][MAXCLOSERELAT][MAXBREAK];
		double scorehap[NBINDIVEA][MAXBREAK];
		{	int groumdtruthdone=0;
			int bestwellphased=0;
			int bestunwellphased=0;
			int64_t score1=0;
			int64_t score2=0;
			int64_t score1tab[MAXBREAK];
			int64_t score2tab[MAXBREAK];			
			for(int  chrtemp1=0;chrtemp1<MAXBREAK;chrtemp1++)		
			{	score1tab[chrtemp1]=0;
				score2tab[chrtemp1]=0;
			}
			int nbonebalanced=0;
		
			int64_t score1averageneg=0;
			int score1nbneg=1;
			int64_t score2averageneg=0;
			int score2nbneg=1;
			int chagemade=1;
			int64_t iter=0;
			int iterbreak=0;
			int64_t maxscorerelattab=0;
			int64_t plus[MAXCLOSERELAT]={0};
			int64_t minus[MAXCLOSERELAT]={0};
			maxscorerelattab=0;
			int64_t maxscorerelattabID=0;
			int64_t score1temp=0;
			 int64_t  score2temp=0;
			printf("%d %d %d %d\n",score1,score2,score1nbneg,score1nbneg);
			int64_t score;
			int64_t averagescore1=1;
			int64_t averagescore2=1;
			if (nbonebalanced*2>nbrelatpihat) nbonebalanced=1+nbonebalanced*2-nbrelatpihat; else nbonebalanced=1+nbrelatpihat-nbonebalanced*2;
			if (bestscoresinceinit<score) 
			{
				lastgenimprove=gen;
				bestscoresinceinit=score;
			};
			int64_t highest=score;
			{	
				bestscoreclustering=score;
				if (bestscoreofbestfitness==0) bestscoreofbestfitness=1;
				char number[100];
				char pathfile[300]; 
				int64_t scorechr1[23];
				int64_t scorechr2[23];
				int64_t scorechr3[23];
				int64_t scorechr4[23];
				int keep[23][4];
				int64_t nbonefromhap[2][23][26229];
				int64_t nbzerofromhap[2][23][26229];
				int64_t nbonequalified[2];
				int64_t nbzeroqualified[2];
				nbonequalified[0]=0;
				nbonequalified[1]=0;
				nbzeroqualified[0]=0;
				nbzeroqualified[1]=0;
				int64_t sumdiff[23];
				int64_t nbcontracditperchr[23];
				for(int nbkickout=1;nbkickout<23;nbkickout++)
				{	keep[nbkickout][ratio]=1;
				};
				{	nbsumerrorbestscore=1;
					int phase1[23][1][1];
					int guessphase1[23][1][1];
					int phase1guessphase1[23][1][1];
					int phase1guessphase2[23][1][1];
					int phase2[23][1][1];
					int guessphase2[23][1][1];
					int phase2guessphase2[23][1][1];
					int phase2guessphase1[23][1][1];
					int phaseguessright[23][1][1];
					int phaseguesswrong[23][1][1];
					int phasenotguessed[23][1][1];
					int nbsum=1;
					int wrong[23];
					int right[23];
					unsigned int x, y, width, height;
					for(int prod=0;prod<1;prod++)
					{	for(int chr=0;chr<23;chr++)
						{	for(int ratio=0;ratio<1;ratio++)
							{	phaseguessright[chr][prod][ratio]=0;
								phaseguesswrong[chr][prod][ratio]=0;
								phasenotguessed[chr][prod][ratio]=0;
								phase1[chr][prod][ratio]=0;
								guessphase1[chr][prod][ratio]=0;
								phase1guessphase1[chr][prod][ratio]=0;
								phase1guessphase2[chr][prod][ratio]=0;
								phase2[chr][prod][ratio]=0;
								guessphase2[chr][prod][ratio]=0;
								phase2guessphase2[chr][prod][ratio]=0;
								phase2guessphase1[chr][prod][ratio]=0;
							};	
						};
					};
					int lenshowchr=height*0.9;
					int widhtshowchr=5;
					int gapshowchry=63;
					int widthshowchry=9;
					int gapshowchrx=5;
					int64_t sumdifftot=0;
					int phaseparent1=-1;
					int phaseparent1tap[23][26229];
					int phaseguesedparent1tap[23][26229];
					for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
					{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
						{	phaseparent1tap[chrtemp1][snp]=0;
							phaseguesedparent1tap[chrtemp1][snp]=0;
						};
					};
					int phaseunknow=0;	
					int bestpur=0;
					double dataperchr[23][110];		
					int nbdivision=2;
					int powerpihat=0;
					int exponentcorp=0;
					double sumoftabcor=0;
					int firstexponent;
					int thirdindice;
					double pihatagainstallchrMPphaseerror[23][MAXPOP][200][2];
					int purcentage;
					int powerpihatDEGREE4=0;
					for(int size=0;size<51;size++)
					{	for(int chr=1;chr<23;chr++)
						{	for(int window=0;window<1;window++)
							{	chrdivider[size][chr][window].segment=0;
							};
						};
					};
					chrdivider[25][1][0].end=2203;chrdivider[25][1][1].start=2203;
					chrdivider[25][1][1].end=7136;chrdivider[25][1][2].start=7136;
					chrdivider[25][1][2].end=14606;chrdivider[25][1][3].start=14606;
					chrdivider[25][1][3].end=18289;chrdivider[25][1][4].start=18289;
					chrdivider[25][1][4].end=nbsnpperchrinfile[1];
					nbchrdivider[25][1]=5;
					nbchrdivider[25][0]=nbchrdivider[25][1];
					chrdivider[25][2][0].end=2097;chrdivider[25][2][1].start=2097;
					chrdivider[25][2][1].end=6555;chrdivider[25][2][2].start=6555;
					chrdivider[25][2][2].end=12181;chrdivider[25][2][3].start=12181;
					chrdivider[25][2][3].end=18340;chrdivider[25][2][4].start=18340;
					chrdivider[25][2][4].end=nbsnpperchrinfile[2];
					nbchrdivider[25][2]=5;
					if (nbchrdivider[25][2]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][2];
					chrdivider[25][3][0].end=4737;chrdivider[25][3][1].start=4737;
					chrdivider[25][3][1].end=9553;chrdivider[25][3][2].start=9553;
					chrdivider[25][3][2].end=14838;chrdivider[25][3][3].start=14838;
					chrdivider[25][3][3].end=18519;chrdivider[25][3][4].start=18519;
					chrdivider[25][3][4].end=nbsnpperchrinfile[3];
					nbchrdivider[25][3]=5;
					if (nbchrdivider[25][3]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][3];
					chrdivider[25][4][0].end=3754;chrdivider[25][4][1].start=3754;
					chrdivider[25][4][1].end=6615;chrdivider[25][4][2].start=6615;
					chrdivider[25][4][2].end=9749;chrdivider[25][4][3].start=9749;
					chrdivider[25][4][3].end=16124;chrdivider[25][4][4].start=16124;
					chrdivider[25][4][4].end=nbsnpperchrinfile[4];
					nbchrdivider[25][4]=5;
					if (nbchrdivider[25][4]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][4];
					chrdivider[25][5][0].end=6946;chrdivider[25][5][1].start=6946;
					chrdivider[25][5][1].end=11357;chrdivider[25][5][2].start=11357;
					chrdivider[25][5][2].end=15071;chrdivider[25][5][3].start=15071;
					chrdivider[25][5][3].end=nbsnpperchrinfile[5];
					nbchrdivider[25][5]=4;
					if (nbchrdivider[25][5]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][5];
					chrdivider[25][6][0].end=3976;chrdivider[25][6][1].start=3976;
					chrdivider[25][6][1].end=11907;chrdivider[25][6][2].start=11907;
					chrdivider[25][6][2].end=nbsnpperchrinfile[6];
					nbchrdivider[25][6]=3;
					if (nbchrdivider[25][6]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][6];
					chrdivider[25][7][0].end=3082;chrdivider[25][7][1].start=3082;
					chrdivider[25][7][1].end=7188;chrdivider[25][7][2].start=7188;
					chrdivider[25][7][2].end=12684;chrdivider[25][7][3].start=12684;
					chrdivider[25][7][3].end=nbsnpperchrinfile[7];
					nbchrdivider[25][7]=4;
					if (nbchrdivider[25][7]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][7];
					chrdivider[25][8][0].end=7271;chrdivider[25][8][1].start=7271;
					chrdivider[25][8][1].end=12465;chrdivider[25][8][2].start=12465;
					chrdivider[25][8][2].end=nbsnpperchrinfile[8];
					nbchrdivider[25][8]=3;
					if (nbchrdivider[25][8]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][8];
					chrdivider[25][9][0].end=3384;chrdivider[25][9][1].start=3384;
					chrdivider[25][9][1].end=7020;chrdivider[25][9][2].start=7020;
					chrdivider[25][9][2].end=9803;chrdivider[25][9][3].start=9803;
					chrdivider[25][9][3].end=nbsnpperchrinfile[9];
					nbchrdivider[25][9]=4;
					if (nbchrdivider[25][9]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][9];
					chrdivider[25][10][0].end=2407;chrdivider[25][10][1].start=2407;
					chrdivider[25][10][1].end=8424;chrdivider[25][10][2].start=8424;
					chrdivider[25][10][2].end=12288;chrdivider[25][10][3].start=12288;
					chrdivider[25][10][3].end=nbsnpperchrinfile[10];
					nbchrdivider[25][10]=4;
					if (nbchrdivider[25][10]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][10];
					chrdivider[25][11][0].end=2953;chrdivider[25][11][1].start=2953;
					chrdivider[25][11][1].end=5344;chrdivider[25][11][2].start=5344;
					chrdivider[25][11][2].end=nbsnpperchrinfile[11];
					nbchrdivider[25][11]=3;
					if (nbchrdivider[25][11]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][11];
					chrdivider[25][12][0].end=5937;chrdivider[25][12][1].start=5937;
					chrdivider[25][12][1].end=9107;chrdivider[25][12][2].start=9107;
					chrdivider[25][12][2].end=14609;chrdivider[25][12][3].start=14609;
					chrdivider[25][12][3].end=nbsnpperchrinfile[12];
					nbchrdivider[25][12]=4;
					if (nbchrdivider[25][12]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][12];
					chrdivider[25][13][0].end=2524;chrdivider[25][13][1].start=2524;
					chrdivider[25][13][1].end=6229;chrdivider[25][13][2].start=6229;
					chrdivider[25][13][2].end=nbsnpperchrinfile[13];
					nbchrdivider[25][13]=3;
					if (nbchrdivider[25][13]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][13];
					chrdivider[25][14][0].end=3228;chrdivider[25][14][1].start=3228;
					chrdivider[25][14][1].end=6437;chrdivider[25][14][2].start=6437;
					chrdivider[25][14][2].end=nbsnpperchrinfile[14];
					nbchrdivider[25][14]=3;
					if (nbchrdivider[25][14]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][14];
					chrdivider[25][15][0].end=3288;chrdivider[25][15][1].start=3288;
					chrdivider[25][15][1].end=6043;chrdivider[25][15][2].start=6043;
					chrdivider[25][15][2].end=8861;chrdivider[25][15][3].start=8861;
					chrdivider[25][15][3].end=nbsnpperchrinfile[15];
					nbchrdivider[25][15]=4;
					if (nbchrdivider[25][15]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][15];
					chrdivider[25][16][0].end=3860;chrdivider[25][16][1].start=3860;
					chrdivider[25][16][1].end=nbsnpperchrinfile[16];
					nbchrdivider[25][16]=2;
					if (nbchrdivider[25][16]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][16];
					chrdivider[25][17][0].end=2196;chrdivider[25][17][1].start=2196;
					chrdivider[25][17][1].end=4676;chrdivider[25][17][2].start=4676;
					chrdivider[25][17][2].end=8643;chrdivider[25][17][3].start=8643;
					chrdivider[25][17][3].end=nbsnpperchrinfile[17];
					nbchrdivider[25][17]=4;
					if (nbchrdivider[25][17]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][17];
					chrdivider[25][18][0].end=3734;chrdivider[25][18][1].start=3734;
					chrdivider[25][18][1].end=7591;chrdivider[25][18][2].start=7591;
					chrdivider[25][18][2].end=nbsnpperchrinfile[18];
					nbchrdivider[25][18]=3;
					if (nbchrdivider[25][18]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][18];
					chrdivider[25][19][0].end=2982;chrdivider[25][19][1].start=2982;
					chrdivider[25][19][1].end=6035;chrdivider[25][19][2].start=6035;
					chrdivider[25][19][2].end=nbsnpperchrinfile[19];
					nbchrdivider[25][19]=3;
					if (nbchrdivider[25][19]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][19];
					chrdivider[25][20][0].end=2054;chrdivider[25][20][1].start=2054;
					chrdivider[25][20][1].end=5072;chrdivider[25][20][2].start=5072;
					chrdivider[25][20][2].end=nbsnpperchrinfile[20];
					nbchrdivider[25][20]=3;
					if (nbchrdivider[25][20]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][20];
					chrdivider[25][21][0].end=3025;chrdivider[25][21][1].start=3025;
					chrdivider[25][21][1].end=nbsnpperchrinfile[21];
					nbchrdivider[25][21]=2;
					if (nbchrdivider[25][21]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][21];
					chrdivider[25][22][0].end=1615;chrdivider[25][22][1].start=1615;
					chrdivider[25][22][1].end=nbsnpperchrinfile[22];
					nbchrdivider[25][22]=2;
					if (nbchrdivider[25][22]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][22];
					int nbdivisionDEGREE4=0; 
{	char segwithav[23][NSNPPERCHR];
						double seuilpihatcorrection=0.022;
						printf("%d %f\n",placefirttoconsider,pihatagainstall[placefirttoconsider]);
						for(int  relattocompare=0;relattocompare<20;relattocompare++)	
						{	int countsnp=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchrinfile[chrtemp1];chrdividerrun++)
								{	segwithav[chrtemp1][chrdividerrun]=0; 
								};
							};
							int IDrelattocompare=bestpihatagainstallID[placefirttoconsider+relattocompare];
							if (pihatagainstall[IDrelattocompare]>seuilpihatcorrection	)
							{	
								#pragma omp parallel for 	
								for(int  chrtemp1=22;chrtemp1>0;chrtemp1--)		
								{	printf("Start correcting chr %d\n",chrtemp1);
									int nbindivmatchseg[MAXPOP][4];
									for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++)
									{	nbindivmatchseg[relat][0]=0;
										nbindivmatchseg[relat][1]=0;
										nbindivmatchseg[relat][2]=0;
										nbindivmatchseg[relat][3]=0;
									};
									int nbsegmentofthislength[NSNPPERCHR]={0};
									int typesegment=-1;
									int lasttypesegment=-1;
									int endlastsegment[4]={-1};
									double ratetoconsiderseg=0.0001;
									double ratetoconsidersegPE=ratetoconsiderseg*2; 
									ratetoconsiderseg=0.0002;   
									ratetoconsidersegPE=2*ratetoconsiderseg; 
									int phaseerrorpossible[4]={0};
									int breaknubercm=25;
									for(int snp=0;snp<nbsnpperchrinfile[chrtemp1];snp++)
									{	
										int snpvalue0=genomeoffpss[0][snp][chrtemp1];		
										int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
										int snpmod4=(((snp%4)*2));
										int snpdiv4=snp/4;
										for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++) if (relat!=ID && pihatagainstall[relat]<seuilpihat[0] )
										{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
											if ((snpvalue0&1)==(snpvalue1&1))
											{	if (nbindivmatchseg[relat][0]<10) nbindivmatchseg[relat][0]++; else 
												nbsegmentofthislength[++nbindivmatchseg[relat][0]]++;
											} else 
											{	if (nbindivmatchseg[relat][0]<10) nbindivmatchseg[relat][0]=0; else 
												{	for(int lenght=10;lenght<nbindivmatchseg[relat][0]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][0]=0;
												}
											};
											if ((snpvalue0&1)==(snpvalue1>>1))
											{	if (nbindivmatchseg[relat][1]<10) nbindivmatchseg[relat][1]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][1]]++;
											} else 
											{	if (nbindivmatchseg[relat][1]<10) nbindivmatchseg[relat][1]=0; else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][1]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][1]=0;
												};
											};
											if ((snpvalue0>>1)==(snpvalue1&1))
											{	if (nbindivmatchseg[relat][2]<10) nbindivmatchseg[relat][2]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][2]]++;
											} else 
											{	if (nbindivmatchseg[relat][2]<10) nbindivmatchseg[relat][2]=0;
												else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][2]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][2]=0;
												}
											};
											if ((snpvalue0>>1)==(snpvalue1>>1))
											{	if (nbindivmatchseg[relat][3]<10) nbindivmatchseg[relat][3]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][3]]++;
											} else 
											{	if (nbindivmatchseg[relat][3]<10) nbindivmatchseg[relat][3]=0;	
												else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][3]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][3]=0;	
												};
											};
										};
										int relat=IDrelattocompare;
										if (nbindivmatchseg[relat][0]==0) phaseerrorpossible[0]=0;
										if (nbindivmatchseg[relat][1]==0) phaseerrorpossible[1]=0;
										if (nbindivmatchseg[relat][2]==0) phaseerrorpossible[2]=0;
										if (nbindivmatchseg[relat][3]==0) phaseerrorpossible[3]=0;
										int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
										int lenght;
										if (nbindivmatchseg[relat][0]>10 && 
												nbsegmentofthislength[nbindivmatchseg[relat][0]]<(phaseerrorpossible[0]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg) && 
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][1] &&
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][2] &&
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 1 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][0], nbsegmentofthislength[nbindivmatchseg[relat][0]]);
											lenght=nbindivmatchseg[relat][0];
											typesegment=1;
											for(int spnrun=snp-nbindivmatchseg[relat][0];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
											endlastsegment[0]=snp;	
										} else if (nbindivmatchseg[relat][1]>10 && nbsegmentofthislength[nbindivmatchseg[relat][1]]<(phaseerrorpossible[1]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg) &&
													nbindivmatchseg[relat][1]>nbindivmatchseg[relat][2] &&
													nbindivmatchseg[relat][1]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 3 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][1], nbsegmentofthislength[nbindivmatchseg[relat][1]]);
											lenght=nbindivmatchseg[relat][1];
											typesegment=3;
											endlastsegment[1]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][1];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
										} else if (nbindivmatchseg[relat][2]>10 && nbsegmentofthislength[nbindivmatchseg[relat][2]]<(phaseerrorpossible[2]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg) &&
													nbindivmatchseg[relat][2]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 2 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][2], nbsegmentofthislength[nbindivmatchseg[relat][2]]);
											lenght=nbindivmatchseg[relat][2];
											typesegment=2;
											endlastsegment[2]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][2];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=2;
										} else if (nbindivmatchseg[relat][3]>10 && nbsegmentofthislength[nbindivmatchseg[relat][3]]<(phaseerrorpossible[3]?NbIndiv*ratetoconsidersegPE:NbIndiv*ratetoconsiderseg)
												)
										{	printf("chr %d segment type 4 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][3], nbsegmentofthislength[nbindivmatchseg[relat][3]]);
											lenght=nbindivmatchseg[relat][3];
											typesegment=4;
											endlastsegment[3]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][3];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=2;
										};
										if (endlastsegment[0]==snp-1 || endlastsegment[1]==snp-1 || endlastsegment[2]==snp-1 || endlastsegment[3]==snp-1)
										{	if ((snpvalue1&1)!=(snpvalue1>>1))
											{	phaseerrorpossible[0]=1;
												phaseerrorpossible[1]=1;
												phaseerrorpossible[2]=1;
												phaseerrorpossible[3]=1;
											}
										}
										if (typesegment>-1 && lasttypesegment>-1 && (typesegment&1)!=(lasttypesegment&1))
										{	int end;
											if (typesegment==1) end=endlastsegment[0];
											if (typesegment==2) end=endlastsegment[2];
											if (typesegment==3) end=endlastsegment[1];
											if (typesegment==4) end=endlastsegment[3];
											printf("change from %d to %d\n",(end+snp-lenght)/2,nbsnpperchrinfile[chrtemp1]);
											for(int snprun=(snp-lenght);snprun<nbsnpperchrinfile[chrtemp1];snprun++)		
											{
												genomeoffpss[0][snprun][chrtemp1]=((genomeoffpss[0][snprun][chrtemp1]&1)<<1)+(genomeoffpss[0][snprun][chrtemp1]>>1);
											};
											snp=end;
											for(int snprun=end;snprun<nbsnpperchrinfile[chrtemp1];snprun++)		
											{	segwithav[chrtemp1][snprun]=0;
											};
											typesegment=-1;
											for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++) if (pihatagainstall[relat]<seuilpihat[0] )
											{	nbindivmatchseg[relat][0]=0;
												nbindivmatchseg[relat][1]=0;
												nbindivmatchseg[relat][2]=0;
												nbindivmatchseg[relat][3]=0;
											};
											for(int snprun=0;snprun<nbsnpperchr[chrtemp1];snprun++)
											{	nbsegmentofthislength[snprun]=0;
											};
										} else if (typesegment>-1)
										{	lasttypesegment=typesegment;
										};									
									};		
								};
							};
						};
						int breaknubercm=25;
						int breaknubercmloop=25;
						{	for(int prod=0;prod<1;prod++)
							{	for(int chr=0;chr<23;chr++)
								{	for(int ratio=0;ratio<1;ratio++)
									{	phaseguessright[chr][prod][ratio]=0;
										phaseguesswrong[chr][prod][ratio]=0;
										phasenotguessed[chr][prod][ratio]=0;
										phase1[chr][prod][ratio]=0;
										guessphase1[chr][prod][ratio]=0;
										phase1guessphase1[chr][prod][ratio]=0;
										phase1guessphase2[chr][prod][ratio]=0;
										phase2[chr][prod][ratio]=0;
										guessphase2[chr][prod][ratio]=0;
										phase2guessphase2[chr][prod][ratio]=0;
										phase2guessphase1[chr][prod][ratio]=0;
									};	
								};
							};
							for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
							{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
								{	phaseparent1tap[chrtemp1][snp]=0;
									phaseguesedparent1tap[chrtemp1][snp]=0;
								};
							};							
							char segwithav[23][NSNPPERCHR];
							int countsnp=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchr[chrtemp1];chrdividerrun++)
								{	segwithav[chrtemp1][chrdividerrun]=0; 
								};
							};
							for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
							{	phaseguessright[chrtemp1][0][0]=0;
								phaseguesswrong[chrtemp1][0][0]=0;	
								wrong[chrtemp1]=0;
								right[chrtemp1]=0;
								int phaseparent1=-1;
								int phaseparent1temp=-1;
								int nbsnperror=0;
								int64_t nbonesnp=1;
								int64_t nbzerosnp=0;
								for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
								{	int marker=genomeoffpss[0][snp][chrtemp1];
									int markerp1=genomeoffpss[1][snp][chrtemp1];
									int snperrofind=0;
									if (marker==0 || marker==3)
									{	if (markerp1==3-marker) snperrofind=1;
									};
									if (snperrofind==0)
									{	int phase=-1;
										if (marker!=0 && marker!=3)
										{	if (markerp1==0) phase=marker;
											else if (markerp1==3) phase=3-marker;	
										};
										if (phase>-1)
										{	if (phaseparent1==-1) 
											{	phaseparent1=phase;
												phaseparent1temp=phase;
												printf("phase %d detected at snp %d\n",phase,snp);
											} else if (phaseparent1!=phase) 
											{	if (phaseparent1temp==phase)
												{	printf("phase change at snp %d: %d %d %d\n",snp,phaseparent1temp,phaseparent1,phase);
													phaseparent1=phase;
													phaseparent1temp=phase;
												} else 
												{	phaseparent1temp=phase;
												};
											} else 
											{	phaseparent1temp=phaseparent1;
											}
										};
									} else 
									{	nbsnperror++;
									};
									phaseparent1tap[chrtemp1][snp]=phaseparent1;								
									if (phaseparent1>-1)
									{	if (nbonesnp>1*nbzerosnp && phaseparent1==1) 
										{	phase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase1[ratio][0][0]++; else guessphase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase1guessphase1[ratio][0][0]++; else  phase1guessphase2[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguessright[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=1;
												}
												else 
												{	phaseguesswrong[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=0;
												};
											};											
											right[chrtemp1]++;
										}
										else if (nbonesnp>1*nbzerosnp && phaseparent1==2)
										{	phase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase1[ratio][0][0]++; else guessphase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase2guessphase1[ratio][0][0]++; else phase2guessphase2[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguesswrong[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=1;
												}
												else 
												{	phaseguessright[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=0;
												};
											};
											wrong[chrtemp1]++;
										}
										else if (nbonesnp*1<nbzerosnp && phaseparent1==2)
										{	phase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase2[ratio][0][0]++; else guessphase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase2guessphase2[ratio][0][0]++; else phase2guessphase1[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguessright[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=0;
												}	
												else 
												{	phaseguesswrong[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=1;
												}
											};
											right[chrtemp1]++;
										}
										else if (nbonesnp*1<nbzerosnp && phaseparent1==1) 
										{	phase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase2[ratio][0][0]++; else guessphase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase1guessphase2[ratio][0][0]++; else phase1guessphase1[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguesswrong[chrtemp1][0][0]++;	
													phaseguesedparent1tap[chrtemp1][snp]=0;
												}
												else 
												{	phaseguessright[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=1;
												};	
											};
											wrong[chrtemp1]++;
										}
										else 
										{	for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) phasenotguessed[ratio][0][0]++;
												else phasenotguessed[ratio][0][0]++;
											};
										};			
									} else 
									{	phaseunknow++;
									};
								};
								bestpur=bestpur+(wrong[chrtemp1]>right[chrtemp1]?wrong[chrtemp1]:right[chrtemp1]);
								sumdifftot=sumdifftot+sumdiff[chrtemp1];
								printf("%d %d %d\n",chrtemp1,phaseguessright[chrtemp1][0][0],phaseguesswrong[chrtemp1][0][0]);
								dataperchr[chrtemp1][9]=(phaseguessright[chrtemp1][0][0]>phaseguesswrong[chrtemp1][0][0]?phaseguessright[chrtemp1][0][0]:phaseguesswrong[chrtemp1][0][0])/
									(1+phaseguessright[chrtemp1][0][0]+phaseguesswrong[chrtemp1][0][0]);
								printf("%d %d %d\n",6,phaseguessright[6][0][0],phaseguesswrong[6][0][0]);
							};
							int64_t purcentage=0;
							ratio=2;
							printf("%d %d %d\n",6,phaseguessright[6][0][0],phaseguesswrong[6][0][0]);
							purcentage=(int64_t) 10000*
													(1+(phaseguessright[ratio][0][0]>phaseguesswrong[ratio][0][0]?
															phaseguessright[ratio][0][0]:
															phaseguesswrong[ratio][0][0]))/ 
													(phaseguesswrong[ratio][0][0]+phaseguessright[ratio][0][0]+1.0);
							if (purcentage>bestphasing) bestphasing=purcentage;
							int nbphaseerror=0;
							double pihatthrehjold[23];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	pihatthrehjold[chrtemp1]=0.10;
							};
							int nbphaseright=0;
							int nbphasewrong=0;
							int tabnbphaseright[23][25];
							int tabnbphasewrong[23][25];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	tabnbphaseright[chrtemp1][chrdividerrun]=0;
									tabnbphasewrong[chrtemp1][chrdividerrun]=0;
								};
							};
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	int nbphaseright=0;
									int nbphasewrong=0;
									for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
									{	if (1==phaseparent1tap[chrtemp1][snp]-1) nbphaseright++; else nbphasewrong++;
									};
									printf("%d %d %d %d\n",chrtemp1,chrdividerrun,nbphaseright,nbphasewrong);
									tabnbphaseright[chrtemp1][chrdividerrun]=nbphaseright;
									tabnbphasewrong[chrtemp1][chrdividerrun]=nbphasewrong;
								};	
							}			
							double pihatagainstallchrMP[MAXPOP][2]={0};
							 #pragma omp parallel for
							for(int relat=0;relat<MAXPOP;relat++)
							{	pihatagainstallchrMP[relat][0]=0;
								pihatagainstallchrMP[relat][1]=0;
							}
							int sizebin=450;
							#pragma omp parallel for					
							for(int relat=0;relat<MAXPOP;relat++)
							{	pihatagainstallchrMPphaseerror[0][relat][0][0]=0;
								pihatagainstallchrMPphaseerror[0][relat][0][1]=0;
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][0]*4;chrdividerrun++)
								{	pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][0]=0;
									pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][1]=0;
								};
							}
							int hetsnp=0;
							double exponentinterchr=1;
							for(int relat=0;relat<MAXPOP;relat++)
							{	
							}
							for(int relat=0;relat<MAXPOP;relat++)
							{	if (pihatagainstall[relat]>0)  printf("PIHAT %d %f %f\n",relat,pihatagainstall[relat],seuilpihat[0]);
							}
							float firstexponent=2;
							int powersnppihat=3;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	
								{	hetsnp=0;
									#pragma omp parallel for 				
									for(int relat=0;relat<MAXPOP;relat++)
									{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1]*4;chrdividerrun++)
										{	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][0]=0;
											pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][1]=0;
										};
									};
									for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
									{	double temp[MAXPOP][4];
										for(int64_t relat=0;relat<(int64_t) MAXPOP;relat++) 
										{	temp[relat][0]=0;
											temp[relat][1]=0;
											temp[relat][2]=0;
											temp[relat][3]=0;
										};
										for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
										{	int snpvalue0=genomeoffpss[0][snp][chrtemp1];
											if (snpvalue0!=0 && snpvalue0!=3 )
											{
												hetsnp++;
												double maffloat=1.0*MAF[snp][chrtemp1]/MAXPOP/2;
												int parent0indiv0=(phaseguesedparent1tap[chrtemp1][snp]==1)?(snpvalue0>>1):(snpvalue0&1);
												int parent1indiv0=(phaseguesedparent1tap[chrtemp1][snp]==1)?(snpvalue0&1):(snpvalue0>>1);
												double pcontribu10=((parent0indiv0)-1.0*maffloat);
												double pcontribu11=((parent1indiv0)-1.0*maffloat);
												double maffloatdiviseur=(maffloat)*(2-maffloat*2)/3*2;
												double pconttab[4][2];
												pconttab[0][0]=pow(((0)-maffloat)/maffloatdiviseur*pcontribu10,1.0/powersnppihat);
												pconttab[0][1]=pow(((0)-maffloat)/maffloatdiviseur*pcontribu11,1.0/powersnppihat);
												pconttab[1][0]=pow(((1)-maffloat)/maffloatdiviseur*pcontribu10,1.0/powersnppihat);
												pconttab[1][1]=pow(((1)-maffloat)/maffloatdiviseur*pcontribu11,1.0/powersnppihat); 
												int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
												int snpmod4=(((snp%4)*2));
												int snpdiv4=snp/4;
												printf("%d %f %d %d %f %f %f %f %f %d %d %d \n",snpvalue0,maffloat,parent0indiv0,parent1indiv0,maffloatdiviseur,pconttab[0][0],pconttab[0][1],pconttab[1][0],pconttab[1][1],
															nbsnpperchrby4,snpmod4,snpdiv4);
												#pragma omp parallel for
												for(int64_t relat=0;relat<(int64_t) NbIndiv;relat++) 
												{	if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )
													{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
														if (pconttab[snpvalue1&1][0]>0) temp[relat][0]=temp[relat][0]+(pconttab[snpvalue1&1][0]);
														else temp[relat][0]=0;
														if (temp[relat][0]<0) temp[relat][0]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]<temp[relat][0]) 
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]=temp[relat][0];
														}
														if (pconttab[snpvalue1>>1][0]>0) temp[relat][1]=temp[relat][1]+(pconttab[snpvalue1>>1][0]);
														else temp[relat][1]=0;
														if (temp[relat][1]<0) temp[relat][1]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]<temp[relat][1])
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]=temp[relat][1];
														};
														if (pconttab[snpvalue1&1][1]>0) temp[relat][2]=temp[relat][2]+(pconttab[snpvalue1&1][1]);
														else temp[relat][2]=0;
														if (temp[relat][2]<0) temp[relat][2]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]<temp[relat][2]) 
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]=temp[relat][2];
														}
														if (pconttab[snpvalue1>>1][1]>0) temp[relat][3]=temp[relat][3]+(pconttab[snpvalue1>>1][1]);
														else temp[relat][3]=0;
														if (temp[relat][3]<0) temp[relat][3]=0;
														else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]<temp[relat][3])
														{	
															pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]=temp[relat][3];
														};
													} else 
													{	printf("%d %f %f\n",relat,pihatagainstall[relat],pihatagainstall2[relat]);
													};
												};
											};
										}; printf("%d %d %d %f %f \n",chrtemp1,0,chrdividerrun,
																	pihatagainstallchrMPphaseerror[chrtemp1][0][chrdividerrun*4][0],
																	pihatagainstallchrMPphaseerror[chrtemp1][0][chrdividerrun*4+2][0]);
										printf("%d %d %d %f %f \n",chrtemp1,0,chrdividerrun,
																	pihatagainstallchrMPphaseerror[chrtemp1][1][chrdividerrun*4][0],
																	pihatagainstallchrMPphaseerror[chrtemp1][1][chrdividerrun*4+2][0]);
										#pragma omp parallel for
										for(int relat=0;relat<NbIndiv;relat++) 
										{	if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 )
											{	int value=0;
												{	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]=
																					(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]>pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][value]?
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]:
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][value]);
													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]=
																					(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]>pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][value]?
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]:
																													pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][value]);
												};
											};
										};
										printf("%d %d %d %f %f \n",chrtemp1,0,chrdividerrun,
																	pihatagainstallchrMPphaseerror[chrtemp1][0][chrdividerrun*4][0],
																	pihatagainstallchrMPphaseerror[chrtemp1][0][chrdividerrun*4+2][0]);
										printf("%d %d %d %f %f \n",chrtemp1,0,chrdividerrun,
																	pihatagainstallchrMPphaseerror[chrtemp1][1][chrdividerrun*4][0],
																	pihatagainstallchrMPphaseerror[chrtemp1][1][chrdividerrun*4+2][0]);
									};
								};	
							};
							float thirdindice=3;
							double corglobmax=0;						
							int chrmax1;						
							int chrmax2;						
							int dividmax1;						
							int dividmax2;						
							double sumall[23][25][2][2];
							double sumallsquare[23][25][2][2];
							#pragma omp parallel for 
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int value=0;value<2;value++)
									{	sumall[chrtemp1][chrdividerrun][0][value]=0;
										sumall[chrtemp1][chrdividerrun][1][value]=0;
										sumallsquare[chrtemp1][chrdividerrun][0][value]=0;
										sumallsquare[chrtemp1][chrdividerrun][1][value]=0;
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	sumall[chrtemp1][chrdividerrun][0][value]=sumall[chrtemp1][chrdividerrun][0][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]),firstexponent)*
													(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]>0?1:-1);
											sumall[chrtemp1][chrdividerrun][1][value]=sumall[chrtemp1][chrdividerrun][1][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]),firstexponent)*
													(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]>0?1:-1);
											sumallsquare[chrtemp1][chrdividerrun][0][value]=sumallsquare[chrtemp1][chrdividerrun][0][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]),firstexponent*2);
											sumallsquare[chrtemp1][chrdividerrun][1][value]=sumallsquare[chrtemp1][chrdividerrun][1][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]),firstexponent*2);
										};
									};
								};
							};
							double allcor[23][MAXNBDIVISOR][23][MAXNBDIVISOR];
							double allseg[23][MAXNBDIVISOR][23][MAXNBDIVISOR];
							#pragma omp parallel for		
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
									{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) 
										{	allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=0;
										};
									};
								};
							};			
							int countseg[23][MAXNBDIVISOR];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	countseg[chrtemp1][chrdividerrun]=1;
								};
							};
							float seuil1=0.5;
							#pragma omp parallel for		
							for(int relat=0;relat<NbIndiv;relat++) if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )
							{	
								#pragma omp parallel for		
								for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
								{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
									{	
										double thresholdseg=seuil1*(chrdivider[breaknubercm][chrtemp1][chrdividerrun].end-chrdivider[breaknubercm][chrtemp1][chrdividerrun].start);				
										if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]>thresholdseg || pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]>thresholdseg)
										{	if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]<thresholdseg || pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]<thresholdseg)
											{	countseg[chrtemp1][chrdividerrun]++;
											};
										};
									};
								};
							};
							double sumallseg=0;
							float penalty=0.7;
							double coefforcor=1.0/4;
							#pragma omp parallel for		
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
									{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2)
										{	double sumproduct; 
											sumproduct=0;
											double nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
											{	nbelem++;
												sumproduct=sumproduct+
															pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
															pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
											};
											double  cor1=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
											double sumproduct1=sumproduct;
											sumproduct=0;
											nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
											};
											double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
											sumproduct=0;
											 nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent)*
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
											};
											double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
											sumproduct=0;
											 nbelem=0;
											for(int relat=0;relat<NbIndiv;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 ) 
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent)*	
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
											};
											double  cor4=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
											double corglob=(pow(fabs(cor1),thirdindice-1)*cor1 -
																							pow(fabs(cor2),thirdindice-1)*cor2+ 
																							pow(fabs(cor3),thirdindice-1)*cor3-
																							pow(fabs(cor4),thirdindice-1)*cor4);
											if (corglob<0 && chrtemp1==chrtemp2) corglob=corglob*penalty;
											allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=corglob;
											allcor[chrtemp2][chrdividerrun2][chrtemp1][chrdividerrun]=corglob;
											if (fabs(corglob)>fabs(corglobmax))
											{	corglobmax=corglob;
												chrmax1=chrtemp1;						
												chrmax2=chrtemp2;						
												dividmax1=chrdividerrun;						
												dividmax2=chrdividerrun2;
												printf("%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f  \n",chrtemp1,chrtemp2,chrdividerrun,chrdividerrun2,cor1,cor2,cor3,cor4,
																							sumproduct,nbelem,
																							sumallsquare[chrtemp1][chrdividerrun][0][0],sumallsquare[chrtemp2][chrdividerrun2][0][0],
																							sumall[chrtemp1][chrdividerrun][0][0],sumall[chrtemp2][chrdividerrun2][0][0],
																							sumproduct1,allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2],corglobmax);
											};											
										};
									};
								};
							};
							printf("MAX %d %d %d %d %f\n",chrmax1,dividmax1,chrmax2,dividmax2,corglobmax);
							int group[23][MAXNBDIVISOR];
							int havemerged[23][MAXNBDIVISOR];
							int nbingroup[23][MAXNBDIVISOR];
							#pragma omp parallel for 	
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	havemerged[chrtemp1][chrdividerrun]=0;
									nbingroup[chrtemp1][chrdividerrun]=1;
									group[chrtemp1][chrdividerrun]=chrtemp1+chrdividerrun*23; 
								};
							};	
							int nbmerge=0;
							float mergingexponent[2];
							mergingexponent[0]=2;
							mergingexponent[1]=5; 
							double sumcorglobal=0;
							int summingroup=0;
							int summaxgroup=0;
							int sumgroup=0;
							int nbpair=0;
							do 
							{	group[chrmax2][dividmax2]=chrmax1+dividmax1*23;
								printf("MERGE chr %d chunk %d to chr %d chnk %d :%d %d %d %d corglobmax %f \n",chrmax2,dividmax2,chrmax1,dividmax1,tabnbphaseright[chrmax1][dividmax1],
											tabnbphasewrong[chrmax1][dividmax1],tabnbphaseright[chrmax2][dividmax2],tabnbphasewrong[chrmax2][dividmax2],corglobmax);
								tappointdec[nbmerge].nbgroup1=nbingroup[chrmax1][dividmax1];
								tappointdec[nbmerge].nbgroup2=nbingroup[chrmax2][dividmax2];
								tappointdec[nbmerge].cor=fabs(corglobmax);
								nbmerge++;
								sumcorglobal=sumcorglobal+fabs(corglobmax);
								sumgroup=sumgroup+nbingroup[chrmax1][dividmax1]+nbingroup[chrmax2][dividmax2];
								summingroup=summingroup+((nbingroup[chrmax1][dividmax1]<nbingroup[chrmax2][dividmax2])?nbingroup[chrmax1][dividmax1]:nbingroup[chrmax2][dividmax2]);
								summaxgroup=summaxgroup+((nbingroup[chrmax1][dividmax1]>nbingroup[chrmax2][dividmax2])?nbingroup[chrmax1][dividmax1]:nbingroup[chrmax2][dividmax2]);
								nbingroup[chrmax1][dividmax1]=nbingroup[chrmax1][dividmax1]+nbingroup[chrmax2][dividmax2];
								if (corglobmax>0)
								{	havemerged[chrmax2][dividmax2]=1;
									for(int  chrtemp3=1;chrtemp3<23;chrtemp3++)	
									{	for(int chrdividerrun3=0;chrdividerrun3<nbchrdivider[breaknubercm][chrtemp3];chrdividerrun3++) 
										{	
										};
									};
									if (nbdivisionDEGREE4==1) 
										countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==0) countseg[chrmax1][dividmax1]=(countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==2) countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2]+
																							((countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2]);
									#pragma omp parallel for		
									for(int relat=0;relat<NbIndiv;relat++) 
										if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 )										
									{	
										for(int value=0;value<2;value++)
										{	double p1=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
											double p2=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
											double p3=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]>0?1:-1);
											double p4=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]=pow(fabs(p1+p3),1.0/mergingexponent[value])*(p1+p3>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]=pow(fabs(p2+p4),1.0/mergingexponent[value])*(p2+p4>0?1:-1);
										};
									}	
								} else
								{	havemerged[chrmax2][dividmax2]=-1;
									#pragma omp parallel for
									for(int  chrtemp3=1;chrtemp3<23;chrtemp3++)		
									{	for(int chrdividerrun3=0;chrdividerrun3<nbchrdivider[breaknubercm][chrtemp3];chrdividerrun3++)
										{
										};
									};
									if (nbdivisionDEGREE4==2) 
										countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==0) countseg[chrmax1][dividmax1]=(countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==1) countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2]+
																							((countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2]);
									#pragma omp parallel for		
									for(int relat=0;relat<NbIndiv;relat++)
										if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
									{	for(int value=0;value<2;value++)
										{
											double p1=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
											double p2=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
											double p3=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]>0?1:-1);
											double p4=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]=pow(fabs(p1+p4),1.0/mergingexponent[value])*(p1+p4>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]=pow(fabs(p2+p3),1.0/mergingexponent[value])*(p2+p3>0?1:-1);
										};
									};					
								};
								sumall[chrmax1][dividmax1][0][0]=0;
								sumall[chrmax1][dividmax1][1][0]=0;
								sumallsquare[chrmax1][dividmax1][0][0]=0;
								sumallsquare[chrmax1][dividmax1][1][0]=0;
								sumall[chrmax1][dividmax1][0][1]=0;
								sumall[chrmax1][dividmax1][1][1]=0;
								sumallsquare[chrmax1][dividmax1][0][1]=0;
								sumallsquare[chrmax1][dividmax1][1][1]=0;
								for(int relat=0;relat<NbIndiv;relat++) 
									if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 ) 
								{	for(int value=0;value<2;value++)
									{	sumall[chrmax1][dividmax1][0][value]=sumall[chrmax1][dividmax1][0][value]+
											pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*(0.3+0.7))*
												(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
										sumall[chrmax1][dividmax1][1][value]=sumall[chrmax1][dividmax1][1][value]+
											pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*(0.3+0.7))*
											(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
										sumallsquare[chrmax1][dividmax1][0][value]=sumallsquare[chrmax1][dividmax1][0][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*2);
										sumallsquare[chrmax1][dividmax1][1][value]=sumallsquare[chrmax1][dividmax1][1][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*2);
									};
								};
								int chrtemp1=chrmax1;
								int chrdividerrun=dividmax1;
								corglobmax=0;	
								double highestnew=0;		
								for(int  chrtemp2=1;chrtemp2<23;chrtemp2++)		
								{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (havemerged[chrtemp2][chrdividerrun2]==0 && (chrtemp1!=chrtemp2 || chrdividerrun!=chrdividerrun2))
									{	double sumproduct; 
										sumproduct=0;
										int nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  )											
										{	nbelem++;
											sumproduct=sumproduct+
														pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
														pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
										};
										double  cor1=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
										sumproduct=0;
										sumproduct=0;
										nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
										};
										double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
										sumproduct=0;
										sumproduct=0;
										 nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
										};
										double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
										sumproduct=0;
										sumproduct=0;
										 nbelem=0;
										for(int relat=0;relat<NbIndiv;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011  ) 
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent*(0.3+0.7));
										};
										double  cor4=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
										sumproduct=0;
										float exponantnumer=1;
										double corglob=(pow(fabs(cor1),thirdindice-1)*cor1 -
																						pow(fabs(cor2),thirdindice-1)*cor2+ 
																						pow(fabs(cor3),thirdindice-1)*cor3-
																						pow(fabs(cor4),thirdindice-1)*cor4)
																				*pow((nbingroup[chrtemp1][chrdividerrun]+nbingroup[chrtemp2][chrdividerrun2]),exponantnumer);
										if (highestnew<fabs(corglob)) 
										{	highestnew=fabs(corglob); 
											printf("%f",corglob);
											printf("NB1 %d MB2 %d \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2]);  
										};
										allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=corglob;										
										allcor[chrtemp2][chrdividerrun2][chrtemp1][chrdividerrun]=corglob;										
									};
								};
								nbpair=0;
								for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
								{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++) if (havemerged[chrtemp1][chrdividerrun]==0)
									{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
										{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (havemerged[chrtemp2][chrdividerrun2]==0 && (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2))
											{	nbpair++;
												double corglob=allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2];
												if (fabs(corglob)>fabs(corglobmax))
												{
													corglobmax=corglob;
													chrmax1=chrtemp1;						
													chrmax2=chrtemp2;						
													dividmax1=chrdividerrun;						
													dividmax2=chrdividerrun2;
													printf("NB1 %d MB2 %d c %f \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2],corglob);  
													if (nbingroup[chrtemp1][chrdividerrun]==1 || nbingroup[chrtemp2][chrdividerrun2]==1)
													{	int chrone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrtemp2:chrtemp1;
														int chrdividerrunone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrdividerrun2:chrdividerrun;
														int chrnotone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrtemp2:chrtemp1;
														int chrdividerrunnotone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrdividerrun2:chrdividerrun;
														printf("1 is %d %d\n",chrnotone,chrdividerrunnotone);
														int numberattached=0;
														for(int chrdividerrunrun=0;chrdividerrunrun<nbchrdivider[breaknubercm][chrone];chrdividerrunrun++) 
														{
															if (chrdividerrunrun==chrdividerrunone)
															{	
															} else
															{	if (havemerged[chrone][chrdividerrunrun]==0)
																{
																} else 
																{	int hrdividerloop=chrdividerrunrun;	
																	int chrloop=chrone;
																	int phaseknown=havemerged[chrloop][hrdividerloop];
																	while (iter<1000 && havemerged[chrloop][hrdividerloop]!=0)
																	{	
																		iter++;
																		phaseknown=phaseknown*(havemerged[chrloop][hrdividerloop]!=0?havemerged[chrloop][hrdividerloop]:1);
																		int savechr=chrloop;
																		chrloop=group[chrloop][hrdividerloop]%23;
																		hrdividerloop=group[savechr][hrdividerloop]/23;
																	};
																	if 	(chrloop==chrnotone && hrdividerloop==chrdividerrunnotone)
																	{	printf("window n%d is attached to %d %d with phaseknown %d\n", chrdividerrunrun,chrnotone,chrdividerrunnotone,phaseknown);
																		numberattached=numberattached+phaseknown;
																	} else 
																	{	printf("window n%d is not attached to %d %d\n", chrdividerrunrun,chrnotone,chrdividerrunnotone);
																	};																		
																};
															};
														};
														if (corglobmax*numberattached<0) corglobmax=corglobmax*penalty;
													};
												};	
											};
										};
									};
								};
								printf("Nb pairs: %d\n",nbpair);
							} while (nbpair>0 && nbmerge<22*19-1);
							 nbphaseright=0;
							 nbphasewrong=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		 
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	int chr=chrtemp1;
									int chrdivi=chrdividerrun;
									int phaseknown=1;
									int iter=0;
									while (iter<1000 && havemerged[chr][chrdivi]!=0)
									{	
										iter++;
										phaseknown=phaseknown*(havemerged[chr][chrdivi]!=0?havemerged[chr][chrdivi]:1);
										int savechr=chr;
										chr=group[chr][chrdivi]%23;
										chrdivi=group[savechr][chrdivi]/23;
									} 
									if (iter>900) exit(0);
									printf("pahse chr %d div %d is %d\n",chrtemp1,chrdividerrun,phaseknown);
									chrdivider[breaknubercm][chrtemp1][chrdividerrun].phasing=phaseknown;
									for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
									{	if ((phaseknown+3)/2==2) 
										{	if (genomeoffpss[0][snp][chrtemp1]==1) 
											{	genomeoffpss[0][snp][chrtemp1]=2;
												*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
													(*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
														(~(3<<((snp%4)*2))))
															|	((2)<<((snp%4)*2));
											}
											else if (genomeoffpss[0][snp][chrtemp1]==2) 
											{	genomeoffpss[0][snp][chrtemp1]=1;
												*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
													(*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
														(~(3<<((snp%4)*2))))
															|	((1)<<((snp%4)*2));
											};
										};
										if ((phaseknown+3)/2==phaseparent1tap[chrtemp1][snp]) 
										{	nbphaseright++;
											chrdivider[breaknubercm][0][0].nbright++;
											chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbright++;
										} else 
										{	nbphasewrong++;
											chrdivider[breaknubercm][0][0].nbwrong++;
											chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbwrong++;
										};
									};
									printf("chr %d div %d nbphaseright %d nbphasewrong %d\n",chrtemp1,chrdividerrun,nbphaseright,nbphasewrong);
								};
								printf("chr %d nbphaseright %d nbphasewrong %d ratio %f \n",chrtemp1,nbphaseright,nbphasewrong,
																							(nbphaseright>nbphasewrong?nbphaseright:nbphasewrong)/(1.0+nbphaseright+nbphasewrong));
							};
							purcentage=(int64_t) 10000*
											(nbphaseright>nbphasewrong?nbphaseright:nbphasewrong)/(nbphaseright+nbphasewrong>0?nbphaseright+nbphasewrong:1);						
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)	
							{	dataperchr[chrtemp1][52+nbdivisionDEGREE4-18+3*powerpihatDEGREE4+breaknubercmloop-23]=purcentage;		
							};
						}
					};
				};					
			};	
		};
	} while (gen<MAXGEN);	
	return (0);
}
	int writeoutput(int chr, char pathoutput[] )
{	FILE *fp; 
    char filename[200]; 
    strcpy(filename,pathoutput);
	char number[100];
    sprintf(number, "%d", chr);
	strcat(filename,number);
	strcat(filename,".ped");
	printf("The file to write the data is %s\n", filename); 
	fp = fopen(filename, "w"); 
    if (fp == NULL) 
    {   printf("Could not open file %s\n", filename); 
        return (0); 
    };
	int32_t person;
	for(person=0;person<NbIndiv;person++)
	{	fprintf(fp,	"%d %d 0 0 -9 0 ",person,person);
		for(int snp=0;snp<nbsnpperchrinfile[chr];snp++)
		{	int geno=(*((genomes[chr]+(unsigned long long) person*(nbsnpperchr[chr]/4+((nbsnpperchr[chr]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			geno=rand()%4; 
			fprintf(fp, "%d %d ", geno%2, geno/2);
		};
		fprintf(fp,"\n");
	};
	fclose (fp);
	return (0);
}
	int run(int argc, char* argv[])
	{
		clock_t begin0 = clock();
		printdetail = 1;
		int chrtodo1 = 1;
		char PathInput[200];
		char PathOutput[200];
		char PathParentInfo[200] = "";
		char PathListIndiv[200] = "";
		for (int input = 1; input < argc; input++)
		{
			if (strncmp(argv[input], "-NbIndiv", strlen("-NbIndiv")) == 0 && input < argc - 1)
				NbIndiv = atoi(argv[++input]);
			else if (strncmp(argv[input], "-PathInput", strlen("-PathInput")) == 0 && input < argc - 1)
				strcpy(PathInput, argv[++input]);
			else if (strncmp(argv[input], "-PathOutput", strlen("-PathOutput")) == 0 && input < argc - 1)
				strcpy(PathOutput, argv[++input]);
			else if (strncmp(argv[input], "-PathParentInfo", strlen("-PathParentInfo")) == 0 && input < argc - 1)
				strcpy(PathParentInfo, argv[++input]);
			else if (strncmp(argv[input], "-ListIndiv", strlen("-ListIndiv")) == 0 && input < argc - 1)
				strcpy(PathListIndiv, argv[++input]);
		}
		if (NbIndiv == 0)
		{
			printf("ERROR: Number of indivudals is zero or undefined\n");
			exit(0);
		}
		if (NbIndiv > NBINDIVMAX)
		{
			printf("ERROR: Number of indivudals is higher than %d\n", NBINDIVMAX);
			exit(0);
		}
		nbsnpperchr[0] = 330005;
		nbsnpperchr[1] = 26229;
		nbsnpperchr[2] = 26210;
		nbsnpperchr[3] = 22209;
		nbsnpperchr[4] = 20690;
		nbsnpperchr[5] = 19027;
		nbsnpperchr[6] = 18418;
		nbsnpperchr[7] = 18367;
		nbsnpperchr[8] = 16283;
		nbsnpperchr[9] = 14990;
		nbsnpperchr[10] = 16494;
		nbsnpperchr[11] = 15818;
		nbsnpperchr[12] = 16008;
		nbsnpperchr[13] = 11510;
		nbsnpperchr[14] = 10804;
		nbsnpperchr[15] = 10884;
		nbsnpperchr[16] = 12195;
		nbsnpperchr[17] = 11486;
		nbsnpperchr[18] = 10222;
		nbsnpperchr[19] = 9806;
		nbsnpperchr[20] = 8985;
		nbsnpperchr[21] = 5227;
		nbsnpperchr[22] = 5882;
		for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
		{
			chrdivider[25][chrtemp1][0].start = 0;
			chrdivider[25][chrtemp1][0].end = nbsnpperchr[chrtemp1];
			nbchrdivider[25][chrtemp1] = 1;
		}
		srand(0);
		for (int method = 0; method < 4; method++)
		{
			for (int res = 0; res < 6; res++)
			{
				result[method][res] = 0;
			}
		}
		if (printdetail) printf("read trio done \n");
		if (printdetail) printf("read ID done \n");
		printf("%d\n", sizeof(genomes));
		clock_t step1 = clock();
		float elapsed_secs1 = (float)(step1 - begin0);
		printf("\nTotal time 1:%f cpu click so %f seconds\n", elapsed_secs1, elapsed_secs1 / CLOCKS_PER_SEC);
		clock_t step11 = clock();
		float elapsed_secs11 = (float)(step11 - step1);
		printf("\nTotal time 1:%f cpu click so %f seconds\n", elapsed_secs11, elapsed_secs11 / CLOCKS_PER_SEC);
		for (int IDrun = 0; IDrun < 12; IDrun++)
		{
			distrigametic[IDrun] = 0;
		}
		for (int IDrun = 0; IDrun < 23; IDrun++)
		{
			distrigametickeep[IDrun] = 0;
		}
		if (strlen(PathParentInfo) > 0)
		{
			if (readParentInfo(PathParentInfo) != 0)
			{
				printf("ERROR: Failed to read parent info file. Exiting.\n");
				exit(1);
			}
			printf("Using parent information from file: %s\n", PathParentInfo);
		}
		else
		{
			printf("No parent info file provided. All individuals will be processed without known parents.\n");
		}
		
		if (strlen(PathListIndiv) > 0)
		{
			if (readListIndiv(PathListIndiv) != 0)
			{
				printf("ERROR: Failed to read individual list file. Exiting.\n");
				exit(1);
			}
			printf("Processing only %d individuals from the list\n", nbIndivToProcess);
		}
		else
		{
			printf("Processing all %d individuals\n", NbIndiv);
		}
		
		for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
		{
			genomes[chrtemp1] = (unsigned char*)calloc((unsigned long long)(1 + nbsnpperchr[chrtemp1] / 4) * NbIndiv, sizeof(char));
			readgenomelocal(PathInput, chrtemp1, 100 + chrtemp1, 105, genomes[chrtemp1]);
		}
	
		if (useListIndiv)
		{
			for (int idx = 0; idx < nbIndivToProcess; idx++)
			{
				int indiv = listIndivToProcess[idx];
				printf("start individual: %d (%d/%d)\n", indiv, idx + 1, nbIndivToProcess); 
				findrelative(indiv,
					indiv,
					indiv,
					indiv,
					0,
					2,
					0,
					PathInput);
				int parent1 = -1;
				int parent2 = -1;
				if (strlen(PathParentInfo) > 0)
				{
					getParentInfo(indiv, &parent1, &parent2);
				}
				if (parent1 == -1 || parent2 == -1 || parent1 == indiv || parent2 == indiv)
				{
					if (parent1 == -1 || parent2 == -1)
					{
						printf("Individual %d: parents unknown, using predict_phasing_without_parents\n", indiv);
					}
					else
					{
						printf("Individual %d: invalid parent IDs (parent1=%d, parent2=%d), using predict_phasing_without_parents\n",
							indiv, parent1, parent2);
					}
					predict_phasing_without_parents(indiv,
						chrtodo1,
						0,
						indiv,
						0,
						0,
						indiv,
						indiv,
						0,
						2,
						0,
						PathOutput);
				}
				else
				{
					printf("Individual %d: parents known (parent1=%d, parent2=%d), using predict_phasing_with_parents_providing_GT\n",
						indiv, parent1, parent2);
					predict_phasing_with_parents_providing_GT(indiv,
						chrtodo1,
						0,
						indiv,
						0,
						0,
						parent1,
						parent2,
						0,
						2,
						0,
						PathOutput);
				}
				for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
				{
					writeoutput(chrtemp1, PathOutput);
				}
			}
		}
		else
		{
			for (int indiv = 0; indiv < NbIndiv; indiv++)
			{
				printf("start individual: %d\n", indiv); 
				findrelative(indiv,
					indiv,
					indiv,
					indiv,
					0,
					2,
					0,
					PathInput);
				int parent1 = -1;
				int parent2 = -1;
				if (strlen(PathParentInfo) > 0)
				{
					getParentInfo(indiv, &parent1, &parent2);
				}
				if (parent1 == -1 || parent2 == -1 || parent1 == indiv || parent2 == indiv)
				{
					if (parent1 == -1 || parent2 == -1)
					{
						printf("Individual %d: parents unknown, using predict_phasing_without_parents\n", indiv);
					}
					else
					{
						printf("Individual %d: invalid parent IDs (parent1=%d, parent2=%d), using predict_phasing_without_parents\n",
							indiv, parent1, parent2);
					}
					predict_phasing_without_parents(indiv,
						chrtodo1,
						0,
						indiv,
						0,
						0,
						indiv,
						indiv,
						0,
						2,
						0,
						PathOutput);
				}
				else
				{
					printf("Individual %d: parents known (parent1=%d, parent2=%d), using predict_phasing_with_parents_providing_GT\n",
						indiv, parent1, parent2);
					predict_phasing_with_parents_providing_GT(indiv,
						chrtodo1,
						0,
						indiv,
						0,
						0,
						parent1,
						parent2,
						0,
						2,
						0,
						PathOutput);
				}
				for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
				{
					writeoutput(chrtemp1, PathOutput);
				}
			}
		}
		return 0;
	}
};

PhasingProgram g_program;

int readgenomelocal(char pathfile[], int chr, int run, int step, unsigned char* gentomodify)
{
	return g_program.readgenomelocal(pathfile, chr, run, step, gentomodify);
}

int readParentInfo(const char* filename)
{
	return g_program.readParentInfo(filename);
}

int getParentInfo(int indivID, int* parent1, int* parent2)
{
	return g_program.getParentInfo(indivID, parent1, parent2);
}

int readListIndiv(const char* filename)
{
	return g_program.readListIndiv(filename);
}

int findrelative(int ID, int numtrio, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[])
{
	return g_program.findrelative(ID, numtrio, IDp1loop, IDp2loop, lenminseg, version, gentostart, pathresult);
}

int predict_phasing_without_parents(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[])
{
	return g_program.predict_phasing_without_parents(ID, chr1, chr2, numtrio, limitnbsnp, run, IDp1loop, IDp2loop, lenminseg, version, gentostart, pathresult);
}

int predict_phasing_with_parents_providing_GT(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[])
{
	return g_program.predict_phasing_with_parents_providing_GT(ID, chr1, chr2, numtrio, limitnbsnp, run, IDp1loop, IDp2loop, lenminseg, version, gentostart, pathresult);
}

int writeoutput(int chr, char pathoutput[])
{
	return g_program.writeoutput(chr, pathoutput);
}

int main(int argc, char* argv[])
{
	return g_program.run(argc, argv);
}
