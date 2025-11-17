#ifndef TYPES_H
#define TYPES_H

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
{
	int chi;
	unsigned char haptohap;
} resultphaing;

typedef struct
{
	int ID;
	unsigned char coef;
} indivclose;

typedef struct
{
	int ID;
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
{
	double cor;
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
{
	int ID;
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

#endif

