#ifndef FILEIO_H
#define FILEIO_H

#include "../types.h"

// HAP file functions
int readgenomelocal(char pathfile[], int chr, int run, int step, unsigned char* gentomodify, int NbIndiv, int nbsnpperchr[], int nbsnpperchrinfile[]);

// PED file functions
int readpedfile(const char* filename, int chr, unsigned char* genomes[], int nbsnpperchr[], int* nbsnpperchrinfile, int* NbIndiv, int maxIndiv);
int writepedfile(int chr, const char* pathoutput, unsigned char* genomes[], int nbsnpperchr[], int nbsnpperchrinfile[], int NbIndiv);

// Parent and individual list functions
int readParentInfo(const char* filename, structparentinfo parentinfo[], int* nbparentinfo, int maxIndiv);
int getParentInfo(int indivID, int* parent1, int* parent2, structparentinfo parentinfo[], int nbparentinfo);
int readListIndiv(const char* filename, int listIndivToProcess[], int* nbIndivToProcess, int* useListIndiv, int NbIndiv, int maxIndiv);

// Legacy output function (writes PED format)
int writeoutput(int chr, char pathoutput[], unsigned char* genomes[], int nbsnpperchr[], int nbsnpperchrinfile[], int NbIndiv);

#endif

