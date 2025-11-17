#ifndef PHASING_H
#define PHASING_H

#include "../types.h"

int findrelative(int ID, int numtrio, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[],
	int NbIndiv, int nbsnpperchr[], int nbsnpperchrinfile[], int MAF[][23], unsigned char* genomes[], 
	float pihatagainstall[], int bestpihatagainstallID[], int* IDbestpihat, int* IDbestpihat2, 
	float bestpihat[], int* placefirttoconsider, float seuilpihat[]);

int predict_phasing_without_parents(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, 
	int lenminseg, int version, int gentostart, char pathresult[], int NbIndiv, int nbsnpperchr[], int nbsnpperchrinfile[], 
	int MAF[][23], unsigned char* genomes[], int genomeoffpss[][NSNPPERCHR][23], float seuilpihat[], 
	typechrdivider chrdivider[][23][20], int nbchrdivider[][23], int* IDjob, int nbbreak, int nbrelatpihat,
	float pihatagainstall[], int bestpihatagainstallID[], float pihatagainstall2[], pointdecision tappointdec[],
	int64_t relatpihatchr[][23][2], int* placefirttoconsider);

int predict_phasing_with_parents_providing_GT(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, 
	int lenminseg, int version, int gentostart, char pathresult[], int NbIndiv, int nbsnpperchr[], int nbsnpperchrinfile[], 
	int MAF[][23], unsigned char* genomes[], int genomeoffpss[][NSNPPERCHR][23], float seuilpihat[], 
	typechrdivider chrdivider[][23][20], int nbchrdivider[][23], int* IDjob, int nbbreak, int nbrelatpihat,
	float pihatagainstall[], int bestpihatagainstallID[], float pihatagainstall2[], pointdecision tappointdec[],
	int64_t relatpihatchr[][23][2], int* placefirttoconsider);

#endif

