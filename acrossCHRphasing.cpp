/***************************************************************************************************************************************************
*
*                                 This program will read segment and output condense files
*                                 
****************************************************************************************************************************************************/
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef __linux__
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sys/syscall.h>
#endif
#include <omp.h>

#define NMAX 664 
#define NCHR 23
#define NSNPPERCHR 26230
#define MAXSNP 330005
#define MAXPOP 435188
#define MAXCLOSERELAT 6000
#define MAXCLOSERELATTEMP 450000
#define MAXID 6026457
#define nbbyteforsegment 7
#define MAXSEGINDIV 10000
#define MAXNBPAIR  152252 
#define MAXNBDIVISOR 25

#define MAXNBTRIO 978
#define NBINDIVEA 1
#define INCREMENTLOOP 1

int bestpihatagainstallID[100];
int IDbestpihat;

float pihatagainstall[MAXPOP];
float pihatagainstall2[MAXPOP];
int nbchrdivider[51][23];

int sumtotindivcount=0;
int placefirttoconsider=0;

typedef struct 
{	int start;
	int end;
	int phasing;
	int nbright;
	int nbwrong;
	int segment;
} typechrdivider;

#define MAXSEGPERCHR 25  /* With minCM=22, hotspots yield ~15 segments max; window uses 5. Matches MAXNBDIVISOR. */
typechrdivider chrdivider[50][23][MAXSEGPERCHR];

#define MAXHOTSPOTS 500
typedef struct { int chr; double cM; double cMperMb; int bp; } hotspot_t;
static int fill_chrdivider_from_hotspots(hotspot_t *all, int ntotal);

float seuilpihat[3];

char pathHap[512] = {0};      /* Path for .haps files (e.g. data/chr -> chr1.haps) */
char pathHotspots[512] = {0}; /* Path to hotspot file (Chr cM cMperMb bp). If set, must exist. */
char pathMap[512] = {0};     /* Path prefix for map files (e.g. /path/SNPCHR -> SNPCHR1.map). Col3=cM Col4=bp. */
char pathResult[512] = {0};  /* Output directory */
char pathwindow[512] = {0};  /* (Legacy) Path to window bounds file. Used only if no hotspot/map. */
int useBinary = 0;          /* When 1, read .bin (raw packed) instead of text .haps */
int minCM = 22;              /* Minimum chunk length in cM for hotspot segments (R: MIN.CM) */
int posOffspring = -1;    /* 0-based position of individual to phase in hap file */
int posParent1 = -1;      /* 0-based position of parent 1 in hap file (-1 = no parent 1) */
int posParent2 = -1;      /* 0-based position of parent 2 in hap file (-1 = no parent 2) */
int nbsnpperchr[23];
int n_individuals = 0;  /* Detected from hap file header, or 0 if using fallback */
int MAF[NSNPPERCHR][23]={0};
double allele_freq[NSNPPERCHR][23]={0};  /* p_k = allele frequency for pi-hat formula */

int genomeoffpss[3][NSNPPERCHR][23];	

unsigned char * genomes[23];

typedef struct
{	double cor;
	int nbgroup1;
	int nbgroup2;

} pointdecision;

pointdecision tappointdec[MAXNBDIVISOR*22];

int nbbreak=0;
#define MAXBREAK 1000

uint64_t * tempload;

/* hotspot-based segment helpers (R logic) */
static int hotspot_cmp_cMperMb(const void *a, const void *b) {
	const hotspot_t *ha = (const hotspot_t *)a, *hb = (const hotspot_t *)b;
	return (ha->cMperMb > hb->cMperMb) ? -1 : ((ha->cMperMb < hb->cMperMb) ? 1 : 0);
}
static int read_hotspot_file(const char *path, hotspot_t *out, int *nout) {
	FILE *f = fopen(path, "r");
	if (!f) return -1;
	char line[1024];
	int n = 0;
	if (fgets(line, sizeof(line), f)); /* skip header */
	while (n < MAXHOTSPOTS && fgets(line, sizeof(line), f)) {
		char chrName[32];
		double cM, cMperMb;
		int bp;
		if (sscanf(line, "%31s %lf %lf %d", chrName, &cM, &cMperMb, &bp) >= 4) {
			int chr = 0;
			if (sscanf(chrName, "chr%d", &chr) != 1)
				sscanf(chrName, "%d", &chr);
			if (chr >= 1 && chr <= NCHR) {
				out[n].chr = chr;
				out[n].cM = cM;
				out[n].cMperMb = cMperMb;
				out[n].bp = bp;
				n++;
			}
		}
	}
	fclose(f);
	*nout = n;
	return 0;
}
/* Map format: chr rsid genpos bp ref alt  (tab/space; col3=cM col4=bp). Tries .map then .TXT */
static FILE* open_map_file(int chr, const char *pathPrefix) {
	char path[512];
	snprintf(path, sizeof(path), "%s%d.map", pathPrefix, chr);
	FILE *f = fopen(path, "r");
	if (f) return f;
	snprintf(path, sizeof(path), "%s%d.TXT", pathPrefix, chr);
	return fopen(path, "r");
}
static int read_map_bp(int chr, const char *pathPrefix, int *bp_out, int max_snps) {
	FILE *f = open_map_file(chr, pathPrefix);
	if (!f) return -1;
	int n = 0;
	char line[512];
	char rsid[128];
	double genpos;
	int bp;
	while (n < max_snps && fgets(line, sizeof(line), f)) {
		if (sscanf(line, "%*d %127s %lf %d", rsid, &genpos, &bp) >= 3) {
			bp_out[n++] = bp;
		}
	}
	fclose(f);
	return n;
}
/* Read map with cM and bp. Returns n SNPs, or -1 on error. */
static int read_map_full(int chr, const char *pathPrefix, double *cM_out, int *bp_out, int max_snps) {
	FILE *f = open_map_file(chr, pathPrefix);
	if (!f) return -1;
	int n = 0;
	char line[512];
	char rsid[128];
	while (n < max_snps && fgets(line, sizeof(line), f)) {
		if (sscanf(line, "%*d %127s %lf %d", rsid, &cM_out[n], &bp_out[n]) >= 3) n++;
	}
	fclose(f);
	return n;
}
/* Count SNPs in map file (for binary mode). Returns -1 on error. */
static int count_map_snps(int chr, const char *pathPrefix) {
	FILE *f = open_map_file(chr, pathPrefix);
	if (!f) return -1;
	int n = 0;
	char line[512];
	char rsid[128];
	double genpos;
	int bp;
	while (fgets(line, sizeof(line), f)) {
		if (sscanf(line, "%*d %127s %lf %d", rsid, &genpos, &bp) >= 3) n++;
	}
	fclose(f);
	return n;
}

/* Count SNPs in text .haps file (one line per SNP). Returns -1 on error. Uses buffered read for speed on large files. */
static int count_haps_snps(int chr, const char *pathHapDir) {
	char pathfile[512];
	snprintf(pathfile, sizeof(pathfile), "%s%d.haps", pathHapDir, chr);
	FILE *f = fopen(pathfile, "r");
	if (!f) return -1;
	int n = 0;
	char buf[65536];
	size_t nr;
	unsigned long chunks = 0;
	while ((nr = fread(buf, 1, sizeof(buf), f)) > 0) {
		for (size_t i = 0; i < nr; i++)
			if (buf[i] == '\n') n++;
		chunks++;
		if (chunks % 10000 == 0) { printf("."); fflush(stdout); }  /* ~640MB per dot */
	}
	fclose(f);
	return n;
}

/* Populate nbsnpperchr before loading HAP files, for segment boundary setup.
 * Uses: useBinary -> BINARY_NBSNPPERCHR; pathMap -> count_map_snps; else -> count_haps_snps. */
static int init_nbsnpperchr_before_haps(void) {
	printf("[window setup] Step 1: Getting SNP counts per chromosome (nbsnpperchr)\n"); fflush(stdout);
	if (useBinary) {
		printf("[window setup]   Using useBinary=1 -> BINARY_NBSNPPERCHR constants\n");
		static const int BINARY_NBSNPPERCHR[23] = {
			330005, 26229, 26210, 22209, 20690, 19027, 18418, 18367, 16283, 14990,
			16494, 15818, 16008, 11510, 10804, 10884, 12195, 11486, 10222, 9806,
			8985, 5227, 5882
		};
		for (int chr = 1; chr <= 22; chr++) {
			nbsnpperchr[chr] = BINARY_NBSNPPERCHR[chr];
			printf("[window setup]   chr%d n_snps=%d\n", chr, nbsnpperchr[chr]);
		}
		printf("[window setup] Step 1 done: nbsnpperchr set for chr 1-22\n");
		return 0;
	}
	if (pathMap[0]) {
		printf("[window setup]   Using pathMap=%s -> count_map_snps per chr\n", pathMap); fflush(stdout);
		for (int chr = 1; chr <= 22; chr++) {
			printf("[window setup]   chr %d/22...", chr); fflush(stdout);
			int n = count_map_snps(chr, pathMap);
			if (n <= 0 || n > MAXSNP) {
				printf("EXIT(1): chr%d map file invalid (n_snps=%d). Check pathMap.\n", chr, n);
				return -1;
			}
			nbsnpperchr[chr] = n;
			printf(" n_snps=%d\n", n); fflush(stdout);
		}
		printf("[window setup] Step 1 done: nbsnpperchr set for chr 1-22\n"); fflush(stdout);
		return 0;
	}
	printf("[window setup]   Using pathHap=%s -> count_haps_snps (line count) per chr\n", pathHap); fflush(stdout);
	for (int chr = 1; chr <= 22; chr++) {
		printf("[window setup]   Scanning chr %d/22...", chr); fflush(stdout);
		int n = count_haps_snps(chr, pathHap);
		if (n <= 0 || n > MAXSNP) {
			printf(" EXIT(1): n_snps=%d invalid. Check pathHap.\n", n); fflush(stdout);
			return -1;
		}
		nbsnpperchr[chr] = n;
		printf(" done n_snps=%d\n", n); fflush(stdout);
	}
	printf("[window setup] Step 1 done: nbsnpperchr set for chr 1-22\n"); fflush(stdout);
	return 0;
}

/* Setup segment boundaries (chrdivider, nbchrdivider) before loading HAP. Requires nbsnpperchr to be set. */
static int setup_segment_boundaries(void) {
	printf("[window setup] Step 2: Initializing chrdivider (zeroing segment field)\n"); fflush(stdout);
	for (int size = 0; size < 51; size++)
		for (int chr = 1; chr < 23; chr++)
			chrdivider[size][chr][0].segment = 0;
	printf("[window setup]   chrdivider[0..50][1..22][0].segment = 0\n"); fflush(stdout);

	const int B = 27;
	if (pathHotspots[0] && pathMap[0]) {
		printf("[window setup] Step 3: Using HOTSPOT mode (pathHotspots=%s pathMap=%s minCM=%d)\n", pathHotspots, pathMap, minCM);
		hotspot_t *all = (hotspot_t *)malloc(MAXHOTSPOTS * sizeof(hotspot_t));
		if (!all) { printf("EXIT(1): Failed to allocate hotspot buffer\n"); return -1; }
		int ntotal = 0;
		if (read_hotspot_file(pathHotspots, all, &ntotal) != 0) {
			printf("EXIT(1): Hotspot file not found or invalid: %s\n", pathHotspots);
			fprintf(stderr, "ERROR: Cannot open hotspot file %s. Check pathHotspots.\n", pathHotspots);
			free(all);
			return -1;
		}
		printf("[window setup]   Read %d hotspots from file\n", ntotal);
		if (fill_chrdivider_from_hotspots(all, ntotal) != 0) {
			printf("EXIT(1): Failed to build segments. Check pathMap=%s and map file format.\n", pathMap);
			free(all);
			return -1;
		}
		free(all);
		printf("[window setup] Step 4: Segment boundaries built from hotspots. Per-chr results:\n");
		for (int chr = 1; chr <= 22; chr++) {
			printf("[window setup]   chr%d: %d segments", chr, nbchrdivider[B][chr]);
			for (int s = 0; s < nbchrdivider[B][chr]; s++)
				printf(" [%d,%d)", chrdivider[B][chr][s].start, chrdivider[B][chr][s].end);
			printf("\n");
		}
		printf("[window setup]   max segments across chr: %d\n", nbchrdivider[B][0]);
	} else if (pathwindow[0]) {
		printf("[window setup] Step 3: Using WINDOW mode (pathwindow=%s)\n", pathwindow); fflush(stdout);
		int CHR27_BOUNDS[23][5];
		FILE *fb = fopen(pathwindow, "r");
		if (!fb) {
			printf("EXIT(1): Cannot open window file %s. Check pathwindow.\n", pathwindow);
			fprintf(stderr, "ERROR: Cannot open %s (window coordinates). Check pathwindow.\n", pathwindow);
			return -1;
		}
		for (int c = 0; c < 23; c++)
			CHR27_BOUNDS[c][0] = CHR27_BOUNDS[c][1] = CHR27_BOUNDS[c][2] = CHR27_BOUNDS[c][3] = CHR27_BOUNDS[c][4] = 0;
		char line[256];
		printf("[window setup]   Reading window bounds (5 boundaries per chr):\n");
		for (int chr = 0; chr < 23; chr++) {
			if (!fgets(line, sizeof(line), fb)) {
				printf("EXIT(1): %s too short (need 23 lines)\n", pathwindow);
				fclose(fb);
				return -1;
			}
			while (line[0] == '#' || line[0] == '\n')
				{ if (!fgets(line, sizeof(line), fb)) { printf("EXIT(1): %s read error (chr %d)\n", pathwindow, chr); fclose(fb); return -1; } }
			int dummy, b0, b1, b2, b3, b4;
			if (sscanf(line, "%d %d %d %d %d %d", &dummy, &b0, &b1, &b2, &b3, &b4) >= 5) {
				CHR27_BOUNDS[chr][0]=b0; CHR27_BOUNDS[chr][1]=b1; CHR27_BOUNDS[chr][2]=b2; CHR27_BOUNDS[chr][3]=b3; CHR27_BOUNDS[chr][4]=b4;
			} else if (sscanf(line, "%d %d %d %d %d", &b0, &b1, &b2, &b3, &b4) >= 5) {
				CHR27_BOUNDS[chr][0]=b0; CHR27_BOUNDS[chr][1]=b1; CHR27_BOUNDS[chr][2]=b2; CHR27_BOUNDS[chr][3]=b3; CHR27_BOUNDS[chr][4]=b4;
			}
			if (chr >= 1 && chr <= 22)
				printf("[window setup]   chr%d raw bounds: %d %d %d %d %d\n", chr, b0, b1, b2, b3, b4);
		}
		fclose(fb);
		nbchrdivider[B][0] = 0;
		printf("[window setup] Step 4: Converting bounds to chrdivider[B][chr] segments:\n");
		for (int chr = 1; chr <= 22; chr++) {
			int start = 0, n = 0;
			for (int i = 0; i < 5; i++) {
				int b = CHR27_BOUNDS[chr][i];
				if (b < 0) {
					chrdivider[B][chr][n].start = start;
					chrdivider[B][chr][n].end = nbsnpperchr[chr];
					n++;
					break;
				}
				chrdivider[B][chr][n].start = start;
				chrdivider[B][chr][n].end = b;
				start = b;
				n++;
			}
			nbchrdivider[B][chr] = n;
			if (n > nbchrdivider[B][0]) nbchrdivider[B][0] = n;
			printf("[window setup]   chr%d: %d segments", chr, n);
			for (int s = 0; s < n; s++)
				printf(" [%d,%d)", chrdivider[B][chr][s].start, chrdivider[B][chr][s].end);
			printf(" (nbsnpperchr=%d)\n", nbsnpperchr[chr]);
		}
		printf("[window setup]   max segments across chr: %d\n", nbchrdivider[B][0]);
	} else {
		printf("EXIT(1): Segment boundaries required. Provide pathwindow= OR (pathHotspots= + pathMap=).\n");
		return -1;
	}
	printf("[window setup] Step 5: Window sizes definition complete (B=27, breaknubercm=27)\n");
	return 0;
}

/* Compute hotspots from map: cMperMb = (delta_cM/delta_bp)*1e6 between consecutive SNPs. */
static int compute_hotspots_from_map(hotspot_t *out, int *nout) {
	*nout = 0;
	double *cM = (double *)malloc(NSNPPERCHR * sizeof(double));
	int *bp = (int *)malloc(NSNPPERCHR * sizeof(int));
	if (!cM || !bp) { free(cM); free(bp); return -1; }
	for (int chr = 1; chr <= 22; chr++) {
		int n = read_map_full(chr, pathMap, cM, bp, NSNPPERCHR);
		if (n < 2) continue;
		for (int i = 0; i < n - 1 && *nout < MAXHOTSPOTS; i++) {
			double dcM = cM[i+1] - cM[i];
			int dbp = bp[i+1] - bp[i];
			if (dbp <= 0) continue;
			double cMperMb = (dcM * 1e6) / (double)dbp;
			out[*nout].chr = chr;
			out[*nout].cM = cM[i+1];
			out[*nout].cMperMb = cMperMb;
			out[*nout].bp = bp[i+1];
			(*nout)++;
		}
	}
	free(cM);
	free(bp);
	return 0;
}

/* Compute chrdivider[B] from hotspots (R logic): top 10 per chr, merge if gap < minCM, map bp->SNP index. */
static int fill_chrdivider_from_hotspots(hotspot_t *all, int ntotal) {
	int *map_bp = (int *)malloc(NSNPPERCHR * sizeof(int));
	if (!map_bp) return -1;

	const int B = 27;
	nbchrdivider[B][0] = 0;

	for (int chr = 1; chr <= 22; chr++) {
		/* Filter by chr, take top 10 by cMperMb */
		hotspot_t sub[12];
		int nsub = 0;
		for (int i = 0; i < ntotal && nsub < 12; i++) {
			if (all[i].chr != chr) continue;
			sub[nsub++] = all[i];
		}
		if (nsub == 0) {
			chrdivider[B][chr][0].start = 0;
			chrdivider[B][chr][0].end = nbsnpperchr[chr];
			chrdivider[B][chr][0].segment = 0;
			nbchrdivider[B][chr] = 1;
			if (1 > nbchrdivider[B][0]) nbchrdivider[B][0] = 1;
			continue;
		}
		qsort(sub, nsub, sizeof(hotspot_t), hotspot_cmp_cMperMb);
		nsub = (nsub > 10) ? 10 : nsub;

		int n_map = read_map_bp(chr, pathMap, map_bp, nbsnpperchr[chr]);
		if (n_map <= 0) {
			free(map_bp);
			return -1;
		}

		/* Build bp boundaries from top hotspots, sorted by cM */
		double cMlist[12];
		int bplist[12];
		int nb = 0;
		for (int i = 0; i < nsub; i++) {
			cMlist[nb] = sub[i].cM;
			bplist[nb] = sub[i].bp;
			nb++;
		}
		/* Sort by cM for merging */
		for (int i = 0; i < nb; i++)
			for (int j = i + 1; j < nb; j++)
				if (cMlist[i] > cMlist[j]) {
					double tc = cMlist[i]; cMlist[i] = cMlist[j]; cMlist[j] = tc;
					int tb = bplist[i]; bplist[i] = bplist[j]; bplist[j] = tb;
				}

		/* Add chr start and end */
		double cM0 = 0.0, cMend = cMlist[nb - 1] + 1.0;
		int bp0 = (n_map > 0) ? map_bp[0] : 0;
		int bpend = (n_map > 0) ? map_bp[n_map - 1] + 1 : 0;

		/* Merge boundaries with gap < minCM */
		double merged_cM[24];
		int merged_bp[24];
		int nm = 0;
		merged_cM[nm] = cM0;
		merged_bp[nm] = bp0;
		nm++;
		for (int i = 0; i < nb; i++) {
			if (cMlist[i] - merged_cM[nm - 1] >= minCM) {
				merged_cM[nm] = cMlist[i];
				merged_bp[nm] = bplist[i];
				nm++;
			}
		}
		if (cMend - merged_cM[nm - 1] >= minCM || nm == 1) {
			merged_cM[nm] = cMend;
			merged_bp[nm] = bpend;
			nm++;
		}

		/* Convert bp boundaries to SNP indices; fill chrdivider */
		int nseg = 0;
		int start = 0;
		for (int i = 1; i < nm && nseg < MAXSEGPERCHR; i++) {
			int end_bp = merged_bp[i];
			int end_ix = 0;
			for (int j = 0; j < n_map; j++) {
				if (map_bp[j] < end_bp) end_ix = j + 1;
			}
			if (end_ix > n_map) end_ix = n_map;
			if (end_ix <= start) end_ix = start + 1;
			chrdivider[B][chr][nseg].start = start;
			chrdivider[B][chr][nseg].end = end_ix;
			chrdivider[B][chr][nseg].segment = 0;
			start = end_ix;
			nseg++;
			if (start >= n_map) break;
		}
		if (start < n_map && nseg < MAXSEGPERCHR) {
			chrdivider[B][chr][nseg].start = start;
			chrdivider[B][chr][nseg].end = n_map;
			chrdivider[B][chr][nseg].segment = 0;
			nseg++;
		}
		if (nseg >= MAXSEGPERCHR && start < n_map) {
			printf("EXIT(1): chr%d exceeded MAXSEGPERCHR=%d. Increase MAXSEGPERCHR.\n", chr, MAXSEGPERCHR);
			free(map_bp);
			return -1;
		}
		nbchrdivider[B][chr] = nseg;
		if (nseg > nbchrdivider[B][0]) nbchrdivider[B][0] = nseg;

		if (nseg > 15) {
			printf("[window setup] WARNING chr%d: nseg=%d (nm=%d) - with minCM=%d expect <=~15. Check hotspot/map consistency.\n", chr, nseg, nm, minCM);
		}
	}

	free(map_bp);
	return 0;
}

int distrigametic[12];	
int distrigametickeep[23];	
clock_t step1;

unsigned char * buffer;

/* Reads HAP data: text .haps or binary .bin.
* Text: one line per SNP: chr rsid position ref alt allele1_0 allele2_0 ... ; converts to packed 2-bit genotypes.
* Binary: raw packed buffer (2 bits per SNP, row-major). Uses fixed nbsnpperchr per chr.
* Writes into *out_buf (allocates and assigns). Sets nbsnpperchr[chr] and n_individuals (global). Exits on error. */
void readgenomelocal(int chr, const char *pathHapDir, unsigned char **out_buf)
{	FILE * fileresult;
	char pathfile[512];
	char rsid[128];

	if (useBinary) {
		/* Binary: {pathHapDir}{chr}.ped2.bin - raw packed genotypes. Uses fixed nbsnpperchr per chr. */
		static const int BINARY_NBSNPPERCHR[23] = {
			330005, 26229, 26210, 22209, 20690, 19027, 18418, 18367, 16283, 14990,
			16494, 15818, 16008, 11510, 10804, 10884, 12195, 11486, 10222, 9806,
			8985, 5227, 5882
		};
		int n_snps = (chr >= 1 && chr <= 22) ? BINARY_NBSNPPERCHR[chr] : 0;
		snprintf(pathfile, sizeof(pathfile), "%s%d.bin", pathHapDir, chr);
		if ((fileresult = fopen(pathfile, "rb")) == NULL) {
			printf("EXIT(1): Cannot open binary file %s. Check pathHap and useBinary.\n", pathfile);
			exit(1);
		}
		fseek(fileresult, 0, SEEK_END);
		uint64_t lSize = (uint64_t)ftell(fileresult);
		rewind(fileresult);
		if (n_snps <= 0 || n_snps > MAXSNP) {
			fclose(fileresult);
			printf("EXIT(1): chr%d invalid n_snps=%d for binary mode.\n", chr, n_snps);
			exit(1);
		}
		uint64_t bytes_per_indiv = (n_snps/4) + ((n_snps%4)>0);
		uint64_t n_indiv = lSize / bytes_per_indiv;
		if (n_indiv * bytes_per_indiv != lSize || n_indiv <= 0 || n_indiv > MAXPOP) {
			fclose(fileresult);
			printf("EXIT(1): chr%d binary file size %" PRIu64 " inconsistent (n_snps=%d). Expected multiple of %" PRIu64 ".\n", chr, lSize, n_snps, bytes_per_indiv);
			exit(1);
		}
		size_t buf_size = (size_t)lSize;
		*out_buf = (unsigned char*)malloc(buf_size);
		if (!*out_buf) {
			fclose(fileresult);
			printf("EXIT(1): chr%d memory allocation failed (%zu bytes).\n", chr, buf_size);
			exit(1);
		}
		size_t result = fread(*out_buf, 1, buf_size, fileresult);
		fclose(fileresult);
		if (result != buf_size) {
			free(*out_buf);
			*out_buf = NULL;
			printf("EXIT(1): chr%d binary read failed (read %zu, expected %zu).\n", chr, result, buf_size);
			exit(1);
		}
		nbsnpperchr[chr] = n_snps;
		n_individuals = (int)n_indiv;
		printf("chr%d: binary loaded n_individuals=%d n_snps=%d\n", chr, n_individuals, nbsnpperchr[chr]);
		return;
	}

	/* Text mode */
	snprintf(pathfile, sizeof(pathfile), "%s%d.haps", pathHapDir, chr);
	if ((fileresult = fopen(pathfile, "r")) == NULL) {
		printf("EXIT(1): Cannot open file %s. Check pathHap.\n", pathfile);
		fprintf(stderr, "ERROR: Cannot open file %s (file not found or permission denied). Check pathHap.\n", pathfile);
		exit(1);
	}
	{	char firstline[512];
		if (fgets(firstline, sizeof(firstline), fileresult)) {
			int truncated = (strlen(firstline) >= 256);
			firstline[256] = '\0';
			printf("  chr%d: first SNP start: %s%s\n", chr, firstline, truncated ? "..." : "");
		}
		rewind(fileresult);
	}

	printf("  chr%d: scanning file...\n", chr);
	fflush(stdout);
	/* First pass: count tokens on first line for n_individuals, then count lines for n_snps */
	int n_snps = 1;
	int n_indiv = 0;
	{	int tok_count = 0, in_token = 0, c;
		while ((c = fgetc(fileresult)) != EOF && c != '\n') {
			if (c == ' ' || c == '\t') in_token = 0;
			else { if (!in_token) { tok_count++; in_token = 1; } }
		}
		if (tok_count >= 6) n_indiv = (tok_count - 6) / 2;  /* 6 fixed: NA chr rsid pos ref alt */
		else if (tok_count >= 5) n_indiv = (tok_count - 5) / 2;  /* 5 fixed: chr rsid pos ref alt */
		else if (tok_count >= 4) n_indiv = (tok_count - 4) / 2;  /* 4 fixed: rsid pos ref alt (no chr) */
	}
	{	int c;
		while ((c = fgetc(fileresult)) != EOF)
		if (c == '\n') n_snps++;
	}
	if (n_indiv <= 0 || n_indiv > MAXPOP) {
		fclose(fileresult);
		printf("EXIT(1): chr%d invalid n_individuals=%d. Check HAP format.\n", chr, n_indiv);
		fprintf(stderr, "ERROR: chr%d - invalid n_individuals=%d (expected 1-%d). Check HAP format: chr rsid pos ref alt allele1_0 allele2_0 ...\n", chr, n_indiv, MAXPOP);
		exit(1);
	}
	if (n_snps <= 0 || n_snps > MAXSNP) {
		fclose(fileresult);
		printf("EXIT(1): chr%d invalid n_snps=%d.\n", chr, n_snps);
		fprintf(stderr, "ERROR: chr%d - invalid n_snps=%d (expected 1-%d).\n", chr, n_snps, MAXSNP);
		exit(1);
	}

	nbsnpperchr[chr] = n_snps;
	n_individuals = n_indiv;
	printf("chr%d: n_individuals=%d n_snps=%d\n", chr, n_individuals, nbsnpperchr[chr]);

	/* Allocate buffer: individuals x SNPs, 2 bits per SNP */
	uint64_t bytes_per_indiv = (n_snps/4) + ((n_snps%4)>0);
	size_t buf_size = (size_t)((uint64_t)n_indiv * bytes_per_indiv);
	*out_buf = (unsigned char*)calloc(1, buf_size);
	if (!*out_buf) {
		fclose(fileresult);
		printf("EXIT(1): chr%d memory allocation failed (%zu bytes).\n", chr, (size_t)buf_size);
		fprintf(stderr, "ERROR: chr%d - memory allocation failed (out of memory for %zu bytes).\n", chr, (size_t)buf_size);
		exit(1);
	}

	/* Second pass: read directly from file. Format: chr rsid pos ref alt a1 a2 a1 a2 ... */
	rewind(fileresult);
	long long cnt0 = 0, cnt1 = 0, cntCH = 0, cntHet = 0, cntRH = 0;
	int report_step = (n_snps >= 100000) ? 10000 : (n_snps >= 10000) ? 2000 : (n_snps >= 1000) ? 500 : 100;
	int last_pct = -1;
	for (int snp = 0; snp < n_snps; snp++)
	{	int dummy_chr, dummy_pos;
		char ref, alt;
		char first_col[32] = {0};
		if (fscanf(fileresult, " %31s", first_col) != 1) break;
		if (strcmp(first_col, "NA") == 0) {
			/* After NA: either "chr rsid pos ref alt" (chr=1-22) or "rsid rsid pos ref alt" (no chr) */
			char second_col[128] = {0};
			if (fscanf(fileresult, " %127s", second_col) != 1) break;
			dummy_chr = 0;
			if (second_col[0] >= '0' && second_col[0] <= '9')
				dummy_chr = atoi(second_col);
			else if (strncmp(second_col, "chr", 3) == 0)
				dummy_chr = atoi(second_col + 3);
			if (dummy_chr >= 1 && dummy_chr <= 22) {
				if (fscanf(fileresult, " %127s %d %c %c", rsid, &dummy_pos, &ref, &alt) < 4) break;
			} else {
				strncpy(rsid, second_col, 127); rsid[127] = '\0';
				if (fscanf(fileresult, " %*s %d %c %c", &dummy_pos, &ref, &alt) < 3) break;  /* skip duplicate rsid */
			}
		} else {
			int maybe_chr = 0;
			if (first_col[0] >= '0' && first_col[0] <= '9')
				maybe_chr = atoi(first_col);
			else if (strncmp(first_col, "chr", 3) == 0)
				maybe_chr = atoi(first_col + 3);
			if (maybe_chr >= 1 && maybe_chr <= 22) {
				/* Chr column present: chr rsid pos ref alt */
				if (fscanf(fileresult, " %127s %d %c %c", rsid, &dummy_pos, &ref, &alt) < 4) break;
			} else {
				/* No chr column: rsid pos ref alt (first_col is rsid) */
				strncpy(rsid, first_col, 127); rsid[127] = '\0';
				if (fscanf(fileresult, " %d %c %c", &dummy_pos, &ref, &alt) < 3) break;
			}
		}
		long long snp0 = 0, snp1 = 0, snpCH = 0, snpHet = 0, snpRH = 0;
		for (int indiv = 0; indiv < n_indiv; indiv++)
		{	float fa1, fa2;
			if (fscanf(fileresult, " %f %f", &fa1, &fa2) != 2) break;
			int a1 = (int)(fa1 + 0.5f), a2 = (int)(fa2 + 0.5f);
			if (a1 == 0) { cnt0++; snp0++; } else { cnt1++; snp1++; }
			if (a2 == 0) { cnt0++; snp0++; } else { cnt1++; snp1++; }
			int geno = (a1 == 0 && a2 == 0) ? 0 : (a1 == 0 && a2 == 1) ? 1 : (a1 == 1 && a2 == 0) ? 2 : (a1 == 1 && a2 == 1) ? 3 : 0;
			if (geno == 0) { cntCH++; snpCH++; } else if (geno == 1 || geno == 2) { cntHet++; snpHet++; } else { cntRH++; snpRH++; }
			size_t byte_idx = (size_t)indiv * bytes_per_indiv + (size_t)(snp/4);
			int shift = (snp % 4) * 2;
			(*out_buf)[byte_idx] |= (unsigned char)((geno & 3) << shift);
		}
		if ((snp + 1) % report_step == 0 || snp == 0 || snp == n_snps - 1) {
			int pct = (int)(100.0 * (snp + 1) / n_snps);
			if (pct != last_pct || snp == 0 || snp == n_snps - 1) {
				printf("  chr%d: %d/%d SNPs (%d%%)\r", chr, snp + 1, n_snps, pct);
				fflush(stdout);
				last_pct = pct;
			}
		}
	}
	printf("  chr%d: done %d SNPs, %lld alleles=0, %lld alleles=1 | CH=%lld Het=%lld RH=%lld\n", chr, n_snps, cnt0, cnt1, cntCH, cntHet, cntRH);
	fclose(fileresult);
}

int nbrelatpihat=0;

/* Compute allele frequencies and pi-hat between individual ID and all others.
* pi-hat formula: hat{pi}(f,i) = (1/m) * sum_k [ (x_kf - p_k)*(x_ki - p_k) / (2*p_k*(1-p_k)) ]
* where x = genotype count (0,1,2) of alternate allele, p_k = allele frequency at SNP k. */
int calculatepihat(int ID)
{	float seuilpihat=0.33;
	int nind = n_individuals > 0 ? n_individuals : MAXPOP;

	printf("\n[calculatepihat] Computing pi-hat for individual ID=%d (vs %d individuals)\n", ID, nind);

	/* Step 1: Compute allele frequency p_k for all SNPs (chr 1..22) */
	printf("[calculatepihat] Step 1: Computing allele frequencies p_k per SNP across chromosomes 1-22...\n");
	for(int chr=1;chr<23;chr++)
	{	int nsp = nbsnpperchr[chr];
		int nsp4 = (nsp/4) + ((nsp%4)>0);
		for(int snp=0;snp<nsp;snp++)
		{	int count_alt = 0, n_valid = 0;
			for(int j=0;j<nind;j++)
			{	int g = (genomes[chr] ? ((*((genomes[chr] + (unsigned long long)j*nsp4) + snp/4)>>((snp%4)*2))&3) : 3);
				/* g: 0=hom ref, 1,2=het, 3=hom alt (all valid). Count alt alleles: (g&1)+(g>>1) = 0,1,1,2 */
				count_alt += (g&1) + (g>>1);
				n_valid += 2;
			}
			double p = (n_valid > 0) ? (double)count_alt / (double)n_valid : 0.5;
			double maf = (p < 0.5) ? p : (1.0 - p);
			MAF[snp][chr] = (int)(maf * n_valid + 0.5);
			allele_freq[snp][chr] = p;  /* p_k for pi-hat */
		}
	}
	int hist[11] = {0};  /* MAF bins: [0,0.05], (0.05,0.1], ..., (0.45,0.5] */
	for(int chr=1;chr<23;chr++)
	for(int snp=0;snp<nbsnpperchr[chr];snp++) {
		double p = allele_freq[snp][chr];
		double maf = (p < 0.5) ? p : (1.0 - p);
		int b = (int)(maf * 20.0);  /* 0..10 for 0 to 0.5 */
		if (b > 10) b = 10;
		hist[b]++;
	}
	printf("\n--- Allele frequency (MAF) distribution ---\n");
	printf("MAF range       Count\n");
	for(int b=0;b<10;b++)
	printf("%.2f-%.2f    %8d\n", b*0.05, (b+1)*0.05, hist[b]);
	printf("0.50            %8d\n", hist[10]);
	printf("-------------------------------------------\n");
	printf("Allele frequencies computed for all chromosomes\n");

	/* Step 2: Compute pi-hat between ID and every other individual using formula */
	for(int relat=0;relat<nind;relat++) pihatagainstall[relat]=0;
	for(int relat=0;relat<nind;relat++) pihatagainstall2[relat]=0;

	for(int relat=0;relat<nind;relat++)
	{	double sum_term = 0.0;
		int m_count = 0;
		for(int chr=1;chr<23;chr++)
		{	int nsp = nbsnpperchr[chr];
			int nsp4 = (nsp/4) + ((nsp%4)>0);
			for(int snp=0;snp<nsp;snp++)
			{	int g0 = (*((genomes[chr] + (unsigned long long)ID*nsp4) + snp/4)>>((snp%4)*2))&3;
				int g1 = (*((genomes[chr] + (unsigned long long)relat*nsp4) + snp/4)>>((snp%4)*2))&3;
				double p_k = allele_freq[snp][chr];
				/* x = allele count (0,1,2): g=0->0, g=1,2->1, g=3->2 */
				double x_kf = (double)((g0&1) + (g0>>1));
				double x_ki = (double)((g1&1) + (g1>>1));
				double denom = 2.0 * p_k * (1.0 - p_k);
				if (denom <= 0.0) continue;
				sum_term += (x_kf - 2.0*p_k) * (x_ki - 2.0*p_k) / denom;
				m_count++;
			}
		}
		float pihat_val = 0.0f;
		if (m_count > 0)
		{	pihat_val = (float)(sum_term / (double)m_count);

		}
		pihatagainstall[relat] = pihat_val;
		/* pihatagainstall2 kept for compatibility (IBS0 proportion - now set to 0 if unused) */
		pihatagainstall2[relat] = 0.0f;
	}

	/* Sort by pi-hat to get best matches for bestpihatagainstall */
	float bestpihatagainstall[100];
	for(int i=0;i<100;i++) bestpihatagainstall[i]=-1;
	for(int i=0;i<100;i++) bestpihatagainstallID[i]=-1;
	nbrelatpihat = 0;
	for(int relat=0;relat<nind;relat++)
	{	float v = pihatagainstall[relat];  /* includes ID vs self = 1.0 */
		for(int i=0;i<100;i++)
		{	/* Insert if v is better than slot i, or slot is empty and v>=0 */
			if (v > bestpihatagainstall[i] )
			{	for(int j=99;j>i;j--) { bestpihatagainstall[j]=bestpihatagainstall[j-1]; bestpihatagainstallID[j]=bestpihatagainstallID[j-1]; }
				bestpihatagainstall[i]=v; bestpihatagainstallID[i]=relat;
				break;
			}
		}
	}
	for(int i=0;i<100 && bestpihatagainstallID[i]>=0;i++) nbrelatpihat++;

	printf("\n--- Top 100 highest pi-hat (most related individuals) ---\n");
	printf("%4s  %10s  %12s\n", "Rank", "Individual", "pi-hat");
	printf("----  ----------  ------------\n");
	for(int i=0;i<100 && bestpihatagainstallID[i]>=0;i++)
	printf("%4d  %10d  %12.6f\n", i+1, bestpihatagainstallID[i], bestpihatagainstall[i]);
	printf("--------------------------------------------------------\n\n");

	int compteur=0;
	do {compteur++;} while (compteur<100 && bestpihatagainstall[compteur]>seuilpihat);
	placefirttoconsider=compteur;
	IDbestpihat=(compteur<100 && bestpihatagainstallID[compteur]>=0)?bestpihatagainstallID[compteur]:-1;
	
	return 0;
}

#ifdef __GNUC__
__attribute__((noinline))
#endif
int acrosschrphasing(int ID,int IDp1loop,int IDp2loop,char pathresult[], const char *pathHapDir, const char *pathwindow,
	int (*phaseparent1tap)[NSNPPERCHR], int (*phaseguesedparent1tap)[NSNPPERCHR],
	double (*pihatagainstallchrMPphaseerror)[MAXPOP][200][2])
{
	printf("Start main phasing program\n");
	fflush(stdout);
	int nind = n_individuals > 0 ? n_individuals : MAXPOP;
	seuilpihat[0]=0.33;
	seuilpihat[1]=0.33;
	seuilpihat[2]=0.33;
	/************************************************************************************************ */
	/*									load genome and initiate focal and parents is any			  */
	/************************************************************************************************ */
	int ratio=2;
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)
	{	/* genomes[] already loaded by main(); readgenomelocal called only there */
		for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)		
		{	genomeoffpss[0][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			genomeoffpss[1][snp][chrtemp1]=(IDp1loop>=0)?((*((genomes[chrtemp1]+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3):3;
			genomeoffpss[2][snp][chrtemp1]=(IDp2loop>=0)?((*((genomes[chrtemp1]+(unsigned long long) IDp2loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3):3;
		};	
	};
	clock_t step2 = clock();
	float elapsed_secs11 = (float)(step2-step1);
	printf("\nTotal time loading focal and parents:%f cpu click so %f seconds\n",elapsed_secs11,elapsed_secs11/CLOCKS_PER_SEC );
	/************************************************************************************************ */
	/*									iniatial variables											  */
	/************************************************************************************************ */
	int64_t iter=0;
	char pathfile[300];
	int phase1[23];
	int guessphase1[23];
	int phase1guessphase1[23];
	int phase1guessphase2[23];
	int phase2[23];
	int guessphase2[23];
	int phase2guessphase2[23];
	int phase2guessphase1[23];
	int phaseguessright[23];
	int phaseguesswrong[23];
	int phasenotguessed[23];
	int wrong[23];
	int right[23];
	for (int chr = 0; chr < 23; chr++) {
		phaseguessright[chr] = 0;
		phaseguesswrong[chr] = 0;
		phasenotguessed[chr] = 0;
		phase1[chr] = 0;
		guessphase1[chr] = 0;
		phase1guessphase1[chr] = 0;
		phase1guessphase2[chr] = 0;
		phase2[chr] = 0;
		guessphase2[chr] = 0;
		phase2guessphase2[chr] = 0;
		phase2guessphase1[chr] = 0;
	}
 
	/* phaseparent1tap, phaseguesedparent1tap, pihatagainstallchrMPphaseerror passed from main() */
	int bestpur=0;
	double (*dataperchr)[110] = (double (*)[110])calloc(23, sizeof(double[110]));
	if (!dataperchr) { printf("RETURN(1): Failed to alloc dataperchr\n"); return 1; }
	int powerpihatDEGREE4=0;
	
	/* Segment boundaries (chrdivider, nbchrdivider) already set up in main() before HAP load */
	int breaknubercm=27;

	///////////******************************************************************************************************************************************************************************
	///////////                                    start main loop                                                                                                                    *
	///////////******************************************************************************************************************************************************************************			

	float seuiltocorrect=0.022;

	char (*segwithav)[NSNPPERCHR] = (char (*)[NSNPPERCHR])calloc(23, NSNPPERCHR);
	if (!segwithav) { printf("RETURN(1): Failed to alloc segwithav\n"); return 1; }
	double seuilpihatcorrection[6]={seuiltocorrect,seuiltocorrect,seuiltocorrect,seuiltocorrect,seuiltocorrect,seuiltocorrect};
	printf("there are %d individuals with a pihat higher than %f\n",placefirttoconsider,seuilpihat[0]);
///////////******************************************************************************************************************************************************************************
///////////                                 correct focal individual                                                                                                           *
///////////******************************************************************************************************************************************************************************			

	for(int relattocompare=0;relattocompare<20;relattocompare++)	
	{	int countsnp=0;
		for(int chrtemp1=1;chrtemp1<23;chrtemp1++)		
		{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchr[chrtemp1];chrdividerrun++)
			{	segwithav[chrtemp1][chrdividerrun]=0; 
			};
		};

		int IDrelattocompare=bestpihatagainstallID[placefirttoconsider+relattocompare];
		printf("Individual %d with a pihat of %f is the individual to consider with a rank of %d in the list of pihats \n",IDrelattocompare,pihatagainstall[IDrelattocompare],placefirttoconsider+relattocompare);
		if (pihatagainstall[IDrelattocompare]>seuiltocorrect )
		{		
			printf("Individual %d pihat of %f is higher than %f\n",IDrelattocompare,pihatagainstall[IDrelattocompare],seuiltocorrect);
//#pragma omp parallel for 	
			for(int chrtemp1=22;chrtemp1>0;chrtemp1--)		
			{	printf("Start correcting chr %d with %d nb of SNPs\n",chrtemp1,nbsnpperchr[chrtemp1]);
				int (*nbindivmatchseg)[4] = (int (*)[4])calloc(nind, sizeof(int[4]));
				if (!nbindivmatchseg) { printf("OMP thread: failed to alloc nbindivmatchseg\n"); continue; }
				for(int64_t relat=0;relat<(int64_t) nind;relat++)
				{	nbindivmatchseg[relat][0]=0;
					nbindivmatchseg[relat][1]=0;
					nbindivmatchseg[relat][2]=0;
					nbindivmatchseg[relat][3]=0;
				};
				int *nbsegmentofthislength = (int *)calloc(NSNPPERCHR + 2, sizeof(int));
				if (!nbsegmentofthislength) { printf("OMP thread: failed to alloc nbsegmentofthislength\n"); free(nbindivmatchseg); continue; }
				for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
				{	nbsegmentofthislength[snp]=0;
				};
				int typesegment=-1;
				int lasttypesegment=-1;
				int endlastsegment[4]={-1};
				double ratetoconsiderseg=0.0001;
				double ratetoconsidersegPE=ratetoconsiderseg*2; 
				ratetoconsiderseg=0.000012; 
				ratetoconsidersegPE=ratetoconsiderseg;  
				int phaseerrorpossible[4]={0};

				for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
				{//	printf("SNP	%d strat",snp);
					int snpvalue0=genomeoffpss[0][snp][chrtemp1];		
					int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
					int snpmod4=(((snp%4)*2));
					int snpdiv4=snp/4;
					for(int64_t relat=0;relat<(int64_t) nind;relat++) if (relat!=ID && pihatagainstall[relat]<seuilpihat[0])
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
				//	printf("SNP	%d check for all \n",snp);
					int relat=IDrelattocompare;
					if (nbindivmatchseg[relat][0]==0) phaseerrorpossible[0]=0;
					if (nbindivmatchseg[relat][1]==0) phaseerrorpossible[1]=0;
					if (nbindivmatchseg[relat][2]==0) phaseerrorpossible[2]=0;
					if (nbindivmatchseg[relat][3]==0) phaseerrorpossible[3]=0;

					int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
					int lenght;
					if (nbindivmatchseg[relat][0]>10 && 
					nbsegmentofthislength[nbindivmatchseg[relat][0]]<(phaseerrorpossible[0]?nind*ratetoconsidersegPE:nind*ratetoconsiderseg) && 
					nbindivmatchseg[relat][0]>nbindivmatchseg[relat][1] &&
					nbindivmatchseg[relat][0]>nbindivmatchseg[relat][2] &&
					nbindivmatchseg[relat][0]>nbindivmatchseg[relat][3] 
					)
					{	//printf("chr %d segment type 1 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][0], nbsegmentofthislength[nbindivmatchseg[relat][0]]);
						lenght=nbindivmatchseg[relat][0];
						typesegment=1;
						for(int spnrun=snp-nbindivmatchseg[relat][0];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
						endlastsegment[0]=snp;	
					} else if (nbindivmatchseg[relat][1]>10 && nbsegmentofthislength[nbindivmatchseg[relat][1]]<(phaseerrorpossible[1]?nind*ratetoconsidersegPE:nind*ratetoconsiderseg) &&
					nbindivmatchseg[relat][1]>nbindivmatchseg[relat][2] &&
					nbindivmatchseg[relat][1]>nbindivmatchseg[relat][3] 
					)
					{//	printf("chr %d segment type 3 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][1], nbsegmentofthislength[nbindivmatchseg[relat][1]]);
						lenght=nbindivmatchseg[relat][1];
						typesegment=3;
						endlastsegment[1]=snp;
						for(int spnrun=snp-nbindivmatchseg[relat][1];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
					} else if (nbindivmatchseg[relat][2]>10 && nbsegmentofthislength[nbindivmatchseg[relat][2]]<(phaseerrorpossible[2]?nind*ratetoconsidersegPE:nind*ratetoconsiderseg) &&
					nbindivmatchseg[relat][2]>nbindivmatchseg[relat][3] 
					)
					{	//printf("chr %d segment type 2 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][2], nbsegmentofthislength[nbindivmatchseg[relat][2]]);
						lenght=nbindivmatchseg[relat][2];
						typesegment=2;
						endlastsegment[2]=snp;
						for(int spnrun=snp-nbindivmatchseg[relat][2];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=2;

					} else if (nbindivmatchseg[relat][3]>10 && nbsegmentofthislength[nbindivmatchseg[relat][3]]<(phaseerrorpossible[3]?nind*ratetoconsidersegPE:nind*ratetoconsiderseg)
					)
					{	//printf("chr %d segment type 4 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][3], nbsegmentofthislength[nbindivmatchseg[relat][3]]);
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
				//	printf("SNP	%d check for indiv %d\n",snp,IDrelattocompare);
					
					if (typesegment>-1 && lasttypesegment>-1 && (typesegment&1)!=(lasttypesegment&1))
					{	int end;
						if (typesegment==1) end=endlastsegment[0];
						if (typesegment==2) end=endlastsegment[2];
						if (typesegment==3) end=endlastsegment[1];
						if (typesegment==4) end=endlastsegment[3];
					//	printf("change from %d to %d\n",(end+snp-lenght)/2,nbsnpperchr[chrtemp1]);
						//	nbphasecorrect++;
						for(int snprun=(snp-lenght);snprun<nbsnpperchr[chrtemp1];snprun++)		
						{//	nbsnpcorrect++;
							genomeoffpss[0][snprun][chrtemp1]=((genomeoffpss[0][snprun][chrtemp1]&1)<<1)+(genomeoffpss[0][snprun][chrtemp1]>>1);
						};
						snp=end;
						for(int snprun=end;snprun<nbsnpperchr[chrtemp1];snprun++)		
						{	segwithav[chrtemp1][snprun]=0;
						};

						typesegment=-1;
						for(int64_t relat=0;relat<(int64_t) nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] )
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
				printf("End correcting chr %d\n",chrtemp1);
				free(nbsegmentofthislength);
				free(nbindivmatchseg);
				
				
			};
			printf("Correcting with individual %d over\n",IDrelattocompare);
		};
	};
	printf("Correcting over\n");
	
//	free(segwithav);
	int breaknubercmloop=breaknubercm;
	
///////////******************************************************************************************************************************************************************************
///////////                                check parents                                                                                                     *
///////////******************************************************************************************************************************************************************************			

	for(int chr=0;chr<23;chr++)
	{	phaseguessright[chr]=0;
		phaseguesswrong[chr]=0;
		phasenotguessed[chr]=0;
		phase1[chr]=0;
		guessphase1[chr]=0;
		phase1guessphase1[chr]=0;
		phase1guessphase2[chr]=0;
		phase2[chr]=0;
		guessphase2[chr]=0;
		phase2guessphase2[chr]=0;
		phase2guessphase1[chr]=0;
	};

	for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
	{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
		{	phaseparent1tap[chrtemp1][snp]=0;
			phaseguesedparent1tap[chrtemp1][snp]=0;
		};
	};

	breaknubercm=breaknubercmloop;

//	char (*segwithav)[NSNPPERCHR] = (char (*)[NSNPPERCHR])calloc(23, NSNPPERCHR);
	if (!segwithav) { printf("RETURN(1): Failed to alloc segwithav\n"); return 1; } else {printf("array segwithav successfully alocated\n");};
	int countsnp=0;
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchr[chrtemp1];chrdividerrun++)
		{	segwithav[chrtemp1][chrdividerrun]=0; 
		};
	};

	if (IDp1loop>=0) for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
	{	printf("Checking GT for chr %d\n",chrtemp1);
		phaseguessright[chrtemp1]=0;
		phaseguesswrong[chrtemp1]=0;	
		wrong[chrtemp1]=0;
		right[chrtemp1]=0;
		int phaseparent1=-1;
		int phaseparent1temp=-1;
		int nbsnperror=0;
		int64_t nbonesnp=1;//rand();
		int64_t nbzerosnp=0;//rand();
		//	printf("%d %d %d\n",6,phaseguessright[6],phaseguesswrong[6]);
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
				//	 printf("snp %d: %d %d %d\n",snp,phaseparent1temp,phaseparent1,phase);
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
				{	phase1[ratio]++;
					guessphase1[ratio]++;
					phase1guessphase1[ratio]++;
					phaseguessright[chrtemp1]++;
					phaseguesedparent1tap[chrtemp1][snp]=1;
					right[chrtemp1]++;
				}
				else if (nbonesnp>1*nbzerosnp && phaseparent1==2)
				{	phase2[ratio]++;
					guessphase1[ratio]++;
					phase2guessphase1[ratio]++;
					phaseguesswrong[chrtemp1]++;
					phaseguesedparent1tap[chrtemp1][snp]=1;
					wrong[chrtemp1]++;
				}
				else if (nbonesnp*1<nbzerosnp && phaseparent1==2)
				{	phase2[ratio]++;
					guessphase2[ratio]++;
					phase2guessphase2[ratio]++;
					phaseguessright[chrtemp1]++;
					phaseguesedparent1tap[chrtemp1][snp]=0;
					right[chrtemp1]++;
				}
				else if (nbonesnp*1<nbzerosnp && phaseparent1==1)
				{	phase1[ratio]++;
					guessphase2[ratio]++;
					phase1guessphase2[ratio]++;
					phaseguesswrong[chrtemp1]++;
					phaseguesedparent1tap[chrtemp1][snp]=0;
					wrong[chrtemp1]++;
				}
				else 
				{	phasenotguessed[ratio]++;
				};			
			} else 
			{
			};
		};
		bestpur=bestpur+(wrong[chrtemp1]>right[chrtemp1]?wrong[chrtemp1]:right[chrtemp1]);
		printf("on chr %d from parent nb snp on each phase: %d %d\n",chrtemp1,phaseguessright[chrtemp1],phaseguesswrong[chrtemp1]);
		dataperchr[chrtemp1][9]=(phaseguessright[chrtemp1]>phaseguesswrong[chrtemp1]?phaseguessright[chrtemp1]:phaseguesswrong[chrtemp1])/(1+phaseguessright[chrtemp1]+phaseguesswrong[chrtemp1]);
	};
	int64_t purcentage=0;
	ratio=2;
	purcentage=(int64_t) 10000*
		(1+(phaseguessright[ratio]>phaseguesswrong[ratio]?
		phaseguessright[ratio]:
		phaseguesswrong[ratio]))/ 
		(phaseguesswrong[ratio]+phaseguessright[ratio]+1.0);
//	printf("FROM PARENT pourcentage before %d %f\n",purcentage,(float) bestpur/330005); 

	int nbphaseright=0;
	int nbphasewrong=0;
	int tabnbphaseright[23][25];
	int tabnbphasewrong[23][25];
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	//for(int  snpphaseerror=0;snpphaseerror<nbsnpperchr[chrtemp1];snpphaseerror=snpphaseerror+100)		
		for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
		{	tabnbphaseright[chrtemp1][chrdividerrun]=0;
			tabnbphasewrong[chrtemp1][chrdividerrun]=0;
		};
	};
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	//for(int  snpphaseerror=0;snpphaseerror<nbsnpperchr[chrtemp1];snpphaseerror=snpphaseerror+100)		
		for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
		{	int nbphaseright=0;
			int nbphasewrong=0;
			for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
			{	if (1==phaseparent1tap[chrtemp1][snp]-1) nbphaseright++; else nbphasewrong++;
			};
			printf("on chr %d window number %d nb SNPs in phase 1 %d nb SNPs in phase 2 %d\n",chrtemp1,chrdividerrun,nbphaseright,nbphasewrong);
			tabnbphaseright[chrtemp1][chrdividerrun]=nbphaseright;
			tabnbphasewrong[chrtemp1][chrdividerrun]=nbphasewrong;
		};	
	}	
	
///////////******************************************************************************************************************************************************************************
///////////                              start phasing by calculating lambda                                                                                            *
///////////******************************************************************************************************************************************************************************			



#pragma omp parallel for					
	for(int relat=0;relat<nind;relat++)// if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
	{	pihatagainstallchrMPphaseerror[0][relat][0][0]=0;
		pihatagainstallchrMPphaseerror[0][relat][0][1]=0;
		for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][0]*4;chrdividerrun++)
		{	pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][0]=0;
			pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][1]=0;
		};
	}
	int hetsnp=0;
	double exponentinterchr=1;
	

	float firstexponent=2;//2+nbdivisionDEGREE4 ;//4.5+.5*nbdivision;//+exponentcorp*.1; DONE

	///exponent=1.5;

	int powersnppihat=5*2;
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	//for(int  snpphaseerror=0;snpphaseerror<nbsnpperchr[chrtemp1];snpphaseerror=snpphaseerror+100)		
		{	hetsnp=0;
			//chrdivider=nbsnpperchr[chrtemp1]/sizebin;
#pragma omp parallel for 				
			for(int relat=0;relat<nind;relat++)// if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
			{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1]*4;chrdividerrun++)
				{	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][0]=0;
					pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][1]=0;
				};
				//	pihatagainstallchrMPphaseerror[chrtemp1][1][relat][0]=0;
			};
			//	int  chrtemp2=chrtemp1;
			for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
			{	double (*temp)[4] = (double (*)[4])calloc(nind, sizeof(double[4]));
				if (!temp) { printf("Failed to alloc temp\n"); continue; }
				for(int64_t relat=0;relat<(int64_t) nind;relat++) 
				{	temp[relat][0]=0;
					temp[relat][1]=0;
					temp[relat][2]=0;
					temp[relat][3]=0;
				};
				for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
				{	int snpvalue0=genomeoffpss[0][snp][chrtemp1];
					if (snpvalue0!=0 && snpvalue0!=3 )//&& phaseparent1tap[chrtemp1][snp]!=-1)
					{//	if (snp<snpphaseerror) snpvalue0=3-snpvalue0;
						hetsnp++;
						double maffloat=(n_individuals>0)?(1.0*MAF[snp][chrtemp1]/(2*n_individuals)):(1.0*MAF[snp][chrtemp1]/nind/2);
						//if (phaseguesedparent1tap[chrtemp1][snp]==0) printf("STOP\n");
						int parent0indiv0=(phaseguesedparent1tap[chrtemp1][snp]==1)?(snpvalue0>>1):(snpvalue0&1);
						int parent1indiv0=(phaseguesedparent1tap[chrtemp1][snp]==1)?(snpvalue0&1):(snpvalue0>>1);
						//	printf("%d %d\n",phaseguesedparent1tap[chrtemp1][snp],snpvalue0);
						double pcontribu10=((parent0indiv0)-1.0*maffloat);
						double pcontribu11=((parent1indiv0)-1.0*maffloat);
						double maffloatdiviseur=(maffloat)*(2-maffloat*2)/3*2;
						double pconttab[4][2];
						pconttab[0][0]=pow(((0)-maffloat)*pcontribu10,1.0/powersnppihat);
						pconttab[0][1]=pow(((0)-maffloat)*pcontribu11,1.0/powersnppihat);
						pconttab[1][0]=pow(((1)-maffloat)*pcontribu10,1.0/powersnppihat);
						pconttab[1][1]=pow(((1)-maffloat)*pcontribu11,1.0/powersnppihat); 
						int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
						int snpmod4=(((snp%4)*2));
						int snpdiv4=snp/4;
#pragma omp parallel for
						for(int64_t relat=0;relat<(int64_t) nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.009 )
						{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
							if (pconttab[snpvalue1&1][0]>0) temp[relat][0]=temp[relat][0]+(pconttab[snpvalue1&1][0]);
							else temp[relat][0]=0;
							if (temp[relat][0]<0) temp[relat][0]=0;
							else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]<temp[relat][0]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]=temp[relat][0];
							if (pconttab[snpvalue1>>1][0]>0) temp[relat][1]=temp[relat][1]+(pconttab[snpvalue1>>1][0]);
							else temp[relat][1]=0;
							if (temp[relat][1]<0) temp[relat][1]=0;
							else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]<temp[relat][1]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]=temp[relat][1];
							if (pconttab[snpvalue1&1][1]>0) temp[relat][2]=temp[relat][2]+(pconttab[snpvalue1&1][1]);
							else temp[relat][2]=0;
							if (temp[relat][2]<0) temp[relat][2]=0;
							else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]<temp[relat][2]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]=temp[relat][2];
							if (pconttab[snpvalue1>>1][1]>0) temp[relat][3]=temp[relat][3]+(pconttab[snpvalue1>>1][1]);
							else temp[relat][3]=0;
							if (temp[relat][3]<0) temp[relat][3]=0;
							else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]<temp[relat][3]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]=temp[relat][3];

						};
					};
				};
				for(int relat=0;relat<nind;relat++) 
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
				
				free(temp);
			};
		};	
	};
	///////////******************************************************************************************************************************************************************************
///////////                              calculating firt correlation                                                                                    *
///////////******************************************************************************************************************************************************************************			


	float thirdindice=3;
	thirdindice=1;
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
				for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
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
	double (*allcor)[MAXNBDIVISOR][23][MAXNBDIVISOR] = (double (*)[MAXNBDIVISOR][23][MAXNBDIVISOR])calloc(23, sizeof(double[MAXNBDIVISOR][23][MAXNBDIVISOR]));
	double (*allseg)[MAXNBDIVISOR][23][MAXNBDIVISOR] = (double (*)[MAXNBDIVISOR][23][MAXNBDIVISOR])calloc(23, sizeof(double[MAXNBDIVISOR][23][MAXNBDIVISOR]));
	if (!allcor || !allseg) { printf("RETURN(1): Failed to alloc allcor/allseg\n"); free(dataperchr); return 1; }
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
	

	float seuil1=0.5;

	float penalty=0.75;
#pragma omp parallel for
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
		{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
			{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2)
				{	double sumproduct; 
					sumproduct=0;
					double nbelem=0;
					for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
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
					for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
					{	nbelem++;
						sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
						pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
					};
					double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
					sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
					((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
					sumproduct=0;
					nbelem=0;
					for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
					{	nbelem++;
						sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent)*
						pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
					};
					double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
					sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
					((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
					sumproduct=0;
					nbelem=0;
					for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
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
						/*printf("%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f  \n",chrtemp1,chrtemp2,chrdividerrun,chrdividerrun2,cor1,cor2,cor3,cor4,
							sumproduct,nbelem,
							sumallsquare[chrtemp1][chrdividerrun][0][0],sumallsquare[chrtemp2][chrdividerrun2][0][0],
							sumall[chrtemp1][chrdividerrun][0][0],sumall[chrtemp2][chrdividerrun2][0][0],
							sumproduct1,allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2],corglobmax);*/
					};											
				};
			};
		};
	};
///////////******************************************************************************************************************************************************************************
///////////                              start merging process                                                                            *
///////////******************************************************************************************************************************************************************************			
	printf("\nStart Merging process\n\n");
	//printf("MAX %d %d %d %d %f\n",chrmax1,dividmax1,chrmax2,dividmax2,corglobmax);
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
	float mergingexponent[2];//nbdivisionDEGREE4
	mergingexponent[0]=2;

	mergingexponent[0]=3;//TRY

	mergingexponent[1]=5; 
	int nbpair=0;
	do 
	{	group[chrmax2][dividmax2]=chrmax1+dividmax1*23;
		printf("Merging chr %d window %d to chr %d window %d based on lambda %f \n",chrmax2,dividmax2,chrmax1,dividmax1,corglobmax);
		printf ("merging info: %d %d\n",havemerged[chrmax1][dividmax1],havemerged[chrmax2][dividmax2]);
		tappointdec[nbmerge].nbgroup1=nbingroup[chrmax1][dividmax1];
		tappointdec[nbmerge].nbgroup2=nbingroup[chrmax2][dividmax2];
		tappointdec[nbmerge].cor=fabs(corglobmax);

		nbmerge++;

		nbingroup[chrmax1][dividmax1]=nbingroup[chrmax1][dividmax1]+nbingroup[chrmax2][dividmax2];
		if (corglobmax>0)
		{	havemerged[chrmax2][dividmax2]=1;
			


#pragma omp parallel for		
			for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
			{	//printf("%f %f %f %f\n",pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2],
				//		pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4],pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2]);
				for(int value=0;value<2;value++)
				{	double p1=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
					double p2=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
					double p3=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]>0?1:-1);
					double p4=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]>0?1:-1);

					pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]=pow(fabs(p1+p3),1.0/mergingexponent[value])*(p1+p3>0?1:-1);
					pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]=pow(fabs(p2+p4),1.0/mergingexponent[value])*(p2+p4>0?1:-1);
					//	printf("%f %f\n",pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2]);
				};
			}

		} else
		{	havemerged[chrmax2][dividmax2]=-1;


#pragma omp parallel for		
			for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
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
///////////******************************************************************************************************************************************************************************
///////////                           recalculating correlation                                                                        *
///////////******************************************************************************************************************************************************************************			

		sumall[chrmax1][dividmax1][0][0]=0;
		sumall[chrmax1][dividmax1][1][0]=0;
		sumallsquare[chrmax1][dividmax1][0][0]=0;
		sumallsquare[chrmax1][dividmax1][1][0]=0;
		sumall[chrmax1][dividmax1][0][1]=0;
		sumall[chrmax1][dividmax1][1][1]=0;
		sumallsquare[chrmax1][dividmax1][0][1]=0;
		sumallsquare[chrmax1][dividmax1][1][1]=0;
		for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
		{	for(int value=0;value<2;value++)
			{	sumall[chrmax1][dividmax1][0][value]=sumall[chrmax1][dividmax1][0][value]+
				pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*(0.3+0.7))*
				(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
				sumall[chrmax1][dividmax1][1][value]=sumall[chrmax1][dividmax1][1][value]+
				pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*(0.3+0.7))*
				(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
				sumallsquare[chrmax1][dividmax1][0][value]=sumallsquare[chrmax1][dividmax1][0][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*2);
				sumallsquare[chrmax1][dividmax1][1][value]=sumallsquare[chrmax1][dividmax1][1][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*2);
				//	printf("%f %f\n",sumall[chrmax1][dividmax1][0],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4]);
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
				for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
				{	nbelem++;
					sumproduct=sumproduct+
					pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
					pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
				};
				double  cor1=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
				sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
				((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
				sumproduct=0;
				//	if (fabs(cor1)<fabs(cor11)*coefforcor) cor1=cor11;
				sumproduct=0;
				nbelem=0;
				for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
				{	nbelem++;
					sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent*(0.3+0.7))*
					pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
				};
				double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
				sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
				((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
				sumproduct=0;
				//		if (fabs(cor2)<fabs(cor22)*coefforcor) cor2=cor22;

				sumproduct=0;
				nbelem=0;
				for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
				{	nbelem++;
					sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
					pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
				};
				double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
				sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
				((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));

				sumproduct=0; 

				//		if (fabs(cor3)<fabs(cor33)*coefforcor) cor3=cor33;				
				sumproduct=0;
				nbelem=0;
				for(int relat=0;relat<nind;relat++) if (pihatagainstall[relat]<seuilpihat[0] ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
				{	nbelem++;
					sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
					pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent*(0.3+0.7));
				};
				double  cor4=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
				sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
				((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
				sumproduct=0;
				//		if (fabs(cor4)<fabs(cor44)*coefforcor) cor4=cor44;				

				float exponantnumer=1;
				exponantnumer=0.5;
				double corglob=(pow(fabs(cor1),thirdindice-1)*cor1 -
				pow(fabs(cor2),thirdindice-1)*cor2+ 
				pow(fabs(cor3),thirdindice-1)*cor3-
				pow(fabs(cor4),thirdindice-1)*cor4)
				*pow((nbingroup[chrtemp1][chrdividerrun]+nbingroup[chrtemp2][chrdividerrun2]),exponantnumer);
				if (highestnew<fabs(corglob))
				{	highestnew=fabs(corglob);
					//printf("%f",corglob);
					//printf("NB1 %d MB2 %d \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2]);  
				};
				allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=corglob;										
				allcor[chrtemp2][chrdividerrun2][chrtemp1][chrdividerrun]=corglob;										
			};
		};
		nbpair=0;
///////////******************************************************************************************************************************************************************************
///////////                          evaluating if merging create phase errore                                                          *
///////////******************************************************************************************************************************************************************************			

		for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
		{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++) if (havemerged[chrtemp1][chrdividerrun]==0)
			{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
				{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (havemerged[chrtemp2][chrdividerrun2]==0 && (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2))
					{	printf("Remaining pairs: %d %d %d %d %d %d \n",chrtemp1,chrdividerrun,chrtemp2,chrdividerrun2,havemerged[chrtemp1][chrdividerrun],havemerged[chrtemp2][chrdividerrun2]);
						nbpair++;
						double corglob=allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2];
						printf("%f %f\n",corglob,corglobmax);
						if (fabs(corglob)>=fabs(corglobmax))
						{//	printf("%d %d %d %d %f\n",chrtemp1,chrtemp2,chrdividerrun,chrdividerrun2,corglob	);

						//	printf("NB1 %d MB2 %d c %f \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2],corglob);  
							if (nbingroup[chrtemp1][chrdividerrun]==1 || nbingroup[chrtemp2][chrdividerrun2]==1)
							{	int chrone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrtemp2:chrtemp1;
								int chrdividerrunone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrdividerrun2:chrdividerrun;
								int chrnotone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrtemp2:chrtemp1;
								int chrdividerrunnotone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrdividerrun2:chrdividerrun;

							//	printf("1 is %d %d\n",chrnotone,chrdividerrunnotone);
								int numberattached=0;
								for(int chrdividerrunrun=0;chrdividerrunrun<nbchrdivider[breaknubercm][chrone];chrdividerrunrun++) 
								{//	printf("checking window %d\n",chrdividerrunrun);

									if (chrdividerrunrun==chrdividerrunone)
									{	//printf("window n %d is consider window\n", chrdividerrunrun);
									} else
									{	if (havemerged[chrone][chrdividerrunrun]==0)
										{//	printf("window n%d hasn t merge\n", chrdividerrunrun);
										} else 
										{	int hrdividerloop=chrdividerrunrun;	
											int chrloop=chrone;
											int phaseknown=havemerged[chrloop][hrdividerloop];
											while (iter<1000 && havemerged[chrloop][hrdividerloop]!=0)
											{	//printf("%d %d %d %d\n", phaseknown,chr,chrdivi,havemerged[chr][chrdivi]);	

												iter++;
												phaseknown=phaseknown*(havemerged[chrloop][hrdividerloop]!=0?havemerged[chrloop][hrdividerloop]:1);
												int savechr=chrloop;
												chrloop=group[chrloop][hrdividerloop]%23;
												hrdividerloop=group[savechr][hrdividerloop]/23;

											};
											if 	(chrloop==chrnotone && hrdividerloop==chrdividerrunnotone)
											{	printf("chr %d window n %d is attached to chr %d window %d with phaseknown %d\n",chrloop, chrdividerrunrun,chrnotone,chrdividerrunnotone,phaseknown);
												numberattached=numberattached+phaseknown;
											} else 
											{	printf("chr %d window n %d is not attached to chr %d window %d\n", chrloop,chrdividerrunrun,chrnotone,chrdividerrunnotone);
											};																		
										}; 
									};
								};
								if (corglobmax*numberattached<0) corglobmax=corglobmax*penalty;
								printf("%f %f %d\n",corglob,corglobmax,numberattached);
						
								if (fabs(corglob)>=fabs(corglobmax))
								{	corglobmax=corglob;
									chrmax1=chrtemp1;						
									chrmax2=chrtemp2;						
									dividmax1=chrdividerrun;						
									dividmax2=chrdividerrun2;
								}
							} else 	
							{	corglobmax=corglob;
								chrmax1=chrtemp1;						
								chrmax2=chrtemp2;						
								dividmax1=chrdividerrun;						
								dividmax2=chrdividerrun2;
							};
						};	
					};
				};
			};
		};
		printf("Nb pairs: %d\n",nbpair);
		for(int  chr=1;chr<23;chr++)		
		{	//printf("%d ",chr);
			for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chr];chrdividerrun++)
			{//	printf("%d ", havemerged[chr][chrdividerrun]+1);
			};
			//printf("\n");
		};
		for(int  chr=1;chr<23;chr++)		
		{	//printf("%d ",chr);
			for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chr];chrdividerrun++)
			{	//printf("%d %d ", (group[chr][chrdividerrun]%23)+9,10+(group[chr][chrdividerrun]/23));
			};
			//	printf("\n");
		};

		//	exit(0);
	} while (nbpair>0 && nbmerge<22*19-1);

	nbphaseright=0;
	nbphasewrong=0;
///////////******************************************************************************************************************************************************************************
///////////                          evaluating merging                                             *
///////////******************************************************************************************************************************************************************************			
	int phase1perchr[23]={0};
	int phase2perchr[23]={0};
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		 
	{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
		{	int chr=chrtemp1;
			int chrdivi=chrdividerrun;
			int phaseknown=1;
			int iter=0;
			while (iter<1000 && havemerged[chr][chrdivi]!=0)
			{	//printf("%d %d %d %d\n", phaseknown,chr,chrdivi,havemerged[chr][chrdivi]);	
				iter++;
				phaseknown=phaseknown*(havemerged[chr][chrdivi]!=0?havemerged[chr][chrdivi]:1);
				int savechr=chr;
				chr=group[chr][chrdivi]%23;
				chrdivi=group[savechr][chrdivi]/23;

			} 
			if (iter>900) { printf("EXIT(0): iter>900 (safety limit)\n"); exit(0); }
		//	printf("pahse chr %d div %d is %d\n",chrtemp1,chrdividerrun,phaseknown);
			chrdivider[breaknubercm][chrtemp1][chrdividerrun].phasing=phaseknown;
			
			for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
			{	if ((phaseknown+3)/2==2) 
				{	if (genomeoffpss[0][snp][chrtemp1]==1) genomeoffpss[0][snp][chrtemp1]=2;
					else if (genomeoffpss[0][snp][chrtemp1]==2) genomeoffpss[0][snp][chrtemp1]=1;
				};
				if ((phaseknown+3)/2==phaseparent1tap[chrtemp1][snp]) 
				{	nbphaseright++;
					phase1perchr[chrtemp1]++;
					chrdivider[breaknubercm][0][0].nbright++;
					chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbright++;
				} else 
				{	nbphasewrong++;
					phase2perchr[chrtemp1]++;
					chrdivider[breaknubercm][0][0].nbwrong++;
					chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbwrong++;

				};

			};
			//printf("chr %d div %d nbphaseright %d nbphasewrong %d\n",chrtemp1,chrdividerrun,nbphaseright,nbphasewrong);
		};
		if (IDp1loop>-1) printf("chr %d nb SNPs in phase 1: %d nb SNPs in phase 2: %d ratio: %f \n",chrtemp1,phase1perchr[chrtemp1],phase2perchr[chrtemp1],
			(phase1perchr[chrtemp1]>phase2perchr[chrtemp1]?phase1perchr[chrtemp1]:phase2perchr[chrtemp1])/(1.0+phase1perchr[chrtemp1]+phase2perchr[chrtemp1]));
	};
	printf("\nResult of merging process\n\n");

	float phasingscore=(int64_t) 100.0*(nbphaseright>nbphasewrong?nbphaseright:nbphasewrong)/(nbphaseright+nbphasewrong>0?0.0+nbphaseright+nbphasewrong:1.0);	
	if (IDp1loop>-1) printf("For all chr nb SNPs in phase 1: %d nb SNPs in phase 2: %d leading to an across chr phasing score of %f\n",nbphaseright,nbphasewrong, phasingscore);					
	free(allcor);
	free(allseg);

	free(segwithav);
	strcpy(pathfile,pathresult);
	/*	sprintf(number, "%d", numtrio);
	strcat(pathfile, number);
	strcat(pathfile, "/");
	sprintf(number, "%d", ID);
	strcat(pathfile, number);*/
	///////////******************************************************************************************************************************************************************************
///////////                          save result                                                         *
///////////******************************************************************************************************************************************************************************			
	free(dataperchr);
	printf("Saving genome of focal with haplotype 1 detected from parent 1 and haplotype 2 detected from parent 2 to file: %s\n", pathfile);
	FILE * filegenomeoffpss;
	if ((filegenomeoffpss = fopen(pathfile, "w")) == NULL) 
		{	printf("Cannot create output file: %s\n", pathfile);
			return (1);
		};
		fprintf(filegenomeoffpss, "#chr snp allelevalue\n");
		for(int chr=1; chr<23; chr++)
	{	for(int snp=0; snp<nbsnpperchr[chr]; snp++)
		{	fprintf(filegenomeoffpss, "%d %d %d\n", chr, snp, genomeoffpss[0][snp][chr]);
		};
	};
	fclose(filegenomeoffpss);
	
	printf("\nAcrosschrphasing completed successfully.\n\n");
	
	return (0);
}

static void parse_arg(const char *arg, char *pathresult)
{
	const char *eq = strchr(arg, '=');
	if (!eq || eq == arg) return;
	const char *val = eq + 1;

	if (strncmp(arg, "pathresult=", 11) == 0 || strncmp(arg, "pathResult=", 11) == 0) {
		strncpy(pathResult, val, sizeof(pathResult)-1); pathResult[sizeof(pathResult)-1]=0;
		strncpy(pathresult, val, 511); pathresult[511]=0;
		return;
	}
	if (strncmp(arg, "pathHap=", 8) == 0) { strncpy(pathHap, val, sizeof(pathHap)-1); pathHap[sizeof(pathHap)-1]=0; return; }
	if (strncmp(arg, "pathwindow=", 11) == 0) { strncpy(pathwindow, val, sizeof(pathwindow)-1); pathwindow[sizeof(pathwindow)-1]=0; return; }
	if (strncmp(arg, "pathHotspots=", 13) == 0) { strncpy(pathHotspots, val, sizeof(pathHotspots)-1); pathHotspots[sizeof(pathHotspots)-1]=0; return; }
	if (strncmp(arg, "pathMap=", 8) == 0) { strncpy(pathMap, val, sizeof(pathMap)-1); pathMap[sizeof(pathMap)-1]=0; return; }
	if (strncmp(arg, "minCM=", 6) == 0) { minCM = atoi(val); if (minCM <= 0) minCM = 22; return; }
	if (strncmp(arg, "useBinary=", 10) == 0) { useBinary = (atoi(val) != 0); return; }
	if (strncmp(arg, "posOffspring=", 13) == 0) { posOffspring = atoi(val); return; }
	if (strncmp(arg, "posParent1=", 11) == 0) { posParent1 = atoi(val); return; }
	if (strncmp(arg, "posParent2=", 11) == 0) { posParent2 = atoi(val); return; }
}


int main(int argc, char *argv[])
{
	printf("BUILD: heap-alloc fix 2025-02-23\n");
	fflush(stdout);
	clock_t begin0 = clock();
	sumtotindivcount = MAXPOP;

	static const char DEFAULT_PATH[] = "./";

	char pathresult[512];
	strcpy(pathresult, DEFAULT_PATH);

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "help") == 0 || strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("Usage: %s [key=value ...]\n", argv[0]);
			printf("Four path inputs:\n");
			printf("  pathHap=...       Path for .haps files (e.g. data/chr -> chr1.haps)\n");
			printf("  pathHotspots=...  Path to hotspot file (Chr cM cMperMb bp). If set, must exist.\n");
			printf("  pathMap=...       Path prefix for map files (e.g. /path/SNPCHR -> SNPCHR1.map). Col3=cM Col4=bp.\n");
			printf("  pathResult=...    Output directory\n");
			printf("Logic: pathwindow= for pre-computed bounds; OR pathHotspots= + pathMap= for hotspot-based segments (both required).\n");
			printf("  useBinary=1       Read .bin (raw packed) instead of text .haps. Uses fixed n_snps per chr.\n");
			printf("Other: posOffspring=N posParent1=N posParent2=N minCM=N (default 22)\n");
			printf("Example (binary): pathHap=/path/pedfile123chr useBinary=1 pathwindow=... pathResult=./out posOffspring=0\n");
			printf("Example (hotspot): pathHap=data/chr pathHotspots=hotspots.txt pathMap=map/SNPCHR pathResult=./out posOffspring=42\n");
			printf("Example (window): pathHap=data/chr pathwindow=data/recombinaisonWindows.txt pathResult=./out posOffspring=42\n");
			printf("RETURN(0): help printed.\n");
			return 0;
		}
		parse_arg(argv[i], pathresult);
	}
	if (pathResult[0] == '\0') strncpy(pathResult, pathresult, sizeof(pathResult)-1);

	printf("pathHap=%s pathHotspots=%s pathMap=%s pathResult=%s pathwindow=%s useBinary=%d minCM=%d posOffspring=%d posParent1=%d posParent2=%d\n",
		pathHap, pathHotspots, pathMap, pathResult, pathwindow, useBinary, minCM, posOffspring, posParent1, posParent2);


	srand(0);
	if (pathHap[0] == '\0' || posOffspring < 0) {
		printf("RETURN(1): pathHap= and posOffspring= are required.\n");
		printf("Run with --help for usage.\n");
		return 1;
	}
	if (!pathwindow[0] && !(pathHotspots[0] && pathMap[0])) {
		printf("RETURN(1): Provide pathwindow= OR (pathHotspots= + pathMap=) for segment boundaries.\n");
		printf("Run with --help for usage.\n");
		return 1;
	}
	if (pathHotspots[0] && pathMap[0] == '\0') {
		printf("RETURN(1): pathHotspots= requires pathMap= for SNP indexing.\n");
		printf("Run with --help for usage.\n");
		return 1;
	}
	int ID = posOffspring;
	int IDp1 = posParent1;
	int IDp2 = posParent2;
	 step1 = clock();
	for(int IDrun=0;IDrun<12;IDrun++)// if (ID<10) MAXNBTRIO
	{	distrigametic[IDrun]=0; 
	}
	for(int IDrun=0;IDrun<23;IDrun++)// if (ID<10) MAXNBTRIO
	{	distrigametickeep[IDrun]=0;
	}

	/* Define window/segment sizes before loading HAP (needs nbsnpperchr from map, binary constants, or haps line count) */
	printf("Setting up segment boundaries (before HAP load)...\n");
	setvbuf(stdout, NULL, _IONBF, 0);  /* Unbuffered so progression shows under SLURM/redirect */
	fflush(stdout);
	if (init_nbsnpperchr_before_haps() != 0) return 1;
	if (setup_segment_boundaries() != 0) return 1;
	printf("Segment boundaries ready.\n");
	fflush(stdout);

	/* Load genomes before calculatepihat (which needs them for MAF and pi-hat) */
	printf("Loading HAP files (chr 22 down to 1)...\n");
	fflush(stdout);
	for(int chrtemp1=22;chrtemp1>=1;chrtemp1--)
	{			printf("Loading chr %d/22...\n", chrtemp1);
		fflush(stdout);
		readgenomelocal(chrtemp1, pathHap, &genomes[chrtemp1]);	
		int nsp = nbsnpperchr[chrtemp1];
		int nsp4 = (nsp/4) + ((nsp%4)>0);
		long long cnt0 = 0, cnt1 = 0;
		for (int j = 0; j < n_individuals; j++)
		for (int s = 0; s < nsp; s++) {
			int g = (*((genomes[chrtemp1] + (unsigned long long)j*nsp4) + s/4)>>((s%4)*2))&3;
			int alt = (g&1) + (g>>1);
			cnt0 += 2 - alt;
			cnt1 += alt;
		}
		printf("genomes[chr%d]: n_indiv=%d n_snp=%d alleles_0=%lld alleles_1=%lld\n", chrtemp1, n_individuals, nsp, cnt0, cnt1);
	}
	if (n_individuals <= 0) { printf("RETURN(1): No individuals detected from hap files\n"); return 1; }
	clock_t step2 = clock();
	float elapsed_secs11 = (float)(step2-step1);
	printf("\nTotal time loading data:%f cpu click so %f seconds\n",elapsed_secs11,elapsed_secs11/CLOCKS_PER_SEC );

	calculatepihat(ID);
	clock_t step3= clock();
	elapsed_secs11 = (float)(step3-step1);
	printf("\nTotal time calculate MAF and Pihats:%f cpu click so %f seconds\n",elapsed_secs11,elapsed_secs11/CLOCKS_PER_SEC );
	fflush(stdout);
	printf("Allocating phasing buffers in main()...\n");
	fflush(stdout);
	int (*phaseparent1tap)[NSNPPERCHR] = (int (*)[NSNPPERCHR])calloc(23, sizeof(int[NSNPPERCHR]));
	int (*phaseguesedparent1tap)[NSNPPERCHR] = (int (*)[NSNPPERCHR])calloc(23, sizeof(int[NSNPPERCHR]));
	double (*pihatagainstallchrMPphaseerror)[MAXPOP][200][2] = (double (*)[MAXPOP][200][2])calloc(23, sizeof(double[MAXPOP][200][2]));
	if (!phaseparent1tap || !phaseguesedparent1tap) {
		printf("RETURN(1): Failed to allocate phase arrays in main\n");
		free(phaseparent1tap); free(phaseguesedparent1tap); free(pihatagainstallchrMPphaseerror);
		return 1;
	}
	if (!pihatagainstallchrMPphaseerror) {
		printf("RETURN(1): Failed to allocate pihatagainstallchrMPphaseerror (~32 GB) in main\n");
		free(phaseparent1tap); free(phaseguesedparent1tap); free(pihatagainstallchrMPphaseerror);
		return 1;
	}
	
	
	printf("Calling acrosschrphasing...\n");
	fflush(stdout);
	int phasing_ret = acrosschrphasing(ID, IDp1, IDp2, pathresult, pathHap, pathwindow,
		phaseparent1tap, phaseguesedparent1tap, pihatagainstallchrMPphaseerror);
	free(phaseparent1tap);
	free(phaseguesedparent1tap);
	free(pihatagainstallchrMPphaseerror);
	clock_t step4= clock();
	if (phasing_ret != 0) return phasing_ret;
	float elapsed_secs2 = (float)(step4-step1);
	printf("\nTotal time phasing:%f cpu click so %f seconds\n",elapsed_secs2,elapsed_secs2/CLOCKS_PER_SEC );

	printf("RETURN(0): main completed successfully.\n");
	return 0;
}
