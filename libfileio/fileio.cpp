#include "fileio.h"
// Note: readinteger.h is already included through types.h via fileio.h

int readgenomelocal(char pathfile[], int chr, int run, int step, unsigned char* gentomodify, int NbIndiv, int nbsnpperchr[], int nbsnpperchrinfile[])
{
	// Suppress unused parameter warnings (these parameters are kept for API compatibility)
	(void)run;
	(void)step;
	
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
	// Read Relate .haps format
	// Format: Chromosome SNP_ID Position Ancestral_allele Alternative_allele allele1 allele2 ...
	// Reference: https://myersgroup.github.io/relate/input_data.html
	char line[1000000];  // Large buffer for line reading
	int snp = 0;
	int debug_snp0_alt_count = 0;  // Debug: count alternate alleles for first SNP
	printf("Reading file: %s corresponding to chromosome %d (Relate .haps format)\n", pathfilechr, chr);
	
	while (fgets(line, sizeof(line), endfile) != NULL)
	{
		// Skip empty lines
		if (line[0] == '\n' || line[0] == '\r' || line[0] == '\0')
			continue;
		
		// Parse the line manually to avoid strtok issues with long lines
		// Format: Chromosome SNP_ID Position Ancestral_allele Alternative_allele allele1 allele2 ...
		char* p = line;
		
		// Skip whitespace at start
		while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')
			p++;
		if (*p == '\0')
			continue;
		
		// Skip first 5 columns (Chromosome, SNP_ID, Position, Ancestral_allele, Alternative_allele)
		for (int col = 0; col < 5; col++)
		{
			// Skip current token
			while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r')
				p++;
			// Skip whitespace
			while (*p == ' ' || *p == '\t')
				p++;
			if (*p == '\0' || *p == '\n' || *p == '\r')
			{
				printf("ERROR: Not enough columns at SNP %d (expected at least 5 metadata columns)\n", snp);
				return 1;
			}
		}
		
		// Now read alleles for each individual (two alleles per individual)
		for (int relat = 0; relat < NbIndiv; relat++)
		{
			// Read first allele
			if (*p == '\0' || *p == '\n' || *p == '\r')
			{
				printf("ERROR: Not enough alleles for individual %d at SNP %d\n", relat, snp);
				return 1;
			}
			int allele1 = (*p == '1') ? 1 : 0;
			p++;
			
			// Skip whitespace
			while (*p == ' ' || *p == '\t')
				p++;
			
			// Read second allele
			if (*p == '\0' || *p == '\n' || *p == '\r')
			{
				printf("ERROR: Not enough alleles for individual %d at SNP %d\n", relat, snp);
				return 1;
			}
			int allele2 = (*p == '1') ? 1 : 0;
			p++;
			
			// Skip whitespace before next pair
			while (*p == ' ' || *p == '\t')
				p++;
			
			// Encode genotype: 0=00, 1=01, 2=10, 3=11
			int geno = (allele1 << 1) | allele2;
			// Debug: count alternate alleles for first SNP
			if (snp == 0 && chr == 1)
			{
				debug_snp0_alt_count += allele1 + allele2;
				
			}
			// Store in compressed format (2 bits per genotype)
			*(gentomodify + (unsigned long long)(((((unsigned long long)relat) * (nbsnpperchr[chr] / 4 + ((nbsnpperchr[chr] % 4) > 0))) + snp / 4))) =
				(*(gentomodify + (unsigned long long)(((((unsigned long long)relat) * (nbsnpperchr[chr] / 4 + ((nbsnpperchr[chr] % 4) > 0))) + snp / 4))) &
				(~(3 << ((snp % 4) * 2)))) |
				((geno) << ((snp % 4) * 2));
		}
		
		snp++;
	}
	
	nbsnpperchrinfile[chr] = snp;
	printf("%d SNPs in file and %d individuals\n", nbsnpperchrinfile[chr], NbIndiv);
	// Debug: print alternate allele count for first SNP
	if (chr == 1 && debug_snp0_alt_count > 0)
	{
		printf("READ SNP 0: total_alt_alleles=%d/%d (expected ~%d-%d for MAF=0.05-0.5)\n", 
			debug_snp0_alt_count, NbIndiv * 2, (int)(NbIndiv * 2 * 0.05), (int)(NbIndiv * 2 * 0.5));
	}
	fclose(endfile);
	return 0;
}

int readParentInfo(const char* filename, structparentinfo parentinfo[], int* nbparentinfo, int maxIndiv)
{
	FILE* fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("ERROR: Could not open parent info file %s\n", filename);
		return 1;
	}
	
	*nbparentinfo = 0;
	char line[256];
	int lineNum = 0;
	
	while (fgets(line, sizeof(line), fp) != NULL && *nbparentinfo < maxIndiv)
	{
		lineNum++;
		if (line[0] == '\n' || line[0] == '\r')
			continue;
		
		int id, p1, p2;
		if (sscanf(line, "%d %d %d", &id, &p1, &p2) == 3)
		{
			if (id < 0 || id >= maxIndiv)
			{
				printf("WARNING: Individual ID %d on line %d is out of range [0, %d). Skipping.\n",
					id, lineNum, maxIndiv);
				continue;
			}
			parentinfo[*nbparentinfo].ID = id;
			parentinfo[*nbparentinfo].parent1 = p1;
			parentinfo[*nbparentinfo].parent2 = p2;
			(*nbparentinfo)++;
		}
		else
		{
			printf("WARNING: Invalid format on line %d. Expected: ID parent1 parent2\n", lineNum);
		}
	}
	
	fclose(fp);
	printf("Successfully loaded %d parent information entries from file %s\n", *nbparentinfo, filename);
	return 0;
}

int getParentInfo(int indivID, int* parent1, int* parent2, structparentinfo parentinfo[], int nbparentinfo)
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

int readListIndiv(const char* filename, int listIndivToProcess[], int* nbIndivToProcess, int* useListIndiv, int NbIndiv, int maxIndiv)
{
	FILE* fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("ERROR: Could not open file %s for reading individual list\n", filename);
		return 1;
	}
	
	*nbIndivToProcess = 0;
	char line[256];
	int lineNum = 0;
	
	while (fgets(line, sizeof(line), fp) != NULL && *nbIndivToProcess < maxIndiv)
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
		
		listIndivToProcess[*nbIndivToProcess] = indivID;
		(*nbIndivToProcess)++;
	}
	
	fclose(fp);
	
	if (*nbIndivToProcess == 0)
	{
		printf("ERROR: No valid individual IDs found in file %s\n", filename);
		return 1;
	}
	
	printf("Successfully loaded %d individuals from file %s\n", *nbIndivToProcess, filename);
	*useListIndiv = 1;
	return 0;
}

/**
 * @brief Read a PED file and load genotype data into the genomes array
 * @param filename Path to the PED file
 * @param chr Chromosome number (1-22)
 * @param genomes Array of genome data pointers (one per chromosome)
 * @param nbsnpperchr Array storing number of SNPs per chromosome
 * @param nbsnpperchrinfile Output: number of SNPs read from file
 * @param NbIndiv Output: number of individuals read from file
 * @param maxIndiv Maximum number of individuals to read
 * @return 0 on success, 1 on error
 * 
 * PED file format:
 * Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype Allele1_SNP1 Allele2_SNP1 Allele1_SNP2 Allele2_SNP2 ...
 * Each line represents one individual
 * Alleles are encoded as '0', '1', or missing data markers
 */
int readpedfile(const char* filename, int chr, unsigned char* genomes[], int nbsnpperchr[], int* nbsnpperchrinfile, int* NbIndiv, int maxIndiv)
{
	FILE* fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("ERROR: Could not open PED file %s for reading\n", filename);
		return 1;
	}
	
	char line[1000000];  // Large buffer for PED file lines
	int lineNum = 0;
	*NbIndiv = 0;
	*nbsnpperchrinfile = 0;
	int first_line = 1;
	
	printf("Reading PED file: %s (chromosome %d)\n", filename, chr);
	
	while (fgets(line, sizeof(line), fp) != NULL && *NbIndiv < maxIndiv)
	{
		lineNum++;
		
		// Skip empty lines and comments
		if (line[0] == '\n' || line[0] == '\r' || line[0] == '#')
			continue;
		
		// Parse the line
		char* token = strtok(line, " \t\n\r");
		if (token == NULL)
			continue;
		
		// Skip first 6 columns: Family_ID, Individual_ID, Father_ID, Mother_ID, Sex, Phenotype
		for (int i = 0; i < 5; i++)
		{
			token = strtok(NULL, " \t\n\r");
			if (token == NULL)
			{
				printf("WARNING: Line %d has insufficient columns, skipping\n", lineNum);
				break;
			}
		}
		
		if (token == NULL)
			continue;
		
		// Count SNPs on first line
		if (first_line)
		{
			// Make a copy of the line for counting (strtok modifies the string)
			char line_copy[1000000];
			strncpy(line_copy, line, sizeof(line_copy) - 1);
			line_copy[sizeof(line_copy) - 1] = '\0';
			
			// Count tokens after the first 6 columns using the copy
			char* temp_token = strtok(line_copy, " \t\n\r");
			for (int i = 0; i < 6 && temp_token != NULL; i++)
				temp_token = strtok(NULL, " \t\n\r");
			
			int snp_count = 0;
			while (temp_token != NULL)
			{
				snp_count++;
				temp_token = strtok(NULL, " \t\n\r");
			}
			*nbsnpperchrinfile = snp_count / 2;  // Two alleles per SNP
			printf("Detected %d SNPs in PED file\n", *nbsnpperchrinfile);
			
			// Now tokenize the original line for processing
			token = strtok(line, " \t\n\r");
			for (int i = 0; i < 6 && token != NULL; i++)
				token = strtok(NULL, " \t\n\r");
			first_line = 0;
		}
		
		// Read genotype data for this individual
		int snp = 0;
		while (token != NULL && snp < *nbsnpperchrinfile && snp < nbsnpperchr[chr])
		{
			// Read first allele
			int allele1 = -1;
			if (token[0] >= '0' && token[0] <= '9')
				allele1 = token[0] - '0';
			else if (token[0] == '-' || token[0] == 'N' || token[0] == 'n')
				allele1 = 0;  // Treat missing as 0
			
			// Read second allele
			token = strtok(NULL, " \t\n\r");
			int allele2 = -1;
			if (token != NULL)
			{
				if (token[0] >= '0' && token[0] <= '9')
					allele2 = token[0] - '0';
				else if (token[0] == '-' || token[0] == 'N' || token[0] == 'n')
					allele2 = 0;  // Treat missing as 0
			}
			
			if (allele1 >= 0 && allele2 >= 0)
			{
				// Encode genotype: 0=00, 1=01, 2=10, 3=11
				int geno = (allele1 << 1) | allele2;
				
				// Store in compressed format (2 bits per genotype)
				unsigned long long offset = (unsigned long long)(*NbIndiv) * (nbsnpperchr[chr] / 4 + ((nbsnpperchr[chr] % 4) > 0)) + snp / 4;
				genomes[chr][offset] = (genomes[chr][offset] & (~(3 << ((snp % 4) * 2)))) | (geno << ((snp % 4) * 2));
			}
			
			snp++;
			token = strtok(NULL, " \t\n\r");
		}
		
		(*NbIndiv)++;
		
		// Progress indicator
		if ((*NbIndiv) % 100 == 0)
		{
			printf("  Progress: %d individuals read\n", *NbIndiv);
		}
	}
	
	fclose(fp);
	printf("Successfully read %d individuals and %d SNPs from PED file %s\n", 
		*NbIndiv, *nbsnpperchrinfile, filename);
	return 0;
}

/**
 * @brief Write genotype data to a PED file
 * @param chr Chromosome number (1-22)
 * @param pathoutput Base path for output file
 * @param genomes Array of genome data pointers
 * @param nbsnpperchr Array storing number of SNPs per chromosome
 * @param nbsnpperchrinfile Number of SNPs to write
 * @param NbIndiv Number of individuals to write
 * @return 0 on success, 1 on error
 * 
 * PED file format:
 * Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype Allele1_SNP1 Allele2_SNP1 ...
 */
int writepedfile(int chr, const char* pathoutput, unsigned char* genomes[], int nbsnpperchr[], int nbsnpperchrinfile[], int NbIndiv)
{
	if (chr < 0 || chr >= 23)
	{
		printf("ERROR: writepedfile received invalid chromosome index %d\n", chr);
		return 1;
	}
	if (genomes[chr] == NULL)
	{
		printf("WARNING: Genome buffer is null for chromosome %d, skipping write.\n", chr);
		return 0;
	}
	if (nbsnpperchrinfile[chr] <= 0)
	{
		printf("WARNING: No SNPs loaded for chromosome %d (nbsnpperchrinfile=%d). Skipping PED write for this chromosome.\n",
			chr, nbsnpperchrinfile[chr]);
		return 0;
	}
	
	char filename[300];
	strcpy(filename, pathoutput);
	char number[100];
	sprintf(number, "%d", chr);
	strcat(filename, number);
	strcat(filename, ".ped");
	
	printf("Writing PED file: %s (chromosome %d, %d individuals, %d SNPs)\n", 
		filename, chr, NbIndiv, nbsnpperchrinfile[chr]);
	
	FILE* fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("ERROR: Could not open file %s for writing\n", filename);
		return 1;
	}
	
	// Write each individual
	for (int person = 0; person < NbIndiv; person++)
	{
		// Write PED header: Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype
		// Using person ID for both Family_ID and Individual_ID
		fprintf(fp, "%d %d 0 0 -9 0", person, person);
		
		// Write genotype data for all SNPs
		for (int snp = 0; snp < nbsnpperchrinfile[chr]; snp++)
		{
			// Extract genotype from compressed format
			unsigned long long offset = (unsigned long long)person * (nbsnpperchr[chr] / 4 + ((nbsnpperchr[chr] % 4) > 0)) + snp / 4;
			int geno = (genomes[chr][offset] >> ((snp % 4) * 2)) & 3;
			
			// Convert to allele pair
			int allele1 = (geno >> 1) & 1;
			int allele2 = geno & 1;
			
			// Write alleles
			fprintf(fp, " %d %d", allele1, allele2);
		}
		
		fprintf(fp, "\n");
		
		// Progress indicator
		if ((person + 1) % 100 == 0)
		{
			printf("  Progress: %d/%d individuals written\n", person + 1, NbIndiv);
		}
	}
	
	fclose(fp);
	printf("Successfully wrote PED file: %s\n", filename);
	return 0;
}

/**
 * @brief Legacy function for writing output (writes PED format)
 * @deprecated Use writepedfile() instead
 */
int writeoutput(int chr, char pathoutput[], unsigned char* genomes[], int nbsnpperchr[], int nbsnpperchrinfile[], int NbIndiv)
{
	return writepedfile(chr, pathoutput, genomes, nbsnpperchr, nbsnpperchrinfile, NbIndiv);
}

