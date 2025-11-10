/**
 * @file GenerateurDonneesAleatoires.cpp
 * @brief Générateur de données d'entrée aléatoires pour le programme de phasing
 * 
 * Ce générateur crée des fichiers HAP au format attendu par ProgramPhasing.cpp
 * avec les caractéristiques suivantes :
 * - 1000 individus
 * - 1000 SNPs par chromosome
 * - Chromosomes 1 à 22 (autosomaux)
 * 
 * Format de sortie HAP :
 * rsID position allele0 allele1 individual1_allele0 individual1_allele1 ...
 * 
 * @author Générateur de données
 * @date 2024
 */

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <iomanip>
#include <vector>
#include <sstream>
#include <cstdlib>

// ============================================================================
// CONSTANTES DE CONFIGURATION
// ============================================================================

const int NB_INDIVIDUS = 1000;           // Nombre d'individus à générer
const int NB_SNPS_PAR_CHROMOSOME = 1000;  // Nombre de SNPs par chromosome
const int CHROMOSOME_MIN = 1;             // Premier chromosome (autosomique)
const int CHROMOSOME_MAX = 22;            // Dernier chromosome (autosomique)
const int POSITION_MIN = 1000;            // Position minimale (en bases)
const int POSITION_MAX = 250000000;       // Position maximale (taille approximative du plus grand chromosome)

// Allèles possibles (nucléotides)
const char ALLELES[] = {'A', 'T', 'G', 'C'};
const int NB_ALLELES = 4;

// ============================================================================
// CLASSE GÉNÉRATEUR DE DONNÉES
// ============================================================================

class GenerateurDonneesAleatoires {
private:
    std::mt19937 generator;                    // Générateur Mersenne Twister
    std::uniform_int_distribution<int> dist_pos;  // Distribution pour les positions
    std::uniform_int_distribution<int> dist_allele; // Distribution pour les allèles
    std::uniform_real_distribution<double> dist_uniforme; // Distribution uniforme [0,1]
    
    /**
     * @brief Génère un identifiant rsID aléatoire
     * @param chromosome Numéro du chromosome
     * @param snpIndex Index du SNP (0-based)
     * @return String au format "rsXXXXXX"
     */
    std::string genererRsID(int chromosome, int snpIndex) {
        int rsNumber = chromosome * 100000 + snpIndex;
        std::ostringstream oss;
        oss << "rs" << std::setw(6) << std::setfill('0') << rsNumber;
        return oss.str();
    }
    
    /**
     * @brief Génère une position génomique aléatoire
     * @return Position en bases
     */
    int genererPosition() {
        return dist_pos(generator);
    }
    
    /**
     * @brief Génère deux allèles différents (allèle de référence et allèle alternatif)
     * @param allele0 Référence pour l'allèle 0 (sortie)
     * @param allele1 Référence pour l'allèle 1 (sortie)
     */
    void genererAlleles(char& allele0, char& allele1) {
        int idx0 = dist_allele(generator);
        int idx1 = dist_allele(generator);
        
        // S'assurer que les deux allèles sont différents
        while (idx1 == idx0) {
            idx1 = dist_allele(generator);
        }
        
        allele0 = ALLELES[idx0];
        allele1 = ALLELES[idx1];
    }
    
    /**
     * @brief Génère un génotype pour un individu selon la fréquence allélique
     * @param maf Minor Allele Frequency (fréquence de l'allèle mineur)
     * @return Génotype encodé : 0 (homozygote référence), 1 (hétérozygote), 2 (homozygote alternatif)
     */
    int genererGenotype(double maf) {
        double random = dist_uniforme(generator);
        
        // Calcul des probabilités selon le modèle Hardy-Weinberg
        // p = 1 - maf (fréquence allèle majeur)
        // q = maf (fréquence allèle mineur)
        // P(homozygote majeur) = p^2
        // P(hétérozygote) = 2pq
        // P(homozygote mineur) = q^2
        
        double p = 1.0 - maf;
        double p_homo_majeur = p * p;
        double p_heterozygote = 2.0 * p * maf;
        
        if (random < p_homo_majeur) {
            return 0; // Homozygote référence (0/0)
        } else if (random < p_homo_majeur + p_heterozygote) {
            return 1; // Hétérozygote (0/1)
        } else {
            return 2; // Homozygote alternatif (1/1)
        }
    }
    
    /**
     * @brief Convertit un génotype en deux allèles pour le format HAP
     * @param genotype Génotype encodé (0, 1, ou 2)
     * @param allele0 Allèle 0 (sortie)
     * @param allele1 Allèle 1 (sortie)
     */
    void genotypeVersAlleles(int genotype, int& allele0, int& allele1) {
        switch (genotype) {
            case 0: // Homozygote référence
                allele0 = 0;
                allele1 = 0;
                break;
            case 1: // Hétérozygote
                allele0 = 0;
                allele1 = 1;
                break;
            case 2: // Homozygote alternatif
                allele0 = 1;
                allele1 = 1;
                break;
            default:
                allele0 = 0;
                allele1 = 0;
        }
    }

public:
    /**
     * @brief Constructeur - initialise les générateurs aléatoires
     * @param seed Graine pour le générateur aléatoire (optionnel)
     */
    GenerateurDonneesAleatoires(unsigned int seed = std::random_device{}()) 
        : generator(seed),
          dist_pos(POSITION_MIN, POSITION_MAX),
          dist_allele(0, NB_ALLELES - 1),
          dist_uniforme(0.0, 1.0) {
    }
    
    /**
     * @brief Génère un fichier HAP pour un chromosome donné
     * @param chromosome Numéro du chromosome (1-22)
     * @param cheminSortie Chemin de base pour les fichiers de sortie
     * @param genererMAF Si true, génère aussi un fichier MAF
     * @return true si succès, false sinon
     */
    bool genererChromosome(int chromosome, const std::string& cheminSortie, bool genererMAF = true) {
        if (chromosome < CHROMOSOME_MIN || chromosome > CHROMOSOME_MAX) {
            std::cerr << "Erreur : Chromosome invalide (" << chromosome << ")" << std::endl;
            return false;
        }
        
        // Créer le nom du fichier HAP
        std::ostringstream oss_hap;
        oss_hap << cheminSortie << chromosome << ".hap";
        std::string nomFichierHAP = oss_hap.str();
        
        // Ouvrir le fichier HAP en écriture
        std::ofstream fichierHAP(nomFichierHAP);
        if (!fichierHAP.is_open()) {
            std::cerr << "Erreur : Impossible d'ouvrir le fichier " << nomFichierHAP << std::endl;
            return false;
        }
        
        // Vecteur pour stocker les MAF si nécessaire
        std::vector<double> maf_values;
        if (genererMAF) {
            maf_values.reserve(NB_SNPS_PAR_CHROMOSOME);
        }
        
        std::cout << "Génération du chromosome " << chromosome << "..." << std::endl;
        
        // Générer chaque SNP
        for (int snp = 0; snp < NB_SNPS_PAR_CHROMOSOME; snp++) {
            // Générer les informations du SNP
            std::string rsID = genererRsID(chromosome, snp);
            int position = genererPosition();
            char allele0, allele1;
            genererAlleles(allele0, allele1);
            
            // Générer une MAF aléatoire (entre 0.01 et 0.5 pour être réaliste)
            double maf = 0.01 + dist_uniforme(generator) * 0.49;
            if (genererMAF) {
                maf_values.push_back(maf);
            }
            
            // Écrire l'en-tête du SNP : rsID position allele0 allele1
            fichierHAP << rsID << " " << position << " " << allele0 << " " << allele1;
            
            // Générer les génotypes pour chaque individu
            for (int indiv = 0; indiv < NB_INDIVIDUS; indiv++) {
                int genotype = genererGenotype(maf);
                int a0, a1;
                genotypeVersAlleles(genotype, a0, a1);
                fichierHAP << " " << a0 << " " << a1;
            }
            
            fichierHAP << std::endl;
            
            // Afficher la progression tous les 100 SNPs
            if ((snp + 1) % 100 == 0) {
                std::cout << "  SNP " << (snp + 1) << "/" << NB_SNPS_PAR_CHROMOSOME << std::endl;
            }
        }
        
        fichierHAP.close();
        std::cout << "Fichier HAP créé : " << nomFichierHAP << std::endl;
        
        // Générer le fichier MAF si demandé
        if (genererMAF) {
            std::ostringstream oss_maf;
            oss_maf << cheminSortie << "MAF.txt";
            std::string nomFichierMAF = oss_maf.str();
            
            std::ofstream fichierMAF(nomFichierMAF, std::ios::app); // Mode append pour tous les chromosomes
            if (fichierMAF.is_open()) {
                for (int snp = 0; snp < NB_SNPS_PAR_CHROMOSOME; snp++) {
                    fichierMAF << chromosome << " " << snp << " " 
                               << std::fixed << std::setprecision(6) << maf_values[snp] << std::endl;
                }
                fichierMAF.close();
            }
        }
        
        return true;
    }
    
    /**
     * @brief Génère tous les chromosomes (1-22)
     * @param cheminSortie Chemin de base pour les fichiers de sortie
     * @param genererMAF Si true, génère aussi un fichier MAF
     * @return true si succès, false sinon
     */
    bool genererTousChromosomes(const std::string& cheminSortie, bool genererMAF = true) {
        // Supprimer le fichier MAF existant s'il existe
        if (genererMAF) {
            std::ostringstream oss_maf;
            oss_maf << cheminSortie << "MAF.txt";
            std::remove(oss_maf.str().c_str());
        }
        
        bool succes = true;
        for (int chr = CHROMOSOME_MIN; chr <= CHROMOSOME_MAX; chr++) {
            if (!genererChromosome(chr, cheminSortie, genererMAF)) {
                succes = false;
            }
        }
        
        return succes;
    }
};

// ============================================================================
// FONCTION PRINCIPALE
// ============================================================================

int main(int argc, char* argv[]) {
    std::cout << "========================================" << std::endl;
    std::cout << "Générateur de Données Aléatoires" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Configuration :" << std::endl;
    std::cout << "  - Nombre d'individus : " << NB_INDIVIDUS << std::endl;
    std::cout << "  - SNPs par chromosome : " << NB_SNPS_PAR_CHROMOSOME << std::endl;
    std::cout << "  - Chromosomes : " << CHROMOSOME_MIN << "-" << CHROMOSOME_MAX << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Déterminer le chemin de sortie
    std::string cheminSortie = "./donnees_aleatoires/";
    if (argc > 1) {
        cheminSortie = argv[1];
        // S'assurer que le chemin se termine par un séparateur
        if (cheminSortie.back() != '/' && cheminSortie.back() != '\\') {
            cheminSortie += "/";
        }
    }
    
    // Créer le répertoire de sortie si nécessaire
    #ifdef _WIN32
        std::string commande = "if not exist \"" + cheminSortie + "\" mkdir \"" + cheminSortie + "\"";
    #else
        std::string commande = "mkdir -p \"" + cheminSortie + "\"";
    #endif
    system(commande.c_str());
    
    // Initialiser le générateur
    unsigned int seed = std::random_device{}();
    if (argc > 2) {
        seed = std::stoul(argv[2]);
        std::cout << "Graine aléatoire spécifiée : " << seed << std::endl;
    } else {
        std::cout << "Graine aléatoire générée : " << seed << std::endl;
    }
    
    GenerateurDonneesAleatoires generateur(seed);
    
    // Générer tous les chromosomes
    std::cout << std::endl << "Début de la génération..." << std::endl;
    bool succes = generateur.genererTousChromosomes(cheminSortie, true);
    
    if (succes) {
        std::cout << std::endl << "========================================" << std::endl;
        std::cout << "Génération terminée avec succès !" << std::endl;
        std::cout << "Fichiers créés dans : " << cheminSortie << std::endl;
        std::cout << "  - 1.hap à 22.hap (fichiers HAP)" << std::endl;
        std::cout << "  - MAF.txt (fréquences alléliques mineures)" << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } else {
        std::cerr << std::endl << "Erreur lors de la génération !" << std::endl;
        return 1;
    }
}

