# Générateur de Données Aléatoires

## Description

Ce générateur crée des fichiers de données d'entrée aléatoires au format HAP pour le programme de phasing `ProgramPhasing.cpp`. Il génère des données génétiques simulées avec les caractéristiques suivantes :

- **1000 individus** : Chaque fichier contient les génotypes de 1000 individus
- **1000 SNPs par chromosome** : Chaque chromosome contient exactement 1000 variants
- **22 chromosomes autosomiques** : Génère les chromosomes 1 à 22
- **Format HAP** : Compatible avec le format d'entrée attendu par le programme de phasing
- **Fichier MAF** : Génère également un fichier de fréquences alléliques mineures (MAF)

## Format des Données Générées

### Format HAP

Chaque fichier `<chromosome>.hap` contient une ligne par SNP au format :

```
rsID position allele0 allele1 individual1_allele0 individual1_allele1 individual2_allele0 individual2_allele1 ...
```

**Exemple :**
```
rs100000 123456 A T 0 1 0 0 1 1 ...
rs100001 234567 G C 0 0 0 1 1 1 ...
```

**Détails :**
- `rsID` : Identifiant du SNP (format rsXXXXXX)
- `position` : Position génomique en bases
- `allele0` : Allèle de référence (A, T, G, ou C)
- `allele1` : Allèle alternatif (A, T, G, ou C, différent de allele0)
- Pour chaque individu : deux valeurs (0 ou 1) représentant les deux allèles
  - `0` = allèle de référence
  - `1` = allèle alternatif

### Format MAF

Le fichier `MAF.txt` contient les fréquences alléliques mineures pour chaque SNP :

```
chromosome snp_index maf_value
```

**Exemple :**
```
1 0 0.234567
1 1 0.156789
2 0 0.345678
...
```

## Compilation

### Prérequis

- Compilateur C++11 ou supérieur (g++, clang++, ou MSVC)
- Bibliothèque standard C++

### Commande de Compilation

**Linux/Mac :**
```bash
g++ -std=c++11 -O2 -o GenerateurDonneesAleatoires GenerateurDonneesAleatoires.cpp
```

**Windows (MinGW/MSVC) :**
```bash
g++ -std=c++11 -O2 -o GenerateurDonneesAleatoires.exe GenerateurDonneesAleatoires.cpp
```

**Ou avec cl.exe (Visual Studio) :**
```bash
cl /EHsc /O2 GenerateurDonneesAleatoires.cpp /Fe:GenerateurDonneesAleatoires.exe
```

## Utilisation

### Utilisation de Base

Génère les données dans le répertoire `./donnees_aleatoires/` :

```bash
./GenerateurDonneesAleatoires
```

### Spécifier le Répertoire de Sortie

```bash
./GenerateurDonneesAleatoires /chemin/vers/sortie/
```

### Spécifier une Graine Aléatoire

Pour obtenir des résultats reproductibles :

```bash
./GenerateurDonneesAleatoires ./donnees_aleatoires/ 12345
```

### Structure des Fichiers Générés

Après exécution, vous obtiendrez :

```
donnees_aleatoires/
├── 1.hap          # Chromosome 1 (1000 SNPs, 1000 individus)
├── 2.hap          # Chromosome 2 (1000 SNPs, 1000 individus)
├── ...
├── 22.hap         # Chromosome 22 (1000 SNPs, 1000 individus)
└── MAF.txt        # Fréquences alléliques mineures
```

## Caractéristiques des Données Générées

### Distribution des Génotypes

Les génotypes sont générés selon le modèle d'équilibre de Hardy-Weinberg :

- **Homozygote référence (0/0)** : Probabilité = (1 - MAF)²
- **Hétérozygote (0/1)** : Probabilité = 2 × (1 - MAF) × MAF
- **Homozygote alternatif (1/1)** : Probabilité = MAF²

### Fréquences Alléliques

- Les MAF sont générées aléatoirement entre **0.01** et **0.5**
- Distribution uniforme pour une variété réaliste de fréquences

### Positions Génomiques

- Positions générées aléatoirement entre **1000** et **250,000,000** bases
- Distribution uniforme (non ordonnée, ce qui est acceptable pour la simulation)

## Utilisation avec ProgramPhasing

Une fois les données générées, vous pouvez les utiliser avec le programme de phasing :

```bash
./ProgramPhasing \
  -NbIndiv 1000 \
  -PathInput ./donnees_aleatoires/ \
  -PathOutput ./resultats_phasing/ \
  -PathMAF ./donnees_aleatoires/MAF.txt \
  -Verbose 1
```

## Paramètres Modifiables

Pour modifier les paramètres de génération, éditez les constantes dans `GenerateurDonneesAleatoires.cpp` :

```cpp
const int NB_INDIVIDUS = 1000;           // Nombre d'individus
const int NB_SNPS_PAR_CHROMOSOME = 1000;  // SNPs par chromosome
const int CHROMOSOME_MIN = 1;             // Premier chromosome
const int CHROMOSOME_MAX = 22;            // Dernier chromosome
const int POSITION_MIN = 1000;            // Position minimale
const int POSITION_MAX = 250000000;       // Position maximale
```

## Notes Techniques

- **Générateur aléatoire** : Utilise Mersenne Twister (std::mt19937) pour une bonne qualité aléatoire
- **Performance** : Génération rapide grâce à l'écriture directe dans les fichiers
- **Mémoire** : Utilisation mémoire minimale (génération ligne par ligne)
- **Reproductibilité** : Utilisez une graine fixe pour obtenir les mêmes données

## Exemple Complet

```bash
# 1. Compiler le générateur
g++ -std=c++11 -O2 -o GenerateurDonneesAleatoires GenerateurDonneesAleatoires.cpp

# 2. Générer les données
./GenerateurDonneesAleatoires ./mes_donnees/ 42

# 3. Vérifier les fichiers générés
ls -lh ./mes_donnees/

# 4. Utiliser avec le programme de phasing
./ProgramPhasing -NbIndiv 1000 -PathInput ./mes_donnees/ -PathOutput ./resultats/
```

## Dépannage

### Erreur : Impossible d'ouvrir le fichier
- Vérifiez que le répertoire de sortie existe ou peut être créé
- Vérifiez les permissions d'écriture

### Fichiers incomplets
- Vérifiez l'espace disque disponible
- Les fichiers peuvent être volumineux (~2-4 MB par chromosome)

### Données non reproductibles
- Utilisez toujours la même graine aléatoire pour obtenir les mêmes résultats

## Auteur

Générateur créé pour le projet de phasing génétique.

