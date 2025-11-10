# Guide pour transférer le projet sur GitHub

## Prérequis

### 1. Installer Git

**Windows :**
- Téléchargez Git depuis : https://git-scm.com/download/win
- Installez-le avec les options par défaut
- Redémarrez votre terminal après l'installation

**Vérification :**
```bash
git --version
```

### 2. Configurer Git (première fois uniquement)

```bash
git config --global user.name "Votre Nom"
git config --global user.email "votre.email@example.com"
```

## Étapes pour transférer sur GitHub

### Étape 1 : Créer un dépôt sur GitHub

1. Allez sur https://github.com
2. Connectez-vous à votre compte (ou créez-en un)
3. Cliquez sur le bouton "+" en haut à droite
4. Sélectionnez "New repository"
5. Donnez un nom à votre dépôt (ex: "ProgramPhasing")
6. Choisissez Public ou Private
7. **NE PAS** cocher "Initialize this repository with a README"
8. Cliquez sur "Create repository"

### Étape 2 : Initialiser Git dans votre projet

Ouvrez un terminal PowerShell dans le répertoire du projet :

```powershell
cd "d:\tranfett D\Fileforgithub"
```

Initialisez le dépôt Git :

```bash
git init
```

### Étape 3 : Ajouter tous les fichiers

```bash
git add .
```

Vérifiez les fichiers qui seront ajoutés :

```bash
git status
```

### Étape 4 : Faire le premier commit

```bash
git commit -m "Initial commit: Programme de phasing génétique avec générateur de données"
```

### Étape 5 : Ajouter le dépôt distant GitHub

Remplacez `VOTRE_USERNAME` et `NOM_DEPOT` par vos valeurs :

```bash
git remote add origin https://github.com/VOTRE_USERNAME/NOM_DEPOT.git
```

Exemple :
```bash
git remote add origin https://github.com/johndoe/ProgramPhasing.git
```

### Étape 6 : Pousser le code sur GitHub

```bash
git branch -M main
git push -u origin main
```

Vous devrez entrer vos identifiants GitHub.

## Commandes Git utiles

### Voir l'état du dépôt
```bash
git status
```

### Ajouter des fichiers modifiés
```bash
git add .
git commit -m "Description des modifications"
git push
```

### Voir l'historique des commits
```bash
git log
```

### Cloner un dépôt (pour récupérer sur un autre ordinateur)
```bash
git clone https://github.com/VOTRE_USERNAME/NOM_DEPOT.git
```

## Note sur les fichiers de données

Le fichier `.gitignore` exclut le dossier `donnees_aleatoires/` car ces fichiers sont volumineux et peuvent être régénérés.

Si vous voulez inclure les fichiers de données générées, modifiez `.gitignore` et supprimez ces lignes :
```
donnees_aleatoires/
*.hap
MAF.txt
```

## Authentification GitHub

Si vous avez des problèmes d'authentification lors du `git push`, vous pouvez :

1. **Utiliser un Personal Access Token** :
   - Allez dans GitHub Settings > Developer settings > Personal access tokens > Tokens (classic)
   - Créez un nouveau token avec les permissions `repo`
   - Utilisez ce token comme mot de passe lors du push

2. **Utiliser GitHub CLI** :
   ```bash
   gh auth login
   ```

## Structure du projet sur GitHub

Une fois poussé, votre dépôt GitHub contiendra :

```
ProgramPhasing/
├── .gitignore
├── README.md
├── README_Generateur.md
├── README_ARCHITECTURE.md
├── README_MODULAR.md
├── USER_MANUAL.md
├── ProgramPhasing.cpp
├── GenerateurDonneesAleatoires.cpp
├── GenererDonnees.ps1
├── Makefile
├── include/
│   └── ...
├── src/
│   └── ...
└── ...
```

Les fichiers dans `donnees_aleatoires/` ne seront pas inclus (voir `.gitignore`).

