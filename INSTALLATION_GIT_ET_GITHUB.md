# Guide d'Installation Git et Création du Dépôt GitHub

## Étape 1 : Installer Git

### Windows

1. **Télécharger Git**
   - Allez sur : https://git-scm.com/download/win
   - Téléchargez la dernière version pour Windows
   - Exécutez l'installateur

2. **Configuration de l'installation**
   - Acceptez les options par défaut
   - Choisissez votre éditeur préféré (ou gardez Vim)
   - Sélectionnez "Git from the command line and also from 3rd-party software"
   - Acceptez les autres options par défaut

3. **Vérifier l'installation**
   - Ouvrez un nouveau PowerShell ou CMD
   - Exécutez : `git --version`
   - Vous devriez voir la version de Git installée

### Configurer Git (première fois)

Ouvrez PowerShell et exécutez :

```powershell
git config --global user.name "Votre Nom"
git config --global user.email "votre.email@example.com"
```

## Étape 2 : Créer un compte GitHub (si vous n'en avez pas)

1. Allez sur https://github.com
2. Cliquez sur "Sign up"
3. Créez votre compte

## Étape 3 : Créer un Personal Access Token (pour l'automatisation)

### Option A : Créer un token pour automatiser la création du dépôt

1. Allez sur : https://github.com/settings/tokens
2. Cliquez sur "Generate new token (classic)"
3. Donnez un nom au token (ex: "Création dépôt ProgramPhasing")
4. Sélectionnez les permissions :
   - ✅ **repo** (toutes les cases)
5. Cliquez sur "Generate token"
6. **COPIEZ LE TOKEN** (vous ne pourrez plus le voir après !)

### Option B : Utiliser la création manuelle (plus simple)

Si vous préférez créer le dépôt manuellement, vous n'avez pas besoin de token.

## Étape 4 : Créer le dépôt GitHub

### Méthode 1 : Automatique (avec token)

1. **Exécutez le script PowerShell** :
   ```powershell
   powershell -ExecutionPolicy Bypass -File .\CreerDepotGitHub.ps1
   ```

2. **Suivez les instructions** :
   - Entrez le nom du dépôt (ex: "ProgramPhasing")
   - Entrez votre nom d'utilisateur GitHub
   - Entrez votre token GitHub (quand demandé)
   - Le script créera le dépôt et poussera le code automatiquement

### Méthode 2 : Manuelle (sans token)

1. **Créer le dépôt sur GitHub** :
   - Allez sur : https://github.com/new
   - Entrez le nom du dépôt (ex: "ProgramPhasing")
   - Ajoutez une description (optionnel)
   - Choisissez Public ou Private
   - **NE COCHEZ PAS** "Initialize this repository with a README"
   - Cliquez sur "Create repository"

2. **Exécutez le script PowerShell** :
   ```powershell
   powershell -ExecutionPolicy Bypass -File .\CreerDepotGitHub.ps1
   ```

3. **Quand demandé** :
   - Répondez "N" à la question sur le token
   - Répondez "O" quand le script demande si vous avez créé le dépôt
   - Entrez l'URL de votre dépôt (ex: `https://github.com/votre-username/ProgramPhasing.git`)

### Méthode 3 : Commandes manuelles

Si vous préférez faire tout manuellement :

```powershell
# 1. Initialiser Git
git init

# 2. Ajouter tous les fichiers
git add .

# 3. Faire le commit initial
git commit -m "Initial commit: Programme de phasing génétique avec générateur de données"

# 4. Renommer la branche en main
git branch -M main

# 5. Créer le dépôt sur GitHub (via le site web)
# Allez sur https://github.com/new et créez le dépôt

# 6. Ajouter le dépôt distant (remplacez par votre URL)
git remote add origin https://github.com/VOTRE_USERNAME/NOM_DEPOT.git

# 7. Pousser le code
git push -u origin main
```

## Étape 5 : Vérifier

Une fois le code poussé, allez sur votre dépôt GitHub :
- URL : `https://github.com/VOTRE_USERNAME/NOM_DEPOT`
- Vous devriez voir tous vos fichiers

## Dépannage

### Erreur : "Git n'est pas installé"
- Installez Git depuis https://git-scm.com/download/win
- Redémarrez PowerShell après l'installation

### Erreur : "Please tell me who you are"
- Configurez Git avec :
  ```powershell
  git config --global user.name "Votre Nom"
  git config --global user.email "votre.email@example.com"
  ```

### Erreur : "Authentication failed"
- Vérifiez vos identifiants GitHub
- Utilisez un Personal Access Token au lieu du mot de passe
- Pour créer un token : https://github.com/settings/tokens

### Erreur : "Repository not found"
- Vérifiez que le dépôt existe sur GitHub
- Vérifiez l'URL du dépôt distant : `git remote -v`
- Vérifiez que vous avez les permissions sur le dépôt

## Commandes Git utiles

```powershell
# Voir l'état du dépôt
git status

# Voir les fichiers modifiés
git diff

# Ajouter des fichiers modifiés
git add .
git commit -m "Description des modifications"
git push

# Voir l'historique
git log

# Voir les dépôts distants
git remote -v
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
├── GUIDE_GITHUB.md
├── ProgramPhasing.cpp
├── GenerateurDonneesAleatoires.cpp
├── GenererDonnees.ps1
├── CreerDepotGitHub.ps1
├── TransfererGitHub.ps1
├── Makefile
├── include/
│   └── ...
├── src/
│   └── ...
└── ...
```

**Note** : Le dossier `donnees_aleatoires/` ne sera pas inclus (voir `.gitignore`).

## Prochaines étapes

Une fois le projet sur GitHub :

1. **Ajoutez une description** sur la page du dépôt GitHub
2. **Ajoutez des topics** (tags) pour faciliter la recherche
3. **Créez un fichier LICENSE** si vous voulez partager le code
4. **Ajoutez des collaborateurs** si nécessaire
5. **Configurez GitHub Actions** pour l'intégration continue (optionnel)

## Support

Si vous rencontrez des problèmes :
- Consultez la documentation Git : https://git-scm.com/doc
- Consultez la documentation GitHub : https://docs.github.com
- Vérifiez le fichier `GUIDE_GITHUB.md` pour plus de détails

