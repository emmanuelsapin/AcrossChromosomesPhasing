# Guide de Transfert vers GitHub - Dépôt Spécifique

## Dépôt cible
**https://github.com/emmanuelsapin/AcrossChromosomesPhasing**

Ce dépôt existe déjà mais est vide. Nous allons y transférer tout votre projet.

## Étapes à suivre

### Étape 1 : Installer Git

1. **Télécharger Git pour Windows**
   - Allez sur : https://git-scm.com/download/win
   - Téléchargez la dernière version
   - Exécutez l'installateur

2. **Configuration de l'installation**
   - Acceptez les options par défaut
   - Choisissez votre éditeur préféré
   - Sélectionnez "Git from the command line and also from 3rd-party software"
   - Acceptez les autres options par défaut

3. **Redémarrer PowerShell**
   - Fermez et rouvrez PowerShell après l'installation

### Étape 2 : Configurer Git (première fois)

Ouvrez PowerShell dans le répertoire du projet et exécutez :

```powershell
git config --global user.name "Emmanuel Sapin"
git config --global user.email "votre.email@example.com"
```

*(Remplacez par votre nom et email réels)*

### Étape 3 : Transférer le projet

Une fois Git installé et configuré, exécutez simplement :

```powershell
powershell -ExecutionPolicy Bypass -File .\TransfererVersGitHub.ps1
```

Le script va :
1. ✅ Vérifier que Git est installé
2. ✅ Initialiser le dépôt Git (si nécessaire)
3. ✅ Ajouter tous les fichiers
4. ✅ Créer un commit initial
5. ✅ Configurer le dépôt distant vers votre GitHub
6. ✅ Pousser le code sur GitHub

### Étape 4 : Authentification GitHub

Lors du `git push`, vous devrez vous authentifier. Deux options :

#### Option A : Personal Access Token (Recommandé)

1. **Créer un token** :
   - Allez sur : https://github.com/settings/tokens
   - Cliquez sur "Generate new token (classic)"
   - Donnez un nom (ex: "Transfert AcrossChromosomesPhasing")
   - Sélectionnez la permission **repo** (toutes les cases)
   - Cliquez sur "Generate token"
   - **COPIEZ LE TOKEN** (vous ne pourrez plus le voir !)

2. **Lors du push** :
   - Username : `emmanuelsapin`
   - Password : collez votre token

#### Option B : GitHub CLI (Alternative)

Si vous préférez utiliser GitHub CLI :

```powershell
# Installer GitHub CLI (optionnel)
# Puis :
gh auth login
```

### Étape 5 : Vérifier

Une fois le transfert réussi, allez sur :
**https://github.com/emmanuelsapin/AcrossChromosomesPhasing**

Vous devriez voir tous vos fichiers !

## Commandes manuelles (si le script ne fonctionne pas)

Si vous préférez faire les étapes manuellement :

```powershell
# 1. Initialiser Git
git init

# 2. Ajouter tous les fichiers
git add .

# 3. Faire le commit initial
git commit -m "Initial commit: Programme de phasing génétique avec générateur de données"

# 4. Renommer la branche en main
git branch -M main

# 5. Ajouter le dépôt distant
git remote add origin https://github.com/emmanuelsapin/AcrossChromosomesPhasing.git

# 6. Pousser le code
git push -u origin main
```

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
- Vérifiez que vous utilisez un Personal Access Token (pas votre mot de passe GitHub)
- Créez un nouveau token si nécessaire : https://github.com/settings/tokens

### Erreur : "Permission denied"
- Vérifiez que vous avez les droits d'écriture sur le dépôt
- Vérifiez que vous êtes connecté avec le bon compte GitHub (emmanuelsapin)

### Erreur : "Repository not found"
- Vérifiez que le dépôt existe : https://github.com/emmanuelsapin/AcrossChromosomesPhasing
- Vérifiez que vous avez les permissions sur ce dépôt

## Fichiers qui seront transférés

Le projet suivant sera transféré :

```
AcrossChromosomesPhasing/
├── .gitignore
├── README.md
├── README_Generateur.md
├── README_ARCHITECTURE.md
├── README_MODULAR.md
├── USER_MANUAL.md
├── GUIDE_GITHUB.md
├── INSTALLATION_GIT_ET_GITHUB.md
├── TRANSFERT_GITHUB_SPECIFIQUE.md
├── ProgramPhasing.cpp
├── GenerateurDonneesAleatoires.cpp
├── GenererDonnees.ps1
├── CreerDepotGitHub.ps1
├── TransfererGitHub.ps1
├── TransfererVersGitHub.ps1
├── Makefile
├── include/
│   └── ...
├── src/
│   └── ...
└── ...
```

**Note** : Le dossier `donnees_aleatoires/` ne sera pas inclus (voir `.gitignore`).

## Après le transfert

Une fois le code sur GitHub, vous pouvez :

1. **Ajouter une description** sur la page du dépôt
2. **Ajouter des topics** (tags) pour faciliter la recherche
3. **Créer un fichier LICENSE** si vous voulez partager le code
4. **Configurer GitHub Pages** pour la documentation (optionnel)
5. **Ajouter des collaborateurs** si nécessaire

## Support

Si vous rencontrez des problèmes :
- Consultez la documentation Git : https://git-scm.com/doc
- Consultez la documentation GitHub : https://docs.github.com
- Vérifiez les autres fichiers de documentation dans le projet

