# Script PowerShell pour transférer le projet vers GitHub avec authentification par token
# IMPORTANT: GitHub n'accepte plus les mots de passe - il faut utiliser un Personal Access Token

$REPO_URL = "https://github.com/emmanuelsapin/AcrossChromosomesPhasing.git"
$USERNAME = "emmanuelsapin"

Write-Host "========================================"
Write-Host "Transfert vers GitHub"
Write-Host "Dépôt : $REPO_URL"
Write-Host "========================================"
Write-Host ""
Write-Host "IMPORTANT: GitHub nécessite un Personal Access Token" -ForegroundColor Yellow
Write-Host "Les mots de passe ne fonctionnent plus depuis août 2021" -ForegroundColor Yellow
Write-Host ""

# Vérifier si Git est installé
try {
    $gitVersion = git --version 2>&1
    Write-Host "[OK] Git est installé : $gitVersion" -ForegroundColor Green
} catch {
    Write-Host "[ERREUR] Git n'est pas installé !" -ForegroundColor Red
    Write-Host "Veuillez installer Git depuis : https://git-scm.com/download/win" -ForegroundColor Yellow
    exit 1
}

# Vérifier la configuration Git
$userName = git config --global user.name 2>&1
$userEmail = git config --global user.email 2>&1

if ([string]::IsNullOrWhiteSpace($userName) -or [string]::IsNullOrWhiteSpace($userEmail)) {
    Write-Host "[ATTENTION] Git n'est pas configuré !" -ForegroundColor Yellow
    Write-Host ""
    $configurer = Read-Host "Voulez-vous configurer Git maintenant ? (O/N)"
    if ($configurer -eq "O" -or $configurer -eq "o") {
        $nom = Read-Host "Entrez votre nom"
        $email = Read-Host "Entrez votre email"
        git config --global user.name $nom
        git config --global user.email $email
        Write-Host "[OK] Git configuré" -ForegroundColor Green
    } else {
        Write-Host "[ERREUR] Git doit être configuré pour continuer" -ForegroundColor Red
        exit 1
    }
}

# Initialiser Git si nécessaire
if (-not (Test-Path .git)) {
    Write-Host "[INFO] Initialisation du dépôt Git..." -ForegroundColor Cyan
    git init
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ERREUR] Échec de l'initialisation Git" -ForegroundColor Red
        exit 1
    }
    Write-Host "[OK] Dépôt Git initialisé" -ForegroundColor Green
} else {
    Write-Host "[INFO] Dépôt Git déjà initialisé" -ForegroundColor Cyan
}

# Vérifier si .gitignore existe
if (-not (Test-Path .gitignore)) {
    Write-Host "[ATTENTION] Le fichier .gitignore n'existe pas !" -ForegroundColor Yellow
    Write-Host "Création d'un .gitignore par défaut..." -ForegroundColor Yellow
    @"
# Fichiers compilés
*.exe
*.o
*.obj
*.out
ProgramPhasing
GenerateurDonneesAleatoires

# Données générées
donnees_aleatoires/
*.hap
MAF.txt

# Fichiers temporaires
*.tmp
*.swp
.DS_Store
"@ | Out-File -FilePath .gitignore -Encoding UTF8
    Write-Host "[OK] .gitignore créé" -ForegroundColor Green
}

# Ajouter tous les fichiers
Write-Host ""
Write-Host "[INFO] Ajout des fichiers au dépôt..." -ForegroundColor Cyan
git add .

# Vérifier s'il y a des changements à committer
$hasChanges = git diff --cached --quiet 2>&1
if ($LASTEXITCODE -ne 0) {
    $messageCommit = "Initial commit: Programme de phasing génétique avec générateur de données"
    Write-Host ""
    Write-Host "[INFO] Création du commit..." -ForegroundColor Cyan
    git commit -m $messageCommit
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ERREUR] Échec du commit" -ForegroundColor Red
        exit 1
    }
    Write-Host "[OK] Commit créé avec succès" -ForegroundColor Green
} else {
    Write-Host "[INFO] Aucun changement à committer" -ForegroundColor Yellow
}

# Renommer la branche en main si nécessaire
$currentBranch = git branch --show-current 2>&1
if ($currentBranch -ne "main" -and $currentBranch -ne "master") {
    Write-Host "[INFO] Renommage de la branche en 'main'..." -ForegroundColor Cyan
    git branch -M main
}

# Configurer le remote
$remoteExists = git remote get-url origin 2>&1
if ($LASTEXITCODE -eq 0) {
    if ($remoteExists -ne $REPO_URL) {
        Write-Host "[INFO] Mise à jour du dépôt distant..." -ForegroundColor Cyan
        git remote remove origin
        git remote add origin $REPO_URL
    }
} else {
    git remote add origin $REPO_URL
}
Write-Host "[OK] Dépôt distant configuré : $REPO_URL" -ForegroundColor Green

# Demander le token
Write-Host ""
Write-Host "========================================"
Write-Host "Authentification GitHub"
Write-Host "========================================"
Write-Host ""
Write-Host "Pour créer un Personal Access Token :" -ForegroundColor Cyan
Write-Host "  1. Allez sur : https://github.com/settings/tokens" -ForegroundColor White
Write-Host "  2. Cliquez sur 'Generate new token (classic)'" -ForegroundColor White
Write-Host "  3. Donnez un nom (ex: 'Transfert AcrossChromosomesPhasing')" -ForegroundColor White
Write-Host "  4. Sélectionnez la permission 'repo' (toutes les cases)" -ForegroundColor White
Write-Host "  5. Cliquez sur 'Generate token'" -ForegroundColor White
Write-Host "  6. COPIEZ LE TOKEN (vous ne pourrez plus le voir après !)" -ForegroundColor Yellow
Write-Host ""

$token = Read-Host "Entrez votre Personal Access Token GitHub" -AsSecureString
$BSTR = [System.Runtime.InteropServices.Marshal]::SecureStringToBSTR($token)
$tokenPlain = [System.Runtime.InteropServices.Marshal]::PtrToStringAuto($BSTR)

if ([string]::IsNullOrWhiteSpace($tokenPlain)) {
    Write-Host "[ERREUR] Le token est requis" -ForegroundColor Red
    exit 1
}

# Configurer l'URL avec le token pour le push
$repoWithToken = $REPO_URL -replace "https://", "https://$USERNAME`:$tokenPlain@"
$repoWithToken = $repoWithToken -replace "github.com", "github.com"

# Pousser le code
Write-Host ""
Write-Host "[INFO] Poussée du code sur GitHub..." -ForegroundColor Cyan

# Utiliser l'URL avec token pour le push
$env:GIT_TERMINAL_PROMPT = "0"
git push -u $repoWithToken main

if ($LASTEXITCODE -eq 0) {
    Write-Host ""
    Write-Host "========================================"
    Write-Host "[SUCCÈS] Projet transféré sur GitHub !" -ForegroundColor Green
    Write-Host "========================================"
    Write-Host ""
    Write-Host "Votre projet est maintenant disponible sur :" -ForegroundColor Cyan
    Write-Host "https://github.com/emmanuelsapin/AcrossChromosomesPhasing" -ForegroundColor White
    Write-Host ""
} else {
    Write-Host ""
    Write-Host "[ERREUR] Échec de la poussée sur GitHub" -ForegroundColor Red
    Write-Host ""
    Write-Host "Vérifiez :" -ForegroundColor Yellow
    Write-Host "  1. Que le token est correct et a la permission 'repo'" -ForegroundColor Yellow
    Write-Host "  2. Que le dépôt existe et que vous avez les permissions" -ForegroundColor Yellow
    Write-Host "  3. Que vous êtes connecté avec le bon compte GitHub" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Vous pouvez réessayer manuellement avec :" -ForegroundColor Yellow
    Write-Host "  git push -u origin main" -ForegroundColor White
    Write-Host "  (utilisez votre token comme mot de passe)" -ForegroundColor White
}

# Nettoyer le token de la mémoire
$tokenPlain = $null
[System.GC]::Collect()

