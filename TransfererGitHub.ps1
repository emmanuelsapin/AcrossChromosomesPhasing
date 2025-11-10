# Script PowerShell pour transférer le projet sur GitHub
# Ce script automatise les étapes de base pour pousser le code sur GitHub

Write-Host "========================================"
Write-Host "Transfert du projet sur GitHub"
Write-Host "========================================"
Write-Host ""

# Vérifier si Git est installé
try {
    $gitVersion = git --version 2>&1
    Write-Host "[OK] Git est installé : $gitVersion" -ForegroundColor Green
} catch {
    Write-Host "[ERREUR] Git n'est pas installé !" -ForegroundColor Red
    Write-Host "Veuillez installer Git depuis : https://git-scm.com/download/win" -ForegroundColor Yellow
    Write-Host "Puis relancez ce script." -ForegroundColor Yellow
    exit 1
}

# Vérifier si un dépôt Git existe déjà
if (Test-Path .git) {
    Write-Host "[INFO] Un dépôt Git existe déjà dans ce répertoire" -ForegroundColor Yellow
    $response = Read-Host "Voulez-vous continuer ? (O/N)"
    if ($response -ne "O" -and $response -ne "o") {
        exit 0
    }
} else {
    Write-Host "[INFO] Initialisation du dépôt Git..." -ForegroundColor Cyan
    git init
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ERREUR] Échec de l'initialisation Git" -ForegroundColor Red
        exit 1
    }
    Write-Host "[OK] Dépôt Git initialisé" -ForegroundColor Green
}

# Vérifier si .gitignore existe
if (-not (Test-Path .gitignore)) {
    Write-Host "[ERREUR] Le fichier .gitignore n'existe pas !" -ForegroundColor Red
    Write-Host "Veuillez créer un fichier .gitignore avant de continuer." -ForegroundColor Yellow
    exit 1
} else {
    Write-Host "[OK] Fichier .gitignore trouvé" -ForegroundColor Green
}

# Ajouter tous les fichiers
Write-Host ""
Write-Host "[INFO] Ajout des fichiers au dépôt..." -ForegroundColor Cyan
git add .

# Afficher le statut
Write-Host ""
Write-Host "[INFO] Statut du dépôt :" -ForegroundColor Cyan
git status

# Demander confirmation pour le commit
Write-Host ""
$messageCommit = Read-Host "Entrez le message de commit (ou appuyez sur Entrée pour utiliser le message par défaut)"
if ([string]::IsNullOrWhiteSpace($messageCommit)) {
    $messageCommit = "Initial commit: Programme de phasing génétique avec générateur de données"
}

Write-Host ""
Write-Host "[INFO] Création du commit..." -ForegroundColor Cyan
git commit -m $messageCommit

if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERREUR] Échec du commit" -ForegroundColor Red
    Write-Host "Vérifiez que vous avez configuré Git avec :" -ForegroundColor Yellow
    Write-Host "  git config --global user.name 'Votre Nom'" -ForegroundColor Yellow
    Write-Host "  git config --global user.email 'votre.email@example.com'" -ForegroundColor Yellow
    exit 1
}

Write-Host "[OK] Commit créé avec succès" -ForegroundColor Green

# Demander l'URL du dépôt GitHub
Write-Host ""
Write-Host "========================================"
Write-Host "Configuration du dépôt distant GitHub"
Write-Host "========================================"
Write-Host ""
Write-Host "Pour continuer, vous devez avoir créé un dépôt sur GitHub." -ForegroundColor Yellow
Write-Host "Si vous ne l'avez pas encore fait :" -ForegroundColor Yellow
Write-Host "  1. Allez sur https://github.com" -ForegroundColor Yellow
Write-Host "  2. Créez un nouveau dépôt (sans initialiser avec README)" -ForegroundColor Yellow
Write-Host "  3. Copiez l'URL du dépôt" -ForegroundColor Yellow
Write-Host ""

$urlDepot = Read-Host "Entrez l'URL de votre dépôt GitHub (ex: https://github.com/username/repo.git)"

if ([string]::IsNullOrWhiteSpace($urlDepot)) {
    Write-Host "[INFO] Aucune URL fournie. Vous pouvez ajouter le dépôt distant plus tard avec :" -ForegroundColor Yellow
    Write-Host "  git remote add origin <URL>" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Pour pousser le code plus tard :" -ForegroundColor Yellow
    Write-Host "  git branch -M main" -ForegroundColor Yellow
    Write-Host "  git push -u origin main" -ForegroundColor Yellow
    exit 0
}

# Vérifier si le remote existe déjà
$remoteExists = git remote get-url origin 2>&1
if ($LASTEXITCODE -eq 0) {
    Write-Host "[INFO] Un dépôt distant 'origin' existe déjà : $remoteExists" -ForegroundColor Yellow
    $response = Read-Host "Voulez-vous le remplacer ? (O/N)"
    if ($response -eq "O" -or $response -eq "o") {
        git remote remove origin
        git remote add origin $urlDepot
        Write-Host "[OK] Dépôt distant mis à jour" -ForegroundColor Green
    }
} else {
    git remote add origin $urlDepot
    if ($LASTEXITCODE -eq 0) {
        Write-Host "[OK] Dépôt distant ajouté" -ForegroundColor Green
    } else {
        Write-Host "[ERREUR] Échec de l'ajout du dépôt distant" -ForegroundColor Red
        exit 1
    }
}

# Renommer la branche en main
Write-Host ""
Write-Host "[INFO] Configuration de la branche principale..." -ForegroundColor Cyan
git branch -M main

# Pousser le code
Write-Host ""
Write-Host "[INFO] Poussée du code sur GitHub..." -ForegroundColor Cyan
Write-Host "Vous devrez peut-être entrer vos identifiants GitHub." -ForegroundColor Yellow
Write-Host ""

git push -u origin main

if ($LASTEXITCODE -eq 0) {
    Write-Host ""
    Write-Host "========================================"
    Write-Host "[SUCCÈS] Projet transféré sur GitHub !" -ForegroundColor Green
    Write-Host "========================================"
    Write-Host ""
    Write-Host "Votre projet est maintenant disponible sur :" -ForegroundColor Cyan
    Write-Host $urlDepot -ForegroundColor White
} else {
    Write-Host ""
    Write-Host "[ERREUR] Échec de la poussée sur GitHub" -ForegroundColor Red
    Write-Host ""
    Write-Host "Solutions possibles :" -ForegroundColor Yellow
    Write-Host "  1. Vérifiez que le dépôt existe sur GitHub" -ForegroundColor Yellow
    Write-Host "  2. Vérifiez vos identifiants GitHub" -ForegroundColor Yellow
    Write-Host "  3. Utilisez un Personal Access Token si nécessaire" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Vous pouvez réessayer manuellement avec :" -ForegroundColor Yellow
    Write-Host "  git push -u origin main" -ForegroundColor White
}

