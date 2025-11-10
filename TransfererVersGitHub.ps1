# Script PowerShell pour transférer le projet vers https://github.com/emmanuelsapin/AcrossChromosomesPhasing

$REPO_URL = "https://github.com/emmanuelsapin/AcrossChromosomesPhasing.git"

Write-Host "========================================"
Write-Host "Transfert vers GitHub"
Write-Host "Dépôt : $REPO_URL"
Write-Host "========================================"
Write-Host ""

# Vérifier si Git est installé
try {
    $gitVersion = git --version 2>&1
    Write-Host "[OK] Git est installé : $gitVersion" -ForegroundColor Green
} catch {
    Write-Host "[ERREUR] Git n'est pas installé !" -ForegroundColor Red
    Write-Host ""
    Write-Host "Veuillez installer Git depuis : https://git-scm.com/download/win" -ForegroundColor Yellow
    Write-Host "Puis relancez ce script." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Après l'installation, configurez Git avec :" -ForegroundColor Cyan
    Write-Host "  git config --global user.name 'Votre Nom'" -ForegroundColor White
    Write-Host "  git config --global user.email 'votre.email@example.com'" -ForegroundColor White
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
    # Créer un .gitignore basique
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

# Afficher le statut
Write-Host ""
Write-Host "[INFO] Statut du dépôt :" -ForegroundColor Cyan
git status --short

# Vérifier s'il y a des changements à committer
$hasChanges = git diff --cached --quiet 2>&1
if ($LASTEXITCODE -ne 0) {
    # Il y a des changements à committer
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

# Vérifier si le remote existe déjà
$remoteExists = git remote get-url origin 2>&1
if ($LASTEXITCODE -eq 0) {
    if ($remoteExists -ne $REPO_URL) {
        Write-Host "[INFO] Un dépôt distant 'origin' existe déjà : $remoteExists" -ForegroundColor Yellow
        $response = Read-Host "Voulez-vous le remplacer par $REPO_URL ? (O/N)"
        if ($response -eq "O" -or $response -eq "o") {
            git remote remove origin
            git remote add origin $REPO_URL
            Write-Host "[OK] Dépôt distant mis à jour" -ForegroundColor Green
        } else {
            Write-Host "[INFO] Conservation du dépôt distant existant" -ForegroundColor Yellow
            $REPO_URL = $remoteExists
        }
    } else {
        Write-Host "[OK] Dépôt distant déjà configuré : $REPO_URL" -ForegroundColor Green
    }
} else {
    git remote add origin $REPO_URL
    if ($LASTEXITCODE -eq 0) {
        Write-Host "[OK] Dépôt distant ajouté : $REPO_URL" -ForegroundColor Green
    } else {
        Write-Host "[ERREUR] Échec de l'ajout du dépôt distant" -ForegroundColor Red
        exit 1
    }
}

# Pousser le code
Write-Host ""
Write-Host "[INFO] Poussée du code sur GitHub..." -ForegroundColor Cyan
Write-Host "Vous devrez peut-être entrer vos identifiants GitHub." -ForegroundColor Yellow
Write-Host "Si vous utilisez l'authentification par token, utilisez votre nom d'utilisateur GitHub" -ForegroundColor Yellow
Write-Host "et votre Personal Access Token comme mot de passe." -ForegroundColor Yellow
Write-Host ""

# Essayer de pousser
git push -u origin main

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
    Write-Host "Solutions possibles :" -ForegroundColor Yellow
    Write-Host "  1. Vérifiez que vous avez les permissions d'écriture sur le dépôt" -ForegroundColor Yellow
    Write-Host "  2. Vérifiez vos identifiants GitHub" -ForegroundColor Yellow
    Write-Host "  3. Utilisez un Personal Access Token :" -ForegroundColor Yellow
    Write-Host "     - Allez sur https://github.com/settings/tokens" -ForegroundColor White
    Write-Host "     - Créez un token avec la permission 'repo'" -ForegroundColor White
    Write-Host "     - Utilisez le token comme mot de passe lors du push" -ForegroundColor White
    Write-Host ""
    Write-Host "Vous pouvez réessayer manuellement avec :" -ForegroundColor Yellow
    Write-Host "  git push -u origin main" -ForegroundColor White
    Write-Host ""
    exit 1
}

