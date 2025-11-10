# Script PowerShell pour transférer le projet vers GitHub avec le token fourni
# Dépôt : https://github.com/emmanuelsapin/AcrossChromosomesPhasing

$REPO_URL = "https://github.com/emmanuelsapin/AcrossChromosomesPhasing.git"
$USERNAME = "emmanuelsapin"
# Le token sera lu depuis une variable d'environnement ou demandé
$TOKEN = $env:GITHUB_TOKEN
if ([string]::IsNullOrWhiteSpace($TOKEN)) {
    # Si pas de variable d'environnement, utiliser le token fourni en paramètre
    if ($args.Count -gt 0) {
        $TOKEN = $args[0]
    } else {
        Write-Host "[INFO] Token GitHub non trouvé" -ForegroundColor Yellow
        Write-Host "Utilisez: `$env:GITHUB_TOKEN='votre_token' ; .\TransfererMaintenant.ps1" -ForegroundColor Yellow
        exit 1
    }
}

# Ajouter Git au PATH si nécessaire
$gitPaths = @(
    "C:\Program Files\Git\cmd",
    "C:\Program Files\Git\bin",
    "C:\Program Files (x86)\Git\cmd",
    "C:\Program Files (x86)\Git\bin"
)

foreach ($path in $gitPaths) {
    if (Test-Path $path) {
        if ($env:Path -notlike "*$path*") {
            $env:Path += ";$path"
        }
    }
}

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
    Write-Host "  git config --global user.name 'Emmanuel Sapin'" -ForegroundColor White
    Write-Host "  git config --global user.email 'votre.email@example.com'" -ForegroundColor White
    exit 1
}

# Vérifier la configuration Git
$userName = git config --global user.name 2>&1
$userEmail = git config --global user.email 2>&1

if ([string]::IsNullOrWhiteSpace($userName) -or [string]::IsNullOrWhiteSpace($userEmail)) {
    Write-Host "[INFO] Configuration automatique de Git..." -ForegroundColor Cyan
    git config --global user.name "Emmanuel Sapin"
    git config --global user.email "emmanuelsapin@users.noreply.github.com"
    Write-Host "[OK] Git configuré automatiquement" -ForegroundColor Green
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

# Afficher le statut
Write-Host ""
Write-Host "[INFO] Statut du dépôt :" -ForegroundColor Cyan
git status --short

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
Write-Host "[INFO] Branche actuelle : $currentBranch" -ForegroundColor Cyan
if ($currentBranch -ne "main") {
    Write-Host "[INFO] Renommage de la branche en 'main'..." -ForegroundColor Cyan
    git branch -M main
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ATTENTION] Impossible de renommer la branche, utilisation de la branche actuelle" -ForegroundColor Yellow
        $currentBranch = git branch --show-current 2>&1
    } else {
        $currentBranch = "main"
    }
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

# Pousser le code
Write-Host ""
Write-Host "[INFO] Poussée du code sur GitHub..." -ForegroundColor Cyan
Write-Host "Utilisation du token fourni..." -ForegroundColor Yellow

# Utiliser un script credential helper personnalisé
$credentialScript = @"
#!/bin/sh
echo username=$USERNAME
echo password=$TOKEN
"@

$scriptPath = "$env:TEMP\git-credential-helper.sh"
$credentialScript | Out-File -FilePath $scriptPath -Encoding ASCII -Force

# Configurer git pour utiliser ce script
git config --global credential.helper "!f() { echo 'username=$USERNAME'; echo 'password=$TOKEN'; }; f"

# Configurer le remote (sans token dans l'URL)
git remote set-url origin $REPO_URL

# Désactiver les prompts Git
$env:GIT_TERMINAL_PROMPT = "0"
$env:GIT_ASKPASS = "echo"

# Pousser avec la branche actuelle
$branchToPush = git branch --show-current 2>&1
Write-Host "[INFO] Poussée de la branche : $branchToPush" -ForegroundColor Cyan

# Pousser vers origin
git push -u origin $branchToPush

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
    Write-Host "  1. Que le token est correct et n'a pas expiré" -ForegroundColor Yellow
    Write-Host "  2. Que le dépôt existe et que vous avez les permissions" -ForegroundColor Yellow
    Write-Host "  3. Que vous êtes connecté avec le bon compte GitHub" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Vous pouvez réessayer manuellement avec :" -ForegroundColor Yellow
    Write-Host "  git push -u origin main" -ForegroundColor White
    Write-Host "  (utilisez votre token comme mot de passe)" -ForegroundColor White
}

# Nettoyer le token de la mémoire
$TOKEN = $null
[System.GC]::Collect()

