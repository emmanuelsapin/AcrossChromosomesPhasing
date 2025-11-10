# Script PowerShell pour créer un dépôt GitHub et transférer le projet
# Ce script peut créer le dépôt via l'API GitHub ou guider la création manuelle

param(
    [string]$NomDepot = "",
    [string]$Description = "Programme de phasing génétique avec générateur de données",
    [string]$TokenGitHub = "",
    [string]$Username = "",
    [switch]$Prive = $false
)

Write-Host "========================================"
Write-Host "Création d'un dépôt GitHub"
Write-Host "========================================"
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

# Demander le nom du dépôt si non fourni
if ([string]::IsNullOrWhiteSpace($NomDepot)) {
    $NomDepot = Read-Host "Entrez le nom du dépôt GitHub (ex: ProgramPhasing)"
}

if ([string]::IsNullOrWhiteSpace($NomDepot)) {
    Write-Host "[ERREUR] Le nom du dépôt est requis" -ForegroundColor Red
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
    exit 1
}

# Ajouter tous les fichiers
Write-Host ""
Write-Host "[INFO] Ajout des fichiers au dépôt..." -ForegroundColor Cyan
git add .

# Faire le commit initial
$messageCommit = "Initial commit: Programme de phasing génétique avec générateur de données"
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

# Renommer la branche en main
git branch -M main

# Essayer de créer le dépôt via l'API GitHub
$creerViaAPI = $false
if (-not [string]::IsNullOrWhiteSpace($TokenGitHub)) {
    $creerViaAPI = $true
} else {
    Write-Host ""
    Write-Host "========================================"
    Write-Host "Création du dépôt GitHub"
    Write-Host "========================================"
    Write-Host ""
    Write-Host "Pour créer le dépôt automatiquement, vous avez besoin d'un Personal Access Token GitHub." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Pour créer un token :" -ForegroundColor Cyan
    Write-Host "  1. Allez sur https://github.com/settings/tokens" -ForegroundColor White
    Write-Host "  2. Cliquez sur 'Generate new token (classic)'" -ForegroundColor White
    Write-Host "  3. Donnez un nom (ex: 'Création dépôt ProgramPhasing')" -ForegroundColor White
    Write-Host "  4. Sélectionnez la permission 'repo' (toutes les cases)" -ForegroundColor White
    Write-Host "  5. Cliquez sur 'Generate token'" -ForegroundColor White
    Write-Host "  6. Copiez le token (vous ne pourrez plus le voir après)" -ForegroundColor White
    Write-Host ""
    
    $response = Read-Host "Avez-vous un token GitHub ? (O/N)"
    if ($response -eq "O" -or $response -eq "o") {
        $TokenGitHub = Read-Host "Entrez votre token GitHub" -AsSecureString
        $BSTR = [System.Runtime.InteropServices.Marshal]::SecureStringToBSTR($TokenGitHub)
        $TokenGitHub = [System.Runtime.InteropServices.Marshal]::PtrToStringAuto($BSTR)
        $creerViaAPI = $true
    }
}

# Créer le dépôt via l'API GitHub
if ($creerViaAPI) {
    Write-Host ""
    Write-Host "[INFO] Création du dépôt via l'API GitHub..." -ForegroundColor Cyan
    
    # Demander le username si non fourni
    if ([string]::IsNullOrWhiteSpace($Username)) {
        $Username = Read-Host "Entrez votre nom d'utilisateur GitHub"
    }
    
    if ([string]::IsNullOrWhiteSpace($Username)) {
        Write-Host "[ERREUR] Le nom d'utilisateur GitHub est requis" -ForegroundColor Red
        $creerViaAPI = $false
    } else {
        # Préparer les données JSON
        $body = @{
            name = $NomDepot
            description = $Description
            private = $Prive.IsPresent
            auto_init = $false
        } | ConvertTo-Json
        
        # Headers avec authentification
        $headers = @{
            "Authorization" = "token $TokenGitHub"
            "Accept" = "application/vnd.github.v3+json"
        }
        
        # Créer le dépôt
        try {
            $response = Invoke-RestMethod -Uri "https://api.github.com/user/repos" -Method Post -Headers $headers -Body $body -ContentType "application/json"
            $urlDepot = $response.clone_url
            
            Write-Host "[OK] Dépôt créé avec succès : $urlDepot" -ForegroundColor Green
            
            # Ajouter le remote et pousser
            git remote add origin $urlDepot
            Write-Host "[OK] Dépôt distant ajouté" -ForegroundColor Green
            
            Write-Host ""
            Write-Host "[INFO] Poussée du code sur GitHub..." -ForegroundColor Cyan
            git push -u origin main
            
            if ($LASTEXITCODE -eq 0) {
                Write-Host ""
                Write-Host "========================================"
                Write-Host "[SUCCÈS] Projet transféré sur GitHub !" -ForegroundColor Green
                Write-Host "========================================"
                Write-Host ""
                Write-Host "Votre projet est maintenant disponible sur :" -ForegroundColor Cyan
                Write-Host $response.html_url -ForegroundColor White
                Write-Host ""
                exit 0
            } else {
                Write-Host "[ERREUR] Échec de la poussée sur GitHub" -ForegroundColor Red
            }
        } catch {
            Write-Host "[ERREUR] Échec de la création du dépôt via l'API" -ForegroundColor Red
            Write-Host "Erreur : $($_.Exception.Message)" -ForegroundColor Yellow
            $creerViaAPI = $false
        }
    }
}

# Si la création via API a échoué, guider l'utilisateur pour la création manuelle
if (-not $creerViaAPI) {
    Write-Host ""
    Write-Host "========================================"
    Write-Host "Création manuelle du dépôt GitHub"
    Write-Host "========================================"
    Write-Host ""
    Write-Host "Veuillez créer le dépôt manuellement sur GitHub :" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "1. Allez sur https://github.com/new" -ForegroundColor Cyan
    Write-Host "2. Entrez le nom du dépôt : $NomDepot" -ForegroundColor White
    Write-Host "3. Ajoutez une description (optionnel)" -ForegroundColor White
    Write-Host "4. Choisissez Public ou Private" -ForegroundColor White
    Write-Host "5. NE COCHEZ PAS 'Initialize this repository with a README'" -ForegroundColor Yellow
    Write-Host "6. Cliquez sur 'Create repository'" -ForegroundColor White
    Write-Host ""
    
    $continuer = Read-Host "Avez-vous créé le dépôt sur GitHub ? (O/N)"
    if ($continuer -eq "O" -or $continuer -eq "o") {
        # Demander l'URL du dépôt
        Write-Host ""
        $urlDepot = Read-Host "Entrez l'URL de votre dépôt GitHub (ex: https://github.com/username/repo.git)"
        
        if (-not [string]::IsNullOrWhiteSpace($urlDepot)) {
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
                Write-Host "[OK] Dépôt distant ajouté" -ForegroundColor Green
            }
            
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
        } else {
            Write-Host "[INFO] Aucune URL fournie. Vous pouvez ajouter le dépôt distant plus tard avec :" -ForegroundColor Yellow
            Write-Host "  git remote add origin <URL>" -ForegroundColor White
            Write-Host "  git push -u origin main" -ForegroundColor White
        }
    } else {
        Write-Host "[INFO] Vous pouvez créer le dépôt plus tard et utiliser :" -ForegroundColor Yellow
        Write-Host "  git remote add origin <URL>" -ForegroundColor White
        Write-Host "  git push -u origin main" -ForegroundColor White
    }
}

