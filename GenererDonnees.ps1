# Script PowerShell pour générer les données aléatoires
# Générateur de données HAP avec 1000 individus et 1000 SNPs par chromosome

$NB_INDIVIDUS = 1000
$NB_SNPS_PAR_CHROMOSOME = 1000
$CHROMOSOME_MIN = 1
$CHROMOSOME_MAX = 22
$POSITION_MIN = 1000
$POSITION_MAX = 250000000

$ALLELES = @('A', 'T', 'G', 'C')

# Fonction pour générer un rsID
function Generer-RsID {
    param($chromosome, $snpIndex)
    $rsNumber = $chromosome * 100000 + $snpIndex
    return "rs$($rsNumber.ToString('000000'))"
}

# Fonction pour générer deux allèles différents
function Generer-Alleles {
    $allele0 = Get-Random -InputObject $ALLELES
    $allele1 = Get-Random -InputObject $ALLELES
    while ($allele1 -eq $allele0) {
        $allele1 = Get-Random -InputObject $ALLELES
    }
    return $allele0, $allele1
}

# Fonction pour générer un génotype selon Hardy-Weinberg
function Generer-Genotype {
    param($maf)
    $random = Get-Random -Minimum 0.0 -Maximum 1.0
    $p = 1.0 - $maf
    $p_homo_majeur = $p * $p
    $p_heterozygote = 2.0 * $p * $maf
    
    if ($random -lt $p_homo_majeur) {
        return 0  # Homozygote référence
    } elseif ($random -lt ($p_homo_majeur + $p_heterozygote)) {
        return 1  # Hétérozygote
    } else {
        return 2  # Homozygote alternatif
    }
}

# Fonction pour convertir génotype en allèles
function Genotype-VersAlleles {
    param($genotype)
    switch ($genotype) {
        0 { return 0, 0 }
        1 { return 0, 1 }
        2 { return 1, 1 }
        default { return 0, 0 }
    }
}

# Fonction pour générer un chromosome
function Generer-Chromosome {
    param($chromosome, $cheminSortie, $genererMAF = $true)
    
    $nomFichierHAP = Join-Path $cheminSortie "$chromosome.hap"
    $mafValues = @()
    
    Write-Host "Génération du chromosome $chromosome..."
    
    $file = [System.IO.StreamWriter]::new($nomFichierHAP)
    
    for ($snp = 0; $snp -lt $NB_SNPS_PAR_CHROMOSOME; $snp++) {
        $rsID = Generer-RsID $chromosome $snp
        $position = Get-Random -Minimum $POSITION_MIN -Maximum $POSITION_MAX
        $allele0, $allele1 = Generer-Alleles
        
        # Générer MAF entre 0.01 et 0.5
        $maf = 0.01 + (Get-Random -Minimum 0.0 -Maximum 1.0) * 0.49
        $mafValues += $maf
        
        # Écrire l'en-tête du SNP
        $file.Write("$rsID $position $allele0 $allele1")
        
        # Générer les génotypes pour chaque individu
        for ($indiv = 0; $indiv -lt $NB_INDIVIDUS; $indiv++) {
            $genotype = Generer-Genotype $maf
            $a0, $a1 = Genotype-VersAlleles $genotype
            $file.Write(" $a0 $a1")
        }
        
        $file.WriteLine()
        
        if (($snp + 1) % 100 -eq 0) {
            Write-Host "  SNP $($snp + 1)/$NB_SNPS_PAR_CHROMOSOME"
        }
    }
    
    $file.Close()
    Write-Host "Fichier HAP créé : $nomFichierHAP"
    
    # Générer le fichier MAF
    if ($genererMAF) {
        $nomFichierMAF = Join-Path $cheminSortie "MAF.txt"
        $mafFile = [System.IO.StreamWriter]::new($nomFichierMAF, $true)  # Append mode
        for ($snp = 0; $snp -lt $NB_SNPS_PAR_CHROMOSOME; $snp++) {
            $mafFile.WriteLine("$chromosome $snp $($mafValues[$snp].ToString('F6'))")
        }
        $mafFile.Close()
    }
}

# Fonction principale
function Main {
    Write-Host "========================================"
    Write-Host "Générateur de Données Aléatoires"
    Write-Host "========================================"
    Write-Host "Configuration :"
    Write-Host "  - Nombre d'individus : $NB_INDIVIDUS"
    Write-Host "  - SNPs par chromosome : $NB_SNPS_PAR_CHROMOSOME"
    Write-Host "  - Chromosomes : $CHROMOSOME_MIN-$CHROMOSOME_MAX"
    Write-Host "========================================"
    
    $cheminSortie = ".\donnees_aleatoires\"
    if ($args.Count -gt 0) {
        $cheminSortie = $args[0]
        if (-not $cheminSortie.EndsWith('\') -and -not $cheminSortie.EndsWith('/')) {
            $cheminSortie += '\'
        }
    }
    
    # Créer le répertoire de sortie
    New-Item -ItemType Directory -Force -Path $cheminSortie | Out-Null
    
    # Supprimer le fichier MAF existant s'il existe
    $nomFichierMAF = Join-Path $cheminSortie "MAF.txt"
    if (Test-Path $nomFichierMAF) {
        Remove-Item $nomFichierMAF
    }
    
    Write-Host "`nDébut de la génération..."
    
    for ($chr = $CHROMOSOME_MIN; $chr -le $CHROMOSOME_MAX; $chr++) {
        Generer-Chromosome $chr $cheminSortie $true
    }
    
    Write-Host "`n========================================"
    Write-Host "Génération terminée avec succès !"
    Write-Host "Fichiers créés dans : $cheminSortie"
    Write-Host "  - 1.hap à 22.hap (fichiers HAP)"
    Write-Host "  - MAF.txt (fréquences alléliques mineures)"
    Write-Host "========================================"
}

Main

