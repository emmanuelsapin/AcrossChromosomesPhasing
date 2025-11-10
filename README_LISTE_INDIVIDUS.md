# Fonctionnalité : Liste d'Individus pour le Phasing

## Description

Cette fonctionnalité permet de spécifier une liste d'individus sur lesquels appliquer la méthode de phasing, au lieu de traiter tous les individus du dataset.

## Utilisation

### Format du fichier de liste

Le fichier de liste doit contenir un ID d'individu par ligne. Les IDs sont indexés à partir de 0.

**Format :**
```
0
5
10
15
20
```

**Caractéristiques :**
- Un ID par ligne
- Les lignes vides sont ignorées
- Les commentaires (lignes commençant par `#`) sont ignorés
- Les IDs doivent être entre 0 et (NbIndiv - 1)

**Exemple de fichier :**
```
# Liste d'individus à traiter
# IDs : 0, 5, 10, 15, 20

0
5
10
15
20
```

### Commande

Ajoutez le paramètre `-ListIndiv` suivi du chemin vers le fichier de liste :

```bash
./ProgramPhasing \
  -NbIndiv 1000 \
  -PathInput ./donnees_aleatoires/ \
  -PathOutput ./resultats/ \
  -ListIndiv ./liste_individus.txt
```

### Exemple complet

1. **Créer un fichier de liste** (`liste_individus.txt`) :
```
0
10
20
30
40
50
```

2. **Exécuter le programme** :
```bash
./ProgramPhasing \
  -NbIndiv 1000 \
  -PathInput ./donnees_aleatoires/ \
  -PathOutput ./resultats/ \
  -ListIndiv ./liste_individus.txt
```

3. **Résultat** :
   - Seuls les individus 0, 10, 20, 30, 40, et 50 seront traités
   - Les fichiers de sortie ne contiendront que ces individus
   - Le temps de traitement sera réduit proportionnellement

## Comportement

### Avec `-ListIndiv` spécifié

- Seuls les individus listés dans le fichier sont traités
- Les fichiers de sortie ne contiennent que ces individus
- Les IDs invalides (hors de la plage [0, NbIndiv)) sont ignorés avec un avertissement
- Si aucun ID valide n'est trouvé, le programme s'arrête avec une erreur

### Sans `-ListIndiv` spécifié

- Comportement par défaut : tous les individus sont traités
- Aucun changement par rapport à la version précédente

## Validation

Le programme valide automatiquement :
- ✅ Que le fichier existe et peut être lu
- ✅ Que les IDs sont dans la plage valide [0, NbIndiv)
- ✅ Qu'au moins un ID valide est présent dans le fichier

## Messages d'erreur

### Fichier introuvable
```
ERROR: Could not open file ./liste_individus.txt for reading individual list
ERROR: Failed to read individual list file. Exiting.
```

### ID invalide
```
WARNING: Individual ID 5000 on line 5 is out of range [0, 1000). Skipping.
```

### Aucun ID valide
```
ERROR: No valid individual IDs found in file ./liste_individus.txt
ERROR: Failed to read individual list file. Exiting.
```

## Avantages

1. **Performance** : Traiter seulement un sous-ensemble d'individus réduit le temps de calcul
2. **Flexibilité** : Permet de tester sur un petit échantillon avant de traiter tout le dataset
3. **Reproductibilité** : Facilite la reproduction d'analyses sur des individus spécifiques
4. **Debugging** : Utile pour déboguer sur un petit nombre d'individus

## Notes techniques

- La fonction `readListIndiv()` lit le fichier et stocke les IDs dans un vecteur
- La boucle principale de traitement est modifiée pour itérer uniquement sur les individus de la liste
- La fonction `writeoutput()` est modifiée pour n'écrire que les individus de la liste dans les fichiers de sortie
- Les données génomiques de tous les individus sont toujours chargées en mémoire (nécessaire pour les calculs de relations)

## Exemple de fichier de liste

Voir `exemple_liste_individus.txt` pour un exemple complet.

