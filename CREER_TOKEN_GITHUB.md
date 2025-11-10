# Comment créer un Personal Access Token GitHub

## ⚠️ IMPORTANT : Sécurité

**GitHub n'accepte plus les mots de passe pour l'authentification Git depuis août 2021.**

Vous devez utiliser un **Personal Access Token** (PAT) à la place.

## Étapes pour créer un token

### 1. Accéder aux paramètres des tokens

1. Allez sur GitHub : https://github.com
2. Connectez-vous avec votre compte (`emmanuelsapin`)
3. Cliquez sur votre photo de profil (en haut à droite)
4. Cliquez sur **Settings**
5. Dans le menu de gauche, allez dans **Developer settings**
6. Cliquez sur **Personal access tokens**
7. Cliquez sur **Tokens (classic)**

### 2. Créer un nouveau token

1. Cliquez sur **Generate new token**
2. Sélectionnez **Generate new token (classic)**

### 3. Configurer le token

1. **Note** : Donnez un nom descriptif
   - Exemple : `Transfert AcrossChromosomesPhasing`
   - Ou : `Accès Git - ProgramPhasing`

2. **Expiration** : Choisissez une durée
   - Pour un usage ponctuel : 30 jours
   - Pour un usage régulier : 90 jours ou No expiration

3. **Permissions** : Sélectionnez les permissions nécessaires
   - ✅ **repo** (toutes les cases)
     - Cela donne accès aux dépôts privés et publics
     - Permet de lire et écrire

4. Cliquez sur **Generate token** (en bas de la page)

### 4. Copier le token

⚠️ **ATTENTION** : Le token s'affiche UNE SEULE FOIS !

1. **COPIEZ LE TOKEN IMMÉDIATEMENT**
2. Collez-le dans un endroit sûr (temporairement)
3. Vous ne pourrez plus le voir après avoir quitté la page

Le token ressemble à : `ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx`

## Utiliser le token

### Option 1 : Avec le script PowerShell

1. Exécutez le script :
   ```powershell
   powershell -ExecutionPolicy Bypass -File .\TransfererAvecToken.ps1
   ```

2. Quand le script demande le token, collez votre token

### Option 2 : Manuellement avec Git

Lors d'un `git push`, utilisez :
- **Username** : `emmanuelsapin`
- **Password** : collez votre token (pas votre mot de passe GitHub)

### Option 3 : Stocker le token (optionnel)

Pour éviter de le retaper à chaque fois :

```powershell
# Windows Credential Manager (recommandé)
git config --global credential.helper wincred

# Puis lors du premier push, entrez le token
# Il sera stocké de manière sécurisée
```

## Sécurité du token

### ⚠️ Bonnes pratiques

1. **Ne partagez jamais votre token**
   - Ne le commitez pas dans Git
   - Ne le partagez pas publiquement
   - Ne l'envoyez pas par email

2. **Utilisez des tokens avec des permissions minimales**
   - Donnez seulement les permissions nécessaires
   - Pour Git : seulement `repo`

3. **Révoquez les tokens inutilisés**
   - Allez sur https://github.com/settings/tokens
   - Supprimez les tokens que vous n'utilisez plus

4. **Utilisez des tokens avec expiration**
   - Évitez "No expiration" sauf si nécessaire
   - Renouvelez régulièrement vos tokens

5. **Si un token est compromis**
   - Révoquez-le immédiatement
   - Créez-en un nouveau

## Dépannage

### Erreur : "Authentication failed"

- Vérifiez que vous utilisez le token (pas le mot de passe)
- Vérifiez que le token n'a pas expiré
- Vérifiez que le token a la permission `repo`

### Erreur : "Permission denied"

- Vérifiez que vous avez les droits d'écriture sur le dépôt
- Vérifiez que vous utilisez le bon compte GitHub

### Le token ne fonctionne plus

- Vérifiez s'il a expiré
- Créez un nouveau token si nécessaire

## Alternative : GitHub CLI

Si vous préférez, vous pouvez utiliser GitHub CLI :

```powershell
# Installer GitHub CLI
# Puis :
gh auth login
```

Cela gère l'authentification automatiquement.

## Ressources

- Documentation GitHub sur les tokens : https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token
- Guide de sécurité GitHub : https://docs.github.com/en/authentication

