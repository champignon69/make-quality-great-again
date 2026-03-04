# MAKE QUALITY GREAT AGAIN

**Module open source d'autoqualification de GEMAUT**

Ce dépôt contient un script Python (`make_quality_great_again.py`) qui
calcule un **masque de qualité** à partir d'une **différence DSM/DTM**
(par exemple **MNS - MNT**).

L'algorithme :

1.  découpe l'image en **tuiles**
2.  traite chaque tuile **en parallèle**
3.  **recoud les résultats**
4.  **interpole les trous (NoData)**
5.  applique un **lissage final**

------------------------------------------------------------------------

# Fonctionnalités principales

-   Découpage en **tuiles avec recouvrement** pour traiter de grands
    rasters
-   **Traitement parallèle (multiprocessing)** sur plusieurs CPU
-   Calcul d'un **masque de qualité par fenêtre glissante** (percentile
    sur la distribution locale)
-   **Assemblage automatique** des tuiles avec gestion du recouvrement
-   **Interpolation des trous (NoData)** avec plusieurs méthodes :
    -   `griddata`
    -   `idw` (IDW vectorisé et parallélisé)
    -   `idw_old`
    -   `window`
    -   `linearnd`
    -   `fast`
    -   `hybrid` (équivalent à `xingng -FB:2:C:50,1:1:50:1`)
-   **Moyenne glissante finale** pour lisser le masque

------------------------------------------------------------------------

# Dépendances

-   Python 3
-   numpy
-   rasterio
-   scipy
-   tqdm
-   osgeo (GDAL)

### Installation typique

``` bash
pip install numpy rasterio scipy tqdm
sudo apt install gdal-bin python3-gdal
```

(adapter selon votre OS)

------------------------------------------------------------------------

# Utilisation

### Exécution de base

``` bash
python3 make_quality_great_again.py \
-diff chemin/vers/diff_MNS_MNT.tif \
-out chemin/vers/masque_qualite.tif \
-RepTra chemin/vers/repertoire_temporaire \
-cpu 8 \
-interp hybrid
```

------------------------------------------------------------------------

# Arguments principaux

  Argument    Description
  ----------- --------------------------------------------
  `-diff`     Raster d'entrée (différence DSM/DTM)
  `-out`      Raster de sortie (masque de qualité final)
  `-RepTra`   Répertoire de travail temporaire
  `-cpu`      Nombre de processus parallèles

Le répertoire temporaire contient :

-   les **tuiles**
-   les **fichiers intermédiaires**

------------------------------------------------------------------------

# Paramètres optionnels

  Argument      Description                                     Défaut
  ------------- ----------------------------------------------- -----------------
  `-no`         Valeur NoData                                   -9999
  `-per`        Percentile utilisé pour le calcul local         0.05
  `-demiwinl`   Demi-taille fenêtre ligne                       50
  `-demiwinc`   Demi-taille fenêtre colonne                     50
  `-tile`       Taille des tuiles                               500 px
  `-pad`        Recouvrement entre tuiles                       50 px
  `-winavg`     Fenêtre moyenne glissante finale                50
  `-interp`     Méthode interpolation NoData                    voir ci-dessous
  `-clean`      Nettoie le répertoire temporaire au lancement   option

------------------------------------------------------------------------

# Méthodes d'interpolation disponibles

    griddata
    idw
    idw_old
    window
    linearnd
    fast
    hybrid

------------------------------------------------------------------------

# Méthode hybride (`-interp hybrid`)

Implémente en Python un équivalent de :

    xingng -FB:2:C:50,1:1:50:1 -EM=-9999

### Principe

1.  Identification des **trous connexes (zones NoData)**
2.  Détection des **pixels de bord du trou** (connexité 4)
3.  Pour chaque trou :

#### Étape 1 --- Calcul d'une constante `V_calc`

-   récupération des pixels de bord
-   suppression des **50 % plus petites valeurs**
-   prise du **minimum des restantes**

#### Étape 2 --- Interpolation locale

Interpolation **IDW sur les pixels de bord**

-   rayon : **50 pixels**
-   poids : **1**

#### Étape 3 --- Combinaison des deux

Selon la **distance au bord du trou** :

  Position         Valeur utilisée
  ---------------- ----------------------
  proche du bord   interpolation locale
  centre du trou   constante `V_calc`

Cette méthode permet :

-   un **remplissage robuste des grands trous**
-   tout en restant **relativement rapide**

------------------------------------------------------------------------

# Organisation du traitement

1.  Lecture des métadonnées de l'image\
    `GetInfo`

2.  Calcul du nombre de tuiles\
    `CalculNombreDallesXY`

3.  Découpage en tuiles\
    `MakeDecoupage`

4.  Calcul du masque par tuile (parallèle)\
    `DoParallel`

5.  Assemblage des tuiles\
    `Make_Assemblage_FINAL`

6.  Interpolation des NoData\
    `interpolate_nodata_*`

(selon `-interp`)

7.  Lissage final\
    `apply_moving_average`

------------------------------------------------------------------------

# Licence et contributions

Code ouvert pour **usage interne et expérimental autour de GEMAUT**.

Les contributions sont bienvenues :

-   nouvelles méthodes d'interpolation
-   optimisation CPU / mémoire
-   amélioration de la CLI
-   documentation

Merci de **documenter clairement toute nouvelle option ou méthode**.
