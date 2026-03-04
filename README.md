# MAKE QUALITY GREAT AGAIN

**code open source module d'autoqualification de GEMAUT**

Ce dépôt contient un script Python (`make_quality_great_again.py`) qui calcule un **masque de qualité** à partir d'une différence DSM/DTM (par exemple `MNS - MNT`).  

L'algorithme découpe l'image en tuiles, traite chaque tuile en parallèle, recolle les résultats, bouche les trous (valeurs NoData) puis applique un lissage final.

---

## Fonctionnalités principales

- **Découpage en tuiles** avec recouvrement pour traiter de grands rasters.
- **Traitement parallèle** (multiprocessing) sur plusieurs CPU.
- **Calcul de masque de qualité** par fenêtre glissante (percentile sur la distribution locale).
- **Assemblage des tuiles** avec gestion des zones de recouvrement.
- **Bouchage des trous (NoData)** avec plusieurs méthodes d’interpolation :  
- `griddata`  
- `idw` (IDW vectorisé et parallélisé)  
- `idw_old`  
- `window`  
- `linearnd`  
- `fast`  
- `hybrid` (méthode hybride équivalente à `xingng -FB:2:C:50,1:1:50:1`)
- **Moyenne glissante finale** pour lisser le masque.

---

## Dépendances- Python 3
- `numpy`
- `rasterio`
- `scipy`
- `tqdm`
- `osgeo` (GDAL)

Installation typique (exemple) : 
pip install numpy rasterio scipy tqdm
sudo apt install gdal-bin python3-gdal  # ou équivalent suivant l’OS

Utilisation
Exécution de base :
python3 make_quality_great_again.py \  
    -diff chemin/vers/diff_MNS_MNT.tif \  
    -out chemin/vers/masque_qualite.tif \  
    -RepTra chemin/vers/repertoire_temporaire \  
    -cpu 8 \  
    -interp hybrid

Arguments principaux
-diff : raster d’entrée (différence DSM/DTM).
-out : raster de sortie (masque de qualité final).
-RepTra : répertoire de travail temporaire (contient les tuiles et fichiers intermédiaires).
-cpu : nombre de processus parallèles.

Paramètres optionnels utiles
-no : valeur NoData (par défaut -9999).
-per : percentile utilisé pour le calcul local (défaut 0.05 = 5 %).
-demiwinl / -demiwinc : demi-taille de la fenêtre d’analyse (défaut 50).
-tile : taille des tuiles (défaut 500 pixels).
-pad : recouvrement entre tuiles (défaut 50 pixels).
-winavg : taille de la fenêtre de moyenne glissante finale (défaut 50).
-interp : méthode d’interpolation des NoData :
griddata, idw, idw_old, window, linearnd, fast, hybrid.
-clean : si présent, nettoie le répertoire temporaire au lancement.

Méthode d’interpolation hybride (-interp hybrid)

La méthode hybrid implémente en Python un équivalent de :
xingng -FB:2:C:50,1:1:50:1 -EM=-9999

Principe :
Identification des trous connexes (zones NoData).
Détection des pixels de bord du trou (connexité 4).
Pour chaque trou :
calcul d’une constante V_calc à partir des pixels de bord
(on retire 50 % des plus petites valeurs, puis on prend le minimum des restantes),
interpolation locale IDW sur les pixels de bord dans un rayon de 50 pixels (poids=1),
combinaison des deux selon la distance au bord :
près du bord → interpolation locale,
au centre du trou → constante V_calc.

Cette méthode remplit mieux les grands trous tout en restant raisonnablement rapide.

Organisation du traitement
Lecture des métadonnées de l’image (GetInfo).
Calcul du nombre de tuiles en X/Y (CalculNombreDallesXY).
Découpage en tuiles (MakeDecoupage).
Calcul du masque de qualité par tuile en parallèle (DoParallel).
Assemblage horizontal puis vertical des tuiles (Make_Assemblage_FINAL).
Interpolation des NoData (interpolate_nodata_*, selon -interp).
Moyenne glissante finale (apply_moving_average).

Licence et contributions
Code ouvert pour usage interne et expérimental autour de GEMAUT.
Contributions, issues et suggestions bienvenues via GitHub.
Merci de documenter clairement les ajouts de méthodes d’interpolation ou d’options en ligne de commande.


