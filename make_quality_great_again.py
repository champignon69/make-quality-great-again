#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, glob, shutil, datetime, time, math, re
#from lxml import etree
#from datetime import date
#from osgeo import ogr
#from shapely.geometry import Point, LineString, Polygon, MultiPolygon
#from shapefile import Shapefile
#from shapefile.const import SHPT_POINT, SHPT_POLYGON, SHPT_POLYLINE
import shutil
import argparse
from osgeo import ogr, gdal
#from geojson import Point, Polygon, MultiPolygon, Feature, FeatureCollection, dump
#from shapely.geometry import Polygon
# from rasterstats import zonal_stats  # Non utilisé actuellement
#
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.windows import from_bounds
import imageio
import scipy as sp
import scipy.misc as sm
import scipy.sparse.linalg as ssl
from scipy.ndimage import generic_filter, uniform_filter
from scipy.interpolate import LinearNDInterpolator

import time

from multiprocessing import Process,Pool,cpu_count
import concurrent.futures

from tqdm import tqdm

import signal

# Chemin vers xingng - peut être défini via variable d'environnement XINGNG_PATH
# ou utilise le chemin par défaut
# _DEFAULT_XINGNG_PATH = '/Volumes/ALI_Serveur/DEPLOIEMENT/bin_linux/xingng'

_DEFAULT_XINGNG_PATH = os.path.expanduser('~/xingng')

def get_xingng_path():
    """Récupère le chemin vers xingng et vérifie qu'il est accessible."""
    chem_xing = os.environ.get('XINGNG_PATH', _DEFAULT_XINGNG_PATH)
    chem_xing = os.path.expanduser(chem_xing)  # Assure l'expansion du ~ même via variable d'env

    # Vérification que xingng existe et est exécutable
    if not os.path.isfile(chem_xing):
        raise FileNotFoundError(
            f"xingng non trouvé à: {chem_xing}.\n"
            f"Définissez la variable d'environnement XINGNG_PATH ou modifiez le chemin par défaut dans le code."
        )
    if not os.access(chem_xing, os.X_OK):
        raise PermissionError(f"xingng n'est pas exécutable: {chem_xing}")
    
    return chem_xing

# Outils OTB (non utilisés actuellement dans le code, commentés pour référence)
# otbcli_ConnectedComponentSegmentation='/home/OTB/OTB-7.4.0-Linux64/bin/otbcli_ConnectedComponentSegmentation'
# otbcli_BinaryMorphologicalOperation='/home/OTB/OTB-7.4.0-Linux64/bin/otbcli_BinaryMorphologicalOperation'
 
############################################################################################################################		 
def GetValue(listInfo,chaine):
	
	for l in listInfo:
		lig = l.strip()
		if lig.find(chaine) != -1:
			return lig.split()[-1]
			
############################################################################################################################		 
############################################################################################################################		 
def GetInfo(chem_xing, cheminTIF):
	"""
	Récupère les métadonnées d'une image GeoTIFF.
	Version open source utilisant rasterio (chem_xing n'est plus utilisé mais conservé pour compatibilité).
	
	Returns:
		[PasX, PasY, Projection, X_0, X_1, Y_0, Y_1, phasage, NbreCol, NbreLig, GModel, GRaster]
	"""
	# Ouvrir l'image avec rasterio
	with rasterio.open(cheminTIF, 'r') as src:
		# Pas en X et Y (résolution)
		transform = src.transform
		PasX = abs(transform[0])  # pixel width
		PasY = abs(transform[4])  # pixel height (généralement négatif, on prend la valeur absolue)
		
		# Bounds (positions)
		bounds = src.bounds
		X_0 = bounds.left   # GAUCHE
		X_1 = bounds.right   # DROITE
		Y_0 = bounds.bottom  # BAS
		Y_1 = bounds.top     # HAUT
		
		# Nombre de colonnes et lignes
		NbreCol = float(src.width)
		NbreLig = float(src.height)
		
		# Code EPSG / Projection
		if src.crs is not None:
			Projection = int(src.crs.to_epsg()) if src.crs.to_epsg() is not None else -1
		else:
			Projection = -1
		
		# GTModelTypeGeoKey et GTRasterTypeGeoKey
		# Ces clés sont dans les tags GeoTIFF, généralement dans les tags de la bande
		GModel = -1
		GRaster = -1
		# Essayer d'abord les tags du dataset
		if hasattr(src, 'tags') and src.tags():
			tags = src.tags()
			if 'GTModelTypeGeoKey' in tags:
				try:
					GModel = int(tags['GTModelTypeGeoKey'])
				except (ValueError, TypeError):
					pass
			if 'GTRasterTypeGeoKey' in tags:
				try:
					GRaster = int(tags['GTRasterTypeGeoKey'])
				except (ValueError, TypeError):
					pass
		# Si pas trouvé, essayer les tags de la première bande
		if (GModel == -1 or GRaster == -1) and hasattr(src, 'tags'):
			try:
				band_tags = src.tags(1)
				if band_tags:
					if 'GTModelTypeGeoKey' in band_tags and GModel == -1:
						try:
							GModel = int(band_tags['GTModelTypeGeoKey'])
						except (ValueError, TypeError):
							pass
					if 'GTRasterTypeGeoKey' in band_tags and GRaster == -1:
						try:
							GRaster = int(band_tags['GTRasterTypeGeoKey'])
						except (ValueError, TypeError):
							pass
			except:
				pass
		
		# Phasage : déterminer si c'est "HG" (haut-gauche) ou "CP" (centre pixel)
		phasage = "HG"  # Par défaut
		if hasattr(src, 'tags') and src.tags():
			# Vérifier les tags pour le phasage
			tags = src.tags()
			if 'AREA_OR_POINT' in tags and tags['AREA_OR_POINT'] == 'Point':
				phasage = "CP"
	
	return [PasX, PasY, Projection, X_0, X_1, Y_0, Y_1, phasage, NbreCol, NbreLig, GModel, GRaster] 
	
#############################################################################################################################	
def read_as_2D_float(filename,no_data):
	
	with rasterio.open(filename, 'r') as dataset:
		data= dataset.read(1).astype(np.float32)
		data[data==no_data] = np.nan
		
	return data

#############################################################################################################################	
def get_pied_histo(histo, seuil_pied):
	for i, elt in enumerate(histo):
		if elt > seuil_pied:
			return i

#############################################################################################################################	
def get_haut_histo(histo, seuil_haut):
	for i, elt in enumerate(histo):
		if elt > seuil_haut:
			return i

#############################################################################################################################	
def calculate_cdf_percent(pixel_array, percentile):
	# Exclude the no-data value and calculate the 5% value of the CDF
	pixel_array = pixel_array[pixel_array != -9999]
	if len(pixel_array) == 0:
		return -9999  # Return no-data value if the window only contains no-data values
	sorted_pixels = np.sort(pixel_array)
	index_5_percent = int(np.ceil(percentile * len(sorted_pixels))) - 1
	return sorted_pixels[max(0, index_5_percent)]
	
#############################################################################################################################	
# def process_image(image, percentile):
# 	# Pad image to handle the borders
# 	padded_image = np.pad(image, 50, mode='constant', constant_values=-9999)
	
# 	# Use generic_filter from scipy.ndimage to apply the function over a 101x101 window
# 	result = generic_filter(padded_image, lambda x: calculate_cdf_percent(x, percentile), size=(101, 101), mode='constant', cval=-9999)
	
# 	# Crop the padded area off the result
# 	return result[50:-50, 50:-50]
	
#############################################################################################################################	
def process_image(image,dl,no_data,percentile):
	# Pad image to handle the borders
	padded_image = np.pad(image, dl, mode='constant', constant_values=no_data)
	# Use generic_filter from scipy.ndimage to apply the function over a 101x101 window
	result = generic_filter(padded_image, lambda x: calculate_cdf_percent(x, percentile), size=(2*dl+1, 2*dl+1), mode='constant', cval=no_data)
	# Crop the padded area off the result
	return result[dl:-dl, dl:-dl]
	
#############################################################################################################################	
def save_ABSOLUTE_image_with_same_geometry(image, output_filename, src_filename):
    # Calculer la valeur absolue de l'image
    abs_image = np.abs(image)
    
    # Ouvrez l'image source pour lire sa géométrie
    with rasterio.open(src_filename) as src:
        metadata = src.meta.copy()  # Copiez les métadonnées de l'image source
        
    # Mettez à jour les métadonnées avec les nouvelles dimensions si nécessaire
    metadata['height'], metadata['width'] = abs_image.shape
    metadata['dtype'] = abs_image.dtype  # Assurez-vous que le type de données correspond à l'image de sortie
    
    # Utilisez lechem_xings métadonnées copiées pour écrire l'image dans un fichier .tif
    with rasterio.open(output_filename, 'w', **metadata) as dst:
        dst.write(abs_image, 1)  # Écrit l'image dans la première bande en assumant qu'il s'agit d'une image à une seule bande
        
#############################################################################################################################	
def save_image_with_same_geometry(image, output_filename, src_filename):
	# Ouvrez l'image source pour lire sa géométrie
	with rasterio.open(src_filename) as src:
		metadata = src.meta.copy()  # Copiez les métadonnées de l'image source
		
	# Mettez à jour les métadonnées avec les nouvelles dimensions si nécessaire
	metadata['height'], metadata['width'] = image.shape
	metadata['dtype'] = image.dtype  # Assurez-vous que le type de données correspond à l'image de sortie
	
	# Utilisez les métadonnées copiées pour écrire l'image dans un fichier .tif
	with rasterio.open(output_filename, 'w', **metadata) as dst:
		dst.write(image, 1)  # Écrit l'image dans la première bande en assumant qu'il s'agit d'une image à une seule bande

#############################################################################################################################	
def interpolate_nodata_with_linearnd(chem_in, chem_out, no_data=-9999):
	"""
	Interpole les pixels nodata (valeur -9999) en utilisant LinearNDInterpolator de scipy.
	
	Args:
		chem_in: Chemin vers l'image d'entrée avec des pixels nodata
		chem_out: Chemin vers l'image de sortie avec les pixels nodata interpolés
		no_data: Valeur nodata (par défaut -9999)
	"""
	# Supprimer le fichier de sortie s'il existe déjà pour garantir l'écrasement
	if os.path.exists(chem_out):
		try:
			os.remove(chem_out)
		except OSError:
			pass  # Ignorer les erreurs si le fichier est verrouillé ou n'existe plus
	
	# Supprimer aussi les fichiers auxiliaires (.aux.xml) s'ils existent
	chem_out_aux = chem_out + '.aux.xml'
	if os.path.exists(chem_out_aux):
		try:
			os.remove(chem_out_aux)
		except OSError:
			pass
	
	# Lire l'image
	with rasterio.open(chem_in, 'r') as src:
		data = src.read(1).astype(np.float32)
		metadata = src.meta.copy()
	
	# Identifier les pixels valides (non nodata) et les pixels nodata
	mask_valid = (data != no_data) & ~np.isnan(data)
	mask_nodata = ~mask_valid
	
	# Si aucun pixel nodata, copier simplement l'image
	if not np.any(mask_nodata):
		print("Aucun pixel nodata à interpoler, copie de l'image originale.")
		with rasterio.open(chem_out, 'w', **metadata) as dst:
			dst.write(data, 1)
		return
	
	# Créer les coordonnées (ligne, colonne) pour tous les pixels
	rows, cols = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]), indexing='ij')
	
	# Coordonnées des pixels valides
	points_valid = np.column_stack([rows[mask_valid], cols[mask_valid]])
	values_valid = data[mask_valid]
	
	# Coordonnées des pixels nodata à interpoler
	points_nodata = np.column_stack([rows[mask_nodata], cols[mask_nodata]])
	
	# Créer l'interpolateur
	interpolator = LinearNDInterpolator(points_valid, values_valid)
	
	# Interpoler les valeurs pour les pixels nodata
	interpolated_values = interpolator(points_nodata)
	
	# Créer une copie de l'image et remplacer les pixels nodata par les valeurs interpolées
	result = data.copy()
	result[mask_nodata] = interpolated_values
	
	# Gérer les cas où l'interpolation peut retourner NaN (hors du domaine convexe)
	# Dans ce cas, on garde la valeur nodata originale
	mask_nan_interp = np.isnan(interpolated_values)
	if np.any(mask_nan_interp):
		print(f"Attention: {np.sum(mask_nan_interp)} pixels nodata n'ont pas pu être interpolés (hors du domaine convexe).")
		# Remplacer les NaN par la valeur nodata
		result[mask_nodata][mask_nan_interp] = no_data
	
	# Sauvegarder l'image résultante
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)
	
	print(f"Interpolation terminée: {np.sum(mask_nodata)} pixels nodata traités.")

#############################################################################################################################	
def apply_moving_average(chem_in, chem_out, window_size=50, no_data=-9999):
	"""
	Applique une moyenne sur une fenêtre glissante à l'image.
	
	Args:
		chem_in: Chemin vers l'image d'entrée
		chem_out: Chemin vers l'image de sortie
		window_size: Taille de la fenêtre glissante (par défaut 50x50)
		no_data: Valeur nodata (par défaut -9999)
	"""
	# Supprimer le fichier de sortie s'il existe déjà pour garantir l'écrasement
	if os.path.exists(chem_out):
		try:
			os.remove(chem_out)
		except OSError:
			pass
	
	# Supprimer aussi les fichiers auxiliaires (.aux.xml) s'ils existent
	chem_out_aux = chem_out + '.aux.xml'
	if os.path.exists(chem_out_aux):
		try:
			os.remove(chem_out_aux)
		except OSError:
			pass
	
	# Lire l'image
	with rasterio.open(chem_in, 'r') as src:
		data = src.read(1).astype(np.float32)
		metadata = src.meta.copy()
	
	# Créer un masque pour les pixels nodata
	mask_nodata = (data == no_data) | np.isnan(data)
	
	# Créer une copie des données en float64 pour les calculs
	data_float = data.astype(np.float64)
	
	# Remplacer les nodata par 0 pour le calcul de la somme
	data_sum = data_float.copy()
	data_sum[mask_nodata] = 0.0
	
	# Créer un masque de poids (1 pour valide, 0 pour nodata)
	weights = (~mask_nodata).astype(np.float64)
	
	# Calculer la somme et le nombre de pixels valides dans chaque fenêtre
	sum_window = uniform_filter(data_sum, size=window_size, mode='constant', cval=0.0)
	count_window = uniform_filter(weights, size=window_size, mode='constant', cval=0.0)
	
	# Calculer la moyenne uniquement là où il y a des pixels valides
	# Éviter la division par zéro
	result = np.where(count_window > 0, sum_window / count_window, no_data)
	
	# Conserver les nodata originaux (si un pixel était nodata, il reste nodata)
	result[mask_nodata] = no_data
	
	# Convertir en float32 pour la sauvegarde
	result = result.astype(np.float32)
	
	# Sauvegarder l'image résultante
	metadata['dtype'] = result.dtype
	with rasterio.open(chem_out, 'w', **metadata) as dst:
		dst.write(result, 1)
	
	print(f"Moyenne sur fenêtre glissante {window_size}x{window_size} terminée.")

#################################################################################################### 
### calcule le nombre de dalles en X et en Y en fonction des paramètres de chantier
def CalculNombreDallesXY_NEW(NbColonnes,NbLignes,dalleX,dalleY):
				
	### Calcul du nombre de dalles en X		
	if (NbColonnes%dalleX == 0):
		NbreDalleX=NbColonnes/dalleX
	else: 
		NbreDalleX=int(NbColonnes/dalleX)+1
		
	### Calcul du nombre de dalles en Y		
	if (NbLignes%dalleY == 0):
		NbreDalleY=NbLignes/dalleY
	else: 
		NbreDalleY=int(NbLignes/dalleY)+1
	
	return (NbreDalleX,NbreDalleY)  
	
#################################################################################################### 
### calcule le nombre de dalles en X et en Y en fonction des paramètres de chantier
def CalculNombreDallesXY(NbColonnes,NbLignes,Taille_dalle,Recouv_entre_dalles):
		
	### pour calculer le nombre de dalles en X et en Y		
	NumX=NbColonnes-Taille_dalle
	NumY=NbLignes-Taille_dalle
	Denom=Taille_dalle-Recouv_entre_dalles
		
	### Calcul du nombre de dalles en X		
	if (NumX%Denom == 0):
		NbreDalleX=int(NumX/Denom+1)
	else: 
		NbreDalleX=int((NumX/Denom)+1)+1
		
	### Calcul du nombre de dalles en Y		
	if (NumY%Denom == 0):
		NbreDalleY=int(NumY/Denom+1)
	else: 
		NbreDalleY=int((NumY/Denom)+1)+1
		
	return (NbreDalleX,NbreDalleY) 
	
#################################################################################################### 
def calculate_overlap_bounds(src1_path, src2_path):
	"""Calcule les bounds de recouvrement entre deux dalles (intersection)."""
	with rasterio.open(src1_path) as src1, rasterio.open(src2_path) as src2:
		bounds1 = src1.bounds
		bounds2 = src2.bounds
		
		# Intersection des bounds
		left = max(bounds1.left, bounds2.left)
		right = min(bounds1.right, bounds2.right)
		bottom = max(bounds1.bottom, bounds2.bottom)
		top = min(bounds1.top, bounds2.top)
		
		# Vérifier qu'il y a bien un recouvrement
		if left >= right or bottom >= top:
			return None
			
		return (left, bottom, right, top)  # X_0, Y_0, X_1, Y_1

#################################################################################################### 
def create_weight_image_horizontal(overlap_bounds, target_shape):
	"""Crée une image de poids basée sur C/NC (colonne/nombre de colonnes).
	Retourne un array numpy avec des valeurs entre 0 (gauche) et 1 (droite).
	
	Args:
		overlap_bounds: (left, bottom, right, top) - bounds de recouvrement
		target_shape: (height, width) - forme de l'array cible
	"""
	height, width = target_shape
	
	# Créer un array de poids basé sur la position en colonne
	# C/NC : numéro de colonne / nombre de colonnes
	col_indices = np.arange(width, dtype=np.float32)
	weight_array = col_indices / (width - 1) if width > 1 else np.ones(width)
	
	# Étendre sur toutes les lignes
	weights = np.tile(weight_array, (height, 1))
	
	return weights

#################################################################################################### 
def create_weight_image_vertical(overlap_bounds, target_shape):
	"""Crée une image de poids basée sur L/NL (ligne/nombre de lignes).
	Retourne un array numpy avec des valeurs entre 0 (haut) et 1 (bas).
	
	Args:
		overlap_bounds: (left, bottom, right, top) - bounds de recouvrement
		target_shape: (height, width) - forme de l'array cible
	"""
	height, width = target_shape
	
	# Créer un array de poids basé sur la position en ligne
	# L/NL : numéro de ligne / nombre de lignes
	row_indices = np.arange(height, dtype=np.float32)
	weight_array = row_indices / (height - 1) if height > 1 else np.ones(height)
	
	# Étendre sur toutes les colonnes
	weights = np.tile(weight_array.reshape(-1, 1), (1, width))
	
	return weights

#################################################################################################### 
def weighted_blend_overlap(dalle1_path, dalle2_path, weight1, weight2, overlap_bounds, output_path):
	"""Fait la moyenne pondérée de deux dalles dans la zone de recouvrement.
	Formule: I1*I2 + I3*I4 où I1=dalle1, I2=poids1, I3=dalle2, I4=poids2
	L'image de sortie a les bounds définis par overlap_bounds (équivalent à -cg:)."""
	with rasterio.open(dalle1_path) as src1, rasterio.open(dalle2_path) as src2:
		# Calculer les fenêtres de recouvrement pour chaque dalle
		window1 = src1.window(*overlap_bounds)
		window2 = src2.window(*overlap_bounds)
		
		# Lire les données dans la zone de recouvrement
		data1 = src1.read(1, window=window1).astype(np.float32)
		data2 = src2.read(1, window=window2).astype(np.float32)
		
		# Vérifier que les poids ont la bonne taille (devrait être le cas maintenant)
		if weight1.shape != data1.shape or weight2.shape != data1.shape:
			raise ValueError(f"Taille des poids incompatible: poids1={weight1.shape}, poids2={weight2.shape}, données={data1.shape}")
		
		# Gérer les no-data
		no_data1 = src1.nodata if src1.nodata is not None else -9999
		no_data2 = src2.nodata if src2.nodata is not None else -9999
		
		# Masques pour les valeurs valides
		valid1 = (data1 != no_data1) & ~np.isnan(data1)
		valid2 = (data2 != no_data2) & ~np.isnan(data2)
		
		# Calculer la moyenne pondérée
		result = np.full_like(data1, no_data1, dtype=np.float32)
		
		# Cas où les deux valeurs sont valides
		both_valid = valid1 & valid2
		result[both_valid] = data1[both_valid] * weight1[both_valid] + data2[both_valid] * weight2[both_valid]
		
		# Cas où seule la dalle1 est valide
		only1 = valid1 & ~valid2
		result[only1] = data1[only1]
		
		# Cas où seule la dalle2 est valide
		only2 = valid2 & ~valid1
		result[only2] = data2[only2]
		
		# Créer le transform pour l'image de sortie avec les bounds de recouvrement
		# overlap_bounds = (left, bottom, right, top)
		left, bottom, right, top = overlap_bounds
		height, width = result.shape
		
		# Calculer la résolution
		pixel_size_x = (right - left) / width
		pixel_size_y = (top - bottom) / height
		
		# Créer le transform (coin haut-gauche)
		from rasterio.transform import Affine
		transform = Affine(pixel_size_x, 0.0, left,
						   0.0, -pixel_size_y, top)
		
		# Métadonnées pour l'image de sortie
		metadata = src1.meta.copy()
		metadata.update({
			'height': height,
			'width': width,
			'transform': transform,
			'dtype': result.dtype,
			'nodata': no_data1,
			'compress': 'lzw'
		})
		
		# Écrire l'image de sortie
		with rasterio.open(output_path, 'w', **metadata) as dst:
			dst.write(result, 1)


def assemble_tiles_and_overlaps(masks, overlaps, output_path):

    # --- 1) MOSAÏQUE DES TUILES PRINCIPALES ---
    src_masks = [rasterio.open(p) for p in masks]
    nodata = src_masks[0].nodata if src_masks[0].nodata is not None else -9999

    mosaic, mosaic_transform = merge(
        src_masks,
        nodata=nodata,
        method="last"
    )

    meta = src_masks[0].meta.copy()
    meta.update({
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": mosaic_transform,
        "nodata": nodata,
        "compress": "lzw"
    })

    # Convertir la mosaïque en édition locale
    final = mosaic.astype(np.float32)

    # --- 2) APPLICATION DES PATCHS DE RECOUVREMENT ---
    for patch_path in overlaps:
        with rasterio.open(patch_path) as patch:

            # Fenêtre d'insertion (patch → mosaïque)
            win = from_bounds(
                patch.bounds.left,
                patch.bounds.bottom,
                patch.bounds.right,
                patch.bounds.top,
                transform=mosaic_transform
            )
            # Rasterio peut retourner des offsets/taille flottants : on les arrondit pour indexer numpy
            win = win.round_offsets().round_lengths()
            row_off, col_off = int(win.row_off), int(win.col_off)
            height, width = int(win.height), int(win.width)

            # Lecture du patch à la taille exacte de la fenêtre
            patch_arr = patch.read(1, out_shape=(height, width))

            # Extraction du bloc correspondant dans la mosaïque
            final_block = final[0, row_off:row_off+height,
                                   col_off:col_off+width]

            # Pixels valides = non nodata
            valid = patch_arr != nodata

            # Remplacement
            final_block[valid] = patch_arr[valid]

            # Réécriture dans la mosaïque
            final[0, row_off:row_off+height,
                     col_off:col_off+width] = final_block

    # --- 3) ENREGISTREMENT ---
    with rasterio.open(output_path, "w", **meta) as dst:
        dst.write(final)

    # Fermeture des tuiles principales
    for s in src_masks:
        s.close()


def assemble_lines_and_overlaps(lines, overlaps, output_path):
    """
    Variante verticale : assemble les lignes mosaïquées, puis applique les patches verticaux.
    Comportement nodata : les pixels nodata des patches n'écrasent pas les pixels valides existants.
    """
    src_lines = [rasterio.open(p) for p in lines]
    nodata = src_lines[0].nodata if src_lines[0].nodata is not None else -9999

    mosaic, mosaic_transform = merge(
        src_lines,
        nodata=nodata,
        method="last"
    )

    meta = src_lines[0].meta.copy()
    meta.update({
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": mosaic_transform,
        "nodata": nodata,
        "compress": "lzw"
    })

    final = mosaic.astype(np.float32)

    for patch_path in overlaps:
        with rasterio.open(patch_path) as patch:
            win = from_bounds(
                patch.bounds.left,
                patch.bounds.bottom,
                patch.bounds.right,
                patch.bounds.top,
                transform=mosaic_transform
            )
            win = win.round_offsets().round_lengths()
            row_off, col_off = int(win.row_off), int(win.col_off)
            height, width = int(win.height), int(win.width)

            patch_arr = patch.read(1, out_shape=(height, width))
            final_block = final[0, row_off:row_off+height,
                                   col_off:col_off+width]
            valid = patch_arr != nodata
            final_block[valid] = patch_arr[valid]
            final[0, row_off:row_off+height,
                     col_off:col_off+width] = final_block

    with rasterio.open(output_path, "w", **meta) as dst:
        dst.write(final)

    for s in src_lines:
        s.close()


#################################################################################################### 
def assemble_horizontal_OLD(image_paths, output_path):
	"""Assemble des images de gauche à droite. Les pixels suivants écrasent les précédents en cas de recouvrement.
	Les pixels NoData ne remplacent jamais les pixels valides grâce à skip_empty=True."""
	from rasterio.merge import merge
	
	# Trier les images par leur position X (left) pour garantir l'ordre de gauche à droite
	image_bounds = []
	for path in image_paths:
		with rasterio.open(path) as src:
			image_bounds.append((src.bounds.left, path))
	
	# Trier par position X (left)
	image_bounds.sort(key=lambda x: x[0])
	sorted_paths = [path for _, path in image_bounds]
	
	# Récupérer la valeur nodata de la première image (on suppose qu'elles sont identiques)
	with rasterio.open(sorted_paths[0]) as first:
		nodata_val = first.nodata if first.nodata is not None else -9999
	
	# Ouvrir toutes les images dans l'ordre trié
	srcs = [rasterio.open(path) for path in sorted_paths]
	
	try:
		# Utiliser merge pour assembler (method='last' pour que les pixels suivants écrasent les précédents)
		# skip_empty=True empêche les pixels NoData d'écraser les pixels valides
		# L'ordre dans la liste détermine la priorité : les dernières images écrasent les premières
		mosaic, out_trans = merge(
			srcs,
			method='last',
			nodata=nodata_val,
			skip_empty=True  # Ne pas remplacer les pixels valides par du NoData
		)
		
		# Métadonnées de sortie
		out_meta = srcs[0].meta.copy()
		out_meta.update({
			'driver': 'GTiff',
			'height': mosaic.shape[1],
			'width': mosaic.shape[2],
			'transform': out_trans,
			'compress': 'lzw',
			'nodata': nodata_val
		})
		
		# Écrire l'image assemblée
		with rasterio.open(output_path, 'w', **out_meta) as dst:
			dst.write(mosaic)
	finally:
		# Fermer toutes les sources
		for src in srcs:
			src.close()

#################################################################################################### 	
def Make_Assemblage_FINAL(chem_out, chem_xing, NbreDalleX, NbreDalleY, RepTra):
	"""
	Assemble les dalles en utilisant rasterio au lieu de xing.
	Note: chem_xing n'est plus utilisé mais conservé pour compatibilité de signature.
	"""
	####################################################################################################################################################
	## on raboute tout d'abord 2 dalles côte à côte (en X)	   #########################################################################################
	####################################################################################################################################################
	
	# Boucle avec barre de progression pour le raboutage des dalles côte à côte
	for y in tqdm(range(NbreDalleY), desc="Raboutage en colonnes"):
		for x in tqdm(range(NbreDalleX-1), desc="Progression en Colonne", leave=False):
			
			## Nom de la dalle courante	
			chem_MASK_QUALITY_dalle_xy = os.path.join(RepTra, "Dalle_%s_%s" % (x, y), "MASK_%s_%s.tif" % (x, y))
			chem_MASK_QUALITY_dalle_xy_droite = os.path.join(RepTra, "Dalle_%s_%s" % (x+1, y), "MASK_%s_%s.tif" % (x+1, y))
			
			# Calculer les bounds de recouvrement
			overlap_bounds = calculate_overlap_bounds(chem_MASK_QUALITY_dalle_xy, chem_MASK_QUALITY_dalle_xy_droite)
			
			if overlap_bounds is None:
				print(f"Attention: Pas de recouvrement entre dalle ({x},{y}) et ({x+1},{y})")
				continue
			
			# Lire la shape des données pour créer les poids de la bonne taille
			with rasterio.open(chem_MASK_QUALITY_dalle_xy_droite) as src:
				window = src.window(*overlap_bounds)
				target_shape = (int(window.height), int(window.width))
			
			# Créer les images de poids (horizontal: C/NC)
			weight_droite = create_weight_image_horizontal(overlap_bounds, target_shape)
			weight_gauche = 1.0 - weight_droite
			
			# Faire la moyenne pondérée dans la zone de recouvrement
			chem_reconstruction = os.path.join(RepTra, 'reconstruction_dalle_%s_%s_%s.tif' % (x, x+1, y))
			weighted_blend_overlap(
				chem_MASK_QUALITY_dalle_xy,
				chem_MASK_QUALITY_dalle_xy_droite,
				weight_gauche,
				weight_droite,
				overlap_bounds,
				chem_reconstruction
			)
			
	####################################################################################################################################################		
	## on raboute toutes les dalles sur une même rangée		#########################################################################################
	####################################################################################################################################################
	print("Raboutage en ligne NbreDalleY = ", NbreDalleY)			

	## on réassemble tout	
	for y in tqdm(range(NbreDalleY), desc="Raboutage en ligne"):
		
		# Construire la liste des images à assembler: dalles originales + images de transition
		image_paths = []
		overlaps = []
		# Ajouter les dalles originales
		for x in range(NbreDalleX):
			image_paths.append(os.path.join(RepTra, "Dalle_%s_%s" % (x, y), "MASK_%s_%s.tif" % (x, y)))
		
		# Ajouter les images de transition
		for x in range(NbreDalleX-1):
			overlaps.append(os.path.join(RepTra, 'reconstruction_dalle_%s_%s_%s.tif' % (x, x+1, y)))
		
		# Assembler de gauche à droite
		chem_final_tmp = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)

		print("image_paths = ", image_paths)
		print("chem_final_tmp = ", chem_final_tmp)
		print("overlaps = ", overlaps)

		try:
			assemble_tiles_and_overlaps(image_paths, overlaps, chem_final_tmp)
		except Exception as e:
			print("ERREUR dans assemble_horizontal:", type(e), e)
			import traceback
			traceback.print_exc()
			print("Vérifiez que la fonction assemble_horizontal est bien importée/définie, qu'il n'y a pas d'erreur de nom de variable ou d'accès aux fichiers ci-dessus.")
			print("Voici la liste des images à assembler, pour vérification des accès fichiers :")
			for ip in image_paths:
				print("  - ", ip, "-->", os.path.exists(ip))
			print("chem_final_tmp =", chem_final_tmp, "--> dossier existe ?", os.path.exists(os.path.dirname(chem_final_tmp)))
			raise  # relancer l'exception pour arrêt si debug

		print("FIN")
	####################################################################################################################################################
	## on assemble les rangées entre elles		 #####################################################################################################
	####################################################################################################################################################
	print('===============================================================')
	## on réassemble tout	
	for y in tqdm(range(NbreDalleY-1), desc="Assemblage final"):
			
		print("Hello")	
		chem_dalle_y = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)
		chem_dalle_y_dessous = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % (y+1))
		
		# Calculer les bounds de recouvrement
		overlap_bounds = calculate_overlap_bounds(chem_dalle_y, chem_dalle_y_dessous)
		
		if overlap_bounds is None:
			print(f"Attention: Pas de recouvrement entre ligne {y} et {y+1}")
			continue
		
		# Lire la shape des données pour créer les poids de la bonne taille
		with rasterio.open(chem_dalle_y_dessous) as src:
			window = src.window(*overlap_bounds)
			target_shape = (int(window.height), int(window.width))
		
		# Créer les images de poids (vertical: L/NL)
		weight_bas = create_weight_image_vertical(overlap_bounds, target_shape)
		weight_haut = 1.0 - weight_bas
		
		# Faire la moyenne pondérée dans la zone de recouvrement
		chem_reconstruction = os.path.join(RepTra, 'reconstruction_dalle_%s_%s.tif' % (y, y+1))
		weighted_blend_overlap(
			chem_dalle_y,
			chem_dalle_y_dessous,
			weight_haut,
			weight_bas,
			overlap_bounds,
			chem_reconstruction
		)
			
	####################################################################################################################################################
	## Assemblage final de toutes les lignes		 #####################################################################################################
	####################################################################################################################################################
	
	# Assemblage final vertical : lignes mosaïquées + patches verticaux
	image_paths = []
	overlaps = []
	
	# Ajouter les lignes assemblées
	for y in range(NbreDalleY):
		chem_dalle_y = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)
		image_paths.append(chem_dalle_y)
		
	# Ajouter les images de transition verticales
	for y in range(NbreDalleY-1):
		overlaps.append(os.path.join(RepTra, 'reconstruction_dalle_%s_%s.tif' % (y, y+1)))
	
	# Assemblage final (vertical)
	assemble_lines_and_overlaps(image_paths, overlaps, chem_out)

	
#################################################################################################### 
def crop_tile(args):
	"""
	Fonction pour découper une dalle d'une image (utilisée en parallèle).
	
	Args:
		args: tuple (chem_in, Chem_decoup, ligmin, ligmax, colmin, colmax)
	"""
	chem_in, Chem_decoup, ligmin, ligmax, colmin, colmax = args
	
	# Ouvrir l'image source
	with rasterio.open(chem_in, 'r') as src:
		# Calculer les dimensions de la fenêtre
		# ligmin inclus, ligmax exclusif -> height = ligmax - ligmin
		# colmin inclus, colmax exclusif -> width = colmax - colmin
		height = ligmax - ligmin
		width = colmax - colmin
		
		# Tronquer si la fenêtre dépasse les limites de l'image
		max_row = src.height
		max_col = src.width
		
		# Si la fenêtre est complètement en dehors de l'image, on ne fait rien
		if ligmin >= max_row or colmin >= max_col or ligmax <= 0 or colmax <= 0:
			print(f"Attention: Fenêtre complètement hors limites pour {Chem_decoup}")
			return
		
		# Ajuster les offsets pour qu'ils soient >= 0
		row_off = max(0, ligmin)
		col_off = max(0, colmin)
		
		# Ajuster la taille si la fenêtre dépasse les limites
		if row_off + height > max_row:
			height = max_row - row_off
		if col_off + width > max_col:
			width = max_col - col_off
		
		# Vérifier que la fenêtre est valide après ajustement
		if height <= 0 or width <= 0:
			print(f"Attention: Fenêtre invalide pour {Chem_decoup} (height={height}, width={width})")
			return
		
		# Créer la fenêtre de découpage
		window = rasterio.windows.Window(col_off=col_off, row_off=row_off, width=width, height=height)
		
		# Lire les données dans la fenêtre
		data = src.read(1, window=window)
		
		# Calculer le nouveau transform pour cette fenêtre
		transform = rasterio.windows.transform(window, src.transform)
		
		# Copier les métadonnées et les mettre à jour
		metadata = src.meta.copy()
		metadata.update({
			'height': height,
			'width': width,
			'transform': transform,
			'compress': 'lzw'
		})
		
		# Écrire la dalle découpée
		with rasterio.open(Chem_decoup, 'w', **metadata) as dst:
			dst.write(data, 1)

#################################################################################################### 
def MakeDecoupage(chem_in, RepTra, NbreDalleX, NbreDalleY, iTailleparcelle, iTailleRecouvrement, iNbreCPU):
	"""
	Découpe l'image en dalles en utilisant rasterio (version open source).
	"""
	tasks = []
	
	for x in range(NbreDalleX):
		for y in range(NbreDalleY):
			#créer le nom du répertoire
			RepDalleXY=os.path.join(RepTra,"Dalle_%s_%s"%(x,y))
			#créer le répertoire, s'il n'existe pas déjà
			if not os.path.isdir(RepDalleXY): os.mkdir(RepDalleXY)
								
			#Détermination de col_min, col_max, lig_min, lig_max pour la dalle XY
			colminDalleXY=x*(iTailleparcelle-iTailleRecouvrement)
			colmaxDalleXY=x*(iTailleparcelle-iTailleRecouvrement)+iTailleparcelle
			ligminDalleXY=y*(iTailleparcelle-iTailleRecouvrement)
			ligmaxDalleXY=y*(iTailleparcelle-iTailleRecouvrement)+iTailleparcelle
				
			#fichier out mns
			Chem_decoup=os.path.join(RepDalleXY,"IN_%s_%s.tif"%(x,y))
			
			# Préparer les arguments pour la fonction de découpage
			args = (chem_in, Chem_decoup, ligminDalleXY, ligmaxDalleXY, colminDalleXY, colmaxDalleXY)
			tasks.append(args)
			
	# Initialize the pool
	with Pool(processes=iNbreCPU, initializer=init_worker) as pool:
		results = list(tqdm(pool.imap_unordered(crop_tile, tasks), total=len(tasks), desc="Découpage en parallèle des dalles"))
		
	return
			
#################################################################################################### 
def diff_2_mask_quality(args):
	chem_in, chem_out, dl, no_data, percentile = args
	#print("chem_in    >>> ",chem_in)
	#print("chem_out   >>> ",chem_out)
	#print("dl         >>> ",dl)
	#print("no_data    >>> ",no_data)
	#print("percentile >>> ",percentile)
	
	data_in = read_as_2D_float(chem_in, no_data)
	result = process_image(data_in, dl, no_data, percentile)
	save_ABSOLUTE_image_with_same_geometry(result, chem_out, chem_in)
	return

#################################################################################################### 
# def diff_2_mask_quality_BIS(chem_in, chem_out, dl, no_data, percentile):

# 	data_in = read_as_2D_float(chem_in,no_data)
# 	result = process_image(data_in,dl,no_data)
# 	save_ABSOLUTE_image_with_same_geometry(result, chem_out, chem_in)
# 	return
	
#################################################################################################### 
def init_worker():
	signal.signal(signal.SIGINT, signal.SIG_IGN)
	
#################################################################################################### 
def DoParallel(RepTra, NbreDalleX, NbreDalleY, dl, no_data, percentile, iNbreCPU):
		
	tasks = []
	
	## lancement sur chaque dalle
	for x in range(NbreDalleX):
		for y in range(NbreDalleY):
			#créer le nom du répertoire
			RepDalleXY=os.path.join(RepTra,"Dalle_%s_%s"%(x,y))
			#fichier out mns
			chem_in_dalle=os.path.join(RepDalleXY,"IN_%s_%s.tif"%(x,y))
			chem_in_dalle_NEG=os.path.join(RepDalleXY,"IN_%s_%s_NEG.tif"%(x,y))
			chem_out_dalle=os.path.join(RepDalleXY,"MASK_%s_%s.tif"%(x,y))
			#
			cmd_xing="%s -i %s -e'I1.1>0?%s:I1.1' -o %s -n:" %(chem_xing,chem_in_dalle,no_data,chem_in_dalle_NEG)
			print('cmd_xing > ',cmd_xing)
			os.system(cmd_xing)
			#
			args=(chem_in_dalle_NEG, chem_out_dalle, dl, no_data, percentile)
			tasks.append(args)

	# Initialize the pool
	with Pool(processes=iNbreCPU, initializer=init_worker) as pool:
		# Use tqdm to show the progress bar
		results = list(tqdm(pool.imap_unordered(diff_2_mask_quality, tasks), total=len(tasks),desc="Calcul des masques de qualité par dalle en //"))
				

#################################################################################################### 
def DoParallel_OLD(RepTra, NbreDalleX, NbreDalleY, dl, no_data, percentile, iNbreCPU):
		
	tasks = []
	
	## lancement sur chaque dalle
	for x in range(NbreDalleX):
		for y in range(NbreDalleY):
			#créer le nom du répertoire
			RepDalleXY=os.path.join(RepTra,"Dalle_%s_%s"%(x,y))
			#fichier out mns
			chem_in_dalle=os.path.join(RepDalleXY,"IN_%s_%s.tif"%(x,y))
			chem_out_dalle=os.path.join(RepDalleXY,"MASK_%s_%s.tif"%(x,y))
			#
			args=(chem_in_dalle, chem_out_dalle, dl, no_data, percentile)
			tasks.append(args)
		
	# Initialize the pool
	with Pool(processes=iNbreCPU, initializer=init_worker) as pool:
		# Use tqdm to show the progress bar
		results = list(tqdm(pool.imap_unordered(diff_2_mask_quality, tasks), total=len(tasks),desc="Calcul des masques de qualité par dalle en //"))
				
#############################################################################################################################	
#############################################################################################################################	
#############################################################################################################################					
if __name__ == '__main__':

	try:
		#
		start_time = time.time()
		#
		parser = argparse.ArgumentParser(description='MAKE QUALITY GREAT AGAIN ALL ZONE - Version Spatialisée & Parallélisée')
		parser.add_argument("-diff", type=str, help="Différence DSM/DTM en entrée")
		parser.add_argument("-out", type=str, help="Masque de Qualité en sortie")
		parser.add_argument("-no", type=int, default=-9999, help="Valeur de No Data")
		#parser.add_argument("-bin", type=float, default=0.1, help="Valeur de bin")
		parser.add_argument("-per", type=float, default=0.05, help="Valeur de percentile")
		#parser.add_argument("-count", type=int, default=100, help="Nombre de points minimal pour que la statistique soit valide")
		parser.add_argument("-demiwinl", type=int, default=50, help="Demie-taille en ligne de la fenêtre d'analyse")
		parser.add_argument("-demiwinc", type=int, default=50, help="Demie-taille en colonne de la fenêtre d'analyse")
		#parser.add_argument("-offl", type=int, default=25, help="Décalage en ligne entre 2 analyses / Facteur de sous-échantillonnage en ligne")
		#parser.add_argument("-offc", type=int, default=25, help="Décalage en colonne entre 2 analyses / Facteur de sous-échantillonnage en colonne")
		parser.add_argument("-tile", type=int, default=500, help="Tile / Taille de la tuile")
		parser.add_argument("-pad", type=int, default=50, help="Pad / Recouvrement entre tuiles")	
		parser.add_argument("-RepTra", type=str, help="Répertoire de Travail")
		parser.add_argument("-cpu", type=int, help="Nombre de CPU diponibles")
		parser.add_argument("-winavg", type=int, default=50, help="Taille de la fenêtre glissante pour la moyenne (par défaut 50x50)")
		
		args = parser.parse_args(sys.argv[1:])
		#
		# Initialisation du chemin vers xingng (vérifie l'existence et les permissions)
		chem_xing = get_xingng_path()
		print(f"Utilisation de xingng: {chem_xing}")
		#
		chem_in=args.diff
		chem_out=args.out
		no_data=args.no
		#bin_step=args.bin
		percentile=args.per
		#count_valid=args.count
		demiwinsize_lig=args.demiwinl
		demiwinsize_col=args.demiwinc
		#offset_lig=args.offl
		#offset_col=args.offc
		RepTra=args.RepTra
		iNbreCPU=args.cpu
		#print(iNbreCPU)
		iTailleparcelle=args.tile
		iTailleRecouvrement=args.pad
		winavg_size=args.winavg
		dl=demiwinsize_lig
		dc=demiwinsize_col
		
		### Découpage
		infos=GetInfo(chem_xing, chem_in)
		
		NbreCol=infos[8]
		NbreLig=infos[9]
				
		NombreDallesXY=CalculNombreDallesXY(NbreCol,NbreLig,iTailleparcelle,iTailleRecouvrement)
		NbreDalleX=NombreDallesXY[0]
		NbreDalleY=NombreDallesXY[1]
		
		#Decoupage en parallèle
		MakeDecoupage(chem_in, RepTra, NbreDalleX, NbreDalleY, iTailleparcelle, iTailleRecouvrement, iNbreCPU)
		
		#Traitement/Calcul en parallèle
		DoParallel(RepTra, NbreDalleX, NbreDalleY, dl, no_data, percentile, iNbreCPU)
		
		#### Assemblage final - avec xingng = A REMPLACER !
		# Le fichier d'assemblage est temporaire (fini par _tmp.tif)
		chem_out_tmp = chem_out.replace('.tif', '_tmp.tif')
		Make_Assemblage_FINAL(chem_out_tmp, chem_xing, NbreDalleX, NbreDalleY, RepTra)
		
		#### Post-traitement 2: Bouchage des zones avec no data (interpolation avec LinearNDInterpolator)
		chem_out_clean_bouchage = chem_out.replace('.tif', '_clean_bouchage.tif')
		print("Post-traitement 2: Interpolation des pixels nodata avec LinearNDInterpolator...")
		interpolate_nodata_with_linearnd(chem_out_tmp, chem_out_clean_bouchage, no_data)


		
		#### Post-traitement 3: Moyenne sur fenêtre glissante (open source)
		# Le fichier final utilise le nom spécifié dans --out
		print(f"Post-traitement 3: Moyenne fenêtre glissante {winavg_size}x{winavg_size}...")
		apply_moving_average(chem_out_clean_bouchage, chem_out, window_size=winavg_size, no_data=no_data)
		
		print("Fichier final généré: %s" % chem_out)
		
		print(time.time() - start_time, "seconds")
		
		print('FIN')
				  
	except (RuntimeError, TypeError, NameError):
		print ("ERREUR: ", NameError)


