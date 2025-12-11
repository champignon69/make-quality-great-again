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
import imageio
import scipy as sp
import scipy.misc as sm
import scipy.sparse.linalg as ssl
from scipy.ndimage import generic_filter

import time

from multiprocessing import Process,Pool,cpu_count
import concurrent.futures

from tqdm import tqdm

import signal

# Chemin vers xingng - peut être défini via variable d'environnement XINGNG_PATH
# ou utilise le chemin par défaut
_DEFAULT_XINGNG_PATH = '/Volumes/ALI_Serveur/DEPLOIEMENT/bin_linux/xingng'

def get_xingng_path():
    """Récupère le chemin vers xingng et vérifie qu'il est accessible."""
    chem_xing = os.environ.get('XINGNG_PATH', _DEFAULT_XINGNG_PATH)
    
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
	
	### récupérer les infos dans un flux
	cmdinfo="%s -i %s  -n:stdout" % (chem_xing,cheminTIF)
	pipe= os.popen(cmdinfo)
	listInfo = pipe.readlines()
	pipe.close()
	
	### initialisation des variables
	PasX, PasY, projection, X_0, X_1, Y_0, Y_1 = 0, 0, 0, 0, 0, 0, 0 
	phasage="HG"
	
	### récupération des données
	PasX=float(GetValue(listInfo,"pas en X"))
	PasY=float(GetValue(listInfo,"pas en Y"))
	X_0=float(GetValue(listInfo,"position en X (GAUCHE)"))
	Y_1=float(GetValue(listInfo,"position en Y (HAUT)"))
	X_1=float(GetValue(listInfo,"position en X (DROITE)"))
	Y_0=float(GetValue(listInfo,"position en Y (BAS)")) 
	NbreCol=float(GetValue(listInfo,"nombre de colonnes"))
	NbreLig=float(GetValue(listInfo,"nombre de lignes"))
	if GetValue(listInfo,"GTModelTypeGeoKey"): GModel=int(GetValue(listInfo,"GTModelTypeGeoKey"))	  
	else: GModel=-1
	if GetValue(listInfo,"GTRasterTypeGeoKey"): GRaster=int(GetValue(listInfo,"GTRasterTypeGeoKey"))
	else: GRaster=-1
	if GetValue(listInfo,"code EPSG"): Projection=int(GetValue(listInfo,"code EPSG"))
	else: Projection=-1
	### traitement différent pour le phasage
	#print('GetValue(listInfo,"CENTRE") >> ',GetValue(listInfo,"CENTRE"))
	if  (not GetValue(listInfo,"CENTRE PIXEL")):
		phasage="HG"
	else: phasage="CP"
	#print('phasage >> ',phasage)
		
	return [PasX,PasY,Projection,X_0,X_1,Y_0,Y_1,phasage,NbreCol,NbreLig,GModel,GRaster] 
	
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
def calculate_cdf_percent(pixel_array, percentile, no_data):
	# Exclude the no-data value and NaN, then calculate the percentile value of the CDF
	# Filtrer les valeurs invalides (no-data et NaN)
	# ~np.isnan(pixel_array) permet d'exclure des pixels qui auraient la valeur NaN ('Not a Number')
	# Même si on code le no_data en entrée et en sortie avec la valeur de l'utilisateur,
	# certaines opérations numériques (division par zéro, etc.) peuvent générer des NaN à l'intérieur du traitement/intermédiaires.
	# On continue donc à filtrer les NaN pour ne jamais les considérer comme valeurs valides dans les fenêtres de calcul de percentile.
	valid_mask = (pixel_array != no_data) & ~np.isnan(pixel_array)
	pixel_array = pixel_array[valid_mask]
	if len(pixel_array) == 0:
		return no_data  # Return no-data value if the window only contains no-data values
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
	# Convertir les NaN en no_data pour uniformiser
	image = np.where(np.isnan(image), no_data, image)
	
	# Pad image to handle the borders
	padded_image = np.pad(image, dl, mode='constant', constant_values=no_data)
	
	# Use generic_filter from scipy.ndimage to apply the function over a (2*dl+1)x(2*dl+1) window
	# calculate_cdf_percent filtre automatiquement les no_data
	result = generic_filter(padded_image, lambda x: calculate_cdf_percent(x, percentile, no_data), size=(2*dl+1, 2*dl+1), mode='constant', cval=no_data)
	
	# Crop the padded area off the result
	result_cropped = result[dl:-dl, dl:-dl]
	
	# S'assurer que les valeurs invalides sont bien à no_data (déjà fait par calculate_cdf_percent, mais on double-vérifie)
	result_cropped = np.where(np.isnan(result_cropped), no_data, result_cropped)
	
	return result_cropped
	
#############################################################################################################################	
def save_ABSOLUTE_image_with_same_geometry(image, output_filename, src_filename, no_data):
    # Calculer la valeur absolue de l'image
    abs_image = np.abs(image)
    
    # Ouvrez l'image source pour lire sa géométrie
    with rasterio.open(src_filename) as src:
        metadata = src.meta.copy()  # Copiez les métadonnées de l'image source
        
    # Convertir les NaN en no_data
    abs_image = np.where(np.isnan(abs_image), no_data, abs_image)
    
    # Mettez à jour les métadonnées avec les nouvelles dimensions si nécessaire
    metadata['height'], metadata['width'] = abs_image.shape
    metadata['dtype'] = abs_image.dtype  # Assurez-vous que le type de données correspond à l'image de sortie
    metadata['nodata'] = no_data  # Définir explicitement le nodata à la valeur de l'utilisateur
    
    # Utilisez les métadonnées copiées pour écrire l'image dans un fichier .tif
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
def weighted_blend_overlap(dalle1_path, dalle2_path, weight1, weight2, overlap_bounds, output_path, no_data):
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
		
		# Gérer les no-data - normaliser à la valeur de l'utilisateur
		no_data1 = src1.nodata if src1.nodata is not None else no_data
		no_data2 = src2.nodata if src2.nodata is not None else no_data
		no_data_out = no_data  # Utiliser la valeur de l'utilisateur comme nodata de sortie
		
		# Convertir les nodata des deux images en no_data_out pour uniformiser
		data1 = np.where((data1 == no_data1) | np.isnan(data1), no_data_out, data1)
		data2 = np.where((data2 == no_data2) | np.isnan(data2), no_data_out, data2)
		
		# Masques pour les valeurs valides (exclure no_data et NaN)
		valid1 = (data1 != no_data_out) & ~np.isnan(data1)
		valid2 = (data2 != no_data_out) & ~np.isnan(data2)
		
		# Calculer la moyenne pondérée
		result = np.full_like(data1, no_data_out, dtype=np.float32)
		
		# Cas où les deux valeurs sont valides : moyenne pondérée
		both_valid = valid1 & valid2
		result[both_valid] = data1[both_valid] * weight1[both_valid] + data2[both_valid] * weight2[both_valid]
		
		# Cas où seule la dalle1 est valide
		only1 = valid1 & ~valid2
		result[only1] = data1[only1]
		
		# Cas où seule la dalle2 est valide
		only2 = valid2 & ~valid1
		result[only2] = data2[only2]
		
		# S'assurer qu'il n'y a pas de NaN dans le résultat
		result = np.where(np.isnan(result), no_data_out, result)
		
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
			'nodata': no_data_out,
			'compress': 'lzw'
		})
		
		# Écrire l'image de sortie
		with rasterio.open(output_path, 'w', **metadata) as dst:
			dst.write(result, 1)

#################################################################################################### 
def assemble_horizontal(image_paths, output_path, no_data):
	"""Assemble des images de gauche à droite. 
	Les pixels suivants écrasent les précédents en cas de recouvrement.
	IMPORTANT: Les pixels no-data ne remplacent PAS les pixels valides existants."""
	from rasterio.merge import merge
	from rasterio.warp import reproject, Resampling
	
	# Trier les images par leur position X (left) pour garantir l'ordre de gauche à droite
	image_bounds = []
	for path in image_paths:
		with rasterio.open(path) as src:
			image_bounds.append((src.bounds.left, path))
	
	# Trier par position X (left)
	image_bounds.sort(key=lambda x: x[0])
	sorted_paths = [path for _, path in image_bounds]
	
	# Ouvrir toutes les images dans l'ordre trié
	srcs = [rasterio.open(path) for path in sorted_paths]
	
	try:
		# Calculer les bounds de toutes les images
		all_bounds = [src.bounds for src in srcs]
		minx = min(b.left for b in all_bounds)
		miny = min(b.bottom for b in all_bounds)
		maxx = max(b.right for b in all_bounds)
		maxy = max(b.top for b in all_bounds)
		
		# Utiliser le transform et la résolution de la première image
		first_src = srcs[0]
		transform = first_src.transform
		width = int((maxx - minx) / abs(transform.a))
		height = int((maxy - miny) / abs(transform.e))
		
		# Créer le transform final
		from rasterio.transform import from_bounds
		out_trans = from_bounds(minx, miny, maxx, maxy, width, height)
		
		# Créer un array pour la mosaïque, initialisé avec no-data
		no_data_value = no_data
		mosaic = np.full((1, height, width), no_data_value, dtype=np.float32)
		
		# Assembler les images une par une, de gauche à droite
		# Seuls les pixels VALIDES remplacent les précédents
		for src in srcs:
			# Lire les données de cette image
			data = src.read(1).astype(np.float32)
			
			# Normaliser les no-data
			src_nodata = src.nodata if src.nodata is not None else no_data_value
			data = np.where((data == src_nodata) | np.isnan(data), no_data_value, data)
			
			# Calculer la fenêtre de cette image dans la mosaïque
			window = rasterio.windows.from_bounds(
				src.bounds.left, src.bounds.bottom, src.bounds.right, src.bounds.top,
				out_trans
			)
			
			# Calculer les indices (arrondir pour éviter les problèmes d'alignement)
			row_start = max(0, int(round(window.row_off)))
			row_end = min(height, row_start + data.shape[0])
			col_start = max(0, int(round(window.col_off)))
			col_end = min(width, col_start + data.shape[1])
			
			# Ajuster data si nécessaire
			data_rows = row_end - row_start
			data_cols = col_end - col_start
			
			if data_rows > 0 and data_cols > 0 and data_rows <= data.shape[0] and data_cols <= data.shape[1]:
				# Extraire la partie de data qui correspond
				data_cropped = data[:data_rows, :data_cols]
				
				# Extraire la région correspondante dans la mosaïque
				mosaic_region = mosaic[0, row_start:row_end, col_start:col_end]
				
				# Masque pour les pixels valides de la nouvelle image
				valid_new = (data_cropped != no_data_value) & ~np.isnan(data_cropped)
				
				# Masque pour les pixels no-data existants dans la mosaïque
				existing_nodata = (mosaic_region == no_data_value) | np.isnan(mosaic_region)
				
				# Masque pour les pixels valides existants dans la mosaïque
				existing_valid = ~existing_nodata
				
				# Remplacer seulement si :
				# - Le nouveau pixel est valide (remplace toujours, même un pixel valide existant), OU
				# - Le pixel existant est no-data (on peut le remplacer même par no-data)
				# Mais on ne remplace JAMAIS un pixel valide par un no-data
				replace_mask = valid_new | (existing_nodata & ~valid_new)
				
				# Appliquer le remplacement
				mosaic_region[replace_mask] = data_cropped[replace_mask]
				mosaic[0, row_start:row_end, col_start:col_end] = mosaic_region
		
		# Métadonnées de sortie
		out_meta = first_src.meta.copy()
		out_meta.update({
			'driver': 'GTiff',
			'height': height,
			'width': width,
			'transform': out_trans,
			'nodata': no_data_value,
			'compress': 'lzw'
		})
		
		# Écrire l'image assemblée
		with rasterio.open(output_path, 'w', **out_meta) as dst:
			dst.write(mosaic)
	finally:
		# Fermer toutes les sources
		for src in srcs:
			src.close()

#################################################################################################### 	
def Make_Assemblage_FINAL(chem_out, chem_xing, NbreDalleX, NbreDalleY, RepTra, no_data):
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
				chem_reconstruction,
				no_data
			)
			
	####################################################################################################################################################		
	## on raboute toutes les dalles sur une même rangée		#########################################################################################
	####################################################################################################################################################
				
	## on réassemble tout	
	for y in tqdm(range(NbreDalleY), desc="Raboutage en ligne"):
		
		# Construire la liste des images à assembler: dalles originales + images de transition
		image_paths = []
		
		# Ajouter les dalles originales
		for x in range(NbreDalleX):
			image_paths.append(os.path.join(RepTra, "Dalle_%s_%s" % (x, y), "MASK_%s_%s.tif" % (x, y)))
		
		# Ajouter les images de transition
		for x in range(NbreDalleX-1):
			image_paths.append(os.path.join(RepTra, 'reconstruction_dalle_%s_%s_%s.tif' % (x, x+1, y)))
		
		# Assembler de gauche à droite
		chem_final_tmp = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)
		assemble_horizontal(image_paths, chem_final_tmp, no_data)
			
	####################################################################################################################################################
	## on assemble les rangées entre elles		 #####################################################################################################
	####################################################################################################################################################
				
	## on réassemble tout	
	for y in tqdm(range(NbreDalleY-1), desc="Assemblage final"):
			
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
			chem_reconstruction,
			no_data
		)
			
	####################################################################################################################################################
	## Assemblage final de toutes les lignes		 #####################################################################################################
	####################################################################################################################################################
	
	# Construire la liste des images à assembler: lignes assemblées + images de transition verticales
	image_paths = []
	
	# Ajouter les lignes assemblées
	for y in range(NbreDalleY):
		chem_dalle_y = os.path.join(RepTra, 'reconstruction_dalle_%s.tif' % y)
		image_paths.append(chem_dalle_y)
		
	# Ajouter les images de transition verticales
	for y in range(NbreDalleY-1):
		image_paths.append(os.path.join(RepTra, 'reconstruction_dalle_%s_%s.tif' % (y, y+1)))
	
	# Assemblage final
	assemble_horizontal(image_paths, chem_out, no_data)

# #################################################################################################### 
# def MakeDecoupage_OLD(chem_in, RepTra, NbreDalleX, NbreDalleY, iTailleparcelle, iTailleRecouvrement, iNbreCPU):
	
	# p=Pool(iNbreCPU)
	
	# for x in range(NbreDalleX):
		# for y in range(NbreDalleY):
			# #créer le nom du répertoire
			# RepDalleXY=os.path.join(RepTra,"Dalle_%s_%s"%(x,y))
			# #print(RepDalleXY)
			# #créer le répertoire, s'il n'existe pas déjà
			# if not os.path.isdir(RepDalleXY): os.mkdir(RepDalleXY)
								
			# #Détermination de col_min, col_max, lig_min, lig_max pour la dalle XY
			# colminDalleXY=x*(iTailleparcelle-iTailleRecouvrement)
			# colmaxDalleXY=x*(iTailleparcelle-iTailleRecouvrement)+iTailleparcelle
			# ligminDalleXY=y*(iTailleparcelle-iTailleRecouvrement)
			# ligmaxDalleXY=y*(iTailleparcelle-iTailleRecouvrement)+iTailleparcelle
				
			# #fichier out mns
			# Chem_decoup=os.path.join(RepDalleXY,"IN_%s_%s.tif"%(x,y))
			# #print(Chem_decoup)
				
			# ### crop MNS
			# cmdCROP = "%s -i %s -ci:%s:%s:%s:%s -o %s -n:" % (chem_xing, chem_in, ligminDalleXY, ligmaxDalleXY, colminDalleXY, colmaxDalleXY, Chem_decoup)
			# #print(cmdCROP)
			# #os.system(cmdCROP)
			# p.apply_async(os.system,[cmdCROP])
		
	# p.close(); p.join(); p.terminate()
	
	# return
	
#################################################################################################### 
def MakeDecoupage(chem_in, RepTra, NbreDalleX, NbreDalleY, iTailleparcelle, iTailleRecouvrement, iNbreCPU):
	
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
			#print(Chem_decoup)
				
			### crop MNS
			cmdCROP = "%s -i %s -ci:%s:%s:%s:%s -o %s -n:" % (chem_xing, chem_in, ligminDalleXY, ligmaxDalleXY, colminDalleXY, colmaxDalleXY, Chem_decoup)
			tasks.append(cmdCROP)
			
	# Initialize the pool
	with Pool(processes=iNbreCPU, initializer=init_worker) as pool:
		results = list(tqdm(pool.imap_unordered(os.system, tasks), total=len(tasks), desc="Découpage en parallèle des dalles"))
		
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
	save_ABSOLUTE_image_with_same_geometry(result, chem_out, chem_in, no_data)
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
		
		#### Assemblage final - avec rasterio (remplace xingng)
		# Le fichier d'assemblage est temporaire (fini par _tmp.tif)
		chem_out_tmp = chem_out.replace('.tif', '_tmp.tif')
		Make_Assemblage_FINAL(chem_out_tmp, chem_xing, NbreDalleX, NbreDalleY, RepTra, no_data)
		
		#### Post-traitement 1: Remplacement des valeurs nodata par -9999
		chem_out_clean = chem_out.replace('.tif', '_clean.tif')
		cmd_clean = "%s -i %s -e'I1.1!=I1.1?-9999:I1.1' -o %s" % (chem_xing, chem_out_tmp, chem_out_clean)
		print("Post-traitement 1: Remplacement des nodata par -9999...")
		os.system(cmd_clean)
		
		#### Post-traitement 2: Bouchage des zones avec no data
		chem_out_clean_bouchage = chem_out.replace('.tif', '_clean_bouchage.tif')
		cmd_bouchage = "%s -i %s -FB:2:C:50,1:1:50:1 -o %s -EM=-9999" % (chem_xing, chem_out_clean, chem_out_clean_bouchage)
		print("Post-traitement 2: Bouchage des zones nodata...")
		os.system(cmd_bouchage)
		
		#### Post-traitement 3: Moyenne sur fenêtre glissante 50x50
		# Le fichier final utilise le nom spécifié dans --out
		cmd_moyenne = "%s -i %s -Xm:50:50 -o %s" % (chem_xing, chem_out_clean_bouchage, chem_out)
		print("Post-traitement 3: Moyenne fenêtre glissante 50x50...")
		os.system(cmd_moyenne)
		
		print("Fichier final généré: %s" % chem_out)
		
		print(time.time() - start_time, "seconds")
		
		print('FIN')
				  
	except (RuntimeError, TypeError, NameError):
		print ("ERREUR: ", NameError)


