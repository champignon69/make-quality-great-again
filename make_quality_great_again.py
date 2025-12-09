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
def Make_Assemblage_FINAL(chem_out, chem_xing, NbreDalleX, NbreDalleY, RepTra):
	
	####################################################################################################################################################
	####################################################################################################################################################		
	## on assemble tout d'abord ligne par ligne	   #####################################################################################################
	####################################################################################################################################################
	####################################################################################################################################################
	
	####################################################################################################################################################		
	## on raboute tout d'abord 2 dalles côte à côté (en X)	   #########################################################################################
	####################################################################################################################################################	

		
	# Boucle avec barre de progression pour le raboutage des dalles côte à côte
	for y in tqdm(range(NbreDalleY), desc="Raboutage en colonnes"):
		for x in tqdm(range(NbreDalleX-1), desc="Progression en Colonne", leave=False):
			
			## Nom de la dalle courante	
			chem_MASK_QUALITY_dalle_xy=os.path.join(RepTra,"Dalle_%s_%s"%(x,y),"MASK_%s_%s.tif"%(x,y))
			#print(chem_MASK_QUALITY_dalle_xy)
			chem_MASK_QUALITY_dalle_xy_droite=os.path.join(RepTra,"Dalle_%s_%s"%(x+1,y),"MASK_%s_%s.tif"%(x+1,y))
			#print(chem_MASK_QUALITY_dalle_xy_droite)
				
			## IMAGE GAUCHE / IMAGE DROITE <> DALLE_0_0 / DALLE_1_0 !
			chem_tmp=os.path.join(RepTra,'diff_tmp.tif')
			cmd_1="%s -i %s %s -X- -o %s -n:" %(chem_xing,chem_MASK_QUALITY_dalle_xy,chem_MASK_QUALITY_dalle_xy_droite,chem_tmp)
			#print(cmd_1)
			os.system(cmd_1)
				
			cmd_2="%s -i %s -e'C/NC' -o  %s -tf -n:" %(chem_xing,chem_tmp,os.path.join(RepTra,'poids_HG_HD_pour_HD.tif'))
			#print(cmd_2)
			os.system(cmd_2)
			
			cmd_3="%s -i %s -e'1-I1' -o  %s -n:" %(chem_xing,os.path.join(RepTra,'poids_HG_HD_pour_HD.tif'),os.path.join(RepTra,'poids_HG_HD_pour_HG.tif'))
			#print(cmd_3)
			os.system(cmd_3)
				
			liste_info_tmp=GetInfo(chem_xing, os.path.join(RepTra,'poids_HG_HD_pour_HD.tif'))
				
			cmd_final="%s -i %s %s %s %s -e'I1*I2+I3*I4' -o %s -cg:%s:%s:%s:%s -n:" %(chem_xing,chem_MASK_QUALITY_dalle_xy,
																					os.path.join(RepTra,'poids_HG_HD_pour_HG.tif'),
																					chem_MASK_QUALITY_dalle_xy_droite, 
																					os.path.join(RepTra,'poids_HG_HD_pour_HD.tif'),
																					os.path.join(RepTra,'reconstruction_dalle_%s_%s_%s.tif' %(x,x+1,y)),
																					liste_info_tmp[3],
																					liste_info_tmp[4],
																					liste_info_tmp[5],
																					liste_info_tmp[6])
																																									
			#print(cmd_final)
			os.system(cmd_final)
			
	####################################################################################################################################################		
	## on raboute toutes les dalles sur une même rangée		#########################################################################################
	####################################################################################################################################################
				
	## on réassemble tout	
	## for y in range(NbreDalleY):
	for y in tqdm(range(NbreDalleY), desc="Raboutage en ligne"):
		
		### initialisation ligne de commande
		cmd_assemblage_final_par_ligne="%s -i " %chem_xing
		
		for x in range(NbreDalleX):		
			cmd_assemblage_final_par_ligne=cmd_assemblage_final_par_ligne+" "+os.path.join(RepTra,"Dalle_%s_%s"%(x,y),"MASK_%s_%s.tif"%(x,y))
								
		for x in range(NbreDalleX-1):
			cmd_assemblage_final_par_ligne=cmd_assemblage_final_par_ligne+" "+os.path.join(RepTra,'reconstruction_dalle_%s_%s_%s.tif' %(x,x+1,y))
				
		#
		chem_final_tmp=os.path.join(RepTra,'reconstruction_dalle_%s.tif' %y) 
		str_tmp=" -a -o %s -n:" %chem_final_tmp
		cmd_assemblage_final_par_ligne+=str_tmp
			
		#print(cmd_assemblage_final_par_ligne)
		os.system(cmd_assemblage_final_par_ligne)	
			
	####################################################################################################################################################
	####################################################################################################################################################				
	## on assemble les rangées entre elles		 #####################################################################################################
	####################################################################################################################################################
	####################################################################################################################################################
				
	## on réassemble tout	
	## for y in range(NbreDalleY-1):
	for y in tqdm(range(NbreDalleY-1), desc="Assemblage final"):
			
		chem_dalle_y=os.path.join(RepTra,'reconstruction_dalle_%s.tif' %y)			 
		chem_dalle_y_dessous=os.path.join(RepTra,'reconstruction_dalle_%s.tif' %(y+1))
			
		cmd_1="%s -i %s %s -X- -o %s -n:" %(chem_xing,chem_dalle_y,chem_dalle_y_dessous,os.path.join(RepTra,'diff_tmp.tif'))
		#print(cmd_1)
		os.system(cmd_1)
			
		cmd_2="%s -i %s -e'L/NL' -o %s -tf -n:" %(chem_xing,os.path.join(RepTra,'diff_tmp.tif'),os.path.join(RepTra,'poids_HG_BG_pour_BG.tif'))
		#print(cmd_2)
		os.system(cmd_2)
			
		cmd_3="%s -i %s -e'1-I1' -o %s -n:" %(chem_xing,os.path.join(RepTra,'poids_HG_BG_pour_BG.tif'),os.path.join(RepTra,'poids_HG_BG_pour_HG.tif'))
		#print(cmd_3)
		os.system(cmd_3)			
			
		liste_info_tmp=GetInfo(chem_xing, os.path.join(RepTra,'poids_HG_BG_pour_BG.tif'))
			
		cmd_final="%s -i %s %s %s %s -e'I1*I2+I3*I4' -o %s -cg:%s:%s:%s:%s -n:" %(chem_xing,chem_dalle_y,
																				os.path.join(RepTra,'poids_HG_BG_pour_HG.tif'),
																				chem_dalle_y_dessous,os.path.join(RepTra,'poids_HG_BG_pour_BG.tif'),
																				os.path.join(RepTra,
																				'reconstruction_dalle_%s_%s.tif' %(y,y+1)),
																				liste_info_tmp[3],
																				liste_info_tmp[4],
																				liste_info_tmp[5],
																				liste_info_tmp[6])
		#print(cmd_final)
		os.system(cmd_final)
			
	### initialisation ligne de commande
	cmd_assemblage_final="%s -i " %chem_xing
			
	## on réassemble tout	
	for y in range(NbreDalleY):
		chem_dalle_y=os.path.join(RepTra,'reconstruction_dalle_%s.tif' %y)	 
		cmd_assemblage_final=cmd_assemblage_final+" "+chem_dalle_y
			
	## on réassemble tout	
	for y in range(NbreDalleY-1):
		cmd_assemblage_final=cmd_assemblage_final+" "+os.path.join(RepTra,'reconstruction_dalle_%s_%s.tif' %(y,y+1))
		
	str_tmp=" -a -o %s -n:" %chem_out
	cmd_assemblage_final+=str_tmp
			
	#print(cmd_assemblage_final)
	os.system(cmd_assemblage_final)   
		
	# MASK_
	# MASK_
	# MASK_

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
		
		#### Assemblage final - avec xingng = A REMPLACER !
		# Le fichier d'assemblage est temporaire (fini par _tmp.tif)
		chem_out_tmp = chem_out.replace('.tif', '_tmp.tif')
		Make_Assemblage_FINAL(chem_out_tmp, chem_xing, NbreDalleX, NbreDalleY, RepTra)
		
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


