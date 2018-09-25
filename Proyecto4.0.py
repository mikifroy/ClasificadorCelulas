# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 17:30:08 2018

@author: Froylan Miguel Diaz Quintano.
"""


from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from skimage import io, util
from skimage.color import rgb2gray
from skimage import measure
import cv2
import os
from skimage.filters import sobel
from skimage.feature import local_binary_pattern
from sklearn.metrics import confusion_matrix








 #Obtener directorio actual + Directorio de las imágenes en blanco y negro
DatosEntrada1 = os.getcwd() + "\\"  + "hela" + "\\" + "actin" 
DatosEntrada2 = os.getcwd() + "\\" + "hela" + "\\" + "dna"
DatosEntrada3 = os.getcwd() + "\\" + "hela" + "\\" + "endosome"
DatosEntrada4 = os.getcwd() + "\\" + "hela" + "\\" + "er"
DatosEntrada5 = os.getcwd() + "\\" +  "hela" + "\\" + "golgia"
DatosEntrada6 = os.getcwd() + "\\" +  "hela" + "\\" + "golgpp"
DatosEntrada7 = os.getcwd() + "\\" +  "hela" + "\\" + "lysosome"
DatosEntrada8 = os.getcwd() + "\\" +  "hela" + "\\" + "microtubules"
DatosEntrada9 = os.getcwd() + "\\" +  "hela" + "\\" + "mitochondria"    
DatosEntrada10 = os.getcwd() + "\\" + "hela" + "\\" + "nucleolus"
    
ARCHIVOS1 = [path.as_posix() for path in Path(DatosEntrada1).iterdir()]
ARCHIVOS2 = [path.as_posix() for path in Path(DatosEntrada2).iterdir()]
ARCHIVOS3 = [path.as_posix() for path in Path(DatosEntrada3).iterdir()]
ARCHIVOS4 = [path.as_posix() for path in Path(DatosEntrada4).iterdir()]
ARCHIVOS5 = [path.as_posix() for path in Path(DatosEntrada5).iterdir()]
ARCHIVOS6 = [path.as_posix() for path in Path(DatosEntrada6).iterdir()]
ARCHIVOS7 = [path.as_posix() for path in Path(DatosEntrada7).iterdir()]
ARCHIVOS8 = [path.as_posix() for path in Path(DatosEntrada8).iterdir()]
ARCHIVOS9 = [path.as_posix() for path in Path(DatosEntrada9).iterdir()]
ARCHIVOS10 = [path.as_posix() for path in Path(DatosEntrada10).iterdir()]
    
    
    
CLEANING1 = ARCHIVOS1[1:]
CLEANING2 = ARCHIVOS2[1:] 
CLEANING3 = ARCHIVOS3[1:] 
CLEANING4 = ARCHIVOS4[1:] 
CLEANING5 = ARCHIVOS5[1:] 
CLEANING6 = ARCHIVOS6[1:] 
CLEANING7 = ARCHIVOS7[1:] 
CLEANING8 = ARCHIVOS8[1:] 
CLEANING9 = ARCHIVOS9[1:] 
CLEANING10 = ARCHIVOS10[1:] 
    
ACTIN_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 1] for sample in CLEANING1 if "actin" in sample])
DNA_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 2] for sample in CLEANING2 if "DNA" in sample])
ENDOSOME_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 3] for sample in CLEANING3 if "endosome" in sample])
ER_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 4] for sample in CLEANING4 if "ER" in sample])
GOLGIA_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 5] for sample in CLEANING5 if "golgia" in sample])
GOLGPP_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 6] for sample in CLEANING6 if "golgpp" in sample])
LYSOSOME_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 7] for sample in CLEANING7 if "lysosome" in sample])
MICROTUBULES_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 8] for sample in CLEANING8 if "microtubules" in sample])
MITOCHONDRIA_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 9] for sample in CLEANING9 if "mitochondria" in sample])
NUCLEOLUS_CELLS = np.array([[cv2.imread(sample,cv2.IMREAD_GRAYSCALE), 10] for sample in CLEANING10 if "nucleolus" in sample])
#Crea una lista con los nombre y directorios de todos los archivos
ARCHIVOS = np.vstack((ACTIN_CELLS,DNA_CELLS,ENDOSOME_CELLS,ER_CELLS,GOLGIA_CELLS,GOLGPP_CELLS,LYSOSOME_CELLS,MICROTUBULES_CELLS,MITOCHONDRIA_CELLS,NUCLEOLUS_CELLS))
    
  
 #Conjunto de imágenes
 
 


CELULAS = ARCHIVOS[:,0]

#os.mkdir("lbp")
    
#Conjunto de etiquetas de las imágenes de 1-10
y = np.asarray(ARCHIVOS[:,1], dtype="int")
x = np.asarray(CELULAS)
 


histogram = []   


    
contador = 0
    
    
for CELULA in CELULAS:
        
           
      

    #se muestran las imagenes.   
    #if contador == 0:
     #   contador += 1
      #  continue
    contador +=1
    #print("Imagen #{}".format(contador))
    #plt.imshow(CELULA, cmap='gray')
    #plt.show()
        
        
        #calculando el lbp de cada imagen  
    lbp_radius = 1
    lbp_n_points = 8 * lbp_radius
    lbp_METHOD = 'default'

        #mostramos el lbp
    image = util.img_as_ubyte(CELULA)
    lbp= local_binary_pattern(image, lbp_n_points, lbp_radius, lbp_METHOD)
    
    
    #plt.imshow(lbp, cmap="gray")
    #plt.show()
        #print(lbp)
        
    
        
        #mostramos el histograma
    #plt.hist(lbp.ravel(),256,[0,256])
    #plt.ylim([0,4000])
        
        #guardar en conjunto los histogramas
                   
    hist, bins = np.histogram(lbp.ravel(),256,[0,256])
    histogram.append(hist)
    np_histogram = np.array(histogram)
    #print(np_histogram)

    #plt.show()
    #print(histogram)
    
from sklearn.cross_validation import train_test_split
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size = 0.20, random_state=0)


from sklearn.svm import SVC
clasificador = SVC(kernel="rbf", random_state=0)
clasificador.fit(X_train, y_train)

# Predecir en el conjunto de prueba
# Vector de predicciones
y_pred = clasificador.predict(X_test)

print("Matriz de Confusión")
cm = confusion_matrix(y_test,y_pred)


# Imprimir la matriz de confusión
print(cm)
# Calcular y presentar en pantalla la exactitud o precisión del clasificador
#print("Accuracy: ", accuracy_score(y_test, y_pred))











        
    

    


