# local Gradient Code
# Recibe una imagen de entrada y devuelve una imagen de salida resultado de aplicar el descriptor LGC
import cv2
import numpy as np

def codigo_binario_LGC(matriz):
    codigo = 0
    if ((matriz[0][0]) - (matriz[2][0])) >= 0:
        codigo = codigo + 2**7
    
    if ((matriz[0][1]) - (matriz[2][1])) >= 0:
        codigo = codigo + 2**6

    if ((matriz[0][2]) -( matriz[2][2])) >= 0:
        codigo = codigo + 2**5

    if ((matriz[0][0]) - (matriz[0][2])) >= 0:
        codigo = codigo + 2**4

    if ((matriz[1][0]) - (matriz[1][2])) >= 0:
        codigo = codigo + 2**3

    if ((matriz[2][0]) - (matriz[2][2])) >= 0:
        codigo = codigo + 2**2

    if ((matriz[0][0]) - (matriz[2][2])) >= 0:
        codigo = codigo + 2**1

    if ((matriz[0][2]) - (matriz[2][0])) >= 0:
        codigo = codigo + 2**0
    
    return codigo

def Local_Gradient_Code(imagen, n_pixeles=8, radio=1):
    LGC_Imagen = np.copy(imagen)
    x,y = LGC_Imagen.shape
    
    for i in range(x):
        LGC_Imagen[i][0] = 0
        LGC_Imagen[i][x-1] = 0
    
    for j in range(y):
        LGC_Imagen[0][j] = 0
        LGC_Imagen[y-1][j] = 0
    
    for i in range(x-2):
        for j in range(y-2):
            LGC_Imagen[i+1][j+1] = codigo_binario_LGC(imagen[i:i+3, j:j+3])
        
    return LGC_Imagen