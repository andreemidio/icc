## Bibliotecas utilizadas
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from math import sqrt
from math import exp

## Declaração de variáveis
fcorte = 60
u = 5
k = 5
Vmin = 1500
Vmax = 5500
Nx = 383    # tamanho no eixo x 
Nz = 121    # tamanho do eixo z (profundidade)

## Definição de funções

# Importação da matriz Marmousi do arquivo binário
def marmousi():
    data = np.fromfile("marmousi_vp_383x121.bin", dtype=np.float32).reshape((383, 121)).T
    plt.imshow(data)
    return data

# Cálculo do h para calcular o DT
def dtcalculoh(Vmin, k, fcorte):
    h = Vmin / (k * fcorte)
    return h
    
# Cálculo do DT
def dtcalculo(u):
    h = dtcalculoh(Vmin, k, fcorte)
    DT = h / (u * Vmax)
    return DT

# Cálculo do valor da função
def fonte(fcorte, n):
    # Cálculo da função fonte proposta por CUNHA (1997)
    # Derivada segunda da gaussiana
    DT = dtcalculo(u) 
    TF = 2 * sqrt(pi) / fcorte  # Período da função Gaussiana
    fc = fcorte / (3. * sqrt(pi))  # Frequência central

    t = ((n - 1) * DT - TF)

    fonte = (-exp(-pi * (pi * fc * t) ** 2) * (1. - 2. * pi * (pi * fc * t) * (pi * fc * t)))

    return fonte

# Cálculo dos passos de tempo totais
def Nf(fcorte, DT):
    Nf = 4 * sqrt(pi) / (fcorte * DT)
    return Nf

## Criando matrizes zeradas
P1 = np.zeros((Nz, Nx))
P2 = np.zeros((Nz, Nx))
P3 = np.zeros((Nz, Nx))
C = np.zeros((Nz, Nx))
A = np.zeros((Nz, Nx))

## Tranferência do valor da matriz data para marmousitrix
marmousitrix = marmousi() 

## Obtenção dos valores de DT e h fora das funções
DT = dtcalculo(u)
h = dtcalculoh(Vmin, k, fcorte)

## Cálculo das matrizes C e A
for i in range(Nz):
    for k in range(Nx):
        C[i, k] = - (marmousitrix[i, k] * (DT/h)) ** 2
print(C)
print(DT)
for i in range(Nz):
    for k in range(Nx):
        A[i, k] = marmousitrix[i, k] * (DT/h)
print(A)

## Obtenção dos passos de tempo totais
Nf = int(Nf(fcorte, DT))

## Computação do campo de pressão
#a = 1
#for n in range(Nf + 2):

    # Termo fonte
    #if n <= Nf + 1:
       # P2[3, int((Nx / 2)) + 1] = P2[3, int((Nx / 2) + 1)] + fonte(60, n)

    # Cálculo do Campo no interior do modelo
    #for k in range(3, Nx - 2):
       # for i in range(3, Nz - 2):
           # P3[i, k] = C[i, k] * (P2[i + 2, k] + P2[i - 2, k] + P2[i, k + 2] + P2[i, k - 2]
                               #   - 16 * (P2[i + 1, k] + P2[i - 1, k] + P2[i, k + 1] + P2[i, k - 1]) + 60 * (
                                 #     P2[i, k])) + 2 * (P2[i, k]) - P1[i, k]
            
    # Atualização do campo de onda
   # P2 = P3
   # P1 = P2

   # print (a)
   # a = a + 1
    
#plt.imshow(P3)
