## Bibliotecas utilizadas
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from math import sqrt
from math import exp


## Declaração de variáveis
frequenciaDeCorte = 60      # frequência da onda
u = 5                       # constante de malha (mi)
k = 5                       # número máximo de amostra de comprimento de onda
velocidadeMinima = 1500     # velocidade mínima de onda
velocidadeMaxima = 5500     # velocidade máxima de onda
Nx = 383                    # tamanho no eixo x 
Nz = 121                    # tamanho do eixo z (profundidade)
a = 1                       # contador do loop
posicaoFonteX = 192         # posição da fonte no eixo x (x - 1)
posicaoFonteZ = 2           # posição da fonte no eixo z (z -1)

## Definição de funções

# Importação da matriz Marmousi do arquivo binário
def marmousi():
    data = np.fromfile("marmousi_vp_383x121.bin", dtype=np.float32).reshape((383, 121)).T
    plt.imshow(data)
    return data

# Cálculo do h para calcular o DT
def dtcalculoh(velocidadeMinima, k, frequenciaDeCorte):
    h = velocidadeMinima / (k * frequenciaDeCorte)
    return h
    
# Cálculo do DT
def dtcalculo(u):
    h = dtcalculoh(velocidadeMinima, k, frequenciaDeCorte)
    variacaoTempo = h / (u * velocidadeMaxima)
    return variacaoTempo

# Cálculo do valor da função
def fonte(frequenciaDeCorte, n):
    # Cálculo da função fonte proposta por CUNHA (1997)
    # Derivada segunda da gaussiana
    variacaoTempo = dtcalculo(u) 
    TF = 2 * sqrt(pi) / frequenciaDeCorte       # Período da função Gaussiana
    fc = frequenciaDeCorte / (3. * sqrt(pi))    # Frequência central

    t = ((n - 1) * variacaoTempo - TF)

    fonte = (-exp(-pi * (pi * fc * t) ** 2) * (1. - 2. * pi * (pi * fc * t) * (pi * fc * t)))

    return fonte

# Obtenção dos passos de tempo totais
def ntotal(Nx, Nz, h, velocidadeMaxima, velocidadeMinima, variacaoTempo):
    variacaoEspaco = sqrt(((Nx * h) ** 2) + ((Nz * h) ** 2))
    velocidadeMedia = (velocidadeMaxima + velocidadeMinima) / 2
    tempo = variacaoEspaco / velocidadeMedia
    ntotal = int(tempo / variacaoTempo)
    return ntotal

# Cálculo dos passos de tempo totais
def Nf(frequenciaDeCorte, variacaoTempo):
    Nf = 4 * sqrt(pi) / (frequenciaDeCorte * variacaoTempo)
    return Nf

## Criando matrizes zeradas
P1 = np.zeros((Nz, Nx))
P2 = np.zeros((Nz, Nx))
P3 = np.zeros((Nz, Nx))
C = np.zeros((Nz, Nx))
A = np.zeros((Nz, Nx))

## Tranferência do valor da matriz data para marmousitrix
marmousitrix = marmousi() 

## Obtenção dos valores de DT, Nf e h fora das funções
variacaoTempo = dtcalculo(u)
h = dtcalculoh(velocidadeMinima, k, frequenciaDeCorte)
Nf = int(Nf(frequenciaDeCorte, variacaoTempo))
ntotal = ntotal(Nx, Nz, h, velocidadeMaxima, velocidadeMinima, variacaoTempo)

## Cálculo das matrizes C e A
for i in range(Nz):
    for k in range(Nx):
        C[i, k] = - (marmousitrix[i, k] * (variacaoTempo/h)) ** 2

for i in range(Nz):
    for k in range(Nx):
        A[i, k] = marmousitrix[i, k] * (variacaoTempo/h)

## Computação do campo de pressão
for n in range(1, ntotal + 2):
    
    # Termo fonte
    if n <= Nf:
        P2[posicaoFonteZ, posicaoFonteX] = P2[posicaoFonteZ, posicaoFonteX] + fonte(frequenciaDeCorte, n)
           
    # Cálculo do Campo no interior do modelo
    for k in range(3, Nx - 2):
        for i in range(3, Nz - 2):
            P3[i, k] = C[i, k] * (P2[i + 2, k] + P2[i - 2, k] + P2[i, k + 2] + P2[i, k - 2]
                                  - 16 * (P2[i + 1, k] + P2[i - 1, k] + P2[i, k + 1] + P2[i, k - 1]) + 60 * (
                                      P2[i, k])) + 2 * (P2[i, k]) - P1[i, k]
            
    # Atualização do campo de onda
    P2 = P3
    P1 = P2
        
    print (a)
    a += 1
    
plt.imshow(P3)
