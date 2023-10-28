import numpy as np
import matplotlib.pyplot as plt

#############Lorenzo Pipa#################
###############27/10/23###################
###########Brasil, Espirito Santo#########

#  Problema de distribuição de temperatura ao longo de um barramento, foi utilizado
#o Metodo de Diferenças finitas para solução 2D do problema.
#  A tranferência de calor ao longo do barramento por condução, e convecção na su-
#percie superior, as laterais esquerda e direita são mantidas a 40°C e 10°C os res-
#tantes das superficies não tem transferência de calor.

##############Constantes do MDF################
inter_max = 80000  #Nº de interações máxima
nvol_x = 100      #Nº de volumes em x
nvol_y = 50       #Nº de volumes em y
tol = 10**-6      #Tolerância mínima
tolmax = 1/tol    #Tolerância máxima

##############Dados de entrada#################
L = 0.1   #Comprimento [m]
H = 0.05  #Altura [m]
t = 0.01  #Espesura [m]

T_inf = 0 + 273   #Temperatura do ar [K]
T_dir = 10 + 273  #Temperatura da parede direita [K]
T_esq = 40 + 273  #Temperatura da parede esquerda [K]

k_cond = 20  #Coeficiente de condução da liga [W/m.K]
h_conv = 75  #Coeficiente de convecção do ar [W/K.m^2]

Q = 10**6    #Taxa de geração de calor [W/m^3]

############Tamanho dos volumes de controle############
dx = L/nvol_x   #Comprimento da coordenada do volume de controle em x [m]
dy = H/nvol_y   #Comprimento da coordenada do volume de controle em x [m]
dx2 = dx*dx
dy2 = dy*dy
aux_dx = dx*h_conv/k_cond
T_q = Q*dx2*dy2/k_cond

#############Matriz da Temperatura###############
T = np.zeros((nvol_y+1,nvol_x+1)) #Matriz zero Nº Volumes em y+1 por Nº Volumes em x+1
T[:,0] = T_esq                    #Temperatura nos pontos da esquerda
T[:,-1] = T_dir                   #Temperatura nos pontos da direita

#####################Malha########################
x = np.zeros(nvol_x+1)       #Vetor x
x[0] = 0
for i in range(1, nvol_x+1): #Preenchimento do vetor x
    x[i] = x[i-1] + dx

y = np.zeros(nvol_y+1)       #Vetor y
y[0] = 0
for j in range(1, nvol_y+1): #Preenchimento do vetor y
    y[j] = y[j-1] + dy

#################Processo interativo################
inter = 1
erro = 1
while (erro >= tol) and (inter <= inter_max) and  (erro <= tolmax):
    soma = 0
    ###Pontos que não estão na borda
    for i in range(1, nvol_x):
        for j in range(1, nvol_y):
            T_old = T[j,i]
            T[j,i] = (dy2*T[j,i+1] + dy2*T[j,i-1]+ dx2*T[j+1,i] + dx2*T[j-1,i]+ T_q)/(2*(dx2+dy2))
            soma = soma + (T[j,i] - T_old)**2
    ###Superficie superior
    for i in range(1, nvol_x):
        T_old =  T[nvol_y,i]
        T[nvol_y,i] = (T[nvol_y-1,i]+aux_dx*T_inf)/(1+aux_dx)
        soma = soma + (T[nvol_y,i] - T_old)**2
    ###Superficie inferior
    for i in range(1, nvol_x):
        T_old =  T[0,i]
        T[0,i] = (dy2*T[0,i+1] + dy2*T[0,i-1]+ dx2*T[1,i])/((dx2+2*dy2))
        soma = soma + (T[0,i] - T_old)**2
    erro = soma**.5
    inter += 1
    print('Interaçao: ',inter,'; Erro: ',erro)

# Criação do gráfico 2D da temperatura
plt.figure(figsize=(12, 5))  # Define o tamanho da figura
plt.contourf(x, y, T, 100, cmap='coolwarm', vmin = T.min(), vmax = T.max())  # Plota o gráfico de contorno
plt.colorbar(label='Temperatura [K]')  # Adiciona uma barra de cores
plt.title('Distribuição de Temperatura ao Longo do Barramento')  # Título do gráfico
plt.xlabel('Comprimento [m]')  # Rótulo do eixo x
plt.ylabel('Altura [m]')  # Rótulo do eixo y
plt.grid(True)  # Ativa as grades no gráfico
plt.show()

print('Temperatura máxima ao longo do barramento:',round(T.max()-273, 1),'°C.')

dQ = h_conv*t*dx*(T[-1,:]-T_inf) #Perda de calor ao longo da superficie superior

#Gráfico da Perda de Calor na superficie
plt.figure(figsize=(10, 5))
plt.plot(x, dQ)
plt.ylabel('$Calor_{perdido}$ [W/m]', fontsize=14)
plt.xlabel('$Comprimento$ [m]' , fontsize=14)
plt.grid(linestyle='--')
plt.ylim(0, 0.06)
plt.xlim(0, 0.1)
plt.title('Distribuição perda de calor ao longo da superficie')
plt.show()