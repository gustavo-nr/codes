#%%
import numpy as np
import matplotlib.pyplot as plt
from time import time

def main():
    exercicio = str(input("Digite i para rodar o item a); ii, para o b); e iii, para o c): "))
    modo = int(input("Digite 0 para algoritmo QR sem deslocamento espectral ou 1 para QR com deslocamento: "))
    if modo == 0:
        algoritmo = 'QR sem deslocamento espectral'
    elif modo == 1:
        algoritmo = 'QR com deslocamento espectral'
    if exercicio == 'i':
        n = int(input("Digite a dimensão da matriz a ser trabalhada: "))
        [kn, An, Vn] = QR_deslocamento('a',n,modo)
        print("Algoritmo: ", algoritmo)
        print("A matriz analisada tem dimensão:", "\n", n, "\n")
        print("Os autovalores determinados são ", "\n", An, "\n")
        print("A matriz contendo os autovetores é","\n", Vn,"\n")
        print("Foram necessárias ", kn," iterações")
        analitico = int(input("Digite 1 caso queira os valores analíticos de autovalores e autovetores: "))
        if analitico == 1:
            print("Os autovalores analíticos são ", "\n", eigenvalue(n), "\n")
            print("A matriz contendo os autovetores analíticos é","\n", normalize(eigenvector(n)),"\n")
    elif exercicio == 'ii':
        caso = int(input("Digite o caso a ser analisado (1, 2 ou 3): "))
        save = int(input("Caso queira salvar os gráficos, digite 1: "))
        [kn, An, Vn] = mola('b',caso,save,modo)
        print("Algoritmo: ", algoritmo)
        print("As frequências do sistema são ", "\n", An, "\n")
        print("Os modos de vibração são","\n", Vn,"\n")
        print("Foram necessárias ", kn," iterações para convergir")
    elif exercicio == 'iii':
        caso = int(input("Digite o caso a ser analisado (1, 2 ou 3): "))
        save = int(input("Caso queira salvar os gráficos, digite 1: "))
        [kn, An, Vn] = mola('c',caso,save,modo)
        print("Algoritmo: ", algoritmo)
        print("As frequências do sistema são ", "\n", An, "\n")
        print("Os modos de vibração são","\n", Vn,"\n")
        print("Foram necessárias ", kn," iterações para convergir")

#Função para criar os vetores associados a cada exercício
def matriz(exercicio,n):
    if exercicio == 'a':
        alpha = 2
        beta = -1
        #vetor da diagonal principal
        a = np.full(n,alpha).tolist()
        #vetores das subdiagonais
        b = np.full(n-1,beta).tolist()
        c = np.full(n-1,beta).tolist()
        return a, b, c

    elif exercicio == 'b':
        #valor da massa
        m = 2
        #vetor da diagonal principal
        a = np.full(n,0).tolist()
        #vetores das subdiagonais
        b = np.full(n-1,0).tolist()
        c = np.full(n-1,0).tolist()
        for i in range(1,n+1,1):
            a[i-1] = 1/m * (k_b(i) + k_b(i+1))
        for j in range(1,n,1):
            b[j-1] = - 1/m * k_b(j+1)
            c[j-1] = - 1/m * k_b(j+1)
        return a, b, c
    
    elif exercicio == 'c':
        #valor da massa
        m = 2
        #vetor da diagonal principal
        a = np.full(n,0).tolist()
        #vetores das subdiagonais
        b = np.full(n-1,0).tolist()
        c = np.full(n-1,0).tolist()
        for i in range(1,n+1,1):
            a[i-1] = 1/m * (k_c(i) + k_c(i+1))
        for j in range(1,n,1):
            b[j-1] = - 1/m * k_c(j+1)
            c[j-1] = - 1/m * k_c(j+1)
        return a, b, c

#Função sinal, que retorna o sinal associado a um valor
def sgn(x):
  if x>= 0:
    return 1
  else:
    return -1

#Função que cria a matriz de rotação de dimensão n
def Qf(n,j,cj,sj):
  Q = identidade(n)
  #elementos com seno e cosseno
  Q[j][j] = cj
  Q[j][j+1] = -sj
  Q[j+1][j] = sj
  Q[j+1][j+1] = cj
  return Q

#Função para transpor uma matriz
def transpose(M):
  n = len(M)
  N = np.zeros(shape=(n,n))
  for i in range(0,n,1):
    for j in range(0,n,1):
      N[i][j] = M[j][i]
  return N

#Função que cria uma matriz identidade de dimensão n
def identidade(n):
  I = np.zeros(shape=(n,n))
  for i in range(0,n):
    I[i][i] = 1
  return I

#Função que calcula os autovalores analíticos para o item a)
def eigenvalue(n):
  A = np.zeros(shape=(n,1))
  for j in range(1,n+1):
    A[n-j][0] = 2*(1-np.cos(j*np.pi/(n+1)))
  return A

#Função que calcula os autovetores analíticos para o item b)
def eigenvector(n):
  V = np.zeros(shape=(n,n))
  for j in range(1,len(V)+1,1):
    for i in range (1,len(V)+1,1):
      V[i-1][n-j] = np.sin(i*j*np.pi/(n+1))
  return V

#Função que normaliza os autovetores dispostos em uma matriz
def normalize(M):
  v_soma = np.zeros(len(M))
  for j in range(1,len(M)+1,1):
    for i in range(1,len(M)+1,1):
      v_soma[i-1] = M[i-1][j-1]
    norm = np.sqrt(v_soma @ v_soma)
    M[:,j-1] = 1/norm * M[:,j-1]

  return M

#Função que calcula o algoritmo QR
#exercicio = item a), b) ou c)
#n = dimensão da matriz
#desloc = parâmetro para determinar se é com ou sem deslocamento
def QR_deslocamento(exercicio,n,desloc):
    a, b, c = matriz(exercicio,n)

    #listas para guardar cossenos e senos do deslocamento
    Ck = [] 
    Sk = [] 

    #listas para guardar temporariamente o resultado da multiplicação
    a_temp = [] 
    b_temp = []
    c_temp = []    

    #inicialização de valores e da matriz V
    V = identidade(n)
    epsilon = 10**-6
    beta_n1 = 100
    k = 0
    muk = 0
    #Matriz de autovalores
    Lambda = np.zeros(shape=(1,n))

    for m in range(n-1,0,-1):
        beta_n1 = abs(b[m-1])
        #checar condição de parada
        while beta_n1 > epsilon:
            #Heurística de Wilkinson  
            if k > 0 and desloc == 1:
                dk = (a[m-1]-a[m])/2
                muk = a[m] + dk-np.sign(dk)*(dk**2 + (b[m-1])**2)**(1/2)
                a = a-muk*np.ones(m+1)

        #Determinação de R
            for i in range(len(a)-1):
            #Cálculo dos senos e cossenos
                if abs(a[i]) > abs(b[i]):
                    tau = -b[i]/a[i]
                    ci = 1/(1+tau**2)**(1/2)
                    si = ci*tau
                    Ck.append(ci)
                    Sk.append(si) 
                else:
                    tau = -a[i]/b[i]
                    si = 1/(1+tau**2)**(1/2)
                    ci = si*tau
                    Ck.append(ci)
                    Sk.append(si)

                #Resultado da multiplicação da matriz Q por a, sem lidar com os 0s das matrizes
                a_temp.append([a[i] * ci - b[i] * si, c[i] * si + a[i+1] * ci])
                b_temp.append(a[i] * si + b[i] * ci)

                #Verifica se foi atingida a borda da matriz (fim do vetor c)
                if i == n-2:
                    c_temp.append([c[i] * ci - a[i+1] * si, 0])
                else:
                    c_temp.append([c[i] * ci - a[i+1] * si, c[i+1] * ci])
                
                #Substituição dos valores nos vetores originais
                a[i] = a_temp[i][0]
                a[i+1] = a_temp[i][1]
                b[i] = b_temp[i]
                c[i] = c_temp[i][0]
                if i < n-2:
                    c[i+1] = c_temp[i][1]

                #Cálculo da matriz de rotação de Givens e V
                Qi = Qf(n,i,ci,si)
                V = V @ transpose(Qi)

            for i in range(len(a)-1):
                a[i] = a[i]*Ck[i] - c[i]*Sk[i]
                b[i] = -a[i+1]*Sk[i]
                c[i] = b[i]
                a[i+1] = a[i+1]*Ck[i]
            

            beta_n1 = abs(b[m-1])

            k += 1
            #Restart dos vetores temporários para próxima iteração
            a_temp = []
            b_temp = []
            c_temp = []
            Ck = []
            Sk = [] 
            a = a + muk*np.ones(m+1)
        #Salvamento do autovalor e redução da dimensão dos vetores
        Lambda[0][m] = a[m]
        a = a[:m]
        b = b[:m]
        c = c[:m]
    Lambda[0][0] = a[0]
    return k, Lambda, V

#Função para calcular os coeficientes das molas no caso b)
def k_b(i):
  k = 40 + 2*i

  return k

#Função para calcular os coeficientes das molas no caso c)
def k_c(i):
  k = 40 + 2*(-1)**i

  return k

#Função para calcular o contexto
def mola(item,case,save,desloc):
    #Vetor de instantes de tempo
    h = 0.025
    tf = 10
    t = np.arange(0,tf+h,h)

    #dimensão da matriz
    if item == 'b':
        dim = 5
    elif item =='c':
        dim = 10

    #Cálculo dos autovetores e autovalores
    [kn, An, Vn] = QR_deslocamento(item,dim,desloc)

    if item == 'b':
        X0 = [[-2,-3,-1,-3,-1],[1,10,-4,3,-2],[Vn[0][0],Vn[1][0],Vn[2][0],Vn[3][0],Vn[4][0]]]
    elif item == 'c':
        X0 = [[-2,-3,-1,-3,-1,-2,-3,-1,-3,-1],
        [1,10,-4,3,-2,1,10,-4,3,-2],
        [Vn[0][0],Vn[1][0],Vn[2][0],Vn[3][0],Vn[4][0],Vn[5][0],Vn[6][0],Vn[7][0],Vn[8][0],Vn[9][0]]]

    n = dim
    N = len(t)
    #matriz de soluções
    y = np.zeros(shape=(n,N))
    #condições iniciais
    y0 = transpose(Vn) @ X0[case-1]
    caso = str(case)
    #cálculo da solução
    for i in range(n):
        y[i] = y0[i] * np.cos(np.sqrt(An[0][i])*t)
    #Retorno às variáveis x
    x = Vn @ y

    #plots
    for i in range(n):
        passo=str(i+1)
        plt.figure(i)
        plt.plot(t,x[i],'b')
        title='Posição do bloco '+passo+' em função do tempo'
        plt.title(title)
        plt.grid()
        plt.xlabel('Tempo (s)')
        plt.ylabel('Posição (m)')
        if save == 1:
            label='x_'+passo+'_caso_'+caso
            plt.savefig(label)
    plt.show()
    omega = np.sqrt(An)
    return kn, omega, Vn
main()