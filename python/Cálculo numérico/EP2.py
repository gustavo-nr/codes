#%%
import numpy as np
from numpy.core.fromnumeric import argmin

def main():
    exercicio = str(input("Digite 'i' para rodar a Tarefa 1; digite 'ii', para a Tarefa 2: "))
    if exercicio == 'i':
        teste = str(input("Digite 'a' para analisar a matriz do item a); digite 'b', para a matriz do item b): "))
        teste = 'input-' + teste
        k, Ln, Vn, Av, Lv, dif, A = tarefa_1(teste)
        print("A matriz analisada é: \n", A, "\n")
        print("Foram necessárias ", k, " iterações para convergência")
        print("Os autovalores são: \n", Ln, "\n")
        print("Os autovetores são: \n", Vn, "\n")
        print("A multiplicação da matriz A pelos autovetores é: \n", Av, "\n")
        print("A multiplicação da matriz de autovetores pelos autovalores é: \n", Lv, "\n")
        print("Verificação de ortogonalidade (multiplicação de V pela transposta): \n", Vn @ np.transpose(Vn), "\n")
        analitico = str(input("Digite 'y' para mostrar os autovalores analíticos: "))
    elif exercicio == 'ii':
        k, omega, V, low_val, low_vec, K, M = tarefa_2('input-c')
        print("As 5 menores frequências de vibração são: \n", low_val, "\n")
        print("e os respectivos 5 modos de vibração são: \n", low_vec, "\n")
        np.savetxt("matriz_K.csv", low_vec, delimiter=" , ")


#leitura da matriz  dos arquivos input-a e input-b
def matriz(file):
    matriz= np.loadtxt(file, dtype = 'float', delimiter = ' ',skiprows=1)
    #print(matriz)
    return matriz

def matriz_c(file):
    props = np.loadtxt(file, dtype = 'float', delimiter = ' ',max_rows=2)
    matriz = np.loadtxt(file, dtype = 'float', delimiter = ' ',skiprows=2)
    return props, matriz

#função que retorna a matriz tridiagonal após aplicações sucessivas de Householder
def householder(A):
    n = len(A)

    #inicialização da matriz H transposta
    Ht = np.identity(n)
    #cálculo dos vetores w_barra
    for i in range(1,n-1):
        #a_barra -> primeira coluna da submatriz das últimas n-1 linhas e colunas de A
        a_col = A[i:n,i-1]
        
        #Sinal do primeiro elemento
        delta = a_col[0] / abs(a_col[0])
        
        #vetor e
        e = np.zeros(shape=(n-i))
        e[0] = 1
        
        #vetor w_barra
        w_col = a_col + delta * np.sqrt(np.inner(a_col,a_col)) * e

        #matriz auxiliar que armazenará, momentaneamente, as modificações de A
        A_aux = np.copy(A)

        #multiplicação de H_w_i por A, usando apenas os vetores
        for j in range(i-1,n):
            w = w_col
            #vetor que contém a coluna j da matriz A
            x = A[i:,j]

            #aplicação da transformação à coluna
            x_new = x - 2 * ( np.inner(w,x) ) / ( np.inner(w,w) ) * w

            #armazenamento na matriz auxiliar
            A_aux[i:,j] = x_new

        #aplicação das modificações à matriz original
        A = A_aux

        #aplicação da propriedade de simetria 
        A[i-1,i:] = A[i:,i-1]

        A_aux = np.copy(A)

        #procedimento análogo, mas para as linhas de A
        for l in range(i,n):
            w = w_col
            x = A[l,i:]
            x_new = x - 2 * ( np.inner(x,w) ) / ( np.inner(w,w) ) * np.transpose(w)
            A_aux[l,i:] = x_new
        A = A_aux

        Ht_aux = np.copy(Ht)

        #procedimento análogo, mas para as linhas de Ht
        for l in range(n):
            w = w_col
            x = Ht[l,i:]
            x_new = x - 2 * ( np.inner(x,w) ) / ( np.inner(w,w) ) * np.transpose(w)
            Ht_aux[l,i:] = x_new

        Ht = Ht_aux
    return A, Ht

#função que retorna vetores contendo a diagonal e subdiagonal
def vectors(A):
    n = len(A)
    d_pri = np.zeros(n)  #diagonal principal
    d_sub = np.zeros(n-1) #diagonais abaixo e acima da principal
    for i in range(0,n):
        d_pri[i] = A[i,i]
    for j in range(0,n-1):
        d_sub[j] = A[j+1,j]
    
    return d_pri, d_sub

#função que cria a matriz de rotação de dimensão n
def Qf(n,j,cj,sj):
    Q = np.identity(n)
    #elementos com seno e cosseno
    Q[j][j] = cj
    Q[j][j+1] = -sj
    Q[j+1][j] = sj
    Q[j+1][j+1] = cj
    return Q

#função do EP1 com deslocamento espectral
def QR_deslocamento(M):
    #obtenção da matriz tridiagonal
    A, Ht = householder(M)

    #obtenção dos vetores contendo a diagonal e subdiagonal
    a, b = vectors(A)
    c = np.copy(b)

    n = len(a)
    #listas para guardar cossenos e senos do deslocamento
    Ck = [] 
    Sk = [] 

    #listas para guardar temporariamente o resultado da multiplicação
    a_temp = [] 
    b_temp = []
    c_temp = []    

    #inicialização da matriz V como a matriz H transposta
    V = Ht

    #inicialização de valores
    epsilon = 10**-6
    beta_n1 = 100
    k = 0
    muk = 0

    #matriz de autovalores
    Lambda = np.zeros(shape=(1,n))

    for m in range(n-1,0,-1):
        beta_n1 = abs(b[m-1])
        #checar condição de parada
        while beta_n1 > epsilon:
            #Heurística de Wilkinson  
            if k > 0:
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
                V = V @ np.transpose(Qi)

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


def lambda_analitico(n):
    autos_an = []
    for i in range(1,n+1):
        auto_an = 1/2 * ( 1 - np.cos(((2 * i - 1) * np.pi)/(2 * n + 1))) ** (-1)
        autos_an.append(auto_an)
    return autos_an

#função que recebe as matrizes de propriedades e retorna K e M
def K_M(M1,M2):
    #propriedades das barras
    rho = M1[1][0]
    A   = M1[1][1]
    E   = M1[1][2] * 10 ** 9

    #conversão de graus para radianos
    deg = np.pi /180

    #criação das matrizes K e M
    q = len(M2)
    K = np.zeros(shape=(q,q))
    M = np.zeros(shape=(q,q))

    #preenchimento das matrizes
    for m in range(len(M2)):
        #número dos nós em que a barra está
        i = int(M2[m][0])
        j = int(M2[m][1])

        #direção e comprimento da barra
        th = M2[m][2] * deg
        L =  M2[m][3]

        #coeficiente que multiplica a matriz de rigidez, além do cos e sen do ângulo
        coef = E * A / L
        c = np.cos(th)
        s = np.sin(th)

        #parcela de contribuição da massa para cada nó
        mi = 0.5 * rho * A * L
        
        #definição da matriz de rigidez
        B = [[c**2, c*s, -c**2, -c*s], 
             [c*s, s**2, -c*s, -s**2], 
             [-c**2, -c*s, c**2, c*s], 
             [-c*s, -s**2, c*s, s**2]]
        B = np.array(B)
        K_ij = coef * B

        #vetor com os coeficientes das posições a substituir K
        pos = np.array([2*i-1,2*i,2*j-1,2*j]) - 1

        #montagem da matriz K
        for p in range(len(K_ij)):
            for q in range(len(K_ij)):
                K[pos[p]][pos[q]] = K[pos[p]][pos[q]] + K_ij[p][q]

        #montagem da matriz M
        M[2*i-2][2*i-2] += mi
        M[2*i-1][2*i-1] += mi
        M[2*j-2][2*j-2] += mi
        M[2*j-1][2*j-1] += mi

    #redução das dimensões das matriz K e M de 28 para 24 
    K = cut(K,4)
    M = cut(M,4)
    return K, M

#função para reduzir a dimensão de uma matriz
def cut(M,n):
    N = np.identity(len(M)-n)
    for i in range(0,len(M)-n):
        for j in range(0,len(M)-n):
            N[i][j] = M[i][j]
    return N

#função para obtenção da matriz simétrica definida positiva K~
def K_sim(K,M):
    M_12 = np.copy(M)

    #cálculo da matriz M^(-1/2)
    for i in range(len(M)):
        M_12[i][i] = 1 / (M[i][i])**(1/2)

    #cálculo da matriz K~
    K_s = M_12 @ K @ M_12
    return K_s

#função para calcular os autovetores e autovalores associados ao sistema de treliças
def trelica(M1,M2):
    #obtenção das matriz K e M
    K, M = K_M(M1, M2)

    #obtenção da matriz simétrica definida positiva K~
    K_s = K_sim(K, M)

    #número de iterações, autovalores e autovetores
    k, Lambda, V = QR_deslocamento(K_s)

    return k, Lambda, V, K, M

def tarefa_1(file):
    A = matriz(file)
    [k,Ln,Vn] = QR_deslocamento(A)
    
    n = len(A)
    Av = np.zeros(shape=(n,n))
    Lv = np.zeros(shape=(n,n))
    
    for i in range(0,len(A)):
        Av[:][i] = A @ Vn[:,i]
    
    for j in range(0,len(A)):
        Lv[:][j] = Ln[0,j] * Vn[:,j]
    
    dif = Av - Lv
    
    return k, Ln, Vn, Av, Lv, dif, A

def tarefa_2(file):
    M1, M2 = matriz_c(file)
    k, Lambda, V, K, M = trelica(M1, M2)

    omega = np.sqrt(Lambda)

    L = np.copy(omega)
    low_val = np.zeros(shape=(1,5))
    low_vec = np.zeros(shape=(24,5))
    for i in range(5):
        indice = int(np.argmin(L))
        autovalor = np.min(L)

        low_val[0][i] = autovalor
        low_vec[:,i] = V[:,indice]

        L[0][indice] = 10 ** 9
    return k, omega, V, low_val, low_vec, K, M
main()
