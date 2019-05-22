# Arquivo com as implementações dos algoritmos usados nas fatorações LU com pivoteamento parcial e completo e com combinação com pivoteamento de Markowitz para aumentar a espasidade.

"""`LU_pivo_parcial(A)

O algoritmo calcula a fatoração LU usando pivoteamento parcial, para quaisquer matrizes n x n. 

"""
function LU_pivo_parcial(A)
    n = size(A)[1]
    
    ### Cálculo dos fatores:
    p = zeros(n)
    for i in 1:n
        p[i] = i # para guardar os índices das linhas que foram trocadas
    end
    
    for k in 1:n-1
        pv = abs(A[k, k])
        r = k
        for i in k+1:n
            if abs(A[i, k]) > pv
                pv = abs(A[i, k])
                r = i
            end
        end
        
        if pv == 0
            println("Matriz singular")
            break # parar; a matriz A é singular
        end
        
        if r != k
            aux = p[k]
            p[k] = p[r]
            p[r] = aux
            for j in 1:n
                aux = A[k, j]
                A[k, j] = A[r, j]
                A[r, j] = aux
            end
        end
        
        for i in k+1:n
            m = A[i, k] / A[k, k]
            A[i, k] = m
            for j in k+1:n
                A[i, j] = A[i, j] - m * A[k, j]
            end
        end
    end
    
    return A, p
end

"""`LU_pivo_completo(A)

O algoritmo calcula a fatoração LU usando pivoteamento completo (total), para quaisquer matrizes n x n. 

"""
function LU_pivo_completo(A)
    n = size(A)[1]
    
    ### Cálculo dos fatores:
    p = zeros(n)
    q = zeros(n)
    for i in 1:n
        p[i] = i # para guardar os índices das linhas que foram trocadas
        q[i] = i # para guardar os índices das colunas que foram trocadas
    end
    
    for k in 1:n-1
        pv = abs(A[k, k])
        r = k
        t = k
        for i in k+1:n
            for j in k:n
                if abs(A[i, k]) > pv
                    pv = abs(A[i, k])
                    r = i
                end
                if abs(A[i, j]) > pv
                    pv = abs(A[i, j])
                    t = j
                end
            end
        end
        
        if pv == 0
            println("Matriz singular")
            break # parar; a matriz A é singular
        end
        
        if r != k
            aux = p[k]
            p[k] = p[r]
            p[r] = aux
            for j in 1:n
                aux = A[k, j]
                A[k, j] = A[r, j]
                A[r, j] = aux
            end
        end
        if t != k
            aux = q[k]
            q[k] = q[t]
            q[t] = aux
            for j in 1:n
                aux = A[j, k]
                A[j, k] = A[j, t]
                A[j, t] = aux
            end
        end
        
        for i in k+1:n
            m = A[i, k] / A[k, k]
            A[i, k] = m
            for j in k+1:n
                A[i, j] = A[i, j] - m * A[k, j]
            end
        end
    end
    
    return A, p, q
end

# LU com pivoteamenteamento parcial e completo numa única função usando tpv
"""`LU_pivo(A; tpv)

O algoritmo calcula a fatoração LU usando pivoteamento parcial (tpv = 0) ou completo (total) (tpv = 1), para quaisquer matrizes n x n.
Sendo tpv o tipo de pivoteamento, ou seja, tpv = 0 pivoteamento parcial e tpv = 1 pivoteamento completo (total).

"""
function LU_pivo(A; tpv = 0)
    n = size(A)[1]
    
    ### Cálculo dos fatores:
    p = zeros(n)
    q = zeros(n)
    for i in 1:n
        p[i] = i # para guardar os índices das linhas que foram trocadas
        q[i] = i # para guardar os índices das colunas que foram trocadas
    end
    
    for k in 1:n-1
        pv = abs(A[k, k])
        r = k
        t = k
        for i in k+1:n
            for j in k:n
                if abs(A[i, k]) > pv && (tpv == 0)
                    pv = abs(A[i, k])
                    r = i
                end
                if (abs(A[i, j]) > pv) && (tpv == 1)
                    pv = abs(A[i, j])
                    t = j
                end
            end
        end
        
        if pv == 0
            println("Matriz singular")
            break # parar; a matriz A é singular
        end
        
        if r != k
            aux = p[k]
            p[k] = p[r]
            p[r] = aux
            for j in 1:n
                aux = A[k, j]
                A[k, j] = A[r, j]
                A[r, j] = aux
            end
        end
        if (t != k) && (tpv == 1)
            aux = q[k]
            q[k] = q[t]
            q[t] = aux
            for j in 1:n
                aux = A[j, k]
                A[j, k] = A[j, t]
                A[j, t] = aux
            end
        end
        
        for i in k+1:n
            # println("A[k, k]: ", A[k, k])
            m = A[i, k] / A[k, k]
            # println("m: ", m)
            A[i, k] = m
            for j in k+1:n
                A[i, j] = A[i, j] - m * A[k, j]
            end
        end
    end
    
    if tpv == 1
        return A, p, q
    else
        return A, p
    end
end

"""`LU_resol_sist(A, p, b)

O algoritmo a partir da fatoração LU, A no formato LU, resolve o sistema Ax = b, sendo A uma matriz n x n. 

Primeiramente resolvemos c = Pb, depois Ly = c e por último Ux = y.
"""
function LU_resol_sist(A, p, b)
    n = size(A)[1]
    
    ### Resolução dos sistemas triangulares:
    # c = Pb
    c = zeros(n)
    for i in 1:n
        r = Int(p[i])
        c[i] = b[r]
    end
    
    # Ly = c
    y = zeros(n)
    for i in 1:n
        soma = 0
        for j in 1:i-1
            soma += A[i, j] * y[j]
        end
        y[i] = c[i] - soma
    end
    
    # Ux = y
    x = zeros(n)
    for i in n:-1:1
        soma = 0
        for j in i+1:n
            soma += A[i,j] * x[j]
        end
        x[i] = (y[i] - soma) / A[i, i]
    end
    
    return x
end

"""`sep_LU(A)

O algoritmo a partir A no formato LU, separa a matriz A em duas matrizes L e U.
"""
function sep_LU(A)
    n = size(A)[1]
    
    ### Separando a matriz A nas matrizes L e U:
    L = tril(A, -1) + diagm(ones(n), 0)
    U = triu(A)
    
    return L, U
end

# LU com pivoteamento completo combinado com pivoteamento de Markowitz para aumentar esparsidade.
"""`LU_pivo_mark(A)

O algoritmo calcula a fatoração LU usando pivoteamento completo (total), para quaisquer matrizes n x n.
Combinando com o pivoteamento de Markowitz.

"""
function LU_pivo_mark(A; u = 1e-4)
    n = size(A)[1]
    
    ### Cálculo dos fatores:
    p = zeros(n)
    q = zeros(n)
    for i in 1:n
        p[i] = i # para guardar os índices das linhas que foram trocadas
        q[i] = i # para guardar os índices das colunas que foram trocadas
    end
    
    pv = 0
    for k in 1:n # 1:n-1
        cont = n^2
        pv = abs(A[k, k])
        r = k # linha
        t = k # coluna
        
        ### Minimiza o número de elementos não nulo na linha i e coluna j sujeito ao pivoteamento de Markowitz:
        for i in k:n # k+1:n
            for j in k:n # nova 1:i-1 # k:n-1
                (nzlin, nzcol) = zeros_mat(A, i, j)
                if cont >= (nzlin * nzcol) # nnz_n <= nnz
                    if abs(A[i, j]) >= u * abs(A[j, j])
                        pv = abs(A[i, j])
                        cont = nzlin * nzcol
                        r = i
                        t = j
                    end
                end
            end
        end

        if pv == 0
            break # parar; a matriz A é singular
        end
        
        if r != k
            aux = p[k]
            p[k] = p[r]
            p[r] = aux
            for j in 1:n
                aux = A[k, j]
                A[k, j] = A[r, j]
                A[r, j] = aux
            end
        end
        if t != k
            aux = q[k]
            q[k] = q[t]
            q[t] = aux
            for j in 1:n
                aux = A[j, k]
                A[j, k] = A[j, t]
                A[j, t] = aux
            end
        end
        
        for i in k+1:n
            m = A[i, k] / A[k, k]
            A[i, k] = m
            for j in k+1:n
                A[i, j] = A[i, j] - m * A[k, j]
            end
        end
    end
    
    if pv == 0
        println("Matriz singular")
    else
        return A, p, q
    end
    
end
