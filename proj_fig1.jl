# Projeto Numérica I - LU com pivoteamento parcial e completo para aumentar a esparsidade
#                      Comandos usados para gerar os dois gráficos que temoos na figura 1 do relatório do projeto.

using Plots
pyplot(size=(500,400)) # size=(400,300)
using LinearAlgebra
using SparseArrays
    
include("proj_alg.jl") # proj_alg.jl - esse arquivo tem os algoritmos implementados

function main()

        # matriz flecha-n
        lmf = zeros(61)
        umf = zeros(61)
        for (i, n) in enumerate([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300])
            A = zeros(n, n)
            for i in 1:n
                for j in 1:n
                    if i == 1
                        A[i, j] = 1.0
                    elseif j == 1
                        A[i, j] = 1.0
                    elseif i == j
                        A[i, j] = 1.0
                    end
                end
            end
            A[1, 1] = n;
            (A, p, q) = LU_pivo_mark(A; u = 1e-1);
            (L, U) = sep_LU(A);
            lmf[i] = esp(L)
            umf[i] = esp(U)
        end

        # matriz diagonal
        ld = zeros(61)
        ud = zeros(61)
        for (i, n) in enumerate([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300])
            A = zeros(n, n)
            for i in 1:n
                for j in 1:n
                    if i == j
                        A[i, j] = i
                    end
                end
            end
            (A, p, q) = LU_pivo_mark(A; u = 1e-1);
            (L, U) = sep_LU(A);
            ld[i] = esp(L)
            ud[i] = esp(U)
        end

        # matriz triangular superior
        lts = zeros(61)
        uts = zeros(61)
        for (i, n) in enumerate([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300])
            A = zeros(n, n)
            for i in 1:n
                for j in 1:n
                    if i <= j
                        A[i, j] = j
                    end
                end
            end
            (A, p, q) = LU_pivo_mark(A; u = 1e-1);
            (L, U) = sep_LU(A);
            lts[i] = esp(L)
            uts[i] = esp(U)
        end

        # matriz triangular inferior
        lti = zeros(61)
        uti = zeros(61)
        for (i, n) in enumerate([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300])
            A = zeros(n, n)
            for i in 1:n
                for j in 1:n
                    if i >= j
                        A[i, j] = i
                    end
                end
            end
            (A, p, q) = LU_pivo_mark(A; u = 1e-1);
            (L, U) = sep_LU(A);
            lti[i] = esp(L)
            uti[i] = esp(U)
        end

        # matriz tridiagonal
        lt = zeros(61)
        ut = zeros(61)
        for (i, n) in enumerate([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300])
            A = zeros(n, n)
            for i in 1:n
                for j in 1:n
                    if (i == j)
                        A[i, j] = i
                    elseif (j == i + 1) || (j == i - 1)
                        A[i, j] = i - j
                    end
                end
            end
            (A, p, q) = LU_pivo_mark(A; u = 1e-1);
            (L, U) = sep_LU(A);
            lt[i] = esp(L)
            ut[i] = esp(U)
        end

        n = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300]
        plot(title="Ordem da matriz x esparsidade dos fator L", xlabel="n", ylab="esparsidade")
        plot!(n, lmf, c=:green, ms=3, lw=1.5, label="matriz flecha-n")
        plot!(n, ld, c=:blue, ms=3, lw=1.5, label="matriz diagonal")
        plot!(n, lts, c=:black, ms=3, lw=1.5, label="matriz triangular superior")
        plot!(n, lti, c=:red, ms=3, lw=1.5, label="matriz triangular inferior")
        plot!(n, lt, c=:orange, ms=3, lw=1.5, label="matriz matriz tridiagonal")

        n = n = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300]
        plot(title="Ordem da matriz x esparsidade dos fator U", xlabel="n", ylab="esparsidade")
        plot!(n, umf, c=:green, ms=3, lw=1.5, label="matriz flecha-n")
        plot!(n, ud, c=:blue, ms=3, lw=1.5, label="matriz diagonal")
        plot!(n, uts, c=:black, ms=3, lw=1.5, label="matriz triangular superior")
        plot!(n, uti, c=:red, ms=3, lw=1.5, label="matriz triangular inferior")
        plot!(n, ut, c=:orange, ms=3, lw=1.5, label="matriz matriz tridiagonal")
    
end

main()
