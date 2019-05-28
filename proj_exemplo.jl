# Projeto Numérica I - LU com pivoteamento parcial e completo para aumentar a esparsidade
#                      Exemplo usado no texto do projeto.

using Plots
pyplot(size=(500,400)) # size=(400,300)
using LinearAlgebra
using SparseArrays

function main()
    A = [2.0 0.0 2.0 0.0 0.0; 3.0 1.0 4.0 0.0 0.0; 0.0 0.0 -3.0 2.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 4.0 2.0 2.0 1.0]
    println("Matriz A: ")
    println(A)
    println(" ")
    
    println("Esparsidade usando LU com pivoteamento completo.")
    (A, p, q) = LU_pivo(A; tpv = 1);
    (L, U) = sep_LU(A);
    println("Esparsidade da Matriz L: ", esp(L))
    println("Esparsidade da Matriz U: ", esp(U))
    println(" ")
    
    println("Esparsidade usando LU esparso.")
    (A, p, q) = LU_pivo_mark(A);
    (L, U) = sep_LU(A);
    println("Esparsidade da Matriz L: ", esp(L))
    println("Esparsidade da Matriz U: ", esp(U))
end

main()
