# Projeto Numérica I - LU com pivoteamento parcial e completo para aumentar a esparsidade

using CSV, DataFrames
using Plots
pyplot(size=(300,200)) # size=(400,300)
    
include("proj_alg.jl") # proj_alg.jl - esse arquivo tem os algoritmos implementados

function main()

    A = readcsv("west0479.csv")
    size(A)

    A = readcsv("west0479.csv")
    (A, p) = LU_pivo(A; tpv = 0);
    (L, U) = sep_LU(A);
    l = @layout [a b]
    plot(spy(L, leg = false, title = "L"), spy(U, leg=false, title = "U"), layout = l) # para LU com tpv = 0, correto
    
    A = readcsv("west0479.csv")
    (A, p) = LU_pivo_mark(A; tpv = 0);
    (L, U) = sep_LU(A);
    l = @layout [a b]
    plot(spy(L, leg = false, title = "L"), spy(U, leg=false, title = "U"), layout = l) # para LU com tpv = 0 usando pivoteamento de Markowitz, correto
    
    A = readcsv("west0479.csv")
    (A, p, q) = LU_pivo(A; tpv = 1);
    (L, U) = sep_LU(A);
    plot(spy(L, leg = false, title = "L"), spy(U, leg=false, title = "U"), layout = l) # para LU com tpv = 0, não está dando correto
    
end

main()