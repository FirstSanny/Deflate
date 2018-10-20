module Solver
export solve

using LinearAlgebra
using IterativeSolvers

    # Solves the System Ax=b with de Deflation-Operator R
    function solve(A, b, R, useM3)
        n = size(A,1)

        P = transpose(R)

        Q = []
        try
            Q = P*inv(R*A*P)*R
        catch
            println("---pinv zur Berechnung von Q benutzt")
            Q = P*pinv(R*A*P)*R
        end

        if(useM3)
            InvDiagA = []
            try
                InvDiagA = inv(Diagonal(A))
            catch
                println("---pinv zur Berechnung der Inversen der Diagonalen von A benutzt")
                InvDiagA = pinv(Diagonal(A))
            end
            PN = InvDiagA*(I - A * Q) + Q
        else
            PN = I - A * Q + Q
        end


        IterativeSolvers.gmres(A, b; tol = 1/100000000, Pl = LinearAlgebra.lu(PN, check = false), log = true)
    end
end
