module Solver
export solve

using IterativeSolvers

    # Solves the System Ax=b with de Deflation-Operator R
    function solve(A, b, R, useM3)
        n = size(A,1)

        P = transpose(R)

        Q = P*inv(R*A*P)*R

        if(useM3)
            PN = inv(Diagonal(A))*(I - A * Q) + Q
        else
            PN = I - A * Q + Q
        end


        gmres(A, b; tol = 1/100000000, Pl = lufact(PN), log = true)   
    end
end