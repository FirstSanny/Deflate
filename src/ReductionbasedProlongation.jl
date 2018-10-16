module ReductionbasedProlongation
export prolongationReduction

using LinearAlgebra

    function computeTauI(A, F, i)
        aij = 0
        for j ∈ F
            aij = aij + abs(A[i,j])
        end

        if(isempty(F))
            return abs(A[i,i])
        end

        return abs(A[i,i]) / aij
    end

    function Adj(A, j, n)
        Adj = Set()
        for k = 1:n
            if(A[j,k] != 0)
                push!(Adj,k)
            end
        end
        return Adj
    end

    function greedyCoarsing(A)
        tau = 1

        n = size(A,1)
        U = Set(collect(1:n))
        F = Set()
        C = Set()
        TauI = zeros(n)

        for i = 1:n
            TauI[i] = computeTauI(A,F,i)
        end

        for i = 1:n
            if(TauI[i] ≥ tau)
                F = ∪(F, i)
                U = setdiff(U, i)
            end
        end

        while(!isempty(U))
            min = Inf
            j = -1
            for i = 1:n
                if(!in(i,U))
                    continue
                end
                if(TauI[i] ≤ min)
                    min = TauI[i]
                    j = i
                end
            end
            C = ∪(C, j)
            U = setdiff(U,j)
            for i ∈ ∩(U, Adj(A,j,n))
                TauI[i] = computeTauI(A,F,i)
                F = ∪(F, i)
                U = setdiff(U, i)
            end
        end
        return F, C, TauI
    end

    function swapColumns(A, i, j)
        if(i == j)
            return A
        end

        m, n = size(A)
        if (1 <= i <= n) && (1 <= j <= n)
            for k = 1:m
                  @inbounds A[k,i],A[k,j] = A[k,j],A[k,i]
            end
            return A
        else
            throw(BoundsError())
        end
    end

    function swapRows(A, i, j)
        if(i == j)
            return A
        end

        m, n = size(A)
        if (1 <= i <= n) && (1 <= j <= n)
            for k = 1:n
                  @inbounds A[i,k],A[j,k] = A[j,k],A[i,k]
            end
            return A
        else
            throw(BoundsError())
        end
    end

    function prolongationReduction(A)
        n = size(A,1)

        FSet, CSet, TauI = greedyCoarsing(A)
        F = collect(FSet)
        C = collect(CSet)
        # Swaps for F
        f = length(F)
        for i in 1:f
            swapRows(A,i,F[i])
            swapColumns(A,i,F[i])
        end

        # Swaps for C
        for i in f+1:n
            # was already swaped
            if in(C[i-f],1:f)
                if(in(i,F) && find(x->x==i,F)[1] == F[C[i-f]])
                    continue
                end
                swapRows(A,i,F[C[i-f]])
                swapColumns(A,i,F[C[i-f]])
            else
                #normal swap
                swapRows(A,i,C[i-f])
                swapColumns(A,i,C[i-f])
            end
        end

        D1 = -1*inv(Diagonal(A[1:f,1:f]))
        TopMatrix = hcat(D1, A[1:f,(f+1):n])
        return Matrix(transpose(vcat(TopMatrix, hcat(zeros(n-f,f), Matrix{Int64}(I, n-f,n-f)))))
    end

end
