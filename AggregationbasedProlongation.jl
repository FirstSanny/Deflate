module AggregationbasedProlongation
export prolongationAggregation

    function prolongationAggregation(A)
        # Init
        n = size(A,1)
        U = Set(collect(1:n))
        W = Set()
        V = Set()
        tau = 0.25
        q = 0

        # Start Computing
        while true
            q = q + 1

            # find u
            u = -1
            for localU ∈ U
                if(!in(localU,W))
                    u = localU
                    #println("Knoten $u in Iteration $q")
                    break
                end
            end
            # Break if no u can be found
            if(u == -1)
                break
            end


            # Compute NiCTau
            NiCTau = Set()
            for j = 1:n
                if( (A[u,j] != 0) & ( sqrt(A[u,j]^2/abs(A[u,u]*A[j,j])) ≥ tau ))
                    push!(NiCTau, j)
                end
            end
            NiCTau = setdiff(NiCTau, W)
            #println("Nahe Nachbarn $NiCTau in Iteration $q")

            # Compute Vq
            Vq = Set()
            if(length(NiCTau) == 2)
                #println("Nur nahe Nachbarn in Iteration $q")
                push!(NiCTau, u)
                Vq = NiCTau
                push!(V, Vq)
            else
                # Compute NiDTau
                NiDTau = Set()
                for i = 1:n

                    isDistantNeighbor = 1
                    for j in NiCTau
                        if( (A[i,j] == 0) | ( sqrt(A[i,j]^2/abs(A[i,i]*A[j,j])) < tau ))
                            isDistantNeighbor = 0
                        end
                    end
                    if(isDistantNeighbor == 1)
                        push!(NiDTau, i)
                    end
                end
                #println("Distante Nachbarn $NiCTau in Iteration $q")
                push!(NiCTau, u)
                Vq = ∪(NiCTau, NiDTau)
                push!(V, Vq)
            end

            # Set W and U
            W = ∪(W, Vq)
            #println("W in Iteration $q ist $W")
            U = setdiff(U,Vq)
            #println("U in Iteration $q ist $U")

            # Break if there are no more entries in U
            if(length(U) == 0)
                break
            end
        end

        #println("V Vektoren $V")

        Z = zeros(q,n)
        # Compute Interpolation-Matrix
        for i = 1:q
            Vq = pop!(V)
            #println("Vq in Iteration $j ist $Vq")
            for j = 1:n
                if(in(j, Vq))
                    Z[i,j] = 1
                end
            end
        end
        return Z

    end

end
