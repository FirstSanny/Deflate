module simpleProlongations
export prolongation1
export prolongation2

    function  prolongation1(n)
        r = floor(Int, n/2)
        R = zeros(r,n)
        i = 1
        j = 1
        while i ≤ r
            R[i,j] = 1
            i = i + 1
            j = j + 2
        end
        return R
    end

    function  prolongation2(n)
        r = floor(Int, n/2)
        R = zeros(r,n)
        i = 1
        j = 1
        while i ≤ r
            R[i,j] = 1
            if(j+1 ≤ n)
                R[i,j+1] = 2
            end
            if(j+2 ≤ n)
                R[i,j+2] = 1
            end
            i = i + 1
            j = j + 2
        end
        return 1/4*R
    end

end
