module Tester
export testAllDiagonal
export testOneDiagonal

    using DataFrames
    using Plots
    using MatrixMarket
    using TimerOutputs

    include("SimpleProlongations.jl")
    include("AggregationbasedProlongation.jl")
    include("ReductionbasedProlongation.jl")
    include("Solver.jl")

    function getA()
        MatrixMarket.mmread("matrices/bcsstk23.mtx")
    end

    function getTestMatrices()
            return ["bcsstm26", "bcspwr07", "bcspwr09", "bcsstk23", "bcsstk24", "bcsstk25", "dwt__234", "dwt__992", "gemat11", "gr_30_30", "jgl009", "lap___25", "mbeause", "nnc666", "nos7", "s3dkt3m2", "saylr3", "watt__1"];
    end

    function testAllDiagonal()
            const to = TimerOutput()
            for testmatrix in getTestMatrices()
                    try
                            testOneDiagonal(testmatrix, to)
                    catch
                            println(string("Fehler bei der Matrix ",testmatrix))
                    end
            end
            print(to)
    end

    function testOneDiagonal(testmatrix, to)
            println(string("Teste die Matrix ",testmatrix))

            A = MatrixMarket.mmread(string("matrices/", testmatrix, ".mtx"))
            println(string("---Matrix ",testmatrix, " eingelesen"))
            n = size(A,1)
            b = ones(Float64, n,1)
            useM3 = true

            @timeit to string(testmatrix) x1, history1 = Solver.solve(A, b, SimpleProlongations.prolongation1(n), useM3)
            println(string("---Berechnung fertig"))
            println(string("---", history1))
            history1data = history1.data[:resnorm]

            pyplot()
            plot(history1data, xlabel="iterations", ylabel="res-norm", label=testmatrix, yscale = :log10)

            savefig(string("../auswertungen/diagonalRestriction/", testmatrix, ".png"))
            println(string("---Plot gespeichert"))

    end

    function test()
        A = getA()
        n = size(A,1)
        b = ones(Float64, n,1)
        useM3 = true
        # solve
        @time x1, history1 = Solver.solve(getA(), b, SimpleProlongations.prolongation1(n), useM3)
        # @time x2, history2 = Solver.solve(getA(), b, SimpleProlongations.prolongation2(n), useM3)
        # @time x3, history3 = Solver.solve(getA(), b, AggregationbasedProlongation.prolongationAggregation(getA()), useM3)
        # @time x4, history4 = Solver.solve(getA(), b, ReductionbasedProlongation.prolongationReduction(getA()), useM3)

        # plot
        #gr()
        pyplot()
        history1data = history1.data[:resnorm]
        # history2data = history2.data[:resnorm]
        # history3data = history3.data[:resnorm]
        # history4data = history4.data[:resnorm]

        println(string("res-norm: ", history1data[end]))
        # println(string("res-norm: ", history2data[end]))
        # println(string("res-norm: ", history3data[end]))
        # println(string("res-norm: ", history4data[end]))

        fig = plot(history1data, xlabel="iterations", ylabel="res-norm", label="simple prolongation", yscale = :log10)
        # plot!(history2data, label="prolongation with weights")
        # plot!(history3data, label="aggregation-based prolongation")
        # plot!(history4data, label="reduction-based prolongation")
        #savefig("prolongations.png")
    end
end
