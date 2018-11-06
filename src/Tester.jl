module Tester
export testAll
export testOneForRes

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
            return ["bcsstm26", "bcspwr07", "bcspwr09", "bcsstk23", "bcsstk24", "bcsstk25", "dwt__234", "dwt__992", "gemat11", "gr_30_30", "jgl009", "mbeause", "nnc666", "nos7", "saylr3", "watt__1"];
    end

    function testAll()
            converged = Dict()
            history = Dict()

            converged["diagonalRestriction"] = Dict()
            converged["weightedDiagonalRestriction"] = Dict()
            converged["aggregationRestriction"] = Dict()
            converged["reductionRestriction"] = Dict()

            converged["diagonalRestriction"][true] = Dict()
            converged["weightedDiagonalRestriction"][true] = Dict()
            converged["aggregationRestriction"][true] = Dict()
            converged["reductionRestriction"][true] = Dict()

            converged["diagonalRestriction"][false] = Dict()
            converged["weightedDiagonalRestriction"][false] = Dict()
            converged["aggregationRestriction"][false] = Dict()
            converged["reductionRestriction"][false] = Dict()

            history["diagonalRestriction"] = Dict()
            history["weightedDiagonalRestriction"] = Dict()
            history["aggregationRestriction"] = Dict()
            history["reductionRestriction"] = Dict()

            history["diagonalRestriction"][true] = Dict()
            history["weightedDiagonalRestriction"][true] = Dict()
            history["aggregationRestriction"][true] = Dict()
            history["reductionRestriction"][true] = Dict()

            history["diagonalRestriction"][false] = Dict()
            history["weightedDiagonalRestriction"][false] = Dict()
            history["aggregationRestriction"][false] = Dict()
            history["reductionRestriction"][false] = Dict()

            testOneForRes(converged, history, "diagonalRestriction")
            testOneForRes(converged, history, "weightedDiagonalRestriction")
            testOneForRes(converged, history, "aggregationRestriction")
            testOneForRes(converged, history, "reductionRestriction")

            testOneForRes(converged, history, "diagonalRestriction", false)
            testOneForRes(converged, history, "weightedDiagonalRestriction", false)
            testOneForRes(converged, history, "aggregationRestriction", false)
            testOneForRes(converged, history, "reductionRestriction", false)

            printForOneM(history, true)
            printForOneM(history, false)

            printM3PN("diagonalRestriction", history)
            printM3PN("weightedDiagonalRestriction", history)
            printM3PN("aggregationRestriction", history)
            printM3PN("reductionRestriction", history)
    end

    function testOneForRes(converged, history, method, useM3 = true)
            to = TimerOutput()
            toR = TimerOutput()
            if(useM3)
                        M = "M3"
            else
                        M = "PN"
            end
            println("Methode $method mit $M")

            localconverged = converged[method][useM3]
            for testmatrix in getTestMatrices()
                    try
                            testOne(testmatrix, method, to, toR, localconverged, history[method][useM3], useM3)
                    catch
                            println("Fehler bei der Matrix $testmatrix")
                    end
            end

            mkpath("../auswertungen/$method/$M")

            open("../auswertungen/$method/$M/converged.txt", "w") do f
                    conv = "";
                    for k in keys(localconverged)
                        dict = localconverged[k]
                        conv = "$conv $k => $dict \n"
                    end
                    write(f, "Konvergiert?\n $conv")
            end

            open("../auswertungen/$method/$M/timeSolve.txt", "w") do f
                    write(f, "$to")
            end

            open("../auswertungen/$method/$M/timeRestriction.txt", "w") do f
                    write(f, "$toR")
            end
    end

    function testOne(testmatrix, method, to=TimerOutput(), toR=TimerOutput(), converged=Dict(), history=Dict(), useM3 = true)
            println("+++ Teste die Matrix $testmatrix")

            A = MatrixMarket.mmread("matrices/$testmatrix.mtx")
            println("+++ --- Matrix $testmatrix eingelesen")

            n = size(A,1)
            b = ones(Float64, n,1)

            if(method == "diagonalRestriction")
                        @timeit toR string(testmatrix) R = SimpleProlongations.prolongation1(n)
            elseif(method == "weightedDiagonalRestriction")
                        @timeit toR string(testmatrix) R = SimpleProlongations.prolongation2(n)
            elseif(method == "aggregationRestriction")
                        @timeit toR string(testmatrix) R = AggregationbasedProlongation.prolongationAggregation(A)
            else
                        @timeit toR string(testmatrix) R = ReductionbasedProlongation.prolongationReduction(A)
            end

            println("+++ --- Restriktionsmatrix berechnet")
            @timeit to string(testmatrix) x1, history1 = Solver.solve(A, b, R, useM3)
            println("+++ --- Berechnung fertig")
            println("+++ --- $history1")
            converged[testmatrix] = "$history1"
            history[testmatrix] = history1.data[:resnorm]

    end


            function printForOneM(history, useM3)
                   print("diagonalRestriction", history, true, false, false, false, useM3)
                   print("weightedDiagonalRestriction", history, false, true, false, false, useM3)
                   print("aggregationRestriction", history, false, false, true, false, useM3)
                   print("reductionRestriction", history, false, false, false, true, useM3)

                   print("simpleRestrictions", history, true, true, false, false, useM3)
                   print("complexRestrictions", history, false, false, true, true, useM3)

                   print("allRestrictions", history, true, true, true, true, useM3)
            end

    function print(dir, history, useSimple1=false, useSimple2=false, useAggreg=false, useReduc=false, useM3=true)
                if(useM3)
                           M = "M3"
                else
                           M = "PN"
                end
                println("Plot für $dir mit $M")
                for testmatrix in getTestMatrices()
                           println("+++ Plot mit $testmatrix")
                           pyplot()
                           try
                                       if useSimple1
                                                   plot(history["diagonalRestriction"][useM3][testmatrix], xlabel="iterations", ylabel="res-norm", label="diagonalRestriction", yscale = :log10)
                                                   if(useSimple2)
                                                               plot!(history["weightedDiagonalRestriction"][useM3][testmatrix], label="weightedDiagonalRestriction")
                                                   end
                                                   if(useAggreg)
                                                               plot!(history["aggregationRestriction"][useM3][testmatrix], label="aggregationRestriction")
                                                   end
                                                   if(useReduc)
                                                               plot!(history["reductionRestriction"][useM3][testmatrix], label="reductionRestriction")
                                                   end
                                       elseif useSimple2
                                                   plot(history["weightedDiagonalRestriction"][useM3][testmatrix], xlabel="iterations", ylabel="res-norm", label="weightedDiagonalRestriction", yscale = :log10)
                                       elseif useAggreg
                                                   plot(history["aggregationRestriction"][useM3][testmatrix], xlabel="iterations", ylabel="res-norm", label="aggregationRestriction", yscale = :log10)
                                                   if(useReduc)
                                                               plot!(history["reductionRestriction"][useM3][testmatrix], label="reductionRestriction")
                                                   end
                                       elseif useReduc
                                                   plot(history["reductionRestriction"][useM3][testmatrix], xlabel="iterations", ylabel="res-norm", label="reductionRestriction", yscale = :log10)
                                       end

                                       mkpath("../auswertungen/$dir/$M")
                                       savefig("../auswertungen/$dir/$M/$testmatrix.png")

                           catch
                                     println("+++ --- Fehler")
                           end
                end

    end


    function printM3PN(dir, history)
                println("Plot für $dir mit beiden Vorkondtionierern")
                for testmatrix in getTestMatrices()
                           println("+++ Plot mit $testmatrix")
                           pyplot()
                           try
                                       if dir == "diagonalRestriction"
                                           plot(history["diagonalRestriction"][true][testmatrix], xlabel="iterations", ylabel="res-norm", label="M3", yscale = :log10)
                                           plot!(history["diagonalRestriction"][false][testmatrix], label="PN")
                                       elseif dir == "weightedDiagonalRestriction"
                                           plot(history["weightedDiagonalRestriction"][true][testmatrix], xlabel="iterations", ylabel="res-norm", label="M3", yscale = :log10)
                                           plot!(history["weightedDiagonalRestriction"][false][testmatrix], label="PN")
                                       elseif dir == "aggregationRestriction"
                                           plot(history["aggregationRestriction"][true][testmatrix], xlabel="iterations", ylabel="res-norm", label="M3", yscale = :log10)
                                           plot!(history["aggregationRestriction"][false][testmatrix], label="PN")
                                       else
                                           plot(history["reductionRestriction"][true][testmatrix], xlabel="iterations", ylabel="res-norm", label="M3", yscale = :log10)
                                           plot!(history["reductionRestriction"][false][testmatrix], label="PN")
                                       end

                                       mkpath("../auswertungen/$dir/both")
                                       savefig("../auswertungen/$dir/both/$testmatrix.png")

                           catch
                                     println("+++ --- Fehler")
                           end
                end

    end

end
