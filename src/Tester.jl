using DataFrames
using Plots
using MatrixMarket

include("SimpleProlongations.jl")
include("AggregationbasedProlongation.jl")
include("ReductionbasedProlongation.jl")
include("Solver.jl")

function getA()
    # mmread("matrices/bcsstm26.mtx")
    # mmread("matrices/nos7.mtx")
    # mmread("matrices/pde225.mtx")
    # mmread("matrices/gr_30_30.mtx")
    # MatrixMarket.mmread("matrices/bcsstk25.mtx")
    MatrixMarket.mmread("matrices/fidap005.mtx")
end

A = getA()
n = size(A,1)
b = ones(Float64, n,1)
useM3 = false
# solve
@time x1, history1 = Solver.solve(getA(), b, SimpleProlongations.prolongation1(n), useM3)
@time x2, history2 = Solver.solve(getA(), b, SimpleProlongations.prolongation2(n), useM3)
@time x3, history3 = Solver.solve(getA(), b, AggregationbasedProlongation.prolongationAggregation(getA()), useM3)
@time x4, history4 = Solver.solve(getA(), b, ReductionbasedProlongation.prolongationReduction(getA()), useM3)

# plot
#gr()
pyplot()
history1data = history1.data[:resnorm]
history2data = history2.data[:resnorm]
history3data = history3.data[:resnorm]
history4data = history4.data[:resnorm]

println(string("res-norm: ", history1data[end]))
println(string("res-norm: ", history2data[end]))
println(string("res-norm: ", history3data[end]))
println(string("res-norm: ", history4data[end]))

fig = plot(history1data, xlabel="iterations", ylabel="res-norm", label="simple prolongation", yscale = :log10)
plot!(history2data, label="prolongation with weights")
plot!(history3data, label="aggregation-based prolongation")
plot!(history4data, label="reduction-based prolongation")
#savefig("prolongations.png")
