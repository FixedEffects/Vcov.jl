using Test
using Vcov
using Vcov: completecases, materialize, simple, robust, cluster, names, nclusters, invcrossmodelmatrix
using GroupedArrays
using LinearAlgebra
using StatsAPI: vcov
using Random
const tests = [
    "estimators"
]

printstyled("Running tests:\n", color=:blue, bold=true)

for test in tests
    include("$test.jl")
    println("\033[1m\033[32mPASSED\033[0m: $(test)")
end
