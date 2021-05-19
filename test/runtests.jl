using Test
using Vcov: GroupedArray, group, factorize!,
    completecases, materialize, simple, robust, cluster, names, nclusters

import Base: ==

const tests = [
    "GroupedArray",
    "estimators"
]

==(x::GroupedArray{N}, y::GroupedArray{N}) where N = x.refs == y.refs && x.n == y.n

printstyled("Running tests:\n", color=:blue, bold=true)

for test in tests
    include("$test.jl")
    println("\033[1m\033[32mPASSED\033[0m: $(test)")
end
