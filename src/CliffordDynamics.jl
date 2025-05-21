module CliffordDynamics

using QuantumClifford
using LinearAlgebra
using Distributions

# Write your package code here.
include("Types.jl")
include("GenerateOperators.jl")
include("Lattice_10_3a.jl")
include("Lattice_Honeycomb.jl")
include("Simulate2D.jl")
include("Observables.jl")
include("Simulate3D.jl")


export create_lattice_honeycomb

end
