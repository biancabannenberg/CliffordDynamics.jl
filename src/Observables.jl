# include("Lattice_10_3a.jl")
# include("Lattice_Honeycomb.jl")
# include("Simulate2D.jl")
# using Distributions
# using DataFrames

### Definition of the Entanglement Entropy and Tripartite mutual information###
# lat = create_lattice_10_3a(4)
#function that calculates entanglement entropy of subsystem of length l
#firstly one need to define the subsystem of length l and total systemsize L (ss_10_3_c(l,L))
#state: MixedDestabilizer
function calc_EE(state, l, lattice::AbstractLattice; cut=:z)
    return entanglement_entropy(state, lattice.subsystem(l, cut=cut), Val(:rref))
end
# st= initialize(lat)
# calc_EE(st, 1, lat, cut = :a1)

#function that calculates entanglement for given subdivistion lcols
function S_A(state, lcols, lattice; cut=:z)
    return entanglement_entropy(state, lattice.subsystem(lcols, cut=cut), Val(:rref))
end

# S_A(st, 2:4, lat, cut = :a1)

function entanglement_arc(state, lattice; cut = :a1)
    S_arc = Array{Float64}(undef, lattice.L+1)

    for (i, l) in enumerate(0:lattice.L)
        S_arc[i] = S_A(state, 1:l, lattice, cut = cut )
    end

    return S_arc
end

# entanglement_arc(st, lat, cut = :a1)
#function that imoplements the tripartite mutual information
#in this definition: system partet in 4 subsystems going from 1:L/4, L/4+1:L/2, L/2+1:3L/4, 3L/4:L
# then: I= SA + SB + SC - SAB - SAC - SBA + SAB
function TMI(state, lattice; cut=:z)
    L = lattice.L
    if L % 4 != 0
        error("L must be a multiple of 4")
    end
    SA = S_A(state, 1:L÷4, lattice, cut=cut)
    SB = S_A(state, L÷4+1:L÷2, lattice, cut=cut)
    SC = S_A(state, L÷2+1:3L÷4, lattice, cut=cut)
    SAB = S_A(state, 1:L÷2, lattice, cut=cut)
    SBC = S_A(state, L÷4+1:3L÷4, lattice, cut=cut)
    SAC = S_A(state, vcat(1:L÷4, L÷2+1:3L÷4), lattice, cut=cut)
    SABC = S_A(state, 1:3L÷4, lattice, cut=cut)
    return SA + SB + SC - SAB - SBC - SAC + SABC
end
# TMI(st, lat, cut = :a1)
create_lattice_honeycomb(16)