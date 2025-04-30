using QuantumClifford
using Distributions
using DataFrames

####### nearest neighbors ######
#the (10,3)c lattice has a basis of 6 sites
L_unitcell = 4

### from LatticePhysics.jl the prepared unit cell and connections
# the lattice vectors
# a1 = [1, 0, 0]
# a2 = [0.5, 0.5, -0.5]
# a3 = [0.5, 0.5, 0.5]
# # Basis Definition
# basis = Array{Float64, 1}[
#     [1/8., 1/8.,  1/8.],
#     [5/8., 3/8., -1/8.],
#     [3/8., 1/8., -1/8.],
#     [7/8., 3/8.,  1/8.]
# ]
# connections of sites in unitcell:
#! [1; 3] XX
#! [2; 4] XX

#connections ourtside of unitcell
#? [4; 1; (1, 0, 0)] ZZ
#? [1; 4; (-1, 0, 0)] ZZ
#? [3; 2; (0, 1, -1)] ZZ
#? [2; 3; (0, -1, 1)] ZZ

#* [2; 1; (0, 1, 0)] YY
#* [1; 2; (0, -1, 0)] YY
#* [3; 4; (0, 0, -1)] YY
#* [4; 3; (0, 0, 1)] YY

#at first define lattice structure by defining the neighbors of a lattice site i in an L-dim system 

#function that returns index of qubit at (x,y,z)
#n : n-th qubit in basis
function index_xyz(x::Int, y::Int, z::Int, n::Int, L::Int)
    return L_unitcell * (mod1(x, L) - 1) + L_unitcell * L * (mod1(y, L) - 1) + L_unitcell * L * L * (mod1(z, L) - 1) + n
end

#function that returns subsystem of a 10_3_c lattice
#general structure of subsystem functions:
#    L_unitcell * (x-1) + L_unitcell * L * (y-1) + L_unitcell * L^2 * (y-1)   --> iterate through all x,y and z with x in l+1:L, and y,z in 1:L
function ss_10_3_a(lcols::Int, L::Int; cut=:y)
    N = L^3 * L_unitcell
    if cut == :a1
        return ([L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l for j in lcols+1:L for k in 1:L for m in 1:L for l in 1:L_unitcell])
    elseif cut == :a2
        return ([L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l for j in 1:L for k in lcols+1:L for m in 1:L for l in 1:L_unitcell])
    elseif cut == :a3
        return ([L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l for j in 1:L for k in 1:L for m in lcols+1:L for l in 1:L_unitcell])
    elseif cut == :x1
        if lcols == 1:L
            return ([1:N...])
        else
            return ([[(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                      for j in 1:L           # iterate over the x-direction
                      for k in 1:L           # iterate over the y-direction
                      for m in lcols+1     # iterate over the z-direction
                      for l in [3; 2; 4]  # select atoms 2 through 4 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in lcols+2:L  # iterate over the z-direction
                 for l in 1:L_unitcell  # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in 1:1      # iterate over the z-direction
                 for l in 1:1            # select atoms 2 through 6 of the unit cell
                ]...])
        end
    elseif cut == :x2
        if lcols == 1:L
            return ([1:N...])
        else
            return ([[(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                      for j in 1:L           # iterate over the x-direction
                      for k in 1:L           # iterate over the y-direction
                      for m in lcols+1      # iterate over the z-direction
                      for l in 4:4             # select atoms 4 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in lcols+2:L  # iterate over the z-direction
                 for l in 1:L_unitcell  # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in 1:1      # iterate over the z-direction
                 for l in [1; 3; 2]            # select atoms 2 through 6 of the unit cell
                ]...])
        end
    elseif cut == :xz
        if lcols == 1:L
            return ([1:N...])
        else
            return ([[(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                      for j in 1:L           # iterate over the x-direction
                      for k in 1:L           # iterate over the y-direction
                      for m in lcols+1   # iterate over the z-direction
                      for l in 2:L_unitcell           # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in lcols+2:L    # iterate over the z-direction
                 for l in 1:L_unitcell           # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in 1:1      # iterate over the z-direction
                 for l in 1:1            # select atoms 2 through 6 of the unit cell
                ]...])
        end
    else
        error("cut must be :a1, :a2, :a3, :x1, :x2 or :xz")
    end
end

function ss_10_3_a(lcols, L::Int; cut=:y)
    N = L^3 * L_unitcell
    if cut == :a1 ## cuts through blue bond == ZZ
        return ([L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l for j in lcols for k in 1:L for m in 1:L for l in 1:L_unitcell])
    elseif cut == :a2 ## cuts through green+blue bond == YY+ZZ
        return ([L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l for j in 1:L for k in lcols for m in 1:L for l in 1:L_unitcell])
    elseif cut == :a3 ## cuts through green bond == YY
        return ([L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l for j in 1:L for k in 1:L for m in lcols for l in 1:L_unitcell])
    elseif cut == :x1
        if lcols == 1:L
            return ([1:N...])
        else
            return ([[(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                      for j in 1:L           # iterate over the x-direction
                      for k in 1:L           # iterate over the y-direction
                      for m in lcols[1]      # iterate over the z-direction
                      for l in [3; 2; 4]  # select atoms 2 through 4 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in lcols[2:end]  # iterate over the z-direction
                 for l in 1:L_unitcell  # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in mod1(lcols[end] + 1, L)      # iterate over the z-direction
                 for l in 1:1            # select atoms 2 through 6 of the unit cell
                ]...])
        end
    elseif cut == :x2
        if lcols == 1:L
            return ([1:N...])
        else
            return ([[(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                      for j in 1:L           # iterate over the x-direction
                      for k in 1:L           # iterate over the y-direction
                      for m in lcols[1]      # iterate over the z-direction
                      for l in 4:4             # select atoms 4 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in lcols[2:end]  # iterate over the z-direction
                 for l in 1:L_unitcell  # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in mod1(lcols[end] + 1, L)      # iterate over the z-direction
                 for l in [1; 3; 2]            # select atoms 2 through 6 of the unit cell
                ]...])
        end
    elseif cut == :xz
        if lcols == 1:L
            return ([1:N...])
        else
            return ([[(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                      for j in 1:L           # iterate over the x-direction
                      for k in 1:L           # iterate over the y-direction
                      for m in lcols[1]   # iterate over the z-direction
                      for l in 2:L_unitcell           # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in lcols[2:end]    # iterate over the z-direction
                 for l in 1:L_unitcell           # select atoms 2 through 6 of the unit cell
                ]...
                [(L_unitcell * (j - 1) + L_unitcell * L * (k - 1) + L_unitcell * L^2 * (m - 1) + l)
                 for j in 1:L           # iterate over the x-direction
                 for k in 1:L           # iterate over the y-direction
                 for m in mod1(lcols[end] + 1, L)      # iterate over the z-direction
                 for l in 1            # select atoms 2 through 6 of the unit cell
                ]...])
        end
    else
        error("cut must be :a1, :a2, :a3, :x1, :x2 or :xz")
    end
end
# println(ss_10_3_c(3, 4))

#function that calculates entanglement entropy of subsystem of length l
#firstly one need to define the subsystem of length l and total systemsize L (ss_10_3_c(l,L))
#state: MixedDestabilizer
function calc_EE(state, l, L; cut=:z)
    return entanglement_entropy(state, ss_10_3_a(l, L, cut=cut), Val(:rref))
end

#function that calculates entanglement for given subdivistion lcols
function S_A(state, lcols, L; cut=:z)
    return entanglement_entropy(state, ss_10_3_a(lcols, L, cut=cut), Val(:rref))
end

#function that imoplements the tripartite mutual information
#in this definition: system partet in 4 subsystems going from 1:L/4, L/4+1:L/2, L/2+1:3L/4, 3L/4:L
# then: I= SA + SB + SC - SAB - SAC - SBA + SAB
function TMI(state, L; cut=:z)
    if L % 4 != 0
        error("L must be a multiple of 4")
    end
    SA = S_A(state, 1:L÷4, L, cut=cut)
    SB = S_A(state, L÷4+1:L÷2, L, cut=cut)
    SC = S_A(state, L÷2+1:3L÷4, L, cut=cut)
    SAB = S_A(state, 1:L÷2, L, cut=cut)
    SBC = S_A(state, L÷4+1:3L÷4, L, cut=cut)
    SAC = S_A(state, vcat(1:L÷4, L÷2+1:3L÷4), L, cut=cut)
    SABC = S_A(state, 1:3L÷4, L, cut=cut)
    return SA + SB + SC - SAB - SBC - SAC + SABC
end


# functions to generate XX, YY and ZZ operators on two given qubits i,j
function genXX(N, i::Int, j::Int)
    Xarr = falses(N)
    Zarr = falses(N)
    Xarr[i] = true
    Xarr[j] = true
    return PauliOperator(0x00, Xarr, Zarr)
end

function genYY(N, i, j)
    Xarr = falses(N)
    Zarr = falses(N)
    Zarr[i] = true
    Zarr[j] = true
    Xarr[i] = true
    Xarr[j] = true
    return PauliOperator(0x00, Xarr, Zarr)
end

function genZZ(N, i, j)
    Xarr = falses(N)
    Zarr = falses(N)
    Zarr[i] = true
    Zarr[j] = true
    return PauliOperator(0x00, Xarr, Zarr)
end

#generates loop-operator of the 10_3_c lattice for given 
#ps::operators ("X","Y","Z")
#is::indizes of operators part of the loop
function genLoop(N, ps, is)
    Xarr = falses(N)
    Zarr = falses(N)
    for (p, i) in zip(ps, is)
        if p == "X"
            Xarr[i] = true
        elseif p == "Z"
            Zarr[i] = true
        elseif p == "Y"
            Xarr[i] = true
            Zarr[i] = true
        else
            error("Invalid Pauli operator")
        end
    end
    return PauliOperator(0x00, Xarr, Zarr)
end


function create_operators(L)
    #variables needed consitently
    #N: number of qubits
    #L: sidelength of the lattice

    N = L^3 * L_unitcell


    # construct, ZZ, XX, YY, W, loop1, loop2
    XX = Array{PauliOperator}(undef, L, L, L, 2)
    for xi in 1:L
        for yi in 1:L
            for zi in 1:L
                XX[xi, yi, zi, 1] = genXX(N, index_xyz(xi, yi, zi, 1, L), index_xyz(xi, yi, zi, 3, L))
                XX[xi, yi, zi, 2] = genXX(N, index_xyz(xi, yi, zi, 2, L), index_xyz(xi, yi, zi, 4, L))
            end
        end
    end

    YY = Array{PauliOperator}(undef, L, L, L, 2)
    for xi in 1:L
        for yi in 1:L
            for zi in 1:L
                YY[xi, yi, zi, 1] = genYY(N, index_xyz(xi, yi, zi, 2, L), index_xyz(xi, yi + 1, zi, 1, L))
                YY[xi, yi, zi, 2] = genYY(N, index_xyz(xi, yi, zi, 3, L), index_xyz(xi, yi, zi - 1, 4, L))
            end
        end
    end

    ZZ = Array{PauliOperator}(undef, L, L, L, 2)
    for xi in 1:L
        for yi in 1:L
            for zi in 1:L
                ZZ[xi, yi, zi, 1] = genZZ(N, index_xyz(xi, yi, zi, 4, L), index_xyz(xi + 1, yi, zi, 1, L))
                ZZ[xi, yi, zi, 2] = genZZ(N, index_xyz(xi, yi, zi, 3, L), index_xyz(xi, yi + 1, zi - 1, 2, L))
            end
        end
    end


    W = Array{PauliOperator}(undef, L, L, L, 4)
    for xi in 1:L
        for yi in 1:L
            for zi in 1:L
                #at first define the 10 operator plaquettes in the order defined in the pictures provided
                #plaquette 1
                is = []
                ps = []
                push!(is, index_xyz(xi, yi, zi, 1, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi, yi, zi, 3, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi, yi + 1, zi - 1, 2, L))
                push!(ps, "X")
                push!(is, index_xyz(xi, yi + 2, zi - 1, 1, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi + 2, zi - 1, 4, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi - 1, yi + 2, zi - 1, 2, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi - 1, yi + 1, zi, 3, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi - 1, yi + 1, zi, 1, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi - 1, yi, zi, 2, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi - 1, yi, zi, 4, L))
                push!(ps, "Y")
                W[xi, yi, zi, 1] = genLoop(N, ps, is)


                #plaquette 2
                is = []
                ps = []
                push!(is, index_xyz(xi, yi, zi, 1, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi, zi, 3, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi, zi - 1, 4, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi, zi - 1, 2, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi, yi - 1, zi, 3, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi, yi - 1, zi, 1, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi - 2, zi, 2, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi - 2, zi, 4, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi - 2, zi + 1, 3, L))
                push!(ps, "X")
                push!(is, index_xyz(xi, yi - 1, zi, 2, L))
                push!(ps, "X")
                W[xi, yi, zi, 2] = genLoop(N, ps, is)


                #plaquette 3
                is = []
                ps = []
                push!(is, index_xyz(xi, yi, zi, 1, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi, yi, zi, 3, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi, yi + 1, zi - 1, 2, L))
                push!(ps, "X")
                push!(is, index_xyz(xi, yi + 2, zi - 1, 1, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi + 2, zi - 1, 4, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi + 2, zi, 3, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi - 1, yi + 2, zi, 1, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi - 1, yi + 1, zi, 2, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi, zi + 1, 3, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi, zi, 4, L))
                push!(ps, "X")
                W[xi, yi, zi, 3] = genLoop(N, ps, is)


                #plaquette 4
                is = []
                ps = []
                push!(is, index_xyz(xi, yi, zi, 2, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi, zi, 4, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi, zi + 1, 3, L))
                push!(ps, "X")
                push!(is, index_xyz(xi, yi + 1, zi, 2, L))
                push!(ps, "X")
                push!(is, index_xyz(xi, yi + 2, zi, 1, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi + 2, zi, 4, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi - 1, yi + 2, zi, 2, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi - 1, yi + 1, zi + 1, 3, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi + 1, zi, 4, L))
                push!(ps, "X")
                push!(is, index_xyz(xi, yi + 1, zi, 1, L))
                push!(ps, "X")
                W[xi, yi, zi, 4] = genLoop(N, ps, is)
            end
        end
    end

    is = []
    ps = [] ### not redundant
    for xi in 1:L
        push!(is, index_xyz(xi, 1, 1, 1, L))
        push!(ps, "Y")
        push!(is, index_xyz(xi, 1, 1, 3, L))
        push!(ps, "Y")
        push!(is, index_xyz(xi, 2, 0, 2, L))
        push!(ps, "Y")
        push!(is, index_xyz(xi, 2, 0, 4, L))
        push!(ps, "Z")
        push!(is, index_xyz(xi, 2, 1, 3, L))
        push!(ps, "Z")
        push!(is, index_xyz(xi, 2, 1, 1, L))
        push!(ps, "Z")
        push!(is, index_xyz(xi, 1, 1, 2, L))
        push!(ps, "Z")
        push!(is, index_xyz(xi, 1, 1, 4, L))
        push!(ps, "Y")
    end
    loop1 = genLoop(N, ps, is)

    is = []
    ps = []
    for yi in 1:L
        push!(is, index_xyz(1, yi, 1, 3, L))
        push!(ps, "X")
        push!(is, index_xyz(1, yi + 1, 0, 2, L))
        push!(ps, "Y")
        push!(is, index_xyz(1, yi + 1, 0, 4, L))
        push!(ps, "Z")
    end
    loop2 = genLoop(N, ps, is)

    is = []
    ps = []
    for zi in 1:L
        push!(is, index_xyz(1, 1, zi, 4, L))
        push!(ps, "Z")
        push!(is, index_xyz(1, 1, zi + 1, 3, L))
        push!(ps, "X")
        push!(is, index_xyz(1, 2, zi, 2, L))
        push!(ps, "Y")
        push!(is, index_xyz(1, 2, zi, 4, L))
        push!(ps, "Z")
        push!(is, index_xyz(1, 2, zi + 1, 3, L))
        push!(ps, "Z")
        push!(is, index_xyz(1, 2, zi + 1, 1, L))
        push!(ps, "Z")
        push!(is, index_xyz(1, 1, zi + 1, 2, L))
        push!(ps, "Z")
    end
    loop3 = genLoop(N, ps, is)

    # loop1,loop2, loop3, W,
    return XX, YY, ZZ, W, loop1, loop2, loop3# [XX..., YY..., ZZ..., W..., loop1, loop2, loop3]

end

# # check if all operators commute, if they dont, rework the operators!!!
# x, y, z, o, l1, l2, l3 = create_operators(2)
# for (i, o1) in enumerate(y)
#     for (j, o2) in enumerate(o)
#         if comm(o1, o2) != 0x00
#             println("ANTICOMMUTE: ", i, " o2: ", j)
#         end
#     end
#     if comm(o1, l1) != 0x00
#         println("ANTICOMMUTE: ", i, " o2: l1")
#     end
#     if comm(o1, l2) != 0x00
#         println("ANTICOMMUTE: ", i, " o2: l2")
#     end
#     if comm(o1, l3) != 0x00
#         println("ANTICOMMUTE: ", i, " o2: l3")
#     end
# end


#initializes state of Length L with complete tableau of stabilizers
function initialize(L::Int; keep_result=false, phases=false)
    #L length of system
    #N number of qubits
    local N = L^3 * L_unitcell

    XX, YY, ZZ, o, l1, l2, l3 = create_operators(L)

    state = MixedDestabilizer(zero(Stabilizer, N))

    # for x in eachindex(XX)
    #     project!(state, XX[x], keep_result=keep_result, phases=phases)
    # end
    # for y in eachindex(YY)
    #     project!(state, YY[y], keep_result=keep_result, phases=phases)
    # end
    for z in eachindex(ZZ)
        project!(state, ZZ[z], keep_result=keep_result, phases=phases)
    end

    for i in eachindex(o)
        project!(state, o[i], keep_result=keep_result, phases=phases)
    end

    project!(state, l1, keep_result=keep_result, phases=phases)
    project!(state, l2, keep_result=keep_result, phases=phases)
    project!(state, l3, keep_result=keep_result, phases=phases) #this loop is not needed for initialization if all x,y,z are measured

    if QuantumClifford.trusted_rank(state) != N
        println("Tableau does not have sufficient Rank: $(QuantumClifford.trusted_rank(state)) / $(N)")
        # return state
    else
        return state
    end
end

# initialize(4)


# #initialize function with already computed operators o
# function initialize(L::Int,o::Array; keep_result = false, phases = false)
#     #L length of system
#     #N number of qubits
#     local N = L^3 * L_unitcell

#     state = MixedDestabilizer(zero(Stabilizer,N))

#     for i in eachindex(o)
#         project!(state, o[i], keep_result=keep_result, phases=phases)
#     end

#     if QuantumClifford.trusted_rank(state) !== N
#         println("Tableau does not have sufficient Rank")
#         return nothing
#     end

#     return state
# end


# function simulate(L::Int, p; cut=:y, iterations::Int=25, thermalization_steps=20, keep_result=false, phases=false)
#     N = L^3 * L_unitcell

#     #calculate operators
#     XX, YY, ZZ = create_operators(L)

#     #probability distribution used : p=[p(XX), p(YY), p(ZZ)]
#     dist = Distributions.Categorical([p[1], p[2], p[3]])

#     #initialize sum of entanglement entropy, mutual information and tripartite mutual information 
#     #:a1, :a2, :a3, :x1, :x2 or :xz
#     if cut == :all
#         EE_a1 = 0
#         EE_a2 = 0
#         EE_a3 = 0
#         EE_xz = 0
#         EE_x1 = 0
#         EE_x2 = 0
#         Tmi_a1 = 0
#         Tmi_a2 = 0
#         Tmi_a3 = 0
#         Tmi_xz = 0
#         Tmi_x1 = 0
#         Tmi_x2 = 0
#     else
#         EE = 0
#         #MI = 0
#         Tmi = 0
#     end


#     for _ in 1:iterations
#         #first initialize state
#         state = initialize(L)

#         #thermalization steps
#         for _ in 1:thermalization_steps
#             #iterate randomly through all qubits
#             for _ in 1:N
#                 xi = rand(1:L)
#                 yi = rand(1:L)
#                 zi = rand(1:L)
#                 si = rand(1:Int(L_unitcell / 2))
#                 ope = rand(dist)
#                 if ope == 1 #do X operation
#                     project!(state, XX[xi, yi, zi, si], keep_result=keep_result, phases=phases)
#                 elseif ope == 2 #do Y operation
#                     project!(state, YY[xi, yi, zi, si], keep_result=keep_result, phases=phases)
#                 else #do Z operation
#                     project!(state, ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
#                 end
#             end
#         end

#         if cut == :all
#             EE_a1 += calc_EE(state, Int(L // 2), L, cut=:a1)
#             EE_a2 += calc_EE(state, Int(L // 2), L, cut=:a2)
#             EE_a3 += calc_EE(state, Int(L // 2), L, cut=:a3)
#             EE_x1 += calc_EE(state, Int(L // 2), L, cut=:x1)
#             EE_x2 += calc_EE(state, Int(L // 2), L, cut=:x2)
#             EE_xz += calc_EE(state, Int(L // 2), L, cut=:xz)
#             Tmi_a1 += TMI(state, L, cut=:a1)
#             Tmi_a2 += TMI(state, L, cut=:a2)
#             Tmi_a3 += TMI(state, L, cut=:a3)   
#             Tmi_x1 += TMI(state, L, cut=:x1)
#             Tmi_x2 += TMI(state, L, cut=:x2)
#             Tmi_xz += TMI(state, L, cut=:xz)
#         else
#             EE += calc_EE(state, Int(L // 2), L, cut=cut)
#             #MI += I_i_j(state, 1, L//2, L//2+1, L)
#             Tmi += TMI(state, L, cut=cut)
#         end


#         GC.gc()
#     end

#     if cut == :all
#         return EE_a1 / iterations, EE_a2 / iterations, EE_a3 / iterations, EE_x1 / iterations, EE_x2 / iterations, EE_xz / iterations,
#                Tmi_a1 / iterations, Tmi_a2 / iterations, Tmi_a3 / iterations, Tmi_x1 / iterations, Tmi_x2 / iterations, Tmi_xz / iterations
#     else
#         return EE / iterations, Tmi / iterations
#     end
# end

function simulate(L::Int, p; cut=[:y], iterations::Int=25, thermalization_steps=20, keep_result=false, phases=false)
    N = L^3 * L_unitcell

    #calculate operators
    XX, YY, ZZ = create_operators(L)

    #probability distribution used : p=[p(XX), p(YY), p(ZZ)]
    dist = Distributions.Categorical([p[1], p[2], p[3]])

    #initialize sum of entanglement entropy, mutual information and tripartite mutual information 
    #:a1, :a2, :a3, :x1, :x2 or :xz
    data = Dict(vcat(["EE_" .* String.(cut) .=> 0.]..., ["TMI_" .* String.(cut) .=> 0.]...))

    for _ in 1:iterations
        #first initialize state
        state = initialize(L)

        #thermalization steps
        for _ in 1:thermalization_steps
            #iterate randomly through all qubits
            for _ in 1:N
                xi = rand(1:L)
                yi = rand(1:L)
                zi = rand(1:L)
                si = rand(1:Int(L_unitcell / 2))
                ope = rand(dist)
                if ope == 1 #do X operation
                    project!(state, XX[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                elseif ope == 2 #do Y operation
                    project!(state, YY[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                else #do Z operation
                    project!(state, ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                end
            end
        end

        for (i, c) in enumerate(cut)
            data["EE_" * String(c)] += calc_EE(state, Int(L // 2), L, cut= c)/iterations
            data["TMI_" * String(c)] += TMI(state, L, cut=c)/iterations
        end


        GC.gc()
    end

    return (data)
end
# using BenchmarkTools
# @benchmark simulate(8, [1//3, 1//3, 1//3]; cut=[:a1, :a2], iterations=1, thermalization_steps=30)

# a = [:a, :x, :y, :z]
# b = Dict(vcat(["EE_" .* String.(a) .=> 0]..., ["TMI_" .* String.(a) .=> 0]...))
# b = Dict(
#     [("EE_" * string(k) => 0) for k in keys(a)]...,
#     [("TMI_" * string(k) => 0) for k in keys(a)]...
# )
# keys(a)
# ("EE_" .*String.(keys(a)))

# zeros(10,10)

"""
   Function that performs the Monte Carlo simulation for the Compass model
    
    parameters:
    p: array of probabilities for X, Y and Z operations
    L: system size
    keep_result: keep the result of the measurement, default is false
    phases: keep the phases of the measurement, default is false
    iterations: number of iterations, default is 2000
    thermalization: number of thermalization steps, default is 500
    between_measurements: number of steps between measurements, default is 200
    newstarts: number of new starts, default is 1

    returns:
    EE: array of entanglement entropy of the state
    Tmi: array of tripartite mutual information of the state
"""
function main(p, L; cut=:y, keep_result=false, phases=false, iterations=20, thermalization=20, between_measurements=0, newstarts=1)
    XX, YY, ZZ = create_operators(L)
    N = L^3 * L_unitcell
    dist = Distributions.Categorical([p[1], p[2], p[3]])
    EE = []
    Tmi = []

    for _ in 1:newstarts
        state = initialize(L)

        for _ in 1:thermalization
            for _ in 1:N
                xi = rand(1:L)
                yi = rand(1:L)
                zi = rand(1:L)
                si = rand(1:Int(L_unitcell / 2))
                ope = rand(dist)
                if ope == 1 #do X operation
                    project!(state, XX[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                elseif ope == 2 #do Y operation
                    project!(state, YY[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                else #do Z operation
                    project!(state, ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                end
            end
        end

        for _ in 1:iterations
            for _ in 1:N
                xi = rand(1:L)
                yi = rand(1:L)
                zi = rand(1:L)
                si = rand(1:Int(L_unitcell / 2))
                ope = rand(dist)
                if ope == 1 #do X operation
                    project!(state, XX[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                elseif ope == 2 #do Y operation
                    project!(state, YY[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                else #do Z operation
                    project!(state, ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                end
            end
            push!(EE, calc_EE(state, Int(L // 2), L, cut=cut))
            #MI += I_i_j(state, 1, L//2, L//2+1, L)
            push!(Tmi, TMI(state, L, cut=cut))

            # push!(EE, calc_EE(state, L//2, L))
            # push!(Tmi, TMI(state, L))

            for _ in 1:between_measurements
                for _ in 1:N
                    xi = rand(1:L)
                    yi = rand(1:L)
                    zi = rand(1:L)
                    si = rand(1:Int(L_unitcell / 2))
                    ope = rand(dist)
                    if ope == 1 #do X operation
                        project!(state, XX[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                    elseif ope == 2 #do Y operation
                        project!(state, YY[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                    else #do Z operation
                        project!(state, ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                    end
                end
            end
            GC.gc()
        end
    end
    #println(QuantumClifford.trusted_rank(state))
    return EE, Tmi
end

# XX, YY, ZZ, o = create_operators(4)
# state = initialize(4, o)
# project!(state, ZZ[3,2,3,2])
# using BenchmarkTools
# @benchmark simulate(4, [1//3, 1//3, 1//3])
# calc_EE(state,2,4)

# Int(4//2)

# function main(args)
#     id = parse(Int, args[1])
#     label = args[2]
#     N = 30
#     pxs = []
#     pys = []
#     pzs = []
#     for xi in 0:N
#         for yi in xi:N-xi
#         zi = N - xi - yi
#         push!(pxs, xi//N)
#         push!(pys, yi//N)
#         push!(pzs, zi//N)
#         end
#     end
#     px, py, pz = collect(zip(pxs, pys, pzs))[id]
#     L = 8
#     EE, Itri = simulate(L, [px, py, pz])
#     outfile = open("./10_3_c/data_$label.dat", "w")
#     println(outfile, "$px $py $pz $EE $Itwo $Itri")
#     close(outfile)
# end