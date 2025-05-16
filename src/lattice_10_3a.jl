include("Types.jl")
include("GenerateOperators.jl")

### Definition of the (10,3)a Lattice ###
struct Lattice10_3a <: AbstractLattice
    L::Int
    unitcell::Int
    XX::Array{PauliOperator, 4}
    YY::Array{PauliOperator, 4}
    ZZ::Array{PauliOperator, 4}
    W::Array{PauliOperator, 4}
    loop1::PauliOperator
    loop2::PauliOperator
    loop3::PauliOperator
    index_fn::Function
    subsystem_fn::Function
end

####### My definition of the lattice ######
#the (10,3)c lattice has a basis of 6 sites
# L_unitcell = 4

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

## function that creates the Operators of the 10_3a lattice
function create_operators_10_3a(L)
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

    ## constructs the Wilson loops
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

    return XX, YY, ZZ, W, loop1, loop2, loop3

end

function create_lattice_10_3a(L::Int)::Lattice10_3a
    unitcell = 4
    XX, YY, ZZ, W, loop1, loop2, loop3 = create_operators_10_3a(L)

    # Define index function
    index_fn = (x, y, z, n) -> begin
        Lmod = L
        return unitcell * (mod1(x, Lmod) - 1) +
               unitcell * Lmod * (mod1(y, Lmod) - 1) +
               unitcell * Lmod * Lmod * (mod1(z, Lmod) - 1) + n
    end

    # Define subsystem function (placeholder)
    subsystem_fn = (lcols, L; cut=nothing) -> begin
        return get_subsystem_10_3_a(lcols, L; cut=cut)
    end

    return Lattice10_3a(L, unitcell, XX, YY, ZZ, W, loop1, loop2, loop3, index_fn, subsystem_fn)
end

create_lattice_10_3a(4)

xx,yy,zz,w,l1,l2,l3 = create_operators_10_3a(4)
typeof(l1)