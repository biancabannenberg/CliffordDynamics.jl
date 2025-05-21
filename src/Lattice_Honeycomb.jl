using QuantumClifford
include("Types.jl")
include("GenerateOperators.jl")

### Definition of the (10,3)a Lattice ###
struct LatticeHoneycomb <: AbstractLattice
    L::Int
    unitcell::Int
    N::Int
    XX::Array{PauliOperator,2}
    YY::Array{PauliOperator,2}
    ZZ::Array{PauliOperator,2}
    W::Array{PauliOperator,2}
    loop1::PauliOperator
    loop2::PauliOperator
    index_xyz::Function
    subsystem::Function
end

####### My definition of the lattice ######
#the honeycomb lattice has a basis of 2 sites
# unitcell = 2

### from LatticePhysics.jl the prepared unit cell and connections
# the lattice vectors
# todo: needs to be added similar to the 10_3a lattice

# connections of sites in unitcell:
#? [1;2] ZZ
#? [2;1] ZZ

#connections ourtside of unitcell
#! [2; 1; ( 1, 0)] XX
#! [1; 2; (-1, 0)] XX

#* [2; 1; (0, 1)] YY
#* [1; 2; (0,-1)] YY


function create_lattice_honeycomb(L::Int)::LatticeHoneycomb
    unitcell = 2
    N = L^2 * unitcell

    # Define index function
    index_xyz = (x, y, n, Lmod) -> begin
        return unitcell * (mod1(x, Lmod) - 1) +
               unitcell * Lmod * (mod1(y, Lmod) - 1) + n
    end

    # function that creates the Operators of the 10_3a lattice
    create_operators_honeycomb = (L) -> begin
        #variables needed consitently
        #N: number of qubits
        #L: sidelength of the lattice

        XX = Matrix{PauliOperator}(undef, L, L)
        for xi in 1:L
            for yi in 1:L
                XX[xi, yi] = genXX(N, index_xyz(xi, yi, 1, L), index_xyz(xi - 1, yi, 2, L))
            end
        end

        YY = Matrix{PauliOperator}(undef, L, L)
        for xi in 1:L
            for yi in 1:L
                YY[xi, yi] = genYY(N, index_xyz(xi, yi, 1, L), index_xyz(xi, yi - 1, 2, L))
            end
        end

        ZZ = Matrix{PauliOperator}(undef, L, L)
        for xi in 1:L
            for yi in 1:L
                ZZ[xi, yi] = genZZ(N, index_xyz(xi, yi, 1, L), index_xyz(xi, yi, 2, L))
            end
        end

        W = Matrix{PauliOperator}(undef, L, L)
        for xi in 1:L
            for yi in 1:L
                is = []
                ps = []
                push!(is, index_xyz(xi, yi, 1, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi, yi - 1, 2, L))
                push!(ps, "X")
                push!(is, index_xyz(xi, yi - 1, 1, L))
                push!(ps, "Y")
                push!(is, index_xyz(xi - 1, yi - 1, 2, L))
                push!(ps, "Z")
                push!(is, index_xyz(xi - 1, yi, 1, L))
                push!(ps, "X")
                push!(is, index_xyz(xi - 1, yi, 2, L))
                push!(ps, "Y")
                W[xi, yi] = genLoop(N, ps, is)
            end
        end

        is = []
        ps = []
        for xi in 1:L
            push!(is, index_xyz(xi, 1, 1, L))
            push!(ps, "Y")
            push!(is, index_xyz(xi, 1, 2, L))
            push!(ps, "Y")
        end
        loop1 = genLoop(N, ps, is)

        is = []
        ps = []
        for yi in 1:L
            push!(is, index_xyz(1, yi, 1, L))
            push!(ps, "X")
            push!(is, index_xyz(1, yi, 2, L))
            push!(ps, "X")
        end
        loop2 = genLoop(N, ps, is)


        # loop1,loop2, loop3, W,
        return XX, YY, ZZ, W, loop1, loop2
    end

    XX, YY, ZZ, W, loop1, loop2 = create_operators_honeycomb(L)

    # Define subsystem function (placeholder)
    subsystem = (cols; cut = :y) -> begin
        lcols = 1:L
        if typeof(cols) == Int
            lcols = cols+1:L
            if lcols == 0
                return []
            end
        elseif typeof(cols) == UnitRange{Int64}
            lcols = cols
            if lcols == 1:0
                return []
            end
        end

        if cut == :a1 ##cuts along a1
            return ([unitcell * (j - 1) + unitcell * L * (k - 1) + l for j in lcols for k in 1:L for l in 1:unitcell])
        elseif cut == :a2 ## cuts along a2
            return ([unitcell * (j - 1) + unitcell * L * (k - 1) + l for j in 1:L for k in lcols for l in 1:unitcell])
        elseif cut == :y
            is = []
            for r in lcols
                for x in 1:L
                    for y in 1:L
                        for s in 1:unitcell
                            if mod1(x + y + s, L) == mod1(r, L)
                                push!(is, unitcell * (x - 1) + unitcell * L * (y - 1) + s)
                            end
                        end
                    end
                end
            end
            return is
        else
            error("Invalid cut type")
        end
    end

    return LatticeHoneycomb(L, unitcell, N, XX, YY, ZZ, W, loop1, loop2, index_xyz, subsystem)
end