using QuantumClifford
include("Types.jl")
include("GenerateOperators.jl")

### Definition of the (10,3)a Lattice ###
struct Lattice10_3a <: AbstractLattice
    L::Int
    unitcell::Int
    N::Int
    XX::Array{PauliOperator,4}
    YY::Array{PauliOperator,4}
    ZZ::Array{PauliOperator,4}
    W::Array{PauliOperator,4}
    loop1::PauliOperator
    loop2::PauliOperator
    loop3::PauliOperator
    index_xyz::Function
    subsystem::Function
end

####### My definition of the lattice ######
#the (10,3)c lattice has a basis of 6 sites
# unitcell = 4

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


function create_lattice_10_3a(L::Int)::Lattice10_3a
    unitcell = 4
    N = L^3 * unitcell

    # Define index function
    index_xyz = (x, y, z, n, Lmod) -> begin
        return unitcell * (mod1(x, Lmod) - 1) +
               unitcell * Lmod * (mod1(y, Lmod) - 1) +
               unitcell * Lmod * Lmod * (mod1(z, Lmod) - 1) + n
    end

    # function that creates the Operators of the 10_3a lattice
    create_operators_10_3a = (L) -> begin
        #variables needed consitently
        #N: number of qubits
        #L: sidelength of the lattice

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

    XX, YY, ZZ, W, loop1, loop2, loop3 = create_operators_10_3a(L)

    # Define subsystem function (placeholder)
    subsystem = (cols; cut=:y) -> begin
        lcols = 1:1
        if typeof(cols) == Int
            lcols = cols+1:L
        elseif typeof(cols) == UnitRange
            lcols = cols
        end
        if cut == :a1 ## cuts through blue bond == ZZ
            return ([unitcell * (j - 1) + unitcell * L * (k - 1) + unitcell * L^2 * (m - 1) + l for j in lcols for k in 1:L for m in 1:L for l in 1:unitcell])
        elseif cut == :a2 ## cuts through green+blue bond == YY+ZZ
            return ([unitcell * (j - 1) + unitcell * L * (k - 1) + unitcell * L^2 * (m - 1) + l for j in 1:L for k in lcols for m in 1:L for l in 1:unitcell])
        elseif cut == :a3 ## cuts through green bond == YY
            return ([unitcell * (j - 1) + unitcell * L * (k - 1) + unitcell * L^2 * (m - 1) + l for j in 1:L for k in 1:L for m in lcols for l in 1:unitcell])
        elseif cut == :xz
            if lcols == 1:L
                return ([1:N...])
            else
                return ([[(unitcell * (j - 1) + unitcell * L * (k - 1) + unitcell * L^2 * (m - 1) + l)
                          for j in 1:L           # iterate over the x-direction
                          for k in 1:L           # iterate over the y-direction
                          for m in lcols[1]   # iterate over the z-direction
                          for l in 2:unitcell           # select atoms 2 through 6 of the unit cell
                    ]...
                    [(unitcell * (j - 1) + unitcell * L * (k - 1) + unitcell * L^2 * (m - 1) + l)
                     for j in 1:L           # iterate over the x-direction
                     for k in 1:L           # iterate over the y-direction
                     for m in lcols[2:end]    # iterate over the z-direction
                     for l in 1:unitcell           # select atoms 2 through 6 of the unit cell
                    ]...
                    [(unitcell * (j - 1) + unitcell * L * (k - 1) + unitcell * L^2 * (m - 1) + l)
                     for j in 1:L           # iterate over the x-direction
                     for k in 1:L           # iterate over the y-direction
                     for m in mod1(lcols[end] + 1, L)      # iterate over the z-direction
                     for l in 1            # select atoms 2 through 6 of the unit cell
                    ]...])
            end
        elseif cut == :xtry
            global todo = []
            if lcols == 1:L
                return ([1:N...])
            elseif lcols == 1:L÷4
                todo = [(1, 1, 1, 3)]
            elseif lcols == L÷4+1:L÷2
                todo = [(1, 1, 1, 2)]
            elseif lcols == L÷2+1:3L÷4
                todo = [(1, 1, 1, 4)]
            elseif lcols == 3L÷4+1:L
                todo = [(1, 1, 1, 1)]
            elseif lcols == 1:L÷2
                todo = [(1, 1, 1, 3), (1, 1, 1, 2)]
            elseif lcols == vcat(1:L÷4, L÷2+1:3L÷4)
                todo = [(1, 1, 1, 3), (1, 1, 1, 4)]
            elseif lcols == 1:3L÷4
                todo = [(1, 1, 1, 3), (1, 1, 1, 2), (1, 1, 1, 4)]
            end
            is = []
            j = 1
    
            while todo != []
                # println("todo: ", todo)
                if index_xyz(todo[1][1], todo[1][2], todo[1][3], todo[1][4], L) in is
                else
                    push!(is, index_xyz(todo[1][1], todo[1][2], todo[1][3], todo[1][4], L))
                end
    
                if todo[1][4] == 1
                    if index_xyz(todo[1][1] - 1, todo[1][2], todo[1][3], 4, L) in is
                        nothing
                    else
                        #first blue ZZ connection
                        push!(todo, (todo[1][1] - 1, todo[1][2], todo[1][3], 4))
                    end
    
                    if index_xyz(todo[1][1], todo[1][2] - 1, todo[1][3], 2, L) in is
                        nothing
                    else
                        #second green YY connection
                        push!(todo, (todo[1][1], todo[1][2] - 1, todo[1][3], 2))
                    end
                end
                if todo[1][4] == 2
                    #first blue ZZ connection
                    if index_xyz(todo[1][1], todo[1][2] - 1, todo[1][3] + 1, 3, L) in is
                        nothing
                    else
                        push!(todo, (todo[1][1], todo[1][2] - 1, todo[1][3] + 1, 3))
                    end
    
                    #second green YY connection
                    if index_xyz(todo[1][1], todo[1][2] + 1, todo[1][3], 1, L) in is
                        nothing
                    else
                        push!(todo, (todo[1][1], todo[1][2] + 1, todo[1][3], 1))
                    end
                end
                if todo[1][4] == 3
                    #first blue ZZ connection
                    if index_xyz(todo[1][1], todo[1][2] + 1, todo[1][3] - 1, 2, L) in is
                        nothing
                    else
                        push!(todo, (todo[1][1], todo[1][2] + 1, todo[1][3] - 1, 2))
                    end
    
                    #second green YY connection
                    if index_xyz(todo[1][1], todo[1][2], todo[1][3] - 1, 4, L) in is
                        nothing
                    else
                        push!(todo, (todo[1][1], todo[1][2], todo[1][3] - 1, 4))
                    end
                end
                if todo[1][4] == 4
                    #first blue ZZ connection
                    if index_xyz(todo[1][1] + 1, todo[1][2], todo[1][3], 1, L) in is
                        nothing
                    else
                        push!(todo, (todo[1][1] + 1, todo[1][2], todo[1][3], 1))
                    end
    
                    #second green YY connection
                    if index_xyz(todo[1][1], todo[1][2], todo[1][3] + 1, 3, L) in is
                        nothing
                    else
                        push!(todo, (todo[1][1], todo[1][2], todo[1][3] + 1, 3))
                    end
                end
                # if todo[1] == (1,1,1,3)
                #     break
                # end
                
                remove!(todo, todo[1])
                ## break condition
                j += 1
                if j > 256
                    println(j)
                    break
                end
    
            end
            return is 
        else
            error("cut must be :a1, :a2, :a3, :x1, :x2 or :xz")
        end
    end 

    return Lattice10_3a(L, unitcell, N, XX, YY, ZZ, W, loop1, loop2, loop3, index_xyz, subsystem)
end
# l = create_lattice_10_3a(16)

# l.subsystem(1:1; cut=:a1)