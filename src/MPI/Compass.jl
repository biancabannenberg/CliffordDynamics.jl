using QuantumClifford
using Distributions

L_unitcell = 1

function index_xy(x::Int, y::Int, L::Int)
    return mod1(x,L) + L * (mod1(y,L)-1)
end

#function that returns subsystem of a 10_3_c lattice
#general structure of subsystem functions:
#    L_unitcell * (x-1) + L_unitcell * L * (y-1) + L_unitcell * L^2 * (y-1)   --> iterate through all x,y and z with x in l+1:L, and y,z in 1:L
function ss_compass(lcols::Rational{Int}, L::Int, dir::String)
    if dir == "x"
        return ([L_unitcell*(j-1) + L_unitcell*L*(k-1) + l for j in lcols+1:L for k in 1:L for l in 1:L_unitcell])
    elseif dir == "y"
        return ([L_unitcell*(j-1) + L_unitcell*L*(k-1) + l for j in 1:L for k in lcols+1:L for l in 1:L_unitcell])
    else    
        error("dir must be x or y")
    end
end

function ss_compass(lcols, L::Int, dir::String)
    if dir == "x"
        return ([L_unitcell*(j-1) + L_unitcell*L*(k-1) + l for j in lcols for k in 1:L for l in 1:L_unitcell])
    elseif dir == "y"
        return ([L_unitcell*(j-1) + L_unitcell*L*(k-1) + l for j in 1:L for k in lcols for l in 1:L_unitcell])
    else
        error("dir must be x or y")
    end
end

#function that calculates entanglement entropy of subsystem of length l
#firstly one need to define the subsystem of length l and total systemsize L (ss_10_3_c(l,L))
#state: MixedDestabilizer
function calc_EE(state, lcols, L, dir::String)
    return entanglement_entropy(state, ss_compass(lcols,L, dir), Val(:rref))
end

#function that calculates entanglement for given subdivistion lcols
function S_A(state, lcols, L, dir::String)
    return entanglement_entropy(state, ss_compass(lcols,L, dir), Val(:rref))
end

#function that imoplements the tripartite mutual information
#in this definition: system partet in 4 subsystems going from 1:L/4, L/4+1:L/2, L/2+1:3L/4, 3L/4:L
# then: I= SA + SB + SC - SAB - SAC - SBA + SAB
function TMI(state, L, dir::String)
    if L % 4 != 0
        error("L must be a multiple of 4")
    end
    SA = S_A(state, 1:L÷4, L, dir)
    SB = S_A(state, L÷4+1:L÷2, L, dir)
    SC = S_A(state, L÷2+1:3L÷4, L, dir)
    SAB = S_A(state, 1:L÷2, L, dir)
    SBC = S_A(state, L÷4+1:3L÷4, L, dir)
    SAC = S_A(state, vcat(1:L÷4, L÷2+1:3L÷4), L, dir)
    SABC = S_A(state, 1:3L÷4, L, dir)
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

function create_operators(L::Int)
    N = L^2 * L_unitcell

    #construct XX,YY,ZZ operators
    XX = Array{PauliOperator}(undef, L,L)
    for xi in 1:L
        for yi in 1:L
            XX[xi,yi] = genXX(N, index_xy(xi,yi, L), index_xy(xi+1,yi, L))
        end
    end

    YY = Array{PauliOperator}(undef, L,L)
    for xi in 1:L
        for yi in 1:L
            YY[xi,yi] = genYY(N, index_xy(xi,yi, L), index_xy(xi,yi+1, L))
        end
    end

    if mod(L,2) != 0
        error("L must be even")
    end
    L_2 = Int(L/2)
    ZZ = Array{PauliOperator}(undef, L,L_2)
    for xi in 1:L
        for yi in 1:L_2
            if mod(xi,2) != 0 
                ZZ[xi,yi] = genZZ(N, index_xy(xi,1+(yi-1)*2, L), index_xy(xi+1,1+(yi-1)*2+1, L))
            end
            if mod(xi,2) == 0
                ZZ[xi,yi] = genZZ(N, index_xy(xi,1+(yi-1)*2, L), index_xy(xi+1,1+(yi-1)*2-1, L))
            end
        end
    end

    return XX, YY, ZZ  
end


#tried all Z measurements but did not work
function Zs(N)
    Xarr = falses(N)
    Zarr = trues(N)
    return PauliOperator(0x00, Xarr, Zarr)
end

function Y_line(L)
    N = L^2
    Xarr = falses(N)
    Zarr = falses(N)
    for i in 1:L
        Xarr[i] = true
        Zarr[i] = true
    end
    return PauliOperator(0x00, Xarr, Zarr)
end

function initialize_state(L; method="prepared", basis = :Z)
    if method == "prepared"
        N = L^2

        state = MixedDestabilizer(zero(Stabilizer, L^2))
        x, y, z = create_operators(L)
        yl = Y_line(L)
        zs = Zs(N)
        
        println(QuantumClifford.trusted_rank(state))

        #project all XX operators
        for i in 1:L
            for j in 1:L
                project!(state, x[i, j]; keep_result=false, phases=false)
            end
        end
        println(QuantumClifford.trusted_rank(state))
        #project all YY operators
        for k in 1:L
            for l in 1:L
                project!(state, y[k, l]; keep_result=false, phases=false)
            end
        end
        println(QuantumClifford.trusted_rank(state))

        #project a line of Y operators
        #project!(state, yl)
        project!(state,yl)
        println(QuantumClifford.trusted_rank(state))

        #if not full rank, throw an error
        if QuantumClifford.trusted_rank(state) != N
            error("Rank is not N")
        end

        return state
    elseif method == "pure" ## method used to be called mixed ->> now pure as for pure state 
        return MixedDestabilizer(one(Stabilizer, L^2; basis = basis))
    elseif method == "mixed" ## method used to be called zero ->> now mixed as for mixed state
        return MixedDestabilizer(zero(Stabilizer, L^2))
    else 
        error("Invalid method")
    end
end

"""
   Function that performs the Monte Carlo simulation for the Compass model
    
    parameters:
    p: array of probabilities for X, Y and Z operations
    L: system size
    initialize: method to initialize the state, default is "mixed"
                methods: "mixed", "pure", "prepared"
    keep_result: keep the result of the measurement, default is false
    phases: keep the phases of the measurement, default is false
    iterations: number of iterations, default is 2000
    thermalization: number of thermalization steps, default is 500
    between_measurements: number of steps between measurements, default is 200

    returns:
    Sx: array of entanglement entropy in x direction
    Sy: array of entanglement entropy in y direction
    Tx: array of tripartite mutual information in x direction
    Ty: array of tripartite mutual information in y direction
    r: array of system rank of the state
"""
function main(p, L;initialize = "mixed", basis = :Z, keep_result=false, phases=false, iterations = 2000, thermalization = 500, between_measurements = 200, newstarts = 1, prep_cat = false) 
    X,Y,Z = create_operators(L)
    N = L^2
    dist  = Distributions.Categorical([p[1],p[2],p[3]])
    Sx = []
    Sy = []
    Tx = []
    Ty = []
    r = []
    
    for _ in 1:newstarts
        state = initialize_state(L, method=initialize, basis = basis)
        
        if p[1] > p[2] && prep_cat
            for i in 1:length(X)
                project!(state, X[i], keep_result=false, phases=false)
            end 
        elseif p[1] < p[2] && prep_cat
            for i in 1:length(Y)
                project!(state, Y[i], keep_result=false, phases=false)
            end
        end
        # for i in 1:L
        #     for j in 1:L
        #         project!(state, Z[i, j]; keep_result=false, phases=false)
        #     end
        # end
        for i in 1:thermalization
            for j in 1:N
                xi = rand(1:L)
                yi = rand(1:L)

                ope = rand(dist)

                if ope == 1 #do X operation
                    project!(state, X[xi,yi], keep_result=keep_result, phases=phases)
                elseif ope == 2 #do Y operation
                    project!(state, Y[xi,yi], keep_result=keep_result, phases=phases)
                else #do Z operation
                    yi = rand(1:Int(L/2))
                    project!(state, Z[xi,yi], keep_result=keep_result, phases=phases)
                end
            end
        end

        for _ in 1:iterations
            for _ in 1:N
                xi = rand(1:L)
                yi = rand(1:L)

                ope = rand(dist)

                if ope == 1 #do X operation
                    project!(state, X[xi,yi], keep_result=false, phases=false)
                elseif ope == 2 #do Y operation
                    project!(state, Y[xi,yi], keep_result=false, phases=false)
                else #do Z operation
                    yi = rand(1:Int(L/2))
                    project!(state, Z[xi,yi], keep_result=false, phases=false)
                end
            end
            push!(Sx, calc_EE(state, 1:Int(L/2), L, "x"))
            push!(Sy, calc_EE(state, 1:Int(L/2), L, "y"))
            push!(Tx, TMI(state, L, "x"))
            push!(Ty, TMI(state, L, "y"))
            push!(r, QuantumClifford.trusted_rank(state))
            for _ in 1:between_measurements
                for _ in 1:N
                    xi = rand(1:L)
                    yi = rand(1:L)

                    ope = rand(dist)

                    if ope == 1 #do X operation
                        project!(state, X[xi,yi], keep_result=false, phases=false)
                    elseif ope == 2 #do Y operation
                        project!(state, Y[xi,yi], keep_result=false, phases=false)
                    else #do Z operation
                        yi = rand(1:Int(L/2))
                        project!(state, Z[xi,yi], keep_result=false, phases=false)
                    end
                end
            end
            GC.gc()
        end
    end
    #println(QuantumClifford.trusted_rank(state))
    return Sx, Sy, Tx, Ty, r
end


"""
   Function that calculates the Entanglement Arc of a given systemsize L and p distribution fo the Compass Model
   in performing the Monte Carlo simulation 
    
    parameters:
    p: array of probabilities for X, Y and Z operations
    L: system size
    n: stepsize of the entanglement arc, default is 1
    initialize: method to initialize the state, default is "mixed"
                methods: "mixed", "pure", "prepared"
    keep_result: keep the result of the measurement, default is false
    phases: keep the phases of the measurement, default is false
    iterations: number of iterations, default is 2000
    thermalization: number of thermalization steps, default is 500
    between_measurements: number of steps between measurements, default is 200
    newstarts: number of newstarts of the simulation, default is 1
    prep_cat: prepare the state in a cat state, default is false    

    returns:
    Sx: array of entanglement entropy in x direction
    Sy: array of entanglement entropy in y direction
    Tx: array of tripartite mutual information in x direction
    Ty: array of tripartite mutual information in y direction
    r: array of system rank of the state
"""    
function entanglement_arc(p, L;n = 1, initialize = "mixed", basis = :Z, keep_result=false, phases=false, iterations = 2000, thermalization = 500, between_measurements = 200, newstarts = 1, prep_cat = false) 
    X,Y,Z = create_operators(L)
    N = L^2
    dist  = Distributions.Categorical([p[1],p[2],p[3]])
    Sx = zeros(L+1)
    Sy = zeros(L+1)
    # r = []
    
    for _ in 1:newstarts
        state = initialize_state(L, method=initialize, basis = basis)
        
        if p[1] > p[2] && prep_cat
            for i in 1:length(X)
                project!(state, X[i], keep_result=false, phases=false)
            end 
        elseif p[1] < p[2] && prep_cat
            for i in 1:length(Y)
                project!(state, Y[i], keep_result=false, phases=false)
            end
        end
        for _ in 1:thermalization
            for _ in 1:N
                xi = rand(1:L)
                yi = rand(1:L)

                ope = rand(dist)

                if ope == 1 #do X operation
                    project!(state, X[xi,yi], keep_result=keep_result, phases=phases)
                elseif ope == 2 #do Y operation
                    project!(state, Y[xi,yi], keep_result=keep_result, phases=phases)
                else #do Z operation
                    yi = rand(1:Int(L/2))
                    project!(state, Z[xi,yi], keep_result=keep_result, phases=phases)
                end
            end
        end

        for _ in 1:iterations
            for _ in 1:N
                xi = rand(1:L)
                yi = rand(1:L)

                ope = rand(dist)

                if ope == 1 #do X operation
                    project!(state, X[xi,yi], keep_result=false, phases=false)
                elseif ope == 2 #do Y operation
                    project!(state, Y[xi,yi], keep_result=false, phases=false)
                else #do Z operation
                    yi = rand(1:Int(L/2))
                    project!(state, Z[xi,yi], keep_result=false, phases=false)
                end
            end
            ## calculate entanglement arc
            for  (i,l) in enumerate(0:n:L)
                Sx[i] += (calc_EE(state, 1:l, L, "x"))./newstarts
                Sy[i] += (calc_EE(state, 1:l, L, "y"))./newstarts
            end
            # push!(r, QuantumClifford.trusted_rank(state))
            for _ in 1:between_measurements
                for _ in 1:N
                    xi = rand(1:L)
                    yi = rand(1:L)

                    ope = rand(dist)

                    if ope == 1 #do X operation
                        project!(state, X[xi,yi], keep_result=false, phases=false)
                    elseif ope == 2 #do Y operation
                        project!(state, Y[xi,yi], keep_result=false, phases=false)
                    else #do Z operation
                        yi = rand(1:Int(L/2))
                        project!(state, Z[xi,yi], keep_result=false, phases=false)
                    end
                end
            end
            GC.gc()
        end
    end
    #println(QuantumClifford.trusted_rank(state))
    return Sx, Sy  #, Tx, Ty, r
end