## compass_initalizations_terminal.jl
#using Plots
using ProgressMeter
using LinearAlgebra
using DelaunayTriangulation
using DelimitedFiles
using MPI
include("Compass.jl")

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

function initialize_state(L; method="prepared")
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
    elseif method == "mixed"
        return MixedDestabilizer(one(Stabilizer, L^2))
    elseif method == "zero"
        return MixedDestabilizer(zero(Stabilizer, L^2))
    else 
        error("Invalid method")
    end
end


function main2(p, L; initialize = "mixed", checkrank = false, keep_result=false, phases=false, iterations = 2000, thermalization = 500, between_measurements = 200) 
    X,Y,Z = create_operators(L)
    N = L^2
    dist  = Distributions.Categorical([p[1],p[2],p[3]])
    Sx = []
    Sy = []
    Tx = []
    Ty = []
    r = []

    #if initialize == "mixed"
    #    state = MixedDestabilizer(one(Stabilizer, N))
    #elseif initialize == "zero"
    #    state = MixedDestabilizer(zero(Stabilizer, N))
    #end

    state = initialize_state(L, method=initialize)

    
    #for i in length(X)
    #    project!(state, X[i], keep_result=false, phases=false)
    #end

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
        if checkrank == true
            push!(r, QuantumClifford.trusted_rank(state))
        end
        #if check == true
        #    println("trusted Rank = ", QuantumClifford.trusted_rank(state))
        #end
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
    #println(QuantumClifford.trusted_rank(state))
    return Sx,Sy,Tx,Ty, r#, state
end

## calculate averages of multiple observables
L = 12
method = "zero"
thermalization = 500
iterations = 2000
between_measurements = 200
Sx_avg = []
Sy_avg = []
Tx_avg = []
Ty_avg = []
r_avg = []
Sx_var = []
Sy_var = []
Tx_var = []
Ty_var = []
r_var = []

px = 1:-0.01:0
py = 1 .-px
pz = 0.0

@showprogress   for (i,p) in enumerate(px)
    Sx, Sy, Tx, Ty, r = main2([p, 1-p, pz], L; initialize = method, checkrank = false, iterations = iterations, thermalization = thermalization, between_measurements = between_measurements)
    open("data/initializations_L$(L)_"*(method)*"_p_meanvar_SxSyTxTyr.txt", "a") do io
        writedlm(io, [p mean(Sx) var(Sx) mean(Sy) var(Sy) mean(Tx) var(Tx) mean(Ty) var(Ty) mean(r) var(r)])
    end
end
