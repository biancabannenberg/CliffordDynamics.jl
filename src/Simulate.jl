# include("Lattice_10_3a.jl")

#initializes state of Length L with complete tableau of stabilizers
function initialize(l::AbstractLattice; keep_result=false, phases=false)
    #L length of system
    #N number of qubits
    state = MixedDestabilizer(zero(Stabilizer, l.N))

    # for x in eachindex(l.XX)
    #     project!(state, l.XX[x], keep_result=keep_result, phases=phases)
    # end
    # for y in eachindex(l.YY)
    #     project!(state, l.YY[y], keep_result=keep_result, phases=phases)
    # end
    for z in eachindex(l.ZZ)
        project!(state, l.ZZ[z], keep_result=keep_result, phases=phases)
    end

    for i in eachindex(l.W)
        project!(state, l.W[i], keep_result=keep_result, phases=phases)
    end

    project!(state, l.loop1, keep_result=keep_result, phases=phases)
    project!(state, l.loop2, keep_result=keep_result, phases=phases)
    project!(state, l.loop3, keep_result=keep_result, phases=phases) #this loop is not needed for initialization if all x,y,z are measured

    if QuantumClifford.trusted_rank(state) != l.N
        println("Tableau does not have sufficient Rank: $(QuantumClifford.trusted_rank(state)) / $(l.N)")
        # return state
    else
        return state
    end
end

lat = create_lattice_10_3a(4)
using BenchmarkTools
@benchmark (initialize(lat))

function simulate(l::AbstractLattice, p; cut=[:y], iterations::Int=25, thermalization_steps=20, keep_result=false, phases=false)
    #probability distribution used : p=[p(XX), p(YY), p(ZZ)]
    dist = Distributions.Categorical([p[1], p[2], p[3]])

    #initialize sum of entanglement entropy, mutual information and tripartite mutual information 
    #:a1, :a2, :a3, :x1, :x2 or :xz
    data = Dict(vcat(["EE_" .* String.(cut) .=> 0.0]..., ["TMI_" .* String.(cut) .=> 0.0]...))

    for _ in 1:iterations
        #first initialize state
        state = initialize(l)

        #thermalization steps
        for _ in 1:thermalization_steps
            #iterate randomly through all qubits
            for _ in 1:l.N
                xi = rand(1:l.L)
                yi = rand(1:l.L)
                zi = rand(1:l.L)
                si = rand(1:Int(l.unitcell / 2))
                ope = rand(dist)
                if ope == 1 #do X operation
                    project!(state, l.XX[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                elseif ope == 2 #do Y operation
                    project!(state, l.YY[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                else #do Z operation
                    project!(state, l.ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
                end
            end
        end

        for (_, c) in enumerate(cut)
            # data["EE_"*String(c)] += calc_EE(state, Int(L // 2), L, cut=c) / iterations
            data["TMI_"*String(c)] += TMI(state, l.L, cut=c) / iterations
        end

    end

    return DataFrame(data)
end
