# include("Lattice_Honeycomb.jl")
# include("Observables.jl")
using Distributions
using DataFrames

#initializes state of Length L with complete tableau of stabilizers
function initialize2D(l::AbstractLattice; keep_result=false, phases=false)
    #L length of system
    #N number of qubits
    state = MixedDestabilizer(zero(Stabilizer, l.N))

    for x in eachindex(l.XX)
        project!(state, l.XX[x], keep_result=keep_result, phases=phases)
    end
    for y in eachindex(l.YY)
        project!(state, l.YY[y], keep_result=keep_result, phases=phases)
    end
    for z in eachindex(l.ZZ)
        project!(state, l.ZZ[z], keep_result=keep_result, phases=phases)
    end

    for i in eachindex(l.W)
        project!(state, l.W[i], keep_result=keep_result, phases=phases)
    end

    project!(state, l.loop1, keep_result=keep_result, phases=phases)
    project!(state, l.loop2, keep_result=keep_result, phases=phases)

    if QuantumClifford.trusted_rank(state) != l.N
        println("Tableau does not have sufficient Rank: $(QuantumClifford.trusted_rank(state)) / $(l.N)")
        # return state
    else
        return state
    end
end

# lat = create_lattice_10_3a(4)
# # using BenchmarkTools
# st = (initialize(lat))

function simulate2D(l::AbstractLattice, p; cut=[:y], iterations::Int=25, thermalization_steps=20, keep_result=false, phases=false, arc=false)
    #probability distribution used : p=[p(XX), p(YY), p(ZZ)]
    dist = Distributions.Categorical([p[1], p[2], p[3]])

    #initialize sum of entanglement entropy, mutual information and tripartite mutual information 
    #:a1, :a2, :a3, :x1, :x2 or :xz
    data = Dict(vcat(["EE_" .* String.(cut) .=> 0.0]..., ["TMI_" .* String.(cut) .=> 0.0]...))
    if arc
        names = "arc_" .* String.(cut)  # ergibt ["arc_a1", "arc_a2", ...]
        arc_data = DataFrame([name => zeros(Float32, l.L+1) for name in names])
    end


    state = initialize2D(l)
    for j in 1:iterations
        #first initialize state
        state = initialize(l)

        #thermalization steps
        for _ in 1:thermalization_steps
            #iterate randomly through all qubits
            for _ in 1:l.N
                xi = rand(1:l.L)
                yi = rand(1:l.L)
                ope = rand(dist)
                if Int(l.unitcell / 2) == 1
                    if ope == 1 #do X operation
                        project!(state, l.XX[xi, yi], keep_result=keep_result, phases=phases)
                    elseif ope == 2 #do Y operation
                        project!(state, l.YY[xi, yi], keep_result=keep_result, phases=phases)
                    else #do Z operation
                        project!(state, l.ZZ[xi, yi], keep_result=keep_result, phases=phases)
                    end
                else
                    si = rand(1:Int(l.unitcell / 2))
                    if ope == 1 #do X operation
                        project!(state, l.XX[xi, yi, si], keep_result=keep_result, phases=phases)
                    elseif ope == 2 #do Y operation
                        project!(state, l.YY[xi, yi, si], keep_result=keep_result, phases=phases)
                    else #do Z operation
                        project!(state, l.ZZ[xi, yi, si], keep_result=keep_result, phases=phases)
                    end
                end
            end
        end

        for (_, c) in enumerate(cut)
            data["EE_"*String(c)] += calc_EE(state, Int(l.L // 2), l, cut=c) / iterations
            data["TMI_"*String(c)] += TMI(state, l, cut=c) / iterations
            if arc
                arc_data[!, "arc_"*String(c)] += entanglement_arc(state, l, cut=c) / iterations
            end
        end
    end

    if arc
        return DataFrame(data), arc_data
    else
        return DataFrame(data)
    end

end