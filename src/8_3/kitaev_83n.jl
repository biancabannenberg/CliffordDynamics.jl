include("../Lattices.jl")
include("../GenerateOperators.jl")
abstract type AchtDreiN <: Lattice end

struct AchtDreiN_lattice <: AchtDreiN
  size ::Int
  L_unitcell::Int
  params::Vector{Any}
  index::Int64
  thermalization_steps::Int64
  iterations::Int64
  bewteen_measurements::Int64
end

function xyz2i1(x::Int, y::Int, z::Int, L::Int)
  return 16*(mod1(x, L) -1) + 16*L*(mod1(y, L)-1) + 16*L*L*(mod1(z, L)-1) + 1
end

function xyzs2i(x::Int, y::Int, z::Int, s::Int, L::Int)
  return 16*(mod1(x, L) -1) + 16*L*(mod1(y, L)-1) + 16*L*L*(mod1(z, L)-1) + s
end

function calc_EE(state, l, L)
  subsystem = [16*(x-1) + 16*L*(y-1) + 16*L*L*(z-1) + s for x in l+1:L for y in 1:L for z in 1:L for s in 1:16]
  return entanglement_entropy(state, subsystem, Val(:rref))
end

function S_A(state, lcols, L)
  subsystem = [16*(x-1) + 16*L*(y-1) + 16*L*L*(z-1) + s for x in lcols for y in 1:L for z in 1:L for s in 1:16]
  return entanglement_entropy(state, subsystem, Val(:rref))
end


function IABC(state, L)
  SA = S_A(state, 1:L÷4, L)
  SB = S_A(state, L÷4+1:L÷2, L)
  SC = S_A(state, L÷2+1:3L÷4, L)
  SAB = S_A(state, 1:L÷2, L)
  SBC = S_A(state, L÷4+1:3L÷4, L)
  SAC = S_A(state, vcat(1:L÷4, L÷2+1:3L÷4), L)
  SABC = S_A(state, 1:3L÷4, L)
  return SA + SB + SC - SAB - SBC - SAC + SABC
end

function I_i_j(state, i, j, L)
  SA = S_A(state, [i], L)
  SB = S_A(state, [j], L)
  SAB = S_A(state, [i, j], L)
  return SA + SB - SAB
end

function initialize_state(L::Int)
  # initialize operators
  XX = Array{PauliOperator}(undef, L, L, L, 8)
  for xi in 1:L
    for yi in 1:L
      for zi in 1:L
        XX[xi, yi, zi, 1] = genXX(N, xyzs2i(xi, yi, zi, 1, L), xyzs2i(xi, yi, zi, 10, L))
        XX[xi, yi, zi, 2] = genXX(N, xyzs2i(xi, yi, zi, 2, L), xyzs2i(xi, yi, zi, 16, L))
        XX[xi, yi, zi, 3] = genXX(N, xyzs2i(xi, yi, zi, 3, L), xyzs2i(xi, yi, zi, 11, L))
        XX[xi, yi, zi, 4] = genXX(N, xyzs2i(xi, yi, zi, 4, L), xyzs2i(xi, yi, zi, 13, L))
        XX[xi, yi, zi, 5] = genXX(N, xyzs2i(xi, yi, zi, 5, L), xyzs2i(xi, yi + 1, zi, 15, L))
        XX[xi, yi, zi, 6] = genXX(N, xyzs2i(xi, yi, zi, 6, L), xyzs2i(xi, yi, zi, 14, L))
        XX[xi, yi, zi, 7] = genXX(N, xyzs2i(xi, yi, zi, 7, L), xyzs2i(xi - 1, yi, zi, 12, L))
        XX[xi, yi, zi, 8] = genXX(N, xyzs2i(xi, yi, zi, 8, L), xyzs2i(xi, yi, zi, 9, L))
      end
    end
  end

  YY = Array{PauliOperator}(undef, L, L, L, 8)
  for xi in 1:L
    for yi in 1:L
      for zi in 1:L
        YY[xi, yi, zi, 1] = genYY(N, xyzs2i(xi, yi, zi, 1, L), xyzs2i(xi, yi, zi, 15, L))
        YY[xi, yi, zi, 2] = genYY(N, xyzs2i(xi, yi, zi, 2, L), xyzs2i(xi, yi, zi, 10, L))
        YY[xi, yi, zi, 3] = genYY(N, xyzs2i(xi, yi, zi, 3, L), xyzs2i(xi, yi, zi, 12, L))
        YY[xi, yi, zi, 4] = genYY(N, xyzs2i(xi, yi, zi, 4, L), xyzs2i(xi, yi + 1, zi, 16, L))
        YY[xi, yi, zi, 5] = genYY(N, xyzs2i(xi, yi, zi, 5, L), xyzs2i(xi, yi, zi, 13, L))
        YY[xi, yi, zi, 6] = genYY(N, xyzs2i(xi, yi, zi, 6, L), xyzs2i(xi, yi, zi, 9, L))
        YY[xi, yi, zi, 7] = genYY(N, xyzs2i(xi, yi, zi, 7, L), xyzs2i(xi, yi, zi, 14, L))
        YY[xi, yi, zi, 8] = genYY(N, xyzs2i(xi, yi, zi, 8, L), xyzs2i(xi - 1, yi, zi, 11, L))
      end
    end
  end

  ZZ = Array{PauliOperator}(undef, L, L, L, 8)
  for xi in 1:L
    for yi in 1:L
      for zi in 1:L
        ZZ[xi, yi, zi, 1] = genZZ(N, xyzs2i(xi, yi, zi, 1, L), xyzs2i(xi, yi, zi, 9, L))
        ZZ[xi, yi, zi, 2] = genZZ(N, xyzs2i(xi, yi, zi, 2, L), xyzs2i(xi, yi, zi, 11, L))
        ZZ[xi, yi, zi, 3] = genZZ(N, xyzs2i(xi, yi, zi, 3, L), xyzs2i(xi, yi, zi + 1, 15, L))
        ZZ[xi, yi, zi, 4] = genZZ(N, xyzs2i(xi, yi, zi, 4, L), xyzs2i(xi, yi, zi, 12, L))
        ZZ[xi, yi, zi, 5] = genZZ(N, xyzs2i(xi, yi, zi, 5, L), xyzs2i(xi, yi, zi, 14, L))
        ZZ[xi, yi, zi, 6] = genZZ(N, xyzs2i(xi, yi, zi, 6, L), xyzs2i(xi - 1, yi, zi + 1, 16, L))
        ZZ[xi, yi, zi, 7] = genZZ(N, xyzs2i(xi, yi, zi, 7, L), xyzs2i(xi, yi + 1, zi - 1, 10, L))
        ZZ[xi, yi, zi, 8] = genZZ(N, xyzs2i(xi, yi, zi, 8, L), xyzs2i(xi, yi, zi - 1, 13, L))
      end
    end
  end



  W = Array{PauliOperator}(undef, L, L, L, 12)
  for xi in 1:L
    for yi in 1:L
      for zi in 1:L
        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 1, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 10, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 2, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 16, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 6, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 14, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 7, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi - 1, 12, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi - 1, 3, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 15, L)); push!(ps, "X")
        W[xi, yi, zi, 1] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 11, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 3, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 12, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi, 7, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 10, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 1, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 15, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 5, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 13, L)); push!(ps, "X")
        push!(is, xyzs2i(xi + 1, yi, zi, 8, L)); push!(ps, "X")
        W[xi, yi, zi, 2] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 14, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 6, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 9, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 8, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi - 1, 13, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi - 1, 4, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 16, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 2, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 10, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 7, L)); push!(ps, "X")
        W[xi, yi, zi, 3] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 4, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 13, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 5, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi + 1, zi, 15, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 3, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 11, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 8, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 9, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 6, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi + 1, zi, 16, L)); push!(ps, "X")
        W[xi, yi, zi, 4] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 1, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 10, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 2, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 11, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 3, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 12, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 4, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 13, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 5, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 14, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 6, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 9, L)); push!(ps, "X")
        W[xi, yi, zi, 5] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi - 1, 12, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 7, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 14, L)); push!(ps, "X")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 5, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 15, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 1, L)); push!(ps, "X")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 9, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 8, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 11, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 2, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 16, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi - 1, 4, L)); push!(ps, "X")
        W[xi, yi, zi, 6] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 12, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 3, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 11, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi, 8, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi, 9, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi, 6, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi, 14, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi + 1, yi, zi, 7, L)); push!(ps, "Z")
        W[xi, yi, zi, 7] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 1, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 15, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi - 1, zi, 5, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi - 1, zi, 13, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi - 1, zi, 4, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 16, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 2, L)); push!(ps, "Z")
        push!(is, xyzs2i(xi, yi, zi, 10, L)); push!(ps, "Z")
        W[xi, yi, zi, 8] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 8, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 9, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 1, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 15, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi - 1, 3, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi - 1, 12, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi - 1, 4, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi - 1, 13, L)); push!(ps, "Y")
        W[xi, yi, zi, 9] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi, yi, zi, 7, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 14, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 5, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi + 1, zi, 15, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 3, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 11, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 2, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi + 1, zi - 1, 10, L)); push!(ps, "X")
        W[xi, yi, zi, 10] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi + 1, yi, zi, 8, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 11, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi, zi, 2, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 16, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 6, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 14, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 5, L)); push!(ps, "X")
        push!(is, xyzs2i(xi + 1, yi, zi - 1, 13, L)); push!(ps, "X")
        W[xi, yi, zi, 11] = genLoop(N, ps, is)

        is = []
        ps = []
        push!(is, xyzs2i(xi + 1, yi, zi, 7, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 12, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi, yi, zi, 4, L)); push!(ps, "X")
        push!(is, xyzs2i(xi, yi + 1, zi, 16, L)); push!(ps, "X")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 6, L)); push!(ps, "X")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 9, L)); push!(ps, "X")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 1, L)); push!(ps, "Y")
        push!(is, xyzs2i(xi + 1, yi + 1, zi - 1, 10, L)); push!(ps, "Y")
        W[xi, yi, zi, 12] = genLoop(N, ps, is)
      end
    end
  end


  is = []
  ps = []
  for x in 1:L
    push!(is, xyzs2i(x, 1, 1, 8, L)); push!(ps, "Z")
    push!(is, xyzs2i(x, 1, 1, 9, L)); push!(ps, "Y")
    push!(is, xyzs2i(x, 1, 1, 1, L)); push!(ps, "Y")
    push!(is, xyzs2i(x, 1, 1, 10, L)); push!(ps, "Z")
    push!(is, xyzs2i(x, 1, 1, 2, L)); push!(ps, "X")
    push!(is, xyzs2i(x, 1, 1, 11, L)); push!(ps, "X")
  end
  loop1 = genLoop(N, ps, is)

  is = []
  ps = []
  for y in 1:L
    push!(is, xyzs2i(1, y, 1, 15, L)); push!(ps, "Z")
    push!(is, xyzs2i(1, y, 1, 1, L)); push!(ps, "X")
    push!(is, xyzs2i(1, y, 1, 9, L)); push!(ps, "X")
    push!(is, xyzs2i(1, y, 1, 6, L)); push!(ps, "Z")
    push!(is, xyzs2i(1, y, 1, 14, L)); push!(ps, "Y")
    push!(is, xyzs2i(1, y, 1, 5, L)); push!(ps, "Y")
  end
  loop2 = genLoop(N, ps, is)

  is = []
  ps = []
  for z in 1:L
    push!(is, xyzs2i(1, 1, z - 1, 12, L)); push!(ps, "Z")
    push!(is, xyzs2i(2, 1, z - 1, 7, L)); push!(ps, "Z")
    push!(is, xyzs2i(2, 1, z - 1, 14, L)); push!(ps, "X")
    push!(is, xyzs2i(2, 1, z - 1, 5, L)); push!(ps, "X")
    push!(is, xyzs2i(2, 1, z - 1, 13, L)); push!(ps, "X")
    push!(is, xyzs2i(2, 1, z, 8, L)); push!(ps, "X")
    push!(is, xyzs2i(1, 1, z, 11, L)); push!(ps, "Z")
    push!(is, xyzs2i(1, 1, z, 3, L)); push!(ps, "Z")
  end
  loop3 = genLoop(N, ps, is)
end



function kitaev_83n(px, py, pz, L, keep_result=true, phases = true)
  Nuc = L*L*L
  N = 16*Nuc
  Nsample = 25

  dist = Categorical([px, py, pz]) #X, Y, Z

  



  maxEEsum = 0
  misum = 0
  Itrisum = 0


  for _ in 1:Nsample
    state = MixedDestabilizer(zero(Stabilizer, N))

    # zz measure
    for xi in 1:L
      for yi in 1:L
        for zi in 1:L
          for si in 1:8
            project!(state, ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
          end
        end
      end
    end

    # W operator
    for xi in 1:L
      for yi in 1:L
        for zi in 1:L
          for si in 1:12
            project!(state, W[xi, yi, zi, si], keep_result=keep_result, phases=phases)
          end
        end
      end
    end

    # loop operators
    project!(state, loop1, keep_result=keep_result, phases=phases)
    project!(state, loop2, keep_result=keep_result, phases=phases)
    project!(state, loop3, keep_result=keep_result, phases=phases)

    #println(N - QuantumClifford.trusted_rank(state))

    for _ in 1:20
      for _ in 1:16*L^3
        xi = rand(1:L)
        yi = rand(1:L)
        zi = rand(1:L)
        si = rand(1:8)
        ope = rand(dist)
        if ope == 1 # X
          project!(state, XX[xi, yi, zi, si], keep_result=keep_result, phases=phases)
        elseif ope == 2 # Y
          project!(state, YY[xi, yi, zi, si], keep_result=keep_result, phases=phases)
        else # Z
          project!(state, ZZ[xi, yi, zi, si], keep_result=keep_result, phases=phases)
        end
      end
    end

    
    maxEEsum += calc_EE(state, L÷2, L)
    misum += I_i_j(state, 1, L÷2+1, L)
    Itrisum += IABC(state, L)

    GC.gc()
  end

  return maxEEsum/Nsample, misum/Nsample, Itrisum/Nsample

end

function main(args)
  id = parse(Int, args[1])
  label = args[2]
  N = 30
  pxs = []
  pys = []
  pzs = []
  for xi in 0:N
    for yi in xi:N-xi
      zi = N - xi - yi
      push!(pxs, xi/N)
      push!(pys, yi/N)
      push!(pzs, zi/N)
    end
  end
  px, py, pz = collect(zip(pxs, pys, pzs))[id]
  L = 8
  EE, Itwo, Itri = kitaev_83n(px, py, pz, L, false, false)
  outfile = open("./data/kitaev_83n_$label.dat", "w")
  println(outfile, "$px $py $pz $EE $Itwo $Itri")
  close(outfile)
end


if abspath(PROGRAM_FILE) == @__FILE__
  main(ARGS)
end
