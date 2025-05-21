using QuantumClifford
using Distributions


function xys2i(x::Int, y::Int, s::Int, L::Int)
  return 2*(mod1(x, L)-1) + 2*L*(mod1(y, L)-1) + s
end

function calc_EE(state, l, L)
  subsystem = [2*(j-1) + 2*L*(k-1) + l for j in l+1:L for k in 1:L for l in 1:2]
  return entanglement_entropy(state, subsystem, Val(:rref))
end

function S_A(state, lcols, L)
  subsystem = [2*(j-1) + 2*L*(k-1) + l for j in lcols for k in 1:L for l in 1:2]
  return entanglement_entropy(state, subsystem, Val(:rref))
end

function ss(lcols,L)
  return [2*(j-1) + 2*L*(k-1) + l for j in lcols for k in 1:L for l in 1:2]
end

print(ss(2, 4))

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

function genXX(N, i, j)
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

function honeycomb(p, L, keep_result=true, phases = true)
  Nuc = L*L
  N = 2*Nuc
  Nsample = 25
  dist = Categorical([(1 - p)/2, p, (1 - p)/2])

  # construct, ZZ, XX, YY, W, loop1, loop2
  XX = Matrix{PauliOperator}(undef, L, L)
  for xi in 1:L
    for yi in 1:L
      XX[xi, yi] = genXX(N, xys2i(xi, yi, 1, L), xys2i(xi - 1, yi, 2, L))
    end
  end

  YY = Matrix{PauliOperator}(undef, L, L)
  for xi in 1:L
    for yi in 1:L
      YY[xi, yi] = genYY(N, xys2i(xi, yi, 1, L), xys2i(xi, yi -1 , 2, L))
    end
  end

  ZZ = Matrix{PauliOperator}(undef, L, L)
  for xi in 1:L
    for yi in 1:L
      ZZ[xi, yi] = genZZ(N, xys2i(xi, yi, 1, L), xys2i(xi, yi, 2, L))
    end
  end

 W = Matrix{PauliOperator}(undef, L, L)
  for xi in 1:L
    for yi in 1:L
      is = []
      ps = []
      push!(is, xys2i(xi, yi, 1, L)); push!(ps, "Z")
      push!(is, xys2i(xi, yi - 1, 2, L)); push!(ps, "X")
      push!(is, xys2i(xi, yi - 1, 1, L)); push!(ps, "Y")
      push!(is, xys2i(xi - 1, yi - 1, 2, L)); push!(ps, "Z")
      push!(is, xys2i(xi - 1, yi, 1, L)); push!(ps, "X")
      push!(is, xys2i(xi - 1, yi, 2, L)); push!(ps, "Y")
      W[xi, yi] = genLoop(N, ps, is)
    end
  end

#=  my own try at the W operator
  W = Matrix{PauliOperator}(undef, L, L)
  for xi in 1:L
    for yi in 1:L
      is = []
      ps = []
      push!(is, xys2i(xi, yi, 1, L)); push!(ps, "Z")
      push!(is, xys2i(xi, yi - 1, 2, L)); push!(ps, "X")
      push!(is, xys2i(xi+1, yi - 1, 1, L)); push!(ps, "Y")
      push!(is, xys2i(xi + 1, yi - 1, 2, L)); push!(ps, "Z")
      push!(is, xys2i(xi + 1, yi, 1, L)); push!(ps, "X")
      push!(is, xys2i(xi, yi, 2, L)); push!(ps, "Y")
      W[xi, yi] = genLoop(N, ps, is)
    end
  end=#

  is = []
  ps = []
  for xi in 1:L
    push!(is, xys2i(xi, 1, 1, L)); push!(ps, "Y")
    push!(is, xys2i(xi, 1, 2, L)); push!(ps, "Y")
  end
  loop1 = genLoop(N, ps, is)

  is = []
  ps = []
  for yi in 1:L
    push!(is, xys2i(1, yi, 1, L)); push!(ps, "X")
    push!(is, xys2i(1, yi, 2, L)); push!(ps, "X")
  end
  loop2 = genLoop(N, ps, is)


  EE_sum = 0
  I3_sum = 0


  for _ in 1:Nsample
    state = MixedDestabilizer(zero(Stabilizer, N))

    for xi in 1:L
      for yi in 1:L
        project!(state, XX[xi, yi], keep_result=keep_result, phases=phases)
      end
    end

    for xi in 1:L
      for yi in 1:L
        project!(state, W[xi, yi], keep_result=keep_result, phases=phases)
      end
    end

    project!(state, loop1, keep_result=keep_result, phases=phases)
    project!(state, loop2, keep_result=keep_result, phases=phases)

    if QuantumClifford.trusted_rank(state) != N
      println(N - QuantumClifford.trusted_rank(state))
      error("rank error")
    end


    for _ in 1:100
      for _ in 1:N
        xi = rand(1:L)
        yi = rand(1:L)
        ope = rand(dist)
        if ope == 1 # X
          project!(state, XX[xi, yi], keep_result=keep_result, phases=phases)
        elseif ope == 2 # Y
          project!(state, YY[xi, yi], keep_result=keep_result, phases=phases)
        else # Z
          project!(state, ZZ[xi, yi], keep_result=keep_result, phases=phases)
        end
      end
    end

    
    EE_sum += calc_EE(state, L÷2, L)
    I3_sum += IABC(state, L)

    GC.gc()
  end

  return EE_sum/Nsample, I3_sum/Nsample

end

function main(args)
  id = parse(Int, args[1])
  savepath = args[2]

  L = 12
  p = id/120

  EE, I3 = honeycomb(p, L, false, false)

  outfile = open(savepath, "w")
  println(outfile, "$p $EE $I3")
  close(outfile)
end





#main(["", "try"])
function main2(ps, path)
  #id = parse(Int, args[1])
  savepath = path#args[2]

  L = 12
  #p = id/120

  outfile = open(savepath, "w")
  Threads.@threads for p in ps
    EE, I3 = honeycomb(p, L, false, false)

    println(outfile, "$p $EE $I3")
  end
  close(outfile)
end

main2(0:0.1:1,"trytreads")
12*12
#main2(0:0.01:1,"try2")


if abspath(PROGRAM_FILE) == @__FILE__
  main(ARGS)
end
