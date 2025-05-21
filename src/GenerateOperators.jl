### This file generates the Pauli operators for a system of N qubits.
### The Pauli operators are represented as bit arrays for X and Z operators.

### The functions genXX, genYY, and genZZ generate the corresponding Pauli operators
### for two qubits i and j. 
### The indizes of the qubits depend on the chosen lattice structure, and can be taken from the ... file.

### The function genLoop generates a Pauli operator for a list of qubits and their corresponding operators.


function genXX(N::Int, i::Int, j::Int)
    Xarr = falses(N)
    Zarr = falses(N)
    Xarr[i] = true
    Xarr[j] = true
    return PauliOperator(0x00, Xarr, Zarr)
  end
  
  function genYY(N::Int, i::Int, j::Int)
    Xarr = falses(N)
    Zarr = falses(N)
    Zarr[i] = true
    Zarr[j] = true
    Xarr[i] = true
    Xarr[j] = true
    return PauliOperator(0x00, Xarr, Zarr)
  end
  
  function genZZ(N::Int, i::Int, j::Int)
    Xarr = falses(N)
    Zarr = falses(N)
    Zarr[i] = true
    Zarr[j] = true
    return PauliOperator(0x00, Xarr, Zarr)
  end
  
  function genLoop(N::Int, ps, is)
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