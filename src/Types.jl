abstract type AbstractLattice end
 


### Definition of the (10,3)b Lattice ###
struct Lattice10_3_b <: AbstractLattice
    L::Int
    unitcell::Int
    XX::Array{PauliOperator, 2}
    YY::Array{PauliOperator, 2}
    ZZ::Array{PauliOperator, 2}
    W::Array{PauliOperator, 4}
    loop1::PauliOperator
    loop2::PauliOperator
    loop3::PauliOperator
    index_fn::Function
    subsystem_fn::Function
end

### Definition of the (10,3)c Lattice ###
struct Lattice10_3_c <: AbstractLattice
    L::Int
    unitcell::Int
    XX::Array{PauliOperator, 2}
    YY::Array{PauliOperator, 2}
    ZZ::Array{PauliOperator, 2}
    W::Array{PauliOperator, 4}
    loop1::PauliOperator
    loop2::PauliOperator
    loop3::PauliOperator
    index_fn::Function
    subsystem_fn::Function
end