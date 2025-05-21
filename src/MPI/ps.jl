using MPI

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
l = 8

send_buf = nothing
recv_buf = Vector{Float64}(undef, l)

if rank == 0
    px = collect(range(0,1,step=0.05))
    s = length(px)
    
    send_buf = reshape(collect(range(0,1,length = l*size)), l, size)
    print("Original data : ", send_buf, "\n")
end

MPI.Scatter!(send_buf, recv_buf, comm; root=0)
print("I got this on rank $(rank): \n $(recv_buf)\n")

siz = 5
s = 7



# reshape(collect(range(0,1,length = siz*siz)), siz, siz)