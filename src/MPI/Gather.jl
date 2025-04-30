using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

send_buf = nothing

if rank == 0
    send_buf = collect(1:size) .* 100
    print("Original array on rank 0: \n $(send_buf) \n")
end

v = MPI.Scatter(send_buf, Int, comm; root=0)
print("I got this on rank $(rank): \n $(v) \n")

v = v*v
recv_buf = MPI.Gather(v, comm; root=0)
if rank == 0
    print("New array on rank 0: \n $(recv_buf) \n")
end