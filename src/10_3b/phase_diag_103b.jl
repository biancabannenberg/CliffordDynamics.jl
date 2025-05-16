# this simulation sweeps the nishimori line
# it calculates the expectation values as well as the mututal information between 
# two quibits in the middle row at 1/4 L and at 3/4 L

##### MPI HEAD #####
using MPI

MPI.Init()
COMM = MPI.COMM_WORLD

RANK = MPI.Comm_rank(COMM)
SIZE = MPI.Comm_size(COMM)

println("Hello, I am ", lpad(RANK, 5), " of $(SIZE). I am running on $(gethostname()).")
MPI.Barrier(COMM)
####################

### load packages and code 
#
#
#
using Dates
using LinearAlgebra
using CSV
using Statistics
using DataFrames
include("10_3_b.jl")

#probabilities
pxs = []
pys = []
pzs = []

for z in range(0,1,step=0.05)
    for x in range(0,round(1-z, digits = 2),step=0.02)
        y = round(1 - x - z, digits=2)
        push!(pxs, x)
        push!(pys, y)
        push!(pzs, z)
    end
end

ps = collect(zip(pxs,pys,pzs))

# define the parameters of the simulation 
L = 8
iterations = 200 
thermalization = 20
cut = :all
# iterations = 1
# between_measurements = 0

#parameters = Dict("L" => L, "method" => method, "basis" => basis, "newstarts" => newstarts, "thermalization" => thermalization, "iterations" => iterations, "between_measurements" => between_measurements)
# header = "#" * "L: $(L), method: $(method), basis: $(basis), CAT: $(prep_cat), newstarts: $(newstarts), iterations: $(iterations), thermalization: $(thermalization), between_measurements: $(between_measurements)"#join(["$(key): $(value)" for (key, value) in parameters], ", ")
header = "#" * "L: $(L), cut: $(cut), iterations: $(iterations), thermalization: $(thermalization)" #join(["$(key): $(value)" for (key, value) in parameters], ", ")

n_samples = length(ps)

# create directories to save data
# if method == "pure"
#     if prep_cat
#         mkpath("data/pxpy/"*method*String(basis)*"/CAT/$(L)/$(newstarts)")
#         pth = "data/pxpy/"*method*String(basis)*"/CAT/$(L)/$(newstarts)/" 
#     else
#         mkpath("data/pxpy/"*method*String(basis)*"/$(L)/$(newstarts)")
#         pth = "data/pxpy/"*method*String(basis)*"/$(L)/$(newstarts)/" 
#     end
# else
#     mkpath("data/pxpy/"*method*"/$(L)/$(newstarts)")
#     pth = "data/pxpy/"*method*"/$(L)/$(newstarts)/" 
# end
# mkpath("data/phaseDiagram/$(L)/$(string(cut))/$(iterations)/")
pth = "data/phaseDiagram/$(L)/$(string(cut))/$(iterations)/" 
mkpath(pth)


timestr() = Dates.format(now(), dateformat"HH:MM:SS") * " --- "

MPI.Barrier(COMM)
RANK == 0 && @info timestr() * "Loaded packages and code on all workers."


if RANK == 0 ##### MANAGER #####
        
    ### do stuff that only the manager should do (eg. define set parameters, ...)
    
    
    
    ### define monitoring function 
    function monitor(t0, jobsdone, n)
        # get runtime in seconds
        runtime = (now() - t0).value / 10^3

        # get estimated time to completion
        total_runtime = runtime / jobsdone * n
        time_done = Dates.format(t0 + Dates.Second(round(Int, total_runtime)), dateformat"dd.mm. - HH:MM:SS")

        @info timestr() * "workers: $(SIZE - 1) --- Samples: " * lpad((jobsdone), 6) * " / " * rpad(n, 6) * " --- Runtime: " * lpad((round(Int, runtime)), 7) * " / " * rpad((round(Int, total_runtime)), 7) * " seconds. Estimated completion time: " * time_done
    end

else ##### WORKER #####
    # # define function to run on all workers
    # function main(task)
    #     # do stuff

    #     # save stuff
    # end
    function main1(task)
        p = ps[task]
        println("Worker $(RANK) is working on sample $(task) / $(n_samples) and p = $p.")
        #do stuff
        try
            if cut == :all
                EE_a1, TMI_a1, EE_a2, Tmi_a2, EE_a3, TMI_a3, EE_z1, TMI_z1, EE_z2, TMI_z2, EE_y, TMI_y = simulate(L, p;cut = cut, iterations=iterations, thermalization_steps=thermalization)
                # df = DataFrame(Sx = Sx,Sy= Sy,Tx= Tx,Ty= Ty,Rank_sys= r)
                
                #save stuff
                # open(pth*"$(p).csv", "w") do io
                #     println(io, header)
                #     println(io, join(names(df), ","))
                #     CSV.write(io, df, append = true)                
                # end

                return (p, EE_a1, TMI_a1, EE_a2, Tmi_a2, EE_a3, TMI_a3, EE_z1, TMI_z1, EE_z2, TMI_z2, EE_y, TMI_y)
            else
                EE, TMI = simulate(L, p;cut = cut, iterations=iterations, thermalization_steps=thermalization)
                # df = DataFrame(Sx = Sx,Sy= Sy,Tx= Tx,Ty= Ty,Rank_sys= r)
                
                #save stuff
                # open(pth*"$(p).csv", "w") do io
                #     println(io, header)
                #     println(io, join(names(df), ","))
                #     CSV.write(io, df, append = true)                
                # end

                return (p, EE, TMI)
            end

        catch 
            println("Worker $(RANK) failed on sample $(task) and p = $p.")
        end
    end        
end

MPI.Barrier(COMM)
RANK == 0 && @info timestr() *  "Loaded main function on all workers."
#von mir einfefügt
motherpath = nothing
if RANK != 0
    motherpath = nothing
end
motherpath = MPI.bcast(motherpath, COMM)
RANK == 0 && @info timestr() *  "Broadcasted motherpath to all workers."
RANK == 0 && @info timestr() *  "Starting simulation."
RANK == 0 && @info timestr() *  "L = $(L), iterations = $(iterations), thermalization = $(thermalization), cut = $(cut), n_samples = $(n_samples)."




##### MPI TAIL #####
"""
    manager(p::Int, n::Int)

The manager distributes n jobs to p-1 workers. Assumes that n >= p-1.
"""

function manager(p::Int, n::Int; t0::DateTime)
    results = []  # Store (p -> mean) pairs
    # used to reduce the number of monitor outputs
    t1 = now()
    min_monitor_interval = 10 # seconds
    
    for j in 1:p-1
        # send the job number to worker j
        MPI.send(j, COMM; dest = j)
    end
    
    jobssent = p-1  # sent already p-1 jobs
    jobsdone = 0    # number of jobs that have been done
    jobsdone_previous_iteration = 0 # used for monitoring progress
    done = 0        # number of workers that have finished
    while done < p-1
        for i in 1:p-1
            messageSent = MPI.Iprobe(COMM; source = i)
            if messageSent
                # report == 1 -> success 
                # report == 0 -> failure
                report, res = MPI.recv(COMM; source = i)

                if report == 1
                    push!(results, res)
                    jobsdone += 1
                end

                # # If report is a tuple (p, mean_value), store it
                # if report isa Tuple{Float64, Float64}
                #     p_val, mean_val = report
                #     results[p_val] = mean_val
                #     println("Manager received mean: $mean_val for p=$p_val")
                # end

                # report == 1 && (jobssent += 1; jobsdone += 1)      # only count successful jobs

                if jobssent >= n
                    MPI.send(-1, COMM; dest = i)
                    done += 1
                else
                    jobssent += 1
                    MPI.send(jobssent, COMM; dest = i)
                end
            end
        end


        # print progress
        if (jobsdone > jobsdone_previous_iteration && (now() - t1).value / 10^3 > min_monitor_interval) || jobsdone == n
            monitor(t0, jobsdone, n)
            jobsdone_previous_iteration = jobsdone
            t1 = now()
        end
    end


    if cut == :all
        # Save collected means to file
        open(pth*"mean.csv", "w") do io
            println(io, header)
            println(io, "#DATAFRAME.p = eval.(Meta.parse.(DATAFRAME.p)) to convert String into Tuple")
            println(io, "p,EE_a1,TMI_a1,EE_a2,Tmi_a2,EE_a3,TMI_a3,EE_z1,TMI_z1,EE_z2,TMI_z2,EE_y,TMI_y")
            CSV.write(io, DataFrame(results), append = true)
        end
    else
        # Save collected means to file
        open(pth*"mean.csv", "w") do io
            println(io, header)
            println(io, "#DATAFRAME.p = eval.(Meta.parse.(DATAFRAME.p)) to convert String into Tuple")
            println(io, "p,EE,TMI")
            CSV.write(io, DataFrame(results), append = true)
        end
    end
    println("Manager saved mean values to summary file.")
end

# a = collect(1:5)
# b = collect(10:18)
# c = collect(19:27)

# println(mean.([a, b, c]))

"""
    worker(i::Int)

The i'th worker receives a number. The worker terminates if the number is -1,
otherwise it sends a report to the manager.

`report == 1` -> job was successful
`report == 0` -> job failed
"""
function worker(i::Int)

    while true
        # receive a task number
        task = MPI.recv(COMM; source = 0)
        
        # terminate if task number is -1
        task == -1 && break
        
        # initialize report variable
        report = -1
        
        means = nothing

        # do the computation
        try 
            means = main1(task)
            report = 1
        catch
            report = 0
            # ps = nothing
            # means = nothing
        end

        # send result to manager
        MPI.send((report, means), COMM; dest = 0)
    end
end

if RANK == 0
    manager(SIZE, n_samples; t0 = now())
else
    worker(RANK)
end
MPI.Barrier(COMM)
# println(S)
RANK == 0 && (@info timestr() * "All done. :)")
####################