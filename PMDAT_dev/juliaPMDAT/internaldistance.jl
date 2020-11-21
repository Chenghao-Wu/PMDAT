using Chemfiles
using LinearAlgebra
using Printf
using ProgressMeter
#using Threads

export R_internal


# using the equation of https://en.wikipedia.org/wiki/Gyration_tensor to calculate the gyration tensor \langle Rg^2 \rangle

function R_internal(trj; output::String="", selection::String="",n_block::Int64=1,blocksize::Int64=10 ,exponentialbase::Float64=1, timeunit::Float64=0.03, num_molecule::Int64=1)

    if selection!=""
        Selection_=Selection(selection)
    else
        Selection_=Selection("atoms: all")
    end

    if output==""
        output="rg.dat"
    end

    println("# Mean Square Radius of Gyration");
    println("# For atoms "*selection);

    trajectory=trj
    frame = read(trajectory)
    natoms::Int64 = size(frame)
    nsteps::Int64 = length(trajectory)
    println("# There are $(size(frame)) atoms in the frame")
    #positions(frame)
    
    expii=0
    displacement_limit=0

    NumberTimeGaps=n_block   +   blocksize
    NumberSteps=n_block   *   blocksize+1

    timeii  =   1
    block_starttime=0
    TimeList=zeros(Float64,NumberSteps)
    for blockii =1:n_block
        for expstepii =1:blocksize
            timeii+=1
            if exponentialbase^(expstepii-1) <= expstepii
                TimeList[timeii] = block_starttime+expstepii*timeunit
            else
                TimeList[timeii] = block_starttime+floor(exponentialbase^(expstepii-1))*timeunit
            end
        end
        block_starttime = TimeList[timeii]
    end

    #=
    calculations
    =#
    selected=convert(Array{Int64,1},evaluate(Selection_,frame)).+1
    n_selectedatoms=length(selected)
    n_atominmolecule::Int64=n_selectedatoms/num_molecule

    R_internal = zeros(Float64,n_atominmolecule)
    p = Progress(nsteps)

    R_internal = Array{Array{Float64,1}}(undef,nsteps)

    @views for timeii=0:nsteps-1

        distance=zeros(Float64,n_atominmolecule)

        frame = read_step(trajectory,timeii)
        selected=convert(Array{Int64,1},evaluate(Selection_,frame)).+1
        n_selectedatoms=length(selected)
        n_atominmolecule=n_selectedatoms/num_molecule
        PositionArray=positions(frame)[:,selected]
        for moleculeii = 1:num_molecule
            atomindexi=(moleculeii-1)*n_atominmolecule+1
            atomindexj=moleculeii*n_atominmolecule
            SpeciesPosArray=PositionArray[:,atomindexi:atomindexj]
            
            for deltaindex=1:n_atominmolecule-1
                index::Int64=n_atominmolecule-deltaindex
                for atomii=1:index
                    atom_this=SpeciesPosArray[:,atomii]
                    atom_that=SpeciesPosArray[:,atomii+deltaindex]
                    distance[deltaindex+1]+=sum((atom_that.-atom_this).^2)/index/deltaindex
                end
            end
        end
        R_internal[timeii+1]=distance/num_molecule
        next!(p)
    end

    R_internal=(sum(R_internal,dims=1)/nsteps)[1]

    io = open(output, "w")
    for stepii=1:n_atominmolecule
        index=stepii
        distanceii=@sprintf("%.4f", R_internal[stepii])
        write(io, "$index"*" "*"$distanceii"*"\n")
    end
    close(io)

    return nothing
end
