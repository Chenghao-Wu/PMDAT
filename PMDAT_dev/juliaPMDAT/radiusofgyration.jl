using Chemfiles
using LinearAlgebra
using Printf
using ProgressMeter
#using Threads

export Rg


# using the equation of https://en.wikipedia.org/wiki/Gyration_tensor to calculate the gyration tensor \langle Rg^2 \rangle

function Rg(trj; output::String="", selection::String="",n_block::Int64=10,blocksize::Int64=10 ,exponentialbase::Float64=1.2, timeunit::Float64=0.03, num_molecule::Int64=1)

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

    Rg = Array{Array{Float64,2}}(undef,nsteps)
    p = Progress(NumberSteps)

    for timeii=0:nsteps-1
        frame = read_step(trajectory,timeii)
        selected=convert(Array{Int64,1},evaluate(Selection_,frame)).+1
        n_selectedatoms::Int64=length(selected)
        n_atominmolecule::Int64=n_selectedatoms/num_molecule

        PositionArray=positions(frame)[:,selected]

        FrameResult=Array{Array{Float64,2}}(undef,num_molecule)

        for moleculeii = 1:num_molecule
            atomindexi::Int64=(moleculeii-1)*n_atominmolecule+1
            atomindexj::Int64=(moleculeii)*(n_atominmolecule)
            SpeciesPosArray=PositionArray[:,atomindexi:atomindexj]
            species_result=Array{Array{Float64,2}}(undef,n_atominmolecule)
            for atomii=1:n_atominmolecule
                atomii_array=SpeciesPosArray[:,atomii]
                tensor=(SpeciesPosArray.-atomii_array)*(SpeciesPosArray.-atomii_array)'
                species_result[atomii]=tensor
            end
            FrameResult[moleculeii]=sum(species_result,dims=1)[1]./((n_atominmolecule^2)*2)
        end
        Rg[timeii+1]=sum(FrameResult,dims=1)[1]./num_molecule
        next!(p)
    end

    io = open(output, "w")
    for stepii=1:nsteps
        time=TimeList[stepii]
        timegap=@sprintf("%.4f", time)
        xx=@sprintf("%.4f", Rg[stepii][1])
        yy=@sprintf("%.4f", Rg[stepii][5])
        zz=@sprintf("%.4f", Rg[stepii][9])
        xy=@sprintf("%.4f", Rg[stepii][2])
        xz=@sprintf("%.4f", Rg[stepii][3])
        yz=@sprintf("%.4f", Rg[stepii][6])
        write(io, "$timegap"*" "*"$xx"*" "*"$yy"*" "*"$zz"*" "*"$xy"*" "*"$xz"*" "*"$yz"*"\n")
    end
    close(io)

    return nothing
end
