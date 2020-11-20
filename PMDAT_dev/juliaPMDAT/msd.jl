using Chemfiles
using LinearAlgebra
using Printf
using ProgressMeter

export msd

function msd(trj;output::String="",n_block::Int64=10,blocksize::Int64=10,selection::String="",exponentialbase::Float64=1.2,timeunit::Float64=0.03)

    if selection!=""
        selection=Selection(selection)
    end

    if output==""
        output="msd.dat"
    else
        output=output
    end

    println("# Mean Square Deviation in");
    println("# For atoms "*selection);

    trajectory=trj
    frame = read(trajectory)
    natoms = size(frame)
    nsteps = length(trajectory)
    println("# There are $(size(frame)) atoms in the frame")
    #positions(frame)
    
    expii=0
    displacement_limit=0

    NumberTimeGaps=n_block   +   blocksize
    NumberSteps=n_block   *   blocksize+1
    
    

    #=
    create time list
    =#
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
    create timegap list for calcaulting time correlation functions
    =#
    TimeGapTable=zeros(Float64,NumberTimeGaps)
    for timeii = 1:blocksize
        TimeGapTable[timeii]=TimeList[timeii]-TimeList[1]
    end

    for timeii = 1:n_block
        TimeGapTable[timeii+blocksize]=TimeList[blocksize*timeii+1]-TimeList[1]
    end

    TimeGapTable[NumberTimeGaps]=TimeList[NumberSteps]-TimeList[1]
    
    #=
    calculations
    =#

    msd = zeros(Float64,NumberTimeGaps)
    weighting = zeros(Float64,NumberTimeGaps)

    p = Progress(NumberTimeGaps)

    for timegapii=0:blocksize-1
        displacement_count=0
        for blockii=1:n_block

            thisii = blocksize*(blockii-1)
            nextii = thisii+timegapii
            frame_this = read_step(trajectory,thisii)
            #matched = evaluate(selection,frame)
            thisii_array = positions(frame_this)

            frame_next = read_step(trajectory,nextii)
            #matched = evaluate(selection,frame)
            nextii_array = positions(frame_next)
            
            delta_array=thisii_array.-nextii_array
            squared_distance=sum(dot(delta_array,delta_array))/natoms
        
            weighting[timegapii+1]+=1
            msd[timegapii+1]+=squared_distance
            

            displacement_count+=1
            if displacement_count == displacement_limit
                break
            end
            next!(p)
        end
    end

    for timegapii=blocksize:NumberTimeGaps-1
        displacement_count=0
        block_timegapii = timegapii - blocksize + 1
        for blockii=1:n_block-block_timegapii
            thisii = blocksize*(blockii-1)+expii
            nextii = thisii + blocksize*block_timegapii
            
            frame_this = read_step(trajectory,thisii)
            #matched = evaluate(selection,frame)
            thisii_array = positions(frame_this)

            frame_next = read_step(trajectory,nextii)
            #matched = evaluate(selection,frame)
            nextii_array = positions(frame_next)

            delta_array=thisii_array.-nextii_array
            squared_distance=sum(dot(delta_array,delta_array))/natoms
            weighting[timegapii+1]+=1
            msd[timegapii+1]+=squared_distance

            displacement_count+=1
            if displacement_count == displacement_limit
                break
            end
            next!(p)
        end
        if displacement_count == displacement_limit
            break
        end
    end

    frame_this = read_step(trajectory,1)
    #matched = evaluate(selection,frame)
    thisii_array = positions(frame_this)

    frame_next = read_step(trajectory,NumberSteps-1)
    #matched = evaluate(selection,frame)
    nextii_array = positions(frame_next)

    delta_array=thisii_array.-nextii_array
    squared_distance=sum(dot(delta_array,delta_array))/natoms

    weighting[NumberTimeGaps]+=1
    msd[NumberTimeGaps]+=squared_distance
    for stepii=1:NumberTimeGaps
        msd[stepii]=msd[stepii]/weighting[stepii]
    end

    io = open(output, "w")
    for stepii=1:NumberTimeGaps
        timegap=@sprintf("%.4f", TimeGapTable[stepii])
        data=@sprintf("%.4f", msd[stepii])
        write(io, "$timegap"*" "*"$data"*"\n")
    end
    close(io)

    return nothing
end
