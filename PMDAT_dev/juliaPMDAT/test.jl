using Chemfiles
include("./PMDAT.jl");using .PMDAT

trj=Trajectory("../trajectory_PAXNa100_000.gen.lammpstrj")

#@time msd(trj,output="test.msd",n_block=10,blocksize=90,timeunit=0.03,exponentialbase=1.2)
Rg(trj,output="PAXNa100.Rg",selection="",n_block=1,blocksize=63,timeunit=0.03,exponentialbase=1.2,num_molecule=100)