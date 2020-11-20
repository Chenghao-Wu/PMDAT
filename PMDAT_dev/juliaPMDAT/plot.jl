using DelimitedFiles
using Plots

data=readdlm("PAXNa100.Rg", ' ', Float64)
Rg=data[:,4]+data[:,2]+data[:,3]

plot(data[:,1],Rg)