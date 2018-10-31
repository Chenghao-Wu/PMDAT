#!/usr/bin/env python3

import sys
sys.path.append('/home/zwu/Dropbox/code/research/')
import SMMSAT

filename = "/home/zwu/trajectory_amdat.custom"
#LogFilename="/home/zwu/50PS160.xml"
DataBuild=SMMSAT.LAMMPSReader(filename)
#DataBuild.set_ReadLog("mass")
sys=SMMSAT.System(DataBuild)
sys.set_ExponentialTimeScheme(10,38,1.2,0.005)
sys.set_Species("polymer",1687,[1],[20]) # follow the order of trajectory file
app=SMMSAT.Application(sys)

#multibody=SMMSAT.MultiBody_List(sys)
#multibody.create_MultiBodyList("1",1,"centroid","species_atomlist","polymer",["S",0,"R",1,"S",2,"S",3,"S",4,"R",5,"S",6])
#multibody.combine_multibody_lists("1")

#msd_multibody=SMMSAT.mean_squared_displacement(multibody,"msd_multibody.dat")
#app.add(msd_multibody)
#isfs_multibody=SMMSAT.intermediate_scattering_function(multibody,"isfs_multibody.dat",26,"xyz",12.4385)
#app.add(isfs_multibody)
#rdf_multibody=SMMSAT.radial_distribution_function(multibody,"rdf_multibody.dat",1000,100)
#app.add(rdf_multibody)

#ChainEnd=SMMSAT.MultiBody_List(sys)
#ChainEnd.create_MultiBodyList("chainend",1,"centroid","species_atomlist","polymer",["S",0,"S",159])
#ChainEnd.combine_multibody_lists("chainend")
#EndEndDsitance=SMMSAT.end_end_distance(ChainEnd,"end_end_distance")
#app.add(EndEndDsitance)

#ChainEnd=SMMSAT.MultiBody_List(sys)
#ChainEnd.create_MultiBodyList("chainend",1,"centroid","species_atomlist","polymer",["S",0,"S",159])
#ChainEnd.combine_multibody_lists("chainend")
#ChainEndVectorAutocorr=SMMSAT.vector_autocorrelation_function(ChainEnd,"chain_end_vector_autocorrelation.dat","xyz")
#app.add(ChainEndVectorAutocorr)

#BafList=SMMSAT.MultiBody_List(sys)
#BafList.create_MultiBodyList("bond1",1,"centroid","species_atomlist","polymer",[1,0,1,1])
#BafList.create_MultiBodyList("bond2",1,"centroid","species_atomlist","polymer",[1,1,1,2])
#BafList.create_MultiBodyList("bond3",1,"centroid","species_atomlist","polymer",[1,2,1,3])
#BafList.create_MultiBodyList("bond4",1,"centroid","species_atomlist","polymer",[1,3,1,4])
#BafList.create_MultiBodyList("bond5",1,"centroid","species_atomlist","polymer",[1,4,1,5])
#BafList.create_MultiBodyList("bond6",1,"centroid","species_atomlist","polymer",[1,5,1,6])
#BafList.create_MultiBodyList("bond7",1,"centroid","species_atomlist","polymer",[1,6,1,7])
#BafList.create_MultiBodyList("bond8",1,"centroid","species_atomlist","polymer",[1,7,1,8])
#BafList.create_MultiBodyList("bond9",1,"centroid","species_atomlist","polymer",[1,8,1,9])
#BafList.create_MultiBodyList("bond10",1,"centroid","species_atomlist","polymer",[1,9,1,10])
#BafList.create_MultiBodyList("bond11",1,"centroid","species_atomlist","polymer",[1,11,1,12])
#BafList.create_MultiBodyList("bond12",1,"centroid","species_atomlist","polymer",[1,12,1,13])
#BafList.create_MultiBodyList("bond13",1,"centroid","species_atomlist","polymer",[1,13,1,14])
#BafList.create_MultiBodyList("bond14",1,"centroid","species_atomlist","polymer",[1,14,1,15])
#BafList.create_MultiBodyList("bond15",1,"centroid","species_atomlist","polymer",[1,15,1,16])
#BafList.create_MultiBodyList("bond16",1,"centroid","species_atomlist","polymer",[1,16,1,17])
#BafList.create_MultiBodyList("bond17",1,"centroid","species_atomlist","polymer",[1,17,1,18])
#BafList.create_MultiBodyList("bond18",1,"centroid","species_atomlist","polymer",[1,18,1,19])
#BafList.combine_multibody_lists("bond1","bond2","bond3","bond4","bond5","bond6","bond7","bond8","bond9","bond10","bond11","bond12","bond13","bond14","bond15","bond16","bond17","bond18")
#baf=SMMSAT.bond_autocorrelation_function(BafList,"bond_autocorrelation_function.dat","xyz")
#app.add(baf)

OrderParameterList=SMMSAT.MultiBody_List(sys)
OrderParameterList.create_MultiBodyList("bond1",1,"centroid","species_atomlist","polymer",[1,0,1,19])
OrderParameterList.combine_multibody_lists("bond1")
OrderParameter=SMMSAT.order_parameter(OrderParameterList,"order_parameter.dat","z",100)
app.add(OrderParameter)

#List_gyrationtensor=SMMSAT.CreateList(sys)
#List_gyrationtensor.create_list("species","polymer")
#GyrationTensor=SMMSAT.gyration_tensor(List_gyrationtensor,"gyration_tensor.dat")
#app.add(GyrationTensor)

#List_=SMMSAT.CreateList(sys)
#List_.create_list("all")

#msd=SMMSAT.mean_squared_displacement(multibody,"msd.dat")
#app.add(msd)
#isfs=SMMSAT.intermediate_scattering_function(List_,"isfs.dat",26,"xyz",12.4385)
#app.add(isfs)
#rdf=SMMSAT.radial_distribution_function(List_,"rdf.dat",1000,100)
#app.add(rdf)

app.run()