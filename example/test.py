#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 14:17:39 2019

@author: zwu
"""
import sys
sys.path.append('/home/zwu/Dropbox/code/research')
import SMMSAT

Reader=SMMSAT.LAMMPSReader("/home/zwu/Dropbox/code/research/SMMSAT/example/PE260_58_frames.custom",start=10)
sys=SMMSAT.System(Reader)
sys.set_LinearTimeScheme(48,0.5)
sys.set_Species("polymer",33,[1],[260])
sys.set_Species("test",20,[1],[260])
app=SMMSAT.Application(sys)

#test all kinds of atom list
all_list=SMMSAT.AtomList(sys,"all")
type_system_list=SMMSAT.AtomList(sys,"type_system",1)
type_species_list=SMMSAT.AtomList(sys,"type_species","polymer",1)
species_list=SMMSAT.AtomList(sys,"species","polymer")
indexatom_list=SMMSAT.AtomList(sys,"index_atom","polymer",np.arange(120,140))


#test multibody list 
ete_list=SMMSAT.MultiBodyList(sys)
ete_list.create_MultiBodyList("ete",1,"species_atomlist","polymer",[1,0,1,259])
ete_list.combine_multibody_lists("ete")

#test center of mass chain list
chain_list=SMMSAT.MultiBodyList(sys)
chain_list.create_MultiBodyList("chain",1,"species_atomlist","polymer",[1,0,1,1,1,2,1,3,1,4,1,5,1,6,1,259])
chain_list.combine_multibody_lists("chain")

#test end end vector autocorrelation function
ete_list=SMMSAT.MultiBodyList(sys)
ete_list.create_MultiBodyList("chain",1,"species_atomlist","polymer",[1,0,1,259])
ete_list.combine_multibody_lists("chain")
eeacf=SMMSAT.vector_autocorrelation_function(ete_list,"test_eeacf","xyz")
app.add(eeacf)

#test end to end distance and distribution of it 
ete=SMMSAT.end_end_distance(type_species_list,"test_ete",calc_Dist=True)
app.add(ete)

#test mean squared internal distance
msid=SMMSAT.mean_squared_internal_distance(type_species_list,"test_msid")
app.add(msid)

#test gyration tensor
rg=SMMSAT.gyration_tensor(type_species_list,"test_rg")
app.add(rg)

#test radial distribution function
rdf=SMMSAT.radial_distribution_function(all_list,"test_rdf",500,20)
app.add(rdf)

#test mean squared displacement for inner beads
msd_inner=SMMSAT.mean_squared_displacement(indexatom_list,"test_msd_inner","xyz")
app.add(msd_inner)

#test mean squared displacement for all beads
msd=SMMSAT.mean_squared_displacement(all_list,"test_msd","xyz")
app.add(msd)

#test self-part intermediate scattering function
isfs=SMMSAT.intermediate_scattering_function(all_list,"test_isfs",26,"xyz",20)
app.add(isfs)

app.run()