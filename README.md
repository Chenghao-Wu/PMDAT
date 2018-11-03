# Soft Matter Molecular Simulation Analysis Toolkit(SMMSAT)

## Overview
### Capabilities
SMMSAT is a toolkit of analysis methods for molecular simulation of Soft Matter(Polymer, liquid etc.). Currently, it can read file format such as .lammpstrj, .xyz, .xml, the format of trajectory file may affect commands of analysis method in input file. The most stable file format at present is .lammpstrj output by LAMMPS CUSTOM command ( with pbc image). SMMSAT enables selection of user-defined sets of particles, which is called " multibody" in SMMSAT and therefore analysis of properties of the multibody's center of mass such static structure, vector correlation, center of mass motion etc. SMMSAT is presently in intensive developing and updating within adding new analysis methods, speed optimization, and stability testing. If you met any questions or bugs or have additional functionality requests, please contact Zhenghao w415146142@gmail.com.
### General Concepts and Terminology
SMMSAT is written in python and accelerated by cython. SMMSAT can be run in terminal and also in ipython(jupyter notebook) as normal python package. To run analysis, the user need to provied some metadata of the system and scripts specifying the analysis method( see Section "Usage"). 

## Prerequisites
* It is fully tested within Anaconda5.3(python 3.7)
* At least: cython, numpy, pandas, pathlib.

## Installing
1. pip install SMMSAT or git clone https://github.com/Chenghao-Wu/SMMSAT.git
2. compile cython code, open terminal and run commands:
```bash
cd ./SMMSAT/cython_func
python setup.py build_ext --inplace
```

## Usage
Scripts ro run SMMSAT are always divided into three parts: 1. Import Package; 2. System Block; 3. Analysis Block
### Import Package
```python
import sys
sys.path.append('.../directory_of_SMMSAT/')
import SMMSAT
```
### System Block
1. < trajectory file reader >
2. < timescheme > (see Subsection "TimeScheme")
3. < additional lines required by trajectory format > (see Subsection "Trajectory File Reader")

Example:
```python
filename = "trajectory.lammpstrj"
DataBuild=SMMSAT.LAMMPSReader(filename)
sys=SMMSAT.System(DataBuild)
sys.set_ExponentialTimeScheme(10,38,1.2,0.005) #number of block: 10; block size: 38; exponent base: 1.2; timestep: 0.005
sys.set_Species("polymer",1687,[1],[20]) # name of speicies: polymer; number of Species( polymer): 1687; 
#atom type in Species( polymer): 1; number of atom 1: 20 in each Species( polymer) 
app=SMMSAT.Application(sys)
```
#### Trajectory File Reader
Currently support: .lammpstrj; .xyz; .xml
1. .lammpstrj
* SMMSAT is presently able to recognize the format of custom file output by LAMMPS CUSTOM command. Custom styles understood by SMMSAT includes wrapped position(x,y and z), unwrapped position(xu,yu,zu),image of pbc box(ix, iy and iz), velocity(v_x, v_y and v_z). At least one type of position needs to be input, and image of pbc box and velocity is according to analysis method. The most common custom format is which includes the sorted of atom index and the image of periodic boundary. The user can output this format by LAMMPS command: 
    ```
    dump ID group-ID custom N trajectory.custom type x y z ix iy iz
    dump_modify commandid first no sort id
    ```

* Additional lines in system block for custom:
    * first species ```system.set_Species(user-defined species name 1,number of species,[ type ],[ number of atoms for each type ])```
    * second species  ```system.set_Species(user-defined species name 2,number of species,[ type ],[ number of atoms for each type ])```

2. .xyz
* xyz is a common trajectory file format produced by simulation packages such as LAMMPS. It is a text file containing the positions of each particle at each time frame sequentially. Each frame is headed by two lines. The first reads “atoms” and the second contains the number of atoms in that frame. These are followed by a list of all atoms, where the first column is the atom type index, and the second, third, and fourth are the x, y, and z coordinates of that atom at that time step. Note that SMMSAT assumes that the atoms are in the same order in each frame; if the trajectory does not conform to this rule results will generally be incorrect. Moreover, because of lack of images of periodic boundary, the dynamical properties calculated from .xyz trajectory are not reliable and only the static structure properties are able to be correctly analyzed.
* Additional lines in system block for xyz:
    * first species ```system.set_Species(user-defined species name 1,number of species,[ type ],[ number of atoms for each type ])```
    * second species  ```system.set_Species(user-defined species name 2,number of species,[ type ],[ number of atoms for each type ])```

3. .xml
* will be added in the future

#### TimeScheme
Currently support: linear timeschmeme, exponential timeschmeme
1. Linear Timescheme
* format: 
```system.set_LinearTimeScheme(<number_frames>, <timestep_between_frames>)```
* linear timescheme is common in molecular simulation especially for analyzing static or struture preporties

2. Exponential Timescheme
* format:
```system.set_ExponentialTimeScheme(<number_blocks>,<block_size>,<>exponential base>,<time unit>)```
* exponential timeschmeme is suitable for analyzing dynamics of polymer which is usually exponential developing with time, and pseudocode are given in following fomula:

        timelist=[]*(block_size*number_block+1)
        block_starttime=0;
        for(blockii=0;blockii<number_exponentials;blockii++)
        {
            for(expii=1;expii<=block_size;expii++)
            {
                timeii++;
                if(power(exponential_base,expii-1) <= expii)
                {
                    timelist[timeii] = block_starttime+expii*time_unit;
                }
                else
                {
                    timelist[timeii] = block_starttime+floor(pow(exponential_base,expii-1))*time_unit;
                }
            }
            block_starttime = timelist[timeii];
        }
        }
    Also there is an example LAMMPS input file yielding a trajectory file corresponding to this time scheme in ./SMMSAT/example.

### Analysis Block
The third part of scripts is analysis block which is to specify the particle set and the analysis method on it.
#### Selecting particle sets for analysis
1. CreateList
* create a list of particles for analysis based on the input parameter such as species, atom type and location within species
    ```
    ListName=SMMSAT.CreateList(system)
    ListName.create_list( <keyword>, <arguments>)
    ```
        keyword         arguments
        
        all             
        * select all atoms in the system
        
        type_system     <atom_type>
        * select atoms of <atom type> in the system
        
        type_species    <species_name> <atom_type>
        * select atoms of <atom type> in <species_name> in the system
        
        species         <species_name>
        * select all atoms in <species_name> in the system
2. Multibody_List
* create a list of "multibody"(particle set) for analysis such as center of mass motion and struture etc according to input parameters
    ```
    MultibodyListName=SMMSAT.MultiBody_List(system)
    MultibodyListName.create_MultiBodyList(<name_multibody>,<type>,<multibodymethod>,<keyword>,<argument>,<atom_list>)
    ```
    < name_multibody >: user-defined name for this multibody
    < multibodymethod >: currently supports only "centroid"( calculate center of mass and vector of the particle's set)

        keyword         arguments
        
        species_atomlist  <species_name>            
        * select atoms of <species_name> to form the multibody
    < atom_list >: [atomtype, atom index,atomtype, atom index...]

#### Analyzing trajectories of particle sets
Every trajectory analysis method in SMMSAT needs a List of Particle sets and specific arguments for this method. 
1. mean_squared_displacement
* calculate the mean squared displacement in 3 dimensions
* msd = SMMSAT.mean_squared_displacement(< ListName >, < FileName >)
2. intermediate_sacttering_function
* calculate the self-part intermediatescattering function
* isfs = SMMSAT.intermediate_sacttering_function(< ListName >, < FileName >, < WaveVectorIndex >, < AnalysisGeometry >, < MaxLengthScale >, < fullblock=0 >)
    * < WaveVectorIndex >: specify the wave vectors to be calculated
    * < AnalysisGeometry >: xyz, xy, xz, yz, x, y, and z. This chooses which dimensions in k-space to include in the calculation of the intermediate scattering function. xyz computes the full radial three dimensional isf, xy, yz, and xz calculate two-dimensional in-plane radial isfs, and x, y, and z compute one-dimensional isfs.
    * < MaxLengthScale >: determines the longest distance which will be decomposed into inverse space. If a given distance is larger than smallest box size, the smallest box size will be taken as MaxLengthScale
    * < fullblock >: an optional argument that can be either 0 or 1. The default is 0, in which case time spacings spanning multiple blocks use only the first time of each block. A setting of 1 specifies that it should use all block times for times spanning blocks. This may result in substantial computational time increases but offers the possibility of modestly increased data quality at very long times.
3. bond_autocorrelation_function
* calculate the bond orientational autocorrelation function
* baf = SMMSAT.bond_autocorrelation_function(< ListName >, < FileName >, < AnalysisGeometry >)
    * < AnalysisGeometry >: xyz, xy, xz, yz. This chooses which dimensions in real space to include in the calculation of the bond autocorrelation function. xyz computes the full radial three dimensional baf, xy, yz, and xz calculate two-dimensional in-plane radial baf.
4. radial_distribution_function
* calculate the radial distribution function
* rdf = SMMSAT.radial_distribution_function(< ListName >, < FileName >, < NumberBins >, < MaxLengthScale >)
    * < NumberBins >: number of bins to be used to do histogram calculation
    * < MaxLengthScale >: maximum length scale for rdf
5. gyration_tensor
* calculate gyration tensor
* gyration_tensor = SMMSAT.gyration_tensor(< ListName >, < FileName >)
    * < ListName >: here < List > should be the molecule species list e.g.: 
    ```python
    List_gyrationtensor=SMMSAT.CreateList(sys)
    List_gyrationtensor.create_list("species","polymer")
    ```
6. end_end_distance
* calculate end to end distance
* end_end_distance = SMMSAT.end_end_distance(< ListName >, < FileName >)
    * < ListName >: here < List > should be the multibody list of chain end atoms e.g.:
    ```python
    ChainEnds=SMMSAT.MultiBody_List(sys)
    ChainEnds.create_MultiBodyList("chainend",1,"centroid","species_atomlist","polymer",["S",0,"S",159])
    ChainEnds.combine_multibody_lists("chainend")
    ```
7. order_parameter
* calculate structure orientation order paramter 
* order_parameter = SMMSAT.order_parameter(< ListName >, < FileName >, < AnalysisGeometry >)
    *  < ListName >: here < List > should be the multibody list of atoms within species to obtain the vector list
    *  < AnalysisGeometry >: x, y and z. This chooses which direction in real space to calculate the relation orientation order parameter. x computes $\vec{r} * (1,0,0)$, y computes $\vec{r} * (0,1,0)$ and z computes $\vec{r} * (0,0,1)$