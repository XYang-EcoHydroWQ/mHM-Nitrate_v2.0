## mHM-Nitrate_v2.0 ##
**The mHM-Nitrate model** is a fully distributed, flexibly designed nitrate model at catchment scale. The model is developed based on the hydrological platform of mesoscale Hydrological Model (mHM, v5.5) (https://www.ufz.de/mhm). Descriptions of nitrate transport and removal processes are mainly introduced from Hydrological predictions for the Environment (HYPE) (http://www.smhi.net/hype/wiki/doku.php?id=start), with several improvements for a better model representation under changing anthropogenic condtions.

Compared to the first version (https://git.ufz.de/yangx/mHM-Nitrate), this version (**v2.0**) has three main modifications: (1) implemented a new regionalization approach for in-stream autotraphic nitrate uptake (Yang et al., Water Res. 2019); (2) improved the crop growth descriptions for winter crops (i.e., witner crops only take up considerable nutrients after the "emerge date"); and (3) added a new input information (i.e., "geoformation.txt" in ./input/water_quality, specifying the permeability of each geologic unit) for a better representation of initial soil nitrate conditions. Detailed information has been provided in **Documentation_WQinputs.pdf** and has been documented within the source code in specific locations.

**This repository contains the nitrate module source code of the mHM-Nitrate model (v2.0)**. Brief instructions of how to compile/run the model are given below (**Steps to run the model**). The case study data of the Selke catchment (456 km2) in central Germany, are provided in https://git.ufz.de/yangx/mHM-Nitrate. 

The presented source code is developed based on the mHM version 5.5. The version based on mHM v5.7 is also available under request to the author Xiaoqiang Yang (xiaoqiang.yang@ufz.de) and Michael Rode (michael.rode@ufz.de).

## Citations ##
Please also refer to the following publications:

Yang, X., Jomaa, S., Zink, M., Fleckenstein, J. H., Borchardt, D., & Rode, M. (2018), A new fully distributed model of nitrate transport and removal at catchment scale, Water Resources Research, 54(8), 5856-5877. https://doi.org/10.1029/2017WR022380.

Yang, X., Jomaa, S., Büttner, O., & Rode, M. (2019), Autotrophic nitrate uptake in river networks: A modeling approach using continuous high-frequency data, Water Research, 157, 258-268. https://doi.org/10.1016/j.watres.2019.02.059.

## Steps to run the model ##
1. Download the mHM (v5.5) hydrological code under the link:https://github.com/mhm-ufz/mhm/releases. 
2. Download the water quality module code and configuration files from this repository.
3. Combine those files and replace the original mHM files where needed*.
4. Prepare the input data:
   - Follow steps in the mHM manual for hydrological inputs (provided associated with the mHM source code);
   - "Documentation_WQinputs.pdf" for nitrate inputs and initial conditions.
5. Set **processCase(11) = 1** in mhm.nml to activate the nitrate simulation.
6. Run the model with "./mhm-nitrate" (or compile again using the "make" command).

*Please note that here we provide several code/configuration files of the mHM v5.5 (./src/mHM/* and *.nml). To run the nitrate simulation, these files have to be modified and replaced.   