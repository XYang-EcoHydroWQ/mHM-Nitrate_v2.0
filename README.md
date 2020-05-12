## mHM-Nitrate_v2.0 ##
The mHM-Nitrate model (v2.0) is a fully distributed, flexibly designed nitrate model at catchment scale. The model is developed based on the hydrological platform of mesoscale Hydrological Model (mHM, v5.5) (https://www.ufz.de/mhm). Descriptions of nitrate transport and removal processes are mainly introduced from Hydrological predictions for the Environment (HYPE) (http://www.smhi.net/hype/wiki/doku.php?id=start), with several improvements for a better model representation under changing anthropogenic condtions.

This repository contains the nitrate module source code of the mHM-Nitrate model (v2.0). Compared to the first version (https://git.ufz.de/yangx/mHM-Nitrate), modifications are mainly: (1) implemented a new regionalization approach for in-stream autotraphic nitrate uptake (Yang et al., Water Res. 2019); (2) improved the crop growth descriptions for winter crops (i.e., witner crops only take up considerable nutrients after the emerge date). Detailed information has been documented within the source code in specific locations.

The presented source code is developed based on the mHM version 5.5 (the version based on mHM v5.7 is also available under request to the author Xiaoqiang Yang (xiaoqiang.yang@ufz.de)).

## Citations ##
Please refer to the following publications:

Yang, X., Jomaa, S., Büttner, O., & Rode, M. (2019), Autotrophic nitrate uptake in river networks: A modeling approach using continuous high-frequency data, Water Research, 157, 258-268. https://doi.org/10.1016/j.watres.2019.02.059.

Yang, X., Jomaa, S., Zink, M., Fleckenstein, J. H., Borchardt, D., & Rode, M. (2018), A new fully distributed model of nitrate transport and removal at catchment scale, Water Resources Research, 54(8), 5856-5877. https://doi.org/10.1029/2017WR022380.

## Steps to run the model ##
1. Download the mHM (v5.5) hydrological code  under the link:https://github.com/mhm-ufz/mhm/releases. 
2. Download the water quality module code and configuration files from here.
3. Combine those files and replace the original mHM files where needed.
4. Prepare the input data:
   - mHM manual for hydrological inputs ()
   - "WQM-InputDescription.pdf" for nitrate inputs).
4. Run the model with "./mhm-nitrate" (or compile again using "make" command).
