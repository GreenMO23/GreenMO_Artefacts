# GreenMO_Artefacts
Artefact submission for Mobicom'23 paper #121, "GreenMO". Github repository link: [https://github.com/GreenMO23/GreenMO_Artefacts](https://github.com/GreenMO23/GreenMO_Artefacts)

## Overview and General Instructions
We have included the data, supporting data analysis and plotting codes for each of the following result figures (Referenced to the accepted paper included with the submission): 

- **Figs 9c,9d**: Trace Level SINR Comparisons of GreenMO with DBF and HBF (Digital and Hybrid Beamforming)
- **Figs 10a,b**: Goodput and SINR performance of GreenMO, how GreenMO guarantees robust goodput and improved SINRs by leveraging more antennas
- **Fig 11a**: Comparison of achieved capacity by GreenMO compared to interference-free FDMA baseline 
- **Fig 11b**: 4 antenna GreenMO's goodput results with 2, 3, 4 users
- **Fig 12b**: Ablation study, of GreenMO's virtual RF chains versus DBF's physical RF chain
- **Fig 12c**: Ablation study, GreenMO's BABF antenna selection approach vs random antenna configs
- **Fig 12d**: Ablation study, Synchronized vs Unsynchronized users and the impact on SINRs
- **Fig 13a**: 64 antenna, 8 users simulation of GreenMO and comaprison with different HBF, DBF architectures
- **Fig 14**: Simulated SE and EE comparison of GreenMO with today's Massive MIMO base stations

The codes for each of these are in separate folders of the [github repo](https://github.com/GreenMO23/GreenMO_Artefacts), and have their own readme's detailing how to run the codes. Please refer to the per-folder readmes for more details. To run results for a particular figure, you need to add the folder to matlab path, **!!and!!** remove other folders included beforehand in the path. Some of the figures share simialr named codes but used in different ways, and to keep simplicity and make the artefact evaluation process easier, we have demarcated different codes in different folders, and these need to be added to path to generate a particular result.

## Data and Required Software
In order to reproduce the results in the paper (Figs 9-14), you would need to download the radio traces collected by us (iq-data) from the GreenMO PCB and the test setup of WARP/USRP SDRs as detailed in the paper. Total space required for the data is **20.4 GB**. 

Please use ["FDT"](http://monalisa.cern.ch/FDT/download.html) to make the data transfer faster (usually takes about 5 minutes). You would need to setup Java JDK on your machine to run FDT (usually pre-installed on linux, for windows follow this [link](https://www.oracle.com/java/technologies/downloads/#jdk20-windows)). Open the command window and navigate to location where you have the "fdt.jar" file and run the following command to get the data required for artefact evaluation of GreenMO:


>  java -jar fdt.jar -p 54323 -pull -r -c 137.110.115.35 -r -d ./ /home/comsol_account/GreenMO_IQ_Data

Alternately, you can download from google drive: [link](https://drive.google.com/file/d/1eQUOTNBp6XQSZd0G7gYAEmD84t4tqREl/view?usp=sharing). There might be some issues with local firewalls which may prevent binding to the FDT port. However, you would need to be patient, since the zip file downloaded from google is about 18 GB and depending on your internet connection it can take anywhere between 30-60 mins to download the zip file.

To run the codes and generate the plots from this downloaded data, you would require a standard Matlab installation, with no requirement of a third party library.

