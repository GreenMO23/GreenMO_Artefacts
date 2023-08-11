# GreenMO_Artefacts: Fig 11a

## What this Figure is:
The figure shows the network capacity as number of antennas in GreenMO are increased from 4-> 8, while serving 4 interfering users in same frequency band, as compared to single antenna FDMA baseline which serves 4 users in non-interfering frequency bands. We see that as number of antennas increase, GreenMO starts approaching the network capacity of FDMA baseline, and sometimes even gets better capacity than FDMA, since it does use more number of antennas.

## How to reproduce the results in this figure?
Make sure you have the associated data downloaded. The iq-data for this figure is in the `capacity_result` subfolder. Please remove any other folder from Matlab path and add Fig11a folder as well as subfolders to matlab path. 

Then run `plot_fig11a_result.m`. It should generate Fig. 11a of the paper. You can follow along the comments in the matlab script if you want to better understand the code flow.