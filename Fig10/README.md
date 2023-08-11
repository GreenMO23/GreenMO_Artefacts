# GreenMO_Artefacts: Fig 10

## What this Figure is:
The figure shows the SINR and Goodput performace of GreenMO when it is serving 4 interfering users, as number of antennas in the PCB prototype are increased from 4 -> 8. The number of virtual RF chains (created from the expanded bandwidth) is fixed to 4. The key implication is that with more number of antennas, GreenMO serves the 4 users reliably, without increasing power consumption (since the power consumption of adding antennas in GreenMO is just the power in RF switches which is minimal compared to RF chains). A detailed numerical analysis is in Table 1.

## How to reproduce the results in this figure?
Make sure you have the associated data downloaded. The iq-data for this figure is in the `goodput_evals` subfolder. Please remove any other folder from Matlab path and add Fig10 folder as well as subfolders to matlab path. 

Then run `plot_fig10_result.m`. It should generate Fig. 10 of the paper. You can follow along the comments in the matlab script if you want to better understand the code flow.