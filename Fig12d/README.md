# GreenMO_Artefacts: Fig 12d

## What this Figure is:

Fig. 12d shows that since GreenMO performs a spatial interference cancellation procedure, it is agnostic of the 4 interfering users being synchronized or not. For this experiment, we collect 50 i-q traces with the users being synchronized and unsynchronized, and plot the result in a CDF to see the performance is largely identical.

## How to reproduce the results in this figure?
Make sure you have the associated data downloaded. The iq-data for this figure is in the `synch_vs_unsynch` subfolder. Please remove any other folder from Matlab path and add Fig12d folder as well as it's subfolders to matlab path. 

Then run `plot_fig12d.m`. This generate Fig. 12d of the paper. You can follow along the comments in the matlab script if you want to better understand the code flow.