# GreenMO_Artefacts: Fig 12c

## What this Figure is:

Fig. 12c shows how GreenMO's BABF method of antenna selection outperforms a on-off configuration chosen randomly. We collect 120 data points with 4 interfering users having fixed position supported with random configurations and 10 data points with BABF (since BABF is deterministic algorithm it does not need a lot of trials to average out). We see BABF significantly outperforming the random choices since BABF tries to prioritize one user per virtual RF chain that helps the MIMO processing atop the chains. 

## How to reproduce the results in this figure?
Make sure you have the associated data downloaded. The iq-data for this figure is in the `random_antenna_configs_vs_babf` subfolder. Please remove any other folder from Matlab path and add Fig12c folder as well as it's subfolders to matlab path. 

Then run `plot_fig_12c.m`. This generate Fig. 12c of the paper. You can follow along the comments in the matlab script if you want to better understand the code flow.