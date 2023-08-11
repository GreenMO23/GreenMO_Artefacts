# GreenMO_Artefacts: Fig 11b, 12b

## What these Figures are:

Fig. 11b shows goodput evaluations with 4 antenna version of GreenMO, to show how GreenMO can still improve the goodput while being power efficient when number of antennas reduce. Commercial devices like smartphones/laptops can furnish 4 antennas to facilittate these results.

Fig. 12b shows one to one comparison with GreenMO's virtual RF chains and standard digital beamformer's physical RF chains. To enable a fair comparison, we fix one antenna per RF chain in case of virtual RF chains (even though it has a many to many connection with number of antennas), since a physical RF chain connects only to a single antenna. Fig. 12b shows that virtual RF chains perform with about the same SINR as a physical RF chain.

## How to reproduce the results in this figure?
Make sure you have the associated data downloaded. The iq-data for this figure is in the `PCB_MMSE_virt_vs_physical` subfolder. Please remove any other folder from Matlab path and add Fig11b12b folder as well as it's subfolders to matlab path. 

Then run `plot_fig11b.m` and `plot_fig12b.m` (they share the same iq-data files). These will seprately generate Fig. 11b and Fig. 12b of the paper. You can follow along the comments in the matlab script if you want to better understand the code flow.