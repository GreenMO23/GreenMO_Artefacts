# GreenMO_Artefacts: Fig 9c, 9d

## What these Figures are:

Fig. 9c/9d shows SINR comaprisons of GreenMO (8 antenna) versus Digital and Hybrid beamforming baselines. Since it is difficult to create different architectures, we do these comparisons via trace level emulations, where we collect channels of each individual antenna, and then emulate performance of digital beamformer by implementing digital combining atop these channels, and for hybrid beamformers, we map the channels to physical RF chains by assuming ideal phase shifter analog network. We see that given same complexity (same number of hardware components/digitized bandwidth), GreenMO outperforms the baselines and comes close to the full digital and fully connected HBF oracle baselines.

## How to reproduce the results in this figure?
Make sure you have the associated data downloaded. The iq-data for this figure is in the `trace_evals` subfolder. Please remove any other folder from Matlab path and add Fig9c9d folder as well as it's subfolders to matlab path. 

Then run `plot_fig9.m` with `config='dbf'` and `config='hbf'`. These will seprately generate Fig. 9c and Fig. 9d of the paper. You can follow along the comments in the matlab script if you want to better understand the code flow.