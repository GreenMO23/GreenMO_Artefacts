# GreenMO_Artefacts: Fig 13a, 14

## What these Figures are:

Fig. 13a shows a simulation result with 8 interfering users and upto 64 antennas in a DBF/HBF/GreenMO architecture to observe of GreenMO's perfromance trends scale with more users/antennas. 8 users, 64 antennas is a common test baseline in modern 5G Massive MIMO settings. We see similar performance trends to hardware measurements: "Given same complexity (same number of hardware components/digitized bandwidth), GreenMO outperforms the baselines and comes close to the full digital and fully connected HBF oracle baselines."

Fig. 14 shows how GreenMO breaks the Spectral vs Energy Efficiency (SE vs EE) tradeoff plaguing the wireless networks today. This is because as we improve SE, we add more antennas which comes at a power cost that reduced EE. In GreenMO, this power cost is made minimal, and this helps scale GreenMO architecture to very high levels of SE without observing a dip in EE. In contrast, with power measurements of existing technology, Massive MIMO implemented via digital beamformers tend to become energy inefficient starting 64 antennas mark. 

## How to reproduce the results in this figure?

Please remove any other folder from Matlab path and add Fig13a14 folder as well as it's subfolders to matlab path. Since this is a simulation, it does not need i-q data, however the simulation takes a large amount of time (as it simulates upto 256 antenna channels via ray-tracing), and hence a pre-run data has been included in the `Pre_run_data` folder to help reproduce the results. It may also need a large memory to simulate Fig. 14, and if requested, we can provide access to a server cluster which can run the simulation for the reviewer. Please contact us via hotcrp and follow up on this if needed.

Run `sim_runs/Plot_Fig13a.m` and `sim_runs/Plot_Fig14,m`. These will seprately generate Fig. 13a and Fig. 14 of the paper. You can follow along the comments in the matlab script if you want to better understand the code flow.