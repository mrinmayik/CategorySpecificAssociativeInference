# Category Specific Associative Inference in Memory

This repository accompanies the paper titled *Scene selective medial temporal lobe subregions are recruited for the integration of non-scene stimuli* published in the Journal of Cognitive Neuroscience (doi.org/10.1162/jocn.a.73). The scripts contained here can be used to replicate the figures and results in the manuscript. You can download the scripts to your machine, or clone this repository to your machine.

The dataset accompanying these scripts is available on OpenNeuro: doi.org/10.18112/openneuro.ds006039.v1.0.2

### Information about scripts:
1. Initialise.R: This script sets up paths and functions that used for the analyses. Update the `scripts_path` and `base_path` in this script to the appropriate paths on your machine.
2. AnalyseBehaviouralData.R: Code in this script will generate the behavioural and eye-tracking results from the manuscript. Update the path in the `setwd` command at the top of the script to point to the folder containing all the scripts.
3. AnalyseBehaviouralData.R: Code in this script will generate the analysis of beta estimates from the manuscript. Update the path in the `setwd` command at the top of the script to point to the folder containing all the scripts.
