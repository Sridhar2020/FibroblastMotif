# FibroblastMotif
Code to generate data for the different Myocyte-FibroblastMotifs

Folders Steep and Shallow contain sub-folders for each of the three coupled Fibroblast-Myocyte motifs for parameters Steep and Shallow in the TP06 model.
(Note that some of the files are common across the different motifs).

For each motif, the programs simulates 4 pacing periods (T=300,400,500,600), 3 Delays ( tau = 0, 10, 25 ms) and 9 values of G_loc (strong) and G_Long (weak) conductance values.

To run the programs first run MakeFolder.py followed by Run_script.py

The programs here are for the case where the Fibroblasts have a resting membrane potential of -24.5 mV. The resting membrane potential of the fibroblast can be changed bys hifting the gating variable voltage dependence of the time dependent potassium current as described in Jacquemet et al (Am. J. Physiol. Circ. Physiol. 294, H2040–H2052 ,2008).




  
