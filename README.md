# EuclideanDecoders

MATLAB functions to decode stimuli from extracellular recordings of neuronal population responses. Decoders based on Euclidean Distance of vectorized responses.

For single neurons, assumes input is spiketimes of neurons in the format: data {modulation frequency}{trial} For populations: data{neuron}{modulation frequency}{trial}

Differences between “full spiketrain,” “rate-only,” and “phase-only” decoders described in: 
Malone et al., 2015. Diverse cortical codes for scene segmentation in primate auditory cortex. J Neurophys. 113(7):2934-52.

Summary: each trial is represented as a vector, which can be the activity of a single neuron over time, or some combination or concatenation of a population response. For each modulation frequency, the average response over trials is calculated. Then for each trial, the euclidean distance is calculated between that trial and the mean of each modulation frequency. As a method of cross-validation, the trial being decoded is removed from the trial-average. The decoding "decision" is taken as the modulation frequency that has the shortest euclidean distance from the trial being decoded. Confusion matricies are constructed to compare the Actual vs Decoded modulation frequency.

Written by Melissa Jane Runfeldt, Spring 2015, in the lab of Brian Malone at UCSF
