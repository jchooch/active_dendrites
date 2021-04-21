# Neuron Model with Active Dendrites

This code implements a two-compartment neuron model with active dendritic spike propagation.
This model was described/inspired by:
**[Larkum et al. 2004](https://pubmed.ncbi.nlm.nih.gov/15115747/) â€” Top-down dendritic input increases the gain of layer 5 pyramidal neurons**.

## TODO

- Figure out all the units.
- Figure out what experiments we need to run (parameters, inputs).
- Add time-dependent applied currents.
    - E.g. a stepped, noisy current.
- Plot currents over time.
- Plot fI curves for different applied currents.
- Calculate CV of ISI.
- Plot ISIs over time.
- See figures to plot in slides: 
    https://docs.google.com/presentation/d/1P3cYdprwlWrXsO4hjGdz60px-bmWu7g9qO2HRPGbbUo/edit#slide=id.gd35cab493c_0_1

"When current was 
injected at the soma, the CV was small and decreased further
with increasing current injection. This is consistent with the
interpretation that at low currents, the CV is mainly determined 
by the input current fluctuations (?), while at high input
currents, CV is mainly determined by the mean (µ) of the input
current."

"For distal current injection, the CV was much higher, 
particularly near threshold where the CV reached values >1. 
APs came in bursts, which are mainly determined by the
intrinsic properties of the cell rather than the fluctuations in
the input current."

