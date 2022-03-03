# code-for-gains-calculator
This repository accompanies the upcoming papers to be submitted to ACP regarding the intercomparison between different stratospheric SO₂ injection latitudes 
in different climate models.
References will be added when the paper is available for discussion.

The file control_gains_calculator.m calculates is given in input zonal mean values for stratospheric AOD and surface temperature 
for a control simulation (no SO2 injection) and 4 simulations considering an injection of SO2 in the stratosphere at four different latitudes (30N, 15N, 15S and 30S.

It then maps the changes in T and AOD on the three Legendere polynomials (T0, T1 and T2 for T and L0, L1 and L2 for AOD). 
It calculates the sensitivity of the given climate model (i.e., how much a change in l0 causes a change in T0, and so on), and determines how T0, T1, and T2 will 
change over time in the background run (how much temperature change the SO2 injections will need to offset).

Then it computes the feedforward gains needed to offset the expected changes in temperature computed in the previous section.
Finally, it computes the matrix "M" for the model, which relates injections at each of the individual locations to changes in l0, l1, and l2.

The code also contains a function to calculate annual means.

Further information is given in the paper.
