# Beta-Gamma Simple Monte Carlo Simulation

This repo will show how to implement Monte Carlo particle simulation for a simple case. In this repo, the code can only be used to perform Monte Carlo simulation inside a cylindrical medium resembling the Liquid Scintillator, as been shown here:

<p align="center">
  <img src="https://github.com/Feryalchemist/Beta-Gamma_MonteCarlo/blob/main/Problem_Description.png?raw=true" alt="Monte Carlo Problem Description"/>
</p>

We assign the upper lid as the detector or the place in which passed radiations will be counted. In this system, we assume that the detector is very efficient and the limitation of detection only comes from the geometry and beta mean free path. The beta mean free path inside this medium is assumed by user input to simplify the simulation of charged particle simulation. 

[H2O_CS](https://github.com/Feryalchemist/Beta-Gamma_MonteCarlo/blob/main/H2O_CS) file is an example that was obtained from this [database](https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html). It can be changed by other material by defining the file path of the new table.
