# Closed-Loop Neuropathic Pain Treatment (CLNPT)
---
## Releases
* v1.0: Software version that is used for the paper submitted to Scientific Reports (at April 17, 2024)

## Folders
* `FOS_MPC_func`: functions belonging to the FOS_MPC_main.m
* `RhoPrescottModel`: 2D model of Rho and Prescott (based on noise)
* `MapBasedModel`: Map-Based model (phenomenological approach)
* `KovalskyModel`: Kovalsky model (modified Hodgkin-Huxley model)
* `figures_paper`: Exported figures that are in the paper


## Main functions
Application of the FOS-MPC using the different models.

### Using the Rho and Prescott model
* `FOS_MPC_rho_main.m`: Inhibiting subthreshold oscillations and ectopic discharge using FOS_MPC.
* `FOS_MPC_rho_main_constr_pulse.m`: Same as above, but application of a biphasic pulse as constraint.
* `FOS_MPC_rho_main_constr_charge.m`: Same as above, but now the integral of the stimulation pulse should equal 0.

### Using the Map-Based model
* `FOS_MPC_map_main.m`: Inhibiting subthreshold oscillations and ectopic discharge using FOS_MPC.

### Using the Kovalsky model (modified Hodgkin-Huxley model)
* `FOS_MPC_HH_main.m`: Inhibiting subthreshold oscillations and ectopic discharge using FOS_MPC.

