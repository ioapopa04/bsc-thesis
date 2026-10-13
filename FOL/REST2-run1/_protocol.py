# md-tools-helper-sha256:93a2a752bc7dbd7fb6b59d0b98aaaefec9b883ea6f755e7e63c661eb9d6342d9
#!/usr/bin/env python
"""REST2: 8 states, tau 0.0 to 0.4, NVT at 300.0 K.

Written by the generated REST2.py at run time. Every state is thermostatted at the SAME
temperature and differs only by Hamiltonian: this is Hamiltonian scaling, not temperature REMD,
and an exchange never rescales velocities.

    tau ladder : [0.0, 0.057143, 0.114286, 0.171429, 0.228571, 0.285714, 0.342857, 0.4]
    (1-tau)^2 on solute-solute terms, (1-tau) on solute-environment terms; amide omega,
    aromatic ring, double bond and improper torsions left unscaled
"""
from md_tools.remd import REST2Protocol

protocol = REST2Protocol(
    tau=[0.0, 0.057143, 0.114286, 0.171429, 0.228571, 0.285714, 0.342857, 0.4],
    temperature_k=300.0,
    timestep_fs=2.0,
    exchange_interval_ps=2.0,
    whole_output_interval_ps=100.0,
    solute_output_interval_ps=2.0,
    checkpoint_interval_ps=100.0,
    number_of_exchanges=25000,
    friction_per_ps=1.0,
    equilibration_ps=100.0,
    random_seed=1,
    hydrogen_mass_amu=None,
    cv_interval_steps=None,
    platform=None,
    precision=None,
)
