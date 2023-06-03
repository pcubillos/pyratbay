# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'run',
]

from . import opacity as op
from . import tools as pt
from . import Pyrat, Atmosphere


def run(cfile, run_step=None, with_log=True):
    """
    Pyrat Bay initialization driver.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    run_step: String
        DEPRECATED
    with_log: Bool
        Flag to save screen outputs to file (True) or not (False)
        (e.g., to prevent overwritting log of a previous run).
    """
    inputs, log = pt.parse(cfile, with_log)
    runmode = inputs.runmode

    # TBD: deprecate run_step
    if run_step == 'dry':
        return inputs

    # Call lineread:
    if runmode == 'tli':
        if inputs.tlifile is None:
            log.error('Undefined TLI file (tlifile)')
        op.make_tli(
            inputs.dblist, inputs.pflist, inputs.dbtype,
            inputs.tlifile[0], inputs.wllow, inputs.wlhigh,
            inputs.wlunits, log,
        )
        return

    # Initialize atmosphere:
    if runmode == 'atmosphere':
        return Atmosphere(inputs, log)

    # Initialize pyrat object:
    pyrat = Pyrat(inputs, log)
    # Stop and return if requested:
    if run_step == 'init':
        return pyrat

    # Calculate extinction-coefficient file if requested:
    if runmode == 'opacity':
        pyrat.compute_opacity()
        return pyrat

    # Compute spectrum and return pyrat object if requested:
    if runmode == "spectrum":
        pyrat.run()
        return pyrat

    if runmode == 'radeq':
        pyrat.radiative_equilibrium()
        return pyrat

    if runmode == 'mcmc':
        pyrat.retrieval()
        return pyrat
