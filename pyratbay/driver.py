# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'run',
]

from . import opacity as op
from . import tools as pt
from . import Pyrat, Atmosphere


def run(cfile, with_log=True):
    """
    Pyrat Bay command-line-interface run

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    with_log: Bool
        Flag to save screen outputs to file (True) or not (False)
        (e.g., to prevent overwritting log of a previous run).
    """
    inputs, log = pt.parse(cfile, with_log)
    runmode = inputs.runmode

    # Call lineread
    if runmode == 'tli':
        if inputs.tlifile is None:
            log.error('Undefined TLI file (tlifile)')
        # return wl variabels to their original units
        inputs.wl_low /= pt.u(inputs.wlunits)
        inputs.wl_high /= pt.u(inputs.wlunits)
        op.make_tli(
            inputs.dblist, inputs.pflist, inputs.dbtype,
            inputs.tlifile[0], inputs.wl_low, inputs.wl_high,
            inputs.wlunits, log,
        )
        return

    # Initialize and run atmosphere
    if runmode == 'atmosphere':
        return Atmosphere(inputs, log=log)

    # Initialize pyrat and execute calculations
    pyrat = Pyrat(inputs, log)
    if runmode == 'opacity':
        pyrat.compute_opacity()

    if runmode == "spectrum":
        pyrat.run()

    if runmode == 'radeq':
        pyrat.radiative_equilibrium()

    if runmode == 'retrieval':
        pyrat.retrieval()

    return pyrat
