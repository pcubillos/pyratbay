# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import subprocess
import pytest
import re

from conftest import make_config

import pyratbay as pb
import pyratbay.io as io
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')

# Short run will trigger warning when evaluating GR tests, ignore it.
@pytest.mark.filterwarnings("ignore: divide by zero encountered")
def test_mcmc_transmission(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)


@pytest.mark.skip(reason='TBD')
def test_mcmc_emission(tmp_path):
    reset = {
        'rt_path': 'emission',
        'kurucz': f'{ROOT}tests/inputs/mock_fp00k0odfnew.pck',
        'log_gstar': '4.5',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset=reset,
        remove=[param],
    )
    error = re.escape(undefined_mcmc[param])
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)

