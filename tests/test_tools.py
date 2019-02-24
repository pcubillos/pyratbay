import os
import sys
import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay.tools as pt


def test_path():
    assert pt.path('file.txt')   == "./file.txt"
    assert pt.path('./file.txt') == "./file.txt"
    assert pt.path('/home/user/file.txt') == "/home/user/file.txt"


@pytest.mark.parametrize('data',
    [[False, True, True, False],
     [0,1,1,0],
     (False, True, True, False),
     np.array([False, True, True, False])])
def test_ifirst_type(data):
    assert pt.ifirst(data) == 1


@pytest.mark.parametrize('data',
    [[False, True, True, False],
     [0,1,1,0],
     (False, True, True, False),
     np.array([False, True, True, False])])
def test_ilast_type(data):
    assert pt.ilast(data) == 2


def test_wrap():
    info = []
    text = "Pyrat atmospheric model\n"
    pt.wrap(info, text)
    assert info == ['Pyrat atmospheric model']
    text = "Pressure = 1.0 bar\nTemperature = 1000.0 K"
    pt.wrap(info, text, indent=2)
    assert info == ['Pyrat atmospheric model',
                    '  Pressure = 1.0 bar',
                    '  Temperature = 1000.0 K']
