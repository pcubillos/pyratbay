import os
import sys
#import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay.tools as pt


def test_path():
    assert pt.path('file.txt')   == "./file.txt"
    assert pt.path('./file.txt') == "./file.txt"
    assert pt.path('/home/user/file.txt') == "/home/user/file.txt"
