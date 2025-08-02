# tests/test_data_analysis_utils.py

import pytest
import pandas as pd
import numpy as np
from src import data_analysis_utils as dau


def test_save_resistance_temperature():
    data = pd.DataFrame({
        'T': [10, 20, 30],
        'R': [100, 200, 300]
    })
    result = dau.save_resistance_temperature(data, 'T', 'R')

    assert np.allclose(result['T'], [10, 20, 30])
    assert np.allclose(result['R'], [100, 200, 300])
    expected_rho = (0.05 / 1.221) * data['R']
    expected_sigma = 1 / expected_rho
    assert np.allclose(result['rho'], expected_rho)
    assert np.allclose(result['sigma'], expected_sigma)


def test_save_magnetoresistance():
    data = pd.DataFrame({
        'T': [10, 20],
        'R': [100, 200],
        'F': [1, 2]
    })
    result = dau.save_magnetoresistance(data, 'T', 'R', 'F', rate=10)

    assert 'Field' in result.columns
    assert np.allclose(result['Field'], 10 * data['F'])


def test_save_magnetoresistance_edc():
    data = pd.DataFrame({
        'T': [4.2],
        'V': [0.5],
        'F': [20000],
        'I': [2.0]  # not used from here
    })
    I_test = 2.0
    result = dau.save_magnetoresistance_edc(data, 'T', 'V', 'F', I=I_test)

    expected_r = 0.5 / (I_test / 1000)
    expected_field = 20000 / 10000

    assert np.isclose(result['R (Ohns)'][0], expected_r)
    assert np.isclose(result['Field (T)'][0], expected_field)
