#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:47:16 2018

@author: Shane Nichols

This test module test that that all primitive subclasses of WaveformABC can be
initialized with their default arguments, the waveform shape is correct, and
that the duration agrees with the length of the waveform and the sample rate.
Basic use of unary, binary, and array composites are also tested.
"""

import waveform as wf
import numpy as np
from waveform._waveform import _WaveformABC, _UnaryComposite
import pytest


@pytest.fixture
def primitives():
    return [cls() for cls in _WaveformABC.__subclasses__() 
                if repr(cls).find("_primitives") != -1]

@pytest.fixture
def unary():
    waveform_object = wf.Ramp(height=10)
    return [
        wf.Map(waveform_object, lambda x: np.exp(x)),
        wf.Shutter(waveform_object),
        wf.PulseRate(waveform_object)
    ]

def test_wf_size(primitives):
    for wfObj in primitives:
        assert wfObj.waveform.ndim == 1

def test_wf_size_unary(unary):
    for wfObj in unary:
        assert wfObj.waveform.ndim == 1

def test_wf_duration(primitives):
    for wfObj in primitives:
        assert wfObj.duration == len(wfObj.waveform) * wfObj.dt

def test_wf_duration_unary(unary):
    for wfObj in unary:
        assert wfObj.duration == len(wfObj.waveform) * wfObj.dt

def test_wf_array_duration(primitives, unary):
    wfArray = wf.Array(primitives + unary)
    assert wfArray.duration == sum(obj.duration for obj in wfArray)

def test_binary_composition():
    a = wf.Step(offset=1)
    b = wf.SinWave(offset=2)
    ab = wf.Array(a, b)
    aa = wf.Array(a, a)
    assert np.array_equal(ab.waveform + aa.waveform, (ab + aa).waveform)
    assert np.array_equal(ab.waveform - aa.waveform, (ab - aa).waveform)
    assert np.array_equal(ab.waveform * aa.waveform, (ab * aa).waveform)
    assert np.array_equal(ab.waveform / aa.waveform, (ab / aa).waveform)
    assert np.array_equal(5 + aa.waveform, (5 + aa).waveform)
    assert np.array_equal(5 - aa.waveform, (5 - aa).waveform)
    assert np.array_equal(5 * aa.waveform, (5 * aa).waveform)
    assert np.array_equal(5 / aa.waveform, (5 / aa).waveform)
    assert np.array_equal(ab.waveform + 5, (ab + 5).waveform)
    assert np.array_equal(ab.waveform - 5, (ab - 5).waveform)
    assert np.array_equal(ab.waveform * 5, (ab * 5).waveform)
    assert np.array_equal(ab.waveform / 5, (ab / 5).waveform)

