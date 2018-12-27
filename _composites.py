#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 09:28:37 2018

@author: Shane Nichols
"""

import numpy as np
from waveform._waveform import _WaveformABC, _UnaryComposite, _BinaryComposite, \
    _Add, _Mul


class Array(list, _WaveformABC):
    # There is probably a better base class for this. 'list' works...
    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]  # to allow passing of interables
        super().__init__(args)

    def __getattr__(self, attr):
        return list(map(lambda x: getattr(x, attr), self))

    def __setattr__(self, name, value):
        for obj in self:
            setattr(obj, name, value)

    # override binary operations defined in 'list'
    def __radd__(self, other):
        return _Add(self, other)

    def __add__(self, other):
        return _Add(self, other)
        
    def __mul__(self, other):
        return _Mul(self, other)

    def __rmul__(self, other):
        return _Mul(self, other)

    def setattrs(self, **kwargs):
        # allows setting multiple attributes at once
        for obj in self:
            obj.setattrs(**kwargs)

    @property
    def waveform(self):
        # putting this method prevents caching Arrays
        return self.get_waveform()

    def get_waveform(self):
        vals = map(lambda x: getattr(x, 'waveform'), self)
        return np.hstack(list(vals))

    def get_duration(self):
        return sum(map(lambda x: getattr(x, 't'), self))

    def get_times(self):
        # times property is overloaded here to return concatenated array
        times = [obj.get_times() for obj in self]
        for i in range(len(times) - 1):
            times[i+1] += (times[i][-1] + times[i+1][1])
        return np.concatenate(times).ravel()


###############################################################################
class Shutter(_UnaryComposite):

    def __init__(self, waveform_object, level=0, tbefore=5, tafter=5):
        super().__init__(waveform_object)
        self.level = level
        self.tbefore = tbefore
        self.tafter = tafter

    def get_waveform(self):
        dt = self.waveform_object.get_homogeneous_dt()
        pts_before = self.tbefore / 1000 / dt
        pts_after = self.tafter / 1000 / dt
        return self.expandlogical(self.waveform_object.waveform > self.level,
                pts_before, pts_after)

    @staticmethod
    def expandlogical(array, ptsBefore, ptsAfter):
        array = array.astype(bool)
        ptsBefore = int(ptsBefore)
        ptsAfter = int(ptsAfter)
        sz = len(array)
        array = np.concatenate(
                (np.array([array[0]]), np.diff(array.astype(np.int8))))
        ris = array == 1
        fal = array == -1
        b4 = sum(ris[:ptsBefore])
        ris = np.concatenate(
                (ris[ptsBefore:], np.zeros(ptsBefore, dtype=bool)))
        fal = np.concatenate(
                (np.zeros(ptsAfter, dtype=bool), fal[0:(sz-ptsAfter)]))
        ris = ris.astype(int)
        ris[0] += b4
        return np.cumsum(ris - fal.astype(int)).astype(bool)


###############################################################################
class Map(_UnaryComposite):

    def __init__(self, waveform_object, function):
        super().__init__(waveform_object)
        self.function = function

    def get_waveform(self):
        return self.function(self.waveform_object.waveform)


###############################################################################
class PulseRate(_UnaryComposite):
    '''
    PulseRate Waveform object. Creates pulse train where the instantanteous 
    pulse rate is given by an arbitrary rate function. An absolute value will
    be applied to the rate function as 'negative rates' are not meaningful.
    The waveform attribute of the waveform_object defines the rate function.
    Properties and defaults:

    waveform_object   Waveform object defining the rate function
    initial_val = 0.5 Number between 0 and 1; With 1 first pulse is at t = 0
    ton = 0.001       Duration of each pulse
    height = 1        Height of pulses above offset
    offset = 0        Offset value from 0
    '''
    def __init__(self, waveform_object, initial_val=0.5, ton=0.001, height=1,
                 offset=0):
        super().__init__(waveform_object)
        self.initial_val = initial_val
        self.ton = ton
        self.height = height
        self.offset = offset

    def get_waveform(self):
        fx = self.waveform_object.waveform
        # integral of fx
        dt = self.waveform_object.get_homogeneous_dt()
        s = self.initial_val + np.cumsum(np.abs(fx)) * dt
        pulse_inds = np.nonzero(np.diff(np.insert(s - np.mod(s, 1), 0, 0)))[0]
        on_samples = int(self.ton / dt) + 1
        pulse_inds = pulse_inds[pulse_inds < (len(fx) - on_samples)]
        wf = np.zeros(int(self.get_duration() / dt))
        for ind in pulse_inds:
            wf[ind:(ind+on_samples)] = self.height
        if self.offset:
            wf = wf + self.offset
        return wf
    