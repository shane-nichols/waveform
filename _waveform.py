#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 09:28:37 2018

@author: Shane Nichols

This test module test that that all primitive subclasses of WaveformABC can be
initialized with their default arguments, the waveform shape is correct, and
that the duration agrees with the length of the waveform and the sample rate.
Basic use of unary, binary, and array composites are also tested.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from abc import ABC, abstractmethod


class WaveformCache:
    '''
    This class defines a singleton dictionary used to cache waveforms. 
    Only the abstract Waveform classes need to interface with the cache 
    (the classes in this file) by modifying of the 'waveform' property 
    method. A cache could be implemented as an instance property owned 
    by each object, but this approach is somewhat more flexible. For 
    example, the entire cache can be emptied by calling 'clear'.
    In the current implementation, all classes except Array are cached.
    '''  
    __instance = None
    @staticmethod
    def getInstance():
        if WaveformCache.__instance == None:
            WaveformCache()
        return WaveformCache.__instance

    def __init__(self):
        if WaveformCache.__instance != None:
            raise Exception("This class is a singleton!")
        else:
            WaveformCache.__instance = self

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __delitem__(self, key):
        if key in self:
            del self.__dict__[key]

    def __contains__(self, key):
        return key in self.__dict__

    def __getitem__(self, key):
            return self.__dict__[key]

    def clear(self):
        self.__dict__.clear()


class _WaveformABC(ABC):
    # get instance of singleton cache object
    _cache = WaveformCache.getInstance()

    def __del__(self):
        # remove waveform from cache upon destruction
        del self._cache[id(self)]

    def __setattr__(self, key, value):
        # if an object's data changes, then it's waveform is 
        # removed from the cache
        object.__setattr__(self, key, value)
        del self._cache[id(self)]

    def setattrs(self, **kwargs):
        # allows setting multiple attributes at once
        for k, v in kwargs.items():
            object.__setattr__(self, k, v)
        del self._cache[id(self)] 

    # primitive subclasses that do not directly use a time property
    # will override this getter and will not implement a setter
    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, t):
        self._t = t

    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, dt):
        self._dt = dt

    @property
    def waveform(self):
        if id(self) not in self._cache:
            self._cache[id(self)] = self.get_waveform()
        return self._cache[id(self)]

    @property
    def times(self):
        return self.get_times()

    @property
    def duration(self):
        return self.get_duration()       

    # default implementations of property methods
    @abstractmethod
    def get_waveform(self):
        pass

    def get_times(self):
        return np.arange(0, self.get_duration(), self.dt)

    def get_duration(self):
        return self.t

    # binary compositions
    def __add__(self, other):
        return _Add(self, other)

    def __radd__(self, other):
        return _Add(self, other)

    def __sub__(self, other):
        return _Sub(self, other)

    def __rsub__(self, other):
        return _rSub(self, other)

    def __mul__(self, other):
        return _Mul(self, other)

    def __rmul__(self, other):
        return _Mul(self, other)

    def __truediv__(self, other):
        return _Div(self, other)

    def __rtruediv__(self, other):
        return _rDiv(self, other)

    # utility methods
    def get_homogeneous_dt(self):
        dts = self.dt
        if isinstance(dts, list):
            if all(x == dts[0] for x in dts):
                return dts[0]
            else:
                raise ValueError('Operation not implemented for Waveform Array with  \
                                  heterogeneous sample rates')
        else:
            return dts

    def clear_cache(self):
        self._cache.clear()

    def plot(self, *args, **kwargs):
        if args and isinstance(args[0], Axes):
            ax = args[0]
            p = 1
        else:
            plt.figure()
            ax = plt.axes()
            p = 0
        ax.plot(self.times, self.waveform, *args[p:], **kwargs)
        ax.set_xlabel('Time (sec)')
        return ax


# DERIVED ABSTRACT CLASSES
class _UnaryComposite(_WaveformABC):
    def __init__(self, waveform_object):
        self.waveform_object = waveform_object

    @property
    def t(self):
        return self.waveform_object.t

    @property
    def dt(self):
        return self.waveform_object.dt

    @property
    def waveform(self):
        if id(self.waveform_object) not in self._cache or id(self) not in self._cache:
            self._cache[id(self)] = self.get_waveform()
        return self._cache[id(self)]

    @abstractmethod
    def get_waveform(self):
        pass

    def get_duration(self):
        return self.waveform_object.get_duration()

    def get_times(self):
        return self.waveform_object.get_times()


class _BinaryComposite(_WaveformABC):
    def __init__(self, obj1, obj2):
        if isinstance(obj2, _WaveformABC):
            if np.array_equal(obj1.times, obj2.times):
                self.isscalar = False
        elif np.isscalar(obj2):
            self.isscalar = True
        else:
            raise ValueError('Composed waveforms must have identical time \
                            arrays or one operand must be scalar numeric')
        self.obj1 = obj1
        self.obj2 = obj2

    @property
    def t(self):
        return self.obj1.t

    @property
    def dt(self):
        return self.obj1.dt

    @property
    def waveform(self):
        if id(self.obj1) not in self._cache or \
                id(self.obj2) not in self._cache or \
                id(self) not in self._cache:
            self._cache[id(self)] = self.get_waveform()
        return self._cache[id(self)]

    @abstractmethod
    def get_waveform(self):
        pass

    def get_duration(self):
        return self.obj1.get_duration()

    def get_times(self):
        return self.obj1.get_times()


# OVERLOADED BUILTINS. Defined here instead of in _composites.py 
# to obviate dealing with circular import statements
class _Add(_BinaryComposite):

    def __init__(self, obj1, obj2):
        super().__init__(obj1, obj2)

    def get_waveform(self):
        if self.isscalar:
            return self.obj1.waveform + self.obj2
        else:
            return self.obj1.waveform + self.obj2.waveform

class _Sub(_BinaryComposite):

    def __init__(self, obj1, obj2):
        super().__init__(obj1, obj2)

    def get_waveform(self):
        if self.isscalar:
            return self.obj1.waveform - self.obj2
        else:
            return self.obj1.waveform - self.obj2.waveform

class _rSub(_BinaryComposite):

    def __init__(self, obj1, obj2):
        super().__init__(obj1, obj2)

    def get_waveform(self):
        if self.isscalar:
            return self.obj2 - self.obj1.waveform
        else:
            return self.obj2.waveform - self.obj1.waveform

class _Div(_BinaryComposite):

    def __init__(self, obj1, obj2):
        super().__init__(obj1, obj2)

    def get_waveform(self):
        if self.isscalar:
            return self.obj1.waveform / self.obj2
        else:
            return self.obj1.waveform / self.obj2.waveform

class _rDiv(_BinaryComposite):

    def __init__(self, obj1, obj2):
        super().__init__(obj1, obj2)

    def get_waveform(self):
        if self.isscalar:
            return self.obj2 / self.obj1.waveform
        else:
            return self.obj2.waveform / self.obj1.waveform 

class _Mul(_BinaryComposite):

    def __init__(self, obj1, obj2):
        super().__init__(obj1, obj2)

    def get_waveform(self):
        if self.isscalar:
            return self.obj1.waveform * self.obj2
        else:
            return self.obj1.waveform * self.obj2.waveform
