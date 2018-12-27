#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 09:28:37 2018
@author: Shane Nichols

Waveform primitive objects. These objects must implement three attributes:
get_waveform(), A method that returns a numpy 1D array
dt,             The sample period of the waveform
t,              The duration of the waveform. If the duration is needed to
                    compute the waveform, then 't' is just a data attribute.
                    If the value of 't' is composed from other attributes, 
                    then the class should implement a @property getter for 
                    't' that returns the duration (e.g., see Step and Ramp).
                    Be careful in how this getter is implemented!
"""
import numpy as np
from waveform._waveform import _WaveformABC


class Step(_WaveformABC):
    '''Step Waveform object. Creates a single step.
    Properties and defaults:
        dt = 0.0001,     Sampling period, in sec
        tBefore = 0.1    Duration at offset before ramp, in seconds
        tStep = 0.8      Duration of ramp, in seconds
        tAfter = 0.1     Duration at offset after ramp, in seconds
        offset = 0       Offset from zero
        height = 1       End value of the ramp, relative to offset '''

    def __init__(self, tBefore=0.1, tStep=0.8, tAfter=0.1, offset=0, height=1, 
                 dt=0.0001):
        self.tBefore = tBefore
        self.tStep = tStep
        self.tAfter = tAfter
        self.offset = offset
        self.height = height
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        iBefore = self.offset + np.zeros(np.int(self.tBefore / self.dt))
        iAfter = self.offset + np.zeros(np.int(self.tAfter / self.dt))
        iStep = (self.offset + self.height) + \
            np.zeros(np.int(self.tStep / self.dt))
        wf = np.concatenate((iBefore, iStep, iAfter))
        return wf

    @property
    def t(self):
        before = int(self.tBefore / self.dt)
        after = int(self.tAfter / self.dt)
        step = int(self.tStep / self.dt)
        return (before + after + step) * self.dt


###############################################################################
class Ramp(_WaveformABC):
    '''
    Ramp Waveform object. Creates a single ramp.
    Properties and defaults:
    dt = 0.0001,     Sampling period, in sec
    tBefore = 0.1    Duration at offset before ramp, in seconds
    tRamp = 0.8      Duration of ramp, in seconds
    tAfter = 0.1     Duration at offset after ramp, in seconds
    offset = 0       offset from zero
    height = 1       End value of the ramp, relative to offset

    '''
    def __init__(self, tBefore=0.1, tRamp=0.8, tAfter=0.1, offset=0, height=1,
                 dt=0.0001):
        self.tBefore = tBefore
        self.tRamp = tRamp
        self.tAfter = tAfter
        self.offset = offset
        self.height = height
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        iBefore = self.offset + np.zeros(int(self.tBefore / self.dt))
        iAfter = self.offset + np.zeros(int(self.tAfter / self.dt))
        iRamp = np.linspace(
                self.offset, self.offset + self.height, int(self.tRamp / self.dt))
        wf = np.concatenate((iBefore, iRamp, iAfter))
        return wf

    @property
    def t(self):
        before = int(self.tBefore / self.dt)
        after = int(self.tAfter / self.dt)
        ramp = int(self.tRamp / self.dt)
        return (before + after + ramp) * self.dt


###############################################################################
class SinModulated(_WaveformABC):
    '''
    SinModulated Waveform object. Creates a sine wave in which the frequency
    is modulated by another sine wave.
    Properties and defaults:
    dt = 0.0001     Sampling period, in sec
    t = 1           Total time, in seconds
    fcarrier = 50   Center pulse frequency, in Hz
    fmod = 1        Modulation frequency, in Hz
    amod = 20       Modulation amplitude, in Hz
    height = 1      Height of pulses above offset
    offset = 0      Offset value from 0
    '''

    def __init__(self, fcarrier=50, fmod=1, amod=20, height=1, offset=0,
                 t=1, dt=0.0001):
        self.fcarrier = fcarrier
        self.fmod = fmod
        self.amod = amod
        self.height = height
        self.offset = offset
        self.t = t
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        x = np.arange(0, self.t, self.dt)
        wf = self.offset + self.height * (
            np.sin(2 * np.pi * x * self.fcarrier -
                   self.amod / self.fmod * np.sin(2 * np.pi * x * self.fmod)))
        return wf


###############################################################################
class SquareWave(_WaveformABC):
    '''SquareWave Waveform object. When the phase is zero, the signal starts at
    a rising edge for any duty cycle.
    Properties and defaults:
    dt = 0.0001     Sampling period, in sec
    t = 1           Total time, in seconds
    f = 50          Frequency, in Hz
    phase = 0       Phase offset, in radians
    duty = 0.5      Duty cycle (fraction of on vs off)
    height = 1      Height of pulses above offset
    offset = 0      Offset value from 0 '''

    def __init__(self, f=50, phase=0, height=1, duty=0.5, offset=0, t=1,
                 dt=0.0001):
        self.f = f
        self.phase = phase
        self.height = height
        self.duty = duty
        self.offset = offset
        self.t = t
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        x = np.arange(0, self.t, self.dt)
        phasei = (1/2 - self.duty) * np.pi
        threshold = np.sin(phasei)
        ''' a tiny number is added to the phase so transitions are not
        exactly zero, otherwise N rising and N falling will not be
        equal for duty = 0.5 '''
        wf = self.offset + self.height * (
            np.sin(2*np.pi*x*self.f + phasei + self.phase + 0.000000001)
            > threshold)
        return wf


###############################################################################
class Constant(_WaveformABC):
    '''Constant Waveform object. Creates a constant waveform.
    Properties and defaults:
        dt = 0.0001,     Sampling period, in sec
        t = 1            Duration, in sec
        value = 1        Value of the waveform '''

    def __init__(self, value=1, t=1, dt=0.0001):
        self.value = value
        self.t = t
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        wf = self.value + np.zeros(int(self.t / self.dt))
        self.t = len(wf) * self.dt
        return wf


###############################################################################
class SinWave(_WaveformABC):
    '''SinWave Waveform object. Simple Sine wave.
    For a cosine, set 'phase = pi/2'.
    This waveform can also be created via Function().
    Properties and defaults:
        dt = 0.0001     Sampling period, in sec
        T = 1           Total time, in seconds
        f = 50          Frequency, in Hz
        phase = 0       Phase offset, in radians
        height = 1      Height of pulses above offset
        offset = 0      Offset value from 0'''

    def __init__(self, f=50, phase=0, height=1, offset=0, t=1, dt=0.0001):
        self.f = f
        self.phase = phase
        self.height = height
        self.offset = offset
        self.t = t
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        x = np.arange(0, self.t, self.dt)
        wf = self.offset + self.height * \
            np.sin(2 * np.pi * x * self.f + self.phase)
        return wf


###############################################################################
class TriangleWave(_WaveformABC):
    '''TriangleWave Waveform object.
    For a cosine, set 'phase = pi/2'.
    Properties and defaults:
        dt = 0.0001     Sampling period, in sec
        t = 1           Total time, in seconds
        f = 50          Frequency, in Hz
        phase = 0       Phase offset, in radians
        height = 1      Height of pulses above offset
        offset = 0      Offset value from 0'''

    def __init__(self, f=50, phase=0, height=1, offset=0, t=1, dt=0.0001):
        self.f = f
        self.phase = phase
        self.height = height
        self.offset = offset
        self.t = t
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        x = np.arange(0, self.t, self.dt)
        wf = x * self.f + self.phase / (2 * np.pi)
        wf = self.offset + (2 * self.height) \
            * (0.5 - np.abs((np.mod(wf, 1) - 0.5)))
        return wf


###############################################################################
class PulseModulated(_WaveformABC):
    '''PulseModulated Waveform object. Creates pulse train of constant duty
    cycle (mean stimulus) and sin-modulated frequency.
    Properties and defaults:
        dt = 0.0001     Sampling period, in sec
        t = 1           Total time, in seconds
        fcarrier = 50   Center pulse frequency, in Hz
        fmod = 1        Modulation frequency, in Hz
        amod = 20       Modulation amplitude, in Hz
        duty = 0.5      Duty cycle (fraction of on vs off)
        height = 1      Height of pulses above offset
        offset = 0      Offset value from 0'''

    def __init__(self, fcarrier=50, fmod=1, amod=20, duty=0.5, height=1,
                 offset=0, dt=0.0001, t=1):
        self.fcarrier = fcarrier
        self.fmod = fmod
        self.amod = amod
        self.duty = duty
        self.height = height
        self.offset = offset
        self.dt = dt
        self.t = t
        # super().__init__()

    def get_waveform(self):
        x = np.arange(0, self.t * 2 * np.pi, self.dt * 2 * np.pi)
        phase = (1/2 - self.duty) * np.pi
        threshold = np.sin(phase)
        wf = self.offset + self.height * (
            np.sin(x * self.fcarrier + phase -
                   self.amod / self.fmod * np.sin(x * self.fmod)) > threshold)
        return wf


###############################################################################
class SinVariance(_WaveformABC):
    '''SinVariance Waveform object. Creates a waveform with constant mean and
    sin modulated variance. The variance is uniformly distributed and
    fluctuates over a characteristic timescale
    Properties and defaults:
        dt = 0.0001     Sampling period, sec
        t = 1           Total time of the stimulus, sec
        tau = 0.003     Characteristic timescale of fluctuations, sec
        I0 = 1          Mean value
        sigma0 = 1      Central variance
        dSigma = 0.5    Variance modulation amplitude
        sigmaF = 0.2    Variance modulation frequency, Hz'''

    def __init__(self, tau=0.003, I0=1, sigma0=1, dSigma=0.5, sigmaF=0.2,
                 t=1, dt=0.0001):
        self.tau = tau
        self.I0 = I0
        self.sigma0 = sigma0
        self.dSigma = dSigma
        self.sigmaF = sigmaF
        self.t = t
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        N = int(self.t / self.dt)
        wf = self.sigma0 * (1 + self.dSigma * np.sin(
                2 * np.pi * self.sigmaF * np.arange(0, N) * self.dt))
        wf = np.sqrt(2 * wf ** 2 * self.dt / self.tau) * np.random.randn(N)
        wf[0] = self.I0  # initialize at mean
        for n in range(1, N):
            wf[n] = wf[n-1] + (self.I0 - wf[n-1]) / self.tau * self.dt + wf[n]
        return wf


###############################################################################
class StepRamp(_WaveformABC):
    '''StepRamp Waveform object. Creates a set of evenly spaced steps
    of linearlly increasing height. An off-time can be placed between steps.
    Properties and defaults:
        dt = 0.0001,            Sampling period, in sec
        tOn = 0.1,              Duration of pulse, in seconds
        tOff = 0.1,             Duration of time between pulses, in seconds
        heightFirst = 0.2,      Amplitude of the first pulse relative to offset
        heightLast = 1,         Amplitude of the last pulse relative to offset
        offset = 0,             Offset value from 0
        Nsteps = 5,             Number of steps, postive integer'''

    def __init__(self, tOn=0.1, tOff=0.1, heightFirst=0.2, heightLast=1,
                 offset=0, Nsteps=5, dt=0.0001):
        self.tOn = tOn
        self.tOff = tOff
        self.heightFirst = heightFirst
        self.heightLast = heightLast
        self.offset = offset
        self.Nsteps = Nsteps
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        N_on = int(self.tOn / self.dt)
        N_off = int(self.tOff / self.dt)
        N = N_on + N_off
        wf = np.empty(self.Nsteps * N + N_off)
        amps = np.linspace(self.heightFirst + self.offset,
                           self.heightLast + self.offset,
                           self.Nsteps)
        for i in range(self.Nsteps):
            wf[i*N:i*N+N_off] = self.offset
            wf[i*N+N_off:(i+1)*N] = amps[i]
        wf[N*self.Nsteps::] = self.offset
        return wf

    @property
    def t(self):
        N_on = int(self.tOn / self.dt)
        N_off = int(self.tOff / self.dt)
        return (self.Nsteps * (N_on + N_off) + N_off) * self.dt


###############################################################################
class General(_WaveformABC):
    '''General Waveform object. Allows one to load any waveform into a Waveform
    object. Properties:
        dt,    Sampling period, in sec
        wf,    A 1D numpy array of a numeric class'''

    def __init__(self, wf=np.array([]), dt=0):
        self._waveform = wf
        self.dt = dt
        # super().__init__()

    @property
    def waveform(self):
        # prevent caching (no point here)
        return self._waveform

    @property
    def t(self):
        return len(self._waveform) * self.dt

    def get_waveform(self):
        return self._waveform


###############################################################################
class SinAmpMod(_WaveformABC):
    '''SinAmpMod Waveform object. Creates a amplitude modulated sin wave.
    Properties and defaults:
        dt = 0.0001     Sampling period, in sec
        t = 1           Total time, in seconds
        fcarrier = 50   Carrier frequency, in Hz
        fmod = 1        Amplitude modulation frequency, in Hz
        amod = 1       Modulation amplitude
        height = 1      Height of waveform above offset
        offset = 0      Offset value from 0'''

    def __init__(self, fcarrier=50, fmod=1, amod=1, height=1, offset=0,
                 t=1, dt=0.0001):
        self.fcarrier = fcarrier
        self.fmod = fmod
        self.amod = amod
        self.height = height
        self.offset = offset
        self.t = t
        self.dt = dt
        # super().__init__()

    def get_waveform(self):
        x = np.arange(0, self.t * 2 * np.pi, self.dt * 2 * np.pi)
        wf = self.offset + self.height * \
            (np.sin(x * self.fcarrier) + 1) / 2 * \
            (1 - self.amod * (np.sin(x * self.fmod) + 1) / 2)
        return wf


###############################################################################
class Function(_WaveformABC):
    '''Function Waveform object. Creates a waveform from a lambda function
    Properties and defaults:
        dt = 0.0001     Sampling period, in sec
        t = 1           Total time, in seconds
        function        function with signature lambda t: ...'''

    def __init__(self, function=lambda t: t, t=1, dt=0.0001):
        self.function = function
        self.t = t
        self.dt = dt

    def get_waveform(self):
        return self.function(self.get_times())


__all__ = [cls.__name__ for cls in _WaveformABC.__subclasses__() 
           if not cls.__name__.startswith("_")] 
""