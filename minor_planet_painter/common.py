#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
from datetime import datetime, timedelta
from astropy.time import Time
from astroquery.jplhorizons import Horizons

## Path to common.py
BASE = os.path.dirname(os.path.abspath(__file__))
## Path to DATA
DATA = os.path.normpath(os.path.join(BASE, "../data"))
## Path to MPCORB
MPCORB = os.path.normpath(os.path.join(DATA, "MPCORB.DAT"))


def jd2utc(jd):
    """Convert jd to utc in '%Y-%m-%dT%H:%M:%S.%f'.

    Parameter
    ---------
    jd : float
        julian day

    Return
    ------
    t_utc : str
        time in utc
    """
    t = Time(str(jd), format='jd', scale='utc')
    t_utc = t.datetime
    return t_utc


def utc2jd(utc):
    """Convert utc to jd in '%Y-%m-%dT%H:%M:%S.%f'.

    Parameter
    ---------
    utc : str
        time in utc

    Return
    ------
    t_jd : float
        julian day
    """
    t = Time(str(utc), format='isot', scale='utc')
    t_jd = t.jd
    return t_jd


def get_planet_positions(t_utc):
    """Get positions of planets.

    Parameter
    ---------
    t_utc : str
        time in utc

    Return
    ------
    positions : dict
        positions of planets
    """
    planets = {
        "Mercury": 1,
        "Venus": 2,
        "Earth": 3,
        "Mars": 4,
        "Jupiter": 5,
        "Saturn": 6,
        "Uranus": 7,
        "Neptune": 8,
        "Pluto": 9
    }
    positions = {}
    jd = utc2jd(t_utc)
    for name, pid in planets.items():
        obj = Horizons(id=pid, location="@0", epochs=jd)
        vec = obj.vectors()
        x = float(vec['x'][0])
        y = float(vec['y'][0])
        z = float(vec['z'][0])
        positions[name] = (x, y, z)
    return positions


def get_planet_orbits(t_utc):
    """Get orbits of planets.

    Parameter
    ---------
    t_utc : str
        time in utc

    Return
    ------
    orbits : dict
        orbits of planets
    """

    # To datetime
    t_utc = datetime.fromisoformat(t_utc)

    planet_params = {
        "Mercury":  {"id": 1, "days": 120,     "step": "1d"},
        "Venus":    {"id": 2, "days": 250,     "step": "1d"},
        "Earth":    {"id": 3, "days": 370,     "step": "1d"},
        "Mars":     {"id": 4, "days": 750,     "step": "2d"},
        "Jupiter":  {"id": 5, "days": 4400,    "step": "10d"},
        "Saturn":   {"id": 6, "days": 11000,   "step": "20d"},
        "Uranus":   {"id": 7, "days": 31000,   "step": "50d"},
        "Neptune":  {"id": 8, "days": 60000,   "step": "100d"},
        "Pluto":    {"id": 9, "days": 90000,   "step": "150d"},
    }
    orbits = {}
    for name, params in planet_params.items():
        start_date = (t_utc - timedelta(days=params["days"]/2)).strftime("%Y-%m-%d")
        stop_date = (t_utc + timedelta(days=params["days"]/2)).strftime("%Y-%m-%d")

        obj = Horizons(id=params["id"], location="@0", epochs={
            'start': start_date,
            'stop': stop_date,
            'step': params["step"]
        })
        vec = obj.vectors()
        x = np.array(vec['x'], dtype=float)
        y = np.array(vec['y'], dtype=float)
        orbits[name] = (x, y)

    return orbits


def solve_kepler_eq(M, e, max_iter=10):
    """Solve Kepler equation.

    See
    https://en.wikipedia.org/wiki/Kepler%27s_equation

    Parameters
    ----------
    M : float
        mean anomaly
    e : float
        eccentricity

    Return
    ------
    E : float
        eccentric anomaly
    """
    E = M.copy()
    for _ in range(max_iter):
        E = E - (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
    return E


def base36_to_int(s):
    """
    Convert a base-36 string to an integer.
    Base-36 uses digits ``0-9`` and letters ``A-Z`` (case-insensitive) to represent
    values from 0 to 35. For example, ``"10"`` in base-36 equals 36 in decimal.

    Parameters
    ----------
    s : str
        The base-36 encoded string. Case-insensitive.

    Return
    ------
    value : int
        The integer value corresponding to the base-36 input.
    """
    s = s.upper()
    value = 0
    for c in s:
        if '0' <= c <= '9':
            digit = ord(c) - ord('0')
        elif 'A' <= c <= 'Z':
            digit = ord(c) - ord('A') + 10
        else:
            raise ValueError(f"Invalid base36 character: {c}")
        value = value*36 + digit
    return value


def mpcepoch2jd(epoch_code):
    """Convert MPC epoch to jd.

    See https://www.minorplanetcenter.net/iau/info/PackedDates.html

    Parameter
    ---------
    epoch_code : str
        MPC epoch in packed format

    Return
    ------
    date_jd : float
        epoch in jd
    """
    epoch_code = epoch_code.strip()
    try:
        return float(epoch_code)
    except ValueError:
        pass

    if '.' in epoch_code:
        integer_part, fractional_part = epoch_code.split('.')
    else:
        integer_part, fractional_part = epoch_code, '0'

    c = integer_part[0].upper()
    if not ('I' <= c <= 'Z'):
        raise ValueError(f"Unknown century code: {c}")
    century = 1800 + 100*(ord(c) - ord('I'))

    yy = int(integer_part[1:3])
    year = century + yy

    day_of_year_int = base36_to_int(integer_part[3:5])
    day_of_year = day_of_year_int + float('0.' + fractional_part)

    date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
    date = date.isoformat()
    date_jd = utc2jd(date)
    return date_jd


# Not used.
def update_mean_anomaly_longterm(M, n, dt_days, step_days=365.25):
    """Alternative to M_t = M + n*dt_days.
    """
    M_t = np.zeros_like(M)
    for idx in range(len(M)):
        dt = dt_days[idx]
        M0 = M[idx]
        n0 = n[idx]

        Nsteps = int(abs(dt)//step_days)
        delta_sign = 1 if dt > 0 else -1
        M_acc = M0
        for _ in range(Nsteps):
            M_acc += n0 * step_days * delta_sign
        remainder = dt - Nsteps*step_days*delta_sign
        M_acc += n0 * remainder
        M_t[idx] = M_acc
    return np.mod(M_t, 2*np.pi)
