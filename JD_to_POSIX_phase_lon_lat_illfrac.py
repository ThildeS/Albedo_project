#!/usr/bin/env python3
"""
Given either:
  (a) a NEF filename  -> read EXIF capture time, compute JD & lunar info
  (b) a UTC Julian Date (number) -> compute lunar info directly

Print one line:
  <JD_f15.7> <YYYY-MM-DDTHH:MM:SS[.ffffff]Z> <phase_angle_deg> <sublunar_lon_deg> <sublunar_lat_deg> <illuminated_fraction>

Examples:
  python JD_to_phase_illfrac_lon_lat.py /path/to/image.NEF
  python JD_to_phase_illfrac_lon_lat.py 2460200.9802545

Deps:
  - exiftool (for NEF)
  - pip install skyfield jplephem
"""

import sys
import json
import re
import subprocess
import math
from datetime import datetime, timedelta, timezone
from typing import Optional

# =========================
# === (Unchanged) EXIF -> JD helpers
# =========================

def _run_exiftool_json(path: str) -> dict:
    """Call exiftool for specific tags and return JSON dict."""
    proc = subprocess.run(
        [
            "exiftool",
            "-j",  # JSON output
            "-n",  # numeric where applicable
            "-DateTimeOriginal",
            "-SubSecTimeOriginal",
            "-OffsetTimeOriginal",
            path,
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    data = json.loads(proc.stdout)
    if not data:
        raise RuntimeError(f"exiftool returned no data for {path}")
    return data[0]

def _parse_exif_datetime(dt_str: str) -> datetime:
    """Parse EXIF 'YYYY:MM:DD HH:MM:SS' to naive datetime."""
    return datetime.strptime(dt_str, "%Y:%m:%d %H:%M:%S")

def _parse_subsec(subsec: Optional[str]) -> float:
    """
    Convert SubSecTimeOriginal (digits) to fractional seconds.
    e.g., '123' -> 0.123 seconds; None -> 0.0
    """
    if not subsec:
        return 0.0
    s = str(subsec).strip()
    if not s.isdigit():
        s = "".join(ch for ch in s if ch.isdigit())
        if not s:
            return 0.0
    return float("0." + s)

def _parse_offset_minutes(offset_str: Optional[str]) -> Optional[int]:
    """
    Parse OffsetTimeOriginal like '+01:00' or '-0530' to minutes.
    Returns None if missing/unparsable.
    """
    if not offset_str:
        return None
    s = offset_str.strip()
    m = re.fullmatch(r"([+-])(\d{2}):?(\d{2})?", s)
    if not m:
        return None
    sign = -1 if m.group(1) == "-" else 1
    hh = int(m.group(2))
    mm = int(m.group(3)) if m.group(3) else 0
    return sign * (hh * 60 + mm)

def _to_julian_date(dt_utc: datetime) -> float:
    """
    Compute Julian Date (UTC) from a datetime.
    Algorithm per Meeus, valid for Gregorian dates.
    """
    if dt_utc.tzinfo is not None:
        dt_utc = dt_utc.astimezone(timezone.utc).replace(tzinfo=None)

    Y = dt_utc.year
    M = dt_utc.month
    D = dt_utc.day
    h = dt_utc.hour
    m = dt_utc.minute
    s = dt_utc.second + dt_utc.microsecond / 1_000_000.0

    frac_day = (h + (m + s / 60.0) / 60.0) / 24.0

    if M <= 2:
        Y -= 1
        M += 12

    A = Y // 100
    B = 2 - A + (A // 4)

    jd = int(365.25 * (Y + 4716)) \
         + int(30.6001 * (M + 1)) \
         + D + B - 1524.5 + frac_day
    return jd

def nef_jd_string(filepath: str) -> str:
    """
    Read EXIF capture time from NEF and return Julian Date as f15.7 string.
    Uses DateTimeOriginal + SubSecTimeOriginal, and applies OffsetTimeOriginal
    if present to convert local capture time to UTC. If DateTimeOriginal is
    missing, raises an error.
    """
    tags = _run_exiftool_json(filepath)

    dto = tags.get("DateTimeOriginal")
    if not dto:
        raise RuntimeError(f"DateTimeOriginal not found in {filepath}")

    dt = _parse_exif_datetime(dto)
    subsec = _parse_subsec(tags.get("SubSecTimeOriginal"))
    dt = dt + timedelta(seconds=subsec)

    # EXIF local time -> UTC by subtracting the offset (if present)
    off_min = _parse_offset_minutes(tags.get("OffsetTimeOriginal"))
    if off_min is None:
        # Assume EXIF time is already UTC if no offset is provided
        dt_utc = dt
    else:
        dt_utc = dt - timedelta(minutes=off_min)

    jd = _to_julian_date(dt_utc)
    return f"{jd:15.7f}"

# Helper to get UTC datetime from NEF (keeps above functions unchanged)
def nef_datetime_utc(filepath: str) -> datetime:
    tags = _run_exiftool_json(filepath)
    dto = tags.get("DateTimeOriginal")
    if not dto:
        raise RuntimeError(f"DateTimeOriginal not found in {filepath}")
    dt = _parse_exif_datetime(dto)
    subsec = _parse_subsec(tags.get("SubSecTimeOriginal"))
    dt = dt + timedelta(seconds=subsec)
    off_min = _parse_offset_minutes(tags.get("OffsetTimeOriginal"))
    if off_min is None:
        return dt.replace(tzinfo=timezone.utc)  # assume UTC if no offset tag
    else:
        return (dt - timedelta(minutes=off_min)).replace(tzinfo=timezone.utc)

# ================================
# === (Unchanged) JD -> phase/LL
# ================================

from skyfield.api import load
from skyfield.framelib import ecliptic_frame

def jd_to_datetime_utc(jd_utc: float) -> datetime:
    """Convert UTC JD to aware datetime in UTC (Meeus)."""
    J = jd_utc + 0.5
    Z = int(math.floor(J))
    F = J - Z
    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - alpha // 4
    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)
    day_float = B - D - int(30.6001 * E) + F
    day = int(day_float)
    frac = day_float - day
    if E < 14:
        month = E - 1
    else:
        month = E - 13
    if month > 2:
        year = C - 4716
    else:
        year = C - 4715
    hours = frac * 24.0
    h = int(hours)
    minutes = (hours - h) * 60.0
    m = int(minutes)
    seconds = (minutes - m) * 60.0
    s = int(seconds)
    us = int(round((seconds - s) * 1_000_000))
    if us >= 1_000_000:
        us -= 1_000_000
        s += 1
    if s >= 60:
        s -= 60
        m += 1
    if m >= 60:
        m -= 60
        h += 1
    return datetime(year, month, day, h, m, s, us, tzinfo=timezone.utc)

def wrap_deg_pm180(x: float) -> float:
    """Wrap degrees to [-180, +180)."""
    y = (x + 180.0) % 360.0 - 180.0
    if y == -180.0:
        y = 180.0
    return y

def cosd(x_deg: float) -> float:
    return math.cos(math.radians(x_deg))

def lunar_phase_and_sublunar_from_jd(jd_utc: float):
    """
    Compute lunar phase and sub-lunar geodetic lat/lon from a UTC Julian Date.

    Returns a dict with:
      illuminated_fraction, phase_angle_deg, elongation_deg, lunar_age_days,
      sublunar_lat_deg, sublunar_lon_deg
    """
    ts = load.timescale()
    dt = jd_to_datetime_utc(jd_utc)
    t = ts.utc(dt)

    eph = load('de421.bsp')   # cached after first download
    earth = eph['earth']
    moon  = eph['moon']
    sun   = eph['sun']

    e = earth.at(t)
    m_app = e.observe(moon).apparent()
    s_app = e.observe(sun).apparent()

    elong_deg = m_app.separation_from(s_app).degrees
    illum = 0.5 * (1.0 - cosd(elong_deg))
    phase_angle = 180.0 - elong_deg

    m_lon, m_lat, m_dist = m_app.frame_latlon(ecliptic_frame)
    s_lon, s_lat, s_dist = s_app.frame_latlon(ecliptic_frame)
    synodic_angle = (m_lon.degrees - s_lon.degrees) % 360.0
    SYNODIC_MONTH = 29.530588853
    _ = synodic_angle / 360.0 * SYNODIC_MONTH  # lunar_age_days (not printed here)

    ra_h, dec_deg, _ = m_app.radec()
    gast_hours = t.gast
    sublon_deg = wrap_deg_pm180((ra_h.hours - gast_hours) * 15.0)

    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    psi = math.radians(dec_deg.degrees)
    phi = math.atan(math.tan(psi) / (1.0 - e2))
    sublat_deg = math.degrees(phi)

    return {
        "illuminated_fraction": illum,
        "phase_angle_deg": phase_angle,
        "elongation_deg": elong_deg,
        "sublunar_lat_deg": sublat_deg,
        "sublunar_lon_deg": sublon_deg,
    }

# ========================
# === Glue / formatting ===
# ========================

def _looks_like_number(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False

def _iso_z(dt_utc: datetime) -> str:
    """Format UTC datetime as ISO-8601 with Z; drop .ffffff if zero."""
    s = dt_utc.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    # strip trailing zeros in fractional part (and the dot if needed)
    if ".000000Z" in s:
        return s.replace(".000000Z", "Z")
    # Otherwise remove trailing zeros but keep at least 1 digit if any
    head, tail = s.split("Z")[0], "Z"
    if "." in head:
        pre, frac = head.split(".")
        frac = frac.rstrip("0")
        head = pre if not frac else f"{pre}.{frac}"
    return head + tail

# ========================
# === CLI ===
# ========================

def main(argv=None) -> int:
    argv = sys.argv if argv is None else argv
    if len(argv) != 2:
        sys.stderr.write("Usage: JD_to_phase_illfrac_lon_lat.py <NEF_path | UTC_Julian_Date>\n")
        return 2

    arg = argv[1]

    try:
        if _looks_like_number(arg):
            # Treat as UTC JD
            jd_val = float(arg)
            jd_str = f"{jd_val:15.7f}"
            dt_utc = jd_to_datetime_utc(jd_val)
        else:
            # Treat as NEF path -> extract JD & datetime via EXIF
            jd_str = nef_jd_string(arg)       # f15.7
            jd_val = float(jd_str)
            dt_utc = nef_datetime_utc(arg)    # aware UTC datetime

        out = lunar_phase_and_sublunar_from_jd(jd_val)

    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"exiftool error: {e}\n{e.stderr}\n")
        return 1
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        return 1

    # One pipe-friendly line:
    # <JD_f15.7> <ISO8601_UTC> <phase_angle_deg> <sublunar_lon_deg> <sublunar_lat_deg> <illuminated_fraction>
    print(
        f"{jd_str} "
        f"{_iso_z(dt_utc)} "
        f"{out['phase_angle_deg']:.6f} "
        f"{out['sublunar_lon_deg']:.6f} "
        f"{out['sublunar_lat_deg']:.6f} "
        f"{out['illuminated_fraction']:.6f}"
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

