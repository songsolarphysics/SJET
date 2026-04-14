"""
Geometric Parameter Extraction for Solar Jets
==============================================
Part of SJET: An Interactive Solar Jet Extraction Tool

This module provides:
  - Quadratic Bézier curve fitting for jet axis modelling
  - Circular-region-based start/end point identification
  - Mask boundary width measurement at ten cross-sectional locations
  - Gaussian FWHM width measurement on the original image intensity
  - Manual control point override for complex or C-shaped jet morphologies

Reference:
    Tan et al. 2026, RAS Techniques and Instruments
    https://github.com/songsolarphysics/SJET
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from datetime import datetime
import sunpy.map
import astropy.units as u
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
import math
from scipy.optimize import curve_fit


# ============================================================================
# Configuration
# ============================================================================

CMAP_DEFAULT = 'sdoaia171'
DEFAULT_LINE_WIDTH = 3
DEFAULT_EXTENSION_RATIO = 0.5
DEFAULT_DPI = 300
MM_PER_PIXEL = 0.435  # Mm/pixel — adjust for your instrument


# ============================================================================
# Bézier Curve Functions
# ============================================================================

def bernstein_poly(i, n, t):
    """Bernstein basis polynomial."""
    return math.comb(n, i) * (t**i) * ((1 - t)**(n - i))


def bezier_curve(points, num_points=100):
    """
    Evaluate a Bézier curve defined by a set of control points.

    Parameters
    ----------
    points : array-like, shape (n+1, 2)
        Control points.
    num_points : int
        Number of points to evaluate on the curve.

    Returns
    -------
    curve_points : ndarray, shape (num_points, 2)
    """
    n = len(points) - 1
    t_values = np.linspace(0, 1, num_points)
    curve_points = np.zeros((num_points, 2))
    for i in range(n + 1):
        curve_points += np.outer(bernstein_poly(i, n, t_values), points[i])
    return curve_points


def fit_bezier_through_points(p0, p1, p2):
    """
    Compute the quadratic Bézier control point that passes through p1.

    The resulting curve passes through the geometric centre of the jet mask.
    """
    control_point = 2 * np.array(p1) - 0.5 * (np.array(p0) + np.array(p2))
    return control_point


def arc_length_parameterization(p0, control, p2, num_sample=1000):
    """
    Compute cumulative arc length along a quadratic Bézier curve.

    Returns
    -------
    cumulative_lengths : ndarray
    t_values : ndarray
    total_length : float
    """
    t_values = np.linspace(0, 1, num_sample)
    points = np.zeros((num_sample, 2))
    for i, t in enumerate(t_values):
        points[i] = ((1 - t)**2 * np.array(p0)
                     + 2 * (1 - t) * t * np.array(control)
                     + t**2 * np.array(p2))
    cumulative_lengths = np.zeros(num_sample)
    for i in range(1, num_sample):
        cumulative_lengths[i] = cumulative_lengths[i - 1] + np.sqrt(
            np.sum((points[i] - points[i - 1])**2))
    return cumulative_lengths, t_values, cumulative_lengths[-1]


def evaluate_bezier_curve_uniform_spacing(p0, control, p2, num_points):
    """
    Evaluate a quadratic Bézier curve with uniform arc-length spacing.
    """
    cumulative_lengths, t_values, total_length = arc_length_parameterization(
        p0, control, p2)
    uniform_t = np.interp(
        np.linspace(0, total_length, num_points),
        cumulative_lengths, t_values)
    pts = np.zeros((num_points, 2))
    for i, t in enumerate(uniform_t):
        pts[i] = ((1 - t)**2 * np.array(p0)
                  + 2 * (1 - t) * t * np.array(control)
                  + t**2 * np.array(p2))
    return pts


def evaluate_bezier_curve(p0, control, p2, num_points=100):
    """Evaluate a quadratic Bézier curve at uniformly spaced t values."""
    t = np.linspace(0, 1, num_points)
    pts = np.zeros((num_points, 2))
    for i, ti in enumerate(t):
        pts[i] = ((1 - ti)**2 * np.array(p0)
                  + 2 * (1 - ti) * ti * np.array(control)
                  + ti**2 * np.array(p2))
    return pts


# ============================================================================
# Gaussian FWHM Measurement
# ============================================================================

def _gaussian_1d(x, amplitude, center, sigma, background):
    """1-D Gaussian with a constant background baseline."""
    return background + amplitude * np.exp(-0.5 * ((x - center) / sigma)**2)


def _extract_perpendicular_profile(data, point, tangent, half_width=30):
    """
    Extract an intensity profile perpendicular to the jet axis.

    Parameters
    ----------
    data : 2-D ndarray
        Original image data (row = y, col = x).
    point : array-like, shape (2,)
        Sampling point on the curve, in (row, col) order.
    tangent : array-like, shape (2,)
        Unit tangent vector at the sampling point, in (row, col) order.
    half_width : int
        Number of pixels to sample on each side of the axis.

    Returns
    -------
    profile : 1-D ndarray
        Intensity profile of length 2*half_width + 1.
    offsets : 1-D ndarray
        Pixel offsets (negative = left, positive = right).
    """
    normal = np.array([-tangent[1], tangent[0]])
    offsets = np.arange(-half_width, half_width + 1, dtype=float)
    profile = np.full(len(offsets), np.nan)
    h, w = data.shape
    for k, off in enumerate(offsets):
        r = int(round(point[0] + off * normal[0]))
        c = int(round(point[1] + off * normal[1]))
        if 0 <= r < h and 0 <= c < w:
            profile[k] = data[r, c]
    return profile, offsets


def _fit_gaussian_fwhm(profile, offsets):
    """
    Fit a 1-D Gaussian to an intensity profile and return the FWHM in pixels.

    Returns NaN if the fit fails or there are fewer than 5 valid data points.
    """
    valid = np.isfinite(profile)
    if valid.sum() < 5:
        return np.nan, None
    x = offsets[valid]
    y = profile[valid]
    bg0  = np.percentile(y, 10)
    amp0 = np.max(y) - bg0
    try:
        popt, _ = curve_fit(
            _gaussian_1d, x, y,
            p0=[amp0, x[np.argmax(y)], (x[-1] - x[0]) / 6.0, bg0],
            bounds=([0, x[0], 0.5, 0],
                    [np.inf, x[-1], len(x), np.inf]),
            maxfev=2000)
        return 2.3548 * abs(popt[2]), popt  # FWHM = 2√(2 ln2) · σ
    except RuntimeError:
        return np.nan, None


# ============================================================================
# Main Analysis Function
# ============================================================================

def analyze_jet_circular_regions(
    file_path,
    data=None,
    visualize=True,
    save_results=False,
    output_file=None,
    control_point_override=None,
    fwhm_half_width=30
):
    """
    Analyse a solar jet binary mask using the circular-region method.

    The algorithm identifies the two most distant pixels in the mask,
    constructs equal-radius circular counting regions centred on each,
    and designates the end with more enclosed mask pixels as the jet base
    (start point). A quadratic Bézier curve is then fitted through the
    geometric centre of the mask to model the jet axis.

    Parameters
    ----------
    file_path : str
        Path to the binary mask FITS file.
    data : 2-D ndarray, optional
        Original intensity image (e.g. sunpy_map.data) used for Gaussian
        FWHM computation. If None, FWHM measurement is skipped.
    visualize : bool
        Whether to display the analysis results.
    save_results : bool
        Whether to save output files.
    output_file : str, optional
        Output file path.
    control_point_override : array-like (row, col) or None
        Manual Bézier control point. If provided, replaces the automatic
        geometric-centre-based estimate. Recommended for strongly curved,
        C-shaped, or force-merged jet structures where the automatic control
        point may deviate from the true jet spine. Both the automatic and
        manually adjusted control points are recorded in the output.
    fwhm_half_width : int
        Half-width (pixels) of the cross-sectional profile for Gaussian
        FWHM fitting. Default: 30.

    Returns
    -------
    dict with keys:
        length                : straight-line distance from start to end (px)
        curve_length          : arc length of the fitted Bézier axis (px)
        average_width         : mean boundary-based width at 10 locations (px)
        average_width_by_area : area-based width = mask area / curve length (px)
        start_point           : (row, col) of the jet base
        center_point          : (row, col) of the mask geometric centre
        end_point             : (row, col) of the jet tip
        control_point_auto    : automatically computed control point
        control_point_used    : control point actually used for fitting
        control_point_manual  : True if a manual override was applied
        curve_points          : ndarray of Bézier curve coordinates
        widths                : list of boundary widths at 10 locations
        width_positions       : indices of the 10 sampling locations
        fwhm_values           : list of per-location FWHM values (px)
        fwhm_mean             : mean Gaussian FWHM (px)
        fwhm_std              : standard deviation of FWHM across locations (px)
        rotation_angle_deg    : axis deflection angle from base to tip (degrees)
        mask_area             : total number of True pixels in the mask
    """
    # ------------------------------------------------------------------
    # 1. Load binary mask
    # ------------------------------------------------------------------
    with fits.open(file_path) as hdul:
        binary_mask = hdul[0].data.astype(bool)

    mask_coords = np.column_stack(np.where(binary_mask))
    mask_area   = int(np.sum(binary_mask))

    if len(mask_coords) == 0:
        print("Warning: mask contains no pixels.")
        return None

    # ------------------------------------------------------------------
    # 2. Find the most distant pixel pair
    # ------------------------------------------------------------------
    max_dist_sq = 0
    point1_idx = point2_idx = 0

    if len(mask_coords) < 1000:
        for i in range(len(mask_coords)):
            for j in range(i + 1, len(mask_coords)):
                d = np.sum((mask_coords[i] - mask_coords[j])**2)
                if d > max_dist_sq:
                    max_dist_sq = d
                    point1_idx, point2_idx = i, j
    else:
        # Deterministic uniform sampling for large masks (>1000 px)
        sample_size = min(1000, len(mask_coords))
        step = max(1, len(mask_coords) // sample_size)
        sample_indices = np.arange(0, len(mask_coords), step)[:sample_size]
        sampled = mask_coords[sample_indices]
        for i in range(len(sample_indices)):
            for j in range(i + 1, len(sample_indices)):
                d = np.sum((sampled[i] - sampled[j])**2)
                if d > max_dist_sq:
                    max_dist_sq = d
                    point1_idx = sample_indices[i]
                    point2_idx = sample_indices[j]

    max_distance = np.sqrt(max_dist_sq)
    point1 = mask_coords[point1_idx]
    point2 = mask_coords[point2_idx]

    # ------------------------------------------------------------------
    # 3. Circular counting regions — determine start / end points
    # ------------------------------------------------------------------
    radius = max_distance / 2
    count1 = np.sum(
        np.sqrt(np.sum((mask_coords - point1)**2, axis=1)) <= radius)
    count2 = np.sum(
        np.sqrt(np.sum((mask_coords - point2)**2, axis=1)) <= radius)

    if count1 >= count2:
        start_point, end_point = point1, point2
    else:
        start_point, end_point = point2, point1

    center_point = np.mean(mask_coords, axis=0).astype(int)

    # ------------------------------------------------------------------
    # 4. Bézier control point (with optional manual override)
    # ------------------------------------------------------------------
    control_point_auto = 2 * center_point - start_point / 2 - end_point / 2

    if control_point_override is not None:
        control_point_used   = np.array(control_point_override, dtype=float)
        control_point_manual = True
        print(f"Using manual control point: {control_point_used}  "
              f"(automatic: {control_point_auto})")
    else:
        control_point_used   = control_point_auto.copy()
        control_point_manual = False

    # ------------------------------------------------------------------
    # 5. Fit Bézier curve and compute arc length
    # ------------------------------------------------------------------
    ctrl_pts     = np.array([start_point, control_point_used, end_point])
    curve_points = bezier_curve(ctrl_pts, num_points=100)

    curve_length = float(sum(
        np.sqrt(np.sum((curve_points[i] - curve_points[i - 1])**2))
        for i in range(1, len(curve_points))))

    average_width_by_area = mask_area / curve_length if curve_length > 0 else 0.0

    # ------------------------------------------------------------------
    # 6. Mask boundary width at 10 cross-sections
    # ------------------------------------------------------------------
    num_width_samples    = 10
    width_sample_indices = np.linspace(
        0, len(curve_points) - 1, num_width_samples).astype(int)
    widths      = []
    width_lines = []

    for idx in width_sample_indices:
        pt = curve_points[idx].astype(int)
        if not (0 <= pt[0] < binary_mask.shape[0] and
                0 <= pt[1] < binary_mask.shape[1]):
            widths.append(0)
            continue

        if 0 < idx < len(curve_points) - 1:
            tan = curve_points[idx + 1] - curve_points[idx - 1]
        elif idx == 0:
            tan = curve_points[1] - curve_points[0]
        else:
            tan = curve_points[-1] - curve_points[-2]

        tan = tan / np.linalg.norm(tan)
        nor = np.array([-tan[1], tan[0]])

        left_edge = right_edge = None
        for dist in range(1, 50):
            lp = (pt + dist * nor).astype(int)
            rp = (pt - dist * nor).astype(int)
            if (0 <= lp[0] < binary_mask.shape[0] and
                    0 <= lp[1] < binary_mask.shape[1] and
                    not binary_mask[lp[0], lp[1]] and
                    left_edge is None):
                left_edge = dist
            if (0 <= rp[0] < binary_mask.shape[0] and
                    0 <= rp[1] < binary_mask.shape[1] and
                    not binary_mask[rp[0], rp[1]] and
                    right_edge is None):
                right_edge = dist
            if left_edge is not None and right_edge is not None:
                break

        if left_edge is not None and right_edge is not None:
            w = left_edge + right_edge
            widths.append(w)
            width_lines.append((pt - (w / 2) * nor,
                                 pt + (w / 2) * nor, w))
        else:
            widths.append(0)

    avg_width = (float(np.mean([w for w in widths if w > 0]))
                 if any(w > 0 for w in widths) else 0.0)

    # ------------------------------------------------------------------
    # 7. Gaussian FWHM width measurement
    # ------------------------------------------------------------------
    fwhm_values = []

    if data is not None:
        for idx in width_sample_indices:
            pt = curve_points[idx]
            if 0 < idx < len(curve_points) - 1:
                tan = curve_points[idx + 1] - curve_points[idx - 1]
            elif idx == 0:
                tan = curve_points[1] - curve_points[0]
            else:
                tan = curve_points[-1] - curve_points[-2]

            norm_t = np.linalg.norm(tan)
            if norm_t == 0:
                fwhm_values.append(np.nan)
                continue

            profile, offsets = _extract_perpendicular_profile(
                data, pt, tan / norm_t, half_width=fwhm_half_width)
            fwhm, _ = _fit_gaussian_fwhm(profile, offsets)
            fwhm_values.append(fwhm)

        valid_fwhm = [f for f in fwhm_values if np.isfinite(f)]
        fwhm_mean  = float(np.mean(valid_fwhm))           if valid_fwhm else np.nan
        fwhm_std   = float(np.std(valid_fwhm))            if len(valid_fwhm) > 1 else np.nan
    else:
        fwhm_mean = fwhm_std = np.nan
        print("Note: no original image provided (data=None); "
              "Gaussian FWHM measurement skipped.")

    # ------------------------------------------------------------------
    # 8. Axis deflection angle
    # ------------------------------------------------------------------
    s_tan = curve_points[1] - curve_points[0]
    s_tan /= np.linalg.norm(s_tan)
    e_tan = curve_points[-1] - curve_points[-2]
    e_tan /= np.linalg.norm(e_tan)

    rotation_angle_deg = float(np.degrees(
        np.arccos(np.clip(np.dot(s_tan, e_tan), -1.0, 1.0))))
    if np.cross(np.append(s_tan, 0), np.append(e_tan, 0))[2] < 0:
        rotation_angle_deg = -rotation_angle_deg

    # ------------------------------------------------------------------
    # 9. Collect results
    # ------------------------------------------------------------------
    results = {
        'length':                float(np.sqrt(np.sum((end_point - start_point)**2))),
        'curve_length':          curve_length,
        'average_width':         avg_width,
        'average_width_by_area': average_width_by_area,
        'start_point':           start_point,
        'center_point':          center_point,
        'end_point':             end_point,
        'control_point_auto':    control_point_auto,
        'control_point_used':    control_point_used,
        'control_point_manual':  control_point_manual,
        'curve_points':          curve_points,
        'widths':                widths,
        'width_positions':       width_sample_indices,
        'width_lines':           width_lines,
        'avg_width':             avg_width,
        'fwhm_values':           fwhm_values,
        'fwhm_mean':             fwhm_mean,
        'fwhm_std':              fwhm_std,
        'rotation_angle_deg':    rotation_angle_deg,
        'mask_area':             mask_area,
        'start_tangent':         s_tan,
        'end_tangent':           e_tan,
        'radius':                radius,
        'point1':                point1,
        'point2':                point2,
        'point1_count':          int(count1),
        'point2_count':          int(count2),
    }

    if visualize:
        _visualize_jet_analysis(binary_mask, results, data, fwhm_half_width)

    # ------------------------------------------------------------------
    # 10. Print summary
    # ------------------------------------------------------------------
    print("\n===== Jet Analysis Results =====")
    print(f"Start point:      ({start_point[0]}, {start_point[1]})")
    print(f"Centre:           ({center_point[0]}, {center_point[1]})")
    print(f"End point:        ({end_point[0]}, {end_point[1]})")
    cp_flag = "[manual]" if control_point_manual else "[auto]"
    print(f"Control point:    ({control_point_used[0]:.1f}, "
          f"{control_point_used[1]:.1f})  {cp_flag}")
    print(f"Boundary width:   {avg_width:.2f} px")
    print(f"Area-based width: {average_width_by_area:.2f} px")
    if np.isfinite(fwhm_mean):
        print(f"Gaussian FWHM:    {fwhm_mean:.2f} ± {fwhm_std:.2f} px")
    else:
        print("Gaussian FWHM:    unavailable")
    print(f"Rotation angle:   {rotation_angle_deg:.2f}°")
    print("================================\n")

    return results


# ============================================================================
# Visualisation
# ============================================================================

def _visualize_jet_analysis(binary_mask, results, data=None, fwhm_half_width=30):
    """
    Display the jet analysis results in a two- or three-panel figure.

    Panel 1: Binary mask with Bézier axis and mask boundary widths (cyan).
    Panel 2: Original image with Gaussian FWHM widths overlaid (yellow).
             [Only shown if data is provided.]
    Panel 3: Example intensity cross-section with Gaussian fit.
             [Only shown if data is provided.]
    """
    num_width_samples = 10
    ncols = 3 if data is not None else 2
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 6))

    cp      = results['curve_points']
    sp      = results['start_point']
    ctr     = results['center_point']
    ep      = results['end_point']
    cp_auto = results['control_point_auto']
    cp_used = results['control_point_used']
    wsi     = results['width_positions']
    p1, p2  = results['point1'], results['point2']
    radius  = results['radius']

    # Panel 1 — mask + axis
    ax = axes[0]
    ax.imshow(binary_mask, cmap='gray', origin='lower', alpha=0.5)
    ax.plot(cp[:, 1], cp[:, 0], 'r-', lw=2, label='Jet Axis')
    ax.plot(sp[1], sp[0], 'go', ms=10, label='Start (Base)')
    ax.plot(ctr[1], ctr[0], 'yo', ms=10, label='Centre')
    ax.plot(ep[1], ep[0], 'bo', ms=10, label='End (Tip)')
    ax.plot(cp_auto[1], cp_auto[0], 's', color='gray', ms=8,
            alpha=0.6, label='Control (auto)')
    ax.plot(cp_used[1], cp_used[0], 'mo', ms=8,
            label='Control (manual)' if results['control_point_manual']
            else 'Control (used)')

    for c, col in [(p1, 'g'), (p2, 'b')]:
        ax.add_patch(plt.Circle((c[1], c[0]), radius,
                                color=col, fill=False,
                                linestyle='--', alpha=0.4))

    for i, idx in enumerate(wsi):
        if i < len(results['widths']) and results['widths'][i] > 0:
            pt = cp[idx]
            tan = (cp[idx + 1] - cp[idx - 1] if 0 < idx < len(cp) - 1
                   else (cp[1] - cp[0] if idx == 0 else cp[-1] - cp[-2]))
            tan = tan / np.linalg.norm(tan)
            nor = np.array([-tan[1], tan[0]])
            w = results['widths'][i]
            s = pt - (w / 2) * nor
            e = pt + (w / 2) * nor
            ax.plot([s[1], e[1]], [s[0], e[0]], 'c-', lw=1)

    ax.set_title('Mask-based Jet Analysis\n(cyan = boundary width)')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_xlabel('Solar-X (px)')
    ax.set_ylabel('Solar-Y (px)')

    if data is not None:
        fwhm_values = results.get('fwhm_values', [])
        fwhm_mean   = results.get('fwhm_mean', np.nan)
        fwhm_std    = results.get('fwhm_std',  np.nan)

        # Panel 2 — original image + FWHM widths
        ax2 = axes[1]
        ax2.imshow(data, cmap='inferno', origin='lower',
                   norm=colors.LogNorm(
                       vmin=np.percentile(data[data > 0], 1),
                       vmax=np.percentile(data, 99.5)))
        ax2.plot(cp[:, 1], cp[:, 0], 'w-', lw=1.5, alpha=0.8)

        for i, idx in enumerate(wsi):
            if i >= len(fwhm_values) or not np.isfinite(fwhm_values[i]):
                continue
            pt = cp[idx]
            tan = (cp[idx + 1] - cp[idx - 1] if 0 < idx < len(cp) - 1
                   else (cp[1] - cp[0] if idx == 0 else cp[-1] - cp[-2]))
            tan = tan / np.linalg.norm(tan)
            nor = np.array([-tan[1], tan[0]])
            fwhm = fwhm_values[i]
            s = pt - (fwhm / 2) * nor
            e = pt + (fwhm / 2) * nor
            ax2.plot([s[1], e[1]], [s[0], e[0]], 'y-', lw=1.5,
                     label=(f'FWHM  mean = {fwhm_mean:.1f} ± '
                            f'{fwhm_std:.1f} px') if i == 0 else '')

        ax2.set_title('Gaussian FWHM Width')
        ax2.legend(fontsize=7, loc='upper right')
        ax2.set_xlabel('Solar-X (px)')
        ax2.set_ylabel('Solar-Y (px)')

        # Panel 3 — example cross-section with Gaussian fit
        ax3 = axes[2]
        mid_i   = len(wsi) // 2
        mid_idx = wsi[mid_i]
        pt = cp[mid_idx]
        tan = (cp[mid_idx + 1] - cp[mid_idx - 1]
               if 0 < mid_idx < len(cp) - 1
               else (cp[1] - cp[0] if mid_idx == 0 else cp[-1] - cp[-2]))
        tan = tan / np.linalg.norm(tan)

        profile, offsets = _extract_perpendicular_profile(
            data, pt, tan, half_width=fwhm_half_width)
        valid = np.isfinite(profile)
        ax3.plot(offsets[valid], profile[valid], 'k.-',
                 label=f'Profile (sample {mid_i + 1}/{num_width_samples})')

        if mid_i < len(fwhm_values) and np.isfinite(fwhm_values[mid_i]):
            _, popt = _fit_gaussian_fwhm(profile, offsets)
            if popt is not None:
                x_fit = np.linspace(offsets[valid][0], offsets[valid][-1], 200)
                fwhm_val = fwhm_values[mid_i]
                ax3.plot(x_fit, _gaussian_1d(x_fit, *popt), 'r-',
                         label=f'Gaussian fit\nFWHM = {fwhm_val:.1f} px')
                ax3.axhline(popt[3] + popt[0] / 2,
                            color='r', linestyle=':', alpha=0.6)
                ax3.axvspan(popt[1] - fwhm_val / 2,
                            popt[1] + fwhm_val / 2,
                            alpha=0.15, color='red')

        ax3.set_title('Example Cross-section')
        ax3.set_xlabel('Offset from axis (px)')
        ax3.set_ylabel('Intensity')
        ax3.legend(fontsize=8)
        ax3.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()
    plt.show()


# ============================================================================
# Intensity Calculation Along Bézier Curve (for time–space diagrams)
# ============================================================================

def calculate_intensity_along_bezier_extended(data, p0, p1, p2, width,
                                               extension_ratio=0,
                                               verbose=False):
    """
    Sample image intensity along a quadratic Bézier axis.

    Parameters
    ----------
    data : 2-D ndarray
        Image data.
    p0, p1, p2 : array-like
        Three control points defining the curve (start, mid, end).
    width : int
        Sampling width perpendicular to the curve (pixels).
    extension_ratio : float
        Fractional extension beyond start and end points. 0 = no extension.
    verbose : bool
        Print curve length information.

    Returns
    -------
    intensities : list of float
    curve_points : ndarray
    control_point : ndarray
    """
    control_point = fit_bezier_through_points(p0, p1, p2)
    _, _, length   = arc_length_parameterization(p0, control_point, p2)
    num_points     = int(length) + 1
    curve_points   = evaluate_bezier_curve_uniform_spacing(
        p0, control_point, p2, num_points)

    if verbose:
        print(f"Bézier axis length: {length:.1f} px")

    intensities = []
    for i in range(len(curve_points)):
        x0c = int(round(curve_points[i][0]))
        y0c = int(round(curve_points[i][1]))
        if 0 < i < len(curve_points) - 1:
            dx = curve_points[i + 1][0] - curve_points[i - 1][0]
            dy = curve_points[i + 1][1] - curve_points[i - 1][1]
        elif i == 0:
            dx = curve_points[1][0] - curve_points[0][0]
            dy = curve_points[1][1] - curve_points[0][1]
        else:
            dx = curve_points[-1][0] - curve_points[-2][0]
            dy = curve_points[-1][1] - curve_points[-2][1]

        norm = np.sqrt(dx**2 + dy**2)
        if norm > 0:
            dx, dy = dx / norm, dy / norm
        perp_dx, perp_dy = -dy, dx

        pts = [(int(round(x0c + w * perp_dx)),
                int(round(y0c + w * perp_dy)))
               for w in range(-(width // 2), width // 2 + 1)
               if 0 <= int(round(y0c + w * perp_dy)) < data.shape[0]
               and 0 <= int(round(x0c + w * perp_dx)) < data.shape[1]]
        if pts:
            intensities.append(float(np.mean([data[y, x] for x, y in pts])))

    return intensities, curve_points, control_point


# ============================================================================
# FITS File Loading Utilities
# ============================================================================

def read_fits_file(file_path):
    """Read the primary HDU data from a FITS file."""
    with fits.open(file_path) as hdul:
        return hdul[0].data


def read_and_trim_fits_files(directory, num_threads=8):
    """
    Load all FITS files in a directory and trim them to the minimum common shape.

    Returns
    -------
    all_data : ndarray, shape (n_files, ny, nx)
    trimmed_files : list of tuples (filename, original_shape, trimmed_shape)
    min_shape : tuple
    """
    fits_files = sorted(f for f in os.listdir(directory)
                        if f.endswith('.fits'))
    paths = [os.path.join(directory, f) for f in fits_files]

    with ThreadPoolExecutor(max_workers=num_threads) as ex:
        data_list = list(ex.map(read_fits_file, paths))

    shapes    = [d.shape for d in data_list]
    min_shape = tuple(min(s[i] for s in shapes) for i in range(len(shapes[0])))
    all_data  = np.empty((len(fits_files), *min_shape), dtype=np.float32)
    trimmed   = []
    lock      = Lock()

    def _proc(args):
        i, d = args
        if d.shape != min_shape:
            with lock:
                trimmed.append((fits_files[i], d.shape, min_shape))
            slices = tuple(slice(0, m) for m in min_shape)
            d = d[slices]
        all_data[i] = d

    with ThreadPoolExecutor(max_workers=num_threads) as ex:
        list(ex.map(_proc, enumerate(data_list)))

    return all_data, trimmed, min_shape
