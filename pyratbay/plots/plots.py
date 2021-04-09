# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'alphatize',
    'spectrum',
    'contribution',
    'temperature',
    'abundance',
    'default_colors',
    ]

import os
from itertools import cycle

from cycler import cycler, Cycler
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_rgb
import numpy as np
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf

from .. import constants as pc
from .. import tools as pt


default_colors = {
    'H2O':"navy",
    'CO2':"red",
    'CO':"limegreen",
    'CH4':"orange",
    'H2':"deepskyblue",
    'He':"seagreen",
    'HCN':"0.7",
    'NH3':"magenta",
    'C2H2':"brown",
    'C2H4':"pink",
    'N2':"gold",
    'H':"olive",
    'TiO':"black",
    'VO':"peru",
    'Na':"darkviolet",
    'K':"cornflowerblue",
    }


def alphatize(colors, alpha, bg='w'):
    """
    Get rgb representation of a color as if it had the specified alpha.

    Parameters
    ----------
    colors: color or iterable of colors
        The color to alphatize.
    alpha: Float
        Alpha value to apply.
    bg: color
        Background color.

    Returns
    -------
    rgb: RGB or list of RGB color arrays
        The RGB representation of the alphatized color (or list of colors).

    Examples
    --------
    >>> import pyrabay.plots as pp
    >>> pp.alphatize('r', 0.5)
    array([1. , 0.5, 0.5])
    >>> pp.alphatize(['r', 'b'], 0.8)
    [array([1. , 0.2, 0.2]), array([0.2, 0.2, 1. ])]
    """
    flatten = False
    if is_color_like(colors):
        colors = [colors]
        flatten = True
    colors = [np.array(to_rgb(color)) for color in colors]
    bg = np.array(to_rgb(bg))

    # https://matplotlib.org/tutorials/colors/colors.html
    rgb = [(1.0-alpha) * bg + alpha*c for c in colors]

    if flatten:
        return rgb[0]
    return rgb


def spectrum(spectrum, wavelength, rt_path,
    data=None, uncert=None, bandwl=None, bandflux=None,
    bandtrans=None, bandidx=None,
    starflux=None, rprs=None, label='model', bounds=None,
    logxticks=None,
    gaussbin=2.0, yran=None, filename=None, fignum=501, axis=None):
    """
    Plot a transmission or emission model spectrum with (optional) data
    points with error bars and band-integrated model.

    Parameters
    ----------
    spectrum: 1D float ndarray
        Planetary spectrum evaluated at wavelength.
    wavelength: 1D float ndarray
        The wavelength of the model in microns.
    rt_path: String
        Radiative-transfer observing geometry (transit, eclipse, or emission).
    data: 1D float ndarray
        Observing data points at each bandwl.
    uncert: 1D float ndarray
        Uncertainties of the data points.
    bandwl: 1D float ndarray
        The mean wavelength for each band/data point.
    bandflux: 1D float ndarray
        Band-integrated model spectrum at each bandwl.
    bandtrans: List of 1D float ndarrays
        Transmission curve for each band.
    bandidx: List of 1D float ndarrays.
        The indices in wavelength for each bandtrans.
    starflux: 1D float ndarray
        Stellar spectrum evaluated at wavelength.
    rprs: Float
        Planet-to-star radius ratio.
    label: String
        Label for spectrum curve.
    bounds: Tuple
        Tuple with -2, -1, +1, and, +2 sigma boundaries of spectrum.
        If not None, plot shaded area between +/-1sigma and +/-2sigma
        boundaries.
    logxticks: 1D float ndarray
        If not None, switch the X-axis scale from linear to log, and set
        the X-axis ticks at the locations given by logxticks.
    gaussbin: Integer
        Standard deviation for Gaussian-kernel smoothing (in number of samples).
    yran: 1D float ndarray
        Figure's Y-axis boundaries.
    filename: String
        If not None, save figure to filename.
    fignum: Integer
        Figure number.
    axis: AxesSubplot instance
        The matplotlib Axes of the figure.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.
    """
    # Plotting setup:
    fs = 14.0
    ms = 6.0
    lw = 1.25

    if axis is None:
        plt.figure(fignum, (8, 5))
        plt.clf()
        ax = plt.subplot(111)
    else:
        ax = axis

    #fscale = {'':1.0, '%':100.0, 'ppt':1e3, 'ppm':1e6}

    spec_kw = {'label':label}
    if bounds is None:
        spec_kw['color'] = 'orange'
    else:
        spec_kw['color'] = 'orangered'


    # Setup according to geometry:
    if rt_path == 'emission':
        fscale = 1.0
        plt.ylabel(r'$F_{\rm p}$ (erg s$^{-1}$ cm$^{-2}$ cm)', fontsize=fs)
    if rt_path == 'eclipse':
        #if starflux is not None and rprs is not None:
        spectrum = spectrum/starflux * rprs**2.0
        if bounds is not None:
            bounds = [bound/starflux * rprs**2.0 for bound in bounds]
        fscale = 1.0 / pc.ppt
        plt.ylabel(r'$F_{\rm p}/F_{\rm s}\ (ppt)$', fontsize=fs)
    elif rt_path == 'transit':
        fscale = 1.0 / pc.percent
        plt.ylabel(r'$(R_{\rm p}/R_{\rm s})^2$ (%)', fontsize=fs)

    gmodel = gaussf(spectrum, gaussbin)
    if bounds is not None:
        gbounds = [gaussf(bound, gaussbin) for bound in bounds]
        ax.fill_between(
            wavelength, fscale*gbounds[0], fscale*gbounds[3],
            facecolor='gold', edgecolor='none')
        ax.fill_between(
            wavelength, fscale*gbounds[1], fscale*gbounds[2],
            facecolor='orange', edgecolor='none')

    # Plot model:
    plt.plot(wavelength, gmodel*fscale, lw=lw, **spec_kw)
    # Plot band-integrated model:
    if bandflux is not None and bandwl is not None:
        plt.plot(
            bandwl, bandflux*fscale, 'o', ms=ms, color='tomato',
            mec='maroon', mew=lw)
    # Plot data:
    if data is not None and uncert is not None and bandwl is not None:
        plt.errorbar(
            bandwl, data*fscale, uncert*fscale, fmt='o', label='data',
            color='blue', ms=ms, elinewidth=lw, capthick=lw, zorder=3)

    if yran is not None:
        ax.set_ylim(np.array(yran))
    yran = ax.get_ylim()

    # Transmission filters:
    if bandtrans is not None and bandidx is not None:
        bandh = 0.06*(yran[1] - yran[0])
        for btrans, bidx in zip(bandtrans, bandidx):
            btrans = bandh * btrans/np.amax(btrans)
            plt.plot(wavelength[bidx], yran[0]+btrans, '0.4', zorder=-10)
        ax.set_ylim(yran)

    if logxticks is not None:
        ax.set_xscale('log')
        plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticks(logxticks)

    ax.tick_params(
        which='both', right=True, top=True, direction='in', labelsize=fs-2)
    plt.xlabel('Wavelength (um)', fontsize=fs)
    plt.legend(loc='best', numpoints=1, fontsize=fs-1)
    plt.xlim(np.amin(wavelength), np.amax(wavelength))
    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename)
    return ax


def contribution(contrib_func, wl, rt_path, pressure, radius, rtop=0,
    filename=None, filters=None, fignum=-21):
    """
    Plot the band-integrated normalized contribution functions
    (emission) or transmittance (transmission).

    Parameters
    ----------
    contrib_func: 2D float ndarray
        Band-integrated contribution functions [nfilters, nlayers].
    wl: 1D float ndarray
        Mean wavelength of the bands in microns.
    rt_path: String
        Radiative-transfer observing geometry (emission or transit).
    pressure: 1D float ndarray
        Layer's pressure array (barye units).
    radius: 1D float ndarray
        Layer's impact parameter array (cm units).
    rtop: Integer
        Index of topmost valid layer.
    filename: String
        Filename of the output figure.
    filters: 1D string ndarray
        Name of the filter bands (optional).
    fignum: Integer
        Figure number.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

    Notes
    -----
    - The dashed lines denote the 0.16 and 0.84 percentiles of the
      cumulative contribution function or the transmittance (i.e.,
      the boundaries of the central 68% of the respective curves).
    - If there are more than 80 filters, this code will thin the
      displayed filter names.
    """
    nfilters = len(wl)
    nlayers = len(pressure)

    wlsort = np.argsort(wl)
    wl = wl[wlsort]
    contrib_func = contrib_func[:,wlsort]
    if filters is not None:
        filters = [filters[i] for i in wlsort]

    press = pressure[rtop:]/pc.bar
    rad = radius[rtop:]/pc.km

    press = pressure[rtop:]/pc.bar
    rad = radius[rtop:]/pc.km
    zz = contrib_func/np.amax(contrib_func)

    is_emission = rt_path == 'emission'
    is_transit = rt_path == 'transit'

    if is_emission:
        yran = np.amax(np.log10(press)), np.amin(np.log10(press))
        xlabel = 'contribution function'
        ylabel = ''
        yright = 0.9
        cbtop  = 0.5
    elif is_transit:
        yran = np.amin(rad), np.amax(rad)
        xlabel = 'transmittance'
        ylabel = 'Impact parameter (km)'
        yright = 0.84
        cbtop  = 0.8
    else:
        print(
            "Invalid radiative-transfer geometry.  "
            "Select from: 'emission' or 'transit'.")
        return

    fs  = 12
    colors = np.asarray(np.linspace(0, 255, nfilters), np.int)
    # 68% percentile boundaries of the central cumulative function:
    lo = 0.5*(1-0.683)
    hi = 1.0 - lo
    # Filter fontsize and thinning:
    ffs = 8.0 + (nfilters<50) + (nfilters<65)
    thin = (nfilters>80) + (nfilters>125) + (nfilters<100) + nfilters//100

    # Colormap and percentile limits:
    z = np.empty((nfilters, nlayers, 4), dtype=float)
    plo = np.zeros(nfilters+1)
    phi = np.zeros(nfilters+1)
    for i in np.arange(nfilters):
        z[i] = plt.cm.rainbow(colors[i])
        z[i,:,-1] = zz[:,i]**(0.5+0.5*(is_transit))
        if is_emission:
            cumul = np.cumsum(zz[:,i])/np.sum(zz[:,i])
            plo[i], phi[i] = press[cumul>lo][0], press[cumul>hi][0]
        elif is_transit:
            plo[i], phi[i] = press[zz[:,i]<lo][0], press[zz[:,i]<hi][0]
    plo[-1] = plo[-2]
    phi[-1] = phi[-2]

    fig = plt.figure(fignum, (8.5, 5))
    plt.clf()
    plt.subplots_adjust(0.105, 0.10, yright, 0.95)
    ax = plt.subplot(111)
    pax = ax.twinx()
    if is_emission:
        ax.imshow(
            z.swapaxes(0,1), aspect='auto',
            extent=[0, nfilters, yran[0], yran[1]],
            origin='upper', interpolation='nearest')
        ax.yaxis.set_visible(False)
        pax.spines['left'].set_visible(True)
        pax.yaxis.set_label_position('left')
        pax.yaxis.set_ticks_position('left')
    elif is_transit:
        ax.imshow(
            z.swapaxes(0,1), aspect='auto',
            extent=[0,nfilters,yran[0],yran[1]],
            origin='upper', interpolation='nearest')
        # Setting the right radius tick labels requires some sorcery:
        fig.canvas.draw()
        ylab = [l.get_text() for l in ax.get_yticklabels()]
        rint = si.interp1d(rad, press, bounds_error=False)
        pticks = rint(ax.get_yticks())
        bounds = np.isfinite(pticks)
        pint = si.interp1d(
            press, np.linspace(yran[1], yran[0], nlayers), bounds_error=False)
        ax.set_yticks(pint(pticks[bounds]))
        ax.set_yticklabels(np.array(ylab)[bounds])

    pax.plot(plo, drawstyle='steps-post', color='0.25', lw=0.75, ls='--')
    pax.plot(phi, drawstyle='steps-post', color='0.25', lw=0.75, ls='--')
    pax.set_ylim(np.amax(press), np.amin(press))
    pax.set_yscale('log')
    pax.set_ylabel(r'Pressure (bar)', fontsize=fs)

    ax.set_xlim(0, nfilters)
    ax.set_ylim(yran)
    ax.set_xticklabels([])
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.set_xlabel(f'Band-averaged {xlabel}', fontsize=fs)

    # Print filter names/wavelengths:
    for i in np.arange(0, nfilters-thin//2, thin):
        fname = f' {wl[i]:5.2f} um '
        # Strip root and file extension:
        if filters is not None:
            fname = (os.path.split(os.path.splitext(filters[i])[0])[1]
                     + ' @' + fname)
        ax.text(
            i+0.1, yran[1], fname, rotation=90, ha='left', va='top',
            fontsize=ffs)

    # Color bar:
    cbar = plt.axes([0.925, 0.10, 0.015, 0.85])
    cz = np.zeros((100, 2, 4), dtype=float)
    cz[:,0,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(is_transit))
    cz[:,1,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(is_transit))
    cbar.imshow(
        cz, aspect='auto', extent=[0, 1, 0, 1],
        origin='lower', interpolation='nearest')
    if is_transit:
        cbar.axhline(0.1585, color='k', lw=1.0, dashes=(2.5,1))
        cbar.axhline(0.8415, color='w', lw=1.0, dashes=(2.5,1))
    cbar.spines['right'].set_visible(True)
    cbar.yaxis.set_label_position('right')
    cbar.yaxis.set_ticks_position('right')
    cbar.set_ylabel(xlabel.capitalize(), fontsize=fs)
    cbar.xaxis.set_visible(False)

    fig.canvas.draw()
    if filename is not None:
        plt.savefig(filename)
    return ax


def temperature(pressure, profiles=None, labels=None, colors=None,
    bounds=None, punits='bar', ax=None, filename=None,
    theme='blue', alpha=[0.8,0.6], fs=13, lw=2.0, fignum=504):
    """
    Plot temperature profiles.

    Parameters
    ----------
    pressure: 1D float ndarray
        The atmospheric pressure profile in barye.
    profiles: iterable of 1D float ndarrays
        Temperature profiles to plot.
    labels: 1D string iterable
        Labels for temperature profiles.
    colors: 1D string iterable.
        Colors for temperature profiles.
    bounds: Tuple
        Tuple with -1sigma, +1sigma, -2sigma, and +2sigma temperature
        boundaries.
        If not None, plot shaded area between +/-1sigma and +/-2sigma
        boundaries.
    punits: String
        Pressure units for output plot (input units are always barye).
    ax: AxesSubplot instance
        If not None, plot into the given axis.
    filename: String
        If not None, save plot to given file name.
    theme: String
        The histograms' color theme for bounds regions.
        Only 'blue' and 'orange' themes are valid at the moment.
    alpha: 2-element float iterable
        Alpha transparency for bounds regions.
    fs: Float
        Labels font sizes.
    lw: Float
        Lines width.
    fignum: Integer
        Figure's number (ignored if axis is not None).

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.
    """
    press = pressure / pt.u(punits)

    if theme == 'blue':
        col1, col2 = 'royalblue', 'royalblue'
    elif theme == 'orange':
        col1, col2 = 'orange', 'gold'
    # alpha != 0 does not work for ps/eps figures:
    alpha1, alpha2 = alpha[:]
    if filename is not None and filename.endswith('ps'):
        fc2 = alphatize(col2, alpha2, 'white')
        fc1 = alphatize(col1, alpha1, fc2)
        alpha1 = alpha2 = 1.0
    else:
        fc1, fc2 = col1, col2

    if profiles is None:
        profiles = []
    if np.ndim(profiles) == 1 and len(profiles) == len(pressure):
        profiles = [profiles]

    if labels is None:
        _labels = [None for _ in profiles]
    else:
        _labels = labels

    if colors is None:
        c = cycle(default_colors.values())
        colors = [next(c) for _ in profiles]

    if ax is None:
        plt.figure(fignum, (7,5))
        plt.clf()
        ax = plt.subplot(111)

    if bounds is not None and len(bounds) == 4:
        low2, high2 = bounds[2:4]
        ax.fill_betweenx(
            press, low2, high2, facecolor=fc2, edgecolor='none', alpha=alpha2)
    if bounds is not None and len(bounds) >= 2:
        low1, high1 = bounds[0:2]
        ax.fill_betweenx(
            press, low1, high1, facecolor=fc1, edgecolor='none', alpha=alpha1)

    for profile, color, label in zip(profiles, colors, _labels):
        plt.plot(profile, press, color, lw=lw, label=label)

    ax.set_ylim(np.amax(press), np.amin(press))
    ax.set_yscale('log')
    plt.xlabel('Temperature (K)', fontsize=fs)
    plt.ylabel(f'Pressure ({punits})', fontsize=fs)
    ax.tick_params(labelsize=fs-2)
    if labels is not None:
        plt.legend(loc='best', fontsize=fs-2)
    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename)
    return ax


def abundance(vol_mix_ratios, pressure, species,
    highlight=None, xlim=None, punits='bar',
    colors=None, dashes=None, filename=None,
    lw=2.0, fignum=505, fs=13, legend_fs=None, ax=None):
    """
    Plot atmospheric volume-mixing-ratio abundances.

    Parameters
    ----------
    vol_mix_ratios: 2D float ndarray
        Atmospheric volume mixing ratios to plot [nlayers,nspecies].
    pressure: 1D float ndarray
        Atmospheric pressure [nlayers], units are given by punits argument.
    species: 1D string iterable
        Atmospheric species names [nspecies].
    highlight: 1D string iterable
        List of species names to highlight.  Non-highlighed species are
        plotted with alpha=0.4, below the highligted species, and are
        not considered to set the default xlim (e.g., might not be shown
        if their abundances are too low).
        If None, all input species are highlighted.
    xlim: 2-element float iterable
        Volume mixing ratio plotting boundaries.
    punits: String
        Pressure units.
    colors: 1D string iterable
        List of colors to use.
        - If len(colors) >= len(species), colors are assigned to each
          species irrespective of highlight.
        - If len(colors) < len(species), the display will cycle the
          color list using solid, long-dashed, short-dashed, and dotted
          line styles (all highlight species being displayed before the rest).
        - If colors == 'default', use pyratbay.plots.default_colors
          dict to assign colors.
        - If colors is None, use matplotlib's default color cycler.
    dashes: 1D dash-sequence iterable
        List of line-styles for each species, irrespective of highlight.
        len(dashes) has to be equal to len(species).
        Alternatively, dashes can by a dash-sequence Cycler.
    filename: String
        If not None, save plot to given file name.
    lw: Float
        Lines width.
    fignum: Integer
        Figure's number (ignored if axis is not None).
    fs: Float
        Labels font sizes.
    legend_fs: Float
        Legend font size.  If legend_fs is None, default to fs-2.
        If legend_fs <= 0, do not plot a legend.
    ax: AxesSubplot instance
        If not None, plot into the given axis.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.plots as pp

    >>> nlayers = 51
    >>> pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers)
    >>> temperature = pa.temperature('isothermal', pressure,  params=1000.0)
    >>> species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 TiO VO H2 H He Na K'.split()
    >>> Q = pa.abundance(pressure, temperature, species, ncpu=3)
    >>> ax = pp.abundance(Q, pressure, species, colors='default',
    >>>     highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())
    """
    if legend_fs is None:
        legend_fs = fs - 2

    if highlight is None:
        highlight = np.copy(species)
    highlight = [spec for spec in species if spec in highlight]
    lowlight  = [spec for spec in species if spec not in highlight]
    sorted_spec = highlight + lowlight

    used_cols = []
    if colors is None:
        colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']

    if len(colors) >= len(species):
        cols = colors
    elif colors == 'default':
        cols = [default_colors[mol] if mol in default_colors
                else None
                for mol in species]
        used_cols = [c for c in default_colors.values() if c in cols]
        remaining_cols = [c for c in default_colors.values() if c not in cols]
        colors = used_cols + remaining_cols
    else:
        cols = [None for _ in species]

    if isinstance(dashes, Cycler):
        dash_cycler = dashes
        dashes = None
    else:
        dash_cycler = cycler(dashes=[(), (8,1.5), (4,1), (1,1)])

    dkws = cycle(dash_cycler * cycler(color=colors))
    for _ in used_cols:
        dkw = next(dkws)

    _dashes = [() for _ in species]
    for i in range(len(species)):
        ispec = list(species).index(sorted_spec[i])
        if cols[ispec] is None:
            dkw = next(dkws)
            cols[ispec] = dkw['color']
            _dashes[ispec] = dkw['dashes']
    if dashes is None or len(dashes) != len(species):
        dashes = _dashes

    press = pressure / pt.u(punits)
    # Plot the results:
    if ax is None:
        plt.figure(fignum, (7,5))
        plt.clf()
        ax = plt.subplot(111)
    for spec in highlight:
        imol = list(species).index(spec)
        ax.loglog(vol_mix_ratios[:,imol], press, label=spec, lw=lw,
            color=cols[imol], dashes=dashes[imol])
    if xlim is None:
        xlim = ax.get_xlim()
    for spec in lowlight:
        imol = list(species).index(spec)
        ax.loglog(
            vol_mix_ratios[:,imol], press, label=spec, lw=lw, zorder=-1,
            color=alphatize(cols[imol],alpha=0.4), dashes=dashes[imol])
    ax.set_xlim(xlim)
    ax.set_ylim(np.amax(press), np.amin(press))
    ax.set_xlabel('Volume mixing ratio', fontsize=fs)
    ax.set_ylabel(f'Pressure ({punits})', fontsize=fs)
    ax.tick_params(labelsize=fs-2)
    if legend_fs > 0:
        ax.legend(loc='best', fontsize=legend_fs)

    if filename is not None:
        plt.savefig(filename)
    return ax

