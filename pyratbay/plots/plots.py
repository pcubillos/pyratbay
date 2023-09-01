# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'alphatize',
    'spectrum',
    'contribution',
    'temperature',
    'abundance',
    'default_colors',
]

from itertools import cycle

from cycler import cycler, Cycler
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_rgb
import numpy as np
from scipy.ndimage import gaussian_filter1d as gaussf
import mc3.plots as mp

from .. import constants as pc
from .. import spectrum as ps
from .. import tools as pt


default_colors = {
    'H2O': 'blue',
    'CO': 'xkcd:green',
    'CO2': 'red',
    'CH4': 'gold',
    'H2': 'indigo',
    'SO2': 'deepskyblue',
    'HCN': '0.6',
    'NH3': 'darkorange',
    'N2': 'darkkhaki',
    'H': 'magenta',
    'TiO': 'black',
    'VO': 'peru',
    'Na': 'darkviolet',
    'K': 'olive',
    'C2H2': 'green',
    'C2H4': 'pink',
    'He': 'dodgerblue',
    'H2S': 'cornflowerblue',
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


def spectrum(
    spectrum, wavelength, rt_path,
    data=None, uncert=None,
    bands_wl0=None, bands_flux=None, bands_response=None, bands_wl=None,
    label='model', bounds=None, logxticks=None,
    resolution=150.0, gaussbin=2.0,
    yran=None, filename=None, fignum=501, axis=None,
    ms=5.0, lw=1.25, fs=14, units=None, dpi=300, theme=None,
    ):
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
        Observing geometry: transit, eclipse, or emission.
    data: 1D float ndarray
        Observing data points at each bands_wl0.
    uncert: 1D float ndarray
        Uncertainties of the data points.
    bands_wl0: 1D float ndarray
        The mean wavelength for each band/data point.
    bands_flux: 1D float ndarray
        Band-integrated model spectrum at each bandwl.
    bands_response: Iterable of 1D float ndarrays
        Transmission response curve for each band.
    bands_wl: Iterable of 1D float ndarrays.
        The wavelength arrasy for each bands_response curve.
    label: String
        Label for spectrum curve.
    bounds: Tuple
        Tuple with -2, -1, +1, and, +2 sigma boundaries of spectrum.
        If not None, plot shaded area between +/-1sigma and +/-2sigma
        boundaries.
    logxticks: 1D float ndarray
        If not None, switch the X-axis scale from linear to log, and set
        the X-axis ticks at the locations given by logxticks.
    resolution: Float
        Binning resolution to display the spectra.
        If None, use gaussbin instead.
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
    ms: Float
        Marker sizes.
    lw: Float
        Line widths.
    fs: Float
        Font sizes.
    units: String
        Flux units. Select from: 'percent', 'ppt', 'ppm', 'none'.
    dpi: Integer
        The resolution in dots per inch for saved files.
    theme: string or mc3.plots.Theme object
        A color theme for the models.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.
    """
    # Plotting setup:
    if units is None:
        if rt_path == 'eclipse':
            units = 'ppm'
        elif rt_path == 'transit':
            units = 'percent'
        else:
            units = 'none'

    flux_scale = 1.0/pt.u(units)
    str_units = f'({units})'
    if str_units == '(none)':
        str_units = ''
    elif str_units == '(percent)':
        str_units = '(%)'

    theme = pt.resolve_theme(theme)
    if theme is None:
        theme = mp.Theme('darkorange')
        theme.light_color = 'gold'
        theme.dark_color = 'maroon'

    # Setup according to geometry:
    if rt_path == 'emission':
        ylabel = r'$F_{\rm p}$ (erg s$^{-1}$ cm$^{-2}$ cm)'
    if rt_path == 'eclipse':
        ylabel = fr'$F_{{\rm p}}/F_{{\rm s}}$ {str_units}'
    elif rt_path == 'transit':
        ylabel = fr'$(R_{{\rm p}}/R_{{\rm s}})^2$ {str_units}'

    # Bin down the spectra
    if resolution is not None:
        wl_min = np.amin(wavelength)
        wl_max = np.amax(wavelength)
        bin_wl = ps.constant_resolution_spectrum(wl_min, wl_max, resolution)
        bin_model = ps.bin_spectrum(bin_wl, wavelength, spectrum)
        if bounds is not None:
            bin_bounds = [
                ps.bin_spectrum(bin_wl, wavelength, bound)
                for bound in bounds
            ]
    else:
        bin_wl = wavelength
        bin_model = gaussf(spectrum, gaussbin)
        if bounds is not None:
            bin_bounds = [gaussf(bound, gaussbin) for bound in bounds]


    # The plot
    if axis is None:
        fig = plt.figure(fignum)
        fig.set_size_inches(8.5, 4.5)
        plt.clf()
        ax = plt.subplot(111)
    else:
        ax = axis

    if bounds is not None:
        ax.fill_between(
            bin_wl, flux_scale*bin_bounds[2], flux_scale*bin_bounds[3],
            facecolor=theme.light_color, edgecolor='none', alpha=0.5,
        )
        ax.fill_between(
            bin_wl, flux_scale*bin_bounds[0], flux_scale*bin_bounds[1],
            facecolor=theme.light_color, edgecolor='none', alpha=0.75,
        )
    plt.plot(
        bin_wl, bin_model*flux_scale, lw=lw, color=theme.color, label=label
    )
    # Plot band-integrated model:
    if bands_flux is not None and bands_wl0 is not None:
        plt.plot(
            bands_wl0, bands_flux*flux_scale,
            ls='', marker='o', ms=ms, mew=lw,
            color=theme.color, mec=theme.dark_color,
        )
    # Plot data:
    if data is not None and uncert is not None and bands_wl0 is not None:
        plt.errorbar(
            bands_wl0, data*flux_scale, uncert*flux_scale,
            fmt='o', label='data',
            mfc='0.45', mec='black', ecolor='0.2',
            ms=ms, elinewidth=lw, capthick=lw, zorder=3,
        )

    if yran is not None:
        ax.set_ylim(np.array(yran))
    yran = ax.get_ylim()

    xmin = np.amin(wavelength)
    xmax = np.amax(wavelength)
    is_log = logxticks is not None
    def color(x, is_log):
        if is_log:
            return np.log(x/xmin) / np.log(xmax/xmin)
        else:
            return (x-xmin) / (xmax-xmin)

    # Transmission filters:
    if bands_response is not None and bands_wl is not None:
        band_height = 0.05*(yran[1] - yran[0])
        for response, wl, wl0 in zip(bands_response, bands_wl, bands_wl0):
            col = plt.cm.viridis_r(color(wl0, is_log))
            btrans = band_height * response/np.amax(response)
            plt.plot(wl, yran[0]+btrans, color=col, lw=1.0, zorder=-10)
        ax.set_ylim(yran)

    if is_log:
        ax.set_xscale('log')
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticks(logxticks)

    ax.tick_params(
        which='both', right=True, top=True, direction='in', labelsize=fs-2,
    )
    ax.set_xlabel('Wavelength (um)', fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.legend(loc='best', numpoints=1, fontsize=fs-1)
    ax.set_xlim(xmin, xmax)
    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, dpi=dpi)
    return ax


def contribution(
        contrib_func, wl, rt_path, pressure,
        filename=None, filters=None, fignum=-21, dpi=300,
    ):
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
    filename: String
        Filename of the output figure.
    filters: 1D string ndarray
        Name of the filter bands (optional).
    fignum: Integer
        Figure number.
    dpi: Integer
        The resolution in dots per inch for saved files.

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

    is_emission = rt_path in pc.emission_rt
    is_transit = rt_path in pc.transmission_rt

    press = pressure / pc.bar
    p_ranges = np.amax(press), np.amin(press)
    log_p_ranges = np.log10(p_ranges)

    if is_emission:
        xlabel = 'contribution functions'
        cbtop  = 0.5
    elif is_transit:
        xlabel = 'transmittance'
        cbtop  = 0.8
    else:
        rt_paths = pc.rt_paths
        print(f"Invalid radiative-transfer geometry. Select from: {rt_paths}")
        return

    fs = 12
    colors = np.asarray(np.linspace(0, 255, nfilters), int)
    # 68% percentile boundaries of the central cumulative function:
    lo = 0.5*(1-0.683)
    hi = 1.0 - lo
    # Filter fontsize and thinning:
    ffs = 8.0 + (nfilters<50) + (nfilters<65)
    thin = (nfilters>80) + (nfilters>125) + (nfilters<100) + nfilters//100

    # Colormap and percentile limits:
    zz = contrib_func / np.amax(contrib_func)
    z = np.empty((nfilters, nlayers, 4), dtype=float)
    plo = np.zeros(nfilters+1)
    phi = np.zeros(nfilters+1)
    for i in range(nfilters):
        z[i] = plt.cm.rainbow(colors[i])
        z[i,:,-1] = zz[:,i]**(0.5+0.5*(is_transit))
        if is_emission:
            cumul = np.cumsum(zz[:,i])/np.sum(zz[:,i])
            plo[i], phi[i] = press[cumul>lo][0], press[cumul>hi][0]
        elif is_transit:
            plo[i], phi[i] = press[zz[:,i]<lo][0], press[zz[:,i]<hi][0]
    plo[-1] = plo[-2]
    phi[-1] = phi[-2]
    log_p_lo = np.log10(plo)
    log_p_hi = np.log10(phi)

    fig = plt.figure(fignum)
    fig.set_size_inches(8.5, 4.5)
    plt.clf()
    plt.subplots_adjust(0.09, 0.10, 0.9, 0.95)
    ax = plt.subplot(111)
    ax.imshow(
        z.swapaxes(0,1),
        aspect='auto',
        extent=[0, nfilters, log_p_ranges[0], log_p_ranges[1]],
        origin='upper',
        interpolation='nearest',
    )
    ax.plot(log_p_lo, drawstyle='steps-post', color='0.25', lw=0.75, ls='--')
    ax.plot(log_p_hi, drawstyle='steps-post', color='0.25', lw=0.75, ls='--')
    ax.yaxis.set_visible(False)
    ax.set_xticklabels([])
    ax.set_xlabel(f'Band-averaged {xlabel}', fontsize=fs)
    ax.tick_params(which='both', top=True, direction='in', labelsize=fs-2)
    ax.set_xlim(0, nfilters)
    ax.set_ylim(log_p_ranges)

    # ax works in log(p) space because imshow only works with linear axes
    # pax shows the pressure axis as intended in log units
    pax = ax.twinx()
    pax.spines['left'].set_visible(True)
    pax.yaxis.set_label_position('left')
    pax.yaxis.set_ticks_position('left')
    pax.set_ylim(p_ranges)
    pax.set_yscale('log')
    pax.set_ylabel(r'Pressure (bar)', fontsize=fs)
    pax.tick_params(which='both', right=True, direction='in', labelsize=fs-2)

    # Bandpass names/wavelengths:
    for i in range(0, nfilters-thin//2, thin):
        fname = f' {wl[i]:5.2f} um '
        # Strip root and file extension:
        if filters is not None:
            fname = str(filters[i]) + ' @' + fname
        ax.text(
            i+0.1, log_p_ranges[1], fname,
            rotation=90, ha='left', va='top', fontsize=ffs,
        )

    # Color bar:
    cbar = plt.axes([0.912, 0.10, 0.020, 0.85])
    cz = np.zeros((100, 2, 4), dtype=float)
    cz[:,0,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(is_transit))
    cz[:,1,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(is_transit))
    cbar.imshow(
        cz, aspect='auto', extent=[0, 1, 0, 1],
        origin='lower', interpolation='nearest',
    )
    if is_transit:
        cbar.axhline(0.1585, color='k', lw=1.0, dashes=(2.5,1))
        cbar.axhline(0.8415, color='w', lw=1.0, dashes=(2.5,1))
    cbar.spines['right'].set_visible(True)
    cbar.yaxis.set_label_position('right')
    cbar.yaxis.set_ticks_position('right')
    cbar.set_ylabel(xlabel.capitalize(), fontsize=fs)
    cbar.xaxis.set_visible(False)
    cbar.axes.tick_params(which='both', direction='in', labelsize=fs-2)

    if filename is not None:
        plt.savefig(filename, dpi=dpi)
    return ax


def temperature(
        pressure, profiles=None, labels=None, colors=None,
        bounds=None, punits='bar', ax=None, filename=None,
        theme='blue', alpha=[0.75,0.5], fs=13, lw=2.0, fignum=504,
        dpi=300,
    ):
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
    theme: string or mc3.plots.Theme object
        A color theme for the profiles and credible regions.
    alpha: 2-element float iterable
        Alpha transparency for bounds regions.
    fs: Float
        Labels font sizes.
    lw: Float
        Lines width.
    fignum: Integer
        Figure's number (ignored if axis is not None).
    dpi: Integer
        The resolution in dots per inch for saved files.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.
    """
    press = pressure / pt.u(punits)

    if profiles is None:
        profiles = []
    if np.ndim(profiles) == 1 and len(profiles) == len(pressure):
        profiles = [profiles]
    nprofiles, nlayers = np.shape(profiles)

    theme = pt.resolve_theme(theme)
    if theme is None:
        theme = mp.THEMES['blue']

    if colors is None and nprofiles <= 2:
        colors = [theme.color, theme.dark_color]
    elif colors is None:
        c = cycle(default_colors.values())
        colors = [next(c) for _ in profiles]

    if labels is None:
        _labels = [None for _ in profiles]
    else:
        _labels = labels

    tighten = ax is None
    if ax is None:
        plt.figure(fignum, (7,5))
        plt.clf()
        ax = plt.subplot(111)

    # Note alpha != 0 does not work for ps/eps figures
    if bounds is not None and len(bounds) == 4:
        ax.fill_betweenx(
            press, bounds[2], bounds[3],
            facecolor=theme.light_color, edgecolor='none', alpha=alpha[1],
        )
    if bounds is not None and len(bounds) >= 2:
        ax.fill_betweenx(
            press, bounds[0], bounds[1],
            facecolor=theme.light_color, edgecolor='none', alpha=alpha[0],
        )

    for profile, color, label in zip(profiles, colors, _labels):
        plt.plot(profile, press, color=color, lw=lw, label=label)

    ax.set_ylim(np.amax(press), np.amin(press))
    ax.set_yscale('log')
    plt.xlabel('Temperature (K)', fontsize=fs)
    plt.ylabel(f'Pressure ({punits})', fontsize=fs)
    ax.tick_params(
        which='both', right=True, top=True, direction='in', labelsize=fs-2,
    )
    if labels is not None:
        plt.legend(loc='best', fontsize=fs-2)
    if tighten:
        plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=dpi)
    return ax


def abundance(
        vol_mix_ratios, pressure, species,
        highlight=None, xlim=None, punits='bar',
        colors=None, dashes=None, filename=None,
        lw=2.0, fignum=505, fs=13, legend_fs=None, ax=None, dpi=300,
    ):
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
    dpi: Integer
        The resolution in dots per inch for saved files.

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
    >>> vmr = pa.chemistry('tea', pressure, temperature, species).vmr
    >>> ax = pp.abundance(
    >>>     vmr, pressure, species, colors='default',
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
        cols = [
            default_colors[mol] if mol in default_colors
            else None
            for mol in species
        ]
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
        ax.loglog(
            vol_mix_ratios[:,imol], press, label=spec, lw=lw,
            color=cols[imol], dashes=dashes[imol],
        )
    if xlim is None:
        xlim = ax.get_xlim()
    for spec in lowlight:
        imol = list(species).index(spec)
        ax.loglog(
            vol_mix_ratios[:,imol], press, label=spec, lw=lw, zorder=-1,
            color=alphatize(cols[imol],alpha=0.4), dashes=dashes[imol],
        )
    ax.set_xlim(xlim)
    ax.set_ylim(np.amax(press), np.amin(press))
    ax.set_xlabel('Volume mixing ratio', fontsize=fs)
    ax.set_ylabel(f'Pressure ({punits})', fontsize=fs)
    ax.tick_params(
        which='both', right=True, top=True, direction='in', labelsize=fs-2,
    )
    if legend_fs > 0:
        ax.legend(loc='best', fontsize=legend_fs)

    if filename is not None:
        plt.savefig(filename, dpi=dpi)
    return ax
