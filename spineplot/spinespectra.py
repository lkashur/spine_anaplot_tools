import numpy as np
import pandas as pd
from style import Style
from variable import Variable
import matplotlib.pyplot as plt

class ConfigException(Exception):
    pass

class SpineSpectra:
    """
    A base class designed to encapsulate spectra for multiple variables
    in an ensemble of samples. Though intended for use with a
    collection of variables, this class is meant to be an abstraction
    representing a single plot instance or plot content. 

    Attributes
    ----------
    _style : str
        The name of the style sheet to use for the plot.
    _variables : list
        The list of Variable objects for the spectra.
    _categories : dict
        A dictionary of the categories for the spectra. This serves as
        a map between the category label in the input TTree and the
        category label for the aggregated data (and therefore what is
        shown in a single legend entry).
    _plotdata : dict
        A dictionary of the data for the spectra. This is a map between
        the category label for the spectra and the histogram data for
        that category.
    """
    def __init__(self, style, variables, categories, colors) -> None:
        """
        Initializes the SpineSpectra object with the given kwargs.

        Parameters
        ----------
        style : str
            The name of the style sheet to use for the plot.
        variables : list
            The list of Variable objects for the spectra.
        categories : dict
            A dictionary of the categories for the spectra. This serves
            as a map between the category label in the input TTree and
            the category label for the aggregated data (and therefore
            what is shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the spectra.
            This serves as a map between the category label for the
            spectra (value in the `_categories` dictionary) and the color
            to use for the histogram. The color can be any valid matplotlib
            color string or a cycle indicator (e.g. 'C0', 'C1', etc.).

        Returns
        -------
        None.
        """
        self._style = style
        self._variables = variables
        self._categories = categories
        self._colors = colors
        self._plotdata = None
        self._binedges = None
        self._onebincount = None

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra1D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra1D object.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        if is_ordinate:
            self._exposure_type = sample._exposure_type
            if self._exposure_type == 'pot':
                self._exposure = sample._exposure_pot
            else:
                self._exposure = sample._exposure_livetime

    def plot(self) -> None:
        """
        Plots the data for the SpineSpectra object. This is a base
        method and contains only common plotting elements. The actual
        plotting is done in the derived classes.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        """
        self._figure = plt.figure()
        self._ax = self._figure.add_subplot()

    def _mark_pot(self) -> None:
        """
        Add the POT information to the plot.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        """
        yrange = self._ax.get_ylim()
        usey = yrange[1] + 0.025*(yrange[1] - yrange[0])
        xrange = self._ax.get_xlim()
        usex = xrange[1] - 0.025*(xrange[1] - xrange[0])
        mag = int(np.floor(np.log10(self._exposure)))
        usepot = self._exposure/10**mag
        s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
        self._ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')

    def _mark_preliminary(self, label) -> None:
        """
        Add a preliminary label to the plot.

        Parameters
        ----------
        label : str
            The label to add to the plot to indicate that the plot is
            preliminary.

        Returns
        -------
        None.
        """
        yrange = self._ax.get_ylim()
        usey = yrange[1] + 0.025*(yrange[1] - yrange[0])
        xrange = self._ax.get_xlim()
        usex = xrange[0] + 0.025*(xrange[1] - xrange[0])
        self._ax.text(x=usex, y=usey, s=label, fontsize=14, color='#d67a11')

    def _mark_title(self, label) -> None:
        yrange = self._ax.get_ylim()
        usey = yrange[1] + 0.025*(yrange[1] - yrange[0])
        xrange = self._ax.get_xlim()
        usex = (xrange[1] + xrange[0])/2
        self._ax.text(x=usex, y=usey, s=label, fontsize=14, color='black')

class SpineSpectra1D(SpineSpectra):
    """
    A class designed to encapsulate a single variable's spectrum for an
    ensemble of samples. The class method add_sample() can be used to
    add a sample to the SpineSpectra

    Attributes
    ----------
    _style : str
        The name of the style sheet to use for the plot.
    _variable : Variable
        The Variable object for the spectrum.
    _categories : dict
        A dictionary of the categories for the spectrum. This serves as
        a map between the category label in the input TTree and the
        category label for the spectrum (and therefore what is shown
        in a single legend entry).
    _colors : dict
        A dictionary of the colors for the categories in the spectrum.
        This serves as a map between the category label for the
        spectrum (value in the `_categories` dictionary) and the color
        to use for the histogram. The color can be any valid matplotlib
        color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
    _plotdata : dict
        A dictionary of the data for the spectrum. This is a map between
        the category label for the spectrum and the histogram data for
        that category.
    """
    def __init__(self, style, variable, categories, colors, category_types) -> None:
        """
        Initializes the SpineSpectra1D.

        Parameters
        ----------
        style : str
            The name of the style sheet to use for the plot.
        variable : Variable
            The Variable object for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This
            serves as a map between the category label in the input
            TTree and the category label for the spectrum (and
            therefore what is shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the color to use for the histogram. The color can be
            any valid matplotlib color string or a cycle indicator
            (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the type of plot to use for the histogram. The type
            should be either 'histogram' or 'scatter' to correspond to
            a stacked histogram or scatter plot, respectively.

        Returns
        -------
        None.
        """
        super().__init__(style, [variable,], categories, colors)
        self._variable = self._variables[0]
        self._category_types = category_types

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra1D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra1D object.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)

        if self._plotdata is None:
            self._plotdata = {}
            self._binedges = {}
            self._onebincount = {}
        data, weights = sample.get_data([self._variable._key,])
        for category, values in data.items():
            values = values[0]
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros(self._variable._nbins)
                self._onebincount[self._categories[category]] = 0
            h = np.histogram(values, bins=self._variable._nbins, range=self._variable._range, weights=weights[category])
            self._onebincount[self._categories[category]] += np.sum(weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

    def plot(self, style, path, name) -> None:
        """
        Plots the data for the SpineSpectra1D object.

        Parameters
        ----------
        style : Style
            The Style object to use for the plot.
        path : str
            The path to the output directory.
        name : str
            The name of the output image file.

        Returns
        -------
        None.
        """
        super().plot()
        self._ax.set_xlabel(self._variable._xlabel)
        self._ax.set_ylabel('Candidates')
        self._ax.set_xlim(*self._variable._range)

        if self._plotdata is not None:
            labels, data = zip(*self._plotdata.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges[l][:-1] + np.diff(self._binedges[l]) / 2 for l in labels]

            histogram_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'histogram']
            scatter_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'scatter']

            denominator = np.sum([self._onebincount[labels[i]] for i in histogram_mask])
            counts = [x for x in self._onebincount.values()]
            if style.get_show_component_number() and style.get_show_component_percentage():
                hlabel = lambda x : f'{np.sum(x):.1f}, {np.sum(x)/denominator:.2%}'
                slabel = lambda x : f'{np.sum(x):.1f}'
                labels = [f'{label} ({hlabel(d) if li in histogram_mask else slabel(d)})' for li, (label, d) in enumerate(zip(labels, counts))]
            elif style.get_show_component_number():
                labels = [f'{label} ({np.sum(d):.1f})' for label, d in zip(labels, counts)]
            elif style.get_show_component_percentage():
                labels = [f'{label} ({np.sum(d)/denominator:.2%})' if li in histogram_mask else label for li, (label, d) in enumerate(zip(labels, counts))]

            if style.get_invert_stack_order():
                reduce = lambda x : [x[i] for i in histogram_mask[::-1]]
            else:
                reduce = lambda x : [x[i] for i in histogram_mask]
            self._ax.hist(reduce(bincenters), weights=reduce(data), bins=self._variable._nbins, range=self._variable._range, histtype='barstacked', label=reduce(labels), color=reduce(colors), stacked=True)

            reduce = lambda x : [x[i] for i in scatter_mask]
            for i, label in enumerate(reduce(labels)):
                self._ax.errorbar(bincenters[scatter_mask[i]], data[scatter_mask[i]], yerr=np.sqrt(data[scatter_mask[i]]), fmt='o', label=label, color=colors[scatter_mask[i]])
        
        if style.get_invert_stack_order():
            h, l = self._ax.get_legend_handles_labels()
            self._ax.legend(h[::-1], l[::-1])
        else:
            self._ax.legend()
        if style.get_mark_pot():
            self._mark_pot()
        if style.get_mark_preliminary() is not None:
            self._mark_preliminary(style.get_mark_preliminary())
        if style.get_mark_title() is not None:
            self._mark_title(style.get_mark_title())
        if style.get_title() is not None:
            self._ax.set_title(style.get_title())
        self._figure.savefig(f'{path}/{name}.png')

class SpineSpectra2D(SpineSpectra):
    """
    A class designed to encapsulate a pair of variables' spectrum for
    an ensemble of samples. The class method add_sample() can be used
    to add a sample to the SpineSpectra2D.

    Attributes
    ----------
    _style : str
        The name of the style sheet to use for the plot.
    _variables : list
        The list of Variable objects for the spectrum.
    _categories : dict
        A dictionary of the categories for the spectrum. This serves as
        a map between the category label in the input TTree and the
        category label for the spectrum (and therefore what is shown
        in a single legend entry).
    _colors : dict
        A dictionary of the colors for the categories in the spectrum.
        This serves as a map between the category label for the
        spectrum (value in the `_categories` dictionary) and the color
        to use for the histogram. The color can be any valid matplotlib
        color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
    _plotdata : dict
        A dictionary of the data for the spectrum. This is a map between
        the category label for the spectrum and the histogram data for
        that category.
    """
    def __init__(self, style, variables, categories, colors, category_types) -> None:
        """
        Initializes the SpineSpectra2D object.

        Parameters
        ----------
        style : str
            The name of the style sheet to use for the plot.
        variables : list
            The list of Variable objects for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This serves
            as a map between the category label in the input TTree and
            the category label for the spectrum (and therefore what is
            shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the spectrum.
            This serves as a map between the category label for the
            spectrum (value in the `_categories` dictionary) and the color
            to use for the histogram. The color can be any valid matplotlib
            color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the spectrum.
            This serves as a map between the category label for the spectrum
            (value in the `_categories` dictionary) and the type of plot to
            use for the histogram. The type should be either 'histogram' or
            'scatter' to correspond to a stacked histogram or scatter plot,
            respectively.

        Returns
        -------
        None.
        """
        super().__init__(style, variables, categories, colors)
        self._category_types = category_types
        self._plotdata_diagonal = None
        self._binedges_diagonal = None

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra2D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra2D object.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)

        if self._plotdata is None:
            self._plotdata = {}
            self._binedges = {}
        if self._plotdata_diagonal is None:
            self._plotdata_diagonal = {}
            self._binedges_diagonal = {}

        data, weights = sample.get_data([self._variables[0]._key, self._variables[1]._key])        
        for category, values in data.items():
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros((self._variables[0]._nbins, self._variables[1]._nbins))
            h = np.histogram2d(values[0], values[1], bins=(self._variables[0]._nbins, self._variables[1]._nbins), range=(self._variables[0]._range, self._variables[1]._range), weights=weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

            if self._categories[category] not in self._plotdata_diagonal:
                self._plotdata_diagonal[self._categories[category]] = np.zeros(self._variables[0]._nbins)
            diag = np.divide(values[1] - values[0], values[0])#, where=values[0] != 0)
            h = np.histogram(diag, bins=self._variables[0]._nbins, range=(-4,4), weights=weights[category])
            self._plotdata_diagonal[self._categories[category]] += h[0]
            self._binedges_diagonal[self._categories[category]] = h[1]
        

    def plot(self, style, path, name) -> None:
        """
        Plots the data for the SpineSpectra2D object.

        Parameters
        ----------
        style : Style
            The Style object to use for the plot.
        path : str
            The path to the output directory
        name : str
            The name of the output image file.
        
        Returns
        -------
        None.
        """
        super().plot()
        self._ax.set_xlabel(self._variables[0]._xlabel)
        self._ax.set_ylabel(self._variables[1]._xlabel)

        if self._plotdata is not None:
            values = np.sum([v for v in self._plotdata.values()], axis=0)
            binedges = self._binedges[list(self._plotdata.keys())[0]]
            self._ax.imshow(values.T, extent=(binedges[0], binedges[-1], binedges[0], binedges[-1]), aspect='auto', origin='lower', cmap='cividis')
            self._figure.savefig(f'{path}/{name}.png')
            
    def plot_diagonal_reduction(self, style, path, name) -> None:
        """
        Plots the data for the SpineSpectra2D object with a diagonal
        reduction.

        Parameters
        ----------
        style : Style
            The Style object to use for the plot.
        path : str
            The path to the output directory.
        name : str
            The name of the output image file.
        
        Returns
        -------
        None.
        """
        super().plot()
        self._ax.set_xlabel(f'{self._variables[1]._xlabel} - {self._variables[0]._xlabel} / {self._variables[0]._xlabel}')
        self._ax.set_ylabel('Candidates')

        if self._plotdata_diagonal is not None:
            labels, data = zip(*self._plotdata_diagonal.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges_diagonal[l][:-1] + np.diff(self._binedges_diagonal[l]) / 2 for l in labels]

            self._ax.hist(bincenters, weights=data, bins=self._variables[0]._nbins, range=(-4,4), histtype='barstacked', label=labels, color=colors, stacked=True)
            self._figure.savefig(f'{path}/{name}_diag.png')
