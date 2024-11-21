import toml
import uproot
from sample import Sample
from spinespectra import SpineSpectra1D
from style import Style
from variable import Variable
from matplotlib import pyplot as plt

class ConfigException(Exception):
    pass

class Analysis:
    """
    The Analysis class is used to containerize the plotting of the
    ensemble of samples for a variety of variables and plotting styles.

    Attributes
    ----------
    _toml_path : str
        The path to the TOML configuration file for the analysis.
    _config : dict
        The configuration dictionary loaded from the TOML file.
    _ordinate_sample : str
        The name of the sample to use as the ordinate for the analysis.
    _category_tree : str
        The label of the branch corresponding to the category label in
        the input TTree.
    """
    def __init__(self, toml_path, rf_path) -> None:
        """
        Initializes the Analysis object with the given kwargs.

        Parameters
        ----------
        toml_path : str
            The path to the TOML configuration file for the analysis.
        rf_path : str
            The path to the ROOT file containing the data.

        Returns
        -------
        None
        """
        self._toml_path = toml_path
        self._config = toml.load(self._toml_path)
        for table in self._config['this_includes']:
            Analysis.handle_include(self._config, table)
        rf = uproot.open(rf_path)

        # Load the categories table
        self._categories = dict()
        for ci, cat in enumerate(self._config['analysis']['category_assignment']):
            self._categories.update({c : self._config['analysis']['category_labels'][ci] for c in cat})
        self._colors = {c : self._config['analysis']['category_colors'][ci] for ci, c in enumerate(self._config['analysis']['category_labels'])}
        self._category_types = {c : self._config['analysis']['category_types'][ci] for ci, c in enumerate(self._config['analysis']['category_labels'])}

        # Initialize the samples
        if 'samples' not in self._config.keys():
            raise ConfigException(f"No samples defined in the TOML file. Please check for a valid sample configuration block in the TOML file ('{toml_path}').")
        self._samples = {name: Sample(name, rf, self._config['analysis']['category_branch'], **self._config['samples'][name]) for name in self._config['samples']}

        # Load the plot styles table
        if 'styles' not in self._config.keys():
            raise ConfigException(f"No plot styles defined in the TOML file. Please check for a valid plot style configuration block (table='styles') in the TOML file ('{toml_path}').")
        self._styles = {name: Style(**self._config['styles'][name]) for name in self._config['styles'].keys()}

        # Load the variable table
        if 'variables' not in self._config.keys():
            raise ConfigException(f"No variables defined in the TOML file. Please check for a valid variable configuration block (table='variables') in the TOML file ('{toml_path}').")
        self._variables = {name: Variable(name, **self._config['variables'][name]) for name in self._config['variables']}

        # Load the plots table and initialize the SpineSpectra objects
        if 'spectra1D' not in self._config.keys():
            raise ConfigException(f"No plots defined in the TOML file. Please check for a valid plot configuration block (table='plots') in the TOML file ('{toml_path}').")
        self._spectra = {name: SpineSpectra1D(v['style'], self._variables[v['variable']], self._categories, self._colors, self._category_types) for name, v in self._config['spectra1D'].items()}
    
    def override_exposure(self, sample_name, exposure, exposure_type='pot') -> None:
        """
        Overrides the exposure for the given sample. This is useful for
        setting the exposure for samples for which the exposure is not
        valid. The exposure type can be either 'pot' or 'livetime'. It
        is not recommended to use this method unless the exposure is
        known to be incorrect.
        
        Parameters
        ----------
        sample_name : str
            The name of the sample to override.
        exposure : float
            The exposure to set for the sample.
        exposure_type : str
            The type of exposure to set. This can be either 'pot' or
            'livetime'. The default is 'pot'.
        
        Returns
        -------
        None.
        """
        if sample_name not in self._samples.keys():
            raise ConfigException(f"Sample '{sample_name}' not found in sample list when attempting to override exposure. Please check the sample configuration block in the TOML file ('{self._toml_path}').")
        self._samples[sample_name].override_exposure(exposure, exposure_type)

    def run(self) -> None:
        """
        Runs the analysis on the samples.

        Returns
        -------
        None.
        """
        if self._config['analysis']['ordinate_sample'] not in self._samples.keys():
            raise ConfigException(f"Ordinate sample '{self._config['analysis']['ordinate_sample']}' not found in sample list. Please check the sample configuration block (table='samples') in the TOML file ('{self._toml_path}').")
        ordinate = self._samples[self._config['analysis']['ordinate_sample']]
        for s in self._samples.values():
            s.set_weight(target=ordinate)

        for name, s in self._spectra.items():
            for sample in self._samples.values():
                s.add_sample(sample)

            with self._styles[s._style] as style:
                s.plot(style, name)

    @classmethod
    def handle_include(self, config, table):
        """
        Handles the inclusion of other configuration files in the main
        configuration file. The include directive may also contain some
        other optional fields which correspond to specific actions:

        - choose: a dictionary of key-value pairs where the key is the
        name of a table and the value is a list of corresponding sub-
        tables to include.

        Parameters
        ----------
        config : dict
            The configuration dictionary to update.
        table : dict
            The block representing the include directive.
        
        Returns
        -------
        None.
        """
        with open(table['file'], 'r') as f:
            if 'choose' in table.keys():
                for key, value in table['choose'].items():
                    c = toml.load(f)
                    config[key] = {v: c[key][v] for v in value}
            else:
                config.update(toml.load(f))