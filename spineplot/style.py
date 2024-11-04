import matplotlib.pyplot as plt

class Style:
    """
    A class designed to encapsulate the configuration of a plotting
    style for a given plot.

    Attributes
    ----------
    _style : str
        The name of the style sheet to use for the plot.
    _show_component_number : bool
        A flag toggling the display of the total number of selected
        interactions per component in the legend.
    _show_component_percentage : bool
        A flag toggling the display of the percentage of selected
        interactions per component in the legend.
    _invert_stack_order : bool
        A flag toggling the inversion of the stack order for the
        components in the histogram.
    """
    def __init__(self, style_sheet, show_component_number, show_component_percentage, invert_stack_order) -> None:
        """
        Initializes the Style object with the given kwargs.

        Parameters
        ----------
        style_sheet : str
            The name of the style sheet to use for the plot.
        show_component_number : bool
            A flag toggling the display of the total number of selected
            interactions per component in the legend.
        show_component_percentage : bool
            A flag toggling the display of the percentage of selected
            interactions per component in the legend.
        invert_stack_order : bool
            A flag toggling the inversion of the stack order for the
            components in the histogram.

        Returns
        -------
        None
        """
        self._style = style_sheet
        self._show_component_number = show_component_number
        self._show_component_percentage = show_component_percentage
        self._invert_stack_order = invert_stack_order

    def __enter__(self):
        """
        Sets the style for the plot.

        Returns
        -------
        self : Style
            The created Style object.
        """
        plt.style.use(self._style)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Resets the style for the plot.

        Parameters
        ----------
        exc_type : type
            The exception type.
        exc_val : Exception
            The exception value.
        exc_tb : traceback
            The traceback.

        Returns
        -------
        None
        """
        plt.style.use('default')

    def get_show_component_number(self) -> bool:
        """
        Returns the value of the show_component_number attribute.

        Returns
        -------
        bool
            The value of the show_component_number attribute.
        """
        return self._show_component_number
    
    def get_show_component_percentage(self) -> bool:
        """
        Returns the value of the show_component_percentage attribute.

        Returns
        -------
        bool
            The value of the show_component_percentage attribute.
        """
        return self._show_component_percentage
    
    def get_invert_stack_order(self) -> bool:
        """
        Returns the value of the invert_stack_order attribute.

        Returns
        -------
        bool
            The value of the invert_stack_order attribute.
        """
        return self._invert_stack_order