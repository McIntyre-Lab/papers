""" Various plotting function """

import numpy as np
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt


class ColorMarker:
    """ A class to handle colors for plotting with with DataFrames. """
    def __init__(self):
        """ Initialize the ColorMarker class.

        Attributes:
            :param list self._colors: A list of the standard python colors.

            :param list self._markers: A list of the standard python markers.

            :type self.colorMarker: list of tuple
            :param self.colorMarker: A list of tuples with (color, marker)
                pairs. There are 114 total combinations.

        """
        # A list of colors
        self._colors = ['k', 'b', 'g', 'c', 'm', 'y']

        # A list of markers
        self._markers = ['o', 's', '^', 'D', 'd', 'h', 'x', '*', '+', 'v', '<', '>', '1', '2', '3', '4', '8', 'p', 'H']

        # Create a list of color|marker pairs
        self.colorMarker = self._get_colors()

    def _get_colors(self):
        """ Get a set of color/marker combinations.

        :rtype: list of tuple
        :returns: A list of tuples containing color|marker pairs. There are a total
            of 114 combinations. Red and white are not used in this color scheme.
            Red is reserved for coloring points beyond a threshold, and white does not
            show up on white backgrounds.

        """
        comb = list()
        for marker in self._markers:
            for color in self._colors:
                comb.append((color, marker))

        return comb

    def groupColorMarker(self, design, groupID):
        """ Adds colors and markers to each group.

        Arguments:
            :type design: pandas.DataFrame
            :param design: A DataFrame with design information.

            :param str groupID: The column name containing group information in
                the design file.

        Returns:
            :rtype: tuple of dict
            :returns: A tuple of dictionaries with sampleID as the key and
                color or marker as value.

        """

        # Get group information from design
        grp = design.groupby(groupID)
        groupIDs = sorted(grp.groups.keys())
        groupNum = len(groupIDs)

        # Pull out the number of color|marker pairs needed for the number of groups
        if groupNum <= 114:
            colorMarker = self.colorMarker[:groupNum]
        else:
            print('Warning - you have more groups than color|marker combinations, will be repeating colors and markers!')
            multiplier = np.ceil(groupNum / 114.0)
            cm = self.colorMarker * multiplier
            colorMarker = cm[:groupNum]

        # Make a DataFrame of color|marker with groups
        dfColorMarker = pd.DataFrame(colorMarker, columns=['color', 'marker'], index=groupIDs)

        # Merge colors to design file
        merged = design.merge(dfColorMarker, how='left', left_on='group', right_index=True)

        # Generate color and marker dictionary
        myDict = merged.drop('group', axis=1).to_dict()
        colorsDict = myDict['color']
        markerDict = myDict['marker']

        return colorsDict, markerDict


def blandAltman(s1, s2, ax=None):
    """ Generate a Bland-Altman plot.

    Arguments:
        :type s1: numpy.array
        :param s1: An array of sample1 data.

        :type s2: numpy.array
        :param s2: An array of sample2 data.

        :param ax:
        :type ax: matplotlib.axes.AxesSubplot

    Returns:
        :returns: If avaiable returns a matplotlib.figure.Figure else adds plot
            to current axis.
        :rtype: matplotlib.figure.Figure

    """

    # Make sure s1 and s2 are numpy arrays
    s1 = np.asarray(s1)
    s2 = np.asarray(s2)

    # Calculate mean and difference
    mean = (s1 + s2) / 2
    diff = s1 - s2

    # make plot if not axis was provided
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    ax.scatter(mean, diff)
    ax.axhline(0, color='r', ls='--', lw=2)
    ax.set_xlabel('Mean')
    ax.set_ylabel('Difference')

    # Try to return the figure if it exists
    try:
        return fig
    except:
        pass


if __name__ == '__main__':
    pass
