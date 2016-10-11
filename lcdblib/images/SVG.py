#!/usr/bin/env python
""" A set of useful tools for handling SVGs in a Jupyter Notebook. """
from IPython.display import SVG, display
from xml.etree import ElementTree
from io import StringIO

class nb_svg(object):
    def __init__(self, svg):
        """ A class for handling SVG drawings and displaying them in the notebook.

        This object can display different layers of an SVG created by Inkscape.
        This allows you to create a set of drawings in a single SVG file and
        then easily display them in a notebook.

        Parameters
        ----------
        svg: str
            Path to an SVG file created with Inkscape.

        """
        # Name spaces
        self.SVG_NS = "http://www.w3.org/2000/svg"
        self.INKSCAPE_NS = 'http://www.inkscape.org/namespaces/inkscape'

        # Use IPython.display.SVG to import SVG file
        self.image = SVG(svg)

        # Use xml.etree.ElementTree to manipulate SVG
        ElementTree.register_namespace('', self.SVG_NS)
        self.ElementTree = ElementTree.ElementTree()
        self.tree = self.ElementTree.parse(StringIO(self.image.data))
        self.layers = self.tree.findall('.//{%s}g' % self.SVG_NS)
    
    def _filter(self, layer):
        """ Turn layers on and off. """
        keep = layer.split('|')
        for l in self.layers:
            if l.attrib.get('{%s}groupmode' % self.INKSCAPE_NS, '') == 'layer':
                if l.attrib['{%s}label' % self.INKSCAPE_NS] in keep:
                    if 'style' in l.attrib:
                        del l.attrib['style']
                else:
                    l.attrib['style'] = 'display:none'
                
        self.image.data = ElementTree.tostring(self.tree, encoding='utf-8', method='xml')
        
    def getLayers(self):
        """ Returns a list of layer Names. """
        names = []
        for l in self.layers:
            if l.attrib.get('{%s}groupmode' % self.INKSCAPE_NS, '') == 'layer':
                names.append(l.attrib['{%s}label' % self.INKSCAPE_NS])
        return names
        
    def display(self, layer=None):
        """ Use IPython display to show a SVG file.

        Specific layers can be show by specifying the name of the layer or
        layers separated by a '|'.

        Parameters
        ----------
        layer: str, None
           Layers to display can be specified as a string separating layers
           names with a '|'. If None then it will display whatever layers were
           visible in the file. 

        Examples
        --------
        >>> s = nv_svg('drawing.svg')
        >>> s.display('One|Two|Three')  # assuming layers are name One, Two, Three
        >>> s.display('One')

        """
        if layer is not None:
            self._filter(layer)
        return display(self.image)


if __name__ == '__main__':
    pass
