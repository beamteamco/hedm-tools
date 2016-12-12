from .crystals import UnitCell, Cubic, Hexagonal
import numpy as np
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
from mpl_toolkits import mplot3d


class HEDMAxes(mplot3d.Axes3D):

    name = 'hedm'

    def __init__(self, *args, **kwds):
        super(HEDMAxes, self).__init__(*args, **kwds)

    def grains3D(self, x, y, z, orient, **kwds):
        """
        Creates a 3D scatter plot of oriented crystal structures.

        Input
        -----
        :x, N-length array-like: x-coordinates
        :y, N-length array-like: y-coordinates
        :z, N-length array-like: z-coordinates
        :orient, Nx3x3 array-like: orientation matrices
        :s, N-length array-like | num: size of the grains. (optional)
        :c, N-length array-like | str: color of the grains.
            Default: 'r'.
        :crystal, str | UnitCell: crystal structure representation.
            'cubic'|'hexagonal'. Default: 'cubic'.
        :cmap, str | Colormap: matplotlib colormap (name or Colormap object).
            Default: 'Spectral'.
        """
        # verify required arguments
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        orient = np.asarray(orient)

        assert x.shape == y.shape, \
            'x shape {xshape:} != y shape {yshape:}'.format(x.shape, y.shape)
        assert y.shape == z.shape, \
            'z shape {zshape:} does not match other dimensions'.format(z.shape)
        assert orient.shape[0] == x.shape[0], \
            '{xsize:}-(3x3) orientation matrices required to match center of ' \
            'mass positions'.format(x.shape[0])

        # get stats
        count = x.shape[0]
        xlim_ = (x.min(), x.max())
        ylim_ = (y.min(), y.max())
        zlim_ = (z.min(), z.max())

        # check and verify optional kwds
        ## size
        s = kwds.get('s', None)
        ## set default value if none specified
        if s is None:
            width = np.max((xlim_[1] - xlim_[0],
                            ylim_[1] - ylim_[0],
                            zlim_[1] - zlim_[0]))
            s = 0.1*width
        s = np.asarray(s)
        s = np.repeat(s, count) if s.shape == () else s

        if len(s) != count:
            msg = 'Invalid size (s), must be a single size or an ' \
            'array-like object the same length as the number of coordinates.'
            raise ValueError(msg)

        ## color
        c = kwds.get('c', 'r') # defaults to red
        c = np.asarray(c)
        c = np.repeat(c, count) if c.shape == () else c

        if len(c) != count:
            msg = 'Invalid color (c), must be a single color or an ' \
            'array-like object the same length as the number of coordinates.'
            raise ValueError(msg)

        ## crystal
        shape = kwds.get('crystal', 'cubic')
        if isinstance(shape, str):
            shape = shape.lower()
        shape = {
            'cubic' : Cubic,
            'hexagonal' : Hexagonal
        }.get(shape, shape)

        if not isinstance(shape(), UnitCell):
            msg = '{xtal:} is not a recognized crystal structure'.format(
                xtal=kwds['crystal'])

        ## colormap
        colormap = kwds.get('cmap', 'Spectral')
        if isinstance(colormap, str):
            colormap = cm.get_cmap(colormap)

        if not isinstance(colormap, Colormap):
            msg = '{cmap:} is not a recognized colormap'.format(cmap=cmap)

        # create a mappable object
        mappable = None
        if np.isreal(c.dtype.type()):
            cmin = kwds.get('vmin', c.min())
            cmax = kwds.get('vmax', c.max())
            mappable = ScalarMappable(
                norm=Normalize(vmin=cmin, vmax=cmax),
                cmap=colormap)

        # create a Path3DCollection object to provide a handle for
        # colorbars, etc.
        if mappable is None:
            collection = self.scatter3D(x, y, z, s=0, c=c)
        else:
            collection = self.scatter3D(x, y, z, s=0, c=c,
                cmap=mappable.cmap,
                norm=mappable.norm)

        # add polygons
        for i in xrange(count):
            poly = shape()
            poly.rotate(orient[i])
            poly.scale(s[i])
            poly.faces = poly + (x[i], y[i], z[i])
            if mappable is None:
                pc = poly.polycollection(facecolors=(c[i],))
            else:
                pc = poly.polycollection(facecolors=(mappable.to_rgba(c[i]),))
            self.add_collection3d(pc)

        # create cubic bouding box to simulate equal aspect ratio
        # from http://stackoverflow.com/questions/13685386/\
        #   matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
        max_range = np.array([xlim_[1] - xlim_[0],
                              ylim_[1] - ylim_[0],
                              zlim_[1] - zlim_[0]]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() \
            + np.mean(xlim_)
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() \
            + np.mean(ylim_)
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() \
            + np.mean(zlim_)
        for xb, yb, zb in zip(Xb, Yb, Zb):
            _ = self.scatter([xb], [yb], [zb], s=0)

        # return the Path3DCollection from the scatter dummy
        return collection
