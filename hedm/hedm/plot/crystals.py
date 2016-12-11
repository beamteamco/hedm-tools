#!/usr/bin/env python

from hedm import quaternion
import numpy as np
from matplotlib.colors import colorConverter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class UnitCell(object):
    def __init__(self, top, bottom, *args, **kwargs):
        """
        Constructs a unit cell object. This base class presents
        a consistent API for subsequently derived types of unit
        cells -- e.g. cubic, tetragonal, etc. -- which follow.

        The vertices from the top face and bottom face are
        sufficient to construct the unit cell. The vertices
        from the top and bottom face should correlate, i.e. an
        edge is shared between the first vertex in the top face
        and the first vertex in the bottom face; the second, and
        the second; etc. The top and the bottom faces should
        "close the loop" on the face, i.e. the first and the
        last vertices should be the same.

        Description

            faces: returns an array of the faces that compose
                the unit cell

            rotate(q): rotates the unit cell by applying the
                specified rotation function. `q` should accept
                a 3-vector, v, and return the rotated 3-vector
                v'

            scale(s): scale the unit cell vertices. If `s` is
                a scalar, then a uniform scaling is applied,
                otherwise, `s = (sx, sy, sz)`.

            aspolycollection(facecolors, alpha): returns the
                UnitCell as a matplotlib.PolyCollection object.
                Facecolors and alpha (both optional) set the
                color of each face and the transparency of the
                object, respectively.

            __add__(vec): operator overload `cell + vec` returns
                a copy of all faces translated by `vec`

            __sub__(vec): operator overload `cell - vec` returns
                a copy of all faces translated by `vec`

        """
        super(UnitCell, self).__init__(*args, **kwargs)
        _faces = [top]
        for ((t1, t2), (b1, b2)) in zip(zip(top[:-1], top[1:]), zip(bottom[:-1], bottom[1:])):
            _faces.append([t1, b1, b2, t2, t1])
        _faces.append(bottom)
        self.faces = _faces

    @property
    def faces(self):
        return self._faces
    @faces.setter
    def faces(self, faces_):
        self._faces = [np.asarray(face) for face in faces_]

    def _rotate(self, qfunc):
        """
        Rotates the faces in UnitCell according to `qfunc(x)`, where `x`
        is a vertex.
        """
        self.faces = [np.apply_along_axis(qfunc, -1, face) for face in self.faces]

    def rotate(self, *rotation):
        """
        Rotate the UnitCell. The rotation is specified as one of:

            axis, angle: axis/angle pair; angle in degrees.
            quaternion: single Quaternion object
            rotation matrix: 3x3 rotation matrix

        Examples

            `ucell.rotate(axis, angle)`
            `ucell.rotate(quaternion)`
            `ucell.rotate(rotation_matrix)`
        """
        # infer argument types based on number/type
        # convert to rotation matrix
        try:
            # two arguments --> axis/angle pair
            axis, angle = rotation
            rotation_matrix = quaternion.to_matrix(quaternion.from_axis_angle(axis, angle))
        except ValueError:
            # invalid arguments passed to rotate
            if len(rotation) != 1:
                raise ValueError('An invalid number of arguments were passed to "rotate".')
            # single argument --> quaternion or rotation matrix
            obj, = rotation # comma because args is a tuple
            if isinstance(obj, quaternion.Quaternion):
                rotation_matrix = quaternion.to_matrix(obj)
            else:
                rotation_matrix = np.asarray(obj)
        # confirm the shape of the rotation matrix
        if rotation_matrix.shape != (3,3):
            raise ValueError('Rotation matrix must be given as a 3x3 array.')
        # perform rotation
        qfunc = lambda x : np.dot(x, rotation_matrix.T)
        return self._rotate(qfunc)

    def scale(self, *s):
        """
        Scale the vertices of the unit cell. If `s` is a single value, scale isotropically;
        if not,
        """
        if len(s) == 1:
            s = s[0]*np.ones(3)
        s = np.asarray(s)
        scale_ = lambda x : np.multiply(s, x)
        self.faces = [np.apply_along_axis(scale_, -1, face) for face in self.faces]

    def __add__(self, vec_):
        op = lambda x : np.add(x, vec_)
        return [np.apply_along_axis(op, -1, face) for face in self.faces]

    def __sub__(self, vec_):
        op = lambda x : np.subtract(x, vec_)
        return [np.apply_along_axis(op, -1, face) for face in self.faces]

    def polycollection(self, facecolors='r', alpha=1.0, linewidths=1):
        """
        Constructs and returns a polygonal representation of the UnitCell.

        Example

        ```
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            poly = ucell.polycollection()
            # uncomment to plot vertices as points
            # ax.scatter(*np.transpose(np.reshape(ucell.faces, (-1, 3))), c='k')
            ax.add_collection3d(poly)
        """
        nfaces = len(self.faces)
        if isinstance(facecolors, str):
            facecolors = (facecolors,)
        if len(facecolors) == 1:
            facecolors = len(self.faces)*facecolors
        if len(facecolors) != nfaces:
            raise ValueError('Either {} or a uniform color must be given '                              'for this UnitCell.'.format(nfaces))
        if isinstance(facecolors[0], str):
            facecolors = [colorConverter.to_rgba(fc, alpha=alpha) for fc in facecolors]
        return Poly3DCollection(
            self.faces,
            facecolors=facecolors,
            linewidths=linewidths,
            alpha=alpha)


class Cubic(UnitCell):
    def __init__(self, *args, **kwargs):
        """
        Construct a cubic unit cell. See `UnitCell` documentation for
        accepted operations.
        """
        top = 0.5*np.array([
            [1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1], [1, 1, 1]])
        bottom = 0.5*np.array([
            [1, 1, -1], [-1, 1, -1], [-1, -1, -1], [1, -1, -1], [1, 1, -1]])
        super(Cubic, self).__init__(top, bottom, *args, **kwargs)


class Hexagonal(UnitCell):
    def __init__(self, c_over_a=1, *args, **kwds):
        """
        Construct a cubic unit cell with an optional c over a ratio
        as the first argument.
        """
        x = np.cos(np.pi/3)
        y = np.sin(np.pi/3)
        h = float(c_over_a)/2.
        top = np.array([
            [1, 0, h], [x, y, h], [-x, y, h], [-1, 0, h],
                [-x, -y, h], [x, -y, h], [1, 0, h]])
        bottom = np.array([
            [1, 0, -h], [x, y, -h], [-x, y, -h], [-1, 0, -h],
                [-x, -y, -h], [x, -y, -h], [1, 0, -h]])
        super(Hexagonal, self).__init__(top, bottom, *args, **kwds)


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Colormap, Normalize
    from matplotlib.colorbar import ColorbarBase
    from mpl_toolkits import mplot3d

    # # ## Test cubes
    #
    # # In[26]:
    #
    # cubes = []
    # for _ in range(2):
    #     c = Cubic()
    #     c.rotate(np.random.random(3), 180*np.random.random())
    #     c.scale(np.random.random()+2)
    #     c.faces = c + (2.*np.random.random(3) - 1)
    #     cubes.append(c)
    #
    #
    # # In[27]:
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # colors = ['r', 'g', 'b']
    # for c in cubes:
    #     color = colors.pop()
    #     polys = c.polycollection(facecolors=color)
    #     ax.add_collection3d(polys)
    # #polys = cubic.polycollection()
    # #ax.scatter(*np.transpose(np.reshape(cubic.faces, (-1, 3))), c='k')
    # #ax.add_collection3d(polys)
    # _ = ax.set_xlim(-2, 2)
    # _ = ax.set_ylim(-2, 2)
    # _ = ax.set_zlim(-2, 2)
    # plt.show()
    #
    #
    # # ## Test hexagonal
    #
    # # In[29]:
    #
    # hexes = []
    # for _ in range(2):
    #     c = Hexagonal()
    #     c.rotate(np.random.random(3), 180*np.random.random())
    #     c.scale(np.random.random()+2)
    #     c.faces = c + (2.*np.random.random(3) - 1)
    #     hexes.append(c)
    #
    #
    # # In[30]:
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # colors = ['r', 'g', 'b']
    # for c in hexes:
    #     color = colors.pop()
    #     polys = c.polycollection(facecolors=color)
    #     ax.add_collection3d(polys)
    # #polys = cubic.polycollection()
    # #ax.scatter(*np.transpose(np.reshape(cubic.faces, (-1, 3))), c='k')
    # #ax.add_collection3d(polys)
    # _ = ax.set_xlim(-2, 2)
    # _ = ax.set_ylim(-2, 2)
    # _ = ax.set_zlim(-2, 2)
    # plt.show()


    # ## Impatience...

    # In[9]:
    # Use a Grains.csv from the Ti7Al data collected by
    # Garrison Hommer in November, 2016 at APS
    grains = np.genfromtxt('Grains.csv',
        skip_header=8, names=True, dtype=None)


    # In[10]:

    print ' '.join(grains.dtype.names)
    print len(grains)


    # In[31]:

    hexes = []
    xyz = grains[np.array(['X', 'Y', 'Z'])]\
        .view(np.float64).reshape((-1, 3))
    orient = grains[np.array(['O11', 'O12', 'O13',
                              'O21', 'O22', 'O23',
                              'O31', 'O32', 'O33'])]\
        .view(np.float64).reshape((-1, 3, 3))
    size = grains['GrainRadius'].view(np.float64)
    size = 80*(size - size.min())/(size.max() - size.min()) + 40
    for i in range(len(grains)):
        c = Hexagonal()
        c.rotate(orient[i])
        c.scale(size[i])
        c.faces = c + xyz[i]
        hexes.append(c)


    # In[43]:

    # set scalar value
    #scalar = grains['Confidence'].view(np.float64)
    exx, eyy, ezz, gxy, gyz, gzx = \
        grains[np.array(['eFab11', 'eFab22', 'eFab33',
                         'eFab12', 'eFab23', 'eFab31'])]\
            .view(np.float64).reshape((-1, 6)).T
    scalar = 2./3.*np.sqrt(3*(exx**2 + eyy**2 + ezz**2)/2. \
        + 3*(gxy**2 + gyz**2 + gzx**2)/4.)
    # select colormap
    cmap = cm.get_cmap('Spectral')
    # create normalizer for colorbar
    norm = Normalize(vmin=scalar.min(), vmax=scalar.max())
    # normalize scalar
    scalar = (scalar - scalar.min())/(scalar.max() - scalar.min())
    # construct figure
    fig = plt.figure(figsize=(12, 10))
    # plot axis
    ax = fig.add_axes([0.05, 0.05, 0.8, 0.9], projection='3d')
    # colorbar axis
    cbax = fig.add_axes([0.85, 0.2, 0.05, 0.7])
    #ax = fig.add_subplot(111, projection='3d')

    # # orthogonal projection
    # # DO NOT USE: While an interesting thought, the Poly3DCollection
    # # does not render properly.
    # def orthogonal_proj(zfront, zback):
    #     a = (zfront+zback)/(zfront-zback)
    #     b = -2*(zfront*zback)/(zfront-zback)
    #     return numpy.array([[1,0,0,0],
    #                         [0,1,0,0],
    #                         [0,0,a,b],
    #                         [0,0,0,zback]])
    # proj3d.persp_transformation = orthogonal_proj

    for i,c in enumerate(hexes):
        polys = c.polycollection(facecolors=(cmap(scalar[i]),))
        ax.add_collection3d(polys)

    ax.set_xlabel('beam direction ($\mathrm{\mu m}$)')
    ax.set_ylabel('inboard-outboard ($\mathrm{\mu m}$)')
    ax.set_zlabel('up ($\mathrm{\mu m}$)')

    cb = ColorbarBase(cbax, cmap=cmap, norm=norm, orientation='vertical')
    cb.set_label('Equivalent strain')

    # Create cubic bounding box to simulate equal aspect ratio
    # from http://stackoverflow.com/questions/13685386/\
    #   matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    X, Y, Z = xyz.T
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    # # Another option: manually set x,y,z limits
    # xmid, ymid, zmid = xyz.mean(axis=0)
    # lower = xyz.min(axis=0)
    # upper = xyz.max(axis=0)
    # halfspan = np.max(upper - lower)/2. # largest span
    # _ = ax.set_xlim(xmid - halfspan, xmid + halfspan)
    # _ = ax.set_ylim(ymid - halfspan, ymid + halfspan)
    # _ = ax.set_zlim(zmid - halfspan, zmid + halfspan)
    plt.show()
