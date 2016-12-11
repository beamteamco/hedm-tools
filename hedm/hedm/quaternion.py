#!/bin/env python
"""
Summary
  Module for handling quaternions.

Contents

Quaternion (class)
  Derived from np.ndarray, the add and multiply features
have been overloaded to perform quaternion addition and
multiplication between quaternions and 3-vectors, e.g. for

$$
    u = (x, y, z)
    q_1 = (w, x, y, z)_1
    q_2 = (w, x, y, z)_2
$$

where $u$, $q_1$ and $q_2$ are a 3D vector and two quaternions,
respectively, the following are acceptable operations:

`uprime = q1*u*inverse(q1)` rotates `u` into `uprime` based on
the rotation encoded in `q1`.

`uprime = q2*q1*u*inverse(q1)*inverse(q2)` performs the rotation
first of `q1` then of `q2`. This is exactly equivalent to:

```
    q3 = q2*q1
    uprime = q3*u*inverse(q3)
```
"""

import numpy as np


class Quaternion(np.ndarray):
    def __new__(cls, input_vec=np.ones(4)):
        """
        Create a new Quaternion object. Except when slicing, which
        is discussed in a moment, this guarantees a 4-element array
        `(w, x, y, z)`, where `w` is the scalar part of the Quaternion
        and `(x, y, z)` is the vector part.

        When slicing, a portion of an existing quaternion may be
        returned. This is useful for normalizing the vector part,
        for example. To create a new Quaternion object, use the
        helper function `asquaternion(vec_)` where `vec_` can be a
        3- or 4- element array-like object (including a slice of a
        quaternion, list, numpy array, tuple, etc.).
        """
        obj = np.asarray(input_vec).view(cls)
        if obj.shape[-1] == 3:
            obj = np.concatenate(((0,), obj))
        if obj.shape[-1] != 4:
            msg = 'Quaternions can only be created from 3- or 4-vectors.'
            raise ValueError(msg)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __add__(self, vec_):
        """
        Quaternion addition between two quaternions or a quaternion and a
        3-vector.
        """
        lhs = self.view(np.ndarray)
        rhs = np.asarray(vec_)
        return Quaternion(lhs + rhs)

    def __radd__(self, vec_):
        return self + vec_

    def __mul__(self, vec_):
        """
        Quaternion multiplication between two quaternions or a quaternion and a
        3-vector.

        $$
            ab = (a_0 b_0 - \vec{a}\cdot\vec{b}; a_0 \vec{b} + b_0 \vec{a} + \vec{a} \times \vec{b})
        $$
        """
        vec_ = Quaternion(vec_)
        a1, b1, c1, d1 = self
        a2, b2, c2, d2 = vec_
        result = Quaternion([
                a1*a2 - b1*b2 - c1*c2 - d1*d2,
                a1*b2 + b1*a2 + c1*d2 - d1*c2,
                a1*c2 - b1*d2 + c1*a2 + d1*b2,
                a1*d2 + b1*c2 - c1*b2 + d1*a2])
        return result


# In[102]:

# implementations from https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
def asquaternion(vec_):
    if isinstance(vec_, Quaternion):
        return vec_
    else:
        return Quaternion(vec_)

def normalize(quaternion):
    """
    Normalize the quaternion. A normalized quaternion is of the form

    $$
        \cos(\theta/2) + \sin(\theta/2) \vec{v}/|\vec{v}|
    $$
    """
    quaternion = asquaternion(quaternion)
    cq = quaternion[0]
    sq = np.sqrt(1. - cq**2)
    quaternion[1:4] /= np.linalg.norm(quaternion[1:4])
    quaternion[1:4] *= sq
    return quaternion

def inverse(quaternion):
    """
    Returns the inverse of the quaternion.
    """
    result  = np.copy(Quaternion(quaternion)).view(np.ndarray)
    result *= (1, -1, -1, -1)
    return asquaternion(result)

def from_matrix(matrix_):
    """
    Sets the quaternion from the rotation matrix.
    """
    if matrix_.shape != (3,3):
        raise ValueError('Rotation matrix must be a (3,3) matrix.')
    xx, xy, xz = matrix_[:,0]
    yx, yy, yz = matrix_[:,1]
    zx, zy, zz = matrix_[:,2]
    K = 1./3.*np.array([
            [xx - yy - zz, yx + xy, zx + xz, yz - zy],
            [yx + xy, yy - xx - zz, zy + yz, zx - xz],
            [zx + xz, zy + yz, zz - xx - yy, xy - yx],
            [yz - zy, zx - xz, xy - yx, xx + yy + zz]])
    # the constructs a quaternion that is a closest fit to the
    # rotation in `matrix_`. If `matrix_` is pure rotation,
    # then `s[0]` is 1.
    U,s,V = np.linalg.svd(K, full_matrices=False)
    if not np.isclose(s[0], 1):
        raise ValueError('Rotation matrix {} is not a pure rotation.'.format(matrix_))
    (x, y, z, w) = U[:,0]
    return Quaternion((w, x, y, z))

def from_axis_angle(axis_, angle_):
    """
    Sets the quaternion from the axis/angle pair. Angle should be in degrees.
    """
    axis_ = np.asarray(axis_)/np.linalg.norm(axis_)
    # Note: quaternion uses half angles because rotation w' <- q w q*
    angle_ = np.radians(angle_)/2.
    cq = np.cos(angle_)
    sq = np.sin(angle_)
    return Quaternion(np.concatenate(((cq,), sq*axis_)))

def from_u_to_v(u_, v_):
    """
    Sets the quaternion to rotation `u` into `v`.
    """
    u_ = np.copy(u_)/np.linalg.norm(u_)
    v_ = np.copy(v_)/np.linalg.norm(v_)
    axis_ = np.cross(u_, v_)
    angle_ = np.degrees(np.arccos(np.dot(u_, v_)))
    return from_axis_angle(axis_, angle_)

def to_matrix(quaternion):
    """
    Converts the quaternion into an equivalent rotation matrix.
    """
    quaternion = asquaternion(quaternion)
    a, b, c, d = normalize(quaternion)
    return np.array([
            [a*a + b*b - c*c - d*d, 2*b*c - 2*a*d, 2*b*d + 2*a*c],
            [2*b*c + 2*a*d, a*a - b*b + c*c - d*d, 2*c*d - 2*a*b],
            [2*b*d - 2*a*c, 2*c*d + 2*a*b, a*a - b*b - c*c + d*d]])

def to_axis_angle(quaternion):
    """
    Converts the quaternion into an equivalent axis/angle pair.
    The angle is given in degrees. The vector is normalized.
    """
    angle_ = np.degrees(2.*np.arccos(quaternion[0]))
    axis_  = np.asarray(quaternion[1:4])
    axis_ /= np.linalg.norm(axis_)
    return (axis_, angle_)


# ## Test quaterion functions

if __name__ =='__main__':
    Rz = lambda q: np.array([
            [ np.cos(q),  np.sin(q), 0],
            [-np.sin(q),  np.cos(q), 0],
            [         0,          0, 1]])

    Rx = lambda q: np.array([
            [ 1,          0,          0],
            [ 0,  np.cos(q),  np.sin(q)],
            [ 0, -np.sin(q),  np.cos(q)]])

    def euler_matrix(phi, theta, psi):
        """
        All angles are in degrees.

        1. Rotation by angle phi about the z-axis
        2. Rotation by angle theta about the former x-axis (now x')
        3. Rotation by angle psi about the former z-axis (now z')
        """
        assert (0 <= theta <= 180.), \
            "Theta must lie between 0 and 180 degrees, inclusively."
        phi = np.radians(phi)
        theta = np.radians(theta)
        psi = np.radians(psi)
        return np.dot(Rz(psi), np.dot(Rx(theta), Rz(phi)))

    # check conversion from rotation matrix
    theta, phi, psi = 180*np.random.random(3)

    # create test vectors
    x = np.random.random(3)
    x /= np.linalg.norm(x)

    u,v = np.random.random((2, 3))
    u /= np.linalg.norm(u)
    v /= np.linalg.norm(v)

    axis, angle = np.random.random(3), 180*np.random.random()
    axis /= np.linalg.norm(axis)

    # create matrices
    R = euler_matrix(theta, phi, psi)
    Qm = from_matrix(R)
    Quv = from_u_to_v(u, v)
    Qaa = from_axis_angle(axis, angle)

    # test matrix rotation
    rxR = np.dot(R, x[:, np.newaxis]).flatten()
    # quaternion rotation
    rxQ = Qm*x*inverse(Qm)

    # verify
    print 'Quaternion rotation matches matrix rotation... ',
    if np.allclose(rxR, np.asarray(rxQ[1:])):
        print "pass"
    else:
        print "fail"

    # check conversion to rotation matrix
    print "Quaternion conversion to matrix...",
    if np.allclose(to_matrix(Qm), R):
        print "pass"
    else:
        print "fail"

    # TODO: test axis/angle rotation
    # create a rotation matrix from an axis/angle pair
    # perform axis/angle rotation
    # perform quaternion rotation
    # compare axis/angle rotation to equivalent quaternion rotation

    # check conversion to axis/angle
    print "Quaternion conversion to axis/angle...",
    caxis, cangle = to_axis_angle(Qaa)
    if np.isclose(cangle, angle) and np.allclose(caxis, axis):
        print "pass"
    else:
        print "fail: ({}, {}) != ({}, {})".format(
            caxis, cangle, axis, angle)
