"""Condition-space interpolation as NICEOPS does it (src/interpolator.cpp,
Database::blend_at), plus the alternatives under test."""
import numpy as np
from scipy.spatial import Delaunay, cKDTree

EPS10 = 10.0 * np.finfo(float).eps


class Conditions:
    def __init__(self, cond, extend=0.1, qhull_options=None):
        self.raw = cond.copy()
        self.lo = cond.min(0)
        self.hi = cond.max(0)
        self.norm = (cond - self.lo) / (self.hi - self.lo)
        # Extended vertex set: every point with any coordinate on 0 or 1 gets one
        # duplicate with those coordinates pushed to -extend / 1 + extend.
        verts, owner = [], []
        for i, p in enumerate(self.norm):
            verts.append(p.copy()); owner.append(i)
            onb = np.any((np.abs(p) < EPS10) | (np.abs(p - 1) < EPS10))
            if onb:
                q = p.copy()
                q[np.abs(p - 1) < EPS10] = 1 + extend
                q[np.abs(p) < EPS10] = -extend
                verts.append(q); owner.append(i)
        self.verts = np.array(verts)
        self.owner = np.array(owner)
        self.tri = Delaunay(self.verts, qhull_options=qhull_options)
        self.tree_ext = cKDTree(self.verts)
        self.tree = cKDTree(self.norm)
        # flight conditions (mu, gamma) and the CT reference copy for each
        fc = np.round(cond[:, :2], 9)
        self.fc, self.fc_of = np.unique(fc, axis=0, return_inverse=True)
        self.fc_of = self.fc_of.ravel()
        fcn = (self.fc - self.lo[:2]) / (self.hi[:2] - self.lo[:2])
        self.fc_tree = cKDTree(fcn)
        cts = np.unique(np.round(cond[:, 2], 12))
        self.ct_ref = cts[cts > 0].min()  # the LF = 1 copy
        self.fc_ref_index = np.full(len(self.fc), -1)
        for i in range(len(cond)):
            if abs(cond[i, 2] - self.ct_ref) < 1e-12:
                self.fc_ref_index[self.fc_of[i]] = i

    def normalize(self, q):
        return np.clip((q - self.lo) / (self.hi - self.lo), 0.0, 1.0)

    def barycentric(self, q):
        """q: (n,3) raw conditions -> (weights (n,4), indices (n,4)), NICEOPS rules."""
        x = self.normalize(q)
        s = self.tri.find_simplex(x)
        n = len(x)
        w = np.zeros((n, 4)); idx = np.zeros((n, 4), int)
        inside = s >= 0
        T = self.tri.transform[s[inside]]
        b = np.einsum('nij,nj->ni', T[:, :3, :], x[inside] - T[:, 3, :])
        raw = np.c_[b, 1 - b.sum(1)]
        verts = self.tri.simplices[s[inside]]
        w[inside] = raw
        idx[inside] = self.owner[verts]
        # a raw weight >= 1 collapses to that single vertex
        one = np.any(raw >= 1.0, axis=1)
        if one.any():
            k = np.argmax(raw[one] >= 1.0, axis=1)
            sub = np.nonzero(inside)[0][one]
            w[sub] = 0.0
            w[sub, 0] = 1.0
            idx[sub, 0] = self.owner[verts[one, k]]
            idx[sub, 1:] = idx[sub, :1]
        out = ~inside
        if out.any():
            _, v = self.tree_ext.query(x[out])
            w[out] = 0.0; w[out, 0] = 1.0
            idx[out] = self.owner[v][:, None]
        return w, idx, x

    def nearest(self, q):
        """Nearest stored condition in the normalized (mu, gamma, CT) space."""
        _, i = self.tree.query(self.normalize(q))
        return i

    def nearest_flight_condition(self, q):
        """Nearest (mu, gamma) sphere, reference CT copy, plus the exact CT law
        the database's own CT copies follow: +20 log10(CT/CT_ref) above LF 1."""
        x = self.normalize(q)
        _, f = self.fc_tree.query(x[:, :2])
        i = self.fc_ref_index[f]
        ct = np.clip(q[:, 2], self.lo[2], self.hi[2])
        corr = np.where(ct > self.ct_ref, 20 * np.log10(np.maximum(ct, 1e-30) / self.ct_ref), 0.0)
        return i, corr


def idw_weights(C, x, w, idx):
    """NORAH2-style weights (D1.5d eqs. 3-10): inverse normalised distance to
    the vertices of the enclosing simplex, to be applied to energies.  Uses the
    NICEOPS simplex and normalisation (no F_fc = 2 on gamma: the axes differ).
    Vertices with zero barycentric weight (padding) are excluded."""
    d = np.linalg.norm(C.norm[idx] - x[:, None, :], axis=2)
    use = w != 0.0
    hit = use & (d < 1e-12)
    with np.errstate(divide='ignore'):
        iw = np.where(use, 1.0 / np.maximum(d, 1e-300), 0.0)
    iw = np.where(hit.any(1)[:, None], hit.astype(float), iw)
    return iw / iw.sum(1)[:, None]


def blend_db(L, w, idx):
    """L: (..., n_cond) levels; one state's 4 weights/indices -> dB-linear blend."""
    return np.tensordot(L[..., idx], w, axes=([-1], [0]))


def blend_energy(L, w, idx):
    E = 10.0 ** (0.1 * L[..., idx])
    return 10.0 * np.log10(np.tensordot(E, w, axes=([-1], [0])))
