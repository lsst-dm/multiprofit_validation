#!/usr/bin/env python
# coding: utf-8

# In[1]:


# imports
import lsst.daf.butler as dafButler
from lsst.meas.extensions.multiprofit.rebuild_coadd_multiband import PatchCoaddRebuilder
from lsst.meas.extensions.multiprofit.plots import plot_blend
from lsst.pipe.base import QuantumGraph
import matplotlib.pyplot as plt
import numpy as np


# In[2]:


# define dataset
load_ci = False
if load_ci:
    skymap = "discrete/ci_imsim/4k"
    collection = "u/dtaranu/DM-42157/fit_ugrizy_merge"
    path_ticket = "/sdf/data/rubin/user/dtaranu/tickets/DM-42157/"
    repo = f"{path_ticket}/ci_imsim/DATA"
    tract = 0
    patch = 24
    matches = {
        "cModel": None,
        "ser_fixedcen": f"{path_ticket}/logs_ci_imsim/fit_src_ser_fixedcen_0_ugrizy.qgraph",
        "ser": f"{path_ticket}/logs_ci_imsim/fit_src_ser_{tract}_ugrizy.qgraph",
    }
else:
    skymap = "DC2"
    collection = "u/dtaranu/DM-42157/fit_ugrizy_merge"
    path_ticket = "/sdf/data/rubin/user/dtaranu/tickets/DM-42157/"
    repo = f"/repo/dc2"
    tract = 3828
    patch = 24
    matches = {
        "cModel": None,
        "ser_fixedcen": f"{path_ticket}/logs_testmed1/fit_src_ser_fixedcen_{tract}_ugrizy.qgraph",
        "ser": f"{path_ticket}/logs_testmed1/fit_src_ser_{tract}_ugrizy.qgraph",
    }

model_ref = "ser_fixedcen"


# In[3]:


# load data
butler = dafButler.Butler(repo, skymap=skymap, collections=[collection])
matches = {
    name: QuantumGraph.loadUri(path) if path is not None else None for name, path in matches.items()
}


# In[4]:


# PatchCoadd


# In[5]:


# make the rebuilder
rebuilder = PatchCoaddRebuilder.from_butler(
    butler=butler,
    skymap=skymap,
    tract=tract,
    patch=patch,
    collection_merged=collection,
    matches=matches,
    model_ref=model_ref,
    format_collection="{run}_match_{name}"
)


# In[6]:


# plot
ra_ref, dec_ref = (x*np.pi/180 for x in (56.53716733, -36.55586800))
dradec = 1.0/3600*np.pi/180

catalog_multi = rebuilder.matches[model_ref].rebuilder.catalog_multi
row_parent = np.where(
    np.hypot(catalog_multi["coord_ra"] - ra_ref, catalog_multi["coord_dec"] - dec_ref) < dradec
)[0][0]

kwargs_parent = dict(Q=5, minimum=-10, rgb_stretch_auto=True)
kwargs_children = dict(Q=5, minimum=-1, rgb_stretch_auto=True)

fig_rgb, ax_rgb, fig_gs, ax_gs = plot_blend(rebuilder, row_parent, kwargs_plot_parent=kwargs_parent, kwargs_plot_children=kwargs_children);


# In[6]:




