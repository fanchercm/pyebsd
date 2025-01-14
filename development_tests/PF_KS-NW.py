import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import pyebsd


T_KS = pyebsd.OR()  # K-S
T_NW = pyebsd.OR(ds=[[1, 1, 0], [1, 0, 0]])  # N-W

proj = [1, 1, 1]

ax = pyebsd.plot_PF(M=T_KS, mfc=[1, 1, 1, 0], mec='r', mew=1,
                    marker='s', ms=5, proj=proj, label='KS')
pyebsd.plot_PF(M=T_NW, mfc=[1, 1, 1, 0], mec='g', mew=1,
               marker='D', ms=5, proj=proj, label='NW', ax=ax)
pyebsd.draw_std_traces(ax, lw=.2)
pyebsd.draw_wulff_net(ax, lw=.2)

ax.legend()
plt.show()
