import sqlite3 as lite
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import pandas.io.sql as psql
import pandas.io.parsers as pp
import rpy2.robjects as robjects
import rpy2.robjects as ro
import pandas.rpy.common as com
import os.path
from PyQt4 import QtGui


def region_area(rc):
    frames= sorted(rc['frame'].unique())
    return np.array([np.sum(rc[rc['frame']==f]['area']) for f in frames])
def region_shape_nematic( rc):
    frames= sorted(rc['frame'].unique())
    av_x= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['area']) for f in frames])
    av_y= np.array([np.sum(rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in frames])
    av_xx= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['area']) for f in frames])
    av_xy= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in frames])
    av_yy= np.array([np.sum(rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in frames])
    ra= region_area(rc)
    m_xx= (av_xx-av_x*av_x/ra)/ra**2.
    m_yy= (av_yy-av_y*av_y/ra)/ra**2.
    m_xy= (av_xy-av_x*av_y/ra)/ra**2.
    s= 0.5*np.log(m_xx*m_yy - m_xy**2)
    Q = np.arcsinh(0.5*np.sqrt((m_xx-m_yy)**2.+(2*m_xy)**2.)/np.exp(s))
    twophi = np.arctan2((2*m_xy),(m_xx-m_yy))
    return Q*np.cos(twophi), Q*np.sin(twophi), s
def region_center( rc):
    frames= sorted(rc['frame'].unique())
    av_x= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['area']) for f in frames])
    av_y= np.array([np.sum(rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in frames])
    ra= region_area(rc)
    av_x, av_y= av_x/ra, av_y/ra
    return av_x, av_y
def region_mean_length_height( rc):
    Q1, Q2, s= region_shape_nematic(rc)
    ra= region_area(rc)
    return np.sqrt(ra)*np.exp(0.5*Q1), np.sqrt(ra)*np.exp(-0.5*Q1)