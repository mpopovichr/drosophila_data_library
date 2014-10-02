from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
import sqlite3 as lite
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas.io.sql as psql
import pandas.io.parsers as pp
import matplotlib.image as mpimg
import scipy.signal as signal
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

def fill_zeros(s,k):
    while len(s)< k:
        s = '0'+s
    return s



avgDeltaQtotDivision= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotDivision.csv')
avgDeltaQtotNoDivision= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotNoDivision.csv')
avgDeltaQtotBlade= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotBlade.csv')

avgDeltaQtotDivision.columns
avgDeltaQtotDivision['QxyThetaU']
(avgDeltaQtotDivision['ThetaU']*(avgDeltaQtotNoDivision['Q_xy.i1']+avgDeltaQtotDivision['Q_xy.i2']))
avgDeltaQtotDivision['ThetaU']*(avgDeltaQtotNoDivision['Q_xy.i1']+avgDeltaQtotDivision['Q_xy.i2'])