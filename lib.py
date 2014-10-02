import sqlite3 as lite
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas.io.sql as psql
import pandas.io.parsers as pp
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com

def connect_DB(DB_path, movie_name):
    return lite.connect(DB_path+movie_name+'/'+movie_name+'.sqlite')

def query_cells(DB_con):
    cells= psql.read_sql('SELECT * FROM cells;', DB_con)


connect_DB('/data/biophys/etournay/DB/', 'WT_25deg_111102')