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
from scipy import optimize

import lib


###### AREA ELONGATION CORRELATION ANALYSIS ##############################

m= {}

name_list= ['WT_25deg_111102',
            'WT_25deg_111103',
            'WT_25deg_120531',
            'WT_25-30deg_130921',
            'MTdp_25deg_140222',
            'WT_sevBdist-25deg_130131',
            'WT_antLinkCut-25deg_131227',
            'HTcdc2_25-30deg_130927',
            'MTcdc2_25-30deg_130919',
            'MTcdc2_25-30deg_130917',
            'MTcdc2_25-30deg_130916',
            'MTcdc2_25deg_130905']

for name in name_list:
    m[name]= lib.Movie(name)
    m[name].cell_area_avg= {}
    m[name].Qxx_avg, m[name].Qxy_avg, m[name].Q_avg= {}, {} ,{}

region= 'blade'
for name in name_list:
    print name
    m[name].load_roiBT()
    m[name].cell_area_avg[region]= lib.region_cells_area_avg(m[name].region_cells(region))
    m[name].Qxx_avg[region], m[name].Qxy_avg[region], m[name].Q_avg[region]= lib.region_cells_shape_avg(m[name].region_cells(region))

global ref_cell_area
ref_cell_area= m['WT_25deg_111102'].cell_area_avg['blade'][0]
print('Done!')

markers = ['o','v','^','D','s']
for region in m['WT_25deg_111102'].regions:
    print region
    plt.figure()
    color_count= 0
    for name in name_list:
        if region in m[name].regions:
            color_frac= 1.0*color_count/len(name_list)
            color_count +=1
            plt.scatter(np.log(m[name].cell_area_avg[region][-1]/ref_cell_area), m[name].Qxx_avg[region][-1],c=plt.cm.hsv(color_frac), s=50, label=name, marker= markers[np.mod(color_count,5)])
    plt.xlabel('Q_xx+Q_yy')
    plt.ylabel('(Q_xx-Q_yy)/2')
    plt.title(region)
    plt.legend(loc=4, fontsize=6)
    plt.savefig('figures/area_elongation_'+region+'_Q_xx-Q_yy.pdf')
    plt.close()

markers = ['o','v','^','D','s']
for region in m['WT_25deg_111102'].regions:
    print region
    plt.figure()
    color_count= 0
    for name in name_list:
        if region in m[name].regions:
            color_frac= 1.0*color_count/len(name_list)
            color_count +=1
            plt.scatter(np.log(m[name].cell_area_avg[region][-1]/ref_cell_area), m[name].Q_avg[region][-1],c=plt.cm.hsv(color_frac), s=50, label=name, marker= markers[np.mod(color_count,5)])
    plt.xlabel('Q_xx+Q_yy')
    plt.ylabel('|Q|')
    plt.title(region)
    plt.legend(loc=4, fontsize=6)
    plt.savefig('figures/area_elongation_'+region+'_Q.pdf')
    plt.close()


markers = ['o','v','^','D','s']
linestyles= ['-','--', '-.', ':']
plt.figure()
color_count= 0
for region in m['WT_25deg_111102'].regions:
    print region
    color_frac= 1.0*color_count/len(name_list)
    color_count +=1
    x,y = [], []
    for name in name_list:
        if region in m[name].regions:
            x.append(np.log(m[name].cell_area_avg[region][-1]/ref_cell_area))
            y.append(m[name].Qxx_avg[region][-1])
    coeff= np.polyfit(x,y,1)
    fit_x= np.linspace(-1, 1, 100)
    pol= np.poly1d(coeff)
    plt.plot(fit_x, pol(fit_x), c=plt.cm.hsv(color_frac), label= region+': '+"{:.4f}".format(coeff[0]), ls= linestyles[np.mod(color_count,4)])
plt.xlabel('Q_xx+Q_yy')
plt.ylabel('(Q_xx-Q_yy)/2')
plt.title('all')
plt.legend(loc=2, fontsize=6)
plt.savefig('figures/area_elongation_all_by_region_Q_xx-Q_yy.pdf')
plt.close()

markers = ['o','v','^','D','s']
linestyles= ['-','--', '-.', ':']
plt.figure()
color_count= 0
for region in m['WT_25deg_111102'].regions:
    print region
    color_frac= 1.0*color_count/len(name_list)
    color_count +=1
    x,y = [], []
    for name in name_list:
        if region in m[name].regions:
            x.append(np.log(m[name].cell_area_avg[region][-1]/ref_cell_area))
            y.append(m[name].Q_avg[region][-1])
    coeff= np.polyfit(x,y,1)
    fit_x= np.linspace(-1, 1, 100)
    pol= np.poly1d(coeff)
    plt.plot(fit_x, pol(fit_x), c=plt.cm.hsv(color_frac), label= region, ls= linestyles[np.mod(color_count,4)])
plt.xlabel('Q_xx+Q_yy')
plt.ylabel('|Q|')
plt.title('all')
plt.legend(loc=2, fontsize=6)
plt.savefig('figures/area_elongation_all_by_region_Q.pdf')
plt.close()

markers = ['o','v','^','D','s']
for name in name_list:
    print name
    plt.figure()
    color_count= 0
    for region in m[name].regions:
        color_frac= 1.0*color_count/len(m[name].regions)
        color_count +=1
        plt.scatter(np.log(m[name].cell_area_avg[region][-1]/ref_cell_area), m[name].Qxx_avg[region][-1],c=plt.cm.hsv(color_frac), s=50, label=region, marker= markers[np.mod(color_count,5)])
    plt.xlabel('Q_xx+Q_yy')
    plt.ylabel('(Q_xx-Q_yy)/2')
    plt.title(region)
    plt.legend(loc=4, fontsize=6)
    plt.savefig('figures/area_elongation_'+name+'_Q_xx-Q_yy.pdf')
    plt.close()

linestyles= ['-','--', '-.', ':']
plt.figure()
color_count= 0
for name in name_list:
    print name
    color_frac= 1.0*color_count/len(name_list)
    color_count +=1
    x,y = [], []
    for region in m[name].regions:
        x.append(np.log(m[name].cell_area_avg[region][-1]/ref_cell_area))
        y.append(m[name].Qxx_avg[region][-1])
    coeff= np.polyfit(x,y,1)
    fit_x= np.linspace(-1, 1, 100)
    pol= np.poly1d(coeff)
    plt.plot(fit_x, pol(fit_x), c=plt.cm.hsv(color_frac), label= name+': '+"{:.4f}".format(coeff[0]), ls= linestyles[np.mod(color_count,4)])
plt.xlabel('Q_xx+Q_yy')
plt.ylabel('(Q_xx-Q_yy)/2')
plt.title('all')
plt.legend(loc=2, fontsize=6)
plt.savefig('figures/area_elongation_all_by_movie_Q_xx-Q_yy.pdf')
plt.close()

x,y= [], []
for name in name_list:
    for region in m[name].regions:
        x.append(np.log(m[name].cell_area_avg[region][-1]/ref_cell_area))
        y.append(m[name].Qxx_avg[region][-1])
plt.scatter(x,y)
coeff= np.polyfit(x,y,1)
fit_x= np.linspace(-1, 1, 100)
pol= np.poly1d(coeff)
plt.plot(fit_x, pol(fit_x), c='red', label= name)
plt.xlabel('Q_xx+Q_yy')
plt.ylabel('(Q_xx-Q_yy)/2')
plt.title('all')
plt.savefig('figures/area_elongation_all_together_fit_Q_xx-Q_yy.pdf')
plt.close()
plt.show()

name= 'WT_25deg_111103'
region= 'blade'
for region in m[name].regions:
    print region
    cells= m[name].region_cells(region)
    last_cells= cells[cells['frame']==m[name].frames[-1]]
    x= np.log(np.array(last_cells['area'])/ref_cell_area)
    y= np.array(last_cells['elong_xx'])
    plt.figure()
    plt.scatter(x,y, linewidths=0)
    coeff= np.polyfit(x,y,1)
    fit_x= np.linspace(-1, 1, 100)
    pol= np.poly1d(coeff)
    plt.plot(fit_x, pol(fit_x), c='red', label= name, linewidth=2)
    plt.title(region)
    plt.savefig('figures/area_elongation_individual_cells_'+name+'_'+region+'_Q_xx-Q_yy.pdf')
    plt.close()

name= 'WT_25deg_111103'
region= 'blade'
for region in m[name].regions:
    print region
    cells= m[name].region_cells(region)
    last_cells= cells[cells['frame']==m[name].frames[-1]]
    x= np.log(np.array(last_cells['area'])/ref_cell_area)
    y1= np.array(last_cells['elong_xx'])
    y2= np.array(last_cells['elong_xy'])
    y= np.sqrt(y1**2+y2**2)
    plt.figure()
    plt.scatter(x,y, linewidths=0)
    coeff= np.polyfit(x,y,1)
    fit_x= np.linspace(-1, 1, 100)
    pol= np.poly1d(coeff)
    plt.plot(fit_x, pol(fit_x), c='red', label= name, linewidth=2)
    plt.title(region)
    plt.savefig('figures/area_elongation_individual_cells_'+name+'_'+region+'_Q.pdf')
    plt.close()
