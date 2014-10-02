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

def fill_zeros(s,n):
    while len(s) < n:
        s= ''.join(('0',s))
    return s

global DB_path
DB_path= '/data/biophys/etournay/DB/'

class Movie:
    def __init__(self, name):
        self.name= name
        self.con= lite.connect(DB_path+name+'/'+name+'.sqlite')
        self.cells_loaded= False
        self.cellinfo_loaded= False
        self.dbonds_loaded= False
        self.Ta_t_loaded= False
        self.triList_loaded= False
        self.roiBT_loaded= False
    def load_cells(self):
        if self.cells_loaded == False:
            print('Loading cells from database...')
            self.cells= psql.read_sql('SELECT * FROM cells;', self.con)
            self.cells_loaded= True
            self.frames= sorted(self.cells['frame'].unique())
    def load_cellinfo(self):
        if self.cellinfo_loaded == False:
            print('Loading cellinfo from database...')
            self.cellinfo= psql.read_sql('SELECT * FROM cellinfo;', self.con)
            self.cellinfo_loaded= True
    def load_Ta_t(self):
        if self.Ta_t_loaded == False:
            print('Loading Ta_t...')
            if not os.path.isfile(DB_path+self.name+'/shear_contrib/Ta_t.csv'):
                print('Converting .RData to .csv...')
                ro.r('load("'+DB_path+self.name+'/shear_contrib/Ta_t.RData")')
                ro.r('write.csv(triList, "'+DB_path+self.name+'/shear_contrib/Ta_t.csv")')
                print('Converted!')
            self.Ta_t= pd.read_csv(DB_path+self.name+'/shear_contrib/Ta_t.csv')
            self.Ta_t_loaded= True
            del self.Ta_t[self.Ta_t.columns[0]]
    def load_roiBT(self):
        if self.roiBT_loaded == False:
            print('Loading roiBT from csv...')
            self.roiBT= pd.read_csv(DB_path+self.name+'/roi_bt/lgRoiSmoothed.csv')
            self.roiBT_loaded= True
            self.regions= self.roiBT['roi'].unique()
    def load_triList(self):
        if self.triList_loaded == False:
            print('Loading triList...')
            if not os.path.isfile(DB_path+self.name+'/shear_contrib/triList.csv'):
                print('Converting .RData to .csv...')
                ro.r('load("'+DB_path+self.name+'/shear_contrib/triList.RData")')
                ro.r('write.csv(triList, "'+DB_path+self.name+'/shear_contrib/triList.csv")')
                print('Converted!')
            self.triList= pd.read_csv(DB_path+self.name+'/shear_contrib/triList.csv')
            self.triList_loaded= True
            del self.triList[self.triList.columns[0]]
    def region_cells(self, roi_name):
        self.load_cells()
        self.load_roiBT()
        if roi_name not in self.regions:
            raise Exception('Region '+roi_name+' is not defined in this movie!')
        else:
            return self.cells[self.cells['cell_id'].isin(self.roiBT[self.roiBT['roi']==roi_name]['cell_id'])]
    def region_area(self, roi_name):
        rc= self.region_cells(roi_name)
        return np.array([np.sum(rc[rc['frame']==f]['area']) for f in self.frames])
    def region_shape_nematic(self, roi_name):
        rc= self.region_cells(roi_name)
        print('Calculating '+roi_name+' shape nematic...')
        av_x= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['area']) for f in self.frames])
        av_y= np.array([np.sum(rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in self.frames])
        av_xx= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['area']) for f in self.frames])
        av_xy= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in self.frames])
        av_yy= np.array([np.sum(rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in self.frames])
        del rc
        ra= self.region_area(roi_name)
        m_xx= (av_xx-av_x*av_x/ra)/ra**2.
        m_yy= (av_yy-av_y*av_y/ra)/ra**2.
        m_xy= (av_xy-av_x*av_y/ra)/ra**2.
        s= 0.5*np.log(m_xx*m_yy - m_xy**2)
        Q = np.arcsinh(0.5*np.sqrt((m_xx-m_yy)**2.+(2*m_xy)**2.)/np.exp(s))
        twophi = np.arctan2((2*m_xy),(m_xx-m_yy))
        return Q*np.cos(twophi), Q*np.sin(twophi), s
    def region_center(self, roi_name):
        rc= self.region_cells(roi_name)
        av_x= np.array([np.sum(rc[rc['frame']==f]['center_x']*rc[rc['frame']==f]['area']) for f in self.frames])
        av_y= np.array([np.sum(rc[rc['frame']==f]['center_y']*rc[rc['frame']==f]['area']) for f in self.frames])
        ra= self.region_area(roi_name)
        av_x, av_y= av_x/ra, av_y/ra
        del rc
        return av_x, av_y
    def region_mean_length_height(self, roi_name):
        Q1, Q2, s= self.region_shape_nematic(roi_name)
        print('Calculating '+roi_name+' mean length and height...')
        return np.sqrt(self.region_area(roi_name))*np.exp(0.5*Q1), np.sqrt(self.region_area(roi_name))*np.exp(-0.5*Q1)


m= Movie('WT_25deg_111102')
m.load_triList()
m.load_Ta_t()
m.load_cells()

m.Ta_t_list= m.Ta_t.merge(m.triList, left_on= 'tri_id', right_on= 'tri_id')

grouped= m.Ta_t_list.groupby('frame')
grouped_Q= grouped['Q_a']


for name, group in grouped_Q:
    print name
    group.hist(bins=np.linspace(0,.9,200), normed=True)
    plt.ylim(0,3.5)
    plt.xlabel(r'$|Q|$')
    plt.savefig('figures/histograms/Q_'+fill_zeros(str(name),4))
    plt.close()
sigma= [group.std() for name, group in grouped_Q]
mean= [group.mean() for name, group in grouped_Q]
for name, group in grouped_Q:
    sigma.append()

compare3= pd.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/compare3_elongation_correlation.csv')
compare3.columns
frames= np.array(compare3['frame'])
crc_xx_lower= np.array(compare3['crc_xx.x'])
area_lower= np.array(compare3['tri_area.x'])
crc_xx_higher= np.array(compare3['crc_xx.y'])
area_higher= np.array(compare3['tri_area.y'])
crc_xx_total= np.array(compare3['crc_xx'])
area_total= np.array(compare3['tri_area.i1'])
def smoothPlot(x, y, *args, **kwargs):
  kernel = np.ones(Nsmooth)/Nsmooth
  plt.plot(np.convolve(x,kernel,'valid'), np.convolve(y,kernel,'valid'), *args, **kwargs)
Nsmooth= 10

plt.figure()
smoothPlot(frames, 13*crc_xx_lower*area_lower/area_total, label= 'lower Q')
smoothPlot(frames, 13*crc_xx_higher*area_higher/area_total, label= 'higher Q')
smoothPlot(frames, 13*crc_xx_total, label= 'total')
plt.legend(loc='best')
plt.savefig('figures/crc_elongation_division.png')
plt.show()



plt.figure()
plt.plot(np.array(sigma))
plt.plot(np.array(mean))
plt.show()
"""
m.
tri1= m.Ta_t[m.Ta_t['frame']==1]
tri1.head()
np.histogram(tri1[])
del m.triList[m.triList.columns[0]]
blade_L, blade_h= m.region_mean_length_height('blade')
hinge_L, hinge_h= m.region_mean_length_height('hinge')
otp_df= pd.DataFrame(pd.Series(blade_L))
otp_df.columns= ['blade_L']
otp_df['blade_h']= pd.Series(blade_h)
otp_df['hinge_L']= pd.Series(hinge_L)
otp_df['hinge_h']= pd.Series(hinge_h)

otp_df.to_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/drosophila_data_library/height_length_data/'+m.name+'.csv')

blade_av_x, blade_av_y= m.region_center('blade')
hinge_av_x, hinge_av_y= m.region_center('hinge')
blade_cells= m.region_cells('blade')
hinge_cells= m.region_cells('hinge')

for frame in m.frames:
    print frame
    im_path= DB_path+'/'+m.name+'/image_data/mutant/tag/segmentationData/frame'+fill_zeros(str(frame),4)+'/original_trafo.png'
    im= plt.imread(im_path)
    y_size, x_size= im.shape[:-1]
    frame_blade_x= blade_cells[blade_cells['frame']==frame]['center_x']
    frame_blade_y= blade_cells[blade_cells['frame']==frame]['center_y']
    frame_hinge_x= hinge_cells[hinge_cells['frame']==frame]['center_x']
    frame_hinge_y= hinge_cells[hinge_cells['frame']==frame]['center_y']
    plt.figure()
    b_rect= patches.Rectangle((blade_av_x[frame]-0.5*blade_L[frame],blade_av_y[frame]-0.5*blade_h[frame]),blade_L[frame],blade_h[frame],fill=None, edgecolor='white', linewidth=2)
    h_rect= patches.Rectangle((hinge_av_x[frame]-0.5*hinge_L[frame],hinge_av_y[frame]-0.5*hinge_h[frame]),hinge_L[frame],hinge_h[frame],fill=None, edgecolor='green', linewidth=2)
    plt.imshow(im)
    plt.gca().add_patch(b_rect)
    plt.gca().add_patch(h_rect)
    plt.scatter(frame_blade_x, frame_blade_y, c='purple', s=15, linewidths=0)
    plt.scatter(frame_hinge_x, frame_hinge_y, c='yellow', s=15, linewidths=0)
    plt.axis('off')
    plt.savefig('figures/blade_size_frame'+fill_zeros(str(frame),4)+'.png', bbox_inches='tight')
    plt.close()

plt.figure()
plt.plot(m.frames, 0.5*Q1)
plt.plot(m.frames, 0.5*Q2)
plt.show()

np.log(a_xx*a_yy-a_xy**2)
a= m.region_mean_length_height('blade')
L2_cells.columns
m.regions
c= np.matrix([[1,2],[3,4]])
c[0,0]
c/3.**2
d= np.arange(10)
c= np.matrix([[d,2*d],[3*d,4*d]])"""