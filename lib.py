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
from sets import Set
import os.path
from PyQt4 import QtGui
import time

global Path
Path= '/data/biophys/etournay/'
#Path= '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/'
global DB_path
DB_path= Path+'DB/'



def smoothPlot(x, y, Nsmooth=10, *args, **kwargs):
  kernel = np.ones(Nsmooth)/Nsmooth
  plt.plot(np.convolve(x,kernel,'valid'), np.convolve(y,kernel,'valid'), *args, **kwargs)
def fill_zeros(s,n):
    while len(s) < n:
        s= ''.join(('0',s))
    return s
def region_area(rc):
    frames= sorted(rc['frame'].unique())
    return np.array([np.sum(rc[rc['frame']==f]['area']) for f in frames])
def region_cells_area_avg(rc):
    frames= sorted(rc['frame'].unique())
    return np.array([rc[rc['frame']==f]['area'].mean() for f in frames])
def region_cells_shape_avg(rc):
    frames= sorted(rc['frame'].unique())
    Q_xx= np.array([rc[rc['frame']==f]['elong_xx'].mean() for f in frames])
    Q_xy= np.array([rc[rc['frame']==f]['elong_xy'].mean() for f in frames])
    rc['Q']= np.sqrt(rc['elong_xx']**2 + rc['elong_xy']**2)
    Q= np.array([rc[rc['frame']==f]['Q'].mean() for f in frames])
    return Q_xx, Q_xy, Q
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
def show_region(m, rc_list, frame):
    im_path= DB_path+'/'+m.name+'/image_data/mutant/tag/segmentationData/frame'+fill_zeros(str(frame),4)+'/original_trafo.png'
    im= plt.imread(im_path)
    plt.figure()
    plt.imshow(im)
    c= 0
    for rc in rc_list:
        rc_frame= rc[rc['frame']==frame]
        plt.scatter(rc_frame['center_x'],rc_frame['center_y'], color= 'red', s=rc_frame['area']/30.)
        c+=1
    plt.show()
    plt.close()
def region_symmetric_difference(ra, rb):
    df= pd.concat([ra,rb])
    df= df.reset_index(drop=True)
    df_gpby= df.groupby(list(df.columns))
    idx= [x[0] for x in df_gpby.groups.values() if len(x) == 1]
    return df.reindex(idx)
def film_region(m, rc, dir_name):
    if not os.path.exists(dir_name): os.makedirs(dir_name)
    for frame in m.frames:
        print frame
        im_path= DB_path+'/'+m.name+'/image_data/mutant/tag/segmentationData/frame'+fill_zeros(str(frame),4)+'/original_trafo.png'
        im= plt.imread(im_path)
        plt.figure()
        plt.imshow(im)
        rc_frame= rc[rc['frame']==frame]
        size= np.sqrt(rc_frame['area'])/15.
        plt.scatter(rc_frame['center_x'], rc_frame['center_y'], color='red', s= size)
        plt.tight_layout()
        plt.savefig(dir_name+'/frame'+fill_zeros(str(frame),4)+'.png')
        plt.close()
def square_difference(a, b):
    return np.sqrt(np.sum((np.array(a)-np.array(b))**2))
def quick_plot(data):
    plt.figure()
    for d in data:
        plt.plot(d)
    plt.show()
def smooth_data(x, NSmooth=10):
    return np.convolve(x, 1.*np.ones(NSmooth)/NSmooth, mode='valid')



class Movie:
    def __init__(self, name):
        self.name= name
        self.con= lite.connect(DB_path+name+'/'+name+'.sqlite')
        self.loaded= Set()
        pupalWings= pp.read_csv(DB_path+'PupalWingMovies.csv', sep='\t')
        if not np.isnan(np.array(pupalWings[pupalWings['nice_name']==name]['time_shift_sec'])[0]):
            self.time= np.array(self.region_deform_tensor('blade')['time_sec'])/3600.+15. \
                   + np.array(pupalWings[pupalWings['nice_name']==name]['time_shift_sec'])[0]/3600.
            self.dt= self.time[1:] - self.time[:-1]
        else:
            print('While loading '+ name+'. No time shift provided for this movie! No time will be loaded.')
    def load_cells(self):
        if not 'cells' in self.loaded:
            print('Loading cells from database...')
            self.cells= psql.read_sql('SELECT * FROM cells WHERE cell_id > 10000;', self.con)
            self.loaded.add('cells')
            self.frames= sorted(self.cells['frame'].unique())
    def RData_to_csv(self, table, file_RData, file_csv):
        if (not os.path.isfile(file_csv)) or (time.ctime(os.path.getmtime(file_RData)) > time.ctime(os.path.getmtime(file_csv))):
            print('Converting .RData to .csv...')
            ro.r('load("'+file_RData+'")')
            ro.r('write.csv('+table+', "'+file_csv+'")')
            print('Converted!')
    def load_cellinfo(self):
        if not 'cellinfo' in self.loaded:
            print('Loading cellinfo from database...')
            self.cellinfo= psql.read_sql('SELECT * FROM cellinfo WHERE cell_id > 10000;', self.con)
            self.loaded.add('cellinfo')
    def load_triTracked(self):
        if not 'triTracked' in self.loaded:
            file_csv=DB_path+self.name+'/shear_contrib/triTracked.csv'
            file_RData=DB_path+self.name+'/shear_contrib/triTracked.RData'
            print('Loading triTracked...')
            self.RData_to_csv('triTracked',file_RData, file_csv)
            print('Done!')
            self.triTracked= pd.read_csv(file_csv)[['frame','cell_a','cell_b','cell_c']]
            self.loaded.add('triTracked')
    def load_triCategories(self):
        if not 'triCategories' in self.loaded:
            file_csv=DB_path+self.name+'/tri_categories/triangleCategories.csv'
            file_RData=DB_path+self.name+'/tri_categories/triangleCategories.RData'
            print('Loading triCategories...')
            self.RData_to_csv('triCategories',file_RData, file_csv)
            print('Done!')
            self.triCategories= pd.read_csv(file_csv)
            self.loaded.add('triCategories')
    def load_Ta_t(self):
        if not 'Ta_t' in self.loaded:
            file_csv= DB_path+self.name+'/shear_contrib/Ta_t.csv'
            file_RData= DB_path+self.name+'/shear_contrib/Ta_t.RData'
            print('Loading Ta_t...')
            self.RData_to_csv('Ta_t',file_RData, file_csv)
            print('Done!')
            self.Ta_t= pd.read_csv(file_csv)
            self.loaded.add('Ta_t')
            del self.Ta_t[self.Ta_t.columns[0]]
    def load_roiBT(self):
        if not 'roiBT' in self.loaded:
            file_csv= DB_path+self.name+'/roi_bt/lgRoiSmoothed.csv'
            file_RData= DB_path+self.name+'/roi_bt/lgRoiSmoothed.RData'
            print('Loading roiBT ...')
            self.RData_to_csv('lgRoiSmoothed', file_RData, file_csv)
            print('Done!')
            self.roiBT= pd.read_csv(DB_path+self.name+'/roi_bt/lgRoiSmoothed.csv')
            self.loaded.add('roiBT')
            self.regions= self.roiBT['roi'].unique()
    def load_triList(self):
        if self.triList_loaded == False:
            file_csv= DB_path+self.name+'/shear_contrib/triList.csv'
            file_RData= DB_path+self.name+'/shear_contrib/triList.RData'
            print('Loading triList...')
            self.RData_to_csv('triList', file_RData, file_csv)
            print('Done!')
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
    def region_deform_tensor(self, roi_name):
        if 'avgDeformTensorsWide.tsv' in os.listdir(DB_path+self.name+'/shear_contrib/'+roi_name):
            df_DB_shear= pp.read_csv(DB_path+self.name+'/shear_contrib/'+roi_name+'/avgDeformTensorsWide.tsv', sep='\t')
        else:
            ro.r('load("'+DB_path+self.name+'/shear_contrib/'+roi_name+'/avgDeformTensorsWide.RData")')
            df_DB_shear= com.load_data('avgDeformTensorsWide')
        return df_DB_shear
    def load_PIV_whole_wing(self, piv_region):
        if not os.path.exists(Path+'PIV/'+self.name):
            print('PIV does not exists for '+self.name)
        else:
            return pp.read_csv(Path+'PIV/'+self.name+'/Segmentation/rotated_pivData_'+piv_region+'.csv')
    def whole_hinge_blade(self):
        ww= self.region_cells('whole_wing')
        HBint= self.region_cells('HBinterface')
        HB_x= np.array([HBint[HBint['frame']==f]['center_x'].mean() for f in self.frames])
        ww_hinge, ww_blade= pd.DataFrame(), pd.DataFrame()
        for frame in self.frames:
            ww_frame= ww[ww['frame']==frame]
            ww_hinge= ww_hinge.append(ww_frame[ww_frame['center_x']<HB_x[frame]])
            ww_blade= ww_blade.append(ww_frame[ww_frame['center_x']>HB_x[frame]])
        return ww_hinge, ww_blade
    def fancy_hinge_blade(self):
        ww= self.region_cells('whole_wing')
        HBint= self.region_cells('HBinterface')
        print('Determining fancy hinge and blade ...')
        f_hinge, f_blade= pd.DataFrame(), pd.DataFrame()
        for frame in self.frames:
            print(str(100.*frame/len(self.frames))+'%')
            ww_frame= ww[ww['frame']==frame]
            HB_frame= HBint[HBint['frame']==frame]
            HB_xmin= HB_frame['center_x'].min()
            HB_xmax= HB_frame['center_x'].max()
            trivial_hinge= ww_frame[ww_frame['center_x']<HB_xmin]
            trivial_blade= ww_frame[ww_frame['center_x']>HB_xmax]
            central_part= ww_frame[(ww_frame['center_x']>HB_xmin)&(ww_frame['center_x']<HB_xmax)]
            coeff= np.polyfit(np.array(HB_frame['center_y']), np.array(HB_frame['center_x']), 4)
            polynomial= np.poly1d(coeff)
            central_part['interface_x']=central_part['center_y'].apply(polynomial)
            central_hinge= central_part[central_part['center_x']<central_part['interface_x']]
            central_blade= central_part[central_part['center_x']>central_part['interface_x']]
            del central_hinge['interface_x']
            del central_blade['interface_x']
            f_hinge= f_hinge.append(trivial_hinge)
            f_hinge= f_hinge.append(central_hinge)
            f_blade= f_blade.append(trivial_blade)
            f_blade= f_blade.append(central_blade)
        return f_hinge, f_blade
    def tri_tracking(self):
        if not self.triTracked:
            self.load_triTracked()
        self.triTracked['tri_hash_a']= (
        self.triTracked['cell_a']*10**10
            + self.triTracked['cell_b']*10**5
            + self.triTracked['cell_c']
            )
        self.triTracked['tri_hash_b']= (
            self.triTracked['cell_b']*10**10
            + self.triTracked['cell_c']*10**5
            + self.triTracked['cell_a']
            )
        self.triTracked['tri_hash_c']= (
            self.triTracked['cell_c']*10**10
            + self.triTracked['cell_a']*10**5
            + self.triTracked['cell_b']
            )
        self.triTracked['tri_hash']= (
            self.triTracked[['tri_hash_a',
                             'tri_hash_b',
                             'tri_hash_c']].min(axis=1)
            )
        self.triTracked= self.triTracked.sort(['tri_hash', 'frame'])
        self.triTracked= self.triTracked.reset_index()
        self.triTracked= self.triTracked[['frame', 'tri_hash','type']]
    def tri_track_first_occ(self, d):
        last_from= [0]
        cd= np.cumsum(d[1:])
        for x in d[1:]:
            if x < 2:
                last_from.append(last_from[-1])
            else:
                last_from.append(cd[len(last_from)-1])
        return np.array(last_from)
    def tri_track_last_occ(self, g, d):
        last_occ=np.append(
            np.array([0 if np.abs(x) < 2 else 1 for x in d[:-1]])*(np.array(g)[:-1]),
            np.array(g)[-1])
        return np.array(pd.Series(last_occ).replace(0,method='bfill'))
    def shear_corrected_dimensions_HB(self):
        self.hinge_fcy, self.blade_fcy= self.fancy_hinge_blade()
        self.blade_L, self.blade_h= region_mean_length_height(self.blade_fcy)
        self.hinge_L, self.hinge_h= region_mean_length_height(self.hinge_fcy)
        self.blade_piv, self.hinge_piv= self.load_PIV_whole_wing('blade_only'), self.load_PIV_whole_wing('hinge_only')
        self.hinge_piv_vxx= np.array(self.hinge_piv.groupby('frame')['Vxx'].mean()*3600.)
        self.hinge_piv_vyy= np.array(self.hinge_piv.groupby('frame')['Vyy'].mean()*3600.)
        self.blade_piv_vxx= np.array(self.blade_piv.groupby('frame')['Vxx'].mean()*3600.)
        self.blade_piv_vyy= np.array(self.blade_piv.groupby('frame')['Vyy'].mean()*3600.)


          
        


