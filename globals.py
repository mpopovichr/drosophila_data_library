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

import lib

m= {}
m['WT']= lib.Movie('WT_25deg_111103')
m['ALC']= lib.Movie('WT_antLinkCut-25deg_131227')
m['DC']= lib.Movie('WT_sevBdist-25deg_130131')
m['DLC']= lib.Movie('WT_distLinkCut-25deg_131226')

pm= 'DLC'

m[pm].hinge_fcy, m[pm].blade_fcy= m[pm].fancy_hinge_blade()
m[pm].hinge_old, m[pm].blade_old= m[pm].whole_hinge_blade()

m[pm].hinge_fcy_L, m[pm].hinge_fcy_h= lib.region_mean_length_height(m[pm].hinge_fcy)
m[pm].blade_fcy_L, m[pm].blade_fcy_h= lib.region_mean_length_height(m[pm].blade_fcy)
m[pm].hinge_old_L, m[pm].hinge_old_h= lib.region_mean_length_height(m[pm].hinge_old)
m[pm].blade_old_L, m[pm].blade_old_h= lib.region_mean_length_height(m[pm].blade_old)
m[pm].hinge_seg_L, m[pm].hinge_seg_h= lib.region_mean_length_height(m[pm].region_cells('hinge'))
m[pm].blade_seg_L, m[pm].blade_seg_h= lib.region_mean_length_height(m[pm].region_cells('blade'))

m[pm].blade_piv, m[pm].hinge_piv= m[pm].load_PIV_whole_wing('blade_only'), m[pm].load_PIV_whole_wing('hinge_only')
m[pm].hinge_piv_vxx= np.array(m[pm].hinge_piv.groupby('frame')['Vxx'].mean()*3600.)
m[pm].hinge_piv_vyy= np.array(m[pm].hinge_piv.groupby('frame')['Vyy'].mean()*3600.)
m[pm].blade_piv_vxx= np.array(m[pm].blade_piv.groupby('frame')['Vxx'].mean()*3600.)
m[pm].blade_piv_vyy= np.array(m[pm].blade_piv.groupby('frame')['Vyy'].mean()*3600.)

m[pm].blade_shear_data, m[pm].hinge_shear_data= m[pm].region_deform_tensor('blade'), m[pm].region_deform_tensor('hinge')
m[pm].dt= m[pm].blade_shear_data['timeInt_sec']/3600.
m[pm].blade_tracked_area, m[pm].hinge_tracked_area= lib.region_area(m[pm].region_cells('blade')), lib.region_area(m[pm].region_cells('hinge'))
m[pm].vkk_blade= (np.array(m[pm].blade_tracked_area[1:])-np.array(m[pm].blade_tracked_area[:-1]))/m[pm].dt/np.array(m[pm].blade_tracked_area[:-1])
m[pm].vkk_hinge= (np.array(m[pm].hinge_tracked_area[1:])-np.array(m[pm].hinge_tracked_area[:-1]))/m[pm].dt/np.array(m[pm].hinge_tracked_area[:-1])

m[pm].hinge_tri_vxx= np.array(m[pm].hinge_shear_data['nu_xx'])/m[pm].dt+0.5*m[pm].vkk_hinge
m[pm].hinge_tri_vyy= -(np.array(m[pm].hinge_shear_data['nu_xx'])/m[pm].dt-0.5*m[pm].vkk_hinge)
m[pm].blade_tri_vxx= np.array(m[pm].blade_shear_data['nu_xx'])/m[pm].dt+0.5*m[pm].vkk_blade
m[pm].blade_tri_vyy= -(np.array(m[pm].blade_shear_data['nu_xx'])/m[pm].dt-0.5*m[pm].vkk_blade)

Nframes= len(m[pm].frames)-1
plt.figure()
plt.plot(m[pm].frames[:Nframes], m[pm].blade_old_h[:Nframes], label='old height', linewidth=1.5)
plt.plot(m[pm].frames[:Nframes], 1450*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '-',label='piv cum. shear', linewidth=1.5)
plt.plot(m[pm].frames[:Nframes], 1370*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '-',label='seg cum. shear', linewidth=1.5)
plt.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[:Nframes], label='seg height', linewidth=1.5)
plt.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
plt.legend(loc='best')
plt.show()

f, ((axWT, axDC),(axALC, axDLC))= plt.subplots(2,2)
pm='WT'
Nframes= len(m[pm].frames)-1
axWT.set_title('wild type')
axWT.plot(m[pm].frames[:Nframes], m[pm].hinge_old_h[:Nframes], label='old height', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[:Nframes], label='seg height', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
pm='DC'
Nframes= len(m[pm].frames)-1
axDC.set_title('distal cut')
axDC.plot(m[pm].frames[:Nframes], m[pm].hinge_old_h[:Nframes], label='old height', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[:Nframes], label='seg height', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
pm='ALC'
Nframes= len(m[pm].frames)-1
axALC.set_title('ant link cut')
axALC.plot(m[pm].frames[:Nframes], m[pm].hinge_old_h[:Nframes], label='old height', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[:Nframes], label='seg height', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axALC.legend(loc='best', fontsize=6)
pm='DLC'
Nframes= len(m[pm].frames)-1
axDLC.set_title('dist link cut')
axDLC.plot(m[pm].frames[:Nframes], m[pm].hinge_old_h[:Nframes], label='old height', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[:Nframes], label='seg height', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].hinge_seg_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].hinge_fcy_h[0]*np.exp(np.cumsum(m[pm].hinge_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
plt.tight_layout()
plt.savefig('figures/hinge_height_plot.png', dpi=2000)
plt.show()

plt.figure()
plt.plot(m.frames, WT_hinge_fcy_h)
plt.show()
len(dt)
len(WT_blade_piv_vxx)

Q1_hinge, Q2_hinge, s_hinge= lib.region_shape_nematic(ww_hinge)
Q1_blade, Q2_blade, s_blade= lib.region_shape_nematic(ww_blade)
area_hinge= m_DlC.region_area(ww_hinge)
area_blade= m_DlC.region_area(ww_blade)
L_hinge, h_hinge= np.sqrt(area_hinge)*np.exp(0.5*Q1_hinge), np.sqrt(area_hinge)*np.exp(-0.5*Q1_hinge)
L_blade, h_blade= np.sqrt(area_blade)*np.exp(0.5*Q1_blade), np.sqrt(area_blade)*np.exp(-0.5*Q1_blade)

blade_vxx= np.array(blade_piv.groupby('frame')['Vxx'].mean()*3600.)
blade_vyy= np.array(blade_piv.groupby('frame')['Vyy'].mean()*3600.)
hinge_vxx= np.array(hinge_piv.groupby('frame')['Vxx'].mean()*3600.)
hinge_vyy= np.array(hinge_piv.groupby('frame')['Vyy'].mean()*3600.)
time= np.array(sorted(blade_piv['time_sec'].unique()))/3600.
dt= time[1:]-time[:-1]

m_DC.blade_L, m_DC.blade_h= m_DC.region_mean_length_height(m_DC.region_cells('blade'))
m_DC.hinge_L, m_DC.hinge_h= m_DC.region_mean_length_height(m_DC.region_cells('hinge'))

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

ax1.plot(time[:230], 1480*np.exp(np.cumsum(blade_vxx[:-1]*dt))[:230], label='PIV')
ax1.plot(time[:230], m_DC.blade_L[:230], label='length')
ax1.legend(loc='best')

ax2.plot(time[:230], 1370*np.exp(np.cumsum(blade_vyy[:-1]*dt))[:230], label='PIV')
ax2.plot(time[:230], m_DC.blade_h[:230], label='height')
ax2.legend(loc='best')

ax3.plot(time[:230], 670*np.exp(np.cumsum(hinge_vxx*dt))[:230], label='PIV')
ax3.plot(time[:230], m_DC.hinge_L[:230], label='length')
ax3.legend(loc='best')

ax4.plot(time[:230], 900*np.exp(np.cumsum(hinge_vyy*dt))[:230], label='PIV')
ax4.plot(time[:230], m_DC.hinge_h[:230], label='height')
ax4.legend(loc='best')

plt.show()

len(L_blade)
f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

ax1.plot(time[:230], 1580*np.exp(np.cumsum(blade_vxx[:-1]*dt))[:230], label='PIV')
ax1.plot(time[:230], L_blade[:230], label='length')
ax1.legend(loc='best')

ax2.plot(time[:230], 1525*np.exp(np.cumsum(blade_vyy[:-1]*dt))[:230], label='PIV')
ax2.plot(time[:230], h_blade[:230], label='height')
ax2.legend(loc='best')

ax3.plot(time[:230], 1550*np.exp(np.cumsum(hinge_vxx*dt))[:230], label='PIV')
ax3.plot(time[:230], L_hinge[:230], label='length')
ax3.legend(loc='best')

ax4.plot(time[:230], 1450*np.exp(np.cumsum(hinge_vyy*dt))[:230], label='PIV')
ax4.plot(time[:230], h_hinge[:230], label='height')
ax4.legend(loc='best')

plt.show()

plt.figure()
plt.plot(time[:230], 2300000*np.exp(np.cumsum(blade_vxx[:-1]*dt))[:230]*np.exp(np.cumsum(blade_vyy[:-1]*dt))[:230], label='PIV')
plt.plot(time[:230], L_blade[:230]*h_blade[:230], label='area')
plt.legend(loc='best')
plt.show()




print m.region_cells('blade').columns
elong_xx= np.array([m.region_cells('blade')[m.region_cells('blade')['frame']==f]['elong_xx'].mean() for f in m.frames])
area= np.array(m.region_cells('blade')[m.region_cells('blade')['frame']==200]['area'])
m.regions
plt.figure()
plt.scatter(np.log(area), elong_xx)
plt.show()

m.region_center()
print m.region_shape_nematic(m.region_cells('blade'))
dir(m)
print region_center(m.region_cells('blade'))
#m.load_triList()
#m.load_Ta_t()
m.load_cells()


av_x_hinge, av_y_hinge= region_center(ww_hinge)
av_x_blade, av_y_blade= region_center(ww_blade)

otp_df= pd.DataFrame(pd.Series(L_blade))
otp_df.columns= ['blade_L']
otp_df['blade_h']= pd.Series(h_blade)
otp_df['hinge_L']= pd.Series(L_hinge)
otp_df['hinge_h']= pd.Series(h_hinge)

otp_df.to_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/drosophila_data_library/height_length_data/DlC_height_length.csv')

#L_hinge_WT, h_hinge_WT= L_hinge, h_hinge
#L_blade_WT, h_blade_WT= L_blade, h_blade
#L_hinge_AC, h_hinge_AC= L_hinge, h_hinge
#L_blade_AC, h_blade_AC= L_blade, h_blade
#L_hinge_DC, h_hinge_DC= L_hinge, h_hinge
#L_blade_DC, h_blade_DC= L_blade, h_blade

f, (ax_L_blade, ax_L_hinge, ax_h_blade, ax_h_hinge)= plt.subplots(4, 1)
ax_L_blade.plot(np.arange(0,len(L_blade_WT)), L_blade_WT, label='WT')
ax_L_blade.plot(np.arange(13,13+len(L_blade_DC)), L_blade_DC, label='DC')
ax_L_blade.plot(np.arange(13*6,13*6+len(L_blade_AC)), L_blade_AC, label='AC')
ax_L_blade.legend(loc='best')
ax_L_blade.set_ylabel('blade L')
ax_L_hinge.plot(np.arange(0,len(L_hinge_WT)), L_hinge_WT, label='WT')
ax_L_hinge.plot(np.arange(13,13+len(L_hinge_DC)), L_hinge_DC, label='DC')
ax_L_hinge.plot(np.arange(13*6,13*6+len(L_hinge_AC)), L_hinge_AC, label='AC')
ax_L_hinge.set_ylabel('hinge L')
ax_h_blade.plot(np.arange(0,len(h_blade_WT)), h_blade_WT, label='WT')
ax_h_blade.plot(np.arange(13,13+len(h_blade_DC)), h_blade_DC, label='DC')
ax_h_blade.plot(np.arange(13*6,13*6+len(h_blade_AC)), h_blade_AC, label='AC')
ax_h_blade.set_ylabel('blade h')
ax_h_hinge.plot(np.arange(0,len(h_hinge_WT)), h_hinge_WT, label='WT')
ax_h_hinge.plot(np.arange(13,13+len(h_hinge_DC)), h_hinge_DC, label='DC')
ax_h_hinge.plot(np.arange(13*6,13*6+len(h_hinge_AC)), h_hinge_AC, label='AC')
ax_h_hinge.set_ylabel('hinge h')
ax_h_hinge.set_xlabel('frames since 16h APF')
f.tight_layout()
plt.savefig('height_length_data/WT_DC_AC_compare.png', dpi=1000)
plt.show()


f, (ax_L, ax_area_blade, ax_area_hinge)= plt.subplots(3, 1)
ax_L.plot(np.arange(0,len(L_blade_WT)), L_blade_WT+L_hinge_WT, label='WT')
ax_L.plot(np.arange(13,13+len(L_blade_DC)), L_blade_DC+L_hinge_DC, label='DC')
ax_L.plot(np.arange(13*6,13*6+len(L_blade_AC)), L_blade_AC+L_hinge_AC, label='AC')
ax_L.set_ylabel('total L')
ax_area_blade.plot(np.arange(0,len(h_blade_WT)), L_blade_WT*h_blade_WT, label='WT')
ax_area_blade.plot(np.arange(13,13+len(h_blade_DC)), L_blade_DC*h_blade_DC, label='DC')
ax_area_blade.plot(np.arange(13*6,13*6+len(h_blade_AC)), L_blade_AC*h_blade_AC, label='AC')
ax_area_blade.set_ylabel('area blade')
ax_area_hinge.plot(np.arange(0,len(h_hinge_WT)), L_hinge_WT*h_hinge_WT, label='WT')
ax_area_hinge.plot(np.arange(13,13+len(h_hinge_DC)), L_hinge_DC*h_hinge_DC, label='DC')
ax_area_hinge.plot(np.arange(13*6,13*6+len(h_hinge_AC)), L_hinge_AC*h_hinge_AC, label='AC')
ax_area_hinge.set_ylabel('area hinge')
ax_area_hinge.legend(loc='best')
ax_area_hinge.set_xlabel('frames since 16h APF')
f.tight_layout()
plt.savefig('height_length_data/WT_DC_AC_totalL_area.png', dpi=1000)
plt.show()

otp_df= pd.DataFrame(pd.Series(L_blade_DC))
otp_df.columns= ['blade_L']
otp_df['blade_h']= pd.Series(h_blade_DC)
otp_df['hinge_L']= pd.Series(L_hinge_DC)
otp_df['hinge_h']= pd.Series(h_hinge_DC)

otp_df.to_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/drosophila_data_library/height_length_data/DC_height_length.csv')


ww.head()
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

blade_av_x, blade_av_y= av_x_blade, av_y_blade
hinge_av_x, hinge_av_y= av_x_hinge, av_y_blade
for frame in m.frames:
    print frame
    im_path= DB_path+'/'+m.name+'/image_data/mutant/tag/segmentationData/frame'+fill_zeros(str(frame),4)+'/original_trafo.png'
    im= plt.imread(im_path)
    y_size, x_size= im.shape[:-1]
    b_rect= patches.Rectangle((blade_av_x[frame]-0.5*L_blade[frame],blade_av_y[frame]-0.5*h_blade[frame]),L_blade[frame],h_blade[frame],fill=None, edgecolor='white', linewidth=2)
    h_rect= patches.Rectangle((hinge_av_x[frame]-0.5*L_hinge[frame],hinge_av_y[frame]-0.5*h_hinge[frame]),L_hinge[frame],h_hinge[frame],fill=None, edgecolor='green', linewidth=2)
    plt.imshow(im)
    frame_blade_x= blade_cells[blade_cells['frame']==frame]['center_x']
    frame_blade_y= blade_cells[blade_cells['frame']==frame]['center_y']
    frame_hinge_x= hinge_cells[hinge_cells['frame']==frame]['center_x']
    frame_hinge_y= hinge_cells[hinge_cells['frame']==frame]['center_y']
    plt.scatter(frame_blade_x, frame_blade_y, c='purple', s=15, linewidths=0)
    plt.scatter(frame_hinge_x, frame_hinge_y, c='yellow', s=15, linewidths=0)
    plt.gca().add_patch(b_rect)
    plt.gca().add_patch(h_rect)
    plt.axis('off')
    plt.savefig('figures/blade_size_frame'+fill_zeros(str(frame),4)+'.png', bbox_inches='tight')
    plt.close()

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