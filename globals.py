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
import time

import lib


name='WT_25deg_111102'
m= lib.Movie(name)

hinge, blade= m.fancy_hinge_blade()
print hinge.columns
hinge_avg_elong_xx= []
for frame in hinge['frame'].unique():
    elong= np.array(hinge[hinge['frame']==frame]['elong_xx'])
    area= np.array(hinge[hinge['frame']==frame]['area'])
    hinge_avg_elong_xx.append(np.sum(elong*area)/np.sum(area))
    

plt.figure()
plt.plot(hinge['frame'].unique(), hinge_avg_elong_xx)
plt.show()

r= m.region_deform_tensor('blade')

r.columns

plt.figure()
plt.plot(lib.smooth_data(m.time), lib.smooth_data(r['ShearT1_xx']))
plt.plot(lib.smooth_data(m.time), lib.smooth_data(r['ShearT1_xy']))
plt.grid()
plt.show()


##### T1 AGE AND TRIANGLE SELECTION ###########################

m= lib.Movie('WT_25deg_111102')
rdt= m.region_deform_tensor('blade')
rdt.columns
time= np.array(rdt['time_sec'])/3600.+ 15.

print test
print test_res
test= m.triTracked[:1000]
tri_hash_list= np.array(sorted(test['tri_hash'].unique()))
test_res= test[['tri_hash', 'frame']].groupby(by='tri_hash').transform(lambda g: tri_track_last_occ(g, np.array(g.shift(-1))-np.array(g)))
test_g= test[test['tri_hash']==tri_hash_list[1]]['frame']
tri_track_last_occ(test_g, np.array(test_g.shift(-1))-np.array(test_g))
test_res= test_res.reset_index()
test['time_dis']= test_res['frame']
test['dis_type']= (
    test[['tri_hash','type']].groupby(by='tri_hash').
    transform(lambda g: g.fillna(method='bfill')).
    fillna('missing')
                   )

tri_track_last_occ(np.array(g.shift(-2))- np.array(g)) - g

m.load_triTracked()
print m.triTracked.columns
m.load_triCategories()

m.triTracked= pd.merge(m.triTracked, m.triCategories, on='tri_id',
                       how='left')
m.triTracked['tri_hash_a']= (
    m.triTracked['cell_a']*10**10
    + m.triTracked['cell_b']*10**5
    + m.triTracked['cell_c']
    )
m.triTracked['tri_hash_b']= (
    m.triTracked['cell_b']*10**10
    + m.triTracked['cell_c']*10**5
    + m.triTracked['cell_a']
    )
m.triTracked['tri_hash_c']= (
    m.triTracked['cell_c']*10**10
    + m.triTracked['cell_a']*10**5
    + m.triTracked['cell_b']
    )
m.triTracked['tri_hash']= (
    m.triTracked[['tri_hash_a',
                  'tri_hash_b',
                  'tri_hash_c']].min(axis=1)
    )
m.triTracked= m.triTracked.sort(['tri_hash', 'frame'])
m.triTracked= m.triTracked.reset_index()
m.triTracked= m.triTracked[['tri_id', 'frame', 'type', 'tri_hash']]
m.triTracked['age']= (
    m.triTracked[['tri_hash', 'frame']]
    .groupby(by='tri_hash')
    .transform(lambda g:
               g - (g.min()
                    + m.tri_track_first_occ(np.array(g)
                                            -np.array(g.shift()))))
    )
m.triTracked['time_dis']= (
    m.triTracked[['tri_hash', 'frame']]
    .groupby(by='tri_hash')
    .transform(lambda g:
               m.tri_track_last_occ(g, np.array(g.shift(-1))
                                  -np.array(g)))
    )

m.triTracked['occ_type']= m.triTracked[['tri_hash', 'type']].groupby(by='tri_hash').transform(lambda g: g.fillna(method='ffill')).fillna('missing')
m.triTracked['dis_type']= m.triTracked[['tri_hash', 'type']].groupby(by='tri_hash').transform(lambda g: g.fillna(method='bfill')).fillna('missing')

m.triTracked.to_csv(lib.DB_path+m.name+'/shear_contrib/triTracked_marko.csv')


test= pp.read_csv('/home/mpopovic/Downloads/tri_age_marko.csv')
late_test=test[test['frame']>50]
test_div= test[test['type']=='cdGain']
test_t1= test[test['type']=='t1_gain']
(1.*len(test_div)+len(test_t1))/len(test)
late_test_div= late_test[late_test['type']=='cdGain']
late_test_t1= late_test[late_test['type']=='t1_gain']
(1.*len(late_test_div)+len(late_test_t1))/len(late_test)


ne_stacked= (test[['tri_id','type',  'age']] != test_mod[['tri_id','type', 'age']]).stack()
changed= ne_stacked[ne_stacked]
changed

#ro.r('load("~/Downloads/111102__triTracked.RData")')
#ro.r('write.csv(triTracked, "~/Downloads/111102__triTracked.csv")')
triTracked_topochanges= pp.read_csv('/data/biophys/etournay/DB/WT_25deg_111102/topochanges/triTracked.csv')
triTracked_holger= pp.read_csv('/home/mpopovic/Downloads/111102__triTracked.csv')
triTracked.sort('tri_id')[['tri_id', 'frame', 'cell_a', 'cell_b', 'cell_c']].head()
triTracked_holger.sort('tri_id')[['tri_id', 'frame', 'cell_a', 'cell_b', 'cell_c']].head()
triTracked_topochanges.sort('tri_id')[['tri_id', 'frame', 'cell_a', 'cell_b', 'cell_c']].head()

test_holger= pd.merge(triTracked_holger, triCategories, on='tri_id', how='left')
test_holger= test_holger.sort(['tri_hash', 'frame'])
test_holger= test_holger.reset_index()
test_holger= test_holger[['tri_id', 'frame', 'type', 'tri_hash']]


grouped_frame_holger= test_holger[['tri_hash', 'frame']].groupby(by='tri_hash')
start= time.time()
test_age_holger= grouped_frame_holger.transform(lambda g: g - (g.min() + process_delta(np.array(g)-np.array(g.shift()))))
end= time.time()
print('Age calculation took: '+str(end-start)+' seconds')
grouped_type_holger=  test_holger[['tri_hash', 'type']].groupby(by='tri_hash')
start= time.time()
test_type_holger= grouped_type_holger.transform(lambda g: g.fillna(method='ffill'))
test_type_holger= test_type_holger.fillna('missing')
end= time.time()
print('Type filling took: '+str(end-start)+' seconds')

test_holger['age']= test_age_holger
test_holger['type']= test_type_holger

comparison= pd.merge(test, test_holger, on='tri_id', suffixes=('', '_h'))

comparison[['tri_id','frame','frame_h', 'type', 'age','min_tri_hash',  'type_h', 'age_h', 'tri_hash']].head()

first_comparison= pd.merge(triTracked[['tri_id', 'frame']], triTracked_holger[['tri_id', 'frame']], on='tri_id')
first_comparison.head()

test['min_tri_hash']= test[['tri_hash_a', 'tri_hash_b', 'tri_hash_c']].min(axis=1)
test= test[['tri_id', 'frame', 'type', 'min_tri_hash', 'tri_hash']]
test= test.sort(['min_tri_hash', 'frame'])
type(grouped)

ro.r('load("'+lib.DB_path+m.name+'/shear_contrib/triList.RData")')
ro.r('write.csv(triList, "'+lib.DB_path+m.name+'/shear_contrib/triList.csv")')
testtriList= pd.read_csv(lib.DB_path+m.name+'/shear_contrib/triList.csv')
testtriList.columns
second_comparison= pd.merge(triTracked[['tri_id','frame']], testtriList[['tri_id','frame']],on='tri_id')
second_comparison.head(20)



T1= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotT1_done.csv')
noT1= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotT1_no.csv')
CD= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotCD_done.csv')
noCD= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotCD_no.csv')
total= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotBlade.csv')
sample= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotSample.csv')
T1_dis= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotT1_dis.csv')
T1_nodis= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotT1_nodis.csv')
CD_dis= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotCD_dis.csv')
CD_nodis= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotCD_nodis.csv')

Ktotal= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotKnown.csv')
KT1= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotKnownT1_done.csv')
KnoT1= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotKnownT1_no.csv')
KT1_dis= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotKnownT1_dis.csv')
KT1_nodis= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotKnownT1_nodis.csv')

T1= T1.sort('frame')
noT1= noT1.sort('frame')
T1_dis= T1_dis.sort('frame')
T1_nodis= T1_nodis.sort('frame')
Ktotal= Ktotal.sort('frame')
KT1= KT1.sort('frame')
KnoT1= KnoT1.sort('frame')
KT1_dis= KT1_dis.sort('frame')
KT1_nodis= KT1_nodis.sort('frame')
CD= CD.sort('frame')
noCD= noCD.sort('frame')
CD_dis= CD_dis.sort('frame')
CD_nodis= CD_nodis.sort('frame')
total= total.sort('frame')

frames_dis= lib.smooth_data(sorted(np.array(T1_dis['frame'])))
frames= lib.smooth_data(sorted(np.array(T1['frame'])))
T1_dis_QxyThetaU= lib.smooth_data(np.array(T1_dis['QxyThetaU']))
T1_dis_area= lib.smooth_data(np.array(T1_dis['tri_area.i1']))
T1_nodis_QxyThetaU= lib.smooth_data(np.array(T1_nodis['QxyThetaU']))
T1_nodis_area= lib.smooth_data(np.array(T1_nodis['tri_area.i1']))
T1_QxyThetaU= lib.smooth_data(np.array(T1['QxyThetaU']))
T1_area= lib.smooth_data(np.array(T1['tri_area.i1']))
noT1_QxyThetaU= lib.smooth_data(np.array(noT1['QxyThetaU']))
noT1_area= lib.smooth_data(np.array(noT1['tri_area.i1']))

CD_QxyThetaU= lib.smooth_data((np.array(CD['QxyThetaU'])))
CD_area= lib.smooth_data(np.array(CD['tri_area.i1']))
noCD_QxyThetaU= lib.smooth_data((np.array(noCD['QxyThetaU'])))
noCD_area= lib.smooth_data(np.array(noCD['tri_area.i1']))
CD_dis_QxyThetaU= lib.smooth_data(np.array(CD_dis['QxyThetaU']))
CD_dis_area= lib.smooth_data(np.array(CD_dis['tri_area.i1']))
CD_nodis_QxyThetaU= lib.smooth_data(np.array(CD_nodis['QxyThetaU']))
CD_nodis_area= lib.smooth_data(np.array(CD_nodis['tri_area.i1']))
total_QxyThetaU= lib.smooth_data(np.array(total['QxyThetaU']))
total_area= lib.smooth_data(np.array(total['tri_area.i1']))

Ktotal_QxyThetaU= lib.smooth_data(np.array(Ktotal['QxyThetaU']))
Ktotal_area= lib.smooth_data(np.array(Ktotal['tri_area.i1']))
KT1_dis_QxyThetaU= lib.smooth_data(np.array(KT1_dis['QxyThetaU']))
KT1_dis_area= lib.smooth_data(np.array(KT1_dis['tri_area.i1']))
KT1_nodis_QxyThetaU= lib.smooth_data(np.array(KT1_nodis['QxyThetaU']))
KT1_nodis_area= lib.smooth_data(np.array(KT1_nodis['tri_area.i1']))
KT1_QxyThetaU= lib.smooth_data(np.array(KT1['QxyThetaU']))
KT1_area= lib.smooth_data(np.array(KT1['tri_area.i1']))
KnoT1_QxyThetaU= lib.smooth_data(np.array(KnoT1['QxyThetaU']))
KnoT1_area= lib.smooth_data(np.array(KnoT1['tri_area.i1']))
#time= np.array(T1_dis['time_sec'])/3600.
dt= time[1:]-time[:-1]
dt= lib.smooth_data(dt)
smooth_time= lib.smooth_data(time[:-1])

len(CD_QxyThetaU)
len(total_area)

plt.figure()
plt.plot(smooth_time, CD_QxyThetaU*CD_area/total_area[1:]/dt, label='recent CD')
plt.plot(smooth_time, total_QxyThetaU[1:]*CD_area/total_area[1:]/dt, label='reference')
plt.plot(smooth_time, noCD_QxyThetaU[1:]*noCD_area[1:]/total_area[1:]/dt, label='other')
plt.plot(smooth_time, total_QxyThetaU[1:]/dt, label='total')
plt.plot(smooth_time,CD_QxyThetaU*CD_area/total_area[1:]/dt+noCD_QxyThetaU[1:]*noCD_area[1:]/total_area[1:]/dt, label='check' )
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'$\langle \omega Q_{xy}\rangle$', fontsize=20)
plt.grid()
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('figures/recent_CD_correlation_contribution.png', dpi=300)
#plt.close()
plt.show()


plt.figure()
plt.plot(smooth_time, T1_QxyThetaU*T1_area/total_area[:-1]/dt, label='recent T1')
plt.plot(smooth_time, total_QxyThetaU[:-1]*T1_area/total_area[:-1]/dt, label='reference')
plt.plot(smooth_time, noT1_QxyThetaU[:-1]*noT1_area[:-1]/total_area[:-1]/dt, label='other')
plt.plot(smooth_time, total_QxyThetaU[:-1]/dt, label='total')
plt.plot(smooth_time,T1_QxyThetaU*T1_area/total_area[:-1]/dt+noT1_QxyThetaU[:-1]*noT1_area[:-1]/total_area[:-1]/dt, label='check' )
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'$\langle \omega Q_{xy}\rangle$', fontsize=20)
plt.grid()
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('figures/recent_T1_correlation_contribution.png', dpi=300)
#plt.close()
plt.show()

plt.figure()
plt.plot(smooth_time, T1_dis_QxyThetaU[:-1]*T1_dis_area[:-1]/total_area[:-1]/dt, label='soon to T1')
plt.plot(smooth_time, total_QxyThetaU[:-1]*T1_dis_area[:-1]/total_area[:-1]/dt, label='reference')
#plt.plot(smooth_time, sample_QxyThetaU*sample_area/total_area/dt, label='sample')
plt.plot(smooth_time, T1_nodis_QxyThetaU[:-1]*T1_nodis_area[:-1]/total_area[:-1]/dt, label='other')
plt.plot(smooth_time, total_QxyThetaU[:-1]/dt, label='total')
plt.plot(smooth_time,T1_dis_QxyThetaU[:-1]*T1_dis_area[:-1]/total_area[:-1]/dt+T1_nodis_QxyThetaU[:-1]*T1_nodis_area[:-1]/total_area[:-1]/dt, label='check' )
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'$\langle \omega Q_{xy}\rangle$', fontsize=20)
plt.grid()
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('figures/soon_to_T1_correlation_contribution.png', dpi=300)
#plt.close()
plt.show()


plt.figure()
plt.plot(smooth_time, CD_dis_QxyThetaU[:-1]*CD_dis_area[:-1]/total_area[:-1]/dt, label='soon to CD')
plt.plot(smooth_time, total_QxyThetaU[:-1]*CD_dis_area[:-1]/total_area[:-1]/dt, label='reference')
#plt.plot(smooth_time, sample_QxyThetaU*sample_area/total_area/dt, label='sample')
plt.plot(smooth_time, CD_nodis_QxyThetaU[:-1]*CD_nodis_area[:-1]/total_area[:-1]/dt, label='other')
plt.plot(smooth_time, total_QxyThetaU[:-1]/dt, label='total')
plt.plot(smooth_time,CD_dis_QxyThetaU[:-1]*CD_dis_area[:-1]/total_area[:-1]/dt+CD_nodis_QxyThetaU[:-1]*CD_nodis_area[:-1]/total_area[:-1]/dt, label='check' )
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'$\langle \omega Q_{xy}\rangle$', fontsize=20)
plt.grid()
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('figures/soon_to_CD_correlation_contribution.png', dpi=300)
#plt.close()
plt.show()


plt.figure()
plt.plot(smooth_time, KT1_QxyThetaU*KT1_area/Ktotal_area[1:]/dt, label='recent T1')
plt.plot(smooth_time, Ktotal_QxyThetaU[1:]*KT1_area/Ktotal_area[1:]/dt, label='reference')
plt.plot(smooth_time, KnoT1_QxyThetaU[1:]*KnoT1_area[1:]/Ktotal_area[1:]/dt, label='other')
plt.plot(smooth_time, Ktotal_QxyThetaU[1:]/dt, label='total')
plt.plot(smooth_time, KT1_QxyThetaU*KT1_area/Ktotal_area[1:]/dt+KnoT1_QxyThetaU[1:]*KnoT1_area[1:]/Ktotal_area[1:]/dt, label='check' )
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'$\langle \omega Q_{xy}\rangle$', fontsize=20)
plt.grid()
plt.legend(loc=4)
plt.tight_layout()
#plt.savefig('figures/recent_T1_correlation_contribution.png', dpi=300)
#plt.close()
plt.show()

plt.figure()
plt.plot(smooth_time, KT1_dis_QxyThetaU[:-1]*KT1_dis_area[:-1]/Ktotal_area[:-1]/dt, label='soon to T1')
plt.plot(smooth_time, Ktotal_QxyThetaU[:-1]*KT1_dis_area[:-1]/Ktotal_area[:-1]/dt, label='reference')
#plt.plot(smooth_time, sample_QxyThetaU*sample_area/total_area/dt, label='sample')
plt.plot(smooth_time, KT1_nodis_QxyThetaU[:-1]*KT1_nodis_area[:-1]/Ktotal_area[:-1]/dt, label='other')
plt.plot(smooth_time, Ktotal_QxyThetaU[:-1]/dt, label='total')
plt.plot(smooth_time, KT1_dis_QxyThetaU[:-1]*KT1_dis_area[:-1]/Ktotal_area[:-1]/dt+KT1_nodis_QxyThetaU[:-1]*KT1_nodis_area[:-1]/Ktotal_area[:-1]/dt, label='check' )
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'$\langle \omega Q_{xy}\rangle$', fontsize=20)
plt.grid()
plt.legend(loc=4)
plt.tight_layout()
#plt.savefig('figures/soon_to_T1_correlation_contribution.png', dpi=300)
#plt.close()
plt.show()




plt.figure()
plt.plot(smooth_time, T1_area, label='recent T1')
#plt.plot(frames, sample_area/total_area, label='sample')
plt.plot(smooth_time, noT1_area[1:], label='other')
plt.plot(smooth_time, total_area[1:], label= 'total')
plt.plot(smooth_time, T1_area+ noT1_area[1:], label='check')
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'area[px]', fontsize=20)
plt.grid()
plt.legend(loc=7)
plt.tight_layout()
#plt.savefig('figures/soon_to_T1_area.png', dpi=300)
plt.show()

plt.figure()
plt.plot(smooth_time, T1_area[:-1], label='soon to T1')
#plt.plot(frames, sample_area/total_area, label='sample')
plt.plot(smooth_time, noT1_area[:-1], label='other')
plt.plot(smooth_time, total_area[:-1], label= 'total')
plt.plot(smooth_time, T1_area[:-1]+ noT1_area[:-1], label='check')
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'area[px]', fontsize=20)
plt.grid()
plt.legend(loc=7)
plt.tight_layout()
plt.savefig('figures/soon_to_T1_area.png', dpi=300)
plt.show()

plt.figure()
smoothPlot(frames, np.array(T1_QxyThetaU), label='T1')
smoothPlot(frames, np.array(noT1_QxyThetaU[:-1]), label='old')
smoothPlot(frames, total_QxyThetaU[:-1], label='total')
plt.grid()
plt.legend(loc=4)
plt.show()

plt.figure()
plt.plot(frames, T1_area, label='T1')
plt.plot(frames, noT1_area[:-1], label= 'old')
plt.plot(frames, total_area[:-1], label='total')
plt.grid()
plt.legend(loc=4)
plt.show()


CD= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotCD_done.csv')
noCD= pp.read_csv('/home/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/Correlations/avgDeltaQtotCD_no.csv')


CD= CD.sort('frame')
noCD= noCD.sort('frame')
total= total.sort('frame')

frames= smooth_data(sorted(np.array(CD['frame'])))
CD_QxyThetaU= smooth_data(np.array(CD['QxyThetaU']))
CD_ThetaU= smooth_data(np.array(CD['ThetaU']))
CD_area= smooth_data(np.array(CD['tri_area.i1']))
CD_Q_xy_i1= smooth_data(np.array(CD['Q_xy.i1']))
CD_Q_xy_ti2= smooth_data(np.array(CD['Q_xy.ti2']))
sample_QxyThetaU= smooth_data(np.array(sample['QxyThetaU'])[:-1])
sample_area= smooth_data(np.array(sample['tri_area.i1'])[:-1])
noCD_QxyThetaU= smooth_data(np.array(noCD['QxyThetaU'])[:-1])
noCD_area= smooth_data(np.array(noCD['tri_area.i1'])[:-1])
total_QxyThetaU= smooth_data(np.array(total['QxyThetaU'])[:-1])
total_area= smooth_data(np.array(total['tri_area.i1'])[:-1])
dt= time[1:]-time[:-1]
dt= smooth_data(dt)
smooth_time= smooth_data(time[:-1])


plt.figure()
plt.plot(smooth_time, CD_QxyThetaU*CD_area/total_area/dt, label='CD')
#plt.plot(smooth_time, sample_QxyThetaU*sample_area/total_area/dt, label='sample')
plt.plot(smooth_time, total_QxyThetaU*CD_area/total_area/dt, label='reference')
plt.plot(smooth_time, noCD_QxyThetaU*noCD_area/total_area/dt, label='old')
plt.plot(smooth_time, total_QxyThetaU/dt, label='total')
plt.plot(smooth_time,CD_QxyThetaU*CD_area/total_area/dt+noCD_QxyThetaU*noCD_area/total_area/dt, label='check' )
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'$\langle \delta \psi Q_{xy}\rangle$', fontsize=20)
plt.grid()
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('figures/CD_correlation_contribution.png', dpi=300)
#plt.close()
plt.show()

plt.figure()
plt.plot(smooth_time, CD_area, label='CD')
#plt.plot(frames, sample_area/total_area, label='sample')
plt.plot(smooth_time, noCD_area, label='old')
plt.xlabel(r'timeAPF[h]',fontsize=20)
plt.ylabel(r'area[px]', fontsize=20)
plt.grid()
plt.legend(loc=7)
plt.tight_layout()
plt.savefig('figures/CD_area.png', dpi=300)
plt.show()

##### AREA ELONGATION CORRELATION PLOT AND MOVIE ######################



#######################################################################

wdf= pp.read_csv(lib.DB_path+'PupalWingMovies.csv', sep='\t')

name= 'MTdp_25deg_140222'
np.isnan(np.array(wdf[wdf['nice_name']==name]['time_shift_sec'])[0])


## SHEAR CORRECTED DIMENSIONS##############

m= {}
m['WT']= lib.Movie('WT_25deg_111103')
m['ALC']= lib.Movie('WT_antLinkCut-25deg_131227')
m['DLC']= lib.Movie('WT_distLinkCut-25deg_131226')
m['dp']= lib.Movie('MTdp_25deg_140222')

m['dp'].time
m['dp'].dt

pm= 'dp'
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

print('Done!')

Nframes= len(m[pm].frames)-1
m[pm].hinge_shear_h= np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].hinge_shear_L= np.exp(np.cumsum(m[pm].hinge_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_h= np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_L= np.exp(np.cumsum(m[pm].blade_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))


pm= 'dp'
fitting_shift= 100
def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_h[fitting_shift:], m[pm].hinge_fcy_h[fitting_shift:-1])
m[pm].beta_hinge_h= optimize.anneal(scaling_func, m[pm].hinge_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_L[fitting_shift:], m[pm].hinge_fcy_L[fitting_shift:-1])
m[pm].beta_hinge_L= optimize.anneal(scaling_func, m[pm].hinge_fcy_L[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_h[fitting_shift:], m[pm].blade_fcy_h[fitting_shift:-1])
m[pm].beta_blade_h= optimize.anneal(scaling_func, m[pm].blade_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_L[fitting_shift:], m[pm].blade_fcy_L[fitting_shift:-1])
m[pm].beta_blade_L= optimize.anneal(scaling_func, m[pm].blade_fcy_L[0])[0]


f, ((axBH, axHH),(axBL, axHL))= plt.subplots(2,2)

axBH.set_ylabel(r'blade height[px]')
axBH.set_xlabel(r'timeAPF[h]')
axBH.plot(m['dp'].time[0]+np.cumsum(m['dp'].dt), m['dp'].beta_blade_h*m['dp'].blade_shear_h, label='dp', color='red')
axBH.plot(m['dp'].time, m['dp'].blade_fcy_h[:-1], color='blue')
axBH.grid()

axHH.set_ylabel(r'hinge height[px]')
axHH.set_xlabel(r'timeAPF[h]')
axHH.plot(m['dp'].time[0]+np.cumsum(m['dp'].dt), m['dp'].beta_hinge_h*m['dp'].hinge_shear_h, label='dp', color='red')
axHH.plot(m['dp'].time, m['dp'].hinge_fcy_h[:-1], color='blue')
axHH.legend(loc='best')
axHH.grid()

axBL.set_ylabel(r'blade length[px]')
axBL.set_xlabel(r'timeAPF[h]')
axBL.plot(m['dp'].time[0]+np.cumsum(m['dp'].dt), m['dp'].beta_blade_L*m['dp'].blade_shear_L, label='dp', color='red')
axBL.plot(m['dp'].time, m['dp'].blade_fcy_L[:-1], color='blue')
axBL.grid()

axHL.set_ylabel(r'hinge_length[px]')
axHL.set_xlabel(r'timeAPF[h]')
axHL.plot(m['dp'].time[0]+np.cumsum(m['dp'].dt), m['dp'].beta_hinge_L*m['dp'].hinge_shear_L, label='dp', color='red')
axHL.plot(m['dp'].time, m['dp'].hinge_fcy_L[:-1], color='blue')
axHL.grid()

f.tight_layout()
#plt.savefig('figures/cumulative_piv_shear_fitted_after_22.png', dpi=1000)
plt.show()

df_otp= pd.DataFrame()
df_otp['dp_blade_shear_h']= m['dp'].beta_blade_h*m['dp'].blade_shear_h
df_otp['dp_hinge_shear_h']= m['dp'].beta_hinge_h*m['dp'].hinge_shear_h
df_otp['dp_blade_shear_L']= m['dp'].beta_blade_L*m['dp'].blade_shear_L
df_otp['dp_hinge_shear_L']= m['dp'].beta_hinge_L*m['dp'].hinge_shear_L

df_otp.to_csv('height_length_data/shear_corrected_dimensions_dumpy.csv')


pm= 'DLC'

m[pm].hinge_fcy, m[pm].blade_fcy= m[pm].fancy_hinge_blade()
m[pm].hinge_old, m[pm].blade_old= m[pm].whole_hinge_blade()

m[pm].hinge_fcy_L, m[pm].hinge_fcy_h= lib.region_mean_length_height(m[pm].hinge_fcy)
m[pm].blade_fcy_L, m[pm].blade_fcy_h= lib.region_mean_length_height(m[pm].blade_fcy)
m[pm].hinge_old_L, m[pm].hinge_old_h= lib.region_mean_length_height(m[pm].hinge_old)
m[pm].blade_old_L, m[pm].blade_old_h= lib.region_mean_length_height(m[pm].blade_old)
m[pm].hinge_seg_L, m[pm].hinge_seg_h= lib.region_mean_length_height(m[pm].region_cells('hinge'))
m[pm].blade_seg_L, m[pm].blade_seg_h= lib.region_mean_length_height(m[pm].region_cells('blade'))

print('Done!')

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

print('Done!')

#pm='WT'
Nframes= len(m[pm].frames)-1
m[pm].hinge_shear_h= np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].hinge_shear_L= np.exp(np.cumsum(m[pm].hinge_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_h= np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_L= np.exp(np.cumsum(m[pm].blade_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))

pm= 'WT'
def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_h[78:], m[pm].hinge_fcy_h[78:-1])
m[pm].beta_hinge_h= optimize.anneal(scaling_func, m[pm].hinge_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_L[78:], m[pm].hinge_fcy_L[78:-1])
m[pm].beta_hinge_L= optimize.anneal(scaling_func, m[pm].hinge_fcy_L[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_h[78:], m[pm].blade_fcy_h[78:-1])
m[pm].beta_blade_h= optimize.anneal(scaling_func, m[pm].blade_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_L[78:], m[pm].blade_fcy_L[78:-1])
m[pm].beta_blade_L= optimize.anneal(scaling_func, m[pm].blade_fcy_L[0])[0]


pm= 'WT'
def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_h[78:], m[pm].hinge_fcy_h[78:-1])
m[pm].beta_hinge_h= optimize.anneal(scaling_func, m[pm].hinge_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_L[78:], m[pm].hinge_fcy_L[78:-1])
m[pm].beta_hinge_L= optimize.anneal(scaling_func, m[pm].hinge_fcy_L[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_h[78:], m[pm].blade_fcy_h[78:-1])
m[pm].beta_blade_h= optimize.anneal(scaling_func, m[pm].blade_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_L[78:], m[pm].blade_fcy_L[78:-1])
m[pm].beta_blade_L= optimize.anneal(scaling_func, m[pm].blade_fcy_L[0])[0]


pm= 'ALC'
Nframes= len(m[pm].frames)-1
m[pm].hinge_shear_h= np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].hinge_shear_L= np.exp(np.cumsum(m[pm].hinge_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_h= np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_L= np.exp(np.cumsum(m[pm].blade_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_h, m[pm].hinge_fcy_h[:-1])
m[pm].beta_hinge_h= optimize.anneal(scaling_func, m[pm].hinge_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_L, m[pm].hinge_fcy_L[:-1])
m[pm].beta_hinge_L= optimize.anneal(scaling_func, m[pm].hinge_fcy_L[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_h, m[pm].blade_fcy_h[:-1])
m[pm].beta_blade_h= optimize.anneal(scaling_func, m[pm].blade_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_L, m[pm].blade_fcy_L[:-1])
m[pm].beta_blade_L= optimize.anneal(scaling_func, m[pm].blade_fcy_L[0])[0]


pm= 'DLC'
Nframes= len(m[pm].frames)-1
m[pm].hinge_shear_h= np.exp(np.cumsum(m[pm].hinge_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].hinge_shear_L= np.exp(np.cumsum(m[pm].hinge_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_h= np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes]))
m[pm].blade_shear_L= np.exp(np.cumsum(m[pm].blade_piv_vxx[:Nframes]*m[pm].dt[:Nframes]))

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_h, m[pm].hinge_fcy_h[:-1])
m[pm].beta_hinge_h= optimize.anneal(scaling_func, m[pm].hinge_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].hinge_shear_L, m[pm].hinge_fcy_L[:-1])
m[pm].beta_hinge_L= optimize.anneal(scaling_func, m[pm].hinge_fcy_L[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_h, m[pm].blade_fcy_h[:-1])
m[pm].beta_blade_h= optimize.anneal(scaling_func, m[pm].blade_fcy_h[0])[0]

def scaling_func(beta):
    return lib.square_difference(beta*m[pm].blade_shear_L, m[pm].blade_fcy_L[:-1])
m[pm].beta_blade_L= optimize.anneal(scaling_func, m[pm].blade_fcy_L[0])[0]



f, ((axBH, axHH),(axBL, axHL))= plt.subplots(2,2)

axBH.set_ylabel(r'blade height[px]')
axBH.set_xlabel(r'timeAPF[h]')
axBH.plot(16.+np.cumsum(m['WT'].dt), m['WT'].beta_blade_h*m['WT'].blade_shear_h, label='WT', color='red')
axBH.plot(22.+np.cumsum(m['ALC'].dt), m['ALC'].beta_blade_h*m['ALC'].blade_shear_h, label='ALC', color='blue')
axBH.plot(22.+np.cumsum(m['DLC'].dt), m['DLC'].beta_blade_h*m['DLC'].blade_shear_h, label='DLC', color='black')
axBH.grid()

axHH.set_ylabel(r'hinge height[px]')
axHH.set_xlabel(r'timeAPF[h]')
axHH.plot(16.+np.cumsum(m['WT'].dt), m['WT'].beta_hinge_h*m['WT'].hinge_shear_h, label='WT', color='red')
axHH.plot(22.+np.cumsum(m['ALC'].dt), m['ALC'].beta_hinge_h*m['ALC'].hinge_shear_h, label='ALC', color='blue')
axHH.plot(22.+np.cumsum(m['DLC'].dt), m['DLC'].beta_hinge_h*m['DLC'].hinge_shear_h, label='DLC', color='black')
axHH.legend(loc='best')
axHH.grid()

axBL.set_ylabel(r'blade length[px]')
axBL.set_xlabel(r'timeAPF[h]')
axBL.plot(16.+np.cumsum(m['WT'].dt), m['WT'].beta_blade_L*m['WT'].blade_shear_L, label='WT', color='red')
axBL.plot(22.+np.cumsum(m['ALC'].dt), m['ALC'].beta_blade_L*m['ALC'].blade_shear_L, label='ALC', color='blue')
axBL.plot(22.+np.cumsum(m['DLC'].dt), m['DLC'].beta_blade_L*m['DLC'].blade_shear_L, label='DLC', color='black')
axBL.grid()

axHL.set_ylabel(r'hinge_length[px]')
axHL.set_xlabel(r'timeAPF[h]')
axHL.plot(16.+np.cumsum(m['WT'].dt), m['WT'].beta_hinge_L*m['WT'].hinge_shear_L, label='WT', color='red')
axHL.plot(22.+np.cumsum(m['ALC'].dt), m['ALC'].beta_hinge_L*m['ALC'].hinge_shear_L, label='ALC', color='blue')
axHL.plot(22.+np.cumsum(m['DLC'].dt), m['DLC'].beta_hinge_L*m['DLC'].hinge_shear_L, label='DLC', color='black')
axHL.grid()

f.tight_layout()
plt.savefig('figures/cumulative_piv_shear_fitted_after_22.png', dpi=1000)
plt.show()

df_otp= pd.DataFrame()
df_otp['WT_blade_shear_h']= m['WT'].beta_blade_h*m['WT'].blade_shear_h
df_otp['ALC_blade_shear_h']= m['ALC'].beta_blade_h*m['ALC'].blade_shear_h
df_otp['DLC_blade_shear_h']= m['DLC'].beta_blade_h*m['DLC'].blade_shear_h
df_otp['WT_hinge_shear_h']= m['WT'].beta_hinge_h*m['WT'].hinge_shear_h
df_otp['ALC_hinge_shear_h']= m['ALC'].beta_hinge_h*m['ALC'].hinge_shear_h
df_otp['DLC_hinge_shear_h']= m['DLC'].beta_hinge_h*m['DLC'].hinge_shear_h
df_otp['WT_blade_shear_L']= m['WT'].beta_blade_L*m['WT'].blade_shear_L
df_otp['ALC_blade_shear_L']= m['ALC'].beta_blade_L*m['ALC'].blade_shear_L
df_otp['DLC_blade_shear_L']= m['DLC'].beta_blade_L*m['DLC'].blade_shear_L
df_otp['WT_hinge_shear_L']= m['WT'].beta_hinge_L*m['WT'].hinge_shear_L
df_otp['ALC_hinge_shear_L']= m['ALC'].beta_hinge_L*m['ALC'].hinge_shear_L
df_otp['DLC_hinge_shear_L']= m['DLC'].beta_hinge_L*m['DLC'].hinge_shear_L

df_otp.to_csv('height_length_data/shear_corrected_dimensions.csv')

#### END SHEAR CORRECTED DIMENSIONS ##############################




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
plt.plot(m[pm].frames, m[pm].blade_fcy_L*m[pm].blade_fcy_h/(m[pm].blade_fcy_L[0]*m[pm].blade_fcy_h[0]), label='old area')
plt.show()

plt.plot(m[pm].frames, m[pm].blade_fcy_L*m[pm].blade_fcy_h, label='fcy area')
plt.plot(m[pm].frames, m[pm].blade_seg_L*m[pm].blade_seg_h, label='seg area')

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

m['WT'].regions


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
