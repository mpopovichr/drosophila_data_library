f, ((axWT, axDC),(axALC, axDLC))= plt.subplots(2,2)
pm='WT'
Nframes= len(m[pm].frames)-1
axWT.set_title('wild type')
axWT.plot(m[pm].frames[:Nframes], m[pm].blade_old_h[:Nframes], label='old height', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[:Nframes], label='seg height', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axWT.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
pm='DC'
Nframes= len(m[pm].frames)-1
axDC.set_title('distal cut')
axDC.plot(m[pm].frames[:Nframes], m[pm].blade_old_h[:Nframes], label='old height', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[:Nframes], label='seg height', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
pm='ALC'
Nframes= len(m[pm].frames)-1
axALC.set_title('ant link cut')
axALC.plot(m[pm].frames[:Nframes], m[pm].blade_old_h[:Nframes], label='old height', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[:Nframes], label='seg height', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axALC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axALC.legend(loc='best', fontsize=6)
pm='DLC'
Nframes= len(m[pm].frames)-1
axDLC.set_title('dist link cut')
axDLC.plot(m[pm].frames[:Nframes], m[pm].blade_old_h[:Nframes], label='old height', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[:Nframes], label='seg height', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[:Nframes], label='fcy height', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].blade_seg_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_piv_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='piv cum. shear', linewidth=1.5)
axDLC.plot(m[pm].frames[:Nframes], m[pm].blade_fcy_h[0]*np.exp(np.cumsum(m[pm].blade_tri_vyy[:Nframes]*m[pm].dt[:Nframes])), '--',label='seg cum. shear', linewidth=1.5)
plt.tight_layout()
plt.savefig('figures/blade_height_plot.png', dpi=2000)
plt.show()





