# http://justingardner.net/doku.php/mrtools/atlas?s[]=benson
# def makeBenson(freesurferSub)
import neuropythy as ny
# freesurferSub = 's0059_7T'
# freesurferSub = 's0033'
# freesurferSub = 's0047'
# freesurferSub = 's0063'
# freesurferSub = 's0059_7T'
# freesurferSub = 's0062'
#freesurferSub = 's0063'
#freesurferSub = 's0077'
#freesurferSub = 's0029'
#freesurferSub = 's0040'
#freesurferSub = 's0049'
#freesurferSub = 's0050'
#freesurferSub = 's0051'
#freesurferSub = 's0055_7T'
#freesurferSub = 's0057_7T'
#freesurferSub = 's0059_7T'
#freesurferSub = 's0060'
#freesurferSub = 's0073'
# freesurferSub = 's0076'
#freesurferSub = 's0086_7T'
#freesurferSub = 's0087_7T'
#freesurferSub = 's0092_3T'
#freesurferSub = 's0088_7T'
#freesurferSub = 's0093_7T'
#freesurferSub = 's0079'
#freesurferSub = 's0080_7T'
#freesurferSub = 's0088_7T'
#freesurferSub = 's0093_7T'
#freesurferSub = 's0074_7T'
#freesurferSub = 's0100_7T'
freesurferSub = 's0077_7T'

sub = ny.freesurfer_subject(freesurferSub)
(lh_retino, rh_retino) = ny.vision.predict_retinotopy(sub)
import scipy.io as sio
sio.savemat('lh_retino.mat',lh_retino)
sio.savemat('rh_retino.mat',rh_retino)
mdl_lh = ny.vision.retinotopy_model('benson17','lh')
with open('lh_labels.txt','w') as f:
    for item in mdl_lh.area_id_to_name:
	   f.write("%s\n" % mdl_lh.area_id_to_name[item])
mdl_rh = ny.vision.retinotopy_model('benson17','rh')
with open('rh_labels.txt','w') as f:
 	for item in mdl_rh.area_id_to_name:          
 	 	f.write("%s\n" % mdl_rh.area_id_to_name[item])