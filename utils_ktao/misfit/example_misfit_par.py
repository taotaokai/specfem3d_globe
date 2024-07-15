
cmt_is_ECEF = True
GPS_ELLPS = 'WGS84'

# seismogram
obs_band_code = ['BHE', 'BHN', 'BHZ']
syn_band_code = ['MXE', 'MXN', 'MXZ']
syn_suffix = 'sem'


left_pad = 50 
right_pad = 50
obs_preevent = 100 

syn_is_grn = False # True in the green,dgreen stages of source inversion

plot_azbin_size = 10
plot_begin_time = -10
plot_end_time = 10
plot_clip_ratio = 2
plot_min_CC0 = 0.8
plot_min_CCmax = 0.9
plot_min_SNR = 10
plot_dist_lim = [0,180]
plot_az0 = 0
plot_align_time = 

misfit_type = 'CC0' # 'dt_CC'
weight_param = 

def make_window_list_P_wave(evdp_km):
    window_list = []
    if evdp_km < 50:
        window_list = [
                {'id':'p,P_Z' , 'phase':'p,P', 'component':'Z', 'time':[-50,50], 'filter':[0.01,0.1,2] , 'dist':[0, 180], 'pre_weight':1},
                {'id':'p,P_R' , 'phase':'p,P', 'component':'R', 'time':[-50,50], 'filter':[0.01,0.1,2] , 'dist':[0, 180], 'pre_weight':1},
                ]
    elif evdp_km < 100:
        window_list = [
                {'id':'p,P_Z' , 'phase':'p,P', 'component':'Z', 'time':[-50,50], 'filter':[0.01,0.1,2] , 'dist':[0, 180], 'pre_weight':1},
                {'id':'pP,sP_Z' , 'phase':'pP,sP', 'component':'Z', 'time':[-50,50], 'filter':[0.01,0.1,2] , 'dist':[0, 180], 'pre_weight':1},
                ]

    return window_list


def make_window_list_S_wave(evdp_km):

def make_window_list_surface_wave(evdp):
