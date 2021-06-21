#!/usr/bin/env python3

#=== read_data.py:
syn_band_code = "BX"
syn_suffix = ".sem.sac"
obs_band_code = "BH"
left_pad = 100
right_pad = 0
obs_preevent = 50
syn_is_grn = False
GPS_ELLPS = "WGS84" # geodesic/ECEF conversion
cmt_is_ECEF = False # True for green function in source inversion

#=== measure_misfit.py
def make_window_list_P_wave(evdp_km):
    """
    Time window for P-wave
    """
    windows = []
    if evdp_km < 150:
        windows = [
            {'id':'p,P,pP_Z', 'phase':'p,P,pP', 'component':'Z', 'time':[-50,50], 'filter':[0.01, 0.2, 2], 'dist':[0, 180], 'pre_weight':1.0},
            {'id':'p,P,pP_Z', 'phase':'p,P,pP', 'component':'R', 'time':[-50,50], 'filter':[0.01, 0.2, 2], 'dist':[0, 180], 'pre_weight':1.0},
        ]
    return windows

def make_window_list_S_wave(evdp_km):
    """Time windows for S wave
    """
    windows = []
    if evdp_km < 150:
        windows = [
            {'id':'s,S,pS,sS_Z', 'phase':'s,S,pS,sS', 'component':'Z', 'time':[-50,50], 'filter':[0.01, 0.2, 2], 'dist':[0, 180], 'pre_weight':1.0},
            {'id':'s,S,pS,sS_R', 'phase':'s,S,pS,sS', 'component':'R', 'time':[-50,50], 'filter':[0.01, 0.2, 2], 'dist':[0, 180], 'pre_weight':1.0},
            {'id':'s,S,sS_T', 'phase':'s,S,sS', 'component':'T', 'time':[-50,50], 'filter':[0.01, 0.2, 2], 'dist':[0, 180], 'pre_weight':1.0},
        ]
    return windows

def make_window_list_surface_wave(evdp_km):
    """Time windows for surface wave
    """
    windows = []
    if evdp_km < 100:
        windows = [
            {'phase':'surface', 'component':'Z', 'time':[-50,50], 'slowness':[25, 40], 'filter':[0.01, 0.02, 2], 'pre_weight':1.0},
            {'phase':'surface', 'component':'R', 'time':[-50,50], 'slowness':[25, 40], 'filter':[0.01, 0.02, 2], 'pre_weight':1.0},
            {'phase':'surface', 'component':'T', 'time':[-50,50], 'slowness':[25, 40], 'filter':[0.01, 0.02, 2], 'pre_weight':1.0},
        ]
    return windows

#misfit.py:measure_adj()
misfit_type = "cc0" # normalized zero-lag xcorrelation ceofficient: cc0, xcorr time shift: ccdt
weight_param = {
    'cc_tshift':[-10,-8, 8,10],
    'SNR':[10,15],
    'CC0':[0.5,0.7],
    #'CCmax':[0.8, 0.9],
    #'dist':[0.5, 1.0], #[unit: degree] exclude short range records (<0.5 deg)
    #'cc_tshift':[-10,-5,5,10],
}

#=== plot_misfit.py
plot_azbin_size = 10
plot_begin_time = -10
plot_end_time = 10
plot_clip_ratio = 1.5
plot_min_CC0 = None
plot_min_CCmax = None
plot_min_SNR = 5
plot_dist_lim = None #[0, 180]
plot_az0 = 0
plot_align_time = False
