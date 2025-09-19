# utils.py
#
# Utility functions for stochprop methods
#
# Philip Blom (pblom@lanl.gov)

import configparser as cnfg
import numpy as np

from importlib.util import find_spec

from scipy.interpolate import interp1d


def eig_wvfrm2json(wvfrm_file, origin_time, c0, ns_option, array_dim, sta, loc, cha):

    infrapy_config = cnfg.ConfigParser()
    infrapy_config.read(find_spec('infrapy').submodule_search_locations[0] + '/resources/default.config')
    fk_dict = dict(infrapy_config.items('FK'))
    det_dict = dict(infrapy_config.items('FD'))

    for key in fk_dict:
        try:
            fk_dict[key] = float(fk_dict[key])
        except:
            pass

    for key in det_dict:
        try:
            det_dict[key] = float(det_dict[key])
        except:
            pass

    # set waveform info
    t1 = origin_time.replace(hour = origin_time.hour - 1, minute=0, second=0, microsecond=0)
    t2 = origin_time.replace(hour = origin_time.hour + 1, minute=0, second=0, microsecond=0)
    with open(wvfrm_file, 'r') as file:
        for line in file:
            if "receiver location" in line:
                line_info = line.split(" ")
                lat, lon = line_info[-3].replace(",", ""), line_info[-2].replace(",", "")
                break
    wvfrm_dicts = []
    for k in range(array_dim):
        wvfrm_dicts = wvfrm_dicts + [{}]
        wvfrm_dicts[-1]['trace id'] = "SY." + sta + str(k + 1) + "." + loc + "." + cha
        wvfrm_dicts[-1]['strattime'] = t1
        wvfrm_dicts[-1]['endtime'] = t2
        wvfrm_dicts[-1]['latitude'] = float(lat)
        wvfrm_dicts[-1]['longitude'] = float(lon)

    # Read arrival info from header
    arrivals = []
    with open(wvfrm_file, 'r') as file:
        in_arrivals = False
        for line in file:
            if "# t [s]	p1 [Pa]" in line:
                break

            if in_arrivals and "#" in line:
                arrivals = arrivals + [[float(item) for item in line.replace("# ","").split(" ")]]

            if "# incl [deg]" in line:
                in_arrivals = True

    arrivals = np.array(arrivals)

    arr_tms = arrivals[:, 5]
    arr_tr_vel = c0 / np.cos(np.radians(arrivals[:, 8]))
    arr_back_az = arrivals[:, 9]

    arr_amplitudes = np.exp((arrivals[:, 10] + arrivals[:, 11]) / 10.0)
    arr_amplitudes = arr_amplitudes / np.sum(arr_amplitudes)

    det_sets = []
    det_indices = np.arange(len(arr_tms))
    while len(det_indices) > 0:
        merge_indices = [det_indices[0]]
        for k in det_indices[1:]:
            if abs(arr_tms[det_indices[0]] - arr_tms[k]) < 10.0:
                merge_indices = merge_indices + [k]

        det_sets = det_sets + [merge_indices]
        det_indices = [m for m in det_indices if m not in merge_indices]

    # Cycle through detection sets and write detection info to InfraPy JSON file format
    wvfrms = np.loadtxt(wvfrm_file)
    wvfrm_tm = wvfrms[:, 0]

    det_cnt = 0

    dets_out = []
    for det_set in det_sets:
        if len(det_set) > 1:
            set_weights = np.array([arr_amplitudes[k] for k in det_set])
            det_back_az = np.average(np.array([arr_back_az[k] for k in det_set]), weights=set_weights)
            det_tr_vel = np.average(np.array([arr_tr_vel[k] for k in det_set]), weights=set_weights)

            det_tms = np.array([arr_tms[k] for k in det_set])
            det_dt = det_tms[np.argmax(set_weights)]
            det_tm = origin_time + det_dt

            det_start = (np.min(det_tms) - det_dt) - 15.0
            det_end = (np.max(det_tms) - det_dt) + 15.0
            
            combined_wvfrm = np.sum(np.array([wvfrms[:, k + 1] for k in det_set]), axis=0)
            det_fstat = 1.0

        else:
            k = det_set[0]
            det_back_az = arr_back_az[k]
            det_tr_vel = arr_tr_vel[k]

            det_dt = arr_tms[k]
            det_tm = origin_time + det_dt

            det_start = -15.0
            det_end = 15.0

            combined_wvfrm = wvfrms[:, det_set[0] + 1]
            det_fstat = 1.0

        # Load IMS noise model and compute SNR to get f-stat
        dt = (wvfrm_tm[1] - wvfrm_tm[0])
        freq = np.fft.rfftfreq(len(wvfrm_tm), d=dt)
        spec = 10.0 * np.log10(abs(np.fft.rfft(combined_wvfrm)) * dt)

        ims_ns = np.loadtxt(find_spec('stochprop').submodule_search_locations[0] + '/resources/IMSNOISE_MIN_MED_MAX.txt')
        if ns_option == "min":
            ns_ref = interp1d(ims_ns[:, 0], 0.5 * (10.0 * ims_ns[:, 1] + 10.0 * np.log10((det_end - det_start) / 1.0)), bounds_error=False, fill_value="extrapolate")
        elif ns_option == "max":
            ns_ref = interp1d(ims_ns[:, 0], 0.5 * (10.0 * ims_ns[:, 3] + 10.0 * np.log10((det_end - det_start) / 1.0)), bounds_error=False, fill_value="extrapolate")
        else:
            ns_ref = interp1d(ims_ns[:, 0], 0.5 * (10.0 * ims_ns[:, 2] + 10.0 * np.log10((det_end - det_start) / 1.0)), bounds_error=False, fill_value="extrapolate")

        freq_mask = np.logical_and(0.005 < freq, freq < 8.0)
        SNR = 10.0 ** (np.max(spec[freq_mask] - ns_ref(freq[freq_mask])) / 10.0)
      
        # Fill out detection dictionary
        tm_mask = np.logical_and(det_start - 30.0 < wvfrm_tm - det_dt, wvfrm_tm - det_dt < 30.0 + det_end)

        if SNR > 1.0:
            dets_out = dets_out + [{}]

            dets_out[-1]['peak f-stat time'] = det_tm
            dets_out[-1]['start/end'] = [[det_start, det_end]]
            dets_out[-1]['f-stat'] = SNR * (array_dim - 1)

            dets_out[-1]['back az'] = det_back_az
            dets_out[-1]['tr vel'] = det_tr_vel

            dets_out[-1]['fk'] = [{}]
            dets_out[-1]['fk'][0]['time'] = [0.0]
            dets_out[-1]['fk'][0]['back az'] = [det_back_az]
            dets_out[-1]['fk'][0]['tr vel'] = [det_tr_vel]
            dets_out[-1]['fk'][0]['f-stat'] = [det_fstat]

            dets_out[-1]['beam'] = [{}]
            dets_out[-1]['beam'][0]['time'] = wvfrm_tm[tm_mask] - det_dt
            dets_out[-1]['beam'][0]['signal'] = combined_wvfrm[tm_mask]
            dets_out[-1]['beam'][0]['resid'] = np.zeros_like(wvfrm_tm)[tm_mask]

            dets_out[-1]['spec'] = [{}]
            dets_out[-1]['spec'][0]['freq'] = freq[freq_mask]
            dets_out[-1]['spec'][0]['signal'] = 10**(spec[freq_mask] / 10.0)
            dets_out[-1]['spec'][0]['resid'] =  10**(ns_ref(freq[freq_mask]) / 10.0)

            dets_out[-1]['wvfrm_info'] = [wvfrm_dicts]
            dets_out[-1]['fk_params'] = [dict(fk_dict)]
            dets_out[-1]['det_params'] = [det_dict]

            pos_SNR = dets_out[-1]['spec'][0]['signal'] > dets_out[-1]['spec'][0]['resid']
            dets_out[-1]['fk_params'][0]["freq_min"] = np.round(dets_out[-1]['spec'][0]['freq'][pos_SNR][0], 3)
            dets_out[-1]['fk_params'][0]["freq_max"] = np.round(dets_out[-1]['spec'][0]['freq'][pos_SNR][-1], 3)

    return dets_out
