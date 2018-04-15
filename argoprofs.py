import os
import glob
import datetime
from netCDF4 import Dataset as ncds
import numpy as np
import gsw
import seawater

# class Argoprof:
#     listofdac = ['aoml', 'bodc', 'coriolis', 'csio', 'csiro', 'incois', 'jma', 'kma', 'kordi', 'meds', 'nmdis']
#     qc_accept_flag = 2
#
#     def __init__(self, datapath = '', Float_ID = '', DAC = None, Cycle = 0):
#         self.datapath = datapath
#         self.Float_ID = Float_ID
#         self.DAC = DAC
#         self.Cycle = Cycle
#         self._getfn()
#         self._getiprof()
#
#         self.coords = {'lon': None, 'lat': None, 'pos_qc': None, 'time': None, 'time_qc': None}
#         self.qc = {'Direction': None, 'Data_mode': None, 'Data_state_indicator': None, 'Profile_pres_qc': None, 'Profile_psal_qc': None, 'Profile_temp_qc': None}
#         self.vars = {'pres': None, 'pres_qc': None, 'pres_adjusted': None, 'pres_adjusted_qc': None, 'pres_adjusted_error': None,
#                      'psal': None, 'psal_qc': None, 'psal_adjusted': None, 'psal_adjusted_qc': None, 'psal_adjusted_error': None,
#                      'temp': None, 'temp_qc': None, 'temp_adjusted': None, 'temp_adjusted_qc': None, 'temp_adjusted_error': None}
#         self.vars['SA'] = None
#         self.vars['CT'] = None
#
#     @property
#     def DAC(self):
#         return self._DAC
#
#     @DAC.setter
#     def DAC(self, value):
#         self._DAC = value
#         if not self._DAC in Argoprof.listofdac:
#             for idac in Argoprof.listofdac:
#                 if os.path.isdir(self.datapath + '/' + idac + '/' + self.Float_ID):
#                     self._DAC = idac
#                     break
#
#     def _getfn(self):
#         self._filename = ''
#         fndir = self.datapath + '/' + self.DAC + '/' + self.Float_ID + '/'
#         fn_prof = fndir + self.Float_ID + '_prof.nc'
#         fn_ind_D = fndir + 'profiles/D' + self.Float_ID + '_' + "{:03d}".format(self.Cycle) + '.nc'
#         fn_ind_R = fndir + 'profiles/R' + self.Float_ID + '_' + "{:03d}".format(self.Cycle) + '.nc'
#         if os.path.isfile(fn_prof):
#             self._filename = fn_prof
#         elif os.path.isfile(fn_ind_D):
#             self._filename = fn_ind_D
#         elif os.path.isfile(fn_ind_R):
#             self._filename = fn_ind_R
#         else:
#             print('Profile file not found! \n' + fn_prof + '\n' + fn_ind_D + '\n' + fn_ind_R + '\nCheck datapath, Float_ID, DAC and cycle.')
#         return
#
#     def _getiprof(self):
#         if self._filename == '':
#             self._iprof = -1
#         elif self._filename[-7:-3] == 'prof':
#             iprof = np.where(ncds(self._filename).variables['CYCLE_NUMBER'][:] == self.Cycle)[0]
#             if iprof.size > 0:
#                 self._iprof = iprof
#             else:
#                 self._iprof = -1
#                 print('Cannot find cycle ' + str(self.Cycle)+ ' in file ' + self._filename)
#         else:
#             self._iprof = 0
#         return
#
#     def readinfo(self, varnm):
#         ncgrp = ncds(self._filename)
#         if ncgrp.variables[varnm].dimensions[0] == 'N_PROF':
#             return ncgrp.variables[varnm][self._iprof, ...]
#         else:
#             return ncgrp.variables[varnm][:]
#
#     def load_coords(self):
#         time_ref = self.readinfo('REFERENCE_DATE_TIME')
#         juld = self.readinfo('JULD')
#         self.coords['time'] = datetime.datetime.strptime(time_ref.tobytes().decode(), '%Y%m%d%H%M%S') + datetime.timedelta(days = juld[0])
#         # juld_loc = self.readinfo('JULD_LOCATION')
#         # self.coords['time'] = datetime.datetime.strptime(time_ref.tobytes().decode(), '%Y%m%d%H%M%S') + datetime.timedelta(days = juld_loc[0])
#         self.coords['time_qc'] = self.readinfo('JULD_QC')
#         self.coords['lon'] = self.readinfo('LONGITUDE')
#         self.coords['lat'] = self.readinfo('LATITUDE')
#         self.coords['pos_qc'] = self.readinfo('POSITION_QC')
#
#     def load_qc(self):
#         self.qc['Direction'] = self.readinfo('DIRECTION')
#         self.qc['Data_mode'] = self.readinfo('DATA_MODE')
#         self.qc['Data_state_indicator'] = self.readinfo('DATA_STATE_INDICATOR')
#         self.qc['Profile_pres_qc'] = self.readinfo('PROFILE_PRES_QC')
#         self.qc['Profile_psal_qc'] = self.readinfo('PROFILE_PSAL_QC')
#         self.qc['Profile_temp_qc'] = self.readinfo('PROFILE_TEMP_QC')
#
#     def load_vars(self):
#         self.vars['pres'] = self.readinfo('PRES')
#         self.vars['psal'] = self.readinfo('PSAL')
#         self.vars['temp'] = self.readinfo('TEMP')
#         self.vars['pres_qc'] = self.readinfo('PRES_QC')
#         self.vars['psal_qc'] = self.readinfo('PSAL_QC')
#         self.vars['temp_qc'] = self.readinfo('TEMP_QC')
#
#     def load_vars_adjusted(self):
#         self.vars['pres_adjusted'] = self.readinfo('PRES_ADJUSTED')
#         self.vars['psal_adjusted'] = self.readinfo('PSAL_ADJUSTED')
#         self.vars['temp_adjusted'] = self.readinfo('TEMP_ADJUSTED')
#         self.vars['pres_adjusted_qc'] = self.readinfo('PRES_ADJUSTED_QC')
#         self.vars['psal_adjusted_qc'] = self.readinfo('PSAL_ADJUSTED_QC')
#         self.vars['temp_adjusted_qc'] = self.readinfo('TEMP_ADJUSTED_QC')
#         self.vars['pres_adjusted_error'] = self.readinfo('PRES_ADJUSTED_ERROR')
#         self.vars['psal_adjusted_error'] = self.readinfo('PSAL_ADJUSTED_ERROR')
#         self.vars['temp_adjusted_error'] = self.readinfo('TEMP_ADJUSTED_ERROR')
#
#     def qc_coords(self):
#         pass
#
#     def qc_vars_adjusted(self):
#         mask_p = self.vars['pres_adjusted_qc'].filled(b'9').astype(np.int) > Argoprof.qc_accept_flag
#         mask_s = self.vars['psal_adjusted_qc'].filled(b'9').astype(np.int) > Argoprof.qc_accept_flag
#         mask_t = self.vars['temp_adjusted_qc'].filled(b'9').astype(np.int) > Argoprof.qc_accept_flag
#         mask = mask_t | mask_s | mask_p
#
#         self.vars['pres_adjusted'].mask = mask
#         self.vars['psal_adjusted'].mask = mask
#         self.vars['temp_adjusted'].mask = mask
#         self.vars['pres_adjusted_error'].mask = mask
#         self.vars['psal_adjusted_error'].mask = mask
#         self.vars['temp_adjusted_error'].mask = mask
#
#     @staticmethod
#     def sw_turner(s, t, p):
#         smid = (s[:-1] + s[1:]) / 2
#         tmid = (t[:-1] + t[1:]) / 2
#         pmid = (p[:-1] + p[1:]) / 2
#
#         tu = seawater.eos80.ptmp(s[:-1], t[:-1], p[:-1], pr=pmid)
#         td = seawater.eos80.ptmp(s[1:] , t[1:] , p[1:] , pr=pmid)
#
#         ds = s[:-1] - s[1:]
#         dt = tu - td
#
#         alpha = seawater.eos80.alpha(smid, tmid, pmid, pt=False)
#         beta  = seawater.eos80.beta (smid, tmid, pmid, pt=False)
#
#         Tu = np.arctan2(alpha * dt + beta * ds, alpha * dt - beta * ds) * 180 / np.pi
#         Rsubrho = (alpha * dt) / (beta * ds)
#
#         return Tu, Rsubrho, pmid
#
#     def calc_SACT(self):
#         self.vars['SA'] = gsw.SA_from_SP(self.vars['psal_adjusted'].compressed(), self.vars['pres_adjusted'].compressed(), self.coords['lon'], self.coords['lat'])
#         self.vars['CT'] = gsw.CT_from_t(self.vars['SA'], self.vars['temp_adjusted'].compressed(), self.vars['pres_adjusted'].compressed())
#
#     def calc_n2_gsw(self):
#         if self.vars['SA'] is None or self.vars['CT'] is None:
#             self.calc_SACT()
#         self.vars['N2'], self.vars['pmid'] = gsw.stability.Nsquared(self.vars['SA'], self.vars['CT'], self.vars['pres_adjusted'].compressed())
#
#     def calc_TA_gsw(self):
#         if self.vars['SA'] is None or self.vars['CT'] is None:
#             self.calc_SACT()
#         self.vars['TA'], self.vars['Rsubrho'], self.vars['pmid'] = gsw.stability.Nsquared(self.vars['SA'], self.vars['CT'], self.vars['pres_adjusted'].compressed())
#
#     def calc_n2_sw(self):
#         self.vars['N2'], q, self.vars['pmid'] = seawater.geostrophic.bfrq(self.vars['psal_adjusted'].compressed(), self.vars['temp_adjusted'].compressed(), self.vars['pres_adjusted'].compressed())
#
#     def calc_TA_sw(self):
#         self.vars['TA'], self.vars['Rsubrho'], self.vars['pmid'] = Argoprof.sw_turner(self.vars['psal_adjusted'].compressed(), self.vars['temp_adjusted'].compressed(), self.vars['pres_adjusted'].compressed())

class Argofloat:
    listofdac = ['aoml', 'bodc', 'coriolis', 'csio', 'csiro', 'incois', 'jma', 'kma', 'kordi', 'meds', 'nmdis']
    qc_accept_flag = 1

    def __init__(self, datapath = '', Float_ID = '', DAC = None):
        self.datapath = datapath
        self.Float_ID = Float_ID
        self.DAC = DAC
        self._getfn()

        self.coords = {'lon': None, 'lat': None, 'pos_qc': None, 'time': None, 'time_qc': None, 'cycle': None, 'direction': None}
        self.qc = {'Data_mode': None, 'Data_state_indicator': None, 'Profile_pres_qc': None, 'Profile_psal_qc': None, 'Profile_temp_qc': None}
        self.vars = {'pres': None, 'pres_qc': None, 'pres_adjusted': None, 'pres_adjusted_qc': None, 'pres_adjusted_error': None,
                     'psal': None, 'psal_qc': None, 'psal_adjusted': None, 'psal_adjusted_qc': None, 'psal_adjusted_error': None,
                     'temp': None, 'temp_qc': None, 'temp_adjusted': None, 'temp_adjusted_qc': None, 'temp_adjusted_error': None}
        self.vars['SA'] = None
        self.vars['CT'] = None

    @property
    def DAC(self):
        return self._DAC

    @DAC.setter
    def DAC(self, value):
        self._DAC = value
        if not self._DAC in Argofloat.listofdac:
            for idac in Argofloat.listofdac:
                if os.path.isdir(self.datapath + '/' + idac + '/' + self.Float_ID):
                    self._DAC = idac
                    break

    def _getfn(self):
        fndir = self.datapath + '/' + self.DAC + '/' + self.Float_ID
        fn_prof = fndir + '/' + self.Float_ID + '_prof.nc'
        fn_ind  = glob.glob(fndir + '/profiles' + '/*.nc')
        if os.path.isfile(fn_prof):
            self._filename = [fn_prof]
            self._nprof = ncds(self._filename[0]).dimensions['N_PROF'].size
            self._nfile = 1
        elif len(fn_ind) > 0:
            fn_ind_rt = glob.glob(fndir + '/profiles/' + 'R*.nc')
            for ifl in fn_ind_rt:
                if ifl.replace('/R', '/D') in fn_ind:
                    del fn_ind[fn_ind.index(ifl)]
            self._filename = fn_ind
            self._nprof = len(fn_ind)
            self._nfile = len(fn_ind)
        else:
            self._filename = []
            self._nprof = 0
            self._nfile = 0
            print('Profile file not found! \n' + fn_prof + '\n' + fn_ind)
        return

    def readnc(self, varnm):
        varlist = [[] for _ in varnm]
        vardim = [[] for _ in varnm]
        nlv = 0

        for fn, iprof in zip(self._filename, range(self._nprof)):
            ncgrp = ncds(fn)

            nlv = max(nlv, ncgrp.dimensions['N_LEVELS'].size)

            for iv, ivnm in zip(range(len(varnm)), varnm):
                varlist[iv].append(ncgrp.variables[ivnm][:])
                vardim [iv].append(ncgrp.variables[ivnm].dimensions)
            ncgrp.close()

        var = []
        if len(self._filename) == 1:
            for iv in range(len(varnm)):
                var.append(varlist[iv][0])
        else:
            for iv in range(len(varnm)):
                if 'N_LEVELS' in vardim[iv][0]:
                    var.append(np.ma.masked_values(np.zeros((self._nprof, nlv))+99999., 99999.))
                    for iprof in range(self._nprof):
                        var[iv][iprof, :varlist[iv][iprof].size] = varlist[iv][iprof]
                else:
                    var.append(np.stack(varlist[iv]))
        return var

    def load_coords(self):
        varnm = ['CYCLE_NUMBER', 'DIRECTION', 'POSITION_QC', 'LONGITUDE', 'LATITUDE', 'JULD_QC', 'REFERENCE_DATE_TIME', 'JULD']
        self.coords['cycle'], self.coords['direction'], \
        self.coords['pos_qc'], self.coords['lon'], self.coords['lat'], \
        self.coords['time_qc'], reftime, juld = self.readnc(varnm)

        # self.coords['time'] = np.zeros((self._nprof,), dtype = 'datetime64')
        self.coords['time'] = []
        for iprof in range(self._nprof):
            if self.coords['time_qc'][iprof] > b'1':
                self.coords['time'].append(None)
                continue
            if len(reftime.shape) == 1:
                self.coords['time'].append(datetime.datetime.strptime(reftime.tobytes().decode(), '%Y%m%d%H%M%S') + datetime.timedelta(days = juld[iprof]))
            else:
                self.coords['time'].append(datetime.datetime.strptime(reftime[iprof].tobytes().decode(), '%Y%m%d%H%M%S') + datetime.timedelta(days = juld[iprof][0]))

    def load_qc(self):
        varnm = ['DATA_MODE', 'DATA_STATE_INDICATOR', 'PROFILE_PRES_QC', 'PROFILE_PSAL_QC', 'PROFILE_TEMP_QC']
        self.qc['Data_mode'], self.qc['Data_state_indicator'], self.qc['Profile_pres_qc'], self.qc['Profile_psal_qc'], self.qc['Profile_temp_qc'] = self.readnc(varnm)

    def load_vars(self):
        varnm = ['PRES', 'PSAL', 'TEMP', 'PRES_QC', 'PSAL_QC', 'TEMP_QC']
        self.vars['pres'], self.vars['psal'], self.vars['temp'], self.vars['pres_qc'], self.vars['psal_qc'], self.vars['temp_qc'] = self.readnc(varnm)

    def load_vars_adjusted(self):
        varnm = ['PRES_ADJUSTED', 'PSAL_ADJUSTED', 'TEMP_ADJUSTED', 'PRES_ADJUSTED_QC', 'PSAL_ADJUSTED_QC', 'TEMP_ADJUSTED_QC',
                 'PRES_ADJUSTED_ERROR', 'PSAL_ADJUSTED_ERROR', 'TEMP_ADJUSTED_ERROR']

        self.vars['pres_adjusted'], self.vars['psal_adjusted'], self.vars['temp_adjusted'], \
        self.vars['pres_adjusted_qc'], self.vars['psal_adjusted_qc'], self.vars['temp_adjusted_qc'], \
        self.vars['pres_adjusted_error'], self.vars['psal_adjusted_error'], self.vars['temp_adjusted_error'] = \
        self.readnc(varnm)

        if not np.ma.is_masked(self.vars['pres_adjusted']):
            self.vars['pres_adjusted'] = np.ma.masked_values(self.vars['pres_adjusted'], 99999.)
        if not np.ma.is_masked(self.vars['psal_adjusted']):
            self.vars['psal_adjusted'] = np.ma.masked_values(self.vars['psal_adjusted'], 99999.)
        if not np.ma.is_masked(self.vars['temp_adjusted']):
            self.vars['temp_adjusted'] = np.ma.masked_values(self.vars['temp_adjusted'], 99999.)
        if not np.ma.is_masked(self.vars['pres_adjusted_error']):
            self.vars['pres_adjusted_error'] = np.ma.masked_values(self.vars['pres_adjusted_error'], 99999.)
        if not np.ma.is_masked(self.vars['psal_adjusted_error']):
            self.vars['psal_adjusted_error'] = np.ma.masked_values(self.vars['psal_adjusted_error'], 99999.)
        if not np.ma.is_masked(self.vars['temp_adjusted_error']):
            self.vars['temp_adjusted_error'] = np.ma.masked_values(self.vars['temp_adjusted_error'], 99999.)

    def qc_coords(self):
        pass

    def qc_vars_adjusted(self):
        if np.ma.is_masked(self.vars['pres_adjusted_qc']):
            self.vars['pres_adjusted_qc'].mask[(self.vars['pres_adjusted_qc'] < b'1') | (self.vars['pres_adjusted_qc'] > b'9')] = True
            mask_p = self.vars['pres_adjusted_qc'].filled(b'9').astype(np.int) > Argofloat.qc_accept_flag
        else:
            self.vars['pres_adjusted_qc'][(self.vars['pres_adjusted_qc'] < b'1') | (self.vars['pres_adjusted_qc'] > b'9')] = b'9'
            mask_p = self.vars['pres_adjusted_qc'].astype(np.int) > Argofloat.qc_accept_flag
        if np.ma.is_masked(self.vars['psal_adjusted_qc']):
            self.vars['psal_adjusted_qc'].mask[(self.vars['psal_adjusted_qc'] < b'1') | (self.vars['psal_adjusted_qc'] > b'9')] = True
            mask_s = self.vars['psal_adjusted_qc'].filled(b'9').astype(np.int) > Argofloat.qc_accept_flag
        else:
            self.vars['psal_adjusted_qc'][(self.vars['psal_adjusted_qc'] < b'1') | (self.vars['psal_adjusted_qc'] > b'9')] = b'9'
            mask_s = self.vars['psal_adjusted_qc'].astype(np.int) > Argofloat.qc_accept_flag
        if np.ma.is_masked(self.vars['temp_adjusted_qc']):
            self.vars['temp_adjusted_qc'].mask[(self.vars['temp_adjusted_qc'] < b'1') | (self.vars['temp_adjusted_qc'] > b'9')] = True
            mask_t = self.vars['temp_adjusted_qc'].filled(b'9').astype(np.int) > Argofloat.qc_accept_flag
        else:
            self.vars['temp_adjusted_qc'][(self.vars['temp_adjusted_qc'] < b'1') | (self.vars['temp_adjusted_qc'] > b'9')] = b'9'
            mask_t = self.vars['temp_adjusted_qc'].astype(np.int) > Argofloat.qc_accept_flag
        mask_dm = self.qc['Data_mode'] != b'D'
        mask_time = self.coords['time_qc'] > b'1'
        mask_pos  = self.coords['pos_qc'] > b'1'
        mask = mask_t | mask_s | mask_p | mask_dm[:, None] | mask_time[:, None] | mask_pos[:, None]

        self.vars['pres_adjusted'].mask = mask
        self.vars['psal_adjusted'].mask = mask
        self.vars['temp_adjusted'].mask = mask
        self.vars['pres_adjusted_error'].mask = mask
        self.vars['psal_adjusted_error'].mask = mask
        self.vars['temp_adjusted_error'].mask = mask

        # sort along depth + moving all the invalid to the bottom
        sort_idx = self.vars['pres_adjusted'].filled().argsort(axis = -1)
        sort_idx_c = np.arange(self._nprof)

        self.vars['pres_adjusted'] = self.vars['pres_adjusted'][sort_idx_c[:, None], sort_idx]
        self.vars['psal_adjusted'] = self.vars['psal_adjusted'][sort_idx_c[:, None], sort_idx]
        self.vars['temp_adjusted'] = self.vars['temp_adjusted'][sort_idx_c[:, None], sort_idx]
        self.vars['pres_adjusted_error'] = self.vars['pres_adjusted_error'][sort_idx_c[:, None], sort_idx]
        self.vars['psal_adjusted_error'] = self.vars['psal_adjusted_error'][sort_idx_c[:, None], sort_idx]
        self.vars['temp_adjusted_error'] = self.vars['temp_adjusted_error'][sort_idx_c[:, None], sort_idx]
        self.vars['pres_adjusted_qc'] = self.vars['pres_adjusted_qc'][sort_idx_c[:, None], sort_idx]
        self.vars['psal_adjusted_qc'] = self.vars['psal_adjusted_qc'][sort_idx_c[:, None], sort_idx]
        self.vars['temp_adjusted_qc'] = self.vars['temp_adjusted_qc'][sort_idx_c[:, None], sort_idx]

    @staticmethod
    def sw_turner(s, t, p):
        smid = (s[:-1] + s[1:]) / 2
        tmid = (t[:-1] + t[1:]) / 2
        pmid = (p[:-1] + p[1:]) / 2
        dz = p[:-1] - p[1:]

        tu = seawater.eos80.ptmp(s[:-1], t[:-1], p[:-1], pr=pmid)
        td = seawater.eos80.ptmp(s[1:] , t[1:] , p[1:] , pr=pmid)

        ds = s[:-1] - s[1:]
        dt = tu - td

        alpha = seawater.eos80.alpha(smid, tmid, pmid, pt=False)
        beta  = seawater.eos80.beta (smid, tmid, pmid, pt=False)

        Tu = np.arctan2(alpha * dt + beta * ds, alpha * dt - beta * ds) * 180 / np.pi
        Rsubrho = (alpha * dt) / (beta * ds)
        adTdz = alpha * dt / dz
        bdSdz = beta * ds / dz

        return Tu, Rsubrho, adTdz, bdSdz, pmid

    @staticmethod
    def sw_mld(s, t, p, pd_c):
        rho0 = seawater.dens0(s, t)
        rho0_10m = np.interp(10., p, rho0)

        depm = p[-1]
        for iz in range(1, p.size, 1):
            if p[iz] > 10. and rho0[iz] - rho0_10m > pd_c:
                depm = np.interp(rho0_10m + pd_c, np.array([rho0[iz], rho0[iz-1]]), np.array([p[iz], p[iz-1]]))
                break

        mld = depm
        return mld

    def calc_SACT(self):
        self.vars['SA'] = np.zeros_like(self.vars['psal_adjusted'])
        self.vars['CT'] = np.zeros_like(self.vars['psal_adjusted'])

        for iprof in range(self._nprof):
            s = self.vars['psal_adjusted'][iprof, :]
            p = self.vars['pres_adjusted'][iprof, :]
            t = self.vars['temp_adjusted'][iprof, :]
            lon = self.coords['lon'][iprof]
            lat = self.coords['lat'][iprof]

            self.vars['SA'][iprof, :] = gsw.SA_from_SP(s, p, lon, lat)
            self.vars['CT'][iprof, :] = gsw.CT_from_t(self.vars['SA'][iprof, :], t, p)

    def calc_n2_gsw(self):
        if self.vars['SA'] is None or self.vars['CT'] is None:
            self.calc_SACT()
        self.vars['N2'], self.vars['pmid'] = gsw.stability.Nsquared(self.vars['SA'], self.vars['CT'], self.vars['pres_adjusted'], axis = 1)

    def calc_TA_gsw(self):
        if self.vars['SA'] is None or self.vars['CT'] is None:
            self.calc_SACT()
        self.vars['TA'], self.vars['Rsubrho'], self.vars['pmid'] = gsw.stability.Turner_Rsubrho(self.vars['SA'], self.vars['CT'], self.vars['pres_adjusted'], axis = 1)

    def calc_n2_sw(self):
        self.vars['N2'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
        self.vars['pmid'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])

    def calc_mld(self):
        pd_c = 0.03
        self.vars['mld'] = np.zeros((self._nprof, ))

        for iprof in range(self._nprof):
            s = self.vars['psal_adjusted'][iprof,:].filled(np.nan)
            t = self.vars['temp_adjusted'][iprof,:].filled(np.nan)
            p = self.vars['pres_adjusted'][iprof,:].filled(np.nan)

            if np.all(np.isnan(s)) or np.all(np.isnan(t)) or np.all(np.isnan(p)):
                self.vars['mld'][iprof] = np.nan
            else:
                self.vars['mld'][iprof] = Argofloat.sw_mld(s, t, p, pd_c)

    def calc_TA_sw(self):
        self.vars['TA'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
        self.vars['Rsubrho'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
        self.vars['adtdz'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
        self.vars['bdsdz'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
        self.vars['pmid'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])

        for iprof in range(self._nprof):
            s = self.vars['psal_adjusted'][iprof,:].filled(np.nan)
            t = self.vars['temp_adjusted'][iprof,:].filled(np.nan)
            p = self.vars['pres_adjusted'][iprof,:].filled(np.nan)
            # s = self.vars['psal_adjusted'][iprof,:].compressed()
            # t = self.vars['temp_adjusted'][iprof,:].compressed()
            # p = self.vars['pres_adjusted'][iprof,:].compressed()

            # ta, rsubrho, pmid = Argofloat.sw_turner(s, t, p)
            ta, rsubrho, dtdz, dsdz, pmid = Argofloat.sw_turner(s, t, p)
            self.vars['TA'][iprof,:] = np.squeeze(ta)
            self.vars['Rsubrho'][iprof,:] = np.squeeze(rsubrho)
            self.vars['adtdz'][iprof,:] = np.squeeze(dtdz)
            self.vars['bdsdz'][iprof,:] = np.squeeze(dsdz)
            self.vars['pmid'][iprof,:] = np.squeeze(pmid)

    def calc_rho0_sw(self):
        s = self.vars['psal_adjusted'].filled(np.nan)
        t = self.vars['temp_adjusted'].filled(np.nan)
        self.vars['rho0'] = seawater.dens0(s, t)
