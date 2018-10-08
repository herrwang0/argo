import os
import glob
import calendar
import datetime
from netCDF4 import Dataset as ncds
import numpy as np
import gsw
import seawater as sw
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from ncep.utils import grid as ncep_grid
import matlab.engine

class Argoprof:
    qc_accept_flag = [1, 2, 5, 8] # can be modified from the outside
    def __init__(self):
        self._nlv = 0
        self.coords = {'lon': None, 'lat': None, 'pos_qc': None, 'time': None, 'time_qc': None, 'cycle': None, 'direction': None}
        self.qc     = {'Data_mode': None, 'Data_state_indicator': None, 'Profile_pres_qc': None, 'Profile_psal_qc': None, 'Profile_temp_qc': None}
        self.vars   = {'pres': None, 'pres_qc': None, 'pres_adjusted': None, 'pres_adjusted_qc': None, 'pres_adjusted_error': None,
                       'psal': None, 'psal_qc': None, 'psal_adjusted': None, 'psal_adjusted_qc': None, 'psal_adjusted_error': None,
                       'temp': None, 'temp_qc': None, 'temp_adjusted': None, 'temp_adjusted_qc': None, 'temp_adjusted_error': None}

    def qc_coords(self):
        mask_time = np.full(self.coords['time_qc'].shape, True)
        mask_pos  = np.full(self.coords['pos_qc' ].shape, True)
        for flg in self.qc_accept_flag:
            mask_time = mask_time & ~(self.coords['time_qc'].astype(np.int) == flg)
            mask_pos  = mask_time & ~(self.coords['pos_qc' ].astype(np.int) == flg)
        return mask_pos, mask_time

    def qc_vars(self, varnm):
        # QC mask1: remove none numbers. QC mask2: in qc_accept_flag
        mask_nan = self.vars[varnm + '_qc'].mask | (self.vars[varnm + '_qc'] < b'1') | (self.vars[varnm + '_qc'] > b'9')

        mask_qc = np.full(self.vars[varnm].shape, True)
        for flg in self.qc_accept_flag:
            mask_qc = mask_qc & ~(self.vars[varnm + '_qc'].filled(b'9').astype(np.int) == flg)
        return mask_nan | mask_qc

    def sort_lev(self, ref):
        # sort along depth + moving all the invalid to the bottom
        sort_idx = self.vars[ref].filled().argsort(axis = -1) # filled makes the masked element to maked_value, in contrast with data (which uses the true value of the element, even if it is masked)
        return sort_idx

    def calc_rho0(self, varsuffix = '_adjusted', toolbox = 'sw'):
        if toolbox == 'sw':
            self.vars['rho0'] = sw.dens0(self.vars['psal_adjusted'], self.vars['temp_adjusted'])
        elif toolbox == 'gsw': #
            pass
        else: # warning here
            pass

    def bin_idx(self, lonb, latb):
        mask_pos, mask_time = self.qc_coords()
        mask_coords = mask_pos | mask_time

        ix = np.zeros_like(mask_coords, dtype = np.int64)
        iy = np.zeros_like(mask_coords, dtype = np.int64)

        lonp = np.copy(self.coords['lon'])
        latp = np.copy(self.coords['lat'])

        lonp[lonp <= lonb.max() - 360] += 360
        lonp[lonp >= lonb.min() + 360] -= 360

        for iprof in range(mask_coords.size):
            if mask_coords[iprof]: continue

            ix[iprof] = np.argmax(lonb > lonp[iprof]) - 1
            if   np.floor(latp[iprof]) <  latb[ 0]:
                iy[iprof] = 0
            elif np.floor(latp[iprof]) >= latb[-1]:
                iy[iprof] = latb.size - 1
            else:
                iy[iprof] = np.argmax(latb > latp[iprof]) - 1
        return ix, iy

class Argofloat(Argoprof):
    listofdac = ['aoml', 'bodc', 'coriolis', 'csio', 'csiro', 'incois', 'jma', 'kma', 'kordi', 'meds', 'nmdis']

    def __init__(self, datapath = '', Float_ID = '', DAC = None):
        Argoprof.__init__(self)

        self.datapath = datapath
        self.Float_ID = Float_ID
        self.DAC = DAC
        self._getfn()

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

    # Generic netCDF reading function. For _prof.nc file, the function simply loads the corresponding variables in the varnm list.
    # For individual profile files, the function goes through each file and concatenate them into one with the maximum level.
    # All arrays are masked.
    def readnc(self, varnm):
        varlist = [None for ivar in varnm]
        if self._nfile == 1 and self._filename[0][-8:-3] == '_prof':
            ncgrp = ncds(self._filename[0])
            self._nlv = ncgrp.dimensions['N_LEVELS'].size
            for ivar, iv in zip(varnm, range(len(varnm))):
                var = ncgrp.variables[ivar]
                if len(var.shape) == 1:
                    if var.dimensions[0] == 'N_PROF':
                        varlist[iv] = var[:][..., None]
                    else:
                        varlist[iv] = np.tile(var[:], [self._nprof, 1])
                else:
                    if np.ma.isMaskedArray(var[:]):
                        var_masked = var[:]
                    else:
                        if var.dtype.type is np.bytes_:
                            var_masked = np.ma.masked_where(var[:] == b'', var[:])
                            var_masked.fill_value = b' '
                        else:  # need a more generic
                            var_masked = np.ma.masked_values(var[:], 99999.)
                    varlist[iv] = var_masked
            ncgrp.close()
        else:
            # Calculate maximum level number on the first run
            if self._nlv == 0:
                for fn in self._filename:
                    ncgrp = ncds(fn)
                    self._nlv = max(self._nlv, ncgrp.dimensions['N_LEVELS'].size)
                    ncgrp.close()

            for fn, iprof in zip(self._filename, range(self._nprof)):
                ncgrp = ncds(fn)
                for ivar, iv in zip(varnm, range(len(varnm))):
                    var = ncgrp.variables[ivar]

                    if iprof == 0:
                        if 'N_LEVELS' in var.dimensions:
                            if var.dtype.type is np.bytes_:
                                varlist[iv] = np.ma.masked_where(np.zeros((self._nprof, self._nlv), dtype = var.dtype) == b'', np.zeros((self._nprof, self._nlv), dtype = var.dtype))
                                varlist[iv].fill_value = b' '
                            else:
                                varlist[iv] = np.ma.masked_values(np.zeros((self._nprof, self._nlv))+99999., 99999.)
                        else:
                            varlist[iv] = np.zeros((self._nprof, var[:].size), dtype = var.dtype)

                    varlist[iv][iprof, :var[:].size] = var[:]
                ncgrp.close()
        return varlist

    def load_coords(self):
        varnm = ['CYCLE_NUMBER', 'DIRECTION', 'POSITION_QC', 'LONGITUDE', 'LATITUDE', 'JULD_QC', 'REFERENCE_DATE_TIME', 'JULD']
        self.coords['cycle'], self.coords['direction'], \
        self.coords['pos_qc'], self.coords['lon'], self.coords['lat'], \
        self.coords['time_qc'], reftime, juld = self.readnc(varnm)

        # self.coords['time'] = np.zeros((self._nprof,), dtype = 'datetime64')
        self.coords['time'] = []
        for iprof in range(self._nprof):
            # if self.coords['time_qc'][iprof] > b'1':
            #     self.coords['time'].append(None)
            #     continue
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

        self.vars['pres_adjusted'      ], self.vars['psal_adjusted'      ], self.vars['temp_adjusted'      ], \
        self.vars['pres_adjusted_qc'   ], self.vars['psal_adjusted_qc'   ], self.vars['temp_adjusted_qc'   ], \
        self.vars['pres_adjusted_error'], self.vars['psal_adjusted_error'], self.vars['temp_adjusted_error'] = \
        self.readnc(varnm)

        # if not np.ma.is_masked(self.vars['pres_adjusted']):
        #     self.vars['pres_adjusted'] = np.ma.masked_values(self.vars['pres_adjusted'], 99999.)
        # if not np.ma.is_masked(self.vars['psal_adjusted']):
        #     self.vars['psal_adjusted'] = np.ma.masked_values(self.vars['psal_adjusted'], 99999.)
        # if not np.ma.is_masked(self.vars['temp_adjusted']):
        #     self.vars['temp_adjusted'] = np.ma.masked_values(self.vars['temp_adjusted'], 99999.)
        # if not np.ma.is_masked(self.vars['pres_adjusted_error']):
        #     self.vars['pres_adjusted_error'] = np.ma.masked_values(self.vars['pres_adjusted_error'], 99999.)
        # if not np.ma.is_masked(self.vars['psal_adjusted_error']):
        #     self.vars['psal_adjusted_error'] = np.ma.masked_values(self.vars['psal_adjusted_error'], 99999.)
        # if not np.ma.is_masked(self.vars['temp_adjusted_error']):
        #     self.vars['temp_adjusted_error'] = np.ma.masked_values(self.vars['temp_adjusted_error'], 99999.)

        # self.vars['pres_adjusted'].mask[self.vars['pres_adjusted'] == 99999.] = True
        # self.vars['psal_adjusted'].mask[self.vars['psal_adjusted'] == 99999.] = True
        # self.vars['temp_adjusted'].mask[self.vars['temp_adjusted'] == 99999.] = True
        # self.vars['pres_adjusted_error'].mask[self.vars['pres_adjusted_error'] == 99999.] = True
        # self.vars['psal_adjusted_error'].mask[self.vars['psal_adjusted_error'] == 99999.] = True
        # self.vars['temp_adjusted_error'].mask[self.vars['temp_adjusted_error'] == 99999.] = True

    def sort_time(self):
        sort_idx = sorted(range(len(self.coords['time'])), key=self.coords['time'].__getitem__)
        sort_idx = np.array(sort_idx)
        return sort_idx

    def qcsort_vars(self):
        pass

    def qcsort_vars_adjusted(self):
        mask_pos, mask_time = self.qc_coords()
        maks_pres   = self.qc_vars('pres_adjusted')
        maks_psal   = self.qc_vars('psal_adjusted')
        maks_temp   = self.qc_vars('temp_adjusted')
        mask_dm     = self.qc['Data_mode'] != b'D'

        mask = mask_pos | mask_time | maks_pres | maks_psal | maks_temp | mask_dm

        self.vars['pres_adjusted'].mask = mask
        self.vars['psal_adjusted'].mask = mask
        self.vars['temp_adjusted'].mask = mask
        self.vars['pres_adjusted_error'].mask = mask
        self.vars['psal_adjusted_error'].mask = mask
        self.vars['temp_adjusted_error'].mask = mask

        sort_l = self.sort_lev('pres_adjusted')
        sort_t = self.sort_time()

        self.vars['pres_adjusted'] = self.vars['pres_adjusted'][np.arange(self._nprof)[:, None], sort_l][sort_t, :]
        self.vars['psal_adjusted'] = self.vars['psal_adjusted'][np.arange(self._nprof)[:, None], sort_l][sort_t, :]
        self.vars['temp_adjusted'] = self.vars['temp_adjusted'][np.arange(self._nprof)[:, None], sort_l][sort_t, :]

        self.vars['pres_adjusted_error'] = self.vars['pres_adjusted_error'][np.arange(self._nprof)[:, None], sort_l][sort_t, :]
        self.vars['psal_adjusted_error'] = self.vars['psal_adjusted_error'][np.arange(self._nprof)[:, None], sort_l][sort_t, :]
        self.vars['temp_adjusted_error'] = self.vars['temp_adjusted_error'][np.arange(self._nprof)[:, None], sort_l][sort_t, :]

        # self.sort_vars('pres_adjusted', 'pres_adjusted_qc')
        # self.sort_vars('pres_adjusted', 'psal_adjusted_qc')
        # self.sort_vars('pres_adjusted', 'temp_adjusted_qc')

    # def calc_ncep_coords(self):
    #     mask_pos, mask_time = self.qc_coords()
    #     mask_coords = mask_pos | mask_time
    #     ncep_ll = ncep_grid()
    #
    #     self.ncep_cell = dict.fromkeys(['ix', 'iy', 'it'])
    #     self.ncep_cell['ix'] = np.zeros((self._nprof, 1), dtype = np.int64)
    #     self.ncep_cell['iy'] = np.zeros((self._nprof, 1), dtype = np.int64)
    #     self.ncep_cell['it'] = np.zeros((self._nprof, 2), dtype = np.int64)
    #
    #     for iprof in range(self._nprof):
    #         if mask_coords[iprof]: continue
    #
    #         ix = np.argmax(ncep_ll.lonb > self.coords['lon'][iprof]) - 1
    #         if   np.floor(self.coords['lat'][iprof]) >= ncep_ll.latb[0]:
    #             iy = 0
    #         elif np.floor(self.coords['lat'][iprof]) <= ncep_ll.latb[-1]:
    #             iy = ncep_ll.lat.size - 1
    #         else:
    #             iy = np.argmax(ncep_ll.latb < self.coords['lat'][iprof]) - 1
    #
    #         dt  = self.coords['time'][iprof]
    #         it  = int((dt - datetime.datetime(dt.year, 1, 1)).total_seconds()/3600/6)
    #
    #         self.ncep_cell['ix'][iprof] = int(ix)
    #         self.ncep_cell['iy'][iprof] = int(iy)
    #         self.ncep_cell['it'][iprof, 0] = int(dt.year)
    #         self.ncep_cell['it'][iprof, 1] = it

    def calc_mld_HT09(self, matlabtoolspath = '/Users/hewang/research/argo/matlabtools/', varsuffix = '_adjusted'):
        eng = matlab.engine.start_matlab()
        eng.addpath(matlabtoolspath)
        self.mld = {'Mixedtp': np.zeros((self._nprof, )), 'Mixedt_ta': np.zeros((self._nprof, )), 'Mixedt_sa': np.zeros((self._nprof, )), 'Mixedt_da': np.zeros((self._nprof, )), \
                    'Mixedsp': np.zeros((self._nprof, )), \
                    'Mixeddp': np.zeros((self._nprof, )), 'Mixedd_ta': np.zeros((self._nprof, )), 'Mixedd_sa': np.zeros((self._nprof, )), 'Mixedd_da': np.zeros((self._nprof, )), \
                    'Mldepthptmpp': np.zeros((self._nprof, )), 'Mldepthptmp_ta': np.zeros((self._nprof, )), 'Mldepthptmp_sa': np.zeros((self._nprof, )), 'Mldepthptmp_da': np.zeros((self._nprof, )), \
                    'Mldepthdensp': np.zeros((self._nprof, )), 'Mldepthdens_ta': np.zeros((self._nprof, )), 'Mldepthdens_sa': np.zeros((self._nprof, )), 'Mldepthdens_da': np.zeros((self._nprof, )), \
                    'Gtmldp': np.zeros((self._nprof, )), 'Gdmldp': np.zeros((self._nprof, )), \
                    'Tanalysis': np.zeros((self._nprof, )), 'Sanalysis': np.zeros((self._nprof, )), 'Danalysis': np.zeros((self._nprof, ))}
        for ipf in range(self._nprof):
            if self.vars['pres'+varsuffix][ipf, :].compressed().shape[0] <= 1: continue

            pres = self.vars['pres'+varsuffix][ipf, :].compressed()
            sal  = self.vars['psal'+varsuffix][ipf, :].compressed()
            temp = self.vars['temp'+varsuffix][ipf, :].compressed()
            pden = self.vars['rho0'][ipf, :].compressed()

            pres = matlab.double(pres.tolist())
            sal  = matlab.double(sal.tolist())
            temp = matlab.double(temp.tolist())
            pden = matlab.double(pden.tolist())

            self.mld['Mixedtp'][ipf,], self.mld['Mixedt_ta'][ipf,], self.mld['Mixedt_sa'][ipf,], self.mld['Mixedt_da'][ipf,], \
            self.mld['Mixedsp'][ipf,], \
            self.mld['Mixeddp'][ipf,], self.mld['Mixedd_ta'][ipf,], self.mld['Mixedd_sa'][ipf,], self.mld['Mixedd_da'][ipf,], \
            self.mld['Mldepthptmpp'][ipf,], self.mld['Mldepthptmp_ta'][ipf,], self.mld['Mldepthptmp_sa'][ipf,], self.mld['Mldepthptmp_da'][ipf,], \
            self.mld['Mldepthdensp'][ipf,], self.mld['Mldepthdens_ta'][ipf,], self.mld['Mldepthdens_sa'][ipf,], self.mld['Mldepthdens_da'][ipf,], \
            self.mld['Gtmldp'][ipf,], self.mld['Gdmldp'][ipf,], \
            self.mld['Tanalysis'][ipf,], self.mld['Sanalysis'][ipf,], self.mld['Danalysis'][ipf,] = eng.calc_mld_prof(pres, sal, temp, pden, nargout=22)

            eng.quit()

    def calc_ncep_idx(self):
        mask_pos, mask_time = self.qc_coords()
        mask_coords = mask_pos | mask_time
        ncep_ll = ncep_grid()

        self.ncep_cell = dict.fromkeys(['ix', 'iy', 'it'])
        self.ncep_cell['ix'], self.ncep_cell['iy'] = self.bin_idx(ncep_ll.lonb, ncep_ll.latb)
        self.ncep_cell['it'] = np.zeros((self._nprof, 2), dtype = np.int64)

        for iprof in range(self._nprof):
            if mask_coords[iprof]: continue

            dt  = self.coords['time'][iprof]
            it  = int((dt - datetime.datetime(dt.year, 1, 1)).total_seconds()/3600/6)

            self.ncep_cell['it'][iprof, 0] = int(dt.year)
            self.ncep_cell['it'][iprof, 1] = it

    def calc_ncep_sfc(self, ncepvars = ['lhtfl.sfc', 'nlwrs.sfc', 'nswrs.sfc', 'shtfl.sfc', 'uwnd.10m', 'vwnd.10m'], datapath = None, ave = None):
        varnm = [var.split('.')[0] for var in ncepvars]
        try:
            self.ncep_sfc
        except AttributeError:
            self.ncep_sfc = {}
        for key in varnm:
            self.ncep_sfc[key] = np.zeros((self._nprof, 1))

        for iprof in range(self._nprof):
            ix = self.ncep_cell['ix'][iprof]
            iy = self.ncep_cell['iy'][iprof]
            it = self.ncep_cell['it'][iprof, :]
            if ix == 0 and iy == 0 and np.all(it==0): continue

            if ave == 'day':
                itst = it[1] - it[1]%4
                ited = it[1] - it[1]%4 + 4
            else:
                itst = it[1]
                ited = it[1] + 1
            for vn_f, vn in zip(ncepvars, varnm):
                fn = datapath + vn_f + '.gauss.' + "{:d}".format(it[0]) + '.nc'
                ncgrp = ncds(fn)
                self.ncep_sfc[vn][iprof] = ncgrp.variables[vn][itst:ited, iy, ix].mean()
                ncgrp.close()

    def calc_bulkTA_mld(self, varsuffix = '_adjusted', mld = 'Mixedtp', mlt = 'Mixedd_ta', mls = 'Mixedd_sa', dep = 50.):
        presnm = 'pres' + varsuffix
        psalnm = 'psal' + varsuffix
        tempnm = 'temp' + varsuffix

        self.bulkTA_mld = dict.fromkeys(['adTdz', 'bdSdz', 'TA'])
        self.bulkTA_mld['adTdz'] = np.zeros((self._nprof, 1))
        self.bulkTA_mld['bdSdz'] = np.zeros((self._nprof, 1))
        self.bulkTA_mld['TA'   ] = np.zeros((self._nprof, 1))

        for iprof in range(self._nprof):
            mask_sec = (self.vars[presnm][iprof, :] > self.mld['Mixedtp'][iprof] & \
                        self.vars[presnm][iprof, :] < self.mld['Mixedtp'][iprof] + dep)
            mask = mask_sec.filled(fill_value=False)
            if mask is None or np.all(~mask): continue

            t = np.array([self.mld[mlt][iprof]  , self.vars[tempnm][iprof, mask].mean()])
            s = np.array([self.mld[mls][iprof]  , self.vars[psalnm][iprof, mask].mean()])
            p = np.array([self.mld[mld][iprof]/2, self.mld[mld][iprof] + dep/2])
            if np.any(np.isnan(t) | np.isnan(s) | np.isnan(p)):continue

            Tu, Rsubrho, adTdz, bdSdz, _ = Argofloat.sw_turner(s, t, p)
            self.bulkTA_mld['adTdz'][iprof] = adTdz
            self.bulkTA_mld['bdSdz'][iprof] = bdSdz
            self.bulkTA_mld['TA'   ][iprof] = Tu

# class Argofloat(Argofloat):

# class Argofloat_GSW(Argofloat):
#     def __init__(self, datapath = '', Float_ID = '', DAC = None):
#         Argofloat.__init__(self, datapath = datapath, Float_ID = Float_ID, DAC = DAC):

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

    @staticmethod
    def sw_turner(s, t, p):
        smid = (s[:-1] + s[1:]) / 2
        tmid = (t[:-1] + t[1:]) / 2
        pmid = (p[:-1] + p[1:]) / 2
        dz = -(p[:-1] - p[1:])

        tu = sw.eos80.ptmp(s[:-1], t[:-1], p[:-1], pr=pmid)
        td = sw.eos80.ptmp(s[1:] , t[1:] , p[1:] , pr=pmid)

        ds = s[:-1] - s[1:]
        dt = tu - td

        alpha = sw.eos80.alpha(smid, tmid, pmid, pt=False)
        beta  = sw.eos80.beta (smid, tmid, pmid, pt=False)

        Tu = np.arctan2(alpha * dt + beta * ds, alpha * dt - beta * ds) * 180 / np.pi
        Rsubrho = (alpha * dt) / (beta * ds)
        adTdz = alpha * dt / dz
        bdSdz = beta * ds / dz

        return Tu, Rsubrho, adTdz, bdSdz, pmid

    @staticmethod
    def sw_mld(s, t, p, pd_c):
        rho0 = sw.dens0(s, t)
        rho0_10m = np.interp(10., p, rho0)

        depm = p[-1]
        for iz in range(1, p.size, 1):
            if p[iz] > 10. and rho0[iz] - rho0_10m > pd_c:
                depm = np.interp(rho0_10m + pd_c, np.array([rho0[iz], rho0[iz-1]]), np.array([p[iz], p[iz-1]]))
                break

        mld = depm
        return mld

    def calc_mld(self):
        self.calc_mld_dt()

    def calc_TA(self):
        self.calc_TA_sw()

    def calc_n2(self):
        self.calc_n2_sw()

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
        self.vars['pv'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
        self.vars['pmid'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])

        for iprof in range(self._nprof):
            s = self.vars['psal_adjusted'][iprof,:].filled(np.nan)
            t = self.vars['temp_adjusted'][iprof,:].filled(np.nan)
            p = self.vars['pres_adjusted'][iprof,:].filled(np.nan)

            n2, q, pmid = sw.geostrophic.bfrq(s, t, p, lat=self.vars['lat'][iprof])
            self.vars['N2'][iprof,:] = np.squeeze(n2)
            self.vars['pv'][iprof,:] = np.squeeze(q)
            self.vars['pmid'][iprof,:] = np.squeeze(pmid)

    def calc_mld_dt(self):
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
            ta[np.isnan(ta)] = 99999.0; ta = np.ma.masked_values(ta, 99999.0)
            rsubrho[np.isnan(rsubrho)] = 99999.0; rsubrho = np.ma.masked_values(rsubrho, 99999.0)
            dtdz[np.isnan(dtdz)] = 99999.0; dtdz = np.ma.masked_values(dtdz, 99999.0)
            dsdz[np.isnan(dsdz)] = 99999.0; dsdz = np.ma.masked_values(dsdz, 99999.0)
            pmid[np.isnan(pmid)] = 99999.0; pmid = np.ma.masked_values(pmid, 99999.0)

            self.vars['TA'][iprof,:] = np.squeeze(ta)
            self.vars['Rsubrho'][iprof,:] = np.squeeze(rsubrho)
            self.vars['adtdz'][iprof,:] = np.squeeze(dtdz)
            self.vars['bdsdz'][iprof,:] = np.squeeze(dsdz)
            self.vars['pmid'][iprof,:] = np.squeeze(pmid)

    def plot_track(self, omitqc_list = ['pos', 'time', 'mode']):
        maskall = {}
        maskall['pos'], maskall['time'] = self.qc_coords()
        maskall['mode'] = self.qc['Data_mode'] != b'D'

        mask = np.full_like(maskall['mode'], False)
        for var in omitqc_list:
            try:
                mask = mask | maskall[var]
            except AttributeError:
                print('Invalid omited field: ' + var)

        lon = self.coords['lon'][~mask]
        lat = self.coords['lat'][~mask]
        proj = ccrs.Mercator()

        xticks = list(range(35, 85, 10))
        yticks = list(range(-10, 40, 10))
        fig, ax = plt.subplots(subplot_kw={'projection': proj})
        ax.set_extent((38,80,-10,31), crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, edgecolor='k', linewidth=0.2)
        ax.gridlines(linewidth=0.2, xlocs = np.arange(35, 95, 10), ylocs = np.arange(-10, 50, 10))
        ax.set_xticks(xticks, crs=ccrs.PlateCarree())
        ax.set_yticks(yticks, crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.plot(lon, lat, transform=ccrs.PlateCarree())
        ax.plot(lon[ 0], lat[ 0], 'r*', transform=ccrs.PlateCarree())
        ax.plot(lon[-1], lat[-1], 'k*', transform=ccrs.PlateCarree())

        ax.plot(72, 29, 'r*', transform=ccrs.PlateCarree())
        ax.plot(72, 27, 'k*', transform=ccrs.PlateCarree())
        ax.text(74, 29, 'start', transform=ccrs.PlateCarree(), va='center')
        ax.text(74, 27, 'end', transform=ccrs.PlateCarree(), va='center')

        return fig, ax

    def plot_cs(self, varnm, lvs = None, cmap = 'viridis', depmax = 1000):
        mask_dm = self.qc['Data_mode'] != b'D'
        mask_time = self.coords['time_qc'] > b'1'
        mask_pos  = self.coords['pos_qc'] > b'1'
        mask = mask_time | mask_dm
        maskb = np.r_[mask[0]+mask[1], mask[:-1]+mask[1:], mask[-2]+mask[-1]]

        if varnm in ['TA', 'N2', 'rho', 'adtdz', 'bdsdz', 'Rsubrho']:
            presb = self.vars['pres_adjusted']
        else:
            presb = np.c_[self.vars['pres_adjusted'][:, 0] - (self.vars['pres_adjusted'][:, 1]-self.vars['pres_adjusted'][:, 0])/2, (self.vars['pres_adjusted'][:, :-1]+self.vars['pres_adjusted'][:, 1:])/2, self.vars['pres_adjusted'][:, -1] + (self.vars['pres_adjusted'][:, -1] - self.vars['pres_adjusted'][:, -2])/2]

        # presb = presb.filled(np.nan)
        presb = presb.filled(99999.)
        # timeb = [calendar.timegm(self.coords['time'][ip].timetuple())-0.5 for ip in range(self._nprof)]
        # timeb.append(calendar.timegm(self.coords['time'][-1].timetuple())+0.5)
        time = mpl.dates.date2num(self.coords['time'])
        timeb = np.r_[time[0] - (time[1] - time[0])/2, (time[:-1]+time[1:])/2, time[-1] + (time[-1] - time[-2])/2]

        years = mpl.dates.YearLocator()   # every year
        months = mpl.dates.MonthLocator()  # every month
        timeFmt = mpl.dates.DateFormatter('%Y-%m')

        norm = mpl.colors.BoundaryNorm(lvs, 256)

        var = self.vars[varnm].filled(np.nan)

        fig, ax = plt.subplots()
        for ip in range(self._nprof):
            pmsh = ax.pcolormesh(timeb[ip:ip+2], presb[ip, :], var[ip, :][:, None], norm = norm, cmap = cmap)
        # ax.set_xlim(timeb[0], timeb[-1])
        ax.set_xlim(timeb[~maskb][0], timeb[~maskb][-1])
        ax.set_ylim(0, depmax)
        ax.invert_yaxis()

        ax.xaxis.set_major_locator(months)
        # ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_major_formatter(timeFmt)

        fig.autofmt_xdate()
        fig.colorbar(pmsh, ticks=lvs[::5]);
        return fig, ax

# class floatdict(Argofloat):
#     def __init__(self):
#         allfloats = {}
# iflt = 0
# exc_flt = []
# # exc_flt = ['2901857', '2901856', '2901855', '2901846', '2901847', '2901852', '2901854', '2901853', '2901859', '2901860']
# for idac in Argofloat.listofdac:
#     if os.path.isdir(datapath + '/' + idac):
#         for fltid in os.listdir(datapath + '/' + idac):
#             print(fltid)
#             if fltid in exc_flt:
#
#         Argofloat.__init__(self, datapath = '', Float_ID = '', DAC = None)
    # def

# class Argofloat:
#     listofdac = ['aoml', 'bodc', 'coriolis', 'csio', 'csiro', 'incois', 'jma', 'kma', 'kordi', 'meds', 'nmdis']
#     qc_accept_flag = 1
#
#     def __init__(self, datapath = '', Float_ID = '', DAC = None):
#         self.datapath = datapath
#         self.Float_ID = Float_ID
#         self.DAC = DAC
#         self._getfn()
#
#         self.coords = {'lon': None, 'lat': None, 'pos_qc': None, 'time': None, 'time_qc': None, 'cycle': None, 'direction': None}
#         self.qc = {'Data_mode': None, 'Data_state_indicator': None, 'Profile_pres_qc': None, 'Profile_psal_qc': None, 'Profile_temp_qc': None}
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
#         if not self._DAC in Argofloat.listofdac:
#             for idac in Argofloat.listofdac:
#                 if os.path.isdir(self.datapath + '/' + idac + '/' + self.Float_ID):
#                     self._DAC = idac
#                     break
#
#     def _getfn(self):
#         fndir = self.datapath + '/' + self.DAC + '/' + self.Float_ID
#         fn_prof = fndir + '/' + self.Float_ID + '_prof.nc'
#         fn_ind  = glob.glob(fndir + '/profiles' + '/*.nc')
#         if os.path.isfile(fn_prof):
#             self._filename = [fn_prof]
#             self._nprof = ncds(self._filename[0]).dimensions['N_PROF'].size
#             self._nfile = 1
#         elif len(fn_ind) > 0:
#             fn_ind_rt = glob.glob(fndir + '/profiles/' + 'R*.nc')
#             for ifl in fn_ind_rt:
#                 if ifl.replace('/R', '/D') in fn_ind:
#                     del fn_ind[fn_ind.index(ifl)]
#             self._filename = fn_ind
#             self._nprof = len(fn_ind)
#             self._nfile = len(fn_ind)
#         else:
#             self._filename = []
#             self._nprof = 0
#             self._nfile = 0
#             print('Profile file not found! \n' + fn_prof + '\n' + fn_ind)
#         return
#
#     def readnc(self, varnm):
#         varlist = [[] for _ in varnm]
#         vardim = [[] for _ in varnm]
#         nlv = 0
#
#         for fn, iprof in zip(self._filename, range(self._nprof)):
#             ncgrp = ncds(fn)
#
#             nlv = max(nlv, ncgrp.dimensions['N_LEVELS'].size)
#
#             for iv, ivnm in zip(range(len(varnm)), varnm):
#                 varlist[iv].append(ncgrp.variables[ivnm][:])
#                 vardim [iv].append(ncgrp.variables[ivnm].dimensions)
#             ncgrp.close()
#
#         var = []
#         if len(self._filename) == 1:
#             for iv in range(len(varnm)):
#                 var.append(varlist[iv][0])
#         else:
#             for iv in range(len(varnm)):
#                 if 'N_LEVELS' in vardim[iv][0]:
#                     var.append(np.ma.masked_values(np.zeros((self._nprof, nlv))+99999., 99999.))
#                     for iprof in range(self._nprof):
#                         var[iv][iprof, :varlist[iv][iprof].size] = varlist[iv][iprof]
#                 else:
#                     var.append(np.stack(varlist[iv]))
#         return var
#
#     def load_coords(self):
#         varnm = ['CYCLE_NUMBER', 'DIRECTION', 'POSITION_QC', 'LONGITUDE', 'LATITUDE', 'JULD_QC', 'REFERENCE_DATE_TIME', 'JULD']
#         self.coords['cycle'], self.coords['direction'], \
#         self.coords['pos_qc'], self.coords['lon'], self.coords['lat'], \
#         self.coords['time_qc'], reftime, juld = self.readnc(varnm)
#
#         # self.coords['time'] = np.zeros((self._nprof,), dtype = 'datetime64')
#         self.coords['time'] = []
#         for iprof in range(self._nprof):
#             if self.coords['time_qc'][iprof] > b'1':
#                 self.coords['time'].append(None)
#                 continue
#             if len(reftime.shape) == 1:
#                 self.coords['time'].append(datetime.datetime.strptime(reftime.tobytes().decode(), '%Y%m%d%H%M%S') + datetime.timedelta(days = juld[iprof]))
#             else:
#                 self.coords['time'].append(datetime.datetime.strptime(reftime[iprof].tobytes().decode(), '%Y%m%d%H%M%S') + datetime.timedelta(days = juld[iprof][0]))
#
#     def load_qc(self):
#         varnm = ['DATA_MODE', 'DATA_STATE_INDICATOR', 'PROFILE_PRES_QC', 'PROFILE_PSAL_QC', 'PROFILE_TEMP_QC']
#         self.qc['Data_mode'], self.qc['Data_state_indicator'], self.qc['Profile_pres_qc'], self.qc['Profile_psal_qc'], self.qc['Profile_temp_qc'] = self.readnc(varnm)
#
#     def load_vars(self):
#         varnm = ['PRES', 'PSAL', 'TEMP', 'PRES_QC', 'PSAL_QC', 'TEMP_QC']
#         self.vars['pres'], self.vars['psal'], self.vars['temp'], self.vars['pres_qc'], self.vars['psal_qc'], self.vars['temp_qc'] = self.readnc(varnm)
#
#     def load_vars_adjusted(self):
#         varnm = ['PRES_ADJUSTED', 'PSAL_ADJUSTED', 'TEMP_ADJUSTED', 'PRES_ADJUSTED_QC', 'PSAL_ADJUSTED_QC', 'TEMP_ADJUSTED_QC',
#                  'PRES_ADJUSTED_ERROR', 'PSAL_ADJUSTED_ERROR', 'TEMP_ADJUSTED_ERROR']
#
#         self.vars['pres_adjusted'], self.vars['psal_adjusted'], self.vars['temp_adjusted'], \
#         self.vars['pres_adjusted_qc'], self.vars['psal_adjusted_qc'], self.vars['temp_adjusted_qc'], \
#         self.vars['pres_adjusted_error'], self.vars['psal_adjusted_error'], self.vars['temp_adjusted_error'] = \
#         self.readnc(varnm)
#
#         if not np.ma.is_masked(self.vars['pres_adjusted']):
#             self.vars['pres_adjusted'] = np.ma.masked_values(self.vars['pres_adjusted'], 99999.)
#         if not np.ma.is_masked(self.vars['psal_adjusted']):
#             self.vars['psal_adjusted'] = np.ma.masked_values(self.vars['psal_adjusted'], 99999.)
#         if not np.ma.is_masked(self.vars['temp_adjusted']):
#             self.vars['temp_adjusted'] = np.ma.masked_values(self.vars['temp_adjusted'], 99999.)
#         if not np.ma.is_masked(self.vars['pres_adjusted_error']):
#             self.vars['pres_adjusted_error'] = np.ma.masked_values(self.vars['pres_adjusted_error'], 99999.)
#         if not np.ma.is_masked(self.vars['psal_adjusted_error']):
#             self.vars['psal_adjusted_error'] = np.ma.masked_values(self.vars['psal_adjusted_error'], 99999.)
#         if not np.ma.is_masked(self.vars['temp_adjusted_error']):
#             self.vars['temp_adjusted_error'] = np.ma.masked_values(self.vars['temp_adjusted_error'], 99999.)
#
#         self.vars['pres_adjusted'].mask[self.vars['pres_adjusted'] == 99999.] = True
#         self.vars['psal_adjusted'].mask[self.vars['psal_adjusted'] == 99999.] = True
#         self.vars['temp_adjusted'].mask[self.vars['temp_adjusted'] == 99999.] = True
#         self.vars['pres_adjusted_error'].mask[self.vars['pres_adjusted_error'] == 99999.] = True
#         self.vars['psal_adjusted_error'].mask[self.vars['psal_adjusted_error'] == 99999.] = True
#         self.vars['temp_adjusted_error'].mask[self.vars['temp_adjusted_error'] == 99999.] = True
#
#     def qc_coords(self):
#         pass
#
#     def qc_vars_adjusted(self):
#         if np.ma.is_masked(self.vars['pres_adjusted_qc']):
#             self.vars['pres_adjusted_qc'].mask[(self.vars['pres_adjusted_qc'] < b'1') | (self.vars['pres_adjusted_qc'] > b'9')] = True
#             mask_p = (self.vars['pres_adjusted_qc'].filled(b'9').astype(np.int) > Argofloat.qc_accept_flag) | (self.vars['pres_adjusted'].mask)
#         else:
#             self.vars['pres_adjusted_qc'][(self.vars['pres_adjusted_qc'] < b'1') | (self.vars['pres_adjusted_qc'] > b'9')] = b'9'
#             mask_p = (self.vars['pres_adjusted_qc'].astype(np.int) > Argofloat.qc_accept_flag) | (self.vars['pres_adjusted'].mask)
#         if np.ma.is_masked(self.vars['psal_adjusted_qc']):
#             self.vars['psal_adjusted_qc'].mask[(self.vars['psal_adjusted_qc'] < b'1') | (self.vars['psal_adjusted_qc'] > b'9')] = True
#             mask_s = (self.vars['psal_adjusted_qc'].filled(b'9').astype(np.int) > Argofloat.qc_accept_flag) | (self.vars['psal_adjusted'].mask)
#         else:
#             self.vars['psal_adjusted_qc'][(self.vars['psal_adjusted_qc'] < b'1') | (self.vars['psal_adjusted_qc'] > b'9')] = b'9'
#             mask_s = (self.vars['psal_adjusted_qc'].astype(np.int) > Argofloat.qc_accept_flag) | (self.vars['psal_adjusted'].mask)
#         if np.ma.is_masked(self.vars['temp_adjusted_qc']):
#             self.vars['temp_adjusted_qc'].mask[(self.vars['temp_adjusted_qc'] < b'1') | (self.vars['temp_adjusted_qc'] > b'9')] = True
#             mask_t = (self.vars['temp_adjusted_qc'].filled(b'9').astype(np.int) > Argofloat.qc_accept_flag) | (self.vars['temp_adjusted'].mask)
#         else:
#             self.vars['temp_adjusted_qc'][(self.vars['temp_adjusted_qc'] < b'1') | (self.vars['temp_adjusted_qc'] > b'9')] = b'9'
#             mask_t = (self.vars['temp_adjusted_qc'].astype(np.int) > Argofloat.qc_accept_flag) | (self.vars['temp_adjusted'].mask)
#         mask_dm = self.qc['Data_mode'] != b'D'
#         mask_time = self.coords['time_qc'] > b'1'
#         mask_pos  = self.coords['pos_qc'] > b'1'
#         mask = mask_t | mask_s | mask_p | mask_dm[:, None] | mask_time[:, None] | mask_pos[:, None]
#
#         self.vars['pres_adjusted'].mask = mask
#         self.vars['psal_adjusted'].mask = mask
#         self.vars['temp_adjusted'].mask = mask
#         self.vars['pres_adjusted_error'].mask = mask
#         self.vars['psal_adjusted_error'].mask = mask
#         self.vars['temp_adjusted_error'].mask = mask
#
#         # sort along depth + moving all the invalid to the bottom
#         sort_idx = self.vars['pres_adjusted'].filled().argsort(axis = -1) # filled makes the masked element to maked_value, in contrast with data (which uses the true value of the element, even if it is masked)
#         sort_idx_c = np.arange(self._nprof)
#
#         self.vars['pres_adjusted'] = self.vars['pres_adjusted'][sort_idx_c[:, None], sort_idx]
#         self.vars['psal_adjusted'] = self.vars['psal_adjusted'][sort_idx_c[:, None], sort_idx]
#         self.vars['temp_adjusted'] = self.vars['temp_adjusted'][sort_idx_c[:, None], sort_idx]
#         self.vars['pres_adjusted_error'] = self.vars['pres_adjusted_error'][sort_idx_c[:, None], sort_idx]
#         self.vars['psal_adjusted_error'] = self.vars['psal_adjusted_error'][sort_idx_c[:, None], sort_idx]
#         self.vars['temp_adjusted_error'] = self.vars['temp_adjusted_error'][sort_idx_c[:, None], sort_idx]
#         self.vars['pres_adjusted_qc'] = self.vars['pres_adjusted_qc'][sort_idx_c[:, None], sort_idx]
#         self.vars['psal_adjusted_qc'] = self.vars['psal_adjusted_qc'][sort_idx_c[:, None], sort_idx]
#         self.vars['temp_adjusted_qc'] = self.vars['temp_adjusted_qc'][sort_idx_c[:, None], sort_idx]
#
#     @staticmethod
#     def sw_turner(s, t, p):
#         smid = (s[:-1] + s[1:]) / 2
#         tmid = (t[:-1] + t[1:]) / 2
#         pmid = (p[:-1] + p[1:]) / 2
#         dz = -(p[:-1] - p[1:])
#
#         tu = sw.eos80.ptmp(s[:-1], t[:-1], p[:-1], pr=pmid)
#         td = sw.eos80.ptmp(s[1:] , t[1:] , p[1:] , pr=pmid)
#
#         ds = s[:-1] - s[1:]
#         dt = tu - td
#
#         alpha = sw.eos80.alpha(smid, tmid, pmid, pt=False)
#         beta  = sw.eos80.beta (smid, tmid, pmid, pt=False)
#
#         Tu = np.arctan2(alpha * dt + beta * ds, alpha * dt - beta * ds) * 180 / np.pi
#         Rsubrho = (alpha * dt) / (beta * ds)
#         adTdz = alpha * dt / dz
#         bdSdz = beta * ds / dz
#
#         return Tu, Rsubrho, adTdz, bdSdz, pmid
#
#     @staticmethod
#     def sw_mld(s, t, p, pd_c):
#         rho0 = sw.dens0(s, t)
#         rho0_10m = np.interp(10., p, rho0)
#
#         depm = p[-1]
#         for iz in range(1, p.size, 1):
#             if p[iz] > 10. and rho0[iz] - rho0_10m > pd_c:
#                 depm = np.interp(rho0_10m + pd_c, np.array([rho0[iz], rho0[iz-1]]), np.array([p[iz], p[iz-1]]))
#                 break
#
#         mld = depm
#         return mld
#
#     def calc_mld(self):
#         self.calc_mld_dt()
#
#     def calc_TA(self):
#         self.calc_TA_sw()
#
#     def calc_rho0(self):
#         self.calc_rho0_sw()
#
#     def calc_n2(self):
#         self.calc_n2_sw()
#
#     def calc_SACT(self):
#         self.vars['SA'] = np.zeros_like(self.vars['psal_adjusted'])
#         self.vars['CT'] = np.zeros_like(self.vars['psal_adjusted'])
#
#         for iprof in range(self._nprof):
#             s = self.vars['psal_adjusted'][iprof, :]
#             p = self.vars['pres_adjusted'][iprof, :]
#             t = self.vars['temp_adjusted'][iprof, :]
#             lon = self.coords['lon'][iprof]
#             lat = self.coords['lat'][iprof]
#
#             self.vars['SA'][iprof, :] = gsw.SA_from_SP(s, p, lon, lat)
#             self.vars['CT'][iprof, :] = gsw.CT_from_t(self.vars['SA'][iprof, :], t, p)
#
#     def calc_n2_gsw(self):
#         if self.vars['SA'] is None or self.vars['CT'] is None:
#             self.calc_SACT()
#         self.vars['N2'], self.vars['pmid'] = gsw.stability.Nsquared(self.vars['SA'], self.vars['CT'], self.vars['pres_adjusted'], axis = 1)
#
#     def calc_TA_gsw(self):
#         if self.vars['SA'] is None or self.vars['CT'] is None:
#             self.calc_SACT()
#         self.vars['TA'], self.vars['Rsubrho'], self.vars['pmid'] = gsw.stability.Turner_Rsubrho(self.vars['SA'], self.vars['CT'], self.vars['pres_adjusted'], axis = 1)
#
#     def calc_n2_sw(self):
#         self.vars['N2'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#         self.vars['pv'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#         self.vars['pmid'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#
#         for iprof in range(self._nprof):
#             s = self.vars['psal_adjusted'][iprof,:].filled(np.nan)
#             t = self.vars['temp_adjusted'][iprof,:].filled(np.nan)
#             p = self.vars['pres_adjusted'][iprof,:].filled(np.nan)
#
#             n2, q, pmid = sw.geostrophic.bfrq(s, t, p, lat=self.vars['lat'][iprof])
#             self.vars['N2'][iprof,:] = np.squeeze(n2)
#             self.vars['pv'][iprof,:] = np.squeeze(q)
#             self.vars['pmid'][iprof,:] = np.squeeze(pmid)
#
#     def calc_mld_dt(self):
#         pd_c = 0.03
#         self.vars['mld'] = np.zeros((self._nprof, ))
#
#         for iprof in range(self._nprof):
#             s = self.vars['psal_adjusted'][iprof,:].filled(np.nan)
#             t = self.vars['temp_adjusted'][iprof,:].filled(np.nan)
#             p = self.vars['pres_adjusted'][iprof,:].filled(np.nan)
#
#             if np.all(np.isnan(s)) or np.all(np.isnan(t)) or np.all(np.isnan(p)):
#                 self.vars['mld'][iprof] = np.nan
#             else:
#                 self.vars['mld'][iprof] = Argofloat.sw_mld(s, t, p, pd_c)
#
#     def calc_TA_sw(self):
#         self.vars['TA'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#         self.vars['Rsubrho'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#         self.vars['adtdz'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#         self.vars['bdsdz'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#         self.vars['pmid'] = np.zeros_like(self.vars['psal_adjusted'][:, :-1]+self.vars['psal_adjusted'][:, 1:])
#
#         for iprof in range(self._nprof):
#             s = self.vars['psal_adjusted'][iprof,:].filled(np.nan)
#             t = self.vars['temp_adjusted'][iprof,:].filled(np.nan)
#             p = self.vars['pres_adjusted'][iprof,:].filled(np.nan)
#             # s = self.vars['psal_adjusted'][iprof,:].compressed()
#             # t = self.vars['temp_adjusted'][iprof,:].compressed()
#             # p = self.vars['pres_adjusted'][iprof,:].compressed()
#
#             # ta, rsubrho, pmid = Argofloat.sw_turner(s, t, p)
#             ta, rsubrho, dtdz, dsdz, pmid = Argofloat.sw_turner(s, t, p)
#             ta[np.isnan(ta)] = 99999.0; ta = np.ma.masked_values(ta, 99999.0)
#             rsubrho[np.isnan(rsubrho)] = 99999.0; rsubrho = np.ma.masked_values(rsubrho, 99999.0)
#             dtdz[np.isnan(dtdz)] = 99999.0; dtdz = np.ma.masked_values(dtdz, 99999.0)
#             dsdz[np.isnan(dsdz)] = 99999.0; dsdz = np.ma.masked_values(dsdz, 99999.0)
#             pmid[np.isnan(pmid)] = 99999.0; pmid = np.ma.masked_values(pmid, 99999.0)
#
#             self.vars['TA'][iprof,:] = np.squeeze(ta)
#             self.vars['Rsubrho'][iprof,:] = np.squeeze(rsubrho)
#             self.vars['adtdz'][iprof,:] = np.squeeze(dtdz)
#             self.vars['bdsdz'][iprof,:] = np.squeeze(dsdz)
#             self.vars['pmid'][iprof,:] = np.squeeze(pmid)
#
#     def calc_rho0_sw(self):
#         s = self.vars['psal_adjusted'].filled(np.nan)
#         t = self.vars['temp_adjusted'].filled(np.nan)
#         r = sw.dens0(s, t)
#         r[np.isnan(r)] = 99999.0
#         self.vars['rho0'] = np.ma.masked_values(r, 99999.0)
#
#     def plot_track(self):
#         mask_dm = self.qc['Data_mode'] != b'D'
#         mask_time = self.coords['time_qc'] > b'1'
#         mask_pos  = self.coords['pos_qc'] > b'1'
#         mask = mask_pos | mask_dm
#
#         lon = self.coords['lon'][~mask]
#         lat = self.coords['lat'][~mask]
#         proj = ccrs.Mercator()
#
#         xticks = list(range(35, 85, 10))
#         yticks = list(range(-10, 40, 10))
#         fig, ax = plt.subplots(subplot_kw={'projection': proj})
#         ax.set_extent((38,80,-10,31), crs=ccrs.PlateCarree())
#         ax.add_feature(cfeature.LAND, edgecolor='k', linewidth=0.2)
#         ax.gridlines(linewidth=0.2, xlocs = np.arange(35, 95, 10), ylocs = np.arange(-10, 50, 10))
#         ax.set_xticks(xticks, crs=ccrs.PlateCarree())
#         ax.set_yticks(yticks, crs=ccrs.PlateCarree())
#         ax.xaxis.set_major_formatter(LongitudeFormatter())
#         ax.yaxis.set_major_formatter(LatitudeFormatter())
#         ax.plot(lon, lat, transform=ccrs.PlateCarree())
#         ax.plot(lon[ 0], lat[ 0], 'r*', transform=ccrs.PlateCarree())
#         ax.plot(lon[-1], lat[-1], 'k*', transform=ccrs.PlateCarree())
#
#         ax.plot(72, 29, 'r*', transform=ccrs.PlateCarree())
#         ax.plot(72, 27, 'k*', transform=ccrs.PlateCarree())
#         ax.text(74, 29, 'start', transform=ccrs.PlateCarree(), va='center')
#         ax.text(74, 27, 'end', transform=ccrs.PlateCarree(), va='center')
#
#         return fig, ax
#
#     def plot_cs(self, varnm, lvs = None, cmap = 'viridis', depmax = 1000):
#         mask_dm = self.qc['Data_mode'] != b'D'
#         mask_time = self.coords['time_qc'] > b'1'
#         mask_pos  = self.coords['pos_qc'] > b'1'
#         mask = mask_time | mask_dm
#         maskb = np.r_[mask[0]+mask[1], mask[:-1]+mask[1:], mask[-2]+mask[-1]]
#
#         if varnm in ['TA', 'N2', 'rho', 'adtdz', 'bdsdz', 'Rsubrho']:
#             presb = self.vars['pres_adjusted']
#         else:
#             presb = np.c_[self.vars['pres_adjusted'][:, 0] - (self.vars['pres_adjusted'][:, 1]-self.vars['pres_adjusted'][:, 0])/2, (self.vars['pres_adjusted'][:, :-1]+self.vars['pres_adjusted'][:, 1:])/2, self.vars['pres_adjusted'][:, -1] + (self.vars['pres_adjusted'][:, -1] - self.vars['pres_adjusted'][:, -2])/2]
#
#         # presb = presb.filled(np.nan)
#         presb = presb.filled(99999.)
#         # timeb = [calendar.timegm(self.coords['time'][ip].timetuple())-0.5 for ip in range(self._nprof)]
#         # timeb.append(calendar.timegm(self.coords['time'][-1].timetuple())+0.5)
#         time = mpl.dates.date2num(self.coords['time'])
#         timeb = np.r_[time[0] - (time[1] - time[0])/2, (time[:-1]+time[1:])/2, time[-1] + (time[-1] - time[-2])/2]
#
#         years = mpl.dates.YearLocator()   # every year
#         months = mpl.dates.MonthLocator()  # every month
#         timeFmt = mpl.dates.DateFormatter('%Y-%m')
#
#         norm = mpl.colors.BoundaryNorm(lvs, 256)
#
#         var = self.vars[varnm].filled(np.nan)
#
#         fig, ax = plt.subplots()
#         for ip in range(self._nprof):
#             pmsh = ax.pcolormesh(timeb[ip:ip+2], presb[ip, :], var[ip, :][:, None], norm = norm, cmap = cmap)
#         # ax.set_xlim(timeb[0], timeb[-1])
#         ax.set_xlim(timeb[~maskb][0], timeb[~maskb][-1])
#         ax.set_ylim(0, depmax)
#         ax.invert_yaxis()
#
#         ax.xaxis.set_major_locator(months)
#         # ax.xaxis.set_minor_locator(months)
#         ax.xaxis.set_major_formatter(timeFmt)
#
#         fig.autofmt_xdate()
#         fig.colorbar(pmsh, ticks=lvs[::5]);
#         return fig, ax
