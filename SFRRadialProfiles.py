#!/usr/bin/python
#
# Lacerda@Corrego - 18/Nov/2015
#
import sys
import time
import numpy as np
import argparse as ap
import CALIFAUtils as C
from CALIFAUtils.lines import Lines
from CALIFAUtils.scripts import calc_xY
from CALIFAUtils.scripts import calc_SFR
from pystarlight.util import redenninglaws
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase
from CALIFAUtils.scripts import radialProfileWeighted

def my_morf(morf_in = None):
    mtype = {
        'Sa' : 0,
        'Sab' : 0,
        'Sb' : 1,
        'Sbc' : 2,
        'Sc' : 3,
        'Scd' : 3,
        'Sd' : 4,
        'Sdm' : 4,
        'Sm' : 4,
        'Ir' : 4,
        'E0' : -2,
        'E1' : -2,
        'E2' : -2,
        'E3' : -2,
        'E4' : -2,
        'E5' : -2,
        'E6' : -2,
        'E7' : -2,
        'S0' : -1,
        'S0a' : -1,
    }
    morf_out = mtype[morf_in]
    return morf_out  

def calc_O3N2(Hb_obs, O3_obs, Ha_obs, N2_obs, mask_zones = None, tau_V = None, correct = False):
    if mask_zones is not None:
        mask = mask_zones
    else:
        mask = np.zeros_like(Hb_obs, dtype = np.bool_)
    Hb = np.ma.masked_array(Hb_obs, mask = mask)
    O3 = np.ma.masked_array(O3_obs, mask = mask)
    Ha = np.ma.masked_array(Ha_obs, mask = mask)
    N2 = np.ma.masked_array(N2_obs, mask = mask)
    if correct is True:
        tau_V_m = np.ma.masked_array(tau_V, mask = mask)
        from pystarlight.util import redenninglaws
        q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
        Hb *= np.ma.exp(q[0] * tau_V_m) 
        O3 *= np.ma.exp(q[1] * tau_V_m) 
        Ha *= np.ma.exp(q[2] * tau_V_m) 
        N2 *= np.ma.exp(q[3] * tau_V_m)
    O3Hb = np.ma.log10(O3/Hb)
    N2Ha = np.ma.log10(N2/Ha)
    O3N2 = np.ma.log10(O3 * Ha / (N2 * Hb))
    return O3, Hb, N2, Ha, O3Hb, N2Ha, O3N2

def parser_args(args_str):
    paths = C.paths
    default_args = {
        'debug' : False,
        'underS06' : False,
        'noplot' : False,
        'whanSF' : None,
        'hdf5' : None,
        'weiradprof' : True,
        'rgbcuts' : False,
        'filter_residual' : False,
        'gasprop' : False,
        'nolinecuts' : False,
        'output' : None,
        'v_run' : -1,
        'minpopx' : np.finfo(np.float_).min,
        'mintauv' : np.finfo(np.float_).min,
        'mintauvneb' : np.finfo(np.float_).min,
        'maxtauvneberr' : np.finfo(np.float_).max,
        'rbinini' : 0.0,
        'rbinfin' : 2.0,
        'rbinstep' : 0.1,
        'gals_filename' : paths.califa_work_dir + 'listv20_q050.d15a.txt',
        'minEWHb' : np.finfo(np.float_).min,
    }
    
    parser = ap.ArgumentParser(description = '%s' % args_str)
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['hdf5'])
    parser.add_argument('--noplot', 
                        action = 'store_true',
                        default = default_args['noplot'])
    parser.add_argument('--gasprop', '-G',
                        action = 'store_true',
                        default = default_args['gasprop'])
    parser.add_argument('--nolinecuts' ,
                        action = 'store_true',
                        default = default_args['nolinecuts'])
    parser.add_argument('--filter_residual', '-R',
                        action = 'store_true',
                        default = default_args['filter_residual'])
    parser.add_argument('--weiradprof', '-W',
                        action = 'store_true',
                        default = default_args['weiradprof'])
    parser.add_argument('--v_run',
                        metavar = 'INT',
                        type = int,
                        default = default_args['v_run'])
    parser.add_argument('--gals_filename', '-L',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['gals_filename'])
    parser.add_argument('--output', '-o',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['output'])
    parser.add_argument('--underS06',
                        action = 'store_true',
                        default = default_args['underS06'])
    parser.add_argument('--whanSF',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['whanSF'])
    parser.add_argument('--rgbcuts',
                        action = 'store_true',
                        default = default_args['rgbcuts'])
    parser.add_argument('--minpopx',
                        help = 'Negative to disable mask in popx',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minpopx'])
    parser.add_argument('--mintauv',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['mintauv'])
    parser.add_argument('--mintauvneb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['mintauvneb'])
    parser.add_argument('--minEWHb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minEWHb'])
    parser.add_argument('--maxtauvneberr',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['maxtauvneberr'])
    parser.add_argument('--rbinini',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['rbinini'])
    parser.add_argument('--rbinfin',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['rbinfin'])
    parser.add_argument('--rbinstep',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['rbinstep'])
    

    return parser.parse_args()

def print_args(args):
    # dummy function to print a dictionary.
    for k, v in args.__dict__.iteritems():
        print k, v 

def verify_files(K, califaID, EL = True, GP = True):
    if K is None:
        print '<<< %s galaxy: miss files' % califaID
        return 0, False
    if EL == True and K.EL is None:
        print '<<< %s galaxy: miss EmLines files' % califaID
        return 1, False
        if K.EL.flux[0, :].sum() == 0.:
            print '<<< %s EmLines FITS problem' % califaID
            return 2, False
    if GP is True and K.GP._hdulist is None:
        print '<<< %s galaxy: miss gasprop file' % califaID
        return 2, False
    # Problem in FITS file
    return 0, True       

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse arguments 
    args = parser_args(sys.argv[0])
    print_args(args)    
    
    Zsun = 0.019

    # Creating radial bins.
    Rbin__r = np.arange(args.rbinini, args.rbinfin + args.rbinstep, args.rbinstep)
    RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
    NRbins = len(RbinCenter__r)
    RColor = [ 'r', 'y', 'b', 'k' ]
    RRange = [  .5, 1., 1.5, 2.  ]
    Rbin_oneHLR = [1. - args.rbinstep, 1. + args.rbinstep]
    
    # Reading galaxies file,
    #gals, _ = C.sort_gals(gals = args.gals_filename, order = 1)

    gals = np.array(['K0195', 'K0201', 'K0272', 'K0136', 'K0197', 'K0017', 'K0708',
                'K0090', 'K0631', 'K0128', 'K0780', 'K0067', 'K0160', 'K0923',
                'K0035', 'K0098', 'K0859', 'K0816', 'K0815', 'K0092', 'K0806',
                'K0602', 'K0018', 'K0210', 'K0782', 'K0051', 'K0279', 'K0903',
                'K0127', 'K0864', 'K0846', 'K0900', 'K0341', 'K0004', 'K0044',
                'K0829', 'K0171', 'K0589', 'K0911', 'K0588', 'K0814', 'K0781',
                'K0835', 'K0881', 'K0112', 'K0387', 'K0845', 'K0076', 'K0612',
                'K0888', 'K0068', 'K0851', 'K0832', 'K0893', 'K0318', 'K0101',
                'K0840', 'K0870', 'K0705', 'K0138', 'K0139', 'K0633', 'K0037',
                'K0391', 'K0046', 'K0121', 'K0744', 'K0162', 'K0080', 'K0339',
                'K0912', 'K0047', 'K0860', 'K0908', 'K0093', 'K0844', 'K0787',
                'K0103', 'K0703', 'K0055', 'K0085', 'K0059', 'K0072', 'K0118',
                'K0173', 'K0134', 'K0916', 'K0919', 'K0050', 'K0281', 'K0917',
                'K0822', 'K0865', 'K0087', 'K0875', 'K0170', 'K0826', 'K0063',
                'K0562', 'K0096', 'K0479', 'K0099', 'K0613', 'K0778', 'K0189',
                'K0874', 'K0607', 'K0119', 'K0883', 'K0858', 'K0592', 'K0673',
                'K0811', 'K0274', 'K0100', 'K0091', 'K0867', 'K0061', 'K0105',
                'K0024', 'K0502', 'K0653', 'K0077', 'K0818', 'K0131', 'K0872',
                'K0057', 'K0838', 'K0111', 'K0186', 'K0174', 'K0634', 'K0933',
                'K0036', 'K0049', 'K0936', 'K0163', 'K0032', 'K0220', 'K0314',
                'K0889', 'K0132', 'K0026', 'K0169', 'K0029', 'K0078', 'K0075',
                'K0066', 'K0319', 'K0135', 'K0902', 'K0083', 'K0020', 'K0386',
                'K0853', 'K0007', 'K0863', 'K0809', 'K0123', 'K0156', 'K0925',
                'K0833', 'K0850', 'K0894', 'K0168', 'K0194', 'K0038', 'K0651',
                'K0219', 'K0791', 'K0381', 'K0062', 'K0886', 'K0364', 'K0624',
                'K0113', 'K0177', 'K0740', 'K0663', 'K0311', 'K0932', 'K0019',
                'K0783', 'K0910', 'K0797', 'K0664', 'K0185', 'K0518', 'K0924',
                'K0414', 'K0013', 'K0804', 'K0326', 'K0615', 'K0672', 'K0837',
                'K0914', 'K0856', 'K0848', 'K0569', 'K0676', 'K0915', 'K0065',
                'K0842', 'K0146', 'K0176', 'K0115', 'K0043', 'K0001', 'K0021',
                'K0824', 'K0192', 'K0023', 'K0097', 'K0890', 'K0188', 'K0868',
                'K0086', 'K0180', 'K0153', 'K0054', 'K0133', 'K0191', 'K0107',
                'K0307', 'K0610', 'K0025', 'K0164', 'K0073', 'K0807', 'K0854',
                'K0834', 'K0821', 'K0871', 'K0151', 'K0070', 'K0010', 'K0830',
                'K0102', 'K0278', 'K0869', 'K0280', 'K0009', 'K0684', 'K0665',
                'K0126', 'K0873', 'K0927', 'K0774', 'K0931', 'K0789', 'K0297',
                'K0934', 'K0383', 'K0825', 'K0630', 'K0715', 'K0152', 'K0768',
                'K0907', 'K0823', 'K0041', 'K0388', 'K0005', 'K0769', 'K0208',
                'K0831', 'K0748', 'K0196', 'K0130', 'K0659', 'K0611', 'K0798',
                'K0880', 'K0184', 'K0879', 'K0820', 'K0754', 'K0309', 'K0861',
                'K0137', 'K0436', 'K0116', 'K0110', 'K0088', 'K0500', 'K0275',
                'K0124', 'K0181', 'K0476', 'K0779', 'K0580', 'K0437', 'K0183',
                'K0901', 'K0042', 'K0094', 'K0147', 'K0108', 'K0052', 'K0898',
                'K0203', 'K0813', 'K0190', 'K0876', 'K0277', 'K0810', 'K0515',
                'K0008', 'K0489', 'K0109', 'K0608', 'K0929', 'K0714', 'K0904',
                'K0896', 'K0149', 'K0122', 'K0857', 'K0165', 'K0002', 'K0178',
                'K0028', 'K0140', 'K0764', 'K0887', 'K0849', 'K0141', 'K0652',
                'K0895', 'K0885', 'K0878', 'K0711', 'K0361', 'K0843', 'K0758',
                'K0713', 'K0069', 'K0852', 'K0089', 'K0697', 'K0827', 'K0095',
                'K0158', 'K0143', 'K0805', 'K0117', 'K0040', 'K0486', 'K0921',
                'K0355', 'K0841', 'K0817', 'K0260', 'K0166', 'K0603', 'K0884',
                'K0148', 'K0039', 'K0930', 'K0707', 'K0125', 'K0935', 'K0016',
                'K0891', 'K0828', 'K0609', 'K0548', 'K0129', 'K0836', 'K0012',
                'K0581', 'K0232', 'K0866', 'K0003', 'K0909', 'K0205', 'K0144',
                'K0906', 'K0159', 'K0030', 'K0187', 'K0157', 'K0053', 'K0775',
                'K0273', 'K0045', 'K0084', 'K0614', 'K0031', 'K0071', 'K0204',
                'K0081', 'K0033', 'K0034', 'K0862', 'K0231', 'K0226', 'K0182',
                'K0027', 'K0312', 'K0060', 'K0606', 'K0528', 'K0150', 'K0306',
                'K0657', 'K0161', 'K0353', 'K0749', 'K0209', 'K0058', 'K0014',
                'K0179', 'K0475', 'K0937'], 
                dtype='|S5')

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # gals = np.array(['K0195', 'K0201', 'K0272', 'K0136', 'K0197', 'K0017', 'K0708'],
    #             dtype='|S5')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    N_gals = len(gals)
    maxGals = None
    if args.debug:
        maxGals = 10
        if N_gals > maxGals:
            N_gals = maxGals

    # SFR-time-scale array (index __T)
    base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
    #tSF__T = np.asarray(base.ageBase)
    #tSF__T = np.array([0.032 , 0.3 , 1.5, 14.2]) * 1.e9
    tSF__T = np.array([1, 3.162, 10, 100]) * 1e7
    N_T = len(tSF__T)

    # SYN
    _mask__Tg = [[] for _ in xrange(N_T)]
    _SFR__Tg = [[] for _ in xrange(N_T)]
    _SFRSD__Tg = [[] for _ in xrange(N_T)]
    aSFRSD__Trg = np.ma.masked_all((N_T, NRbins, N_gals), dtype = np.float_)
    McorSD__Trg = np.ma.masked_all((N_T, NRbins, N_gals), dtype = np.float_)
                                    
    # EmLines
    _mask_eml__Tg = [[] for _ in xrange(N_T)]
    _SFR_Ha__Tg = [[] for _ in xrange(N_T)]
    _SFRSD_Ha__Tg = [[] for _ in xrange(N_T)]
    aSFRSD_Ha__Trg = np.ma.masked_all((N_T, NRbins, N_gals), dtype = np.float_)
    EW_Ha__rg = np.ma.masked_all((NRbins, N_gals), dtype = np.float_)
    F_int_Ha__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    F_int_Hb__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    F_int_O3__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    F_int_N2__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    
    #Sebas
    SFR_int__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.float_)
    SFR_Ha_int__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.float_)
    SFR_Ha_int_masked__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.float_)
    O3N2M13_int__g = np.ma.masked_all((N_gals), dtype = np.float_)
    O3N2M13_int_masked__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.float_)

    morf__g = np.ma.masked_all((N_gals), dtype = np.int_)
    N_zones__g = np.ma.masked_all((N_gals), dtype = np.int_)
    N_zones_notmasked__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.int_)

    # automatically read PyCASSO and EmLines data cubes.
    EL = True
    baseCode = 'Bgsd6e'
    for iGal, K in C.loop_cubes(gals.tolist(), imax = maxGals, EL = EL, GP = args.gasprop, v_run = args.v_run, baseCode = baseCode):        
    #for iGal in xrange(len(gals)):
        t_init_gal = time.clock()
        califaID = gals[iGal] 
        
        sit, verify = verify_files(K, califaID, EL = EL, GP = args.gasprop)
        
        if verify is not True:
            print '<<< ', califaID, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue
        tipos, tipo, tipo_m, tipo_p = C.get_morfologia(califaID)
        
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # P = C.CALIFAPaths(baseCode = 'Bzca6e')
        # eml_file = P.get_emlines_file(califaID).replace('Bzca6e', 'Bgsd6e')
        # K.loadEmLinesDataCube(eml_file)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        
        morf__g[iGal] = my_morf(tipos)
        print '>>> Doing' , iGal , califaID , 'hubtyp=', tipo,'(',tipos, ',',morf__g[iGal],')|  Nzones=' , K.N_zone

        N_zones__g[iGal] = K.N_zone

        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)

        # Saving for later :D
        zone_distance_HLR = K.zoneDistance_HLR
        zone_area_pc2 = K.zoneArea_pc2
                    
        #######################
        ### RESID.EML MASKS ###
        #######################
        mask_eml__z = np.zeros((K.N_zone), dtype = np.bool_)
        Hb_central_wl = '4861'
        O3_central_wl = '5007'
        Ha_central_wl = '6563'
        N2_central_wl = '6583'
        lines_central_wl = [ 
            Hb_central_wl,
            O3_central_wl,
            Ha_central_wl,
            N2_central_wl,
        ]
        i_Hb = K.EL.lines.index(Hb_central_wl)
        i_O3 = K.EL.lines.index(O3_central_wl)
        i_Ha = K.EL.lines.index(Ha_central_wl)
        i_N2 = K.EL.lines.index(N2_central_wl)
        minSNR = 3.
        mask_lines_dict__L = {}
        if args.nolinecuts is True:
            for l in lines_central_wl:
                mask_lines_dict__L = np.zeros((K.N_zone), dtype = np.bool_)
        else:
            for l in lines_central_wl:
                C.debug_var(args.debug, l = l)
                mask_lines_dict__L[l] = K.EL._setMaskLineFluxNeg(l)
                mask_lines_dict__L[l] |= K.EL._setMaskLineSNR(l, minSNR)
        if args.rgbcuts is True:
            for l in lines_central_wl:
                if args.gasprop is True:
                    pos = K.GP._dlcons[l]['pos']
                    sigma = K.GP._dlcons[l]['sigma']
                    snr = K.GP._dlcons[l]['SN']
                    if snr < minSNR: snr = minSNR
                else:
                    pos, sigma, snr = 3.0, 3.0, 3.0
                mask_lines_dict__L[l] = K.EL._setMaskLineFluxNeg(l)
                mask_lines_dict__L[l] |= K.EL._setMaskLineDisplacement(l, pos)
                mask_lines_dict__L[l] |= K.EL._setMaskLineSigma(l, sigma)
                mask_lines_dict__L[l] |= K.EL._setMaskLineSNR(l, snr)
        mask_tau_V_neb__z = np.less(K.EL.tau_V_neb__z, args.mintauvneb)
        mask_tau_V_neb__z = np.ma.masked_array(mask_tau_V_neb__z, dtype = np.bool_, fill_value = True)
        mask_tau_V_neb__z = mask_tau_V_neb__z.data
        mask_tau_V_neb_err__z = np.greater(K.EL.tau_V_neb_err__z, args.maxtauvneberr)
        mask_tau_V_neb_err__z = np.ma.masked_array(mask_tau_V_neb_err__z, dtype = np.bool_, fill_value = True)
        mask_tau_V_neb_err__z = mask_tau_V_neb_err__z.data
        mask_EW_Hb__z = np.less(K.EL.EW[i_Hb], args.minEWHb)
        mask_EW_Hb__z = np.ma.masked_array(mask_EW_Hb__z, dtype = np.bool_, fill_value = True)
        mask_EW_Hb__z = mask_EW_Hb__z.data
        mask_bpt__z = np.zeros((K.N_zone), dtype = np.bool_)
        if args.underS06:
            L = Lines()
            N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
            O3Hb = np.ma.log10(K.EL.O3_obs__z / K.EL.Hb_obs__z)
            mask_BPT__z = ~(L.belowlinebpt('S06', N2Ha, O3Hb))
        mask_whan__z = np.zeros((K.N_zone), dtype = np.bool_)
        if args.whanSF:
            N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
            WHa = K.EL.EW[i_Ha, :]
            mask_whan__z = np.bitwise_or(np.less(WHa, 3.), np.greater(N2Ha, -0.4))
        mask_eml__z = np.zeros(K.N_zone, dtype = np.bool_)
        #C.debug_var(True, type_mask = type(mask_eml__z))
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_lines_dict__L[Hb_central_wl])
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_lines_dict__L[O3_central_wl])
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_lines_dict__L[Ha_central_wl])
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_lines_dict__L[N2_central_wl])
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_EW_Hb__z)
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_bpt__z)
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_whan__z)
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_tau_V_neb__z)
        mask_eml__z = np.bitwise_or(mask_eml__z, mask_tau_V_neb_err__z)
        #C.debug_var(True, type_mask = type(mask_eml__z), mask_eml = mask_eml__z, Nmask_eml = mask_eml__z.sum())
        #######################
        ### STARLIGHT MASKS ###
        #######################
        mask__Tz = np.zeros((N_T, K.N_zone), dtype = np.bool_)
        mask_syn__Tz = np.zeros((N_T, K.N_zone), dtype = np.bool_)
        for iT, tSF in enumerate(tSF__T):
            tau_V__z = np.ma.masked_array(K.tau_V__z)
            mask_tau_V = np.less(tau_V__z, args.mintauv) 
            x_Y__z, integrated_x_Y = calc_xY(K, tSF)
            mask_popx = np.less(x_Y__z, args.minpopx)
            mask_residual = np.zeros(K.N_zone, dtype = np.bool_)
            if args.filter_residual is True:
                mask_residual = ~(K.filterResidual(w2 = 4600))
            mask_syn__Tz[iT] = np.bitwise_or(np.bitwise_or(mask_tau_V, mask_popx), mask_residual)
            mask__Tz[iT] = np.bitwise_or(mask_syn__Tz[iT], mask_eml__z)
            N_zones_notmasked__Tg[iT, iGal] = K.N_zone - np.asarray(mask__Tz[iT], dtype = np.int).sum()
            #C.debug_var(True, iT = iT, type_mask = type(mask_syn__Tz[iT]), mask_syn = mask_syn__Tz[iT], Nmask_syn = mask_syn__Tz[iT].sum())
            C.debug_var(True, iGal = iGal, N_zones_notmasked = N_zones_notmasked__Tg[iT, iGal])
        #######################
        #######################
        #######################

        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        Mcor__z = K.Mcor__z
        McorSD__z = K.Mcor__z / K.zoneArea_pc2
        Mcor_GAL = K.Mcor_tot.sum()
        McorSD_GAL = K.McorSD__yx.mean()
        McorSD__r = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale=K.HLR_pix)
        
        # Composition by StarForming time scale
        for iT, tSF in enumerate(tSF__T):
            mask = mask__Tz[iT]
            _mask__Tg[iT].append(mask__Tz[iT])
            Mcor__z = np.ma.masked_array(K.Mcor__z)
            aux = calc_SFR(K, tSF)
            SFR__z = np.ma.masked_array(aux[0])
            SFRSD__z = np.ma.masked_array(aux[1])
            
            SFR__z[mask] = np.ma.masked
            SFRSD__z[mask] = np.ma.masked
            Mcor__z[mask] = np.ma.masked
            
            _SFR__Tg[iT].append(SFR__z)
            _SFRSD__Tg[iT].append(SFRSD__z)
            
            # Radial Profiles:
            #x_Y__r = K.zoneToRad(x_Y__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            McorSD__Trg[iT, :, iGal] = K.zoneToRad(Mcor__z, Rbin__r, rad_scale=K.HLR_pix)
            aSFRSD__Trg[iT, :, iGal] = K.zoneToRad(SFRSD__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            SFR_int__Tg[iT, iGal] = SFR__z.sum()
        
            ########## tau_V_neb #########
            mask_tauVNeb_aux__z = np.bitwise_or(mask_lines_dict__L[Ha_central_wl], mask_lines_dict__L[Hb_central_wl])
            mask_tauVNeb_aux__z = np.bitwise_or(mask_tauVNeb_aux__z, mask_tau_V_neb__z)
            mask_tauVNeb_aux__z = np.bitwise_or(mask_tauVNeb_aux__z, mask_tau_V_neb_err__z)
            tau_V_neb__z = np.ma.masked_array(K.EL.tau_V_neb__z, mask = mask_tauVNeb_aux__z)
            tau_V_neb_err__z = np.ma.masked_array(K.EL.tau_V_neb_err__z, mask = mask_tauVNeb_aux__z)
            tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
            tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z, extensive = False)
            tau_V_neb__r = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
            tau_V_neb_err__r = K.radialProfile(tau_V_neb_err__yx, Rbin__r, rad_scale = K.HLR_pix)
     
            ########### EW ###########
            EW_Ha__z = np.ma.masked_array(K.EL.EW[i_Ha, :], mask = mask_lines_dict__L[Ha_central_wl])
            EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
            EW_Ha__rg[:, iGal] = K.radialProfile(EW_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
     
            #### intrinsic Ha Lum ####
            q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
            L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
            L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux) / L_sun        
            L_obs_Ha__z = np.ma.masked_array(L_obs__Lz[i_Ha, :], mask = mask_lines_dict__L[Ha_central_wl])
            L_obs_Hb__z = np.ma.masked_array(L_obs__Lz[i_Hb, :], mask = mask_lines_dict__L[Hb_central_wl])
            L_obs_Ha_err__z = np.ma.masked_array(L_obs_err__Lz[i_Ha, :], mask = mask_lines_dict__L[Ha_central_wl])
            L_obs_Hb_err__z = np.ma.masked_array(L_obs_err__Lz[i_Hb, :], mask = mask_lines_dict__L[Hb_central_wl])
            L_obs_HaHb__z = L_obs_Ha__z / L_obs_Hb__z
            # L_int_Ha__Lz intrinsic Ha luminosity 
            eHa = np.ma.exp(q[2] * tau_V_neb__z)
            # For the zones where I don't have values for tau_V_neb I don't correct the Lum_Ha
            L_int_Ha__z = np.where(~mask_tauVNeb_aux__z, L_obs_Ha__z * eHa, L_obs_Ha__z)
            L_int_Ha__z = np.ma.masked_array(L_int_Ha__z, mask = mask_tauVNeb_aux__z)
     
            # L_int_Ha_err__Lz intrinsic Ha luminosity propagated error
            qq = q[2] / (q[0] - q[2])
            a = L_obs_Ha_err__z
            b = qq * L_obs_HaHb__z * L_obs_Hb_err__z
            L_int_Ha_err__z = np.where(~mask_tauVNeb_aux__z, L_obs_Ha_err__z, eHa * np.sqrt(a ** 2.0 + b ** 2.0))
            L_int_Ha_err__z = np.ma.masked_array(L_int_Ha_err__z, mask = mask_tauVNeb_aux__z)
                     
            ###### OTH BPT LINES #####
            F_obs_Ha__z = np.ma.masked_array(K.EL.flux[i_Ha, :], mask = mask_lines_dict__L[Ha_central_wl])
            F_obs_Hb__z = np.ma.masked_array(K.EL.flux[i_Hb, :], mask = mask_lines_dict__L[Hb_central_wl])
            F_obs_O3__z = np.ma.masked_array(K.EL.flux[i_O3, :], mask = mask_lines_dict__L[O3_central_wl])
            F_obs_N2__z = np.ma.masked_array(K.EL.flux[i_N2, :], mask = mask_lines_dict__L[N2_central_wl])
            eHb = np.ma.exp(q[0] * tau_V_neb__z)
            eO3 = np.ma.exp(q[1] * tau_V_neb__z)
            eHa = np.ma.exp(q[2] * tau_V_neb__z)
            eN2 = np.ma.exp(q[3] * tau_V_neb__z)
            F_int_Ha__z = np.where(~mask_tauVNeb_aux__z, F_obs_Ha__z * eHa, F_obs_Ha__z)
            F_int_Hb__z = np.where(~mask_tauVNeb_aux__z, F_obs_Hb__z * eHb, F_obs_Hb__z)
            F_int_O3__z = np.where(~mask_tauVNeb_aux__z, F_obs_O3__z * eO3, F_obs_O3__z)
            F_int_N2__z = np.where(~mask_tauVNeb_aux__z, F_obs_N2__z * eN2, F_obs_N2__z)
            F_int_Ha__rg[:, iGal] = K.zoneToRad(F_int_Ha__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            F_int_Hb__rg[:, iGal] = K.zoneToRad(F_int_Hb__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            F_int_O3__rg[:, iGal] = K.zoneToRad(F_int_O3__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            F_int_N2__rg[:, iGal] = K.zoneToRad(F_int_N2__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
             
            #### SFR and SigmaSFR ####
            # 3.13 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter
            SFR_Ha__z = np.ma.masked_array(3.13 * L_int_Ha__z.data / (1.e8), mask = L_int_Ha__z.mask)
            SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
            SFR_Ha_int__Tg[iT, iGal] = SFR_Ha__z.sum() 
            SFR_Ha_int_masked__Tg[iT, iGal] = np.ma.masked_array(3.13 * L_int_Ha__z.data / (1.e8), mask = mask).sum()

            _mask_eml__Tg[iT].append(SFR_Ha__z.mask)
            _SFR_Ha__Tg[iT].append(SFR_Ha__z)
            _SFRSD_Ha__Tg[iT].append(SFRSD_Ha__z)
             
            SFRSD_Ha__yx = K.zoneToYX(np.ma.masked_array(SFRSD_Ha__z, mask = mask__Tz[iT]), extensive = False)
            aSFRSD_Ha__Trg[iT, :, iGal] = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
            ####################################################
            ####################################################
            ####################################################
            O3, Hb, N2, Ha, _, _, _ = calc_O3N2(F_obs_Hb__z, F_obs_O3__z, F_obs_Ha__z, F_obs_N2__z, mask, tau_V = tau_V_neb__z, correct = True)
            O3N2M13_int_masked__Tg[iT, iGal] = 8.533 - 0.214 * np.ma.log10(O3.sum() * Ha.sum()/ (N2.sum() * Hb.sum()))
        #O3N2M13_int__g[iGal] = K.GP.EMPAB.integrated_O_O3N2_M13
        
 
        #K.EL.close()
        #K.close()
        #del K
        print 'time per galaxy: %s %.2f' % (califaID, time.clock() - t_init_gal)
    SFR__Tg = []
    SFRSD__Tg = []
    SFR_Ha__Tg = []
    SFRSD_Ha__Tg = []
    for iT in xrange(N_T):
        auxMask = np.hstack(_mask__Tg[iT])
        aux = np.hstack(_SFR__Tg[iT])
        SFR__Tg.append(np.ma.masked_array(aux, mask = auxMask)) 
        aux = np.hstack(_SFRSD__Tg[iT])
        SFRSD__Tg.append(np.ma.masked_array(aux, mask = auxMask)) 
        
        auxMask = np.hstack(_mask_eml__Tg[iT])
        aux = np.hstack(_SFR_Ha__Tg[iT])
        SFR_Ha__Tg.append(np.ma.masked_array(aux, mask = auxMask)) 
        aux = np.hstack(_SFRSD_Ha__Tg[iT])
        SFRSD_Ha__Tg.append(np.ma.masked_array(aux, mask = auxMask)) 
                
    print 'total time: %.2f' % (time.clock() - t_init_prog)

    if args.hdf5:
        t_init_hdf5 = time.clock() 
        import h5py
        filename = args.hdf5
        h5 = h5py.File(filename, 'w')
        D = {}
        D['HEADER/command_line'] = sys.argv
        D['HEADER/RbinIni'] = args.rbinini
        D['HEADER/RbinFin'] = args.rbinfin
        D['HEADER/RbinStep'] = args.rbinstep
        D['HEADER/xOkMin'] = args.minpopx
        D['HEADER/tauVOkMin'] = args.mintauv
        D['HEADER/tauVNebOkMin'] = args.mintauvneb
        D['HEADER/tauVNebErrMax'] = args.maxtauvneberr
        D['HEADER/RbinCenter__r'] = RbinCenter__r
        D['HEADER/Rbin__r'] = Rbin__r
        D['HEADER/tSF__T'] = tSF__T
        D['HEADER/N_rbins'] = NRbins
        D['HEADER/N_gals'] = N_gals
        D['HEADER/N_T'] = N_T
        D['HEADER/gals'] = gals,
        D['data/N_zones_notmasked__Tg'] = N_zones_notmasked__Tg.data
        D['mask/N_zones_notmasked__Tg'] = N_zones_notmasked__Tg.mask
        D['data/morf__g'] = morf__g.data
        D['mask/morf__g'] = morf__g.mask
        D['data/aSFRSD__Trg'] = aSFRSD__Trg.data
        D['mask/aSFRSD__Trg'] = aSFRSD__Trg.mask
        D['data/McorSD__Trg'] = McorSD__Trg.data
        D['mask/McorSD__Trg'] = McorSD__Trg.mask
        D['data/aSFRSD_Ha__Trg'] = aSFRSD_Ha__Trg.data
        D['mask/aSFRSD_Ha__Trg'] = aSFRSD_Ha__Trg.mask
        D['data/SFR_Ha_int__Tg'] = SFR_Ha_int__Tg.data
        D['mask/SFR_Ha_int__Tg'] = SFR_Ha_int__Tg.mask
        D['data/SFR_Ha_int_masked__Tg'] = SFR_Ha_int_masked__Tg.data
        D['mask/SFR_Ha_int_masked__Tg'] = SFR_Ha_int_masked__Tg.mask
        D['data/SFR_int__Tg'] = SFR_int__Tg.data
        D['mask/SFR_int__Tg'] = SFR_int__Tg.mask
        D['data/O3N2M13_int__g'] = O3N2M13_int__g
        D['data/O3N2M13_int_masked__Tg'] = O3N2M13_int_masked__Tg.data
        D['mask/O3N2M13_int_masked__Tg'] = O3N2M13_int_masked__Tg.mask
        for iT in xrange(N_T):
            D['data/SFR__Tg/%d'%iT] = SFR__Tg[iT].data
            D['mask/SFR__Tg/%d'%iT] = SFR__Tg[iT].mask
            D['data/SFRSD__Tg/%d'%iT] = SFRSD__Tg[iT].data
            D['mask/SFRSD__Tg/%d'%iT] = SFRSD__Tg[iT].mask
            D['data/SFR_Ha__Tg/%d'%iT] = SFR_Ha__Tg[iT].data
            D['mask/SFR_Ha__Tg/%d'%iT] = SFR_Ha__Tg[iT].mask
            D['data/SFRSD_Ha__Tg/%d'%iT] = SFRSD_Ha__Tg[iT].data
            D['mask/SFRSD_Ha__Tg/%d'%iT] = SFRSD_Ha__Tg[iT].mask
        
        for k in D.keys():
            try:
                h5.create_dataset(k, data = D[k], compression = 'gzip', compression_opts = 4)
            except TypeError:
                h5.create_dataset(k, data = D[k])
        h5.close()
        print 'time hdf5: %.2f' % (time.clock() - t_init_hdf5)

    