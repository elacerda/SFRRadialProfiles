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

def gimmedatrealmtfkingmorf(notrealmtfkingmorf = None):
    mtype = {
        'E0' : 0,
        'E1' : 0,
        'E2' : 0,
        'E3' : 0,
        'E4' : 0,
        'E5' : 0,
        'E6' : 0,
        'E7' : 0,
        'S0' : 1,
        'S0a' : 1,
        'Sa' : 2,
        'Sab' : 2,
        'Sb' : 3,
        'Sbc' : 4,
        'Sc' : 5,
        'Scd' : 5,
        'Sd' : 6,
        'Sdm' : -1,
        'Sm' : -1,
        'Ir' : -1,
    }
    realmtfkingmof = mtype[notrealmtfkingmorf]
    return realmtfkingmof   


def parser_args(args_str):
    paths = C.paths
    default_args = {
        'debug' : False,
        'underS06' : False,
        'noplot' : False,
        'whanSF' : None,
        'weiradprof' : True,
        'minpopx' : 0.05,
        'mintauv' : 0.05,
        'mintauvneb' : 0.05,
        'maxtauvneberr' : 999.,
        'rbinini' : 0.0,
        'rbinfin' : 2.0,
        'rbinstep' : 0.1,
        'gals_filename' : paths.califa_work_dir + 'listv20_q050.d15a.txt',
        'rgbcuts' : False,
        'filter_residual' : False,
        'gasprop' : False,
        'v_run' : -1,
        'nolinecuts' : False,
        'output' : None,
    }
    
    parser = ap.ArgumentParser(description = '%s' % args_str)
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
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
    # gals = np.array(['K0037', 'K0046', 'K0047', 'K0050', 'K0055', 'K0059', 'K0072',
    #    'K0080', 'K0085', 'K0087', 'K0093', 'K0096', 'K0099', 'K0103',
    #    'K0104', 'K0118', 'K0119', 'K0121', 'K0134', 'K0138', 'K0139',
    #    'K0162', 'K0170', 'K0173', 'K0189', 'K0281', 'K0339', 'K0340',
    #    'K0360', 'K0391', 'K0479', 'K0562', 'K0607', 'K0613', 'K0703',
    #    'K0705', 'K0744', 'K0778', 'K0787', 'K0822', 'K0826', 'K0844',
    #    'K0858', 'K0860', 'K0865', 'K0870', 'K0874', 'K0875', 'K0883',
    #    'K0908', 'K0912', 'K0916', 'K0917', 'K0919'], dtype='|S5')
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
    tSF__T = np.array([1, 3.2, 10, 100]) * 1e7
    N_T = len(tSF__T)

    # Z-time-scale array (index __U).
    tZ__U = np.array([1.0 , 2.0 , 5.0 , 8.0 , 11.3 , 14.2]) * 1.e9
    N_U = len(tZ__U)

    # SYN
    aSFRSD__Trg = np.ma.masked_all((N_T, NRbins, N_gals), dtype = np.float_)
    McorSD__Trg = np.ma.masked_all((N_T, NRbins, N_gals), dtype = np.float_)
                                    
    # EmLines
    aSFRSD_Ha__Trg = np.ma.empty((N_T, NRbins, N_gals), dtype = np.float_)
    EW_Ha__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    F_int_Ha__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    F_int_Hb__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    F_int_O3__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)
    F_int_N2__rg = np.ma.empty((NRbins, N_gals), dtype = np.float_)

    morf__g = np.ma.empty((N_gals), dtype = np.int_)
    N_zones__g = np.ma.empty((N_gals), dtype = np.int_)
    N_zones_notmasked__Tg = np.ma.zeros((N_T, N_gals), dtype = np.int_)

    # automatically read PyCASSO and EmLines data cubes.
    for iGal, K in C.loop_cubes(gals.tolist(), imax = maxGals, EL = True, GP = args.gasprop, v_run = args.v_run):        
    #for iGal in xrange(len(gals)):
        t_init_gal = time.clock()
        califaID = gals[iGal] 
        
        sit, verify = verify_files(K, califaID, EL = True, GP = args.gasprop)
        
        if verify is not True:
            print '<<< ', califaID, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue
        tipos, tipo, tipo_m, tipo_p = C.get_morfologia(califaID)
        
        morf__g[iGal] = gimmedatrealmtfkingmorf(tipos)
        print '>>> Doing' , iGal , califaID , 'hubtyp=', tipo,'(',tipos, ',',morf__g[iGal],')|  Nzones=' , K.N_zone

        N_zones__g[iGal] = K.N_zone

        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)

        # Saving for later :D
        zone_distance_HLR = K.zoneDistance_HLR
        zone_area_pc2 = K.zoneArea_pc2
        
        maskNotOk__Tz = np.zeros((N_T, K.N_zone), dtype = np.bool_)
        
        for iT, tSF in enumerate(tSF__T):
            tau_V__z = np.ma.masked_array(K.tau_V__z)
            x_Y__z, integrated_x_Y = calc_xY(K, tSF)
            if args.filter_residual is True:
                maskNotOk__Tz[iT] |= ~(K.filterResidual(w2 = 4600))
                C.debug_var(args.debug, pref = '    >>>',
                    tSF = '%.2fMyr' % (tSF / 1e6),
                    NnotOkResidual = (~(K.filterResidual())).sum()
                )
            if args.minpopx >= 0.:
                # Compute xOk "raw" image
                maskNotOk__Tz[iT] |= (x_Y__z < args.minpopx) 
                maskNotOk__Tz[iT] |= (tau_V__z < args.mintauv)
                C.debug_var(args.debug, pref = '    >>>',
                    tSF = '%.2fMyr' % (tSF / 1e6),
                    NnotOkxY = (x_Y__z < args.minpopx).sum(),
                    NnotOktauV = (tau_V__z < args.mintauv).sum(),
                    NnotOk = maskNotOk__Tz[iT].sum(),
                )
            
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
            
            ###### MASK EmLines ######
            maskOkLine = {}
            if args.nolinecuts is True:
                for l in lines_central_wl:
                    maskOkLine[l] = np.ones((K.N_zone), dtype = np.bool)
            else:
                if args.rgbcuts is True:
                    for l in lines_central_wl:
                        if args.gasprop is True:
                            pos = K.GP._dlcons[l]['pos']
                            sigma = K.GP._dlcons[l]['sigma']
                            snr = K.GP._dlcons[l]['SN']
                            if snr < 3.0: snr = 3.0
                        else:
                            pos, sigma, snr = 3.0, 3.0, 3.0
                        C.debug_var(args.debug, pref = l, pos = pos, sigma = sigma, snr = snr)
                        maskOkLine[l] = K.EL._setMaskLineFluxNeg(l) 
                        maskOkLine[l] |= K.EL._setMaskLineDisplacement(l, pos)
                        maskOkLine[l] |= K.EL._setMaskLineSigma(l, sigma)
                        maskOkLine[l] |= K.EL._setMaskLineSNR(l, snr)
                        maskOkLine[l] = ~(maskOkLine[l])
                        #C.debug_var(args.debug, maskOkLine = maskOkLine[l])
                else:
                    for l in lines_central_wl:
                        C.debug_var(args.debug, l = l)
                        minSNR = 3
                        maskOkLine[l] = K.EL._setMaskLineFluxNeg(l)
                        maskOkLine[l] |= K.EL._setMaskLineSNR(l, minSNR)
                        maskOkLine[l] = ~(maskOkLine[l])
                        #C.debug_var(args.debug, maskOkLine = maskOkLine[l])
                
            maskOkLines__z = np.bitwise_and(maskOkLine[Hb_central_wl], maskOkLine[O3_central_wl])
            maskOkLines__z = np.bitwise_and(maskOkLines__z, maskOkLine[Ha_central_wl])
            maskOkLines__z = np.bitwise_and(maskOkLines__z, maskOkLine[N2_central_wl])  
    
            maskOkBPT__z = np.ones((K.N_zone), dtype = np.bool)
            maskOkwhan__z = np.ones((K.N_zone), dtype = np.bool)
            if args.underS06:
                L = Lines()
                N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
                O3Hb = np.ma.log10(K.EL.O3_obs__z / K.EL.Hb_obs__z)
                maskOkBPT__z = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
            if args.whanSF is not None:
                N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
                WHa = K.EL.EW[i_Ha, :]
                maskOkwhan__z = np.greater(WHa, args.whanSF)
                #maskOkwhan__z = np.bitwise_and(np.greater(WHa, 3.), np.less_equal(N2Ha, -0.4))
            ##########################
            
            ########## tau_V_neb #########
            maskOkTauVNeb__z = np.ones((K.N_zone), dtype = np.bool)
            if args.mintauvneb >= 0:
                maskOkTauVNeb__z = (K.EL.tau_V_neb__z >= args.mintauvneb) & (K.EL.tau_V_neb_err__z <= args.maxtauvneberr)
            #maskOkTauVNeb__z &= maskOkLine[Ha_central_wl]
            #maskOkTauVNeb__z &= maskOkLine[Hb_central_wl]
            #maskOkTauVNeb__z &= maskOkBPT__z
            #maskOkTauVNeb__z &= maskOkwhan__z

            maskNotOk__Tz[iT] |= ~maskOkLines__z
            maskNotOk__Tz[iT] |= ~maskOkTauVNeb__z
            maskNotOk__Tz[iT] |= ~maskOkBPT__z
            maskNotOk__Tz[iT] |= ~maskOkwhan__z

            C.debug_var(args.debug, 
                NOkWhan = maskOkwhan__z.sum(),
                NOkHb = maskOkLine[Hb_central_wl].sum(),
                NOkO3 = maskOkLine[O3_central_wl].sum(),
                NOkHa = maskOkLine[Ha_central_wl].sum(),
                NOkN2 = maskOkLine[N2_central_wl].sum(),
                NOkBPT = maskOkBPT__z.sum(),
                NOkminTauVNeb = (K.EL.tau_V_neb__z >= args.mintauvneb).sum(),
                NOkmaxTauVNebErr = (K.EL.tau_V_neb_err__z <= args.maxtauvneberr).sum(),
                NOkTauVNeb = maskOkTauVNeb__z.sum(),
                NOkHaHb = (maskOkLine[Ha_central_wl] & maskOkLine[Hb_central_wl]).sum(),
                NOkO3N2 = (maskOkLine[O3_central_wl] & maskOkLine[N2_central_wl]).sum(),
                NOkLines = (maskOkLines__z).sum(),
            )
            
            N_zones_notmasked__Tg[iT, iGal] = K.N_zone - maskNotOk__Tz[iT].sum() 

        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        Mcor__z = K.Mcor__z
        McorSD__z = K.Mcor__z / K.zoneArea_pc2
        Mcor_GAL = K.Mcor_tot.sum()
        McorSD_GAL = K.McorSD__yx.mean()
        McorSD__r = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale=K.HLR_pix)
        
        # Composition by StarForming time scale
        for iT, tSF in enumerate(tSF__T):
            Mcor__z = np.ma.masked_array(K.Mcor__z)
            aux = calc_SFR(K, tSF)
            SFR__z = np.ma.masked_array(aux[0])
            SFRSD__z = np.ma.masked_array(aux[1])
            
            SFR__z[maskNotOk__Tz[iT]] = np.ma.masked
            SFRSD__z[maskNotOk__Tz[iT]] = np.ma.masked
            Mcor__z[maskNotOk__Tz[iT]] = np.ma.masked

            # Radial Profiles:
            #x_Y__r = K.zoneToRad(x_Y__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            McorSD__Trg[iT, :, iGal] = K.zoneToRad(Mcor__z, Rbin__r, rad_scale=K.HLR_pix)
            aSFRSD__Trg[iT, :, iGal] = K.zoneToRad(SFRSD__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
        
        ########## tau_V_neb #########
        tau_V_neb__z = np.ma.masked_array(K.EL.tau_V_neb__z, mask = ~maskOkTauVNeb__z)
        tau_V_neb_err__z = np.ma.masked_array(K.EL.tau_V_neb_err__z, mask = ~maskOkTauVNeb__z)
        tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
        tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z, extensive = False)
        tau_V_neb__r = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
        tau_V_neb_err__r = K.radialProfile(tau_V_neb_err__yx, Rbin__r, rad_scale = K.HLR_pix)

        ########### EW ###########
        EW_Ha__z = np.ma.masked_array(K.EL.EW[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
        EW_Ha__rg[:, iGal] = K.radialProfile(EW_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)

        #### intrinsic Ha Lum ####
        q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
        L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
        L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux) / L_sun        
        L_obs_Ha__z = np.ma.masked_array(L_obs__Lz[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        L_obs_Hb__z = np.ma.masked_array(L_obs__Lz[i_Hb, :], mask = ~(maskOkLine[Hb_central_wl]))
        L_obs_Ha_err__z = np.ma.masked_array(L_obs_err__Lz[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        L_obs_Hb_err__z = np.ma.masked_array(L_obs_err__Lz[i_Hb, :], mask = ~(maskOkLine[Hb_central_wl]))
        L_obs_HaHb__z = L_obs_Ha__z / L_obs_Hb__z
        # L_int_Ha__Lz intrinsic Ha luminosity 
        eHa = np.ma.exp(q[2] * tau_V_neb__z)
        # For the zones where I don't have values for tau_V_neb I don't correct the Lum_Ha
        L_int_Ha__z = np.where(maskOkTauVNeb__z, L_obs_Ha__z * eHa, L_obs_Ha__z)
        L_int_Ha__z = np.ma.masked_array(L_int_Ha__z, mask = ~maskOkTauVNeb__z)

        # L_int_Ha_err__Lz intrinsic Ha luminosity propagated error
        qq = q[2] / (q[0] - q[2])
        a = L_obs_Ha_err__z
        b = qq * L_obs_HaHb__z * L_obs_Hb_err__z
        L_int_Ha_err__z = np.where(maskOkTauVNeb__z == True, L_obs_Ha_err__z, eHa * np.sqrt(a ** 2.0 + b ** 2.0))
        L_int_Ha_err__z = np.ma.masked_array(L_int_Ha_err__z, mask = ~maskOkTauVNeb__z)
                
        ###### OTH BPT LINES #####
        F_obs_Ha__z = np.ma.masked_array(K.EL.flux[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        F_obs_Hb__z = np.ma.masked_array(K.EL.flux[i_Hb, :], mask = ~(maskOkLine[Hb_central_wl]))
        F_obs_O3__z = np.ma.masked_array(K.EL.flux[i_O3, :], mask = ~(maskOkLine[O3_central_wl]))
        F_obs_N2__z = np.ma.masked_array(K.EL.flux[i_N2, :], mask = ~(maskOkLine[N2_central_wl]))
        eHb = np.ma.exp(q[0] * tau_V_neb__z)
        eO3 = np.ma.exp(q[1] * tau_V_neb__z)
        eHa = np.ma.exp(q[2] * tau_V_neb__z)
        eN2 = np.ma.exp(q[3] * tau_V_neb__z)
        F_int_Ha__z = np.where(maskOkTauVNeb__z, F_obs_Ha__z * eHa, F_obs_Ha__z)
        F_int_Hb__z = np.where(maskOkTauVNeb__z, F_obs_Hb__z * eHb, F_obs_Hb__z)
        F_int_O3__z = np.where(maskOkTauVNeb__z, F_obs_O3__z * eO3, F_obs_O3__z)
        F_int_N2__z = np.where(maskOkTauVNeb__z, F_obs_N2__z * eN2, F_obs_N2__z)
        F_int_Ha__rg[:, iGal] = K.zoneToRad(F_int_Ha__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
        F_int_Hb__rg[:, iGal] = K.zoneToRad(F_int_Hb__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
        F_int_O3__rg[:, iGal] = K.zoneToRad(F_int_O3__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
        F_int_N2__rg[:, iGal] = K.zoneToRad(F_int_N2__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
        
        #### SFR and SigmaSFR ####
        # 3.13 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter
        SFR_Ha__z = np.ma.masked_array(3.13 * L_int_Ha__z.data / (1.e8), mask = L_int_Ha__z.mask)
        SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
        
        for iT in xrange(N_T):
            SFRSD_Ha__yx = K.zoneToYX(np.ma.masked_array(SFRSD_Ha__z, mask = maskNotOk__Tz[iT]), extensive = False)
            aSFRSD_Ha__Trg[iT, :, iGal] = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
        ####################################################
        ####################################################
        ####################################################
        #K.EL.close()
        #K.close()
        #del K
        print 'time per galaxy: %s %.2f' % (califaID, time.clock() - t_init_gal)

    print 'total time: %.2f' % (time.clock() - t_init_prog)

    colortipo = ['brown', 'red', 'orange', 'green', '#00D0C9', '#0076C9', 'blue']
    Ntype = len(colortipo)
    mtypes = [ 0, 1, 2, 3, 4, 5, 6 ]
    mtype_labels = [ 'E', 'S0', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd' ]
    halfbinstep = np.diff(mtypes)[0]/2.
    tickpos = np.linspace(mtypes[0] + halfbinstep, mtypes[-1] - halfbinstep, Ntype)
            
    aSFRSD__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.float_)
    aSFRSD_Ha__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.float_)
    nbinelem_aSFRSD__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.int_)
    nbinelem_aSFRSD_Ha__Ttr = np.empty((N_T, Ntype, NRbins), dtype = np.int_)
    perc_aSFRSD__pTtr = np.ma.masked_all((3, N_T, Ntype, NRbins), dtype = np.float_)
    perc_aSFRSD_Ha__pTtr = np.ma.masked_all((3, N_T, Ntype, NRbins), dtype = np.float_)

    Ngals__t = np.empty((Ntype), dtype = np.int_)
    for it in mtypes:
        mask_morf__g = (morf__g == it)
        Ngals__t[it] = np.sum(mask_morf__g)
        for iT in range(N_T):

            # G means masked g
            aSFRSD__rG = aSFRSD__Trg[iT][:, mask_morf__g]
            aSFRSD_Ha__rG = aSFRSD_Ha__Trg[iT][:, mask_morf__g]

            for ir in np.where(RbinCenter__r <= 2 + args.rbinstep)[0]:
                for ig in xrange(Ngals__t[it]):
                    if aSFRSD__rG[ir].mask[ig]:
                        aSFRSD__rG[ir, ig] = 0
                    if aSFRSD_Ha__rG[ir].mask[ig]:    
                        aSFRSD_Ha__rG[ir, ig] = 0.
                    
            aSFRSD__Ttr[iT][it] = aSFRSD__rG.mean(axis = 1)
            nbinelem_aSFRSD__Ttr[iT][it] = aSFRSD__rG.count(axis = 1)
            aSFRSD_Ha__Ttr[iT][it] = aSFRSD_Ha__rG.mean(axis = 1)
            nbinelem_aSFRSD_Ha__Ttr[iT][it] = aSFRSD_Ha__rG.count(axis = 1)
            
            for ir in xrange(NRbins):
                if nbinelem_aSFRSD__Ttr[iT, it, ir] > 2:
                    perc_aSFRSD__pTtr[:, iT, it, ir] = np.percentile(aSFRSD__rG[ir].compressed(), [16, 50, 84])

                if nbinelem_aSFRSD_Ha__Ttr[iT, it, ir] > 2:
                    perc_aSFRSD_Ha__pTtr[:, iT, it, ir] = np.percentile(aSFRSD__rG[ir].compressed(), [16, 50, 84])

    if not args.noplot:                    
        import matplotlib as mpl
        from matplotlib import pyplot as plt
        from matplotlib.pyplot import MultipleLocator
        from CALIFAUtils.plots import plot_text_ax
        
        iT = 1
        ylim = (-3, 0)
        f = plt.figure()
        NRows = 2
        NCols = 8
        page_size_inches = [12, 8]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        cmap = mpl.colors.ListedColormap(colortipo)
        suptitle = 'Ngals:%d  Nzones:%d  NRbins:%d  tSF:%.2f Myr' % ( np.sum(morf__g >= 0), N_zones_notmasked__Tg[iT].sum(), aSFRSD__Trg[iT].count(), tSF__T[iT]/1e6)
        f.suptitle(suptitle, fontsize = 15) 
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        for it in mtypes:
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
        ax.set_ylabel(r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(.2))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(.2))
        plt.setp(ax.get_xticklabels(), visible = False)
        ax.grid(which = 'major')
        ax.set_xlim(0,3)
        ax.set_ylim(ylim)
        
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        for it in mtypes:
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD_Ha__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
        ax.set_ylabel(r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$')
        ax.set_xlabel(r'R [HLR]')
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(.2))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(.2))
        ax.grid(which = 'major')
        ax.set_xlim(0,3)
        ax.set_ylim(ylim)
        
        for it in mtypes:
            ax = plt.subplot2grid(grid_shape, loc = (0, it + 1))
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD__pTtr[0, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD__pTtr[1, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD__pTtr[2, iT, it, :] * 1e6), 'k--')
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(.2))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(.2))
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
            plot_text_ax(ax, '%s (%d)' % (mtype_labels[it], Ngals__t[it]), 0.95, 0.95, 15, 'top', 'right', colortipo[it])
            ax.grid(which = 'major')
            ax.set_xlim(0,3)
            ax.set_ylim(ylim)
    
            ax = plt.subplot2grid(grid_shape, loc = (1, it + 1))
            ax.plot(RbinCenter__r, np.ma.log10(aSFRSD_Ha__Ttr[iT][it] * 1e6), c = colortipo[it], label = mtype_labels[it])
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD_Ha__pTtr[0, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD_Ha__pTtr[1, iT, it, :] * 1e6), 'k--')
            ax.plot(RbinCenter__r, np.ma.log10(perc_aSFRSD_Ha__pTtr[2, iT, it, :] * 1e6), 'k--')
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(.2))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(.2))
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid(which = 'major')
            ax.set_xlim(0,3)
            ax.set_ylim(ylim)
            
            C.debug_var(
                True, 
                tipoM = mtypes[it], 
                Nelemperbin = nbinelem_aSFRSD__Ttr[iT][it],
                Nelemperbin_Ha = nbinelem_aSFRSD_Ha__Ttr[iT][it],
            )
        
        f.subplots_adjust(bottom = 0.15, hspace = 0, wspace = 0, right = 0.95, left = 0.07)
        if args.output is not None:
            f.savefig(args.output)
        else:
            f.savefig('SFRSDRadialProfiles_%.2f.pdf' % args.minpopx)
        
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        # ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        # ax.set_ylabel(r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$')
        # ax.set_xlabel(r'R [HLR]')
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
     
        #cb = f.colorbar(sc, ticks = tickpos)
        #cb.ax.set_yticklabels(mtype_labels)
        #cb.ax.yaxis.set_ticks_position('none')    
         
    