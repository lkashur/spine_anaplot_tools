#if IS_NUMI
    #define SELCUT all_0mu2gamma_numi_cut
    #define FLASHCUT cuts::flash_cut_numi
    #define BEAMDIR 0.39431672, 0.04210058, 0.91800973

#else
    #define SELCUT all_0mu2gamma_bnb_cut
    #define FLASHCUT cuts::flash_cut_bnb
    #define BEAMDIR 0,0,1

#endif
