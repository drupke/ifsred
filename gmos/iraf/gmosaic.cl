# Copyright(c) 2002-2013 Association of Universities for Research in Astronomy, Inc.

procedure gmosaic (inimages)

# Mosaic the 3 GMOS ccds
# geotran gdata$g24173[1][33:2080,1:4608] chip1 "" "" xshift=-2.50 yshift=-1.58 \
# xrotation=-0.0027827 yrotation=-0.0027827
#
# imcopy gdata$g24173[2][33:2080,1:4608] chip2
# geotran gdata$g24173[3][1:2048,1:4608] chip3 "" "" xshift=3.3723 yshift=-1.89 \
# xrotation=-0.056477 yrotation=-0.056477
# followed by imtile

# Version Feb 28, 2002  IJ,BM  v1.3 release
#         Apr 26, 2002  IJ     handle raw images with no complaints, WCS warning
#                              optional cleaning of OBSMODE=IMAGE
#         Jun 21, 2002  IJ     fixed warning for missing WCS
#         Aug 26, 2002  KGB    added parameter to specify interpolation type, default as before
#         Aug 26, 2002  IJ     changed parameter name to geointer, parameter encoding
#         Aug 27, 2002  IJ     allowed values on geointer, geointer in output header
#         Sept 20, 2002 IJ     v1.4 release
#         Oct 14, 2002  BM     generalize to GMOS-S
#         Jan 16, 2003  BM     alternate amps 5,6 for 6 amp mode
#         Jan 24, 2003  ML     updated numbers for GMOS-S 
#         Feb 6, 2003   IJ     fix the binx2 bug - I hope; enabled binx2 use for GMOS-S
#         Feb 14, 2003  IJ     fixed handling of WCS.
#         Feb 19, 2003  IJ     tweaked GMOS-N ccd3 shift with 0.5pix
#         March 4, 2003 IJ     fixed rotation for Xbin!=Ybin case
#                              fix rotation correction for multiple images processed
#                              merged in new GMOS-S values from BM/ML, corrected lxbshift[2,3]
#         Mar 19, 2003  BM     default pix value in chip gaps is -1, generalize mbpm for GMOS-N,S
#         Mar 20, 2003  BM     changed lxbshift[2,3] for GMOS-S to get chip-gap separation correct (checked with images)
#         Apr 10, 2003  BM     generalized to handle original GMOS-S geometry and geometry with new CCD1
#         May 9, 2003   BM     swap order for amps 1,2 for GMOS-S in 6 amp mode with new CCD
#         May 9, 2003   IJ     clean-up logfile output for rotations
#         May 15, 2003  BM     new GMOS-S numbers, 6amp mode for GMOS-N and GMOS-S
#         Jun 2, 2003   IJ     transparent handling of old BPM names for GMOS-N. Released as bugfix.
#         Jun 18, 2003  BM     more mods for new best amp ordering for GMOS-S 6amp mode
#         Jun 19, 2003  BM     bug fix
#         Aug 26, 2003  KL     IRAF2.12 - new parameters
#                                 hedit: addonly;
#                                 imstat: nclip, lsigma, usigma, cache
#         Dec 05, 2003  BM     New yoffsets for GMOS-S CCD1,3 from IFU anal.
#         Apr 07, 2004  BM     New rotations and xoffsets for GMOS-S CCDs
#         Jun 25, 2004  BM     Add status=0 at start, otherwise can get stuck
#                              at status=1
#
#    ************************************************************************
#    **** Gemini developers:
#    ****   See CVS and package revision notes for more recent modifications.
#    ************************************************************************

char    inimages    {prompt="GMOS images to mosaic"}                 # OLDP-1-input-primary-single-prefix=m
char    outimages   {"" ,prompt="Output images or list"}            # OLDP-1-output
char    outpref     {"m", prompt="Prefix for output images"}        # OLDP-4
bool    fl_paste    {no, prompt="Paste images instead of mosaic"}   # OLDP-3
bool    fl_vardq    {no, prompt="Propagate the variance and data quality planes"}       # OLDP-2 
bool    fl_fixpix   {no, prompt="Interpolate across chip gaps"}     # OLDP-3
bool    fl_clean    {yes, prompt="Clean imaging data outside imaging field"}            # OLDP-3
char    geointer    {"linear", min="linear|nearest|poly3|poly5|spline3|sinc", prompt="Interpolant to use with geotran"} # OLDP-2
char    gap         {"default", prompt="Gap between the CCDs in unbinned pixels"}       # OLDP-3
char    bpmfile     {"gmos$data/chipgaps.dat", prompt="Info on location of chip gaps"}  # OLDP-3
char    statsec     {"default", prompt="Statistics section for cleaning"}               # OLDP-3
char    obsmode     {"IMAGE", prompt="Value of key_obsmode for imaging data"}           # OLDP-3
char    sci_ext     {"SCI", prompt="Science extension(s) to mosaic, use '' for raw data"}   # OLDP-3
char    var_ext     {"VAR", prompt="Variance extension(s) to mosaic"}                   # OLDP-3
char    dq_ext      {"DQ", prompt="Data quality extension(s) to mosaic"}                # OLDP-3
char    mdf_ext     {"MDF", prompt="Mask definition file extension name"}               # OLDP-3
char    key_detsec  {"DETSEC", prompt="Header keyword for detector section"}            # OLDP-3
char    key_ccdsec  {"CCDSEC", prompt="Header keyword for CCD section"}            # OLDP-3
char    key_datsec  {"DATASEC", prompt="Header keyword for data section"}               # OLDP-3
char    key_ccdsum  {"CCDSUM", prompt="Header keyword for CCD binning"}                 # OLDP-3
char    key_obsmode {"OBSMODE", prompt="Header keyword for observing mode"}             # OLDP-3
char    logfile     {"", prompt="Logfile"}                          # OLDP-1
bool    fl_real     {no, prompt="Convert file to real before transforming"}
bool    verbose     {yes, prompt="Verbose"}                         # OLDP-4
int     status      {0, prompt="Exit status (0=good)"}              # OLDP-4
struct  *scanfile   {"", prompt="Internal use only"}                # OLDP-4

begin

    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outpref = ""
    char    l_logfile = ""
    char    l_bpmfile = ""
    char    l_geointer = ""
    char    l_statsec = ""
    char    l_obsmode = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_mdf_ext = ""
    char    l_key_detsec = ""
    char    l_key_ccdsec = ""
    char    l_key_datsec = ""
    char    l_key_ccdsum = ""
    char    l_key_obsmode = ""
    char    l_gap
    bool    l_fl_paste, l_fl_vardq, l_fl_fixpix, l_fl_clean
    bool    l_fl_real, l_verbose
    
    char    filename, inimg[200], outimg[200]   #maxfiles=200
    char    tmpfile, tmpinimg, tmpoutimg, tmpimg, tmptile, tmplog, tmpmdf
    char    tmpbpmtext, tmpbpm, tmpgtiledimg
    char    tmpchipsci[3], tmpchipvar[3], tmpchipdq[3]
    char    tmpgeoinsci[3], tmpgeoinvar[3], tmpgeoindq[3]
    char    imglist, outputstr
    char    sciext, varext, dqext
    char    pastesci, pastevar, pastedq, pixtype
    char    dampswap, paramstr, logmsg, errmsg
    char    keyfound, currentimg, real_currentimg, groups
    char    dateobs, detsec, datasec
    char    thisstatsec
    int     maxfiles, Xbingap           #constants defined in 'Default values'
    int     defstatsec[4]               #constants defined in 'Default values'
    int     nimages, noutimages, ngood, ndetroi
    int     Xbin, Ybin, thisimggap, nsci, inst, iccd
    int     ampsperchip, nctile, nltile
    int     pos, jstart, jstop
    int     i, j, jj, k                 # i=input img ; j=extension ; k=chip
    int     junk, tmpint, gapvalue, gap_to_add
    real    geoxshift[3], geoyshift[3], medsky, crpix1
    real    orig_sci_repvalue, orig_var_repvalue, orig_dq_repvalue
    real    sci_repvalue, var_repvalue, dq_repvalue
    real    croi_scivalue, croi_varvalue, croi_dqvalue
    bool    cleanthis, fixthis, isbinned
    bool    success, run_gtile
    bool    debug = no
    bool    ishamamatsu = no
    struct  l_struct
    int     chipstart, chipstop, cnt, number_of_ampsperccd[3]
    int     number_of_ccds, namps
    int     tx1,tx2,ty1,ty2, ampsorder[12], crpix1_adjust
    int     detx1[12], detx2[12], dety1[12], dety2[12]
    int     datx1[12], datx2[12], daty1[12], daty2[12]
    int     ccdx1[12], ccdx2[12], ccdy1[12], ccdy2[12]
    char    ccdsec, correctextn, dettype, detector

    # var[instrument, set, chip]
    #   instrument: North = 1; South = 2
    # chip set:   
    #    original = 1; blue set = 2 (GMOS-S only and new GMOS-N EEV2 CCDs; 
    #        hamamatsu upgrades = 3
    # chip: 3-amp mode: chip 1 = 1; chip 2 = 2; chip 3 = 3
    #       6-amp mode: amp1-2 => chip 1; amp3-4 => chip 2; amp5-6 => chip 3
    #       12-amp mode: amp1-4 => chip 1; amp5-8 => chip 2; amp9-12 => chip 3
    # 
    real    lxshift[2,3,3], lyshift[2,3,3], lrotation[2,3,3], lxbshift[2,3,3]
    real    lxrot[3], lyrot[3]
    int     refextn, refextn_naxis1

    ############################################
    ### Setting constants for chips geometry ###
    # Kept for prosperity - MS
    #    int     defrefextn[3] # 1-> 3 amps ; 2-> 6 amps; 3-> 12 amps
    #defrefextn[1] = 2       # extn 2 is the fixed one in 3-amps mode
    #defrefextn[2] = 3       # extn 3 is the fixed one in 6-amps mode (x4 too)
    #defrefextn[3] = 6       # extn 6 is the fixed one in 12-amps mode #hcode

    ####################################################
    # Transformation: GMOS-NORTH - Old CCDs
    #
    #   Feb 14, 2003
    #       Chip 1:
    #           lxshift[1,1,1]      =   -2.50
    #           lyshift[1,1,1]      =   -1.58
    #           lrotation[1,1,1]    =   -0.004 
    #       Chip 2:
    #           lxshift[1,1,2]      =   0.0
    #           lyshift[1,1,2]      =   0.0
    #           lrotation[1,1,2]    =   0.0
    #       Chip 3: (Feb 14, 2003)
    #           lxshift[1,1,3]      =   3.8723
    #           lyshift[1,1,3]      =   -1.86
    #           lrotation[1,1,3]    =   -0.046
    #
    #
    # Corrections for BINNED pixels: (still unbinned pixel values in shift)
    #
    #       Chip 1:
    #           lxbshift[1,1,1] =   -3.5
    #       Chip 2:
    #           lxbshift[1,1,2] =   0.0
    #       Chip 3:
    #           lxbshift[1,1,3] =   4.8723
    #
    #
    # Previous iterations:
    #       Chip 1:
    #           lxshift[1]=-2.50 ; lyshift[1]=-1.58 ; lrotation[1]=-0.0027827
    #           lxshift[1]=-2.50 ; lyshift[1]=-1.62 ; lrotation[1]=-0.004
    #               (IFU-optimized)
    #       Chip 3:
    #           lxshift[3]=3.3723 ; lyshift[3]=-1.89 ; lrotation[3]=-0.056477
    #           lxshift[3]=3.3723 ; lyshift[3]=-1.99 ; lrotation[3]=-0.054
    #               (IFU-optimized)
    #           lxshift[3]=3.3723 ; lyshift[3]=-2.05 ; lrotation[3]=-0.056 
    #           lxshift[3]=3.3723 ; lyshift[3]=-1.89 ; lrotation[3]=-0.048 
    #           lxshift[3]=3.3723 ; lyshift[3]=-1.86 ; lrotation[3]=-0.046
    #
    #       Binned:
    #           lxbshift[1]=-3.5 ; lxbshift[2]=0 ; lxbshift[3]=4.3723
    
    # GMOS-N - Unbinned - Chip 1
    lxshift[1,1,1]      =   -2.50
    lyshift[1,1,1]      =   -1.58
    lrotation[1,1,1]    =   -0.004
    
    # GMOS-N - Unbinned - Chip 2
    lxshift[1,1,2]      =   0.0
    lyshift[1,1,2]      =   0.0
    lrotation[1,1,2]    =   0.0
    
    # GMOS-N - Unbinned - Chip 3
    lxshift[1,1,3]      =   3.8723
    lyshift[1,1,3]      =   -1.86
    lrotation[1,1,3]    =   -0.046
    
    # GMOS-N - Binned - Chips 1,2,3
    lxbshift[1,1,1] = -3.5
    lxbshift[1,1,2] = 0.0
    lxbshift[1,1,3] = 4.8723
 
    ####################################################
    ## Transformation: GMOS-NORTH - New e2vDD CCDs
    ##M Needs checking/updating

    # GMOS-N - Unbinned - Chip 1 
    lxshift[1,2,1]      =   -2.7000 
    lyshift[1,2,1]      =   -0.7490
    lrotation[1,2,1]    =   -0.0090
    
    # GMOS-N - Unbinned - Chip 2 
    lxshift[1,2,2]      =    0.0 
    lyshift[1,2,2]      =    0.0 
    lrotation[1,2,2]    =    0.0 
    
    # GMOS-N - Unbinned - Chip 3 
    lxshift[1,2,3]      =    2.8014
    lyshift[1,2,3]      =    2.0500
    lrotation[1,2,3]    =   -0.0030
    
    # GMOS-N - Binned - Chips 1,2,3 
    lxbshift[1,2,1]     =   -3.7000
    lxbshift[1,2,2]     =    0.0 
    lxbshift[1,2,3]     =    3.8014 
   
    ####################################################
    ## Transformation: GMOS-NORTH - New CCDs (Hamamatsu)
    ##M Needs checking/updating

    # GMOS-N - Unbinned - Chip 1 New CCDs #hcode 
    lxshift[1,3,1]      =   -2.50 #hcode 
    lyshift[1,3,1]      =   -1.58 #hcode 
    lrotation[1,3,1]    =   -0.004 #hcode 
    
    # GMOS-N - Unbinned - Chip 2 New CCDs #hcode 
    lxshift[1,3,2]      =   0.0 #hcode 
    lyshift[1,3,2]      =   0.0 #hcode 
    lrotation[1,3,2]    =   0.0 #hcode 
    
    # GMOS-N - Unbinned - Chip 3 New CCDs #hcode 
    lxshift[1,3,3]      =   3.8723 #hcode 
    lyshift[1,3,3]      =   -1.86 #hcode 
    lrotation[1,3,3]    =   -0.046 #hcode 
    
    # GMOS-N - Binned - Chips 1,2,3 New CCDs #hcode 
    lxbshift[1,3,1] = -3.5 #hcode 
    lxbshift[1,3,2] = 0.0 #hcode 
    lxbshift[1,3,3] = 4.8723 #hcode 
       
    ####################################################
    # Transformation: GMOS-SOUTH - OLD CCDs
    #
    #       Chip 1:
    #           lxshift[2,1,1]      =   -1.44
    #           lyshift[2,1,1]      =   5.46
    #           lrotation[2,1,1]    =   -0.01
    #       Chip 2:
    #           lxshift[2,1,2]      =   0.0
    #           lyshift[2,1,2]      =   0.0
    #           lrotation[2,1,2]    =   0.0
    #       Chip 3:
    #           lxshift[2,1,3]      =   7.53
    #           lyshift[2,1,3]      =   9.57
    #           lrotation[2,1,3]    =   0.02
    #
    # Corrections for BINNED pixels: (still unbinned pixel values in shift)
    #
    #       Chip 1:
    #           lxbshift[1,1,1] =   -2.44
    #       Chip 2:
    #           lxbshift[1,1,2] =   0.0
    #       Chip 3:
    #           lxbshift[1,1,3] =   7.53
    #
    # Previous iterations:
    #       Chip 1:
    #           lxshift[1]=5.07 ; lyshift[1]=5.49 ; lrotation[1]=0.0
    #       Chip 3:
    #           lxshift[3]=-0.81 ; lyshift[3]=9.37 ; lrotation[3]=0.02
    #
    #       Binned:
    #           lxbshift[2,1]=4.07 ; lxbshift[2,2]=0 ; lxbshift[2,3]=-1.81

    # GMOS-S - Unbinned - Chip 1 - Old CCDs
    lxshift[2,1,1]      =   -1.44
    lyshift[2,1,1]      =   5.46
    lrotation[2,1,1]    =   -0.01

    # GMOS-S - Unbinned - Chip 2 - Old CCDs
    lxshift[2,1,2]      =   0.0
    lyshift[2,1,2]      =   0.0
    lrotation[2,1,2]    =   0.0

    # GMOS-S - Unbinned - Chip 3 - Old CCDs
    lxshift[2,1,3]      =   7.53
    lyshift[2,1,3]      =   9.57
    lrotation[2,1,3]    =   0.02
    
    # GMOS-S - Binned - Chips 1,2,3 - Old CCDs
    lxbshift[2,1,1] = -2.44
    lxbshift[2,1,2] = 0.0
    lxbshift[2,1,3] = 7.53

    ####################################################
    # Transformation: GMOS-SOUTH - NEW CCDs
    #
    #   April 7, 2004:
    #       Chip 1:
    #           lxshift[2,2,1]      =   -1.49
    #           lyshift[2,2,1]      =   -0.22
    #           lrotation[2,2,1]    =   0.011
    #       Chip 2:
    #           lxshift[2,2,2]      =   0.0
    #           lyshift[2,2,2]      =   0.0
    #           lrotation[2,2,2]    =   0.0
    #       Chip 3:
    #           lxshift[2,2,3]      =   4.31
    #           lyshift[2,2,3]          2.04
    #           lrotation[2,2,3]    =   0.012
    #
    # Corrections for BINNED pixels: (still unbinned pixel values in shift)
    #
    #       Chip 1:
    #           lxbshift[2,2,1]     =   -2.49
    #       Chip 2:
    #           lxbshift[2,2,2]     =   0.0
    #       Chip 3:
    #           lxbshift[2,2,3]     =   5.31
    #
    # Previous iterations:
    #       Chip 1:
    #           lxshift[1]=-2.43 ; lyshift[1]=0.0 ; lrotation[1]=-0.01
    #           lxshift[1]=-0.49 ; lyshift[1]=-0.10 ; lrotation[1]=0.012
    #               (May 15, 2003)
    #       Chip 3:
    #           lxshift[3]=2.43 ; lyshift[3]=2.12 ; lrotation[3]=0.0
    #           lxshift[3]=2.31 ; lyshift[3]=1.91 ; lrotation[3]=0.015
    #               (May 15, 2003)
    #
    #       Binned:
    #           lxbshift[1]=-3.43 ; lxbshift[2]=0 ; lxbshift[3]=2.43

    # GMOS-S - Unbinned - Chip 1 - New CCDs
    lxshift[2,2,1]      =   -1.49
    lyshift[2,2,1]      =   -0.22
    lrotation[2,2,1]    =   0.011

    # GMOS-S - Unbinned - Chip 2 - New CCDs
    lxshift[2,2,2]      =   0.0
    lyshift[2,2,2]      =   0.0
    lrotation[2,2,2]    =   0.0

    # GMOS-S - Unbinned - Chip 3 - New CCDs
    lxshift[2,2,3]      =   4.31
    lyshift[2,2,3]      =   2.04
    lrotation[2,2,3]    =   0.012

    # GMOS-S - Binned - Chips 1,2,3 - New CCDs
    lxbshift[2,2,1] = -2.49
    lxbshift[2,2,2] = 0.0
    lxbshift[2,2,3] = 5.31

    ########################################################
    # Dates of chip swapping
    #
    #   dampswap:
    #       Date when 'best amp' for new chip 1 was swapped to L
    #       => June 18, 2003
    #
    # !!! This is not actually used by the code anymore as the 
    # !!! chips are sorted according to their data section.
    # !!! The information was kept here for reference only.
    
    dampswap = "2003-06-18"

    ###########################################################
    ###########################################################
    
    # Set up caches
    unlearn geotran         # ??? is this necessary ???
    cache ("gloginit", "gemextn", "tinfo", "gemdate")
    
    # Default values
    maxfiles = 200
    ngood = 0
    status = 0

    # Replacement values used for gaps - these get reset for custom ROIs for
    # display purposes
    orig_sci_repvalue = -1
    orig_var_repvalue = -1
    orig_dq_repvalue = -1

    # Custom ROI replacement values - used to replace the unilluminated pixels
    # of each of the arrays
    croi_scivalue = 0
    croi_varvalue = 0
    croi_dqvalue = 16
    
    # Read in parameter values  (fscan removes blank spaces)
    junk = fscan (inimages, l_inimages)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outpref, l_outpref)
    junk = fscan (gap, l_gap)
    l_fl_paste = fl_paste
    l_fl_vardq = fl_vardq
    l_fl_fixpix = fl_fixpix
    l_fl_clean = fl_clean
    junk = fscan (geointer, l_geointer)
    junk = fscan (bpmfile, l_bpmfile)
    junk = fscan (statsec, l_statsec)
    junk = fscan (obsmode, l_obsmode)
    
    junk = fscan (sci_ext, l_sci_ext)
    junk = fscan (var_ext, l_var_ext)
    junk = fscan (dq_ext, l_dq_ext)
    junk = fscan (mdf_ext, l_mdf_ext)
    
    junk = fscan (key_detsec, l_key_detsec)
    junk = fscan (key_ccdsec, l_key_ccdsec)
    junk = fscan (key_datsec, l_key_datsec)
    junk = fscan (key_ccdsum, l_key_ccdsum)
    junk = fscan (key_obsmode, l_key_obsmode)
    
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    l_fl_real = fl_real
    
    # Temporary files
    tmpfile   = mktemp("tmpfile")
    tmpinimg  = mktemp("tmpinimg")
    tmpoutimg = mktemp("tmpoutimg")
    tmplog    = mktemp("tmplog")
    tmpgtiledimg = mktemp ("tmpgtiledimg")//".fits"
    
    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string. Line delimeter is '\n'.
    
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "fl_paste       = "//fl_paste.p_value//"\n"
    paramstr += "fl_vardq       = "//fl_vardq.p_value//"\n"
    paramstr += "fl_fixpix      = "//fl_fixpix.p_value//"\n"
    paramstr += "fl_clean       = "//fl_clean.p_value//"\n"
    paramstr += "geointer       = "//geointer.p_value//"\n"
    paramstr += "gap            = "//gap.p_value//"\n"
    paramstr += "bpmfile        = "//bpmfile.p_value//"\n"
    paramstr += "statsec        = "//statsec.p_value//"\n"
    paramstr += "obsmode        = "//obsmode.p_value//"\n"
    paramstr += "sci_ext        = "//sci_ext.p_value//"\n"
    paramstr += "var_ext        = "//var_ext.p_value//"\n"
    paramstr += "dq_ext         = "//dq_ext.p_value//"\n"
    paramstr += "key_detsec     = "//key_detsec.p_value//"\n"
    paramstr += "key_ccdsec     = "//key_ccdsec.p_value//"\n"
    paramstr += "key_datsec     = "//key_datsec.p_value//"\n"
    paramstr += "key_ccdsum     = "//key_ccdsum.p_value//"\n"
    paramstr += "key_obsmode    = "//key_obsmode.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value//"\n"
    paramstr += "fl_real        = "//fl_real.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    
    gloginit (l_logfile, "gmosaic", "gmos", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Load up the array of input file names.
    # Make sure the input files exist and that they are MEF files
    gemextn (l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, glogpars="", verbose=l_verbose)

    gemextn ("@"//tmpfile, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpinimg, logfile=l_logfile, glogpars="", verbose=l_verbose)

    nimages = gemextn.count
    delete (tmpfile, verify-, >& "dev$null")
    
    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maxfiles)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" images were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input images defined."
            status = 121
        } else if (nimages > maxfiles) {
            errmsg = "Maximum number of input images ["//str(maxfiles)//"] \
                has been exceeded."
            status = 121
        }
        
        glogprint (l_logfile, "gmosaic", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    }
    
    ######################################
    # Load up the image names to inimg[]
    i = 0
    scanfile = tmpinimg
    while (fscan (scanfile, filename) != EOF) {
        # check if image already mosaiced
        keyfound = ""
        hselect (filename//"[0]", "GMOSAIC", yes) | scan(keyfound)
        if (keyfound != "") {
            errmsg = "Image "//filename//" has already been gmosaic'ed"
            status = 121
            glogprint (l_logfile, "gmosaic", "status", type="error",
                errno=status, str=errmsg, verbose=l_verbose)
        } else {
            i += 1
            inimg[i] = filename
        }
    }
    scanfile = ""
    if (status != 0)
        goto clean
    
    if (i != nimages) {
        errmsg = "Error while counting the input images."
        status = 99
        glogprint (l_logfile, "gmosaic", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
    }
    
    ######################################
    # Check the input images for the required SCI, VAR, DQ names
    if ("" == l_sci_ext) {          # assume user asking to process RAW DATA
        if (l_fl_vardq) {
            status = 121
            errmsg = "No science extension name given indicating raw data. \
                Forcing fl_vardq to 'no'"
            glogprint (l_logfile, "gmosaic", "status", type="warning",
                errno=status, str=errmsg, verbose=l_verbose)
                
            l_fl_vardq = no
            status = 0          # warning issued, continuing.
        }
        sciext = ""
    } else {                        # working on PREPARED DATA
    
        # checking for img[l_sci_ext], img[l_var_ext], img[l_dq_ext]
        imglist = inimg[1]//"["//l_sci_ext//"]"
        for (i=2; i<=nimages; i=i+1) {
            imglist = imglist//","//inimg[i]//"["//l_sci_ext//"]"
        }        
        if (l_fl_vardq) {
            for (i=1; i<=nimages; i=i+1) {
                imglist = imglist//","//inimg[i]//"["//l_var_ext//"]"
                imglist = imglist//","//inimg[i]//"["//l_dq_ext//"]"
            }
        }
        
        gemextn (imglist, check="exists", process="none", index="",
            extname="", extversion="", ikparams="", omit="", replace="",
            outfile="dev$null", logfile=l_logfile, glogpars="", 
            verbose=l_verbose)

        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" missing "//l_sci_ext//", "//\
                l_var_ext//" and/or "//l_dq_ext//" extensions"
            status = 101
            glogprint (l_logfile, "gmosaic", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
        
        # set sciext, varext, and dqext to make the
        # rest of the code a bit more friendly (this is because we
        # need to allow for l_sci_ext="")
        # The extra "," allows "img["//sciext//"1]" whether l_sci_ext
        # is empty or not.
        
        sciext = l_sci_ext//","
        varext = l_var_ext//","
        dqext = l_dq_ext//","
    }
    
    
    ########################################
    # Load up the array of output file names
    #  Make they don't exist
    
    if (l_outimages != "")
        outputstr = l_outimages
    else if (l_outpref != "") {
        gemextn("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, glogpars="",
            verbose=l_verbose)
        outputstr = l_outpref//"@"//tmpoutimg
    } else {
        errmsg = "Neither output image name nor output prefix is defined"
        status = 121
        glogprint (l_logfile, "gmosaic", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    }
    
    gemextn(outputstr, check="", process="none", index="", extname="",
        extversion="", ikparam="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, glogpars="", verbose=l_verbose)
    delete (tmpoutimg, verify-, >& "dev$null")
    gemextn ("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, glogpars="", verbose=l_verbose)
    noutimages = gemextn.count
    delete (tmpfile, verify-, >& "dev$null")
    
    if ((gemextn.fail_count > 0) || (noutimages == 0) || \
        (noutimages != nimages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" image(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "Maximum number of output images exceeded: " \
                //str(maxfiles)
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "Different number of input ("//nimages//") and output \
                images ("//noutimages//")."
            status = 121
        }
        
        glogprint (l_logfile, "gmosaic", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
 
    } else {
        scanfile=tmpoutimg
        i=0
        while (fscan(scanfile, filename) != EOF) {
            i += 1
            outimg[i] = filename//".fits"
        }
        scanfile=""
        if (i != noutimages) {
            errmsg = "Error while counting the output images"
            status = 99
            glogprint (l_logfile, "gmosaic", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
    }
    
    
    ####################################
    #  Check on bpmfile if fl_fixpix=yes
    if (l_fl_fixpix && (no == access(l_bpmfile))) {
        errmsg = "Unable to find 'bpmfile'"
        status = 101
        glogprint (l_logfile, "gmosaic", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean      
    }

    
    #############################################################
    #############################################################
    
    # Start the work.  Loop through the input images.
    for (i=1; i<=nimages; i=i+1) {
        # define temp imgs used with this loop
        tmpimg = mktemp("tmpimg")   # do not delete until end of the loop
        tmpmdf = mktemp("tmpmdf")   
        tmpbpmtext = mktemp("tmpbpmtext")   # used only in fixpix block
        tmpbpm = mktemp("tmpbpm")   # used only in fixpix block
        
        # reset defaults
        currentimg = inimg[i]
        cleanthis = l_fl_clean
        fixthis = l_fl_fixpix

        sci_repvalue = orig_sci_repvalue
        var_repvalue = orig_var_repvalue
        dq_repvalue = orig_dq_repvalue
        
        glogprint (l_logfile, "gmosaic", "task", type="string",
            str="Input: "//currentimg//"  Output: "//outimg[i],
            verbose=l_verbose)
        
        
        #########################################################
        # Get all required info from the image
        
        # Set ishamamatsu flag
        imgets (currentimg//"[0]", "DETTYPE", >& "dev$null")
        if (imgets.value == "S10892-01")
            ishamamatsu = yes #hcode
        else
            ishamamatsu = no

        # Set gapvalue
        if (l_gap == "default") {
            if (ishamamatsu)
                gapvalue = 37 #hcode ##M Needs checking
            else
                gapvalue = 37
        } else
            gapvalue = int(l_gap)
        
        # Set statistic regions for use when cleaning the image
        if (ishamamatsu) { ##M Needs checking
            Xbingap = 36 #hcode 
            defstatsec[1]=2150 ; defstatsec[2]=3970 #hcode
            defstatsec[3]=100  ; defstatsec[4]=4000 #hcode
        } else {
            Xbingap = 36
            defstatsec[1]=2150 ; defstatsec[2]=3970
            defstatsec[3]=100  ; defstatsec[4]=4400
        }

        # Set instrument north or south 
        keyfound = ""
        hselect (currentimg//"[0]", "INSTRUME", yes) | scan(keyfound)
        if (keyfound == "GMOS-S")
            inst = 2
        else
            inst = 1    # assume GMOS-N
        
        l_struct = ""
        # Set the detector type here
        # GMOS-N
        if (inst == 1) {
            hselect (currentimg//"[0]", "DETTYPE", yes) | scan(l_struct)
            if (l_struct == "SDSU II CCD") { #Current EEV CCDs
                iccd = 1
            } else if (l_struct == "SDSU II e2v DD CCD42-90") { 
                # New e2vDD CCDs
                iccd = 2
            } else if (l_struct == "S10892-01") { # Hamamatsu CCDs
                iccd = 3
            }
        } else if (inst == 2) { # GMOS-S
            hselect (currentimg//"[0]", "DETECTOR", yes) | scan(l_struct)
            if (l_struct == "GMOS + Blue1 + new CCD1") {
                iccd = 2 # GMOS-S blue CCD upgrades
            } else {
                iccd = 1 # Original CCDs
            } 
        }
        
        # get observation date
        dateobs = ""
        hselect (currentimg//"[0]", "DATE-OBS", yes) | scan(dateobs)
        if (dateobs == "") {
            status = 131
            errmsg = "DATE-OBS keyword not found in "//currentimg
            glogprint (l_logfile, "gmosaic", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
        
        # Get binning
        Xbin=0 ; Ybin=0
#        hselect (currentimg//"["//sciext//"1]", l_key_ccdsum, yes) |\
#            translit ("STDIN", '"', """, delete+, collapse-) | \
#            scan (Xbin, Ybin)
        # To make PyRAF-compatible, we need to get rid of the translit.
        # translit will not delete the double-quotes for some reason.
        hselect (currentimg//"["//sciext//"1]", l_key_ccdsum, yes) |\
            scanf ('"%d %d"', Xbin, Ybin)
        if ((Xbin <= 0) || (Ybin <= 0)) {
            errmsg = "Unable to determine binning for "//currentimg
            status = 131
            glogprint (l_logfile, "gmosaic", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
        
        if (Xbin > 1 || Ybin > 1)
            isbinned = yes
        else
            isbinned = no
        
        if (Xbin > 1)
            thisimggap = nint(Xbingap/Xbin)
        else
            thisimggap = gapvalue          # Xbin = 1
        
        # Check if cleaning can be done.
        if (cleanthis) {
            # check against binning
            if ( (Xbin != Ybin) || Xbin > 2) {
                errmsg = "Setting fl_clean=no for "//currentimg//". \
                    No cleaning for bin>2 or Xbin!=Ybin"
                status = 121
                glogprint (l_logfile, "gmosaic", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
            }
       
            # check against OBSMODE
            keyfound = ""
            hselect (currentimg//"[0]", l_key_obsmode, yes) | scan (keyfound)
            if (keyfound == "") {
                errmsg = "Setting fl_clean=no for "//currentimg//". \
                    Cannot determine OBSMODE."
                status = 131
                glogprint (l_logfile, "gmosaic", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
            } else if (keyfound != l_obsmode) {
                cleanthis = no
            }
            
            # set flag, reset status
            if (status != 0) {
                cleanthis = no
                status = 0      # warning issued, continuing
            }
        }
        
        # fl_fixpix and fl_clean are mutually exclusive, clean takes precedence
        if (cleanthis && l_fl_fixpix) {
            errmsg = "Setting fl_fixpix = no for this image. Gaps cannot be \
                cleaned and interpolated."
            status = 121

            glogprint (l_logfile, "gmosaic", "status", type="warning",
                errno=status, str=errmsg, verbose=l_verbose)

            fixthis=no
            status = 0          # warning issued, continuing
        }


        ####
        # Custom ROI - Handling
        ndetroi = 1
        keypar (currentimg//"[0]", "DETNROI", silent+)
        if (!keypar.found) {
            errmsg = "DETNROI keyword not found - assuming non-custom ROI \
                image"
            status = 131
            glogprint (l_logfile, "gmosaic", "status", type="warning",
                errno=status, str=errmsg, verbose=l_verbose)
        } else {
            ndetroi = int(keypar.value)
        }

        run_gtile = no
        if (ndetroi > 1) {
            # Check to see if it's been GTILE'd already
            keypar (currentimg//"[0]", "GTILE", silent+)
            if (keypar.found) {
                # Check the type of tiling
                keypar (currentimg//"[0]", "GTILETYP", silent+)
                if (keypar.found) {
                    # 3 possible values - ROI APPENDED|ARRAY|DETECTOR
                    if (keypar.value == "ROI APPENDED") {
                        # Needs to be gtiled again
                        run_gtile = yes
                    } else if (keypar.value == "DETECTOR") {
                        errmsg = currentimg//" has been GTILE'd into \
                            the full detector"
                        status = 121
        
                        glogprint (l_logfile, "gmosaic", "status", \
                            type="error", errno=status, str=errmsg, verbose+)
                        goto clean
                    }
                } else {
                    errmsg = "Unable to determine type of tilling performed \
                        by GTILE"
                    status = 121
    
                    glogprint (l_logfile, "gmosaic", "status", \
                        type="error", errno=status, str=errmsg, verbose+)
                    goto clean
                }
            } else {
                run_gtile = yes
            }
        } # End of custom ROI checking

        if (run_gtile) {

            # Reset the replacement value for the chip gaps
            sci_repvalue = croi_scivalue
            var_repvalue = croi_varvalue
            dq_repvalue = croi_dqvalue

            glogprint (l_logfile, "gmosaic", "task", type="string", \
                str="Custom ROIs found", verbose=l_verbose)

            glogprint (l_logfile, "gmosaic", "status", type="fork", \
                fork="forward", child="gtile", verbose=l_verbose )

            gtile (currentimg, outimages=tmpgtiledimg, out_ccds="all", \
                ret_roi=yes, req_roi=0, fl_stats_only=no, fl_tile_det=no, \
                fl_app_rois=no, fl_pad=no, sci_padval=0., var_padval=0., \
                dq_padval=16., sci_fakeval=croi_scivalue, \
                var_fakeval=croi_varvalue., \
                dq_fakeval=croi_dqvalue, chipgap="default", \
                sci_ext=l_sci_ext, \
                var_ext=l_var_ext, dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, \
                key_detsec=l_key_detsec, key_ccdsec=l_key_ccdsec, \
                key_datasec=l_key_datsec, key_biassec="BIASSEC", \
                key_ccdsum=l_key_ccdsum, rawpath="", logfile=l_logfile, \
                fl_verbose=no)

            glogprint (l_logfile, "gmosaic", "status", type="fork", \
                fork="backward", child="gtile", verbose=l_verbose )

            if (gtile.status != 0) {
                errmsg = "GTILE returned a non-zero status"
                status = 121
    
                glogprint (l_logfile, "gmosaic", "status", \
                    type="error", errno=status, str=errmsg, verbose+)
                goto clean
            }

            currentimg = tmpgtiledimg
        }

        # Custom ROI - Handling END
        ####
        
        # Get number of science extensions
        nsci = 0
        hselect (currentimg//"[0]", "NSCIEXT", yes) | scan (nsci)
        if (nsci == 0) {    # NSCIEXT not in header => raw, count images
            fxhead (currentimg) | match ("IMAGE", "STDIN", stop-) | \
                count ("STDIN") | scan (nsci)
        }

        # Set the image's pixel type
        keyfound=""
        hselect (currentimg//"["//sciext//"1]", "i_pixtype", yes) |\
            scan (keyfound)
        if (keyfound == "11")
            pixtype = "l"   #long int
        else
            pixtype = "r"   #real
        
        # Convert to real to avoid rounding errors in geotran
        real_currentimg="r_"//currentimg
        
        if (l_fl_real && pixtype == "l") {
            if ("" == l_sci_ext) {
                gemarith (currentimg, "*", 1, result=real_currentimg, \
                    sci_ext="SCI", var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    mdf_ext="MDF", fl_vardq=l_fl_vardq, dims="default", \
                    intype="default", outtype="real", refim="default", \
                    rangecheck+, verbose-, lastout="", logfile=l_logfile)
            } else {
                gemarith (currentimg, "*", 1, result=real_currentimg, \
                    sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    mdf_ext="MDF", fl_vardq=l_fl_vardq, dims="default", \
                    intype="default", outtype="real", refim="default", \
                    rangecheck+, verbose-, lastout="", logfile=l_logfile)
            }
            currentimg=real_currentimg
            pixtype = "r"
        }
        
        # Get the MDF out of there for further use (do that before
        #  messing with 6-amp data since currentimg is changed in that
        #  block)
        
        tinfo (currentimg//".fits["//l_mdf_ext//"]", ttout-, >& "dev$null")
        if (tinfo.tbltype == "fits") {
            glogprint (l_logfile, "gmosaic", "task", type="string",
                str="Found a MDF, copying for later use.", verbose=l_verbose)
            tcopy (currentimg//".fits["//l_mdf_ext//"]", tmpmdf//".fits",
                verbose=no)
        }
        
        #################################################
        #  Must deal with n-amps
        #  Basically combine the amps to get full chips

        if (debug) {
            print ("    DEBUG:")
            print ("    DEBUG: instrument="//inst)
            print ("    DEBUG: detector="//iccd)
            print ("    DEBUG: date-obs="//dateobs)
            print ("    DEBUG: binning (x,y)=("//Xbin//","//Ybin//")")
            print ("    DEBUG: gap="//thisimggap)
            print ("    DEBUG: fixing activated="//fixthis)
            print ("    DEBUG: cleaning activated="//cleanthis)
            print ("    DEBUG: number of science extension="//nsci)
            print ("    DEBUG:")
        }
        
        # Get detector sections  (tell me where you come
        # from I'll figure out where you go)
        for (j=1; j<=nsci; j=j+1) {
            hselect (currentimg//"["//sciext//j//"]", l_key_detsec,
                yes) | scan (detsec)
            junk = fscanf (detsec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
            detx1[j] = tx1
            detx2[j] = tx2
            dety1[j] = ty1
            dety2[j] = ty2
            
            if (debug) {
                printf ("    DEBUG: DATASEC %d %d %d %d\n", 
                    detx1[j], detx2[j], dety1[j], dety2[j])
            }
            
            hselect (currentimg//"["//sciext//j//"]", l_key_datsec,
                yes) | scan (detsec)
            junk = fscanf (detsec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
            datx1[j] = tx1
            datx2[j] = tx2
            daty1[j] = ty1
            daty2[j] = ty2
            
            if (debug) {
                printf ("    DEBUG: DETSEC %d %d %d %d\n", 
                    datx1[j], datx2[j], daty1[j], daty2[j])
            }
            
            hselect (currentimg//"["//sciext//j//"]", l_key_ccdsec,
                yes) | scan (ccdsec)
            junk = fscanf (ccdsec, "[%d:%d,%d:%d]", tx1, tx2, ty1, ty2)
            ccdx1[j] = tx1
            ccdx2[j] = tx2
            ccdy1[j] = ty1
            ccdy2[j] = ty2
            
            if (debug) {
                printf ("    DEBUG: CCDSEC %d %d %d %d\n", 
                    ccdx1[j], ccdx2[j], ccdy1[j], ccdy2[j])
            }
            
            #initialize ampsorder to current non sorted amps order
            ampsorder[j] = j
        }
        
        # Now use that detsec and ccdsec info to figure out where each ext go.
        # Let's sort out the extensions then
        for (j=2; j<=nsci; j=j+1) {
            for (jj=1; jj<j; jj=jj+1) {
                if (detx1[j] < detx1[jj]) {
                    tmpint = ampsorder[jj]
                    ampsorder[jj] = ampsorder[j]
                    ampsorder[j] = tmpint
                }
            }
        }

        # Initiate number_of_ampsperccd
        for (j=1;j<=3;j+=1) {
            number_of_ampsperccd[j] = 0
        }

        # Determine the number of CCDs used. Makes the assumption it's either 
        # all 3 CCDs or all of CCD2 or a stamp. - MS
        # initiate counter number_of_ccds - there is always one - MS
        number_of_ccds = 1 
        for (j=1; j<=nsci; j+=1) {
            # Count the number of times a CCD edge ( x1 == 1 ) is found
            if (ccdx1[ampsorder[j]] == 1 && j != 1 && number_of_ccds == 1) {
                # First edge of CCD2
                # Calculate the number of amps on CCD1
                number_of_ampsperccd[1] = j - 1
                number_of_ccds += 1
            } else if (ccdx1[ampsorder[j]] == 1 && j != 1) { 
                # First edge of CCD3
                # Calculate the number of amps on CCD2
                number_of_ampsperccd[2] = j - 1 - number_of_ampsperccd[1]
                number_of_ccds += 1
            }

            if ( j==nsci && number_of_ccds == 3) {
                # Calculate the number of amps on ccd3
                number_of_ampsperccd[3] = nsci - ( number_of_ampsperccd[2] +\
                    number_of_ampsperccd[1])
            } else if ( j==nsci && number_of_ccds == 1 ) { 
                # CCD2 only / Stamp images
                # Calculate the number of amps on CCD2
                number_of_ampsperccd[2] = nsci
            }
        } # End of for (j=1; j<=nsci; j+=1) loop

        # Set reference extention - MS
        refextn = ampsorder[nint(number_of_ampsperccd[2] / 2.0) + \
            number_of_ampsperccd[1]]

        if (debug) {
            print ("Reference extension is: "//refextn)
        }

        # Obtain refextn's DATASEC Xaxis length; for updating WCS at the end 
        # - MS
        refextn_naxis1 = datx2[refextn] - (datx1[refextn] - 1)

        # Adjust CRPIX1 if required:
        # This catches all cases; so if trimmed CRPIX doesn't get updated
        # as DATASEC X1 will be 1 - MS
        if (datx1[refextn] == 1) {
            crpix1_adjust = 0
        } else {
            crpix1_adjust = (datx1[refextn] - 1)
        }

        # Set the range of ccds to loop over; used throughout code - MS
        if (number_of_ccds==1){
            chipstart = 2
            chipstop = 2
        } else {
            chipstart = 1
            chipstop = 3
        }

        if (debug) {
            printf ("    DEBUG: Amps order =")
            for (j=1; j<=nsci; j=j+1) {
                printf (" %d", ampsorder[j])
            }
            printf ("\n")
        }

        # Now let's tile the adjacent extensions.  This
        # will create _one_ data file with number_of_ccds image extensions
        # representing the newly tiled chips
        
        # Copy phu to tmpimg
        fxcopy (currentimg, tmpimg//".fits", groups="0", new_file=yes,
            verbose=no)            
        
        # Tile and append to tmpimg
        for (k=chipstart; k<=chipstop; k=k+1) {

            # tiling configuration
            nctile = number_of_ampsperccd[k]
            nltile = 1

            # Initiate parameters and temporary file names - MS
            pastesci="" ; pastevar="" ; pastedq=""
            tmptile = mktemp("tmptile")
            
            # Copy phu of current image to tmptile - the output file from
            # tiling the amps into one CCD - MS 
            fxcopy (currentimg, tmptile//".fits", groups="0", new_file=yes,
                verbose=no)

            # Set the extensions to tile together. - MS
            if (chipstart == chipstop) { # Only one CCD
                jstart = 1
            } else { # More than one CCD
                jstart = k + (k-1)*(number_of_ampsperccd[k]-1)
            }
            jstop = jstart + number_of_ampsperccd[k]

            for (j=jstart; j<jstop; j=j+1) {
                jj = ampsorder[j]
                datasec=""
                hselect (currentimg//"["//sciext//jj//"]", l_key_datsec,
                    yes) | scan(datasec)
                if (j==jstart) {
                    pastesci = currentimg//"["//sciext//jj//"]"//datasec
                    if (l_fl_vardq) {
                        pastevar = currentimg//"["//varext//jj//"]" \
                            //datasec
                        pastedq = currentimg//"["//dqext//jj//"]"//datasec
                   }
                } else {
                    pastesci = pastesci//","//\
                        currentimg//"["//sciext//jj//"]"//datasec
                     if (l_fl_vardq) {
                        pastevar = pastevar//","//\
                            currentimg//"["//varext//jj//"]"//datasec
                        pastedq = pastedq//","//\
                            currentimg//"["//dqext//jj//"]"//datasec
                   }
                }
            }
            
            if (debug) {
                print ("    DEBUG: Pasting "//pastesci)
                if (l_fl_vardq) {
                    print ("    DEBUG: Pasting "//pastevar)
                    print ("    DEBUG: Pasting "//pastedq)
                }
            }

            # If there is only one amp on a CCD there is no need to tile it
            # here. If that is the case the current image just gets copied 
            # to tmptile in the else to this if - MS
            if (number_of_ampsperccd[k] != 1) {
                imtile (pastesci, tmptile//"[1,append]", nctile, nltile, \
                    ncoverlap=0, trim_section="", missing_input="", \
                    start_tile="ll", row_order=yes, raster_order=no, \
                    median_section="", subtract=no, ncols=INDEF, \
                    nlines=INDEF, nloverlap=-1, opixtype=pixtype, \
                    ovalue=-1, verbose=l_verbose, >& tmplog)
            
                if (l_fl_vardq) {
                    imtile (pastevar, tmptile//"[2,append]", nctile, \ 
                        nltile, ncoverlap=0, trim_section="", \
                        missing_input="", start_tile="ll", row_order=yes, \
                        raster_order=no, median_section="", subtract=no, \
                        ncols=INDEF, nlines=INDEF, nloverlap=-1, \
                        opixtype=pixtype, ovalue=-1, verbose=l_verbose, \
                        >>& tmplog)
                    imtile (pastedq, tmptile//"[3,append]", nctile, nltile,
                        ncoverlap=0, trim_section="", missing_input="",
                        start_tile="ll", row_order=yes, raster_order=no,
                        median_section="", subtract=no, ncols=INDEF,
                        nlines=INDEF, nloverlap=-1, opixtype=pixtype, 
                        ovalue=-1, verbose=l_verbose, >>& tmplog)
                }
            } else { # Not tiled: just the data section is copied for each
                     # extension - MS
                imcopy (pastesci, tmptile//"[1,append]", verbose-, >>& tmplog)

                if (l_fl_vardq) {
                    imcopy (pastevar, tmptile//"[2,append]", verbose-, \
                        >>& tmplog)
                    imcopy (pastedq, tmptile//"[3,append]", verbose-, \
                        >>& tmplog)
                }
            }

            glogprint (l_logfile, "gmosaic", "engineering", type="file",
                str=tmplog, verbose=l_verbose)
            delete (tmplog, verify-, >& "dev$null")
            
            # Build the new datasec value
            tx2 = 0 ; ty2 = 0
            hselect (tmptile//"[1]", "i_naxis1,i_naxis2", yes) | \
                scan (tx2, ty2)
            datasec = "[1:"//str(tx2)//",1:"//str(ty2)//"]"
            
            # Set the headers for this tiled chip
            gemhedit (tmptile//"[1]", l_key_datsec, datasec, "")
            gemhedit (tmptile//"[1]", "EXTNAME", l_sci_ext, "")
            gemhedit (tmptile//"[1]", "EXTVER", k, "")

            if (l_fl_vardq) {
                gemhedit (tmptile//"[2]", l_key_datsec, datasec, "")
                gemhedit (tmptile//"[2]", "EXTNAME", l_var_ext, "")
                gemhedit (tmptile//"[2]", "EXTVER", k, "")
                    
                gemhedit (tmptile//"[3]", l_key_datsec, datasec, "")
                gemhedit (tmptile//"[3]", "EXTNAME", l_dq_ext, "")
                gemhedit (tmptile//"[3]", "EXTVER", k, "")
            }
                
            # Add this tile to the tmpimg
            if (l_fl_vardq) {
                pos = 3*(k - 1)
                groups="1,2,3"
            } else {
                pos = 1*(k - 1)
                groups="1"
            }

            # Reset pos if only 1 ccd - MS
            if (chipstart == chipstop) {
                pos = 0
            }

            fxinsert (tmptile, tmpimg//".fits["//pos//"]", groups=groups,
                verbose-)

            # clean up
            imdelete (tmptile, verify-, >& "dev$null")
        } # End of k loop to tile extensions for each ccd
        
        # assign tmp output to 'currentimg'.  Nothing else
        # special needs to be done for 6amps after that.
        currentimg = tmpimg
        
        #################################################
        # Bookkeeping, setting, etc. (keep separated from
        #  real work, even if it means two for-loops through
        #  the extensions; the sanity of the programmer is 
        #  at stake)
        
        for (k=chipstart; k<=chipstop; k=k+1) {
            if (l_fl_paste == no) {      # ie. will transform with geotran
            
                # Correct the rotation and the shift for the binning
                if (isbinned) {
                    lxrot[k] = lrotation[inst,iccd,k] * Xbin/Ybin
                    lyrot[k] = lrotation[inst,iccd,k] * Ybin/Xbin
                    
                    if (Xbin == 1)      #can be isbinned if Ybin!=1
                        geoxshift[k] = lxshift[inst,iccd,k]
                    else
                        geoxshift[k] = lxbshift[inst,iccd,k]/Xbin

                    geoyshift[k] = lyshift[inst,iccd,k]/Ybin
                        
                } else {
                    lxrot[k] = lrotation[inst,iccd,k]
                    lyrot[k] = lrotation[inst,iccd,k]
                    
                    geoxshift[k] = lxshift[inst,iccd,k]
                    geoyshift[k] = lyshift[inst,iccd,k]
                }
                
                if (l_fl_paste == no) {
                    printf ("Setting rotation Xrot[%d]=%9.6f, \
                        Yrot[%d]=%9.6f\n", k, lxrot[k], k, lyrot[k]) | \
                        scan(l_struct)
                    glogprint (l_logfile, "gmosaic", "engineering", 
                        type="string", str=l_struct, verbose=l_verbose)
                }
                
            }
            
            # Set chip names.  if paste then use currentimg,
            # otherwise we need tmp names for the geotran outputs.
            datasec=""

            # A way of minimising the code changes for it to work on RAW
            # CCD2 amplifiers only data - MS
            if (ishamamatsu && chipstart==chipstop && sciext=="" ) {
                 correctextn="1" # hardcoded and will therefore not work
                                # with vardq+ as this is for raw data
            } else {
                correctextn = sciext//k
            }

            hselect (currentimg//"["//correctextn//"]", l_key_datsec, yes) |\
                scan(datasec)
            if (l_fl_paste) {
                tmpchipsci[k] = currentimg//"["//correctextn//"]"//datasec
                if (l_fl_vardq) {
                    tmpchipvar[k] = currentimg//"["//varext//k//"]"//datasec
                    tmpchipdq[k] = currentimg//"["//dqext//k//"]"//datasec
                }
            } else {
                tmpgeoinsci[k] = currentimg//"["//correctextn//"]"//datasec
                tmpchipsci[k] = mktemp("tmpchip")
                if (l_fl_vardq) {
                    tmpgeoinvar[k] = currentimg//"["//varext//k//"]"//datasec
                    tmpchipvar[k] = mktemp("tmpchip")
                    tmpgeoindq[k] = currentimg//"["//dqext//k//"]"//datasec
                    tmpchipdq[k] = mktemp("tmpchip")
                }
            }
            
        }
        
        if (debug) {
            printf ("    DEBUG: sci_ext->")
            for (k=chipstart; k<=chipstop; k=k+1) {
                printf ("  %s", tmpchipsci[k])
            }
            printf ("\n")
            if (l_fl_vardq) {
                printf ("    DEBUG: var_ext->")
                for (k=chipstart; k<=chipstop; k=k+1) {
                    printf ("  %s", tmpchipvar[k]) 
                }
                printf ("\n")
                printf ("    DEBUG: dq_ext->")
                for (k=chipstart; k<=chipstop; k=k+1) {
                    printf ("  %s", tmpchipdq[k]) 
                }
                printf ("\n")
            }
        }


        #################################################
        # Ready to do work
        
        #################################################
        # Do geotran, if full transformation requested.
        if (!l_fl_paste) {
            # Only if more than one CCD
            if (number_of_ccds != 1) {
                for (k=chipstart; k<=chipstop; k=k+1) { 
                    geotran (tmpgeoinsci[k], tmpchipsci[k], "", "",
                        xshift=geoxshift[k],  yshift=geoyshift[k],
                        xrotation=lxrot[k], yrotation=lyrot[k], xmag=1, ymag=1,
                        xmin=INDEF, xmax=INDEF, ymin=INDEF, ymax=INDEF, 
                        ncols=INDEF, nlines=INDEF, fluxconserve=yes, 
                        nxblock=2048, nyblock=2048, interpolant=l_geointer, 
                        boundary="constant", constant=0, verbose=l_verbose, 
                        >& tmplog)

                    if (l_fl_vardq) {
                        geotran (tmpgeoinvar[k], tmpchipvar[k], "", "",
                            xshift=geoxshift[k],  yshift=geoyshift[k],
                            xrotation=lxrot[k], yrotation=lyrot[k], xmag=1, 
                            ymag=1, xmin=INDEF, xmax=INDEF, ymin=INDEF, 
                            ymax=INDEF, ncols=INDEF, nlines=INDEF, 
                            fluxconserve=yes, nxblock=2048, nyblock=2048, 
                            interpolant=l_geointer, boundary="constant", 
                            constant=0, verbose=l_verbose, >>& tmplog)

                        geotran (tmpgeoindq[k], tmpchipdq[k], "", "",
                            xshift=geoxshift[k],  yshift=geoyshift[k],
                            xrotation=lxrot[k], yrotation=lyrot[k], xmag=1, 
                            ymag=1, xmin=INDEF, xmax=INDEF, ymin=INDEF, 
                            ymax=INDEF,ncols=INDEF, nlines=INDEF, 
                            fluxconserve=yes, nxblock=2048, nyblock=2048, 
                            interpolant="nearest", boundary="constant", 
                            constant=0, verbose=l_verbose, >>& tmplog)
                    }
                    glogprint (l_logfile, "gmosaic", "engineering", 
                        type="file", str=tmplog, verbose=l_verbose)
                    delete (tmplog, verify-, >& "dev$null")
                } # End of loop over chips
            } else {
                # Warn user that not using geotran as only central CCD is used
                print ("GMOSAIC: Warning: Not using geotran as only data from \
                    the central CCD is being mosaiced", >& tmplog)
                glogprint (l_logfile, "gmosaic", "engineering", 
                    type="file", str=tmplog, verbose=l_verbose)
                delete (tmplog, verify-, >& "dev$null")

                # Copy files to filenames used in the next section
                imcopy (tmpgeoinsci[k], tmpchipsci[k], verbose+)
                if (l_fl_vardq) {
                    imcopy (tmpgeoinvar[k], tmpchipvar[k], verbose+)
                    imcopy (tmpgeoindq[k], tmpchipdq[k], verbose+)
                }
            } # End of if loop number_of_ccds != 1
        } # end of !fl_paste loop
        
        ###################################################
        # Tile it up!

        # Only tile if there are more than one CCD used
        if (number_of_ccds != 1) {
            # imtile input file list
            pastesci = tmpchipsci[chipstart]
            for (k=chipstart+1; k<=chipstop; k=k+1) {
                pastesci = pastesci//","//tmpchipsci[k]
            }

            if (l_fl_vardq) {
                pastevar = tmpchipvar[chipstart]
                pastedq = tmpchipdq[chipstart]
                for (k=chipstart+1; k<=chipstop; k=k+1) {
                    pastevar = pastevar//","//tmpchipvar[k]
                    pastedq = pastedq//","//tmpchipdq[k]
                }
            }

            # Copy PHU to new output file
            fxcopy (currentimg, outimg[i], groups="0", new_file=yes,
                verbose=no)

            # tiling configuration
            nctile = number_of_ccds
            nltile = 1
    
            # Do the tiling
            imtile (pastesci, outimg[i]//"[1,append]", nctile, nltile,
                ncoverlap=-thisimggap, trim_section="", missing_input="",
                start_tile="ll", row_order=yes, raster_order=no, 
                median_section="",
                subtract=no, ncols=INDEF, nlines=INDEF, nloverlap=-1,
                opixtype=pixtype, ovalue=sci_repvalue, verbose=l_verbose, \
                >& tmplog)

            if (l_fl_vardq) {
                imtile (pastevar, outimg[i]//"[2,append]", nctile, nltile,
                    ncoverlap=-thisimggap, trim_section="", missing_input="",
                    start_tile="ll", row_order=yes, raster_order=no, 
                    median_section="", subtract=no, ncols=INDEF, nlines=INDEF,
                    nloverlap=-1, opixtype=pixtype, ovalue=var_repvalue, 
                    verbose=l_verbose,
                    >>& tmplog)
                imtile (pastedq, outimg[i]//"[3,append]", nctile, nltile,
                    ncoverlap=-thisimggap, trim_section="", missing_input="",
                    start_tile="ll", row_order=yes, raster_order=no,
                    median_section="", subtract=no, ncols=INDEF, nlines=INDEF,
                    nloverlap=-1, opixtype=pixtype, ovalue=dq_repvalue, 
                    verbose=l_verbose,
                    >>& tmplog)
            }
            glogprint (l_logfile, "gmosaic", "engineering", type="file",
                str=tmplog, verbose=l_verbose)
            delete (tmplog, verify-, >& "dev$null")

        } else {

            if (!fl_vardq) {
                fxcopy (currentimg, outimg[i], groups="0-1", new_file=yes,
                     verbose=no)
            } else {
                fxcopy (currentimg, outimg[i], groups="0-3", new_file=yes,
                     verbose=no)
            }                    
        }

        # Name the extensions of the tiled image
        if (l_sci_ext == "")     # Presumably a raw frame, set name to SCI
            gemhedit (outimg[i]//"[1]", "EXTNAME", "SCI", "")
        else
            gemhedit (outimg[i]//"[1]", "EXTNAME", l_sci_ext, "")
        gemhedit (outimg[i]//"[1]", "EXTVER", 1, "")
        if (l_fl_vardq) {
            gemhedit (outimg[i]//"[2]", "EXTNAME", l_var_ext, "")
            gemhedit (outimg[i]//"[2]", "EXTVER", 1, "")
            gemhedit (outimg[i]//"[3]", "EXTNAME", l_dq_ext, "")
            gemhedit (outimg[i]//"[3]", "EXTVER", 1, "")
        }
        
        ####################################################
        # Fix pixel if requested
        
        if (fixthis) {
            glogprint (l_logfile, "gmosaic", "task", type="string",
                str="Interpolating across chip gaps using "//l_bpmfile,
                verbose=l_verbose)
                
            if (ishamamatsu) { # New chipgaps in pixels for Hamamatsu CCDs -MS
                fields (l_bpmfile,"1,2,3,4",lines="5,6", >& tmpbpmtext) 
            } else { # Old chipgaps in pixels for EEV2 CCDs - MS
                fields (l_bpmfile,"1,2,3,4",lines="2,3", >& tmpbpmtext) 
            }
            if (isbinned) {
                if (Xbin>1) {
                    tcalc (tmpbpmtext, "c1", "int((c1-1.)/"//str(Xbin)// \
                        "+0.5)")
                    tcalc (tmpbpmtext, "c2", "int((c2-1.)/"//str(Xbin)// \
                        "+0.5)")
                }
                tcalc (tmpbpmtext, "c4", "int(c4/"//str(Ybin)//")")
            }
            badpiximage (tmpbpmtext, outimg[i]//"["//sciext//"1]", \
                tmpbpm//".pl", goodvalue=0, badvalue=1)

            proto.fixpix (outimg[i]//"["//sciext//"1]", tmpbpm//".pl", \
                linterp="INDEF", cinterp="INDEF", verbose=no, pixels=no, \
                >>& "dev$null")

            if (l_fl_vardq) {
                # The imexpr call below wasn't working for variance 
                # just used same bpm for same location of chip gaps
                # in the above fixpix call - JH 12.17.07
                #imexpr("(b>0) ? "//str(0)//" : a",
                #    outimg[i]//"["//varext//"1,overwrite]",
                #    outimg[i]//"["//varext//"1]", tmpbpm, verbose-)

                # The imexpr call below is randomly failing due to an "extname 
                # and/or extver value not found ('DQ,1')" error. However, the 
                # temporary file that is left behind when the task fails at 
                # this point is in tact as has a valid DQ extension. Putting 
                # the fixpix call after the imexpr call seems to have solved 
                # the problem. EH

                imexpr("(b>0) ? "//str(1)//" : a", \
                    outimg[i]//"["//dqext//"1,overwrite]", \
                    outimg[i]//"["//dqext//"1]", tmpbpm//".pl", \
                    outtype="short", verbose-)  

                proto.fixpix (outimg[i]//"["//varext//"1]", tmpbpm//".pl", \
                    linterp="INDEF", cinterp="INDEF", verbose=no, \
                    pixels=no, >>& "dev$null")

            }

            delete (tmpbpmtext, verify-, >& "dev$null")
            imdelete (tmpbpm//".pl", verify-, >>& "dev$null")
        }
        
        ####################################################
        # Clean data if requested
        
        if (cleanthis) {
            # Build the BPM name
            tmpbpm = "gmos$data/gmos-"
            if (inst == 1) {
                tmpbpm = tmpbpm//"n_"
            } else {
                tmpbpm = tmpbpm//"s_"
            }

            tmpbpm = tmpbpm//"bpm_"
            
            keypar (outimg[i]//"[0]", "DETTYPE", silent+)
            dettype = str(keypar.value)
            if (dettype == "SDSU II CCD") { #Current EEV CCDs
                # EEV
                tmpbpm = tmpbpm//"EEV"
            } else if (dettype == "SDSU II e2v DD CCD42-90") { 
                # e2vDD
                tmpbpm = tmpbpm//"e2v"
            } else if (dettype == "S10892-01") { # Hamamatsu CCDs
                # Hamamatsu
                tmpbpm = tmpbpm//"HAM"
            }

            tmpbpm = tmpbpm//"_"//str(Xbin)//str(Ybin)//"_"

            keypar (outimg[i]//"[0]", "NAMPS", silent+)
            namps = int(keypar.value) * 3
            tmpbpm = tmpbpm//str(namps)//"amp_"

            keypar (outimg[i]//"[0]", "DETECTOR", silent+)
            detector = str(keypar.value)
            if (detector == "GMOS + Blue1 + new CCD1") {
                # GMOS-S; newer EEV CCDs
                tmpbpm = tmpbpm//"v2"
            } else {
                iccd = 1 # Original CCDs
                tmpbpm = tmpbpm//"v1"
            }

            tmpbpm = tmpbpm//"_mosaic.pl"

            if (l_statsec == "default") {
                thisstatsec = "["//str(int(defstatsec[1]/Xbin))//":"//
                    str(int(defstatsec[2]/Xbin))//","//
                    str(int(defstatsec[3]/Ybin))//":"//
                    str(int(defstatsec[4]/Ybin))//"]"
            } else {
                thisstatsec = l_statsec
            }

            imstat (outimg[i]//"["//sciext//"1]"//thisstatsec, fields="midpt",
                lower=INDEF, upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF,
                binwidth=0.1, format-, cache-) | scan (medsky)

            glogprint (l_logfile, "gmosaic", "task", type="string",
                str="Cleaning areas outside imaging field, cleaning \
                value = "//str(medsky), verbose=l_verbose)
            imexpr ("(b>0) ? "//str(medsky)//" : a",
                outimg[i]//"["//sciext//"1,overwrite]",
                outimg[i]//"["//sciext//"1]",
                tmpbpm, verbose-)
            if (l_fl_vardq) {
                imexpr ("(b>0) ? "//str(0)//" : a",
                    outimg[i]//"["//varext//"1,overwrite]",
                    outimg[i]//"["//varext//"1]",
                    tmpbpm, verbose-)
                imexpr ("(b>0) ? "//str(1)//" : a",
                    outimg[i]//"["//dqext//"1,overwrite]",
                    outimg[i]//"["//dqext//"1]",
                    tmpbpm, verbose-)
                
            }
            gemhedit (outimg[i]//"[0]", "GMOSCLVA", medsky, \
                "GMOSAIC cleaning value")
        }
        
        ####################################################
        # Headers
        
        # Get the WCS from the reference extension in the original image
        #
        # CTYPE1  = 'RA---TAN'           / R.A. in tangent plane projection
        # CRPIX1  =     1044.60927497064 / Ref pix of axis 1
        # CRVAL1  =     308.611677174582 / RA at Ref pix in decimal degrees
        # CTYPE2  = 'DEC--TAN'           / DEC. in tangent plane projection
        # CRPIX2  =     2289.53839638225 / Ref pix of axis 2
        # CRVAL2  =     28.2797249736673 / DEC at Ref pix in decimal degrees
        # CD1_1   = -2.0204567874974E-05 / WCS matrix element 1 1
        # CD1_2   =  1.8157912081056E-08 / WCS matrix element 1 2
        # CD2_1   = 2.39303378045856E-08 / WCS matrix element 2 1
        # CD2_2   = 2.02198294404948E-05 / WCS matrix element 2 2
        # CCDSUM  = '2 2     '           / CCD sum

        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CTYPE", "STDIN", stop-, > tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CRPIX", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CRVAL", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CD1_1", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CD1_2", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CD2_1", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CD2_2", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("CCDSUM", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("GAIN", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("RDNOISE", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("RADECSYS", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("EQUINOX", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("MJD-OBS", "STDIN", stop-, >> tmpfile)
        imhead (inimg[i]//"["//sciext//refextn//"]", long+) |\
            match ("DISPAXIS", "STDIN", stop-, >> tmpfile)
        
        mkheader (outimg[i]//"[0]", tmpfile, append+, verbose-)
        mkheader (outimg[i]//"["//sciext//"1]", tmpfile, append+, verbose-)
        if (l_fl_vardq) {
            mkheader (outimg[i]//"["//varext//"1]", tmpfile, append+, verbose-)
            mkheader (outimg[i]//"["//dqext//"1]", tmpfile, append+, verbose-)
        }
        delete (tmpfile, verify-, >& "dev$null")

        keyfound = ""
        hselect (outimg[i]//"["//sciext//"1]", "CRPIX1", yes) | scan (keyfound)
        if (keyfound == "") {
            errmsg = "Image "//outimg[i]//" has no WCS"
            status = 131
            glogprint (l_logfile, "gmosaic", "status", type="warning",
                errno=status, str=errmsg, verbose=l_verbose)
            status = 0      # warning issued, continuing
        } else {
            # Update the WCS crpix1 value
            if (debug) {
                print ("CRPIX1 of REFEXTN: "//real(keyfound)//\
                    "\nChip gap value used "//real(thisimggap)//\
                    "\nTotal length of amplifiers upto but not "//\
                    "including the REFEXTN: "//\
                    real(nint(((number_of_ampsperccd[2] / 2.0) - 1) + \
                        number_of_ampsperccd[1]) * refextn_naxis1)//\
                    "\nCRPIX1 adjustment value: "//real(crpix1_adjust))
            }

            # NOTE: CRPIX 1 is relative to the pixels in the image (including 
            #       BIASSEC in all data, so add crpix1_adjust - MS

            if (number_of_ccds == 1) {
                gap_to_add = 0
            } else {
                gap_to_add = thisimggap
            }

            crpix1 = real(keyfound) + real(gap_to_add) + \
                real(nint(((number_of_ampsperccd[2] / 2.0) - 1) + \
                    number_of_ampsperccd[1]) * refextn_naxis1) - \
                real(crpix1_adjust)

            gemhedit (outimg[i]//"[0]", "CRPIX1", crpix1, "")
            gemhedit (outimg[i]//"["//sciext//"1]", "CRPIX1", crpix1, "")
            if (l_fl_vardq) {
                gemhedit (outimg[i]//"["//varext//"1]", "CRPIX1", crpix1, "")
                gemhedit (outimg[i]//"["//dqext//"1]", "CRPIX1", crpix1, "")
            }
        }

        # Done with the WCS. now with the final edits to the PHU
        gemdate ()
        gemhedit (outimg[i]//"[0]", "GMOSAIC", gemdate.outdate, \
            "UT Time stamp for GMOSAIC")
        gemhedit (outimg[i]//"[0]", "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI")
        if (!l_fl_paste)
            gemhedit (outimg[i]//"[0]", "GMSINTER", l_geointer, \
                "Interpolant used by GMOSAIC")
        
        # Add MDF if it exists
        if (access(tmpmdf//".fits")) {
            fxinsert (tmpmdf//"[1]", outimg[i]//"[0]", groups="", verbose-)
            delete (tmpmdf//".fits", verify-, >& "dev$null")
        }
        
        # Count the extensions and update PHU
        gemhedit (outimg[i]//"[0]", "NSCIEXT", 1, \
            "Number of science extensions")
        gemextn (outimg[i], check="", process="expand", index="1-",
            extname="", extver="", ikparam="", omit="", replace="",
            outfile="dev$null", logfile="", glogpars="", verbose=no)
        gemhedit (outimg[i]//"[0]", "NEXTEND", gemextn.count, \
            "Number of extensions", delete-)
       ####################################################
        # Done with this image, cleaning and moving on to
        # the next one.

        delete (tmpmdf//".fits", verify-, >& "dev$null")
        delete (tmpgtiledimg, verify-, >& "dev$null")
        delete (real_currentimg//".fits", verify-, >& "dev$null")
        imdelete (tmpimg, verify-, >& "dev$null")
        for (k=chipstart; k<=chipstop; k=k+1) {
            imdelete (tmpchipsci[k], verify-, >& "dev$null")
            imdelete (tmpchipvar[k], verify-, >& "dev$null")
            imdelete (tmpchipdq[k], verify-, >& "dev$null")
        }
        
        ngood = ngood+1
    }  #end for-loop through the input images


clean:
    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")
    delete (tmpinimg, verify-, >& "dev$null")
    delete (tmpoutimg, verify-, >& "dev$null")
    delete (tmpgtiledimg, verify-, >& "dev$null")

exit:
    if (status == 0) {
        logmsg = "All "//str(nimages)//" images successfully mosaiced."
        success = yes
    } else {
        logmsg = str(ngood)//" of "//str(nimages)//" successfully mosaiced."
        success = no
    }
    glogprint (l_logfile, "gmosaic", "status", type="string",
        str=logmsg, verbose=l_verbose)
    glogclose (l_logfile, "gmosaic", fl_success=success, verbose=l_verbose)
   
exitnow:
    ;  
    
end
