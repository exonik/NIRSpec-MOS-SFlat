# Utility functions

import numpy as np
import pandas as pd


# LXML
import lxml
from lxml import etree


# Plot functions
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# JWST Pipeline
import jwst
from jwst import datamodels
from jwst.pipeline import calwebb_spec2
from jwst.assign_wcs import nirspec
from jwst.assign_wcs import AssignWcsStep


# Astropy and other spectroscopy utilities
import astropy
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
import astropy.io.fits as fits
from astropy.timeseries import LombScargle
from astropy.utils.data import download_file





def list2np(in_list):
    l2 = [int(v) for v in in_list] # convert to np array of long ints
    pids = np.array(l2)
    return pids




def display_characters(num):
    # Returns a char for a given number <=26
    import string
    i = 0
    for char in string.ascii_lowercase:
        i += 1
        if (i == num):
            letter_id = char
            #print(letter_id)
            #print(char, end=' ')
    return letter_id


#print(display_characters(2))




def jwst_activity_naming_nn(anbr):
    # Computes the MAST activity number with
    # base 36 using actual activity number
    if anbr <= 9:
        jnbr = anbr

    if ((anbr >= 10) & ( anbr <= 35)):
        #print(anbr-10+1)
        jnbr = display_characters(anbr-10+1)

    if (anbr > 35):
        #print('Outside range for naming activity!!!')
        idxx = anbr-35 - 1
        b = 10 + np.arange(90)
        jnbr = b[idxx]

    return jnbr




# test
# a = 1+np.arange(41)
#for i in range(len(a)):
#    print( str(a[i]).zfill(2), ' ', jwst_activity_naming_nn(a[i]) )
#






def readcol(path_with_file):
    # Reads one-column txt
    # Input: file with a path
    # Output: list
    # N.Nikolov, STScI, 2023

    f = open(path_with_file, "r")
    lines = f.readlines()

    colws = []

    for x in lines:
        colws.append(x.split(' ')[0])
    f.close()

    col1 = []

    for sub in colws:
        col1.append(sub.replace("\n", ""))
    #print(col1)

    return col1








# Add functions that parse the individual sections of the APT



def nn_parse_apt(apt_pid, path_to_file, tagname_of_node_to_parse):
    # Reads a SINGLE APT file and parses Information
    # given a user provided node
    # Input: APT file, path and
    # Output: df with the contents
    # N.Nikolov, STScI, 2024

    web = "{http://www.stsci.edu/JWST/APT}"

    tree = etree.parse(path_to_file)
    root = tree.getroot()

    # Define a data frame to keep all results
    df_cols = [
    "ProposalID",
    "Label",
    "ObsNumber",
    "TargetID",
    "Target_Number",
    "Label2",
    "Instrument",
    "CoordinatedParallel",
    "OperatingMode",
    "Subarray",
    "MsaConfigFile",
    "ReadoutPattern",
    "Groups",
    "Integrations",
    "EtcID",
    "Lamp",
    "Grating",
    "ActivityNum",
    "ObsVisit",
    "ScienceDuration"]

    # Define an empty list to append list of dictionaries
    rows = []


    for node in root:
        #print(node.tag)

        #-------------------  2  ---------------------
        # Look inside the Targets
        if (node.tag == web+'DataRequests'):

            #print("Hello!")
            #{http://www.stsci.edu/JWST/APT}ProposalInformation
            #{http://www.stsci.edu/JWST/APT}VisitStatuses
            #{http://www.stsci.edu/JWST/APT}Targets
            #{http://www.stsci.edu/JWST/APT}PureParallelSlots
            #{http://www.stsci.edu/JWST/APT}DataRequests
            #{http://www.stsci.edu/JWST/APT}LinkingRequirements
            #{http://www.stsci.edu/JWST/APT}ProposalHistory
            #{http://www.stsci.edu/JWST/APT}ToolData

            # Loop the tags to identify Title
            for element in node:
                #print('OG', element.tag)

                # Enter the Observation Group
                if (element.tag == web+'ObservationGroup'):
                    dr_label = 'NONE'
                    obs_target = 'NONE'
                    obs_number = -1
                    obs_target_id = -1
                    #print('dr_label', dr_label)

                    # Show the contents of the Observation Group
                    #print("Total number of children: ", len(element))
                    for elem2 in element:
                        #print('Start of TAG: ', elem2.tag)

                        # Check the contents of "Label"
                        if (elem2.tag == web+'Label'):
                            dr_label = elem2.text
                            #print(dr_label)

                        # Check the contents of "Comments"
                        if (elem2.tag == web+'Comments'):
                            dr_comm = elem2.text
                            #print(dr_comm)

                        # Check the contents of "Observation"
                        if (elem2.tag == web+'Observation'):
                            #dr_autotarget = elem2.attrib.get("AutoTarget")
                            #print("Tname: ", elem2.tag)
                            #print(elem2.attrib.get("AutoTarget"))
                            #print("Start of an observation!")



                            obs_sci_dur = -1
                            obs_coord_par = False
                            obs_number = -1
                            obs_visit = -1
                            obs_target = None

                            for elem3 in elem2:
                                #print(elem3.tag)


                                # Print Coordinated Parallel information
                                if (elem3.tag == web+'CoordinatedParallel'):
                                    obs_coord_par = elem3.text
                                    #print('Coordinated Parallel', obs_coord_par)


                                if (elem3.tag == web+'Number'):
                                    obs_number = elem3.text
                                    #print('Observation number: ', obs_number)



                                if (elem3.tag == web+'TargetID'):
                                    obs_target = elem3.text
                                    #print('obs_target ', obs_target)#, type(obs_target))


                                if (elem3.tag == web+'Label'):
                                    obs_label2 = elem3.text
                                    #print('Obs. Label: ', obs_label2)

                                if (elem3.tag == web+'Instrument'):
                                    obs_instrument = elem3.text
                                    #print('Obs. Instrument: ', obs_instrument)



                                if (elem3.tag == web+'Template'):
                                    #print(elem3[0][0].text)

                                    i_templ = 0
                                    for templ in elem3:
                                        #print(i_templ, ' template')
                                        i_templ = i_templ + 1
                                        #print(templ.tag)
                                        #print(templ.text)
                                        #print(templ[0].tag)
                                        templ_suffix = '/Template'
                                        obs_templ = templ.tag
                                        #print('obs_templ', obs_templ)


                                        # Define initial value for variables to be searched
                                        #obs_subarray = -1
                                        #obs_rop      = -1
                                        #obs_groups   = -1
                                        #obs_opt_elem = 'X'
                                        obs_subarray = 'X'
                                        #obs_msa_config_file = 'X'
                                        obs_rot_pattern = 'X'
                                        obs_groups = 'X'
                                        obs_integr = 'X'
                                        obs_etcid = 'X'
                                        obs_lamp = 'X'
                                        obs_grating = 'X'

                                        # Search only for Lamp exposures
                                        web3 = web[0:30] + templ_suffix + '/NirspecInternalLamp}'
                                        case3 = web3 + 'NirspecInternalLamp'

                                        if (templ.tag == case3):
                                            #print('1111')
                                            # Loop he contents of lamp exposures
                                            for templ_attr in templ:
                                                if (templ_attr.tag == web3+'ExposureList'):
                                                    # Enter the exposure list
                                                    # count the number of exposures
                                                    count_exposure = 0
                                                    for exp_id in templ_attr:
                                                        count_exposure += 1
                                                        #print("exp_id: ", exp_id.tag)

                                                        # Loop the exposure list and obtain information
                                                        for exp_cnts in exp_id:
                                                            #print("exp_cnts:", exp_cnts.tag, exp_cnts.text)

                                                            if (exp_cnts.tag == web3+'OperatingMode'):
                                                                obs_op_mode = exp_cnts.text

                                                            if (exp_cnts.tag == web3+'Subarray'):
                                                                obs_subarray = exp_cnts.text

                                                            if (exp_cnts.tag == web3+'MsaConfigFile'):
                                                                obs_msa_config_file = exp_cnts.text
                                                                #print(obs_msa_config_file)

                                                            if (exp_cnts.tag == web3+'ReadoutPattern'):
                                                                obs_rot_pattern = exp_cnts.text

                                                            if (exp_cnts.tag == web3+'Groups'):
                                                                obs_groups = exp_cnts.text

                                                            if (exp_cnts.tag == web3+'Integrations'):
                                                                obs_integr = exp_cnts.text

                                                            if (exp_cnts.tag == web3+'EtcId'):
                                                                obs_etcid = exp_cnts.text

                                                            if (exp_cnts.tag == web3+'Lamp'):
                                                                obs_lamp = exp_cnts.text

                                                            if (exp_cnts.tag == web3+'Grating'):
                                                                obs_grating = exp_cnts.text

                                                        # Append list of dictionaries
                                                        rows.append({
                                                        "ProposalID":           apt_pid,
                                                        "Label":                dr_label,
                                                        "ObsNumber":            obs_number,
                                                        "TargetID":             obs_target,
                                                        "Target_Number":        obs_target_id,
                                                        "Label2":               obs_label2,
                                                        "Instrument":           obs_instrument,
                                                        "CoordinatedParallel":  obs_coord_par,
                                                        "OperatingMode":        obs_op_mode,
                                                        "Subarray":             obs_subarray,
                                                        "MsaConfigFile":        obs_msa_config_file,
                                                        "ReadoutPattern":       obs_rot_pattern,
                                                        "Groups":               obs_groups,
                                                        "Integrations":         obs_integr,
                                                        "EtcID":                obs_etcid,
                                                        "Lamp":                 obs_lamp,
                                                        "Grating":              obs_grating,
                                                        "ActivityNum":          count_exposure,
                                                        "ObsVisit":             obs_visit,
                                                        "ScienceDuration":      obs_sci_dur})

                                                        #print("Obs Numb: ", obs_number)
                                                        #print("Obs Visit: ", obs_visit)



                                if (elem3.tag == web+'Visit'):
                                    #print("Echo!")
                                    obs_visit = elem3.attrib.get("Number")
                                    #print('Observation visit: ', obs_visit)
                                    for d1 in rows:
                                        if (d1['ObsNumber'] == obs_number):
                                            d1['ObsVisit'] = obs_visit

                                # Print the Science duration
                                if (elem3.tag == web+'ScienceDuration'):
                                    obs_sci_dur = elem3.text
                                    #print('Science Duration: ', obs_sci_dur)
                                    for d2 in rows:
                                        if (d2['ObsNumber'] == obs_number):
                                            d2['ScienceDuration'] = obs_sci_dur





    out_df = pd.DataFrame(rows, columns = df_cols)
    return out_df












# SFlat functions


def plot_mos_msa(im1, im2, msa_im, im1_nm, im2_nm, msa_nm, im_lo, im_hi, msa_lo, msa_hi, fig_sz, path_and_name, grating_id):
    # Plot the data
    # Check toward the end here: https://matplotlib.org/stable/users/explain/axes/constrainedlayout_guide.html

    fig = plt.figure(constrained_layout=True, figsize=fig_sz)
    #(subfig_l, subfig_r) = fig.subfigures(nrows=1, ncols=2)
    (subfig_l, subfig_m, subfig_r) = fig.subfigures(nrows=1, ncols=3)

    # Plot the first two figures
    #axes_l = subfig_l.subplots(nrows=1, ncols=2, sharey=True)
    axes_l = subfig_l.subplots(nrows=1, ncols=1, sharey=True)
    im = axes_l.imshow(im1,
                   interpolation='None',
                   aspect='auto',
                   cmap='inferno',
                   origin='lower',
                   vmin=im_lo,
                   vmax=im_hi)
    axes_l.set_xlabel('x-column', fontsize=15)
    axes_l.set_ylabel('y-row', fontsize=15)
    axes_l.set_box_aspect(1)
    axes_l.set_title(im1_nm + '\n fmin=' + str(round(im_lo,1)) + ', fmax=' + str(round(im_hi,1)))

    axes_m = subfig_m.subplots(nrows=1, ncols=1, sharey=True)
    im = axes_m.imshow(im2,
                   interpolation='None',
                   aspect='auto',
                   cmap='inferno',
                   origin='lower',
                   vmin=im_lo,
                   vmax=im_hi)
    axes_m.set_xlabel('x-column', fontsize=15)
    axes_m.set_ylabel('y-row', fontsize=15)
    axes_m.set_box_aspect(1)
    axes_m.set_title(im2_nm + '\n ' + grating_id)


    axes_r = subfig_r.subplots(nrows=1, ncols=1, sharey=True)
    #for ax in axes_r:
    #    im = ax.imshow((msa_image), vmin=f_lo_2, vmax=f_hi_2)
    im = axes_r.imshow(msa_im,
                       interpolation='None',
                       aspect='auto',
                       cmap='inferno',
                       origin='lower',
                       vmin=msa_lo,
                       vmax=msa_hi)
    axes_r.set_xlabel('x-column', fontsize=15)
    axes_r.set_ylabel('y-row', fontsize=15)
    axes_r.set_title(msa_nm + '\n fmin=' + str(round(msa_lo,1)) + ', fmax=' + str(round(msa_hi,4)))
    axes_r.set_box_aspect(1/1)

    # Save the figure
    #plt.savefig(path_and_name,bbox_inches='tight',pad_inches=0.25,dpi=(200))

    pp = PdfPages(path_and_name)
    plt.savefig(pp,format='pdf',bbox_inches='tight',pad_inches=0.25)
    plt.close()
    pp.close()




def shutter_correct(msa_shutter):
    #----------------------------------
    # JWST NIRSpec MOS convert the MSA_SHUTTER
    # to match the data
    # First version: 2023
    # N.Nikolov 2024
    # STScI - Baltimore, USA
    #----------------------------------
    # Conda environment: sflat_tools
    #----------------------------------
    msa_shutter = np.rot90(msa_shutter, k=3)
    # after rotation the array is x:0->730 and y:0:342
    #------------------------------------------------------
    '''
    plt.figure(figsize= (20,20))
    f_lo_2 = -1. # np.nanmin(image) # lower limit for flux
    f_hi_2 =  1. # np.nanmax(image) # upper limit for flux
    im2_name = ' '

    plt.imshow(msa_shutter,
               interpolation='None',
               cmap='inferno',
               origin='lower',
               clim=(f_lo_2, f_hi_2))
    plt.xlabel('x-column', fontsize=15)
    plt.ylabel('y-row', fontsize=15)
    cb1 = plt.colorbar(label=r'Flag: '+str(f_lo_2)+' to '+str(f_hi_2))
    plt.title('MSA shutter image:\n' +im2_name, fontsize=12)

    plt.show()
    '''
    #------------------------------------------------------
    Q2 = msa_shutter[171:   , 0  :365]
    Q4 = msa_shutter[  0:171, 0  :365]
    Q1 = msa_shutter[171:   , 365:   ]
    Q3 = msa_shutter[  0:171, 365:   ]

    # Each quadrant is also flipped upside down
    Q2 = np.flipud(Q2)
    Q4 = np.flipud(Q4)
    Q1 = np.flipud(Q1)
    Q3 = np.flipud(Q3)

    Q12 = np.concatenate((Q2, Q1), axis=0) # axis=0 puts 1 above 2; axis=1 1 next to 2
    Q34 = np.concatenate((Q4, Q3), axis=0)
    new_image = np.concatenate((Q34, Q12), axis=1)

    return new_image



def find_shutter_coords(msa_image):
    # Convert to a function for long slit
    msa_xlength = len(msa_image[0,:])
    msa_ylength = len(msa_image[:,0])
    #print(msa_xlength, msa_ylength)

    idx = np.where(msa_image != 0)
    n_idx = len(idx[0])
    idx_xx = idx[1]
    idx_yy = idx[0]
    #print("idx0: ", idx[0])
    #print("idx1: ", idx[1])
    #print("NE 0 len: ", len(idx))

    msa_shutter_quadrant = np.zeros(n_idx, dtype=int) #msa_ylength
    msa_shutter_row      = np.zeros(n_idx, dtype=int)
    msa_shutter_col      = np.zeros(n_idx, dtype=int)
    msa_shutter_row_Qcrd = np.zeros(n_idx, dtype=int)

    for i in range(n_idx):
        msa_shutter_row[i] = msa_ylength - idx_yy[i]
        msa_shutter_col[i] = msa_xlength - idx_xx[i]

        if (idx_yy[i] > 170):
            msa_shutter_row_Qcrd[i] = np.abs(342 - idx_yy[i])
        else:
            msa_shutter_row_Qcrd[i] = np.abs(171 - idx_yy[i])

        # Q4
        if ((idx_yy[i] >= 0 and idx_yy[i] <= 170) and (idx_xx[i] >= 0 and idx_xx[i] <= 364)):
            msa_shutter_quadrant[i] = 4

        #Q3
        if ((idx_yy[i] >= 171 and idx_yy[i] < msa_ylength) and (idx_xx[i] >= 0 and idx_xx[i] <= 364)):
            msa_shutter_quadrant[i] = 3

        # Q2
        if ((idx_yy[i] >= 0 and idx_yy[i] <= 170) and (idx_xx[i] >= 365 and idx_xx[i] < msa_xlength)):
            msa_shutter_quadrant[i] = 2

        #Q3
        if ((idx_yy[i] >= 171 and idx_yy[i] < msa_ylength) and (idx_xx[i] >= 365 and idx_xx[i] < msa_xlength)):
            msa_shutter_quadrant[i] = 1

    return msa_shutter_quadrant, msa_shutter_row, msa_shutter_col, msa_shutter_row_Qcrd




def print_meta(dm_im):
    print(dm_im.meta.filename,
          str(dm_im.meta.exposure.count+1).zfill(2), # Exposure count in vist (Starts from 0 in the header!!!!)
          str(dm_im.meta.exposure.ngroups).zfill(2),
          dm_im.meta.exposure.type, # can be NRS_LAMP, NRS_AUTOWAVE or NRS_AUTOFLAT
          dm_im.meta.instrument.lamp_mode, # can be MSASPEC, IFU, FIXEDSLIT, BRIGHTOBJ, GRATING-ONLY, NONE, IMAGE
          dm_im.meta.instrument.lamp_state,
          dm_im.meta.observation.date,
          dm_im.meta.observation.time,
          dm_im.meta.instrument.filter,
          str(dm_im.meta.instrument.grating).zfill(00000000),
          dm_im.meta.instrument.msa_configuration_id)




def msa_modify_v2(path_to_msa_file, msa_file_id, prog_id, use_subset_of_slitlets, nsl_start, nsl, disperser):
    # Modifies existing MSA file
    # written by NNikolov - STScI 2024
    #use_subset_of_slitlets = True
    #nsl_start = 18
    #nsl = 3

    msa_quad_len = 171

    hdul = fits.open(path_to_msa_file + msa_file_id)
    #hdul.info()
    #print(hdul['SHUTTER_INFO'].data[3][3])
    #print(hdul['SHUTTER_INFO'].data)


    c1 = hdul['SHUTTER_INFO'].data['source_id']
    n_slitlets = len(c1)
    sl = np.linspace(1, n_slitlets, n_slitlets)
    # convert from floats to integers
    sl = sl.astype(int)
    #print(sl)

    # Modify the shutter info
    hdul['SHUTTER_INFO'].data['source_id'] = sl
    hdul['SHUTTER_INFO'].data['slitlet_id'] = sl
    hdul['SHUTTER_INFO'].data['background'] = 'N'
    hdul['SHUTTER_INFO'].data['estimated_source_in_shutter_x'] = 0.5
    hdul['SHUTTER_INFO'].data['estimated_source_in_shutter_y'] = 0.5
    hdul['SHUTTER_INFO'].data['primary_source'] = 'Y'
    #print(hdul['SHUTTER_INFO'].data['source_id'])

    # Modify the source info
    tabcol21 = fits.Column(name='program', format='I', array = [prog_id]*n_slitlets)
    tabcol22 = fits.Column(name='source_id', format='I', array = sl)
    #in_list = [prog_id+'_'+str(i+1) for i in range(n_slitlets)]
    #tabcol23 = fits.Column(name='source_name', format='10A', array = [prog_id + '_1']*n_slitlets)
    tabcol23 = fits.Column(name='source_name', format='10A', array = [prog_id+'_'+str(i+1) for i in range(n_slitlets)])
    tabcol24 = fits.Column(name='alias', format='5A', array = ['SFlat']*n_slitlets)
    tabcol25 = fits.Column(name='ra', format='E', array=sl*0)
    tabcol26 = fits.Column(name='dec', format='E', array=sl*0)
    tabcol27 = fits.Column(name='preimage_id', format='8A', array=['NO_PREIM']*n_slitlets)
    tabcol28 = fits.Column(name='stellarity', format='E', array=sl*0) # 0 -- fully extended; 1 -- a perfect point source

    hdul['SOURCE_INFO'] = fits.BinTableHDU.from_columns([tabcol21, tabcol22, tabcol23, tabcol24,  tabcol25, tabcol26, tabcol27, tabcol28], name='SOURCE_INFO')

    if use_subset_of_slitlets:

        if (disperser != 'PRISM'):
            hdul['SHUTTER_INFO'].data = hdul['SHUTTER_INFO'].data[nsl_start:-1*nsl_start]
            hdul['SOURCE_INFO'].data = hdul['SOURCE_INFO'].data[nsl_start:-1*nsl_start]

        if (disperser == 'PRISM'):
            Q1_sh = hdul['SHUTTER_INFO'].data[0:msa_quad_len]
            Q2_sh = hdul['SHUTTER_INFO'].data[msa_quad_len:2*msa_quad_len]
            Q3_sh = hdul['SHUTTER_INFO'].data[2*msa_quad_len:3*msa_quad_len]
            Q4_sh = hdul['SHUTTER_INFO'].data[3*msa_quad_len:]
            # print('Q1_sh: ', Q1_sh)

            Q1_sh = Q1_sh[nsl_start:]
            Q2_sh = Q2_sh[:-1*nsl_start]
            Q3_sh = Q3_sh[nsl_start:]
            Q4_sh = Q4_sh[:-1*nsl_start]

            Q_sh = np.concatenate([Q1_sh, Q2_sh, Q3_sh, Q4_sh], axis=0)
            #print(Q_sh)
            hdul['SHUTTER_INFO'] = fits.BinTableHDU(data = Q_sh, name='SHUTTER_INFO')

            Q1_sr = hdul['SOURCE_INFO'].data[0:msa_quad_len]
            Q2_sr = hdul['SOURCE_INFO'].data[msa_quad_len:2*msa_quad_len]
            Q3_sr = hdul['SOURCE_INFO'].data[2*msa_quad_len:3*msa_quad_len]
            Q4_sr = hdul['SOURCE_INFO'].data[3*msa_quad_len:]
            # print('Q1_sh: ', Q1_sh)

            Q1_sr = Q1_sr[nsl_start:]
            Q2_sr = Q2_sr[:-1*nsl_start]
            Q3_sr = Q3_sr[nsl_start:]
            Q4_sr = Q4_sr[:-1*nsl_start]

            Q_sr = np.concatenate([Q1_sr, Q2_sr, Q3_sr, Q4_sr], axis=0)
            #print(Q_sr)

            hdul['SOURCE_INFO'] = fits.BinTableHDU(data = Q_sr, name='SOURCE_INFO')

    #print('SHUTTER_INFO:')
    #print(hdul['SHUTTER_INFO'].data)
    #print('SOURCE_INFO:')
    #print(hdul['SOURCE_INFO'].data)

    hdul.writeto(path_to_msa_file + msa_file_id, overwrite=True)
    hdul.close()
    return




def msa_modify_v1(path_to_msa_file, msa_file_id, prog_id, use_subset_of_slitlets, nsl_start, nsl):
    # Modifies existing MSA file
    # written by NNikolov - STScI 2024
    #use_subset_of_slitlets = True
    #nsl_start = 18
    #nsl = 3

    hdul = fits.open(path_to_msa_file + msa_file_id)
    #hdul.info()
    #print(hdul['SHUTTER_INFO'].data[3][3])
    #print(hdul['SHUTTER_INFO'].data)


    c1 = hdul['SHUTTER_INFO'].data['source_id']
    n_slitlets = len(c1)
    sl = np.linspace(1, n_slitlets, n_slitlets)
    # convert from floats to integers
    sl = sl.astype(int)
    #print(sl)

    # Modify the shutter info
    hdul['SHUTTER_INFO'].data['source_id'] = sl
    hdul['SHUTTER_INFO'].data['slitlet_id'] = sl
    hdul['SHUTTER_INFO'].data['background'] = 'N'
    hdul['SHUTTER_INFO'].data['estimated_source_in_shutter_x'] = 0.5
    hdul['SHUTTER_INFO'].data['estimated_source_in_shutter_y'] = 0.5
    hdul['SHUTTER_INFO'].data['primary_source'] = 'Y'
    #print(hdul['SHUTTER_INFO'].data['source_id'])

    # Modify the source info
    tabcol21 = fits.Column(name='program', format='I', array = [prog_id]*n_slitlets)
    tabcol22 = fits.Column(name='source_id', format='I', array = sl)
    #in_list = [prog_id+'_'+str(i+1) for i in range(n_slitlets)]
    #tabcol23 = fits.Column(name='source_name', format='10A', array = [prog_id + '_1']*n_slitlets)
    tabcol23 = fits.Column(name='source_name', format='10A', array = [prog_id+'_'+str(i+1) for i in range(n_slitlets)])
    tabcol24 = fits.Column(name='alias', format='5A', array = ['SFlat']*n_slitlets)
    tabcol25 = fits.Column(name='ra', format='E', array=sl*0)
    tabcol26 = fits.Column(name='dec', format='E', array=sl*0)
    tabcol27 = fits.Column(name='preimage_id', format='8A', array=['NO_PREIM']*n_slitlets)
    tabcol28 = fits.Column(name='stellarity', format='E', array=sl*0) # 0 -- fully extended; 1 -- a perfect point source

    hdul['SOURCE_INFO'] = fits.BinTableHDU.from_columns([tabcol21, tabcol22, tabcol23, tabcol24,  tabcol25, tabcol26, tabcol27, tabcol28], name='SOURCE_INFO')

    if use_subset_of_slitlets:
        #hdul['SHUTTER_INFO'].data = hdul['SHUTTER_INFO'].data[nsl_start:nsl_start+nsl]
        #hdul['SOURCE_INFO'].data = hdul['SOURCE_INFO'].data[nsl_start:nsl_start+nsl]
        hdul['SHUTTER_INFO'].data = hdul['SHUTTER_INFO'].data[nsl_start:-1*nsl_start]
        hdul['SOURCE_INFO'].data = hdul['SOURCE_INFO'].data[nsl_start:-1*nsl_start]
        # Modify the data table to have only 3 rows
        #hdul['SHUTTER_INFO'].data = hdul['SHUTTER_INFO'].data[0:3]  # first lines
        #hdul['SHUTTER_INFO'].data = hdul['SHUTTER_INFO'].data[-3:]  # last few lines
        #hdul['SHUTTER_INFO'].data = hdul['SHUTTER_INFO'].data[18:21]  # agter line 17

    #print('SHUTTER_INFO:')
    #print(hdul['SHUTTER_INFO'].data)
    #print('SOURCE_INFO:')
    #print(hdul['SOURCE_INFO'].data)

    hdul.writeto(path_to_msa_file + msa_file_id, overwrite=True)
    hdul.close()
    return






def populate_msa_meta_full(path_to_msa_file, msa_file_id, print_individual_columns, print_msa_shutter_info_before_after, msa_configuration_id):
    #----------------------------------
    # Populate MSA metadata
    # Using the shutter image derive shutter info for long slit S-flat lamps
    # and populate the source info table too
    # First version: 2024
    # N.Nikolov 2024
    # STScI - Baltimore, USA
    #----------------------------------
    # Conda environment: jwst.1.12.5
    #----------------------------------

    #print_individual_columns = 'no'
    #print_msa_shutter_info_before_after = 'yes'

    # Open the MSA file as a data model
    msa = datamodels.open(path_to_msa_file + msa_file_id)

    # Print the full contents
    if (print_individual_columns == 'yes'):
        msa.info(max_rows=None)

    # Obtain the shutter image extension
    msa_im = msa.extra_fits.SHUTTER_IMAGE.data

    # Reorient the shutter image to match the data
    msa_im = shutter_correct(msa_im)

    # Obtain the table with the shutter info and print
    t = Table(msa.extra_fits.SHUTTER_INFO.data)

    if (print_msa_shutter_info_before_after == 'yes'):
        print('')
        print('Input MSA:')
        print("info:")
        t.info()
        print(t)


    # Compute the horizontal and vertical sizes of the MSA image
    msa_xlen = len(msa_im[0,:])
    msa_ylen = len(msa_im[:,0])

    if (print_individual_columns == 'yes'):
        print(msa_xlen, msa_ylen)


    # Obtain the indeces of the slits (these are different from 0)
    idx_slitlets = np.where(msa_im != 0)

    # Compute the number of sluts
    n_slitlets = len(idx_slitlets[0])



    # Prepare the data to populate the individual columns of the table

    # 01: Assign slitled ID of an integer
    msa_slitlet_id = 1+np.arange(n_slitlets)

    if (print_individual_columns == 'yes'):
        print('msa_slitlet_id', msa_slitlet_id)


    # 02: Obtain the metafile ID and set it in an array of the same size as the array for slitlets
    msa_metad_id = msa_slitlet_id * 0 + msa_configuration_id # datamodels.open(nrs_file_with_path).meta.instrument.msa_configuration_id

    if (print_individual_columns == 'yes'):
        #print(datamodels.open(path_to_msa_file+nrs1_file+'_rate.fits').meta.instrument.msa_metadata_id)
        print('msa_metad_id', msa_metad_id)


    # 03: Compute shutter, row, column and quadrant
    msa_shutter_quadrant, msa_shutter_row, msa_shutter_col, msa_shutter_row_quad_coord = find_shutter_coords(msa_im)

    if (print_individual_columns == 'yes'):
        print('quadrant IDs:', msa_shutter_quadrant)
        print('shuttrow IDs:', msa_shutter_row)
        print('shuttcol IDs:', msa_shutter_col)
        print('shuttrow Q IDs:', msa_shutter_row_quad_coord)


    # source_id: need to be an integer number
    msa_source_id = msa_slitlet_id #0*msa_shutter_col + 1111

    if (print_individual_columns == 'yes'):
        print('msa_source_id:', msa_source_id)


    # background    str1, N for long slits
    backg = ["N"]*len(msa_shutter_row)

    if (print_individual_columns == 'yes'):
        print('Background: ', backg)
        #print(backg[0])


    # shutter_state    str6, OPEN for all long slits
    shutter_state = ['OPEN']*len(msa_shutter_row)

    if (print_individual_columns == 'yes'):
        print('Shutter state: ', shutter_state)


    # estimated_source_in_shutter_x float32, 0.5 for the center of the slitlet
    # estimated_source_in_shutter_y float32, 0.5 for the center of the slitlet
    estimated_source_in_shutter_x = 0.5 + 0.*msa_shutter_quadrant
    estimated_source_in_shutter_y = 0.5 + 0.*msa_shutter_quadrant

    if (print_individual_columns == 'yes'):
        print('estimated_source_in_shutter_x', estimated_source_in_shutter_x)
        print('estimated_source_in_shutter_y', estimated_source_in_shutter_y)


    # dither_point_index   int16, 1 for no dithers
    dither_point_index = 1 + 0*msa_shutter_quadrant

    if (print_individual_columns == 'yes'):
        print('dither_point_index', dither_point_index)


    # primary_source    str1, all long slits are primary sources
    primary_source = ["Y"]*len(msa_shutter_row)

    if (print_individual_columns == 'yes'):
        print('primary_source: ', primary_source)




    # Open the MSA file as FITS, modify and save the contents

    # Open the existing MSA file and update its contents
    hdul = fits.open(path_to_msa_file + msa_file_id)

    #idx_use = 'good' # This works great!
    #idx_use = 'top'
    idx_use = 'all'

    #idx_use = 'bottom'

    if (idx_use == 'good'):
        n1=53 #18
        n2=56 #20

    if (idx_use == 'top'):
        #n1=15
        #n2=25
        n1 = 161
        n2 = 171

    if (idx_use == 'all'):
        #n1=15
        #n2=25
        n1 = 15
        n2 = 171


    # construct a new table for the SHUTTER_INFO
    tabcol1 = fits.Column(name='slitlet_id', format='I', array=msa_slitlet_id[n1:n2])
    tabcol2 = fits.Column(name='msa_metadata_id', format='I', array=msa_metad_id[n1:n2])
    tabcol3 = fits.Column(name='shutter_quadrant', format='I', array=msa_shutter_quadrant[n1:n2])
    tabcol4 = fits.Column(name='shutter_row', format='I', array=msa_shutter_col[n1:n2])
    tabcol5 = fits.Column(name='shutter_column', format='I', array=msa_shutter_row_quad_coord[n1:n2])
    tabcol6 = fits.Column(name='source_id', format='I', array=msa_source_id[n1:n2])
    tabcol7 = fits.Column(name='background', format='A', array=backg[n1:n2])
    tabcol8 = fits.Column(name='shutter_state', format='4A', array=shutter_state[n1:n2])
    tabcol9 = fits.Column(name='estimated_source_in_shutter_x', format='E', array=estimated_source_in_shutter_x[n1:n2])
    tabcol10 = fits.Column(name='estimated_source_in_shutter_y', format='E', array=estimated_source_in_shutter_y[n1:n2])
    tabcol11 = fits.Column(name='dither_point_index', format='I', array=dither_point_index[n1:n2])
    tabcol12 = fits.Column(name='primary_source', format='1A', array=primary_source[n1:n2])
    hdul['SHUTTER_INFO'] = fits.BinTableHDU.from_columns([tabcol1, tabcol2, tabcol3, tabcol4, tabcol5,
                                                          tabcol6, tabcol7, tabcol8, tabcol9, tabcol10,
                                                          tabcol11, tabcol12], name='SHUTTER_INFO')




    # construct a new table for the SOURCE_INFO
    # PROGRAM, SOURCE_ID, SOURCE_NAME, ALIAS, RA, DEC, PREIMAGE_ID, STELLARITY
    ms = datamodels.open(path_to_msa_file + msa_file_id)
    t2 = Table(msa.extra_fits.SOURCE_INFO.data)
    #print("SOURCE_INFO:")
    #t2.info()

    n_sources = n_slitlets #2
    arr_pid = 0*np.arange(n_sources)+int(msa.meta.observation.visit_id[1:5])
    #print(int(msa.meta.observation.visit_id[1:5]))
    #print(arr_pid)
    tabcol21 = fits.Column(name='program', format='I', array=arr_pid[n1:n2])

    arr_source_id = msa_slitlet_id #0*np.arange(n_sources)+1111
    tabcol22 = fits.Column(name='source_id', format='I', array=arr_source_id[n1:n2])

    arr_source_name = [str(msa.meta.observation.visit_id[1:5])+'_'+str(int(msa_source_id[0]))]*n_sources
    # ["N"]*len(msa_shutter_row)
    tabcol23 = fits.Column(name='source_name', format='10A', array=arr_source_name[n1:n2])

    arr_alias = ['SFlat']*n_sources
    tabcol24 = fits.Column(name='alias', format='5A', array=arr_alias[n1:n2])
    #print(arr_alias[n1:n2])

    arr_ra = 0.+0.*np.arange(n_sources)
    tabcol25 = fits.Column(name='ra', format='E', array=arr_ra[n1:n2])

    arr_dec = 0.+0.*np.arange(n_sources)
    tabcol26 = fits.Column(name='dec', format='E', array=arr_dec[n1:n2])

    arr_preim_id = ['NO_PREIM']*n_sources
    tabcol27 = fits.Column(name='preimage_id', format='8A', array=arr_preim_id[n1:n2])

    arr_stellarity = 0.+0.*np.arange(n_sources)
    tabcol28 = fits.Column(name='stellarity', format='E', array=arr_stellarity[n1:n2]) # 0 -- fully extended; 1 -- a perfect point source
    #hdul['SOURCE_INFO'] = fits.BinTableHDU.from_columns([tabcol21, tabcol22, tabcol23], name='SOURCE_INFO')
    #hdul['SOURCE_INFO'] = fits.BinTableHDU.from_columns([tabcol21, tabcol22, tabcol23], name='SOURCE_INFO')

    hdul['SOURCE_INFO'] = fits.BinTableHDU.from_columns([tabcol21, tabcol22, tabcol23, tabcol24, tabcol25, tabcol26, tabcol27, tabcol28], name='SOURCE_INFO')

    #    SOURCE_INFO:
    #<Table length=0>
    #    name     dtype
    #----------- -------
    #    program   int32
    #  source_id   int32
    #source_name   str20
    #      alias   str31
    #         ra float64
    #        dec float64
    #preimage_id   str30
    # stellarity float64



    # save to the same file (this will overwrite the existing file)
    print('Saving MSA modifications to file ...')
    hdul.writeto(path_to_msa_file + msa_file_id, overwrite=True)
    hdul.close()


    # Check the modified file
    if (print_msa_shutter_info_before_after == 'yes'):
        print('')
        print('Checking new MSA ...')
        print('Output MSA:')
        msa = datamodels.open(path_to_msa_file + msa_file_id)
        t = Table(msa.extra_fits.SHUTTER_INFO.data)
        print("SHUTTER_INFO:")
        print(t)
        t.info()

        print(" ")
        t = Table(msa.extra_fits.SOURCE_INFO.data)
        print("SOURCE_INFO:")
        print(t)
        t.info()



    if (print_individual_columns == 'yes'):
        print('')
        print('MSA datamodel contents (FULL): ')
        msa.info(max_rows=None)





def msa_ls_detstat(shutter_image):
    # Finds which detector data has the long slit
    # The PRISM SFlats has two slits
    # The high resolution gratings have a single slit

    row_data_low  = shutter_image[85,:]
    row_data_high = shutter_image[250,:]

    #print(row_data_low)

    #idx_low  = np.where(row_data_low < 0.)[0]
    #idx_high = np.where(row_data_high < 0.)[0]

    idx_low  = np.where(row_data_low < 0.)[0]
    idx_high = np.where(row_data_high < 0.)[0]
    #if ((len(idx_low) == 1) & ((len(idx_high) == 1)):

    if ((len(idx_low) == 1)):
        if ((idx_low < 365) & (idx_high < 365)):
            detector_status = 'nrs1'
            #print('Detector status: ', detector_status)
        if ((idx_low >= 365) & (idx_high >= 365)):
            detector_status = 'nrs2'
            #print('Detector status: ', detector_status)

    if (len(idx_low) > 1):
        detector_status = 'nrs1nrs2'

    #print(idx_low)
    #print(idx_high)

    return detector_status





def nn_print_header(fits_with_path):
    # Prints the full header of a FITS file given a path with a file
    with fits.open(fits_with_path) as hdulist:  # this is like the "hdulist = fits.open('test1.fits')"
        print('')
        print('')
        hdulist.info()
        print('')
        print('')
        for hdu in hdulist:
            print(repr(hdu.header)) # repr returns a printable representation
    # After leaving the "with" the file closes.




def nn_find_header_keyword(fits_with_path, key_word):
    # Prints the value of a FITS header, given a path with the file and a keyword
    with fits.open(fits_with_path) as hdulist:  # this is like the "hdulist = fits.open('test1.fits')"
        hdr = hdulist[0].header
        #print('Header ('+key_word+'): ', hdr[key_word])
        flag = 0
        for cnts in hdr:
            if (cnts == key_word):
                flag = 1
            else:
                flag += 0
        if (flag == 0):
            print(key_word, ': Not found in header!')
        else:
            print(key_word+': ', hdr[key_word])
