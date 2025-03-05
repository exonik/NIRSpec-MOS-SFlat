# NIRSpec S-Flat calweb

# Author: Nikolay Nikolov, AURA Associate Scientist, NIRSpec branch
# Pipeline Version**: 1.15.1
# Last Update**: November 12, 2024

# Conda environment (laptop): jwst.1.15.1
# Conda environment (dlsci122): jwst_1.12.5


import sflat_red_utils

# Python math functions
import numpy as np

# Time Utilities
import time as tt
# from time import gmtime, strftime
import datetime
from datetime import datetime
from datetime import date

# Directory and file services e.g., unpacking
import os
# import shutil
import subprocess

# Plot functions
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

# Astropy and other spectroscopy utilities
# import astropy
from astropy.io import fits
# from astropy.visualization import astropy_mpl_style
# from astropy.table import Table
# from astropy.time import Time
# from astropy.io import fits
# import astropy.io.fits as fits
# from astropy.timeseries import LombScargle
# from astropy.utils.data import download_file

# Astroquerry
# import astroquery
from astroquery.mast import Observations


# Multiprocessing functions
import multiprocessing as mp

# JWST Pipeline
import jwst
from jwst import datamodels
# from jwst.pipeline import calwebb_spec2
# from jwst.assign_wcs import nirspec
# from jwst.assign_wcs import AssignWcsStep

from jwst.pipeline import Spec2Pipeline


# Pandas
import pandas as pd

import crds
import json

from jwst.assign_wcs.util import NoDataOnDetectorError



prog_id = '1122'
disperser = 'PRISM'
# obs = [1]
# -1: ALL; 1,2,3: observation numbers

# apt_list = [1122]#, 1485]

work_dir = './'
apt_dir = work_dir #+ 'APT/'
dat_dir = work_dir + 'DATA/'
stage2_dir = work_dir + 'stage_2/'
plot_leakcal_dir = work_dir + 'plot_leakcal/'
plot_rate_dir = work_dir + 'plot_rate/'

# -------------------------------
do_plot_and_correct_data = 'yes'
#use_subset_list = 'no'
# -------------------------------
do_run_calwebb = 'yes'
# -------------------------------


# ------------------------------------------
# Time the execution time -> START
time_start = tt.time()
now = datetime.now()
start_time_formatted = now.strftime("%Y-%m-%d %H:%M:%S")
print("Run time (sta): ", start_time_formatted)
# ------------------------------------------


print(' ')
print("JWST Calibration Pipeline Version={}".format(jwst.__version__))
print("Current Operational CRDS Context = {}".format(crds.get_default_context()))
print("Number of processors: ", mp.cpu_count())
print("CRDS directory: ", os.environ['CRDS_PATH'])
print(' ')

# ------------------------------------------------
# STEP 1.0: Identify and remove single exposure
# observations and observations without leakcals
# Subtract the leakcal for each observations
# ------------------------------------------------
print('Identifying leakcals and reducing the data list ...')

# ---
# Filter the observations that have just one exposure
# It will be either an LS without leakcal or a leakcal without an LS
# In both cases wont be used
# ---

# Read the data frame
filtered_df = pd.read_csv(apt_dir+prog_id+'_'+disperser+'.csv')

# Count the occurrences of each value in 'ObsNumber'
ObsNum_counts = filtered_df['ObsNumber'].value_counts()

# Filter the rows where the value occurs more than once
red_df = filtered_df[filtered_df['ObsNumber'].map(ObsNum_counts) > 1]

# Update the data frame
filtered_df = red_df
print('')


# ---
# Filter observations that dont have leakcals
# ---

# Function to check if 'LS-' exposure exists
def has_LS(x):
    return any(val.startswith('LS-') for val in x['MsaConfigFile'])

# Group by observation number and filter groups that have exposures called 'ALLCLOSED' and 'LS-*' variation
red_df = filtered_df.groupby('ObsNumber').filter(lambda x: 'ALLCLOSED' in x['MsaConfigFile'].values and has_LS(x))

filtered_df = red_df

#print(filtered_df)
print(filtered_df[['ProposalID', 'ObsNumber', 'OperatingMode', 'MsaConfigFile', 'Groups', 'Integrations', 'Lamp', 'Grating', 'ActivityNum', 'DataFileNameNRS1', 'DataFileNameNRS2']])

print('')


# ---
# Check if the data directory exists and create it if not present
# ---
isExist_DATA = os.path.exists(dat_dir)

if isExist_DATA == False:
    # create a new empty DATA directory to store the data files
    print("Creating a new DATA directory ...")
    # create the results directory
    os.mkdir(dat_dir)

# ---
# Check if the stage_2 directory exists and create it if not present
# ---
isExist_stage2 = os.path.exists(stage2_dir)

if isExist_stage2 == False:
    # create a new empty stage_2 directory to store the data files
    print("Creating a new stage_2 directory ...")
    # create the results directory
    os.mkdir(stage2_dir)

# ---
# Check if the plot_leakcal directory exists and create it if not present
# ---
isplot_leakcal = os.path.exists(plot_leakcal_dir)

if isplot_leakcal == False:
    # create a new empty plot_leakcal directory to store the plot files from leakcal correction
    print("Creating a new plot_leakcal directory ...")
    os.mkdir(plot_leakcal_dir)

# ---
# Check if the plot_rate directory exists and create it if not present
# ---
isplot_rate = os.path.exists(plot_rate_dir)

if isplot_rate == False:
    # create a new empty plot_rate directory to store the plot files from leakcal correction
    print("Creating a new plot_rate directory ...")
    os.mkdir(plot_rate_dir)


# ---
# Loop the data frame and group by observation; download the ALLCLOSED and LS- frames per exposure and subtract the LEAKCAL
# ---

# Define the MAST dir to intiate the download
mast_dir = 'mast:jwst/product'

for obs, group in filtered_df.groupby('ObsNumber'):
    leakcal = group[group['MsaConfigFile'] == 'ALLCLOSED']
    LS = group[group['MsaConfigFile'].str.startswith('LS-')]

    # LeakCals
    print(f"Observation {obs}:")
    print('LeakCals root names: ', leakcal['DataFileNameNRS1'].values, leakcal['DataFileNameNRS2'].values)

    # NRS1
    mast_path = os.path.join(mast_dir, leakcal['DataFileNameNRS1'].values[0] + '_rate.fits')
    Observations.download_file(mast_path, cache=False)

    # NRS2
    mast_path = os.path.join(mast_dir, leakcal['DataFileNameNRS2'].values[0] + '_rate.fits')
    Observations.download_file(mast_path, cache=False)

    # Read the nrs1/2 LeakCal files
    with fits.open(work_dir + leakcal['DataFileNameNRS1'].values[0] + '_rate.fits') as hdul_leakcal_nrs1:
        print(hdul_leakcal_nrs1.info())
        im_leakcal_nrs1 = hdul_leakcal_nrs1[1].data
        err_leakcal_nrs1 = hdul_leakcal_nrs1[2].data

    with fits.open(work_dir + leakcal['DataFileNameNRS2'].values[0] + '_rate.fits') as hdul_leakcal_nrs2:
        print(hdul_leakcal_nrs2.info())
        im_leakcal_nrs2 = hdul_leakcal_nrs2[1].data
        err_leakcal_nrs2 = hdul_leakcal_nrs2[2].data

    # LS
    for ls_nrs in LS['DataFileNameNRS1'].values:
        print('LS root names: ', ls_nrs[:-4]+'nrs1', ls_nrs[:-4]+'nrs2')

        # Download the LS NRS1 files
        mast_path = os.path.join(mast_dir, ls_nrs[:-4]+'nrs1_rate.fits')
        Observations.download_file(mast_path, cache=False)

        # Download the LS NRS2 files
        mast_path = os.path.join(mast_dir, ls_nrs[:-4]+'nrs2_rate.fits')
        Observations.download_file(mast_path, cache=False)

        # Copy temporarily at tests to check difference
        #subprocess.run(["cp " + work_dir + ls_nrs[:-4] + 'nrs1_rate.fits ' + dat_dir], shell=True)
        #subprocess.run(["cp " + work_dir + ls_nrs[:-4] + 'nrs2_rate.fits '  + dat_dir], shell=True)

        # Read the nrs1/2 LS files
        with fits.open(work_dir + ls_nrs[:-4]+'nrs1_rate.fits', mode="update") as hdul_LS_nrs1:
            print(hdul_LS_nrs1.info())
            im_LS_nrs1 = hdul_LS_nrs1[1].data.copy()

            # Compute the difference in-place
            hdul_LS_nrs1[1].data -= im_leakcal_nrs1.copy()
            # Propagate the error of difference into the error array of the corrected for leakcal FITS
            hdul_LS_nrs1[2].data = np.sqrt(hdul_LS_nrs1[2].data**2 + err_leakcal_nrs1**2)
            hdul_LS_nrs1.flush()



        with fits.open(work_dir + ls_nrs[:-4]+'nrs2_rate.fits', mode="update") as hdul_LS_nrs2:
            print(hdul_LS_nrs2.info())
            im_LS_nrs2 = hdul_LS_nrs2[1].data.copy()

            # Compute the difference in-place
            hdul_LS_nrs2[1].data -= im_leakcal_nrs2.copy()
            hdul_LS_nrs2[2].data = np.sqrt(hdul_LS_nrs2[2].data**2 + err_leakcal_nrs2**2)
            hdul_LS_nrs2.flush()

        im_LS_diff_nrs1 = im_LS_nrs1 - im_leakcal_nrs1
        im_LS_diff_nrs2 = im_LS_nrs2 - im_leakcal_nrs2



        # Check what MSA file is used and download it as well
        msa_12 = datamodels.open(work_dir + ls_nrs[:-4]+'nrs1_rate.fits').meta.instrument.msa_metadata_file
        # print(f"{filtered_df['DataFileNameNRS1'][j]}_rate.fits uses:", msa_12)

        # Download the MSA file
        mast_path = os.path.join(mast_dir, msa_12)
        Observations.download_file(mast_path, cache=False)



        # ------------------------------------------------
        # Plot Original, LeakCal and Subtracted NRS 1/2
        # ------------------------------------------------
        fig_sz = (20, 10)

        f_lo_1 = -10.  # median_image - delta_im_flux
        f_hi_1 = +10. #median_image + delta_im_flux

        f_lo_2 = -10. #np.nanmin(msa_im)-0.5
        f_hi_2 = +10. # np.nanmax(msa_im)+0.5

        # Plot the NRS1 and NRS2 images along with the MSA shutter image and save to pdf
        sflat_red_utils.plot_mos_msa(
        im_LS_nrs1,
        im_leakcal_nrs1,
        im_LS_diff_nrs1,
        'Original rate: ' + ls_nrs[:-4]+'nrs1_rate.fits',
        'Leackcal rate: ' + leakcal['DataFileNameNRS1'].values[0] + '_rate.fits',
        'Difference: Original - Leackcal',
        f_lo_1, f_hi_1,
        f_lo_2, f_hi_2,
        fig_sz,
        work_dir + ls_nrs[:-4]+'nrs1_plot_leakcal.pdf',
        'Disperser: ' + disperser)

        sflat_red_utils.plot_mos_msa(
        im_LS_nrs2,
        im_leakcal_nrs2,
        im_LS_diff_nrs2,
        'Original rate: ' + ls_nrs[:-4]+'nrs2_rate.fits',
        'Leackcal rate: ' + leakcal['DataFileNameNRS2'].values[0] + '_rate.fits',
        'Difference: Original - Leackcal',
        f_lo_1, f_hi_1,
        f_lo_2, f_hi_2,
        fig_sz,
        work_dir + ls_nrs[:-4]+'nrs2_plot_leakcal.pdf',
        'Disperser: ' + disperser)

        # Move the data in the DATA directory
        subprocess.run(["mv " + work_dir + ls_nrs[:-4] + 'nrs1_rate.fits ' + dat_dir], shell=True)
        subprocess.run(["mv " + work_dir + ls_nrs[:-4] + 'nrs2_rate.fits '  + dat_dir], shell=True)
        subprocess.run(["mv " + work_dir + msa_12 + ' ' + dat_dir], shell=True)

    subprocess.run(["mv " + leakcal['DataFileNameNRS1'].values[0] + '_rate.fits ' + dat_dir], shell=True)
    subprocess.run(["mv " + leakcal['DataFileNameNRS2'].values[0] + '_rate.fits ' + dat_dir], shell=True)
    print('')


subprocess.run(["cd " + work_dir + " ; mv *plot_leakcal.pdf " + plot_leakcal_dir], shell=True)


# Remove the LEAKCALS from the data frame and proceed with the plotting and pipeline run
new_df = filtered_df[filtered_df['MsaConfigFile'] != 'ALLCLOSED']
print('new_dataframe:')
print(new_df[['ProposalID', 'ObsNumber', 'OperatingMode', 'MsaConfigFile', 'Groups', 'Integrations', 'Lamp', 'Grating', 'ActivityNum', 'DataFileNameNRS1', 'DataFileNameNRS2']])
print('')

# Save the filtered data frame
new_df.to_csv(work_dir + prog_id+'_LS_only.csv', index=False)






# Old code, where filtering is needed
# Load the list of files to be processed for the desired observations
# df = pd.read_csv(apt_dir + 'APT_'+prog_id+'_Filtered.csv')


# if (prog_id == '1122'):
#    if (len(obs) >= 1):
#        filtered_df = df[df['ObsNumber'].isin(obs)]
#
#    if ((len(obs) == 1) & (obs[0] == -1)):
#        filtered_df = df



# New Code with pre-filtered dfs
#filtered_df = pd.read_csv(apt_dir+prog_id+'_'+disperser+'.csv')
filtered_df = pd.read_csv(apt_dir+prog_id+'_LS_only.csv')

nrs1_file = filtered_df['DataFileNameNRS1'].values
nrs2_file = filtered_df['DataFileNameNRS2'].values

obs_num_id = filtered_df['ObsNumber'].values

msa_config_file = filtered_df['MsaConfigFile'].values
lamp_id = filtered_df['Lamp'].values
grating_id = filtered_df['Grating'].values
rop_id = filtered_df['ReadoutPattern'].values
groups_id = filtered_df['Groups'].values
activity_id = filtered_df['ActivityNum'].values


# Use only a subset of files for testing
#if (use_subset_list == 'yes'):
#    nrs1_file = nrs1_file[0:2]
#    nrs2_file = nrs2_file[0:2]


print(' ')
print('Number of files to be processed: ', len(nrs1_file))
print('Files names to be processed: ', nrs1_file)
print(' ')







# ------------------------------------------------
# STEP 1: Download and plot rates and shutter
# ------------------------------------------------
if (do_plot_and_correct_data == 'yes'):
    # Produce the metadata and plot the data along with the shutter info image
    # This step needs to be done just once after each data
    # download to make sure the metadata is not corrected twice

    # Define the MAST dir to intiate the download
    #mast_dir = 'mast:jwst/product'

    for i in range(len(nrs1_file)):

        # ---------- A -----------
        # Download and move data
        # ------------------------

        #mast_path = os.path.join(mast_dir, nrs1_file[i] + '_rate.fits')
        #Observations.download_file(mast_path, cache=False)

        # Download the NRS2 file
        #mast_path = os.path.join(mast_dir, nrs2_file[i] + '_rate.fits')
        #Observations.download_file(mast_path, cache=False)

        # Check what MSA file is used and download it as well
        #msa_12 = datamodels.open(work_dir + nrs1_file[i] + '_rate.fits').meta.instrument.msa_metadata_file
        # print(f"{filtered_df['DataFileNameNRS1'][j]}_rate.fits uses:", msa_12)

        # Download the MSA file
        #mast_path = os.path.join(mast_dir, msa_12)
        #Observations.download_file(mast_path, cache=False)

        # Move the data in the DATA directory
        #subprocess.run(["mv " + work_dir + nrs1_file[i] + '_rate.fits ' + dat_dir], shell=True)
        #subprocess.run(["mv " + work_dir + nrs2_file[i] + '_rate.fits ' + dat_dir], shell=True)
        #subprocess.run(["mv " + work_dir + msa_12 + ' ' + dat_dir], shell=True)

        # ---------- B -----------
        # Plot NRS1,2 and SHUTTER
        # ------------------------

        # print(nrs1_file[i])

        # im1_dm = datamodels.open(dat_dir + nrs1_file[i] + '_rate.fits')
        # im1_dm.info(max_rows=None)
        # grating_ id = im1_dm.meta.instrument.grarting
        # instrument_ id = im1_dm.meta.instrument
        # lamp_state = im1_dm.meta.instrument.lamp_state

        image1 = datamodels.open(dat_dir + nrs1_file[i] + '_rate.fits').data
        image2 = datamodels.open(dat_dir + nrs2_file[i] + '_rate.fits').data
        msa_id = datamodels.open(dat_dir + nrs1_file[i] + '_rate.fits').meta.instrument.msa_metadata_file

        # print(msa_id)

        # Obtain the MSA shutter image
        msa = datamodels.open(dat_dir + msa_id)
        # msa.info(max_rows=None)
        msa_im = msa.extra_fits.SHUTTER_IMAGE.data

        # ---------------------------------------
        # Correct the shutter image orientation
        # for plotting purposes and populating
        # the metadata
        # ---------------------------------------
        msa_im = sflat_red_utils.shutter_correct(msa_im)
        # this is used for plotting here
        # ---------------------------------------

        fig_sz = (20, 10)

        # Obtain the PID
        # print(nrs1_file[3:7])
        pid_from_name = nrs1_file[3:7]
        # print(pid_from_name)

        # Obtain the lamp_state value
        flat_state = datamodels.open(dat_dir + nrs1_file[i] + '_rate.fits').meta.instrument.lamp_state
        # if (flat_state[0:4] == 'FLAT'):

        # Find out which detector images the long slit spectra
        det_stat = sflat_red_utils.msa_ls_detstat(msa_im)
        if (det_stat == 'nrs1'):
            median_image = np.nanmedian(image1[1250, :])
        if (det_stat == 'nrs2'):
            median_image = np.nanmedian(image2[1250, :])
        if (det_stat == 'nrs1nrs2'):
            # print('len image1[0,:]: ', len(image1[0,:]))
            # print('len image1[:,0]: ', len(image1[:,0]))
            # print('Image 1:', image1[1250,:])
            # print('Image 2:', image2[1250,:])
            median_image = np.nanmedian(np.concatenate((image1[1250, :], image2[1250, :])))

        delta_im_flux = 500.

        f_lo_1 = -10.  # median_image - delta_im_flux
        f_hi_1 = median_image + delta_im_flux

        f_lo_2 = np.nanmin(msa_im)-0.5
        f_hi_2 = np.nanmax(msa_im)+0.5

        # Plot the NRS1 and NRS2 images along with the MSA shutter image and save to pdf
        sflat_red_utils.plot_mos_msa(image1, image2, msa_im, nrs1_file[i] + '_rate.fits, activ.=' + str(activity_id[i]) , nrs2_file[i] + '_rate.fits ',
        msa_id, f_lo_1, f_hi_1, f_lo_2, f_hi_2, fig_sz,
        stage2_dir + nrs1_file[i][0:-4]+'plot.pdf', grating_id[i] + ', ' + lamp_id[i] + ', ' + msa_config_file[i] + ', ' + rop_id[i] + ', gr=' + str(groups_id[i]))

        # sflat_red_utils.print_meta(datamodels.open(dat_dir + nrs1_file[i] + '_rate.fits'))
        # sflat_red_utils.print_meta(datamodels.open(dat_dir + nrs2_file[i] + '_rate.fits'))


    subprocess.run(["cd " + stage2_dir + " ; mv *plot.pdf ../plot_rate/"], shell=True)


# ------------------------------------------------
# STEP 2: Populate the metadata and run calwebb
# ------------------------------------------------
if (do_run_calwebb == 'yes'):

    for i in range(len(nrs1_file)):

        # ---------- A -----------
        # Populate the metadata
        # ------------------------
        msa_id = datamodels.open(dat_dir + nrs1_file[i] + '_rate.fits').meta.instrument.msa_metadata_file
        #msa_configuration_id = datamodels.open(dat_dir + nrs1_file[i] + '_rate.fits').meta.instrument.msa_configuration_id
        # print("MSA ID: ", msa_id)

        #print_individual_columns = 'no'
        #print_msa_shutter_info_before_after = 'yes'

        subprocess.run(["cp " + dat_dir + nrs1_file[i] + '_rate.fits ' + work_dir], shell=True)
        subprocess.run(["cp " + dat_dir + nrs2_file[i] + '_rate.fits ' + work_dir], shell=True)
        subprocess.run(["cp " + dat_dir + msa_id + ' ' + work_dir], shell=True)
        #sflat_red_utils.populate_msa_meta_full(work_dir, msa_id, print_individual_columns, print_msa_shutter_info_before_after, msa_configuration_id)

        # Modify metafile
        path_to_msa_file = work_dir
        msa_file_id = msa_id
        use_subset_of_slitlets = True
        nsl_start = 17 #17 #169 # 27 # 169 would give only 4 slitlets
        nsl = 3
        sflat_red_utils.msa_modify_v2(path_to_msa_file, msa_file_id, prog_id, use_subset_of_slitlets, nsl_start, nsl, disperser)

        # ---------- B -----------
        # Identify NRS? behind LS
        # ------------------------
        #msa_id = datamodels.open(work_dir + nrs1_file[i] + '_rate.fits').meta.instrument.msa_metadata_file
        #msa = datamodels.open(work_dir + msa_id)
        # msa.info(max_rows=None)
        #msa_im = msa.extra_fits.SHUTTER_IMAGE.data
        #msa_im = sflat_red_utils.shutter_correct(msa_im)
        #det_stat = sflat_red_utils.msa_ls_detstat(msa_im)

        #if (det_stat == 'nrs1'):
        #    rate_file = [work_dir + nrs1_file[i] + '_rate.fits']
        #if (det_stat == 'nrs2'):
        #    rate_file = [work_dir + nrs2_file[i] + '_rate.fits']
        #if (det_stat == 'nrs1nrs2'):
        #    rate_file = [work_dir + nrs1_file[i] + '_rate.fits', work_dir + nrs2_file[i] + '_rate.fits']

        # ---------- C -----------
        # Update extraction reference file
        # ------------------------

        # We modify the parameter reference file of the extract_1d step to specify the width of the extraction aperture.
        # This reference file is located in the CRDS_PATH directory.

        json_name = 'jwst_nirspec_extract1d_0003.json'
        json_dir = os.environ['CRDS_PATH'] + '/references/jwst/nirspec/'

        # open file in read-mode
        with open(json_dir + json_name, 'r') as file:
            # read JSON data
            data = json.load(file)
            # data['apertures'][0]['extract_width'] = 4 # 22 extraction aperture radius in pixels
            data['apertures'][0].pop("extract_width")
            # 22 extraction aperture radius in pixels
            data['apertures'][0]['ystart'] = 4#2.5
            # 22 extraction aperture radius in pixels
            data['apertures'][0]['ystop'] = 5# 5.5

        newData = json.dumps(data, indent=0)
        # print(json_name[:-5])
        # add the suffix '_bots' to distinguish the file from the default version.
        with open(json_dir + json_name[:-5] + '_MOS_SFLAT_NN.json', 'w') as file:
            # write
            file.write(newData)


    rate_file = np.concatenate((nrs1_file, nrs2_file))
    rate_file.sort()
    #print(rate_file)

    for j in range(len(rate_file)):
        rate_file[j] = work_dir + rate_file[j] + '_rate.fits'
        # ---------- D -----------
        # Update rate header
        # ------------------------
        print(rate_file[j])
        data_model = datamodels.open(rate_file[j])

        # Obtain the MSA file name
        # msa_conf_id = data_model.meta.instrument.msa_configuration_id
        # print(' msa_conf_id FILE: ', msa_conf_id)

        # Update the value of the FITS file header using the datamodel information
        #fits.setval(rate_file[j], 'MSAMETID', value=msa_conf_id)

        # Update the value of the FITS file header exp_type
        # fits.setval(rate_file[j], 'EXP_TYPE', value='LAMP_FLAT')

        # Open the FITS file and check if the header has been updated
        hdul = fits.open(rate_file[j])
        hdr = hdul[0].header  # the primary HDU header
        print('Header (after change): MSAMETID', hdr['MSAMETID'])
        print('Header (after change): EXP_TYPE', hdr['EXP_TYPE'])

        # ---------- E -----------
        # Run calwebb
        # ------------------------

        run_stage2_mode = 'all' # all or step

        if (run_stage2_mode == 'all'):
            # Prepare to run the pipeline
            fits.setval(rate_file[j], 'EXP_TYPE', value='NRS_MSASPEC')

            param_stage2 = {
                'assign_wcs': {'skip': False, 'save_results': True},  # Assigns wavelength solution for spectroscopic data
                'imprint_subtract': {'skip': True},  # MOS and IFU
                'bkg_subtract': {'skip': True},
                'msa_flagging': {'skip': False, 'save_results': False},  # MSA and IFU; assign alway OPEN and CLOSED shutters
                'extract_2d': {'skip': False, 'save_results': True},  # Obtains 2D cutouts; MOS and FS and WFSS;
                'srctype': {'skip': False, 'save_results': False, 'source_type': 'EXTENDED'},  # Assigns source type: 'POINT' or 'EXTENDED'   !!!! Check the intermediate product to see why Barshadow is skipped
                'master_background_mos': {'skip': True},  # MOS master bckg
                'wavecorr': {'skip': True},  # FS, MOS
                'flat_field': {'skip': False, 'save_results': False},  # Fore- and Spectroscopic- F-, and S-flats are skipped, but D-flat is applied
                'pathloss': {'skip': True},  # need to be skipped; loss of signal along the optical path
                'barshadow': {'skip': False, 'save_results': False, 'source_type': 'EXTENDED'},  # MSA and MOS only for standard target; not for MOS internal flat
                'photom': {'skip': True},  # We also skip this step, because BOTS exoplanet observations are relative; applies flux calibration to surface brightness
                'pixel_replace': {'skip': False, 'save_results': True, 'n_adjacent_cols': 15, 'algorithm': "fit_profile"},   #, 'algorithm':'mingrad'; resamples the 2D image based on the WCS and distortion information
                'resample_spec': {'skip': False, 'save_results': True},  # ONLY if multiple files exist; this would combine multi-dithered exposures
                'extract_1d': {'skip': False, 'save_results': True, 'override_extract1d': f'{json_dir + json_name[:-5]}_MOS_SFLAT_NN.json', 'apply_apcorr': False, 'subtract_background': False, 'use_source_posn': False}}

            start_pp = tt.time()

            # From Chris; no need to identify the which NRS1/2 file has the data
            try:
                Spec2Pipeline.call(rate_file[j], save_results=True, steps=param_stage2, output_dir=stage2_dir)
            except NoDataOnDetectorError:
                continue

            #Spec2Pipeline.call(rate_file[j], save_results=True,
            #                   steps=param_stage2, output_dir=stage2_dir)
            end_pp = tt.time()
            print("Run time on one image: ", round(end_pp-start_pp,1)/60.0, " min")

# ------------------------------------------
# Step 0: time the execution time -> END
time_end = tt.time()
dtime = (time_end-time_start)/(60.)
end_time_formatted = now.strftime("%Y-%m-%d %H:%M:%S")

print("Run time (sta): ", start_time_formatted)
print("Run time (end): ", end_time_formatted)

print("Time for execution: ", round(dtime, 1), ' minutes')
print("Done!")
print(" ")
# ------------------------------------------
