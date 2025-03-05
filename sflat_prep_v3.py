# NIRSpec S-Flat APT prepare

# Author: Nikolay Nikolov, AURA Associate Scientist, NIRSpec branch
# Pipeline Version**: 1.15.1
# Last Update**: November 13, 2024
# First Version**: February 1, 2024

# Conda environment (laptop): jwst.1.15.1

import sflat_red_utils

# Math and numpy
import math
import numpy as np

# Directory and file services e.g., unpacking
import os
import shutil
import subprocess

# Import requests, multithreading and multiprocessing
import requests
import multiprocessing

# LXML
from lxml import etree

# Pandas
import pandas as pd

apt_list = [1122]  #, 1485, 4457]#, 1485]
work_dir = './'
#work_dir = '/Users/nnikolov/Documents/Functional/NIRSpec/S-FLAT/20240208/' # laptop
#work_dir = '/home/nnikolov/MOS_SFlat/20240214/' # dlsci122

apt_dir  = work_dir + 'APT/'

do_step_download = True
do_parse_apt = True
do_filter_dfs = True
filter_mode = 'ALL' # Available options: 'ALL', 'GRATINGS', 'PRISM'
do_step_organize = True



if do_step_download:
    # -------------------------------------
    # Step 1: Create directory to download and store all APT files
    # -------------------------------------
    # Define APT directory

    # Check if APT directory exists and create it, if not
    isExist_APT = os.path.exists(apt_dir)

    if isExist_APT == False:
        # create a new empty APT directory to store apt and vsr files
        print("Creating a new APT directory ...")
        os.mkdir(apt_dir) # create the results directory




    prog_id = sflat_red_utils.list2np(apt_list)

    # convert the np array of numbers to an np array of strings
    prog_id = prog_id.astype(str)

    print(prog_id[0])



    apt_urls = []

    for i in range(len(prog_id)):
        apt_wl="https://www.stsci.edu/jwst/phase2-public/"+prog_id[i]+".aptx"
        apt_urls.append(apt_wl)

        response = requests.get(apt_urls[i], timeout=None)

        if response.status_code == 200:
            print("Success! Downloading file " + apt_urls[i][-9:-5] + "_APT.zip ...")
            open(apt_dir + apt_urls[i][-9:-5]+"_APT.zip", 'wb').write(response.content)

            # Unpack the zip file
            shutil.unpack_archive(apt_dir + apt_urls[i][-9:-5]+"_APT.zip", apt_dir)

            # Check if the APT zip file is there and delete it
            isZIPexist = os.path.exists(apt_dir + apt_urls[i][-9:-5]+"_APT.zip")

            # remove if the directory exists
            if isZIPexist == True:
                #print("Deleting APT ZIP ...")
                os.remove(apt_dir + apt_urls[i][-9:-5]+"_APT.zip")
                os.rename(apt_dir + apt_urls[i][-9:-5]+".xml", apt_dir + apt_urls[i][-9:-5]+"_APT.xml")
        else:
            print(f"Failed to retrieve {apt_urls[-9:-5]}_APT.xml. Status code: {response.status_code}")

    # Create lists of all APT and VSR files
    subprocess.run(["cd " + work_dir + "APT/ ; ls -1 *APT.xml > apt_list.txt"], shell=True)









if do_parse_apt:
    # -------------------------------------
    # Step 2: Parse all APT files and save as df
    # -------------------------------------
    prog_id = sflat_red_utils.list2np(apt_list)
    prog_id = prog_id.astype(str)


    for i in range(len(prog_id)):
        print("Parsing PID: ", prog_id[i])
        infile = apt_dir + prog_id[i] + "_APT.xml"


        # APT: Data Requests (observations info)
        tagname_of_node_to_parse = 'DataRequests'
        df_dr = sflat_red_utils.nn_parse_apt(prog_id[i], apt_dir + prog_id[i] + "_APT.xml", tagname_of_node_to_parse)
        df_dr.to_csv(apt_dir + 'APT_'+prog_id[i]+'.csv', index=False)









if do_filter_dfs:
    # -------------------------------------
    # Step 2: Read the data frames of the
    # individual PIDs and filter non-FLAT
    # exposures; define file names to be
    # ready for download
    # -------------------------------------
    prog_id = sflat_red_utils.list2np(apt_list)
    prog_id = prog_id.astype(str)

    for i in range(len(prog_id)):

        # Read each csv in a data frame
        df = pd.read_csv(apt_dir + 'APT_'+prog_id[i]+'.csv')

        # Filter the data frame with respect to needed quantities; keep ony FLAT*
        if (filter_mode == 'ALL'):
            filtered_df = df[ \
            (df['Lamp'] != 'REF') & \
            (df['Lamp'] != 'LINE1') & \
            (df['Lamp'] != 'LINE2') & \
            (df['Lamp'] != 'LINE3') & \
            (df['Lamp'] != 'LINE4') & \
            (df['Lamp'] != 'TEST') & \
            (df['OperatingMode'] == 'MSASPEC')]# & \
            #(df['MsaConfigFile'] != 'ALLCLOSED')]


        if (filter_mode == 'GRATINGS'):
            filtered_df = df[ \
            (df['Lamp'] != 'REF') & \
            (df['Lamp'] != 'LINE1') & \
            (df['Lamp'] != 'LINE2') & \
            (df['Lamp'] != 'LINE3') & \
            (df['Lamp'] != 'LINE4') & \
            (df['Lamp'] != 'TEST') & \
            (df['OperatingMode'] == 'MSASPEC') & \
            (df['MsaConfigFile'].str.len() < 17)] # this filters all LS PRISM configurations i.e., LS on both detectors #             (df['MsaConfigFile'] != 'ALLCLOSED') & \

        if (filter_mode == 'PRISM'):
            filtered_df = df[ \
            (df['Lamp'] != 'REF') & \
            (df['Lamp'] != 'LINE1') & \
            (df['Lamp'] != 'LINE2') & \
            (df['Lamp'] != 'LINE3') & \
            (df['Lamp'] != 'LINE4') & \
            (df['Lamp'] != 'TEST') & \
            (df['OperatingMode'] == 'MSASPEC') & \
            (df['MsaConfigFile'].str.len() == 26)] # this filters all GRATING configurations i.e., LS on one detector  #             (df['MsaConfigFile'] != 'ALLCLOSED') & \

        #print(filtered_df['MsaConfigFile'].values)

        # Obs 7, 8, 9, 10 are not on the MAST
        # Obs 17, 18, 19 were not executed
        # Remove these Observations
        new_df = filtered_df[ \
        (filtered_df['ObsNumber'] != 7) & \
        (filtered_df['ObsNumber'] != 8) & \
        (filtered_df['ObsNumber'] != 9) & \
        (filtered_df['ObsNumber'] != 10) & \
        (filtered_df['ObsNumber'] != 17) & \
        (filtered_df['ObsNumber'] != 18) & \
        (filtered_df['ObsNumber'] != 19) & \
        (filtered_df['ObsNumber'] != 25)]

        #a = new_df['OperatingMode'].values
        #for ik in range(len(a)):
        #    print(a[ik])

        # Add a new column with the file name for easier download
        # For details, check here: https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
        # 'jw01122002001_0210a_00001_nrs1'

        # Add column for the filename of NRS1
        new_df['DataFileNameNRS1'] = 'jw' + \
        new_df['ProposalID'].apply(lambda x: str(x).zfill(5)) + \
        new_df['ObsNumber'].apply(lambda x: str(x).zfill(3)) + \
        new_df['ObsVisit'].apply(lambda x: str(x).zfill(3)) + \
        '_021' + \
        new_df['ActivityNum'].apply(lambda x: str(sflat_red_utils.jwst_activity_naming_nn(x)).zfill(2)) + '_00001_nrs1'\

        # Add column for the filename of NRS2
        new_df['DataFileNameNRS2'] = 'jw' + \
        new_df['ProposalID'].apply(lambda x: str(x).zfill(5)) + \
        new_df['ObsNumber'].apply(lambda x: str(x).zfill(3)) + \
        new_df['ObsVisit'].apply(lambda x: str(x).zfill(3)) + \
        '_021' + \
        new_df['ActivityNum'].apply(lambda x: str(sflat_red_utils.jwst_activity_naming_nn(x)).zfill(2)) + '_00001_nrs2'\

        # Save the filtered data frame
        new_df.to_csv(apt_dir + 'APT_'+prog_id[i]+'_Filtered.csv', index=False)




if do_step_organize:
    # -------------------------------------
    # Read the filtered data frames
    # and organize
    # -------------------------------------
    prog_id = sflat_red_utils.list2np(apt_list)
    prog_id = prog_id.astype(str)

    calweb_code = 'sflat_calweb_v5.py'
    calweb_utils = 'sflat_red_utils.py'

    for i in range(len(prog_id)):

        # Read each csv in a data frame
        df = pd.read_csv(apt_dir + 'APT_'+prog_id[i]+'_Filtered.csv')

        all_dispersers = df['Grating'].values
        unique_dispersers = np.unique(all_dispersers)
        # print(unique_dispersers)

        # Obtsin data frames only for a given disperser
        for disperser in unique_dispersers:
            #print(disperser)

            # PRISM
            if (disperser == 'PRISM'):
                df_disperser = df[(df['Grating'] == str(disperser)) & (df['MsaConfigFile'].str.len() < 28)]

                df_disperser.to_csv(apt_dir + prog_id[i]+'_' + disperser + '.csv', index=False)

                # Check if a directory for the prism exists
                disp_dir = work_dir + disperser
                isExist_csv = os.path.exists(disp_dir)
                if isExist_csv == False:
                    os.mkdir(disp_dir)
                    subprocess.run(["cp " + apt_dir + prog_id[i] + '_' + disperser + '.csv '+ disp_dir + '/'], shell=True)
                    subprocess.run(["cp " + work_dir + calweb_code + ' ' + disp_dir + '/'], shell=True)
                    subprocess.run(["cp " + work_dir + calweb_utils + ' ' + disp_dir + '/'], shell=True)
                    subprocess.run(["mv " + disp_dir + '/' + calweb_code + ' ' + disp_dir + '/' + disperser + '_' + calweb_code], shell=True)

            if (disperser != 'PRISM'):
                df_disperser = df[(df['Grating'] == str(disperser)) & (df['MsaConfigFile'].str.len() < 17)]

                # G140M&H
                if ((disperser == 'G140M') or (disperser == 'G140H')):

                    if ((disperser == 'G140M') or (disperser == 'G140H')):
                        all_lamps = df_disperser['Lamp'].values
                        unique_lamps = np.unique(all_lamps)

                        for lamp_id in unique_lamps:
                            # print(disperser, ' ', str(lamp_id))
                            df_lamp = df_disperser[df_disperser['Lamp'] == lamp_id]
                            df_lamp.to_csv(apt_dir + prog_id[i] + '_' + disperser + '_' + lamp_id + '.csv', index=False)

                            # Check if a directory for the prism exists
                            disp_dir = work_dir + disperser + '_' + lamp_id
                            isExist_csv = os.path.exists(disp_dir)
                            if isExist_csv == False:
                                os.mkdir(disp_dir)
                                subprocess.run(["cp " + apt_dir + prog_id[i] + '_' + disperser + '_' + lamp_id + '.csv '+ disp_dir + '/'], shell=True)
                                subprocess.run(["cp " + work_dir + calweb_code + ' ' + disp_dir + '/'], shell=True)
                                subprocess.run(["cp " + work_dir + calweb_utils + ' ' + disp_dir + '/'], shell=True)
                                subprocess.run(["mv " + disp_dir + '/' + calweb_code + ' ' + disp_dir + '/' + disperser + '_' + lamp_id + '_' + calweb_code], shell=True)

            # G235H&M and G395H&M
            if ((disperser != 'G140M') & (disperser != 'G140H')):
                df_disperser.to_csv(apt_dir + prog_id[i] + '_' + disperser + '.csv', index=False)

                # Check if a directory for the prism exists
                disp_dir = work_dir + disperser
                isExist_csv = os.path.exists(disp_dir)
                if isExist_csv == False:
                    os.mkdir(disp_dir)
                    subprocess.run(["cp " + apt_dir + prog_id[i] + '_' + disperser + '.csv '+ disp_dir + '/'], shell=True)
                    subprocess.run(["cp " + work_dir + calweb_code + ' ' + disp_dir + '/'], shell=True)
                    subprocess.run(["cp " + work_dir + calweb_utils + ' ' + disp_dir + '/'], shell=True)
                    subprocess.run(["mv " + disp_dir + '/' + calweb_code + ' ' + disp_dir + '/' + disperser + '_' + calweb_code], shell=True)

print('Done!')
