#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# This file is part of the BHC Complexity Toolkit.
#
# The BHC Complexity Toolkit is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The BHC Complexity Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the BHC Complexity Toolkit.  If not, 
# see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
# Copyright 2019, Mark D. Flood
#
# Author: Mark D. Flood
# Last revision: 16-Jul-2019
# -----------------------------------------------------------------------------

__all__ = ['make_dirs', 
           'make_dir_nic',
           'make_dir_fdiccb',
           'make_dir_fdicsod',
           'make_dir_fdicfail',
           'download_data', 
           'download_data_nic',
           'download_data_fdiccb',
           'download_data_fdicsod',
           'download_data_fdicfail',
           'main',
           ]
__version__ = '0.3'
__author__ = 'Mark D. Flood'


import os
import logging
import time
import urllib.request

import progressbar as pb

import bhc_util as UTIL
import bhc_data as DATA

LOG = UTIL.log_config(__file__.split(os.path.sep)[-1].split('.')[0])



def make_dirs(config):
    """Creates local data directories for all requested data downloads
    
    Four regulatory data downloads are possible:
        
        * NIC -- Bank holding company (BHC) structure data, based on the 
                 Federal Reserve's FR Y-10 data collection
        * FDIC CB -- History of banking institutions, including their
                 FDIC Cert identifiers and RSSD identifiers for their
                 regulatory high holders. Note that this sample does not
                 include all BHC subsidiaries, just the regulated 
                 depository institutions
        * FDIC SoD -- Summary of deposits data, including banks' and
                 their branches' FDIC Cert identifiers and RSSD identifiers 
                 for their regulatory high holders. Note that this sample 
                 does not include all BHC subsidiaries, just the regulated 
                 depository institutions
        * FDIC Fail -- FDIC historical list of failed (insured) depository
                 institutions, including their resolution costs
                 
    This covers the first step of the data-download process: creating 
    local target directories. Whether to process each of
    the three downloads is controlled by the "fetch" variables in the 
    configuration: nic_fetch, fdiccb_fetch, and fdicfail_fetch. 
    """
    if ('TRUE'==config['www2dat']['nic_fetch'].upper()):
        LOG.debug('Creating NIC directory: ')
        make_dir_nic(config)
    if ('TRUE'==config['www2dat']['fdiccb_fetch'].upper()):
        LOG.debug('Creating FDIC CB directory: ')
        make_dir_fdiccb(config)
    if ('TRUE'==config['www2dat']['fdicsod_fetch'].upper()):
        LOG.debug('Creating FDIC Summary of Deposits directory: ')
        make_dir_fdicsod(config)
    if ('TRUE'==config['www2dat']['fdicfail_fetch'].upper()):
        LOG.debug('Creating FDIC Fails directory: ')
        make_dir_fdicfail(config)

def make_dir_nic(config):
    nic_path = UTIL.resolve_dir_nic(config['www2dat']['nic_dir'], 
      config['www2dat']['nic_subdir'])
    os.makedirs(nic_path, exist_ok=True)
    LOG.info('NIC path: '+nic_path)

def make_dir_fdiccb(config):
    fdiccb_dir = config['www2dat']['fdiccb_dir']
    os.makedirs(fdiccb_dir, exist_ok=True)
    LOG.info('FDIC CB dir: '+fdiccb_dir)

def make_dir_fdicsod(config):
    fdicsod_dir = config['www2dat']['fdicsod_dir']
    os.makedirs(fdicsod_dir, exist_ok=True)
    LOG.info('FDIC Summary of Deposits dir: '+fdicsod_dir)

def make_dir_fdicfail(config):
    fdicfail_dir = config['www2dat']['fdicfail_dir']
    os.makedirs(fdicfail_dir, exist_ok=True)
    LOG.info('FDIC Fail dir: '+fdicfail_dir)


    
def download_data(config):
    """Downloads all requested data to local data directories 
    
    Four regulatory data downloads are possible:
        
        * NIC -- Bank holding company (BHC) structure data, based on the 
                 Federal Reserve's FR Y-10 data collection
        * FDIC CB -- History of banking institutions, including their
                 FDIC Cert identifiers and RSSD identifiers for their
                 regulatory high holders. Note that this sample does not
                 include all BHC subsidiaries, just the regulated 
                 depository institutions
        * FDIC SoD -- Summary of deposits data, including banks' and
                 their branches' FDIC Cert identifiers and RSSD identifiers 
                 for their regulatory high holders. Note that this sample 
                 does not include all BHC subsidiaries, just the regulated 
                 depository institutions
        * FDIC Fail -- FDIC historical list of failed (insured) depository
                 institutions, including their resolution costs
                 
    This covers the second step of the data-download process: actually 
    pulling down the data from the Internet. Whether to process each of
    the three downloads is controlled by the "fetch" variables in the 
    configuration: nic_fetch, fdiccb_fetch, and fdicfail_fetch. 
    """
    if ('TRUE'==config['www2dat']['nic_fetch'].upper()):
        LOG.debug('Downloading NIC data: ')
        download_data_nic(config)
    if ('TRUE'==config['www2dat']['fdiccb_fetch'].upper()):
        LOG.debug('Downloading FDIC CB data: ')
        download_data_fdiccb(config)
    if ('TRUE'==config['www2dat']['fdicsod_fetch'].upper()):
        LOG.debug('Downloading FDIC Summary of Deposits data: ')
        download_data_fdicsod(config)
    if ('TRUE'==config['www2dat']['fdicfail_fetch'].upper()):
        LOG.debug('Downloading FDIC Failures data: ')
        download_data_fdicfail(config)

def download_data_nic(config):
    nic_format = config['www2dat']['nic_format'].upper()
    tgtdir = UTIL.resolve_dir_nic(config['www2dat']['nic_dir'], 
      config['www2dat']['nic_subdir'])
    S = config['www2dat']
    if ('XML'==nic_format or 'BOTH'==nic_format):
        downloads = [
            (S['nic_xml_attributesactive'],S['nic_xmltgt_attributesactive']),
            (S['nic_xml_attributesbranch'],S['nic_xmltgt_attributesbranch']),
            (S['nic_xml_attributesclosed'],S['nic_xmltgt_attributesclosed']),
            (S['nic_xml_relationships'],   S['nic_xmltgt_relationships']),
            (S['nic_xml_transformations'], S['nic_xmltgt_transformations']) ]
        download_urllist(config, downloads, tgtdir)
    if ('CSV'==nic_format or 'BOTH'==nic_format):
        downloads = [
            (S['nic_csv_attributesactive'],S['nic_csvtgt_attributesactive']),
            (S['nic_csv_attributesbranch'],S['nic_csvtgt_attributesbranch']),
            (S['nic_csv_attributesclosed'],S['nic_csvtgt_attributesclosed']),
            (S['nic_csv_relationships'],   S['nic_csvtgt_relationships']),
            (S['nic_csv_transformations'], S['nic_csvtgt_transformations']) ]
        download_urllist(config, downloads, tgtdir)

def download_data_fdiccb(config):
    tgtdir = config['www2dat']['fdiccb_dir']
    cfgsection = config['www2dat']
    downloads = [
        (cfgsection['fdiccb_csv_8487'],cfgsection['fdiccb_csvtgt_8487']),
        (cfgsection['fdiccb_csv_8891'],cfgsection['fdiccb_csvtgt_8891']),
        (cfgsection['fdiccb_csv_9296'],cfgsection['fdiccb_csvtgt_9296']),
        (cfgsection['fdiccb_csv_9702'],cfgsection['fdiccb_csvtgt_9702']),
        (cfgsection['fdiccb_csv_0309'],cfgsection['fdiccb_csvtgt_0309']),
        (cfgsection['fdiccb_csv_1016'],cfgsection['fdiccb_csvtgt_1016']),
        (cfgsection['fdiccb_csv_1719'],cfgsection['fdiccb_csvtgt_1719']) ]
    download_urllist(config, downloads, tgtdir)

def download_data_fdicsod(config):
    tgtdir = config['www2dat']['fdicsod_dir']
    urllist = eval(config['www2dat']['fdicsod_url_filelist'])
    filelist = eval(config['www2dat']['fdicsod_ziptgt_filelist'])
    downloads = []
    # Assemble pairings of download URLs with normalized local filenames
    for i,url in enumerate(urllist):
        dnld = (urllist[i], filelist[i])
        downloads.append(dnld)
#    cfgsection = config['www2dat']
#    downloads = [
#        (cfgsection['fdicsod_csv_2018'],cfgsection['fdicsod_csvtgt_2018']),
#        (cfgsection['fdicsod_csv_2017'],cfgsection['fdicsod_csvtgt_2017']),
#        (cfgsection['fdicsod_csv_2016'],cfgsection['fdicsod_csvtgt_2016']),
#        (cfgsection['fdicsod_csv_2015'],cfgsection['fdicsod_csvtgt_2015']),
#        (cfgsection['fdicsod_csv_2014'],cfgsection['fdicsod_csvtgt_2014']),
#        (cfgsection['fdicsod_csv_2013'],cfgsection['fdicsod_csvtgt_2013']),
#        (cfgsection['fdicsod_csv_2012'],cfgsection['fdicsod_csvtgt_2012']),
#        (cfgsection['fdicsod_csv_2011'],cfgsection['fdicsod_csvtgt_2011']),
#        (cfgsection['fdicsod_csv_2010'],cfgsection['fdicsod_csvtgt_2010']),
#        (cfgsection['fdicsod_csv_2009'],cfgsection['fdicsod_csvtgt_2009']),
#        (cfgsection['fdicsod_csv_2008'],cfgsection['fdicsod_csvtgt_2008']),
#        (cfgsection['fdicsod_csv_2007'],cfgsection['fdicsod_csvtgt_2007']),
#        (cfgsection['fdicsod_csv_2006'],cfgsection['fdicsod_csvtgt_2006']),
#        (cfgsection['fdicsod_csv_2005'],cfgsection['fdicsod_csvtgt_2005']),
#        (cfgsection['fdicsod_csv_2004'],cfgsection['fdicsod_csvtgt_2004']),
#        (cfgsection['fdicsod_csv_2003'],cfgsection['fdicsod_csvtgt_2003']),
#        (cfgsection['fdicsod_csv_2002'],cfgsection['fdicsod_csvtgt_2002']),
#        (cfgsection['fdicsod_csv_2001'],cfgsection['fdicsod_csvtgt_2001']),
#        (cfgsection['fdicsod_csv_2000'],cfgsection['fdicsod_csvtgt_2000']),
#        (cfgsection['fdicsod_csv_1999'],cfgsection['fdicsod_csvtgt_1999']),
#        (cfgsection['fdicsod_csv_1998'],cfgsection['fdicsod_csvtgt_1998']),
#        (cfgsection['fdicsod_csv_1997'],cfgsection['fdicsod_csvtgt_1997']),
#        (cfgsection['fdicsod_csv_1996'],cfgsection['fdicsod_csvtgt_1996']),
#        (cfgsection['fdicsod_csv_1995'],cfgsection['fdicsod_csvtgt_1995']),
#        (cfgsection['fdicsod_csv_1994'],cfgsection['fdicsod_csvtgt_1994']) ]
    download_urllist(config, downloads, tgtdir)

def download_data_fdicfail(config):
    tgtdir = config['www2dat']['fdicfail_dir']
    cfgsection = config['www2dat']
    downloads = [
        (cfgsection['fdicfail_csv_curr'],cfgsection['fdicfail_csvtgt_curr']) ]
    download_urllist(config, downloads, tgtdir)



def download_urllist(config, downloads, tgtdir):
    sleep_interval = int(config['www2dat']['sleep_interval'])
    download_count = 0
    if ('TRUE'==config['DEFAULT']['progressbars'].upper()):
        for dld in pb.progressbar(downloads, redirect_stdout=True):
            download_count = download_count + 1
            url = dld[0]
            tgt = os.path.join(tgtdir, dld[1])
            if (os.path.isfile(tgt)):
                LOG.warning('Skipping download, target already exists: '+tgt)
            else:
                LOG.info('Downloading: '+url+' to: '+tgt)
                try:
                    urllib.request.urlretrieve(url, tgt)
                    if (download_count<len(downloads)):
                        time.sleep(sleep_interval)
                    DATA.test_zipfile_integrity(tgt)
                except (ConnectionError) as ce:
                    LOG.error('Download failed for: '+tgt+' from: '+url+'. '+
                      'You may need to download this file manually. '+str(ce))
    else:
        for dld in downloads:
            url = dld[0]
            tgt = os.path.join(tgtdir, dld[1])
            if (os.path.isfile(tgt)):
                LOG.warning('Skipping download, target already exists: '+tgt)
            else:
                LOG.info('Downloading: '+url+' to: '+tgt)
                try:
                    urllib.request.urlretrieve(url, tgt)
                    if (download_count<len(downloads)):
                        time.sleep(sleep_interval)
                    DATA.test_zipfile_integrity(tgt)
                except (ConnectionError) as ce:
                    LOG.error('Download failed for: '+tgt+' from: '+url+'. '+
                      'You may need to download this file manually. '+str(ce))



def main(argv=None):
    """A main function for command line execution
    
    This function parses the command line, loads the configuration, and 
    invokes the local functions:
        
         * make_dirs(config)
         * download_data(config)
         
    Parameters
    ----------
    argv : dict
        The collection of arguments submitted on the command line
    """
    config = UTIL.parse_command_line(argv, __file__)
    try:
        make_dirs(config)
        download_data(config)
    except Exception as e:
        logging.exception("message")
    LOG.info('**** Processing complete ****')
    
# This tests whether the module is being run from the command line
if __name__ == "__main__":
    main()