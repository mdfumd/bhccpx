#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Copyright (c) 2020 University of Maryland # All rights reserved.
# -----------------------------------------------------------------------------
#
# This file is part of the BHC Complexity Toolkit.
#
# The toolkit is part of a larger research project undertaken by:
#   * Mark D. Flood, U. of Maryland
#   * Dror Kenett, John Hopkins U. and London School of Economics
#   * Robin Lumsdaine, American U., Erasmus U., and Tinbergen Inst.
#   * Jonathan K. Simon, U. of Iowa
#
# Users of the software should cite the following research paper:
#
#   * M. Flood, D. Kenett, R. Lumsdaine, and J. Simon
#     The Complexity of Bank Holding Companies: A Topological Approach
#     Journal of Banking and Finance, 2020, forthcoming
#     https://doi.org/10.1016/j.jbankfin.2020.105789
#
#     Abstract:
#     We develop metrics to assess the complexity of a bank holding
#     company (BHC), based on its ownership structure. Large BHCs have
#     intricate ownership hierarchies involving hundreds or even thousands
#     of legal entities that may contribute to increased operational risk
#     and greater opacity. Our measures are mathematically grounded,
#     intuitive, and easy to implement. They may be particularly useful
#     in the context of resolution, where regulators often face significant
#     time pressure and coordination challenges. We use regulatory filing
#     data from the Federal Reserve to validate the measures, demonstrating
#     that they provide a useful complement to balance sheet information in
#     assessing BHC complexity. Notably, the proposed measures are highly
#     correlated with existing complexity indicators that are not based on
#     organizational structure and are less correlated with size than these
#     existing complexity measures. We show that the proposed measures
#     provide additional explanatory power for the regulatory indicators,
#     even after controlling for size.
#
# A preprint is available at:
#
#   * The Complexity of Bank Holding Companies: A Topological Approach
#     https://ssrn.com/abstract=3031726
#
# -----------------------------------------------------------------------------

import os
import logging
import time
import urllib.request

import progressbar as pb

import bhc_util as UTIL
import bhc_data as DATA


__all__ = ['make_dirs', 
           'make_dir_cache',
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
__version__ = '0.6'
__author__ = 'Mark D. Flood'


MODNAME = __file__.split(os.path.sep)[-1].split('.')[0]
LOG = UTIL.log_config(MODNAME)


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
    LOG.debug('Creating cache directory: ')
    make_dir_cache(config)
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

def make_dir_cache(config):
    cache_dir = config['DEFAULT']['cachedir']
    os.makedirs(cache_dir, exist_ok=True)
    LOG.info('Cache dir: '+cache_dir)

def make_dir_nic(config):
    nic_path = DATA.resolve_dir_nic(config['www2dat']['nic_dir'], 
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
    tgtdir = DATA.resolve_dir_nic(config['www2dat']['nic_dir'], 
      config['www2dat']['nic_subdir'])
    if ('XML'==nic_format or 'BOTH'==nic_format):
        urllist = eval(config['www2dat']['nic_xml_urllist'])
        filelist = eval(config['www2dat']['nic_xmlzip_filelist'])
        download_urllist(config, urllist, filelist, tgtdir)
    if ('CSV'==nic_format or 'BOTH'==nic_format):
        urllist = eval(config['www2dat']['nic_csv_urllist'])
        filelist = eval(config['www2dat']['nic_csvzip_filelist'])
        download_urllist(config, urllist, filelist, tgtdir)

def download_data_fdiccb(config):
    tgtdir = config['www2dat']['fdiccb_dir']
    urllist = eval(config['www2dat']['fdiccb_urllist'])
    filelist = eval(config['www2dat']['fdiccb_zip_filelist'])
    download_urllist(config, urllist, filelist, tgtdir)

def download_data_fdicsod(config):
    tgtdir = config['www2dat']['fdicsod_dir']
    urllist = eval(config['www2dat']['fdicsod_urllist'])
    filelist = eval(config['www2dat']['fdicsod_zip_filelist'])
    download_urllist(config, urllist, filelist, tgtdir)

def download_data_fdicfail(config):
    tgtdir = config['www2dat']['fdicfail_dir']
    urllist = eval(config['www2dat']['fdicfail_urllist'])
    filelist = eval(config['www2dat']['fdicfail_csv_filelist'])
    download_urllist(config, urllist, filelist, tgtdir, testzip=False)



def download_urllist(config, urllist, filelist, tgtdir, testzip=True):
    sleep_interval = int(config['www2dat']['sleep_interval'])
    # Assemble pairings of download URLs with normalized local filenames
    downloads = []
    for i,url in enumerate(urllist):
        dnld = (urllist[i], filelist[i])
        downloads.append(dnld)
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
                    if (testzip):
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
                    if (testzip):
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