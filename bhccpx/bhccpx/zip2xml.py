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
import zipfile as zf
import progressbar as pb

import bhc_util as UTIL
import bhc_data as DATA

LOG = UTIL.log_config(__file__.split(os.path.sep)[-1].split('.')[0])

    
def unzip_data(config):
    """Unzips downloaded data within the local data directories 
    
    Three regulatory data downloads require unzipping:
        
        * NIC -- Bank holding company (BHC) structure data, based on the 
                 Federal Reserve's FR Y-10 data collection
        * FDIC SoD -- Summary of deposits data, including banks' and
                 their branches' FDIC Cert identifiers and RSSD identifiers 
                 for their regulatory high holders. Note that this sample 
                 does not include all BHC subsidiaries, just the regulated 
                 depository institutions
        * FDIC CB -- History of banking institutions, including their
                 FDIC Cert identifiers and RSSD identifiers for their
                 regulatory high holders. Note that this sample does not
                 include all BHC subsidiaries, just the regulated 
                 depository institutions
                 
    This is the final step of the data download process. 
    """
    if ('TRUE'==config['zip2xml']['nic_unzip'].upper()):
        LOG.debug('Unzipping NIC data: ')
        unzip_data_nic(config)
    if ('TRUE'==config['zip2xml']['fdicsod_unzip'].upper()):
        LOG.debug('Unzipping FDIC Summary of Deposits data: ')
        unzip_data_fdicsod(config)
    if ('TRUE'==config['zip2xml']['fdiccb_unzip'].upper()):
        LOG.debug('Unzipping FDIC CB data: ')
        unzip_data_fdiccb(config)


def unzip_data_nic(config):
    """Unzips the NIC files, as specified in the configuration
    """
    tgtdir = DATA.resolve_dir_nic(config['zip2xml']['nic_dir'], 
      config['zip2xml']['nic_subdir'])
    nic_format = config['zip2xml']['nic_format'].upper()
    pbar = ('TRUE'==config['DEFAULT']['progressbars'].upper())
    force = ('TRUE'==config['zip2xml']['force'].upper())
    if ('XML'==nic_format or 'BOTH'==nic_format):
        zipfiles = [
            config['zip2xml']['nic_xmltgt_attributesactive'],
            config['zip2xml']['nic_xmltgt_attributesbranch'],
            config['zip2xml']['nic_xmltgt_attributesclosed'],
            config['zip2xml']['nic_xmltgt_relationships'],
            config['zip2xml']['nic_xmltgt_transformations'] ]
        unzip_filelist(tgtdir, zipfiles, force, pbar)
    if ('CSV'==nic_format or 'BOTH'==nic_format):
        zipfiles = [
            config['zip2xml']['nic_csvtgt_attributesactive'],
            config['zip2xml']['nic_csvtgt_attributesbranch'],
            config['zip2xml']['nic_csvtgt_attributesclosed'],
            config['zip2xml']['nic_csvtgt_relationships'],
            config['zip2xml']['nic_csvtgt_transformations'] ]
        normfiles = [
            config['zip2xml']['nic_csvnrm_attributesactive'],
            config['zip2xml']['nic_csvnrm_attributesbranch'],
            config['zip2xml']['nic_csvnrm_attributesclosed'],
            config['zip2xml']['nic_csvnrm_relationships'],
            config['zip2xml']['nic_csvnrm_transformations'] ]
        unzip_filelist(tgtdir, zipfiles, force, pbar, normfiles)
        scrub_headers_nic(tgtdir, pbar, normfiles)



def unzip_data_fdiccb(config):
    """Unzips the FDIC CB files, as specified in the configuration
    """
    tgtdir = config['zip2xml']['fdiccb_dir']
    zipfiles = [
        config['zip2xml']['fdiccb_csvtgt_8487'],
        config['zip2xml']['fdiccb_csvtgt_8891'],
        config['zip2xml']['fdiccb_csvtgt_9296'],
        config['zip2xml']['fdiccb_csvtgt_9702'],
        config['zip2xml']['fdiccb_csvtgt_0309'],
        config['zip2xml']['fdiccb_csvtgt_1016'],
        config['zip2xml']['fdiccb_csvtgt_1719'] ]
    normfiles = [
        config['zip2xml']['fdiccb_csvnrm_8487'],
        config['zip2xml']['fdiccb_csvnrm_8891'],
        config['zip2xml']['fdiccb_csvnrm_9296'],
        config['zip2xml']['fdiccb_csvnrm_9702'],
        config['zip2xml']['fdiccb_csvnrm_0309'],
        config['zip2xml']['fdiccb_csvnrm_1016'],
        config['zip2xml']['fdiccb_csvnrm_1719'] ]
    force = ('TRUE'==config['zip2xml']['force'].upper())
    pbar = ('TRUE'==config['DEFAULT']['progressbars'].upper())
    unzip_filelist(tgtdir, zipfiles, force, pbar, normfiles)



def unzip_data_fdicsod(config):
    """Unzips the FDIC SoD files, as specified in the configuration
    """
    tgtdir = config['zip2xml']['fdicsod_dir']
    zipfiles = eval(config['zip2xml']['fdicsod_ziptgt_filelist'])
#    zipfiles = [
#        config['zip2xml']['fdicsod_csvtgt_2018'],
#        config['zip2xml']['fdicsod_csvtgt_2017'],
#        config['zip2xml']['fdicsod_csvtgt_2016'],
#        config['zip2xml']['fdicsod_csvtgt_2015'],
#        config['zip2xml']['fdicsod_csvtgt_2014'],
#        config['zip2xml']['fdicsod_csvtgt_2013'],
#        config['zip2xml']['fdicsod_csvtgt_2012'],
#        config['zip2xml']['fdicsod_csvtgt_2011'],
#        config['zip2xml']['fdicsod_csvtgt_2010'],
#        config['zip2xml']['fdicsod_csvtgt_2009'],
#        config['zip2xml']['fdicsod_csvtgt_2008'],
#        config['zip2xml']['fdicsod_csvtgt_2007'],
#        config['zip2xml']['fdicsod_csvtgt_2006'],
#        config['zip2xml']['fdicsod_csvtgt_2005'],
#        config['zip2xml']['fdicsod_csvtgt_2004'],
#        config['zip2xml']['fdicsod_csvtgt_2003'],
#        config['zip2xml']['fdicsod_csvtgt_2002'],
#        config['zip2xml']['fdicsod_csvtgt_2001'],
#        config['zip2xml']['fdicsod_csvtgt_2000'],
#        config['zip2xml']['fdicsod_csvtgt_1999'],
#        config['zip2xml']['fdicsod_csvtgt_1998'],
#        config['zip2xml']['fdicsod_csvtgt_1997'],
#        config['zip2xml']['fdicsod_csvtgt_1996'],
#        config['zip2xml']['fdicsod_csvtgt_1995'],
#        config['zip2xml']['fdicsod_csvtgt_1994'] ]
    force = ('TRUE'==config['zip2xml']['force'].upper())
    pbar = ('TRUE'==config['DEFAULT']['progressbars'].upper())
    unzip_filelist(tgtdir, zipfiles, force, pbar)



def scrub_headers_nic(tgtdir, pbar, normfiles):
    """Scrubs the header row for each of a list of CSV files
    
    The CSV files available for download from the FFIEC site include a 
    special character ('#') at the start of the header row, making the 
    first column name incompatible with the field names in the XML 
    download. This routine scrubs the CSV header row to remove the 
    special character, so that the CSV and XML files are consistent.
    
    Parameters
    ----------
    tgtdir : directory name
        Local file system location for the CSV files, and also where the
        unzipped contents will land
    pbar : bool
        Whether to display a progress bar in the console
    normfiles : list
        List of file names to be scrubbed
    """
    if (pbar):
        for i in pb.progressbar(range(len(normfiles)), redirect_stdout=True):
            normf = normfiles[i]
            normfilepath = os.path.join(tgtdir, normf)
            LOG.debug('Scrubbing header row: '+normfilepath)
            DATA.sed(normfilepath, '#', '', N=1)
    else:
        for i in range(len(normfiles)):
            normf = normfiles[i]
            normfilepath = os.path.join(tgtdir, normf)
            LOG.debug('Scrubbing header row: '+normfilepath)
            DATA.sed(normfilepath, '#', '', N=1)



def unzip_filelist(tgtdir, zipfiles, force, pbar, normfiles=None):
    """Unzips each of a list of zip files in a target directory
    
    Parameters
    ----------
    zipfiles : list
        List of file names to be unzipped
    tgtdir : directory name
        Local file system location for the zip files, and also where the
        unzipped contents will land
    force: bool
        Whether to force unzipping a given file, even if the target file
        already exists
    pbar : bool
        Whether to display a progress bar in the console
    normfiles : list
        List of normalized filenames to apply to the unzipped files
    """
    if (pbar):
        for i in pb.progressbar(range(len(zipfiles)), redirect_stdout=True):
            zipf = zipfiles[i]
            zipfilepath = os.path.join(tgtdir, zipf)
            zipfile = zf.ZipFile(zipfilepath, 'r')
            # Sort names, because each FDIC SoD zip has multiple files
            tgtfilename = sorted(zipfile.namelist())[0]
            tgtfilepath = os.path.join(tgtdir, tgtfilename)
            if not(None==normfiles):
                nrmfilepath = os.path.join(tgtdir, normfiles[i])
                if (os.path.isfile(nrmfilepath) and not(force)):
                    LOG.warning('Skipping unzip, target already exists: '+
                                zipfilepath)
                else:
                    LOG.info('Unzipping: '+zipfilepath+
                             ' to: '+nrmfilepath)
                    zipfile.extract(tgtfilename, tgtdir)
                    LOG.debug('Renaming '+tgtfilepath+' to '+nrmfilepath)
                    os.rename(tgtfilepath, nrmfilepath)
            else:
                if (os.path.isfile(tgtfilepath) and not(force)):
                    LOG.warning('Skipping unzip, target already exists: '+
                                zipfilepath)
                else:
                    LOG.info('Unzipping: '+zipfilepath+' to: '+tgtfilepath)
                    zipfile.extract(tgtfilename, tgtdir)
    else:
        for i in range(len(zipfiles)):
            zipf = zipfiles[i]
            zipfilepath = os.path.join(tgtdir, zipf)
            zipfile = zf.ZipFile(zipfilepath, 'r')
            # Sort names, because each FDIC SoD zip has multiple files
            tgtfilename = sorted(zipfile.namelist())[0]
            tgtfilepath = os.path.join(tgtdir, tgtfilename)
            if not(None==normfiles):
                nrmfilepath = os.path.join(tgtdir, normfiles[i])
                if (os.path.isfile(nrmfilepath) and not(force)):
                    LOG.warning('Skipping unzip, target already exists: '+
                                zipfilepath)
                else:
                    LOG.info('Unzipping: '+zipfilepath+
                             ' to: '+nrmfilepath)
                    zipfile.extract(tgtfilename, tgtdir)
                    LOG.debug('Renaming '+tgtfilepath+' to '+nrmfilepath)
                    os.rename(tgtfilepath, nrmfilepath)
            else:
                if (os.path.isfile(tgtfilepath) and not(force)):
                    LOG.warning('Skipping unzip, target already exists: '+
                                zipfilepath)
                else:
                    LOG.info('Unzipping: '+zipfilepath+
                             ' to: '+tgtfilepath)
                    zipfile.extract(tgtfilename, tgtdir)



def main(argv=None):
    """A main function for command line execution
    
    This function parses the command line, loads the configuration, and 
    invokes the local function:
        
         * unzip_data(config)
         
    Parameters
    ----------
    argv : dict
        The collection of arguments submitted on the command line
    """
    config = UTIL.parse_command_line(argv, __file__)
    try:
        unzip_data(config)
    except Exception as e:
        logging.exception("message")
    LOG.info('**** Processing complete ****')    
    
# This tests whether the module is being run from the command line
if __name__ == "__main__":
    main()