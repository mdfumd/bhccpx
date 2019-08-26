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
# Last revision: 22-Jun-2019
# -----------------------------------------------------------------------------

__all__ = ['resolve_dir_nic',
           'ATTcsv2df',
           'RELcsv2df',
           'FDICCBcsv2df',
           'FDICSoDcsv2df',
           'FDICFailcsv2df',
           'fetch_data',
           'create_NIC',
           'create_FDICSoD',
           'create_FDICCB',
           'create_FDICFail',
           'NIC_highholders',
           ]
__version__ = '0.3'
__author__ = 'Mark D. Flood'


import os
import time
import logging
import re
import shutil
import tempfile

import zipfile as zf
import numpy as np 
import pandas as pd
import _pickle as pik 

import bhc_util as UTIL

LOG = logging.getLogger(__file__.split(os.path.sep)[-1].split('.')[0])


# Standard prefixes to use for cached objects
case_NIC = 'NIC'
case_FDICSoD = 'FDICSoD'
case_FDICCB = 'FDICCB'
case_FDICFail = 'FDICFail'
case_IDmaps = 'IDmaps'
case_BankSys = 'BankSys'
def fetch_data(spec, case):
    data = None
    asofdate = spec['asofdate']
    cachedir = spec['cachedir']
    cachefilename = f'{case}_{normalize_cache_date(asofdate, case)}.pik'
    cachepath = os.path.join(cachedir, cachefilename)
    # Clear the cache, if requested
    if (case_NIC==case and spec['nic_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    elif (case_FDICCB==case and spec['fdiccb_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    elif (case_FDICSoD==case and spec['fdicsod_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    elif (case_FDICFail==case and spec['fdicfail_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    elif (case_IDmaps==case and spec['idmaps_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    elif (case_BankSys==case and spec['banksys_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    # Retrieve from cache, if available
    if (os.path.isfile(cachepath)):
        f = open(cachepath, 'rb')
        data = pik.load(f)
        f.close()
    # If not available, create the cache file
    else:
        if (case_NIC==case):
            indir = spec['nic_dir']
            fA = spec['nic_attributesactive']
            fB = spec['nic_attributesbranch']
            fC = spec['nic_attributesclosed']
            fREL = spec['nic_relationships']
            sep = spec['nic_sep']
            filt = spec['nic_filterasof']
            data = create_NIC(indir, fA, fB, fC, fREL, asofdate, sep, filt)
        elif (case_FDICCB==case):
            indir = spec['fdiccb_dir']
            filelist = spec['fdiccb_filelist']
            data = create_FDICCB(indir, filelist, asofdate)
        elif (case_FDICSoD==case):
            indir = spec['fdicsod_dir']
            filelist = spec['fdicsod_filelist']
            data = create_FDICSoD(indir, filelist, asofdate)
        elif (case_FDICFail==case):
            indir = spec['fdicfail_dir']
            filelist = spec['fdicfail_file']
            data = create_FDICFail(indir, filelist, asofdate)
        elif (case_IDmaps==case):
            LOG.error(f'Cached IDmaps not found at {cachepath}')
        elif (case_BankSys==case):
            LOG.error(f'Cached BankSys not found at {cachepath}')
        if (data is not None):
            f = open(cachepath, 'wb')
            pik.dump(data, f)
            f.close()
    return data

def cache_data(spec, case, data):
    asofdate = spec['asofdate']
    cachedir = spec['cachedir']
    cachefilename = f'{case}_{normalize_cache_date(asofdate, case)}.pik'
    cachepath = os.path.join(cachedir, cachefilename)
    # Clear the cache, if requested
    if (case_IDmaps==case and spec['idmaps_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    elif (case_BankSys==case and spec['banksys_clearcache']):
        if os.path.isfile(cachepath): os.remove(cachepath)
    elif (case_IDmaps==case):
        LOG.error(f'Cached IDmaps not found at {cachepath}')
    elif (case_BankSys==case):
        LOG.error(f'Cached BankSys not found at {cachepath}')
    if (data is not None):
        f = open(cachepath, 'wb')
        pik.dump(data, f)
        f.close()
    return data

# =============================================================================
#   Cache management
# =============================================================================
def clear_cache(cachedir, case, YQ0, YQ1):
    """Clears the cache of any pickled files for chosen prefix and date range
    """
    LOG.debug(f'Clearing cache of {case}_*.pik files for range: {YQ0}-{YQ1}')
    asof_list = UTIL.assemble_asofs(YQ0, YQ1)
    for asofdate in asof_list:
        filename = f'{case}_{asofdate}.pik'
        filepath = os.path.join(cachedir, filename)
        if os.path.isfile(filepath):
            os.remove(filepath)
            


def resolve_dir_nic(nic_dir, nic_subdir):
    if ('NONE'==nic_subdir.upper()):
        today = int(time.strftime("%Y%m%d"))
        qtrend = UTIL.rcnt_qtrend(today)
        nic_subdir = UTIL.stringify_qtrend(qtrend)
        LOG.debug('Generated NIC subdir: '+nic_subdir)
    else:
        nic_subdir = nic_subdir
        LOG.debug('Provided NIC subdir: '+nic_subdir)
    nic_path = os.path.join(nic_dir, nic_subdir)
    LOG.debug('Resolved NIC directory: '+nic_path)
    return nic_path


def normalize_cache_date(asofdate, case):
    norm_date = 10000101
    if (case_NIC==case):
        norm_date = UTIL.rcnt_qtrend(asofdate)
    elif (case_FDICSoD==case):
        norm_date = UTIL.rcnt_qtrend(asofdate)
    elif (case_FDICCB==case):
        norm_date = UTIL.rcnt_qtrend(asofdate)
    elif (case_FDICFail==case):
        norm_date = UTIL.rcnt_qtrend(asofdate)
    elif (case_IDmaps==case):
        norm_date = UTIL.rcnt_qtrend(asofdate)
    elif (case_BankSys==case):
        norm_date = UTIL.rcnt_qtrend(asofdate)
    else:
        LOG.error('Invalid cache case: '+case)
    return norm_date
    
# =============================================================================
#   NIC data
# =============================================================================
# Mnemonic indices for the NIC data list
IDX_Attributes = 0
IDX_Relationships = 1
IDX_HighHolder = 2
IDX_Entities = 3
IDX_Parents = 4
IDX_Offspring = 5
def create_NIC(indir, fA, fB, fC, fRel, asofdate, sep=',', filter_asof=True):
    """Parses the NIC data and creates a list of DataFrames for the asofdate

    Assembles the NIC data for given asofdate into a single list object.
    The returned list contains pointers to six objects, indexed as follows:
        
       0. IDX_Attributes - NIC attributes, combining active, branch, and closed
       1. IDX_Relationships - NIC relationships
       2. IDX_HighHolder - List of the high holder BHCs in the system
       3. IDX_Entities - Map (dict)
       4. IDX_Parents       = 4
       5. IDX_Offspring     = 5
       
    """
    data = []
    # First, populate DATA with the raw info from the CSV files:    
#    csvfilepathA = os.path.join(indir, file_attA)
#    csvfilepathB = os.path.join(indir, file_attB)
#    csvfilepathC = os.path.join(indir, file_attC)
#    ATTdf_a = ATTcsv2df(csvfilepathA, asofdate, 'A')
#    ATTdf_b = ATTcsv2df(csvfilepathB, asofdate, 'B')
#    ATTdf_c = ATTcsv2df(csvfilepathC, asofdate, 'C')
#    ATTdf = makeATTs(indir, file_attA, file_attB, file_attC, asofdate)
#def makeATTs(indir, fA, fB, fC, asofdate, sep=',', filter_asof=False):
    csvfilepathA = os.path.join(indir, fA)
    csvfilepathB = os.path.join(indir, fB)
    csvfilepathC = os.path.join(indir, fC)
    print(f'Active: {csvfilepathA}, Branch: {csvfilepathB}, Closed: {csvfilepathC}')
    ATTdf_a = ATTcsv2df(csvfilepathA, asofdate, 'A', sep, filter_asof)
    ATTdf_b = ATTcsv2df(csvfilepathB, asofdate, 'B', sep, filter_asof)
    ATTdf_c = ATTcsv2df(csvfilepathC, asofdate, 'C', sep, filter_asof)
    print(f'{len(ATTdf_a)} Active, {len(ATTdf_b)} Branch, {len(ATTdf_c)} Closed')
    ATTdf = pd.concat([ATTdf_a, ATTdf_b, ATTdf_c])
    data.insert(IDX_Attributes, ATTdf)
    csvfilepathR = os.path.join(indir, fRel)
    RELdf = RELcsv2df(csvfilepathR, asofdate)
    data.insert(IDX_Relationships, RELdf)
    # Then, add derived structures based on the Relationships data:    
    derived_data = NIC_highholders(RELdf, asofdate)
    data.insert(IDX_HighHolder, derived_data[0])
    data.insert(IDX_Entities, derived_data[1])
    data.insert(IDX_Parents, derived_data[2])
    data.insert(IDX_Offspring, derived_data[3])
    return data

def NIC_highholders(RELdf, asofdate):
    """Creates four key derived objects from a NIC relationships dataframe
    
    The derived objects are:
      * entities is the set of all NIC nodes appearing the relationships, 
        either as parents or offspring
      * parents is a dictionary, keyed by individual node_ids, of the 
        set of immediate parents of each node
      * offspring is a dictionary, keyed by individual node_ids, of the 
        set of immediate offspring of each node
      * high_holders is the set of all high-holder entities, defined as any 
        node with no immediate parent 
    
    Note that a high-holder node will have an entry in the parents dict, but
    this entry will point to an empty set (high holders have no parents)
    """
#    (ID_RSSD_PARENT, ID_RSSD_OFFSPRING, DT_START, DT_END) = REL_IDcols(RELdf)
    # Create some containers for derived structures
    parents = {}     # Dictionary of immediate parents (a set) for each node
    offspring = {}   # Dictionary of immediate children (a set) for each node
    entities = set()
    high_holders = set()
    # Loop through Relationships to assemble entities, parents, and offspring
    for row in RELdf.iterrows():
        date0 = int(row[1]['DT_START'])
        date1 = int(row[1]['DT_END'])
        rssd_par = row[1]['ID_RSSD_PARENT']
        rssd_off = row[1]['ID_RSSD_OFFSPRING']
#    for row in RELdf.iterrows():
#        date0 = int(row[0][DT_START])
#        date1 = int(row[0][DT_END])
#        rssd_par = row[0][ID_RSSD_PARENT]
#        rssd_off = row[0][ID_RSSD_OFFSPRING]
        if (asofdate < date0 or asofdate > date1):
            LOG.warning(f'Bad asofdate: {asofdate}; '+
                'out of bounds: {rssd_par}, {rssd_off}, {date0}, {date1}')
            continue   
        entities.add(rssd_par)
        try:
            offspring[rssd_par].add(rssd_off)
        except KeyError:
            offspring[rssd_par] = set()
            offspring[rssd_par].add(rssd_off)
        entities.add(rssd_off)
        try:
            parents[rssd_off].add(rssd_par)
        except KeyError:
            parents[rssd_off] = set()
            parents[rssd_off].add(rssd_par)
    # Filter entities to find the high_holders
    for ent in entities:
        try:
            len(parents[ent])      # Count the parents, if they exist
        except KeyError:
            high_holders.add(ent)  # High holders are those w/zero parents
    return high_holders, entities, parents, offspring


def ATTcsv2df(csvfile, asofdate, nicsource, sep=',', filter_asof=False):
    """# Reads NIC attributes from a CSV file
    
    Note that NIC downloads start in XML format; you must convert this from
    XML to tab-delimited CSV before using this function.
    The contents of the CSV file are converted to appropriate primitive types
    and stored in a Pandas dataframe, which is returned. 
    The returned dataframe is indexed on one field: ID_RSSD.
    This function also adds a 'NICsource' column (not in the source CSV) 
    to the dataframe, which indicates the nature (A/B/C) of the NIC source.

    Parameters
    ----------
    csvfile : str
        Location of the CSV input file - an open, readable pointer to a 
        tab-delimited CSV file containing a NIC attributes download
    asofdate : int
        The as-of date (in form form YYYYMMDD) to seek in the NIC data. 
        When filtering, NIC observations are excluded if 
        DT_END < asofdate or DT_START > asofdate 
    nicsource : a single character that indicating the nature of the source:
        * 'A' indicates an "active" or going-concern node
        * 'B' indicates a "branch" of an active node; not a distinct entity
        * 'C' indicates a "closed" or "inactive" node
    """
    DTYPES_ATT = {
        'ACT_PRIM_CD':object, 
        'AUTH_REG_DIST_FRS':np.int8, 
        'BHC_IND':np.int8, 
        'BNK_TYPE_ANALYS_CD':np.int8, 
        'BROAD_REG_CD':np.int8, 
        'CHTR_AUTH_CD':np.int8, 
        'CHTR_TYPE_CD':np.int16, 
        'CITY':object, 
        'CNSRVTR_CD':np.int8, 
        'CNTRY_CD':np.int32, 
        'CNTRY_INC_CD':np.int32, 
        'CNTRY_INC_NM':object, 
        'CNTRY_NM':object, 
        'COUNTY_CD':np.int32, 
        'DIST_FRS':np.int8, 
        'DOMESTIC_IND':object, 
        'DT_END':np.int32, 
        'DT_EXIST_CMNC':np.int32, 
        'DT_EXIST_TERM':np.int32, 
        'DT_INSUR':np.int32, 
        'DT_OPEN':np.int32, 
        'DT_START':np.int32, 
        'D_DT_END':object, 
        'D_DT_EXIST_CMNC':object, 
        'D_DT_EXIST_TERM':object, 
        'D_DT_INSUR':object, 
        'D_DT_OPEN':object, 
        'D_DT_START':object, 
        'ENTITY_TYPE':object, 
        'EST_TYPE_CD':np.int8, 
        'FBO_4C9_IND':np.int8, 
        'FHC_IND':np.int8, 
        'FISC_YREND_MMDD':np.int16, 
        'FNCL_SUB_HOLDER':np.int8, 
        'FNCL_SUB_IND':np.int8, 
        'FUNC_REG':np.int8, 
        'IBA_GRNDFTHR_IND':np.int8, 
        'IBF_IND':np.int8, 
        'ID_ABA_PRIM':np.int32, 
        'ID_CUSIP':object, 
        'ID_FDIC_CERT':np.int32, 
        'ID_LEI':object, 
        'ID_NCUA':np.int32, 
        'ID_OCC':np.int32, 
        'ID_RSSD':np.int32, 
        'ID_RSSD_HD_OFF':np.int32, 
        'ID_TAX':np.int32, 
        'ID_THRIFT':np.int32, 
        'ID_THRIFT_HC':object, 
        'INSUR_PRI_CD':np.int8, 
        'MBR_FHLBS_IND':bool,        # Boolean
        'MBR_FRS_IND':bool,          # Boolean
        'MJR_OWN_MNRTY':np.int8, 
        'NM_LGL':object, 
        'NM_SHORT':object, 
        'NM_SRCH_CD':np.int32, 
        'ORG_TYPE_CD':np.int8, 
        'PLACE_CD':np.int32, 
        'PRIM_FED_REG':object, 
        'PROV_REGION':object, 
        'REASON_TERM_CD':np.int8, 
        'SEC_RPTG_STATUS':np.int8, 
        'SLHC_IND':bool,             # Boolean
        'SLHC_TYPE_IND':np.int8,
        'STATE_ABBR_NM':object, 
        'STATE_CD':np.int8, 
        'STATE_HOME_CD':np.int8, 
        'STATE_INC_ABBR_NM':object, 
        'STATE_INC_CD':np.int8, 
        'STREET_LINE1':object, 
        'STREET_LINE2':object, 
        'URL':object, 
        'ZIP_CD':object
    }
    ATTdf = pd.read_csv(csvfile, dtype=DTYPES_ATT, sep=sep)
    # Create new column to serve as the primary key
    ATTdf['rssd'] = ATTdf['ID_RSSD']
    if (filter_asof):
        ATTdf = ATTdf[ATTdf.DT_END >= asofdate]
        ATTdf = ATTdf[ATTdf.DT_OPEN <= asofdate]
    ATTdf.insert(len(ATTdf.columns), 'NICsource', nicsource, allow_duplicates=True)
    ATTdf.reset_index(inplace=True)
    ATTdf.set_index(['rssd'], append=True, inplace=True)
    ATTdf.sort_index(inplace=True)
    return ATTdf



def RELcsv2df(csvfile, asofdate, sep=',', filter_asof=True):
    """Reads a NIC relationships CSV file into a Pandas dataframe
    
    NIC downloads typically start in XML format, although more recent
    (since 2018) vintages are also available as CSV downloads. 
    You must convert an XML download to tab-delimited CSV before using 
    this function. Upon reading, the contents of the CSV file are converted 
    to appropriate primitive types and stored in a Pandas dataframe, 
    which RELcsv2df returns. 
    The returned dataframe is indexed (and sorted) on four new fields:
      * rssd_par is a copy of the ID_RSSD_PARENT field
      * rssd_off is a copy of the ID_RSSD_OFFSPRING field
      * dt0 is a copy of the DT_START field
      * dt1 is a copy of the DT_END field

    Parameters
    ----------    
    csvfile : str
        Filename of a tab- or comma-delimited CSV file containing 
        a NIC relationships table
    asofdate : int
        The as-of date to include; an integer value of the form YYYYMMDD
    sep : str
        The field delimiter in the CSV file
    filter_asof : bool
        Whether to filter records based on the as-of date; if True, only
        relationships that were 'live' (asofdate between DT_START and 
        DT_END) are included
        
    Returns
    -------
    RELdf : Pandas datframe
        Dataframe of the CSV contents, perhaps filterd by as-of date
    """
    DTYPES_REL = {
        'CTRL_IND':np.int8, 
        'DT_RELN_EST':object, 
        'DT_START':np.int32, 
        'DT_END':np.int32, 
        'D_DT_RELN_EST':object, 
        'D_DT_START':object, 
        'D_DT_END':object, 
        'EQUITY_IND':np.int8, 
        'FC_IND':np.int8, 
        'ID_RSSD_OFFSPRING':np.int32, 
        'ID_RSSD_PARENT':np.int32, 
        'MB_COST':np.float64, 
        'OTHER_BASIS_IND':np.int8, 
        'PCT_EQUITY':np.float64, 
        'PCT_EQUITY_BRACKET':object, 
        'PCT_EQUITY_FORMAT':object, 
        'PCT_OTHER':np.float64, 
        'REASON_ROW_CRTD':np.int8, 
        'REASON_TERM_RELN':np.int8, 
        'REGK_INV':np.int8, 
        'REG_IND':np.int8, 
        'RELN_LVL':np.int8
    }
    RELdf = pd.read_csv(csvfile, dtype=DTYPES_REL, sep=sep)
    # Create new columns to serve as the (compound) primary key
    RELdf['rssd_par'] = RELdf['ID_RSSD_PARENT']
    RELdf['rssd_off'] = RELdf['ID_RSSD_PARENT']
    RELdf['dt0'] = RELdf['DT_START']
    RELdf['dt1'] = RELdf['DT_END']
    if (filter_asof):
        RELdf = RELdf[RELdf.DT_START <= asofdate]
        RELdf = RELdf[RELdf.DT_END >= asofdate]
    RELdf.reset_index(inplace=True)
    RELdf.set_index(['rssd_par', 'rssd_off', 'dt0', 'dt1'], append=True, inplace=True)
    RELdf.sort_index(inplace=True)
    return RELdf
    


# =============================================================================
#   FDIC Failures data
# =============================================================================
def create_FDICFail(indir, filename, asofdate):
    LOG.info(f'Processing for as-of date: {asofdate}')
    data = None
    filepath = os.path.join(indir, filename)
    LOG.debug(f'Path to FDIC Failure data input: {filepath}')
    data = FDICFailcsv2df(filepath, asofdate)
    LOG.debug(f'FDIC Failure data has {len(data)} obs')
    return data
        
def FDICFailcsv2df(csvfile):
    DTYPES_FAIL = {
        'CERT':np.int32, 
        'CHCLASS1':object, 
        'CITYST':object, 
        'COST':object, 
        'FAILDATE':object, 
        'FIN':np.int32, 
        'NAME':object, 
        'QBFASSET':np.int32, 
        'QBFDEP':np.int32, 
        'RESTYPE':object, 
        'RESTYPE1':object, 
        'SAVR':object, 
    }
    faildate_parser = lambda x: pd.datetime.strptime(x, '%m/%d/%y')
    nans = {'COST': [''], 'SAVR': ['***']}
    FAILdf = pd.read_csv(csvfile, dtype=DTYPES_FAIL, sep=',', 
      parse_dates=['FAILDATE'], date_parser=faildate_parser, na_values=nans)
    FAILdf['COST'] = pd.to_numeric(FAILdf.COST)
    FAILdf['cert'] = pd.to_numeric(FAILdf.CERT)
    FAILdf.reset_index(inplace=True)
    FAILdf.set_index(['cert'], inplace=True)
    return FAILdf



# =============================================================================
#   FDIC Summary of Deposits data
# =============================================================================
def create_FDICSoD(indir, filelist, asofdate):
    LOG.info(f'Processing for as-of date: {asofdate}')
    data = None
    # Find the correct input file to use, by parsing the filenames
    asof_yyyy = int(asofdate/10000)
    asof_mmdd = asofdate - asof_yyyy*10000
    maxyyyy = -1
    minyyyy = 9999
    max_file = ''
    min_file = ''
    if (asof_mmdd < 630):
        asof_yyyy = asof_yyyy - 1
    filename = None
    for fn in filelist:
        LOG.debug(f'Checking file: {fn}')
        DD = re.sub('[A-Za-z_]*([0-9]+)\\.csv', '\\1', fn)
        print(fn, DD)
        file_yyyy = int(DD)
        if (file_yyyy < minyyyy):
            minyyyy = min(minyyyy,file_yyyy)
            min_file = fn
        if (file_yyyy > maxyyyy):
            maxyyyy = max(maxyyyy,file_yyyy)
            max_file = fn
        maxyyyy = max(maxyyyy,file_yyyy)
        LOG.debug(f'Checking DD: {DD}, for {file_yyyy}')
        if (asof_yyyy==file_yyyy):
            filename = fn
            break
    if (filename is None):
        LOG.warning(f'Bad asofdate: {asofdate}; out of file range')
        if (asof_yyyy >= maxyyyy):
            filename = max_file
        elif (asof_yyyy <= minyyyy):
            filename = min_file
        else:
            LOG.error(f'No FDIC SoD file for: {asofdate}')
        LOG.info(f'Matched file: {filename} for asofdate={asofdate}')
    filepath = os.path.join(indir, filename)
    LOG.debug(f'Path to FDIC SoD input: {filepath}')
    data = FDICSoDcsv2df(filepath)
    LOG.debug(f'FDIC SoD has {len(data)} obs')
    return data

def FDICSoDcsv2df(csvfile):
    DTYPES_FDICSoD = {
        'YEAR':np.int16,
        'CERT':np.int32,
        'BRNUM':np.int32,
        'UNINUMBR':np.int32,
        'NAMEFULL':object,
        'ADDRESBR':object,
        'CITYBR':object,
        'CNTYNAMB':object,
        'STALPBR':object,
        'ZIPBR':object,
        'BRCENM':object,
        'CONSOLD':object,
        'BRSERTYP':np.int8,
        'DEPSUMBR':np.float64,
        'BKMO':np.int8,
        'CBSA_DIV_NAMB':object,
        'CITY2BR':object,
        'CNTRYNAB':object,
        'CNTYNUMB':np.int16,
        'CSABR':np.int16,
        'CSANAMBR':object,
        'DIVISIONB':np.int32,
        'MSABR':np.int32,
        'MSANAMB':object,
        'METROBR':np.int8,
        'MICROBR':np.int8,
        'NAMEBR':object,
        'NECTABR':np.int8,
        'NECNAMB':object,
        'PLACENUM':np.int32,
        'SIMS_ACQUIRED_DATE':object,
        'SIMS_ESTABLISHED_DATE':object,
        'SIMS_LATITUDE':np.float32,
        'SIMS_LONGITUDE':np.float32,
        'SIMS_DESCRIPTION':object,
        'SIMS_PROJECTION':object,
        'STCNTYBR':np.float32,
        'STNAMEBR':object,
        'STNUMBR':np.int8,
        'HCTMULT':object,
        'RSSDHCR':np.int32,
        'NAMEHCR':object,
        'CITYHCR':object,
        'STALPHCR':object,
        'RSSDID':np.int32,
        'UNIT':object,
        'ADDRESS':object,
        'CITY':object,
        'STALP':object,
        'ZIP':object,
        'ASSET':np.float64,
        'BKCLASS':object,
        'CALL':object,
        'CHARTER':object,
        'CHRTAGNN':object,
        'CHRTAGNT':object,
        'CLCODE':np.int8,
        'CNTRYNA':object,
        'DENOVO':object,
        'DEPDOM':np.float64,
        'DEPSUM':np.float64,
        'DOCKET':np.int32,
        'ESCROW':np.float64,
        'FDICDBS':np.int8,
        'FDICNAME':object,
        'FED':np.int8,
        'FEDNAME':object,
        'INSAGNT1':object,
        'INSURED':object,
        'INSBRDD':np.float64,
        'INSBRTS':np.float64,
        'OCCDIST':np.int8,
        'OCCNAME':object,
        'REGAGNT':object,
        'SPECGRP':np.int8,
        'SPECDESC':object,
        'STCNTY':np.int32,
        'STNAME':object,
        'USA':np.int8,
    }
    # By default, Pandas parses 'NA' as missing so no special treatment 
    # is needed for the FDIC CB files (which also use 'NA')
    FDICSoDdf = pd.read_csv(csvfile, dtype=DTYPES_FDICSoD, encoding='latin_1', sep=',', thousands=r',')
    # Create new columns to serve as the primary key
    FDICSoDdf['cert'] = FDICSoDdf['CERT']
    FDICSoDdf['brnum'] = FDICSoDdf['BRNUM']
    FDICSoDdf.reset_index(inplace=True)
    FDICSoDdf.set_index(['cert','brnum'], append=True, inplace=True)
    FDICSoDdf.sort_index(inplace=True)
    return FDICSoDdf
    


# =============================================================================
#   FDIC Community Banking data
# =============================================================================
def create_FDICCB(indir, filelist, asofdate):
    LOG.info(f'Processing for as-of date: {asofdate}')
    data = None
    # Find the correct input file to use, by parsing the filenames
    asof_recent = UTIL.rcnt_qtrend(asofdate)
    filename = None
    for fn in filelist:
        LOG.info('Checking file: '+fn)
        DD = re.sub('[A-Za-z_]*([0-9]+)\\-([0-9]+)\\.csv', '\\1_\\2', fn)
        LOG.info('Checking DD: '+DD)
        yyyymmdd0 = int(DD[0:4]+'0331')
        yyyymmdd1 = int(DD[5:9]+'1231')
        if (asof_recent>=yyyymmdd0 and asof_recent<=yyyymmdd1):
            filename = fn
            break
    if (filename is None):
        LOG.error('No FDIC CB file for: '+str(asofdate))
    filepath = os.path.join(indir, filename)
    LOG.debug(f'Path to FDIC CB input: {filepath}')
    data = FDICCBcsv2df(filepath, asofdate)
    LOG.debug(f'FDIC CB has {len(data)} obs')
    return data


    
def FDICCBcsv2df(csvfile, asofdate, filter_asof=True):
    DTYPES_FDICCB = {
        'Id':np.int32,
        'NameHCR':object,
        'RSSDHCR':np.int32,
        'Excl_Specialty':np.int8,
        'Excl_Foreign':np.int8,
        'Excl_Nocore':np.int8,
        'Excl_Noloans':np.int8,
        'Excl_Credit':np.int8,
        'CERT':np.int32,
        'NameFull':object,
        'Address':object,
        'City':object,
        'Stalp':object,
        'Zip':object,
        'Stcnty':np.int32,
        'Year':np.int16,
        'Callym':np.int32,
        'CB':np.int8,
        'Asset':np.float64,
        'HCAsset':np.float64,
        'Size':object,
        'ForAsset':object,
        'LoanToAsset':object,
        'CoreRatio':object,
        'Max_Deposits':np.float64,
        'Office_Count':np.float32,
        'Unique_Metros':np.float32,
        'State_Count':np.float32,
        'Business_Line':object,
#        'Meets_LTA':np.float32,
#        'Meets_Core':np.float32,
#        'Meets_Deposits':np.float32,
#        'Meets_Offices':np.float32,
#        'Meets_LrgMSA':np.float32,
#        'Meets_States':np.float32,
        'Run_Date':np.int32,
    }
    # By default, Pandas parses 'NA' as missing so no special treatment 
    # is needed for the FDIC CB files (which also use 'NA')
    FDICCBdf = pd.read_csv(csvfile, dtype=DTYPES_FDICCB, sep=',', thousands=r',')
    # Create new columns to serve as the (compound) primary key
    FDICCBdf['cert'] = FDICCBdf['CERT']
    FDICCBdf['calldate'] = FDICCBdf['Callym']
    if (filter_asof):
        asofyyyymm = int(UTIL.rcnt_qtrend(asofdate)/100)
        FDICCBdf = FDICCBdf[FDICCBdf.calldate == asofyyyymm]
    FDICCBdf.reset_index(inplace=True)
    FDICCBdf.set_index(['cert', 'calldate'], append=True, inplace=True)
    FDICCBdf.sort_index(inplace=True)
    return FDICCBdf
    
        

        
#def fetch_NIC(outdir, asofdate, indir=None, fA=None, fB=None, fC=None, fREL=None):
#    NIC = None
#    datafilename = 'DATA_'+str(asofdate)+'.pik'
#    datafilepath = os.path.join(outdir, datafilename)
#    nonefiles = (None==indir or None==fA or None==fB or None==fC or None==fREL)
#    if (os.path.isfile(datafilepath)):
#        f = open(datafilepath, 'rb')
#        NIC = pik.load(f)
#        f.close()
#    elif (not(nonefiles)):
#        NIC = makeDATA(indir, fA, fB, fC, fREL, asofdate)
#        f = open(datafilepath, 'wb')
#        pik.dump(NIC, f)
#        f.close()
#    return NIC


    
    
#def fetch(cachedir, filename):
#    data = None
#    datafilepath = os.path.join(cachedir, filename)
#    if (os.path.isfile(datafilepath)):
#        f = open(datafilepath, 'rb')
#        data = pik.load(f)
#        f.close()
#    elif (spec is not None):
#        case = spec['type']
#        asofdate = spec['asofdate']
#        if ('NIC'==case):
#            fA = spec['file_ATTRIBUTES_ACTIVE']
#            fB = spec['file_ATTRIBUTES_BRANCH']
#            fC = spec['file_ATTRIBUTES_CLOSED']
#            fREL = spec['file_RELATIONSHIPS']
#            data = makeDATA(cachedir, fA, fB, fC, fREL, asofdate)
#        elif ('FDICCB'==case):
#            pass
#        elif ('FDICSoD'==case):
#            pass
#        elif ('FDICFail'==case):
#            pass
#        f = open(datafilepath, 'wb')
#        pik.dump(data, f)
#        f.close()
#    return data



## A convenience function to look up and return the column number for the 
## four columns composing the index in the relationships dataframe. 
## See the function RELcsv2df for further details. 
#def REL_IDcols(RELdf):
#    # Get the column numbers to dereference the values packed in the multiindex
#    ID_RSSD_PARENT = RELdf.index.names.index('ID_RSSD_PARENT')
#    ID_RSSD_OFFSPRING = RELdf.index.names.index('ID_RSSD_OFFSPRING')
#    DT_START = RELdf.index.names.index('DT_START')
#    DT_END = RELdf.index.names.index('DT_END')
#    return ID_RSSD_PARENT, ID_RSSD_OFFSPRING, DT_START, DT_END



#def maps_rssd_cert(DATA):
#    rssd2cert = dict()
#    cert2rssd = dict()
#    ATTdf = DATA[IDX_Attributes]
#    ATTdf = ATTdf[ATTdf.ID_FDIC_CERT > 0]
#    for idx,row in ATTdf.iterrows():
#        rssd = idx
#        cert = row['ID_FDIC_CERT']
#        rssd2cert[rssd] = cert
#        cert2rssd[cert] = rssd
#    return (rssd2cert, cert2rssd)
#
#
#def augment_FAILdf(FAILdf, outdir, dataasof):
#    FAILdf.sort_values(by=['FAILDATE'], inplace=True)
#    DATA = fetch_NIC(outdir, dataasof)
#    ATTdf = DATA[IDX_Attributes]
#    ATTdf = ATTdf[ATTdf.ID_FDIC_CERT > 0]
#    (rssd2cert, cert2rssd) = maps_rssd_cert(DATA)
#    FAILdf2 = FAILdf.copy(deep=True)
#    FAILdf2['RSSD']=-1
#    FAILdf2['RSSD_HH']=-1
#    FAILdf2['ENTITY_TYPE']=''
#    FAILdf2['CNTRY_NM']=''
#    FAILdf2['STATE_ABBR_NM']=''
#    rcntasof = -1
#    banksys = None
#    for idx,row in FAILdf.iterrows():
#        failasof = FAILdf.loc[idx]['FAILDATE']
#        failasof = failasof.year*10000 + failasof.month*100 + failasof.day
#        if (rcntasof != UTIL.rcnt_qtrend(failasof)):
#            rcntasof = UTIL.rcnt_qtrend(failasof)
#            sysfilename = 'NIC_'+'_'+str(rcntasof)+'.pik'
#            sysfilepath = os.path.join(outdir, sysfilename)
#            f = open(sysfilepath, 'rb')
#            banksys = pik.load(f)
#            f.close()
#        cert = FAILdf.loc[idx]['CERT']
#        rssd = cert2rssd[cert]
#        NICdict = ATTdf.loc[rssd].to_dict()
#        FAILdf2.loc[cert,('RSSD')] = rssd
#        FAILdf2.loc[cert,('RSSD_HH')] = rssd
##        FAILdf2.loc[cert,('DT_OPEN')] = NICdict['DT_OPEN']
##        FAILdf2.loc[cert,('DT_START')] = NICdict['DT_START']
##        FAILdf2.loc[cert,('DT_END')] = NICdict['DT_END']
#        FAILdf2.loc[cert,('ENTITY_TYPE')] = NICdict['ENTITY_TYPE']
##        FAILdf2.loc[cert,('CNTRY_CD')] = NICdict['CNTRY_CD']
#        FAILdf2.loc[cert,('COUNTRY')] = NICdict['CNTRY_NM'].strip()
#        FAILdf2.loc[cert,('STATE')] = NICdict['STATE_ABBR_NM']
#    return FAILdf2



#def fetch_banksys(sysfilepath, csvfilepath, asofdate):
##    DATA = None
##    datafilename = 'DATA_'+str(asofdate)+'.pik'
##    datafilepath = os.path.join(outdir, datafilename)
##    nonefiles = (None==indir or None==fA or None==fB or None==fC or None==fREL)
##    if (os.path.isfile(datafilepath)):
##        f = open(datafilepath, 'rb')
##        DATA = pik.load(f)
##        f.close()
##    elif (not(nonefiles)):
##        DATA = makeDATA(indir, fA, fB, fC, fREL, asofdate)
##        f = open(datafilepath, 'wb')
##        pik.dump(DATA, f)
##        f.close()
##    return DATA
#
#    BankSys = None
#    if os.path.isfile(sysfilepath):
##        if (veryverbose): print('FOUND: Banking system file path:   ', sysfilepath)
#        f = open(sysfilepath, 'rb')
#        BankSys = pik.load(f)
#        f.close()
#    else:
##        if (veryverbose): print('CREATING: Banking system file path:', sysfilepath, asofdate)
#        BankSys = nx.DiGraph()
##        if (veryverbose): print('CSV file path:', csvfilepath, asofdate)
#        RELdf = UTIL.RELcsv2df(csvfilepath, asofdate)
#        (ID_RSSD_PARENT, ID_RSSD_OFFSPRING, DT_START, DT_END) = UTIL.REL_IDcols(RELdf)
#        for row in RELdf.iterrows():
#            date0 = int(row[0][DT_START])
#            date1 = int(row[0][DT_END])
#            rssd_par = row[0][ID_RSSD_PARENT]
#            rssd_off = row[0][ID_RSSD_OFFSPRING]
#            if (asofdate < date0 or asofdate > date1):
#                if (verbose): print('ASOFDATE,', asofdate, 'out of bounds:', rssd_par, rssd_off, date0, date1)
#                continue   
#            BankSys.add_edge(rssd_par, rssd_off)
#        f = open(sysfilepath, 'wb')
#        pik.dump(BankSys, f)
#        f.close()
##    if (veryverbose): print('System as of '+str(asofdate)+' has', BankSys.number_of_nodes(), 'nodes and', BankSys.number_of_edges(), 'edges')
#    return BankSys

# =============================================================================
#   Utility methods
# =============================================================================
def test_zipfile_integrity(zipfile):
    """Tests the integrity of a local zip file
    
    This function runs zf.ZipFile(zipfile) to test for a bad file. Then runs
    zipfile.testzip() to check again. In either case, this function logs an
    error if an invalid file is detected.
    
    Parameters
    ----------
    zipfile : path
        Full system path to the zip file to be tested
    """
    try:
        zfile = zf.ZipFile(zipfile)
    except zf.BadZipfile as bzf:
        LOG.error('Corrupt zip file: '+zipfile, ': '+str(bzf))
    test = zfile.testzip()
    if not(None==test):
        LOG.error('Corrupt zip file: '+zipfile)
    else:
        LOG.debug('Valid zip file: '+zipfile)

def sed(textfile, pattern, replace, N=0):
    """Replaces instances of a pattern in a text file

    Iterates through the lines of the file, replacing each instance of
    the regular expression pattern. Details on Python handling of 
    regular expressions appear in the 
    `regular expression docs <https://docs.python.org/3/library/re.html>`_.

    Parameters
    ----------
    pattern : str
        Regular expression to match
    replace : str
        Regular expression replacement string
    textfile : str
        Filename of input to be modified
    N : int
        Number of lines in the textfile to treat. If N < 1 (the default),
        then a global replace occurs.
    """
    cpattern = re.compile(pattern)
    count = 0
    global_replace = (N < 1)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        with open(textfile) as src:
            for line in src:
                line_revised = cpattern.sub(replace, line)
                tmp.write(line_revised)
                if (line_revised != line):
                    count = count + 1
                if not(global_replace):
                    if (count > N):
                        break
            try:
                tmp.writelines(src.readlines())
                # Copy file attributes (permissions) to the temporary file
                shutil.copystat(textfile, tmp.name)
                shutil.move(tmp.name, textfile)
            except (Exception) as e:
                LOG.error('File update failed, count='+str(count)+' '+str(e))
                raise e

        



if __name__ == "__main__":
    import doctest
    doctest.testmod()