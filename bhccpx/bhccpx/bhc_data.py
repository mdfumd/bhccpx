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
# Last revision: 5-Sep-2019
# -----------------------------------------------------------------------------

import os
import time
import re
import shutil
import tempfile

import graphviz as gv
import zipfile as zf
import networkx as nx
import numpy as np 
import pandas as pd
import _pickle as pik 

import bhc_util as UTIL

__all__ = ['resolve_dir_nic',
           'fetch_data_nic',
           'fetch_data_fdicsod',
           'fetch_data_fdiccb',
           'fetch_data_fdicfail',
           'fetch_data_banksys',
           ]
__version__ = '0.5'
__author__ = 'Mark D. Flood'


MODNAME = __file__.split(os.path.sep)[-1].split('.')[0]
LOG = UTIL.log_config(MODNAME)

# =============================================================================
#   Cache management
# =============================================================================
def fetch_data_banksys(config, asofdate, sect):
    spec = cache_spec(config, sect)
    data = None
    cachepath = build_cachepath(spec, 'BankSys', asofdate)
    if (spec['banksys_clearcache'] and os.path.isfile(cachepath)): 
        os.remove(cachepath)
    data = cache_retrieve(cachepath)
    if (data is None):
        LOG.debug(f'Cache file not found: {cachepath}; creating')
        # Callback to the banksys constructor method
        data = make_banksys(config, asofdate, sect)
        cache_create(cachepath, data)
    return data


def fetch_data_fdicfail(config, sect):
    spec = cache_spec(config, sect)
    data = None
    cachepath = os.path.join(spec['cachedir'], 'FDICFail.pik')
    if (spec['fdicfail_clearcache'] and os.path.isfile(cachepath)): 
        os.remove(cachepath)
    data = cache_retrieve(cachepath)
    if (data is None):
        LOG.debug(f'Cache file not found: {cachepath}; creating')
        indir = spec['fdicfail_dir']
        filelist = spec['fdicfail_file']
        data = create_FDICFail(indir, filelist)
        cache_create(cachepath, data)
    return data


def fetch_data_fdiccb(config, asofdate, sect):
    spec = cache_spec(config, sect)
    data = None
    cachepath = build_cachepath(spec, 'FDICCB', asofdate)
    if (spec['fdiccb_clearcache'] and os.path.isfile(cachepath)): 
        os.remove(cachepath)
    data = cache_retrieve(cachepath)
    if (data is None):
        LOG.debug(f'Cache file not found: {cachepath}; creating')
        indir = spec['fdiccb_dir']
        filelist = spec['fdiccb_filelist']
        data = create_FDICCB(indir, filelist, asofdate)
        cache_create(cachepath, data)
    return data


def fetch_data_fdicsod(config, asofdate, sect):
    spec = cache_spec(config, sect)
    data = None
    cachepath = build_cachepath(spec,'FDICSoD', asofdate, quarterly_data=False)
    if (spec['fdicsod_clearcache'] and os.path.isfile(cachepath)): 
        os.remove(cachepath)
    data = cache_retrieve(cachepath)
    if (data is None):
        LOG.debug(f'Cache file not found: {cachepath}; creating')
        indir = spec['fdicsod_dir']
        filelist = spec['fdicsod_filelist']
        data = create_FDICSoD(indir, filelist, asofdate)
        cache_create(cachepath, data)
    return data


def fetch_data_nic(config, asofdate, sect):
    spec = cache_spec(config, sect)
    data = None
    cachepath = build_cachepath(spec, 'NIC', asofdate)
    if (spec['nic_clearcache'] and os.path.isfile(cachepath)): 
        LOG.error(f"SHOULD NOT BE HERE: {spec['nic_clearcache']}")
        os.remove(cachepath)
    data = cache_retrieve(cachepath)
    if (data is None):
        LOG.debug(f'Cache file not found: {cachepath}; creating')
        indir = spec['nic_dir']
        fA = spec['nic_attributesactive']
        fB = spec['nic_attributesbranch']
        fC = spec['nic_attributesclosed']
        fREL = spec['nic_relationships']
        sep = spec['nic_sep']
        filt = spec['nic_filterasof']
        data = create_NIC(indir, fA, fB, fC, fREL, asofdate, sep, filt)
        cache_create(cachepath, data)
    return data



def cache_create(cachepath, data):
    """ Creates or replace a cache file with the pickled data object
    """
    if (data is not None):
        LOG.debug(f'Creating cache file: {cachepath}')
        f = open(cachepath, 'wb')
        pik.dump(data, f)
        f.close()
    else:
        LOG.warning(f'Skipping cache creation, data==None: {cachepath}')



def cache_retrieve(cachepath):
    """ Retrieves a given file from the cache, if available
    """
    data = None
    if (os.path.isfile(cachepath)):
        LOG.debug(f'Retrieving file from: {cachepath}')
        f = open(cachepath, 'rb')
        data = pik.load(f)
        f.close()
    else:
        LOG.info(f'Requested cache file does not exist: {cachepath}')
    return data



def cache_spec(config, sectname):
    """Creates a cache specification dictionary from a config object.
    """
    spec = dict()
    sect = config[sectname]
    spec['cachedir'] = sect['cachedir']
    # NIC items
    spec['nic_dir'] = resolve_dir_nic(sect['nic_dir'], sect['nic_subdir'])
    spec['nic_attributesactive'] = eval(sect['nic_csv_filelist'])[0]
    spec['nic_attributesbranch'] = eval(sect['nic_csv_filelist'])[1]
    spec['nic_attributesclosed'] = eval(sect['nic_csv_filelist'])[2]
    spec['nic_relationships'] = eval(sect['nic_csv_filelist'])[3]
    spec['nic_transformations'] = eval(sect['nic_csv_filelist'])[4]
    spec['nic_sep'] = sect['delim']
    spec['nic_filterasof'] = True
    spec['nic_clearcache'] = ('TRUE'==sect['nic_clearcache'].upper())
    # FDIC CB items
    spec['fdiccb_dir'] = sect['fdiccb_dir']
    spec['fdiccb_filelist'] = eval(sect['fdiccb_csv_filelist'])
    spec['fdiccb_clearcache'] = ('TRUE'==sect['fdiccb_clearcache'].upper())
    # FDIC SoD items
    spec['fdicsod_dir'] = sect['fdicsod_dir']
    spec['fdicsod_filelist'] = eval(sect['fdicsod_csv_filelist'])
    spec['fdicsod_clearcache'] = ('TRUE'==sect['fdicsod_clearcache'].upper())
    # FDIC fail items
    spec['fdicfail_dir'] = sect['fdicfail_dir']
    spec['fdicfail_file'] = sect['fdicfail_csv_file']
    spec['fdicfail_clearcache'] = ('TRUE'==sect['fdicfail_clearcache'].upper())
    # BankSys items
    spec['banksys_clearcache'] = ('TRUE'==sect['banksys_clearcache'].upper())
    return spec



def build_cachepath(spec, case, asofdate, quarterly_data=True):
    cachedir = spec['cachedir']
    if (quarterly_data):
        cachefilename = f'{case}_{UTIL.rcnt_qtrend(asofdate)}.pik'
    else:
        # FDIC SoD is annual data, based at the mid-year point
        cachefilename = f'{case}_{UTIL.rcnt_midyear(asofdate)}.pik'
    cachepath = os.path.join(cachedir, cachefilename)
    return cachepath

    

# =============================================================================
#   NIC data
# =============================================================================
# Mnemonic indices for the NIC data list
IDX_Attributes = 0
IDX_Relationships = 1
#IDX _HighHolder = 2
#IDX _Entities = 3
#IDX _Parents = 4
#IDX _Offspring = 5
def create_NIC(indir, fA, fB, fC, fRel, asofdate, sep=',', filter_asof=True):
    """Parses the NIC data and creates a list of DataFrames for the asofdate

    Assembles the NIC data for given asofdate into a single list object.
    The returned list contains pointers to six objects, indexed as follows:
        
       0. IDX_Attributes - NIC attributes, combining active, branch, and closed
       1. IDX_Relationships - NIC relationships      
    """
#       2. IDX _HighHolder - List of the high holder BHCs in the system
#       3. IDX _Entities - Map (dict)
#       4. IDX _Parents       = 4
#       5. IDX _Offspring     = 5
    
    data = []
    csvfilepathA = os.path.join(indir, fA)
    csvfilepathB = os.path.join(indir, fB)
    csvfilepathC = os.path.join(indir, fC)
    ATTdf_a = ATTcsv2df(csvfilepathA, asofdate, 'A', sep, filter_asof)
    ATTdf_b = ATTcsv2df(csvfilepathB, asofdate, 'B', sep, filter_asof)
    ATTdf_c = ATTcsv2df(csvfilepathC, asofdate, 'C', sep, filter_asof)
    ATTdf = pd.concat([ATTdf_a, ATTdf_b, ATTdf_c])
    data.insert(IDX_Attributes, ATTdf)
    csvfilepathR = os.path.join(indir, fRel)
    RELdf = RELcsv2df(csvfilepathR, asofdate, sep)
    data.insert(IDX_Relationships, RELdf)
#    # Then, add derived structures based on the Relationships data:    
#    derived_data = NIC_highholders(RELdf, asofdate)
#    data.insert(IDX _HighHolder, derived_data[0])
#    data.insert(IDX _Entities, derived_data[1])
#    data.insert(IDX _Parents, derived_data[2])
#    data.insert(IDX _Offspring, derived_data[3])
    return data



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
    sep = UTIL.delim_norm(sep)
    LOG.info(f'Reading CSV file {csvfile} with delimiter: {sep}')
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
    which is returned.
    . 
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
    sep = UTIL.delim_norm(sep)
    LOG.info(f'Reading CSV file {csvfile} with delimiter: {sep}')
    RELdf = pd.read_csv(csvfile, dtype=DTYPES_REL, sep=sep)
    # Create new columns to serve as the (compound) primary key
    RELdf['rssd_par'] = RELdf['ID_RSSD_PARENT']
    RELdf['rssd_off'] = RELdf['ID_RSSD_OFFSPRING']
    RELdf['dt0'] = RELdf['DT_START']
    RELdf['dt1'] = RELdf['DT_END']
    if (filter_asof):
        RELdf = RELdf[RELdf.DT_START <= asofdate]
        RELdf = RELdf[RELdf.DT_END >= asofdate]
    RELdf.reset_index(inplace=True)
    RELdf.set_index(['rssd_par', 'rssd_off', 'dt0', 'dt1'], append=True, inplace=True)
    RELdf.sort_index(inplace=True)
    return RELdf
    


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



def make_banksys(config, asofdate, sect):
    """Reads or creates a graph of the full banking system on a given date
    
    The function looks for an existing graph in a pickle file
    located at sysfilepath (for example, .../cachedir/NIC__YYYYMMDD.pik),
    where YYYYMMDD is the asofdate. 
    If this file exists, it is unpackeded from the pickle and returned.
    If the file does not (yet) exist, the NetworkX DiGraph is instead created
    from the relationships data and dumped into a new pickle at sysfilepath.
    The graph is a naked directed graph whose nodes are NIC entities and 
    whose edges point from parent nodes to offspring nodes. 
    The function then returns this digraph (either newly created or unpickled). 
    """
    UTIL.tic()
    trace_logging = ('TRUE'==config[sect]['trace_logging'].upper())
    NICdata = fetch_data_nic(config, asofdate, sect)
    ATTdf = NICdata[IDX_Attributes]
    RELdf = NICdata[IDX_Relationships]
    BankSys = nx.DiGraph()
    LOG.debug(f'Relationships table has {len(RELdf)} obs')
    for row in RELdf.iterrows():
        date0 = int(row[1]['DT_START'])
        date1 = int(row[1]['DT_END'])
        rssd_par = row[1]['ID_RSSD_PARENT']
        rssd_off = row[1]['ID_RSSD_OFFSPRING']
        if (trace_logging):
            LOG.debug(f'Reln: asofdate={asofdate}, rssd_par={rssd_par}'+
                ', rssd_off={rssd_off}, date0={date0}, date1={date1}')
        if (asofdate < date0 or asofdate > date1):
            LOG.warning(f'As-of date: {asofdate} is out of bounds: '+
                f'rssd_par={rssd_par}, rssd_off={rssd_off}, '+
                f'date0={date0}, date1={date1}')
            continue   
        BankSys.add_edge(rssd_par, rssd_off)
    # Adding in the singleton institutions (no edges in Relationships file)
    LOG.info(f'System, as of {asofdate}, has {BankSys.number_of_nodes()} ' +
              f'nodes and {BankSys.number_of_edges()} edges')
    ATTdf = ATTdf[ATTdf.DT_END >= asofdate]
    ATTdf = ATTdf[ATTdf.DT_OPEN <= asofdate]
    nodes_BankSys = set(BankSys.nodes)
    nodes_ATTdf = set(ATTdf['ID_RSSD'].unique())
    nodes_new = nodes_ATTdf.difference(nodes_BankSys)
    BankSys.add_nodes_from(nodes_new)
    LOG.info(f'PROFILE building BankSys {asofdate} took {UTIL.toc()} secs.')
    LOG.debug(f'System, as of {asofdate}, has {BankSys.number_of_nodes()} ' +
              f'nodes and {BankSys.number_of_edges()} edges; '+
              f'{len(nodes_new)} added, out of {len(nodes_ATTdf)} candidates')
    return BankSys
    


# =============================================================================
#   FDIC Failures data
# =============================================================================
def create_FDICFail(indir, filename):
    LOG.info(f'Processing')
    data = None
    filepath = os.path.join(indir, filename)
    LOG.debug(f'Path to FDIC Failure data input: {filepath}')
    data = FDICFailcsv2df(filepath)
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
    LOG.info(f'Reading CSV file {csvfile}')
    FAILdf = pd.read_csv(csvfile, dtype=DTYPES_FAIL, sep=',', 
      parse_dates=['FAILDATE'], date_parser=faildate_parser, na_values=nans)
    FAILdf['COST'] = pd.to_numeric(FAILdf.COST)
    FAILdf['cert'] = pd.to_numeric(FAILdf.CERT)
    FAILdf.reset_index(inplace=True)
    FAILdf.set_index(['cert'], inplace=True)
    # Some derived date fields:
    FAILdf['asofdate'] = \
      pd.DatetimeIndex(FAILdf.FAILDATE).year*10000 + \
      pd.DatetimeIndex(FAILdf.FAILDATE).month*100 + \
      pd.DatetimeIndex(FAILdf.FAILDATE).day
    FAILdf['asofqtr'] = \
      FAILdf.apply(lambda row: UTIL.rcnt_qtrend(row.asofdate), axis=1)   
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
#        print(fn, DD)
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
    LOG.info(f'Path to FDIC SoD input: {filepath}')
    data = FDICSoDcsv2df(filepath)
    LOG.info(f'FDIC SoD has {len(data)} obs')
    return data



def FDICSoDcsv2df(csvfile):
    DTYPES_FDICSoD = {
        'YEAR':np.int16,
        'CERT':np.int32,
        'BRNUM':np.int32,
        'UNINUMBR':np.float64,
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
        'PLACENUM':object,
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
        'SPECGRP':np.float64,
        'SPECDESC':object,
        'STCNTY':np.int32,
        'STNAME':object,
        'USA':np.int8,
    }
#    FDICSoD_missings = {
#        'UNINUMBR': ['', ], 
#    }    
    
    # By default, Pandas parses 'NA' as missing so no special treatment 
    # is needed for the FDIC CB files (which also use 'NA')
    LOG.info(f'Reading CSV file {csvfile}')
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
    LOG.info(f'Reading CSV file {csvfile}')
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

        
def log_bhc_svg(BHC, outdir, fileroot, colormap, dim='entity_type', 
                title=None, condense=False):
    """Create an SVG image file representing a BHC. 
    
    The file is stored in the outdir, with the filename: 
    RSSD_<rssd_hh>_<asofdate>.svg. If popup is set to True, then the 
    function will also launch a browser to display the file. 
    """
    warnings = []
    vis = gv.Digraph()
    if (title is not None):
#        dot.attr(label=r'\n'+title)
        vis.node('TITLE', "info", shape="doublecircle", fontsize="20", style="filled", fillcolor="orange", tooltip=title)
    vis.attr('node', fontsize='8')
    vis.attr('node', fixedsize='true')
    if condense:
        # Reasonable svg viewbox=1700x700pt; graph scale=.03,.40; graph translate=4,250
        vis.attr('graph', size='9, 11')
        vis.attr('graph', scale='2, 5')
        vis.attr('node', shape='circle')
        vis.attr('node', width='0.04')
    else:
        vis.attr('node', shape='ellipse')
        vis.attr('node', width='0.7')
        vis.attr('node', height='0.3')
    for N in BHC.nodes():
        NM_LGL = ''
        ENTITY_TYPE = 'ZZZ'
        GEO_JURISD = 'ZZZ'
        attribute_error = True
        try:
            NM_LGL = BHC.node[N]['nm_lgl'].strip()
            ENTITY_TYPE = BHC.node[N]['entity_type']
            GEO_JURISD = BHC.node[N]['GEO_JURISD']
            attribute_error = False
        except KeyError as KE:
            warnings.append(f'Invalid attribute data for RSSD={N} for {fileroot}')
        color_key = ENTITY_TYPE
        if (dim=='GEO_JURISD'):
            color_key = GEO_JURISD.replace(' - 0', '').strip()
        tt = f'[{N}] {ENTITY_TYPE}\\n------------\\n{NM_LGL}\\n------------\\n{GEO_JURISD}'
        if condense:
            lbl = ''
        else:
            lbl = str(N)
        if (attribute_error):
            vis.node('rssd'+str(N), lbl, style="filled", fillcolor="red;.5:green", tooltip=tt)
        else:
            fc = colormap[color_key]
            vis.node('rssd'+str(N), lbl, style="filled", fillcolor=fc, tooltip=tt)
    for E in BHC.edges():
        src = 'rssd' + str(E[0])
        tgt = 'rssd' + str(E[1])
        Vs = BHC.node[E[0]]
        Vt = BHC.node[E[1]]
        col = 'red'
        if (dim not in Vs) or (dim not in Vt):
            col='green'
        elif Vs[dim]==Vt[dim]:
            col='black'
        vis.edge(src, tgt, arrowsize='0.3', color=col)
    svg_path = os.path.join(outdir, f'{fileroot}')
    vis.render(filename=svg_path, format='svg')

def log_bhcq_svg(BHC, outdir, fileroot, colormap, dim='entity_type', title=None):
    """Create an SVG image file representing a BHC quotient graph. 
    
    The file is stored in the outdir, with the filename: 
    RSSD_<rssd_hh>_<asofdate>.svg. If popup is set to True, then the 
    function will also launch a browser to display the file. 
    """
    warnings = []
    vis = gv.Graph()
    if (title is not None):
#        dot.attr(label=r'\n'+title)
        vis.node('TITLE', "info", shape="doublecircle", fontsize="20", style="filled", fillcolor="orange", tooltip=title)
    vis.attr('node', fontsize='8')
    vis.attr('node', fixedsize='true')
#    vis.attr('node', height='0.3')
    vis.attr('node', shape='circle')
    vis.attr('node', width='0.5')
    for N in BHC.nodes():
        NM_LGL = ''
        ENTITY_TYPE = 'ZZZ'
        GEO_JURISD = 'ZZZ'
        attribute_error = True
        try:
            NM_LGL = BHC.node[N]['nm_lgl'].strip()
            ENTITY_TYPE = BHC.node[N]['entity_type']
            GEO_JURISD = BHC.node[N]['GEO_JURISD']
            attribute_error = False
        except KeyError as KE:
            warnings.append(f'Invalid attribute data for RSSD={N} for {fileroot}')
        color_key = ENTITY_TYPE
        if (dim=='GEO_JURISD'):
            color_key = GEO_JURISD.replace(' - 0', '').strip()
        tt = f'[{N}] {ENTITY_TYPE}\\n------------\\n{NM_LGL}\\n------------\\n{GEO_JURISD}'
        if (attribute_error):
#            vis.node('rssd'+str(N), str(N), style="filled", fillcolor="red;.5:green", tooltip=tt)
#            vis.node('rssd'+str(N), str(N), style="filled", fillcolor="red;.5:green", tooltip=tt)
            vis.node('rssd'+str(N), str(N), style="filled", fillcolor="beige", tooltip=tt)
        else:
            fc = colormap[color_key]
            fc = 'beige'
#            vis.node('rssd'+str(N), str(N), style="filled", fillcolor=fc, tooltip=tt)
            vis.node('rssd'+str(N), str(N), style="filled", fillcolor=fc, tooltip=tt)
    for E in BHC.edges():
        src = 'rssd' + str(E[0])
        tgt = 'rssd' + str(E[1])
        Vs = BHC.node[E[0]]
        Vt = BHC.node[E[1]]
        col = 'red'
        if (dim not in Vs) or (dim not in Vt):
            col='green'
            col='black'
        elif Vs[dim]==Vt[dim]:
            col='black'
        vis.edge(src, tgt, arrowsize='0.3', color=col)
    svg_path = os.path.join(outdir, f'{fileroot}')
    vis.render(filename=svg_path, format='svg')


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    
    