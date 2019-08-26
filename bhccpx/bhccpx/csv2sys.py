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

import os
import logging

import networkx as nx
import _pickle as pik
import multiprocessing as mp
import progressbar as pb

import bhc_data as DATA
import bhc_util as UTIL
import bhca

LOG = UTIL.log_config(__file__.split(os.path.sep)[-1].split('.')[0])


MAP_CERT2HCR = 0
MAP_HCR2CERT = 1
MAP_CERT2RSSD = 2
MAP_RSSD2CERT = 3
MAP_RSSD2HCR = 4
MAP_HCR2RSSD = 5
MAP_RSSD2HH = 6
MAP_HH2RSSD = 7


def cache_spec(config, asofdate):
    spec = dict()
    sect = config['csv2sys']
    spec['asofdate'] = asofdate
    spec['cachedir'] = sect['cachedir']
    # NIC items
    spec['nic_dir'] = DATA.resolve_dir_nic(sect['nic_dir'], sect['nic_subdir'])
    spec['nic_attributesactive'] = sect['attributesactive']
    spec['nic_attributesbranch'] = sect['attributesbranch']
    spec['nic_attributesclosed'] = sect['attributesclosed']
    spec['nic_relationships'] = sect['relationships']
    spec['nic_sep'] = sect['delim']
    spec['nic_filterasof'] = True
    spec['nic_clearcache'] = ('TRUE'==sect['clearcache_nic'].upper())
    # FDIC CB items
    spec['fdiccb_dir'] = sect['fdiccb_dir']
    spec['fdiccb_filelist'] = eval(sect['fdiccb_csvnrm_filelist'])
    spec['fdiccb_clearcache'] = ('TRUE'==sect['clearcache_fdiccb'].upper())
    # FDIC SoD items
    spec['fdicsod_dir'] = sect['fdicsod_dir']
    spec['fdicsod_filelist'] = eval(sect['fdicsod_csv_filelist'])
    spec['fdicsod_clearcache'] = ('TRUE'==sect['clearcache_fdicsod'].upper())
    # BankSys items
    spec['banksys_clearcache'] = ('TRUE'==sect['clearcache_banksys'].upper())
    # IDmaps items
    spec['idmaps_clearcache'] = ('TRUE'==sect['clearcache_idmaps'].upper())
    return spec

def make_banksys(config, asofdate, read_data=False):
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
    sysfilename = 'NIC_'+str(asofdate)+'.pik'
    sysfilepath = os.path.join(config['csv2sys']['cachedir'], sysfilename)
    trace_logging = ('TRUE'==config['csv2sys']['trace_logging'].upper())
    spec = cache_spec(config, asofdate)
    BankSys = DATA.fetch_data(spec, DATA.case_BankSys)
    LOG.info(f'Processing for asofdate={asofdate}')
#    if os.path.isfile(sysfilepath):
#        LOG.info('FOUND: Banking system file: '+sysfilepath+
#                 ' for as-of date: '+str(asofdate))
#        if (read_data):
#            f = open(sysfilepath, 'rb')
#            BankSys = pik.load(f)
#            f.close()
#        else:
#            LOG.debug('Not reading: '+sysfilepath+' (read_data==False)')
    if (BankSys is None):
        NICdata = DATA.fetch_data(spec, DATA.case_NIC)
        ATTdf = NICdata[DATA.IDX_Attributes]
        RELdf = NICdata[DATA.IDX_Relationships]
        BankSys = nx.DiGraph()
#        relfilename = config['csv2sys']['relationships']
#        nicdir = DATA.resolve_dir_nic(config['csv2sys']['nic_dir'], 
#                                      config['csv2sys']['nic_subdir'])
#        csvfilepath = os.path.join(nicdir, relfilename)
#        LOG.info('Reading relationshipx file: '+csvfilepath)
#        delim = UTIL.delim_norm(config['csv2sys']['delim'])
#        RELdf = DATA.RELcsv2df(csvfilepath, asofdate, sep=delim)
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
                LOG.warning('As-of date: '+str(asofdate)+' out of bounds: '+
                    'rssd_par='+str(rssd_par)+', rssd_off='+str(rssd_off)+
                    ', date0='+str(date0)+', date1='+str(date1))
                continue   
            BankSys.add_edge(rssd_par, rssd_off)
        # Adding in the singleton institutions (no edges in Relationships file)
        LOG.info('System (pre), asofdate='+str(asofdate)+' has '+
                 str(BankSys.number_of_nodes())+' nodes and '+
                 str(BankSys.number_of_edges())+' edges ')
#        fA = config['csv2sys']['attributesactive']
#        fB = config['csv2sys']['attributesbranch']
#        fC = config['csv2sys']['attributesclosed']
#        ATTdf = DATA.makeATTs(nicdir, fA, fB, fC, asofdate, 
#                              sep=delim, filter_asof=True)
        ATTdf = ATTdf[ATTdf.DT_END >= asofdate]
        ATTdf = ATTdf[ATTdf.DT_OPEN <= asofdate]
        nodes_BankSys = set(BankSys.nodes)
        nodes_ATTdf = set(ATTdf['ID_RSSD'].unique())
        nodes_new = nodes_ATTdf.difference(nodes_BankSys)
        BankSys.add_nodes_from(nodes_new)
        # Storing the new banking system file
        LOG.info(f'CACHING: System file: {sysfilepath}')
        DATA.cache_data(spec, DATA.case_BankSys, BankSys)
#        f = open(sysfilepath, 'wb')
#        pik.dump(BankSys, f)
#        f.close()
#        LOG.debug('CACHED: System file: '+sysfilepath)
        LOG.debug(f'System (pre), asofdate={asofdate} has {BankSys.number_of_nodes()} nodes and {BankSys.number_of_edges()} edges; {len(nodes_new)} added, out of {len(nodes_ATTdf)} candidates')
    return BankSys
    


def build_sys(config):
    """Builds requested network representations of the full banking system
    
    Converts the NIC data into NetworkX objects (DiGraphs), each representing
    the ownership relationships in the banking system at a given as-of date.
    The new network objects are stored in the local cache for faster 
    subsequent retrieval. However, if requested, this function will clear 
    existing versions of the requested network objects from the cache 
    before beginning. 
    
    Parameters
    ----------
    config : ConfigParser object
        The key parameters are:
          * config['csv2sys']['clearcache_nic']
          * config['csv2sys']['cachedir'] 
          * config['csv2sys']['asofdate0']
          * config['csv2sys']['asofdate1']
    """
    # If clearcache, then remove existing banking system pik files and recreate
#    if ('TRUE'==config['csv2sys']['clearcache_nic'].upper()):
#        LOG.warning('Clearing output cache of *.pik files in the range: '+
#            config['csv2sys']['asofdate0']+'-'+config['csv2sys']['asofdate1'])
#        DATA.clear_cache(config['csv2sys']['cachedir'], 'NIC_',
#            config['csv2sys']['asofdate0'], config['csv2sys']['asofdate1'])
    asof_list = UTIL.assemble_asofs(config['csv2sys']['asofdate0'], 
                                    config['csv2sys']['asofdate1'])
    LOG.info('As-of dates for NIC data: '+str(asof_list))
    if (int(config['csv2sys']['parallel']) > 0):
        LOG.warning('Begin parallel processing ('+str(len(asof_list))+
                    ' tasks across '+config['csv2sys']['parallel']+
                    ' cores) of system snapshots')
        pcount = min(int(config['csv2sys']['parallel']), 
                     os.cpu_count(), len(asof_list))
        pool = mp.Pool(processes=pcount)
        for asof in asof_list:
            pool.apply_async(make_banksys, (config, asof))
        pool.close()
        pool.join()
        LOG.info('Complete parallel processing of system snapshots')
    else:
        LOG.info('Begin sequential processing of system snapshots ('+ 
              str(len(asof_list))+' dates)')
        for asof in pb.progressbar(asof_list, redirect_stdout=True):
             make_banksys(config, asof)
        LOG.info('Complete sequential processing of system snapshots')




#def build_hhdata(config):
#    """Caches the FDIC CB data containing high-holder mappings
#    
#    Converts the FDIC CB data into Pandas dataframes, each representing
#    the high-holder relationships at a given as-of date.
#    The new dataframe objects are stored in the local cache for faster 
#    subsequent retrieval. However, if requested, this function will clear 
#    existing versions of the requested dataframes from the cache 
#    before beginning. 
#    
#    Parameters
#    ----------
#    config : ConfigParser object
#        The key parameters are:
#          * config['csv2sys']['clearcache_fdiccb']
#          * config['csv2sys']['cachedir'] 
#          * config['csv2sys']['asofdate0']
#          * config['csv2sys']['asofdate1']
#    """
#    
#    if ('TRUE'==config['csv2sys']['clearcache_fdiccb'].upper()):
#        # Remove cached pik files and recreate
#        LOG.debug('Clearing output cache of *.pik files in the range: '+
#            config['csv2sys']['asofdate0']+'-'+config['csv2sys']['asofdate1'])
#        DATA.clear_cache(config['csv2sys']['cachedir'], 'FDICCB_',
#            config['csv2sys']['asofdate0'], config['csv2sys']['asofdate1'])
#        
#    asof_list = UTIL.assemble_asofs(config['csv2sys']['asofdate0'], 
#                                    config['csv2sys']['asofdate1'])
#    LOG.info('As-of dates for FDIC CB: '+str(asof_list))
#    if (int(config['csv2sys']['parallel']) > 0):
#        # Parallel processing of each requested asofdate
#        LOG.info('Begin parallel processing ('+str(len(asof_list))+
#                    ' tasks across '+config['csv2sys']['parallel']+
#                    ' cores) of high-holders')
#        pcount = min(int(config['csv2sys']['parallel']), 
#                     os.cpu_count(), len(asof_list))
#        pool = mp.Pool(processes=pcount)
#        for asof in asof_list:
#            pool.apply_async(make_fdiccb_data, (config, asof))
#        pool.close()
#        pool.join()
#        LOG.info('Complete parallel processing of high-holders')
#    else:
#        # Sequential processing of each requested asofdate
#        LOG.info('Begin sequential processing of high-holders ('+ 
#              str(len(asof_list))+' dates)')
#        for asof in pb.progressbar(asof_list, redirect_stdout=True):
#             make_fdiccb_data(config, asof)
#        LOG.info('Complete sequential processing of high-holders')
    
#def maps_rssd_cert(data):
#    rssd2cert = dict()
#    cert2rssd = dict()
#    ATTdf = data[DATA.IDX_Attributes]
#    ATTdf = ATTdf[ATTdf.ID_FDIC_CERT > 0]
#    for idx,row in ATTdf.iterrows():
#        rssd = idx
#        cert = row['ID_FDIC_CERT']
#        rssd2cert[rssd] = cert
#        cert2rssd[cert] = rssd
#    return (rssd2cert, cert2rssd)
#
#

def get_idmaps(config, asofdate):
#    id_filename = 'IDmaps_'+str(UTIL.rcnt_qtrend(asofdate))+'.pik'
#    id_filepath = os.path.join(config['csv2sys']['cachedir'], id_filename)
#    LOG.info('Processing IDmaps, asofdate='+str(asofdate))
    spec = cache_spec(config, asofdate)
    IDmaps = DATA.fetch_data(spec, DATA.case_IDmaps)
#    if os.path.isfile(id_filepath):
#        LOG.info('Found ID maps file: '+id_filepath)
#        if (read_data):
#            f = open(id_filepath, 'rb')
#            IDmaps = pik.load(f)
#            f.close()
#        else:
#            LOG.debug('Not reading: '+id_filepath+' (read_data==False)')
    if (IDmaps is None):
        LOG.info('ID maps file not found: {id_filepath}, building instead')
        LOG.info(f'Building BankSys object for {asofdate}')
        BankSys = make_banksys(config, asofdate, read_data=True)
#        fA = config['csv2sys']['attributesactive']
#        fB = config['csv2sys']['attributesbranch']
#        fC = config['csv2sys']['attributesclosed']
#        nicdir = DATA.resolve_dir_nic(config['csv2sys']['nic_dir'], 
#                                      config['csv2sys']['nic_subdir'])
#        delim = UTIL.delim_norm(config['csv2sys']['delim'])
        LOG.debug(f'Assembling cache spec for {asofdate}')
        LOG.info(f'Building NICATTdf object for {asofdate}')
        NICATTdf = DATA.fetch_data(spec, DATA.case_NIC)[DATA.IDX_Attributes]
        LOG.info(f'Building FDICCBdf object for {asofdate}')
        FDICCBdf = DATA.fetch_data(spec, DATA.case_FDICCB)
        LOG.info(f'Building FDICSoDdf object for {asofdate}')
        FDICSoDdf = DATA.fetch_data(spec, DATA.case_FDICSoD)
#        FDICSoDcsv2df(config['csv2sys']['attributesclosed'])
        IDdata = (BankSys, NICATTdf, FDICCBdf, FDICSoDdf)
        # Make all the ID mappings
        LOG.info(f'Building CERT-HCR and CERT-RSSD objects for {asofdate}')
        (CERT2HCR, HCR2CERT, CERT2RSSD, RSSD2CERT) = \
            make_idmaps_CERT_HCR_RSSD(config, asofdate, IDdata)  
#        (CERT2RSSD,RSSD2CERT) = make_idmaps_CERT_RSSD(config, asofdate, IDdata)
        LOG.info(f'Building RSSD-HCR objects for {asofdate}')
        (RSSD2HCR,HCR2RSSD) = \
            make_idmaps_HCR_RSSD(config, asofdate, IDdata, RSSD2CERT, CERT2HCR)
        LOG.info(f'Building HH-RSSD objects for {asofdate}')
        (RSSD2HH,HH2RSSD) = make_idmaps_HH_RSSD(config, asofdate, IDdata)
        # Pack and pickle the ID mappings
        IDmaps = (CERT2HCR, HCR2CERT, CERT2RSSD, RSSD2CERT, 
                  RSSD2HCR, HCR2RSSD, RSSD2HH, HH2RSSD)
#        IDmaps = (CERT2HCR, HCR2CERT, CERT2RSSD, RSSD2CERT, 
#                  RSSD2HCR, HCR2RSSD)
#        IDlens = (len(CERT2HCR), len(HCR2CERT), len(CERT2RSSD), len(RSSD2CERT), 
#                  len(RSSD2HCR), len(HCR2RSSD))
        DATA.cache_data(spec, DATA.case_IDmaps, IDmaps)
#        LOG.info('CACHING: IDmaps file: '+id_filepath)
#        f = open(id_filepath, 'wb')
#        pik.dump(IDmaps, f)
#        f.close()
#        LOG.debug('CACHED: IDmaps file: '+id_filepath)
        IDlens = (len(CERT2HCR), len(HCR2CERT), len(CERT2RSSD), len(RSSD2CERT), 
                  len(RSSD2HCR), len(HCR2RSSD), len(RSSD2HH), len(HH2RSSD))
        LOG.debug(f'IDmaps have {IDlens} obs')
    else:
        LOG.info('ID maps found')
    return IDmaps
        
def make_idmaps_CERT_HCR_RSSD(config, asofdate, IDdata):
    """Scans the NIC and FDIC data for CERT-HCR mappings for the asofdate
    
    Many RSSDs will not have a CERT, because the CERT only applies to 
    FDIC-insured depository institutions, and there are many other 
    uninsured entity types in the NIC. If no CERT is found, the entity
    is ignored by both the CERT2HCR and HCR2CERT mappings. 
    
    Many CERTs will not have an HCR, because there exist free-standing 
    banks that are not a member of any holding company. If a CERT exists, 
    but no matching HCR is found, the entity will appear in the 
    CERT2HCR mapping (pointing to HCR=0), but there will be no 
    corresponding entry in the HCR2CERT mapping.

    """
    LOG.info('Beginning mappings for asofdate='+str(asofdate))
    trace_logging = ('TRUE'==config['csv2sys']['trace_logging'].upper())
    BankSys = IDdata[0]
    NICATTdf =  IDdata[1]
    FDICCBdf =  IDdata[2]
    FDICSoDdf =  IDdata[3]
    CERT2HCR = dict()
    HCR2CERT = dict()
    CERT2RSSD = dict()
    RSSD2CERT = dict()
    processed = set()
    SysU = BankSys.to_undirected()
    for E in BankSys:
        if (trace_logging):
            LOG.debug(f'Processing RSSD={E}, asofdate={asofdate}')
        if (E in processed):
            continue
        # Nodes in the graph component that contains entity E
        BHCnodes = nx.algorithms.components.node_connected_component(SysU, E)
        hcrset = set()
        for rssd in BHCnodes:
            cert_found = False
            hcr_found = False
            # First, check the NIC attributes to look up the RSSD and CERT
            dfQ = f'rssd=={rssd}'
            resultset = NICATTdf.query(dfQ)[['ID_FDIC_CERT', 'ID_RSSD']]
            for i in range(len(resultset)):
                cert = resultset.iat[i,0]
                if (cert > 0):
                    cert_found = True
                    continue
            if (len(resultset)>1):
                LOG.info(f'Multiple matches for RSSD={rssd} in NIC attributes')
            elif (len(resultset)<1):
                LOG.info(f'RSSD={rssd} not in NIC attrs, asofdate={asofdate}')
            # Fall back and try the FDIC SoD data
            if not(cert_found):
                dfQ = f'RSSDID=={rssd}'
                resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDID', 'RSSDHCR']]
                for i in range(len(resultset)):
                    cert = resultset.iat[i,0]
                    if (cert > 0):
                        cert_found = True
                        continue
            # CERT search is over, one way or the other
            if not(cert_found):
                if (trace_logging):
                    LOG.debug(f'No CERT found for RSSD={rssd}')
                processed.add(rssd)
            else:
                # We have a CERT, now look for the matching HCR
                dfQ = f'cert=={cert}'
                resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDID', 'RSSDHCR']]
                for i in range(len(resultset)):
                    hcr = resultset.iat[i,2]
                    if (hcr > 0):
                        hcr_found = True
                        hcrset.add(hcr)
                        continue
                if not(hcr_found):
                    # Last place to look for HCR is in the FDIC CB data
                    resultset = FDICCBdf.query(dfQ)[['CERT', 'RSSDHCR']]
                    for i in range(len(resultset)):
                        hcr = resultset.iat[i,1]
                        if (hcr > 0):
                            hcr_found = True
                            hcrset.add(hcr)
                            LOG.info(f'HCR={hcr} for {rssd} found in CB data')
                            continue
            if (cert_found):
                CERT2RSSD[cert] = rssd
                RSSD2CERT[rssd] = cert
                if (hcr_found):
                    CERT2HCR[cert] = hcr
                    HCR2CERT[hcr] = cert
                else:
                    CERT2HCR[cert] = 0
            else:
                RSSD2CERT[rssd] = 0
            processed.add(rssd)
        if (len(hcrset) > 1):
            LOG.info(f'Multiple HCRs ({hcrset}) for component w/RSSD={rssd}')
    LOG.info(f'Completing {len(processed)} mappings for asofdate={asofdate}')
    LOG.info(f'Completing CERT2HCR:{len(CERT2HCR)}, HCR2CERT:{len(HCR2CERT)}, CERT2RSSD:{len(CERT2RSSD)}, RSSD2CERT:{len(RSSD2CERT)}')
    return (CERT2HCR, HCR2CERT, CERT2RSSD, RSSD2CERT)

#def make_idmaps_CERT_HCR(config, asofdate, IDdata):
#    """Scans the NIC and FDIC CB data for CERT-HCR mappings for the asofdate
#    """
#    LOG.info('Beginning mappings for asofdate='+str(asofdate))
#    BankSys = IDdata[0]
#    NICATTdf =  IDdata[1]
#    FDICCBdf =  IDdata[2]
#    asoffdic = int(UTIL.rcnt_qtrend(asofdate)/100)
#    TRACE = ('TRUE'==config['csv2sys']['trace_logging'].upper())
#    CERT2HCR = dict()
#    HCR2CERT = dict()
#    processed = set()
#    NICATT_icol_ID_FDIC_CERT = NICATTdf.columns.get_loc('ID_FDIC_CERT')
#    FDICCB_icol_RSSDHCR = FDICCBdf.columns.get_loc('RSSDHCR')
#    SysU = BankSys.to_undirected()
#    for E in BankSys:
#        LOG.debug('Processing RSSD='+str(E)+', asofdate='+str(asofdate))
#        if (E in processed):
#            continue
#        # Nodes in the graph component that contains entity E
#        BHCnodes = nx.algorithms.components.node_connected_component(SysU, E)
#        for rssd in BHCnodes:
#            cert = -1
#            try:
#                cert = NICATTdf.loc[rssd].iat[NICATT_icol_ID_FDIC_CERT]
#            except (KeyError) as ke:
#                LOG.info('NIC lacks RSSD='+str(rssd)+
#                         ', asofdate='+str(asofdate))
#                LOG.info('Current progress: '+str(len(processed)))
#            FDICCBrow = None
#            if (cert > 0):
#                try:
#                    FDICCBrow = FDICCBdf.loc[(cert, asoffdic)]
#                except (KeyError) as ke:
#                    if (TRACE):
#                        LOG.debug('Missing CERT='+str(cert)+' in FDIC CB '+
#                          'for RSSD='+str(rssd)+' on asoffdic='+str(asoffdic))
#            if not(FDICCBrow is None):
#                hcr = FDICCBrow.iat[FDICCB_icol_RSSDHCR]
#                CERT2HCR[cert] = hcr
#                if (hcr > 0):
#                    # Not every CERT has an official HCR
#                    idmap_add(HCR2CERT, hcr, cert)
#            processed.add(rssd)
#    LOG.info('Completing mappings for asofdate='+str(asofdate)+
#             ', processed='+str(len(processed)))
#    return (CERT2HCR, HCR2CERT)


#def make_idmaps_CERT_RSSD(config, asofdate, IDdata):
#    """Scans the NIC data for CERT-RSSD mappings for the asofdate
#
#    Many RSSDs will not have a CERT, because the CERT only applies to 
#    FDIC-insured depository institutions, and there are many other 
#    uninsured entity types in the NIC. If no CERT is found, the entity
#    is ignored by both the CERT2HCR and HCR2CERT mappings. 
#    
#    """
#    LOG.info('Beginning mappings for asofdate='+str(asofdate))
#    BankSys = IDdata[0]
#    NICATTdf =  IDdata[1]
#    FDICSoDdf =  IDdata[3]
#    CERT2RSSD = dict()
#    RSSD2CERT = dict()
#    processed = set()
#    SysU = BankSys.to_undirected()
#    for E in BankSys:
#        LOG.debug('Processing RSSD='+str(E)+', asofdate='+str(asofdate))
#        if (E in processed):
#            continue
#        # Nodes in the graph component that contains entity E
#        BHCnodes = nx.algorithms.components.node_connected_component(SysU, E)
#        for rssd in BHCnodes:
#            cert_found = False
#            # First, check the NIC attributes to look up the RSSD and CERT
#            dfQ = f'rssd=={rssd}'
#            resultset = NICATTdf.query(dfQ)[['ID_FDIC_CERT', 'ID_RSSD']]
#            for i in range(len(resultset)):
#                cert = resultset.iat[i,0]
#                if (cert > 0):
#                    cert_found = True
#                    continue
#            if (len(resultset)>1):
#                LOG.info(f'Multiple matches for RSSD={rssd} in NIC attributes')
#            elif (len(resultset)<1):
#                LOG.info(f'RSSD={rssd} not found in NIC attributes')
#            # Fall back and try the FDIC SoD data
#            if not(cert_found):
#                dfQ = f'RSSDID=={rssd}'
#                resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDID', 'RSSDHCR']]
#                for i in range(len(resultset)):
#                    cert = resultset.iat[i,0]
#                    if (cert > 0):
#                        cert_found = True
#                        continue
#            # CERT search is over, one way or the other
#            if (cert_found):
#                CERT2RSSD[cert] = rssd
#                RSSD2CERT[rssd] = cert
#            else:
#                LOG.debug(f'No CERT found for RSSD={rssd}')
#                RSSD2CERT[rssd] = 0
#            processed.add(rssd)
#    LOG.info(f'Completing {len(processed)} mappings for asofdate={asofdate}')
#    return (CERT2RSSD, RSSD2CERT)                  

#def make_idmaps_CERT_RSSD(config, asofdate, IDdata):
#    """Scans the NIC data for CERT-RSSD mappings for the asofdate
#    """
#    LOG.info('Beginning mappings for asofdate='+str(asofdate))
#    BankSys = IDdata[0]
#    NICATTdf =  IDdata[1]
#    CERT2RSSD = dict()
#    RSSD2CERT = dict()
#    processed = set()
#    NICATT_icol_ID_FDIC_CERT = NICATTdf.columns.get_loc('ID_FDIC_CERT')
#    BSU = BankSys.to_undirected()
#    for E in BankSys:
#        LOG.debug('Processing RSSD='+str(E)+', asofdate='+str(asofdate))
#        if (E in processed):
#            continue
#        BHCnodes = nx.algorithms.components.node_connected_component(BSU, E)
#        for rssd in BHCnodes:
#            cert = -1
#            try:
#                cert = NICATTdf.loc[rssd].iat[NICATT_icol_ID_FDIC_CERT]
##                cert = NICATTdf.loc[rssd]['ID_FDIC_CERT']
#            except (KeyError) as ke:
#                LOG.info('NIC lacks RSSD='+str(rssd)+
#                         ', asofdate='+str(asofdate))
#            if (cert > 0):
#                CERT2RSSD[cert] = rssd
#                RSSD2CERT[rssd] = cert
#            else:
#                RSSD2CERT[rssd] = 0
#                LOG.debug('Missing CERT for RSSD='+str(rssd)+
#                          ' in NIC attributes for asofdate='+str(asofdate))
#            processed.add(rssd)
#    LOG.info('Completing mappings for asofdate='+str(asofdate)+
#             ', processed='+str(len(processed)))
#    return (CERT2RSSD, RSSD2CERT)                  

def make_idmaps_HCR_RSSD(config, asofdate, IDdata, RSSD2CERT, CERT2HCR):
    """Builds mappings between RSSD IDs and the associated high holder IDs
    """
    LOG.info(f'Beginning mappings for asofdate={asofdate}')
    BankSys = IDdata[0]
#    asoffdic = int(UTIL.rcnt_qtrend(asofdate)/100)
    trace_logging = ('TRUE'==config['csv2sys']['trace_logging'].upper())
    RSSD2HCR = dict()
    HCR2RSSD = dict()
    processed = set()
    SysU = BankSys.to_undirected()
    BHC_trickyset = set()
    for E in BankSys:
        if (trace_logging):
            LOG.debug(f'Processing RSSD={E}, asofdate={asofdate}')
        LOG.debug(f'Processing RSSD={E}, asofdate={asofdate}')
        if (E in processed):
            continue
        BHCnodes = nx.algorithms.components.node_connected_component(SysU, E)
        BHC = BankSys.subgraph(BHCnodes)
        # Begin with a triage, to set aside the tricky cases
        HCRset = set()
        for rssd in BHCnodes:
            if not(0==RSSD2CERT[rssd]):
                cert = RSSD2CERT[rssd]
                if not(0==CERT2HCR[cert]):
                    hcr = CERT2HCR[cert]
                    HCRset.add(hcr)
        if (0==len(HCRset)):
            for rssd in BHCnodes:
                RSSD2HCR[rssd] = 0
                processed.add(rssd)
        elif (1==len(HCRset)):
            hcr = HCRset.pop()
            (HH, case) = bhca.high_holder_impute(BHC)
            if not(HH==hcr):
                LOG.error(f'Mismatch between imputed HH ({HH}) vs HCR ({hcr})')
            for rssd in BHCnodes:
                RSSD2HCR[rssd] = hcr
                idmap_add(HCR2RSSD, hcr, rssd)
                processed.add(rssd)
        else:
            BHC_trickyset.add(BHC)
    LOG.info(f'Processed {len(processed)} of {len(BankSys)} as simple cases')
    LOG.info(f'There are {len(BHC_trickyset)} components in the tricky set')
    for BHC in BHC_trickyset:
        for E in BHC.nodes:
            if (E in processed):
                continue
            rssd_anc = nx.algorithms.dag.ancestors(BHC, E)
            HCRset = set()
            for rssd in rssd_anc:
                if not(0==RSSD2CERT[rssd]):
                    cert = RSSD2CERT[rssd]
                    if not(0==CERT2HCR[cert]):
                        hcr = CERT2HCR[cert]
                        HCRset.add(hcr)
            if (0==len(HCRset)):
                for rssd in rssd_anc:
                    RSSD2HCR[rssd] = 0
                    processed.add(rssd)
            elif (1==len(HCRset)):
                hcr = HCRset.pop()
                (HH, case) = bhca.high_holder_impute(BHC)
                if not(HH==hcr):
                    LOG.error(f'Mismatch between imputed HH ({HH}) vs HCR ({hcr})')
                for rssd in rssd_anc:
                    RSSD2HCR[rssd] = hcr
                    idmap_add(HCR2RSSD, hcr, rssd)
                    processed.add(rssd)
            else:
                LOG.error(f'Ambiguous HCR: {HCRset} for RSSD={rssd}')
                hcr = HCRset.pop()
                (HH, case) = bhca.high_holder_impute(BHC)
                if not(HH==hcr):
                    LOG.error(f'Mismatch between imputed HH ({HH}) vs HCR ({hcr})')
                for rssd in rssd_anc:
                    RSSD2HCR[rssd] = hcr
                    idmap_add(HCR2RSSD, hcr, rssd)
                    processed.add(rssd)
    LOG.info(f'Completing {len(processed)} mappings for asofdate={asofdate}')
    return (RSSD2HCR, HCR2RSSD)

#def make_idmaps_HCR_RSSD(config, asofdate, IDdata):
#    LOG.info('Beginning mappings for asofdate='+str(asofdate))
#    BankSys = IDdata[0]
#    NICATTdf =  IDdata[1]
#    FDICCBdf =  IDdata[2]
#    asoffdic = int(UTIL.rcnt_qtrend(asofdate)/100)
#    TRACE = ('TRUE'==config['csv2sys']['trace_logging'].upper())
#    RSSD2HCR = dict()
#    HCR2RSSD = dict()
#    processed = set()
#    BSU = BankSys.to_undirected()
#    for E in BankSys:
#        LOG.debug('Processing RSSD='+str(E)+', asofdate='+str(asofdate))
#        if (E in processed):
#            continue
#        BHCnodes = nx.algorithms.components.node_connected_component(BSU, E)
#        BHC = BankSys.subgraph(BHCnodes)
#        for rssd in BHCnodes:
#            climb_tree = True
#            hcr_found = False
#            context = rssd
#            while (climb_tree and not(hcr_found)):
#                # Inspect the context node for an HCR match ...
#                try:
#                    hcr = RSSD2HCR[context]
#                    RSSD2HCR[rssd] = hcr
#                    idmap_add(HCR2RSSD, hcr, rssd)
#                    hcr_found = True
#                    continue
#                except (KeyError) as ke:
#                    pass
#                # No existing match for context node, map instead from NIC/FDIC
#                ctx_cert = -1
#                FDICCBrow = None
#                try:
#                    ctx_cert = NICATTdf.loc[context]['ID_FDIC_CERT']
#                except (KeyError) as ke:
#                    LOG.info('NIC lacks RSSD='+str(context)+
#                             ', asofdate='+str(asofdate))
#                if (ctx_cert > 0):
#                    try:
#                        FDICCBrow = FDICCBdf.loc[(ctx_cert, asoffdic)]
#                    except (KeyError) as ke:
#                        if (TRACE):
#                          LOG.debug('FDIC CB lacks CERT='+str(ctx_cert)+
#                                    'for RSSD='+str(context)+
#                                    ' on asoffdic='+str(asoffdic))
#                if (FDICCBrow is not None):
#                    # We have an FDIC match for the context node
#                    hcr = FDICCBrow['RSSDHCR']
#                    if (hcr > 0):
#                        RSSD2HCR[context] = hcr
#                        RSSD2HCR[rssd] = hcr
#                        idmap_add(HCR2RSSD, hcr, context)
#                        idmap_add(HCR2RSSD, hcr, rssd)
#                        processed.add(context)
#                        hcr_found = True
#                        continue
#                # Still no match; climb the hierarchy and try the parent
#                parentage = BHC.in_edges(context)
#                if (0==len(parentage)):
#                    # Top of the house, nowhere to go
#                    climb_tree = False
#                elif (1==len(parentage)):
#                    context = list(parentage)[0][0]
#                else:
#                    # Ambiguous parentage
#                    climb_tree = False
#            processed.add(rssd)
#    LOG.info('Completing mappings for asofdate='+str(asofdate)+
#             ', processed='+str(len(processed)))
#    return (RSSD2HCR, HCR2RSSD)
#
def make_idmaps_HH_RSSD(config, asofdate, IDdata):
    LOG.info(f'Beginning mappings for asofdate={asofdate}')
    trace_logging = ('TRUE'==config['csv2sys']['trace_logging'].upper())
    BankSys = IDdata[0]
    # Maps from NIC RSSDs to imputed high-holder RSSDs, and vice versa
    RSSD2HH = dict()
    HH2RSSD = dict()
    processed = set()
    SysU = BankSys.to_undirected()
    for E in BankSys:
        if (trace_logging):
            LOG.debug(f'Processing RSSD={E}, asofdate={asofdate}')
        if (E in processed):
            continue
        BHCnodes = nx.algorithms.components.node_connected_component(SysU,E)
        BHC = BankSys.subgraph(BHCnodes)
        (BHCs, JVs) = bhca.high_holder_partition(BHC)
        for hc in BHCs.keys():
            # These are the BHC subgraphs with unambiguous high holder
            for ent in BHCs[hc].nodes:
                RSSD2HH[ent] = hc
                idmap_add(HH2RSSD, hc, ent)
                processed.add(ent)
        for jvk in JVs.keys():
            # These are the joint ventures, where the high holder is ambiguous
            max_hh_size = -1
            max_hh = None
            for i in range(len(jvk)):
                hh = jvk[i]
                bhc = BHCs[hh]
                if (len(bhc) > max_hh_size):
                    max_hh_size = len(bhc)
                    max_hh = hh
            for ent in JVs[jvk][0][0].nodes:
                RSSD2HH[ent] = max_hh 
                idmap_add(HH2RSSD, max_hh, ent)
                processed.add(ent)
    LOG.info(f'Completing {len(processed)} mappings for asofdate={asofdate}')
    return (RSSD2HH, HH2RSSD)

    
def idmap_add(idmap, key, val):
    """Adds value to the set (creating it, if needed) for the key in the idmap
    
    Examples
    --------
    >>> newmap = dict()
    >>> idmap_add(newmap, 1234, 9876)
    >>> sorted(newmap[1234])
    [9876]
    
    >>> idmap_add(newmap, 1234, 7654)
    >>> sorted(newmap[1234])
    [7654, 9876]
    """
    valset = None
    try:
        valset = idmap[key]
    except (KeyError) as ke:
        valset = set()
        idmap[key] = valset
    valset.add(val)
            


def build_idmaps(config):
    """Builds six maps among regulatory identifiers for each quarterly asofdate
    
    This function assembles the following mappings (dictionaries) between 
    regulatory IDs:
        
      * CERT2HCR - For a given FDIC CERT, provide the RSSD of the official
                   regulatory high holder BHC (the HCR, if any)
      * HCR2CERT - For a given regulatory high holder RSSD (the HCR), 
                   provide a list of its subsidiary FDIC CERTs
      * RSSD2HCR - For a given NIC RSSD, provide the RSSD of the official
                   regulatory high holder BHC (the HCR, if any)
      * HCR2RSSD - For a given regulatory high holder RSSD (the HCR), 
                   provide a list of its subsidiary NIC RSSDs
      * RSSD2HH -  For a given NIC RSSD, provide the RSSD of the official
                   regulatory high holder BHC (the HCR, if any); if 
                   no designated HCR exists, impute a high holder from 
                   the firm's structure
      * HH2RSSD -  For a given high holder RSSD (either an official HCR 
                   or an imputed high holder), provide a list of its 
                   subsidiary NIC RSSDs
    
    Parameters
    ----------
    config : ConfigParser object
        The key parameters are:
          * config['csv2sys']['clearcache_fdiccb']
          * config['csv2sys']['cachedir'] 
          * config['csv2sys']['asofdate0']
          * config['csv2sys']['asofdate1']
    """
#    if ('TRUE'==config['csv2sys']['clearcache_idmaps'].upper()):
#        # Remove cached pik files and recreate
#        LOG.warning('Clearing output cache of ID map *.pik files for: '+
#            config['csv2sys']['asofdate0']+'-'+config['csv2sys']['asofdate1'])
#        DATA.clear_cache(config['csv2sys']['cachedir'], 'IDmaps_',
#            config['csv2sys']['asofdate0'], config['csv2sys']['asofdate1'])
    asof_list = UTIL.assemble_asofs(config['csv2sys']['asofdate0'], 
                                    config['csv2sys']['asofdate1'])
    LOG.info(f'Building ID maps for asof_list={asof_list}')
    if (int(config['csv2sys']['parallel']) > 0):
        # Parallel processing of each requested asofdate
        LOG.info('Begin parallel processing ('+str(len(asof_list))+
                    ' tasks across '+config['csv2sys']['parallel']+
                    ' cores) of high-holders')
        pcount = min(int(config['csv2sys']['parallel']), 
                     os.cpu_count(), len(asof_list))
        pool = mp.Pool(processes=pcount)
        for asof in asof_list:
            pool.apply_async(get_idmaps, (config, asof))
        pool.close()
        pool.join()
        LOG.info('Complete parallel processing of high-holders')
    else:
        # Sequential processing of each requested asofdate
        LOG.info('Begin sequential processing of ID maps ('+ 
              str(len(asof_list))+' dates)')
        for asof in pb.progressbar(asof_list, redirect_stdout=True):
             get_idmaps(config, asof)
        LOG.info('Complete sequential processing of ID maps')
    


def main(argv=None):
    """A main function for command line execution
    
    This function parses the command line, loads the configuration, and 
    invokes the local function:
        
         * build_sys(config)
         
    Parameters
    ----------
    argv : dict
        The collection of arguments submitted on the command line
    """
    config = UTIL.parse_command_line(argv, __file__)
    try:
        build_sys(config)
#        build_hhdata(config)
        build_idmaps(config)
    except Exception as e:
        logging.exception("message")
    LOG.info('**** Processing complete ****')
    
# This tests whether the module is being run from the command line
if __name__ == "__main__":
    main()
    
