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
    """Creates a cache specification dictionary from a config object.
    """
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
    # FDIC fail items
    spec['fdicfail_dir'] = sect['fdicfail_dir']
    spec['fdicfail_file'] = sect['fdicfail_csv_file']
    spec['fdicfail_clearcache'] = ('TRUE'==sect['clearcache_fdicfail'].upper())
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
    UTIL.tic()
    sysfilename = 'NIC_'+str(asofdate)+'.pik'
    sysfilepath = os.path.join(config['csv2sys']['cachedir'], sysfilename)
    trace_logging = ('TRUE'==config['csv2sys']['trace_logging'].upper())
    spec = cache_spec(config, asofdate)
    BankSys = DATA.fetch_data(spec, DATA.case_BankSys)
    LOG.info(f'PROFILE fetch_data BankSys {asofdate} took {UTIL.toc()} secs.')

    LOG.info(f'Processing for asofdate={asofdate}')
    if (BankSys is None):
        LOG.info(f'Cache not found for BankSys as of {asofdate}, creating')
        NICdata = DATA.fetch_data(spec, DATA.case_NIC)
        ATTdf = NICdata[DATA.IDX_Attributes]
        RELdf = NICdata[DATA.IDX_Relationships]
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
                LOG.warning('As-of date: '+str(asofdate)+' out of bounds: '+
                    'rssd_par='+str(rssd_par)+', rssd_off='+str(rssd_off)+
                    ', date0='+str(date0)+', date1='+str(date1))
                continue   
            BankSys.add_edge(rssd_par, rssd_off)
        # Adding in the singleton institutions (no edges in Relationships file)
        LOG.info('System (pre), asofdate='+str(asofdate)+' has '+
                 str(BankSys.number_of_nodes())+' nodes and '+
                 str(BankSys.number_of_edges())+' edges ')
        ATTdf = ATTdf[ATTdf.DT_END >= asofdate]
        ATTdf = ATTdf[ATTdf.DT_OPEN <= asofdate]
        nodes_BankSys = set(BankSys.nodes)
        nodes_ATTdf = set(ATTdf['ID_RSSD'].unique())
        nodes_new = nodes_ATTdf.difference(nodes_BankSys)
        BankSys.add_nodes_from(nodes_new)
        LOG.info(f'PROFILE building BankSys {asofdate} took {UTIL.toc()} secs.')

        # Storing the new banking system file
        LOG.info(f'CACHING: System file: {sysfilepath}')
        DATA.cache_data(spec, DATA.case_BankSys, BankSys)
        LOG.info(f'PROFILE cache_data BankSys {asofdate} took {UTIL.toc()} secs.')

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



def get_idmaps(config, asofdate):
        
    spec = cache_spec(config, asofdate)
    IDmaps = DATA.fetch_data(spec, DATA.case_IDmaps)
    
    if (IDmaps is None):
        UTIL.tic()
        
        # Assembling the raw materials
        LOG.info(f'Building BankSys object for {asofdate}')
        BankSys = make_banksys(config, asofdate, read_data=True)
        LOG.info(f'PROFILE {UTIL.toc()} secs: make_banksys {asofdate}')
        
        LOG.debug(f'Assembling cache spec for {asofdate}')
        LOG.info(f'Building NICATTdf object for {asofdate}')
        NICdata = DATA.fetch_data(spec, DATA.case_NIC)
        NICATTdf = NICdata[DATA.IDX_Attributes]
        LOG.info(f'PROFILE {UTIL.toc()} secs: fetch NICdata {asofdate}')
        
#        LOG.info(f'Building FDICCBdf object for {asofdate}')
#        FDICCBdf = DATA.fetch_data(spec, DATA.case_FDICCB)
#        LOG.info(f'PROFILE {UTIL.toc()} secs: fetch FDICCBdf {asofdate}')
        
        LOG.info(f'Building FDICSoDdf object for {asofdate}')
        FDICSoDdf = DATA.fetch_data(spec, DATA.case_FDICSoD)
        FDICSoDdf = FDICSoDdf[FDICSoDdf['BRNUM']==0]        # Main office only
        LOG.info(f'PROFILE {UTIL.toc()} secs. fetch FDICSoDdf {asofdate}')
        
        # Building RSSD2CERT from BankSys, FDIC SoD, and NIC data
        RSSD2CERT = dict.fromkeys(BankSys.nodes, 0)
        dfQ = f'RSSDID > 0'
        resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDID']]
        RSSD2CERTsod = dict(zip(resultset.RSSDID, resultset.CERT))
        dfQ = f'ID_FDIC_CERT > 0'
        resultset = NICATTdf.query(dfQ)[['ID_FDIC_CERT', 'ID_RSSD']]
        RSSD2CERTnic = dict(zip(resultset.ID_RSSD, resultset.ID_FDIC_CERT))
        if (True):
            RSSD2CERT.update(RSSD2CERTsod)
            RSSD2CERT.update(RSSD2CERTnic)
        else:
            RSSD2CERT.update(RSSD2CERTnic)
            RSSD2CERT.update(RSSD2CERTsod)
        LOG.info(f'PROFILE {UTIL.toc()} secs. RSSD2CERT {asofdate}')
        
        # Reversing the mapping:  CERT2RSSD
        CERT2RSSD = dict()
        for k, v in RSSD2CERT.items():
            CERT2RSSD.setdefault(v, set()).add(k)
        LOG.info(f'PROFILE {UTIL.toc()} secs. CERT2RSSD {asofdate}')
        
        # Building CERT2HCR from FDIC SoD data
        dfQ = f'CERT > 0'
        resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDHCR']]
        CERT2HCR = dict(zip(resultset.CERT, resultset.RSSDHCR))
        LOG.info(f'PROFILE {UTIL.toc()} secs. CERT2HCR {asofdate}')
        
        # Reversing the mapping:  HCR2CERT
        HCR2CERT = dict()
        for k, v in CERT2HCR.items():
            HCR2CERT.setdefault(v, set()).add(k)
        LOG.info(f'PROFILE {UTIL.toc()} secs. HCR2CERT {asofdate}')
        
        IDmaps = (CERT2HCR, HCR2CERT, CERT2RSSD, RSSD2CERT)
        DATA.cache_data(spec, DATA.case_IDmaps, IDmaps)
        LOG.info(f'PROFILE {UTIL.toc()} secs: cache_data {asofdate}')
        
        IDlens = (len(CERT2HCR), len(HCR2CERT), len(CERT2RSSD), len(RSSD2CERT))
        LOG.debug(f'IDmaps have {IDlens} obs')
    else:
        LOG.info('ID maps found')
    return IDmaps



def make_fail_cpx(config):
    spec = cache_spec(config)
    FDICFailsdf = DATA.fetch_data(spec, DATA.case_FDICFail)
    FDICFailsdf['LOA'] = FDICFailsdf.COST / FDICFailsdf.QBFASSET
    FDICFailsdf['LOD'] = FDICFailsdf.COST / FDICFailsdf.QBFDEP
    
def get_idmaps_extra(config, asofdate):

#        # Assembling the raw materials
#        LOG.info(f'Building BankSys object for {asofdate}')
#        BankSys = make_banksys(config, asofdate, read_data=True)
#        BankSysU = BankSys.to_undirected()
#        LOG.info(f'PROFILE {UTIL.toc()} secs: make_banksys {asofdate}')
#        
#        LOG.debug(f'Assembling cache spec for {asofdate}')
#        LOG.info(f'Building NICATTdf object for {asofdate}')
#        NICdata = DATA.fetch_data(spec, DATA.case_NIC)
#        NICATTdf = NICdata[DATA.IDX_Attributes]
#        LOG.info(f'PROFILE {UTIL.toc()} secs: fetch NICdata {asofdate}')
        
        RSSD2HCR = dict.fromkeys(BankSys.nodes, -1)
        RSSD2HH = dict.fromkeys(BankSys.nodes, -1)
        HCR2RSSD = {0 : set()}
        HH2RSSD = {0 : set()}
        allHCRs = frozenset(HCR2CERT.keys())
        compct = nx.algorithms.components.number_connected_components(BankSysU)
        LOG.info(f'There are {compct} components, as of {asofdate}')
        comp_counter = 0 
        for comp in nx.connected_components(BankSysU):
            LOG.debug(f'Processing component {comp_counter}, as of {asofdate}')
            comp_counter = comp_counter + 1
            bhc = BankSys.subgraph(comp)
            (hh, case) = bhca.high_holder_impute(bhc)

            # Handle the HH calculations
            if (1==len(comp)):
                # Singleton component
                rssd = next(iter(comp))
                LOG.debug(f'Singleton node: {rssd}, as of {asofdate}')
                RSSD2HH[rssd] = rssd
                HH2RSSD[rssd] = comp
            else:
                # Component with multiple entities
                if (0==case):
                    # There is a unique HH node (entity of indegree==0)
                    RSSD2HH.update(dict.fromkeys(comp, hh))
                    HH2RSSD[hh] = comp
                elif (2==case):
                    # There is no HH node at all -- a bit strange, but possible
                    RSSD2HH.update(dict.fromkeys(comp, hh))
                    HH2RSSD[hh] = comp
                    fname = f'anomaly_{asofdate}_rogues_{hh}'
                    msg = f'Rogue nodes {comp}, as of {asofdate}. '
                    LOG.error(msg + ' SVG output to: '+fname)
                    log_error_svg(config, bhc, NICdata, fname, msg)
                else:
                    # There are multiple HH nodes 
                    indegrees = dict(bhc.in_degree())
                    HHs = set([n for n, d in indegrees.items() if 0==d])
                    eachHHs_desc = dict()
                    HHsize = dict()
                    for hhi in HHs:
                        eachHHs_desc[hhi] = bhca.reachable_nodes(bhc, hhi)
                        HHsize[hhi] = len(eachHHs_desc[hhi])
                    unassigned_nodes = comp.copy()
                    while len(HHsize)>0:
                        maxHH =  max(HHsize, key=lambda key: HHsize[key])
                        maxHH_desc = eachHHs_desc[maxHH]
                        RSSD2HH.update(dict.fromkeys(maxHH_desc, maxHH))
                        HH2RSSD[maxHH] = maxHH_desc
                        [unassigned_nodes.discard(n) for n in maxHH_desc]
                        del eachHHs_desc[maxHH]
                        del HHsize[maxHH]
                    if len(unassigned_nodes)>0:
                        # Rogue nodes (no HH) are still left over...
                        RSSD2HH.update(dict.fromkeys(unassigned_nodes, 0))
                        HH2RSSD[0].union(unassigned_nodes)

            # Handle the HCR calculations
            thisHCRs = frozenset(allHCRs.intersection(comp))
            if (len(thisHCRs)==0):
                fname = f'anomaly_{asofdate}_noHCR_{hh}'
                msg = f'No HCR for BHC with HH={hh}, as of {asofdate}. '
                LOG.debug(msg + ' SVG output to: '+fname)
            elif (len(thisHCRs)>1):
                fname = f'anomaly_{asofdate}_multHCRs_{hh}'
                msg = f'Multiple HCRs {set(thisHCRs)} for HH={hh},'+\
                      f' as of {asofdate}. '
                LOG.error(msg + ' SVG output to: '+fname)
                log_error_svg(config, bhc, NICdata, fname, msg)
                if not(hh in thisHCRs):
                    fname = f'anomaly_{asofdate}_HHmismatchHCR_{hh}'
                    msg = f'HH {hh} not in HCRs {set(thisHCRs)}, '+\
                          f'as of {asofdate}. '
                    LOG.error(msg + ' SVG output to: '+fname)
                    log_error_svg(config, bhc, NICdata, fname, msg)
            
            if (1==len(comp)):
                # Singleton component
                rssd = next(iter(comp))
                LOG.debug(f'Singleton node: {rssd}, as of {asofdate}')
#                RSSD2HH[rssd] = rssd
#                HH2RSSD[rssd] = comp
                if (rssd in allHCRs):
                    RSSD2HCR[rssd] = rssd
                    HCR2RSSD[rssd] = comp
                else:
                    RSSD2HCR[rssd] = 0
                    HCR2RSSD[0].add(rssd)
            elif (0==len(thisHCRs)):
                # Multiple nodes, no HCR: assign all nodes in comp to HCR 0
                RSSD2HCR.update(dict.fromkeys(comp, 0))
                HCR2RSSD[0].union(comp)
            elif (1==len(thisHCRs)):
                # Unique HCR, find its descendants and assign them
                hcr = next(iter(thisHCRs))
                hcr_desc = bhca.reachable_nodes(bhc, hcr)
                RSSD2HCR.update(dict.fromkeys(comp, 0))
                RSSD2HCR.update(dict.fromkeys(hcr_desc, hcr))
                HCR2RSSD[hcr] = hcr_desc
                hcr_complement = comp.difference(hcr_desc)
                if not(0==len(hcr_complement)):
                    fname = f'anomaly_{asofdate}_noHCRancestor_{hh}'
                    msg= f'Node(s) {hcr_complement} are unreachable from '+\
                         f'HCR {hcr}, as of {asofdate}. '
                    LOG.error(msg + ' SVG output to: '+fname)
                    log_error_svg(config, bhc, NICdata, fname, msg)
                    HCR2RSSD[0].union(hcr_complement)
            else:
                # Multiple HCRs - we need to know how they interact
                thisHCRs_desc = set()
                eachHCRs_desc = dict()
                for hcr in thisHCRs:
                    eachHCRs_desc[hcr] = bhca.reachable_nodes(bhc, hcr)
                    thisHCRs_desc.union(eachHCRs_desc[hcr])
                hcrnx = nx.DiGraph()
                for hcr in thisHCRs:
                    # Build a little digraph of only the HCRs in this component
                    hcr_desc = eachHCRs_desc[hcr].copy()
                    hcr_desc.remove(hcr)
                    hcr_subs = hcr_desc.intersection(thisHCRs)
                    hcr_rels = [(hcr,subhcr) for subhcr in hcr_subs]
                    hcrnx.add_edges_from(hcr_rels)
                    # Check the digraph for (messy) HCR cycles
                    hcrcycs = nx.algorithms.cycles.simple_cycles(hcrnx)
                    hcrcyclist = [c for c in hcrcycs]
                    if (len(hcrcyclist) > 0):
                        # TODO: If this ever occurs, we NEED a response
                        fname = f'anomaly_{asofdate}_HCRcycles_{hcr}'
                        msg= f'HCR {hcr} in mutual (HCR) ownership '+\
                             f'cycle(s) {hcrcyclist}, as of {asofdate}. '
                        LOG.error(msg + ' SVG output to: '+fname)
                        log_error_svg(config, bhc, NICdata, fname, msg)
                    # Work bottem up through hcrnx, assigning HCRs
#                    assigned = set()
                    while (hcrnx.number_of_nodes()>0):
                        outdegs = dict(hcrnx.out_degree())
                        OUT0 = set([n for n,d in outdegs.items() if 0==d])
                        for leaf in OUT0:
                            leaf_desc = eachHCRs_desc[leaf].copy()
                            RSSD2HCR.update(dict.fromkeys(leaf_desc, 0))
                            HCR2RSSD[leaf] = leaf_desc
#                            assigned = assigned.union(leaf_desc)
                            hcrnx.remove_node(leaf)
        LOG.info(f'PROFILE {UTIL.toc()} secs. RSSD2HCR {asofdate}')
        # Check whether we mis-assigned a bank
        RSSD2CERT_valid = {k:v for k,v in RSSD2CERT.items() if (v>0)}
        missing_certs = {v for v in RSSD2CERT_valid.values() if v not in CERT2HCR.keys()}
        CERT2HCR.update(dict.fromkeys(missing_certs, -1))
        LOG.info(f'Filtered RSSD2CERT mapping has {len(RSSD2CERT_valid)} elements')
        certs_mismatched = {k:v for k,v in RSSD2CERT_valid.items() if RSSD2HCR[k] != CERT2HCR[v]}
        LOG.info(f'PROFILE {UTIL.toc()} secs. HERE!! {len(certs_mismatched)} for {asofdate}')
        for rssd in certs_mismatched.keys():
            cert = certs_mismatched[rssd]
            try:
                hcr_mapped = RSSD2HCR[rssd]
                hcr_official = CERT2HCR[cert]
                comp = nx.algorithms.components.node_connected_component(BankSysU, rssd)
                bhc = BankSys.subgraph(comp)
                HCRid = hcr_mapped
                if (HCRid<=0):
                    HCRid = hcr_official
                if (HCRid<=0):
                    HCRid = f'T{int(UTIL.currtime())}'
                fname = f'anomaly_{asofdate}_HCRassgn_{HCRid}'
                msg= f'HCR misassigned for CERT={cert} and '+\
                     f'RSSD={rssd}, (HCR mapped={hcr_mapped} vs '+\
                     f'official={hcr_official}) as of {asofdate}. '
                LOG.error(msg + ' SVG output to: '+fname)
#                log_error_svg(config, bhc, NICdata, fname, msg)
            except KeyError as ke:
                comp = nx.algorithms.components.node_connected_component(BankSysU, rssd)
                bhc = BankSys.subgraph(comp)
                fname = f'anomaly_{asofdate}_CERTlack_{cert}'
                msg= f'CERT={cert} in NIC (for RSSD={rssd}) missing '+\
                     f'from FDIC SoD (HCR mapped={hcr_mapped}), '+\
                     f'as of {asofdate}. '
                LOG.error(msg + ' SVG output to: '+fname)
#                log_error_svg(config, bhc, NICdata, fname, msg)
#        LOG.info(f'CERT check reveals {hcrs_agree} agree '+\
#                 f'and {hcrs_disagree} disagree')
        LOG.info(f'PROFILE {UTIL.toc()} secs. checking CERT-HCRs {asofdate}')
                        
            

#def foobar():
#    # Fall back and try the FDIC SoD data
#    if not(cert_found):
#        dfQ = f'RSSDID=={rssd}'
#        resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDID', 'RSSDHCR']]
#        if len(resultset)>0:
#            cert = resultset.iat[0,0]
#            cert_found = True
#            if len(resultset)>1:
#                multcerts = set()
#                for i in range(len(resultset)):
#                    multcerts.add(resultset.iat[i,0])
#                if len(multcerts) > 0:
#                    msg = f'Multiple CERTs {multcerts} for RSSD={rssd} in FDIC SoD, as of {asofdate}'
#                    LOG.error(msg)
#                    BHC = BankSys.subgraph(BHCnodes)
#                    cert0 = sorted(multcerts)[0]
#                    filename = f'anomaly_ALY_{asofdate}_multCERTs_SoD_{cert0}'
#                    log_error_svg(config, BHC, NICdata, filename, msg)
##                for i in range(len(resultset)):
##                    cert = resultset.iat[i,0]
##                    if (cert > 0):
##                        cert_found = True
##                        continue
#    # CERT search is over, one way or the other
#    if not(cert_found):
#        if (trace_logging):
#            LOG.debug(f'No CERT found for RSSD={rssd}')
#    else:
#        # We have a CERT, now look for the matching HCR
#        dfQ = f'cert=={cert} and RSSDHCR>0'
#        resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDID', 'RSSDHCR']]
#        if len(resultset)>0:
#            hcr = resultset.iat[0,2]
#            hcr_found = True
#            if len(resultset)>1:
#                multhcrs = set()
#                for i in range(len(resultset)):
#                    multhcrs.add(resultset.iat[i,2])
#                if len(multhcrs) > 0:
#                    msg = f'Multiple HCRs {multhcrs} for RSSD={rssd} in FDIC SoD, as of {asofdate}'
#                    LOG.error(msg)
#                    BHC = BankSys.subgraph(BHCnodes)
#                    hcr0 = sorted(multhcrs)[0]
#                    filename = f'anomaly_ALY_{asofdate}_multHCRs_SoD_{hcr0}'
#                    log_error_svg(config, BHC, NICdata, filename, msg)
##                for i in range(len(resultset)):
##                    hcr = resultset.iat[i,2]
##                    if (hcr > 0):
##                        hcr_found = True
##                        hcrset.add(hcr)
##                        continue
#        if not(hcr_found):
#            # Last place to look for HCR is in the FDIC CB data
#            dfQ = f'cert=={cert} and RSSDHCR>0'
#            resultset = FDICCBdf.query(dfQ)[['CERT', 'RSSDHCR']]
#            for i in range(len(resultset)):
#                hcr = resultset.iat[i,1]
#                if (hcr > 0):
#                    hcr_found = True
#                    hcrset.add(hcr)
#                    LOG.debug(f'HCR={hcr} for {rssd} found in CB data')
#                    continue
#

        

def log_error_svg(config, BHC, NICdata, filename, msg):
    imagelogdir = config['DEFAULT']['imagelogdir']
    colormap = eval(config['bhc2out']['colormap'])
    extraattributes = eval(config['csv2sys']['extraattributes'])
    add_attributes(BHC, NICdata, attlist=extraattributes)
    DATA.log_bhc_svg(BHC, imagelogdir, filename, colormap, msg)

    
    
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
    


def add_attributes(BHC, NICdata, attlist=None):
    """# Decorates a BHC graph with certain important attibutes
    
    For each node in the BHC, looks up a number of important attributes 
    from the attributes dataframe (in DATA), and attaches them to the node. 
    
    Parameters
    ----------
    BHC : NetworkX DiGraph 
        The bank holding company object to decorate
    NICdata : DataFrame
        A list of key information resources, as assembled by ATTcsv2df()
    attlist : list
        A list of attributes to include from NICdata
    derive_geo : bool
        Whether to derive the 
        
    Returns
    -------
    The decorated BHC graph is returned.

    """
    ATTdf = NICdata[DATA.IDX_Attributes]
#    RELdf = NICdata[DATA.IDX_Relationships]
    if (attlist is None):
        attlist = ['NICsource', 'ENTITY_TYPE', 'CNTRY_NM', 'STATE_ABBR_NM',]
    for node in BHC.nodes():
#        node_id = 
        try:
            ent = ATTdf.query(f'rssd=={int(node)}').iloc[0]
            node_dict = dict()
#            ent = ATTdf.loc[node_id]
            node_dict = {
            'nicsource': ent['NICsource'],
            'entity_type': ent['ENTITY_TYPE'],
            'GEO_JURISD': ent['CNTRY_NM'].strip() +' - '+ ent['STATE_ABBR_NM'],
            }
            # Now add the extra params requested in the config file
#            extras = eval(config['sys2bhc']['extraattributes'])
            for x in attlist:
                node_dict[x.lower()] = ent[x]
        except:
            node_dict = {
            'nicsource': 'XXX',
            'entity_type': 'XXX',
            'GEO_JURISD': 'XXX'
            }
        nx.set_node_attributes(BHC, {node: node_dict})
    return BHC



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
    starttime = UTIL.currtime()
    config = UTIL.parse_command_line(argv, __file__)
    try:
        build_sys(config)
        build_idmaps(config)
    except Exception as e:
        logging.exception("message")
    LOG.info('**** Processing complete ****')
    delta = UTIL.currtime() - starttime
    LOG.info(f'PROFILE total time: {delta} seconds (or {delta/60.0} minutes)')
    
# This tests whether the module is being run from the command line
if __name__ == "__main__":
    main()
    
