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
import logging

import networkx as nx
import numpy as np
import pandas as pd
import csv
import multiprocessing as mp
import progressbar as pb

import bhc_data as DATA
import bhc_util as UTIL
import bhca
import csv2cch as c2c

MODNAME = __file__.split(os.path.sep)[-1].split('.')[0]
LOG = UTIL.log_config(MODNAME)



def extractBHC(config, asofdate, rssd, NICdata=None):
    """Extracts a single BHC graph from a full banking system graph
    
    A function to extract a single BHC graph from a full banking system 
    graph, starting with the BHC's high-holder RSSD ID    
    -- config is a configuration module, containing pointers to key files, etc.
    -- asofdate is the point in (historical) time of the desired BHC snapshot
    -- rssd is the RSSD ID of a particular BHC high-holder to extract
    The function also decorates the BHC graph with a number of important 
    attributes, by calling the add_attibutes().
    """
    BHC = None
#    if (config['sys2out']['indir'] != config['csv2sys']['indir']):
#        LOG.warning(f'csv2sys.indir: {config["csv2sys"]["indir"]} '+
#               f'differs from sys2out_indir: {config["sys2out"]["indir"]}')
#    if (config['sys2out']['outdir'] != config['csv2sys']['outdir']):
#        LOG.warning(f'csv2sys.outdir: {config["csv2sys"]["outdir"]} '+
#               f'differs from sys2out_outdir: {config["sys2out"]["outdir"]}')
    BankSys = DATA.fetch_data_banksys(config, asofdate, MODNAME)
    if (NICdata is None):
        LOG.debug('Fetching DATA (not provided)')
        NICdata = DATA.fetch_data_nic(config, asofdate, MODNAME)
    HHs = NICrel_breakdown(NICdata, asofdate)[0]
#    HHs = NICdata[DATA.IDX_HighHolder]
    if (rssd in HHs):
        BHC = populate_bhc(config, BankSys, NICdata, rssd)
        LOG.debug(f'BHC with high holder={rssd} has {BHC.number_of_nodes()} '+
            f'nodes and {BHC.number_of_edges()} edges, on asofdate={asofdate}')
        bhcfilename = 'NIC_'+str(rssd)+'_'+str(asofdate)+'.pik'
        bhcfilepath = os.path.join(config[MODNAME]['outdir'], bhcfilename)
        nx.write_gpickle(BHC, bhcfilepath)
    else:
        LOG.info(f'RSSD {rssd} missing from high-holder list as of {asofdate}')
    return BHC

def NICrel_breakdown(NICdata, asofdate):
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
    RELdf = NICdata[DATA.IDX_Relationships]
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



def populate_bhc(config, BankSys, NICdata, rssd):
    usebranches = ('TRUE'==config[MODNAME]['usebranches'].upper())
    node_atts = eval(config[MODNAME]['node_attributes'])
    edge_atts = eval(config[MODNAME]['edge_attributes'])
    bhc_entities = nx.algorithms.dag.descendants(BankSys, rssd)
    bhc_entities.add(rssd)    # Include HH in the BHC too
    BHC = BankSys.subgraph(bhc_entities)
    BHC = c2c.add_attributes_node(BHC, NICdata, node_atts)
    BHC = c2c.add_attributes_edge(BHC, NICdata, edge_atts)
    if (not(usebranches)):
        BHC = remove_branches (config, NICdata, BHC)
    BHC = BHC.to_directed()
    return BHC

def remove_branches(config, NICdata, BHC):
    """Duplicates a BHC DiGraph, but omitting the branch node
    
    Copies BHC to a new DiGraph object that is identical to BHC, except
    that it lacks any branches or subsidiaries of branches
    """
    removal_set = set()
    ATTdf = NICdata[DATA.IDX_Attributes]
    # Create a copy of the BHC that is editable
    BHC2 = nx.DiGraph()
    BHC2.add_nodes_from(BHC)
    BHC2.add_edges_from(BHC.edges)
    for node in BHC.nodes():
        node_id = np.int32(node)
        nicsource = 'XXX'
        # Copy any of node's existing attributes from BHC to BHC2
        nx.set_node_attributes(BHC2, {node: BHC.nodes()[node]})
        try:
            ent = ATTdf.loc[node_id]
            nicsource = ent['NICsource']
        except:
            pass
        if ('B'==nicsource):
            branchdescendants = nx.algorithms.dag.descendants(BHC,node)
            removal_set.add(node)
            for b in branchdescendants:
                removal_set.add(b)
    for n in removal_set:
        BHC2.remove_node(n)
    return BHC2


def edge_type_histogram(BHC, dim):
    """Builds a histogram of edge "types" based on the endpoints' node types
    
    Cycles through all edges of the graph, identifying the "type" of each
    endpoint, based on the value of the given (dim) attribute for the node.
    Creates and returns a dictionary of possible parent->child type 
    pairings as keys, and the number of instances of that type in the graph
    as values.
    """
    type_histogram = dict()
    dim_types = nx.get_node_attributes(BHC, dim)
    for e in BHC.edges():
        par = e[0]
        off = e[1]
        type_pair = str(dim_types[par])+"->"+str(dim_types[off])
        if (type_pair not in type_histogram):
            type_histogram[type_pair] = 0
        type_histogram[type_pair] = type_histogram[type_pair] + 1
    return type_histogram


def extract_bhcs_ondate(config, asofdate):
    """# Loads/creates a cached pik file for each BHC on the asofdate
    
    Returns
    -------
    List containing those BHCs
    """
    BHCs = []
    for rssd in eval(config[MODNAME]['bhclist']):
        BHC = extractBHC(config, asofdate, rssd)
        BHCs.append(BHC)
    return BHCs

def make_sample(config):
    """Extracts all BHCs in a given list, for a given list of as-of dates
    """
    asof_list = []    
    for YQ in eval(config[MODNAME]['asoflist']):
        asof_list.append(UTIL.make_asof(YQ)[0])
#    if ("TRUE"==config[MODNAME]['clearcache'].upper()):
#        clear_cache(config[MODNAME]['outdir'], asof_list)
    parallel_cores = int(config[MODNAME]['parallel'])
    if (parallel_cores > 0):
        LOG.warning(f'Begin parallel processing ({len(asof_list)} '+
                    f'tasks across {parallel_cores} cores')
        pcount = min(int(config[MODNAME]['parallel']), 
                     os.cpu_count(), len(asof_list))
        pool = mp.Pool(pcount)
        for asofdate in asof_list:
            pool.apply_async(extract_bhcs_ondate, (config, asofdate))
        pool.close()
        pool.join()
        LOG.info('Parallel processing complete')
    else:
        LOG.info('Beginning sequential processing for each asofdate')
        for asof in pb.progressbar(asof_list, redirect_stdout=True):
            extract_bhcs_ondate(config, asof)
        LOG.info('Sequential processing complete')


# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
        
        
# A dedicated function that produces the summary comparison of complexity
# measures for the Wachovia-Wells Fargo case study. This appears as Table 2
# in the NBER version of the paper
def make_persistent(config):
    if not('TRUE'==config[MODNAME]['make_persistent'].upper()):
        LOG.warning(f'Declined make_persistent, returning')
        return  # NOTE: Alternate return 
    BHCconfs = eval(config[MODNAME]['make_persistent_bhcqtrs'])
    # Loop to extract all of BHC snapshots
    pbar = ('TRUE'==config['DEFAULT']['progressbars'].upper())
    cutincr = int(config[MODNAME]['make_persistent_cutincr'])
    bar = None
    metrics = None
    BHCdict = {}
    idx = 0
    cutseq = [c for c in range(0, 101, cutincr)]
    LOG.info('Creating BHC networks')
    if (pbar): 
        widg = ['BHC-Qtr pairs: ', pb.Percentage(),' ', pb.Bar(),' ', pb.ETA()]
        bar = pb.ProgressBar(max_value=len(BHCconfs), widgets=widg)
        bar.start()
    for rssd, asof in BHCconfs:
        idx = idx + 1
        BHC = extractBHC(config, asof, rssd)
        BHCseq = bhca.edge_persistence_seq(BHC, 'pct_equity', cutseq)
        for b in BHCseq: 
#            print(f'BHCseq is: {type(BHCseq)} ({len(BHCseq)}), starting with {BHCseq[0]}; b={b.number_of_nodes()}')
#            print(f'BHCseq.index(b)={BHCseq.index(b)}')
            cut = cutseq[BHCseq.index(b)]
            key = f'{rssd}_{asof}_{cut}'
            metrics = bhca.complexity_workup(b)
            BHCdict[key] = [rssd, asof, cut] + list(metrics.values())
            if (pbar): bar.update(idx)
    if (pbar): 
        bar.update(bar.max_value)
        bar.finish(end='', dirty=True)
    cols = ['rssd', 'asofdate', 'cut'] + list(metrics.keys())
    persist_metrics = pd.DataFrame.from_dict(BHCdict, orient='index')
    persist_metrics.columns = cols
    persist_metrics['rssd'] = persist_metrics['rssd'].astype(int)
    persist_metrics.sort_values(['rssd', 'asofdate', 'cut'],
                                ascending=[True, True, True], inplace=True)
    if (pbar): 
        bar.update(bar.max_value)
        bar.finish(end='', dirty=True)
    csvdir = config[MODNAME]['make_persistent_outdir']
    csvfile = config[MODNAME]['make_persistent_outfile']
    csvpath = os.path.join(csvdir, csvfile)
    persist_metrics.to_csv(csvpath)
    LOG.warning('*** Processing persistent metrics ***')


# A dedicated function that produces the summary comparison of complexity
# measures for the Wachovia-Wells Fargo case study. This appears as Table 2
# in the NBER version of the paper
def make_wachwells(config):
    if not('TRUE'==config[MODNAME]['make_wachwells'].upper()):
        LOG.warning(f'Declined make_wachwells, returning')
        return  # NOTE: Alternate return 
    BHCconfigs = eval(config[MODNAME]['make_wachwells_bhcqtrs'])
    # Loop to extract all of BHC snapshots defined by BHCconfigs
    pbar = ('TRUE'==config['DEFAULT']['progressbars'].upper())
    bar = None
    metrics = None
    BHCdict = {}
    idx = 0
    LOG.info('Creating BHC networks')
    if (pbar): 
        widg = ['BHC-Qtr pairs: ', pb.Percentage(),' ', pb.Bar(),' ', pb.ETA()]
        bar = pb.ProgressBar(max_value=len(BHCconfigs), widgets=widg)
        bar.start()
    for rssd, asof in BHCconfigs:
        idx = idx + 1
        BHC = extractBHC(config, asof, rssd)
        metrics = bhca.complexity_workup(BHC)
        BHCdict[f'{rssd}_{asof}'] = [rssd, asof] + list(metrics.values())
        if (pbar): bar.update(idx)
    if (pbar): 
        bar.update(bar.max_value)
        bar.finish(end='', dirty=True)
    cols = ['rssd', 'asofdate'] + list(metrics.keys())
    table2 = pd.DataFrame.from_dict(BHCdict, orient='index')
    table2.columns = cols
    table2['rssd'] = table2['rssd'].astype(int)
    table2.sort_values(['rssd','asofdate'],ascending=[True,True],inplace=True)
    if (pbar): 
        bar.update(bar.max_value)
        bar.finish(end='', dirty=True)
    print(table2.iloc[:,2:6])
    print(table2.iloc[:,6:14])
    print(table2.iloc[:,14:22])
    LOG.warning('*** Processing Table2 Complete ***')
    return table2


def make_panel(config):
    """Creates a full panel of complexity measures
    
    Iterates through all BHCs for all quarters in a given list 
    of as-of dates between panel_asofdate0 and panel_asofdate1.
    """
    if not('TRUE'==config[MODNAME]['make_panel'].upper()):
        LOG.warning(f'Declined make_panel, returning')
        return  # NOTE: Alternate return 
    asof_list = UTIL.assemble_asofs(config[MODNAME]['panel_asofdate0'],
                                    config[MODNAME]['panel_asofdate1'])
    parallel_cores = int(config[MODNAME]['parallel'])
    if (parallel_cores > 0):
        LOG.warning(f'Begin parallel processing ({len(asof_list)} '+
                    f'tasks across {parallel_cores} cores')
        pcount = min(parallel_cores, os.cpu_count(), len(asof_list))
        pool = mp.Pool(pcount)
        results = [pool.apply_async(make_panel_asof, (config, asof)) 
                                    for asof in asof_list]
        results = [res.get() for res in results]
        pool.close()
        pool.join()
        LOG.debug('Parallel processing complete')
    else:
        LOG.info('Beginning sequential processing for each asofdate')
        results = []
        for asofdate in pb.progressbar(asof_list, redirect_stdout=True):
            results.append(make_panel_asof(config, asofdate))
        LOG.debug('Sequential processing complete')
    panelfilepath = os.path.join(config[MODNAME]['outdir'], 
                                 config[MODNAME]['panel_filename'])
    with open(panelfilepath, mode='w') as csvfile:
        fields = ['ASOF', 'RSSD'] + eval(config[MODNAME]['panel_metricslist'])
        csvwriter = csv.DictWriter(csvfile, fieldnames=fields)
        csvwriter.writeheader()
        # TODO: NEED TO SORT results BY ASOF AND RSSD BEFORE SAVING TO CSV
        for asof_dict in results:
            asofdate = asof_dict.pop('ASOF')
            for rssd in asof_dict.keys():
                metric_dict = asof_dict[rssd]
                metric_dict['ASOF'] = asofdate
                metric_dict['RSSD'] = rssd
                csvwriter.writerow(metric_dict)
    csvfile.close()
    LOG.info('**** Processing complete ****')


def make_panel_asof(config, asofdate):
    NICdata = DATA.fetch_data_nic(config, asofdate, MODNAME)
    BankSys = DATA.fetch_data_banksys(config, asofdate, MODNAME)
    HHs = eval(config[MODNAME]['panel_bhclist'])
    if (len(HHs) <= 0):
        HHs = NICrel_breakdown(NICdata, asofdate)[0]
    LOG.debug('Identified {len(HHs)} high-holders for {asofdate}')
    BHCs = dict()
    BHCs['ASOF'] = asofdate
    for rssd in HHs:
        BHC = populate_bhc(config, BankSys, NICdata, rssd)
        metrics = bhca.complexity_workup(BHC)
#        if ('TRUE'==config[MODNAME]['test_metrics'].upper()):
#            context = f'ASOF={asofdate}, RSSD={rssd}'
#            bhca.test_metrics(metrics, context)
        BHCs[rssd] = metrics
    return BHCs
        

def make_failscatter(config):
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
    """
#    if not('TRUE'==config[MODNAME]['make_failscatter'].upper()):
#        LOG.warning(f'Declined make_failscatter, returning')
#        return  # NOTE: Alternate return 
#    asof_list = UTIL.assemble_asofs(config[MODNAME]['failscatter_asofdate0'], 
#                                    config[MODNAME]['failscatter_asofdate1'])
#    LOG.info(f'As-of dates {asof_list}')
#    parallel_cores = int(config[MODNAME]['parallel'])
#    if (parallel_cores > 0):
#        LOG.warning(f'Begin parallel processing ({len(asof_list)} '+
#                    f'tasks across {parallel_cores} cores')
#        pcount = min(parallel_cores, os.cpu_count(), len(asof_list))
#        pool = mp.Pool(processes=pcount)
#        for asof in asof_list:
#            pool.apply_async(make_scatter, (config, asof))
#        pool.close()
#        pool.join()
#        LOG.info('Complete parallel processing')
#    else:
#        LOG.info(f'Begin sequential processing, ({len(asof_list)} dates)')
#        for asof in pb.progressbar(asof_list, redirect_stdout=True):
#            make_scatter(config, asof)
#        LOG.info('Complete sequential processing')
#
#
#def make_scatter(config, asofdate):
    if not('TRUE'==config[MODNAME]['make_failscatter'].upper()):
        LOG.warning(f'Declined make_failscatter, returning')
        return  # NOTE: Alternate return 
    FDICFailsdf = DATA.fetch_data_fdicfail(config, MODNAME)
    RSSD2CERT = dict()
    for yyyy in range(1994,2019):
        asof_mid = yyyy*10000+630
#        FDICSoDdf = DATA.fetch_data(spec, DATA.case_FDICSoD)
        FDICSoDdf = DATA.fetch_data_fdicsod(config, asof_mid, MODNAME)
        dfQ = f'RSSDID > 0'
        resultset = FDICSoDdf.query(dfQ)[['CERT', 'RSSDID']]
        RSSD2CERTsod = dict(zip(resultset.RSSDID, resultset.CERT))
        RSSD2CERT.update(RSSD2CERTsod)
    CERT2RSSD = dict()
    for k, v in RSSD2CERT.items():
        CERT2RSSD.setdefault(v, set()).add(k)
    extraatts = eval(config[MODNAME]['node_attributes'])
    metriccols = ['M_BVct_c', 'M_BEct_c', 'M_BCrk_c', 'M_BCmp_c', 
                    'M_EQfxB1_c', 'M_EQhxB1_c', 'M_EQfcB1_c', 'M_EQhcB1_c', 
                    'M_EQecB1_c', 'M_EDHmB1_c', 'M_EDHmM_c', 'M_ENlbl_c', 
                    'M_GQfxB1_c', 'M_GQhxB1_c', 'M_GQfcB1_c', 'M_GQhcB1_c', 
                    'M_GQecB1_c', 'M_GDHmB1_c', 'M_GDHmM_c', 'M_GNlbl_c',]
    FDICFailsdf = pd.concat([FDICFailsdf, pd.DataFrame(columns = metriccols)])
    for idx, row in FDICFailsdf.iterrows():
        cert = FDICFailsdf.loc[idx,'CERT']
        name = FDICFailsdf.loc[idx,'NAME']
        rssd = int(min(CERT2RSSD[FDICFailsdf.loc[idx,'CERT']]))
        FDICFailsdf.loc[idx,'RSSDID'] = rssd
        faildate = UTIL.timestamp2asofdate(FDICFailsdf.loc[idx, 'FAILDATE'])
        asofdate = UTIL.rcnt_qtrend(faildate)
        banksys = DATA.fetch_data_banksys(config, asofdate, MODNAME)
        LOG.debug(f'Processing CERT={cert} RSSD={rssd} NAME={name}')
        nicdat = DATA.fetch_data_nic(config, asofdate, MODNAME)
        try: 
            bhc_c = bhca.containing_component(banksys, rssd)
            c2c.add_attributes_node(bhc_c, nicdat, extraatts)
            metrics_c = bhca.complexity_workup(bhc_c)
            FDICFailsdf.loc[idx,'M_BVct_c'] = metrics_c[bhca.M_BVct]
            FDICFailsdf.loc[idx,'M_BEct_c'] = metrics_c[bhca.M_BEct]
            FDICFailsdf.loc[idx,'M_BCrk_c'] = metrics_c[bhca.M_BCrk]
            FDICFailsdf.loc[idx,'M_BCmp_c'] = metrics_c[bhca.M_BCmp]
            FDICFailsdf.loc[idx,'M_EQfxB1_c'] = metrics_c[bhca.M_EQfxB]
            FDICFailsdf.loc[idx,'M_EQhxB1_c'] = metrics_c[bhca.M_EQhxB]
            FDICFailsdf.loc[idx,'M_EQfcB1_c'] = metrics_c[bhca.M_EQfcB]
            FDICFailsdf.loc[idx,'M_EQhcB1_c'] = metrics_c[bhca.M_EQhcB]
            FDICFailsdf.loc[idx,'M_EQecB1_c'] = metrics_c[bhca.M_EQecB]
            FDICFailsdf.loc[idx,'M_EDHmB1_c'] = metrics_c[bhca.M_EDHmB]
            FDICFailsdf.loc[idx,'M_EDHmM_c'] = metrics_c[bhca.M_EDHmM]
            FDICFailsdf.loc[idx,'M_ENlbl_c'] = metrics_c[bhca.M_ENlbl]
            FDICFailsdf.loc[idx,'M_GQfxB1_c'] = metrics_c[bhca.M_GQfxB]
            FDICFailsdf.loc[idx,'M_GQhxB1_c'] = metrics_c[bhca.M_GQhxB]
            FDICFailsdf.loc[idx,'M_GQfcB1_c'] = metrics_c[bhca.M_GQfcB]
            FDICFailsdf.loc[idx,'M_GQhcB1_c'] = metrics_c[bhca.M_GQhcB]
            FDICFailsdf.loc[idx,'M_GQecB1_c'] = metrics_c[bhca.M_GQecB]
            FDICFailsdf.loc[idx,'M_GDHmB1_c'] = metrics_c[bhca.M_GDHmB]
            FDICFailsdf.loc[idx,'M_GDHmM_c'] = metrics_c[bhca.M_GDHmM]
            FDICFailsdf.loc[idx,'M_GNlbl_c'] = metrics_c[bhca.M_GNlbl]        
            bhc_r = banksys.subgraph(bhca.reachable_nodes(banksys, rssd))
            c2c.add_attributes_node(bhc_r, nicdat, extraatts)
            metrics_r = bhca.complexity_workup(bhc_r)
            FDICFailsdf.loc[idx,'M_BVct_r'] = metrics_r[bhca.M_BVct]
            FDICFailsdf.loc[idx,'M_BEct_r'] = metrics_r[bhca.M_BEct]
            FDICFailsdf.loc[idx,'M_BCrk_r'] = metrics_r[bhca.M_BCrk]
            FDICFailsdf.loc[idx,'M_BCmp_r'] = metrics_r[bhca.M_BCmp]
            FDICFailsdf.loc[idx,'M_EQfxB1_r'] = metrics_r[bhca.M_EQfxB]
            FDICFailsdf.loc[idx,'M_EQhxB1_r'] = metrics_r[bhca.M_EQhxB]
            FDICFailsdf.loc[idx,'M_EQfcB1_r'] = metrics_r[bhca.M_EQfcB]
            FDICFailsdf.loc[idx,'M_EQhcB1_r'] = metrics_r[bhca.M_EQhcB]
            FDICFailsdf.loc[idx,'M_EQecB1_r'] = metrics_r[bhca.M_EQecB]
            FDICFailsdf.loc[idx,'M_EDHmB1_r'] = metrics_r[bhca.M_EDHmB]
            FDICFailsdf.loc[idx,'M_EDHmM_r'] = metrics_r[bhca.M_EDHmM]
            FDICFailsdf.loc[idx,'M_ENlbl_r'] = metrics_r[bhca.M_ENlbl]
            FDICFailsdf.loc[idx,'M_GQfxB1_r'] = metrics_r[bhca.M_GQfxB]
            FDICFailsdf.loc[idx,'M_GQhxB1_r'] = metrics_r[bhca.M_GQhxB]
            FDICFailsdf.loc[idx,'M_GQfcB1_r'] = metrics_r[bhca.M_GQfcB]
            FDICFailsdf.loc[idx,'M_GQhcB1_r'] = metrics_r[bhca.M_GQhcB]
            FDICFailsdf.loc[idx,'M_GQecB1_r'] = metrics_r[bhca.M_GQecB]
            FDICFailsdf.loc[idx,'M_GDHmB1_r'] = metrics_r[bhca.M_GDHmB]
            FDICFailsdf.loc[idx,'M_GDHmM_r'] = metrics_r[bhca.M_GDHmM]
            FDICFailsdf.loc[idx,'M_GNlbl_r'] = metrics_r[bhca.M_GNlbl]
        except (KeyError):
            LOG.warning(f' -- Missing RSSD: {rssd}')
    outfile = 'MyScatter.csv'
    FDICFailsdf.to_csv(outfile, sep='\t')




# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================


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
#        make_sample(config)
        make_wachwells(config)
        make_panel(config)
        make_failscatter(config)
        make_persistent(config)
    except Exception as e:
        logging.exception("message")
    LOG.info('**** Processing complete ****')
    
# This tests whether the module is being run from the command line
if __name__ == "__main__":
    main()
    

## Create an SVG image file representing a BHC. The file is stored in the 
## outdir, with the filename: RSSD_<rssd_hh>_<asofdate>.svg.
## If popup is set to True, then the function will also launch a browser to
## display the file. 
#def makeSVG(BHC, outdir, rssd_hh, asofdate, popup=False):
#    svg_filename = 'RSSD'+'_'+str(rssd_hh)+'_'+str(asofdate)
#    svg_file = outdir + svg_filename
#    colormap = eval(CONFIG[MODNAME]['colormap'])
#    dot = gv.Digraph(comment='RSSD:'+str(rssd_hh), engine='dot')
#    dot.attr('node', fontsize='8')
#    dot.attr('node', fixedsize='true')
#    dot.attr('node', width='0.7')
#    dot.attr('node', height='0.3')
#    for N in BHC.nodes():
#        NM_LGL = ''
#        ENTITY_TYPE = 'ZZZ'
#        GEO_JURISD = 'ZZZ'
#        attribute_error = True
#        try:
#            NM_LGL = BHC.node[N]['nm_lgl'].strip()
#            ENTITY_TYPE = BHC.node[N]['entity_type']
#            GEO_JURISD = BHC.node[N]['GEO_JURISD']
#            attribute_error = False
#        except KeyError as KE:
#            print('WARNING: Invalid attribute data for RSSD='+str(N)+' at asofdate='+str(asofdate))
#        tt = '['+str(N)+']'+' '+ENTITY_TYPE+'\\n' +'------------\\n'+ NM_LGL +'\\n' +'------------\\n'+ GEO_JURISD
#        if (attribute_error):
#            dot.node('rssd'+str(N), str(N), style="filled", fillcolor="red;.5:green", tooltip=tt)
#        else:
#            fc = colormap[ENTITY_TYPE]
#            dot.node('rssd'+str(N), str(N), style="filled", fillcolor=fc, tooltip=tt)
#    for E in BHC.edges():
#        src = 'rssd' + str(E[0])
#        tgt = 'rssd' + str(E[1])
#        Vs = BHC.node[E[0]]
#        Vt = BHC.node[E[1]]
#        col = 'red'
#        if ('entity_type' not in Vs) or ('entity_type'not in Vt):
#            col='green'
#        elif Vs['entity_type']==Vt['entity_type']:
#            col='black'
#        dot.edge(src, tgt, arrowsize='0.3', color=col)
#    dot.render(filename=svg_file, format='svg')
#    dot.save(filename=svg_file+'.dot', directory=outdir)
#    if (popup):
#        os.system("%s %s" % (CONFIG[MODNAME]['browsercmd'], svg_file+'.svg'))


#def make_failscatter(config):
#    """ Create a CSV file correlating complexity with FDIC failure severity 
#    """
#    if not('TRUE'==config[MODNAME]['make_failscatter'].upper()):
#        LOG.warning(f'Declined make_failscatter, returning')
#        return  # NOTE: Alternate return 
##    spec = cache_spec(config)
#    FDICFailsdf = DATA.fetch_data(spec, DATA.case_FDICFail)
#    FDICFailsdf['LOA'] = FDICFailsdf.COST / FDICFailsdf.QBFASSET
#    FDICFailsdf['LOD'] = FDICFailsdf.COST / FDICFailsdf.QBFDEP
#    asof_list = UTIL.assemble_asofs(config[MODNAME]['failscatter_asofdate0'],
#                                    config[MODNAME]['failscatter_asofdate1'])
##    parallel_cores = int(config[MODNAME]['parallel'])
##    if (parallel_cores > 0):
##        LOG.info(f'Starting parallel processing ({parallel_cores} cores) for '+
##                 f'{len(asof_list)} dates (threads may trap process messages)')
##        pcount = min(int(config[MODNAME]['parallel']), os.cpu_count(), len(asof_list))
##        pool = mp.Pool(pcount)
##        results = [pool.apply_async(all_bhc_complex, (config, asof)) for asof in asof_list]
##        results = [res.get() for res in results]
##        pool.join()
##        LOG.debug('Parallel processing complete')
##        pool.close()
##    else:
#    LOG.info(f'Beginning sequential processing for {len(asof_list)} dates')
#    for idx,row in FDICFailsdf.iterrows():
#        cert = idx
#        faildate = UTIL.timestamp2asofdate( row['FAILDATE'])
#        rcntdate = UTIL.rcnt_qtrend(faildate)
#        specd = cache_spec(config, rcntdate)
#        IDmaps = DATA.fetch_data(spec, DATA.case_IDmaps)
##        cert = row['ID_FDIC_CERT']
##        rssd2cert[rssd] = cert
##        cert2rssd[cert] = rssd    
#    results = []
#    for asofdate in pb.progressbar(asof_list, redirect_stdout=True):
#        results.append(all_bhc_complex(config, asofdate))
#    LOG.debug('Sequential processing complete')
#    panelfilepath = os.path.join(config[MODNAME]['outdir'], config['sys2out']['panel_filename'])
#    with open(panelfilepath, mode='w') as csvfile:
#        fields = ['ASOF', 'RSSD'] + eval(config['sys2out']['metric_list'])
#        csvwriter = csv.DictWriter(csvfile, fieldnames=fields)
#        csvwriter.writeheader()
#        # TODO: NEED TO SORT results BY ASOF AND RSSD BEFORE SAVING TO CSV
#        for asof_dict in results:
#            asofdate = asof_dict.pop('ASOF')
#            for rssd in asof_dict.keys():
#                metric_dict = asof_dict[rssd]
#                metric_dict['ASOF'] = asofdate
#                metric_dict['RSSD'] = rssd
#                csvwriter.writerow(metric_dict)
#    csvfile.close()
#    LOG.info('**** Processing complete ****')
#  
   
## Deletes any DATA_* files in the cache corresponding to the dates in asof_list
#def clear_cache(cachedir, asof_list):
#    for asofdate in asof_list:
#        sysfilename = 'DATA_'+str(asofdate)+'.pik'
#        sysfilepath = os.path.join(cachedir, sysfilename)
#        if os.path.isfile(sysfilepath):
#            os.remove(sysfilepath)


     
