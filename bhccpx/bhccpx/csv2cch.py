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

import networkx as nx
import multiprocessing as mp
import progressbar as pb

import bhc_data as DATA
import bhc_util as UTIL
import bhca

__all__ = ['add_attributes_node',
           'build_fdiccb',
           'build_fdicsod',
           'build_nic',
           'build_banksys',
           'main',
           ]
__version__ = '0.5'
__author__ = 'Mark D. Flood'

MODNAME = __file__.split(os.path.sep)[-1].split('.')[0]
LOG = UTIL.log_config(MODNAME)


def add_attributes_node(BHC, NICdata, attlist=None):
    """# Decorates a BHC graph with certain important attibutes
    
    For each node in the BHC, looks up a number of important attributes 
    from the attributes dataframe (in DATA), and attaches them to the node. 
    
    Parameters
    ----------
    BHC : NetworkX DiGraph 
        The bank holding company object to decorate
    NICdata : DataFrame
        A list of key information resources, from ATTcsv2df() and RELcsv2df()
    attlist : list
        A list of node attributes to include from NICdata
        
    Returns
    -------
    The decorated BHC graph is returned.

    """
    ATTdf = NICdata[DATA.IDX_Attributes]
    if (attlist is None):
        attlist = ['NICsource', 'ENTITY_TYPE', 'CNTRY_NM', 'STATE_ABBR_NM',]
    for node in BHC.nodes():
        try:
            ent = ATTdf.query(f'rssd=={int(node)}').iloc[0]
            node_dict = dict()
            node_dict = {
            # GEO_JURISD key is all caps, to distinguish as a derived attribute
            'GEO_JURISD': ent['CNTRY_NM'].strip() +' - '+ ent['STATE_ABBR_NM'],
            }
            # Now add the extra params requested in the config file
            for x in attlist:
                node_dict[x.lower()] = ent[x]
        except (IndexError) as IE:
            # Nodes (RSSDs) in NIC relationships can be no-shows in attributes 
            LOG.debug(f'Node attributes missing for {node}: {type(IE)}: {IE}')
            node_dict = {
            'nicsource': 'XXX',
            'entity_type': 'XXX',
            'GEO_JURISD': 'XXX'
            }
        except (Exception) as E:
            # An unexpected error has occurred
            LOG.error(f'Node attribute error for {node}: {type(E)}: {E}')
            node_dict = {
            'nicsource': 'XXX',
            'entity_type': 'XXX',
            'GEO_JURISD': 'XXX'
            }
        nx.set_node_attributes(BHC, {node: node_dict})
    return BHC


def add_attributes_edge(BHC, NICdata, attlist=None):
    """# Decorates a BHC graph with certain important attibutes
    
    For each edge in the BHC, looks up a number of important attributes 
    from the attributes dataframe (in DATA), and attaches them to the edge. 
    
    Parameters
    ----------
    BHC : NetworkX DiGraph 
        The bank holding company object to decorate
    NICdata : DataFrame
        A list of key information resources, from ATTcsv2df() and RELcsv2df()
    attlist : list
        A list of edge attributes to include from NICdata
        
    Returns
    -------
    The decorated BHC graph is returned.

    """
    RELdf = NICdata[DATA.IDX_Relationships]
    if (attlist is None):
        attlist = ['DT_START', 'DT_END',]
    for edge in BHC.edges():
        par = edge[0]
        off = edge[1]
        edge_dict = dict()
        try:
            rel = RELdf.query(f'rssd_par=={par} and rssd_off=={off}').iloc[0]
            for x in attlist:
                edge_dict[x.lower()] = rel[x]
        except (Exception) as E:
            # An unexpected error has occurred
            LOG.error(f'Edge attribute error for {edge}: {type(E)}: {E}')
        nx.set_edge_attributes(BHC, {edge: edge_dict})
    return BHC


def labeled_bhc(config, asofdate, hh_rssd):
    node_atts = eval(config[MODNAME]['node_attributes'])
    banksys = DATA.fetch_data_banksys(config, asofdate, MODNAME)
    BHCrset = bhca.reachable_nodes(banksys, hh_rssd)
    BHC = banksys.subgraph(BHCrset)
    NICdata = DATA.fetch_data_nic(config, asofdate, MODNAME)
    BHC = add_attributes_node(BHC, NICdata, node_atts)
    return BHC


def build_fdiccb(config):
    """Builds requested data objects for the FDIC CB data
    """
    asof_list = UTIL.assemble_asofs(config[MODNAME]['asofdate0'], 
                                    config[MODNAME]['asofdate1'])
    LOG.info('As-of dates: '+str(asof_list))
    parallel_cores = int(config[MODNAME]['parallel'])
    if (parallel_cores > 0):
        LOG.warning(f'Begin parallel processing ({len(asof_list)} '+
                    f'tasks across {parallel_cores} cores)')
        pcount = min(parallel_cores, os.cpu_count(), len(asof_list))
        pool = mp.Pool(processes=pcount)
        for asof in asof_list:
            pool.apply_async(DATA.fetch_data_fdiccb, (config, asof, MODNAME))
        pool.close()
        pool.join()
        LOG.info('Complete parallel processing')
    else:
        LOG.info(f'Begin sequential processing, ({len(asof_list)} dates)')
        for asof in pb.progressbar(asof_list, redirect_stdout=True):
             DATA.fetch_data_fdiccb(config, asof, MODNAME)
        LOG.info('Complete sequential processing')

    
    
def build_fdicsod(config):
    """Builds requested data objects for the FDIC SoD data
    """
    asof_list = UTIL.assemble_asofs(config[MODNAME]['asofdate0'], 
                                    config[MODNAME]['asofdate1'])
    LOG.info('As-of dates: '+str(asof_list))
    parallel_cores = int(config[MODNAME]['parallel'])
    if (parallel_cores > 0):
        LOG.warning(f'Begin parallel processing ({len(asof_list)} '+
                    f'tasks across {parallel_cores} cores)')
        pcount = min(parallel_cores, os.cpu_count(), len(asof_list))
        pool = mp.Pool(processes=pcount)
        for asof in asof_list:
            pool.apply_async(DATA.fetch_data_fdicsod, (config, asof, MODNAME))
        pool.close()
        pool.join()
        LOG.info('Complete parallel processing')
    else:
        LOG.info(f'Begin sequential processing, ({len(asof_list)} dates)')
        for asof in pb.progressbar(asof_list, redirect_stdout=True):
             DATA.fetch_data_fdicsod(config, asof, MODNAME)
        LOG.info('Complete sequential processing')

    
    
def build_nic(config):
    """Builds requested data objects for the NIC data
    """
    asof_list = UTIL.assemble_asofs(config[MODNAME]['asofdate0'], 
                                    config[MODNAME]['asofdate1'])
    LOG.info('As-of dates: '+str(asof_list))
    parallel_cores = int(config[MODNAME]['parallel'])
    if (parallel_cores > 0):
        LOG.warning(f'Begin parallel processing ({len(asof_list)} '+
                    f'tasks across {parallel_cores} cores)')
        pcount = min(parallel_cores, os.cpu_count(), len(asof_list))
        pool = mp.Pool(processes=pcount)
        for asof in asof_list:
            pool.apply_async(DATA.fetch_data_nic, (config, asof, MODNAME))
        pool.close()
        pool.join()
        LOG.info('Complete parallel processing')
    else:
        LOG.info(f'Begin sequential processing, ({len(asof_list)} dates)')
        for asof in pb.progressbar(asof_list, redirect_stdout=True):
             DATA.fetch_data_nic(config, asof, MODNAME)
        LOG.info('Complete sequential processing')

    
    
def build_banksys(config):
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
          * config['csv2cch']['nic_clearcache']
          * config['csv2cch']['cachedir'] 
          * config['csv2cch']['asofdate0']
          * config['csv2cch']['asofdate1']
    """
    asof_list = UTIL.assemble_asofs(config[MODNAME]['asofdate0'], 
                                    config[MODNAME]['asofdate1'])
    LOG.info('As-of dates: '+str(asof_list))
    parallel_cores = int(config[MODNAME]['parallel'])
    if (parallel_cores > 0):
        LOG.warning(f'Begin parallel processing ({len(asof_list)} '+
                    f'tasks across {parallel_cores} cores)')
        pcount = min(parallel_cores, os.cpu_count(), len(asof_list))
        pool = mp.Pool(processes=pcount)
        for asof in asof_list:
            pool.apply_async(DATA.fetch_data_banksys, (config, asof, MODNAME))
        pool.close()
        pool.join()
        LOG.info('Complete parallel processing')
    else:
        LOG.info(f'Begin sequential processing, ({len(asof_list)} dates)')
        for asof in pb.progressbar(asof_list, redirect_stdout=True):
             DATA.fetch_data_banksys(config, asof, MODNAME)
        LOG.info('Complete sequential processing')

    
    
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
        build_nic(config)
        build_banksys(config)
        build_fdicsod(config)
        build_fdiccb(config)
    except Exception as e:
        logging.exception("message")
    LOG.info('**** Processing complete ****')
    delta = UTIL.currtime() - starttime
    LOG.info(f'PROFILE total time: {delta} seconds (or {delta/60.0} minutes)')
    
# This tests whether the module is being run from the command line
if __name__ == "__main__":
    main()
    
    
