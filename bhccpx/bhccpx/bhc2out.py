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

import os
import csv
import logging

import pandas as pd
import graphviz as gv
import progressbar as pb
import multiprocessing as mp

import csv2sys
import sys2bhc
import bhca

import bhc_data as DATA
import bhc_util as UTIL

LOG = logging.getLogger(__file__.split(os.path.sep)[-1].split('.')[0])



# A dedicated function that produces the summary comparison of complexity
# measures for the Wachovia-Wells Fargo case study. This appears as Table 2
# in the NBER version of the paper
def make_wachwells_comparison(config):
#    (verbose, veryverbose) = UTIL.verbosity(config)
    # BHCconfigs lists all of the (RSSD,asofdate) pairs to process
    # RSSD 1073551 is Wachovia Corp.
    # RSSD 1120754 is Wells Fargo & Co.
    BHCconfigs = [(1120754,20061231),
                  (1073551,20061231),
                  (1120754,20080930), 
                  (1073551,20080930),
                  (1120754,20081231),
                  (1120754,20101231)]
    # Loop to extract all of BHC snapshots defined by BHCconfigs
    pbar = ('TRUE'==config['DEFAULT']['progressbars'].upper())
    bar = None
    metrics = None
    BHCdict = {}
    idx = 0
    LOG.info('Creating BHC networks')
    if (pbar): 
        widg = ['BHC-Quarter pairs: ', pb.Percentage(),' ', pb.Bar(),' ', pb.ETA()]
        bar = pb.ProgressBar(max_value=len(BHCconfigs), widgets=widg)
        bar.start()
    for rssd, asof in BHCconfigs:
        idx = idx + 1
        BHC = sys2bhc.extractBHC(config, asof, rssd)
        metrics = complexity_workup(BHC)
        BHCdict[str(rssd)+'_'+str(asof)] = [rssd, asof] + list(metrics.values())
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


# String constants for metric names
M_BVct = 'Bas_Vertex_count'
M_BEct = 'Bas_Edge_count'
M_BCrk = 'Bas_Cycle_rank'
M_BCmp = 'Bas_Num_CComp'
M_EQfxB = 'Ent_Qfull_B1'
M_EQhxB = 'Ent_Qhetr_B1'
M_EQfcB = 'Ent_Qfcon_B1'
M_EQhcB = 'Ent_Qhcon_B1'
M_EQecB = 'Ent_edgcn_B1'
M_EDHmB = 'Ent_DjHom_B1'
M_EDHmM = 'Ent_DjHom_M'
M_ENlbl = 'Ent_Nlabl'
M_GQfxB = 'Geo_Qfull_B1'
M_GQhxB = 'Geo_Qhetr_B1'
M_GQfcB = 'Geo_Qfcon_B1'
M_GQhcB = 'Geo_Qhcon_B1'
M_GQecB = 'Geo_edgcn_B1'
M_GDHmB = 'Geo_DjHom_B1'
M_GDHmM = 'Geo_DjHom_M'
M_GNlbl = 'Geo_Nlabl'

def complexity_workup(BHC):
    """Calculates a standard set of complexity metrics for a BHC
    
    Most of the metrics involve quotienting the nodes of the BHC graph.
    The nodes are quotiented by entity type and (separately) geographic
    juridiction. 
    
    The following metrics are calculated:
        
    A. Basic metrics
    
      * Bas_Vertex_count = Number of nodes, original BHC graph 
      * Bas_Edge_count   = Number of edges, BHC graph 
      * Bas_Cycle_rank   = Cycle rank (b1), BHC graph 
      * Bas_Num_CComp    = Number of connected components (b0), BHC graph
      
    B. Entity quotients
    
      * Ent_Qfull_B1     = Cycle rank, full entity quotient
      * Ent_Qhetr_B1     = Cycle rank, heterogeneous entity quotient
      * Ent_Qfcon_B1     = Cycle rank, condensed entity quotient
      * Ent_Qhcon_B1     = Cycle rank, heterogeneous condensed entity quotient
      * Ent_edgcn_B1 = 'Ent_edgcn_B1'
      * Ent_DjHom_B1 = 'Ent_DjHom_B1'
      * Ent_DjHom_M = 'Ent_DjHom_M'
      * Ent_Nlabl = 'Ent_Nlabl'
      
    B. Beography quotients
    
      * Geo_Qfull_B1 = 'Geo_Qfull_B1'
      * Geo_Qhetr_B1 = 'Geo_Qhetr_B1'
      * Geo_Qfcon_B1 = 'Geo_Qfcon_B1'
      * Geo_Qhcon_B1 = 'Geo_Qhcon_B1'
      * Geo_edgcn_B1 = 'Geo_edgcn_B1'
      * Geo_DjHom_B1 = 'Geo_DjHom_B1'
      * Geo_DjHom_M = 'Geo_DjHom_M'
      * Geo_Nlabl = 'Geo_Nlabl'

    Parameters
    ----------
    BHC : networkx.DiGraph
        A directed graph representing a bank holding company
    
    Returns
    -------
    int
        Components in the projection of BHC to a simple undirected graph
        

    """
    
    metrics = dict()
    # Basic metrics, using the key constants defined above
    metrics[M_BVct] = BHC.number_of_nodes()
    metrics[M_BEct] = bhca.edge_count(BHC)
    metrics[M_BCrk] = bhca.cycle_rank(BHC)
    metrics[M_BCmp] = bhca.number_of_components(BHC)
    # Quotiented by entity type
    DIMEN = 'entity_type'
    QEF = bhca.get_quotient(BHC, DIMEN, bhca.Q_FULL)
    QEH = bhca.get_quotient(BHC, DIMEN, bhca.Q_HETERO)
    QEFC = bhca.get_quotient(BHC, DIMEN, bhca.Q_FULL_COND)
    QEHC = bhca.get_quotient(BHC, DIMEN, bhca.Q_HETERO_COND)
    CE = bhca.get_contraction(BHC, DIMEN).to_undirected()
    DMHE = bhca.get_disjoint_maximal_homogeneous_subgraphs(BHC, DIMEN)
    metrics[M_EQfxB] = bhca.cycle_rank(QEF)
    metrics[M_EQhxB] = bhca.cycle_rank(QEH)
    metrics[M_EQfcB] = bhca.cycle_rank(QEFC)
    metrics[M_EQhcB] = bhca.cycle_rank(QEHC)
    metrics[M_EQecB] = bhca.cycle_rank(CE)
    metrics[M_EDHmB] = bhca.cycle_rank(DMHE)
    metrics[M_EDHmM] = bhca.number_of_components(DMHE)
    metrics[M_ENlbl] = len(bhca.get_labels(BHC, DIMEN))
    # Quotiented by geographic jurisdiction
    DIMEN = 'GEO_JURISD'
    QGF = bhca.get_quotient(BHC, DIMEN, bhca.Q_FULL)
    QGH = bhca.get_quotient(BHC, DIMEN, bhca.Q_HETERO)
    QGFC = bhca.get_quotient(BHC, DIMEN, bhca.Q_FULL_COND)
    QGHC = bhca.get_quotient(BHC, DIMEN, bhca.Q_HETERO_COND)
    CG = bhca.get_contraction(BHC, DIMEN).to_undirected()
    DMHG = bhca.get_disjoint_maximal_homogeneous_subgraphs(BHC, DIMEN)
    metrics[M_GQfxB] = bhca.cycle_rank(QGF)
    metrics[M_GQhxB] = bhca.cycle_rank(QGH)
    metrics[M_GQfcB] = bhca.cycle_rank(QGFC)
    metrics[M_GQhcB] = bhca.cycle_rank(QGHC)
    metrics[M_GQecB] = bhca.cycle_rank(CG)
    metrics[M_GDHmB] = bhca.cycle_rank(DMHG)
    metrics[M_GDHmM] = bhca.number_of_components(DMHG)
    metrics[M_GNlbl] = len(bhca.get_labels(BHC, DIMEN))
    return metrics


# Calculates a full set of complexity metrics for a BHC, quotienting by
# both entity type and geographic jurisdiction, and returns them as a dict.
def test_metrics(metrics, context):
    # Ensure that the BHC is a single connected component
    if (1 != metrics[M_BCmp]):
        LOG.error('BHC is not completely connected. '+M_BCmp+': '+str(metrics[M_BCmp])+', Context: '+context)
    # Confirm that Equation 3 (Euler-Poincare) holds
    if (metrics[M_BCrk] != metrics[M_BEct] - metrics[M_BVct] + metrics[M_BCmp]):
        LOG.error('Euler-Poincare fails. '+M_BCrk+': '+str(metrics[M_BCrk])+', '+M_BEct+': '+str(metrics[M_BEct])+', '+M_BVct+': '+str(metrics[M_BVct])+', '+M_BCmp+': '+str(metrics[M_BCmp])+', Context: '+context)

    # Confirm that Theorem 3.2 Equation 6 (NBER version) holds -- entity type
    if (metrics[M_EQfxB] != metrics[M_BCrk] + metrics[M_BVct] - metrics[M_ENlbl]):
        LOG.error('Theorem 3.2 fails. '+M_EQfxB+': '+str(metrics[M_EQfxB])+', '+M_BCrk+': '+str(metrics[M_BCrk])+', '+M_BVct+': '+str(metrics[M_BVct])+', '+M_ENlbl+': '+str(metrics[M_ENlbl])+', Context: '+context)
    # Confirm that Corollary 3.4 (NBER version) holds -- entity type
    if (metrics[M_EQhxB] != metrics[M_EDHmM] - metrics[M_ENlbl] + metrics[M_BCrk] - metrics[M_EDHmB]):
        LOG.error('Corollary 3.4 fails. '+M_EQhxB+': '+str(metrics[M_EQhxB])+', '+M_EDHmM+': '+str(metrics[M_EDHmM])+', '+M_ENlbl+': '+str(metrics[M_ENlbl])+', '+M_BCrk+': '+str(metrics[M_BCrk])+', '+M_EDHmB+': '+str(metrics[M_EDHmB])+', Context: '+context)

    # Confirm that Theorem 3.2 Equation 6 (NBER version) holds -- geography
    if (metrics[M_GQfxB] != metrics[M_BCrk] + metrics[M_BVct] - metrics[M_GNlbl]):
        LOG.error('Theorem 3.2 fails (Geographic quotient). '+M_GQfxB+': '+str(metrics[M_GQfxB])+', '+M_BCrk+': '+str(metrics[M_BCrk])+', '+M_BVct+': '+str(metrics[M_BVct])+', '+M_GNlbl+': '+str(metrics[M_GNlbl])+', Context: '+context)
    # Confirm that Corollary 3.4 (NBER version) holds -- geography
    if (metrics[M_GQhxB] != metrics[M_GDHmM] - metrics[M_GNlbl] + metrics[M_BCrk] - metrics[M_GDHmB]):
        print('ERROR: Corollary 3.4 fails (Geographic quotient). '+M_GQhxB+': '+str(metrics[M_GQhxB])+', '+M_GDHmM+': '+str(metrics[M_GDHmM])+', '+M_GNlbl+': '+str(metrics[M_GNlbl])+', '+M_BCrk+': '+str(metrics[M_BCrk])+', '+M_GDHmB+': '+str(metrics[M_GDHmB])+', Context: '+context)

# Create an SVG image file representing a BHC. The file is stored in the 
# outdir, with the filename: RSSD_<rssd_hh>_<asofdate>.svg.
# If popup is set to True, then the function will also launch a browser to
# display the file. 
def makeSVG(BHC, outdir, rssd_hh, asofdate, popup=False):
    svg_filename = 'RSSD'+'_'+str(rssd_hh)+'_'+str(asofdate)
    svg_file = outdir + svg_filename
    colormap = eval(CONFIG['bhc2out']['colormap'])
    dot = gv.Digraph(comment='RSSD:'+str(rssd_hh), engine='dot')
    dot.attr('node', fontsize='8')
    dot.attr('node', fixedsize='true')
    dot.attr('node', width='0.7')
    dot.attr('node', height='0.3')
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
            print('WARNING: Invalid attribute data for RSSD='+str(N)+' at asofdate='+str(asofdate))
        tt = '['+str(N)+']'+' '+ENTITY_TYPE+'\\n' +'------------\\n'+ NM_LGL +'\\n' +'------------\\n'+ GEO_JURISD
        if (attribute_error):
            dot.node('rssd'+str(N), str(N), style="filled", fillcolor="red;.5:green", tooltip=tt)
        else:
            fc = colormap[ENTITY_TYPE]
            dot.node('rssd'+str(N), str(N), style="filled", fillcolor=fc, tooltip=tt)
    for E in BHC.edges():
        src = 'rssd' + str(E[0])
        tgt = 'rssd' + str(E[1])
        Vs = BHC.node[E[0]]
        Vt = BHC.node[E[1]]
        col = 'red'
        if ('entity_type' not in Vs) or ('entity_type'not in Vt):
            col='green'
        elif Vs['entity_type']==Vt['entity_type']:
            col='black'
        dot.edge(src, tgt, arrowsize='0.3', color=col)
    dot.render(filename=svg_file, format='svg')
    dot.save(filename=svg_file+'.dot', directory=outdir)
    if (popup):
        os.system("%s %s" % (CONFIG['bhc2out']['browsercmd'], svg_file+'.svg'))


# Create a full panel of complexity measures for all BHCs for all quarters in 
# the list of as-of dates between asofdate0 and asofdate1.
def make_panel(config):
#    (verbose, veryverbose) = UTIL.verbosity(config)
    if (veryverbose): UTIL.print_config(config, __file__)
    asof_list = UTIL.assemble_asofs(config['bhc2out']['asofdate0'], config['bhc2out']['asofdate1'])
    if (int(config['bhc2out']['parallel']) > 0):
        if (verbose): print('Beginning parallel processing for each asofdate (process messages may be trapped by parallel threads)')
        pcount = min(int(config['bhc2out']['parallel']), os.cpu_count(), len(asof_list))
        pool = mp.Pool(pcount)
        results = [pool.apply_async(all_bhc_complex, (config, asof)) for asof in asof_list]
        results = [res.get() for res in results]
        pool.close()
        pool.join()
        if (veryverbose): print('Parallel processing complete')
    else:
        if (verbose): print('Beginning sequential processing for each asofdate')
        results = []
        for asofdate in pb.progressbar(asof_list, redirect_stdout=True):
            results.append(all_bhc_complex(config, asofdate))
        if (veryverbose): print('Sequential processing complete')
    panelfilepath = os.path.join(config['bhc2out']['outdir'], config['bhc2out']['panel_filename'])
#    print('Finished up...', len(results), type(results), type(results[0]))
#    print('Checking', results[0]['ASOF'], results[16]['ASOF'])
    with open(panelfilepath, mode='w') as csvfile:
        fields = ['ASOF', 'RSSD'] + eval(config['bhc2out']['metric_list'])
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
    if (verbose): print('**** Processing complete ****')


def all_bhc_complex(config, asofdate):
#    (verbose, veryverbose) = UTIL.verbosity(config)
    fA = config['bhc2out']['attributesactive']
    fB = config['bhc2out']['attributesbranch']
    fC = config['bhc2out']['attributesclosed']
    fREL = config['bhc2out']['relationships']
    dats = DATA.makeDATA(config['bhc2out']['indir'], fA, fB, fC, fREL, asofdate)
    BankSys = csv2sys.make_banksys(config, asofdate)
    HHs = eval(config['bhc2out']['bhclist'])
    if (len(HHs) <= 0):
        HHs = dats[DATA.IDX_HighHolder]      # Include all RSSDs when HHs empty
    if (veryverbose): print('Identified', str(len(HHs)), 'high-holders for', asofdate)
    BHCs = dict()
    BHCs['ASOF'] = asofdate
    for rssd in HHs:
        BHC = sys2bhc.populate_bhc(config, BankSys, dats, rssd)
        metrics = complexity_workup(BHC)
        if ('TRUE'==config['bhc2out']['test_metrics'].upper()):
            context = 'ASOF='+str(asofdate)+', RSSD='+str(rssd)
            test_metrics(metrics, context)
        BHCs[rssd] = metrics
    return BHCs
        

def main(argv=None):
    """A main function for command line execution
    
    This function parses the command line, loads the configuration, and 
    invokes the local functions:
        
         * make_panel(config)
         * make_wachwells_comparison(config)
         
    Parameters
    ----------
    argv : dict
        The collection of arguments submitted on the command line
    """
    config = UTIL.parse_command_line(argv, __file__)
    if ('TRUE'==config['bhc2out']['make_panel'].upper()): 
        make_panel(config)
    if ('TRUE'==config['bhc2out']['make_wachwells_comparison'].upper()): 
        make_wachwells_comparison(config)
    LOG.warning('**** Processing complete ****')
    
    
# Test whether the module is invoked from the command line
if __name__ == "__main__":
    # Module unit tests
    import doctest
    doctest.testmod()
    # Invoking the primary functionality
    main()
