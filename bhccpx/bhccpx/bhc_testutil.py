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

import random

import networkx as nx


entity_types=[
    # Holding companies 
    'BHC', # Bank Holding Company
    'FBH', # FBH Foreign Banking Organization as a BHC'
    'FHD', # FHD Financial Holding Company / BHC '
    'FHF', # FHF Financial Holding Company / FBO'
    'IHC', # IHC Intermediate Holding Company'
    'SLHC', # SLHC Savings and Loan Holding Company'
    # Bank-like entities 
    'AGB', # AGB Agreement Corporation - Banking'
    'CPB', # CPB Cooperative Bank'
    'EDB', # EDB Edge Corporation - Banking'
    'FSB', # FSB Federal Savings Bank'
    'FBK', # FBK Foreign Bank'
    'FBO', # FBO Foreign Banking Organization'
    'FCU', # FCU Federal Credit Union'
    'IBK', # IBK Int'l Bank of US Depository - Edge/Trust Co.'
    'NAT', # NAT National Bank'
    'NMB', # NMB Non-member Bank'
    'SAL', # SAL Savings & Loan Association'
    'SCU', # SCU State Credit Union'
    'SMB', # SMB State Member Bank'
    'SSB', # SSB State Savings Bank'
    # Branch-like entities 
    'DBR', # DBR Domestic Branch of a Domestic Bank'
    'EBR', # EBR Edge Corporation - Domestic Branch'
    'IBR', # IBR Foreign Branch of a U.S. Bank'
    'IFB', # IFB Insured Federal Branch of an F'=BO
    'ISB', # ISB Insured State Branch of an FBO
    'PST', # PST Pseudo-TWIG
    'REP', # REP Representative Office
    'TWG', # TWG TWIG
    'UFA', # UFA Uninsured Federal Agency of an FBO
    'UFB', # UFB Uninsured Federal Branch of an FBO
    'USA', # USA Uninsured State Agency of an FBO
    'USB', # USB Uninsured State Branch of an FBO
    # Investment entities 
    'AGI', # AGI Agreement Corporation - Investment
    'EDI', # EDI Edge Corporation - Investment
    'FNC', # FNC Finance Company
    'MTC', # MTC Non-deposit Trust Company - Member
    'NTC', # NTC Non-deposit Trust Company - Non-member
    'NYI', # NYI New York Investment Company
    'SBD', # SBD Securities Broker / Dealer
    # Miscellaneous entity types 
    'DEO', # DEO Domestic Entity Other
    'DPS', # DPS Data Processing Servicer
    'FEO', # FEO Foreign Entity Other
    'INB', # INB Int'l Non-bank Subs of Domestic Entities
    ]

state_abbr_nms = [
    'GA', 'WA', 'IL', 'TX', 'KS', 'MN', 'NY', 'MA', 'NE', 'CA', 'MO', 'WY', 
    'OK', 'WI', 'LA', 'ND', 'NJ', 'AL', 'IA', 'IN', 'MD', 'SC', 'AZ', 'FL',
    'HI', 'OR', 'PA', 'OH', 'MI', 'VA', 'ME', 'DE', 'WV', 'AR', 'CT', 'TN',
    'MS', 'NC', 'DC', 'NH', 'MT', 'UT', 'CO', 'GU', 'KY', 'NM', 'RI', 'NV',
    'SD', 'VT', 'ID', 'AK', 'PR', 'VI', 'MP', 'AS', '0 ',]

cntry_nms = [
    'UNITED STATES',
    'GUAM',
    'PUERTO RICO',
    'VIRGIN ISLANDS OF THE U.S.',
    'FEDERATED STATES OF MICRONESIA',
    'SOUTH AFRICA',
    'JAPAN',
    'LIBERIA',
    'THAILAND',
    'CAYMAN ISLANDS',
    'FRANCE (OTHER)',
    'SWITZERLAND (OTHER)',
    'CURACAO, BONAIRE, SABA, ST. MARTIN & ST.',
    'AUSTRALIA',
    'MEXICO',
    'HONG KONG',
    'BRAZIL',
    'ENGLAND',
    'ARGENTINA',
    'DOMINICAN REPUBLIC',
    'VENEZUELA',
    'SLOVENIA',
    'PANAMA',
    'INDONESIA (OTHER)',
    'PHILIPPINES',
    'TAIWAN',
    'PORTUGAL',
    'EL SALVADOR',
    'AUSTRIA',
    'ARUBA',
    'JAMAICA',
    'COSTA RICA',
    'GERMANY',
    'JERSEY',
    'LUXEMBOURG',
    'BELGIUM',
    'SPAIN',
    'CANADA',
    'BAHAMAS, THE',
    'MOROCCO (OTHER)',
    'HONDURAS',
    'NETHERLANDS',
    'CONGO (KINSHASA)',
    'SINGAPORE',
    'PARAGUAY',
    'COLOMBIA',
    'IRELAND',
    'EGYPT',
    'ZAMBIA',
    'NIGER',
    'CHINA, PEOPLES REPUBLIC OF',
    'NETHERLANDS ANTILLES',
    'TRINIDAD & TOBAGO (OTHER)',
    'KOREA, SOUTH',
    'NIGERIA',
    'ITALY (OTHER)',
    'TURKEY',
    'CHILE',
    'CHANNEL ISLANDS',
    'MACAU',
    'PERU',
    'URUGUAY',
    'MALAYSIA (OTHER)',
    'KENYA',
    'ECUADOR',
    'INDIA (OTHER)',
    'GREECE',
    'ISRAEL',
    'UNITED ARAB EMIRATES (OTHER)',
    'BAHRAIN',
    'DENMARK (OTHER)',
    'IRAN',
    'NEW ZEALAND (OTHER)',
    'NORWAY',
    'PAKISTAN',
    'SWEDEN',
    'SCOTLAND',
    'MACEDONIA (FORMER YUGOSLAV REPUBLIC OF)',
    'DUBAI',
    'SERBIA AND MONTENEGRO(FORMER YUGOSLAVIA)',
    'JORDAN',
    'KUWAIT',
    'QATAR',
    'SAUDI ARABIA',
    'BERMUDA',
    'COMORO ISLANDS',
    'FINLAND',
    'NORTHERN MARIANA ISLANDS',
    'ANGOLA',
    'LEBANON',
    'LIECHTENSTEIN',
    'FRENCH POLYNESIA',
    'TRINIDAD',
    'UGANDA',
    'SRI LANKA',
    'ETHIOPIA (OTHER)',
    'CYPRUS',
    'IRAQ',
    'HAITI',
    'CONGO (BRAZZAVILLE)',
    'RWANDA',
    'MONACO',
    'CAMEROON, UNITED REPUBLIC OF',
    'IVORY COAST',
    'TUNISIA',
    'GABON',
    'SENEGAL',
    'AFGHANISTAN',
    'WALES',
    'NICARAGUA',
    'TONGA',
    'ZIMBABWE',
    'GUATEMALA',
    'MALAWI',
    'NEW CALEDONIA',
    'BOLIVIA',
    'BOSNIA AND HERZEGOVINA',
    'CROATIA',
    'RUSSIA',
    'TURKS & CAICOS ISLANDS',
    'BARBADOS',
    'GIBRALTAR',
    'ANDORRA',
    'HUNGARY',
    'ICELAND',
    'BRITISH VIRGIN ISLANDS',
    'ROMANIA',
    'CZECH REPUBLIC',
    'POLAND',
    'TANZANIA, UNITED REPUBLIC OF',
    'GUYANA',
    'UNITED KINGDOM (OTHER)',
    'KAZAKHSTAN',
    'UKRAINE',
    'BHUTAN',
    'BRITISH WEST INDIES (OTHER)',
    'ESTONIA',
    'MAURITIUS',
    'SLOVAKIA',
    'GUERNSEY',
    'CAMBODIA',
    'TAJIKISTAN',
    'ARMENIA',
    'MONGOLIA',
    'ISLE OF MAN',
    'AZERBAIJAN',
    'BELARUS',
    'AMERICAN SAMOA',
    'MALTA',
    'SAINT LUCIA',
    'OMAN',
    'VIETNAM',
    'BRUNEI',
    'SWAZILAND',
    'ABU DHABI',
    'MARSHALL ISLANDS',
    'PALAU',
    'BANGLADESH',
    'ALGERIA',
    'BULGARIA',
    'NAMIBIA',
    'BOTSWANA',
    'MOZAMBIQUE',
    'LATVIA',
    'UZBEKISTAN',
    'GHANA',
    'LITHUANIA',
    ]



def BHC_simpleDAG(max_node_count=7):
    """Creates a BHC (a DiGraph object) with a simple directed tree structure.
    
    Parameters
    ----------
    max_node_count : int
        The maximum number of nodes in the BHC
        
    Returns
    -------
    BHC : nx.DiGraph
        The completed BHC graph
        
    Examples
    --------
    The default parameter values create a tree with 7 nodes
    
    >>> bhc7 = BHC_simpleDAG()
    >>> bhc7.number_of_nodes()
    7
    >>> bhc7.number_of_edges()
    6
    
    If we shrink the max_node_count by 1 (from 7 to 6) the tree is smaller
    
    >>> bhc3 = BHC_simpleDAG(max_node_count=6)
    >>> bhc3.number_of_nodes()
    3
    >>> bhc3.number_of_edges()
    2
    
    """
    depth = 0
    bfact = 2
    node_count = 1
    while (node_count<=max_node_count):
        depth = depth + 1
        node_count = node_count + pow(bfact,depth)
    BHC = nx.generators.classic.balanced_tree(2, depth-1, create_using=nx.DiGraph)
    return BHC
 


def BHC_systemDAG(num_comps = 2, max_node_count=7):
    """Creates a banking system (a DiGraph object) with several simpleDAGs.
    
    Parameters
    ----------
    num_comps : int, optional
        The number of BHC components in the system
    max_node_count : int, optional
        The maximum number of nodes in a typical BHC component
        
    Returns
    -------
    nx.DiGraph
        The completed system graph
        
    Examples
    --------
    The default parameter values create two trees with seven nodes each
    
    >>> BankSys = BHC_systemDAG()
    >>> nx.number_connected_components(BankSys.to_undirected())
    2
    >>> [len(c) for c in nx.connected_components(BankSys.to_undirected())]
    [7, 7]
    
    """
    BankSys = nx.DiGraph()
    for n in range(num_comps):
        ncnt = BankSys.number_of_nodes()
        bhc = BHC_simpleDAG(max_node_count)
        bhc = nx.relabel.convert_node_labels_to_integers(bhc, first_label=ncnt)
        BankSys.add_nodes_from(bhc.nodes())
        BankSys.add_edges_from(bhc.edges())
    return BankSys
 


def BHC_simpleDAG_plusreverseedge(max_node_count=7, reverse_edges=1):
    """Create a BHC as a DAG, except some edges have matching reverse edges
    
    Parameters
    ----------
    node_count
        The number of nodes in the BHC
    reverse_edges
        The number of reverse_edges to add
        
    Returns
    -------
    nx.DiGraph
        The completed BHC graph
        
    Examples
    --------
    The default creates a tree with 7 nodes, and one extra (reverse) edge
    
    >>> bhc71 = BHC_simpleDAG_plusreverseedge()
    >>> bhc71.number_of_nodes()
    7
    >>> bhc71.number_of_edges()
    7
    
    If we increase the reverse_edges (from 1 to 3) edge count goes up
    
    >>> bhc73 = BHC_simpleDAG_plusreverseedge(reverse_edges=3)
    >>> bhc73.number_of_nodes()
    7
    >>> bhc73.number_of_edges()
    9
    
    """
    BHC = BHC_simpleDAG(max_node_count)
    edges_unreversed = list(BHC.edges)
    for x in range(reverse_edges):
        e = random.choice(edges_unreversed)
        edges_unreversed.remove(e)
        BHC.add_edge(e[1],e[0])
    return BHC




def BHC_attribDAG(max_node_count=7, ent_labels=None, geo_labels=None):
    """Create a BHC as a DAG, adding entity type and geography attributes
    
    The default behavior (if ent_labels not provided),
    is to assign entity labels according to this rule:
    
        * First (root) node has entity type 'BHC'
        * Leaf nodes have entity type 'DEO'
        * All other nodes have entity type 'SMB'
    
    The default behavior (if geo_labels not provided),
    is to assign geographic labels according to this rule:
    
        * First (root) node has jurisdiction 'UNITED STATES - NY'
        * Last node has jurisdiction 'GERMANY'
        * Every (other) fifth node has jurisdiction 'UNITED STATES - CA'
        * Every other node has jurisdiction 'UNITED STATES - NC'
    
    Parameters
    ----------
    max_node_count : int
        The number of nodes in the BHC
    ent_labels : dict
        A mapping of entity type codes, keyed by node number
    geo_labels : dict
        A mapping of geographic jurisdiction codes, keyed by node number
        
    Returns
    -------
    BHC : nx.DiGraph
        The completed BHC graph, with entity and geography labels

    Examples
    --------
    The default creates a tree with 7 nodes, where the root node is a BHC in NY
    
    >>> bhc = TEST.BHC_attribDAG()
    >>> bhc.nodes(data=True)[0]
    {'entity_type': 'BHC', 'GEO_JURISD': 'UNITED STATES - NY'}
    >>> bhc.nodes(data=True)[1]
    {'entity_type': 'SMB', 'GEO_JURISD': 'UNITED STATES - NC'}
    >>> bhc.nodes(data=True)[2]
    {'entity_type': 'SMB', 'GEO_JURISD': 'UNITED STATES - NC'}
    >>> bhc.nodes(data=True)[3]
    {'entity_type': 'DEO', 'GEO_JURISD': 'UNITED STATES - NC'}
    >>> bhc.nodes(data=True)[4]
    {'entity_type': 'DEO', 'GEO_JURISD': 'UNITED STATES - NC'}
    >>> bhc.nodes(data=True)[5]
    {'entity_type': 'DEO', 'GEO_JURISD': 'UNITED STATES - CA'}
    >>> bhc.nodes(data=True)[6]
    {'entity_type': 'DEO', 'GEO_JURISD': 'GERMANY'}
    
    """
    BHC = BHC_simpleDAG(max_node_count)
    if (None==ent_labels):
        ent_labels = dict()
        first = True
        for nod in BHC:
            if (first): 
                ent_labels[nod] = 'BHC'
                first = False
            elif (0==len(list(BHC.successors(nod)))): 
                ent_labels[nod] = 'DEO'
            else: 
                ent_labels[nod] = 'SMB'
    if (None==geo_labels):
        geo_labels = dict()
        first = True
        for nod in BHC:
            if (first): 
                geo_labels[nod] = 'UNITED STATES - NY'
                first = False
            elif (0==nod%(len(BHC)-1)): 
                geo_labels[nod] = 'GERMANY'
            elif (0==nod%5): 
                geo_labels[nod] = 'UNITED STATES - CA'
            else:
                geo_labels[nod] = 'UNITED STATES - NC'
    nx.set_node_attributes(BHC, ent_labels, 'entity_type')
    nx.set_node_attributes(BHC, geo_labels, 'GEO_JURISD')
    return BHC



if __name__ == "__main__":
    import doctest
    doctest.testmod() 
 