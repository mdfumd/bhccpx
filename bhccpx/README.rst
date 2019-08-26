BHCcpx Module Repository
========================

This project gathers tools for analyzing the complexity of bank holding
companies (BHCs), based on the internal ownership and control hierarchies
among a BHC's subsidiaries.

---------------

For additional background, please refer to the associated working paper:

Mark Flood, Dror Kenett, Robin Lumsdaine, and Jonathan Simon
`The Complexity of Bank Holding Companies: A Topological Approach <https://www.nber.org/papers/w23755>`_. 
NBER Working Paper 23755, August 2017

Large bank holding companies (BHCs) are structured into intricate ownership 
hierarchies involving hundreds or even thousands of legal entities. Each 
subsidiary in these hierarchies has its own legal form, assets, liabilities, 
managerial goals, and supervisory authorities. In the event of BHC default 
or insolvency, regulators may need to resolve the BHC and its constituent 
entities. Each entity individually will require some mix of cash infusion, 
outside purchase, consolidation with other subsidiaries, legal guarantees, 
and outright dissolution. The subsidiaries are not resolved in isolation, 
of course, but in the context of resolving the consolidated BHC at the top 
of the hierarchy. The number, diversity, and distribution of subsidiaries 
within the hierarchy can therefore significantly ease or complicate the 
resolution process. We propose a set of related metrics intended to assess 
the complexity of the BHC ownership graph. These proposed metrics focus on 
the graph quotient relative to certain well identified partitions on the 
set of subsidiaries, such as charter type and regulatory jurisdiction. The 
intended measures are mathematically grounded, intuitively sensible, and 
easy to implement. We illustrate the process with a case study of one 
large U.S. BHC. 
