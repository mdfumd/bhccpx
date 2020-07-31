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

import sys
import os
import re
import logging
import progressbar as pb

# The config object defines key user-controlled variables
# If running from the command line, these choices may be overridden below
import bhc_util as UTIL
import bhc_data as DATA

LOG = UTIL.log_config(__file__.split(os.path.sep)[-1].split('.')[0])


def write_head(config, outfile, template):
    """ Writes the header row in the CSV output file.
    """
    LOG.debug('Entering function')
    delim = re.sub('<TAB>', '\t', config['xml2csv']['delim'])
    for k in template:
        outfile.write(k + delim)
    outfile.write('\n')
    
    

def write_elem(config, elem, outfile, template):
    """Writes an individual observation (elem) row in the CSV output file
    """
    delim = re.sub('<TAB>', '\t', config['xml2csv']['delim'])
    keyvals = {key:None for key in template}
    # Carve off the first XML tag as opentag, keeping the rest as elem
    (opentag, elem) = elem.split('>',1)
    # Now extract attributes in opentag as key-value pairs:  '<KEY>value'
    opentag = opentag.replace('<','')
    open_pairs = opentag.split(' ')
    open_pairs.pop(0)
    for i, p in enumerate(open_pairs):
        open_pairs[i] = '<'+p.replace('="', '>').replace('"','')
    # Replace element end tags with the CSV delimiter
    elem_mod = re.sub('</[A-Z_0-9]+>', delim, elem)
    # Strip off meaningless trailing blanks to avoid parsing errors
    elem_mod = elem_mod.rstrip()
    # Split elem into a list of key-value pairs, each of the form: '<KEY>value'
    keyval_pairs = open_pairs + elem_mod.split(delim) 
    # Populate the keyvals dict with values parsed from elem
    for kv in keyval_pairs:
        try:
            key_val = kv.split('>',1)
            key = key_val[0].replace('<','')
            val = key_val[1]
            keyvals[key] = val
        except Exception:
            print('ERROR: Cannot parse key-value pair', kv, elem)
    for k in template:
        val = keyvals[k]
        if (None==val):
            outfile.write(delim)
        else:
            outfile.write(keyvals[k] + delim)
    outfile.write('\n')
    
    
    
def parse_nic_file(config, xmlfilename):
    UTIL.tic()
    
    xmlfiledir = DATA.resolve_dir_nic(config['xml2csv']['nic_dir'], 
      config['xml2csv']['nic_subdir'])
    xmlfilepath = os.path.join(xmlfiledir, xmlfilename)
    chunksize = int(config['xml2csv']['chunksize'])
    
    # Inspect XML file to determine the type of data
    try:
        LOG.info('Parsing: '+xmlfilepath)
        infile = open(xmlfilepath, 'r')
    except (FileNotFoundError) as fnfe:
        LOG.error('Could not open: '+xmlfilepath+', skipping: '+str(fnfe)) 
        return
    chunk = infile.read(chunksize).upper()
    LOG.debug('Next chunk: '+chunk)
    infile.close()
    elemtype = ''
    if (chunk.find('<ATTRIBUTES')>=0):
        elemtype = 'ATTRIBUTES'
        template = eval(config['xml2csv']['attributestemplate'])
    elif (chunk.find('<RELATIONSHIP')>=0):
        elemtype = 'RELATIONSHIP'
        template = eval(config['xml2csv']['relationshipstemplate'])
    elif (chunk.find('<TRANSFORMATION')>=0):
        elemtype = 'TRANSFORMATION'
        template = eval(config['xml2csv']['transformationstemplate'])
    else:
        raise Exception('XML data not recognized as any of: '+
                        'ATTRIBUTES, RELATIONSHIP, TRANSFORMATION')
    LOG.warning('NIC element type: '+elemtype)

    # Open files to read XML and write CSV
    infile = open(xmlfilepath, 'r')
    LOG.warning('Opening for input: '+xmlfilepath)
    outfilename = (os.path.splitext(xmlfilename)[0]+
                   config['xml2csv']['outfileext'])
    outfilepath = os.path.join(config['xml2csv']['outdir'], outfilename)
    outfile = open(outfilepath, 'w', 1)
    LOG.warning('Opening for output: '+outfilepath)
    
    # Read from XML and write to CSV
    pbar = ('TRUE'==config['DEFAULT']['progressbars'].upper())
    needs_csv_head = True
    infileEOF = False
    xml0=''
    xml1='</DATA>'
    chunk = ''
    elem0='<'+elemtype
    elem1='</'+elemtype+'>'
    processed = 0
    xmlfilesize = os.path.getsize(xmlfilepath)
    sys.stdout.flush()
    bar = None
    if (pbar):
        widg = [xmlfilename+': ', pb.Percentage(),' ', pb.Bar(),' ', pb.ETA()]
        bar = pb.ProgressBar(max_value=xmlfilesize, widgets=widg)
        bar.start()
    while not(infileEOF):
        chunk = chunk + infile.read(chunksize)
        processed = processed + min(chunksize, xmlfilesize-processed)
        if ('TRUE'==config['DEFAULT']['progressbars'].upper()):
            bar.update(processed)
        if (len(xml0)<=0):
            idx0 = chunk.upper().find(elem0)
            xml0 = chunk.upper()[0:idx0]
            chunk = chunk.upper()[idx0:len(chunk)]
        # This next block means we have found the end point (elem1) of
        # the XML element (elem), and can proceed to parse and write it out
        while (chunk.upper().find(elem1)>=0):
            endelem = chunk.upper().find(elem1) + len(elem1)
            elem = chunk[0:endelem].upper()
            elem = elem.replace('&AMP;', '&')
            elem = elem.replace('&LT;', '<')
            elem = elem.replace('&GT;', '>')
            chunk = chunk[endelem:len(chunk)]
            if (needs_csv_head):
                write_head(config, outfile, template)
                needs_csv_head = False
            write_elem(config, elem, outfile, template)
        # When we find the XML end tag (xml1), shut things down
        if (chunk.upper().find(xml1)>=0):
            sys.stdout.write('\n' )
            infileEOF = True
            outfile.close()
            LOG.warning('Closing output file: '+ outfilepath)
            infile.close()
            LOG.warning('Closing input file: ', xmlfilepath)
    if (pbar):
        bar.update(bar.max_value)
        pb.utils.streams.flush()
        bar.finish(end='', dirty=True)
        
    LOG.info(f'PROFILE {UTIL.toc()} secs: parse_nic_file {xmlfilename}')


# This default invocation parses all five *.xml files, as
# indentified in the configuration object
def parse_nic(config):
    """Parses NIC XML files into equivalent CSV format
    
    The following files as specified in the configuration are processed:
      * config['xml2csv']['attributesactive']
      * config['xml2csv']['attributesbranch']
      * config['xml2csv']['attributesclosed']
      * config['xml2csv']['relationships']
      * config['xml2csv']['transformations']
    
    Parameters
    ----------
    config : configparser.ConfigParser
        Configuration object holding parameter choices
    """
    xml_input_files = [config['xml2csv']['attributesactive'],
                       config['xml2csv']['attributesbranch'],
                       config['xml2csv']['attributesclosed'],
                       config['xml2csv']['relationships'],
                       config['xml2csv']['transformations'] ]
    if ('TRUE'==config['xml2csv']['perform_xml_parse'].upper()):
        for xfile in xml_input_files:
            parse_nic_file(config, xfile)
    else:
        LOG.warning('Option perform_xml_parse==False; skipping parsing')



def main(argv=None):
    """A main function for command line execution
    
    This function parses the command line, loads the configuration, and 
    invokes the local function:
        
         * parse_nic(config)
         
    Parameters
    ----------
    argv : dict
        The collection of arguments submitted on the command line
    """
    config = UTIL.parse_command_line(argv, __file__)
    try:
        parse_nic(config)
    except Exception as e:
        logging.exception("message")
    LOG.info('**** Processing complete ****')
    
# This tests whether the module is being run from the command line
if __name__ == "__main__":
    main()
