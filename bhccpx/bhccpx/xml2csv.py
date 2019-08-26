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

import sys
import os
import re
import logging
import progressbar as pb

# The config object defines key user-controlled variables
# If running from the command line, these choices may be overridden below
import bhc_util as UTIL

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
    xmlfiledir = UTIL.resolve_dir_nic(config['xml2csv']['nic_dir'], 
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
