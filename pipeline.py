
__author__ = "Katharina Moos"

import gzip
import gc
import argparse
import subprocess

import sys
import os
#sys.path.remove("/opt/rrzk/software/python/Python-3.4.3/lib/python3.4/site-packages/pysam-0.8.3-py3.4-linux-x86_64.egg")
#for d in os.walk("/home/kmoos3/.local"):
    #sys.path.append(d[0])
import pysam

import math
import numpy
import scipy.special

r_to_c_hg19 = {"NC_000001.10": "1", "NC_000002.11": "2", "NC_000003.11": "3", "NC_000004.11": "4", "NC_000005.9": "5",
               "NC_000006.11": "6", "NC_000007.13": "7", "NC_000008.10": "8", "NC_000009.11": "9", "NC_000010.10": "10",
               "NC_000011.9": "11", "NC_000012.11": "12", "NC_000013.10": "13", "NC_000014.8": "14", "NC_000015.9": "15",
               "NC_000016.9": "16", "NC_000017.10": "17", "NC_000018.9": "18", "NC_000019.9": "19", "NC_000020.10": "20",
               "NC_000021.8": "21", "NC_000022.10": "22", "NC_000023.10": "X", "NC_000024.9": "Y", "NC_012920.1": "MT"}
c_to_r_hg19 = {c: r for r, c in r_to_c_hg19.items()}
chr_numbers_hs = {"1": 0, "2": 0, "3": 0, "4": 0, "5": 0, "6": 0, "7": 0, "8": 0, "9": 0, "10": 0, "11": 0, "12": 0,
                  "13": 0, "14": 0, "15": 0, "16": 0, "17": 0, "18": 0, "19": 0, "20": 0, "21": 0, "22": 0,
                  "X": 0, "Y": 0, "MT": 0}

# Helper Functions #####################################################################################################

def average(values):
    """
    :param values: list of values for which the average is to be computed
    :return: the average
    """

    nvals = len(values)
    valsum = sum(values)
    avg = valsum/nvals
    return(avg)

def standard_deviation(values):
    """
    :param values:  list of values for which the standard deviation is to be computed
    :return: the standard deviation
    """

    nvals = len(values)
    valavg = average(values)
    sum_of_squares = 0
    for v in values:
        sum_of_squares += (valavg - v) ** 2
    stdev = math.sqrt(sum_of_squares / (nvals - 1))
    return(stdev)

def median(values):
    """
    :param values: list of values for which the median is to be computed
    :return: the median
    """

    values.sort()   # sorting values in ascending order
    n_values = len(values)
    if n_values % 2 == 1:
        med = values[int((n_values - 1) / 2)]
    else:
        med = (values[int(n_values / 2)] + values[int(n_values / 2 - 1)]) / 2

    return(med)

def str_is_int(s):
    """
    :param s:   string to be checked whether it represents an integer or not
    :return:    True if the input string represents an integer, False otherwise
    """

    try:
        int(s)
        return(True)
    except ValueError:
        return(False)

def write_annotation_statistics(outfile_stream, info_list, bam_file_des):
    """
    :param outfile_stream:      output file stream
    :param info_hash:           list storing annotation information (either multiple lists or a single number)
    :param bam_file_des:        list of BAM file descriptors

    :return: void

    This function writes annotation information to an output file, dependent on whether read count data is recorded.
    """

    if bam_file_des:                       # read counts are recorded
        outfile_stream.write(str(len(info_list[0])) + "\n")     # number of regions with no annotation
        outfile_stream.write("median read counts:")
        for i, c in enumerate(info_list):                       # iterate through BAM lists with region read counts
            outfile_stream.write("\t" + bam_file_des[i])            # write BAM file descriptor
            median_count = 0                                        # initialize BAM file's median region read count
            if c: median_count = median(c)                          # BAM file's median region read count
            outfile_stream.write("\t" + str(median_count))          # write median read count
    else:                                   # read counts are not recorded
        outfile_stream.write(str(info_list[0]))                 # number of regions with no annotation
    return()

def update_region_borders(region, minpos_hash, maxpos_hash, ref_start, ref_end):
    """
    :param region:          genomic region name
    :param minpos_hash:     hash storing region_name -> minimum position covered (0-based, including)
    :param maxpos_hash:     hash storing region_name -> maximum position covered (0-based, excluding)
    :param ref_start:       read starting position
    :param ref_end:         read ending position

    :return:    void

    This function updates the minimum and maximum position of a genomic region covered by reads.
    """

    if minpos_hash[region] > ref_start:
        minpos_hash[region] = ref_start
    if maxpos_hash[region] < ref_end:
        maxpos_hash[region] = ref_end

    return()

def labeling_ntpair(region_strand):
    """
    :param region_strand:   genomic region strand (':': unspecified, '+': reference reverse, '-': reference forward)
    :return:                a list containing the converted (original) and conversion (observed) nucleotide
    """

    if region_strand == "-":
        converted_nt = "A"  # reference forward region: "A" to be converted (labeled)
        conversion_nt = "G"  # reference forward region: "G" to be observed
    elif region_strand == "+":
        converted_nt = "T"  # reference reverse region: "T" to be converted (labeled)
        conversion_nt = "C"  # reference forward region: "C" to be observed
    else:
        converted_nt = None
        conversion_nt = None

    return([converted_nt, conversion_nt])

def rtable_file_generator(rtable_file):
    """
    :param rtable_file:     read table file
    :return:                generator; returning one genomic region's entry in each iteration
    """

    in_file = open(rtable_file, "r")        # opening read table
    line = in_file.readline()               # reading in first line
    gr = line.strip().split("#")[1]         # initialize current genomic region (first genomic region in file)
    line = in_file.readline().strip()       # reading in second line
    gr_entry = {gr: []}                     # initialize current genomic region's entry (empty list)
    table = []                              # initialize current table (empty list)

    while line:                             # iterating through read table file until end
        if line.startswith("#"):                # next genomic region encountered
            gr_entry[gr].append(table)          # storing current conversion table to current genomic region's entry
            yield(gr_entry)                     # returning current genomic region entry
            gr = line.split("#")[1]             # setting next genomic region
            gr_entry = {gr: []}                 # initializing next genomic region's entry
            table = []                          # re-setting table
        elif line.startswith(">SNP"):       # SNP table encountered
            pass                                # nothing to do, everything was prepared in the step before
        elif line.startswith(">CONV"):      # conversion table encountered
            gr_entry[gr].append(table)          # storing current SNP table to current genomic region's entry
            table = []                          # re-setting table
        else:                               # table row
            entries = line.split("\t")          # retrieving single row entries
            # re-transforming entries into integers, Nones and strings
            entries = [entries[0]] + [int(i) if i == '0' or i == '1' else None if i == 'None' else i for i in entries[1:]]
            table.append(entries)               # storing row entries to current table
        line = in_file.readline().strip()   # loading in next line

    gr_entry[gr].append(table)              # finally, adding last conversion table to last genomic region
    yield(gr_entry)                         # finally, returning last genomic region's entry

def table_from_cluster_results(infile, start, end, time_series):
    """
    :param infile:          file containing the clustering results
    :param start:           string indicating start of line in the file to start processing at (including)
    :param end:             string indicating start of line in the file to end processing at (excluding)
    :param time_series:     list containing time points at which measurement were taken

    :return:    a list containing (in the following order):
                    a hash storing row_name -> [transcript count ratios] for the nucleus (sorted by time)
                    a hash storing row_name -> [transcript count ratios] for the cytosol (sorted by time)

    This function retrieves data table entries from a clustering results input file. In detail, it retrieves the
    transcript count ratios of the nucleus as well as the cytosol, sorted by time.
    To find the location of the data table within the file, parameters 'start' and 'end' are to be defined. Parameter
    'start' defines a substring occurring at the start of the line which corresponds to the data table's first row (by
    which the data table can be found in the file), parameter 'end' defines a substring occurring at the start of the
    line which follows up the last row of the data table (first line that does not belong to the data table anymore,
    by which the end of the data table can be found in the file).
    """

    # initializing

    in_file = open(infile, "r")     # opening input file

    ratios_nu = {}      # hash storing row_name -> [ratios] for nucleus (sorted by time)
    ratios_cy = {}      # hash storing row_name -> [ratios] for cytosol (sorted by time)

    column_names = []   # initialize list storing table column names
    row_entries = {}    # initialize hash storing row_name -> [row entries]

    # iterating to input file part where table starts (and retrieving column names)

    for line in in_file:
        if line.startswith(start):
            column_names = line.strip().split("\t")
            break

    # getting time points' and ratios' column indices

    nu_colidx = [0 for i in time_series]    # initializing nuclear ratios time points' column indices
    cy_colidx = [0 for i in time_series]    # initializing cytosolic ratios time points' column indices

    for i, t in enumerate(time_series):     # iterating through time series
        for j, name in enumerate(column_names):  # searching correct column index (correct time and correct ratio)
            if "_" + (str(t) + "min") in name:      # correct time
                    if "nu" in name:                    # nucleus
                        nu_colidx[i] = j                    # assigning correct column index
                    elif "cy" in name:                  # cytosol
                        cy_colidx[i] = j                    # assigning correct column index

    # getting single rows

    for line in in_file:
        if line.startswith(end):
            break
        else:
            fields = line.strip().split("\t")   # getting single row entries
            row_name = fields[0][1:-1]          # getting row name
            row_entries[row_name] = fields[1:]  # assigning row's name its row values

    # getting ratio table values, sorted by time series

    for row_name in row_entries:            # iterating through row entries
        ratios_nu[row_name] = []                # initialize row name's hash entry with empty list (nucleus)
        ratios_cy[row_name] = []                # initialize row name's hash entry with empty list (cytosol)
        for i, t in enumerate(time_series):     # iterating through time series, appending ratios' row values
            nu_idx = nu_colidx[i]                   # nuclear mod/total column index
            cy_idx = cy_colidx[i]                   # cytosolic mod/total column index

            ratios_nu[row_name].append(row_entries[row_name][nu_idx])
            ratios_cy[row_name].append(row_entries[row_name][cy_idx])

    # closing input file and returning

    in_file.close()
    return([ratios_nu, ratios_cy])

def process_cigar_snps_conv(mode, cigar, converted_nt, conversion_nt,
                            read_id, read_seq, read_chr, read_region, ref_seq, ref_start,
                            region_insinfo, regions_convpos, snps={}):
    """
    :param mode:                mode to be executed; choose 'snp_all' to compute SNPs reporting all kind of conversions,
                                choose 'snp_exclusive' to compute SNPs reporting all but labeling-specific conversions,
                                choose 'conv' to compute conversions only, choose 'both' to compute both SNPs (exclu-
                                sive) and conversions
    :param cigar:               read cigar string to be processed
    :param converted_nt:        nucleotide converted due to metabolic labeling
    :param conversion_nt:       nucleotide observed due to metabolic labeling
    :param read_id:             read ID
    :param read_seq:            read sequence
    :param read_chr:            chromosome number the read is mapped to
    :param read_region:         genomic region name the read is mapped to
    :param ref_seq:             nucleotide sequence the read is mapped to
    :param ref_start:           reference starting position of the read's mapping
    :param region_insinfo:      hash storing read_region -> hash of potential conversion positions (0-based)
    :param regions_convpos:     hash storing read_region -> {position -> {read_ids with insertion at that position}}
    :param snps:                hash storing SNPs (will not be considered for conversions) (default: empty hash)

    :return:    a list containing (in the following order) a list storing read SNP information and a list storing read
                conversion information

    This function extracts read SNP and conversion information from a given cigar string. SNPs are denoted as None if
    the position is not covered by the read, 0 in case no SNP is present, and the SNP nucleotide otherwise (deletions
    are denoted with 'D', insertions with 'I' in additional columns at the end of the table). Conversions are stored as
    None in case that the position is not covered by the read, masked by a deletion or SNP correction, 0 in case no
    conversion is present, and as 1 otherwise.
    """

    read_snps = []  # list storing read SNP information
    read_conv = []  # list storing read conversion information
    read_pos = 0    # current read position (0-based)
    ref_pos = 0     # current reference position (0-based)

    if mode == "snp_all":           # MODE: all SNPs ...................................................................

        for c in cigar:     # filling read SNPs and conversions with the cigar string

            operation = c[0]                # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]                # length (stretch) of operation

            if operation == 4: continue     # skip soft-clipped positions

            if operation == 1:      # --- CASE: insertion ---

                ins_position = ref_start + ref_pos                  # insertion position start
                if ins_position in region_insinfo[read_region]:     # storing insertions
                    region_insinfo[read_region][ins_position][read_id] = None
                else:
                    region_insinfo[read_region][ins_position] = {read_id: None}
                read_pos += op_length                           # update read position (ref position does not change)

            elif operation == 2:    # --- CASE: deletion ---

                read_snps.extend(["D"] * op_length)     # storing deletions
                ref_pos += op_length                    # update reference position (read pos. does not change)

            else:                   # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]      # reference sequence nucleotide
                    read_nt = read_seq[readp]   # read nucleotide
                    if ref_nt == read_nt:       # real match
                        read_snps.append(0)         # no SNP recorded
                    else:                       # mismatch
                        read_snps.append(read_nt)   # SNP recorded

                read_pos += op_length   # update read position
                ref_pos += op_length    # update reference position

    elif mode == "snp_exclusive":   # MODE: all SNPs but labeling conversions ..........................................

        for c in cigar:     # filling read SNPs and conversions with the cigar string

            operation = c[0]                # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]                # length (stretch) of operation

            if operation == 4: continue     # skip soft-clipped positions

            if operation == 1:      # --- CASE: insertion ---

                ins_position = ref_start + ref_pos                  # insertion position start
                if ins_position in region_insinfo[read_region]:     # storing insertions
                    region_insinfo[read_region][ins_position][read_id] = None
                else:
                    region_insinfo[read_region][ins_position] = {read_id: None}
                read_pos += op_length                           # update read position (ref position does not change)

            elif operation == 2:    # --- CASE: deletion ---

                read_snps.extend(["D"] * op_length)     # storing deletions
                ref_pos += op_length                    # update reference position (read pos. does not change)

            else:                   # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]      # reference sequence nucleotide
                    read_nt = read_seq[readp]   # read nucleotide

                    if ref_nt == read_nt:       # real match
                        read_snps.append(0)         # no SNP recorded

                    else:                       # mismatch
                        if ref_nt == converted_nt and read_nt == conversion_nt:     # excluding labeling conversions
                            read_snps.append(0)         # no SNP recorded
                        else:
                            read_snps.append(read_nt)   # SNP recorded

                read_pos += op_length       # update read position
                ref_pos += op_length        # update reference position

    elif mode == "conv":            # MODE: conversions ................................................................

        for c in cigar:     # filling read SNPs and conversions with the cigar string

            operation = c[0]                # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]                # length (stretch) of operation
            snp_key_i = read_chr + "_"      # SNP sites key (initialized, incomplete)

            if operation == 4: continue     # skip soft-clipped positions

            if operation == 1:      # --- CASE: insertion ---

                read_pos += op_length                           # update read position (ref position does not change)

            elif operation == 2:    # --- CASE: deletion ---

                ref_nts = ref_seq[ref_pos:ref_pos + op_length]  # getting reference nucleotides
                for p, nt in enumerate(ref_nts):                # iterating through reference nucleotides
                    if nt == converted_nt:                          # conversion is stored as None (masked by deletion)
                        read_conv.append(None)                          # position is stored to position hash
                        regions_convpos[read_region][ref_start + ref_pos + p] = None
                ref_pos += op_length                            # update reference position (read pos. does not change)

            else:                   # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]                          # reference sequence nucleotide
                    read_nt = read_seq[readp]                       # read nucleotide

                    if ref_nt == read_nt:                           # real match
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:     # masked by SNP correction
                                read_conv.append(None)                                  # (stored as None)
                            else:                                                   # not converted
                                read_conv.append(0)                                     # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None

                    else:                                           # mismatch
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:     # masked by SNP correction
                                read_conv.append(None)                                  # (stored as None)
                            elif read_nt == conversion_nt:                          # is converted
                                read_conv.append(1)                                     # (stored as 1)
                            else:                                                   # another SNP
                                read_conv.append(0)                                     # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None

                read_pos += op_length   # update read position
                ref_pos += op_length    # update reference position

    elif mode == "both":            # MODE: both SNPs and conversions ..................................................

        for c in cigar:     # filling read SNPs and conversions with the cigar string

            operation = c[0]                # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]                # length (stretch) of operation
            snp_key_i = read_chr + "_"      # SNP sites key (initialized, incomplete)

            if operation == 4: continue     # skip soft-clipped positions

            if operation == 1:      # --- CASE: insertion ---

                ins_position = ref_start + ref_pos                  # insertion position start
                if ins_position in region_insinfo[read_region]:     # storing insertions
                    region_insinfo[read_region][ins_position][read_id] = None
                else:
                    region_insinfo[read_region][ins_position] = {read_id: None}
                read_pos += op_length                           # update read position (ref position does not change)

            elif operation == 2:    # --- CASE: deletion ---

                ref_nts = ref_seq[ref_pos:ref_pos + op_length]  # getting reference nucleotides
                for p, nt in enumerate(ref_nts):                # iterating through reference nucleotides
                    if nt == converted_nt:                          # conversion is stored as None (masked by deletion)
                        read_conv.append(None)                          # position is stored to position hash
                        regions_convpos[read_region][ref_start + ref_pos + p] = None
                read_snps.extend(["D"] * op_length)             # storing deletions
                ref_pos += op_length                            # update reference position (read pos. does not change)

            else:                   # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]                          # reference sequence nucleotide
                    read_nt = read_seq[readp]                       # read nucleotide

                    if ref_nt == read_nt:                           # real match
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:     # masked by SNP correction
                                read_conv.append(None)                                  # (stored as None)
                            else:                                                   # not converted
                                read_conv.append(0)                                     # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None
                        read_snps.append(0)                                         # no SNP recorded

                    else:                                           # mismatch
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:     # masked by SNP correction
                                read_conv.append(None)                                  # (stored as None)
                            elif read_nt == conversion_nt:                          # is converted
                                read_conv.append(1)                                     # (stored as 1)
                            else:                                                   # another SNP
                                read_conv.append(0)                                     # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None
                        read_snps.append(read_nt)                                   # SNP recorded

                read_pos += op_length   # update read position
                ref_pos += op_length    # update reference position

    return([read_snps, read_conv])  # RETURNING ........................................................................

def process_cigar_ntpairs(ref_seq, read_seq, cigar, chr, ref_start, snp_ai_hash):
    """
    :param ref_seq:         references' nucleotide sequence the read is aligned to
    :param read_seq:        read nucleotide sequence
    :param cigar:           read cigar string
    :param chr:             reference chromosome the read is aligned to
    :param ref_start:       reference position of read alignment start
    :param snp_ai_hash:     hash storing potential SNP and A-to-I positions as ChrNumber_ChrPos

    :return:    a list containing (in the following order):
                    list of positions of insertions (position after insertion event, 1-based)
                    list of positions of deletions (position of deletion event, 1-based)
                    confusion table (hash) assigning each nucleotide pairing its abundance within the read's mapping

    This function extracts read SNP information from a given cigar string. It records positions of insertions (position
    before insertion event, w.r.t. reference, 1-based), positions of deletions (position before deletion event, w.r.t.
    reference, 1-based) and nucleotide pairing counts.
    Nucleotide pairings counts are stored in a hash, where the first nucleotide refers to the reference, the second to
    the read. 'XI' and 'XD' refer to insertions and deletions, respectively, where X is the next reference nucleotide
    after the insertion event or, in the case of deletions, the first nucleotide in the reference deleted. The counts
    refer to single nucleotides inserted or deleted rather than to insertion or deletion events.
    """

    insertions = list()     # initialize list storing insertion positions (positions w.r.t. reference,
                            # position before insertion event, 1-based)
    deletions = list()      # initialize list storing deletion positions (positions w.r.t. reference,
                            # position before deletion event, 1-based)
    conf_table = {"AA": 0, "AT": 0, "AG": 0, "AC": 0, "AN": 0, "AI": 0, "AD": 0,
                  "TA": 0, "TT": 0, "TG": 0, "TC": 0, "TN": 0, "TI": 0, "TD": 0,
                  "GA": 0, "GT": 0, "GG": 0, "GC": 0, "GN": 0, "GI": 0, "GD": 0,
                  "CA": 0, "CT": 0, "CG": 0, "CC": 0, "CN": 0, "CI": 0, "CD": 0,    # initialize confusion table
                  "NA": 0, "NT": 0, "NG": 0, "NC": 0, "NN": 0, "NI": 0, "ND": 0}    # (all counts are 0)

    # computing confusion table

    read_pos = 0        # current read position (0-based)
    ref_pos = 0         # current reference position (0-based)
    for c in cigar:     # filling confusion table with the cigar string

        operation = c[0]        # type of operation is encoded by numbers (m, ins, del)
        op_length = c[1]        # length (stretch) of operation
        snp_key_i = chr + "_"   # SNP and A-to-I sites key (initialized, incomplete)

        if operation == 4: continue     # skip soft-clipped positions

        if operation == 1:              # CASE: insertion
            ref_nt = ref_seq[ref_pos]                   # reference nucleotide
            table_key = ref_nt + "I"                    # generate confusion table key
            conf_table[table_key] += op_length          # incrementing insertion counter
            insertions.extend([ref_pos] * op_length)    # storing insertion positions
            read_pos += op_length                       # update read position (ref position does not change)

        elif operation == 2:            # CASE: deletion
            ref_nt = ref_seq[ref_pos]                   # reference nucleotide
            table_key = ref_nt + "D"                    # generate confusion table key
            conf_table[table_key] += op_length          # incrementing deletion counter
            deletions.extend([ref_pos] * op_length)     # storing deletion positions
            ref_pos += op_length                        # update reference position (read pos. does not change)

        else:                           # CASE: match (real match or mismatch)
            for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):
                if (snp_key_i + str(ref_start + 1 + refp)) in snp_ai_hash:
                    continue                            # skipping SNP / A-to-I sites
                if ref_seq[refp] == read_seq[readp]:    # real match
                    table_key = ref_seq[refp] * 2           # generate confusion table key
                    conf_table[table_key] += 1              # incrementing match counter
                else:                                   # mismatch
                    ref_nt = ref_seq[refp]                  # reference nucleotide
                    table_key = ref_nt + read_seq[readp]    # generate confusion table key
                    conf_table[table_key] += 1              # incrementing mismatch counter

            read_pos += op_length   # update read position
            ref_pos += op_length    # update reference position

    return([insertions, deletions, conf_table])     # returning

def construct_snp_table(reads, read_snpinfo, region_insinfo, minpos, maxpos):
    """
    :param reads:               list of read IDs to be processed
    :param read_snpinfo:        hash storing read_ID -> [starting position (0-based, incl), SNPs per position]
    :param region_insinfo:      hash storing reference_position -> {read_ids with ins at that position}
    :param minpos:              minimum reference position covered by a read
    :param maxpos:              maximum reference position covered by a read

    :return:    a detailed SNP table; nested list, first index corresponding to rows and second to columns:
                    row 1: genomic region's reference sequence positions (0-based), other rows: reads,
                    column 1: read names, other columns: reference sequence positions,
                    entries: None if the position is not covered by the read, 0 if read has no SNP at that position,
                        SNP nucleotide otherwise

    This function computes a detailed SNP table for all reads of a certain genomic region.
    """

    snptable = []       # initializing SNP table

    region_length = maxpos - minpos                                     # genomic region total length
    inspos = list(region_insinfo.keys())                                # genomic region insertion positions
    inspos.sort()                                                       # sorting insertion positions
    snptable.append([None] + list(range(minpos, maxpos)) + inspos)      # SNP table's first row

    for r in reads:     # iterating through genomic region's reads
        r_start = read_snpinfo[r][0]                            # read starting position
        r_entry = [r] + [None] * (r_start - minpos)             # uncovered leading positions are set to None
        r_entry += read_snpinfo[r][1:]                          # adding read's SNP information
        r_entry += [None] * (region_length - len(r_entry) + 1)  # uncovered trailing positions are set to None
        for insp in inspos:                                     # finally adding insertions, too
            if r in region_insinfo[insp]:
                r_entry.append(1)
            else:
                r_entry.append(0)
        snptable.append(r_entry)                                # adding read's entry to SNP table

    return(snptable)    # returning SNP table

def construct_conv_table(reads, read_convinfo, region_convpos):
    """
    :param reads:               list of read IDs to be processed
    :param read_convinfo:       hash storing read_ID -> [starting position (0-based, incl), conv per position]
    :param region_convpos:      list of potential conversion positions (0-based)

    :return:    a detailed conversion table; nested list, first index corresponding to rows and second to columns):
                    row 1: potential reference sequence conversion position, other rows: reads,
                    column 1: read names, other columns: potential reference sequence conversion positions,
                    entries: None if the position is not covered, masked by a deletion or SNP correction, 0 if read
                        has no conversion at that position (may have any other SNP, though), 1 otherwise

    This function computes a detailed conversion table for all reads of a certain genomic region.
    """

    convtable = []      # initializing conversion table

    region_convpos.sort()                           # sorting potential conversion positions
    conv_length = len(region_convpos)               # number of potential conversion positions
    convtable.append([None] + region_convpos)       # conversion table's first row

    for r in reads:     # iterating through genomic region's reads
        r_start = read_convinfo[r][0]                           # read starting positions
        r_entry = [r]                                           # read ID
        for p in region_convpos:                                # uncovered leading positions are set to None
            if p < r_start:
                r_entry.append(None)
            else:
                break
        r_entry += read_convinfo[r][1:]                         # adding read's conversion information
        r_entry += [None] * (conv_length - len(r_entry) + 1)    # uncovered trailing positions are set to None
        convtable.append(r_entry)                               # adding read's entry to conversion table

    return (convtable)  # returning conversion table

def construct_snp_conv_summary_table(mode, snptable, convtable):
    """
    :param mode:            mode to be executed; choose 'snp' for computing a SNP table summary, choose 'conv' for
                            computing a conversion table summary, choose 'both' for computing both summary tables
    :param snptable:        detailed SNP table; nested list, first index corresponding to rows and second to columns:
                                row 1: genomic region's reference sequence positions (0-based), other rows: reads,
                                column 1: read names, other columns: reference sequence positions,
                                entries: None if the position is not covered by the read, 0 if read has no SNP at that
                                    position, SNP nucleotide otherwise
    :param convtable:       detailed conversion table; nested list, first index corresponding to rows and second to
                            columns):
                                row 1: potential reference sequence conversion position, other rows: reads,
                                column 1: read names, other columns: potential reference sequence conversion positions,
                                entries: None if the position is not covered, masked by a deletion or SNP correction,
                                    0 if read has no conversion at that position (may have any other SNP, though),
                                    1 otherwise

    :return:    a list containing the SNP and conversion summary tables; summary tables are nested lists, with first
                index corresponding to rows and second index to columns; first row stores positions (first entry: None),
                second row stores total counts (first entry: 'total_counts'), third row stores event counts (first
                entry: 'event_counts')

    This function computes summary tables of detailed SNP and conversion tables.
    """

    snptable_summary = []
    convtable_summary = []

    if mode in ("snp", "both"):     # SNP summary

        total_counts = [0 for i in snptable[0][1:]]     # initializing list storing total counts
        event_counts = [0 for i in snptable[0][1:]]     # initializing list storing event counts
        for read in snptable[1:]:                       # iterating through reads
            for p, entry in enumerate(read[1:]):            # iterating through read positions
                if entry not in (None, "I"):                    # counting up total counts
                    total_counts[p] += 1
                    if entry != 0:                              # counting up event counts
                        event_counts[p] += 1
        snptable_summary = [snptable[0], ["total_counts"] + total_counts, ["event_counts"] + event_counts]

    if mode in ("conv", "both"):    # conversion summary

        total_counts = [0 for i in convtable[0][1:]]    # initializing list storing total counts
        event_counts = [0 for i in convtable[0][1:]]    # initializing list storing event counts
        for read in convtable[1:]:                      # iterating through reads
            for p, entry in enumerate(read[1:]):            # iterating through read positions
                if entry != None:                               # counting up total counts
                    total_counts[p] += 1
                    if entry == 1:                              # counting up event counts
                        event_counts[p] += 1
        convtable_summary = [convtable[0], ["total_counts"] + total_counts, ["event_counts"] + event_counts]

    return([snptable_summary, convtable_summary])

def construct_new_ce_estimators_table(stable_file, ce_original_file, ce_var_file, ce_const_file, outfile, fn, fp):
    """
    :param stable_file:         input summary table file
    :param ce_original_file:    input conversion efficiency table file (of original efficiency estimation)
    :param ce_var_file:         input conversion efficiency table file (using time-dependent conversion efficiency)
    :param ce_const_file:       input conversion efficiency table file (using time-constant conversion efficiency)
    :param outfile:             output file for overview table
    :param fn:                  false negative rate
    :param fp:                  false positive rate

    :return: void

    This function contains an overview table of conversion efficiency estimated and newly synthesized estimates, with
    columns as the following:
    1) gene name
    2) labeled/total ratio measured
    3) new/total as inferred from labeled/total with time-dependent conversion efficiency
    4) new/total as inferred from labeled/total with time-constant conversion efficiency
    5) new/total as originally estimated
    6) new/total as back-estimated with time-dependent conversion efficiency
    7) new/total as back-estimated with time-constant conversion efficiency
    8) normalized difference of n/t inferred time-dependent and originally estimated
    9) normalized difference of n/t inferred time-constant and originally estimated
    10) normalized difference of n/t time-dependent and originally estimated
    11) normalized difference of n/t time-constant and originally estimated
    12) normalized difference of n/t time-dependent inferred and estimated
    13) normalized difference of n/t time-constant inferred and estimated
    14) quotient of n/t inferred time-dependent and originally estimated
    15) quotient difference of n/t inferred time-constant and originally estimated
    16) quotient difference of n/t time-dependent and originally estimated
    17) quotient difference of n/t time-constant and originally estimated
    18) quotient difference of n/t time-dependent inferred and estimated
    19) quotient difference of n/t time-constant inferred and estimated
    20) total read counts
    21) average conversion positions
    22) conversion efficiency originally estimated

    This function creates an overview table for the comparison of conversion efficiency and new/total ratio estimates
    under different conditions.
    """

    # opening files, loading in data

    out_file = open(outfile, "w")
    stab = load_summary_table(stable_file)
    ce_o = {}
    ce_v = {}
    ce_c = {}

    ce_intab = load_conveff_table(ce_original_file)
    for g in ce_intab:
        if ce_intab[g][0] == None: continue
        ce_o[g] = [ce_intab[g][0][-1], ce_intab[g][1][-1]]
    ce_intab = load_conveff_table(ce_var_file)
    for g in ce_intab:
        if ce_intab[g][0] == None: continue
        ce_v[g] = [ce_intab[g][0][-1], ce_intab[g][1][-1]]
    ce_intab = load_conveff_table(ce_const_file)
    for g in ce_intab:
        if ce_intab[g][0] == None: continue
        ce_c[g] = [ce_intab[g][0][-1], ce_intab[g][1][-1]]
    del ce_intab

    # writing output file header

    out_file.write("gene_name\tlab_ratio\tnew_inf_var\tnew_inf_const\tnew_est_original\tnew_est_var\tnew_est_const\t"
                   "nd_inf_var_est_or\tnd_inf_const_est_or\tnd_est_var_or\tnd_est_const_or\t"
                   "nd_est_inf_var\tnd_est_inf_const\t"
                   "qo_inf_var_est_or\tqo_inf_const_est_or\tqo_est_var_or\tqo_est_const_or\t"
                   "qo_est_inf_var\tqo_est_inf_const\t"
                   "tot\tcp\n")

    # iterating through genes estimated, recording metadata, estimators and estimator differences

    for g in ce_v:

        # measured / estimated values

        for d in stab[g]:
            tot = stab[g][d][1]
            lab = stab[g][d][2]
            cp = stab[g][d][3]

        ce_original = ce_o[g][1]
        new_original = ce_o[g][0]
        new_est_var = ce_v[g][0]
        new_est_const = ce_c[g][0]
        lab_ratio = lab / tot

        # inferred new/total ratios

        ce_v_value = ce_v[g][1]
        ce_c_value = ce_c[g][1]
        p_mod_v = 1 - (1 - ce_v_value * (1 - fn)) ** cp
        p_mod_c = 1 - (1 - ce_c_value * (1 - fn)) ** cp
        p_un = 1 - (1 - fp) ** cp
        new_inf_var = (lab_ratio + p_un - 1) / (p_mod_v + p_un - 1)
        new_inf_const = (lab_ratio + p_un - 1) / (p_mod_c + p_un - 1)

        # normalized differences and quotiens

        nd_inf_var_est_original = abs((new_inf_var - new_original) / new_original)
        nd_inf_const_est_original = abs((new_inf_const - new_original) / new_original)
        nd_est_var_original = abs((new_est_var - new_original) / new_original)
        nd_est_const_original = abs((new_est_const - new_original) / new_original)
        nd_est_inf_var = abs((new_est_var - new_inf_var) / new_inf_var)
        nd_est_inf_const = abs((new_est_const - new_inf_const) / new_inf_const)

        qo_inf_var_est_original = new_inf_var / new_original
        qo_inf_const_est_original = new_inf_const / new_original
        qo_est_var_original = new_est_var / new_original
        qo_est_const_original = new_est_const / new_original
        qo_est_inf_var = new_est_var / new_inf_var
        qo_est_inf_const = new_est_const / new_inf_const

        # writing to output file

        out_file.write(g + "\t" + str(lab_ratio) + "\t" + str(new_inf_var) + "\t" + str(new_inf_const) + "\t" +
                       str(new_original) + "\t" + str(new_est_var) + "\t" + str(new_est_const) + "\t" +
                       str(nd_inf_var_est_original) + "\t" + str(nd_inf_const_est_original) + "\t" +
                       str(nd_est_var_original) + "\t" + str(nd_est_const_original) + "\t" +
                       str(nd_est_inf_var) + "\t" + str(nd_est_inf_const) + "\t" +
                       str(qo_inf_var_est_original) + "\t" + str(qo_inf_const_est_original) + "\t" +
                       str(qo_est_var_original) + "\t" + str(qo_est_const_original) + "\t" +
                       str(qo_est_inf_var) + "\t" + str(qo_est_inf_const) + "\t" +
                       str(tot) + "\t" + str(cp) + "\t" + str(ce_original) + "\n")

    # closing output file and returning

    out_file.close()
    return ()

# Data Storing and Loading #############################################################################################

# --- additional information ---

def load_snps(snpfile, chr_selected):
    """
    :param snpfile:         gzipped vcf file containing SNP positions of the reference
    :param chr_selected:    chromosome number for which SNPs are to be loaded

    :return: a hash storing single-nucleotide exchange SNPs as 'ChrNumber_ChrPos'

    This function loads single-nucleotide exchange SNPs from a gzipped vcf file for a selected chromosome number.
    """

    snps = dict()   # hash storing ChrNumber_ChrPos of SNP sites

    snp_file = gzip.open(snpfile, "rt")     # opening SNP file
    while next(snp_file).startswith("##"):  # skipping meta information (and header with last iteration)
        pass

    for line in snp_file:                   # iterating through SNP table

        fields = line.strip().split("\t")       # getting single vcf fields
        ref_chr = fields[0]                     # reference chromosome number
        if ref_chr != chr_selected:             # (skipping SNPs not from the current reference sequence)
            continue

        ref_pos = int(fields[1])                # reference chromosome position (1-based)
        ref_nt = fields[3]                      # reference nucleotides
        snp_nt = fields[4]                      # read nucleotides

        if (len(snp_nt) > 1) or (len(ref_nt) > 1):      # only considering single nucleotide exchanges
            continue
        else:
            snps[ref_chr + "_" + str(ref_pos)] = 0  # storing SNP site

    return(snps)    # returning snp hash

def load_ai_sites(aifile, chr_selected):
    """
    :param aifile:          REDIportal file containing A-to-I editing positions of the reference
    :param chr_selected:    chromosome number for which A-to-I sites are to be loaded

    :return: a hash storing A-to-I sites as 'ChrNumber_ChrPos'

    This function loads A-to-I sites from a REDIportal file for a selected chromosome number.
    """

    ai_file = open(aifile, "r")     # opening file storing A-to-I sites
    ai_hash = {}                    # hash storing ChrNumber_ChrPos of SNP sites

    for line in ai_file:            # iterating through A-to-I file
        fields = line.strip().split("\t")           # getting single A-to-I site fields
        ref_chr = fields[0][3:]                     # reference sequence name
        if ref_chr != chr_selected: continue        # (skipping AIs not from the selected reference sequence)

        ref_pos = fields[1]                         # reference chromosome position (1-based)
        ai_hash[ref_chr + "_" + str(ref_pos)] = 0   # storing A-to-I site

    return(ai_hash)     # returning A-to-I sites hash

def write_edit_sites(rtable_file, editing_cutoff, outfile, chr_numbers=chr_numbers_hs):
    """
    :param rtable_file:     read table file
    :param editing_cutoff:  maximum editing rate (SNP rate) of a nucleotide position to be tolerated
    :param outfile:         output file to write editing sites to

    :param chr_numbers:     hash storing valid chromosome numbers (default: hs chromosome numbers)

    :return: void

    This function detects editing sites, computing position-wise SNP rates from a read table file and storing all
    positions with a SNP rate above a certain threshold. Editing positions are written to an output file as
    'ChrNumber_ChrPos\n'
    """

    rtable = load_read_table(rtable_file)       # loading read table
    out_file = open(outfile, "w")               # opening output file

    for gr in rtable:                           # iterating through genomic regions

        gr_name_fields = gr.split("_")              # splitting genomic region name into single parts
        for field in gr_name_fields:                # getting genomic region's corresponding chromosome number
            if field.startswith("chr"):
                gr_chr = field[3:]
        if gr_chr not in chr_numbers: continue      # skipping genomic regions not matching a valid chromosome number

        gr_snp_table = rtable[gr][0]                    # genomic region's SNP table

        for i, pos in enumerate(gr_snp_table[0][1:]):   # iterating through genomic region positions

            total = gr_snp_table[1][i + 1]                          # getting total counts
            if total == 0: continue                                 # (skipping uncovered positions)
            events = gr_snp_table[2][i + 1]                         # getting event (SNP) counts
            if (events / total) > editing_cutoff:                   # computing SNP rate
                out_file.write(gr_chr + "_" + str(pos + 1) + "\n")  # storing positions with SNP rate above threshold

    out_file.close()    # closing output file
    return()            # returning

def load_edit_sites(editsites_file, chr_selected):
    """
    :param editsites_file:  file storing editing sites as 'ChrNumber_ChrPos\n'
    :param chr_selected:    chromosome number for which editing sites are to be loaded

    :return: a hash storing editing sites as 'ChrNumber_ChrPos'

    This function loads editing sites from a editing sites file for a selected chromosome number.
    """

    edit_sites_file = open(editsites_file, "r")     # opening editing sites file
    edit_hash = {}                                  # hash storing editing sites as 'ChrNumber_ChrPos'

    for line in edit_sites_file:                    # iterating through editing sites file
        if line.split("_")[0] == chr_selected:          # checking if editing site belongs to selected chromosome
            edit_hash[line.strip()] = 0                     # if so, storing editing site to return hash

    edit_sites_file.close()     # closing input editing sites file
    return(edit_hash)           # returning

def load_gff3_annotation(gff3_file):
    """
    :param gff3_file:   input gzipped GFF3 file (gff.gz)

    :return:    a list containing (in the following order):
                features' hash:
                    feature_ID -> [[[ref_chr, start, end, strand], ...], feature_type, name, database_crossref]
                    starting positions 1-based, including
                    ending positions 1-based, excluding
                parent-child relation hash:
                    feature_ID -> [parent_id, {child_ids}]
                recorded feature types' hash:
                    feature_type -> number of entries with that feature type

    This function loads an annotation in GFF3 file format, storing all feature entries to a hash. Parent-child relations
    amongst the features are written to a second hash. Furthermore, all feature types recorded within the GFF file are
    stored in a separate hash map.
    """

    gff_infile = gzip.open(gff3_file, "rt")     # opening gff.gz file
    feature_hash = {}                           # initialize features' hash
    parent_child_hash = {}                      # initialize parent-child relations' hash
    feature_type_hash = {}                      # initialize feature types' hashmap

    for line in gff_infile:         # iterating through GFF3 file

        if line.startswith("#"): continue   # skipping meta information
        fields = line.strip().split("\t")   # getting entry's single fields
        if fields[2] == "region": continue  # skipping chromosome / scaffold entries

        # --- getting location information ---

        ref_chr = fields[0]                 # reference chromosome or scaffold
        start = int(fields[3])              # starting position (1-based, including)
        end = int(fields[4])                # ending position (1-based, excluding)
        strand = fields[6]                  # strand orientation

        # --- getting attribute information ---

        type = fields[2]                    # feature type
        attributes = fields[8].split(";")   # further attributes
        id = ""                             # initialize feature ID
        parent = ""                         # initialize parent's ID
        name = ""                           # initialize name
        dbxref = ""                         # initialize database cross references
        for a in attributes:
            if a.startswith("ID="):
                id = a[3:]
            elif a.startswith("Parent="):
                parent = a[7:]
            elif a.startswith("Name="):
                name = a[5:]
            elif a.startswith("Dbxref="):
                dbxref = a[7:]

        # --- filling hashes ---

        # filling feature hash

        if id not in feature_hash:          # feature has no hash entry: build new entry
            feature_hash[id] = [[[ref_chr, start, end, strand]], type, name, dbxref]
        else:                               # feature has a hash entry: adding location information
            feature_hash[id][0].append([ref_chr, start, end, strand])

        # filling parent-child hash

        if id not in parent_child_hash:     # feature has no hash entry: build new entry
            parent_child_hash[id] = [parent, {}]
        else:                               # feature has a hash entry (appeared as parent before): add parent
            parent_child_hash[id][0] = parent

        if parent:                          # feature has a parent
            if parent not in parent_child_hash:     # parent has no hash entry: build new entry
                parent_child_hash[parent] = ["", {id: None}]
            else:                                   # parent has a hash entry: add feature id as child
                parent_child_hash[parent][1][id] = None

        # filling feature type hash

        if type not in feature_type_hash:
            feature_type_hash[type] = 1
        else:
            feature_type_hash[type] += 1

    # closing GFF3 file and returning

    gff_infile.close()
    return([feature_hash, parent_child_hash, feature_type_hash])

def load_bed_regions(bedfile):
    """
    :param bedfile:     input BED file for which genomic regions are to be loaded

    :return:    a list of two hashes (for (-) and (+) strand BED regions, respectively) storing BED regions as
                {ref_name -> [[start_1, end_1, region_id_1], [start_2, end_2, region_id_2], [...], ...]}
                (starting position 0 based, including, ending position 0-based, excluding)

    This function loads BED regions, separating them by strand orientation and reference sequence.
    """

    bed_file = open(bedfile)    # opening BED file
    regions_minus = {}          # initializing hash storing (-) BED regions
    regions_plus = {}           # initializing hash storing (+) BED regions

    for line in bed_file:       # iterating though BED file entries

        fields = line.strip().split("\t")   # getting BED file entry's single fields

        name = fields[3]                    # region identifier
        strand = fields[5]                  # region strand orientation
        chr_number = fields[0][3:]          # genomic region chromosome number
        startpos = int(fields[1])           # genomic region starting position, 0-based, including
        endpos = int(fields[2])             # genomic region ending position, 0-based, excluding

        if strand == "-":                   # (-) hash to store region
            selected_hash = regions_minus
        elif strand == "+":                 # (+) hash to store region
            selected_hash = regions_plus
        else:                               # skipping regions with unspecified strand orientation
            continue

        if chr_number not in selected_hash:     # initializing chromosome number's entries
            selected_hash[chr_number] = [[startpos, endpos, name]]
        else:                                   # appending BED region to chromosome number's entries
            selected_hash[chr_number].append([startpos, endpos, name])

    bed_file.close()                        # closing BED file
    return([regions_minus, regions_plus])   # returning

def bed_index(bed_regions, strand, ref_name):
    """
    :param bed_regions:     list of two hashes (for (-) and (+) strand BED regions, respectively) storing BED regions
                            as {ref_name -> [[start_1, end_1, region_id_1], [start_2, end_2, region_id_2], [...], ...]}
                            (starting position 0 based, including, ending position 0-based, excluding)
    :param strand:          strand for which the index should be built (either "-" or "+")
    :param ref_name:        reference name for which the index should be built

    :return:    a hash storing the BED region positions as position -> [region_ids] (positions 0-based)

    This function builds a BEd region index. In detail, for a given strand orientation and reference sequence, a hash
    is created which assigns all positions defined in the BED regions a list with the corresponding region IDs.
    """

    bed_index_hash = {}         # initializing BED index hash
    bed_regions_list = []       # initializing list of BED regions

    if strand == "-":           # CASE: selecting (-) strand region
        if ref_name in bed_regions[0]:  # checking if reference sequence has BED regions on specified strand
            bed_regions_list = bed_regions[0][ref_name]     # selecting those regions
    elif strand == "+":         # CASE: selecting (+) strand regions
        if ref_name in bed_regions[1]:  # checking if reference sequence has BED regions on specified strand
            bed_regions_list = bed_regions[1][ref_name]     # selecting those regions

    for startpos, endpos, id in bed_regions_list:   # iterating through BED regions
        for p in range(startpos, endpos):               # iterating through BED region positions
            if p not in bed_index_hash:                     # storing new BED region position and region ID
                bed_index_hash[p] = [id]
            else:                                           # storing existing BED region position and region ID
                bed_index_hash[p].append(id)

    return(bed_index_hash)      # returning

def write_bed_annotation(bed_annotation_hash, outfile):
    """
    :param bed_annotation_hash:     annotation hash to be written to an output file
                                    a hash assigning each BED genomic region ID a list with the region's annotation as:
                                    {region_ID -> [ [BAM1_descriptor, BAM1_counts, BAM2_descriptor, BAM2_counts, ...],
                                                    [exon_id], [intron_id], [3UTR_id], [5UTR_id],
                                                    {feature_ID -> [    type, name, db_crossref,
                                                                        {sub_feature_ID -> [type, name, db_crossref]}
                                                    ]}
                                    ]}
    :param outfile:                 output file to write annotations to

    :return: void

    ....................................................................................................................

    This function writes annotations to an output file. The file starts with one header-line (marked by a leading '#'
    symbol) and is tab-separated. Columns are in the following order:

    > region identifier
    > first BAM file descriptor
    ...
    > n-th BAM file descriptor
    > exonic region IDs
    > intronic region IDs
    > 3'UTR IDs
    > 5'UTR IDs
    > primary features' IDs
    > primary features' types
    > primary features' names
    > primary features' database cross-references
    > secondary features' IDs
    > secondary features' types
    > secondary features' names
    > secondary features' database cross references

    If annotation information are not available for a region, columns are filled with '.'.
    Within one column ...
    ... multiple exonic, intronic, 3'UTR or 5'UTR IDs are ';' separated (example: intronID1;intronID2)
    ... multiple primary feature information are ';' separated (example: ID1;ID2)
    ... multiple secondary feature information are '*' separated (example: ID1_1*ID1_2;ID2_1*ID2_2*ID2_3)
    """

    out_file = open(outfile, "w")                   # opening output file

    n_bam_counts = len(
        bed_annotation_hash[next(iter(bed_annotation_hash))][0]
    ) / 2                                           # number of BAM files for which region read counts are recorded
    n_bam_counts = int(n_bam_counts)
    bam_des = [
        bed_annotation_hash[next(iter(bed_annotation_hash))][0][i]
        for i in range(0, n_bam_counts*2, 2)]       # BAM file descriptors

    c_idx = [n_bam_counts + i for i in range(12)]   # column indices for all non-read-count information

    out_file.write("#" + "\t".join(["region_id"] + bam_des) + "\texon_id\tintron_id\t3'UTR_id\t5'UTR_id\t" +
                   "primary_ids\tprimary_types\tprimary_names\tprimary_xref\t" +
                   "secondary_ids\tsecondary_types\tsecondary_names\tsecondary_xref\n")     # writing header line

    for r_id in bed_annotation_hash:  # iterating through region identifiers

        r_entry = bed_annotation_hash[r_id]             # region's annotation entry

        columns = [str(r_entry[0][i]) for i in range(1, n_bam_counts*2, 2)] + \
                  ["", "", "", "",
                   "", "", "", "",
                   "", "", "", ""]                      # initializing columns' entries r_entry[1:5] + \

        for i, ann_list in enumerate(r_entry[1:5]):     # iterating through exon, intron, 3'UTR and 5'UTR IDs
            columns[c_idx[i]] = ";".join(ann_list)

        for f in r_entry[5]:                            # iterating through primary features
            columns[c_idx[4]] += f + ";"                    # updating IDs
            columns[c_idx[5]] += r_entry[5][f][0] + ";"     # updating types
            columns[c_idx[6]] += r_entry[5][f][1] + ";"     # updating names
            columns[c_idx[7]] += r_entry[5][f][2] + ";"     # updating cross references
            f_children = r_entry[5][f][3]                               # secondary features' sub-hash
            child_ids = list(f_children.keys())                         # secondary features' IDs
            columns[c_idx[8]] += "*".join(child_ids) + ";"                              # updating IDs
            columns[c_idx[9]] += "*".join(f_children[c][0] for c in child_ids) + ";"    # updating types
            columns[c_idx[10]] += "*".join(f_children[c][1] for c in child_ids) + ";"   # updating names
            columns[c_idx[11]] += "*".join(f_children[c][2] for c in child_ids) + ";"   # updating cross refs

        # writing to output file (removing trailing semi-colons)
        columns[c_idx[4]:] = \
            [i[:-1] for i in columns[n_bam_counts+4:]]              # cutting trailing ';'
        columns = [i if i != "" else "." for i in columns]          # filling missing info
        out_file.write(r_id + "\t" + "\t".join(columns) + "\n")     # writing column entries

    out_file.close()    # closing output file
    return()            # returning

def load_bed_annotation(bed_annotation_file):
    """
    :param bed_annotation_file:     input file storing annotations

    :return:    a hash assigning each BED genomic region ID a list with the region's annotation as:
                {region_ID -> [ [BAM1_descriptor, BAM1_counts, BAM2_descriptor, BAM2_counts, ...],
                                [exon_ids], [intron_ids], [3UTR_ids], [5UTR_ids],
                                {feature_ID -> [    type, name, db_crossref,
                                                    {sub_feature_ID -> [type, name, db_crossref]}
                                ]}
                ]}

    ....................................................................................................................

    This function loads annotations from a file. The file format should be as described in the following:

    Tab-separated, starting with one header-line marked by a leading '#' symbol.
    Columns in the following order:

    > region identifier
    > first BAM file descriptor
    ...
    > n-th BAM file descriptor
    > exonic region ID
    > intronic region ID
    > 3'UTR ID
    > 5'UTR ID
    > primary features' IDs
    > primary features' types
    > primary features' names
    > primary features' database cross-references
    > secondary features' IDs
    > secondary features' types
    > secondary features' names
    > secondary features' database cross references

    If annotation information are not available for a region, columns are filled with '.'.
    Within one column ...
    ... multiple exonic, intronic, 3'UTR or 5'UTR IDs are ';' separated (example: intronID1;intronID2)
    ... multiple primary feature information are ';' separated (example: ID1;ID2)
    ... multiple secondary feature information are '*' separated (example: ID1_1*ID1_2;ID2_1*ID2_2*ID2_3)
    """

    ann_hash = {}                                   # initializing annotation hash to be returned
    ann_file = open(bed_annotation_file)            # opening annotation file

    header = next(ann_file).strip().split("\t")                 # getting header line fields
    n_count_cols = len(header) - 13                             # number of columns containing read count data
    count_des = [header[i] for i in range(1, n_count_cols+1)]   # descriptors for BAM files of read count data
    c_idx = [i + n_count_cols for i in range(1, 13)]            # indices of non-read count columns

    for line in ann_file:       # iterating through annotation file

        # region ID and read counts ....................................................................................

        fields = line.strip().split("\t")                           # getting entry's fields
        fields = [i if i != "." else "" for i in fields]            # converting missing info "." to ""
        r_id = fields[0]                                            # getting region ID

        count_info = []                                             # initialize list with read count information
        for des, c in zip(count_des, fields[1:n_count_cols+1]):     # zipping BAM file descriptors and count data
            count_info += [des, int(c)]

        ann_hash[r_id] = [count_info, [], [], [], [], {}]           # initialize region's hash entry

        # exon, intron, 3'UTR and 5'UTR IDs ............................................................................

        for i, entry in enumerate(fields[c_idx[0]:c_idx[4]]):   # iterate through exon, intron, 3'UTR, 5'UTR entries
            if entry:                                               # check if IDs are recorded
                ids = entry.split(";")                                  # get single IDs
                ann_hash[r_id][i+1] = ids                               # store single IDs

        # primary and secondary feature information ....................................................................

        if fields[c_idx[4]] == "":              # no further annotation information recorded: continue
            continue

        prim_ids = fields[c_idx[4]].split(";")          # getting primary features' IDs
        prim_types = fields[c_idx[5]].split(";")        # primary features' types
        prim_names = fields[c_idx[6]].split(";")        # primary features' names
        prim_xrefs = fields[c_idx[7]].split(";")        # primary features' cross references
        prim_sec_ids = fields[c_idx[8]].split(";")      # primary features' secondary feature IDs
        prim_sec_types = fields[c_idx[9]].split(";")    # primary features' secondary feature tyes
        prim_sec_names = fields[c_idx[10]].split(";")   # primary features' secondary feature names
        prim_sec_xrefs = fields[c_idx[11]].split(";")   # primary features' secondary feature cross references

        for i, f_id in enumerate(prim_ids):     # iterating through primary features

            ann_hash[r_id][5][f_id] = \
                [prim_types[i], prim_names[i], prim_xrefs[i], {}]       # storing primary feature info
            if prim_sec_ids[i] == "":                                   # no secondary annotations recorded: continue
                continue

            sec_ids = prim_sec_ids[i].split("*")                        # secondary features' IDs
            sec_types = prim_sec_types[i].split("*")                    # secondary features' types
            sec_names = prim_sec_names[i].split("*")                    # secondary features' names
            sec_xrefs = prim_sec_xrefs[i].split("*")                    # secondary features' cross references

            for j, c_id in enumerate(sec_ids):      # iterating through secondary features
                ann_hash[r_id][5][f_id][3][c_id] = \
                    [sec_types[j], sec_names[j], sec_xrefs[j]]              # storing secondary feature info

    ann_file.close()        # closing annotation file
    return(ann_hash)        # returning annotation hash

def write_covered_positions(covpos, outfile):
    """
    :param covpos:      a list containing two hashes (one for forward and one for reversely mapped reads, respectively);
                        the hashes assign each reference sequence name another hash which assigns a position (1-based)
                        within the reference the coverage
    :param outfile:     output file to write covered positions to

    :return: void

    This function writes covered positions to an output file. Forwards and reversely mapped reads are listed separately
    (therefore, the same reference sequence name and nucleotide position can occur twice within the file). The output
    file is formatted as described in the following (mapping orientation is indicated by "f" for forwards and "r" for
    reversely mapped reads):

    #reference_sequence_name_1\n
    position_1\tcoverage\tmapping_orientation\n
    position_2\tcoverage\tmapping_orientation\n
    ...
    position_n\tcoverage\tmapping_orientation\n
    #reference_sequence_name_2\n
    ...
    """

    out_file = open(outfile, "w")       # opening output file
    covpos_forward = covpos[0]          # covered positions (forwards mapped reads)
    covpos_reverse = covpos[1]          # covered positions (reversely mapped reads)

    # writing covered positions of forwards mapped reads

    for refname in covpos_forward:          # iterate through reference sequence names
        out_file.write("#" + refname + "\n")        # writing reference sequence name
        refname_covpos = covpos_forward[refname]    # getting reference sequence's covered positions
        for cp in refname_covpos:                   # iterate through covered positions
            out_file.write(str(cp) + "\t" + str(refname_covpos[cp]) + "\t" + "f\n")     # writing covered position

    # writing covered positions of reversely mapped reads

    for refname in covpos_reverse:          # iterate through reference sequence names
        out_file.write("#" + refname + "\n")        # writing reference sequence name
        refname_covpos = covpos_reverse[refname]    # getting reference sequence's covered positions
        for cp in refname_covpos:                   # iterate through covered positions
            out_file.write(str(cp) + "\t" + str(refname_covpos[cp]) + "\t" + "r\n")     # writing covered position

    # closing output file and returning

    out_file.close()
    return()

def load_covered_positions(covpos_file):
    """
    :param covpos_file:     input file storing covered positions

    :return:    a list containing two hashes (one for forward and one for reversely mapped reads, respectively); the
                hashes assign each reference sequence name another hash which assigns a position (1-based) within the
                reference the coverage

    This function loads covered positions from a file. The input file should be formatted as following (listing forwards
    and reversely mapped reads separately, mapping orientation being indicated by "f" and "r", respectively):

    #reference_sequence_name_1\n
    position_1\tcoverage\tmapping_orientation\n
    position_2\tcoverage\tmapping_orientation\n
    ...
    position_n\tcoverage\tmapping_orientation\n
    #reference_sequence_name_2\n
    ...
    """

    covpos_forward = {}     # initialize covered positions of forwards mapped reads
    covpos_reverse = {}     # initialize covered positions of reversely mapped reads

    cp_file = open(covpos_file, "r")        # opening covered positions' file

    ref_name = next(cp_file).strip()[1:]    # getting first reference sequence name
    ref_forward = {}                        # initialize reference's forwards mapped reads covered positions hash
    ref_reverse = {}                        # initialize reference's reversely mapped reads covered positions hash

    for line in cp_file:                    # iterating through covered positions' file
        if line.startswith("#"):                # CASE: next reference sequence encountered
            if ref_forward: covpos_forward[ref_name] = ref_forward      # storing previous reference's forwards hash
            if ref_reverse: covpos_reverse[ref_name] = ref_reverse      # storing previous reference's reverse hash
            ref_name = line.strip()[1:]                                 # getting next reference sequence's name
            ref_forward = {}                                            # resetting forwards hash
            ref_reverse = {}                                            # resetting reverse hash
        else:                                   # CASE: covered position's information encountered
            cp_info = line.strip().split("\t")                  # getting covered position's information
            cp = int(cp_info[0])                                # covered position
            cov = int(cp_info[1])                               # coverage
            orientation = cp_info[2]                            # mapping orientation
            if orientation == "f": ref_forward[cp] = cov        # storing covered position to correct hash (according
            elif orientation == "r": ref_reverse[cp] = cov      # to mapping orientation)

    if ref_forward: covpos_forward[ref_name] = ref_forward      # storing last encountered reference sequence's ...
    if ref_reverse: covpos_reverse[ref_name] = ref_reverse      # ... covered positions

    cp_file.close()                             # closing input file
    return([covpos_forward, covpos_reverse])    # returning

def write_coverage_profile(cov_profile, outfile):
    """
    :param cov_profile:     a list containing two hashes (for forward and reversely mapped reads, respectively); the
                            hashes assign each reference sequence name a list which stores the coverage profile; the
                            smallest list index (0) represents the downstream-most nucleotide position, the largest
                            list index (radius * 2 + 1) the upstream-most nucleotide and the middle list index (radius)
                            the central covered nucleotide position
    :param outfile:         output file to write coverage profile to

    :return: void

    This function writes coverage profiles to an output file in tabular format. Profile window positions are stored as
    rows (the central position being 0) and profiles are stored as columns. The first row of the table contains the
    column names, which are the reference sequence names and mapped reads' orientation of the profiles: 'refname_f' for
    forwards and 'refname_r' for reversely mapped reads. The following rows of the table contain the respective window
    position as rowname and the profile count values as following entries.
    """

    out_file = open(outfile, "w")           # opening output file
    profiles_forwards = cov_profile[0]      # forwards mapped reads' coverage profile
    profiles_reverse = cov_profile[1]       # reversely mapped reads' coverage profile
    refnames_forward = list(profiles_forwards.keys())   # reference names in forwards hash
    refnames_reverse = list(profiles_reverse.keys())    # reference names in reverse hash

    window_size = len(profiles_forwards[next(iter(profiles_forwards))])     # coverage profile window size
    radius = int((window_size - 1) / 2)                                     # coverage profile radius

    out_file.write("\t".join([r + "_f" for r in refnames_forward] +
                             [r + "_r" for r in refnames_reverse]) + "\n")  # writing output file column names

    # writing window positions' counts
    for i in range(0, window_size):         # iterating through window positions
        out_file.write(str(i - radius))         # writing window position

        for ref_name in refnames_forward:                               # iterating through forwards references
            out_file.write("\t" + str(profiles_forwards[ref_name][i]))      # writing reference count
        for ref_name in refnames_reverse:                               # iterating through reverse references
            out_file.write("\t" + str(profiles_reverse[ref_name][i]))       # writing reference count
        out_file.write("\n")                                            # ending row with newline

    out_file.close()    # closing output file
    return()            # returning

def load_coverage_profile(profile_file):
    """
    :param profile_file:    input file storing coverage profile

    :return:    a list containing two hashes (for forward and reversely mapped reads, respectively); the hashes assign
                each reference sequence name a list which stores the coverage profile; the smallest list index (0)
                represents the downstream-most nucleotide position, the largest list index (radius * 2 + 1) the
                upstream-most nucleotide and the middle list index (radius) the central covered nucleotide position

    This function loads a coverage profile from a tabular input file. The input file should be formatted as follows:
    Profile window positions are stored as rows (the central position being 0) and profiles are stored as columns.
    The first row of the table contains the column names, which are the reference sequence names and mapped reads'
    orientation of the profiles: 'refname_f' for forwards and 'refname_r' for reversely mapped reads. The following rows
    of the table contain the respective window position as rowname and the profile count values as following entries.
    """

    # reading out coverage profile reference names and window size .....................................................

    in_file = open(profile_file, "r")               # opening coverage profile file
    window_size = -1                                # initialize window size (-1 for line with column names)
    for line in in_file: window_size += 1           # counting the window size
    in_file.close()                                 # closing coverage profile file again

    # loading coverage profiles ........................................................................................

    # initializing

    in_file = open(profile_file, "r")                       # opening coverage profile file
    ref_names = next(in_file).strip().split("\t")           # getting reference names and mapped reads' orientations
    all_profiles = [[0 for pos in range(0, window_size)]
                    for r in ref_names]                     # initializing empty profiles (list of profile lists)
    forwards_profiles = {}                                  # initialize hash storing forwards profiles
    reverse_profiles = {}                                   # initialize hash storing reverse profiles

    for i, r in enumerate(ref_names):       # iterating through profiles' ref. names and orientations
        rname = r[:-2]                          # getting name
        ori = r[-1]                             # getting orientation
        if ori == "f":                          # CASE: profile is forward
            forwards_profiles[rname] = all_profiles[i]  # forwards hash stores corresponding profile
        elif ori == "r":                        # CASE: profile is reverse
            reverse_profiles[rname] = all_profiles[i]   # reverse hash stores corresponding profile

    for i, line in enumerate(in_file):      # iterating through window positions (input file rows)
        count_values = line.strip().split("\t")[1:]     # getting profiles' count values
        for j, prof in enumerate(all_profiles):         # iterating through profiles
            prof[i] = int(count_values[j])                  # storing count value to profile

    # closing input file and returning .................................................................................

    in_file.close()
    return([forwards_profiles, reverse_profiles])

def write_coverage_convolution(cov_convolution, outfile, append=False):
    """
    :param cov_convolution:     a list storing two sublists, which again store two hashes:
                                sublist 1 stores convolution value hashes, and sublist 2 stores the corresponding
                                    center position hashes
                                hash 1 refers to forwards mapped and hash 2 refers to reversely mapped reads
                                a hash assigns a reference sequence name its list of convolution values or corresponding
                                    center position values, respectively
    :param outfile:     output file to store convolution values to
    :param append:      set to TRUE to append convolution results to the output file (default: FALSE)

    :return: void

    This function stores coverage convolution values and their corresponding center positions to a given output file.
    The output file contains one tab-separated table for each reference sequence and each read orientation. Each table
    is preceded by a header line which indicates the reference name and read orientation as 'refname_f' or 'refname_r'
    for forwards and reversely mapped reads, respectively. The table itself stores center positions in the first column
    and convolution values in the second column.
    """

    mode = "w"                                      # default output file opening mode (writing)
    if append: mode = "wa"                          # 'append' output file opening mode (write and append)
    out_file = open(outfile, mode)                  # opening output file
    convolution_hashes = cov_convolution[0]         # sublist of coverage convolution hashes
    center_hashes = cov_convolution[1]              # sublist of center position hashes

    for i, o in zip([0, 1], ["f", "r"]):    # iterating through read orientation

        convolution_hash = convolution_hashes[i]        # convolution values' hash
        center_hash = center_hashes[i]                  # center position's hash

        for ref in convolution_hash:        # iterating through reference sequences

            out_file.write("#" + ref + "_" + o + "\n")      # writing reference header
            ref_convolution = convolution_hash[ref]         # reference's convolution values
            ref_centers = center_hash[ref]                  # reference's center positions

            for c, v in zip(ref_centers, ref_convolution):  # iterating through center positions and convolution values

                out_file.write(str(c) + "\t" + str(v) + "\n")   # writing center position and convolution value

    out_file.close()        # closing output file
    return()                # returning

def load_coverage_convolution(covconv_file):
    """
    :param covconv_file:    input file storing coverage convolution values

    :return:    a list storing two sublists, which again store two hashes:
                sublist 1 stores convolution value hashes, and sublist 2 stores the corresponding center position hashes
                hash 1 refers to forwards mapped and hash 2 refers to reversely mapped reads
                a hash assigns a reference sequence name its list of convolution values or corresponding center position
                    values, respectively

    This function loads coverage convolution values and their corresponding center positions from a given input file.
    The input file should be formatted as described in the following:
    The file contains one tab-separated table for each reference sequence and each read orientation. Each table is
    preceded by a header line which indicates the reference name and read orientation as 'refname_f' or 'refname_r'
    for forwards and reversely mapped reads, respectively. The table itself stores center positions in the first column
    and convolution values in the second column.
    """

    in_file = open(covconv_file, "r")       # opening input file
    convolution_values = [{}, {}]           # list storing convolution value hashes (for forwards and reverse reads)
    center_positions = [{}, {}]             # list storing center positions (for forwards and reverse mapped reads)

    ref_name = None             # initialize current reference sequence name
    ref_cp = []                 # initialize current reference sequence's center positions
    ref_cv = []                 # initialize current reference sequence's convolution values
    current_cp_hash = None      # initalize current center positions hash
    current_cv_hash = None      # initialize current convolution value hash

    for line in in_file:    # iterating through input file

        if line.startswith("#"):    # CASE: next reference sequence table encountered
            if ref_cv:                  # storing convolution values and center positions for the current reference
                current_cp_hash[ref_name] = ref_cp
                current_cv_hash[ref_name] = ref_cv
            ref_cp = []                                         # re-set current reference's center positions
            ref_cv = []                                         # re-set current reference's convolution values
            ref_name = line.strip()[1:-2]                       # getting next reference name
            orientation = line.strip()[-1]                      # getting next read orientation
            if orientation == "f":                              # read orientation forward: select respective hashes
                current_cp_hash = center_positions[0]
                current_cv_hash = convolution_values[0]
            elif orientation == "r":                            # read orientation reverse: select respective hashes
                current_cp_hash = center_positions[1]
                current_cv_hash = convolution_values[1]

        else:                       # CASE: table entry
            (cp, cv) = line.strip().split("\t")     # getting table entry's center position and convolution value
            ref_cp.append(int(cp))                  # appending center position to list
            ref_cv.append(float(cv))                # appending convolution value to list

    if ref_cv:  # storing last reference sequence's convolution values and center positions
        current_cp_hash[ref_name] = ref_cp
        current_cv_hash[ref_name] = ref_cv

    in_file.close()                                 # closing input file
    return([convolution_values, center_positions])  # returning

# --- derived data tables ---

def write_read_table(r_table, outfile):
    """
    :param r_table:     a hash assigning each genomic region's name a list containing a SNP and a conversion table
                        (implemented as nested lists, the first index corresponding to rows and the second to columns):
                            SNP table: row 1: genomic region's reference sequence positions (0-based), other rows:
                                reads, column 1: read names, other columns: reference sequence positions, entries: None
                                if the position is not covered by the read, 0 if read has no SNP at that position,
                                SNP nucleotide otherwise
                            conversion table: row 1: potential reference sequence conversion position, other rows:
                                reads, column 1: read names, other columns: potential reference sequence conversion
                                positions, entries: None if the position is not covered, masked by a deletion, a SNP or
                                SNP correction, 0 if read has no conversion at that position, 1 otherwise
    :param outfile:     output file to store read tables to

    :return: void

    ....................................................................................................................

    This function writes the contents of a read table hash to a file. The file is formatted as described in the
    following:

    #genomic_region_name\n
    >SNP_TABLE\n
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    >CONV_TABLE\n
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    """

    out_file = open(outfile, "w")   # opening output file
    for gr in r_table:              # iterating through genomic regions
        gr_snps = r_table[gr][0]            # getting genomic region's SNP table
        gr_conv = r_table[gr][1]            # getting genomic region's conversion table
        out_file.write("#" + gr + "\n")     # writing genomic region's name
        out_file.write(">SNP_TABLE\n")      # indicating SNP table listing in the following
        for row in gr_snps:                 # iterating through SNP table rows
            for entry in row[:-1]:              # iterating through row entries
                out_file.write(str(entry) + "\t")   # writing row entries (tab-separated)
            out_file.write(str(row[-1]) + "\n")
        out_file.write(">CONV_TABLE\n")     # indicating conversion table listing in the following
        for row in gr_conv:                 # iterating through conversion table rows
            for entry in row[:-1]:              # iterating through row entries
                out_file.write(str(entry) + "\t")   # writing row entries (tab-separated)
            out_file.write(str(row[-1]) + "\n")
        out_file.write("\n")

    out_file.close()  # closing output file
    return ()  # returning

def load_read_table(infile):
    """
    :param infile:  file storing read tables

    :return: read table; a hash assigning each genomic region's name a list containing a SNP and a conversion table
                        (implemented as nested lists, the first index corresponding to rows and the second to columns):
                            SNP table: row 1: genomic region's reference sequence positions (0-based), other rows:
                                reads, column 1: read names, other columns: reference sequence positions, entries: None
                                if the position is not covered by the read, 0 if read has no SNP at that position,
                                SNP nucleotide otherwise
                            conversion table: row 1: potential reference sequence conversion position, other rows:
                                reads, column 1: read names, other columns: potential reference sequence conversion
                                positions, entries: None if the position is not covered, masked by a deletion, a SNP or
                                SNP correction, 0 if read has no conversion at that position, 1 otherwise

    ....................................................................................................................

    This function loads a read table from a file. The file format should be as described in the following:

    #genomic_region_name\n
    >SNP_TABLE\n
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    >CONV_TABLE
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    """

    r_table = {}  # initializing read table (empty hash)
    in_file = open(infile, "r")  # opening input file

    gr = next(in_file).strip().split("#")[1]    # initialize current genomic region (first genomic region in file)
    r_table[gr] = []                            # initialize current genomic region's entry (empty list)
    table = []                                  # initialize current table (empty list)
    for line in in_file:            # iterating through file
        if line.startswith("#"):        # next genomic region encountered
            r_table[gr].append(table)           # storing current conversion table to current genomic region's entry
            table = []                          # re-setting table
            gr = line.strip().split("#")[1]     # setting next genomic region
            r_table[gr] = []                    # initializing next genomic region's entry
        elif line.startswith(">SNP"):   # SNP table encountered
            pass                                # nothing to do, everything was prepared in the step before
        elif line.startswith(">CONV"):  # conversion table encountered
            r_table[gr].append(table)           # storing current SNP table to current genomic region's entry
            table = []                          # re-setting table
        else:                           # table row
            entries = line.strip().split("\t")  # retrieving single row entries
            # re-transforming entries into integers, Nones and strings
            entries = [int(i) if str_is_int(i) else None if i == 'None' else i for i in entries]
            table.append(entries)               # storing row entries to current table
    r_table[gr].append(table)  # finally, adding last conversion table to last genomic region

    return(r_table)     # returning read table

def write_r_summary(r_summary, outfile):
    """
    :param r_summary:   read summary; a hash assigning each read's ID a list containing (in the following order):
                            reference sequence name mapped to
                            read orientation
                            read length
                            list of positions of T>C conversions observed (1-based)
                            list of positions of insertions (position after insertion event, 1-based)
                            list of positions of deletions (position of deletion event, 1-based)
                            hash assigning a nucleotide pairing its abundance within the read's mapping
    :param outfile:     output file to store read summary to

    :return:            void

    ....................................................................................................................

    This function writes the content of a read summary hash to a file. The file is formatted as described in the
    following:

    #read_id\n
    name_of_reference_sequence_read_is_mapped_to\n
    read_orientation\n
    read_length\n
    tab_delimited_positions_of_TC_conversions\n
    tab_delimited_positions_of_insertions\n
    tab_delimited_positions_of_deletions\n
    tab_delimited_key_value_pairs_for_nucleotide_pairing_abundances\n

    In case that there are no T>C conversions, insertions or deletions, respectively, the corresponding lines will be
    empty. The key-value pairs are written as key:value (nucleotide_pairing:abundance).
    """

    out_file = open(outfile, "w")  # opening output file
    for read_id in r_summary:  # iterating through read IDs in read summary hash
        read_record = r_summary[read_id]        # getting read ID's summary record
        out_file.write("#" + read_id + "\n")    # writing read ID
        out_file.write(read_record[0] + "\n")   # writing reference sequence mapped t
        out_file.write(str(read_record[1]) + "\n")  # writing read length
        out_file.write("\t".join([str(p) for p in read_record[2]]) + "\n")  # writing insertion positions
        out_file.write("\t".join([str(p) for p in read_record[3]]) + "\n")  # writing deletion positions
        out_file.write("\t".join([ntp + ":" + str(read_record[4][ntp]) for ntp in read_record[4]]) + "\n")  # writing
        #   nucleotide pairing abundances
    out_file.close()    # closing output file
    return ()           # returning

def load_r_summary(infile, gene_selection=None):
    """
    :param infile:              file storing read summary
    :param gene_selection:      file storing new-line separated list of genes to load only (default: None)

    :return:        read summary; a hash assigning each read's ID a list containing (in the following order):
                        reference sequence name mapped to
                        read orientation (if stored in the input file)
                        read length
                        list of positions of insertions (position after insertion event, 1-based)
                        list of positions of deletions (position of deletion event, 1-based)
                        hash assigning a nucleotide pairing its abundance within the read's mapping

    ....................................................................................................................

    This function loads a read summary from a file. The file format should be as described in the following:

    #read_id\n
    name_of_reference_sequence_read_is_mapped_to\n
    read_orientation\n                                  (this line is optional)
    read_length\n
    tab_delimited_positions_of_TC_conversions\n
    tab_delimited_positions_of_insertions\n
    tab_delimited_positions_of_deletions\n
    tab_delimited_key_value_pairs_for_nucleotide_pairing_abundances\n

    In case that there are no T>C conversions, insertions or deletions, respectively, the corresponding lines should be
    empty. The key-value pairs should be written as key:value (nucleotide_pairing:abundance).
    """

    in_file = open(infile)  # opening read summary file
    r_summary = {}          # initialize empty read summary hash
    genes = {}              # initialize empty hash with gene names to load
    if gene_selection:      # getting gene names to load
        in_genes = open(gene_selection)
        for line in in_genes:
            genes[line.strip()] = None

    for line in in_file:    # iterating though read summary file
        read_id = line.strip()[1:]                      # getting read ID
        ref_name = next(in_file).strip()                # getting reference sequence name mapped to
        if gene_selection and ref_name not in genes:    # skipping all genes not to be loaded
            next(in_file)
            next(in_file)
            next(in_file)
            next(in_file)
            continue
        read_len = int(next(in_file).strip())           # getting read length
        ins_pos = [int(p) for p in next(in_file).strip().split("\t") if p]  # getting list of insertion positions
        del_pos = [int(p) for p in next(in_file).strip().split("\t") if p]  # getting list of deletion positions
        ntp_abundance = {}                                      # initialize hash of nucleotide pair abundances
        for key_value in next(in_file).strip().split("\t"):     # iterate through key-value pairs
            key_value = key_value.split(":")                        # splitting key-value pair
            ntp_abundance[key_value[0]] = int(key_value[1])         # storing nucleotide pair abundance

        # storing read ID's record to read_summary hash
        r_summary[read_id] = [ref_name, read_len, ins_pos, del_pos, ntp_abundance]

    in_file.close()     # closing read summary file
    return(r_summary)   # returning read summary hash

def write_r_stats(rsum_stats, outfile):
    """
    :param rsum_stats:  read summary statistics; a list containing (in the following order):
                            total number of reads
                            hash: fragment length -> ratio of reads of that length
                            relative position of T>C conversions
                            relative position of insertions
                            relative position of deletions
                            a hash storing for each nucleotide pairing:
                                number of pairings -> % reads with that many pairings
                                % reference nucleotides paired with this read nucleotide (not considering indels)
                                % reads that have at least one such pairing and additionally at least one insertion
                                % reads that have at least one such pairing and additionally at least one deletion
                                relative position of additional insertions
                                relative position of additional deletions
                                hash: read length -> % reads of that length (of reads with at least one such pairing)
    :param outfile:     output file to write read summary statistics to
    :return:            void

    ....................................................................................................................

    This function writes read summary statistics to an output file. In particular, it writes the total number of reads
    recorded, nucleotide pairing statistics regarding the percentage of reads containing at least one such nucleotide
    pairing and the percentage of reference nucleotides being paired with the corresponding read nucleotide.
    """

    # defining nucleotide pairing keys for % reads and % reference nucleotides with one such nucleotide pairing
    read_percent_keys = ["AA", "TT", "CC", "GG", "AT", "AC", "AG", "AI", "AD", "TA", "TC", "TG", "TI", "TD",
                         "CA", "CT", "CG", "CI", "CD", "GA", "GT", "GC", "GI", "GD"]
    nt_percent_keys = ["AA", "TT", "CC", "GG", "AT", "AC", "AG", "TA", "TC", "TG", "CA", "CT", "CG", "GA", "GT", "GC"]

    # opening output file
    out_file = open(outfile, "w")

    # writing total number of reads
    out_file.write("total number of reads:\n" + str(rsum_stats[0]) + "\n")

    # writing % reads with one such nucleotide pairing
    out_file.write("% reads with at least one such nucleotide pairing\n")
    for k in read_percent_keys:
        out_file.write(k + "\t" + str((1 - rsum_stats[4][k][0][0]) * 100) + "\n")

    # writing % reference nucleotides paired with the corresponding read nucleotide
    out_file.write("% reference nucleotides paired with the corresponding read nucleotide\n")
    for k in nt_percent_keys:
        out_file.write(k + "\t" + str(rsum_stats[4][k][1] * 100) + "\n")

    # closing output file and returning
    out_file.close()
    return ()

def write_convcount_table(convcount_table, outfile):
    """
    :param convcount_table:     a hash assigning each genomic region's name another hash, assigning each read ID within
                                the genomic region a list containing the number of potential conversion positions and
                                the number of actual conversions
    :param outfile:             output file to store conversion counts table to

    :return: void

    ....................................................................................................................

    This function writes the content of a conversion counts table to a file. The file is formatted as described in the
    following:

    #genomic_region\n
    >read_id\n
    number_of_potential_conversion_positions\tnumber_of_conversions\n
    """

    out_file = open(outfile, "w")       # opening output file
    for gr in convcount_table:          # iterating through genomic regions in conversion counts table
        out_file.write("#" + gr + "\n")     # writing genomic region
        gr_entry = convcount_table[gr]      # getting genomic region's entry
        for read_id in gr_entry:            # iterating through genomic region's reads
            out_file.write(">" + read_id + "\n")    # writing read ID
            out_file.write(str(gr_entry[read_id][0]) + "\t" + str(gr_entry[read_id][1]) + "\n")     # writing counts

    out_file.close()    # closing output file
    return()            # returning

def load_convcount_table(infile):
    """
    :param infile:  file storing conversion counts table

    :return:    a hash assigning each genomic region's name another hash, assigning each read ID within the genomic
                region a list containing the number of potential conversion positions and the number of actual
                conversions

    ....................................................................................................................

    This function loads a conversion counts table from a file. The file format should be as described in the following:

    #genomic_region\n
    >read_id\n
    number_of_potential_conversion_positions\tnumber_of_conversions\n
    """

    in_file = open(infile)  # opening conversion counts table file
    convcounts_table = {}   # initialize empty conversion counts table hash

    current_gr = next(in_file).strip()[1:]      # storing current genomic region being processed
    gr_entry = {}                               # initialize current genomic region's entry
    current_read = next(in_file).strip()[1:]    # storing current read ID being processed
    for line in in_file:                        # iterating through conversion counts table file
        if line.startswith("#"):                    # new genomic region encountered
            convcounts_table[current_gr] = gr_entry     # storing former genomic region
            current_gr = line.strip()[1:]               # setting new genomic region
            gr_entry = {}                               # re-setting genomic region's entry
        elif line.startswith(">"):                  # new read ID encountered
            current_read = line.strip()[1:]             # setting new read ID
        else:                                       # counts encountered
            gr_entry[current_read] = [int(c) for c in line.strip().split("\t")]     # storing counts
    convcounts_table[current_gr] = gr_entry     # storing last genomic region, too

    in_file.close()             # closing input file
    return(convcounts_table)    # returning

def write_conveff_table(conveff_table, outfile):
    """
    :param conveff_table:   a hash assigning each genomic region's name a list containing two sublists:
                                the first sublist contains the sequence of estimated newly synthesized by total RNA
                                    ratios
                                the second sublist contains the sequence of estimated conversion efficiencies
    :param outfile:         output file to store conversion efficiency table to

    :return: void

    ....................................................................................................................

    This function writes the content of a conversion efficiencies table to a file. The file is formatted as described
    in the following:

    #genomic_region\n
    >new_by_total\n
    tab-separated listing of new-by-total ratio estimations\n
    >conv_eff\n
    tab-separated listing of conversion efficiency estimations\n
    """

    out_file = open(outfile, "w")   # opening output file
    for gr in conveff_table:        # iterating through genomic regions in conversion efficiency table

        out_file.write("#" + gr + "\n")         # writing genomic region
        out_file.write(">new_by_total\n")       # writing new-by-total ratios headline
        if conveff_table[gr][0] == None:
            out_file.write("None\n")
        else:
            for i in conveff_table[gr][0][:-1]:     # iterating through new-by-total entries
                out_file.write(str(i) + "\t")           # writing entry
            out_file.write(str(conveff_table[gr][0][-1]) + "\n")

        out_file.write(">conv_eff\n")           # writing conversion efficiencies headline
        if conveff_table[gr][0] == None:
            out_file.write("None\n")
        else:
            for i in conveff_table[gr][1][:-1]:     # iterating through conversion efficiency entries
                out_file.write(str(i) + "\t")           # writing entry
            out_file.write(str(conveff_table[gr][1][-1]) + "\n")

    out_file.close()    # closing output file
    return()            # returning

def load_conveff_table(infile):
    """
    :param infile:  file storing conversion efficiency table

    :return:    a hash assigning each genomic region's name a list containing two sublists:
                    the first sublist contains the sequence of estimated newly synthesized by total RNA ratios
                    the second sublist contains the sequence of estimated conversion efficiencies

    ....................................................................................................................

    This function loads a conversion efficiency table from a file. The file format should be as described in the
    following:

    #genomic_region\n
    >new_by_total\n
    tab-separated listing of new-by-total ratio estimations\n
    >conv_eff\n
    tab-separated listing of conversion efficiency estimations\n
    """

    in_file = open(infile)  # opening conversion efficiency table file
    conveff_table = {}      # initialize empty conversion efficiency table hash

    current_gr = next(in_file).strip()[1:]      # storing current genomic region being processed
    gr_entry = []                               # initialize current genomic region's entry
    for line in in_file:                        # iterating through conversion efficiency table file
        if line.startswith("#"):                    # new genomic region encountered
            conveff_table[current_gr] = gr_entry        # storing former genomic region
            current_gr = line.strip()[1:]               # setting new genomic region
            gr_entry = []                               # re-setting genomic region's entry
        elif line.startswith(">"):                  # headline for new-by-total ratios or conv. eff. encountered
            pass                                        # nothing to do, everything was prepared in the step before
        else:                                       # single entries encountered
            if line.startswith("None"):
                gr_entry.append(None)                                   # transforming entries to None values
            else:
                entries = [float(i) for i in line.strip().split("\t")]  # transforming entries to floats
                gr_entry.append(entries)
    conveff_table[current_gr] = gr_entry        # storing last genomic region, too

    in_file.close()         # closing input file
    return(conveff_table)   # returning

def write_summary_table(stable, outfile):
    """
    :param stable:      a nested hash, the outer hash assigning each genomic region's name an inner hash assigning the
                        description a list containing library size, total read counts, labeled read counts, conversion
                        efficiency estimation and newly synthesized transcripts ratio estimation:
                        region_name -> description -> [libsize, total, labeled, average potential conversion positions,
                                                       conv. efficiency, newly ratio]
    :param outfile:     output file to write summary matrix to

    :return: void

    This function writes a summary table to a file (according to the following format):
    region_name\tdescription\tlibrary_size\ttotal_reads\tlabeled_reads\taverage_potential_conversion_positions\t
    conv_efficiency\tnewly_synthesized_ratio\n
    """

    out_file = open(outfile, "w")   # opening output file

    for gr in stable:               # iterating through summary table hash, writing entries to output file
        for d in stable[gr]:
            entry = stable[gr][d]
            out_file.write(gr + "\t" + d + "\t" + "\t".join([str(e) for e in entry]) + "\n")

    out_file.close()        # closing output file
    return ()               # returning

def load_summary_table(infile):
    """
    :param infile:  file storing the summary table

    :return:    a nested hash, the outer hash assigning each genomic region's name an inner hash assigning the
                description a list containing library size, total read counts, labeled read counts, conversion
                efficiency estimation and newly synthesized transcripts ratio estimation:
                region_name -> description -> [libsize, total, labeled, average potential conversion positions,
                                               conv. efficiency, newly ratio]

    This function loads a summary table from a file. The file format should be as described in the following:
    region_name\tdescription\tlibrary_size\ttotal_reads\tlabeled_reads\tconv_efficiency\tnewly_synthesized_ratio\n
    """

    stable = {}                     # initializing summary table as empty hash
    in_file = open(infile, "r")     # opening summary table file

    for line in in_file:            # iterating through summary matrix file
        entry = line.strip().split("\t")    # getting a genomic region's entries
        gr_name = entry[0]                  # getting a genomic region's name

        if gr_name in stable:  # genomic region already stored:
            stable[gr_name][entry[1]] = [float(e) for e in entry[2:]]       # adding description's entry
        else:  # genomic region not yet stored:
            stable[entry[0]] = {entry[1]: [float(e) for e in entry[2:]]}    # setting up new genomic region's entry

    in_file.close()     # closing summary matrix file
    return(stable)      # returning summary matrix

# Data (Pre-) Processing ###############################################################################################

# --- file manipulation ---

def subset_fasta(fastafile, bedfile, outfile, refname_to_chr=r_to_c_hg19):
    """
    :param fastafile:   input fasta file from which a subset is to be built
    :param bedfile:     input BED file defining genomic regions to be extracted from the input fasta file
    :param outfile:     output file to store fasta subset to

    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number
                            (default: default hash)

    :return: void

    This function takes a subset from a given fasta file, extracting genomic regions as defined within a given BED file.
    """

    in_fasta = open(fastafile, "r")     # opening input fasta file
    in_bed = open(bedfile, "r")         # opening inpput BED file
    out_fasta = open(outfile, "w")      # opening output fasta file

    bed_regions = {}                    # hash storing BED regions as: chr_number -> {start_pos -> [end_pos, name]}

    # ..................................................................................................................
    # iterating through BED file, storing regions to hash
    # ..................................................................................................................

    for line in in_bed:                 # iterating though BED file entries

        fields = line.strip().split("\t")                   # getting BED file entry's single fields
        chr_number = fields[0][3:]                          # genomic region chromosome number
        startpos = int(fields[1])                           # genomic region starting position, 0-based, including
        endpos = int(fields[2])                             # genomic region ending position, 0-based, excluding
        name = fields[3]                                    # genomic region name

        if chr_number in bed_regions:                       # storing region to existing chromosome number
            bed_regions[chr_number][startpos] = [endpos, name]
        else:                                               # storing region with new chromosome number
            bed_regions[chr_number] = {startpos: [endpos, name]}

    # ..................................................................................................................
    # iterating through fasta file, extracting genomic regions:
    # 1) finding reference for which genomic regions are stored
    # 2) iterating through reference, processing genomic regions assigned to that reference
    # 3) as soon as all genomic regions of a reference were processed, skipping lines until new ref is encountered
    # ..................................................................................................................

    ref_name = ""                   # reference sequence name
    ref_chr = ""                    # reference sequence's assigned chromosome number
    ref_endpos = 0                  # reference sequence's current end position (0-based, excluding)
    ref_region_starts = []          # start positions of genomic regions assigned to current refseq
    skip = False                    # variable indicating whether the current reference sequence is to be skipped

    # iterating through fasta file

    for line in in_fasta:

        line = line.strip()     # removing newline character

        # 1) start of new reference sequence ---------------------------------------------------------------------------

        if line.startswith(">"):

            ref_name = line.split(" ")[0][1:]   # getting reference sequence name
            if ref_name not in refname_to_chr:  # (skipping refs not mappable to a chromosome number)
                skip = True
            else:                               # initializing refseq to process
                ref_chr = refname_to_chr[ref_name]      # setting refseq's chromosome number
                if ref_chr not in bed_regions:          # (skipping refs not in BED file)
                    skip = True
                else:
                    ref_endpos = 0                          # resetting refseq position counter
                    ref_region_starts = list(bed_regions[ref_chr].keys())  # getting genomic region start positions
                    ref_region_starts.sort()                # gr start positions in ascending order
                    skip = False                            # resetting the skip variable

        # 2) processing reference sequence's assigned genomic regions --------------------------------------------------

        elif not skip:

            for gr_start in ref_region_starts:

                while gr_start >= ref_endpos + len(line):       # going to start position
                    ref_endpos += len(line)                         # update reference ending position
                    line = next(in_fasta).strip()                   # get next line

                gr_start_idx = gr_start - ref_endpos            # genomic region in-line start index
                ref_endpos += len(line[:gr_start_idx])          # update reference end position
                line = line[gr_start_idx:]                      # update line to genomic region's start

                gr_name = bed_regions[ref_chr][gr_start][1]     # genomic region name
                gr_end = bed_regions[ref_chr][gr_start][0]      # genomic region end
                gr_seq = ""                                     # initialize genomic region sequence

                while gr_end > ref_endpos + len(line):          # going to end position
                    gr_seq += line                                  # update genomic region sequence
                    ref_endpos += len(line)                         # update reference end position
                    line = next(in_fasta).strip()                   # get next line

                gr_end_idx = gr_end - ref_endpos                # genomic region in-line end idx
                gr_seq += line[:gr_end_idx]                     # genomic region sequence

                seq_lines = math.ceil(len(gr_seq) / 80)         # gr sequence output lines (line length 80)
                out_fasta.write(">" + gr_name + "\n")           # writing genomic region identifier
                for l in range(0, seq_lines, 1):                # writing genomic region sequence
                    idx_start = l * 80                                  # gr seq start idx for output line
                    idx_end = (l + 1) * 80                              # gr seq end idx for output line
                    out_fasta.write(gr_seq[idx_start:idx_end] + "\n")   # writing output line

            skip = True     # as soon as all genomic regions are processed, mark remaining refseq to be skipped

    # ..................................................................................................................
    # closing files and returning
    # ..................................................................................................................

    in_bed.close()
    in_fasta.close()
    out_fasta.close()
    return()

def degen_fasta(fastafile, bedfile, outfile):
    """
    :param fastafile:   input fasta file to be degenerated
    :param bedfile:     input BED file defining genomic regions and their strand orientation
    :param outfile:     output file to store degenerated fasta to

    :return: void

    This function degenerates the sequences from a given fasta file, based on the sequences' strand orientations as
    defined within a given BED file. Nucleotides are degenerated based on the metabolic labeling conversion expected.
    """

    in_fasta = open(fastafile, "r")     # opening input fasta file
    in_bed = open(bedfile, "r")         # opening input BED file
    out_fasta = open(outfile, "w")      # opening output fasta file

    # iterating through BED file, storing genomic regions to hash ------------------------------------------------------

    bed_regions = {}        # hash storing BED regions as: name -> strand
    for line in in_bed:     # iterating though BED file entries

        fields = line.strip().split("\t")   # getting BED file entry's single fields
        name = fields[3]                    # genomic region name
        strand = fields[5]                  # genomic region strand
        bed_regions[name] = strand          # storing genomic region

    # iterating through input fasta file, degenerating sequences -------------------------------------------------------

    current_name = next(in_fasta).strip().split(" ")[0][1:]     # getting first sequence's name
    current_strand = bed_regions[current_name]                  # getting first sequence's strand orientation
    current_seq = ""                                            # initializing currently processed nucleotide sequence

    for line in in_fasta:       # iterating through input fasta file

        if not line.startswith(">"):        # still processing current sequence
            current_seq += line.strip()         # updating current sequence's nucleotide sequence
        else:                               # new sequence encountered
            current_seq = current_seq.upper()   # transforming soft-masked lower-case

            if current_strand == "-" or current_strand == ".":          # --- CASE: --- reference forward or both

                seq1 = ["A" if i == "G" else i for i in current_seq]    # transform Gs to As
                seq1 = "".join(seq1)

                seq_lines = math.ceil(len(current_seq) / 80)        # nucleotide sequence output lines (line length 80)
                out_fasta.write(">" + current_name + "_1degen\n")   # writing sequence identifier (G>A)
                for l in range(0, seq_lines, 1):                    # writing degenerated sequence (G>A)
                    idx_start = l * 80                                  # ntseq start idx for output line
                    idx_end = (l + 1) * 80                              # ntseq end idx for output line
                    out_fasta.write(seq1[idx_start:idx_end] + "\n")     # writing output line

            if current_strand == "+" or current_strand == ".":          # --- CASE: --- reference reverse or both

                seq2 = ["T" if i == "C" else i for i in current_seq]    # transform Cs to Ts
                seq2 = "".join(seq2)

                seq_lines = math.ceil(len(current_seq) / 80)        # nucleotide sequence output lines (line length 80)
                out_fasta.write(">" + current_name + "_2degen\n")   # writing sequence identifier (C>T)
                for l in range(0, seq_lines, 1):                    # writing degenerated sequence (C>T)
                    idx_start = l * 80                                  # ntseq start idx for output line
                    idx_end = (l + 1) * 80                              # ntseq end idx for output line
                    out_fasta.write(seq2[idx_start:idx_end] + "\n")     # writing output line

            current_name = line.strip().split(" ")[0][1:]       # getting next sequence's name
            current_strand = bed_regions[current_name]          # getting strand orientation
            current_seq = ""                                    # resetting nucleotide sequence

    # processing last input fasta sequence aswell ----------------------------------------------------------------------

    current_seq = current_seq.upper()       # transforming soft-masked lower-case

    if current_strand == "-" or current_strand == ".":      # --- CASE: --- reference forward or both

        seq1 = ["A" if i == "G" else i for i in current_seq]    # transform Gs to As
        seq1 = "".join(seq1)

        seq_lines = math.ceil(len(current_seq) / 80)            # nucleotide sequence output lines (line length 80)
        out_fasta.write(">" + current_name + "_1degen\n")       # writing sequence identifier (G>A)
        for l in range(0, seq_lines, 1):                        # writing degenerated sequence (G>A)
            idx_start = l * 80                                      # ntseq start idx for output line
            idx_end = (l + 1) * 80                                  # ntseq end idx for output line
            out_fasta.write(seq1[idx_start:idx_end] + "\n")         # writing output line

    if current_strand == "+" or current_strand == ".":      # --- CASE: --- reference reverse or both

        seq2 = ["T" if i == "C" else i for i in current_seq]    # transform Cs to Ts
        seq2 = "".join(seq2)

        seq_lines = math.ceil(len(current_seq) / 80)            # nucleotide sequence output lines (line length 80)
        out_fasta.write(">" + current_name + "_2degen\n")       # writing sequence identifier (C>T)
        for l in range(0, seq_lines, 1):                        # writing degenerated sequence (C>T)
            idx_start = l * 80                                      # ntseq start idx for output line
            idx_end = (l + 1) * 80                                  # ntseq end idx for output line
            out_fasta.write(seq2[idx_start:idx_end] + "\n")         # writing output line

    # closing files and returning --------------------------------------------------------------------------------------

    in_fasta.close()
    in_bed.close()
    out_fasta.close()
    return()

def degen_fastq(fastqfile, outfile):
    """
    :param fastqfile:   input fastq.gz file to be degenerated
    :param outfile:     output file to store degenerated fastq.gz to

    :return: void

    This function degenerates the reads from a fastq.gz file. Two sets of reads are created, one being C>T and the other
    one G>A degenerated. Read names are adapted, appending "_degen1" for C>T and "_degen2" for G>A reads.
    """

    in_file = gzip.open(fastqfile, "rt")        # opening input fastq.gz file
    out_file = gzip.open(outfile, "wb")         # opening output fastq.gz file

    for line in in_file:    # iterating through fastq file reads

        if line.startswith("@"):                # beginning of a read's entry
            read_description = line.split(" ")          # getting read description
            read_id = read_description[0]               # getting read ID
            read_seq = next(in_file).strip().upper()    # getting read sequence (soft-masked lower-case patched)
            spacer = next(in_file)                      # getting this weird spacer noone knows what its purpose is
            read_qual = next(in_file)                   # getting read quality string

            read_id_1 = read_id + "_degen1\n"           # read ID for degenerated read (1)
            read_id_2 = read_id + "_degen2\n"           # read ID for degenerated read (2)
            read_seq_1 = "".join(["T" if i == "C" else i for i in read_seq]) + "\n"     # degenerating read (1)
            read_seq_2 = "".join(["A" if i == "G" else i for i in read_seq]) + "\n"     # degenerating read (2)

            out_file.write(read_id_1.encode() + read_seq_1.encode() + spacer.encode() +
                           read_qual.encode())          # writing degenerated read (1) to output fastq
            out_file.write(read_id_2.encode() + read_seq_2.encode() + spacer.encode() +
                           read_qual.encode())          # writing degenerated read (2) to output fastq

    in_file.close()     # closing input file
    out_file.close()    # closing output file
    return()            # returning

def bed_overlap_regions(bed_infile, bed_outfile, merge=False):
    """
    :param bed_infile:      input BED file to be processed
    :param bed_outfile:     output BED file
    :param merge:           set to TRUE to also merge BED regions with different strand orientation (default: FALSE)
    :return: void

    This function merges all overlapping regions defined in a BED file.
    """

    bed_in = open(bed_infile, "r")      # opening input BED file
    bed_out = open(bed_outfile, "w")    # opening output BED file

    # ..................................................................................................................
    # SUBROUTINE: keep strand orientation separate
    # ..................................................................................................................

    if not merge:

        fields = next(bed_in).strip().split("\t")   # BED first line
        bed_hash = {fields[0]: {fields[5]: {int(fields[1]): {int(fields[2]): [fields[3], fields[4]]}
                                            }}}     # hash storing chr -> strand -> start -> end -> [name, score]

        # loading in BED file to hash ----------------------------------------------------------------------------------

        for line in bed_in:                 # iterating through input BED file

            fields = line.strip().split("\t")   # getting genomic region's single fields
            chrom = fields[0]                   # chromosome
            start = int(fields[1])              # start position
            end = int(fields[2])                # end position
            name = fields[3]                    # name
            score = fields[4]                   # score
            strand = fields[5]                  # strand

            if chrom in bed_hash:               # chromosome already stored to BAM hash
                if strand in bed_hash[chrom]:           # duplicate strand
                    if start in bed_hash[chrom][strand]:    # duplicate strand and start
                        if end in bed_hash[chrom][strand][start]:   # duplicate strand, start, end (identical regions)
                            pass
                        else:
                            bed_hash[chrom][strand][start][end] = [name, score]
                    else:
                        bed_hash[chrom][strand][start] = {end: [name, score]}
                else:
                    bed_hash[chrom][strand] = {start: {end: [name, score]}}
            else:
                bed_hash[chrom] = {strand: {start: {end: [name, score]}}}

        # iterating through genomic regions, merging overlapping ones --------------------------------------------------

        # iterating through chromosomes

        for chrom in bed_hash:

            chrom_hash = bed_hash[chrom]                                        # getting chromosome's subhash
            chrom_hash_minus = chrom_hash["-"] if "-" in chrom_hash else {}     # getting minus strand subhash
            chrom_hash_plus = chrom_hash["+"] if "+" in chrom_hash else {}      # getting plus strand subhash

            startpos_minus = list(chrom_hash_minus.keys())      # getting chromosome's minus strand start positions
            startpos_minus.sort()                               # sorting minus strand start positions
            startpos_plus = list(chrom_hash_plus.keys())        # getting chromosome's plus strand start positions
            startpos_plus.sort()                                # sorting plus strand start positions

            # processing strand orientations separately

            for strand, o, startpos in [["-", "r", startpos_minus], ["+", "f", startpos_plus]]:

                if not startpos: continue                                   # strand information empty: continue

                current_start = startpos.pop(0)                             # initialize current start position
                endpos = list(chrom_hash[strand][current_start].keys())     # getting assigned end positions
                current_end = max(endpos)                                   # initialize current end position

                if len(endpos) > 1:                                 # multiple assigned end positions (overlap):
                    overlap = True                                      # initialize 'overlap' with 'True'
                else:                                               # else:
                    overlap = False                                     # initialize 'overlap' with 'False'

                # iterating through genomic region's start positions

                for p in startpos:

                    endpos = list(chrom_hash[strand][p].keys())     # getting assigned end positions

                    # --- CASE: ---
                    # current and next genomic region do not overlap

                    if p >= current_end:

                        # CASE: current region is not built from multiple, overlapping regions
                        # ACTION: store current region
                        if not overlap:

                            name = chrom_hash[strand][current_start][current_end][0]    # current region name
                            score = chrom_hash[strand][current_start][current_end][1]   # current region score
                            bed_out.write("\t".join([chrom, str(current_start), str(current_end),
                                                     name, score, strand]) + "\n")      # writing current region

                        # CASE: current region is built from multiple, overlapping regions
                        # ACTION: re-specify and store current region
                        else:

                            name = "_".join([chrom, str(current_start), str(current_end), o])   # define name
                            score = "0"                                                         # set score to 0
                            bed_out.write("\t".join([chrom, str(current_start), str(current_end),
                                                     name, score, strand]) + "\n")      # writing current region

                        # updating current genomic region

                        current_start = p               # update current start position
                        current_end = max(endpos)       # update current end position

                        if len(endpos) > 1:             # multiple assigned end positions (overlap):
                            overlap = True                  # update 'overlap' to 'True'
                        else:                           # else:
                            overlap = False                 # update 'overlap' to 'False'

                    # --- CASE: ---
                    # current and next genomic region do overlap
                    # ACTION: update current region

                    else:

                        current_end = max(endpos + [current_end])   # update current end position
                        overlap = True                              # update 'overlap'

                # writing chromosome's last remaining genomic region ---------------------------------------------------

                if not overlap:     # no multiple, overlapped regions
                    name = chrom_hash[strand][current_start][current_end][0]
                    score = chrom_hash[strand][current_start][current_end][1]
                    bed_out.write("\t".join([chrom, str(current_start), str(current_end), name, score, strand]) + "\n")
                else:               # multiple, overlapped regions
                    name = "_".join([chrom, str(current_start), str(current_end), o])
                    score = "0"
                    bed_out.write("\t".join([chrom, str(current_start), str(current_end), name, score, strand]) + "\n")

    # ..................................................................................................................
    # SUBROUTINE: merge strand orientation
    # ..................................................................................................................

    else:

        fields = next(bed_in).strip().split("\t")       # BED first line
        bed_hash = {fields[0]: {int(fields[1]): {int(fields[2]): {fields[5]: [fields[3], fields[4]]}
                                                 }}}    # hash storing chr -> start -> end -> strand -> [name, score]

        # loading in BED file to hash ----------------------------------------------------------------------------------

        for line in bed_in:                 # iterating through input BED file

            fields = line.strip().split("\t")   # getting genomic region's single fields
            chrom = fields[0]                   # chromosome
            start = int(fields[1])              # start position
            end = int(fields[2])                # end position
            name = fields[3]                    # name
            score = fields[4]                   # score
            strand = fields[5]                  # strand

            if chrom in bed_hash:               # chromosome already stored to BAM hash
                if start in bed_hash[chrom]:            # duplicate start
                    if end in bed_hash[chrom][start]:       # duplicate start and end
                        if strand in bed_hash[chrom][start][end]:   # duplicate start, end, strand (identical regions)
                            pass
                        else:
                            bed_hash[chrom][start][end][strand] = [name, score]
                    else:
                        bed_hash[chrom][start][end] = {strand: [name, score]}
                else:
                    bed_hash[chrom][start] = {end: {strand: [name, score]}}
            else:
                bed_hash[chrom] = {start: {end: {strand: [name, score]}}}

        # iterating through genomic regions, merging overlapping ones --------------------------------------------------

        for chrom in bed_hash:              # iterating through chromosomes
    
            chrom_hash = bed_hash[chrom]            # getting chromosome's subhash
            startpos = list(chrom_hash.keys())      # getting chromosome's start positions
            startpos.sort()                         # sorting start positions
    
            current_start = startpos.pop(0)                     # initialize current start position
            endpos = list(chrom_hash[current_start].keys())     # getting assigned end positions
            current_end = max(endpos)                           # initialize current end position
            current_strands = [k for e in endpos
                               for k in chrom_hash[current_start][e]]   # initialize assigned strands

            if len(endpos) > 1:                                 # multiple assigned end positions (overlap):
                overlap = True                                      # initialize 'overlap' with 'True'
            else:                                               # else:
                overlap = False                                     # initialize 'overlap' with 'False'
    
            for p in startpos:              # iterating through genomic region's start positions
    
                endpos = list(chrom_hash[p].keys())                         # getting assigned end positions
                strands = [k for e in endpos for k in chrom_hash[p][e]]     # getting assigned strands
    
                # --- CASE: ---
                # current and next genomic region do not overlap
    
                if p >= current_end:
    
                    # CASE: current region is not built from multiple, overlapping regions
                    # ACTION: store current region
                    if not overlap:
    
                        name = chrom_hash[current_start][current_end][current_strands[0]][0]    # current region name
                        score = chrom_hash[current_start][current_end][current_strands[0]][1]   # current region score
                        bed_out.write("\t".join([chrom, str(current_start), str(current_end),
                                                 name, score, current_strands[0]]) + "\n")      # writing current region
    
                    # CASE: current region is built from multiple, overlapping regions
                    # ACTION: re-specify and store current region
                    else:
    
                        name = chrom + "_" + str(current_start) + "_" + str(current_end)    # define name
                        score = "0"                                                         # set score to 0
                        strand = current_strands[0]                                         # set strand
                        if "+" in current_strands:                                  # if multiple strand specifications
                            if "-" in current_strands or "." in current_strands:    # are present, set 'strand' to '.'
                                strand = "."
                        elif "-" in current_strands and "." in current_strands:
                            strand = "."
    
                        bed_out.write("\t".join([chrom, str(current_start), str(current_end),
                                                 name, score, strand]) + "\n")              # writing current region
    
                    # updating current genomic region
    
                    current_start = p               # update current start position
                    current_end = max(endpos)       # update current end position
                    current_strands = strands       # update current strands
                    if len(endpos) > 1:             # multiple assigned end positions (overlap):
                        overlap = True                  # update 'overlap' to 'True'
                    else:                           # else:
                        overlap = False                 # update 'overlap' to 'False'
    
                # --- CASE: ---
                # current and next genomic region do overlap
                # ACTION: update current region
    
                else:
    
                    current_end = max(endpos + [current_end])   # update current end position
                    current_strands += strands                  # update current strands
                    overlap = True                              # update 'overlap'
    
            # writing chromosome's last remaining genomic region -------------------------------------------------------
    
            if not overlap:     # no multiple, overlapped regions
                name = chrom_hash[current_start][current_end][current_strands[0]][0]
                score = chrom_hash[current_start][current_end][current_strands[0]][1]
                bed_out.write("\t".join([chrom, str(current_start), str(current_end),
                                         name, score, current_strands[0]]) + "\n")
            else:               # multiple, overlapped regions
                name = chrom + "_" + str(current_start) + "_" + str(current_end)
                score = "0"
                strand = current_strands[0]
                if "+" in current_strands:
                    if "-" in current_strands or "." in current_strands: strand = "."
                elif "-" in current_strands and "." in current_strands: strand = "."
                bed_out.write("\t".join([chrom, str(current_start), str(current_end), name, score, strand]) + "\n")

    # ..................................................................................................................
    # closing files and returning
    # ..................................................................................................................

    bed_in.close()
    bed_out.close()
    return()

def bam_region_to_chromosome(bamfile_in, bamfile_template, bedfile, bamfile_out, suffix=None,
                             chr_to_refname=c_to_r_hg19):
    """
    :param bamfile_in:          input BAM file for which reference names should be changed from regions to chromosomes
    :param bamfile_template:    BAM file containing header template for output BAM file
    :param bedfile:             BED file defining genomic regions
    :param bamfile_out:         output BAM file

    :param suffix:              length of a string suffix appended to the references' names (as originally defined
                                within the BED file) in the input BAM file (default: None)
    :param chr_to_refname:      hash storing which chromosome number refers to which reference sequence name
                                (default: default hash)

    :return: void

    This function modifies the reference names of a given input BAM file, changing genomic region names to chromosome
    names and updating mapping positions correspondingly.
    """

    # initializing -----------------------------------------------------------------------------------------------------

    sam_file_in = bamfile_in + "_tmp.sam"           # input BAM file temporary SAM version (need SAM for read/write)
    sam_file_out = bamfile_out + "_tmp.sam"         # output BAM file temporary SAM version (need SAM for read/write)
    bam_unsorted = bamfile_out + "_unsorted.bam"    # unsoted output BAM file (temporary)

    bam_to_sam = subprocess.Popen("samtools view -h -o " + sam_file_in + " " + bamfile_in,
                                  shell=True)       # converting input BAM file to input SAM file
    bam_to_sam.wait()                               # waiting for conversion to finish

    write_header = subprocess.Popen("samtools view -H -o " + sam_file_out + " " + bamfile_template,
                                    shell=True)     # writing header to output SAM file
    write_header.wait()                             # waiting for header writing to finish

    # loading genomic regions into hash
    # region_name -> [chr_number, starting_pos] ------------------------------------------------------------------------

    bed_in = open(bedfile, "r")     # opening BED file with genomic regions
    genomic_regions = {}            # initialize hash storing genomic region information
    for line in bed_in:             # iterating through BED file (genomic region entries)

        fields = line.strip().split("\t")                   # getting entrie's single fields
        chr_number = fields[0][3:]                          # genomic region chromosome number
        startpos = int(fields[1])                           # genomic region starting position, 0-based, including
        name = fields[3]                                    # genomic region name
        genomic_regions[name] = [chr_number, startpos]      # storing region

    # iterating through input SAM file,
    # writing modified reads to output SAM file ------------------------------------------------------------------------

    if suffix: suffix *= -1                 # initializing index of reference name suffix
    sam_in = open(sam_file_in, "r")         # opening input SAM file
    sam_out = open(sam_file_out, "a")       # opening output SAM file

    for line in sam_in:                     # iterating through input SAM file
        if line.startswith("@"): continue       # skipping header lines

        read_info = line.strip().split("\t")    # getting read information fields
        read_refname = read_info[2][:suffix]    # getting read's reference name (genomic region name)
        read_startpos = int(read_info[3])       # getting read's mapping start position (1-based, including)

        ref_region = genomic_regions[read_refname]  # getting reference genomic region's information
        ref_chr = ref_region[0]                     # genomic region's chromosome number
        ref_start = ref_region[1]                   # genomic region starting position (0-based, including)

        read_refname_updated = chr_to_refname[ref_chr]      # updating reference name (chromosome name)
        read_startpos_updated = ref_start + read_startpos   # updating mapping start position (1-based, including)
        read_info[2] = read_refname_updated                 # modify read information (reference name)
        read_info[3] = str(read_startpos_updated)           # modify read information (mapping position)

        read_modified = "\t".join(read_info) + "\n"     # generating modified read
        sam_out.write(read_modified)                    # writing modified read

    # converting output SAM to output BAM, creating BAM index,
    # deleting temporary files, returning ------------------------------------------------------------------------------

    sam_in.close()              # closing input SAM file
    sam_out.close()             # closing output SAm file

    sam_to_bam = subprocess.Popen("samtools view -h -b -o " + bam_unsorted + " " + sam_file_out,
                                  shell=True)   # converting temporary output SAM to output BAM
    sam_to_bam.wait()                           # waiting for conversion to be finished
    bam_sort = subprocess.Popen("samtools sort -o " + bamfile_out + " " + bam_unsorted,
                                shell=True)     # sorting BAM file (needed for indexing)
    bam_sort.wait()                             # waiting for sorting to be finished
    bam_index = subprocess.Popen("samtools index " + bamfile_out,
                                 shell=True)    # indexing output BAM file
    bam_index.wait()                            # waiting for indexing to be finished

    os.remove(sam_file_in)      # removing temporary input SAM file
    os.remove(sam_file_out)     # removing temporary output SAM file
    os.remove(bam_unsorted)     # removing temporary unsorted BAM file

    return()    # returning

def bam_recover_sequence(bamfile_in, fastq_in, bamfile_out):
    """
    :param bamfile_in:      input BAM file to be restored
    :param fastq_in:        input fastq.gz file containing original reads (not-degenerated reads)
    :param bamfile_out:     output BAM file to store recovered reads to

    :return: void

    This function recovers the read sequence within a BAM file for degenerated reads' mapping. It also checks which
    reads have been mapped twice (both degenerated versions have been mapped) and filters out these reads completely.
    """

    # initializing -----------------------------------------------------------------------------------------------------

    sam_file_in = bamfile_in + "_tmp.sam"       # input BAM file temporary SAM version (need SAM for read/write)
    sam_file_out = bamfile_out + "_tmp.sam"     # output BAM file temporary SAM version (need SAM for read/write)

    bam_to_sam = subprocess.Popen("samtools view -h -o " + sam_file_in + " " + bamfile_in,
                                  shell=True)   # converting input BAM file to input SAM file
    bam_to_sam.wait()                           # waiting for conversion to finish

    write_header = subprocess.Popen("samtools view -H -o " + sam_file_out + " " + bamfile_in,
                                    shell=True)     # writing header to output SAM file
    write_header.wait()                             # waiting for header writing to finish

    # building a hash of all read IDs that were mapped uniquely --------------------------------------------------------
    # (uniquely: only one degenerated version of the read has been mapped, and only once)

    in_bam = pysam.AlignmentFile(bamfile_in, "rb")      # opening input BAM file
    unique_ids = {}                                     # initialize hash storing uniquely mapped read IDs

    for read in in_bam.fetch():         # iterating through reads
        read_id_raw = read.query_name       # getting read ID
        read_id_pure = read_id_raw[:-7]     # getting pure read ID (ID without suffix marking the degenerated version)
        if read_id_pure in unique_ids:      # read ID occurred before: increment count
            unique_ids[read_id_pure] += 1
        else:                               # read ID is new: initialize with 0
            unique_ids[read_id_pure] = 0

    for read_id in list(unique_ids.keys()):     # iterating through read IDs' counts
        if unique_ids[read_id] > 0:                 # any read ID that occurred multiple times is deleted from hash
            del unique_ids[read_id]

    in_bam.close()                      # closing input BAM file again

    # getting original read sequences for uniquely mapped read IDs -----------------------------------------------------

    in_fastq = gzip.open(fastq_in, "rt")    # opening input fastq.gz file
    for line in in_fastq:
        read_id = line.strip().split(" ")[0][1:]        # getting read's ID
        read_sequence = next(in_fastq).strip()[12:]     # getting read's sequence (removing 12, hard-clipped nts)
        if read_id in unique_ids:                       # read is mapped uniquely
            unique_ids[read_id] = read_sequence             # store the read's original sequence
        next(in_fastq)                                  # skipping next line (read entry spacer)
        next(in_fastq)                                  # skipping next line (read entry quality string)

    # recovering read sequences in BAM file ----------------------------------------------------------------------------

    sam_in = open(sam_file_in, "r")         # opening input SAM file
    sam_out = open(sam_file_out, "a")       # opening output SAM file

    for line in sam_in:                     # iterating through input SAM file
        if line.startswith("@"): continue       # skipping header lines

        read_info = line.strip().split("\t")            # getting read information fields
        read_id = read_info[0][:-7]                     # getting read ID (without suffix of the degenerated version)
        if read_id not in unique_ids: continue          # read ID is not mapped uniquely: continue with next read

        read_seq = unique_ids[read_id]                          # initialize read sequence with original sequence
        if len(read_seq) != len(read_info[9]):                  # length of original and mapped sequence disagree:
            continue                                                # discard that read (dunno what happened here)

        read_flag = int(read_info[1])                           # read flag
        binary_flag = "{0:b}".format(read_flag)                 # read flag (binary mode)
        if len(binary_flag) > 4 and binary_flag[-5] == "1":     # flag for read being reverse complemented
            read_seq = read_seq[::-1]                               # reverting read sequence
            read_seq = "".join(["A" if i == "T" else "T" if i == "A"
                                else "G" if i == "C" else "C" if i == "G"
                                else i for i in read_seq])          # complementing read sequence

        read_info[0] = read_id                          # updating read ID (without degenerated version's suffix)
        read_info[9] = read_seq                         # updating read information (read sequence)
        read_modified = "\t".join(read_info) + "\n"     # generating modified read
        sam_out.write(read_modified)                    # writing modified read

    # converting output SAM to output BAM, creating BAM index,
    # deleting temporary files, returning ------------------------------------------------------------------------------

    sam_in.close()      # closing input SAM file
    sam_out.close()     # closing output SAm file

    sam_to_bam = subprocess.Popen("samtools view -h -b -o " + bamfile_out + " " + sam_file_out,
                                  shell=True)   # converting temporary output SAM to output BAM
    sam_to_bam.wait()                           # waiting for conversion to be finished
    bam_index = subprocess.Popen("samtools index " + bamfile_out,
                                 shell=True)    # indexing output BAM file
    bam_index.wait()                            # waiting for indexing to be finished

    os.remove(sam_file_in)      # removing temporary input SAM file
    os.remove(sam_file_out)     # removing temporary output SAm file

    return ()   # returning

def bam_merge(bam_one_in, bam_two_in, bam_out):
    """
    :param bam_one_in:          input BAM file 1 to be merged
    :param bam_two_in:          input BAM file 2 to be merged
    :param bam_out:             output BAM file

    :return: void

    This function merges two BAM files, automatically sorting and creating an index for the output file.
    """

    bam_one_sorted = bam_one_in + "_tmp_sorted.bam"     # sorted input BAM file 1 (temporary)
    bam_one_sort = subprocess.Popen("samtools sort -o " + bam_one_sorted + " " + bam_one_in,
                                    shell=True)         # sorting input BAM file 1
    bam_one_sort.wait()                                 # waiting for sorting to be finished

    bam_two_sorted = bam_two_in + "_tmp_sorted.bam"     # sorted input BAM file 2 (temporary)
    bam_two_sort = subprocess.Popen("samtools sort -o " + bam_two_sorted + " " + bam_one_in,
                                    shell=True)         # sorting input BAM file 2
    bam_two_sort.wait()                                 # waiting for sorting to be finished

    bam_merge = subprocess.Popen("samtools merge " + bam_out + " " + bam_one_sorted + " " + bam_two_sorted,
                                 shell=True)            # merging input files
    bam_merge.wait()                                    # waiting for merging to be finished

    bam_index = subprocess.Popen("samtools index " + bam_out, shell=True)   # indexing output BAM file
    bam_index.wait()                                                        # waiting for indexing to be finished

    os.remove(bam_one_sorted)   # removing temporary sorted BAM file 1
    os.remove(bam_two_sorted)   # removing temporary sorted BAM file 2
    return()                    # returning

def covered_positions_merge(cp_files, output_file):
    """
    :param cp_files:        list of files storing covered positions that are to be merged
    :param output_file:     output file to store merged covered positions to

    :return: void

    This function merges covered positions as stored in the input files given. In detail, the coverages are summed up
    over all input files for each position that is recorded in any of the input files.
    """

    # initialize merged covered positions' hashes with first covered positions hashes given
    [cp_merged_forward, cp_merged_reverse] = load_covered_positions(cp_files.pop(0))

    # iterating through covered positions' files, merging to merged covered positions' hash
    for file in cp_files:

        [cp_hash_forward, cp_hash_reverse] = load_covered_positions(file)   # loading covered positions

        # processing forwards mapped reads' covered positions

        for ref_name in cp_hash_forward:            # iterating through reference sequences
            ref_hash = cp_hash_forward[ref_name]    # reference sequences' covered positions hash

            if ref_name not in cp_merged_forward:           # CASE: reference not in merged cov. pos. hash
                cp_merged_forward[ref_name] = ref_hash          # simply copy reference hash to merged hash

            else:                                           # CASE: reference in merged cov. pos. hash
                ref_merged = cp_merged_forward[ref_name]        # reference sequence's merged hash
                for cp in ref_hash:                             # iterating through covered positions, ...
                    if cp in ref_merged:                        # ... CASE: covered position in merged hash
                        ref_merged[cp] += ref_hash[cp]              # adding coverage to merged hash entry
                    else:                                       # ... CASE: covered position not in merged hash
                        ref_merged[cp] = ref_hash[cp]               # copy coverage to merged hash entry

        # processing reversely mapped reads' covered positions

        for ref_name in cp_hash_reverse:            # iterating through reference sequences
            ref_hash = cp_hash_reverse[ref_name]    # reference sequences' covered positions hash

            if ref_name not in cp_merged_reverse:           # CASE: reference not in merged cov. pos. hash
                cp_merged_reverse[ref_name] = ref_hash          # simply copy reference hash to merged hash

            else:                                           # CASE: reference in merged cov. pos. hash
                ref_merged = cp_merged_reverse[ref_name]        # getting reference sequence's merged hash
                for cp in ref_hash:                             # iterating through covered positions, ...
                    if cp in ref_merged:                        # ... CASE: covered position in merged hash
                        ref_merged[cp] += ref_hash[cp]              # adding coverage to merged hash entry
                    else:                                       # ... CASE: covered position not in merged hash
                        ref_merged[cp] = ref_hash[cp]               # copy coverage to merged hash entry

    # writing merged covered positions to output file and returning
    write_covered_positions([cp_merged_forward, cp_merged_reverse], output_file)
    return()

def cluster_tables(cluster_results, regions_out, clusters_out, time_series, data_description=""):
    """
    :param cluster_results:     file containing the clustering results
    :param regions_out:         output file to write genomic regions' ratios table to
    :param clusters_out:        output file to write clusters' center table to
    :param time_series:         list containing time points at which measurement were taken
    :param data_description:    optional; data description to add to column names in output tables (default: "")

    :return:    void

    ....................................................................................................................

    This function retrieves data tables from a clustering results input file. In detail, it extracts the genomic
    regions' ratios table underlying the clustering process, and the clusters' center table being the result of the
    clustering process. It writes the tables into two distinct output files:

    ....................................................................................................................
    Genomic regions' ratio values
    ....................................................................................................................

    a tab-separated output file in the following format:

    Columns: time points' ratios
        first column:   genomic region names
        second column:  genomic region's cluster
        first half:     modified/total ratios, for nucleus and cytosol, ordered by time series
        second half:    total/library_size ratios, for nucleus and cytosol, ordered by time series
    Rows:
        first row:      genomic region name, cluster and time points itself (four times listing of the time series,
                        for mod/total nucleus, mod/total cytosol, total/lib nucleus, total/lib cytosol)
        next rows:      genomic region's names, clusters and values

    Output file example:

    columns (in the following order):
    region names, region clusters, mod/total nucleus 0min, mod/total nucleus 5min, mod/total cytosol 0min,
    mod/total cytosol 5min, total/lib nucleus 0min, total/lib nucleus 5min, total/lib cytosol 0min,
    total/lib cytosol 5min

    name    cluster     0   5   0   5   0   5   0   5       (time points)
    gr1     6           1   2   1   3   1   1   1   1       (genomic region 1 values)
    gr2     9           3   15  2   3   1   2   1   1       (genomic region 2 values)
    gr3     2           4   4   1   6   2   4   3   6       (genomic region 3 values)

    ....................................................................................................................
    Clusters' center values
    ....................................................................................................................

    a tab-separated output file in the following format:

    Columns: time points' ratios
        first column:   cluster name
        second column:  cluster size (number of genomic regions within that cluster)
        first half:     modified/total ratios, for nucleus and cytosol, ordered by time series
        second half:    total/library_size ratios, for nucleus and cytosol, ordered by time series
    Rows:
        first row:      cluster size and time points itself (four times listing of the time series,
                        for mod/total nucleus, mod/total cytosol, total/lib nucleus, total/lib cytosol)
        next rows:      clusters' size and center values

    Output file example:

    columns (in the following order):
    cluster name, cluster size, mod/total nucleus 0min, mod/total nucleus 5min, mod/total cytosol 0min,
    mod/total cytosol 5min, total/lib nucleus 0min, total/lib nucleus 5min, total/lib cytosol 0min,
    total/lib cytosol 5min

    name    size    0   5   0   5   0   5   0   5       (time points)
    1       120     1   2   1   3   1   1   1   1       (cluster 1 center values)
    2       86      3   15  2   3   1   2   1   1       (cluster 2 center values)
    3       237     4   4   1   6   2   4   3   6       (cluster 3 center values)
    """

    # getting genomic region's assigned clusters

    region_clusters = {}                                # hash storing region_name -> cluster_assigned
    in_file = open(cluster_results, "r")                # opening cluster results file
    while not next(in_file).startswith('"cluster"'):    # jumping to listing of genomic regions' cluster assignments
        pass
    for line in in_file:                                # iterating through listing of region's cluster assignments
        if line.startswith('"centers'): break               # stopping iteration if end of list is reached
        else:
            fields = line.strip().split("\t")               # getting region name and cluster assigned
            region_name = fields[0][1:-1]                   # getting region name
            cluster = fields[1]                             # getting cluster
            region_clusters[region_name] = cluster          # storing region_name -> cluster to hash
    in_file.close()                                     # closing clustering results file

    # getting clusters' sizes (number of genomic regions assigned)

    cluster_sizes = {}                                  # hash storing cluster -> cluster_size
    in_file = open(cluster_results, "r")                # opening cluster results file
    while not next(in_file).startswith('"size"'):       # jumping to listing of cluster sizes
        pass
    for line in in_file:                                # iterating through listing of cluster sizes
        fields = line.strip().split("\t")                   # getting cluster name and cluster size
        cluster = fields[0][1:-1]                           # cluster name
        size = fields[1]                                    # cluster size
        cluster_sizes[cluster] = size                       # storing cluster -> cluster_size
    in_file.close()                                     # closing clustering results file

    # getting genomic region's and cluster center ratios

    region_ratios = table_from_cluster_results(cluster_results, '"', '"cluster"', time_series)          # region ratios
    cluster_ratios = table_from_cluster_results(cluster_results, '"centers', '"totss"', time_series)    # cluster ratios

    # writing results to output files

    regions_outfile = open(regions_out, "w")        # opening regions' ratio values output file
    clusters_outfile = open(clusters_out, "w")      # opening clusters' center values output file

    time_points_nu = "\t".join(["nu_" + data_description + "_" + str(t) for t in time_series])  # nu ratio time series
    time_points_cy = "\t".join(["cy_" + data_description + "_" + str(t) for t in time_series])  # cy ratio time series
    time_points_series = "\t".join([time_points_nu, time_points_cy])        # time series

    regions_outfile.write("name\tcluster\t" + time_points_series + "\n")    # regions: first row
    clusters_outfile.write("name\tsize\t" + time_points_series + "\n")      # clusters: first row

    for gr in region_clusters:      # iterating through genomic regions, writing region ratio values

        gr_entry = gr + "\t" + region_clusters[gr] + "\t"               # genomic region name and cluster assigned
        gr_entry += "\t".join(region_ratios[0][gr]) + "\t" + "\t".join(region_ratios[1][gr]) + "\n"     # ratios
        regions_outfile.write(gr_entry)                                 # writing entry to ouput file

    for c in cluster_sizes:         # iterating through clusters, writing cluster center values

        c_entry = c + "\t" + cluster_sizes[c] + "\t"                    # cluster name and size
        c_entry += "\t".join(cluster_ratios[0][c]) + "\t" + "\t".join(cluster_ratios[1][c]) + "\n"      # ratios
        clusters_outfile.write(c_entry)                                 # writing entry to ouput file

    # closing files and returning

    regions_outfile.close()
    clusters_outfile.close()
    return()

# --- data filtering ---

def bam_ambiguous(fastqfile, bamfile, outfile):
    """
    :param fastqfile:   input fastq.gz file containing original reads that have been mapped
    :param bamfile:     input BAM file containing all aligned reads
    :param outfile:     output file to store ambiguously mapped reads to (fastq.gz format)

    :return: void

    This function filters all ambiguously mapped reads from a BAM file and writes them to a fastq output file.
    """

    # initializing

    sam_sorted = bamfile + "_tmp_sorted.sam"        # file to store sorted BAM file (temporary, deleted in the end)
    samtools_sort = subprocess.Popen("samtools sort " + bamfile + " -n -o " + sam_sorted + " -O SAM",
                                     shell=True)    # sorting BAM file by read name (read ID)
    samtools_sort.wait()                            # waiting for subprocess to finish

    in_fastq = gzip.open(fastqfile, "rt")           # opening input fastq file (unzipped)
    in_sam = open(sam_sorted, "r")                  # opening sorted SAM file
    out_fastq = gzip.open(outfile, "wb")            # opening output fastq file

    ambiguous_ids = {}      # hash storing ambiguously mapped read IDs

    # iterating through BAM sorted file, finding ambiguously mapped reads
    # (all reads occurring more than once have multiple, ambiguous mappings)

    current_id = ""             # current read ID processed
    current_id_count = 0        # counting number of alignments of current read ID
    ambiguous_count = 0         # counting number of ambiguously mapped reads

    for line in in_sam:         # iterating through reads in sorted BAM file

        if line.startswith("@"): continue       # skipping header lines

        r_id = line.split("\t")[0]              # getting read ID

        if r_id == current_id:                  # read ID equals current ID processed:
            current_id_count += 1                   # increment read ID alignment counter
        else:                                   # new read ID encountered:
            if current_id_count > 1:                # former read ID has multiple alignments:
                ambiguous_ids[current_id] = 0           # store former read ID to ambiguous reads' hash (else pass)
                ambiguous_count += 1                    # incrementing ambiguously mapped reads counter
            current_id = r_id                       # resetting current read ID
            current_id_count = 1                    # resetting read ID count

    if current_id_count > 1:   # also processing last read ID in BAM file
        ambiguous_ids[current_id] = 0

    # iterating through input fastq file, writing all ambiguously mapped read to output fastq file

    for line in in_fastq:       # iterating through input fastq file
        if line.startswith("@"):            # beginning of a read's entry (else pass)
            r_id = line.split(" ")[0][1:]                   # getting read ID
            if r_id in ambiguous_ids:                       # checking if read ID is ambiguous
                out_fastq.write(line.encode() + next(in_fastq).encode() + next(in_fastq).encode() +
                                next(in_fastq).encode())    # writing to output fastq

    # removing temporary sorted BAM, closing files and returning

    in_fastq.close()
    in_sam.close()
    os.remove(sam_sorted)
    out_fastq.close()
    print(ambiguous_count)
    return()

def bam_unmapped(fastqfile, bamfile, outfile):
    """
    :param fastqfile:   input fastq.gz file containing original reads that have been mapped
    :param bamfile:     input BAM file containing all aligned reads
    :param outfile:     output file to store unmapped reads to (fastq.gz format)

    :return: void

    This function filters all unmapped reads from a BAM file and writes them to a fastq output file.
    """

    # initializing

    sam_unmapped = bamfile + "_tmp_unmapped.sam"    # file to store unmapped aligned reads (temporary)
    get_unmapped = subprocess.Popen("samtools view -h -f 4 -o " + sam_unmapped + " " + bamfile,
                                    shell=True)     # getting unmapped reads
    get_unmapped.wait()                             # waiting for samtools to finish

    in_fastq = gzip.open(fastqfile, "rt")           # opening input fastq file (unzipped)
    in_sam = open(sam_unmapped, "r")                # opening SAM file with unmapped reads
    out_fastq = gzip.open(outfile, "wb")            # opening output fastq file

    unmapped_ids = {}                       # hash storing unmapped read IDs

    # iterating through BAM file, finding unmapped reads' IDs

    unmapped_count = 0              # counting number of ambiguously mapped reads

    for line in in_sam:             # iterating through reads in BAM file

        if line.startswith("@"): continue   # skipping header lines

        unmapped_count += 1                 # increment unmapped reads' counter
        read_id = line.split("\t")[0]       # getting read ID
        unmapped_ids[read_id] = 0           # storing read's ID

    # iterating through input fastq file, writing all unmapped reads to output fastq file

    for line in in_fastq:           # iterating through input fastq file
        if line.startswith("@"):                    # beginning of a read's entry (else pass)
            r_id = line.split(" ")[0][1:]               # getting read ID
            if r_id in unmapped_ids:                    # checking if read ID is ambiguous
                out_fastq.write(line.encode() + next(in_fastq).encode() + next(in_fastq).encode() +
                                next(in_fastq).encode())    # writing to output fastq

    # removing temporary input SAM file, closing files and returning

    in_fastq.close()
    in_sam.close()
    os.remove(sam_unmapped)
    out_fastq.close()
    print(unmapped_count)
    return ()

def bam_overlapping(bamfile_1, bamfile_2):
    """
    :param bamfile_1:   BAM file 1 containing aligned reads
    :param bamfile_2:   BAM file 2 containing aligned reads

    :return: a list containing read IDs of reads contained in both BAM files

    This function filters all reads contained in both given BAM files and returns the corresponding read IDs.
    """

    spikein_refnames = {"chrS2": None, "chrS4": None, "chrS5": None,
                        "chrS8": None, "chrS9": None, "chrS12": None}

    bam_file_1 = pysam.AlignmentFile(bamfile_1, "rb")                   # opening first BAM file with aligned reads
    bam_file_1_ids = [read.query_name for read in bam_file_1.fetch()]   # getting first BAM file read IDs
    bam_file_1_ids = dict.fromkeys(bam_file_1_ids)                      # (converting list to hash)
    bam_file_1.close()                                                  # closing first BAM file

    bam_file_2 = pysam.AlignmentFile(bamfile_2, "rb")       # opening second BAM file with aligned reads
    overlapping_reads = []                                  # initialize list storing read IDs contained in both files
    for read in bam_file_2.fetch():                         # iterating through second BAM file reads

        # checking if read ID occurs in first BAM file, and if it is not mapped to spike-in sequences
        if (read.query_name in bam_file_1_ids) and (read.reference_name not in spikein_refnames):
            overlapping_reads.append(read.query_name)       # if so, storing to return list

    bam_file_2.close()              # closing second BAM file

    return(overlapping_reads)       # returning list containing read IDs that occur in both BAM files

def primary_annotations(feature_hash, relation_hash):
    """
    :param feature_hash:    a hash storing annotated features in the following format:
                            feature_ID -> [[[ref_chr, start, end, strand], ...], feature_type, name, database_crossref]
                            starting positions 1-based, including
                            ending positions 1-based, excluding
    :param relation_hash:   a hash storing features' parent-child relations in the following format:
                            feature_ID -> [parent_id, {child_ids}]

    :return:    a list containing one hash for each the antisense, sense and unspecified strands' primary annotations
                as:
                feature_ID -> [[[ref_chr, start, end, strand], ...], feature_type, name, database_crossref]
                starting positions 1-based, including
                ending positions 1-based, excluding

    This function collects all primary annotations (all features that don't have a parent feature). Antisense and sense
    annotated features are stored in separate hashes.
    """

    # collecting primary annotations' IDs (features without parent feature)

    primary_feature_ids = [k for k in relation_hash if relation_hash[k][0] == ""]

    # sorting primary annotations by strand orientation

    antisense_annotations = {}  # hash storing positions with antisense annotations
    sense_annotations = {}      # hash storing positios with sense annotations
    other_annotations = {}      # hash storing positions with unspecified strand orientation

    for id in primary_feature_ids:      # iterating through IDs of primary annotations

        id_entry = feature_hash[id]
        id_regions = id_entry[0]        # getting genomic positions of annotated feature
        id_type = id_entry[1]           # getting feature type
        id_name = id_entry[2]           # getting feature name
        id_xref = id_entry[3]           # getting feature ID

        for region in id_regions:       # iterating through genomic positions (feature may be spread across genome)

            strand = region[3]              # getting strand orientation

            # choosing annotation hash depending on strand orientation

            if strand == "-": selected_hash = antisense_annotations
            elif strand == "+": selected_hash = sense_annotations
            else: selected_hash = other_annotations

            # storing primary annotation to corresponding hash

            if id not in selected_hash:
                selected_hash[id] = [[region], id_type, id_name, id_xref]
            else:
                selected_hash[id][0].append(region)

    # returning

    return([antisense_annotations, sense_annotations, other_annotations])

def filter_regions(bamfile, bedfile, refname_to_chr=r_to_c_hg19):
    """
    :param bamfile:         input BAM file containing aligned reads
    :param bedfile:         BED file containing genomic regions for which reads are to be filtered
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate BED chromosome names to the reference sequences) (default: default hash)

    :return: a hash containing all genomic regions and the read IDs overlapping with these regions:
                hash keys:      genomic regions' names
                hash values:    lists containing the region's chromosome number, starting position (0-based, including),
                                ending position (0-based, excluding), strand orientation, and read IDs overlapping

    This function filters reads from a BAM file for regions as defined in a BED file, respecting read mapping
    orientation. Forwards mapped reads are only assigned to antisense (-) regions and vice versa (this is a consequence
    of the 3' sequencing approach)
    """

    # ..................................................................................................................
    # opening files and initializing
    # ..................................................................................................................

    bam_file = pysam.AlignmentFile(bamfile, "rb")       # opening BAM file with aligned reads
    bed_file = open(bedfile)                            # opening BED file with genomic regions defined
    reads_per_region = {}               # initializing hash with genomic regions and reads filtered
                                        #   (region_name -> [chr_number, region_start, region_end, strand, read_ids])

    # reference names are chromosomes: need to find genomic regions in chromosomes

    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]         # getting all reference sequence names
    for i, r in enumerate(refnames):                                # replacing reference sequence names with
        if r in refname_to_chr: refnames[i] = refname_to_chr[r]     # chromosome numbers, if given

    genomic_regions = {}                # initializing hash with genomic regions' nucleotide positions
    for r in refnames:                  #   (chr_number -> {strand -> {nt_position -> region_name}})
        genomic_regions[r] = {"+": {}, "-": {}}

    # ..................................................................................................................
    # storing genomic regions (single nucleotide positions) in nested hash
    # chr_number -> {strand -> {nt_position -> region_name}}
    # ..................................................................................................................

    for line in bed_file:                   # iterating though BED file entries

        fields = line.strip().split("\t")                   # getting BED file entry's single fields
        chr_number = fields[0][3:]                          # genomic region chromosome number
        if chr_number not in genomic_regions: continue      # if chromosome number is unknown, continue

        startpos = int(fields[1])                           # genomic region starting position, 0-based, including
        endpos = int(fields[2])                             # genomic region ending position, 0-based, excluding
        name = fields[3]                                    # genomic region name
        strand = fields[5]                                  # genomic region strand orientation
        all_pos = range(startpos, endpos)                   # all positions of genomic region

        for p in all_pos:                                                   # storing the genomic region's nucleotide
            genomic_regions[chr_number][strand][p] = name                   #   positions to the genomic regions' hash

        reads_per_region[name] = [chr_number, startpos, endpos, strand]     # initialize genomic region information to
                                                                            #  reads-per-region hash

    # ..................................................................................................................
    # filtering reads from BAM file
    # region_name -> [chr_number, region_start, region_end, strand, read_ids]
    # ..................................................................................................................

    for read in bam_file.fetch():   # iterating through reads

        read_id = read.query_name           # read ID
        ref_name = read.reference_name      # name of reference sequence mapped to read
        if ref_name in refname_to_chr:      # replacing reference name with chromosome number, if given
            ref_name = refname_to_chr[ref_name]
        read_o = "+" if read.is_reverse else "-"    # read mapping orientation (is reversed on purpose, due to 3'Seq)
        ref_start = read.reference_start            # reference alignment start, 0-based, including
        ref_end = read.reference_end                # reference alignment end, 0-based, excluding
        ref_positions = range(ref_start, ref_end)   # reference alignment nucleotide positions

        for p in ref_positions:                     # checking if read alignment overlaps with any genomic region
            if p in genomic_regions[ref_name][read_o]:              # overlap was found:
                gregion_name = genomic_regions[ref_name][read_o][p]     # getting genomic region name
                reads_per_region[gregion_name].append(read_id)          # storing read ID to that genomic region's list
                break                                                   # stop iteration; continue with next read

    # ..................................................................................................................
    # returning hash with genomic regions and filtered reads
    # ..................................................................................................................

    return(reads_per_region)

def filter_coverage(reads_per_region, coverage=100):
    """
    :param reads_per_region:    a hash storing genomic regions and read IDs overlapping with these regions:
                                    hash keys:      genomic regions' names
                                    hash values:    lists containing the region's chromosome number, starting position
                                                    (0-based, including), ending position (0-based, excluding), strand
                                                    orientation, and read IDs overlapping
    :param coverage:            minimum coverage a genomic region has to satisfy (default: 100)

    :return:    a hash containing only these genomic regions that hit the minimum coverage

    This function selects genomic regions according to a minimum coverage a region has to satisfy. Coverage here is
    defined as the number of reads assigned to a region.
    """

    return_hash = {}                # initializing return hash
    cov_adjusted = coverage + 4     # minimum length a genomic region's hash entry has to have

    for gr in reads_per_region:         # iterating through genomic regions

        if len(reads_per_region[gr]) >= cov_adjusted:    # checking if region hits minimum coverage
            return_hash[gr] = reads_per_region[gr]       # storing region's input hash entry to output return list

    return(return_hash)             # returning hash with selected genomic regions

def antisense_sense_transcripts(r_summary, reads_per_region):
    """
    :param r_summary:           read summary; a hash assigning each read's ID a list containing (in the following
                                order):
                                    reference sequence name mapped to
                                    read length
                                    list of positions of T>C conversions observed (1-based)
                                    list of positions of insertions (position after insertion event, 1-based)
                                    list of positions of deletions (position of deletion event, 1-based)
                                    hash assigning a nucleotide pairing its abundance within the read's mapping
    :param reads_per_region:    a hash storing genomic regions and the read IDs overlapping with these regions:
                                    hash keys:      genomic regions' names
                                    hash values:    lists containing the region's chromosome number, starting position
                                                    (0-based, including), ending position (0-based, excluding), strand
                                                    orientation, and read IDs overlapping

    :return:    a list containing three read summaries storing antisense, sense and unspecified orientation genomic
                regions' reads
                *** NOTE *** reference sequences in the read summaries are replaced by the name of the corresponding
                genomic region

    This function separates reads from antisense (-) and sense (+) genomic regions. Reads are filtered from a read
    summary, using a hash which assigns read IDs to single genomic regions. For both antisense and sense genomic
    regions, a new read summary is returned.
    """

    # initializing read summary hashes

    ref_antisense = {}
    ref_sense = {}
    ref_unspecified = {}

    # iterating through genomic regions, separating reads by orientation

    for region in reads_per_region:

        region_entry = reads_per_region[region]     # getting genomic region's entry
        region_strand = region_entry[3]             # getting genomic region's strand orientation
        region_reads = region_entry[4:]             # getting genomic region's reads

        if (region_strand == "."):              # CASE: genomic region strand orientation not specified
            for read_id in region_reads:                # iterating through reads
                read_entry = r_summary[read_id]             # getting read's summary entry
                read_entry[0] = region                      # replacing reference sequence with genomic region
                ref_unspecified[read_id] = read_entry       # read to unspecified read summary

        elif (region_strand == "-"):            # CASE: genomic region strand orientation is antisense
            for read_id in region_reads:                # iterating through reads
                read_entry = r_summary[read_id]                 # getting read's summary entry
                read_entry[0] = region                          # replacing reference sequence with genomic region
                ref_antisense[read_id] = r_summary[read_id]     # adding all reads to forward read summary

        elif (region_strand == "+"):            # CASE: genomic region strand orientation is sense
            for read_id in region_reads:                # iterating through reads
                read_entry = r_summary[read_id]             # getting read's summary entry
                read_entry[0] = region                      # replacing reference sequence with genomic region
                ref_sense[read_id] = r_summary[read_id]     # adding all reads to reverse read summary

    # returning

    return([ref_antisense, ref_sense, ref_unspecified])

def filter_spare_reads(bamfile_in, rsummary_files, recorded_reads_bam, spare_reads_bam):
    """
    :param bamfile_in:          input BAM file containing aligned reads
    :param rsummary_files:      list containing one or more read summary files
    :param recorded_reads_out:  output BAM file to store aligned reads recorded within the read summaries
    :param spare_reads_bam:     output BAM file to store aligned reads not recorded within the read summaries

    :return:    void

    This function separates all aligned reads of a given BAM input file by reads which are recorded in any of the read
    summaries specified and those which are not. The two groups of reads are written to separate output BAM files.
    """

    # getting read IDs stored in read summary tables

    rsummary_ids = {}   # initialize hash storing read IDs occurring in the read summaries

    for rsum in rsummary_files:     # iterating through read summaries

        rsum_file = open(rsum, "r")     # opening read summary file
        for line in rsum_file:          # iterating through read summary
            if line.startswith("#"):        # finding read IDs in the read summary file
                read_id = line.strip()[1:]      # extracting read ID
                rsummary_ids[read_id] = None    # storing read ID
        rsum_file.close()               # closing read summary file

    # separating aligned reads from BAM file by read which are and which are not stored in the read summary tables

    bam_in = pysam.AlignmentFile(bamfile_in, "rb")                          # opening input BAM file
    bam_recorded_out = \
        pysam.AlignmentFile(recorded_reads_bam, "wb", template=bam_in)      # opening recorded reads' output BAM file
    bam_spare_out = \
        pysam.AlignmentFile(spare_reads_bam, "wb", template=bam_in)         # opening spare reads' output BAM file

    for read in bam_in.fetch():             # iterating through input BAM file reads

        read_id = read.query_name               # getting read ID
        if read_id in rsummary_ids:             # if read ID occurs in any read summary,
            bam_recorded_out.write(read)            # writing aligned read to recorded reads' output BAM file
        else:                                   # else,
            bam_spare_out.write(read)               # writing aligned read to spare reads' output BAM file

    # closing files and returning

    bam_recorded_out.close()    # closing output BAM file
    bam_spare_out.close()       # closing output BAM file
    bam_in.close()              # closing input BAM file
    return()                    # returning

# Data Analysis ########################################################################################################

# --- coverage ---

def covered_positions(bamfile):
    """
    :param bamfile:     input BAM file containing aligned reads

    :return:    a list containing two hashes (one for forward and one for reversely mapped mRNAs, respectively); the
                hashes assign each reference sequence name another hash which assigns a position (1-based) within the
                reference the coverage (uncovered positions are not stored);

    This function evaluates covered positions for the reads stored in a BAM file. To do so, the position of the 3'-most
    nucleotide of a read is considered only (since these are the suggested sites for the 3'UTR end / pA start).
    """

    bam_file = pysam.AlignmentFile(bamfile, "rb")               # opening BAM for retrieving refseq names
    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]     # getting all reference sequence names
    bam_file.close()                                            # closing BAM file again

    coverage_forward = dict.fromkeys(refnames)      # initializing hash to store ref forward mapped reads' coverage
    coverage_reverse = dict.fromkeys(refnames)      # initializing hash to store ref reverse mapped reads' coverage
    for r in coverage_forward:      # using for loop since dict.fromkeys() assigns the same hash object to all keys
        coverage_forward[r] = {}
        coverage_reverse[r] = {}

    # Read orientation and 3'most nucleotide location:

    # 1) reference sequence is (+, forward) strand, if its sequence corresponds to the 5'-> 3' mRNA's sequence
    # 2) reference sequence is (-, reverse) strand, if its sequence encodes the 3'-> 5' mRNA sequence
    # 3) read is the reverse complement of the mRNA strand

    # read is mapped forward (mRNA is reverse: mRNA orientation is 3'-> 5'): 3'most base is alignment start
    # read is mapped reverse (mRNA is forward: mRNA orientation is 5'-> 3'): 3'most base is alignment end

    bam_file = pysam.AlignmentFile(bamfile, "rb")  # opening BAM file

    for read in bam_file.fetch():   # iterating through reads in BAM file
        aln_ref = read.reference_name       # alignment reference sequence name
        aln_reverse = read.is_reverse       # alignment orientation (TRUE: ref reverse, FALSE: ref forward)
        if aln_reverse:                     # CASE: read is mapped in reference reverse orientation
            aln_end = read.reference_end                    # alignment end on reference sequence, 1-based, including
            if aln_end in coverage_forward[aln_ref]:        # increment coverage counter ...
                coverage_forward[aln_ref][aln_end] += 1
            else:                                           # ... or initialize coverage counter with 1
                coverage_forward[aln_ref][aln_end] = 1
        else:                               # CASE: read is mapped in reference forward orientation
            aln_start = read.reference_start + 1            # alignment start on reference sequence, 1-based, including
            if aln_start in coverage_reverse[aln_ref]:      # increment coverage counter ...
                coverage_reverse[aln_ref][aln_start] += 1
            else:                                           # ... or initialize coverage counter with 1
                coverage_reverse[aln_ref][aln_start] = 1

    # closing BAM file and returning

    bam_file.close()
    return([coverage_forward, coverage_reverse])

def coverage_distribution(covpos_file, outfile):
    """
    :param covpos_file:     input file storing covered positions
    :param outfile:         output file storing coverage frequency

    :return: void

    This function computes the coverage distribution (coverage frequencies) from a given covered positions' input file.
    Coverage distribution is calculated over all covered positions recorded (it is NOT computed for each reference
    sequence and/or read mapping orientation separately).

    The input file should be formatted as following (listing forwards and reversely mapped reads separately, mRNA
    orientation being indicated by "f" and "r", respectively):

    #reference_sequence_name_1\n
    position_1\tcoverage\tmRNA_orientation\n
    position_2\tcoverage\tmRNA_orientation\n
    ...
    position_n\tcoverage\tmRNA_orientation\n
    #reference_sequence_name_2\n

    The output file is tab-separated, storing the coverage value in the first column and the coverage value frequency
    in the second column.
    ...
    """

    cp_hash = load_covered_positions(covpos_file)   # loading in covered positions
    cov_distribution = {}                           # initializing hash storing {coverage_value: frequency}
    out_file = open(outfile, "w")                   # opening output file

    for h in cp_hash:   # iterating through forwards and reverse covered positions' hashes
        for ref_name in h:  # iterating through reference sequences
            refhash = h[ref_name]   # getting sub-hash
            for p in refhash:       # iterating through reference sequence's covered positions
                cov = refhash[p]        # getting coverage
                if cov in cov_distribution:     # coverage value is stored in coverage distribution hash, ...
                    cov_distribution[cov] += 1      # ... increment frequency
                else:                           # coverage value is not stored in coverage distribution hash, ...
                    cov_distribution[cov] = 1       # ... initialize frequency with 1

    cov_values = list(cov_distribution.keys())  # getting all coverage values stored
    min_value = min(cov_values)                 # minimum coverage value stored
    max_value = max(cov_values)                 # maximum coverage value stored

    # iterating through range of coverage values, writing coverage distribution to output file
    for i in range(min_value, max_value + 1):
        if i in cov_distribution:   # if coverage value is stored, writing its frequency
            out_file.write(str(i) + "\t" + str(cov_distribution[i]) + "\n")
        else:                       # if coverage value is not stored, writing frequency 0
            out_file.write(str(i) + "\t0\n")

    out_file.close()    # closing output file
    return()            # returning

def coverage_profile(cov_pos, radius=1000):
    """
    :param cov_pos:     a list containing two hashes (for forward and reversely mapped mRNAs, respectively); the hashes
                        assign each reference sequence name another hash which assigns a position within the reference
                        the coverage
    :param radius:      number of positions both upstream and downstream to be evaluated for the coverage profile
                        (default: 1000)

    :return:    a list containing two hashes (for forward and reversely mapped mRNAs, respectively); the hashes assign
                each reference sequence name a list which stores the coverage profile; the smallest list index (0)
                represents the downstream-most nucleotide position, the largest list index (radius * 2 + 1) the
                upstream-most nucleotide and the middle list index (radius) the central covered nucleotide position

    ....................................................................................................................

    This function computes a coverage profile from a reference sequence's covered positions, for each reference sequence
    and for forward and reversely mapped mRNAs separately. A coverage profile reflects the covered positions within a
    certain window around a central covered position, which is evaluated for each covered position as central position.
    The local profile of a single central position is a binary vector, 0 encoding for another position to not be covered
    and 1 otherwise. Local profiles are summed up to obtain the final profile.
    The window size can be specified by parameter 'radius', which defines the number of positions to be evaluated both
    upstream and downstream of the central position.
    """

    # ..................................................................................................................
    # initialize
    # ..................................................................................................................

    [covpos_forward, covpos_reverse] = cov_pos          # forwards and reversely mapped mRNAs' hashes of covered pos.
    ref_names_forward = list(covpos_forward.keys())     # forwards mapped mRNAs' evaluated reference sequences
    ref_names_reverse = list(covpos_reverse.keys())     # reversely mapped mRNAs' evaluated reference sequences

    # initialize profile hashes with reference sequence names as keys and lists of 0 counts as values (initialize
    # lists of 0 counts with for loop, since dict.fromkeys() assigns one and the same list object to all keys)

    window_size = radius * 2 + 1    # size of nucleotide window to be evaluated (2 * radius + central position)
    profile_forward = dict.fromkeys(ref_names_forward)
    profile_reverse = dict.fromkeys(ref_names_reverse)
    for r in profile_forward:
        profile_forward[r] = [0 for i in range(0, window_size)]
    for r in profile_reverse:
        profile_reverse[r] = [0 for i in range(0, window_size)]

    # ..................................................................................................................
    # calculate coverage profiles
    # for forward and reversely mapped mRNAs separately, for each reference sequence
    # ..................................................................................................................

    ref_names_number = str(len(ref_names_forward) + len(ref_names_reverse))     # number of references to be processed
    progress_counter = 0                                                        # reference progress counter

    for [cp_hash, profile] in zip([covpos_forward, covpos_reverse],
                                  [profile_forward, profile_reverse]):  # get covered positions' and profile hashes

        for ref_name in cp_hash:                                        # iterating through reference sequence names

            # initializing ---------------------------------------------------------------------------------------------

            progress_counter += 1                                                               # increment progress
            print("processing reference " + str(progress_counter) + " of " + ref_names_number)  # print progress

            # initialize
            ref_name_coverage = cp_hash[ref_name]       # coverage hash for the current reference sequence
            positions = list(cp_hash[ref_name].keys())  # covered positions for the current reference sequence
            positions.sort()                            # sorting covered positions in increasing order
            current_pos = positions[0]              # initialize current central position with first covered position
            current_profile = \
                [0 for i in range(0, window_size)]  # initialize current central position's profile with 0 counts

            # only process upstream positions, no downstream positions for the first covered position
            for i, p in enumerate(range(current_pos, current_pos + radius + 1)):
                if p in ref_name_coverage:             # if a position in the window is covered, ...
                    current_profile[radius + i] += 1   # ... update profile (set position's value to 1)

            # update overall profile with current profile
            for i in range(0, window_size):
                profile[ref_name][i] += current_profile[i]

            # iterate through reference sequence's covered positions,
            # constructing coverage profile ----------------------------------------------------------------------------

            for next_pos in positions[1:]:  # iterating through covered positions

                # initialize
                last_pos = current_pos              # storing last position
                current_pos = next_pos              # setting next position as current position
                distance = current_pos - last_pos   # distance between last position and next position

                # construct current profile

                if distance <= 2 * radius:      # CASE: current and last position's profiles overlap

                    recycled_profile = current_profile[distance:]       # recycle profile for last position
                    new_profile = [0 for i in range(distance)]          # initialize remaining profile with 0
                    for i, p in enumerate(range(
                                    last_pos + radius + 1, current_pos + radius + 1)):  # construct remaining profile
                        if p in ref_name_coverage:                      # if a position in the window is covered, ...
                            new_profile[i] += 1                         # ... update profile (set position's value to 1)
                    current_profile = recycled_profile + new_profile    # re-assemble complete profile

                else:                           # CASE: current and last position's profiles don't overlap

                    current_profile = [0 for i in range(0, window_size)]                # resetting current profile
                    for i, p in enumerate(range(
                                    current_pos - radius, current_pos + radius + 1)):   # construct profile
                        if p in ref_name_coverage:          # if a position in the window is covered, ...
                            current_profile[i] += 1         # ... update profile (set position's value to 1)

                # update overall profile with current profile
                for i in range(0, window_size):
                    profile[ref_name][i] += current_profile[i]

    # returning

    return([profile_forward, profile_reverse])

def coverage_convolution(cov_pos, radius=50, eval_length=10 ** 6):
    """
    :param cov_pos:         a list containing two hashes (for forward and reversely mapped mRNAs, respectively); the
                            hashes assign each reference sequence name another hash which assigns a position within the
                            reference the coverage
    :param radius:          nucleotide radius used for convolution (default: 50)
    :param eval_length:     length of a nucleotide stretch to be evaluated at a time (this is to avoid loading in all
                            data at once) (default: 10**6)

    :return:    a list storing two sublists, which again store two hashes:
                sublist 1 stores convolution value hashes, and sublist 2 stores the corresponding center position hashes
                hash 1 refers to forwards mapped and hash 2 refers to reversely mapped mRNAs
                a hash assigns a reference sequence name its list of convolution values or corresponding center position
                    values (1-based), respectively

    ....................................................................................................................

    This function convolves coverage across nucleotide positions using a triangular function. The width of the triangle
    is defined by a nucleotide radius considered around a central position, which is 2 * radius + 1. Additionally to
    the convolved coverage values, the corresponding central positions are recorded.
    """

    window_size = radius * 2 + 1                                # calculating window size

    filter_weights = [i for i in range(1, radius + 1)] + \
                     [radius + 1] + \
                     [i for i in reversed(range(1, radius + 1))]        # initializing filter weights
    filter_integral = sum(filter_weights)                               # getting the integral of the filter weights
    filter_sequence = [i / filter_integral for i in filter_weights]     # normalizing filter weights to integral of 1

    convolution_hash = [{}, {}]         # (for forwards and reversely mapped mRNAs) ref_name -> convolution_values
    center_positions_hash = [{}, {}]    # (for forwards and reversely mapped mRNAs) ref_name -> center_position

    for o, cp_hash in enumerate(cov_pos):       # iterating through mRNA orientation (forwards and reverse)
        for ref in cp_hash:                         # iterating through reference sequences

            # initializing .............................................................................................

            convolution_hash[o][ref] = []           # initialize list of convolution values
            center_positions_hash[o][ref] = []      # initialize list of center positions

            ref_hash = cp_hash[ref]                 # getting reference sequences' covered positions hash
            ref_cp = list(ref_hash.keys())          # getting covered positions
            ref_cp.sort()                           # sorting covered positions

            start_pos = ref_cp[0] - window_size     # start position (append empty window before first covered position)
            end_pos = ref_cp[-1] + 1 + window_size  # end position (append empty window after last covered position)
            pos_range = end_pos - start_pos + 1     # get length of nucleotide sequence from start to end position

            # starting (1-based, including), ending (1-based, excluding) and center positions (1-based)
            # of evaluated nucleotide stretches ........................................................................

            # ... start positions: beginning of the next window evaluated (beginning of last window evaluated + 1,
            # old_start + evaluation_length - window_size + 1)
            eval_starts = [start_pos + i * (eval_length - window_size + 1)
                           for i in range(0, int(math.floor(pos_range / (eval_length - window_size + 1))))]
            # ... end positions: start position + evaluation_length, the very end position is appended manually
            eval_ends = [i + eval_length for i in eval_starts[:-1] if (i + eval_length) < end_pos] + \
                        [end_pos]
            eval_starts = eval_starts[0:len(eval_ends)]
            # ... center positions: window beginning + radius, for each window beginning (from start position to
            # beginning of last window evaluated, from start to end - window_size + 1)
            center_positions = [[p + radius for p in range(start, end - window_size + 1)]
                                for start, end in zip(eval_starts, eval_ends)]

            # processing nucleotide stretches ..........................................................................

            zero_stretch = False        # boolean indicating whether a convolution 0-stretch is encountered

            for n, (start, end) in enumerate(zip(eval_starts, eval_ends)):  # iterating through nucleotide stretches

                # --- filling coverage values --

                coverage_values = [0 for i in range(start, end)]    # initialize coverage values with 0
                for i, pos in enumerate(range(start, end)):         # iterate through nucleotide positions

                    if pos in ref_hash:                                 # if position is covered, ...
                        coverage_values[i] = ref_hash[pos]              # ... setting corresponding coverage value

                # --- convolution ---

                convolution_values = numpy.convolve(coverage_values, filter_sequence, mode="valid")

                # --- filter out convolution 0-values ---

                conv_filtered = []                  # convolution values filtered
                centers_filtered = []               # center positions of convolution values filtered
                centers_all = center_positions[n]   # all center positions (of current nucleotide stretch)

                for i, v in enumerate(convolution_values):      # iterating through convolution values

                    if v == 0:                      # CASE: convolution value is 0
                        if not zero_stretch:            # first 0 value of a 0 stretch:
                            zero_stretch = True                         # marking start of a 0 stretch
                            conv_filtered.append(0)                     # storing first 0 value of a stretch
                            centers_filtered.append(centers_all[i])     # storing corresponding center position
                        else:                           # not the first 0 value of a 0 stretch:
                            pass                                        # omitting all other 0 values

                    else:                           # CASE: convolution value is not 0
                        zero_stretch = False                            # marking a currently non-0 stretch
                        conv_filtered.append(v)                         # storing convolution value
                        centers_filtered.append(centers_all[i])         # storing corresponding center position

                # --- storing ---

                convolution_hash[o][ref] = convolution_hash[o][ref] + conv_filtered
                center_positions_hash[o][ref] = center_positions_hash[o][ref] + centers_filtered

    # returning

    return([convolution_hash, center_positions_hash])

# --- annotation ---

def bed_counts(input_files, bedfile, refname_to_chr=r_to_c_hg19):
    """
    :param input_files:         list of input BAM files
    :param bedfile:             input BED file storing genomic regions
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number (needed to
                                associate BED chromosome names to the reference sequences) (default: default hash)

    :return: a hash assigning each BED region ID a list storing the region's read counts for the single input BAM files

    This function counts the number of reads that map to the genomic regions as defined in a BED file, for one or more
    given BAM files separately.
    """

    # initialize hash storing {region_id: [counts_per_bam]}

    bam_region_counts = {}              # initialize empty hash

    bed_file = open(bedfile, "r")       # opening BED file
    for line in bed_file:               # iterating through BED file regions
        region_id = line.strip().split("\t")[3]     # getting region ID
        bam_region_counts[region_id] = []           # initialize region ID entry
    bed_file.close()                    # closing BED file

    # counting BAM read counts per BED region

    for bam_file in input_files:        # iterate through BAM files
        region_info = filter_regions(
            bam_file, bedfile, refname_to_chr=refname_to_chr)       # getting BAM reads per region
        for rid in bam_region_counts:                               # iterating through region IDs
            if rid not in region_info:                                  # CASE: region ID has no information record
                bam_region_counts[rid].append(0)                            # append counts = 0
            else:                                                       # CASE: region ID has a information record
                r_counts = len(region_info[rid]) - 4                        # counting reads
                bam_region_counts[rid].append(str(r_counts))                # storing read counts

    # returning

    return(bam_region_counts)

def bed_annotation_gff(bedfile, annotation_file, refname_to_chr=r_to_c_hg19):
    """
    :param bedfile:             BED file storing genomic regions
    :param annotation_file:     GFF3 file storing genome annotations
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number
                                (default: default hash)

    :return:    a hash assigning each BED region ID another hash, assigning each primary feature's ID a list storing the
                feature's type, name, database cross references, and a hash storing each secondary feature's (sub-
                feature of the primary feature) type, name and database cross references.
                If no annotation is available, an empty hash is assigned.
                BED_ID -> {
                    feature_ID -> [type, name, db_crossref, {sub_feature_ID -> [type, name, db_crossref]}]
                }

    This function annotates all regions as defined in a given BED file with the help of annotated features as recorded
    in a given GFF3 file. Each region is assigned its associated primary features' IDs, types, names, and database cross
    references. Furthermore, all derived features' (sub-features of the primary feature) IDs, types and names are
    collected aswell.
    """

    # ..................................................................................................................
    # initializing BED annotation hash, storing BED genomic regions in nested hash
    # ..................................................................................................................

    bed_annotation_hash = {}                # initializing hash to store BED annotations to

    bed_file = open(bedfile)                # opening BED file
    for line in bed_file:                   # iterating though BED file entries
        name = line.strip().split("\t")[3]      # region identifier
        bed_annotation_hash[name] = {}          # initializing region's annotation hash entry
    bed_file.close()                        # closing BED file

    # ..................................................................................................................
    # loading in BED regions and annotations
    # ..................................................................................................................

    # loading in BED regions (chr_number -> [[start1, end1, id1], [...], ...])
    # 0-based!

    bed_regions_hash = load_bed_regions(bedfile)

    # loading annotated features from GFF3 annotation file
    # 1-based!

    annotation_tables = load_gff3_annotation(annotation_file)   # loading annotation information
    feature_table = annotation_tables[0]                        # hash storing annotated features
    parent_child_table = annotation_tables[1]                   # hash storing features' parent-child relations

    # filtering for primary annotations and sorting by strand orientation

    prim_annotations = \
        primary_annotations(feature_table, parent_child_table)  # primary features, sorted by strand orientation
    minus_annotations = prim_annotations[0]                     # hash for (-) strand annotated features
    plus_annotations = prim_annotations[1]                      # hash for (+) strand annotated features

    # ..................................................................................................................
    # iterating through BED regions, finding and assigning matching annotations
    # ..................................................................................................................

    for strand, region_hash, ann_hash in [
        ["-", bed_regions_hash[0], minus_annotations],
        ["+", bed_regions_hash[1], plus_annotations]]:      # strand-specific iteration
        for ref_name in region_hash:                            # iterating through reference sequences

            # building index for BED regions {position -> [region_ids]}
            ref_region_index = bed_index(bed_regions_hash, strand, ref_name)

            # iterating through primary annotations --------------------------------------------------------------------

            for f_id in ann_hash:        # iterating through primary annotations (primary features)

                f_entry = ann_hash[f_id]            # primary annotation's entry
                f_regions = f_entry[0]              # primary annotation's regions
                f_type = f_entry[1]                 # primary annotation's type
                f_name = f_entry[2]                 # primary annotation's name
                f_xref = f_entry[3]                 # primary annotation's cross references

                if "match" in f_type: continue      # skipping alignments (non-informative)

                for region in f_regions:            # iterating through annotation's region, checking if region matches

                    ref_chr = region[0]                     # annotation's reference sequence
                    if ref_chr in refname_to_chr:           # assign chromosome number (if possible)
                        ref_chr = refname_to_chr[ref_chr]
                    if ref_chr != ref_name: continue        # skipping annotations from non-matching references

                    for p in range(region[1], region[2]):   # iterating through annotated region positions,
                        p_zb = p - 1                            # (transforming position to 0-based scale)
                        if p_zb in ref_region_index:            # check if positions are found in BED region index

                            assigned_rids = ref_region_index[p_zb]      # getting overlapping BED region IDs
                            for rid in assigned_rids:                   # assigning annotation to BED region IDs
                                bed_annotation_hash[rid][f_id] = [f_type, f_name, f_xref, {}]

                # iterating through secondary annotations --------------------------------------------------------------

                f_children = parent_child_table[f_id][1]    # primary feature ID's child IDs (secondary feature IDs)
                for c_id in f_children:                     # iterating through secondary annotations

                    c_entry = feature_table[c_id]   # secondary feature's entry
                    c_regions = c_entry[0]          # secondary feature's regions

                    for region in c_regions:        # iterating through annotation's region, checking if region matches

                        ref_chr = region[0]                     # annotation's reference sequence
                        if ref_chr in refname_to_chr:           # assign chromosome number (if possible)
                            ref_chr = refname_to_chr[ref_chr]
                        if ref_chr != ref_name: continue        # skipping annotations from non-matching references

                        for p in range(region[1], region[2]):   # iterating through annotated region positions,
                            p_zb = p - 1                            # (transforming position to 0-based scale)
                            if p_zb in ref_region_index:            # check if positions are found in BED region index

                                assigned_rids = ref_region_index[p_zb]      # getting overlapping BED region IDs
                                for rid in assigned_rids:                   # assigning annotation to BED region IDs
                                    bed_annotation_hash[rid][f_id][3][c_id] = c_entry[1:]

            # garbage collection ---------------------------------------------------------------------------------------

            del ref_region_index
            gc.collect()

    # ..................................................................................................................
    # garbage collection and returning
    # ..................................................................................................................

    del bed_regions_hash
    del annotation_tables, feature_table, parent_child_table, prim_annotations, minus_annotations, plus_annotations
    gc.collect()
    return(bed_annotation_hash)

def bed_annotation_bed(bedfile, annotation_file):
    """
    :param bedfile:             BED file storing genomic regions
    :param annotation_file:     BED file storing annotated regions

    :return:    a hash assigning each BED region ID another hash storing all annotated regions' IDs that overlap with
                the BED region

    This function annotates all regions as defined in a given BED file with the help of annotated features as recorded
    in another BED file. Each region is assigned its associated annotated regions' IDs.
    """

    # initializing BED annotation hash .................................................................................

    bed_annotation_hash = {}    # initializing hash to store BED annotations to

    bed_file = open(bedfile)    # opening genomic regions' BED file
    for line in bed_file:       # iterating though file entries
        name = line.strip().split("\t")[3]      # region identifier
        bed_annotation_hash[name] = {}          # initializing region's annotation hash entry
    bed_file.close()            # closing BED file

    # loading in BED regions {chr_number -> [[start1, end1, id1], [...], ...]} .........................................

    bed_regions_hash = load_bed_regions(bedfile)            # BED genomic regions of interest
    annotation_hash = load_bed_regions(annotation_file)     # BED annotated regions

    # iterating through BED genomic regions, searching for annotations .................................................

    for strand, regions_hash, ann_hash in [
        ["-", bed_regions_hash[0], annotation_hash[0]],
        ["+", bed_regions_hash[1], annotation_hash[1]]]:    # iterating through strand orientation
        for ref_name in regions_hash:                           # iterating through reference sequences

            # skipping references sequences without annotations
            if ref_name not in ann_hash: continue

            # building BED indices {position -> [region_ids]}
            ref_regions_index = bed_index(bed_regions_hash, strand, ref_name)
            ref_ann_index = bed_index(annotation_hash, strand, ref_name)

            # iterating through annotated positions, finding matching BED genomic regions
            for p in ref_ann_index:         # iterating through annotated positions
                if p in ref_regions_index:      # finding matching BED genomic regions

                    ann_ids = ref_ann_index[p]      # get annotated regions' IDs
                    r_ids = ref_regions_index[p]    # get genomic regions' IDs
                    for a_id in ann_ids:            # assign annotated region IDs ...
                        for r_id in r_ids:              # ... to genomic region IDs
                            bed_annotation_hash[r_id][a_id] = None

            # garbage collection
            del ref_regions_index, ref_ann_index
            gc.collect()

    # closing BED file, garbage collection and returning ...............................................................

    bed_file.close()
    del bed_regions_hash, annotation_hash
    gc.collect()
    return(bed_annotation_hash)

def annotation_statistics(infile, outfile):
    """
    :param infile:      input file storing genomic region annotation
    :param outfile:     output file to write annotation statistics to

    :return: void

    ....................................................................................................................

    This function calculates annotation statistics from a given input file storing genomic region annotations.

    The input file format should be as described in the following:

    Tab-separated, starting with one header-line marked by a leading '#' symbol.
    Columns in the following order:

    > region identifier
    > first BAM file descriptor
    ...
    > n-th BAM file descriptor
    > exonic region ID
    > intronic region ID
    > 3'UTR ID
    > 5'UTR ID
    > primary features' IDs
    > primary features' types
    > primary features' names
    > primary features' database cross-references
    > secondary features' IDs
    > secondary features' types
    > secondary features' names
    > secondary features' database cross references

    If annotation information are not available for a region, columns are filled with '.'.
    Within one column ...
    ... multiple exonic, intronic, 3'UTR or 5'UTR IDs are ';' separated (example: intronID1;intronID2)
    ... multiple primary feature information are ';' separated (example: ID1;ID2)
    ... multiple secondary feature information are '*' separated (example: ID1_1*ID1_2;ID2_1*ID2_2*ID2_3)
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    # loading annotation
    # {region_ID -> [   [BAM1_descriptor, BAM1_counts, BAM2_descriptor, BAM2_counts, ...],
    #                   [exon_ids], [intron_id], [3UTR_id], [5UTR_id],
    #                   {feature_ID -> [    type, name, db_crossref,
    #                                       {sub_feature_ID -> [type, name, db_crossref]}
    #                                   ]}
    # ]}

    annotation_hash = load_bed_annotation(infile)

    # checking the number of BAM files for which region read counts are recorded

    any_ann_entry = annotation_hash[next(iter(annotation_hash))]            # get any annotation entry
    n_bam_entries = len(any_ann_entry[0])                                   # number of BAM file entries
    n_bam_files = int(n_bam_entries / 2)                                    # number of BAM files
    bam_des = [any_ann_entry[0][i] for i in range(0, n_bam_entries, 2)]     # BAM file descriptors

    # initializing hash storing annotation statistics, assigning each annotation type either:
    #   ... a list storing one list per BAM file with region read counts (if BAM file counts are recorded)
    #   ... a counter for that annotation type
    # annotation types:
    #   'none': no annotation at all
    #   'element': no primary annotation but still some 'exon' (only exon), 'intron' (only intron), 'utr3' (only 3'UTR),
    #       'utr5' (only 5'UTR) or 'others' (any other combination) annotation
    #   'primary': primary annotation and 'exon' (only exon), 'intron' (only intron), 'utr3' (only 3'UTR), 'utr5' (only
    #       5'UTR), 'intron_plus' (intron + any other element), 'exon_3utr' (exon + 3'UTR), 'exon_5utr' (exon + 5'UTR),
    #       'others' (any other combination), 'none' (no additional annotation)
    #   'secondary': secondary annotation feature types

    def ini_val():      # initialize assignment values (list of empty lists or 0)
        return([[] for i in range(n_bam_files)] if n_bam_files else [0])
    stats_hash = {"none": ini_val(),
                  "element": {"exon": ini_val(), "intron": ini_val(), "utr3": ini_val(), "utr5": ini_val(),
                              "others": ini_val()},
                  "primary": {"exon": ini_val(), "intron": ini_val(), "utr3": ini_val(), "utr5": ini_val(),
                              "intron_plus": ini_val(), "exon_utr3": ini_val(), "exon_utr5": ini_val(),
                              "others": ini_val(), "none": ini_val()},
                  "secondary": {"none": ini_val()}}

    # ..................................................................................................................
    # reading annotations and storing region read counts
    # ..................................................................................................................

    for r_id in annotation_hash:    # iterating through region IDs

        rid_ann = annotation_hash[r_id]               # region ID's annotation entry

        rid_counts = [[rid_ann[0][i]] for i in range(1, n_bam_entries, 2)] \
            if n_bam_entries else [1]                 # region ID's BAM file read counts (if recorded, else 1)

        exon, intron, utr3, utr5, prim = [bool(i) for i in rid_ann[1:6]]    # boolean for single annotation types

        # CASE: no primary annotation ----------------------------------------------------------------------------------

        if not prim:

            if exon:                                # --- CASE: --- exon recorded
                if intron or utr3 or utr5:                  # --- CASE: --- any other combination
                    for i, c in enumerate(rid_counts):
                        stats_hash["element"]["others"][i] += c
                else:                                       # --- CASE: --- exon only
                    for i, c in enumerate(rid_counts):
                        stats_hash["element"]["exon"][i] += c

            elif intron:                            # --- CASE: --- intron recorded
                if utr3 or utr5:                            # --- CASE: --- any other combination
                    for i, c in enumerate(rid_counts):
                        stats_hash["element"]["others"][i] += c
                else:                                       # --- CASE: --- intron only
                    for i, c in enumerate(rid_counts):
                        stats_hash["element"]["intron"][i] += c

            elif utr3:                              # --- CASE: --- 3'UTR recorded
                if utr5:                                    # --- CASE: --- any other combination
                    for i, c in enumerate(rid_counts):
                        stats_hash["element"]["others"][i] += c
                else:                                       # --- CASE: --- 3'UTR only
                    for i, c in enumerate(rid_counts):
                        stats_hash["element"]["utr3"][i] += c

            elif utr5:                              # --- CASE: --- 5'UTR only
                for i, c in enumerate(rid_counts):
                    stats_hash["element"]["utr5"][i] += c

            else:                                   # --- CASE: --- no annotation at all
                for i, c in enumerate(rid_counts):
                    stats_hash["none"][i] += c

        # CASE: primary annotation -------------------------------------------------------------------------------------

        else:

            if intron:                              # --- CASE: --- intron recorded
                if exon or utr3 or utr5:                    # --- CASE: --- intron_plus combination
                    for i, c in enumerate(rid_counts):
                        stats_hash["primary"]["intron_plus"][i] += c
                else:                                       # --- CASE: --- intron only
                    for i, c in enumerate(rid_counts):
                        stats_hash["primary"]["intron"][i] += c

            elif exon:                              # --- CASE: --- exon recorded
                if not utr3:
                    if not utr5:                            # --- CASE: --- exon only
                        for i, c in enumerate(rid_counts):
                            stats_hash["primary"]["exon"][i] += c
                    else:                                   # --- CASE: --- exon_utr5 combination
                        for i, c in enumerate(rid_counts):
                            stats_hash["primary"]["exon_utr5"][i] += c
                else:
                    if not utr5:                            # --- CASE: --- exon_utr3 combination
                        for i, c in enumerate(rid_counts):
                            stats_hash["primary"]["exon_utr3"][i] += c
                    else:                                   # --- CASE: --- other combination
                        for i, c in enumerate(rid_counts):
                            stats_hash["primary"]["others"][i] += c

            elif utr3:                              # --- CASE: --- 3'UTR recorded
                if not utr5:                                # --- CASE: --- 3'UTR only
                    for i, c in enumerate(rid_counts):
                        stats_hash["primary"]["utr3"][i] += c
                else:                                       # --- CASE: --- other combination
                    for i, c in enumerate(rid_counts):
                        stats_hash["others"]["utr3"][i] += c

            elif utr5:                              # --- CASE: --- 5'UTR only
                for i, c in enumerate(rid_counts):
                    stats_hash["primary"]["utr5"][i] += c

            else:                                   # --- CASE: --- no additional annotation
                for i, c in enumerate(rid_counts):
                    stats_hash["primary"]["none"][i] += c

            # processing secondary annotations -------------------------------------------------------------------------

            sec_types = {}      # hash storing secondary annotation feature types

            for f_id in rid_ann[5]:             # iterating through primary feature IDs
                for c_id in rid_ann[5][f_id][3]:    # iterating through secondary feature IDs
                    sec_types[rid_ann[5][f_id][3][c_id][0]] = None  # storing secondary feature type

            if not sec_types:   # no secondary features recorded
                for i, c in enumerate(rid_counts):
                    stats_hash["secondary"]["none"][i] += c
            else:               # secondary features recorded
                for t in sec_types:     # iterating through secondary feature types
                    if t not in stats_hash["secondary"]:    # initialize feature type's hash entry (if not yet stored)
                        stats_hash["secondary"][t] = ini_val()
                    for i, c in enumerate(rid_counts):      # updating secondary feature type's hash entry
                        stats_hash["secondary"][t][i] += c

    # ..................................................................................................................
    # writing to output file and returning
    # ..................................................................................................................

    out_file = open(outfile, "w")       # opening output file

    # writing 'no annotation'

    out_file.write("regions with no annotation:\t")
    stats_info = stats_hash["none"]                             # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics
    out_file.write("\n\n")

    # write 'no primary but exon/intron/UTR annotation'

    out_file.write("regions with no primary gene annotation, but an Exon / Intron / 3'UTR / 5'UTR annotation:\n")

    out_file.write("exon only:\t")                              # exon only statistics
    stats_info = stats_hash["element"]["exon"]                  # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\nintron only:\t")                          # intron only statistics
    stats_info = stats_hash["element"]["intron"]                # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\n3'UTR only:\t")                           # 3'UTR only statistics
    stats_info = stats_hash["element"]["utr3"]                  # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\n5'UTR only:\t")                           # 5'UTR only statistics
    stats_info = stats_hash["element"]["utr5"]                  # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\nany combinations:\t")                     # combinations statistics
    stats_info = stats_hash["element"]["others"]                # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\n\n")

    # write 'primary annotation'

    out_file.write("regions with primary gene annotation\n")

    out_file.write("exon only:\t")                              # exon only statistics
    stats_info = stats_hash["primary"]["exon"]                  # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\nexon and 3'UTR:\t")                       # exon and 3'UTR statistics
    stats_info = stats_hash["primary"]["exon_utr3"]             # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\nexon and 5'UTR:\t")                       # exon and 5'UTR statistics
    stats_info = stats_hash["primary"]["exon_utr5"]             # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\nintron only:\t")                          # intron only statistics
    stats_info = stats_hash["primary"]["intron"]                # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\nintron and others:\t")                    # intron and others statistics
    stats_info = stats_hash["primary"]["intron_plus"]           # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\n3'UTR only:\t")                           # 3'UTR only statistics
    stats_info = stats_hash["primary"]["utr3"]                  # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\n5'UTR only:\t")                           # 5'UTR only statistics
    stats_info = stats_hash["primary"]["utr5"]                  # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\n3'UTR and 5'UTR:\t")                      # 3'UTR and 5'UTR statistics
    stats_info = stats_hash["primary"]["others"]                # getting statistics information
    write_annotation_statistics(out_file, stats_info, bam_des)  # writing statistics

    out_file.write("\n\n")

    # write secondary feature types

    out_file.write("feature types assigned to primary annotations:")

    for feature_type in stats_hash["secondary"]:        # iterating through secondary feature types
        out_file.write("\n" + feature_type + ":\t")
        stats_info = stats_hash["secondary"][feature_type]              # getting statistics information
        write_annotation_statistics(out_file, stats_info, bam_des)      # writing statistics

    out_file.write("\n")

    # closing output file and returning

    out_file.close()
    return()

# --- conversions and labeling ---

def editing_counts(r_summary, conversion, outfile):
    """
    :param r_summary:   read summary; a hash assigning each read's ID a list containing (in the following order):
                            reference sequence name mapped to
                            read length
                            list of positions of insertions (position after insertion event, 1-based)
                            list of positions of deletions (position of deletion event, 1-based)
                            hash assigning a nucleotide pairing its abundance within the read's mapping
    :param conversion:  nucleotide conversion to be investigated (nucleotide pairing as 'ReferenceNtReadNt')
    :param outfile:     output file to write editing counts to

    :return: void

    This function counts the number of reads carrying a specified nucleotide conversion per reference sequence (per
    genomic region) and writes the results to an output file (tab-separated, column 1: genomic region name, column 2:
    total number of reads, column 3: number of reads carrying a conversion).
    """

    edit_hash = {}    # initialize hash storing genomic_region -> [total reads, converted reads]
    out_file = open(outfile, "w")

    for read in r_summary:      # iterating through read summary's reads

        read_region = r_summary[read][0]                    # getting read's genomic region
        read_ntcounts = r_summary[read][4][conversion]      # getting read's nucleotide pairing counts
        read_status = 0                                     # storing if read is converted (0 for no, 1 for yes)
        if read_ntcounts > 0:
            read_status = 1

        if read_region in edit_hash:                        # genomic region already recorded: counting up
            edit_hash[read_region][0] += 1
            edit_hash[read_region][1] += read_status
        else:                                               # genomic region not yet recorded: initializing
            edit_hash[read_region] = [1, read_status]

    out_file.write("name\ttotal\tconverted\n")              # writing output file header line
    for gr in edit_hash:                                    # writing genomic regions' counts to output file
        out_file.write(gr + "\t" + str(edit_hash[gr][0]) + "\t" + edit_hash[gr][1] + "\n")

    return()        # returning

def editing_binning(editing_file, outfile_ratios, outfile_binning, coverage=50, resolution=100):
    """
    :param editing_file:        file storing genomic regions' converted read counts (tab separated, columns referring
                                to genomic region names, total counts, and converted counts)
    :param outfile_ratios:      output file to store genomic regions' converted reads ratios to
    :param outfile_binning:     output file to store converted reads ratios' bins to

    :param coverage:            minimum coverage a genomic region has to satisfy (default: 50)
    :param resolution:          binning resolution (number of bins) (default: 100)

    :return: void

    This function bins the editing ratios of a collection of genomic regions and writes both the computed ratios and
    the bins to separate output files (tab-separated, first column referring to gene region name or bin upper border,
    second column referring to converted reads ratios or bin counts, respectively).
    """

    edit_counts = {}                        # initialize hash to store editing counts to
    edit_infile = open(editing_file, "r")   # opening editing counts file
    next(edit_infile)                       # skipping header
    for line in edit_infile:                # iterating through editing infile
        entries = line.strip().split("\t")                                  # getting single entries
        edit_counts[entries[0]] = [float(entries[1]), float(entries[2])]    # storing entries to editing hash
    edit_infile.close()                     # closing editing counts file

    out_ratios = open(outfile_ratios, "w")      # opening ratios output file
    out_bins = open(outfile_binning, "w")       # opening bins output file
    out_ratios.write("name\tratio\n")           # writing header line (ratios)
    out_bins.write("bin\tcount\n")              # writing header line (counts)

    edit_ratios = {}            # initializing hash storing genomic_region -> converted ratio
    for gr in edit_counts:      # iterating through genomic regions, computing ratios
        gr_total = edit_counts[gr][0]               # total read counts
        if gr_total < coverage: continue            # skipping all genomic regions not fulfilling minimum coverage
        gr_converted = edit_counts[gr][1]           # converted read counts
        gr_edit_ratio = gr_converted / gr_total     # converted reads ratio
        edit_ratios[gr] = gr_edit_ratio             # storing converted ratio

    binstep = 1 / resolution                                    # bin step size for converted ratios
    bins = [binstep * i for i in range(1, resolution + 1)]      # bin upper limits
    bin_counts = [0 for i in range(0, resolution)]              # bin counts (list indices refer to bins)

    for gr in edit_ratios:                                      # iterating through ratios
        ratio_bin = int(math.floor(edit_ratios[gr] / binstep))      # computing ratio bin
        if ratio_bin == resolution:                                 # maximum ratio falls into bin (resolution + 1)
            ratio_bin -= 1
        bin_counts[ratio_bin] += 1                                  # incrementing bin count
        out_ratios.write(gr + "\t" + str(edit_ratios[gr]) + "\n")   # writing ratio to ratio output file

    for i, b in enumerate(bins):    # iterating through bins, writing counts to output file
        out_bins.write(str(b) + "\t" + str(bin_counts[i]) + "\n")

    out_ratios.close()      # closing ratio output file
    out_bins.close()        # closing bins output file
    return()                # returning

def snp_pos_binning(read_table_file, outfile, coverage=20, resolution=100):
    """
    :param read_table_file:     file storing read table
    :param outfile:             output file to store position-wise SNP ratios' bins to

    :param coverage:            minimum coverage a position has to satisfy (default: 20)
    :param resolution:          binning resolution (number of bins) (default: 100)

    :return: void

    This function computes position-wise SNP rates for all genomic regions stored within the input SNP table summary
    and bins the SNP rates (according to a given resolution). The bins and bin counts are written to an output file
    (tab-separated, first column referring to bin upper border, second column referring to bin counts).
    """

    rtable = load_read_table(read_table_file)                   # loading in read table
    out_file = open(outfile, "w")                               # opening output file
    binstep = 1 / resolution                                    # bin step size
    bins = [binstep * i for i in range(1, resolution + 1)]      # bin upper limits
    bin_counts = [0 for i in range(0, resolution)]              # bin counts (list indices refer to bins)

    # computing position-wise ratios, counting up bins

    for gr in rtable:                           # iterating though genomic regions
        gr_snp_table = rtable[gr][0]                        # getting genomic region's SNP table summary
        for p, total in enumerate(gr_snp_table[1][1:]):     # iterating through region's positions (skipping row name)

            if total < coverage: continue                       # skipping positions with insufficient coverage

            event = gr_snp_table[2][p + 1]                      # getting event counts (skipping row name)
            ratio = event/total                                 # computing ratio (event_counts/total_counts)
            ratio_bin = int(math.floor(ratio / binstep))        # computing ratio bin
            if ratio_bin == resolution:                         # maximum ratio falls into bin (resolution + 1)
                ratio_bin -= 1
            bin_counts[ratio_bin] += 1                          # incrementing bin count

    # writing bins to output file

    out_file.write("bin\tcount\n")              # writing header
    for b, bin_border in enumerate(bins):       # iterating through bins
        out_file.write(str(bin_border) + "\t" + str(bin_counts[b]) + "\n")  # writing bin counts

    # closing output file and returning

    out_file.close()
    return()

def conversion_counts(conversion_table):
    """
    :param conversion_table:    a nested list storing read conversions:
                                    row 1: potential reference sequence conversion position, other rows: reads,
                                    column 1: read names, other columns: potential reference sequence conversion
                                        positions,
                                    entries: None if the position is not covered, masked by a deletion or SNP
                                        correction, 0 if read has no conversion at that position (may have any other
                                        SNP, though), 1 otherwise

    :return:    a hash assigning each read ID a list containing the total number of potential conversion positions and
                the number of conversions observed

    This function counts the total number of potential conversion positions as well as the number of actual conversions
    of the reads stored within a conversion table.
    """

    read_counts = {}    # hash storing read_id -> [number of potential conversion positions, number of conversions]

    for read_record in conversion_table[1:]:    # skipping first row (only contains reference genome positions)
        read_id = read_record[0]                # getting read ID
        # number of potential conversion positions (all entries unequal None and excluding read name)
        total_counts = sum([1 for i in read_record if i != None]) - 1
        # number of conversions (all entries being 1)
        conv_counts = sum([1 for i in read_record if i == 1])
        # storing counts to hash
        read_counts[read_id] = [total_counts, conv_counts]

    return(read_counts)     # returning reads' counts

def efficiency_estimation(r_conversions, base_error, sequencing_error,
                          fix_efficiency=None, new_ini=None, iterations=100):
    """
    :param r_conversions:       a hash assigning each read ID a list containing the total number of potential conversion
                                positions and the number of conversions observed
    :param base_error:          probability of observing a conversion due to a sequencing error
    :param sequencing_error:    probability that the converted nucleotide is masked by a sequencing error

    :param fix_efficiency:      optional; a value to fix for conversion efficiency (as a consequence, only new/total
                                ratio will be estimated) (default: None
    :param new_ini:             optional; initial value for new/total ratio (default: None)
    :param iterations:          optional; number of iterations to use for efficiency estimation (default: 100)

    :return:    a list containing two sublists: the first sublist contains the sequence of estimated newly synthesized
                by total RNA ratios, the second sublist contains the sequence of estimated conversion efficiencies (the
                final estimations being the last list entries); in case no estimation could be made (labeled read counts
                too low), the list contains two 'None' values instead

    This function estimates the conversion efficiency from a collection of reads, using an EM algorithm that in parallel
    estimates the new/total ratio as well.
    """

    # ..................................................................................................................
    # collecting all combinations: (number of potential conversion positions and number of conversions) read counts
    # hash: number potential conversion pos. -> {number conversions -> number of reads}
    # ..................................................................................................................

    pos_conv_numbers = {}           # initializing hash
    for r in r_conversions:         # iterating through reads
        p = r_conversions[r][0]         # number of potential conversion positions
        c = r_conversions[r][1]         # number of conversions
        if p in pos_conv_numbers:       # number of positions already stored in hash
            if c in pos_conv_numbers[p]:    # number of conversions already stored in hash
                pos_conv_numbers[p][c] += 1     # incrementing number of reads in existing entry
            else:                           # number of conversions not yet stored in hash
                pos_conv_numbers[p][c] = 1      # setting new entry
        else:                           # number of positions not yet stored in hash
            pos_conv_numbers[p] = {c: 1}    # setting new entry

    # ..................................................................................................................
    # initial estimations of 'new' (ratio of newly synthesized RNA transcripts) and 'e' (conversion efficiency)
    # initialize 'new' with the ratio of labeled by total RNA transcripts
    # initialize 'e' with an efficiency estimation excluding base and sequencing errors, or base error
    # ..................................................................................................................

    conv_1 = 0      # number of reads with one conversion
    conv_2 = 0      # number of reads with two conversions
    conv_all = 0    # number of reads with one or more conversions (number of labeled reads)
    for p in pos_conv_numbers:      # iterating over hash storing combination counts, counting up labeled reads
        for c in pos_conv_numbers[p]:   # iterating over conversion numbers
            if c == 1:                      # conversion number 1
                conv_1 += pos_conv_numbers[p][c]    # incrementing counter for one-time converted reads
            elif c == 2:                    # conversion number 2
                conv_2 += pos_conv_numbers[p][c]    # incrementing counter for two-time converted reads
            elif c != 0:
                conv_all += pos_conv_numbers[p][c]  # incrementing counter for multiple-times converted reads
    conv_all += conv_1 + conv_2     # adding 1- and 2-converted read counts to one-or-more-converted read counts

    # initial values for 'e' and 'new'

    e = base_error                              # e: initialize with base error
    if fix_efficiency:                          # e: if given, fixed efficiency
        e = fix_efficiency
    elif (conv_1 != 0 and conv_2 != 0):         # e: else, ratio of 2 and 1-times converted, excluding base errors
        ratio = conv_2/conv_1                           # ratio of 2- and 1-times converted read counts
        npos = list(pos_conv_numbers.keys())[0] + 2     # number of potential conversion positions (anything > 1)
        e = ratio / ((npos - 1) / 2 + ratio)            # initial 'e' estimation

    new = conv_all / len(r_conversions)         # new: labeled/total transcript counts ratio
    if new_ini: new = new_ini                   # new: user-defined value

    # ..................................................................................................................
    # iteratively estimating 'new' and 'e' using EM algorithm
    # ..................................................................................................................

    new_estimation_sequence = [new]     # sequence of 'new' estimations during EM algorithm iterations
    e_estimation_sequence = [e]         # sequence of 'e' estimations during EM algorithm iterations

    for i in range(0, iterations):  # iteratively re-estimating

        # conversion probability at one position for newly synthesized transcripts
        a = e * (1 - sequencing_error - base_error) + base_error

        # P(H|O ; theta) probability of hidden variable (pre-existing or newly synthesized transcript) given observation
        # (read number of conversions and number of potential conversion positions)
        # hash: number of potential conversion positions -> {number of conversions -> P(H | O; theta)}

        prob_h0 = dict.fromkeys(pos_conv_numbers, {})       # initialize hash for H = 0 (pre-existing)
        prob_h1 = dict.fromkeys(pos_conv_numbers, {})       # initialize hash for H = 1 (newly synthesized)

        for p in pos_conv_numbers:      # iterating over number of potential conversion positions
            for c in pos_conv_numbers[p]:   # iterating over number of conversions
                combi = scipy.special.binom(p, c)       # combinations to draw 'c' out of 'p' (binomial coefficient)
                h0 = combi * (base_error ** c) * ((1 - base_error) ** (p - c)) * (1 - new)   # P(O|H)*P(H) for H = 0
                h1 = combi * (a ** c) * ((1 - a) ** (p - c)) * new                           # P(O|H)*P(H) for H = 1
                denominator = h0 + h1                                                        # P(O); denominator
                prob_h0[p][c] = h0 / denominator        # computing and storing P(H|O; theta) for H = 0
                prob_h1[p][c] = h1 / denominator        # computing and storing P(H|O; theta) for H = 1

        # sums over P(H|O; theta) needed for subsequent computations

        w0 = 0      # sum over all P(H=0|O=c; theta)
        w1 = 0      # sum over all P(H=1|O=c; theta)
        v = 0       # sum over all ( P(H=1|O=c; theta) * c )
        z = 0       # sum over all ( P(H=1|O=c; theta) * (p-c) )

        for p in pos_conv_numbers:  # iterating over number of potential conversion positions
            for c in pos_conv_numbers[p]:   # iterating over number of conversions
                counts = pos_conv_numbers[p][c]     # read counts
                h0 = prob_h0[p][c]                  # corresponding P(H=0|O=c)
                h1 = prob_h1[p][c]                  # corresponding P(H=1|O=c)
                w0 += counts * h0           # incrementing w0
                w1 += counts * h1           # incrementing w1
                v += counts * h1 * c        # incrementing v
                z += counts * h1 * (p - c)  # incrementing z

        # re-estimating 'new' and 'e', appending to estimation sequences

        new = w1 / (w0 + w1)
        if fix_efficiency:          # fixed conversion efficiency:
            e = fix_efficiency          # setting e to fixed value
        else:                       # conversion efficiency is to be estimated:
            if (v + z) == 0:            # no conversion efficiency estimation for can be made:
                return ([None, None])       # returning None values
            else:                       # conversion efficiency estimation possible:
                e = 1 / (1 - base_error - sequencing_error) * (v / (v + z) - base_error)

        new_estimation_sequence.append(new)
        e_estimation_sequence.append(e)

    # ..................................................................................................................
    # returning
    # ..................................................................................................................

    return([new_estimation_sequence, e_estimation_sequence])

def efficiency_median(cetable_file, cutoff):
    """
    :param cetable_file:    conversion efficiency table file
    :param cutoff:          minimum value of a conversion efficiency estimator to be considered (this is to filter out
                            estimates of genes with barely any signal, being biased towards 0)

    :return: the median conversion efficiency of the input conversion efficiency table

    This function computes the median conversion efficiency of a set of genes' conversion efficiency estimates as
    stored within an input table.
    """

    ce_table = load_conveff_table(cetable_file)     # loading conversion efficiency table
    ce_values = []                                  # list storing conversion efficiency estimates

    for g in ce_table:                          # iterating through genes listed in conversion efficiency table
        if ce_table[g][0] == None: continue         # skipping genes without efficiency estimate
        ce = ce_table[g][1][-1]                     # getting conversion efficiency
        if ce >= cutoff: ce_values.append(ce)       # collecting estimate if sufficiently large

    m = median(ce_values)       # computing median
    return(m)                   # returning

def efficiency_binning(conveff_table, outfile_ce, outfile_ratio, resolution=100):
    """
    :param conveff_table:   conversion efficiency table
    :param outfile_ce:      output file to store conversion efficiency bins to
    :param outfile_ratio:   output file to store ratios bins to

    :param resolution:      binning resolution (number of bins) (default: 100)

    :return: void

    This function bins the results of a conversion efficiency table (i.e. both the conversion efficiencies and newly
    synthesized by total transcript ratios) and writes the results to separate output files (tab-separated, the first
    column referring to upper bin limits and the second to bin counts, respectively).
    """

    n_vals = len(conveff_table)                                         # getting number of data points
    conveff_list = [conveff_table[gr][1][-1] if conveff_table[gr][1] else -1
                    for gr in conveff_table]                            # getting conversion efficiencies
    ratio_list = [conveff_table[gr][0][-1] if conveff_table[gr][0] else -1
                  for gr in conveff_table]                              # getting newly synthesized ratios
    out_file_ce = open(outfile_ce, "w")                                 # opening conversion efficiency output file
    out_file_ratio = open(outfile_ratio, "w")                           # opening ratio output file

    min_ce = min(conveff_list)                      # minimum conversion efficiency value
    max_ce = max(conveff_list)                      # maximum conversion efficiency value
    ce_span = abs(min_ce) + abs(max_ce)             # conversion efficiency value span
    min_ratio = min(ratio_list)                     # minimum ratio value
    max_ratio = max(ratio_list)                     # maximum ratio value
    ratio_span = abs(min_ratio) + abs(max_ratio)    # ratio value span

    ce_binstep = ce_span / resolution               # bin step size for conversion efficiencies
    ratio_binstep = ratio_span / resolution         # bin step size for ratios
    ce_bins = [min_ce + ce_binstep * i for i in range(1, resolution + 1)]        # bin upper limits for ce
    ratio_bins = [min_ce + ratio_binstep * i for i in range(1, resolution + 1)]  # bin upper limits for ratios

    ce_bin_counts = [0 for i in range(0, resolution)]       # list storing ce bin counts (list indices refer to bins)
    ratio_bin_counts = [0 for i in range(0, resolution)]    # list storing ratio bin counts (list indices refer to bins)

    for i in range(0, n_vals):              # iterating through data points
        ce = conveff_list[i]                                                # conversion efficiency
        r = ratio_list[i]                                                   # newly synthesized ratio
        ce_bin = int(math.floor((ce - min_ce) / ce_binstep))                # conversion efficiency bin
        ratio_bin = int(math.floor((r - min_ratio) / ratio_binstep))        # newly synthesized ratio bin
        if ce_bin == resolution:        # maximum conversion efficiency falls into bin (resolution + 1)
            ce_bin -= 1
        if ratio_bin == resolution:     # maximum ratio falls into bin (resolution + 1)
            ratio_bin -= 1
        ce_bin_counts[ce_bin] += 1                              # incrementing conversion efficiency bin counts
        ratio_bin_counts[ratio_bin] += 1                        # incrementing newly synthesized ratio bon counts

    for i in range(0, resolution):      # iterating through bins, writing bin counts to output file
        out_file_ce.write(str(ce_bins[i]) + "\t" + str(ce_bin_counts[i]) + "\n")
        out_file_ratio.write(str(ratio_bins[i]) + "\t" + str(ratio_bin_counts[i]) + "\n")

    out_file_ce.close()     # closing conversion efficiency output file
    out_file_ratio.close()  # closing ratio output file
    return()                # returning

# Data Table Computations ##############################################################################################

def covered_regions(covconv_file, outfile, radius=50, refname_to_chr=r_to_c_hg19):
    """
    :param covconv_file:    input file storing coverage convolution
    :param outfile:         output BED file to store covered regions
    :param radius:          radius around a coverage peak position to be defined as covered region (default: 50)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number
                            (default: default hash)

    :return: void

    This function defines covered regions and writes them to an output BED file. A covered region is defined based on
    a coverage peak and a nucleotide radius of fixed size around that peak.
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    # loading coverage convolution and corresponding center positions
    (conv_values, center_positions) = load_coverage_convolution(covconv_file)
    # hashes storing reference_name -> [coverage peak positions] for forwards and reversely mapped mRNAs
    peak_positions = [{}, {}]

    out_file = open(outfile, "w")               # opening output BED file

    # ..................................................................................................................
    # finding coverage peaks
    # ..................................................................................................................

    for cv_hash, cp_hash, peak_hash in zip(
            conv_values, center_positions, peak_positions):     # iterating through mRNA orientation
        for ref_name in cv_hash:                                    # iterating through reference sequences

            cv_values = cv_hash[ref_name]       # reference's convolution values
            center_pos = cp_hash[ref_name]      # reference's corresponding center positions
            peak_list = []                      # list storing reference's peak positions

            cv_diffs = numpy.diff(cv_values)    # computing change in convolved coverage to find peaks
            increase = False                    # boolean indicating whether coverage increases
            for i, d in enumerate(cv_diffs):    # iterating through coverage changes

                if d > 0:               # CASE: positive change
                    increase = True         # indicate that coverage increases
                elif d == 0:            # CASE: no change
                    pass                    # nothing to do
                else:                   # CASE: negative change
                    if increase:            # if coverage was increasing before, there is a coverage peak
                        peak_list.append(center_pos[i])     # storing corresponding nucleotide position
                    increase = False        # indicate that coverage is not increasing

            peak_hash[ref_name] = peak_list     # storing reference's peak positions

    # ..................................................................................................................
    # defining and writing covered regions
    # ..................................................................................................................

    for peak_hash, o, s in zip(peak_positions, ["+", "-"],
                               ["f", "r"]):                 # iterating through mRNA orientation
        for ref_name in peak_hash:                              # iterating through reference sequences
            for peak_pos in peak_hash[ref_name]:                    # iterating through peak positions
                chr_number = ref_name                                   # setting chromosome number to reference name,
                if ref_name in refname_to_chr:                          # unless a chromosome number is specified
                    chr_number = "chr" + refname_to_chr[ref_name]
                region_start = peak_pos - radius - 1                    # region start position (0-based, including)
                region_end = peak_pos + radius                          # region end position (0-based, excluding)
                id = "_".join([chr_number, str(region_start), str(region_end), s])  # region identifier
                # writing BED file entry
                out_file.write("\t".join([chr_number, str(region_start), str(region_end), id, "0", o]) + "\n")

    # ..................................................................................................................
    # closing files and returning
    # ..................................................................................................................

    out_file.close()
    return()

def bed_annotations(bed_infile, annotation_file, bed_exons, bed_introns, bed_3utrs, bed_5utrs, outfile,
                    bam_files=None, bam_des=None, refname_to_chr=r_to_c_hg19):
    """
    :param bed_infile:          BED file storing genomic regions
    :param annotation_file:     GFF3 file storing genome annotations
    :param bed_exons:           BED file storing exon regions
    :param bed_introns:         BED file storing intronic regions
    :param bed_3utrs:           BED file storing 3'UTR regions
    :param bed_5utrs:           BED file storing 5'UTR regions
    :param outfile:             output file to store annotation to

    :param bam_files:           optional; list of BAM files for which the BED genomic region read counts are to be
                                calculated (default: None)
    :param bam_des:             optional; list of BAM file descriptions to use within the output file header
                                (default: ascending numbers)
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number
                                (default: default hash)

    :return: void

    ....................................................................................................................

    This function annotates all regions as defined within a given BED file. Annotation includes a region's read counts
    within one or more BAM files, IDs of associated exonic, intronic, 3' and 5' UTR regions (each of which defined in a
    separate BED file), and IDs, types, names, and database cross references of associated primary and secondary
    features.

    The output file is formatted as described in the following:

    Tab-separated, starting with one header-line marked by a leading '#' symbol.
    Columns in the following order:

    > region identifier
    > first BAM file descriptor
    ...
    > n-th BAM file descriptor
    > exonic region ID
    > intronic region ID
    > 3'UTR ID
    > 5'UTR ID
    > primary features' IDs
    > primary features' types
    > primary features' names
    > primary features' database cross-references
    > secondary features' IDs
    > secondary features' types
    > secondary features' names
    > secondary features' database cross references

    If annotation information are not available for a region, columns are filled with '.'.
    Within one column ...
    ... multiple exonic, intronic, 3'UTR or 5'UTR IDs are ';' separated (example: intronID1;intronID2)
    ... multiple primary feature information are ';' separated (example: ID1;ID2)
    ... multiple secondary feature information are '*' separated (example: ID1_1*ID1_2;ID2_1*ID2_2*ID2_3)
    """

    # initializing hash storing all annotations ........................................................................
    # {region_ID -> [   [BAM1_descriptor, BAM1_counts, BAM2_descriptor, BAM2_counts, ...],
    #                   [exon_ids], [intron_id], [3UTR_id], [5UTR_id],
    #                   {feature_ID -> [    type, name, db_crossref,
    #                                       {sub_feature_ID -> [type, name, db_crossref]}
    #                                   ]}
    # ]}

    all_annotations = {}                # initialize hash storing all annotations
    in_file = open(bed_infile, "r")     # opening BED file
    for line in in_file:                # iterating through BED file regions
        region_id = line.strip().split("\t")[3]                 # getting region ID
        all_annotations[region_id] = [[], [], [], [], [], {}]   # initialize region ID entry
    in_file.close()                     # closing BED file

    # getting BED genomic region read counts ...........................................................................

    if bam_files:
        # define BAM file descriptors as increasing numbers (if not given)
        if not bam_des:
            bam_des = [str(i) for i in range(1, len(bam_files) + 1)]
        # getting BED genomic regions' read counts
        region_counts = bed_counts(bam_files, bed_infile,refname_to_chr=refname_to_chr)
        # storing regions' read counts to annotation hash
        for rid in region_counts:
            for des, c in zip(bam_des, region_counts[rid]):
                all_annotations[rid][0] += [des, str(c)]
        # garbage collection
        del region_counts
        gc.collect()

    # getting and merging exon and intron annotations ..................................................................

    exon_annotations = bed_annotation_bed(bed_infile, bed_exons)
    intron_annotations = bed_annotation_bed(bed_infile, bed_introns)

    for r_id in exon_annotations:       # iterating through genomic region IDs

        # getting annotated region IDs (set to "" if no annotation is available)
        exon_ids = [a_id for a_id in exon_annotations[r_id]] if exon_annotations[r_id] else []
        intron_ids = [a_id for a_id in intron_annotations[r_id]] if intron_annotations[r_id] else []

        # merging
        all_annotations[r_id][1:3] = [exon_ids, intron_ids]

    del exon_annotations, intron_annotations    # intermediate garbage collection
    gc.collect()

    # getting and merging 3'UTR and 5'UTR annotations ..................................................................

    utr3_annotations = bed_annotation_bed(bed_infile, bed_3utrs)
    utr5_annotations = bed_annotation_bed(bed_infile, bed_5utrs)

    for r_id in all_annotations:        # iterating through genomic region IDs

        # getting annotated region IDs (set to "" if no annotation is available)
        utr3_ids = [a_id for a_id in utr3_annotations[r_id]] if utr3_annotations[r_id] else []
        utr5_ids = [a_id for a_id in utr5_annotations[r_id]] if utr5_annotations[r_id] else []
        # merging
        all_annotations[r_id][3:5] = [utr3_ids, utr5_ids]

    del utr3_annotations, utr5_annotations
    gc.collect()

    # getting and merging GFF3 annotations .............................................................................

    gff3_annotations = bed_annotation_gff(bed_infile, annotation_file, refname_to_chr=refname_to_chr)

    for r_id in all_annotations:  # iterating through genomic region IDs

        # merging
        all_annotations[r_id][5] = gff3_annotations[r_id]

    del gff3_annotations
    gc.collect()

    # writing annotations and returning ................................................................................

    write_bed_annotation(all_annotations, outfile)
    return()

def read_table(mode, bamfile, reffile, bedfile, snpfile=None, coverage=1, refname_to_chr=r_to_c_hg19):
    """
    :param mode:        mode to be executed; choose 'snp_all' to compute SNP tables reporting all kind of conversions,
                        choose 'snp_exclusive' to compute SNP tables reporting all but labeling-specific conversions,
                        choose 'conv' to compute conversion tables only, use 'snp_all_summary', 'snp_exclusive_summary',
                        'conv_summary' to compute summarized tables that report SNP and total read counts for each
                        position, choose 'both' or 'both_summary' to compute both a SNP (snp_exclusive) and a conversion
                        table
    :param bamfile:     input BAM file containing aligned reads
    :param reffile:     fasta file containing all reference sequences; an index of the same file needs to be located in
                        the same directory, having the same filename with the additional extension '.fai'
    :param bedfile:     BED file containing genomic regions for which reads are to be filtered

    :param snpfile:         gzipped vcf file containing SNP positions of the reference (default: None)
    :param coverage:        minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate the SNP sites / A-to-I sites to the reference sequences) (default: default hash)

    :return:    a hash assigning each genomic region's name a list containing a SNP and a conversion table (implemented
                as nested lists, the first index corresponding to rows and the second to columns):
                    SNP table: row 1: genomic region's reference sequence positions (0-based), other rows: reads,
                        column 1: read names, other columns: reference sequence positions, entries: None if the position
                        is not covered by the read, 0 if read has no SNP at that position, SNP nucleotide otherwise
                    conversion table: row 1: potential reference sequence conversion position, other rows: reads,
                        column 1: read names, other columns: potential reference sequence conversion positions,
                        entries: None if the position is not covered, masked by a deletion or SNP correction, 0 if read
                        has no conversion at that position (may have any other SNP, though), 1 otherwise

    This function records SNPs as well as conversions for each genomic region, resolved by reads. SNPs and conversions
    are stored in separate tables, with rows corresponding to reads and columns to potential SNP or conversion posi-
    tions, respectively, in the reference genome. The first row always contains reference sequence positions, the
    first column always contains read IDs (in case of summary tables, the read IDs are set to 'total_counts' and
    'event_counts'). The first row + first column entry is set to None. SNPs are denoted as None if the position is not
    covered by the read, 0 in case no SNP is present, and the SNP nucleotide otherwise (deletions are denoted with 'D',
    insertions with 'I' in additional columns at the end of the table). Conversions are stored as None in case that the
    position is not covered by the read, masked by a deletion or SNP correction, 0 in case no conversion is present,
    and as 1 otherwise.

    Optionally, conversion positions can be corrected for SNP positions by providing a SNP vcf file. In this case, all
    positions listed to be prone to SNPs will be reported as None in the table.
    Optionally, a minimum coverage can be specified by which genomic regions are filtered; only regions hitting this
    minimum coverage will be processed and included into the read table hash.
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    print("\nINITIALIZING")

    reference = pysam.FastaFile(reffile)                        # opening reference file
    bam_file = pysam.AlignmentFile(bamfile, "rb")               # opening BAM for retrieving refseqs
    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]     # getting all reference sequence names
    n_refs = len(refnames)                                      # getting number of reference sequences
    bam_file.close()                                            # closing BAM file again

    # filtering reads and genomic regions ------------------------------------------------------------------------------

    genomic_regions = filter_regions(           # filtering reads by genomic region
        bamfile, bedfile, refname_to_chr=refname_to_chr)
    for gr in list(genomic_regions.keys()):     # filtering genomic regions by strand orientation
        if genomic_regions[gr][3] == ".":       # (filter out unspecified strand orientations)
            del genomic_regions[gr]
    if coverage:        # filtering genomic regions by coverage (keeping regions hitting minimum coverage threshold)
        genomic_regions = filter_coverage(genomic_regions, coverage=coverage)
    read_regions = {}   # reverting genomic regions' hash: read_ID -> genomic region and strand orientation
    for gr in genomic_regions:
        for read_id in genomic_regions[gr][4:]:
            read_regions[read_id] = [gr, genomic_regions[gr][3]]

    # initialize SNP and conversion table related hashes ---------------------------------------------------------------

    # SNP table related hashes

    regions_minpos = {}     # hash storing region_name -> minimum position covered (0-based, including)
    regions_maxpos = {}     # hash storing region_name -> maximum position covered (0-based, excluding)
    region_insinfo = {}     # hash storing region_name -> {position -> {read_ids with ins at that position}}
    read_snpinfo = {}       # hash storing read_ID -> [starting position (0-based, incl), SNPs per position]

    if mode.startswith("snp") or mode.startswith("both"):
        for gr in genomic_regions:
            regions_minpos[gr] = sys.maxsize
            regions_maxpos[gr] = 0
            region_insinfo[gr] = {}

    # conversion table related hashes

    regions_convpos = {}    # hash storing region_name -> hash of potential conversion positions (0-based)
    snps = {}               # hash storing ChrNumber_ChrPos of reported SNPs
    read_convinfo = {}      # hash storing read_ID -> [starting position (0-based, incl), conv per position]

    if mode.startswith("conv") or mode.startswith("both"):
        for gr in genomic_regions:
            regions_convpos[gr] = {}

    # initializing hash storing read tables

    read_tables = dict.fromkeys([k for k in genomic_regions])

    # ..................................................................................................................
    # collecting reads' SNP and conversion information
    # (per reference sequence, since SNP sites' hash may need too much memory for all eference sequences at once)
    # ..................................................................................................................

    for (i, current_ref) in enumerate(refnames):

        # initializing -------------------------------------------------------------------------------------------------

        print("\nPROCESSING REFERENCE " + current_ref + " (" + str(i + 1) + "/" + str(n_refs) + ")\n")

        bam_file = pysam.AlignmentFile(bamfile, "rb")           # re-opening BAM file
        if refname_to_chr and (current_ref in refname_to_chr):  # reference sequence chromosome number (is set to the
            current_chr = refname_to_chr[current_ref]           # reference sequence name if no number is given by
        else:                                                   # 'refname_to_chr')
            current_chr = current_ref

        # loading reported SNPs' hash (if needed)
        if mode.startswith("conv") or mode.startswith("both"):
            print("loading SNPs ...")
            if snpfile and current_ref in refname_to_chr:
                snps = load_snps(snpfile, current_chr)
            print(" done")

        # iterating through all reads in the BAM file,
        # getting reads' SNPs and conversions per positions ------------------------------------------------------------

        print("computing SNP and conversion tables ...")
        for read in bam_file.fetch():

            # getting read and reference metadata

            read_id = read.query_name                       # read ID
            if read_id not in read_regions: continue            # skipping reads not mapped to any genomic region
            if read.reference_name != current_ref: continue     # skipping reads not mapped to current refseq
            read_region = read_regions[read_id][0]          # read genomic region
            region_strand = read_regions[read_id][1]        # read genomic region's strand (ref forward or reverse)
            if region_strand == ".": continue                   # skipping reads with unspecified strand

            read_seq = read.query_alignment_sequence        # read sequence (without hard/soft-clipped parts)
            ref_start = read.reference_start                # reference alignment start, 0-based, including
            ref_end = read.reference_end                    # reference alignment end, 0-based, excluding
            ref_seq = reference.fetch(                      # reference alignment seq (soft-masked lower-case patched)
                reference=read.reference_name, start=ref_start, end=ref_end).upper()

            update_region_borders(read_region, regions_minpos, regions_maxpos,
                                  ref_start, ref_end)                       # updating genomic regions' min and max pos
            (converted_nt, conversion_nt) = labeling_ntpair(region_strand)  # defining conversion (labeling) nt pair

            # getting and storing read SNP and conversion information

            construction_mode = mode.split("_summary")[0]
            read_snpconv = process_cigar_snps_conv(construction_mode, read.cigartuples, converted_nt, conversion_nt,
                                                   read_id, read_seq, current_chr, read_region, ref_seq, ref_start,
                                                   region_insinfo, regions_convpos, snps=snps)  # read SNP and conv info

            read_snpinfo[read_id] = [ref_start] + read_snpconv[0]
            read_convinfo[read_id] = [ref_start] + read_snpconv[1]

        # closing BAM file (end of file)

        bam_file.close()
        print(" done")

    # ..................................................................................................................
    # constructing genomic regions' read tables and storing them
    # ..................................................................................................................

    print("\nCONSTRUCTING SNP AND CONVERSION TABLES ...")

    for gr in genomic_regions:  # iterating through all genomic regions, building read tables

        gr_snptable = []            # initializing SNP table for current genomic region
        gr_convtable = []           # initializing conversion table for current genomic region

        if mode.startswith("snp") or mode.startswith("both"):       # --- SNP table ---
            gr_snptable = construct_snp_table(genomic_regions[gr][4:], read_snpinfo, region_insinfo[gr],
                                              regions_minpos[gr], regions_maxpos[gr])
        if mode.startswith("conv") or mode.startswith("both"):      # --- conversion table ---
            gr_convtable = construct_conv_table(genomic_regions[gr][4:], read_convinfo,
                                                list(regions_convpos[gr].keys()))
        if mode.endswith("summary"):                                # --- summary tables ---
            construction_mode = mode.split("_")[0]
            table_summaries = construct_snp_conv_summary_table(construction_mode, gr_snptable, gr_convtable)
            gr_snptable = table_summaries[0]
            gr_convtable = table_summaries[1]

        read_tables[gr] = [gr_snptable, gr_convtable]               # adding genomic region's table to return hash

    print(" done")

    # ..................................................................................................................
    # returning hash storing genomic regions' read tables
    # ..................................................................................................................

    return(read_tables)

def read_summary(bamfile, reffile, snpfile=None, aifile=None, editfile=None,
                 refname_to_chr=r_to_c_hg19):
    """
    :param bamfile: input BAM file containing aligned reads
    :param reffile: fasta file containing all reference sequences; an index of the same file needs to be located in the
                    same directory, having the same filename with the additional extension '.fai'

    :param snpfile:         gzipped vcf file containing SNP positions of the reference (default: None)
    :param aifile:          REDIportal file containing A-to-I editing positions of the reference (default: None)
    :param editfile:        editing sites file storing editing sites as 'ChrNumber_ChrPos' (default: None)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate the SNP sites / A-to-I sites to the reference sequences) (default: default hash)

    :return: a hash assigning each read's ID a list containing (in the following order):
                reference sequence name mapped to
                read length
                list of positions of insertions (position after insertion event, 1-based)
                list of positions of deletions (position of deletion event, 1-based)
                confusion table (hash) assigning a nucleotide pairing its abundance within the read's mapping

    ....................................................................................................................

    This function collects read information for all mapped reads within a given BAM file. It is designed for the
    analysis of metabolically labeled RNA transcripts.
    The information retrieved are: reference sequence name mapped to, read length, positions of insertions (position
    before insertion event, w.r.t. reference, 1-based), positions of deletions (position before deletion event, w.r.t.
    reference, 1-based), confusion table listing the number of nucleotide pairings as well as insertions and deletions
    within the read.

    The confusion table is a hash assigning each nucleotide pairing its abundance within the read's mapping. The first
    nucleotide refers to the reference, the second to the read. Pairings with 'N' are also contained (i.e., all 'NX'
    and 'XN' for some standard nucleotide X). 'XI' and 'XD' refer to insertions and deletions, respectively, where X is
    the next reference nucleotide after the insertion event or, in the case of deletions, the first nucleotide in the
    reference deleted. The counts refer to single nucleotides inserted or deleted rather than to insertion or deletion
    events.

    Optionally, confusion table counts (nucleotide pairing counts) can be corrected for SNPs and A-to-I editing sites by
    providing corresponding vcf files. Here, only single nucleotide exchanges are considered. All reference positions
    being reported to contain a SNP or being an A-to-I editing site are excluded from counting (are not comprised in
    confusion table). Also, position-wise SNP rates can be computed from an input read table file, and a maximum editing
    rate (SNP rate) can be defined for a nucleotide position to be considered.

    Read information is returned as hash, assigning a read's ID a list with its information.
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    print("\nINITIALIZING")
    bam_file = pysam.AlignmentFile(bamfile, "rb")                               # opening BAM for retrieving read IDs
    read_information = dict.fromkeys([r.query_name for r in bam_file.fetch()])  # initialize hash keys with read IDs
    bam_file.close()                                                            # closing BAM (end of file)

    bam_file = pysam.AlignmentFile(bamfile, "rb")                               # opening BAM for retrieving refseqs
    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]                     # getting all reference sequence names
    n_refs = len(refnames)                                                      # number of all reference sequences
    bam_file.close()                                                            # closing BAM file again

    reference = pysam.FastaFile(reffile)                                        # opening reference fasta file

    # ..................................................................................................................
    # retrieving read information per reference sequence, since SNP and A-to-I editing sites' hash may need too much
    # memory for all reference sequences at once
    # ..................................................................................................................

    for (i, current_ref) in enumerate(refnames):

        # initializing -------------------------------------------------------------------------------------------------

        print("\nPROCESSING REFERENCE " + current_ref + " (" + str(i + 1) + "/" + str(n_refs) + ")\n")
        bam_file = pysam.AlignmentFile(bamfile, "rb")           # re-opening BAM file
        if refname_to_chr and (current_ref in refname_to_chr):  # reference sequence chromosome number (is set to the
            current_chr = refname_to_chr[current_ref]           # reference sequence name if no number is given by
        else:                                                   # 'refname_to_chr')
            current_chr = current_ref
        snps_ai = dict()                                        # hash storing ChrNumber_ChrPos of SNP and A-to-I sites

        # loading SNPs from SNP vcf file (if file is given)

        if snpfile and current_ref in refname_to_chr:
            print("loading SNPs ...")
            snps_ai = load_snps(snpfile, current_chr)
            print(" done")

        # loading A-to-I editing sites (if file is given)

        if aifile and current_ref in refname_to_chr:
            print("loading A-to-I sites ...")
            ai_hash = load_ai_sites(aifile, current_chr)
            snps_ai.update(ai_hash)
            print(" done")

        # loading nt sites with SNP rates above tolerated threshold (if file is given)

        if editfile:
            print("loading potential editing sites ...")
            edit_hash = load_edit_sites(editfile, current_chr)
            snps_ai.update(edit_hash)
            print(" done")

        # iterating through all reads in the BAM file, getting: --------------------------------------------------------
        # read ID, reference, read length, read T->C / insertion / deletion positions, read confusion table

        print("computing confusion table ...")
        for read in bam_file.fetch():

            # getting read and reference metadata

            ref_name = read.reference_name              # name of reference sequence mapped to read

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue                                # (skipping unmapped and ambiguously mapped reads)
            if ref_name != current_ref: continue        # (skipping reads not mapped to current reference sequence)

            read_id = read.query_name                   # read ID
            read_seq = read.query_alignment_sequence    # read sequence (without hard/soft-clipped parts)
            read_length = read.query_alignment_length   # read length (without hard/soft-clipped parts)
            read_cigar = read.cigartuples               # read cigar string (as tuples)

            ref_start = read.reference_start            # reference alignment start, 0-based, including
            ref_end = read.reference_end                # reference alignment end, 0-based, excluding
            ref_seq = reference.fetch(reference=ref_name, start=ref_start, end=ref_end)     # alignment ref. seq.
            ref_seq = ref_seq.upper()                   # converting soft-masked lower-case nucleotides to upper case

            # computing and storing confusion table

            read_cigar_info = process_cigar_ntpairs(ref_seq, read_seq, read_cigar, current_chr, ref_start,
                                                    snps_ai)
            insertions = read_cigar_info[0]
            deletions = read_cigar_info[1]
            conf_table = read_cigar_info[2]

            read_information[read_id] = [ref_name, read_length, insertions, deletions, conf_table]

        bam_file.close()                                # closing BAM file (end of file)
        print(" done")

    # ..................................................................................................................
    # returning hash storing all read information as read_id -> list of read information
    # ..................................................................................................................

    return(read_information)

def read_summary_statistics(r_summary):
    """
    :param r_summary:   a hash assigning each read's ID a list containing information on that read (as returned by the
                        function 'read_summary'):
                            reference sequence name mapped to
                            read length
                            list of positions of insertions (position before insertion event, 1-based)
                            list of positions of deletions (position before deletion event, 1-based)
                            hash assigning a nucleotide pairing its abundance within the read's mapping

    :return: a list containing (in the following order):
                total number of reads
                hash: fragment length -> ratio of reads of that length
                relative position of insertions
                relative position of deletions
                a hash storing for each nucleotide pairing:
                    number of pairings -> % reads with that many pairings
                    % reference nucleotides paired with this read nucleotide (only considering SNPs, not indels)
                    % reads that have at least one such pairing and additionally contain at least one insertion
                    % reads that have at least one such pairing and additionally contain at least one deletion
                    relative position of additional insertions
                    relative position of additional deletions
                    hash: read length -> % reads of that length (of reads with at least one such pairing)

    ....................................................................................................................

    This function computes a summary statistics on read information collected for all mapped reads within a given BAM
    file as returned by the function 'read_summary'. Just as 'read_summary', it is designed for the analysis of
    metabolically labeled RNA transcripts.

    The statistics includes:

    total number of reads
    fragment length distribution as hash: fragment length -> % reads of that length
    relative position of insertions
    relative position of deletions
    a hash storing for each nucleotide pairing:
        number of pairings distribution as hash: number of pairings -> % reads with that many pairings
        % reference nucleotides being paired with the corresponding read nucleotide (only considering SNPs, not indels)
        % reads that have at least one such pairing and additionally contain at least one insertion
        % reads that have at least one such pairing and additionally contain at least one deletion
        relative position of additional insertions
        relative position of additional deletions
        fragment length distr. (of reads with at least one such pairing) as hash: read length -> % reads of that length

    The percentage of reference nucleotides being paired with the corresponding read nucleotide are computed from the
    total number of reference nucleotides being paired to any other nucleotide, i.e. insertion and deletion events are
    not taken into account here.
    Relative positions of insertions and deletions are computed as the average relative position within an read and
    then averaged over all reads. Position of an inserted or deleted nucleotide is the position of the last matched
    nucleotide before the insertion or deletion event. Each nucleotide is considered individually rather than taking
    into account only the insertion or deletion event; therefore, larger indels have higher impact on an read's average
    relative indel position.
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    total_reads = len(r_summary)            # total number of reads
    read_lengths = {}                       # hash storing read length -> % reads of that length
    ins_positions = 0                       # relative position of insertions
    del_positions = 0                       # relative position of deletions

    # stats to be recorded for each nucleotide pairing:
        # number of pairings -> % reads with that many pairings,
        # % reference nucleotides being paired with the corresponding read nucleotide (doesn't consider indels),
        # % reads that additionally contain insertions, % reads that additionally contain deletions,
        # average relative position of additional insertion, average relative position of additional deletion,
        # read length -> % reads of that length
    ntpair_stats = {
        "AA": [{0: 0}, 0, 0, 0, 0, 0, {}], "AT": [{0: 0}, 0, 0, 0, 0, 0, {}], "AG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "AC": [{0: 0}, 0, 0, 0, 0, 0, {}], "AI": [{0: 0}, 0, 0, 0, 0, 0, {}], "AD": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "TA": [{0: 0}, 0, 0, 0, 0, 0, {}], "TT": [{0: 0}, 0, 0, 0, 0, 0, {}], "TG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "TC": [{0: 0}, 0, 0, 0, 0, 0, {}], "TI": [{0: 0}, 0, 0, 0, 0, 0, {}], "TD": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "GA": [{0: 0}, 0, 0, 0, 0, 0, {}], "GT": [{0: 0}, 0, 0, 0, 0, 0, {}], "GG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "GC": [{0: 0}, 0, 0, 0, 0, 0, {}], "GI": [{0: 0}, 0, 0, 0, 0, 0, {}], "GD": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "CA": [{0: 0}, 0, 0, 0, 0, 0, {}], "CT": [{0: 0}, 0, 0, 0, 0, 0, {}], "CG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "CC": [{0: 0}, 0, 0, 0, 0, 0, {}], "CI": [{0: 0}, 0, 0, 0, 0, 0, {}], "CD": [{0: 0}, 0, 0, 0, 0, 0, {}]}
    nt_pairs = ["AA", "AT", "AG", "AC", "AI", "AD", "TA", "TT", "TG", "TC", "TI", "TD", "GA", "GT", "GG", "GC",
                "GI", "GD", "CA", "CT", "CG", "CC", "CI", "CD"]     # keys for nucleotide pairing statistics hash

    # ..................................................................................................................
    # iterating through read summary, filling read summary statistics
    # ..................................................................................................................

    ins = 0     # total number of insertion affected reads (needed for averaging)
    de = 0      # total number of deletion affected reads (needed for averaging)

    for read_id in r_summary:    # iterating through reads

        read_records = r_summary[read_id]    # all information on read
        read_len = read_records[1]                  # length of read
        read_ins = read_records[2]                  # read insertion positions
        read_del = read_records[3]                  # read deletion positions
        read_conf = read_records[4]                 # read confusion table

        # filling read length distribution -----------------------------------------------------------------------------

        if read_len in read_lengths: read_lengths[read_len] += 1
        else: read_lengths[read_len] = 1

        # counting up insertions and relative insertion positions ------------------------------------------------------

        is_ins = 0          # variable indicating whether insertions are present
        pos_ins = 0         # initialize average relative position of read's insertions
        if read_ins:
            is_ins = 1                                              # tagging that insertions are present in read
            pos_ins = sum(read_ins) / (len(read_ins) * read_len)    # compute average relative pos. of read's insertions

        ins += is_ins               # incrementing insertion reads counter
        ins_positions += pos_ins    # adding average relative position to insertion relative position counter

        # counting up deletions and relative deletion positions --------------------------------------------------------

        is_del = 0          # variable indicating whether insertions are present
        pos_del = 0         # initialize average relative position of read's deletions
        if read_del:
            is_del = 1                                              # tagging that deletions are present in read
            pos_del = sum(read_del) / (len(read_del) * read_len)    # compute average relative pos. of read's deletions

        de += is_del                # incrementing deletion reads counter
        del_positions += pos_del    # adding average relative position to deletion relative position counter

        # filling nucleotide pairing statistics ------------------------------------------------------------------------

        for pair in nt_pairs:       # iterating through nucleotide pairings
            pair_frequency = read_conf[pair]                # getting nt pairing's frequency from confusion table

            # counting up nucleotide pairing frequencies
            if pair_frequency in ntpair_stats[pair][0]:     # if existent: incrementing pairing frequency counter
                ntpair_stats[pair][0][pair_frequency] += 1
            else:                                           # if not yet existent: initialize frequency counter with 1
                ntpair_stats[pair][0][pair_frequency] = 1

            # (following statistics refer to observing at least one nt pairing)
            if pair_frequency > 0:

            # counting up additional insertions and relative insertion positions
                ntpair_stats[pair][2] += is_ins
                ntpair_stats[pair][4] += pos_ins

            # counting up additional deletions and relative deletion position
                ntpair_stats[pair][3] += is_ins
                ntpair_stats[pair][5] += pos_del

            # filling read length distribution
                if read_len in ntpair_stats[pair][6]: ntpair_stats[pair][6][read_len] += 1
                else: ntpair_stats[pair][6][read_len] = 1

    # ..................................................................................................................
    # computing averages and ratios (yet, absolute numbers stored)
    # ..................................................................................................................

    # counting number of reference A, T, C, G nucleotides for % nucleotide pairings
    ref_A = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["AA", "AT", "AC", "AG"]])
    ref_T = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["TA", "TT", "TC", "TG"]])
    ref_C = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["CA", "CT", "CC", "CG"]])
    ref_G = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["GA", "GT", "GC", "GG"]])

    # % reads for fragment length distribution
    for l in read_lengths:
        read_lengths[l] = read_lengths[l] / total_reads
    # average relative insertion position
    if ins: ins_positions /= ins
    # average relative deletion position
    if de: del_positions /= de

    # nucleotide pairing information
    for pair in ntpair_stats:

        # % reference nucleotides paired with corresponding read nucleotide (not considering indels)
        if pair[1] != "I" and pair[1] != "D":
            ref_N = ref_A                                                       # reference nucleotide (A)
            if pair[0] == "T": ref_N = ref_T                                    # reference nucleotide (T)
            elif pair[0] == "C": ref_N = ref_C                                  # reference nucleotide (C)
            elif pair[0] == "G": ref_N = ref_G                                  # reference nucleotide (G)
            pair_counts = sum(
                [f * ntpair_stats[pair][0][f] for f in ntpair_stats[pair][0]])  # nt pair counts
            ntpair_stats[pair][1] = pair_counts / ref_N                         # % ref nt paired with read nt

        # % reads for nt pair frequency distribution
        pair_observed = total_reads - ntpair_stats[pair][0][0]  # number of reads with at least one such nt pair
        if not pair_observed:                                   # if there is no such nt pair observed,
            ntpair_stats[pair][0][0] /= total_reads                 # only adapt 0 - occurrence information
            continue                                                # and continue
        for f in ntpair_stats[pair][0]:
            ntpair_stats[pair][0][f] /= total_reads

        # % reads with additional insertion
        ntpair_stats[pair][2] /= pair_observed
        # % reads with additional deletions
        ntpair_stats[pair][3] /= pair_observed
        # average relative insertion position
        if ntpair_stats[pair][2]: ntpair_stats[pair][4] /= ntpair_stats[pair][2]
        # average relative deletion position
        if ntpair_stats[pair][3]: ntpair_stats[pair][5] /= ntpair_stats[pair][3]
        # % reads for read length distribution
        for l in ntpair_stats[pair][6]:
            ntpair_stats[pair][6][l] /= pair_observed

    # ..................................................................................................................
    # returning
    # ..................................................................................................................

    return([total_reads, read_lengths, ins_positions, del_positions, ntpair_stats])

def conversion_counts_table(mode, input_files):
    """
    :param mode:            type of data provided for computing conversion counts; choose 'table' for providing a
                            read table, choose 'summary' for providing a read summary
    :param input_files:     list containing the file storing the read table (mode == 'table') or the antisense and
                            sense read summary files (mode == 'summary'), respectively

    :return:    a hash assigning each genomic region's name another hash, assigning each read ID within the genomic
                region a list containing the number of potential conversion positions and the number of actual
                conversions

    This function computes the number of potential conversion positions as well as the number of actual conversions of
    the genomic regions' reads stored within the read table file. The counts are stored as nested hash, assigning each
    genomic region's name another hash, assigning each read ID within the genomic region a list containing its counts.
    """

    region_convcounts = {}  # hash storing region_name -> {read_id -> [potential conversion positions, conversions]

    # computing conversion counts from read table

    if mode == "table":

        for gr_entry in rtable_file_generator(input_files[0]):  # iterating through read table file's genomic region
                                                                #  entries
            gr = list(gr_entry.keys())[0]                           # getting genomic region's name
            gr_counts = conversion_counts(gr_entry[gr][1])          # computing genomic region's read conversion counts
            region_convcounts[gr] = gr_counts                       # storing genomic region's conversion counts

    # computing conversion counts from read summary

    elif mode == "summary":

        antisense_summary = load_r_summary(input_files[0])      # loading antisense read summary
        sense_summary = load_r_summary(input_files[1])          # loading sense read summary

        for read in antisense_summary:        # initializing hash storing conversion counts with genomic region names
            region_convcounts[antisense_summary[read][0]] = {}
        for read in sense_summary:
            region_convcounts[sense_summary[read][0]] = {}

        for r in antisense_summary:           # iterating through antisense summary, counting conversions

            read_region = antisense_summary[r][0]               # getting read's region
            read_ntpairs = antisense_summary[r][4]              # getting read's nucleotide pairing table
            conv = read_ntpairs["AG"]                           # conversions
            pot = conv + read_ntpairs["AA"] + read_ntpairs["AT"] + read_ntpairs["AC"]   # potential conv. pos.
            region_convcounts[read_region][r] = [pot, conv]     # storing conversions

        for r in sense_summary:             # iterating through reverse summary, counting conversions

            read_region = sense_summary[r][0]                   # getting read's region
            read_ntpairs = sense_summary[r][4]                  # getting read's nucleotide pairing table
            conv = read_ntpairs["TC"]                           # conversions
            pot = conv + read_ntpairs["TA"] + read_ntpairs["TT"] + read_ntpairs["TG"]   # potential conv. pos.
            region_convcounts[read_region][r] = [pot, conv]     # storing conversions

    return(region_convcounts)       # returning

def conversion_efficiency_table(cctable_file, base_error, sequencing_error,
                                fix_efficiency=None, new_ini=None, coverage=50, iterations=100):
    """
    :param cctable_file:        file storing genomic regions' read conversion counts table
    :param base_error:          probability of observing a conversion due to a sequencing error
    :param sequencing_error:    probability that the converted nucleotide is masked by a sequencing error

    :param fix_efficiency:  optional; a value to fix for conversion efficiency (as a consequence, only new/total ratio
                            will be estimated) (default: None
    :param new_ini:         optional; initial value for new/total ratio (default: None)
    :param coverage:        optional; minimum coverage a genomic region has to satisfy for estimation (default: 50)
    :param iterations:      optional; number of iterations to use for efficiency estimation (default: 100)

    :return:    a hash assigning each genomic region's name a list containing two sublists (the final estimations being
                the last sublist entries, respectively):
                    the first sublist contains the sequence of estimated newly synthesized by total RNA ratios
                    the second sublist contains the sequence of estimated conversion efficiencies

    This function estimates the ratio of newly synthesized transcripts as well as the conversion efficiency for each
    genomic region, based on a given conversion counts table file. Optionally, a minimum coverage threshold can be set
    a genomic region has to fulfill to be considered for estimation (all regions not fulfilling the minimum coverage
    will be assigned a 'None' in the returned hash).
    """

    # loading in genomic regions' read conversion counts from file
    # (hash storing genomic_region -> {read_id -> [potential conversion positions, conversions]})
    cc_table = load_convcount_table(cctable_file)
    # hash storing genomic_region -> [[newly synthesized ratio estimations], [conversion efficiency estimations]]
    estimations = {}

    # iterating through genomic regions, for each region performing estimations
    for genomic_region in cc_table:
        conv_counts = cc_table[genomic_region]          # getting genomic region's read conversion counts
        if len(conv_counts) < coverage:                 # minimum coverage threshold is not fulfilled:
            estimations[genomic_region] = [None, None]      # assigning 'None' to the genomic region
            continue                                        # continue with next genomic region
        est = efficiency_estimation(conv_counts, base_error, sequencing_error, fix_efficiency=fix_efficiency,
                                    new_ini=new_ini, iterations=iterations)     # estimation
        estimations[genomic_region] = est               # storing genomic region's estimations to return hash

    # returning
    return(estimations)

def summary_table(description, libsize, cctable_file, cetable_file=None, coverage=1, refname_to_chr=r_to_c_hg19):
    """
    :param description:         a string describing the data processed in the summary matrix
    :param libsize:             library size of the data processed in the summary matrix
    :param cctable_file:        conversion counts table file

    :param cetable_file:        conversion efficiency table file (default: None)
    :param coverage:            minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number (needed to
                                associate BED chromosome names to the reference sequences) (default: default hash)

    :return:    a nested hash, the outer hash assigning each genomic region's name an inner hash assigning the
                description a list containing library size, total read counts, labeled read counts, conversion
                efficiency estimation and newly synthesized transcripts ratio estimation:
                region_name -> description -> [libsize, total, labeled, average number of potential conversion positions,
                                               conv. efficiency, newly synthesized ratio]

    ....................................................................................................................

    This function records, for each genomic region, the number of total and labeled reads as well as the estimated
    conversion efficiencies and newly synthesized transcripts ratios, given a conversion counts table file and a
    conversion efficiency table file. If no conversion efficiency table file is given, estimated conversion efficiencies
    and newly synthesized ratios are set to -1. A description of the data set processed as well as the total library
    size has to be given, too.

    Optionally, only such regions can be chosen for evaluation that fulfill a certain coverage minimum (as defined by
    parameter 'coverage'). Coverage here is defined as the number of reads associated to a genomic region.
    """

    stable = {}                                         # initializing summary table (empty hash)
    cc_table = load_convcount_table(cctable_file)       # loading in conversion counts table
    if cetable_file:
        ce_table = load_conveff_table(cetable_file)     # if given, loading in conversion efficiency table
    else:
        ce_table = None                                 # else, setting ce_table to None

    # computing read counts per region

    for gr in cc_table:             # iterating through genomic regions

        gr_cc_entry = cc_table[gr]          # getting region's cc table entry
        gr_total = len(cc_table[gr])        # getting genomic region's total read counts
        if gr_total < coverage:             # skip genomic regions not fulfilling minimum coverage
            continue

        gr_convpos = 0                      # total number of potential conversion positions in all reads of the region
        gr_labeled = 0                      # number of labeled reads
        for read_id in gr_cc_entry:         # iterating through reads
            gr_convpos += gr_cc_entry[read_id][0]   # counting up total number of potential conversion positions
            if gr_cc_entry[read_id][1] > 0:         # counting up reads with at least one conversion
                gr_labeled += 1
        gr_cp_avg = gr_convpos / gr_total   # average total number of potential conversion positions

        if ce_table:                        # if given, getting conversion efficiency and newly synthesized ratio
            gr_ce_entry = ce_table[gr]          # getting region's ce table entry
            if gr_ce_entry[1] != None:          # if successfully computed, getting region's estimations
                gr_conveff = gr_ce_entry[1][-1]     # getting region's estimated conversion efficiency
                gr_ratio = gr_ce_entry[0][-1]       # getting region's estimated newly synthesized transcripts ratio
            else:                               # else, setting estimations to -1
                gr_conveff = -1
                gr_ratio = -1
        else:                               # else, setting efficiency and ratio to -1
            gr_conveff = -1                     # setting region's estimated conversion efficiency to -1
            gr_ratio = -1                       # setting region's estimated newly synthesized transcripts ratio to -1

        # storing genomic region's read counts to summary table hash
        stable[gr] = {description: [libsize, gr_total, gr_labeled, gr_cp_avg, gr_conveff, gr_ratio]}

    # returning

    return(stable)

# Wrappers #############################################################################################################

def new_read_table(mode, bamfile, reffile, bedfile, outfile, snpfile=None, coverage=1, refname_to_chr=r_to_c_hg19):
    """
    :param mode:        mode to be executed; choose 'snp_all' to compute SNP tables reporting all kind of conversions,
                        choose 'snp_exclusive' to compute SNP tables reporting all but labeling-specific conversions,
                        choose 'conv' to compute conversion tables only, use 'snp_all_summary', 'snp_exclusive_summary',
                        'conv_summary' to compute summarized tables that report SNP and total read counts for each
                        position, choose 'both' or 'both_summary' to compute both a SNP (snp_exclusive) and a conversion
                        table
    :param bamfile:     input BAM file containing aligned reads
    :param reffile:     fasta file containing all reference sequences; an index of the same file needs to be located in
                        the same directory, having the same filename with the additional extension '.fai'
    :param bedfile:     BED file containing genomic regions for which reads are to be filtered
    :param outfile:     output file to store read table to

    :param snpfile:         gzipped vcf file containing SNP positions of the reference (default: None)
    :param coverage:        minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate the SNP sites / A-to-I sites to the reference sequences) (default: default hash)

    :return: void

    ....................................................................................................................

    This function computes a read table for a BAM file and writes it to an output file. The read table records SNPs as
    well as conversions for each genomic region, resolved by reads. Optionally, conversion positions can be corrected
    for SNP positions by providing a SNP vcf file. In this case, all positions listed to be prone to SNPs will be
    reported as None in the table. Optionally, a minimum coverage can be specified by which genomic regions are
    filtered; only regions hitting this minimum coverage will be processed and included into the read table.
    """

    r_table = read_table(mode, bamfile, reffile, bedfile, snpfile=snpfile,
                         coverage=coverage, refname_to_chr=refname_to_chr)
    write_read_table(r_table, outfile)
    return()

def snp_corrected_rsummary(bamfile, reffile, snpfile, outfile, aifile=None, editfile=None, refname_to_chr=r_to_c_hg19):
    """
    :param bamfile:         input BAM file containing aligned reads
    :param reffile:         fasta file containing all reference sequences; an index of the same file needs to be located
                            in the same directory, having the same filename with the additional extension '.fai'
    :param snpfile:         gzipped vcf file containing SNP positions of the reference
    :param outfile:         output file to store read summary to

    :param aifile:          REDIportal file containing A-to-I editing positions of the reference (default: None)
    :param editfile:        editing sites file storing editing sites as 'ChrNumber_ChrPos' (default: None)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate the SNP sites / A-to-I sites to the reference sequences) (default: default hash)

    :return: void

    ....................................................................................................................

    This function collects read information for all mapped reads within a given BAM file and writes it to an output
    file.
    The information retrieved are: reference sequence name mapped to, read length, positions of insertions (position
    before insertion event, w.r.t. reference, 1-based), positions of deletions (position before deletion event, w.r.t.
    reference, 1-based), confusion table listing the number of nucleotide pairings as well as insertions and deletions
    within the read.
    Confusion table counts (nucleotide pairing counts) are corrected for SNP sites and, optionally, for editing sites.
    """

    r_summary = read_summary(bamfile, reffile, snpfile=snpfile, aifile=aifile, editfile=editfile,
                             refname_to_chr=refname_to_chr)
    write_r_summary(r_summary, outfile)
    return()

def antisense_sense_rsummary(bamfile, rsummary_infile, bedfile, antisense_outfile, sense_outfile,
                             coverage=1, refname_to_chr=r_to_c_hg19):
    """
    :param bamfile:             input BAM file containing aligned reads
    :param rsummary_infile:     file storing the BAM file's read summary
    :param bedfile:             BED file containing genomic regions for which reads are to be filtered and separated
                                into antisense and sense genomic region reads
    :param antisense_outfile:   output file to store antisense genomic regions' reads summary to
    :param sense_outfile:       output file to store sense genomic regions' reads summary to

    :param coverage:            minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number (needed to
                                associate BED chromosome names to the reference sequences) (default: default hash)

    :return: void

    ....................................................................................................................

    This function filters all reads from a BAM file that map to user-defined genomic regions and separates them into
    reads originating from antisense genomic regions and sense genomic regions. Based on the BAM file and its original
    read summary as well as a BED file specifying genomic regions, a read summary for both the antisense and sense
    genomic regions' reads is written.
    """

    genomic_regions = filter_regions(bamfile, bedfile, refname_to_chr=refname_to_chr)
    genomic_regions = filter_coverage(genomic_regions, coverage=coverage)
    original_rsummary = load_r_summary(rsummary_infile)
    antisense_sense_separation = antisense_sense_transcripts(original_rsummary, genomic_regions)
    write_r_summary(antisense_sense_separation[0], antisense_outfile)
    write_r_summary(antisense_sense_separation[1], sense_outfile)
    return()

def rsummary_stats(rsummary_infile, rstats_outfile):
    """
    :param rsummary_infile: read summary input file
    :param rstats_outfile:  read summary statistics output file

    :return: void

    ....................................................................................................................

    This function computes statistics to a given read summary and writes the statistics to an output file.
    """

    rsummary = load_r_summary(rsummary_infile)
    rstats = read_summary_statistics(rsummary)
    write_r_stats(rstats, rstats_outfile)
    return()

def new_convcount_table(mode, input_files, outfile):
    """
    :param mode:            type of data provided for computing conversion counts; choose 'table' for providing a
                            read table, choose 'summary' for providing a read summary
    :param input_files:     file storing the read table (mode == 'table') or the antisense and sense read summary files
                            (mode == 'summary'), respectively
    :param outfile:         output file to store conversion counts table to

    :return: void

    ....................................................................................................................

    This function computes a conversion count table from a given read table file and writes it to an output file.
    """
    cc_table = conversion_counts_table(mode, input_files)
    write_convcount_table(cc_table, outfile)
    return()

def merged_convcount_table(cctable_files, outfile):
    """
    :param cctable_files:   list of conversion count table files
    :param outfile:         output file to store merged conversion count table file to

    :return: void

    ....................................................................................................................

    This function takes multiple conversion count table files and merges them into one. Only genomic regions present in
    all conversion count tables will be considered.
    """

    all_cctables = [load_convcount_table(file) for file in cctable_files]   # loading in all conversion count tables
    merged_cctable = {}                                                     # initializing merged cc table

    for gr in all_cctables[0]:      # iterating through genomic regions

        shared = True                   # variable indicating whether genomic region is present in all cc tables
        for cct in all_cctables[1:]:    # iterating through all conversion counts tables
            if gr not in cct:               # if region is not present in one of them, setting 'shared' to False
                shared = False
                break

        if shared:                      # if region is present in all cc tables, storing region's single entries
            gr_entry = {}                   # initializing regions entry
            for cct in all_cctables:            # iterating through cc tables
                for read_id in cct[gr]:             # iterating through cc table's read entries
                    gr_entry[read_id] = cct[gr][read_id]    # adding read entries to region entry
            merged_cctable[gr] = gr_entry   # storing region to merged conversion counts table

    write_convcount_table(merged_cctable, outfile)      # writing merged conversion counts table
    return()                                            # returning

def new_conveff_table(cctable_file, outfile, base_error, sequencing_error,
                      fix_efficiency=None, new_ini=None, coverage=50, iterations=100):
    """
    :param cctable_file:        input conversion counts table file
    :param outfile:             output file to store conversion efficiency table to
    :param base_error:          probability of observing a conversion due to a sequencing error
    :param sequencing_error:    probability that the converted nucleotide is masked by a sequencing error

    :param fix_efficiency:      optional; a value to fix for conversion efficiency (as a consequence, only new/total
                                ratio will be estimated) (default: None
    :param new_ini:             optional; initial value for new/total ratio (default: None)
    :param coverage:            minimum coverage a genomic region has to satisfy for estimation (default: 50)
    :param iterations:          number of iterations to use for efficiency estimation (default: 100)

    :return: void

    ....................................................................................................................

    This function computes a conversion efficiency table from a given conversion counts table file and writes it to an
    output file. In detail, the ratio of newly synthesized transcripts and the conversion efficiency are estimated for
    each genomic region using an EM algorithm. Optionally, a minimum coverage threshold can be set a genomic region has
    to fulfill to be considered for estimation.
    """

    conveff_table = conversion_efficiency_table(cctable_file, base_error, sequencing_error,
                                                fix_efficiency=fix_efficiency, new_ini=new_ini,
                                                coverage=coverage, iterations=iterations)
    write_conveff_table(conveff_table, outfile)
    return()

def new_summary_table(outfile, description, libsize, cctable_file, cetable_file=None, coverage=1,
                      refname_to_chr=r_to_c_hg19):
    """
    :param outfile:             output file to write summary table to
    :param description:         a string describing the data processed in the summary matrix
    :param libsize:             library size of the data processed in the summary matrix
    :param cctable_file:        conversion counts table file

    :param cetable_file:        conversion efficiency table file (default: None)
    :param coverage:            minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number (needed to
                                associate BED chromosome names to the reference sequences) (default: default hash)

    :return:    a hash assigning each genomic region's name a list containing (in the following order) the description,
                unlabeled read counts, labeled read counts

    ....................................................................................................................

    This function computed a summary table for a given a conversion counts table and conversion efficiency table file
    and writes the table to an output file. If no conversion efficiency table file is given, estimated conversion
    efficiencies and newly synthesized ratios are set to -1. Genomic regions to be included in the matrix can be
    selected by defining a coverage minimum (as defined by parameter 'coverage'). Coverage here is defined as the
    number of reads associated to a genomic region.

    The output file is of tab-delimited format:
    region_name\tdescription\tlibrary_size\ttotal_reads\tlabeled_reads\taverage_potential_conversion_positions\t
    conv_efficiency\tnewly_synthesized_ratio\n
    """

    stable = summary_table(description, libsize, cctable_file, cetable_file=cetable_file,
                           coverage=coverage, refname_to_chr=refname_to_chr)    # computing summary table
    write_summary_table(stable, outfile)                                        # writing summary table
    return()                                                                    # returning

def merged_summary_table(mode, stable_files, outfile):
    """
    :param mode:            mode to be executed: choose 'integrate' for integrating tables' entries, choose 'merge'
                            for merging tables into one matrix
    :param stable_files:    list of summary table files
    :param outfile:         output file to store integrated / merged summary table to

    :return: void

    ....................................................................................................................

    This function takes multiple summary tables and either integrates or merges them into one.

    Integration will integrate all entries that are of equal genomic region name and data description. In detail,
    library size as well as read counts will be summed up, average conversion positions as well as conversion effi-
    ciencies and newly synthesized ratios will be weightedly averaged. Only genomic region name and data description
    combinations present in all input summary tables will be considered.
    Merging will collect genomic regions' entries. Only genomic regions that are present in all input summary tables
    will be collected.

    The integrated / merged summary table will be sorted by genomic regions.
    """

    # initializing

    out_file = open(outfile, "w")                                           # opening output file
    all_stables = [load_summary_table(file) for file in stable_files]       # loading all summary tables
    first_table = all_stables[0]                                            # getting first table
    shared_regions = []                 # list storing genomic regions that occur in all summary tables

    # filtering genomic regions (and data descriptions) that occur in all summary tables

    if mode == "integrate":     # CASE: integrate matrices

        for gr in first_table:              # iterating through regions (of first summary table)

            shared = True                   # storing if region + [descriptions] occurs in all tables (set to True)
            for des in first_table[gr]:     # iterating through region's data descriptions
                for tab in all_stables[1:]:     # iterating through remaining matrices
                    if gr not in tab:                       # checking if region is missing in matrix
                        shared = False                      # 'shared' set to False (doesn't occur in all matrices)
                    elif des not in tab[gr]:                # checking if description is missing in matrix
                        shared = False                      # 'shared' set to False (doesn't occur in all matrices)
            if shared:                      # if region + [descriptions] occurs in all matrices,
                for des in first_table[gr]:    # store region and descriptions to selected regions
                    shared_regions.append([gr, des])

    elif mode == "merge":       # CASE: merge matrices

        for gr in first_table:              # iterating through regions (of first summary matrix)

            shared = True                           # storing if region occurs in all matrices (set to True)
            for tab in all_stables[1:]:             # iterating through remaining matrices
                if gr not in tab:                       # checking if region is missing in matrix
                    shared = False                      # 'shared' set to False (doesn't occur in all matrices)

            if shared:                              # if region occurs in all matrices:
                shared_regions.append(gr)           # store region to selected regions

    # integrating summary tables sorted by region name

    if mode == "integrate":

        for gr_des in shared_regions:       # iterating through shared genomic regions + descriptions
            gr = gr_des[0]                      # genomic region name
            des = gr_des[1]                     # data description
            libsize = 0                         # summed library size (initialized with 0)
            total = 0                           # summed total read counts (initialized with 0)
            labeled = 0                         # summed labeled read counts (initialized with 0)
            convpos = 0                         # weighted summed potential conversion positions
            conveff = 0                         # weighted summed conversion efficiencies (initialized with 0)
            ratios = 0                          # weighted summed newly synthesized ratios (initialized with 0)

            for tab in all_stables:             # iterating through summary tables
                tab_total = tab[gr][des][1]             # table entry's total number of reads
                libsize += tab[gr][des][0]              # summing up library size
                total += tab_total                      # summing up total reads
                labeled += tab[gr][des][2]              # summing up labeled reads
                convpos += tab_total * tab[gr][des][3]  # summing up weighted potential conversion positions
                conveff += tab_total * tab[gr][des][4]  # summing up weighted conversion efficiencies
                ratios += tab_total * tab[gr][des][5]   # summing up weighted newly synthesized ratios

            convpos /= total                    # averaging potential conversion positions
            conveff /= total                    # averaging conversion efficiencies
            ratios /= total                     # averaging newly synthesized ratios

            out_file.write(
                gr + "\t" + des + "\t" +
                "\t".join([str(i) for i in [libsize, total, labeled, convpos,
                                            conveff, ratios]]) + "\n")              # writing entry

    # merging summary matrices sorted by region name

    elif mode == "merge":

        for gr in shared_regions:               # iterating through shared genomic regions
            for tab in all_stables:                 # iterating through summary matrices
                for d in tab[gr]:                       # iterating through genomic region's entries
                    entry = tab[gr][d]                                                      # getting entry
                    out_file.write(
                        gr + "\t" + d + "\t" + "\t".join([str(e) for e in entry]) + "\n")   # writing entry

    # returning

    return()

# Pipelines ############################################################################################################

def preprocess_reads(bamfile, reffile, snpfile, bedfile, rsum_outfile, rsum_antisense_outfile, rsum_sense_outfile,
                     bam_recorded_outfile, bam_spare_outfile, cctable_outfile, editfile=None):

    """
    :param bamfile:     input BAM file containing aligned reads
    :param reffile:     fasta file containing all reference sequences; an index of the same file needs to be located
                        in the same directory, having the same filename with the additional extension '.fai'
    :param snpfile:     gzipped vcf file containing SNP positions of the reference
    :param bedfile:     BED file containing genomic regions for which reads are to be filtered and separated into
                        antisense and sense genomic region reads

    :param rsum_outfile:                output file to store read summary to
    :param rsum_antisense_outfile:      output file to store antisense genomic regions' reads summary to
    :param rsum_sense_outfile:          output file to store sense genomic regions' reads summary to
    :param bam_recorded_outfile:        output BAM file to store aligned reads matched to genomic regions defined
                                        within the BED file
    :param bam_spare_outfile:           output BAM file to store aligned reads not matched to genomic regions defined
                                        within the BED file
    :param cctable_outfile:             output file to store conversion counts table to

    :param editfile:    editing sites file storing editing sites as 'ChrNumber_ChrPos' (default: None)

    :return:    void

    ....................................................................................................................

    This function computes a read summary, filters and separates the read summary according to antisense and sense
    genomic regions as defined in a given BED file, computes read summary statistics, splits the input BAM file into
    reads matched and reads not matched to the BED file regions, and generates a conversion counts table based on the
    reads matched to the BED file regions.
    """

    # computing read summary
    snp_corrected_rsummary(bamfile, reffile, snpfile, rsum_outfile, editfile=editfile)
    gc.collect()
    # filter and split read summary according to BED file regions
    antisense_sense_rsummary(bamfile, rsum_outfile, bedfile, rsum_antisense_outfile, rsum_sense_outfile)
    gc.collect()
    # compute read summary statistics
    rsummary_stats(rsum_outfile, rsum_outfile + "_stats.txt")
    rsummary_stats(rsum_antisense_outfile, rsum_antisense_outfile + "_stats.txt")
    rsummary_stats(rsum_sense_outfile, rsum_sense_outfile + "_stats.txt")
    gc.collect()

    # split BAM file into reads matching and not matching to BED file regions
    filter_spare_reads(bamfile, [rsum_antisense_outfile, rsum_sense_outfile], bam_recorded_outfile, bam_spare_outfile)
    gc.collect()

    # compute conversion counts table for the reads matching to BED file regions
    new_convcount_table("summary", [rsum_antisense_outfile, rsum_sense_outfile], cctable_outfile)

    return()

# Main Method (Executes Wrappers and Pipelines) ########################################################################

if __name__ == "__main__":

    ####################################################################################################################
    # getting docstrings for argument help strings
    ####################################################################################################################

    # bed_annotation (bed_annotations) help strings --------------------------------------------------------------------

    ban_help = bed_annotations.__doc__      # getting docstring
    ban_help = ban_help.split(":param")     # splitting docstring at 'param' tags

    ban_description = ban_help[-1].split("....")[-1].strip()    # bed_annotation description

    ban_infile_help = ban_help[1].split(":")[1].strip()     # input file help
    ban_gff_help = ban_help[2].split(":")[1].strip()        # gff3 file help
    ban_exon_help = ban_help[3].split(":")[1].strip()       # exon file help
    ban_intron_help = ban_help[4].split(":")[1].strip()     # intron file help
    ban_3utr_help = ban_help[5].split(":")[1].strip()       # 3'UTRs file help
    ban_5utr_help = ban_help[6].split(":")[1].strip()       # 5'UTRs file help
    ban_outfile_help = ban_help[7].split(":")[1].strip()    # output file help

    ban_bamfiles_help = ban_help[8].split(":")[1].strip()   # BAM files help
    ban_bamdes_help = ban_help[9].split(":")[1].strip()     # BAM descriptors help

    # bed_annotation_statistics (annotation_statistics) help strings ---------------------------------------------------

    bas_description = \
        "This function calculates basic statistics on an annotation file"   # bed_annotation_statistcs description
    bas_infile_help = "input file storing annotations"                      # input file help
    bas_outfile_help = "output file to write annotation statistics to"      # output file help

    # denovo_editsites (write_edit_sites) help strings -----------------------------------------------------------------

    des_description = "This function detects editing sites, computing position-wise SNP rates from a read table file " \
                      "and storing all positions with a SNP rate above a certain threshold. Editing positions are " \
                      "written to an output file as 'ChrNumber_ChrPos'"     # denovo_editsites description
    des_rtable_help = "read table file"                                     # read table help
    des_cutoff_help = "maximum editing rate (SNP rate) of a nucleotide " \
                      "position to be tolerated"                            # editing cutoff help
    des_outfile_help = "output file to store editing sites to"              # output file help

    # filter_spare_reads (filter_spare_reads) help strings -------------------------------------------------------------

    fsr_description = "This function separates all aligned reads of a given BAM input file by reads which are " \
                      "recorded in any of the read summaries specified and those which are not. The two groups of " \
                      "reads are written to separate output BAM files."     # filter_spare_reads help
    fsr_bamin_help = "input BAM file containing aligned reads"              # input BAM file help
    fsr_bamoutrec_help = "output BAM file to store aligned reads " \
                         "recorded within the read summaries"               # output BAM file (recorded reads) help
    fsr_bamoutspare_help = "output BAM file to store aligned reads not " \
                           "recorded within the read summaries"             # output BAM file (spare reads) help
    fsr_rsum_help = "one or more read summary files"                        # read summary files help

    # covered_positions (covered_positions) help strings ---------------------------------------------------------------

    cpo_description = "This function evaluates covered positions for the reads stored in a BAM file. To do so, the " \
                      "position of the 3'-most nucleotide of a read is considered only (since these are the " \
                      "suggested sites for the 3'UTR end / pA start)."      # covered_positions help
    cpo_bamin_help = "input BAM file containing aligned reads"              # input BAM file help
    cpo_outfile_help = "output file to store covered positions to"          # output file help

    # merge_covered_positions (covered_positions_merge) help strings ---------------------------------------------------

    mcp_description = "This function merges covered positions as stored in the input files given. In detail, the " \
                      "coverages are summed up over all input files for each position that is recorded in any of the " \
                      "input files."                                                # covered_positions_merge help
    mcp_infiles_help = "files storing covered positions that are to be merged"      # input files help
    mcp_outfile_help = "output file to store merged covered positions to"           # output file help

    # coverage_profile (coverage_profile) help strings -----------------------------------------------------------------

    cpr_description = coverage_profile.__doc__.split("....")[-1].strip()    # coverage_profile help
    cpr_infile_help = "file storing covered positions"                      # input file help
    cpr_outfile_help = "output file to store coverage profile"              # output file help
    cpr_radius_help = "number of positions both upstream and downstream to be evaluated for the coverage profile " \
                      "(default: 1000)"                                     # radius help

    # coverage_convolution (coverage_convolution) help strings ---------------------------------------------------------

    ccv_description = coverage_convolution.__doc__.split("....")[-1].strip()            # coverage convolution help
    ccv_infile_help = "input file storing covered positions"                            # input file help
    ccv_outfile_help = "output file to store coverage convolution"                      # output file help
    ccv_radius_help = "nucleotide radius used for convolution (default: 50)"            # radius help
    ccv_eval_help = "length of a nucleotide stretch to be evaluated at a time " \
                    "(this is to avoid loading in all data at once) (default: 10**6)"   # evaluation length help

    # covered_regions (covered_regions) help strings -------------------------------------------------------------------

    cvr_description = "This function defines covered regions and writes them to an output BED file. A covered region " \
                      "is defined based on a coverage peak and a nucleotide radius of " \
                      "fixed size around that peak."                        # covered regions help
    cvr_infile_help = "input file storing coverage convolution"             # input file help
    cvr_outfile_help = "output BED file to store covered regions"           # output file help
    cvr_radius_help = "radius around a coverage peak position to be defined as " \
                      "covered region (default: 50)"                        # radius help

    # denovo_readtable (new_read_table) help strings -------------------------------------------------------------------

    drt_help = new_read_table.__doc__           # getting docstring
    drt_help = drt_help.split(":param")         # splitting docstring at 'param' tags

    drt_description = drt_help[-1].split("....")[-1].strip()    # denovo_readtable description
    drt_mode_help = drt_help[1].split(":")[1].strip()           # mode help
    drt_bam_help = drt_help[2].split(":")[1].strip()            # BAM file help
    drt_ref_help = drt_help[3].split(":")[1].strip()            # reference file help
    drt_bed_help = drt_help[4].split(":")[1].strip()            # BED file help
    drt_out_help = drt_help[5].split(":")[1].strip()            # output file help
    drt_snp_help = drt_help[6].split(":")[1].strip()            # SNP file help
    drt_cov_help = drt_help[7].split(":")[1].strip()            # coverage help

    # denovo_rsummary (snp_corrected_rsummary) help strings ------------------------------------------------------------

    drs_help = snp_corrected_rsummary.__doc__       # getting docstring
    drs_help = drs_help.split(":param")             # splitting docstring at 'param' tags

    drs_description = drs_help[-1].split("....")[-1].strip()    # denovo_rsummary description
    drs_bam_help = drs_help[1].split(":")[1].strip()            # BAM file help
    drs_ref_help = drs_help[2].split(":")[1].strip()            # reference file help
    drs_snpfile_help = drs_help[3].split(":")[1].strip()        # SNP file help
    drs_outfile_help = drs_help[4].split(":")[1].strip()        # output file help
    drs_editfile_help = drs_help[6].split(":")[1].strip()       # edit sites file help

    # filter_rsummary (antisense_sense_rsummary) help strings ----------------------------------------------------------

    frs_help = antisense_sense_rsummary.__doc__     # getting docstring
    frs_help = frs_help.split(":param")             # splitting docstring at 'param' tags

    frs_description = frs_help[-1].split("....")[-1].strip()    # filter_rsummary description
    frs_bam_help = frs_help[1].split(":")[1].strip()            # BAM file help
    frs_rsumfile_help = frs_help[2].split(":")[1].strip()       # read summary file help
    frs_bed_help = frs_help[3].split(":")[1].strip()            # BED file help
    frs_aout_help = frs_help[4].split(":")[1].strip()           # antisense read summary output file help
    frs_sout_help = frs_help[5].split(":")[1].strip()           # sense read summary output file help
    frs_cov_help = frs_help[6].split(":")[1].strip()            # coverage help

    # rsummary_stats (rsummary_stats) help strings ---------------------------------------------------------------------

    rss_help = rsummary_stats.__doc__       # getting docstring
    rss_help = rss_help.split(":param")     # splitting docstring at 'param' tags

    rss_description = rss_help[-1].split("....")[-1].strip()    # rsummary_stats description
    rss_rsumfile_help = rss_help[1].split(":")[1].strip()       # read summary file help
    rss_outfile_help = rss_help[2].split(":")[1].strip()        # summary statistics output file

    # denovo_cctable (new_convcount_table) help strings ----------------------------------------------------------------

    dcc_help = new_convcount_table.__doc__      # getting docstring
    dcc_help = dcc_help.split(":param")         # splitting docstring at 'param' tags

    dcc_description = dcc_help[-1].split("....")[-1].strip()    # denovo_cctable description
    dcc_mode_help = dcc_help[1].split(":")[1].strip()           # mode help
    dcc_infiles_help = dcc_help[2].split(":")[1].strip()        # read information files help
    dcc_out_help = dcc_help[3].split(":")[1].strip()            # output file help

    # merge_cctable (merged_convcount_table) help strings --------------------------------------------------------------

    mcc_help = merged_convcount_table.__doc__   # getting docstring
    mcc_help = mcc_help.split(":param")         # splitting docstring at 'param' tags

    mcc_description = mcc_help[-1].split("....")[-1].strip()    # merged_cctable description
    mcc_cctables_help = mcc_help[1].split(":")[1].strip()       # cc table files help
    mcc_out_help = mcc_help[2].split(":")[1].strip()            # output file help

    # denovo_convefftable (new_conveff_table) help strings -------------------------------------------------------------

    dce_help = new_conveff_table.__doc__    # getting docstring
    dce_help = dce_help.split(":param")     # splitting docstring at 'param' tags

    dce_description = dce_help[-1].split("....")[-1].strip()        # denovo_confefftable description
    dce_cctable_help = dce_help[1].split(":")[1].strip()            # cc table file help
    dce_outfile_help = dce_help[2].split(":")[1].strip()            # output file help
    dce_base_help = dce_help[3].split(":")[1].strip()               # base error help
    dce_seq_help = dce_help[4].split(":")[1].strip()                # sequencing error help
    dce_fixce_help = dce_help[5].split(":")[1].strip()              # fixed conversion efficiency help
    dce_newini_help = dce_help[6].split(":")[1].strip()             # new initial value help
    dce_cov_help = dce_help[7].split(":")[1].strip()                # coverage help
    dce_iter_help = dce_help[8].split(":")[1].strip()               # iterations help

    # denovo_summarytable (new_summary_table) help strings -------------------------------------------------------------

    dst_help = new_summary_table.__doc__        # getting docstring
    dst_help = dst_help.split(":param")         # splitting docstring at 'param' tags

    dst_description = dst_help[-1].split("....")[-1].strip()    # denovo_summarytable description
    dst_out_help = dst_help[1].split(":")[1].strip()            # output file help
    dst_descr_help = dst_help[2].split(":")[1].strip()          # data description help
    dst_libsize_help = dst_help[3].split(":")[1].strip()        # library size help
    dst_cctable_help = dst_help[4].split(":")[1].strip()        # cc table file help
    dst_cetable_help = dst_help[5].split(":")[1].strip()        # ce table file help
    dst_cov_help = dst_help[6].split(":")[1].strip()            # coverage help

    # merge_summarytable (merged_summary_table) help strings -----------------------------------------------------------

    mst_help = merged_summary_table.__doc__         # getting docstring
    mst_help = mst_help.split(":param")             # splitting docstring at 'param' tags

    mst_description = mst_help[-1].split("....")[-1].strip()    # merge_summarytable description
    mst_mode_help = mst_help[1].split(":")[1].strip()           # mode help
    mst_out_help = mst_help[3].split(":")[1].strip()            # output file help
    mst_in_help = mst_help[2].split(":")[1].strip()             # input summary table files help

    # preprocess_reads (preprocess_reads) help strings -----------------------------------------------------------------

    ppr_help = preprocess_reads.__doc__     # getting docstring
    ppr_help = ppr_help.split(":param")     # splitting dicstring at 'param' tags

    ppr_description = ppr_help[-1].split("....")[-1].strip()    # preprocess_reads description
    ppr_bamfile_help = ppr_help[1].split(":")[1].strip()        # BAM file help
    ppr_reffile_help = ppr_help[2].split(":")[1].strip()        # reference file help
    ppr_snpfile_help = ppr_help[3].split(":")[1].strip()        # SNP file help
    ppr_bedfile_help = ppr_help[4].split(":")[1].strip()        # BED file help
    ppr_rsumfile_help = ppr_help[5].split(":")[1].strip()       # read summary file help
    ppr_antifile_help = ppr_help[6].split(":")[1].strip()       # antisense read summary file help
    ppr_sensefile_help = ppr_help[7].split(":")[1].strip()      # sense read summary file help
    ppr_recfile_help = ppr_help[8].split(":")[1].strip()        # recorded reads BAM file help
    ppr_sparefile_help = ppr_help[9].split(":")[1].strip()      # spare reads BAM file help
    ppr_ccfile_help = ppr_help[10].split(":")[1].strip()        # conversion counts file help
    ppr_editfile_help = ppr_help[11].split(":")[1].strip()      # editing sites file help

    ####################################################################################################################
    # defining arguments
    ####################################################################################################################

    # general parser and subparsers ------------------------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description="This module can run different data analyses:\n\n"

                    "Choose 'bed_annotation' to annotate 3'UTRs.\n"
                    "Choose 'bed_annotation_statistics' to compute basic statistics on an annotation file.\n"
                    "Choose 'denovo_editsites' to detect editing sites.\n"
                    "Choose 'filter_spare_reads' to filter BAM file reads for those being recorded in a read summary "
                    "and those which are not.\n"
                    "Choose 'covered_positions' to compute covered positions.\n"
                    "Choose 'merge_covered_positions' to merge covered positions.\n"
                    "Choose 'coverage_profile' to compute a coverage profile.\n"
                    "Choose 'coverage_convolution' to perform a coverage convolution.\n\n"

                    "Choose 'denovo_readtable' to compute a new read table.\n"
                    "Choose 'denovo_rsummary' to compute a new read summary.\n"
                    "Choose 'filter_rsummary' to filter and split an existing read summary for antisense and sense "
                    "genomic regions.\n"
                    "Choose 'rsummary_stats' to compute read summary statistics.\n"
                    "Choose 'denovo_cctable' to compute a new conversion counts table.\n"
                    "Choose 'merge_cctable' to merge multiple conversion counts tables into one table.\n"
                    "Choose 'denovo_convefftable' to compute a new conversion efficiency table.\n"
                    "Choose 'denovo_summarytable' to compute a new summary table.\n"
                    "Choose 'merge_summarytable' to integrate or merge multiple summary tables into one summary "
                    "table.\n\n"
        
                    "Choose 'preprocess_reads' to compute a read summary, filter the read summary for genomic regions, "
                    "compute read summary statistics, filter the BAM file for reads matched and not matched to the "
                    "genomic regions, and compute a conversion counts table based on the reads matched to the regions."
                    "\n\n"

                    "For a more detailed description of the single data analyses, please call the analyses' help "
                    "messages by typing 'confusion_table <analysis_name> -h'.")     # main parser

    subparsers = parser.add_subparsers(dest="analysis", help="data analysis to be performed")       # subparsers

    parser_ban = subparsers.add_parser("bed_annotation", description=ban_description,
                                       help="annotate 3'UTRs")                                          # ban subparser
    parser_bas = subparsers.add_parser("bed_annotation_statistics", description=bas_description,
                                       help="compute annotation statistics")                            # bas subparser
    parser_des = subparsers.add_parser("denovo_editsites", description=des_description,
                                       help="computes a new editing sites file")                        # des subparser
    parser_fsr = subparsers.add_parser("filter_spare_reads", description=fsr_description,
                                       help="filter BAM file reads by whether or not they are recorded within a "
                                            "read summary")                                             # fsr subparser
    parser_cpo = subparsers.add_parser("covered_positions", description=cpo_description,
                                       help="compute covered positions")                                # cpo subparser
    parser_mcp = subparsers.add_parser("merge_covered_positions", description=mcp_description,
                                       help="merge covered positions")                                  # mcp subparser
    parser_cpr = subparsers.add_parser("coverage_profile", description=cpr_description,
                                       help="compute coverage profile")                                 # cpr subparser
    parser_ccv = subparsers.add_parser("coverage_convolution", description=ccv_description,
                                       help="execute coverage convolution")                             # ccv subparser
    parser_cvr = subparsers.add_parser("covered_regions", description=cvr_description,
                                       help="find covered regions")                                     # cvr subparser

    parser_drt = subparsers.add_parser("denovo_readtable", description=drt_description,
                                       help="computes a new read table")                                # drt subparser
    parser_drs = subparsers.add_parser("denovo_rsummary", description=drs_description,
                                       help="compute a new read summary")                               # drs subparser
    parser_frs = subparsers.add_parser("filter_rsummary", description=frs_description,
                                       help="filter an existing read summary for genomic regions")      # frs subparser
    parser_rss = subparsers.add_parser("rsummary_stats", description=rss_description,
                                       help="compute statistics for a read summary")                    # rss subparser
    parser_dcc = subparsers.add_parser("denovo_cctable", description=dcc_description,
                                       help="compute a new conversion counts table")                    # dcc subparser
    parser_mcc = subparsers.add_parser("merge_cctable", description=mcc_description,
                                       help="merge multiple conversion count tables")                   # mcc subparser
    parser_dce = subparsers.add_parser("denovo_convefftable", description=dce_description,
                                       help="compute a new conversion efficiency table")                # dce subparser
    parser_dst = subparsers.add_parser("denovo_summarytable", description=dst_description,
                                       help="compute a new summary table")                              # dst subparser
    parser_mst = subparsers.add_parser("merge_summarytable", description=mst_description,
                                       help="add or merge multiple summary tables")                     # mst subparser

    parser_ppr = subparsers.add_parser("preprocess_reads", description=ppr_description,
                                       help="pre-process reads, applying multiple analysis steps")      # ppr subparser

    # bed_annotation (bed_annotations) arguments -----------------------------------------------------------------------

    parser_ban.add_argument("bedfile", type=str, help=ban_infile_help)                      # input BED file
    parser_ban.add_argument("gff3file", type=str, help=ban_gff_help)                        # GFF3 file
    parser_ban.add_argument("exonfile", type=str, help=ban_exon_help)                       # exon BED file
    parser_ban.add_argument("intronfile", type=str, help=ban_intron_help)                   # intron BED file
    parser_ban.add_argument("utr3file", type=str, help=ban_3utr_help)                       # 3'UTR BED file
    parser_ban.add_argument("utr5file", type=str, help=ban_5utr_help)                       # 5'UTR BED file
    parser_ban.add_argument("outfile", type=str, help=ban_outfile_help)                     # output file
    parser_ban.add_argument("--bamfiles", type=str, nargs="*", help=ban_bamfiles_help)      # BAM files
    parser_ban.add_argument("--bamdes", type=str, nargs="*", help=ban_bamdes_help)          # BAM descriptions

    # bed_annotation_statistics (annotation_statistics) arguments ------------------------------------------------------

    parser_bas.add_argument("infile", type=str)     # input file
    parser_bas.add_argument("outfile", type=str)    # output file

    # denovo_editsites (write_edit_sites) arguments --------------------------------------------------------------------

    parser_des.add_argument("rtablefile", type=str, help=des_rtable_help)       # read table file
    parser_des.add_argument("cutoff", type=float, help=des_cutoff_help)         # editing cutoff
    parser_des.add_argument("outfile", type=str, help=des_outfile_help)         # output file

    # filter_spare_reads (filter_spare_reads) arguments ----------------------------------------------------------------

    parser_fsr.add_argument("bam_in", type=str, help=fsr_bamin_help)                # input BAM file
    parser_fsr.add_argument("bam_out_rec", type=str, help=fsr_bamoutrec_help)       # output (recorded reads) BAM file
    parser_fsr.add_argument("bam_out_spare", type=str, help=fsr_bamoutspare_help)   # output (spare reads) BAM file
    parser_fsr.add_argument("rsum_files", type=str, nargs="*", help=fsr_rsum_help)  # read summary files

    # covered_positions (covered_positions) arguments ------------------------------------------------------------------

    parser_cpo.add_argument("bam_in", type=str, help=cpo_bamin_help)        # input BAM file
    parser_cpo.add_argument("outfile", type=str, help=cpo_outfile_help)     # output file

    # merge_covered_positions (covered_positions_merge) arguments ------------------------------------------------------

    parser_mcp.add_argument("outfile", type=str, help=mcp_outfile_help)                 # output file
    parser_mcp.add_argument("input_files", type=str, nargs="*", help=mcp_infiles_help)  # input files

    # coverage_profile (coverage_profile) arguments --------------------------------------------------------------------

    parser_cpr.add_argument("infile", type=str, help=cpr_infile_help)                   # input file
    parser_cpr.add_argument("outfile", type=str, help=cpr_outfile_help)                 # output file
    parser_cpr.add_argument("--radius", type=int, default=1000, help=cpr_radius_help)   # radius

    # coverage_convolution (coverage_convolution) arguments ------------------------------------------------------------

    parser_ccv.add_argument("infile", type=str, help=ccv_infile_help)                       # input file
    parser_ccv.add_argument("outfile", type=str, help=ccv_outfile_help)                     # output file
    parser_ccv.add_argument("--radius", type=int, help=ccv_radius_help, default=50)         # radius
    parser_ccv.add_argument("--eval_length", type=int, help=ccv_eval_help, default=10**6)   # evaluation length

    # covered_regions (covered_regions) help strings -------------------------------------------------------------------

    parser_cvr.add_argument("infile", type=str, help=cvr_infile_help)                   # input file
    parser_cvr.add_argument("outfile", type=str, help=cvr_outfile_help)                 # output file
    parser_cvr.add_argument("--radius", type=int, help=cvr_radius_help, default=50)     # radius

    # denovo_readtable (new_read_table) arguments ----------------------------------------------------------------------

    parser_drt.add_argument("mode", type=str, help=drt_mode_help)               # mode
    parser_drt.add_argument("bamfile", type=str, help=drt_bam_help)             # BAM file
    parser_drt.add_argument("reffile", type=str, help=drt_ref_help)             # reference file
    parser_drt.add_argument("bedfile", type=str, help=drt_bed_help)             # BED file
    parser_drt.add_argument("outfile", type=str, help=drt_out_help)             # output file
    parser_drt.add_argument("snpfile", type=str, help=drt_snp_help)             # SNP file
    parser_drt.add_argument("--cov", type=int, help=drt_cov_help, default=1)    # coverage

    # denovo_rsummary (snp_corrected_rsummary) arguments ---------------------------------------------------------------

    parser_drs.add_argument("bamfile", type=str, help=drs_bam_help)             # BAM file
    parser_drs.add_argument("reffile", type=str, help=drs_ref_help)             # reference file
    parser_drs.add_argument("snpfile", type=str, help=drs_snpfile_help)         # SNP file
    parser_drs.add_argument("outfile", type=str, help=drs_outfile_help)         # output file
    parser_drs.add_argument("--editfile", type=str, help=drs_editfile_help)     # edit sites file

    # filter_rsummary (antisense_sense_rsummary) arguments -------------------------------------------------------------

    parser_frs.add_argument("bamfile", type=str, help=frs_bam_help)                 # BAM file
    parser_frs.add_argument("rsumfile", type=str, help=frs_rsumfile_help)           # read summary input file
    parser_frs.add_argument("bedfile", type=str, help=frs_bed_help)                 # BED file
    parser_frs.add_argument("antisense_outfile", type=str, help=frs_aout_help)      # antisense output file
    parser_frs.add_argument("sense_outfile", type=str, help=frs_sout_help)          # sense output file
    parser_frs.add_argument("--cov", type=int, help=frs_cov_help, default=1)        # coverage

    # rsummary_stats (rsummary_stats) arguments ------------------------------------------------------------------------

    parser_rss.add_argument("rsumfile", type=str, help=rss_rsumfile_help)       # read summary input file
    parser_rss.add_argument("outfile", type=str, help=rss_outfile_help)         # output file help

    # denovo_cctable (new_convcount_table) arguments -------------------------------------------------------------------

    parser_dcc.add_argument("mode", type=str, help=dcc_mode_help)                           # mode help
    parser_dcc.add_argument("outfile", type=str, help=dcc_out_help)                         # output file
    parser_dcc.add_argument("input_files", type=str, nargs="*", help=dcc_infiles_help)      # read information files

    # merge_cctable (merged_convcount_table) arguments -----------------------------------------------------------------

    parser_mcc.add_argument("outfile", type=str, help=mcc_out_help)                         # output file
    parser_mcc.add_argument("input_files", type=str, nargs="*", help=mcc_cctables_help)     # cc table input files

    # denovo_convefftable (new_conveff_table) arguments ----------------------------------------------------------------

    parser_dce.add_argument("cctable_file", type=str, help=dce_cctable_help)                # conversion counts table
    parser_dce.add_argument("outfile", type=str, help=dce_outfile_help)                     # output file
    parser_dce.add_argument("base_error", type=float, help=dce_base_help)                   # base error
    parser_dce.add_argument("seq_error", type=float, help=dce_seq_help)                     # sequencing error
    parser_dce.add_argument("--fix_ce", type=float, help=dce_fixce_help, default=None)      # fixed conv. eff.
    parser_dce.add_argument("--new_ini", type=float, help=dce_newini_help, default=None)    # new initial value
    parser_dce.add_argument("--cov", type=int, help=dce_cov_help, default=1)                # coverage
    parser_dce.add_argument("--iter", type=int, help=dce_iter_help, default=100)            # number of iterations

    # denovo_summarytable (new_summary_table) arguments ----------------------------------------------------------------

    parser_dst.add_argument("outfile", type=str, help=dst_out_help)                 # output file
    parser_dst.add_argument("description", type=str, help=dst_descr_help)           # data description
    parser_dst.add_argument("libsize", type=str, help=dst_libsize_help)             # library size
    parser_dst.add_argument("cctable_file", type=str, help=dst_cctable_help)        # cc table file
    parser_dst.add_argument("--cetable_file", type=str, help=dst_cetable_help,
                            default=None)                                           # ce table file
    parser_dst.add_argument("--cov", type=int, help=dst_cov_help, default=1)        # coverage

    # merge_summarytable (merged_summary_table) arguments --------------------------------------------------------------

    parser_mst.add_argument("mode", type=str, help=mst_mode_help)                       # mode
    parser_mst.add_argument("outfile", type=str, help=mst_out_help)                     # output file
    parser_mst.add_argument("input_files", type=str, nargs="*", help=mst_in_help)       # summary table input files

    # preprocess_reads (preprocess_reads) arguments --------------------------------------------------------------------

    parser_ppr.add_argument("bamfile", type=str, help=ppr_bamfile_help)         # BAM file
    parser_ppr.add_argument("reffile", type=str, help=ppr_reffile_help)         # reference file
    parser_ppr.add_argument("snpfile", type=str, help=ppr_snpfile_help)         # SNP file
    parser_ppr.add_argument("bedfile", type=str, help=ppr_bedfile_help)         # BED file
    parser_ppr.add_argument("rsumfile", type=str, help=ppr_rsumfile_help)       # read summary file
    parser_ppr.add_argument("antifile", type=str, help=ppr_antifile_help)       # antisense read summary file
    parser_ppr.add_argument("sensefile", type=str, help=ppr_sensefile_help)     # sense read summary file help
    parser_ppr.add_argument("recfile", type=str, help=ppr_reffile_help)         # recorded reads' BAM file
    parser_ppr.add_argument("sparefile", type=str, help=ppr_sparefile_help)     # spare reads' BAM file
    parser_ppr.add_argument("ccfile", type=str, help=ppr_ccfile_help)           # conversion counts file
    parser_ppr.add_argument("--editfile", type=str, help=ppr_editfile_help,
                            default=None)                                       # editing sites file

    # parsing arguments ------------------------------------------------------------------------------------------------

    args = parser.parse_args()

    ####################################################################################################################
    # running functions
    ####################################################################################################################

    # analysis to be performed -----------------------------------------------------------------------------------------

    analysis = args.analysis

    # bed_annotation (bed_annotations) ---------------------------------------------------------------------------------

    if analysis == "bed_annotation":

        bed_annotations(args.bedfile, args.gff3file, args.exonfile, args.intronfile, args.utr3file, args.utr5file,
                        args.outfile, bam_files=args.bamfiles, bam_des=args.bamdes)

    # bed_annotation_statistics (annotation_statistics) ----------------------------------------------------------------

    elif analysis == "bed_annotation_statistics":

        annotation_statistics(args.infile, args.outfile)

    # denovo_editsites (write_edit_sites) ------------------------------------------------------------------------------

    elif analysis == "denovo_editsites":

        write_edit_sites(args.rtablefile, args.cutoff, args.outfile)

    # filter_spare_reads (filter_spare_reads) --------------------------------------------------------------------------

    elif analysis == "filter_spare_reads":

        filter_spare_reads(args.bam_in, args.rsum_files, args.bam_out_rec, args.bam_out_spare)

    # covered_positions (covered_positions) ----------------------------------------------------------------------------

    elif analysis == "covered_positions":

        cpos = covered_positions(args.bam_in)
        write_covered_positions(cpos, args.outfile)

    # merge_covered_positions (covered_positions_merge) ----------------------------------------------------------------

    elif analysis == "merge_covered_positions":

        covered_positions_merge(args.input_files, args.outfile)

    # coverage_profile (coverage_profile) ------------------------------------------------------------------------------

    elif analysis == "coverage_profile":

        cov_hash = load_covered_positions(args.infile)
        cov_profile = coverage_profile(cov_hash, args.radius)
        write_coverage_profile(cov_profile, args.outfile)

    # coverage_convolution (coverage_convolution) ----------------------------------------------------------------------

    elif analysis == "coverage_convolution":

        cov_hash = load_covered_positions(args.infile)
        convolution_hash = coverage_convolution(cov_hash, radius=args.radius, eval_length=args.eval_length)
        write_coverage_convolution(convolution_hash, args.outfile)

    # covered_regions (covered_regions) --------------------------------------------------------------------------------

    elif analysis == "covered_regions":

        covered_regions(args.infile, args.outfile, radius=args.radius)

    # denovo_readtable (new_read_table) --------------------------------------------------------------------------------

    elif analysis == "denovo_readtable":

        new_read_table(args.mode, args.bamfile, args.reffile, args.bedfile, args.outfile,
                       snpfile=args.snpfile, coverage=args.cov)

    # denovo_rsummary (snp_corrected_rsummary) -------------------------------------------------------------------------

    elif analysis == "denovo_rsummary":

        snp_corrected_rsummary(args.bamfile, args.reffile, args.snpfile, args.outfile, editfile=args.editfile)

    # filter_rsummary (antisense_sense_rsummary) arguments -------------------------------------------------------------

    elif analysis == "filter_rsummary":

        antisense_sense_rsummary(args.bamfile, args.rsumfile, args.bedfile, args.antisense_outfile, args.sense_outfile,
                                 coverage=args.cov)

    # rsummary_stats (rsummary_stats) ----------------------------------------------------------------------------------

    elif analysis == "rsummary_stats":

        rsummary_stats(args.rsumfile, args.outfile)

    # denovo_cctable (new_convcount_table) -----------------------------------------------------------------------------

    elif analysis == "denovo_cctable":

        new_convcount_table(args.mode, args.input_files, args.outfile)

    # merge_cctable (merged_convcount_table) ---------------------------------------------------------------------------

    elif analysis == "merge_cctable":

        merged_convcount_table(args.input_files, args.outfile)

    # denovo_convefftable (new_conveff_table) --------------------------------------------------------------------------

    elif analysis == "denovo_convefftable":

        new_conveff_table(args.cctable_file, args.outfile, args.base_error, args.seq_error,
                          fix_efficiency=args.fix_ce, new_ini=args.new_ini, coverage=args.cov, iterations=args.iter)

    # denovo_summarytable (new_summary_table) --------------------------------------------------------------------------

    elif analysis == "denovo_summarytable":

        new_summary_table(args.outfile, args.description, args.libsize, args.cctable_file, args.cetable_file,
                          coverage=args.cov)

    # merge_summarytable (merged_summary_table) ------------------------------------------------------------------------

    elif analysis == "merge_summarytable":

        merged_summary_table(args.mode, args.input_files, args.outfile)

    # preprocess_reads (preprocess_reads) ------------------------------------------------------------------------------

    elif analysis == "preprocess_reads":

        preprocess_reads(args.bamfile, args.reffile, args.snpfile, args.bedfile, args.rsumfile, args.antifile,
                         args.sensefile, args.recfile, args.sparefile, args.ccfile, editfile=args.editfile)


