__author__ = 'jgwall'

import argparse
import gzip
import numpy as np
from scipy.stats import binom_test
from warnings import warn

debug = False
het_calls = ["0/1", "1/0"]  # Not bothering dealing with tertiary allele states ("0/2", etc)
changed_count, total_hets = 0, 0  # Keep track of how many locations got changed

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help = "VCF file to have hets corrected (can be gzipped)")
    parser.add_argument("-o", "--outfile", help="Output VCF file (optionally gzipped)")
    parser.add_argument("-p", "--p-cutoff", type=float, default=0.01, help="P-value cutoff to change het calls. Het calls with p <= this will be changed")
    parser.add_argument("-m", "--min-count", type=int, default=3,
                        help="Minimum count for a het call. Only het calls with fewer than this many reads in the minor allele will be changed")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def main():
    args = parse_args()
    print("Correcting het calls in", args.infile,"based on binomial distribution with p-value cutoff of", args.p_cutoff,"and minimum read count of",args.min_count)
    print("\tNote: this script ignores tertiary allele states for simplicity. Only het calls of '0/1' or '1/0' are corrected.")

    # Get file handles to work with
    IN = get_filehandle(args.infile, "rt")
    OUT = get_filehandle(args.outfile, "wt")

    # Go through VCF and parse each line
    n=0
    for line in IN:

        newline = parse_line(line, p_cutoff = args.p_cutoff, min_count=args.min_count)
        OUT.write(newline)

        # Track number of lines handled
        n += 1
        if n % 100000 == 0:
            print("\tProcessed",n,"lines")
        if debug and n > 1000:
            break

    global changed_count, total_hets
    print("Total", changed_count,"out of",total_hets,"het calls changed (=",changed_count/total_hets,")")

    # Close out file handles
    IN.close()
    OUT.close()


# Parse individual VCF lines
def parse_line(line, p_cutoff, min_count):

    # Return header lines unchanged
    if line.startswith("#"):
        return line

    # Split line into fields and extract info as to where in each field the allele depth is stored
    records = line.strip().split('\t')
    genoID, depthID = get_field_locations(records)

    # Confirm this line actually has the needed info; skip if it doesn't
    if genoID is None:
        warn("Unable to find genotype (GT) field in line", records[:4], "; Skipping")
        return line
    if depthID is None:
        warn("Unable to find allele depth (AD) field in line", records[:4], "; Skipping")
        return line

    # Go through records and correct hets if needed. VCF standards indicate that records start at index 9 (=column 10)
    for i in range(9, len(records), 1):
        records[i] = correct_hets(records[i], genoID=genoID, depthID=depthID, p_cutoff=p_cutoff, min_count=min_count)

    # Reconstitute the line of text and return
    return "\t".join(records) + "\n"


# Take a list of fields (corresponding to a line of a VCF file) and determine where the genotype and allele depth are stored in each record
def get_field_locations(records):
    formatkey = records[8]  # VCF standards set that column 9 (=index 8) is FORMAT
    formatkey = formatkey.split(":")  # Split on delimiter
    genoID = formatkey.index("GT")
    depthID = formatkey.index("AD")
    return genoID, depthID


# Change errant het calls based on binomial probability
def correct_hets(record, genoID, depthID, p_cutoff, min_count, delimiter=":"):
    fields = record.split(delimiter)
    global changed_count, total_hets

    # If this genotype call isn't a het, return it unchanged
    if fields[genoID] not in het_calls:
        #if debug:
            #print("\tRecord",record,"is not heterozygous; returning unchanged")
        return record
    else:
        total_hets+=1

    # Calculate probability of a het call based on these counts
    depths = fields[depthID].split(',')[:2] # Take only first two alleles; above check should filter out any with allele calls outside these
    depths = [int(d) for d in depths] # Convert to numeric
    total_depth = sum(depths)
    p_value = binom_test(x=depths, p=0.5, alternative="two-sided")

    # Debugging help message
    if debug:
        print("\tRecord",record,"has genotype",fields[genoID], ", allele depths",depths, ", total depth", total_depth,
              "and two-sided binomial p-value",p_value)

    # Correct het if needed. If p-value and allele counts too low, just assume it goes to the majority call
    if (p_value < p_cutoff) and min(depths) < min_count:
        majority_allele = np.argmax(depths)
        newgeno = str(majority_allele) + "/" + str(majority_allele)
        fields[genoID] = newgeno
        changed_count +=1 # Keep track of how many changed

        if debug:
            print("\t\tCorrected to genotype call",newgeno,"yielding record", fields)

    # Return values
    return delimiter.join(fields)


# Helper function for oepening files for reading/writing based on file extension
def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    else:
        return open(file,mode)


if __name__ == '__main__': main()