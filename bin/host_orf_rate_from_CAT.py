import csv
import sys
from collections import defaultdict

# yangm@idsse.ac.cn 2024-09-29 17:23 

def parse_classification(classification):
    'Parse a taxonmic ID to judge if it belongs to virus of host.'
    if classification.startswith('1;10239'):  # Virus class IDs start with "1;10239"
        return 'virus'
    elif classification.startswith('1;131567'):  # Host class IDs start with "1;10239"
        return 'host'
    return 'other'

def process_cat_file(input_file, output_file):
    # contig ORF catagory dictionary
    contig_stats = defaultdict(lambda: {'virus': 0, 'host': 0, 'other': 0})
    # Read input file from Contig Annotation Tool (CAT) output file *.ORF2LCA.txt.
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        next(reader, None)
        for row in reader:
            orf_id = row[0]
            # Extract contigID (before the last underscore of the ORF ID)
            contig_id = '_'.join(orf_id.split('_')[:-1])
            category = 'other'
            if len(row) > 2:
                classification = row[2]
                # Parse and categorize to determine whether it's a virus or a cellular organism
                category = parse_classification(classification)
            contig_stats[contig_id][category] += 1
    # output
    with open(output_file, 'w') as outfile:
        outfile.write('Contig_ID\tViral_ORFs\tHost_ORFs\tOther_ORFs\tHost_ORF_Rate(%)\tIs_Virus\n')
        for contig_id, counts in contig_stats.items():
            virus_count = counts['virus']
            host_count = counts['host']
            other_count = counts['other']
            total_count = virus_count + host_count
            host_percentage = (host_count / total_count) * 100 if total_count > 0 else 0
            is_virus = 0 if host_percentage > 40 else 1
            outfile.write(f'{contig_id}\t{virus_count}\t{host_count}\t{other_count}\t{host_percentage:.2f}\t{is_virus}\n')

# process_cat_file function
if __name__=='__main__':
    process_cat_file(sys.argv[1], sys.argv[2])
