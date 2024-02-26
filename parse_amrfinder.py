#!/usr/bin/env python

from argparse import ArgumentParser
import re, os

class GeneRule(object):
    def __init__(self, species, allele, context, drug, expected_pheno):
        self.species = species
        self.allele = allele
        self.context = context
        self.drug = drug
        self.expected_pheno = expected_pheno

class AMRFinderResult(object):
    def __init__(self, allele_id, amr_class, amr_subclass, name=None):
        self.allele_id = allele_id
        self.amr_class = amr_class
        self.amr_subclass = amr_subclass
        # add the name if the AMRFinder report has it, as this is an optional column when running AMRFinder
        if name:
            self.name = name
    def add_to_result(self, expected_pheno, drug, name, context):
        self.expected_pheno = expected_pheno
        self.drug = drug
        self.name = name # update the name value here - will be what's in the result if it's there, otherwise will stay the same
        self.context = context
        return self

def get_arguments():
    parser = ArgumentParser(description='Parse AMRFinderPlus files with organism specific rules.')
    
    parser.add_argument('--reports', nargs='+', type=str, required=True, help='One or more AMRFinderPlus results files (should all belong to the same species).')
    parser.add_argument('--rules', required=True, type=str, help='Organism specific rule set table.')
    parser.add_argument('--output', required=True, type=str, help='Name for output file.')

    return parser.parse_args()

def create_rule_list(rule_infile):
    rule_list = []
    with open(rule_infile, 'r') as rule_file:
        header = 0
        for line in rule_file:
            if header == 0:
                header += 1
            else:
                fields = line.strip().split('\t')
                new_rule = GeneRule(fields[0], fields[1], fields[2], fields[3], fields[4])
                rule_list.append(new_rule)
    return rule_list

def parse_amr_report(report_file):

    amrfinder_report_lines = []
    
    with open(report_file, 'r') as report:
        header = 0
        for line in report:
            fields = line.strip().split('\t')
            # let's keep the column headers saved, in case we have a name column at the start (people might not always use this option when running AMRFinder)
            if header == 0:
                col_headers = fields
                # get the index for the relevant columns
                gene_symbol_col = col_headers.index('Gene symbol')
                class_col = col_headers.index("Class")
                subclass_col = col_headers.index("Subclass")
                try:
                    name_col = col_headers.index("Name")
                except:
                    name_col = None
                element_type_col = col_headers.index("Element type")
                header += 1
            else:
                # only parse the line if the element_type is AMR
                if fields[element_type_col] == "AMR":
                    gene_allele = fields[gene_symbol_col]
                    class_type = fields[class_col]
                    subclass_type = fields[subclass_col]
                    if name_col != None:
                        name_id = fields[name_col]
                        amrfinder_report_lines.append(AMRFinderResult(gene_allele, class_type, subclass_type, name=name_id))
                    else:
                        amrfinder_report_lines.append(AMRFinderResult(gene_allele, class_type, subclass_type))
    
    return amrfinder_report_lines

def determine_rules(amrfinder_report_lines, rule_list, sampleID):

    # this is our list of result line classes that have the info we want
    output_lines = []

    for amrfinder_result in amrfinder_report_lines:
        # grab the allele allele AMRFinder found
        amrfinder_allele = amrfinder_result.allele_id
        # okay so how do we do this? Do we take the allele and check that it's in our allele rule list?
        # how are we dealing with cases where we have like blaSHV*? because we want to include all alleles, which assumes then some kind of regex matching or substring in x?
        # allele will be blaSHV-28, we want to match blaSHV*
        for rule in rule_list:
            # see if there's a matching rule for that allele
            search_value = re.search(rule.allele, amrfinder_allele)
            # this will return a value if there's something, otherwise it will be None and this won't evaluate
            if search_value:
                # so now we've found something, what do we want to extract? we want to extract the expected phenotype to add to the table, for that rule
                expanded_result = amrfinder_result.add_to_result(rule.expected_pheno, rule.drug, sampleID, rule.context)
                # now escape this for loop!!
                break
            # if we're at the final rule, and still no search result then add an empty version, as there is no rule for this allele call
            if (rule_list.index(rule) + 1) == len(rule_list) and not search_value:
                expanded_result = amrfinder_result.add_to_result("", "", sampleID, context="")
        # add it to our new list
        output_lines.append(expanded_result)

    return output_lines

def write_output(output_lines, out_file):

    with open(out_file, "w") as out:
        # add name to the header if it exists, otherwise don't bother
        header = ['Sample', 'Allele', 'Context', 'Org interpretation', 'Drug', 'Class', 'Subclass']
        out.write('\t'.join(header) + '\n')
        for out_line in output_lines:
            # correctly format the wt resistant/susceptible codes to match poster
            if out_line.expected_pheno == 'wt resistant':
                expected_pheno = 'wt (R)'
            elif out_line.expected_pheno == 'wt susceptible':
                expected_pheno = 'wt (S)'
            else:
                expected_pheno = out_line.expected_pheno
            final_line = [out_line.name, out_line.allele_id, out_line.context, expected_pheno, out_line.drug, out_line.amr_class, out_line.amr_subclass]
            out.write('\t'.join(final_line) + '\n')

def main():

    args = get_arguments()

    # parse the organism rule file a list. Each item in the list is an object with the relevant details
    rule_list = create_rule_list(args.rules)

    full_output_entries = []

    for report in args.reports:
        amrfinder_results = parse_amr_report(report)
        # if we have a name value in the amrfinder report, use that
        # otherwise just use the name of the report file in the 'Sample' column
        if amrfinder_results[0].name:
            sampleID = amrfinder_results[0].name
        else:
            sampleID = os.path.basename(report)
        output_entries = determine_rules(amrfinder_results, rule_list, sampleID)
        full_output_entries = full_output_entries + output_entries
        
    # now write out the output into a single file
    write_output(full_output_entries, args.output)

if __name__ == '__main__':
    main()
