#!/usr/bin/env python

from argparse import ArgumentParser
import re, os
from itertools import chain

class GeneRule(object):

    def __init__(self, species, allele, context, drug, expected_pheno):
        self.species = species
        self.allele = allele
        self.context = context
        self.drug = drug
        self.expected_pheno = expected_pheno

class AMRFinderResult(object):

    def __init__(self, allele_id, amr_class, amr_subclass, drug, element_type, element_subtype, scope, sequence_name, protein_onwards, method_onwards, name=None):
        self.allele_id = allele_id
        self.amr_class = amr_class
        self.amr_subclass = amr_subclass
        self.drug = str(drug)
        self.name = name
        self.element_type = element_type
        self.element_subtype = element_subtype
        self.scope = scope
        self.sequence_name = sequence_name
        self.protein_onwards = protein_onwards
        self.method_onwards = method_onwards

        # Assigned during rule processing
        self.expected_pheno = str()
        self.context = str()
    

class OrganismAwareReport(object):
    
    pass

def get_arguments():
    parser = ArgumentParser(description='Parse AMRFinderPlus files with organism specific rules.')
    
    parser.add_argument('--reports', nargs='+', type=str, required=True, help='One or more AMRFinderPlus results files (should all belong to the same species).')
    parser.add_argument('--species', required=False, default='Klebsiella pneumoniae', type=str, help='Species of the genomes in the input files (default is Klebsiella pneumoniae)')
    parser.add_argument('--organism_rules', required=True, type=str, help='Organism-specific rule set table, used to generate the annotated AMR report (all genomes included in a single report).')
    parser.add_argument('--drug_dictionary', required=False, default='./example_dict_kleb/Kleb_local_dict.tsv', help='Path to drug dictionary that matches allele names with drug classes. Current default is the temporary dictionary in this repo, "./example_dict_kleb/Kleb_local_dict.tsv')
    parser.add_argument('--output', required=True, type=str, help='Name for output file.')

    return parser.parse_args()

def create_drug_list(drug_gene_file):

    # initialise dictionary, key=allele, value=drug name
    drug_dict = {}

    with open(drug_gene_file, 'r') as drug_genes:
        header = 0
        for line in drug_genes:
            fields = line.strip().split('\t')
            if header == 0:
                header +=1 # ignore header
            else:
                if fields[0] == "NA":   
                    drug_dict[fields[1]] = str.lower(fields[3])
                else:
                    drug_dict[fields[0]] = str.lower(fields[3])

    return drug_dict

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

def parse_amr_report(report_file, drug_dict):

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
                element_subtype_col = col_headers.index("Element subtype")
                scope_col = col_headers.index("Scope")
                seq_name_col = col_headers.index("Sequence name")
                # we want all the columns from the original AMRFinder output, but it's going to be annoying to save them all as individual values inside the class object
                # so instead let's just get everything after 'name' (or 'Protein identifier' if name not present) up to 'Gene symbol' and save that as one item
                protein_id_col = col_headers.index("Protein identifier")
                # and then let's grab everything after 'Subclass' and save that as it's own single entry
                method_col = col_headers.index("Method")
                header += 1
            else:
                # only parse the line if the element_type is AMR
                if fields[element_type_col] == "AMR":
                    gene_allele = fields[gene_symbol_col]
                    # grab the drug name for this allele from the dictionary, if it exists. If not, leave blank
                    try:
                        gene_drug = drug_dict[gene_allele]
                    except KeyError:
                        gene_drug = ''
                    class_type = fields[class_col]
                    subclass_type = fields[subclass_col]
                    # now get all values from protein ID to (but not including) gene symbol
                    protein_id_onwards = fields[protein_id_col:protein_id_col+5]
                    # and now grab the stuff at the end
                    method_onwards = fields[method_col:]
                    if name_col != None:
                        name_id = fields[name_col]
                        amrfinder_report_lines.append(AMRFinderResult(gene_allele, class_type, subclass_type, gene_drug, fields[element_type_col],
                                                                      fields[element_subtype_col], fields[scope_col], fields[seq_name_col], 
                                                                      protein_id_onwards, method_onwards, name=name_id))
                    else:
                        amrfinder_report_lines.append(AMRFinderResult(gene_allele, class_type, subclass_type, gene_drug, fields[element_type_col],
                                                                      fields[element_subtype_col], fields[scope_col], fields[seq_name_col], 
                                                                      protein_id_onwards, method_onwards))
    
    return amrfinder_report_lines

def determine_rules(amrfinder_report_lines, rule_list, sampleID, species):

    # this is our list of result line classes that have the info we want
    output_lines = []

    # extract the rules relevant to the species
    relevant_rules = [rule for rule in rule_list if rule.species == species]

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
                #expanded_result = amrfinder_result.add_to_result(rule.expected_pheno, rule.drug, sampleID, rule.context)

                amrfinder_result.expected_pheno = rule.expected_pheno
                amrfinder_result.drug = rule.drug
                amrfinder_result.name = sampleID
                amrfinder_result.context = rule.context

                # now escape this for loop!!
                break
            # if we're at the final rule, and still no search result then add an empty version, as there is no rule for this allele call
            if (rule_list.index(rule) + 1) == len(rule_list) and not search_value:
                amrfinder_result.name = sampleID
        # add it to our new list
        output_lines.append(amrfinder_result)

    return output_lines

def remove_nested_lists(line_to_write):
    unnested_list = []
    for item in line_to_write:
        if type(item) == list:
            for value in item:
                unnested_list.append(value)
        else:
            unnested_list.append(item)
    return unnested_list

def write_output(output_lines, out_file, species, version):

    with open(out_file, "w") as out:
        # add name to the header if it exists, otherwise don't bother
        header = ['Name', 'Protein identifier', 'Contig id', 'Start', 'Stop', 'Strand', 'Gene symbol', 'Sequence name', 
                  'Species interpretation', 'Context', 'Org interpretation', 'Drug', 'Scope', 'Element type', 'Element subtype', 'Class', 'Subclass',
                  'Method', 'Target length', 'Reference sequence length', '% Coverage of reference sequence', '% Identity to reference sequence',
                  'Alignment length', 'Accession of closest sequence', 'Name of closest sequence', 'HMM id', 'HMM description']
        out.write('\t'.join(header) + '\n')

        # get the value for species and version for the 'Species interpretation' column
        species_interp = species + '; ' + version

        for out_line in output_lines:
            # correctly format the wt resistant/susceptible codes to match poster
            if out_line.expected_pheno == 'wt resistant':
                expected_pheno = 'wt (R)'
            elif out_line.expected_pheno == 'wt susceptible':
                expected_pheno = 'wt (S)'
            else:
                expected_pheno = out_line.expected_pheno
            final_line = [out_line.name, out_line.protein_onwards, out_line.allele_id, out_line.sequence_name,
                        species_interp, out_line.context, expected_pheno, out_line.drug, 
                        out_line.scope, out_line.element_type, out_line.element_subtype, 
                        out_line.amr_class, out_line.amr_subclass, out_line.method_onwards]
            #TODO: see if there is a faster/more memory efficient way to do this
            final_line = remove_nested_lists(final_line)
            out.write('\t'.join(final_line) + '\n')

def organism_aware_report(output_lines, local_drug_list):

    # get the list of drugs we care about, match to their AMRFinder name if needed
    #key: drug name to write in report, value: AMR drug name, if needed, otherwise empty string
    drugs_to_report = []
    drug_amrfinder_name_conversion = {}
    with open(local_drug_list, 'r') as in_file:
        header = 0
        for line in in_file:
            if header == 0:
                header += 1
            else:
                fields = line.strip().split('\t')
                drugs_to_report.append(fields[0])
                if len(fields) == 2:
                    drug_amrfinder_name_conversion[fields[1]] = fields[0]
    
    for entry in output_lines:
        if entry.drug == '' or entry.drug == 'multiple drug':
            drug_interest = entry.amr_subclass
            # check to see the conversion 
        else:
            drug_interest = entry.drug

    return drug_interest

def main():

    args = get_arguments()

    # parse the drug dictionary
    drug_dict = create_drug_list(args.drug_dictionary)

    # parse the organism rule file a list. Each item in the list is an object with the relevant details
    rule_list = create_rule_list(args.organism_rules)


    # to get to the wt(R) thing we want, we're going to need to, for each gene that has a WT R, make a note that this
    #drug is an expected pheno. And then when we get to a call that doesn't have a rule attached to it, check if it's a drug with a known wt(R) expected pheno? but this will depend on the order things happen in
    # kind of need to draw out the class hierachy I think to make this make more sense

    full_output_entries = []

    for report in args.reports:
        amrfinder_results = parse_amr_report(report, drug_dict)
        # if we have a name value in the amrfinder report, use that
        # otherwise just use the name of the report file in the 'Sample' column
        if amrfinder_results[0].name:
            sampleID = amrfinder_results[0].name
        else:
            sampleID = os.path.basename(report)
        output_entries = determine_rules(amrfinder_results, rule_list, sampleID, args.species)

        # sort the lines
        sorted_entries = dict()
        for entry in output_entries:
            # group calls by drug so we can update the expected phenotype as needed
            if entry.drug not in sorted_entries:
                sorted_entries[entry.drug] = list()
            sorted_entries[entry.drug].append(entry)

        # TODO(JH): check out enums
        # check each drug and update to be wt resistant if there is a core gene present
        for drug_name, entries in sorted_entries.items():
            has_wt_r = False
            # Other rules?
            for entry in entries:
                has_wt_r |= entry.expected_pheno == 'wt resistant'

            # Second pass, update information in other entries as required
            for entry in entries:
                if has_wt_r:
                    entry.expected_pheno = 'wt resistant'

        full_output_entries = full_output_entries + output_entries
        
    # now write out the output into a single file
    #TODO: encode version number correctly 
    write_output(full_output_entries, args.output, args.species, 'v1.1')

if __name__ == '__main__':
    main()
