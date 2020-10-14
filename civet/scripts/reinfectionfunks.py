import csv
import os
from collections import defaultdict

from reportfunk.funks import io_functions as qcfunk

def configure_input_query(config):

    query_metadata = config["query"]
    patient_id_col = config["patient_id_col"]
    path_to_dir = config["outdir"]

    new_query = os.path.join(path_to_dir, "query_with_groups.csv") #this is not actually where we want it - tempdir maybe?

    patient_to_group = {}
    count = 1
    single_line = False
    with open(query_metadata) as f:
        reader = csv.DictReader(f)
        fieldnames=reader.fieldnames
        data = [r for r in reader]
        for seq in data:
            
            if len(fieldnames) == 1:
                patient = "Patient A"
                single_line = True
            else:
                patient = seq[patient_id_col]

            if patient not in patient_to_group.keys():
                patient_to_group[patient] = count
                count += 1

    
    with open(query_metadata) as csv_input:
        with open(new_query, "w") as csv_output:
            reader = csv.DictReader(csv_input)
            data = [r for r in reader]
            
            if single_line:
                fieldnames.append("patient_id")
            fieldnames.append("group")
            
            writer = csv.DictWriter(csv_output, fieldnames = fieldnames)
            writer.writeheader()
            
            for row in data:
                if single_line:
                    row["patient_id"] = patient
                    row["group"] = patient_to_group[patient]
                else:
                    row["group"] = patient_to_group[row[patient_id_col]]
                
                writer.writerow(row)

    config["query"] = new_query

def reinfection_args_to_config(args, config, defaultdict):

    patient_id_col = qcfunk.check_arg_config_default("patient_id_col", args.patient_id_col, config, defaultdict)
    config["patient_id_col"] = patient_id_col

def sort_into_patients(query_dict):

    patient_to_taxa = defaultdict(list)
    patient_to_tree = defaultdict(set)
    patient_list = set()
    
    for tax in query_dict.values():
        
        patient = tax.attribute_dict["patient"]
        patient_to_taxa[patient].append(tax)
        patient_to_tree[patient].add(tax.tree)
        patient_list.add(patient)


    patient_list = list(patient_list)

    return patient_list, patient_to_taxa, patient_to_tree


    
    