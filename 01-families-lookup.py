# Input:
#
# 1. ukb_rel.dat - Kinship file fetched from the UK BioBank Bulk folder. Its contents are described at  https://biobank.ndph.ox.ac.uk/ukb/ukb/docs/ukb_genetic_data_description.txt
# 2. ukbiobank_agesex.csv - Age and sex file for each of the samples in the UK BioBank (fields are: eid,p22001,p34)
# 3. all-ids.csv - A list of all the sample IDs in the UK BioBank
#
# Output:
#
# 1. uk-trios.txt - Tab Delimited file of the trios found: Mother, Father, Offspring, Siblings, Year of birth of the three
# 2. uk-offsprings-samples.txt - Bash array format of the offsprings and parents

import csv

# These samples are ignored as their matching data was not found in the UK BioBank
# or potentially due to other data quality issues.
IGNORED_OFFSPRINGS = [1157301, 1450512, 2801593, 3198967, 4268248]
IGNORED_PARENTS = [1228021, 2023070, 2238694, 2388375, 2417693, 2448536, 2972106, 3016733, 3249750, 3283001, 3894815, 4175152, 4568364, 4887634, 4952183, 5544932, 5839148, 5885321, 5918382]

sexYOB = {}  # Global variable to store sex and Year of Birth (YOB) data for each sample ID.
# Key: sample ID (int), Value: [sex (str), YOB (int)]


def load_sex_yob():
    """Loads sex and year of birth (YOB) data from 'ukbiobank_agesex.csv'.

    The data is stored in the global `sexYOB` dictionary.
    Handles cases where sex or YOB might be missing in the input file.
    """
    global sexYOB
    sexYOB = {}  # Initialize or clear the dictionary
    with open('ukbiobank_agesex.csv', 'r') as fp:
        reader = csv.reader(fp)
        next(reader)  # Skip header row (e.g., 'eid,p22001,p34')
        for row in reader:
            id, sex, YOB = row
            # Store sex as the first character (e.g., 'M' or 'F') or empty string if missing.
            # Store YOB as an integer, or 0 if missing.
            sexYOB[int(id)] = [sex[0] if sex else "", int(YOB) if YOB else 0]  # Modified line: handle empty YOB


def analyze_relatedness(patients_fn, output_fn, samples_fn):
    """
    Analyzes relatedness data from 'ukb_rel.dat' to identify parent-offspring relationships
    and subsequently, trios (mother, father, offspring).

    Args:
        patients_fn (str): Filename of the CSV containing all sample IDs (e.g., 'all-ids.csv').
        output_fn (str): Filename for the output tab-delimited file of trios (e.g., 'uk-trios.txt').
        samples_fn (str): Filename for the output bash array of offsprings and parents (e.g., 'uk-offsprings-samples.txt').
    """
    global sexYOB
    patient_ids = {}  # Dictionary to store potential relatives for each patient ID.
    # Key: patient ID, Value: list of [relative_id, sex, age_diff, type, ibs0, kinship]

    # Load all patient IDs from the specified file.
    with open(patients_fn, 'r') as fp:
        for line in fp:
            stripped = line.strip()
            id_val = int(stripped) if stripped.isnumeric() else 0
            if id_val:
                patient_ids[id_val] = []  # Initialize an empty list for relatives

    # Process the relatedness data file ('ukb_rel.dat').
    # This file contains pairs of related individuals and their genetic similarity metrics.
    with open('ukb_rel.dat', 'r') as fp:
        for line in fp:
            ID1, ID2, HetHet, IBS0, Kinship = line.strip().split()

            if not ID1.isnumeric() or not ID2.isnumeric():
                continue

            id1_int = int(ID1)
            id2_int = int(ID2)
            ibs0_float = float(IBS0)
            kinship_float = float(Kinship)

            # If ID1 is in our list of patients, add ID2 as a potential relative.
            if id1_int in patient_ids:
                add_relative(id1_int, id2_int, ibs0_float, kinship_float, patient_ids, sexYOB)

            # If ID2 is in our list of patients, add ID1 as a potential relative.
            # This ensures relationships are captured regardless of which ID appears first.
            if id2_int in patient_ids:
                add_relative(id2_int, id1_int, ibs0_float, kinship_float, patient_ids, sexYOB)

    # Sort patient_ids by ID for consistent output, though not strictly necessary for logic.
    patient_ids = dict(sorted(patient_ids.items()))

    # Open the output file for trios and write the header.
    with open(output_fn, 'w') as fp:
        fp.write('\t'.join(['Mother', 'Father', 'Offspring', 'Siblings', 'YOB M', 'YOB F', 'YOB C']) + '\n')

        offsprings = {}  # Dictionary to store identified offsprings and their parents.
        # Key: offspring_id, Value: {'M': mother_id, 'F': father_id}

        # Iterate through each patient and their identified relatives.
        for id_val, arr in patient_ids.items():  # id_val is the patient, arr is list of their relatives
            for row in arr:  # row contains [relative_id, sex, age_diff, type, ibs0, kinship]
                # Check if the relationship type is 'parent-offspring'.
                if row[3] == 'parent-offspring':
                    # Determine who is the offspring and who is the parent based on age difference.
                    # age_diff is calculated as relative_YOB - patient_YOB.
                    # If age_diff < 0, relative is younger (offspring), patient is older (parent).
                    if int(row[2]) < 0:  # row[2] is age_diff
                        offspring = id_val
                        parent = row[0]  # row[0] is relative_id
                    else:
                        offspring = row[0]
                        parent = id_val

                    if offspring in IGNORED_OFFSPRINGS or parent in IGNORED_PARENTS:
                        continue

                    # Get the sex of the identified parent.
                    parent_sex = sexYOB[parent][0]

                    # Initialize the offspring in the offsprings dictionary if not already present.
                    if offspring not in offsprings:
                        offsprings[offspring] = {'M': '', 'F': ''}  # M for Mother, F for Father
                    # Assign the parent ID to the corresponding sex (Mother or Father).
                    offsprings[offspring][parent_sex] = parent

        # Dictionary to collect unique parent and offspring IDs for the samples file.
        samples = {'p': {}, 'c': {}}  # 'p' for parents, 'c' for offsprings

        # Iterate through the identified offsprings to find complete trios.
        for id_val, row in offsprings.items():  # id_val is offspring_id, row is {'M': mother_id, 'F': father_id}
            # Check if both a mother and a father have been identified for the offspring.
            if row['M'] and row['F']:
                siblings = 0
                # Count siblings: other offsprings with the same mother and father.
                for id2, row2 in offsprings.items():
                    if id_val != id2 and row2['M'] == row['M'] and row2['F'] == row['F']:
                        siblings += 1

                # Add parents and offspring to the samples list. Using dict keys for uniqueness.
                samples['p'][int(row['M'])] = 1
                samples['p'][int(row['F'])] = 1
                samples['c'][int(id_val)] = 1
                # Write the trio information to the output file.
                fp.write('\t'.join([
                    str(row['M']), str(row['F']), str(id_val), str(siblings),
                    str(sexYOB[row['M']][1]), str(sexYOB[row['F']][1]), str(sexYOB[id_val][1])
                ]) + '\n')

    # Prepare the content for the samples file (bash array format).
    samples['p'] = sorted(list(samples['p'].keys()))
    samples['c'] = sorted(list(samples['c'].keys()))
    samples_str = 'SAMPLES_PARENTS=(' + ' '.join(map(str, samples['p'])) + ')\n' + \
                  'SAMPLES_OFFSPRINGS=(' + ' '.join(map(str, samples['c'])) + ')\n'

    # Write the samples data to the specified file.
    with open(samples_fn, 'w') as file:
        file.write(samples_str)


def add_relative(patient_id, relative_id, ibs0, kinship, patient_ids, sex_yob):
    """
    Determines the type of relationship between a patient and a relative based on
    kinship coefficient and IBS0 (Identity By State 0) score.
    Adds the relative's information to the patient's list if a relevant relationship is found.

    Args:
        patient_id (int): The ID of the patient.
        relative_id (int): The ID of the relative.
        ibs0 (float): IBS0 score between the patient and relative.
        kinship (float): Kinship coefficient between the patient and relative.
        patient_ids (dict): The main dictionary storing patient IDs and their relatives.
        sex_yob (dict): Dictionary containing sex and YOB for all individuals.
    """
    type_val = ''  # Initialize relationship type

    # Determine relationship type based on kinship coefficient.
    if kinship > 0.354:
        type_val = 'twin'  # Or identical individuals
    # For first-degree relatives (parent-offspring or siblings)
    if 0.177 <= kinship <= 0.354:
        # Differentiate between parent-offspring and siblings using IBS0.
        if ibs0 < 0.0012:
            type_val = 'parent-offspring'
        else:
            type_val = 'sibling'

    # If a relevant relationship type is identified:
    if type_val:
        # Calculate age difference if YOB data is available for both individuals.
        if relative_id in sex_yob and patient_id in sex_yob:
            # Age difference = YOB_relative - YOB_patient
            # Positive if relative is younger, negative if relative is older.
            age_diff_val = sex_yob[relative_id][1] - sex_yob[patient_id][1]
            # Ensure YOBs are not 0 (missing) for a meaningful age difference.
            if sex_yob[relative_id][1] == 0 or sex_yob[patient_id][1] == 0:
                age_diff_str = ''  # Indicate missing age diff if YOB is unknown
            else:
                age_diff_str = str(age_diff_val)
        else:
            age_diff_str = ''  # Age difference unknown if YOB data is missing for one or both

        # Append relative's information to the patient's list.
        # [relative_id, sex_of_relative, age_difference_str, relationship_type, ibs0, kinship]
        patient_ids[patient_id].append([
            relative_id,
            sex_yob.get(relative_id, [''])[0],  # Get sex, default to empty string if not found
            age_diff_str,
            type_val,
            ibs0,
            kinship
        ])


# Main execution
if __name__ == "__main__":
    # Load sex and Year of Birth data first, as it's needed by other functions.
    load_sex_yob()
    # Perform the main relatedness analysis and generate output files.
    analyze_relatedness('all-ids.csv', 'uk-trios.txt', 'uk-offsprings-samples.txt')
