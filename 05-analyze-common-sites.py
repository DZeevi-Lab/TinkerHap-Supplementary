import sys
import os
import pickle

# ------------------------
# Load all *.pkl files into a nested dictionary: detailed[sample][method][position] = result
# These files contain per-sample, per-method phasing result dictionaries.
# ------------------------

detailed = {}
for fn in os.listdir('.'):
    if not fn.endswith('.pkl'):
        continue

    method_name = os.path.splitext(fn)[0]  # Extract method name from filename
    with open(fn, "rb") as myfile:
        data = pickle.load(myfile)

    for sample_name in data:
        if not sample_name in detailed:
            detailed[sample_name] = {}
        detailed[sample_name][method_name] = data[sample_name]

# ------------------------
# Decide method set based on whether PacBio flag was passed
# ------------------------
if (len(sys.argv) == 2) and (sys.argv[1] == '--pacbio'):
    methods = ["TinkerHap", "WhatsHap", "HapCUT2"]
else:
    methods = ["TinkerHap+SI", "TinkerHap", "WhatsHap", "HapCUT2", "ShapeIT"]

title = "\t".join(methods)

# ------------------------
# Prepare output file
# ------------------------
list_file = open("common-sites.txt", "w")
list_file.write("Sample\tSites\t"+title+"\n")

# Assign method subsets
common_methods = methods
first_method = common_methods[0]
other_methods = common_methods[1:]

# ------------------------
# Main comparison loop: iterate over each sample
# ------------------------
for sample in detailed:
    print(sample, end="      \r")  # Display progress in terminal

    matches = {}  # Track number of correct phasings per method for shared sites
    missing_method = False

    # Ensure that all methods are present for the current sample
    for method in common_methods:
        if method not in detailed[sample]:
            missing_method = True
            break
        matches[method] = 0

    if missing_method:
        continue

    count = 0  # Number of common positions
    for pos in detailed[sample][first_method]:
        common = True
        pos_matches = {first_method: 1 if detailed[sample][first_method][pos] else 0}

        for method in other_methods:
            if pos not in detailed[sample][method]:
                common = False
                break

            if method not in pos_matches:
                pos_matches[method] = 0

            if detailed[sample][method][pos]:
                pos_matches[method] += 1

        if not common:
            continue

        count += 1
        res = []
        for method in common_methods:
            res.append(str(pos_matches[method]))
            matches[method] += pos_matches[method]

    # ------------------------
    # Calculate and store error rate for each method
    # ------------------------
    output = {}
    for method in matches:
        output[method] = (count - matches[method]) / count if count > 0 else 0

    # Write output line to file
    list_line = ""
    for method in common_methods:
        list_line += "\t" + str(output[method])
    list_file.write(f"{sample}\t{count}{list_line}\n")

# Close output file
list_file.close()
