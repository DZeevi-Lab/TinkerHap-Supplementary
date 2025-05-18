import pysam
import argparse
import os
import time
import sys
import pickle

# Application metadata
APP_NAME = "TruthMaker"
APP_VERSION = "1.0"
APP_DATE = "2024/06/30"
APP_URL = "https://github.com/DZeevi-Lab/"


def isHomozygote(alleles: list) -> bool:
    """Helper function to check if a genotype is homozygous"""
    return (alleles[0] == alleles[1])


def isHeterozygote(alleles: list) -> bool:
    """Helper function to check if a genotype is heterozygous"""
    return (alleles[0] != alleles[1])


def fill_del(allele: str, ref: str) -> str:
    """Pad deletion allele format for normalization"""
    if (allele == ref) and (len(ref) > 1):
        return allele[0]
    else:
        return allele.ljust(1+len(ref)-len(allele), ".")


def fill_dels(alleles: list, ref: str) -> list:
    """Normalize both alleles using fill_del"""
    return [fill_del(alleles[0], ref), fill_del(alleles[1], ref)]


class AnalyzeOffspring:
    """Main class to process parent and offspring VCFs and build truth dataset"""
    vcf_suffix: str = "_24053_0_0.vcf.gz"  # Suffix of VCF files for samples
    vcf_path: str = ""
    tab_path: str = ""
    output_path: str = ""

    def log(self, text: str, print_out: bool = True) -> None:
        """Outputs a text to the log file, and optionally print it to stdout"""
        if not hasattr(self, "time_started"):
            self.time_started = time.time()
        text = str(round(time.time() - self.time_started))+"] " + text
        self.logFile.write(text + "\n")
        if (print_out):
            print(text)

    def phase_by_homozygote(self, parents, offspring, parentIndex) -> list:
        """Infer offspring phase using homozygous parent genotypes"""
        if not isHomozygote(parents[parentIndex]):
            return []
        phased = ["", ""]
        otherParent = 1 if parentIndex == 0 else 0
        gt = parents[parentIndex][0]
        if (offspring[0] == gt) and (offspring[1] in parents[otherParent]):
            phased[otherParent] = offspring[1]
        elif (offspring[1] == gt) and (offspring[0] in parents[otherParent]):
            phased[otherParent] = offspring[0]
        else:
            return []
        phased[parentIndex] = gt
        return phased

    def parse_offspring_vcf(self, offspring_vcf_path: str) -> dict:
        """Parses the VCF file of the offspring and builds truth dataset based on parents"""
        self.truth = {}
        with pysam.VariantFile(offspring_vcf_path) as reader:
            for record in reader:
                if (record.qual == 0):
                    continue
                if (len(record.samples) != 1):
                    self.log("ERROR: Too many samples in file", True)
                    sys.exit(1)
                if (not record.pos in self.parents):
                    continue
                offspring = record.samples[0].alleles
                if (isHomozygote(offspring)):
                    continue
                offspring = fill_dels(offspring, record.ref)
                parents = self.parents[record.pos]
                phased = self.phase_by_homozygote(parents, offspring, 0)
                if (len(phased) == 0):
                    phased = self.phase_by_homozygote(parents, offspring, 1)
                    if (len(phased) == 0):
                        continue
                self.truth[record.pos] = phased
        self.log("  extracted "+str(len(self.truth))+" locations", True)

    def read_parent_vcf(self, vcf_path: str) -> dict:
        """Parses a parent VCF file and returns a dictionary of allele calls"""
        result = {}
        with pysam.VariantFile(vcf_path) as reader:
            for record in reader:
                if (record.qual == 0):
                    continue
                if (len(record.samples) != 1):
                    self.log("ERROR: Too many samples in file", True)
                    sys.exit(1)
                alleles = record.samples[0].alleles
                if (alleles[0] == None) or (alleles[1] == None):
                    continue
                result[record.pos] = [fill_del(alleles[0], record.ref), fill_del(alleles[1], record.ref), record.ref[0]]
        self.log("  extracted "+str(len(result))+" locations", True)
        return result

    def add_parent(self, parent1, parent2, pos) -> None:
        """Adds a variant position to the parent dictionary if it's phase-informative"""
        if ((isHeterozygote(parent1) and isHomozygote(parent2)) or
                (isHeterozygote(parent2) and isHomozygote(parent1)) or
                (isHomozygote(parent2) and isHomozygote(parent1) and (parent2[0] != parent1[0]))):
            self.parents[pos] = [parent1[0:2], parent2[0:2]]

    def parse_parents(self, parent1_vcf_path: str, parent2_vcf_path: str) -> None:
        """Parses both parent VCFs and builds a parent allele dictionary for shared positions"""
        self.parent1 = self.read_parent_vcf(parent1_vcf_path)
        self.parent2 = self.read_parent_vcf(parent2_vcf_path)
        self.parents = {}
        for pos in self.parent1:
            parent1 = self.parent1[pos]
            parent2 = self.parent2[pos] if pos in self.parent2 else [parent1[2], parent1[2]]
            self.add_parent(parent1, parent2, pos)
        for pos in self.parent2:
            if pos in self.parent1:
                continue
            parent2 = self.parent2[pos]
            parent1 = [parent2[2], parent2[2]]
            self.add_parent(parent1, parent2, pos)

    def save_output(self, output_fn) -> None:
        """Saves the truth dataset to a pickle file"""
        self.log("  saving "+str(len(self.truth))+" samples to output "+output_fn, True)
        with open(self.output_path+"/"+output_fn, "wb") as outfile:
            pickle.dump(self.truth, outfile)

    def parse_offspring(self, parent1_id: str, parent2_id: str, offspring_id: str):
        """Full process for a single trio (mother, father, offspring)"""
        self.log(f"Parsing {offspring_id} ({parent1_id},{parent2_id})")
        self.parse_parents(self.vcf_path+"/"+parent1_id+self.vcf_suffix, self.vcf_path+"/"+parent2_id+self.vcf_suffix)
        self.parse_offspring_vcf(self.vcf_path+"/"+offspring_id+self.vcf_suffix)
        self.save_output(offspring_id+".pkl")

    def go(self):
        """Main entry method to iterate through all trios listed in the input file"""
        self.logFile = open("analysis.log", "w")
        with open(self.tab_path, 'r') as file:
            for line in file:
                row = line.strip().split('\t')
                if (len(row) < 3):
                    continue
                parent1_id = row[0]
                parent2_id = row[1]
                offspring_id = row[2]
                if not parent1_id.isdigit():
                    continue
                self.parse_offspring(parent1_id, parent2_id, offspring_id)
        self.logFile.close()


def arg_file_in(x):
    """Validation helper for argparse: check if a file path is valid"""
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError(f"File '{x}' doesn't exist")
    return x


def main():
    parser = argparse.ArgumentParser(description=APP_NAME + " - parents data extractor. Version "+APP_VERSION+" ("+APP_DATE+")",
                                     epilog="For more info visit "+APP_URL)

    parser.add_argument("-i", "--input",
                        dest="tab_path",
                        type=arg_file_in,
                        required=True,
                        help="List in tabbed text format containing at least 3 columns: Mother, Father, Offspring")
    parser.add_argument("-v", "--vcf",
                        dest="vcf_path",
                        type=arg_file_in,
                        required=True,
                        help="Input VCFs path")
    parser.add_argument("-o", "--output",
                        dest="output_path",
                        required=True,
                        help="Output path")

    args = vars(parser.parse_known_args()[0]).items()

    obj = AnalyzeOffspring()
    for param, value in args:
        setattr(obj, param, value)

    obj.go()


if __name__ == "__main__":
    main()
