import pysam
import argparse
import os
import time
import sys
import pickle
import csv
import statistics
import gc

# Metadata constants
APP_NAME = "OffspringTester"
APP_VERSION = "0.1"
APP_DATE = "2024/06/30"
APP_URL = "https://github.com/DZeevi-Lab/"

# ----------------- Utility functions -----------------


def isHomozygote(alleles: list) -> bool:
    """Returns True if the two alleles are identical (homozygous)"""
    return (alleles[0] == alleles[1])


def compare_alleles(alleles1: list, alleles2: list) -> bool:
    """Returns True if both alleles match in order"""
    return (alleles1[0] == alleles2[0]) and (alleles1[1] == alleles2[1])


def compare_alleles_inverted(alleles1: list, alleles2: list) -> bool:
    """Returns True if alleles match in reverse order (inverted phase)"""
    return (alleles1[0] == alleles2[1]) and (alleles1[1] == alleles2[0])


def isSNP(alleles: list) -> bool:
    """Returns True if both alleles are single-nucleotide variants (SNPs)"""
    return (len(alleles[0]) == 1) and (len(alleles[1]) == 1)


def alleles_to_str(alleles: list) -> str:
    """Returns a string representation of the allele pair"""
    return alleles[0] + "|" + alleles[1]


def fill_del(allele: str, ref: str) -> str:
    """Pads deletions with dots if allele is shorter than reference"""
    if (allele == ref) and (len(ref) > 1):
        return allele[0]
    else:
        return allele.ljust(1 + len(ref) - len(allele), ".")


def fill_dels(alleles: list, ref: str) -> list:
    """Applies fill_del to both alleles in the pair"""
    return [fill_del(alleles[0], ref), fill_del(alleles[1], ref)]

# ----------------- Main analysis class -----------------


class AnalyzePhased:
    """
    Class for analyzing phased offspring VCFs against truth datasets.
    Evaluates phasing accuracy and outputs error reports.
    """

    # File paths and configuration
    summaryFiles = {}
    vcf_path: str = ""
    truth_path: str = ""
    output_path: str = ""
    ignore_haplotypes: bool = False
    sample_name: str = ""
    method_name: str = ""
    list_path: str = ""

    # Results containers
    output: dict = {}
    errors: dict = {}

    def log(self, text: str, print_out: bool = False) -> None:
        """Outputs a text to the log file, and optionally prints to stdout"""
        if not hasattr(self, "time_started"):
            self.time_started = time.time()
        text = str(round(time.time() - self.time_started)) + "] " + text
        if self.logFile != None:
            self.logFile.write(text + "\n")

    def phase_by_homozygote(self, parents, child, parentIndex) -> list:
        """Attempts phasing using a homozygous parent allele"""
        if not isHomozygote(parents[parentIndex]):
            return []
        phased = ["", ""]
        otherParent = 1 if parentIndex == 0 else 0
        gt = parents[parentIndex][0]
        if (child[0] == gt) and (child[1] in parents[otherParent]):
            phased[otherParent] = child[1]
        elif (child[1] == gt) and (child[0] in parents[otherParent]):
            phased[otherParent] = child[0]
        else:
            return []
        phased[parentIndex] = gt
        return phased

    def parse_phased_vcf(self) -> dict:
        """Main parsing function to evaluate phasing from VCF against truth"""
        self.output = {'sites': len(self.truth), 'errors': 0, 'errors_snp': 0, 'missing': 0, 'missing_snp': 0}
        self.errors = {"TYPE": ["HT", "POS1", "POS2", "ALLLES1", "ALLELES2", "TRUTH1", "TRUTH2"]}
        haplotypes = {}
        prev_sample = {'ht': None}

        with pysam.VariantFile(self.vcf_path) as reader:
            for record in reader:
                if (record.qual == 0):
                    continue
                if (len(record.samples) != 1):
                    self.log("ERROR: Too many samples in file", True)
                    sys.exit(1)
                if (not record.pos in self.truth) or (not record.samples[0].phased):
                    continue

                phased = record.samples[0].alleles
                if (isHomozygote(phased)):
                    continue
                haplotype = 1 if self.ignore_haplotypes else record.samples[0]['PS']
                this_sample = {'ht': haplotype, 'pos': record.pos, 'alleles': fill_dels(phased, record.ref)}

                if (prev_sample['ht'] == None) or (prev_sample['ht'] != this_sample['ht']):
                    prev_sample = this_sample
                    continue

                self.truth[record.pos].append('found')
                if not haplotype in haplotypes:
                    haplotypes[haplotype] = [0, 0, record.pos, 0]
                snp = isSNP(phased)
                haplotypes[haplotype][1 if snp else 0] += 1
                haplotypes[haplotype][3] = record.pos

                success = ((compare_alleles(this_sample['alleles'], self.truth[this_sample['pos']]) and
                            compare_alleles(prev_sample['alleles'], self.truth[prev_sample['pos']])) or
                           (compare_alleles_inverted(this_sample['alleles'], self.truth[this_sample['pos']]) and
                            compare_alleles_inverted(prev_sample['alleles'], self.truth[prev_sample['pos']])))

                if not self.method_name in self.detailed:
                    self.detailed[self.method_name] = {}
                if not self.sample_name in self.detailed[self.method_name]:
                    self.detailed[self.method_name][self.sample_name] = {}

                detail = 1 if success else 0
                if ((len(this_sample['alleles'][0]) != 1) or
                    (len(this_sample['alleles'][1]) != 1) or
                    (len(prev_sample['alleles'][0]) != 1) or
                        (len(prev_sample['alleles'][1]) != 1)):
                    detail += 0.5
                self.detailed[self.method_name][self.sample_name][record.pos] = detail

                if not success:
                    self.output['errors'] += 1
                    if isSNP(phased):
                        self.output['errors_snp'] += 1
                    self.errors[this_sample['pos']] = [
                        str(prev_sample['ht']),
                        str(prev_sample['pos']),
                        str(this_sample['pos']),
                        alleles_to_str(prev_sample['alleles']),
                        alleles_to_str(this_sample['alleles']),
                        alleles_to_str(self.truth[prev_sample['pos']]),
                        alleles_to_str(self.truth[this_sample['pos']])
                    ]
                prev_sample = this_sample

        haplotypes_single = 0
        haplotypes_sizes = []
        for i in haplotypes:
            haplotypes_sizes.append(1 + haplotypes[i][3] - haplotypes[i][2])
            if haplotypes[i][0] + haplotypes[i][1] == 1:
                haplotypes_single += 1

        for pos in self.truth:
            if (len(self.truth[pos]) == 3):
                continue
            self.output['missing'] += 1
            if isSNP(self.truth[pos]):
                self.output['missing_snp'] += 1
            self.errors[pos] = ['?',
                                str(pos),
                                '',
                                '',
                                '',
                                alleles_to_str(self.truth[pos]),
                                '']

        pfull = round(100 * (self.output['sites']-(self.output['errors']+self.output['missing'])) / self.output['sites'], 1)
        psnp = round(100 * (self.output['sites']-(self.output['errors_snp']+self.output['missing_snp'])) / self.output['sites'], 1)
        self.output['summary'] = f"{self.sample_name} / {self.method_name}: %f={pfull}, %s={str(psnp)}"
        self.output['pfull'] = pfull
        self.output['psnp'] = psnp

        # Output summary line
        headerStr = "\t".join([
            "sample", "method", "sites", "sites_snp", "errors", "errors_snp",
            "missing", "missing_snp", "success %", "success % snp",
            "haplotypes", "haplotypes_single", "haplotypes_size_med"
        ])
        summaryStr = "\t".join([
            self.sample_name,
            self.method_name,
            str(self.output['sites']),
            str(self.truth_snps),
            str(self.output['errors']),
            str(self.output['errors_snp']),
            str(self.output['missing']),
            str(self.output['missing_snp']),
            str(pfull),
            str(psnp),
            str(len(haplotypes)),
            str(haplotypes_single),
            "0" if len(haplotypes_sizes) == 0 else str(statistics.median(haplotypes_sizes))
        ])
        self.log(headerStr, True)
        self.log(summaryStr, True)

        fn = f"summary-{self.method_name}.txt"
        if not fn in self.summaryFiles:
            self.summaryFiles[fn] = open(fn, "a")
            self.summaryFiles[fn].write(headerStr + "\n")
        self.summaryFiles[fn].write(summaryStr + "\n")

    def load_truth(self) -> None:
        """Loads the truth dataset (pickled dictionary of allele lists)"""
        with open(self.truth_path, "rb") as outfile:
            self.truth = pickle.load(outfile)
        self.truth_snps = sum(1 for i in self.truth if len(self.truth[i][0]) == 1 and len(self.truth[i][1]) == 1)

    def save_output(self) -> None:
        """Saves summary and error details to output file"""
        if self.output_path == "":
            return
        with open(self.output_path, "wt") as outfile:
            for key in self.output:
                outfile.write(f"INFO\t{key}\t{self.output[key]}\n")
            for key in self.errors:
                outfile.write(f"ERROR\t{'\t'.join(self.errors[key])}\n")

    def cleanup(self) -> None:
        """Resets objects and forces garbage collection"""
        self.output = None
        self.errors = None
        gc.collect()

    def goSingle(self):
        """Process a single sample based on configured paths"""
        if self.output_path != "":
            self.logFile = open(self.output_path[0:-4]+".log", "w")
        else:
            self.logFile = None

        self.load_truth()
        self.parse_phased_vcf()
        self.save_output()

        if self.logFile != None:
            self.logFile.close()
        self.cleanup()

    def go(self):
        """Main execution method � supports single and multiple sample list processing"""
        self.detailed = {}
        if self.list_path == "":
            self.goSingle()
        else:
            reader = csv.reader(open(self.list_path), delimiter='\t')
            for row in reader:
                self.method_name = row[0]
                self.sample_name = row[1]
                self.vcf_path = row[2]
                self.truth_path = row[3]
                self.ignore_haplotypes = row[5] == '-ih'
                print(self.sample_name, end="\r")
                self.goSingle()

        for fn in self.summaryFiles:
            self.summaryFiles[fn].close()

        for method_name in self.detailed:
            with open(method_name+".pkl", "wb") as myfile:
                pickle.dump(self.detailed[method_name], myfile)

# ----------------- Argparse utility functions -----------------


def arg_file_in(x):
    """For argparse - validates that a file exists"""
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError(f"File '{x}' doesn't exist")
    return x


def arg_file_out(x):
    """For argparse - checks that the file path is writable and does not exist"""
    if os.path.exists(x):
        raise argparse.ArgumentTypeError(f"{x} Already exists")
    try:
        f = open(x, "w")
        f.close()
        os.remove(x)
        return x
    except Exception as error:
        raise argparse.ArgumentTypeError(f"File '{x}' isn't writable: {str(error)}")

# ----------------- Main script entry point -----------------


def main():
    """Parse command-line arguments and run the main evaluation logic"""
    parser = argparse.ArgumentParser(description=APP_NAME + " - parents data extractor. Version "+APP_VERSION+" ("+APP_DATE+")",
                                     epilog="For more info visit "+APP_URL)
    parser.add_argument("-s", "--sample", dest="sample_name", help="Sample name")
    parser.add_argument("-m", "--method", dest="method_name", help="Method name")
    parser.add_argument("-v", "--vcf", dest="vcf_path", type=arg_file_in, help="Child VCF file")
    parser.add_argument("-t", "--truth", dest="truth_path", type=arg_file_in, help="Truth VCF file")
    parser.add_argument("-ih", "--ignore-haplotypes", dest="ignore_haplotypes", default=False, action='store_true', help="Ignore haplotypes (assume 1 big haplotype)")
    parser.add_argument("-o", "--output", dest="output_path", default="", help="Output file")
    parser.add_argument("-l", "--list", dest="list_path", default="", help="List file")

    args = vars(parser.parse_known_args()[0]).items()
    obj = AnalyzePhased()
    for param, value in args:
        setattr(obj, param, value)
    obj.go()


if __name__ == "__main__":
    main()
