"""
Reference free genotype demultiplexing on pooled scRNA-seq
"""

import vcf  # https://pyvcf.readthedocs.io/en/latest/INTRO.html
import numpy as np
import pysam as ps  # http://pysam.readthedocs.io/en/latest/api.html#sam-bam-files
import pandas as pd
import subprocess
import random
import math
import datetime

# TODO: Avoid hard k-means, calculate genotype based on weighted contribution for likelihood of cell coming from genotype
# TODO: Remove the h

BASE_RR_PROB = 0.4995
BASE_RA_PROB = 0.001
BASE_AA_PROB = 0.4995
MIN_POST_PROB = 1e-6
PRECISION = 1e-16



class SNV_data:
    """
    Stores data on each SNV
    """

    def __init__(self, chrom, pos, ref, alt, loglik):
        """
        Parameters:
            chrom (int): chromosome number
            pos (int): position on chromosome
            ref (str): reference base
            alt (str): alternate base
            loglik (tuple(float)): log likelihood of genotypes normalised to most likely allele (RR,RA,AA)
        """
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALT = alt
        self.loglik = loglik
        self.RR, self.RA, self.AA = self._calculate_genotype(loglik)


    def _calculate_genotype(self, loglik):
        """
        Calculates genotype probabilities from normalised log likelihoods in vcf file(by likely allele)

        Example:
            if allele RR is the most likely then:
            likRR = P(RR)/P(RR)
            likRA = P(RA)/P(RR)
            likAA = P(AA)/P(RR)

            [P(RR) + P(RA) + P(AA)]/P(RR) = (likRR + likRA + likAA)

            P(RR) = 1/(likRR + likRA + likAA)
            P(RA) = likRA * P(RR)
            P(AA) = likRA * P(RR)

        Returns:
            RR, RA, AA (tuple(floats)): genotype probabilities
        """
        llikRR, llikRA, llikAA = loglik

        likRR = 10**llikRR
        likRA = 10**llikRA
        likAA = 10**llikAA

        # get probability of the most likely allele that everything was normalised to
        pmax = 1/(likRR + likRA + likAA)

        # Calcuate allele probabilities
        RR = likRR * pmax
        RA = likRA * pmax
        AA = likAA * pmax

        return RR, RA, AA


    def get_all_SNV_pos(all_SNVs):
        """
        Return ordered list of positions of SNVs as chr:pos
        Parameters:
            all_SNVs (list(SNV_data obj)): list of SNV_data objects
        Returns:
            sorted list of SNV unique positions as chr:pos
        """
        all_POS = []
        for entry in all_SNVs:
            pos = str(entry.CHROM) + ':' + str(entry.POS)
            if pos not in all_POS:
                all_POS.append(pos)
        return all_POS


class model_genotype:
    """
    Calculates model genotypes from base call frequencies
    Given an overall error probability 0.01
    For a observed single allele (A,C,G,T)
    P(error) = 0.01/4
    P(hom) = 1 - (3 * P(error))
    P(het) = 0.5 - P(error)

    Example likelihood of hom allele AA given base frequencies F(X):
    L(AA | F, P(error)) = P(hom)^F(A) * P(error)^(F(C)+(F(G))+F(T))

    Example likelihood of het allele CT:
    L(CT | F, P(error)) = P(het)^(F(C)+F(T)) * P(error)^(F(A)+F(G))
    """

    def __init__(self, all_SNVs, base_calls_mtx, barcodes, num=2,
                 model_genotypes=[], assigned=None):
        """
        Parameters:
             all_SNVs: list[SNV_data objects]
             base_calls_mtx: SNV-barcode matrix containing lists of base calls
             num(int): number of individual genotypes
             barcodes: list of cell barcodes
             model_genotypes(list(DataFrame): list of num model genotypes represented in snv-probdistrib DataFrame as <RR,RA,AA>
             assigned(list): lists of cell/barcode assigned to each genotype model
        """
        self.all_SNVs = all_SNVs
        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.num = num
        self.barcodes = barcodes
        self.assigned = assigned
        self.model_genotypes = model_genotypes

        self.assign_cells_llhood = []
        for n in range(num):
            self.assign_cells_llhood.append([])
            self.assign_cells_llhood[n] = pd.DataFrame(
                np.ones((len(self.all_SNVs), 3)),
                index=SNV_data.get_all_SNV_pos(self.all_SNVs),
                columns=['RR', 'RA', 'AA'])

        self.model_llhood = []
        for n in range(num):
            self.model_llhood.append([])
            self.model_llhood[n] = pd.DataFrame(
                np.ones((len(self.all_SNVs), 3)),
                index=SNV_data.get_all_SNV_pos(self.all_SNVs),
                columns=['RR', 'RA', 'AA'])

        if self.model_genotypes == []:
            for n in range(num):
                self.model_genotypes.append([])
                self.model_genotypes[n] = pd.DataFrame(
                    np.zeros((len(self.all_SNVs), 3)),
                    index=SNV_data.get_all_SNV_pos(self.all_SNVs),
                    columns=['RR', 'RA', 'AA'])

    def initialise_cell_assignments(self):
        """ Random initialisation of cell assignments"""
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])
        for barcode in self.barcodes:
            n = random.randint(0, self.num - 1)
            self.assigned[n].append(barcode)
        return self.assigned

    def initialise_model_genotypes(self):
        """ Initialises model from dirichlet distribution based on calculated
        average genotype in vcf file """
        for snv in self.all_SNVs:
            pos = "{}:{}".format(snv.CHROM, snv.POS)
            gt = np.random.dirichlet((snv.RR, snv.RA, snv.AA),
                                     self.num).transpose()
            for n in range(self.num):
                self.model_genotypes[n].loc[pos] = (
                    gt[0][n], gt[1][n], gt[2][n])

    def calculate_model_genotypes(self):
        """
        Calculates a probability distribution <RR, RA, AA> for each SNV position based on counts of bases
        """
        # TODO: How to handle when no snv calls at a position for that model
        coverage = [0, 0]
        no_coverage = [0, 0]
        for m in range(self.num):
            self.model_llhood[m].loc[:] = 1
        for snv in self.all_SNVs:  # for each snv position
            if len(snv.REF) == 1 and len(snv.ALT) == 1:  # skip indels
                pos = "{}:{}".format(snv.CHROM, snv.POS)
                for n in range(
                        self.num):  # for the number of genotypes being modelled
                    cluster = self.assigned[n]  # list of barcodes

                    # get snv base calls for cells assigned to model, returns pandas.Series
                    ref_count = sum(self.ref_bc_mtx.loc[pos, cluster])
                    alt_count = sum(self.alt_bc_mtx.loc[pos, cluster])

                    if ref_count > 0 or alt_count > 0:  # if cells in genotype carry reads at this snv
                        coverage[n] += 1
                        llhood = self.calc_model_likelihood(
                            ref_count, alt_count, pos)
                        self.model_llhood[n].loc[pos] = llhood

                        P_RR_given_D, P_RA_given_D, P_AA_given_D = \
                            self.calc_model_posterior(llhood)
                        self.model_genotypes[n].loc[pos] = (
                            P_RR_given_D, P_RA_given_D, P_AA_given_D)
                    else:
                        no_coverage[n] += 1
                        self.model_genotypes[n].loc[pos] = \
                            (BASE_RR_PROB, BASE_RA_PROB, BASE_AA_PROB)
        print("Coverage for model 1 = {}, no coverage = {}".format(
            coverage[0], no_coverage[0]))
        print("Coverage for model 2 = {}, no coverage = {}".format(
            coverage[1], no_coverage[1]))

    def calc_model_likelihood(self, ref_count, alt_count, pos):
        """
        Calculates log likelihood (P(D|G))
        Parameters:
            ref_count: count of reference base calls
            alt_count: count of alternate base calls

        https://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
        """
        e_rate = 0.01

        log_precision_on_n = math.log2(PRECISION / 3)

        # Sum log, product normal

        # reference count
        lRR = math.log2(1 - e_rate) * ref_count
        lRA = math.log2(0.5 - (e_rate / 3)) * ref_count
        lAA = math.log2(e_rate / 3) * ref_count
        # alternate count
        lRR += math.log2(e_rate / 3) * alt_count
        lRA += math.log2(0.5 - (e_rate / 3)) * alt_count
        lAA += math.log2(1 - e_rate) * alt_count

        maximum = max(lRR, lRA, lAA)
        RR = self.check_log_value(lRR, maximum)
        RA = self.check_log_value(lRA, maximum)
        AA = self.check_log_value(lAA, maximum)

        denom = RR + RA + AA
        P_D_given_RR = RR / denom
        P_D_given_RA = RA / denom
        P_D_given_AA = AA / denom

        return (P_D_given_RR, P_D_given_RA, P_D_given_AA)

    def check_log_value(self, prob, maximum):
        log_precision_on_n = math.log2(PRECISION / 3)
        if (prob - maximum) < log_precision_on_n:
            return 0
        else:
            return math.pow(2, prob - maximum)

    def calc_model_posterior(self, llhood):
        """
        Calculates posterior probability (P(G|D))
        Parameter:
            Log likelihood probabilities (P(D|G))
        """
        P_D_given_RR, P_D_given_RA, P_D_given_AA = llhood

        h = 0.001  # uniform heterozygosity rate
        P_hom = (1 - h) / 4
        P_het = h / 6

        z = (P_D_given_RR * P_hom) + \
            (P_D_given_RA * P_het) + \
            (P_D_given_AA * P_hom)

        P_RR_given_D = (P_D_given_RR * P_hom) / z
        P_RA_given_D = (P_D_given_RA * P_het) / z
        P_AA_given_D = (P_D_given_AA * P_hom) / z

        for prob in [P_RR_given_D, P_RA_given_D, P_AA_given_D]:
            if prob < MIN_POST_PROB:
                prob = MIN_POST_PROB

        return (P_RR_given_D, P_RA_given_D, P_AA_given_D)

    def assign_cells(self):
        """
        Assign cells/barcodes to most likely genotype model
        Calculate likelihood cell/barcode came from individual i (i.e. model genotype 1...n):
        For each cell/barcode c
        Product over all variants v
        Sum over all genotypes <RR,RA,AA>
        Product over all reads in cell c for variant v
        Pr(base_called | G)*P(G)
        where P(G) = set(RR,RA,AA)

        Data required:
        cell assignments : self.assigned = [[barcodes],[barcodes]]
        genotype model : self.model_genotypes
        snv pos : self.model_genotypes[0].index
        self.bc_mtx : SNV-barcode matrix containing lists of base calls

        Issue: Some barcodes in matrix have no base calls as they do not cover snp region
        Due to base_calls_mtx being instantiated from list of all barcodes in bam file (not just those covering snv)
        """
        self.old_assignment = self.assigned
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])

        e = 0.01
        self.old_llhood = self.assign_cells_llhood
        for m in range(self.num):
            self.assign_cells_llhood[m].loc[:] = 1

        for barcode in self.barcodes:  # for each cell

            cell_llhood = []  # store likelihood values of cell/barcode c for model 1 to n
            indices = []  # index of non-zero values in column
            for idx in self.ref_bc_mtx.loc[:, barcode].nonzero():
                indices.extend(idx)
            for idx in self.alt_bc_mtx.loc[:, barcode].nonzero():
                indices.extend(idx)
            for n in range(self.num):  # for each model

                all_var_llhood = 1  # store product over all variants for model n

                for i in indices:  # for each variant with non-zero values in this cell
                    snv = self.all_SNVs[i]
                    pos = "{}:{}".format(snv.CHROM, snv.POS)
                    RR, RA, AA = self.model_genotypes[n].loc[pos]

                    ref_count = self.ref_bc_mtx.loc[pos, barcode]
                    alt_count = self.alt_bc_mtx.loc[pos, barcode]

                    # sum log over all base counts (or product normal)
                    # calculate reference
                    P_D_given_RR = ((1 - e) * RR) ** ref_count
                    P_D_given_RA = ((0.5 - e / 3) * RA) ** ref_count
                    P_D_given_AA = ((e / 3) * AA) ** ref_count

                    # calculate alternate
                    P_D_given_RR *= ((e / 3) * RR) ** alt_count
                    P_D_given_RR *= ((0.5 - e / 3) * RA) ** alt_count
                    P_D_given_AA *= ((1 - e) * AA) ** alt_count

                    # for each snv sum log likelihood for genotype
                    self.assign_cells_llhood[n].loc[pos, 'RR'] *= P_D_given_RR
                    self.assign_cells_llhood[n].loc[pos, 'RA'] *= P_D_given_RA
                    self.assign_cells_llhood[n].loc[pos, 'AA'] *= P_D_given_AA

                    # Sum over all genotypes at variant position
                    P_D_given_G = P_D_given_RR + \
                                  P_D_given_RA + \
                                  P_D_given_AA
                    # Product over all variant positions
                    all_var_llhood *= P_D_given_G

                cell_llhood.append(all_var_llhood)
            # if likelihoood is equal in cases such as:
            # 1. barcode has no coverage over any snv region
            # 2. barcode covers only snv regions where <RR, RA, AA> is same for all model genotypes
            if all(i == cell_llhood[0] for i in
                   cell_llhood):
                pass
            # for n in range(self.num):
            #                     self.assigned[n].append(barcode)
            else:
                n = cell_llhood.index(max(cell_llhood))
                self.assigned[n].append(barcode)
        for n in range(self.num):
            self.assigned[n] = sorted(self.assigned[n])

    def compare_cell_assignments(self):
        dif = []
        for n in range(self.num):
            old = self.old_assignment[n]
            new = self.assigned[n]
            count = 0
            for item in new:
                if item in old:
                    count += 1
            print('{}% of barcodes in new also in old'.format(
                count * 100 / len(new)))

            sm = jaccard_similarity(old, new)
            dif.append(sm)
        return dif


def get_all_barcodes(bc_file):
    """
    Load barcodes from text file
    Parameters:
        bc_file: text file of line separated cell barcodes
    Returns:
        list of cell barcodes
    """
    file = open(bc_file, 'r')
    barcodes = []
    for line in file:
        barcodes.append(line.strip())
    return barcodes


def build_base_calls_matrix(sam_filename, all_SNVs, all_POS, barcodes):
    """
    Build pandas DataFrame
    Parameters:
        sam_filename(str): Path to sam file (0-based positions)
        all_SNVs: list of SNV_data objects
        all_POS(list('chr:pos)): snv positions (1-based positions from vcf file)
        barcodes(list): cell barcodes
    """
    #TODO: fix up index names (in case of duplicate position)

    in_sam = ps.AlignmentFile(sam_filename, 'rb')
    ref_base_calls_mtx = pd.DataFrame(np.zeros((len(all_POS), len(barcodes))),
                                  index=all_POS, columns=barcodes)
    alt_base_calls_mtx = pd.DataFrame(np.zeros((len(all_POS), len(barcodes))),
                                      index=all_POS, columns=barcodes)

    print('Matrix size: Num Pos:', len(all_POS), 'Num barcodes:',len(barcodes))

    all_POS = []
    for entry in all_SNVs:
        pos = str(entry.CHROM) + ':' + str(entry.POS)
        if pos not in all_POS:
            all_POS.append(pos)
    print("Searching positions: 0 - 1000000")
    m = 10000
    for snv in all_SNVs:
        position = str(snv.CHROM) + ':' + str(snv.POS)
        if snv.POS/m >= 1:
            print("Searching positions: {} - {}".format(m,m+10000))
            m += 10000
        for pileupcolumn in in_sam.pileup(snv.CHROM, snv.POS-1, snv.POS, truncate=True):  # pysam uses 0 based positions
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    barcode = get_read_barcodes(sam_filename, snv.CHROM, snv.POS,
                                                pileupread.alignment.query_name)
                    if barcode is not None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if base == snv.REF:
                            ref_base_calls_mtx.loc[position, barcode] += 1
                        if base == snv.ALT:
                            alt_base_calls_mtx.loc[position, barcode] += 1

    return (ref_base_calls_mtx, alt_base_calls_mtx)



def get_read_barcodes(sam_file, chr, pos, readname):
    """
    Get the cell barcode (CB tag information) from sam file
    Pysam does not hold this information therefore requires use of the
    command line functions.
    Parameters:
         sam_file: Absolute location of sam/bam file
         chr: chromosome name in sam
         pos: position of snv
         readname: unique read name in sam file
    """
    loc = "{0}:{1}-{2}".format(chr,int(pos)-1, int(pos)+1)
    cmd = ['samtools', 'view', sam_file,loc]
    result = subprocess.check_output(cmd)
    for entry in str(result).lstrip('b\'').split('\\n'):
        if (entry.split("\\t")[0]) == readname:
            try:
                return (entry.split('CB:Z:')[1].split('\\tUR:Z:')[0])
            except IndexError:
                print('No CB tag for', readname) # can return error if no CB tag attached to read in bam file
        else:
            continue


def run_model(all_SNVs, base_calls_mtx, barcodes, num_models):
    model = model_genotype(all_SNVs, base_calls_mtx, barcodes, num_models, model_genotypes=[], assigned=None)
    dif = [0 for _ in range(num_models)]
    print("\ninitialising model genotypes", datetime.datetime.now().time())
    model.initialise_model_genotypes()
    print("initial cell assignments", datetime.datetime.now().time())
    model.assign_cells()
    length=[]
    for assigned in model.assigned:
        length.append(len(assigned))
    print(length)
    print("Commencing E-M")
    while (any(i < 0.80 for i in dif) == True):
        print("calculating model ", datetime.datetime.now().time())
        model.calculate_model_genotypes()
        print("assigning cells", datetime.datetime.now().time())
        model.assign_cells()
        old_dif = dif
        dif = model.compare_cell_assignments()
        length=[]
        for assigned in model.assigned:
            length.append(len(assigned))
        print(length)
        print("difference:", dif)
        if dif == old_dif and sum(dif) > 0:
            break
    print("final difference:", dif)
    return model


def jaccard_similarity(x,y):
    intersect = len(set.intersection(*[set(x), set(y)]))
    union = len(set.union(*[set(x), set(y)]))
    return intersect/float(union)




def main():
    path = '/Users/Caitlin/Documents/Bioinformatics/Summer_Project_2017/analysis/'
    out_dir = 'data/'
    outfile_A = 'pbmc50mix/model_A{}_chr1.txt'
    outfile_B = 'pbmc50mix/model_B{}_chr1.txt'
    out_csv_ref = 'pbmc50mix/ref_chr1.csv'
    out_csv_alt = 'pbmc50mix/alt_chr1.csv'


    # # 50:50 10k
    # file_v = "pbmc50mix/pbmc50_10k.vcf"
    # file_s = "pbmc50mix/pbmc_50_50_10k.bam"
    # file_bc = "pbmc50mix/bc_sorted_10k.txt"

    # # 50:50 files 100k
    # file_v = "pbmc50mix/pbmc50_100k.vcf"
    # file_s = "pbmc50mix/pbmc50_100k_bc_only.bam"
    # file_bc = "pbmc50mix/barcodes50_100k_sorted.txt"

    # # 50:50 files 500k
    # file_v = "pbmc50mix/pbmc50_500k.vcf"
    # file_s = "pbmc50mix/pbmc50_500k_bc_only.bam"
    # file_bc = "pbmc50mix/barcodes50_500k_sorted.txt"

    # 50:50 chr1
    file_v = "pbmc50mix/pbmc_chr1.vcf"
    file_s = "pbmc50mix/pbmc_chr1_cleaned.bam"
    file_bc = "pbmc50mix/pbmc_chr1_cleaned_bc.txt"

    # # Mixed donor files
    # file_v = "pbmc_A-C_mix/pbmc_30k.vcf"
    # file_s = "pbmc_A-C_mix/pbmc_30k_bc_only.bam"
    # file_bc = "pbmc_A-C_mix/bc_30k_sorted.txt"

    vcf_filename = str(path + out_dir + file_v)
    sam_filename = str(path + out_dir + file_s)
    bc_filename = str(path + out_dir + file_bc)

    in_vcf = vcf.Reader(open(vcf_filename, 'r'))

    vcf_records = []
    for record in in_vcf:
        vcf_records.append(record)

    all_SNVs = []  # list of SNV_data objects
    for record in vcf_records:
        all_SNVs.append(
            SNV_data(record.CHROM, record.POS, record.REF, record.ALT[0],
                     record.samples[0]['GL'][0:3]))

    print("Last entry in vcf record: ", all_SNVs[-1].CHROM, ":",
          all_SNVs[-1].POS, " REF:", all_SNVs[-1].REF, " ALT:",
          all_SNVs[-1].ALT, " <", all_SNVs[-1].RR, ", ", all_SNVs[-1].RA, ", ",
          all_SNVs[-1].AA, ">", sep='')

    print("Starting data collection", datetime.datetime.now().time())
    all_POS = SNV_data.get_all_SNV_pos(all_SNVs)
    barcodes = get_all_barcodes(bc_filename)
    base_calls_mtx = build_base_calls_matrix(sam_filename, all_SNVs, all_POS,
                                             barcodes)
    base_calls_mtx[0].to_csv('{}{}{}'.format(path, out_dir, out_csv_ref))
    base_calls_mtx[1].to_csv('{}{}{}'.format(path, out_dir, out_csv_alt))
    print("Base call matrix finished", datetime.datetime.now().time())


    # show some data entered okay by entering location of last entry in vcf
    # majority of bases should match ref and alt in vcf above
    location = "{}:{}".format(all_SNVs[-1].CHROM, all_SNVs[-1].POS)
    print("RR,RA,AA : <", all_SNVs[-1].RR, all_SNVs[-1].RA, all_SNVs[-1].AA,
          ">")
    alt = 0
    ref = 0

    for base in base_calls_mtx[0].loc[location]:
        ref += base
    for base in base_calls_mtx[1].loc[location]:
        alt += base
    print("Reference count:", ref, "; Alternate count:", alt)


    num_runs = 3
    num_models = 2

    all_models = []
    for run in range(num_runs):
        model = run_model(all_SNVs, base_calls_mtx, barcodes, num_models)
        all_models.append(model)
        print("Finished model {} at {}".format(run + 1, datetime.datetime.now().time()))

    print("\nFinished all models... Comparing models...")

    # Separate models output by multiple runs based on similarity to each other
    for e, model in enumerate(all_models):
        if e == 0:  # first iteration
            model_A = [model.assigned[0]]
            model_B = [model.assigned[1]]
        else:
            # compare current two to first model entered in A and B
            for n in range(2):
                smA = jaccard_similarity(model_A[0], model.assigned[n])
                smB = jaccard_similarity(model_B[0], model.assigned[n])
                print("model_A to model{}[{}] :".format(e + 1, n), smA)
                print("model_B to model{}[{}] :".format(e + 1, n), smB)

                match = (smA, smB).index(max(smA, smB))
                if match == 0:
                    model_A.append(model.assigned[n])
                elif match == 1:
                    model_B.append(model.assigned[n])

    for n in range(len(model_A)):
        for m in range(n + 1, len(model_A)):
            sm = jaccard_similarity(model_A[n], model_A[m])
            print("model_A-{}, model_A-{} :".format(n + 1, m + 1), sm)

    for n in range(len(model_B)):
        for m in range(n + 1, len(model_B)):
            sm = jaccard_similarity(model_B[n], model_B[m])
            print("model_B-{}, model_B-{} :".format(n + 1, m + 1), sm)


    for m, model in enumerate(model_A):
        filename = '{}{}{}'.format(path, out_dir, outfile_A).format(m + 1)
        file = open(filename, "w")
        for barcode in model:
            file.write(barcode + '\n')
        file.close()

    for m, model in enumerate(model_B):
        filename = '{}{}{}'.format(path, out_dir, outfile_B).format(m + 1)
        file = open(filename, "w")
        for barcode in model:
            file.write(barcode + '\n')
        file.close()

if __name__ == "__main__":
    main()
