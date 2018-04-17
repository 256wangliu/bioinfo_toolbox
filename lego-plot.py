#!/usr/bin/env python
# Author: Junko Tsuji <jtsuji@broadinstitute.org>

import sys
import os.path
import fileinput
import warnings
from argparse import ArgumentParser

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# MutSig 2CV precomputed coverage data based on CGA's lego plotter
PRECOMP_EXOME = {
    'AAA':516645, 'AAC':336113, 'AAG':545347, 'AAT':365625,
    'CAA':462172, 'CAC':488457, 'CAG':797766, 'CAT':489731,
    'GAA':549735, 'GAC':357085, 'GAG':540864, 'GAT':397040,
    'TAA':230606, 'TAC':269162, 'TAG':226201, 'TAT':252346,
    'ACA':469526, 'ACC':429764, 'ACG':160802, 'ACT':389408,
    'CCA':664422, 'CCC':477669, 'CCG':235217, 'CCT':616696,
    'GCA':509693, 'GCC':512648, 'GCG':186053, 'GCT':541099,
    'TCA':537050, 'TCC':568802, 'TCG':169936, 'TCT':589050,
    'AGA':577447, 'AGC':538629, 'AGG':601391, 'AGT':389620,
    'CGA':168530, 'CGC':186767, 'CGG':234416, 'CGT':162102,
    'GGA':556167, 'GGC':509431, 'GGG':469219, 'GGT':426880,
    'TGA':536918, 'TGC':514744, 'TGG':661776, 'TGT':473758,
    'ATA':252140, 'ATC':397574, 'ATG':486994, 'ATT':367907,
    'CTA':226652, 'CTC':552864, 'CTG':803001, 'CTT':556516,
    'GTA':267878, 'GTC':359998, 'GTG':486805, 'GTT':338852,
    'TTA':231118, 'TTC':560177, 'TTG':467349, 'TTT':524231
}
PRECOMP_GENOME = {
    'AAA':94564228, 'AAC':36158942, 'AAG':49494329, 'AAT':61919448,
    'CAA':46789183, 'CAC':35927089, 'CAG':48482195, 'CAT':45521663,
    'GAA':48882158, 'GAC':22927567, 'GAG':40130275, 'GAT':33037787,
    'TAA':51852385, 'TAC':28182087, 'TAG':32110875, 'TAT':50214797,
    'ACA':49545058, 'ACC':27998332, 'ACG': 5899505, 'ACT':39750998,
    'CCA':44139802, 'CCC':30001011, 'CCG': 5654254, 'CCT':42486275,
    'GCA':34775000, 'GCC':27270757, 'GCG': 4806689, 'GCT':33579428,
    'TCA':48349332, 'TCC':37043310, 'TCG': 5098333, 'TCT':54501232,
    'AGA':54403725, 'AGC':33559069, 'AGG':42441923, 'AGT':39792917,
    'CGA': 5088120, 'CGC': 4802745, 'CGG': 5651927, 'CGT': 5914345,
    'GGA':37066788, 'GGC':27263586, 'GGG':30032551, 'GGT':28031174,
    'TGA':48352258, 'TGC':34807746, 'TGG':44225804, 'TGT':49728134,
    'ATA':50164185, 'ATC':33013168, 'ATG':45519675, 'ATT':61986409,
    'CTA':32073940, 'CTC':40143640, 'CTG':48528214, 'CTT':49570546,
    'GTA':28197095, 'GTC':22968750, 'GTG':36009505, 'GTT':36288035,
    'TTA':51937650, 'TTC':48943412, 'TTG':46991718, 'TTT':94906144
}

# mutational context
SNP = ['C>T', 'C>A', 'C>G', 'A>G', 'A>C', 'A>T']
CONTEXT = ['T_G', 'T_A', 'T_C', 'T_T',
           'C_G', 'C_A', 'C_C', 'C_T',
           'A_G', 'A_A', 'A_C', 'A_T',
           'G_G', 'G_A', 'G_C', 'G_T']
MUT_CONTEXT = ['T_G.C>T', 'T_A.C>T', 'T_C.C>T', 'T_T.C>T', 'T_G.C>A', 'T_A.C>A',
               'T_C.C>A', 'T_T.C>A', 'T_G.C>G', 'T_A.C>G', 'T_C.C>G', 'T_T.C>G',
               'C_G.C>T', 'C_A.C>T', 'C_C.C>T', 'C_T.C>T', 'C_G.C>A', 'C_A.C>A',
               'C_C.C>A', 'C_T.C>A', 'C_G.C>G', 'C_A.C>G', 'C_C.C>G', 'C_T.C>G',
               'A_G.C>T', 'A_A.C>T', 'A_C.C>T', 'A_T.C>T', 'A_G.C>A', 'A_A.C>A',
               'A_C.C>A', 'A_T.C>A', 'A_G.C>G', 'A_A.C>G', 'A_C.C>G', 'A_T.C>G',
               'G_G.C>T', 'G_A.C>T', 'G_C.C>T', 'G_T.C>T', 'G_G.C>A', 'G_A.C>A',
               'G_C.C>A', 'G_T.C>A', 'G_G.C>G', 'G_A.C>G', 'G_C.C>G', 'G_T.C>G',
               'T_G.A>G', 'T_A.A>G', 'T_C.A>G', 'T_T.A>G', 'T_G.A>C', 'T_A.A>C',
               'T_C.A>C', 'T_T.A>C', 'T_G.A>T', 'T_A.A>T', 'T_C.A>T', 'T_T.A>T',
               'C_G.A>G', 'C_A.A>G', 'C_C.A>G', 'C_T.A>G', 'C_G.A>C', 'C_A.A>C',
               'C_C.A>C', 'C_T.A>C', 'C_G.A>T', 'C_A.A>T', 'C_C.A>T', 'C_T.A>T',
               'A_G.A>G', 'A_A.A>G', 'A_C.A>G', 'A_T.A>G', 'A_G.A>C', 'A_A.A>C',
               'A_C.A>C', 'A_T.A>C', 'A_G.A>T', 'A_A.A>T', 'A_C.A>T', 'A_T.A>T',
               'G_G.A>G', 'G_A.A>G', 'G_C.A>G', 'G_T.A>G', 'G_G.A>C', 'G_A.A>C',
               'G_C.A>C', 'G_T.A>C', 'G_G.A>T', 'G_A.A>T', 'G_C.A>T', 'G_T.A>T']
MUT_SET = set(MUT_CONTEXT) # to access contents faster than list

# reverse complement
REVCOMP = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}


# function to plot lego plot (made based on Justin's code)
def generate_lego_plot(muts, total_muts, output_prefix, plot_title, af_slice, is_rate_per_Mb=False):
    # filter warnings from matplotlib plotting
    warnings.filterwarnings("ignore")

    fig = plt.figure(figsize=(13,8), dpi=300)
    ax = fig.add_subplot(111, projection="3d")

    # coordinates of bottom boxes
    xpos = []
    ypos = []

    # coordinates of top boxes and widths
    dx = [0.8] * 112
    dy = [0.8] * 112
    dz = muts + [0] * 16

    # box colors
    cols = []
    lego_colors = {"yellow":"#FEFD38", "aqua":"#1CB2B1", "red":"#FC0F1C",
                   "green":"#2BCA2E", "blue":"#0C3BC8", "purple":"#7D4DAC",
                   "white":"#FFFFFF"}

    # axis coordinate
    for j in range(1,9):
        for i in range(1,13):
            xpos.append(i)
            ypos.append(j)
            if i <= 4 and j <= 4:
                cols.append(lego_colors["yellow"])
            elif (4 < i and i <= 8) and j <= 4:
                cols.append(lego_colors["aqua"])
            elif 8 < i and j <= 4:
                cols.append(lego_colors["red"])
            elif i <= 4 and 4 < j:
                cols.append(lego_colors["green"])
            elif (4 < i and i <= 8) and 4 < j:
                cols.append(lego_colors["blue"])
            elif 8 < i and 4 < j:
                cols.append(lego_colors["purple"])

    # context label coordinate
    for j in range(9,13):
        for i in range(5,9):
            xpos.append(i)
            ypos.append(j)
            cols.append(lego_colors["white"])

    ax.view_init(elev=60, azim=55)
    ax.set_zlim(bottom=0, top=max([5,max(dz)]))
    ax.set_xlim(left=1, right=12.5)
    ax.set_ylim(bottom=1, top=12)

    # plot boxes
    for idx in range(len(xpos)):
        ax.bar3d(xpos[idx], ypos[idx], 0, dx[idx], dy[idx], dz[idx],
                 linewidth=0.3, color=cols[idx], zsort="max", edgecolor="k")
        ax.collections[idx].set_sort_zpos(idx * 0.1)
        ax.collections[idx].set_facecolors([cols[idx]]*6)

    # get rid of the panes
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # get rid of the spines
    ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # get rid of the ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # add legend
    C_to_T = plt.Rectangle((0,0), 1, 1, fc=lego_colors["yellow"])
    C_to_A = plt.Rectangle((0,0), 1, 1, fc=lego_colors["aqua"])
    C_to_G = plt.Rectangle((0,0), 1, 1, fc=lego_colors["red"])
    A_to_G = plt.Rectangle((0,0), 1, 1, fc=lego_colors["green"])
    A_to_C = plt.Rectangle((0,0), 1, 1, fc=lego_colors["blue"])
    A_to_T = plt.Rectangle((0,0), 1, 1, fc=lego_colors["purple"])
    ax.legend([C_to_T, C_to_A, C_to_G, A_to_G, A_to_C, A_to_T], SNP, loc=(0.8, 0.6))

    zlabel = "Mutation Count"
    if is_rate_per_Mb:
        zlabel = "Mutation Rate / Mb"
    ax.set_zlabel(zlabel)

    # add context text
    for idx, pos in enumerate(range(96,112)):
        ax.text(xpos[pos]+0.45, ypos[pos]+0.45, 0, CONTEXT[idx],
                size="small", ha="center", va="center", zorder=112)  

    # add title and number of mutation
    title = output_prefix
    if plot_title:
        title = plot_title
    if sum(af_slice) < 1 and not plot_title:
        title += " ({} <= AF < {})".format(*af_slice)
    ax.text2D(0.1, 0.94, title, transform=ax.transAxes)
    ax.text2D(0.1, 0.90, "n = {}".format(total_muts), transform=ax.transAxes)

    fig.savefig(output_prefix+".pdf", bbox_inches="tight")


# function to calculate mutation rate per Mb
def compute_rate_per_Mb(muts, coverage):
    # mutational context without alt base
    context = [mc.split(">")[0] for mc in MUT_CONTEXT]
    weight = [0.0 for mc in MUT_CONTEXT]

    # translate 64 3-base contexts into 96 mutational contexts
    for bases in coverage:
        mc = mutational_context(bases)
        if mc not in context:
            bases = "".join([REVCOMP[b] for b in bases[::-1]])
            mc = mutational_context(bases)

        # integrate the counts
        for i, c in enumerate(context):
            if mc != c:
                continue
            weight[i] += coverage[bases]

    # compute mutation per Mb
    for i in range(len(context)):
        muts[i] = round(muts[i] * 10**6 / max(weight[i],1.0), 4)
    return muts


# function to store sequences in FASTA
def get_fasta_sequence(fin):
    ref_fasta = {}
    name = ""
    seq = ""
    for x in open(fin):
        x = x.rstrip().split()[0]
        if x.startswith(">"):
            if name:
                ref_fasta.setdefault(name, seq)
            name = x[1:]
            seq = ""
            continue
        seq += x
    if name:
        ref_fasta.setdefault(name, seq)
    return ref_fasta


# function to return mutational context string
def mutational_context(bases):
    return "{0}_{2}.{1}".format(*tuple(bases))


# function to read input VCF
def read_vcf(fin, include_all, filtered_only, ref_fasta_file):
    ref_fasta = get_fasta_sequence(ref_fasta_file)

    passed = [0 for mc in MUT_CONTEXT]
    filtered = [0 for mc in MUT_CONTEXT]

    for x in fileinput.input(fin):
        # skip VCF header lines
        if x.startswith("#"):
            continue
        x = x.rstrip("\n").split("\t")

        ref_base = x[3].upper()
        alt_base = x[4].upper()

        is_snp = len(ref_base) == 1 and len(alt_base) == 1
        is_valid_base = ref_base in REVCOMP and alt_base in REVCOMP

        if not is_snp or not is_valid_base:
            continue

        chrom = x[0]
        pos = int(x[1])
        bases = ref_fasta[chrom][pos-2:pos+1]
        judgement = x[6]

        mc = mutational_context(bases) + ">" + alt_base

        if mc not in MUT_SET:
            # convert to reverse compliment
            bases = "".join([REVCOMP[b] for b in bases[::-1]])
            mc = mutational_context(bases) + ">" + REVCOMP[alt_base]
            if mc not in MUT_SET:
                sys.stderr("cannot find context for {}".format(mc))
                continue

        index = MUT_CONTEXT.index(mc)
        if judgement == "PASS":
            passed[index] += 1
        else:
            filtered[index] += 1

    # tally passed and filtered if counting all mutations
    if include_all:
        passed = [p+f for p,f in zip(passed, filtered)]
    elif filtered_only:
        passed = filtered

    return passed


# function to read input MAF
def read_maf(fin, af_slice, include_all, filtered_only):
    # 4 required columns for MAF input
    # i_tumor_f = TCGA style | tumor_f = Oncotator style
    header = ["ref_context", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type", "i_tumor_f", "tumor_f"]

    xcol = []
    passed = [0 for mc in MUT_CONTEXT]
    filtered = [0 for mc in MUT_CONTEXT]
    for x in fileinput.input(fin):
        if x.startswith("#"):
            continue
        x = x.rstrip("\n").split("\t")
        if x[0].startswith("Hugo"):
            # extract required column index
            for n in header:
                try:
                    xcol.append(x.index(n))
                except ValueError:
                    pass
            if len(xcol) < 5:
                raise Exception("MAF header needs to have {}, and either {}".format( \
                                ", ".join(header[:4]), " or ".join(header[4:])))
            continue

        if x[xcol[3]] != "SNP":
            continue

        try:
            tumor_af = float(x[xcol[4]])
        except ValueError:
            tumor_af = float(x[xcol[5]])

        # ignore lines with allele fraction not in the specified range
        if tumor_af < af_slice[0] or af_slice[1] <= tumor_af:
            continue

        ref_context = x[xcol[0]].upper()
        ref_base = x[xcol[1]].upper()
        alt_base = x[xcol[2]].upper()

        is_snp = len(ref_base) == 1 and len(alt_base) == 1
        is_valid_base = ref_base in REVCOMP and alt_base in REVCOMP

        if not is_snp or not is_valid_base:
            continue

        # obtain mutational context
        center = len(ref_context)//2
        base_f = ref_context[center-1]
        base_r = ref_context[center+1]

        bases = "".join([base_f, ref_base, base_r])
        mc = mutational_context(bases) + ">" + alt_base
        if mc not in MUT_SET:
            # convert to reverse compliment
            bases = "".join([REVCOMP[b] for b in bases[::-1]])
            mc = mutational_context(bases) + ">" + REVCOMP[alt_base]
            if mc not in MUT_SET:
                sys.stderr("cannot find context for {}".format(mc))
                continue
        index = MUT_CONTEXT.index(mc)
        if 'FAIL' in x:
            filtered[index] += 1
        else:
            passed[index] += 1

    # tally passed and filtered if counting all mutations
    if include_all:
        passed = [p+f for p,f in zip(passed, filtered)]
    elif filtered_only:
        passed = filtered

    return passed


# function to return context weighting factors
def get_context_coverage(user_input, mutsig_genome, mutsig_exome):
    coverage = {}

    if user_input:
        if not os.path.exists(user_input):
            raise Exception("user coverage input {} does not exist".format(user_input))
        # when specified with other coverage options, prioritize user input
        if mutsig_genome or mutsig_exome:
            sys.stderr("prioritizing user coverage input, turning off MutSig coverage")
            args.mutsig_exome = False
            args.mutsig_genome = False
        # read user coverage input
        for x in open(user_input):
            x = x.rstrip().split()
            if len(x) != 2:
                continue
            if x[0] not in PRECOMP_GENOME:
                continue
            coverage.setdefault(x[0], int(x[1]))

        if len(coverage) != 64:
            raise Exception("needs 64 3-base contexts in user coverage input")

    if mutsig_genome and mutsig_exome:
        raise Exception("specify either '--mutsig-genome' or '--mutsig-exome'")

    elif mutsig_genome:
        coverage = PRECOMP_GENOME

    elif mutsig_exome:
        coverage = PRECOMP_EXOME

    return coverage


def main(args):
    # check input variant file path
    if not os.path.exists(args.variants):
        raise Exception("{} does not exist".format(args.variants))

    # check allele fraction range
    if args.af_slice[0] >= args.af_slice[1]:
        raise Exception("left AF value needs to be smaller than right")
    if args.af_slice[1] > 1:
        raise Exception("AF value needs to be from 0 to 1")

    # get coverage for computing rate per Mb
    coverage = get_context_coverage(args.user_coverage, args.mutsig_genome, args.mutsig_exome)

    if args.all_variants and args.filtered_variants:
        sys.stderr("turning off '--filitered-variants' since '--all-variants' is on")
        args.filtered_variants = False

    # load variant data
    if args.format == "maf":
        muts = read_maf(args.variants, args.af_slice, args.all_variants, args.filtered_variants)
    else:
        if not args.ref_fasta:
            raise Exception("reference fasta is required for VCF input")
        if not os.path.exists(args.ref_fasta):
            raise Exception("{} does not exist".format(args.ref_fasta))
        if args.af_slice[1] < 1:
            sys.stderr("AF slice functionality is not yet supported for VCF input, using cutoff of 1.0")
        muts = read_vcf(args.variants, args.all_variants, args.filtered_variants, args.ref_fasta)

    is_rate_per_Mb = False
    total_muts = sum(muts)

    # compute rate per Mb if coverage is supplied
    if coverage:
        muts = compute_rate_per_Mb(muts, coverage)
        is_rate_per_Mb = True
    # draw lego plot
    generate_lego_plot(muts, total_muts, args.output_prefix, args.plot_title, args.af_slice, is_rate_per_Mb)



if __name__ == "__main__":
    description = "Generate lego plots from MAF or VCF"
    parser = ArgumentParser(description=description)

    parser.add_argument("-o", "--output-prefix",
                        default="lego_plot",
                        help="output prefix (default: %(default)s)")
    parser.add_argument("-s", "--ref-fasta",
                        default=None,
                        help="reference genome fasta, required for VCF input")
    parser.add_argument("-t", "--plot-title",
                        default=None,
                        help="lego plot title (default: '--output-prefix')")

    flt_opts = parser.add_argument_group("filter options",
                                         "filter allele fractions and mutations")
    flt_opts.add_argument("--af-slice",
                          nargs=2,
                          type=float,
                          default=[0.0, 1.0],
                          metavar=("AF_LEFT", "AF_RIGHT"),
                          help="slice specific allele fraction range, LEFT <= AF < RIGHT (not supported for VCF)")
    flt_opts.add_argument("--all-variants",
                          action="store_true",
                          help="plot both PASS and discarded mutations (default: PASS only)")
    flt_opts.add_argument("--filtered-variants",
                          action="store_true",
                          help="plot only discarded mutations (default: PASS only)")

    cov_opts = parser.add_argument_group("coverage options",
                                         "coverage setting for calculating mutation rate per Mb")
    cov_opts.add_argument("--mutsig-genome",
                          action="store_true",
                          help="use precomputed genome coverage in MutSig 2CV")
    cov_opts.add_argument("--mutsig-exome",
                          action="store_true",
                          help="use precomputed exome coverage in MutSig 2CV")
    cov_opts.add_argument("--user-coverage",
                          type=str,
                          help="use custom coverage information (tab-separated, 3-base context and the frequency per line)")

    # required arguments
    parser.add_argument("format", 
                        choices=["maf","vcf"],
                        help="variant format")
    parser.add_argument("variants",
                        help="variant file or standard input")

    args = parser.parse_args()

    try:
        main(args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(__file__)
        sys.exit(prog + ": error: " + str(e))
        
