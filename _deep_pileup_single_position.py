import pandas as pd
import numpy as np
import glob
import time
import os
import sys

from tqdm import tqdm as tqdm
import subprocess
import json
import datetime


import pandas as pd
import plotly.graph_objects as go
from optparse import OptionParser


def scanPileup(path, colNo=4, minFreq=0.9, minAbs=2, minHetSnp=0.25):
    headerList = (
        [
            "Cohort_File",
            "Samples",
            "Reads",
            "Rel_A",
            "Rel_C",
            "Rel_G",
            "Rel_T",
            f"SNPs_AF>{int(100*minHetSnp)}_%",
        ]
        + [f"min_{minFreq:.2f}_{n}" for n in ["A", "C", "G", "T"]]
        + [f"min_{minAbs}_{n}" for n in ["A", "C", "G", "T"]]
    )

    output_path_to_overview_file = os.path.join(path, "Overview.tsv")
    with open(output_path_to_overview_file, "w") as f:
        f.write("\t".join(headerList) + "\n")

        paths_to_iterate = glob.glob(os.path.join(path, "pileup*.txt"))
        for file in paths_to_iterate:
            sA = 0
            sC = 0
            sG = 0
            sT = 0
            minA = 0
            minC = 0
            minG = 0
            minT = 0
            minAbsA = 0
            minAbsC = 0
            minAbsG = 0
            minAbsT = 0

            hetPos = 0
            n = 0
            with open(file, "r") as j:
                for line in j:
                    pileup = line.split("\t")[colNo].upper()
                    n += 1
                    A = pileup.count("A")
                    C = pileup.count("C")
                    G = pileup.count("G")
                    T = pileup.count("T")
                    count = float(A + C + G + T)
                    first = False
                    if count == 0:
                        continue
                    for B in (A, C, G, T):
                        if first and B / count >= minHetSnp:
                            hetPos += 1
                            break
                        if not first and B / count >= minHetSnp:
                            first = True

                    # Patients with min Frequency
                    if A / count >= minFreq:
                        minA += 1
                    if C / count >= minFreq:
                        minC += 1
                    if G / count >= minFreq:
                        minG += 1
                    if T / count >= minFreq:
                        minT += 1

                    # Patients with min Absolute Count
                    if A >= minAbs:
                        minAbsA += 1
                    if C >= minAbs:
                        minAbsC += 1
                    if G >= minAbs:
                        minAbsG += 1
                    if T >= minAbs:
                        minAbsT += 1

                    sA += A
                    sC += C
                    sG += G
                    sT += T

                # summerize pileup file
                cnt = sA + sC + sG + sT

                # Calculate frequencies and percentages
                freq_A = sA * 100.0 / cnt
                freq_C = sC * 100.0 / cnt
                freq_G = sG * 100.0 / cnt
                freq_T = sT * 100.0 / cnt
                het_percentage = hetPos * 100.0 / n

                # Calculate minimum allele frequencies and percentages

                output = (
                    [
                        file.split("/")[-1],
                        str(n),
                        str(cnt),
                        freq_A,
                        freq_C,
                        freq_G,
                        freq_T,
                        het_percentage,
                    ]
                    + [
                        "%.4f" % (minA * 100.0 / n),
                        "%.4f" % (minC * 100.0 / n),
                        "%.4f" % (minG * 100.0 / n),
                        "%.4f" % (minT * 100.0 / n),
                    ]
                    + [
                        "%.4f" % (minAbsA * 100.0 / n),
                        "%.4f" % (minAbsC * 100.0 / n),
                        "%.4f" % (minAbsG * 100.0 / n),
                        "%.4f" % (minAbsT * 100.0 / n),
                    ]
                )
            f.write("\t".join([str(i) for i in output]) + "\n")
    return output_path_to_overview_file


def _sort_overview_file(path_to_overview_file: str):
    overview = pd.read_csv(path_to_overview_file, delimiter="\t")

    cohort_dict = {}
    for idx, row in overview.iterrows():
        cohort = row["Cohort_File"].split("_")[-1].replace(".txt", "")
        control_or_tumor = row["Cohort_File"].split("_")[1]
        if cohort not in cohort_dict:
            cohort_dict[cohort] = {}
        cohort_dict[cohort][control_or_tumor] = idx

    line_order = []
    sorted_cohorts = sorted(cohort_dict.keys())
    for cohort in sorted_cohorts:
        line_order.append(cohort_dict[cohort]["tumor"])
        line_order.append(cohort_dict[cohort]["control"])

    overview = overview.iloc[line_order]
    overview.to_csv(path_to_overview_file, index=None, sep="\t")
    return path_to_overview_file


def _plot_af_greater_than_25(
    path_to_overview_file: str,
    gene_name: str,
    chromosome: str,
    position: str,
    only_relevant: bool = True,
):

    data = pd.read_csv(path_to_overview_file, delimiter="\t")

    cohorts = {}
    for idx, row in data.iterrows():
        cohort_file = row["Cohort_File"]
        af = float(row["SNPs_AF>25_%"])

        cohort = cohort_file.split("_")[2]
        if cohort not in cohorts:
            cohorts[cohort] = {"tumor": -1, "control": -1}
        if "control" in cohort_file:
            cohorts[cohort]["control"] = float(af)
        else:
            cohorts[cohort]["tumor"] = float(af)

    num_original_cohorts = len(cohorts.keys())

    if only_relevant:
        cohorts_to_remove = []
        for cohort in cohorts:
            control_value = cohorts[cohort]["control"]
            tumor_value = cohorts[cohort]["tumor"]
            if control_value == 0 and tumor_value == 0:
                cohorts_to_remove.append(cohort)
        for cohort in cohorts_to_remove:
            cohorts.pop(cohort, None)

    x_axis = list(cohorts.keys())

    num_current_cohorts = len(x_axis)

    if len(x_axis) == 0:
        return _plot_af_greater_than_25(path_to_overview_file, False)

    else:

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                x=x_axis,
                y=[cohorts[cohort]["tumor"] for cohort in x_axis],
                mode="markers",
                marker_color="red",
                marker={"size": 10},
                name="Tumor",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=x_axis,
                y=[cohorts[cohort]["control"] for cohort in x_axis],
                mode="markers",
                marker_color="green",
                marker={"symbol": "circle-x-open", "size": 10},
                name="Control",
            )
        )

        fig.update_layout(
            title=f"<b>Patients with a minor allele frequency > 25%</b><br><sup>Gene Name: {gene_name} / Position: chr{chromosome}:{position} / Displaying {num_current_cohorts} of {num_original_cohorts} cohorts",
            xaxis_title="Cohorts",
            yaxis_title="Percent of Patients",
            yaxis_range=[-2, 100],
        )

        fig.update_xaxes(tickangle=45)

        return fig


def _plot_at_least_two_variant_alleles(
    path_to_overview_file: str,
    gene_name: str,
    chromosome: str,
    position: str,
    only_relevant: bool = True,
):
    data = pd.read_csv(path_to_overview_file, delimiter="\t")

    cohorts = {}

    for idx, row in data.iterrows():
        # Get the maximum.
        maximum = np.max(
            [row["min_2_A"], row["min_2_C"], row["min_2_G"], row["min_2_T"]]
        )

        # Get the cohort.
        cohort = row["Cohort_File"].split("_")[2]
        tumor_or_control = row["Cohort_File"].split("_")[1]

        # Get the counts.
        min_2_A = float(row["min_2_A"]) if float(row["min_2_A"]) != maximum else 0
        min_2_C = float(row["min_2_C"]) if float(row["min_2_C"]) != maximum else 0
        min_2_G = float(row["min_2_G"]) if float(row["min_2_G"]) != maximum else 0
        min_2_T = float(row["min_2_T"]) if float(row["min_2_T"]) != maximum else 0

        if cohort not in cohorts:
            cohorts[cohort] = {
                "tumor": {"A": -1, "C": -1, "G": -1, "T": -1},
                "control": {"A": -1, "C": -1, "G": -1, "T": -1},
            }
        cohorts[cohort][tumor_or_control]["A"] = min_2_A
        cohorts[cohort][tumor_or_control]["C"] = min_2_C
        cohorts[cohort][tumor_or_control]["G"] = min_2_G
        cohorts[cohort][tumor_or_control]["T"] = min_2_T

    num_original_cohorts = len(cohorts.keys())

    if only_relevant:
        cohorts_to_remove = []
        for cohort in cohorts:
            tumor_sum = np.sum(
                [
                    cohorts[cohort]["tumor"][nucleotide]
                    for nucleotide in ["A", "C", "G", "T"]
                ]
            )
            control_sum = np.sum(
                [
                    cohorts[cohort]["control"][nucleotide]
                    for nucleotide in ["A", "C", "G", "T"]
                ]
            )
            if (tumor_sum == 0) and (control_sum == 0):
                cohorts_to_remove.append(cohort)

        for cohort in cohorts_to_remove:
            cohorts.pop(cohort, None)

    x_axis = list(cohorts.keys())

    num_current_cohorts = len(x_axis)

    if len(x_axis) == 0:
        return _plot_at_least_two_variant_alleles(path_to_overview_file, False)
    else:

        fig = go.Figure()

        color_dict = {"A": "green", "C": "blue", "G": "orange", "T": "red"}

        for nucleotide in ["A", "C", "T", "G"]:
            fig.add_trace(
                go.Scatter(
                    x=[cohort for cohort in cohorts],
                    y=[cohorts[cohort]["tumor"][nucleotide] for cohort in cohorts],
                    name=f"{nucleotide} (tumor)",
                    mode="markers",
                    marker_color=color_dict[nucleotide],
                    marker={"size": 10},
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=[cohort for cohort in cohorts],
                    y=[cohorts[cohort]["control"][nucleotide] for cohort in cohorts],
                    name=f"{nucleotide} (control)",
                    mode="markers",
                    marker_color=color_dict[nucleotide],
                    marker={"symbol": "circle-x-open", "size": 10},
                )
            )

            fig.update_layout(
                title=f"<b>Patients with at least 2 Variant Alleles</b><br><sup>Gene Name: {gene_name} / Position: chr{chromosome}:{position} / Displaying {num_current_cohorts} of {num_original_cohorts} cohorts",
                xaxis_title="Cohorts",
                yaxis_title="Percent of Affected Patients",
                yaxis_range=[-2, 100],
            )

            fig.update_xaxes(tickangle=45)

        return fig


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--path-to-metadata",
        action="store",
        type="str",
        dest="metadata",
        help="Path to the metadata file. \n",
    )

    parser.add_option(
        "--output-path-to-repository",
        action="store",
        type="str",
        dest="output_path_to_repository",
        help="Output path to the repository \n",
    )

    parser.add_option(
        "--gene",
        action="store",
        type="str",
        dest="gene",
        default=None,
        help="Gene to scan Deep Pileup. If unknown, this can be left blank. \n",
    )

    parser.add_option(
        "--chromosome",
        action="store",
        type="str",
        dest="chromosome",
        help="Chromosome to scan Deep Pileup. \n",
    )

    parser.add_option(
        "--position",
        action="store",
        type="str",
        dest="position",
        help="Genomic position to scan Deep Pileup. \n",
    )

    (options, args) = parser.parse_args()

    path_to_metadata = options.metadata
    output_path_to_repository = options.output_path_to_repository
    gene = str(options.gene) if options.gene is not None else "no_gene_specified"
    chromosome = str(options.chromosome)
    position = int(options.position)

    # Load in metadataframe.
    metadata = pd.read_csv(path_to_metadata)

    subset_of_data_with_bams = metadata[
        (metadata["path_to_control_bam"] != "not_available")
        | (metadata["path_to_tumor_bam"] != "not_available")
    ]
    subset_of_data_with_bams = subset_of_data_with_bams[
        ["pid", "cohort", "path_to_control_bam", "path_to_tumor_bam"]
    ].drop_duplicates()
    subset_of_data_with_bams.reset_index(inplace=True, drop=True)

    # Check if a deep pileup of this already exists. If so, skip.
    gene_directory = os.path.join(output_path_to_repository, gene)
    if not os.path.exists(gene_directory):
        os.mkdir(gene_directory)

    chromosome_and_position_directory = os.path.join(
        gene_directory, f"chr{chromosome}:{position}"
    )
    if not os.path.exists(chromosome_and_position_directory):
        os.mkdir(chromosome_and_position_directory)
    else:
        print(f"{gene} (chr{chromosome}:{position}) already exists. Skipping...")
        sys.exit()

    # Define a dictionary that will contain the cohort, control/tumor, paths of the bams, samtools commands that's being called and the samtools output.
    samtools_dictionary = {}

    """
    samtools_dictionary[cohort]["control"|"tumor"] = {
        "paths": [],
        "commands": [],
        "output": []
    }
    """

    # Iterate through each row within the dataframe and run samtools on all available .bam files.
    for idx, row in tqdm(
        subset_of_data_with_bams.iterrows(),
        total=subset_of_data_with_bams.shape[0],
        desc="Calling .bam files.",
    ):
        cohort = row["cohort"]

        if cohort not in samtools_dictionary:
            samtools_dictionary[cohort] = {
                "tumor": {"paths": [], "commands": [], "output": []},
                "control": {"paths": [], "commands": [], "output": []},
            }

        path_to_control_bam = row["path_to_control_bam"]
        path_to_tumor_bam = row["path_to_tumor_bam"]

        # Check if the _control_ file exists.
        if os.path.exists(path_to_control_bam):
            # Add the path to the dictionary.
            samtools_dictionary[cohort]["control"]["paths"].append(path_to_control_bam)

            # Add the command to the dictionary.
            command = f"module load samtools/1.14; samtools view -b {path_to_control_bam} {chromosome}:{int(position)}-{int(position) + 1} | samtools mpileup --min-BQ 0 - | grep '{position}'"
            samtools_dictionary[cohort]["control"]["commands"].append(command)

            # Add the results after running to the dictionary.
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            samtools_dictionary[cohort]["control"]["output"].append(result.stdout)

        # Check if the _tumor_ file exists.
        if os.path.exists(path_to_tumor_bam):
            # Add the path to the dictionary.
            samtools_dictionary[cohort]["tumor"]["paths"].append(path_to_tumor_bam)

            # Add the command to the dictionary.
            command = f"module load samtools/1.14; samtools view -b {path_to_tumor_bam} {chromosome}:{int(position)}-{int(position) + 1} | samtools mpileup --min-BQ 0 - | grep '{position}'"
            samtools_dictionary[cohort]["tumor"]["commands"].append(command)

            # Add the results after running to the dictionary.
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            samtools_dictionary[cohort]["tumor"]["output"].append(result.stdout)

    # Save the samtools_dictionary .json file.
    current_datetime = datetime.datetime.now().strftime("%d-%m-%Y_%H:%M:%S")
    with open(
        os.path.join(
            chromosome_and_position_directory,
            f"samtools_dictionary_{current_datetime}.json",
        ),
        "w",
    ) as f:
        json.dump(samtools_dictionary, f, indent=4)

    # For each cohort in samtools_dictionary, write control and tumor .json files with all outputs.
    for cohort in samtools_dictionary:
        for tumor_or_control in samtools_dictionary[cohort]:
            with open(
                os.path.join(
                    chromosome_and_position_directory,
                    f"pileup_{tumor_or_control}_{cohort}.txt",
                ),
                "w",
            ) as f:
                for line in samtools_dictionary[cohort][tumor_or_control]["output"]:
                    f.write(line)

    # Create the overview file from the pileups.
    output_path_to_overview_file = scanPileup(chromosome_and_position_directory)

    # Sort the overview file.
    output_path_to_overview_file = _sort_overview_file(output_path_to_overview_file)

    # Plot and save with an AF > 25.
    fig = _plot_af_greater_than_25(
        path_to_overview_file=output_path_to_overview_file,
        gene_name=gene,
        chromosome=chromosome,
        position=position,
        only_relevant=True,
    )
    fig.write_image(
        os.path.join(
            chromosome_and_position_directory, "af_greater_than_25_only_relevant.png"
        )
    )

    fig = _plot_af_greater_than_25(
        path_to_overview_file=output_path_to_overview_file,
        gene_name=gene,
        chromosome=chromosome,
        position=position,
        only_relevant=False,
    )
    fig.write_image(
        os.path.join(chromosome_and_position_directory, "af_greater_than_25_all.png")
    )

    # Plot and save at least two variant alleles.
    fig = _plot_at_least_two_variant_alleles(
        path_to_overview_file=output_path_to_overview_file,
        gene_name=gene,
        chromosome=chromosome,
        position=position,
        only_relevant=True,
    )
    fig.write_image(
        os.path.join(
            chromosome_and_position_directory,
            "at_least_two_variant_alleles_only_relevant.png",
        )
    )

    fig = _plot_at_least_two_variant_alleles(
        path_to_overview_file=output_path_to_overview_file,
        gene_name=gene,
        chromosome=chromosome,
        position=position,
        only_relevant=False,
    )
    fig.write_image(
        os.path.join(
            chromosome_and_position_directory, "at_least_two_variant_alleles_all.png"
        )
    )
