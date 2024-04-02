
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import upsetplot
from matplotlib_venn import venn2, venn3
from sklearn.metrics import pairwise_distances

from ..._utils import convert_vcf_to_json, read_vcf_into_df


class VariantComparator:
    def __init__(self, variant_files: list[dict], bed_file: str = None):
        """
        args:
            variant_files: list of dicts with keys:
                - path: path to the variant file
                - caller: caller type of the variant file
                - mapper: mapper type of the variant file
        """
        self.variant_files = variant_files
        self.color_palette = ("skyblue", "lightcoral", "lightgreen")
        self.variants_df = self._read_variant_files()

    def _read_variant_files(self):
        variants_df = pd.DataFrame()
        for variant_file in self.variant_files:
            # Convert vcf to tsv
            variants_json = convert_vcf_to_json(
                variant_file["path"], variant_file["caller"]
            )
            df = pd.json_normalize(variants_json)

            df["pipeline"] = variant_file["mapper"] + "_" + variant_file["caller"]

            # If caller is vardict, the status of the variant should be StrongSomatic
            if variant_file["caller"] == "vardict":
                df = df[df["STATUS"].str.contains("StrongSomatic")]

            # Set variant_id as CHROM-POS-REF-ALT
            df["variant_id"] = df[["CHROM", "POS", "REF", "ALT"]].apply(
                lambda x: "-".join(x.astype(str)), axis=1
            )

            df.rename(
                columns={
                    "AF": "af",
                    "AD": "ad",
                    "DP": "dp",
                },
                inplace=True,
            )

            # Keep only "variant_id, caller, mapper" columns
            df = df[["variant_id", "af", "pipeline"]]

            variants_df = pd.concat([variants_df, df], ignore_index=True)

        return variants_df

    def draw_upset_plot(self):
        upset_df = self.variants_df.copy()
        upset_df["value"] = True
        upset_df = upset_df.pivot(
            index="variant_id", columns="pipeline", values="value"
        ).replace(pd.NA, False)

        upset_df = pd.merge(
            upset_df,
            self.variants_df[["variant_id", "af"]],
            on="variant_id",
            how="left",
        )

        upset_df = (
            upset_df.groupby(
                ["variant_id"]
                + [f"{i['mapper']}_{i['caller']}" for i in self.variant_files]
            )
            .agg(
                {
                    "af": "mean",
                }
            )
            .reset_index()
        )

        upset_df = upset_df.drop("variant_id", axis=1)

        # Set all callers as index and af as value
        upset_df = upset_df.set_index(list(upset_df.columns[:-1]))

        upset = upsetplot.UpSet(
            upset_df,
            sort_by="-degree",
        )

        # Add catplot which shows allele frequency distribution of each cardinality
        upset.add_catplot(
            value="af",
            kind="strip",
            color="blue",
            size=2.5,
        )

        upset.plot()

    def draw_similarity_plot(self):
        sim_df = self.variants_df.copy()

        sim_df = sim_df.pivot_table(
            index="variant_id", values="pipeline", aggfunc=lambda x: ",".join(x)
        )["pipeline"].str.get_dummies(sep=",")

        jac_sim = 1 - pairwise_distances(sim_df.T, metric="hamming")
        jac_sim = pd.DataFrame(jac_sim, columns=sim_df.columns, index=sim_df.columns)

        sns.heatmap(jac_sim, annot=True, cmap="YlGnBu")

    def draw_venn2_plot(self, pipeline1: str, pipeline2: str):
        venn2(
            subsets=(
                set(
                    self.variants_df[self.variants_df["pipeline"] == pipeline1][
                        "variant_id"
                    ].tolist()
                ),
                set(
                    self.variants_df[self.variants_df["pipeline"] == pipeline2][
                        "variant_id"
                    ].tolist()
                ),
            ),
            set_labels=(pipeline1, pipeline2),
        )

    def draw_venn3_plot(self, pipeline1: str, pipeline2: str, pipeline3: str):
        venn3(
            subsets=(
                set(
                    self.variants_df[self.variants_df["pipeline"] == pipeline1][
                        "variant_id"
                    ].tolist()
                ),
                set(
                    self.variants_df[self.variants_df["pipeline"] == pipeline2][
                        "variant_id"
                    ].tolist()
                ),
                set(
                    self.variants_df[self.variants_df["pipeline"] == pipeline3][
                        "variant_id"
                    ].tolist()
                ),
            ),
            set_labels=(pipeline1, pipeline2, pipeline3),
        )

    def draw_precision_recall_plot(self, truth_vcf: str):
        """
        Scatter plot of precision and recall of each caller where x-axis is precision and y-axis is recall.
        """

        precision_recall_values = self._create_precision_recall_dict(truth_vcf)
        precision_recall_df = pd.DataFrame(precision_recall_values).T
        precision_recall_df.columns = ["precision", "recall"]

        g = sns.scatterplot(
            data=precision_recall_df,
            x="precision",
            y="recall",
            hue=precision_recall_df.index,
            s=150,
        )

    def _create_precision_recall_dict(self, truth_vcf: str):
        """
        Calculate precision and recall of each caller with respect to the truth vcf file.
        """

        truth_df = read_vcf_into_df(truth_vcf)

        truth_df["variant_id"] = truth_df[["CHROM", "POS", "REF", "ALT"]].apply(
            lambda x: "-".join(x.astype(str)), axis=1
        )

        truth_set = set(truth_df["variant_id"].tolist())

        precision_recall_values = {}
        for pipeline in self.variants_df["pipeline"].unique():
            pipeline_set = set(
                self.variants_df[self.variants_df["pipeline"] == pipeline][
                    "variant_id"
                ].tolist()
            )

            precision_value = self._precision(pipeline_set, truth_set)
            recall_value = self._recall(pipeline_set, truth_set)

            precision_recall_values[pipeline] = (precision_value, recall_value)

        return precision_recall_values

    def get_variants(self):
        return self.variants_df

    def create_intersection_bed(self, query: str, output_file: str):
        """
        Create a bed file with the intersection of the variants with +- 50 bp
        bed_format: CHROM, START, END

        args:
            pipelines: list of pipelines
            intersections: list of intersections in format of 'pipeline1 & pipeline2 ~ pipeline3'
            output_file: output file path
        """
        operations = {
            "|": {"precedence": 1, "function": lambda x, y: x.union(y)},
            "&": {"precedence": 2, "function": lambda x, y: x.intersection(y)},
            "~": {"precedence": 3, "function": lambda x, y: x.difference(y)},
        }

        def shunting_yard(query):
            output = []
            stack = []
            for token in query.split():
                if token in operations:
                    while (
                        stack
                        and stack[-1] in operations
                        and operations[token]["precedence"]
                        <= operations[stack[-1]]["precedence"]
                    ):
                        output.append(stack.pop())
                    stack.append(token)
                else:
                    output.append(token)
            while stack:
                output.append(stack.pop())
            return output

        postfix = shunting_yard(query)

        stack = []
        for token in postfix:
            if token in operations:
                y = stack.pop()
                x = stack.pop()

                x = (
                    set(
                        self.variants_df[self.variants_df["pipeline"] == x][
                            "variant_id"
                        ].tolist()
                    )
                    if isinstance(x, str)
                    else x
                )
                y = (
                    set(
                        self.variants_df[self.variants_df["pipeline"] == y][
                            "variant_id"
                        ].tolist()
                    )
                    if isinstance(y, str)
                    else y
                )

                results = operations[token]["function"](x, y)
                stack.append(results)
            else:
                stack.append(token)

        results = stack[0]

        bed_df = pd.DataFrame()
        for variant in results:
            chrom, pos, ref, alt = variant.split("-")
            start = int(pos) - 40
            end = int(pos) + 40
            bed_df = pd.concat(
                [
                    bed_df,
                    pd.DataFrame({"CHROM": [chrom], "START": [start], "END": [end]}),
                ]
            )

        bed_df.to_csv(output_file, sep="\t", index=False, header=False)

    def _precision(self, retrieved_set: set, relevant_set: set) -> float:
        if len(retrieved_set) == 0:
            return 0.0

        correct_findings = self._intersection_size(retrieved_set, relevant_set)
        precision_value = correct_findings / len(retrieved_set)
        return precision_value

    def _recall(self, retrieved_set: set, relevant_set: set) -> float:
        if len(relevant_set) == 0:
            return 0.0

        correct_findings = self._intersection_size(retrieved_set, relevant_set)
        recall_value = correct_findings / len(relevant_set)
        return recall_value

    def _intersection_size(self, a: set, *b: set) -> int:
        return len(a.intersection(*b))
