import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

print("Starting metagenomics downstream analysis...")

# 1. Load the HUMAnN pathway relative abundance table
df = pd.read_csv("merged_pathabundance_relab.tsv", sep="\t")

# Ensure the first column is named 'Pathway'
df = df.rename(columns={df.columns[0]: "Pathway"})

# 2. Data Cleaning & Filtering
# Remove stratified rows (lines containing species breakdown denoted by "|")
df_unstrat = df[~df["Pathway"].str.contains("\\|", regex=True)].copy()

# Remove HUMAnN structural/special rows
df_unstrat = df_unstrat[
    ~df_unstrat["Pathway"].isin(["UNMAPPED", "UNINTEGRATED"])
]

# Set pathway names as row index for processing
df_unstrat = df_unstrat.set_index("Pathway")

# Clean sample names by stripping common suffixes
df_unstrat.columns = (
    df_unstrat.columns
    .str.replace("_merged_Abundance", "", regex=False)
    .str.replace("_Abundance", "", regex=False)
)

print("Samples detected in table:")
print(df_unstrat.columns.tolist())

# 3. Define Experimental Cohort Groups
atlanta = ["SRR15202440", "SRR15202441", "SRR15202442"]
maputo = ["SRR15209203", "SRR15209176", "SRR15209171"]

all_samples = atlanta + maputo
df_unstrat = df_unstrat[all_samples]

# 4. Extract Top 20 Overall Pathways
top20 = df_unstrat.mean(axis=1).sort_values(ascending=False).head(20)
top20.to_csv("top20_pathways.csv", header=["Mean_Relative_Abundance"])

print("\nTop 20 most abundant pathways calculated.")

# 5. Visualizations: Barplot & Heatmap
# Horizontal Bar Plot
plt.figure(figsize=(10, 8))
top20.sort_values().plot(kind="barh", color="teal")
plt.xlabel("Mean Relative Abundance")
plt.ylabel("Pathway")
plt.title("Top 20 Most Abundant Metabolic Pathways")
plt.tight_layout()
plt.savefig("figures/top20_pathways_barplot.png", dpi=300)
plt.close()

# Cluster Heatmap
heatmap_data = df_unstrat.loc[top20.index]
plt.figure(figsize=(12, 10))
sns.heatmap(heatmap_data, cmap="viridis", annot=False)
plt.title("Top 20 Pathway Abundances Across Samples")
plt.xlabel("Sample ID")
plt.ylabel("Pathway")
plt.tight_layout()
plt.savefig("figures/top20_pathways_heatmap.png", dpi=300)
plt.close()

# 6. Group Comparisons (Atlanta vs. Maputo)
group_df = pd.DataFrame({
    "Atlanta_mean": df_unstrat[atlanta].mean(axis=1),
    "Maputo_mean": df_unstrat[maputo].mean(axis=1)
})

group_df["Difference_Maputo_minus_Atlanta"] = (
    group_df["Maputo_mean"] - group_df["Atlanta_mean"]
)

# Safe fold change calculation avoiding division by zero errors
group_df["Fold_change_Maputo_over_Atlanta"] = (
    (group_df["Maputo_mean"] + 1e-9) / (group_df["Atlanta_mean"] + 1e-9)
)

group_df = group_df.sort_values("Difference_Maputo_minus_Atlanta", ascending=False)
group_df.to_csv("atlanta_vs_maputo_pathway_means.csv")

# Identify top 20 differing pathways by absolute magnitude
top_diff = group_df.reindex(
    group_df["Difference_Maputo_minus_Atlanta"].abs()
    .sort_values(ascending=False)
    .head(20)
    .index
)
top_diff.to_csv("top_differing_pathways.csv")

# Horizontal Bar Plot for Differences
plt.figure(figsize=(10, 8))
top_diff["Difference_Maputo_minus_Atlanta"].sort_values().plot(kind="barh", color="coral")
plt.axvline(0, color="black", linewidth=0.8, linestyle="--")
plt.xlabel("Difference in Mean Relative Abundance (Maputo - Atlanta)")
plt.ylabel("Pathway")
plt.title("Top Pathway Differences: Maputo vs. Atlanta")
plt.tight_layout()
plt.savefig("figures/top_differing_pathways_barplot.png", dpi=300)
plt.close()

# 7. Non-parametric Statistical Testing (Mann-Whitney U)
stats_results = []

for pathway in df_unstrat.index:
    atl_vals = df_unstrat.loc[pathway, atlanta]
    map_vals = df_unstrat.loc[pathway, maputo]

    try:
        stat, p = mannwhitneyu(atl_vals, map_vals, alternative="two-sided")
    except ValueError:
        p = None

    stats_results.append({
        "Pathway": pathway,
        "Atlanta_mean": atl_vals.mean(),
        "Maputo_mean": map_vals.mean(),
        "Difference_Maputo_minus_Atlanta": map_vals.mean() - atl_vals.mean(),
        "p_value": p
    })

stats_df = pd.DataFrame(stats_results)
stats_df = stats_df.sort_values("p_value", na_position="last")
stats_df.to_csv("pathway_stats_mannwhitney.csv", index=False)

print("\nAnalysis complete! Figures saved to the 'figures/' directory.")
