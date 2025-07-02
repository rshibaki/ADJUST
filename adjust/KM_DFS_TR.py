import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

# ✅ Google Colab の場合、適切なディレクトリを指定
output_dir = "/Users/rshibaki/Documents/project/ADJUST/Figure"  # 必要に応じて変更してください

# ファイルのパス（Google Driveから直接読み込む場合は変更が必要）
input_path = "/Users/rshibaki/Library/CloudStorage/GoogleDrive-ryota.shibaki@gmail.com/マイドライブ/ADJUST試験/DATA/dataset_250129.xlsx"

# Excelからデータ読み込み
df = pd.read_excel(input_path, sheet_name="dataset", engine="openpyxl")

# 必要なデータの取得
columns_needed = ["DFSMONTHS", "DFS_event", "Lymphoid Compartment", "JAK-STAT signaling", "hypoxia"]
df_selected = df[columns_needed].dropna()

# カテゴリごとのグループ名のマッピング
# group_labels = {
#     "Lymphoid Compartment": {0: "High group", 1: "Low group"},
#     "JAK-STAT signaling": {0: "Low group", 1: "High group"},
#     "hypoxia": {0: "Low group", 1: "High group"}
# }

def plot_km_by_category(ax, category_name, colors=["#CE1C48", "#1180D9"], ci_show=True, ci_alpha=0.1, invert_order=False):
    """
    指定したカテゴリでKaplan-Meier生存曲線をプロットする関数。

    Parameters:
    - ax: matplotlibのAxesオブジェクト
    - category_name: カテゴリのカラム名
    - colors: グループごとの色指定（デフォルトは["blue", "red"]）
    - ci_show: 95%信頼区間を表示するかどうか（デフォルト: False）
    - invert_order: グループの表示順を逆にするか（デフォルト: False）
    """

    unique_groups = sorted(df_selected[category_name].dropna().unique())
    if invert_order:
        unique_groups.reverse()  # グループの順番を逆にする

    risk_times = [0, 6, 12, 18, 24, 30, 36, 42]
    risk_counts_by_group = {group: [] for group in unique_groups}

    kmf = KaplanMeierFitter()

    for i, group in enumerate(unique_groups):
        df_group = df_selected[df_selected[category_name] == group].copy()

        # グループ名を取得
        # label_name = group_labels[category_name].get(group, f"{category_name} {group}")
        label_name = "Low group" if group == 0 else "High group"

        # Kaplan-Meierフィッティング
        kmf.fit(df_group["DFSMONTHS"], event_observed=df_group["DFS_event"], label=label_name)
        kmf.plot_survival_function(ax=ax, color=colors[i], linewidth=2, ci_show=ci_show, ci_alpha=ci_alpha)

        # センサリングデータのプロット
        censor_times = df_group["DFSMONTHS"][df_group["DFS_event"] == 0]
        censor_probs = kmf.survival_function_at_times(censor_times).values
        ax.scatter(censor_times, censor_probs, color=colors[i], marker="|", s=100)

        # Number at Risk の計算
        risk_counts_by_group[group] = [df_group[df_group["DFSMONTHS"] >= t].shape[0] for t in risk_times]

    # 軸ラベルと範囲の設定
    ax.set_xlabel("Months since Registration", fontsize=12)
    ax.set_ylabel("Probability of Disease-free survival", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, df_selected["DFSMONTHS"].max() + 5)

    # 補助線
    for x in [12, 24, 36]:
        ax.axvline(x=x, linestyle="--", color="gray", alpha=0.7)

    # タイトル修正
    title_name = "Hypoxia" if category_name == "hypoxia" else category_name
    ax.set_title(f"{title_name}", fontsize=14)

    # Number at Risk の追加（色分け）
    for i, t in enumerate(risk_times):
        for j, group in enumerate(unique_groups):
            ax.text(t, -0.20 - (j * 0.05), f"{risk_counts_by_group[group][i]}", fontsize=12, color=colors[j], ha="center")

    # Number at Risk のラベル
    ax.text(-2, -0.15, "Number at Risk", fontsize=12, ha="right", color="black")
    for i, group in enumerate(unique_groups):
        ax.text(-2, -0.20 - (i * 0.05), "Low group" if group == 0 else "High group", fontsize=12, ha="right", color=colors[i])

# ✅ KMプロットの作成・保存
fig, axes = plt.subplots(1, 3, figsize=(24, 6))
plot_km_by_category(axes[0], "Lymphoid Compartment", colors=["#1180D9", "#CE1C48"])
plot_km_by_category(axes[1], "JAK-STAT signaling")
plot_km_by_category(axes[2], "hypoxia")
plt.tight_layout()

# ✅ PNGとして保存
output_filename = "KaplanMeier_DFS_TR.png"
output_path = os.path.join(output_dir, output_filename)
plt.savefig(output_path, dpi=300, bbox_inches="tight")

##### PD-L1解析 #####
# PD-L1の<1%, 1-49%, 50%で解析
columns_pdl1 = ["DFSMONTHS", "DFS_event", "PD-L1", "PD-L1_np"]
df_pdl1 = df[columns_pdl1].dropna()

pdl1_3g_counts = df_pdl1["PD-L1"] .value_counts()
pdl1_np_counts = df_pdl1["PD-L1_np"].value_counts()

df_dfs_pdl1_0 = df_pdl1[df_pdl1["PD-L1"] == "<1%"]
df_dfs_pdl1_1 = df_pdl1[df_pdl1["PD-L1"] == "1-49%"]
df_dfs_pdl1_2 = df_pdl1[df_pdl1["PD-L1"] == ">50%"]

km_dfs_pdl1_0 = KaplanMeierFitter()
km_dfs_pdl1_1 = KaplanMeierFitter()
km_dfs_pdl1_2 = KaplanMeierFitter()

km_dfs_pdl1_0.fit(df_dfs_pdl1_0["DFSMONTHS"], event_observed=df_dfs_pdl1_0["DFS_event"])
km_dfs_pdl1_1.fit(df_dfs_pdl1_1["DFSMONTHS"], event_observed=df_dfs_pdl1_1["DFS_event"])
km_dfs_pdl1_2.fit(df_dfs_pdl1_2["DFSMONTHS"], event_observed=df_dfs_pdl1_2["DFS_event"])

plt.rcParams["font.family"] = "Arial"
fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(bottom=0.4)

km_dfs_pdl1_0.plot(ax=ax, ci_show=True, ci_alpha=0.1, linewidth=2, color="#1180D9", label="<1%")
km_dfs_pdl1_1.plot(ax=ax, ci_show=True, ci_alpha=0.1, linewidth=2, color="#008080", label="1-49%")
km_dfs_pdl1_2.plot(ax=ax, ci_show=True, ci_alpha=0.1, linewidth=2, color="#CE1C48", label=">50%")

cencor_time_0 = df_dfs_pdl1_0['DFSMONTHS'][df_dfs_pdl1_0['DFS_event'] == 0]
cencor_probs_0 = km_dfs_pdl1_0.survival_function_at_times(cencor_time_0).values
ax.scatter(cencor_time_0, cencor_probs_0, color="#1180D9", marker='|', s=100)

cencor_time_1 = df_dfs_pdl1_1['DFSMONTHS'][df_dfs_pdl1_1['DFS_event'] == 0]
cencor_probs_1 = km_dfs_pdl1_1.survival_function_at_times(cencor_time_1).values
ax.scatter(cencor_time_1, cencor_probs_1, color="#008080", marker='|', s=100)

cencor_time_2 = df_dfs_pdl1_2['DFSMONTHS'][df_dfs_pdl1_2['DFS_event'] == 0]
cencor_probs_2 = km_dfs_pdl1_2.survival_function_at_times(cencor_time_2).values
ax.scatter(cencor_time_2, cencor_probs_2, color="#CE1C48", marker='|', s=100)

ax.set_xlabel("Months since Registration", fontsize=12, labelpad=8)
ax.set_ylabel("Probability of Disease-free Survival", fontsize=12)
ax.set_ylim(0, 1.05)

# 以前の横軸の最大値を使用
x_max = df_pdl1["DFSMONTHS"].max() # 以前の最大値を使用
ax.set_xlim(0, x_max + 5)
ax.set_xticks(np.arange(0, x_max + 5, 12))
ax.set_yticks(np.arange(0, 1.1, 0.2))
ax.minorticks_on()
ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

# 12, 24, 36ヶ月の補助線
for x in [12, 24, 36]:
    ax.axvline(x=x, linestyle="--", color="gray", alpha=0.7)

# Number at Risk
risk_times = [0, 6, 12, 18, 24, 30, 36, 42]
risk_counts_0 = [df_dfs_pdl1_0 [df_dfs_pdl1_0["DFSMONTHS"] >= t].shape[0] for t in risk_times]
risk_counts_1 = [df_dfs_pdl1_1[df_dfs_pdl1_1["DFSMONTHS"] >= t].shape[0] for t in risk_times]
risk_counts_2 = [df_dfs_pdl1_2 [df_dfs_pdl1_2["DFSMONTHS"] >= t].shape[0] for t in risk_times]

for i, (t, r0, r1, r2) in enumerate(zip(risk_times, risk_counts_0, risk_counts_1, risk_counts_2)):
    ax.text(t, -0.35, f"{r0}", fontsize=12, color="#1180D9", ha="center")
    ax.text(t, -0.40, f"{r1}", fontsize=12, color="#008080", ha="center")
    ax.text(t, -0.45, f"{r2}", fontsize=12, color="#CE1C48", ha="center")

ax.text(-2, -0.30, "Number at risk", fontsize=12, ha="right", color="black")
ax.text(-2, -0.35, "<1%", fontsize=12, ha="right", color="#1180D9")
ax.text(-2, -0.40, "1-49%", fontsize=12, ha="right", color="#008080")
ax.text(-2, -0.45, ">50%", fontsize=12, ha="right", color="#CE1C48")

ax.legend(fontsize=10)

output_filename = "KaplanMeier_DFS_PDL1_3gropu.png"
output_path = os.path.join(output_dir, output_filename)
plt.savefig(output_path, dpi=300, bbox_inches="tight")

##### PD-L1のnega, posiで解析 #####
df_dfs_pdl1_n = df_pdl1[df_pdl1["PD-L1_np"] == "negative"]
df_dfs_pdl1_p = df_pdl1[df_pdl1["PD-L1_np"] == "positive"]

km_dfs_pdl1_n = KaplanMeierFitter()
km_dfs_pdl1_p = KaplanMeierFitter()

km_dfs_pdl1_n.fit(df_dfs_pdl1_n["DFSMONTHS"], event_observed=df_dfs_pdl1_n["DFS_event"])
km_dfs_pdl1_p.fit(df_dfs_pdl1_p["DFSMONTHS"], event_observed=df_dfs_pdl1_p["DFS_event"])

plt.rcParams["font.family"] = "Arial"
fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(bottom=0.4)

km_dfs_pdl1_n.plot(ax=ax, ci_show=True, ci_alpha=0.1, linewidth=2, color="#1180D9", label="negative")
km_dfs_pdl1_p.plot(ax=ax, ci_show=True, ci_alpha=0.1, linewidth=2, color="#CE1C48", label="positive")

cencor_time_n = df_dfs_pdl1_n['DFSMONTHS'][df_dfs_pdl1_n['DFS_event'] == 0]
cencor_probs_n = km_dfs_pdl1_n.survival_function_at_times(cencor_time_n).values
ax.scatter(cencor_time_n, cencor_probs_n, color="#1180D9", marker='|', s=100)

cencor_time_p = df_dfs_pdl1_p['DFSMONTHS'][df_dfs_pdl1_p['DFS_event'] == 0]
cencor_probs_p = km_dfs_pdl1_p.survival_function_at_times(cencor_time_p).values
ax.scatter(cencor_time_p, cencor_probs_p, color="#CE1C48", marker='|', s=100)

ax.set_xlabel("Months since Registration", fontsize=12, labelpad=8)
ax.set_ylabel("Probability of Disease-free Survival", fontsize=12)
ax.set_ylim(0, 1.05)

# 以前の横軸の最大値を使用
x_max = df_pdl1["DFSMONTHS"].max() # 以前の最大値を使用
ax.set_xlim(0, x_max + 5)
ax.set_xticks(np.arange(0, x_max + 5, 12))
ax.set_yticks(np.arange(0, 1.1, 0.2))
ax.minorticks_on()
ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

# 12, 24, 36ヶ月の補助線
for x in [12, 24, 36]:
    ax.axvline(x=x, linestyle="--", color="gray", alpha=0.7)

# Number at Risk
risk_times = [0, 6, 12, 18, 24, 30, 36, 42]
risk_counts_n = [df_dfs_pdl1_n [df_dfs_pdl1_n["DFSMONTHS"] >= t].shape[0] for t in risk_times]
risk_counts_p = [df_dfs_pdl1_p[df_dfs_pdl1_p["DFSMONTHS"] >= t].shape[0] for t in risk_times]

for i, (t, r0, r1) in enumerate(zip(risk_times, risk_counts_n, risk_counts_p)):
    ax.text(t, -0.35, f"{r0}", fontsize=12, color="#1180D9", ha="center")
    ax.text(t, -0.40, f"{r1}", fontsize=12, color="#CE1C48", ha="center")

ax.text(-2, -0.30, "Number at risk", fontsize=12, ha="right", color="black")
ax.text(-2, -0.35, "Negative", fontsize=12, ha="right", color="#1180D9")
ax.text(-2, -0.40, "Positive", fontsize=12, ha="right", color="#CE1C48")

ax.legend(fontsize=10)

output_filename = "KaplanMeier_DFS_PDL1_negaposi.png"
output_path = os.path.join(output_dir, output_filename)
plt.savefig(output_path, dpi=300, bbox_inches="tight")

# ✅ 生存解析の統計情報をTXTで保存
output_filename_txt = "KM_DFS_TR_probabilities.txt"
output_path_txt = os.path.join(output_dir, output_filename_txt)


with open(output_path_txt, "w") as f:
    f.write("Kaplan-Meier Survival Analysis Results\n")
    f.write("=" * 50 + "\n\n")

    time_points = [12, 24, 36, 48]

    for category_name in ["Lymphoid Compartment", "JAK-STAT signaling", "hypoxia"]:
        f.write(f"Category: {category_name}\n")
        f.write("-" * 50 + "\n")

        unique_groups = sorted(df_selected[category_name].dropna().unique())
        kmf = KaplanMeierFitter()

        group_survival_probs = {}
        group_medians = {}
        group_event_counts = {}

        for group in unique_groups:
            df_group = df_selected[df_selected[category_name] == group].copy()
            label_name = group_labels[category_name].get(group, f"{category_name} {group}")

            kmf.fit(df_group["DFSMONTHS"], event_observed=df_group["DFS_event"])
            survival_probs = kmf.survival_function_at_times(time_points)

            f.write(f"{label_name}:\n")
            for time, prob in zip(time_points, survival_probs.values.flatten()):
                f.write(f"  {time} months: {prob:.4f}\n")

            median_survival = kmf.median_survival_time_
            f.write(f"  Median Survival: {median_survival:.2f} months\n\n")

            # 保存用データ
            group_survival_probs[group] = survival_probs
            group_medians[group] = median_survival
            group_event_counts[group] = kmf.event_table.observed.sum()

        # ✅ ログランク検定の実施
        df_group_0 = df_selected[df_selected[category_name] == unique_groups[0]]
        df_group_1 = df_selected[df_selected[category_name] == unique_groups[1]]

        logrank_result = logrank_test(
            df_group_0["DFSMONTHS"], df_group_1["DFSMONTHS"],
            event_observed_A=df_group_0["DFS_event"], event_observed_B=df_group_1["DFS_event"]
        )

        p_value = logrank_result.p_value

        # ✅ Cox比例ハザードモデルでHRを正確に推定
        df_cox = df_selected[[category_name, "DFSMONTHS", "DFS_event"]].copy()
        df_cox["group_binary"] = df_cox[category_name].map({
            unique_groups[0]: 0,
            unique_groups[1]: 1
        })

        cph = CoxPHFitter()
        cph.fit(df_cox[["DFSMONTHS", "DFS_event", "group_binary"]],
                duration_col="DFSMONTHS",
                event_col="DFS_event")

        hr = cph.hazard_ratios_["group_binary"]
        ci_lower = np.exp(cph.summary.loc["group_binary", "coef lower 95%"])
        ci_upper = np.exp(cph.summary.loc["group_binary", "coef upper 95%"])

        # ✅ ログランク検定結果を記載
        f.write("Statistical Analysis:\n")
        f.write(f"  Log-rank test P-value: {p_value:.4f}\n")
        f.write(f"  Hazard Ratio (HR): {hr:.4f}\n")
        f.write(f"  95% CI: {ci_lower:.4f} - {ci_upper:.4f}\n")
        f.write("=" * 50 + "\n\n")
        
    ###PD-L1###
    f.write("Kaplan-Meier Survival Analysis of PD-L1\n")
    f.write("=" * 50 + "\n")
    f.write(f"{pdl1_3g_counts.to_string()}\n")
    f.write(f"{pdl1_np_counts.to_string()}\n")

    ### PD-L1 (3群) 年次生存割合だけ
    f.write("Survival Probabilities (PD-L1 3-group):\n")
    f.write("-" * 50 + "\n")
    km_pdl1_0 = KaplanMeierFitter().fit(df_dfs_pdl1_0["DFSMONTHS"], df_dfs_pdl1_0["DFS_event"])
    km_pdl1_1 = KaplanMeierFitter().fit(df_dfs_pdl1_1["DFSMONTHS"], df_dfs_pdl1_1["DFS_event"])
    km_pdl1_2 = KaplanMeierFitter().fit(df_dfs_pdl1_2["DFSMONTHS"], df_dfs_pdl1_2["DFS_event"])

    for label, km in zip(["<1%", "1-49%", ">50%"], [km_pdl1_0, km_pdl1_1, km_pdl1_2]):
        survival_probs = km.survival_function_at_times(time_points)
        f.write(f"{label}:\n")
        for t, prob in zip(time_points, survival_probs.values.flatten()):
            f.write(f"  {t} months: {prob:.4f}\n")
        f.write("\n")

    ### PD-L1_np (nega vs posi) 年次生存割合 + ログランク・HR
    f.write("Survival Probabilities (PD-L1 negative vs positive):\n")
    f.write("-" * 50 + "\n")
    km_pdl1_n = KaplanMeierFitter().fit(df_dfs_pdl1_n["DFSMONTHS"], df_dfs_pdl1_n["DFS_event"])
    km_pdl1_p = KaplanMeierFitter().fit(df_dfs_pdl1_p["DFSMONTHS"], df_dfs_pdl1_p["DFS_event"])

    for label, km in zip(["Negative", "Positive"], [km_pdl1_n, km_pdl1_p]):
        survival_probs = km.survival_function_at_times(time_points)
        f.write(f"{label}:\n")
        for t, prob in zip(time_points, survival_probs.values.flatten()):
            f.write(f"  {t} months: {prob:.4f}\n")
        f.write("\n")

    # ログランクとHR
    logrank_result_pdl1_np = logrank_test(
        df_dfs_pdl1_n["DFSMONTHS"], df_dfs_pdl1_p["DFSMONTHS"],
        event_observed_A=df_dfs_pdl1_n["DFS_event"], event_observed_B=df_dfs_pdl1_p["DFS_event"]
    )
    p_value_pdl1_np = logrank_result_pdl1_np.p_value

    cph_pdl1 = CoxPHFitter()
    df_pdl1["PD-L1_np_binary"] = df_pdl1["PD-L1_np"].map({"negative": 0, "positive": 1})
    cph_pdl1.fit(df_pdl1[["DFSMONTHS", "DFS_event", "PD-L1_np_binary"]],
                         duration_col="DFSMONTHS",
                         event_col="DFS_event")
    ci_lower = np.exp(cph_pdl1.confidence_intervals_.iloc[0, 0])
    ci_upper = np.exp(cph_pdl1.confidence_intervals_.iloc[0, 1])

    f.write("Statistical Analysis (PD-L1 negative vs positive):\n")
    f.write(f"  Log-rank test P-value: {p_value_pdl1_np:.4f}\n")
    f.write(f"  Hazard Ratio (HR): {cph_pdl1.hazard_ratios_["PD-L1_np_binary"]:.4f}\n")
    f.write(f"  95% CI: {ci_lower:.4f} - {ci_upper:.4f}\n")
    f.write("=" * 50 + "\n\n")


print(f"✅ Kaplan-Meier 曲線を保存しました: {output_path}")
print(f"✅ 生存解析の統計情報を保存しました: {output_path_txt}") 
