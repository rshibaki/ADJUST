import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import median_survival_times

#出力先フォルダを指定
output_dir = "/Users/rshibaki/Documents/project/ADJUST/Figure"

file_path = "/Users/rshibaki/Library/CloudStorage/GoogleDrive-ryota.shibaki@gmail.com/マイドライブ/ADJUST試験/DATA/dataset_250129.xlsx"
df = pd.read_excel(file_path, sheet_name="dataset", engine="openpyxl")

# OS_mainの作図
df_os_main = df[df["STUDY"] == 0][["OSMONTHS", "OS_event"]].dropna()
df_os_main.columns = ["OSMONTHS", "OS_event"]

km_os_main = KaplanMeierFitter()
km_os_main.fit(df_os_main["OSMONTHS"], event_observed=df_os_main["OS_event"])

plt.rcParams["font.family"] = "Arial"
fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(bottom=0.35)

# 生存曲線と95%信頼区間をプロット
km_os_main.plot_survival_function(ax=ax, ci_show=True, linewidth=2, color="blue", legend=False)

# ax.step(np.insert(km_os_main.timeline, 0, 0), np.insert(km_os_main.survival_function_.values.flatten(), 0, 1.0), 
#     where="post", label="Kaplan-Meier Estimate", linewidth=2, color="blue")

censor_times = df_os_main["OSMONTHS"][df_os_main["OS_event"] == 0]
censor_probs = km_os_main.survival_function_at_times(censor_times).values

ax.scatter(censor_times, censor_probs, color='blue', marker='|', s=100, label="Censored")

ax.set_xlabel("Time (months)", fontsize=12, labelpad=8)
ax.set_ylabel("Survival Probability", fontsize=12)
ax.set_ylim(0, 1.05)
ax.set_xlim(0, df_os_main["OSMONTHS"].max() + 5)

ax.set_xticks(np.arange(0, df_os_main["OSMONTHS"].max() + 12, 12))
ax.set_yticks(np.arange(0, 1.1, 0.2))

ax.minorticks_on()
ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

for x in [12, 24, 36]:
    ax.axvline(x=x, linestyle="--", color="gray", alpha=0.7)

for x in [12, 24, 36]:
    if x in km_os_main.timeline:
        survival_prob = km_os_main.survival_function_at_times([x]).values[0]
        ax.text(x, survival_prob + 0.03, f"{survival_prob*100:.1f}%", fontsize=10, color="black", ha="right")

risk_times = [0, 6, 12, 18, 24, 30, 36, 42]
risk_counts = [df_os_main[df_os_main["OSMONTHS"] >= t].shape[0] for t in risk_times]

for i, (t, r) in enumerate(zip(risk_times, risk_counts)):
    ax.text(t, -0.30, f"{r}", fontsize=12, color="black", ha="center")

ax.text(-2, -0.30, "Number at Risk", fontsize=12, ha="right")

#PNG形式でグラフをフォルダに出力
output_filename = "KaplanMeier_OS_main.png"
output_path = os.path.join(output_dir, output_filename)
plt.savefig(output_path, dpi=300, bbox_inches="tight")


# OS integrated analysis (matching 1:1)
# STUDY == 0 のデータ
df_os_integ_0 = df[df["STUDY"] == 0][["OSMONTHS", "OS_event"]].dropna()
df_os_integ_0.columns = ["OSMONTHS", "OS_event"]

# STUDY == 1 のデータ
df_os_integ_1 = df[df["STUDY"] == 1][["OSMONTHS", "OS_event"]].dropna()
df_os_integ_1.columns = ["OSMONTHS", "OS_event"]

# KaplanMeierFitter を使用
km_os_0 = KaplanMeierFitter()
km_os_1 = KaplanMeierFitter()

km_os_0.fit(df_os_integ_0["OSMONTHS"], event_observed=df_os_integ_0["OS_event"])
km_os_1.fit(df_os_integ_1["OSMONTHS"], event_observed=df_os_integ_1["OS_event"])

plt.rcParams["font.family"] = "Arial"
fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(bottom=0.4)

# STUDY == 0 生存曲線と95%信頼区間をプロット
km_os_0.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.1, linewidth=2, color="#CE1C48", label="Atezo+CDDP+VNR")

# STUDY == 1 生存曲線と95%信頼区間をプロット
km_os_1.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.1, linewidth=2, color="#1180D9", label="Matched CDDP+VNR")

# # STUDY == 0 の生存曲線（青）
# ax.step(np.insert(km_os_0.timeline, 0, 0), np.insert(km_os_0.survival_function_.values.flatten(), 0, 1.0), 
#         where="post", label="ADJUST", linewidth=2, color="blue")

# # STUDY == 1 の生存曲線（赤）
# ax.step(np.insert(km_os_1.timeline, 0, 0), np.insert(km_os_1.survival_function_.values.flatten(), 0, 1.0), 
#         where="post", label="IMPACT", linewidth=2, color="red")

# センサリングデータのプロット
censor_times_0 = df_os_integ_0["OSMONTHS"][df_os_integ_0["OS_event"] == 0]
censor_probs_0 = km_os_0.survival_function_at_times(censor_times_0).values
ax.scatter(censor_times_0, censor_probs_0, color='#CE1C48', marker='|', s=100)

censor_times_1 = df_os_integ_1["OSMONTHS"][df_os_integ_1["OS_event"] == 0]
censor_probs_1 = km_os_1.survival_function_at_times(censor_times_1).values
ax.scatter(censor_times_1, censor_probs_1, color='#1180D9', marker='|', s=100)

ax.set_xlabel("Months since Registration", fontsize=12, labelpad=8)
ax.set_ylabel("Probability of Overall Survival", fontsize=12)
ax.set_ylim(0, 1.05)

# 以前の横軸の最大値を使用
x_max = df_os_integ_0["OSMONTHS"].max() # 以前の最大値を使用
ax.set_xlim(0, x_max + 5)

ax.set_xticks(np.arange(0, x_max + 5, 12))
ax.set_yticks(np.arange(0, 1.1, 0.2))

ax.minorticks_on()
ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

# 12, 24, 36, 48ヶ月の補助線
for x in [12, 24, 36, 48]:
    ax.axvline(x=x, linestyle="--", color="gray", alpha=0.7)

# 12, 24, 36, 48ヶ月の生存割合をグラフに表示
for x in [12, 24, 36, 48]:
    if x in km_os_0.timeline:
        survival_prob_0 = km_os_0.survival_function_at_times([x]).values[0]
        ax.text(x, survival_prob_0 + 0.03, f"{survival_prob_0*100:.1f}%", fontsize=10, color="#CE1C48", ha="right")

    if x in km_os_1.timeline:
        survival_prob_1 = km_os_1.survival_function_at_times([x]).values[0]
        ax.text(x, survival_prob_1 + 0.03, f"{survival_prob_1*100:.1f}%", fontsize=10, color="#1180D9", ha="right")

# Number at Risk
risk_times = [0, 6, 12, 18, 24, 30, 36, 42, 48]
risk_counts_0 = [df_os_integ_0[df_os_integ_0["OSMONTHS"] >= t].shape[0] for t in risk_times]
risk_counts_1 = [df_os_integ_1[df_os_integ_1["OSMONTHS"] >= t].shape[0] for t in risk_times]

for i, (t, r0, r1) in enumerate(zip(risk_times, risk_counts_0, risk_counts_1)):
    ax.text(t, -0.35, f"{r0}", fontsize=12, color="#CE1C48", ha="center")
    ax.text(t, -0.40, f"{r1}", fontsize=12, color="#1180D9", ha="center")

ax.text(-2, -0.30, "Number at risk", fontsize=12, ha="right", color="black")
ax.text(-2, -0.35, "Atezo+CDDP+VNR", fontsize=12, ha="right", color="#CE1C48")
ax.text(-2, -0.40, "Matched CDDP+VNR", fontsize=12, ha="right", color="#1180D9")

ax.legend(fontsize=10)

#PNG形式でグラフをフォルダに出力
output_filename = "KaplanMeier_OS_inte.png"
output_path = os.path.join(output_dir, output_filename)
plt.savefig(output_path, dpi=300, bbox_inches="tight")


# OS integrated analysis (matching 1: 1)の生存割合
# 指定した時点の生存確率を取得
time_points = [12, 24, 36, 48]
survival_probabilities0 = km_os_0.survival_function_at_times(time_points)
survival_probabilities1 = km_os_1.survival_function_at_times(time_points)

# 生存中央値と95%信頼区間の計算
median_survival_0 = km_os_0.median_survival_time_
median_survival_1 = km_os_1.median_survival_time_

# 95%信頼区間の計算
ci_0 = km_os_0.confidence_interval_survival_function_
ci_1 = km_os_1.confidence_interval_survival_function_

# ログランク検定の実施
logrank_result = logrank_test(
    df_os_integ_0["OSMONTHS"], df_os_integ_1["OSMONTHS"],
    event_observed_A=df_os_integ_0["OS_event"], event_observed_B=df_os_integ_1["OS_event"]
)

# P値を取得
p_value = logrank_result.p_value

# ハザード比 (HR) の計算
hr = km_os_1.event_table.observed.sum() / km_os_0.event_table.observed.sum()

# 95%信頼区間の計算（カプランマイヤー推定法による近似）
ci_lower = np.exp(np.log(hr) - 1.96 * np.sqrt(1/km_os_1.event_table.observed.sum() + 1/km_os_0.event_table.observed.sum()))
ci_upper = np.exp(np.log(hr) + 1.96 * np.sqrt(1/km_os_1.event_table.observed.sum() + 1/km_os_0.event_table.observed.sum()))

# ✅ 保存フォルダとファイル名を指定
output_filename_txt = "KM_OS_probabilities.txt"
output_path_txt = os.path.join(output_dir, output_filename_txt)

# ✅ テキスト形式で保存
with open(output_path_txt, "w") as f:
    f.write("Survival probabilities\n")
    f.write("="*40 + "\n\n")
    f.write("Time (months)\tSurvival Probability\n")  # ヘッダー
    f.write("ADJUST study\n")
    for time, prob in zip(time_points, survival_probabilities0.values.flatten()):
        f.write(f"{time}\t{prob:.4f}\n")  # タブ区切りで保存
    f.write("IMPACT study\n")
    for time, prob in zip(time_points, survival_probabilities1.values.flatten()):
        f.write(f"{time}\t{prob:.4f}\n")  # タブ区切りで保存
    
    f.write("\nKaplan-Meier 生存解析結果\n")
    f.write("="*40 + "\n\n")

    # ✅ 生存中央値の書き込み
    f.write(f"生存中央値 (Median Survival Time):\n")
    f.write(f"  - ADJUST群 0: {median_survival_0:.2f} ヶ月\n")
    f.write(f"  - IMPACT群 1: {median_survival_1:.2f} ヶ月\n\n")

    # ✅ 95% 信頼区間の書き込み
    f.write("95% 信頼区間 (Confidence Interval):\n")
    f.write(f"  - ADJUST群 0: {ci_0.iloc[:, 0].min():.2f} - {ci_0.iloc[:, 0].max():.2f} ヶ月\n")
    f.write(f"  - IMPACT群 1: {ci_1.iloc[:, 0].min():.2f} - {ci_1.iloc[:, 0].max():.2f} ヶ月\n\n")

    # ✅ ログランク検定のP値の書き込み
    f.write(f"ログランク検定 (Log-rank test):\n")
    f.write(f"  - P値: {p_value:.4f}\n\n")

    # ✅ ハザード比 (HR) と 95% 信頼区間の書き込み
    f.write(f"ハザード比 (Hazard Ratio, HR):\n")
    f.write(f"  - HR: {hr:.2f}\n")
    f.write(f"  - 95% CI: {ci_lower:.2f} - {ci_upper:.2f}\n\n")

print(f"✅ txtデータを保存しました: {output_path_txt}")