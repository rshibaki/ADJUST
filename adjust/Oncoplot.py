import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.pandas2ri
from PIL import Image
import IPython.display as display
import os

# --- rpy2 による R との連携を有効化 ---
rpy2.robjects.pandas2ri.activate()

# ※ 保存先フォルダを指定（存在しない場合は作成）
save_folder = "/Users/rshibaki/Documents/project/ADJUST/Figure/"
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

# ============================================================
# 1. 変異データ（MAF 用）の読み込み・前処理
# ============================================================
file_path = "/Users/rshibaki/Library/CloudStorage/GoogleDrive-ryota.shibaki@gmail.com/マイドライブ/ADJUST試験/DATA/ADJUST_DNA.xlsx"
mutation_df = pd.read_excel(file_path, sheet_name="main_Selected")
print("元のデータ件数:", len(mutation_df))

# 必要なカラムを抽出（"Coding region change" も含む）
columns_needed = ['Patient_ID', 'Chromosome', 'Region', 'Type', 'Reference', 'Allele',
                  'Gene Cards', 'Non-synonymous', 'Coding region change']
mutation_df = mutation_df[columns_needed]

# Chromosome, Region は文字列に変換
mutation_df['Chromosome'] = mutation_df['Chromosome'].astype(str)
mutation_df['Region'] = mutation_df['Region'].astype(str)

# 変異位置の抽出
def extract_positions(region):
    try:
        if ".." in region:
            start, end = region.split("..")
        else:
            start = end = region
        start = ''.join(filter(str.isdigit, start))
        end = ''.join(filter(str.isdigit, end))
        return float(start) if start else np.nan, float(end) if end else np.nan
    except Exception as e:
        print(f"位置抽出エラー: {region} → {e}")
        return np.nan, np.nan

mutation_df[['Start_Position', 'End_Position']] = mutation_df['Region'].apply(lambda x: pd.Series(extract_positions(x)))

# 変異タイプの決定
def determine_variant_type(reference, allele):
    if isinstance(reference, str) and isinstance(allele, str):
        if len(reference) == 1 and len(allele) == 1:
            return "SNP"
        elif len(reference) > len(allele):
            return "DEL"
        elif len(reference) < len(allele):
            return "INS"
    return "OTHER"

mutation_df['Variant_Type'] = mutation_df.apply(
    lambda row: determine_variant_type(row['Reference'], row['Allele']),
    axis=1
)

# 変異分類の決定
def classify_variant(row):
    gene = row['Gene Cards']
    variant_type = row['Type']
    non_syn = str(row['Non-synonymous']).upper() if pd.notnull(row['Non-synonymous']) else ""
    coding_change = str(row['Coding region change']) if pd.notnull(row['Coding region change']) else ""
    
    if variant_type == "SNV" and (("-1" in coding_change) or ("-2" in coding_change) or ("splice" in coding_change.lower())):
        return "Splice_Site"
    if gene in ["WRN", "DICER1"]:
        return "Splice_Site" if variant_type == "SNV" else "Frame_Shift_Del"
    if non_syn in ['YES', 'TRUE']:
        return "Missense_Mutation" if variant_type == "SNV" else "Frame_Shift_Del"
    return "Missense_Mutation" if variant_type == "SNV" else "Frame_Shift_Del"

mutation_df['Variant_Classification'] = mutation_df.apply(classify_variant, axis=1)
mutation_df.fillna("Unknown", inplace=True)

# ※ "Missense_Mutation" は元の名称のままとする

# MAF フォーマットに変換
maf_df = mutation_df.rename(columns={
    'Patient_ID': 'Tumor_Sample_Barcode',
    'Gene Cards': 'Hugo_Symbol',
    'Reference': 'Reference_Allele',
    'Allele': 'Tumor_Seq_Allele2'
})
maf_df = maf_df[['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
                 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type']]

# ============================================================
# 2. TMB/DFS データの読み込みとサンプル順序の決定
# ============================================================
tmb_df = pd.read_excel(file_path, sheet_name="TMB")
tmb_df = tmb_df.drop_duplicates(subset='Patient_ID')
tmb_df['DFS'] = pd.to_numeric(tmb_df['DFS'], errors='coerce')
tmb_df_sorted = tmb_df.sort_values(by='DFS', ascending=False)

sample_order = tmb_df_sorted['Patient_ID'].astype(str).tolist()
print("サンプル順序（DFS 降順）:", sample_order)

# clinicalData 用に整形
clinical_df = tmb_df[['Patient_ID', 'TMB (Mutations/Mb)', 'DFS']].copy()
clinical_df = clinical_df.rename(columns={'Patient_ID': 'Tumor_Sample_Barcode', 'TMB (Mutations/Mb)': 'TMB'})
clinical_df['Tumor_Sample_Barcode'] = clinical_df['Tumor_Sample_Barcode'].astype(str)

# ============================================================
# 3. R (maftools) による Oncoplot 作成（TMB ヒートマップ削除版）
# ============================================================
maf_r_df = rpy2.robjects.pandas2ri.py2rpy(maf_df)
clinical_r_df = rpy2.robjects.pandas2ri.py2rpy(clinical_df)
sample_order_r = robjects.StrVector(sample_order)

r_code = f"""
library(maftools)

plot_oncoplot <- function(maf_data, clinical_data, sample_order) {{
    maf_df <- as.data.frame(maf_data)
    maf_object <- read.maf(maf = maf_df, clinicalData = clinical_data, vc_nonSyn = NULL)
    
    custom_cols <- c("Missense_Mutation" = "#589F99", 
                     "Splice_Site" = "#ADD7F9",
                     "Frame_Shift_Del" = "#F4E06F",
                     "Multi_Hit" = "#C8DB81",
                     "SNP" = "#E69F00", 
                     "DEL" = "#F0E442", 
                     "INS" = "#0072B2", 
                     "OTHER" = "#999999")
    
    output_path <- file.path("{save_folder}", "Oncoplot_main.png")
    png(filename = output_path, width = 1200, height = 800, bg="white")
    par(cex=5.0)
    oncoplot(maf = maf_object, 
             sampleOrder = sample_order,
             showTumorSampleBarcodes = FALSE,
             colors = custom_cols,
             top = 18)
    dev.off()
    print(paste("Oncoplot を保存しました:", output_path))
}}
"""
robjects.r(r_code)
plot_oncoplot = robjects.globalenv['plot_oncoplot']
_ = plot_oncoplot(maf_r_df, clinical_r_df, sample_order_r)
print("Oncoplot の描画が完了しました。")

# ============================================================
# 4. TMB の棒グラフ作成（症例順は oncoplot と同じ）
# ============================================================
r_code_tmb = f"""
plot_tmb_barplot <- function(clinical_data, sample_order) {{
    clinical_df <- as.data.frame(clinical_data)
    clinical_df <- clinical_df[match(sample_order, clinical_df$Tumor_Sample_Barcode), ]
    output_path <- file.path("{save_folder}", "Oncoplot_TMB.png")
    png(filename = output_path, width = 900, height = 120, bg="white")
    par(cex=0.8, mar=c(2, 4, 1, 2))  # 下部余白を確保
    # 棒グラフを描画し、棒の中央位置を取得
    bp <- barplot(clinical_df$TMB, 
                  names.arg = rep("", nrow(clinical_df)),   # 症例番号を非表示
                  las = 2, col = "#F4B1C2", border = NA, 
                  main = "", ylab = "TMB", 
                  ylim = c(0, 40))
    # y = 10 に破線を追加
    abline(h = 10, lty = 2)
    # 各棒の下に TMB の実数値を小数点第一位で表示するために sprintf() を使用
    text(x = bp, y = -2, labels = sprintf("%.1f", clinical_df$TMB), srt = 0, adj = c(0.5, 1), xpd = TRUE, cex = 1.5)
    axis(2, at = seq(0, 40, by = 10), las = 1)
    dev.off()
    print(paste("TMB barplot saved at:", output_path))
}}
"""
robjects.r(r_code_tmb)
plot_tmb_barplot = robjects.globalenv['plot_tmb_barplot']
_ = plot_tmb_barplot(clinical_r_df, sample_order_r)


# ============================================================
# 5. DFS の棒グラフ作成（症例順は oncoplot と同じ）
# ============================================================
r_code_dfs = f"""
plot_dfs_barplot <- function(clinical_data, sample_order) {{
    clinical_df <- as.data.frame(clinical_data)
    clinical_df <- clinical_df[match(sample_order, clinical_df$Tumor_Sample_Barcode), ]
    output_path <- file.path("{save_folder}", "Oncoplot_DFS.png")
    png(filename = output_path, width = 900, height = 120, bg="white")
    par(cex=0.8, mar=c(2, 4, 1, 2))  # 下部余白を確保
    bp <- barplot(clinical_df$DFS, 
                  names.arg = rep("", nrow(clinical_df)),   # 症例番号を非表示
                  las = 2, col = "#2D6DCC", border = NA, 
                  main = "", ylab = "DFS", 
                  ylim = c(0, 48),
                  yaxt = "n")
    axis(2, at = seq(0, 48, by = 12), las = 1)
    dev.off()
    print(paste("DFS barplot saved at:", output_path))
}}
"""
robjects.r(r_code_dfs)
plot_dfs_barplot = robjects.globalenv['plot_dfs_barplot']
_ = plot_dfs_barplot(clinical_r_df, sample_order_r)

