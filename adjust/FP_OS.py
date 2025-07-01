import pandas as pd
from rpy2.robjects import r, globalenv
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# pandas と R の DataFrame 間の自動変換を有効化
pandas2ri.activate()

# R 側のパッケージをインポート
survival   = importr('survival')
forestplot = importr('forestplot')

# 正しいファイルパスを指定
file_path = "/Users/rshibaki/Library/CloudStorage/GoogleDrive-ryota.shibaki@gmail.com/マイドライブ/ADJUST試験/DATA/dataset_250129.xlsx"

# データの読み込み
df = pd.read_excel(file_path, sheet_name="dataset")

# 警告が出ている列（例：Last_visit_DATE）の型を統一（ここでは日付型に変換）
df['Last_visit_DATE'] = pd.to_datetime(df['Last_visit_DATE'], errors='coerce')

# R 側に DataFrame を渡す（オブジェクト名：r_df）
globalenv['r_df'] = df

r('''
### --- 変数のラベル変換 ---
# SEX: 0 -> Male, 1 -> Female
r_df$SEX <- factor(r_df$SEX, levels = c(0,1), labels = c("Male", "Female"))
# Smoking status (smoke): 0 -> Never, 1 -> Current or former
r_df$smoke <- factor(r_df$smoke, levels = c(0,1), labels = c("Never", "Current or former"))
# EGFR mutation type (EGFR): 0 -> exon 19DEL, 1 -> L858R
r_df$EGFR <- factor(r_df$EGFR, levels = c(0,1), labels = c("exon 19DEL", "L858R"))

### ---------------------------
### PNG 出力の設定
### ---------------------------
output_dir = "/Users/rshibaki/Documents/project/ADJUST/Figure"
png(file.path(output_dir, "forestplot_OS.png"), width=800, height=600)
library(survival)
library(forestplot)
library(grid)

### ---------------------------
### 全症例の解析（All patients 行の作成）
### ---------------------------
model_all <- coxph(Surv(OSMONTHS, OS_event) ~ STUDY, data = r_df)
summ_all  <- summary(model_all)
if("STUDY1" %in% rownames(summ_all$coefficients)) {
    hr_all    <- exp(summ_all$coefficients["STUDY1", "coef"])
    lower_all <- summ_all$conf.int["STUDY1", "lower .95"]
    upper_all <- summ_all$conf.int["STUDY1", "upper .95"]
} else if ("STUDY" %in% rownames(summ_all$coefficients)) {
    hr_all    <- exp(summ_all$coefficients["STUDY", "coef"])
    lower_all <- summ_all$conf.int["STUDY", "lower .95"]
    upper_all <- summ_all$conf.int["STUDY", "upper .95"]
} else {
    hr_all    <- NA
    lower_all <- NA
    upper_all <- NA
}
n0_all <- sum(r_df$STUDY == "0")
n1_all <- sum(r_df$STUDY == "1")
total_all <- n0_all + n1_all
# 全症例では各グループの分母は自身なので 100%
count0_all <- paste0(n0_all, "/", n0_all, " (100%)")
count1_all <- paste0(n1_all, "/", n1_all, " (100%)")
all_row <- data.frame(Subgroup="All patients", Level="", HR=hr_all, lower=lower_all, upper=upper_all, 
                      Count0=count0_all, Count1=count1_all, N=total_all, isHeader=FALSE, isSummary=TRUE, stringsAsFactors = FALSE)

### ---------------------------
### データの前処理
### ---------------------------
# STUDY を因子型に変換（baseline が "0" になる）
r_df$STUDY <- as.factor(r_df$STUDY)
# Age を70歳以上／70歳未満の2群に分けるため、AgeGroup 変数を作成
r_df$AgeGroup <- ifelse(r_df$Age >= 70, ">=70", "<70")

### ---------------------------
### サブグループごとの解析
### ---------------------------
# サブグループ変数と表示用の名称の対応表（PS を ECOG PS に変更）
subgroup_map <- list(PS="ECOG PS", SEX="SEX", smoke="Smoking status", EGFR="EGFR mutation type", AgeGroup="Age")
results <- data.frame(Subgroup=character(), Level=character(), HR=numeric(), lower=numeric(), upper=numeric(), 
                      Count0=character(), Count1=character(), N=numeric(), isHeader=logical(), isSummary=logical(), stringsAsFactors=FALSE)

# 各サブグループについて、まずヘッダー行を追加し、その後各水準の解析結果を追加
for (subgroup in names(subgroup_map)) {
  header_label <- subgroup_map[[subgroup]]
  results <- rbind(results, data.frame(Subgroup=header_label, Level="", HR=NA, lower=NA, upper=NA, 
                                       Count0="", Count1="", N=NA, isHeader=TRUE, isSummary=FALSE, stringsAsFactors=FALSE))
  levels_sub <- sort(unique(r_df[[subgroup]]))
  for (lev in levels_sub) {
    subset_data <- r_df[r_df[[subgroup]] == lev, ]
    if(nrow(subset_data) > 0 && length(unique(subset_data$STUDY)) > 1) {
      model <- coxph(Surv(OSMONTHS, OS_event) ~ STUDY, data = subset_data)
      summ <- summary(model)
      if("STUDY1" %in% rownames(summ$coefficients)) {
         hr <- exp(summ$coefficients["STUDY1", "coef"])
         lower <- summ$conf.int["STUDY1", "lower .95"]
         upper <- summ$conf.int["STUDY1", "upper .95"]
      } else if ("STUDY" %in% rownames(summ$coefficients)) {
         hr <- exp(summ$coefficients["STUDY", "coef"])
         lower <- summ$conf.int["STUDY", "lower .95"]
         upper <- summ$conf.int["STUDY", "upper .95"]
      } else {
         next
      }
      n0 <- sum(subset_data$STUDY == "0")
      n1 <- sum(subset_data$STUDY == "1")
      total <- n0 + n1
      p0 <- round(n0 / n0_all * 100, 1)
      p1 <- round(n1 / n1_all * 100, 1)
      count0_str <- paste0(n0, "/", n0_all, " (", p0, "%)")
      count1_str <- paste0(n1, "/", n1_all, " (", p1, "%)")
      results <- rbind(results, data.frame(Subgroup=subgroup_map[[subgroup]], Level=as.character(lev), HR=hr, lower=lower, upper=upper, 
                                           Count0=count0_str, Count1=count1_str, N=total, isHeader=FALSE, isSummary=FALSE, stringsAsFactors=FALSE))
    }
  }
}
results <- rbind(all_row, results)

dev.off()
''')