# 0. Description --------------------------------------------------------------
# Author: [Author Name]
# Last Modified: 2026-04-27 (v03f)
# Goal: Examine the association between CT scan exposure in childhood and
# subsequent malignancy diagnosis using Taiwan's National Health
# Insurance Research Database (NHIRD), corresponding to three study
# designs:
#         Study 1 - Population-level cohort (all children aged 0-18)
#         Study 2 - Sibling/twin matched cohort
#         Study 3 - Appendicitis-restricted cohort
#
# 本版本（v03f）變更：對照 00-basic_rate.R 修正 arrow lazy query 問題。
# ───────────────────────────────────────────────────────────────────────
# 【F-1】§2.3 / §2.4 / §2.5：所有對 icd9cm_1 的 trimws/substr/條件運算
#        全部移到 collect() 之後。arrow lazy query 階段只保留
#        `select` + `filter(id %in% eligible_ids)`，對齊 00-basic_rate.R
#        的慣用法（該腳本所有衍生欄位計算都在 collect 後）。
#        症狀：原版會在某些 arrow 版本因 trimws/substr/if_else kernel
#        不支援或 NA 行為差異而報錯或回傳錯誤結果。
#
# 【F-2】§2.3 加入 ICD-9 平行查詢分支（2000–2015）。
#        原版只用 ICD-10 cancer_code (C00–C97) 對所有年份比對，會漏掉
#        2000–2015 期間 ICD-9-CM (140–208) 確診的惡性腫瘤紀錄，
#        導致 §3.3 prior_malignancy_ids 嚴重低估，
#        部分研究前已確診癌症者沒有被排除。
#
# 【F-3】§2.5 appendicitis_records_raw bind 後加 distinct(ID, FUNC_DATE,
#        SOURCE)，避免同人同次就醫多筆診斷在後續 inner_join 時 cartesian。
#
# 【F-4】§2.6 enrol_relation_ext 的 coalesce 改在 collect 後做，
#        避免 dictionary-encoded vs character 的 type mismatch。
#
# 【F-5】§2.7 ENROL 處理全面對齊 00-basic_rate.R 慣例：
#        (a) 加入 rename_with(tolower) 統一欄位名；
#        (b) arrow 階段只做 select + 簡單 filter，型別轉換移到 collect 後；
#        (c) 欄位名以小寫 id_status / id_roc 操作（對齊 03i 用 ID_STATUS
#            而非 STATUS 的觀察；rename_with(tolower) 後皆為 id_status）；
#        (d) join key (PREM_YM, ID1) 在 R 端統一型別後再 join。
#
# 【F-6】§2.7 ENROL 欄位名按官方手冊（H_NHI_ENROL）修正：
#        原版用了不存在的 INS_ID / INS_AMT 欄位，正確欄位是：
#          INS_ID  → ID1     （手冊序號 7：被保險人身分證字號）
#          INS_AMT → ID1_AMT （手冊序號 13：投保金額）
#        判斷主投保者：id == id1（手冊「注意事項 5」）。
#        同時移除 detect_status_col() 動態偵測，依手冊確認就是 id_status，
#        直接寫死。
#
# 【F-7】find_parquet() 加 suffix 參數（預設 "_10"）：
#        ENROL 月檔沒有 _10 sample group 後綴（檔名為 H_NHI_ENROL{yymm}），
#        原 find_parquet 寫死 _10 後綴會找不到 ENROL 檔。
#        §2.7 三處 find_parquet("ENROL", y, m) 改傳 suffix = ""。
#        OPDTE/OPDTO/IPDTE 不受影響（預設 _10 維持原行為）。
#
# 【F-8】§4.1 修正 ID_S 欄位撞名 bug：
#        base_cohort 已從 pers_info 衍生（已含 ID_S 與 BIRTH_DATE），原版
#        在 §4.1 又做 left_join(pers_info |> select(ID, ID_S))，導致兩邊
#        都有 ID_S → dplyr 自動加 suffix 變成 ID_S.x / ID_S.y → 後續
#        mutate(ID_S == "1") 找不到欄位而報錯。修正：刪除這個重複 join。
#
# 【F-9】§2.x 各建檔區段加 read_rds 開關（對齊 03i 風格）：
#        所有建檔迴圈預設註解掉，腳本起始就 read_rds 中間檔，讓二次跑
#        分析時可以一次跑完不必重建檔。要重新建檔請取消對應區塊註解。
#        例外：§2.5（appendicitis）保留實際執行，因為這個建檔目前在
#        debug 階段（rds 只 3KB，需要看 §2.5 的診斷 print 才能定位
#        哪一步抓 0 筆）。
#        §2.5 也加上完整的 nrow 診斷 print，方便快速看出斷點。
#
# ───────────────────────────────────────────────────────────────────────
# v03e (2026-04-27) 的歷史變更（保留）：
#   (C-1) §2.7 enrol-income-urban：改用「家戶主投保者」proxy
#   (C-2) §4 / §6 主分析 LM/Logistic + Cox 副分析
#   (C-3) §5.2 study2_control 新增 AGE_AT_INDEX ∈ [0, 18] filter
#   (C-5) §2.5 + §6.1 Study 3 暴露窗 ±7 天，主分析限 Group A
#   (I-1) §5.2 SEX_CONCORDANCE 改為 pair-level
#   (I-3) §1.4 闌尾炎 ICD-10 substr(..., 1, 3) 涵蓋 K35/K36/K37
#   (I-5) §4.2.1 / §6.2.1 分開呈現 FIRST_CT_AGE 與 AGE_AT_2016
#   (I-7) §5.2 study2_exposed select 加入 BIRTH_DATE
#   (M-*) 各種小重構（log_progress / find_parquet helper 等）
#
# ⚠ 中間檔結構變更（重跑時會覆寫）：
#   - malignancy-dx.rds：v03f 起涵蓋 ICD-9 (2000–2015) + ICD-10 (2016–2023)
#     兩段；CODE_VERSION 欄位標記紀錄來源。重跑會明顯增加 row 數。
#   - hereditary-records-full.rds：邏輯不變但 collect 後再做 substr，結果
#     應與舊版一致；若有差異代表舊版有 arrow kernel bug。
#   - enrol-income-urban.rds：F-6 修正後內容會大幅改變（從找不到任何
#     主投保者 → 正確抓到 ID1 / ID1_AMT 後 INCOME_Q 涵蓋率應大幅提升）。
#   - enrol-relation-ext.rds：v03f 新增（從原本只在記憶體中改為持久化）。
#   - 其他中間檔不變。

# 1. Libraries and Configs ----------------------------------------------------

## 1.1 Load Libraries ---------------------------------------------------------

library(tidyverse)
library(arrow)
library(rlang)
library(fixest)
library(broom)
library(scales)
library(survival)
library(lubridate)

## 1.2 Constants --------------------------------------------------------------

DATA_PATH           <- "../../data/parquet"
PROCESSED_DATA_PATH <- "../../CleanData/processed-data"
INTERMEDIATE_PATH   <- "data/intermediate"
OUTPUT_TABLE_PATH   <- file.path("outputs", "tables")
OUTPUT_FIGURE_PATH  <- file.path("outputs", "figures")

# 確保資料夾存在
dir.create(INTERMEDIATE_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TABLE_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIGURE_PATH, recursive = TRUE, showWarnings = FALSE)

# Study window
STUDY_START      <- ymd(20160101)
STUDY_END        <- ymd(20231231)
MAX_AGE_AT_INDEX <- 18L
LATENCY_YEARS    <- 2L

# 承保檔篩選閾值（對齊 00-basic_rate.R 與 03i 慣例）
AMT_CUTOFF <- 2000L

options(scipen = 123)

## 1.3 CT Procedure Codes in NHIRD --------------------------------------------
# 頭部型電腦斷層（33067B - 33069B）：
# 33067B 頭部型電腦斷層造影 — 無造影劑
# 33068B 頭部型電腦斷層造影 — 有造影劑
# 33069B 頭部型電腦斷層造影 — 有/無造影劑
# 一般/全身型電腦斷層（33070B - 33072B）：
# 33070B 電腦斷層造影 — 無造影劑
# 33071B 電腦斷層造影 — 有造影劑
# 33072B 電腦斷層造影 — 有/無造影劑
# ⚠ 健保 CT 代碼無部位資訊；Study 3 以非頭部型代碼作為腹部 CT proxy

all_ct_codes           <- c("33070B", "33071B", "33072B", "33067B", "33068B", "33069B")
abdomen_ct_proxy_codes <- c("33070B", "33071B", "33072B")

## 1.4 Malignancy ICD Codes ---------------------------------------------------
# NHIRD 跨 ICD 版本：
#   2000–2015：ICD-9-CM   惡性腫瘤主碼為 140–208（不含 209x 神經內分泌、
#                          210–229 良性，亦不含 230 之後 in situ）
#   2016–2023：ICD-10-CM  惡性腫瘤為 C00–C97
# ⚠ icd9cm_1 欄位左靠右補空白，需在 R 端先 trimws() 再取前 3 碼比對
# ⚠ F-2 修正：原版只比對 C00–C97 對所有年份查詢，ICD-9 段（2000–2015）永遠
#    匹配 0 筆，會導致 §3.3 prior_malignancy_ids 漏掉研究前已確診者。
#    本版分段：年份 ≤ 2015 用 ICD-9 cancer_code_icd9；年份 ≥ 2016 用
#    cancer_code_icd10。
# ⚠ 為簡化處理，ICD-9 cancer code 同樣使用 3 碼前綴比對（140–208）。
#    這涵蓋 ICD-9-CM 標準分類中的 malignant neoplasm 主章。

cancer_code_icd10 <- paste0("C", sprintf("%02d", 0:97))
cancer_code_icd9  <- as.character(140:208)

# 為向後相容（描述用）保留舊變數名 cancer_code = ICD-10
cancer_code <- cancer_code_icd10

## 1.4.1 Appendicitis ICD-10 Codes (I-3 修正) --------------------------------
# 闌尾炎相關 ICD-10-CM 三碼前綴：
#   K35  Acute appendicitis（包含 K35.2x、K35.3x、K35.80–K35.891 等所有子碼）
#   K36  Other appendicitis
#   K37  Unspecified appendicitis
# 原版本僅列 K358/K359/K370–K379 4 碼前綴會漏掉 K352/K353/K36 等 ICD-10-CM
# 標準子碼，故改為 3 碼前綴 substr(..., 1, 3)，覆蓋更廣。

appendicitis_icd_prefix <- c("K35", "K36", "K37")

## 1.5 六都定義（ID1_CITY 前兩碼判別，依醫療機構現況檔代碼簿）----------------
# 六都：
#   01 = 台北市、03 = 台中市（99年起）、05 = 台南市（99年起）、
#   07 = 高雄市（99年起）、31 = 新北市、32 = 桃園市（104年起）
# 99 年以前高雄/台南/台中以不同代碼出現，但研究期間為 2016-2023（民國 105-112），
# 完全落在 104 年後代碼（含 32xx 桃園市）的區間，毋須處理舊版縣市代碼。
# ⚠ 若資料實際落在 99–103 年間，桃園縣仍為 32xx；故 c("01","03","05","07","31","32")
#   對 2016 年後一致成立。

SIX_CITIES_PREFIX <- c("01", "03", "05", "07", "31", "32")

classify_urban <- function(id1_city) {
  # id1_city 為 4 碼字串（例：0101 = 台北市松山區）
  # 取前 2 碼對照六都前綴
  # M-2 修正：先強制轉 character、空字串視為 NA，避免 NA literal 進來時誤判
  id1_city <- as.character(id1_city)
  id1_city[id1_city == "" | id1_city == "NA"] <- NA_character_
  prefix <- substr(id1_city, 1, 2)
  dplyr::case_when(
    is.na(id1_city)               ~ NA_character_,
    prefix %in% SIX_CITIES_PREFIX ~ "Metro",
    TRUE                          ~ "Non-Metro"
  )
}

## 1.6 Helpers ----------------------------------------------------------------

calc_age_years <- function(birth_date, ref_date) {
  floor(time_length(interval(birth_date, ref_date), "year"))
}

# BIRTH_DATE 解析：兼容 pers_info$ID_BIRTHYM 的不同存儲型別
# ⚠ 上游檔案在不同版本中 ID_BIRTHYM 可能是：
#   (a) Date（如 2016-05-01，最常見；03i 即是此型別）
#   (b) character 六碼（如 "201605"）
#   (c) integer 六碼（如 201605L）
# 若直接 paste0(ID_BIRTHYM, "15") 在 (a) 情形會產出 "2016-05-0115"，
# ymd() 全部 parse 失敗 → NA → 整個 eligible_ids 為 0。
# 本 helper 統一三種型別都回傳合法 Date（月中 15 號，避開月初 edge case）。
resolve_birth_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }
  # character / numeric / integer 六碼 YYYYMM → 拼 "15" 後 parse
  ymd(paste0(as.character(x), "15"), quiet = TRUE)
}

# 6-key join vector（OPDTO join OPDTE 共用）
join_keys <- c("fee_ym", "appl_date", "appl_type", "case_type", "seq_no", "hosp_id")

# Join key 標準化：統一型別、去空白
normalise_join_keys <- function(dt, keys = join_keys) {
  dt |> mutate(across(all_of(keys), ~ toupper(trimws(as.character(.)))))
}

# 進度顯示 helper：印出時間戳與訊息
log_progress <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

# M-7 統一處理 H_NHI_ 前綴 fallback
# 月檔：name = "OPDTE", y, m → H_NHI_OPDTE{roc_yr}{mm}_10.parquet（預設 suffix）
#       name = "ENROL", y, m, suffix = "" → H_NHI_ENROL{roc_yr}{mm}.parquet
# 年檔：name = "IPDTE", y     → H_NHI_IPDTE{roc_yr}.parquet
# 若帶 H_NHI_ 前綴的檔案不存在，嘗試不帶前綴；皆不在則回 NA_character_
# ⚠ ENROL 月檔沒有 _10 sample group 後綴，呼叫時必須指定 suffix = ""
find_parquet <- function(name, y, m = NULL, suffix = "_10") {
  base <- if (is.null(m)) {
    sprintf("%s%d.parquet", name, y - 1911)
  } else {
    sprintf("%s%d%02d%s.parquet", name, y - 1911, m, suffix)
  }
  for (prefix in c("H_NHI_", "")) {
    p <- file.path(DATA_PATH, paste0(prefix, base))
    if (file.exists(p)) return(p)
  }
  NA_character_
}

# 分層 RR 計算函式
stratify_rr <- function(data, strat_var) {
  data |>
    group_by({{ strat_var }}, EXPOSED) |>
    summarise(
      N_EVENT = sum(OUTCOME),
      N_TOTAL = n(),
      RISK    = N_EVENT / N_TOTAL,
      .groups = "drop"
    ) |>
    group_by({{ strat_var }}) |>
    summarise(
      RR = RISK[EXPOSED] / RISK[!EXPOSED],
      .groups = "drop"
    )
}

# 2. Load Intermediate Data ---------------------------------------------------
# ⚠ 下方「建檔程式碼區塊」會實際執行（首次跑完後可考慮拆成獨立 00-build 腳本）。
#   每個建檔區塊結尾都會 write_rds() 到 INTERMEDIATE_PATH，
#   之後重跑分析時可直接從 read_rds() 開始，免得重新讀 parquet。

## 2.1 Personal Information ---------------------------------------------------

log_progress("讀取 pers-info.parquet …")

pers_info <- read_parquet(
  file.path(PROCESSED_DATA_PATH, "pers-info.parquet")
)

# 預先計算符合年齡條件的 eligible ID 清單
# 研究期間 2016-2023，納入條件：0-18 歲
# → 最晚出生：2023-12-31（AGE_2023 >= 0）  → 出生年 ≤ 2023
# → 最早出生：2016-01-01 時年齡 <= 18      → 出生年 ≥ 1997
# 因此 eligible_ids 的出生年範圍為 1997–2023
# 這批人在民國89年（西元2000年）時最大為3歲，均已存在於 NHIRD 中，
# 故可安全地從 2000 年起回查惡性腫瘤及遺傳疾患歷史紀錄。
#
# ⚠ 判斷方式（2026-04-21 更新）：
#   改用 "研究期間內任一天年齡 0-18" 的 windowed 判準：
#     BIRTH_DATE <= STUDY_END                              → 研究期內已出生
#     BIRTH_DATE + years(MAX_AGE_AT_INDEX + 1) > STUDY_START → 研究期間尚未滿 19 歲
#   此寫法與 calc_age_years(ymd-month-mid, ...) 不同之處：
#     直接用 BIRTH_DATE + years(n) 避免 1997 年 1 月出生者被 floor 掉成 19 歲
#     邊界 cohort（尤其 1997 年一整年）不會被錯殺。

HISTORY_START_YEAR <- 2000L   # 民國89年，惡性腫瘤/遺傳疾患歷史回查起始年

pers_info_with_birth <- pers_info |>
  mutate(BIRTH_DATE = resolve_birth_date(ID_BIRTHYM))

n_birth_na <- sum(is.na(pers_info_with_birth$BIRTH_DATE))
log_progress(sprintf("pers_info 總人數 = %s；BIRTH_DATE 解析失敗 = %s",
                     format(nrow(pers_info_with_birth), big.mark = ","),
                     format(n_birth_na, big.mark = ",")))
if (n_birth_na == nrow(pers_info_with_birth)) {
  stop("所有 BIRTH_DATE 都是 NA — 請檢查 pers_info$ID_BIRTHYM 的型別與內容")
}

eligible_ids <- pers_info_with_birth |>
  filter(
    !is.na(BIRTH_DATE),
    BIRTH_DATE <= STUDY_END,
    BIRTH_DATE + years(MAX_AGE_AT_INDEX + 1L) > STUDY_START
  ) |>
  pull(ID)

log_progress(sprintf("eligible_ids 筆數 = %s", format(length(eligible_ids), big.mark = ",")))
if (length(eligible_ids) == 0L) {
  stop("eligible_ids 為 0 — 請檢查 pers_info$ID_BIRTHYM 的內容")
}

## 2.2 CT Scan Records (OPDTO only; IPDTO 不可用) -----------------------------
# 【建檔步驟:對齊 03i 風格已註解】
# ⚠ 資料來源:
#   - 門診:OPDTO(含 drug_no)join OPDTE(含 id, func_date)
#   - 本專案無 IPDTO,不抓住院 CT 醫令
#   - OPDTE 月份格式:H_NHI_OPDTE{roc_yr}{mm}_10.parquet
# ⚠ 二次跑分析時直接 read_rds 即可;要重建檔請取消下方註解區塊

# log_progress("===== 建檔 2.2:門診 CT 紀錄(OPDTO join OPDTE)=====")
#
# N_OP <- length(2016:2023) * 12L   # 96 個 slot
# ct_op_list <- vector("list", N_OP)
# idx <- 1L
# for (y in 2016:2023) {
#   log_progress(sprintf("  [CT/OP] year %d", y))
#   for (m in 1:12) {
#     file_o <- find_parquet("OPDTO", y, m)
#     file_e <- find_parquet("OPDTE", y, m)
#
#     if (is.na(file_o) || is.na(file_e)) {
#       idx <- idx + 1L
#       next
#     }
#     # 先篩醫令(小檔),確認有 CT 才讀就醫紀錄(大檔)
#     dt_o <- open_dataset(file_o) |>
#       rename_with(tolower) |>
#       select(fee_ym, appl_date, appl_type, case_type, seq_no, hosp_id, drug_no) |>
#       filter(drug_no %in% all_ct_codes) |>
#       collect() |>
#       rename(order_code = drug_no) |>
#       normalise_join_keys()
#     if (nrow(dt_o) == 0L) {
#       idx <- idx + 1L
#       next
#     }
#     dt_e <- open_dataset(file_e) |>
#       rename_with(tolower) |>
#       select(fee_ym, appl_date, appl_type, case_type, seq_no, hosp_id, id, func_date) |>
#       filter(id %in% eligible_ids) |>
#       collect() |>
#       normalise_join_keys()
#     ct_op_list[[idx]] <- inner_join(dt_e, dt_o, by = join_keys)
#     idx <- idx + 1L
#   }
# }
#
# ct_records_raw <- bind_rows(ct_op_list) |>
#   mutate(
#     func_date = ymd(func_date),
#     source    = "op",
#     ct_type   = if_else(
#       order_code %in% c("33067B", "33068B", "33069B"),
#       "head_type", "body_type"
#     )
#   ) |>
#   select(id, func_date, order_code, source, ct_type) |>
#   rename(ID = id, FUNC_DATE = func_date, ORDER_CODE = order_code,
#          SOURCE = source, CT_TYPE = ct_type)
#
# log_progress(sprintf("CT 紀錄合併完成:%s 筆", format(nrow(ct_records_raw), big.mark = ",")))
# write_rds(ct_records_raw, file.path(INTERMEDIATE_PATH, "ct-records-all.rds"))
# log_progress("已儲存 ct-records-all.rds")

ct_records_raw <- read_rds(file.path(INTERMEDIATE_PATH, "ct-records-all.rds"))
ct_records <- ct_records_raw
log_progress(sprintf("讀入 ct-records-all.rds:%s 筆", format(nrow(ct_records), big.mark = ",")))

## 2.3 Malignancy Diagnosis ---------------------------------------------------
# 【建檔步驟:對齊 03i 風格已註解】
# ⚠ 資料來源:
#   - 診斷碼(icd9cm_1)存在 OPDTE(門診)與 IPDTE(住院)
#   - 2000–2015 ICD-9-CM (140-208);2016–2023 ICD-10-CM (C00-C97)
#   - icd9cm_1 左靠右補空白,在 R 端 collect 後再 trimws + substr
# ⚠ 二次跑分析時直接 read_rds 即可;要重建檔請取消下方註解區塊
#
# log_progress("===== 建檔 2.3:惡性腫瘤診斷(OPDTE + IPDTE)=====")
#
# ICD9_CUTOFF_YEAR <- 2015L
#
# # Helper:給定一筆 collected df 與年份,回傳 cancer dx 紀錄
# # 在 R 端做 trimws + substr 比對(不在 arrow 上)
# filter_cancer <- function(df, year) {
#   if (nrow(df) == 0L) return(df[0, , drop = FALSE])
#   code_set <- if (year <= ICD9_CUTOFF_YEAR) cancer_code_icd9 else cancer_code_icd10
#   df |>
#     mutate(icd3 = substr(trimws(icd9cm_1), 1, 3)) |>
#     filter(icd3 %in% code_set)
# }
#
# N_MAL_OP <- length(HISTORY_START_YEAR:2023) * 12L
# malignancy_op_list <- vector("list", N_MAL_OP)
# idx <- 1L
# for (y in HISTORY_START_YEAR:2023) {
#   log_progress(sprintf("  [Malignancy/OP] year %d", y))
#   for (m in 1:12) {
#     file_e <- find_parquet("OPDTE", y, m)
#     if (is.na(file_e)) {
#       idx <- idx + 1L
#       next
#     }
#     # F-1:arrow 階段只做 select + 簡單 filter,不做任何字串運算
#     df_raw <- open_dataset(file_e) |>
#       rename_with(tolower) |>
#       select(id, func_date, icd9cm_1) |>
#       filter(id %in% eligible_ids) |>
#       collect()
#     # 在 R 端做 ICD trimws/substr/比對與日期 parse
#     df_cancer <- filter_cancer(df_raw, y)
#     if (nrow(df_cancer) > 0L) {
#       malignancy_op_list[[idx]] <- df_cancer |>
#         mutate(
#           func_date     = ymd(func_date),
#           source        = "op",
#           code_version  = if_else(y <= ICD9_CUTOFF_YEAR, "ICD-9", "ICD-10")
#         ) |>
#         select(id, func_date, icd3, source, code_version)
#     }
#     idx <- idx + 1L
#   }
# }
#
# N_MAL_IP <- length(HISTORY_START_YEAR:2023)
# malignancy_ip_list <- vector("list", N_MAL_IP)
# idx <- 1L
# for (y in HISTORY_START_YEAR:2023) {
#   file_e <- find_parquet("IPDTE", y)
#   log_progress(sprintf("  [Malignancy/IP] year %d", y))
#   if (is.na(file_e)) {
#     idx <- idx + 1L
#     next
#   }
#   df_raw <- open_dataset(file_e) |>
#     rename_with(tolower) |>
#     select(id, func_date = in_date, icd9cm_1) |>
#     filter(id %in% eligible_ids) |>
#     collect()
#   df_cancer <- filter_cancer(df_raw, y)
#   if (nrow(df_cancer) > 0L) {
#     malignancy_ip_list[[idx]] <- df_cancer |>
#       mutate(
#         func_date    = ymd(func_date),
#         source       = "ip",
#         code_version = if_else(y <= ICD9_CUTOFF_YEAR, "ICD-9", "ICD-10")
#       ) |>
#       select(id, func_date, icd3, source, code_version)
#   }
#   idx <- idx + 1L
# }
#
# malignancy_raw <- bind_rows(
#   bind_rows(malignancy_op_list),
#   bind_rows(malignancy_ip_list)
# ) |>
#   distinct(id, func_date, icd3, code_version) |>
#   group_by(id) |>
#   arrange(func_date, icd3, .by_group = TRUE) |>
#   slice(1) |>
#   ungroup() |>
#   rename(
#     ID                    = id,
#     FIRST_MALIGNANCY_DATE = func_date,
#     ICD3                  = icd3,
#     CODE_VERSION          = code_version
#   ) |>
#   mutate(
#     icd_num = suppressWarnings(as.integer(substr(ICD3, 2, 3))),
#     MALIGNANCY_TYPE = case_when(
#       CODE_VERSION == "ICD-9"     ~ "icd9_unclassified",
#       icd_num <= 14               ~ "head_and_neck",
#       icd_num <= 26               ~ "intestinal",
#       icd_num <= 39               ~ "chest_and_lung",
#       icd_num <= 41               ~ "bone",
#       icd_num <= 49               ~ "cnt",
#       icd_num <= 58               ~ "breast_and_female",
#       icd_num <= 68               ~ "urinary_and_fertile",
#       icd_num <= 72               ~ "brain_and_cns",
#       icd_num <= 75               ~ "endocrine",
#       icd_num <= 80               ~ "other_malignant_tumor",
#       icd_num <= 97               ~ "lymphoma_and_leukemia",
#       TRUE                        ~ NA_character_
#     )
#   ) |>
#   select(-icd_num)
#
# log_progress(sprintf("惡性腫瘤首次診斷:%s 人(含 ICD-9 段 %s 人,ICD-10 段 %s 人)",
#                      format(nrow(malignancy_raw), big.mark = ","),
#                      format(sum(malignancy_raw$CODE_VERSION == "ICD-9"), big.mark = ","),
#                      format(sum(malignancy_raw$CODE_VERSION == "ICD-10"), big.mark = ",")))
# write_rds(malignancy_raw, file.path(INTERMEDIATE_PATH, "malignancy-dx.rds"))
# log_progress("已儲存 malignancy-dx.rds")

malignancy_raw <- read_rds(file.path(INTERMEDIATE_PATH, "malignancy-dx.rds"))
malignancy_dx  <- malignancy_raw
log_progress(sprintf("讀入 malignancy-dx.rds:%s 人", format(nrow(malignancy_dx), big.mark = ",")))

## 2.4 Hereditary Cancer Exclusion --------------------------------------------
# 【建檔步驟:對齊 03i 風格已註解】
# ⚠ 資料來源:OPDTE / IPDTE 的 icd9cm_1
#   - 2000–2015:ICD-9-CM;2016–2023:ICD-10-CM
#   - icd9cm_1 欄位名稱不因版本而改變,但內容格式不同
# ⚠ 二次跑分析時直接 read_rds 即可;要重建檔請取消下方註解區塊
#
# log_progress("===== 建檔 2.4:遺傳性癌症疾患排除清單(OPDTE + IPDTE)=====")
#
# hereditary_icd9 <- c(
#   "23770", "23771", "23772", "23779",
#   "7595",
#   "7580",
#   "V8401", "V8402", "V8409"
# )
#
# hereditary_icd10_5 <- c(
#   "Q8500", "Q8501", "Q8502", "Q8503", "Q8509",
#   "Q900",  "Q901",  "Q902",  "Q909",
#   "Z1501", "Z1502", "Z1509"
# )
#
# hereditary_icd10_4 <- c("Q851")
#
# hereditary_code_labels <- bind_rows(
#   tibble(
#     CODE    = hereditary_icd9,
#     VERSION = "ICD-9",
#     LABEL   = c("NF 未分類", "NF type 1(von Recklinghausen)", "NF type 2",
#                 "Schwannomatosis / 其他NF",
#                 "結節性硬化症(Tuberous sclerosis)",
#                 "唐氏症(Down syndrome)",
#                 "BRCA1 遺傳易感性", "BRCA2 遺傳易感性", "其他遺傳性癌症易感性")
#   ),
#   tibble(
#     CODE    = c(hereditary_icd10_5, hereditary_icd10_4),
#     VERSION = "ICD-10",
#     LABEL   = c("NF 未分類", "NF type 1", "NF type 2", "Schwannomatosis", "其他NF",
#                 "唐氏症 trisomy 21", "唐氏症 mosaic", "唐氏症 translocation", "唐氏症 未分類",
#                 "BRCA1 遺傳易感性", "BRCA2 遺傳易感性", "其他遺傳性癌症易感性",
#                 "結節性硬化症(Tuberous sclerosis)")
#   )
# )
#
# # F-1:把 hereditary 比對邏輯抽成 R 端 helper(不在 arrow 上做任何字串運算)
# filter_hereditary <- function(df, year) {
#   if (nrow(df) == 0L) return(df[0, c("id", "func_date", "HERED_CODE"), drop = FALSE])
#   df <- df |> mutate(icd_raw = trimws(icd9cm_1))
#   if (year <= ICD9_CUTOFF_YEAR) {
#     df |>
#       filter(icd_raw %in% hereditary_icd9) |>
#       mutate(HERED_CODE = icd_raw) |>
#       select(id, func_date, HERED_CODE)
#   } else {
#     df |>
#       mutate(
#         icd5 = substr(icd_raw, 1, 5),
#         icd4 = substr(icd_raw, 1, 4)
#       ) |>
#       filter(icd5 %in% hereditary_icd10_5 | icd4 %in% hereditary_icd10_4) |>
#       mutate(
#         HERED_CODE = if_else(icd5 %in% hereditary_icd10_5, icd5, icd4)
#       ) |>
#       select(id, func_date, HERED_CODE)
#   }
# }
#
# N_HER_OP <- length(HISTORY_START_YEAR:2023) * 12L
# hereditary_op_list <- vector("list", N_HER_OP)
# idx <- 1L
# for (y in HISTORY_START_YEAR:2023) {
#   log_progress(sprintf("  [Hered/OP] year %d", y))
#   for (m in 1:12) {
#     file_e <- find_parquet("OPDTE", y, m)
#     if (is.na(file_e)) {
#       idx <- idx + 1L
#       next
#     }
#     df_raw <- open_dataset(file_e) |>
#       rename_with(tolower) |>
#       select(id, func_date, icd9cm_1) |>
#       filter(id %in% eligible_ids) |>
#       collect()
#     df_hered <- filter_hereditary(df_raw, y)
#     if (nrow(df_hered) > 0L) {
#       hereditary_op_list[[idx]] <- df_hered |>
#         mutate(
#           func_date = ymd(func_date),
#           source    = "op"
#         )
#     }
#     idx <- idx + 1L
#   }
# }
#
# N_HER_IP <- length(HISTORY_START_YEAR:2023)
# hereditary_ip_list <- vector("list", N_HER_IP)
# idx <- 1L
# for (y in HISTORY_START_YEAR:2023) {
#   file_e <- find_parquet("IPDTE", y)
#   log_progress(sprintf("  [Hered/IP] year %d", y))
#   if (is.na(file_e)) {
#     idx <- idx + 1L
#     next
#   }
#   df_raw <- open_dataset(file_e) |>
#     rename_with(tolower) |>
#     select(id, func_date = in_date, icd9cm_1) |>
#     filter(id %in% eligible_ids) |>
#     collect()
#   df_hered <- filter_hereditary(df_raw, y)
#   if (nrow(df_hered) > 0L) {
#     hereditary_ip_list[[idx]] <- df_hered |>
#       mutate(
#         func_date = ymd(func_date),
#         source    = "ip"
#       )
#   }
#   idx <- idx + 1L
# }
#
# hereditary_records_raw <- bind_rows(
#   bind_rows(hereditary_op_list),
#   bind_rows(hereditary_ip_list)
# ) |>
#   rename(ID = id, HERED_DATE = func_date)
#
# hereditary_exclusion_raw <- hereditary_records_raw |>
#   group_by(ID) |>
#   arrange(HERED_DATE, HERED_CODE, .by_group = TRUE) |>
#   slice(1) |>
#   ungroup() |>
#   select(ID, FIRST_HERED_DATE = HERED_DATE, FIRST_HERED_CODE = HERED_CODE)
#
# log_progress(sprintf("遺傳疾患排除人數:%s", format(nrow(hereditary_exclusion_raw), big.mark = ",")))
#
# write_rds(hereditary_records_raw,
#           file.path(INTERMEDIATE_PATH, "hereditary-records-full.rds"))
# write_rds(hereditary_exclusion_raw,
#           file.path(INTERMEDIATE_PATH, "hereditary-exclusion-ids.rds"))
# log_progress("已儲存 hereditary-records-full.rds 與 hereditary-exclusion-ids.rds")

hereditary_records_raw   <- read_rds(file.path(INTERMEDIATE_PATH, "hereditary-records-full.rds"))
hereditary_exclusion_raw <- read_rds(file.path(INTERMEDIATE_PATH, "hereditary-exclusion-ids.rds"))
hereditary_exclusion <- hereditary_exclusion_raw |> select(ID)
log_progress(sprintf("讀入 hereditary-exclusion-ids.rds:%s 人",
                     format(nrow(hereditary_exclusion), big.mark = ",")))


## 2.5 Appendicitis + Appendectomy Records (Study 3) --------------------------
# 【建檔步驟：實際執行】
# ⚠ 資料來源：
#   - 闌尾炎診斷（icd9cm_1）：OPDTE + IPDTE
#   - 闌尾切除術醫令：僅 OPDTO（本專案無 IPDTO）
#   - 闌尾切除術代碼：74002B、74004B
#
# ⚠ 無 IPDTO 的處理方式（C-5 修正版）：
#   由於 IPDTO 不可用，住院中進行的闌尾切除術無法從醫令端確認。
#   臨床上「因闌尾炎而住院的病人 ≈ 將接受闌尾切除術」雖為合理假設，
#   但此處無法用醫令驗證；此外，住院期間若有 CT 也不會出現在 OPDTO，
#   會嚴重低估 Group B 的 CT 暴露率。
#
#   本版本 appendicitis_appendectomy 的定義：
#     (A) OPDTE 有闌尾炎診斷 AND ±7 天內 OPDTO 有闌尾切除術醫令 → 門診確認組
#     (B) IPDTE 有闌尾炎住院診斷                                → 住院 proxy 組
#
#   【主要分析】Study 3 主分析僅使用 (A)（APPX_SOURCE = "op_op_medorder"）
#   【敏感性分析】(B) 改為敏感性分析（APPX_SOURCE = "ip_diagnosis_proxy"），
#   並在報告中註明：因 IPDTO 缺失，此組 CT 暴露率必然低估。
#
#   日後取得 IPDTO 後，可把 (B) 改為同 (A) 的嚴謹醫令確認邏輯，併入主分析。

log_progress("===== 建檔 2.5：闌尾炎 + 闌尾切除術（僅 OPDTO）=====")

# I-3：闌尾炎 ICD-10 已於 §1.4.1 統一定義為 appendicitis_icd_prefix
# 此處沿用 3 碼前綴比對（substr(..., 1, 3)）
appendectomy_codes <- c("74002B", "74004B")

# -- 闌尾炎診斷（從 OPDTE / IPDTE 取 icd9cm_1）---------------------------------
# F-1：對齊 00-basic_rate.R，arrow 階段只 select + filter(id)，
# substr/trimws 移到 collect 後在 R 端執行

# Helper：在 R 端做 ICD-10 闌尾炎 3 碼前綴比對
filter_appendicitis <- function(df) {
  if (nrow(df) == 0L) return(df[0, , drop = FALSE])
  df |>
    mutate(icd3 = substr(trimws(icd9cm_1), 1, 3)) |>
    filter(icd3 %in% appendicitis_icd_prefix)
}

N_APP_OP <- length(2016:2023) * 12L
appendicitis_op_list <- vector("list", N_APP_OP)
idx <- 1L
for (y in 2016:2023) {
  log_progress(sprintf("  [Appendicitis/OP] year %d", y))
  for (m in 1:12) {
    file_e <- find_parquet("OPDTE", y, m)
    
    if (is.na(file_e)) {
      idx <- idx + 1L
      next
    }
    df_raw <- open_dataset(file_e) |>
      rename_with(tolower) |>
      select(id, func_date, icd9cm_1) |>
      filter(id %in% eligible_ids) |>
      collect()
    df_app <- filter_appendicitis(df_raw)
    if (nrow(df_app) > 0L) {
      appendicitis_op_list[[idx]] <- df_app |>
        mutate(func_date = ymd(func_date), source = "op")
    }
    idx <- idx + 1L
  }
}

N_APP_IP <- length(2016:2023)
appendicitis_ip_list <- vector("list", N_APP_IP)
idx <- 1L
for (y in 2016:2023) {
  file_e <- find_parquet("IPDTE", y)
  log_progress(sprintf("  [Appendicitis/IP] year %d", y))
  
  if (is.na(file_e)) {
    idx <- idx + 1L
    next
  }
  df_raw <- open_dataset(file_e) |>
    rename_with(tolower) |>
    select(id, func_date = in_date, icd9cm_1) |>
    filter(id %in% eligible_ids) |>
    collect()
  df_app <- filter_appendicitis(df_raw)
  if (nrow(df_app) > 0L) {
    appendicitis_ip_list[[idx]] <- df_app |>
      mutate(func_date = ymd(func_date), source = "ip")
  }
  idx <- idx + 1L
}

# F-3：bind 後加 distinct，避免同人同日多筆診斷在後續 inner_join 時 cartesian
appendicitis_records_raw <- bind_rows(
  bind_rows(appendicitis_op_list),
  bind_rows(appendicitis_ip_list)
) |>
  select(ID = id, FUNC_DATE = func_date, SOURCE = source) |>
  distinct(ID, FUNC_DATE, SOURCE)

# -- 闌尾切除術醫令（只從 OPDTO join OPDTE）---------------------------------

N_ECT_OP <- length(2016:2023) * 12L
appendectomy_op_list <- vector("list", N_ECT_OP)
idx <- 1L
for (y in 2016:2023) {
  log_progress(sprintf("  [Appendectomy/OP] year %d", y))
  for (m in 1:12) {
    file_o <- find_parquet("OPDTO", y, m)
    file_e <- find_parquet("OPDTE", y, m)
    
    if (is.na(file_o) || is.na(file_e)) {
      idx <- idx + 1L
      next
    }
    dt_o <- open_dataset(file_o) |>
      rename_with(tolower) |>
      select(all_of(join_keys), drug_no) |>
      filter(drug_no %in% appendectomy_codes) |>
      collect() |>
      rename(order_code = drug_no) |>
      normalise_join_keys()
    if (nrow(dt_o) == 0L) {
      idx <- idx + 1L
      next
    }
    dt_e <- open_dataset(file_e) |>
      rename_with(tolower) |>
      select(all_of(join_keys), id, func_date) |>
      filter(id %in% eligible_ids) |>
      collect() |>
      normalise_join_keys()
    appendectomy_op_list[[idx]] <- inner_join(dt_e, dt_o, by = join_keys)
    idx <- idx + 1L
  }
}

appendectomy_records_raw <- bind_rows(appendectomy_op_list) |>
  select(ID = id, FUNC_DATE = func_date) |>
  mutate(FUNC_DATE = ymd(FUNC_DATE))

# 組別 (A)：門診闌尾炎診斷 + ±7 天內門診闌尾切除術醫令
#   以原本的 inner_join 邏輯確認
appx_group_A <- appendicitis_records_raw |>
  filter(SOURCE == "op") |>
  select(ID, FUNC_DATE) |>
  inner_join(appendectomy_records_raw |> rename(OP_DATE = FUNC_DATE), by = "ID") |>
  filter(abs(as.numeric(OP_DATE - FUNC_DATE)) <= 7) |>
  select(ID, APPENDIX_DATE = FUNC_DATE) |>
  mutate(APPX_SOURCE = "op_op_medorder")

# 組別 (B)：住院闌尾炎（IPDTE）→ 假設為接受闌尾切除術
#   由於住院因闌尾炎幾乎必做闌尾切除，以 IPDTE 入院日為 APPENDIX_DATE
#   ⚠ 因本專案無 IPDTO，此組病人住院期間若有 CT 也無法觀測 → CT 暴露率必然低估
appx_group_B <- appendicitis_records_raw |>
  filter(SOURCE == "ip") |>
  select(ID, APPENDIX_DATE = FUNC_DATE) |>
  mutate(APPX_SOURCE = "ip_diagnosis_proxy")

# 聯集（A + B）：每人取最早一次事件 → 用於敏感性分析
# M-6：slice_min 前 arrange APPX_SOURCE，確保同日時 op_op_medorder 優先於 ip_diagnosis_proxy
appendicitis_appendectomy_full <- bind_rows(appx_group_A, appx_group_B) |>
  group_by(ID) |>
  arrange(APPENDIX_DATE, APPX_SOURCE, .by_group = TRUE) |>
  slice(1) |>
  ungroup()

# 主分析資料：只用 Group A（C-5）
appendicitis_appendectomy_main <- appx_group_A |>
  group_by(ID) |>
  arrange(APPENDIX_DATE, .by_group = TRUE) |>
  slice(1) |>
  ungroup()

# F-7 診斷：印出 §2.5 各步驟筆數，以便確認是哪一步斷掉
# （若主分析的 appendicitis-appendectomy-main.rds 只有 3KB，
#   通常意味著某一步抓到 0 筆。比對下面的數字可以快速定位）
log_progress("===== §2.5 闌尾炎建檔診斷 =====")
log_progress(sprintf("  闌尾炎診斷 (raw bind, distinct)        ：%s 筆",
                     format(nrow(appendicitis_records_raw), big.mark = ",")))
log_progress(sprintf("    └ 來自 OPDTE                          ：%s 筆",
                     format(nrow(appendicitis_records_raw |> filter(SOURCE == "op")),
                            big.mark = ",")))
log_progress(sprintf("    └ 來自 IPDTE                          ：%s 筆",
                     format(nrow(appendicitis_records_raw |> filter(SOURCE == "ip")),
                            big.mark = ",")))
log_progress(sprintf("  闌尾切除術醫令 (OPDTO 74002B/74004B)   ：%s 筆",
                     format(nrow(appendectomy_records_raw), big.mark = ",")))
log_progress(sprintf("  Group A (OP 診斷 + ±7 天內 OP 醫令)    ：%s 筆 / %s 人",
                     format(nrow(appx_group_A), big.mark = ","),
                     format(nrow(appx_group_A |> distinct(ID)), big.mark = ",")))
log_progress(sprintf("  Group B (IP 診斷 proxy)                 ：%s 筆 / %s 人",
                     format(nrow(appx_group_B), big.mark = ","),
                     format(nrow(appx_group_B |> distinct(ID)), big.mark = ",")))
log_progress(sprintf("  Main = appx_group_A 去重後              ：%s 人",
                     format(nrow(appendicitis_appendectomy_main), big.mark = ",")))
log_progress(sprintf("  Full = (A ∪ B) 去重後                   ：%s 人",
                     format(nrow(appendicitis_appendectomy_full), big.mark = ",")))
log_progress("===== §2.5 闌尾炎建檔診斷 結束 =====")

log_progress(sprintf("Study 3 納入人數明細："))
log_progress(sprintf("  (A) OP 診斷 + OP 醫令確認 (主分析)：%s",
                     format(nrow(appx_group_A |> distinct(ID)), big.mark = ",")))
log_progress(sprintf("  (B) IP 住院診斷 proxy   (敏感性)  ：%s",
                     format(nrow(appx_group_B |> distinct(ID)), big.mark = ",")))
log_progress(sprintf("  聯集（A + B，每人最早一次）       ：%s",
                     format(nrow(appendicitis_appendectomy_full), big.mark = ",")))

# 同時寫出兩個版本的 RDS：
#   appendicitis-appendectomy.rds       = 聯集（保留向後相容；含 APPX_SOURCE 欄位）
#   appendicitis-appendectomy-main.rds  = 僅 Group A（C-5 主分析用）
write_rds(appendicitis_appendectomy_full,
          file.path(INTERMEDIATE_PATH, "appendicitis-appendectomy.rds"))
write_rds(appendicitis_appendectomy_main,
          file.path(INTERMEDIATE_PATH, "appendicitis-appendectomy-main.rds"))
log_progress("已儲存 appendicitis-appendectomy.rds 與 appendicitis-appendectomy-main.rds")

# 後續 §6 主分析使用 main；敏感性分析使用 full
appendicitis_appendectomy      <- appendicitis_appendectomy_main
appendicitis_appendectomy_sens <- appendicitis_appendectomy_full

## 2.6 Sibling Relationship File (Study 2) ------------------------------------
# ⚠ 這段保持執行(不註解):
#   1. 建檔很快(只是讀一個 parquet 檔),沒必要 read_rds
#   2. 這個 .rds 不曾在前版產生過,read_rds 會直接掛掉
# 若日後想加快,可在跑完一次後手動把下方 build 區塊註解、改用 read_rds

log_progress("===== 建檔 2.6:手足配對關係 =====")

# F-4:原版在 arrow lazy query 上做 coalesce,當 F_ID/Guess_F 為
# dictionary-encoded 而另一邊是 character 時可能 type mismatch。
# 對齊 00-basic_rate.R:先 select + collect,所有 mutate 在 R 端做。
enrol_relation_ext <- open_dataset(
  file.path(PROCESSED_DATA_PATH, "enrol-relation-extend.parquet")
) |>
  select(ID, F_ID, M_ID, Guess_F, Guess_M) |>
  collect() |>
  mutate(
    EFF_F = coalesce(as.character(F_ID), as.character(Guess_F)),
    EFF_M = coalesce(as.character(M_ID), as.character(Guess_M))
  ) |>
  select(ID, EFF_F, EFF_M)

log_progress(sprintf("enrol_relation_ext 筆數:%s",
                     format(nrow(enrol_relation_ext), big.mark = ",")))
write_rds(enrol_relation_ext, file.path(INTERMEDIATE_PATH, "enrol-relation-ext.rds"))
log_progress("已儲存 enrol-relation-ext.rds")


## 2.7 Annual Enrollment: Urban Status & Income Quartile (C-1 + F-5 + F-6) ----
# 【建檔步驟:對齊 03i 風格已註解;F-6 修正版已跑過一次 rds 為正確欄位】
# ⚠ F-6 修正(依官方手冊 H_NHI_ENROL):
#   官方手冊 NHIRD ENROL 真實欄位是:
#     - ID1     = 被保險人身分證字號(主投保者)
#     - ID1_AMT = 投保金額
#     - ID_STATUS = 投保狀態
#     - ID_ROC = 身分證檢誤
#     - ID1_CITY = 投保單位地區代號
#   原版本誤用 INS_ID / INS_AMT 等不存在的欄位,故全面修正:
#     INS_ID  → ID1     (對齊手冊序號 7)
#     INS_AMT → ID1_AMT (對齊手冊序號 13)
#   判斷主投保者:id == id1(手冊「注意事項 5」)。
#
# ⚠ ENROL 月檔無 _10 sample group 後綴:呼叫 find_parquet 時必須傳
#   suffix = "" 才能找到檔案。
#
# ⚠ 二次跑分析時直接 read_rds 即可;要重建檔請取消下方註解區塊
#
# log_progress("===== 建檔 2.7:承保檔年度家戶縣市 + 薪資四分位(F-6 修正版)=====")
#
# eligible_ids_tbl <- tibble(ID = eligible_ids)
#
# # 兩個收集 list:
# #   child_ins_list:每個小孩 → 該月主投保者 id1
# #   ins_amt_list  :每個主投保者 → 該月薪資與縣市
# N_ENROL <- length(2016:2023) * 12L
# child_ins_list <- vector("list", N_ENROL)
# ins_amt_list   <- vector("list", N_ENROL)
# idx <- 1L
# # Step 2 主投保者篩選用的欄位集合(小寫;對齊手冊真實欄位名)
# ins_amt_select_cols <- c("id", "id1_amt", "id1_city", "prem_ym",
#                          "id_status", "id_roc", "id1")
# # Step 1 小孩篩選用的欄位集合
# child_select_cols   <- c("id", "id1", "prem_ym", "id_roc")
#
# for (y in 2016:2023) {
#   log_progress(sprintf("  [Enrol] year %d", y))
#   for (m in 1:12) {
#     # ⚠ ENROL 月檔沒有 _10 後綴,必須傳 suffix = ""
#     file_en <- find_parquet("ENROL", y, m, suffix = "")
#     if (is.na(file_en)) {
#       idx <- idx + 1L
#       next
#     }
#
#     # ===== Step 1:小孩 → 主投保者 id1 =====
#     child_raw <- open_dataset(file_en) |>
#       rename_with(tolower) |>
#       select(any_of(child_select_cols)) |>
#       filter(id %in% eligible_ids) |>
#       collect()
#     if (nrow(child_raw) > 0L) {
#       child_ins_list[[idx]] <- child_raw |>
#         mutate(
#           id_roc  = as.character(id_roc),
#           id1     = as.character(id1),
#           prem_ym = as.character(prem_ym),
#           year    = suppressWarnings(as.integer(substr(prem_ym, 1, 4)))
#         ) |>
#         filter(id_roc == "0", year == y) |>
#         distinct(id, id1, year, prem_ym)
#     }
#
#     # ===== Step 2:主投保者本人薪資與縣市 =====
#     ins_raw <- open_dataset(file_en) |>
#       rename_with(tolower) |>
#       select(any_of(ins_amt_select_cols)) |>
#       filter(id == id1) |>
#       collect()
#
#     if (nrow(ins_raw) > 0L) {
#       ins_clean <- ins_raw |>
#         mutate(
#           id        = as.character(id),
#           id_status = as.character(id_status),
#           id_roc    = as.character(id_roc),
#           id1_city  = as.character(id1_city),
#           prem_ym   = as.character(prem_ym),
#           year      = suppressWarnings(as.integer(substr(prem_ym, 1, 4)))
#         ) |>
#         filter(
#           !is.na(id1_amt),
#           id1_amt    >= AMT_CUTOFF,
#           id_status  %in% c("1", "2", "3"),
#           id_roc     == "0",
#           year       == y
#         ) |>
#         transmute(id1 = id, id1_amt, id1_city, prem_ym, year)
#       ins_amt_list[[idx]] <- ins_clean
#     }
#
#     idx <- idx + 1L
#   }
# }
#
# child_ins_monthly <- bind_rows(child_ins_list)
# ins_amt_monthly   <- bind_rows(ins_amt_list)
#
# log_progress(sprintf("child_ins_monthly:%s 筆(小孩-月-主投保者)",
#                      format(nrow(child_ins_monthly), big.mark = ",")))
# log_progress(sprintf("ins_amt_monthly  :%s 筆(主投保者-月-薪資)",
#                      format(nrow(ins_amt_monthly), big.mark = ",")))
#
# child_ins_monthly <- child_ins_monthly |>
#   mutate(
#     id1     = as.character(id1),
#     prem_ym = as.character(prem_ym),
#     year    = as.integer(year)
#   )
# ins_amt_monthly <- ins_amt_monthly |>
#   mutate(
#     id1     = as.character(id1),
#     prem_ym = as.character(prem_ym),
#     year    = as.integer(year)
#   )
#
# household_monthly <- child_ins_monthly |>
#   inner_join(ins_amt_monthly, by = c("id1", "year", "prem_ym"))
#
# log_progress(sprintf("household_monthly:%s 筆(小孩-月,含家戶薪資/縣市)",
#                      format(nrow(household_monthly), big.mark = ",")))
#
# annual_amt <- household_monthly |>
#   group_by(id, year) |>
#   summarise(
#     TOTAL_AMT = sum(id1_amt, na.rm = TRUE),
#     N_MONTHS  = dplyr::n(),
#     .groups   = "drop"
#   ) |>
#   rename(ID = id, YEAR = year)
#
# annual_city <- household_monthly |>
#   filter(!is.na(id1_city), id1_city != "") |>
#   group_by(id, year, id1_city) |>
#   summarise(N_MONTHS_CITY = dplyr::n(), .groups = "drop") |>
#   group_by(id, year) |>
#   arrange(desc(N_MONTHS_CITY), id1_city, .by_group = TRUE) |>
#   slice(1) |>
#   ungroup() |>
#   select(ID = id, YEAR = year, MAIN_CITY = id1_city)
#
# annual_amt_quartile <- annual_amt |>
#   group_by(YEAR) |>
#   mutate(
#     INCOME_Q = cut(
#       TOTAL_AMT,
#       breaks         = quantile(TOTAL_AMT, probs = c(0, 0.25, 0.5, 0.75, 1),
#                                 na.rm = TRUE, type = 7),
#       labels         = c("Q1", "Q2", "Q3", "Q4"),
#       include.lowest = TRUE
#     )
#   ) |>
#   ungroup()
#
# enrol_income_urban_raw <- annual_amt_quartile |>
#   left_join(annual_city, by = c("ID", "YEAR")) |>
#   mutate(URBAN_STATUS = classify_urban(MAIN_CITY)) |>
#   select(ID, YEAR, TOTAL_AMT, N_MONTHS, INCOME_Q, MAIN_CITY, URBAN_STATUS)
#
# log_progress(sprintf("enrol_income_urban_raw:%s 筆(人-年)",
#                      format(nrow(enrol_income_urban_raw), big.mark = ",")))
#
# n_eligible_with_income <- enrol_income_urban_raw |>
#   filter(!is.na(INCOME_Q)) |>
#   distinct(ID) |>
#   nrow()
# log_progress(sprintf("eligible_ids 中取得 INCOME_Q 的人數:%s / %s (%.1f%%)",
#                      format(n_eligible_with_income, big.mark = ","),
#                      format(length(eligible_ids), big.mark = ","),
#                      100 * n_eligible_with_income / length(eligible_ids)))
#
# write_rds(enrol_income_urban_raw,
#           file.path(INTERMEDIATE_PATH, "enrol-income-urban.rds"))
# log_progress("已儲存 enrol-income-urban.rds")

enrol_income_urban_raw <- read_rds(file.path(INTERMEDIATE_PATH, "enrol-income-urban.rds"))
enrol_income_urban     <- enrol_income_urban_raw
log_progress(sprintf("讀入 enrol-income-urban.rds:%s 筆(人-年)",
                     format(nrow(enrol_income_urban), big.mark = ",")))

# ----------------------------------------------------------------------------
# 取每人「研究基準年」的縣市 / 薪資四分位，用於主要模型的共變項
#   - 暴露者：使用 INDEX_DATE（= FIRST_CT_DATE）所在年份
#   - 非暴露者：使用 STUDY_START 所在年份（2016）
# 實作策略：先幫每人建「主要分配年」-> 取該年 enrol_income_urban 資料
#   （此步驟會在 3.x 建 cohort 時再做 left_join）
# ----------------------------------------------------------------------------

# 3. Build Analytic Cohorts ---------------------------------------------------

## 3.1 Base Eligibility -------------------------------------------------------
# 條件：2016-2023 年間年齡 0-18 歲，排除遺傳性癌症

log_progress("===== 3.1 Base cohort =====")

base_cohort <- pers_info |>
  anti_join(hereditary_exclusion, by = "ID") |>
  mutate(
    BIRTH_DATE = resolve_birth_date(ID_BIRTHYM),
    AGE_2016   = calc_age_years(BIRTH_DATE, STUDY_START),
    AGE_2023   = calc_age_years(BIRTH_DATE, STUDY_END)
  ) |>
  filter(
    !is.na(BIRTH_DATE),
    BIRTH_DATE <= STUDY_END,
    BIRTH_DATE + years(MAX_AGE_AT_INDEX + 1L) > STUDY_START
  )

log_progress(sprintf("base_cohort：%s", format(nrow(base_cohort), big.mark = ",")))

## 3.2 CT Exposure Summary per Person ----------------------------------------

ct_summary <- ct_records |>
  filter(FUNC_DATE >= STUDY_START, FUNC_DATE <= STUDY_END) |>
  group_by(ID) |>
  summarise(
    FIRST_CT_DATE = min(FUNC_DATE),
    N_CT_TOTAL    = n(),
    CT_TYPES      = paste(sort(unique(CT_TYPE)), collapse = "; "),
    .groups = "drop"
  ) |>
  left_join(
    ct_records |>
      group_by(ID) |>
      arrange(FUNC_DATE, ORDER_CODE, .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      select(ID, FIRST_CT_TYPE = CT_TYPE),
    by = "ID"
  )

## 3.3 Prior Malignancy Exclusion Flag ----------------------------------------
# 定義：確診日期早於研究起始日（STUDY_START = 2016-01-01）

prior_malignancy_ids <- malignancy_dx |>
  semi_join(base_cohort, by = "ID") |>
  filter(FIRST_MALIGNANCY_DATE < STUDY_START) |>
  distinct(ID)

## 3.4 Cohort Flowchart / Attrition Table ------------------------------------

n_pers_info       <- nrow(pers_info)
n_age_eligible    <- length(eligible_ids)
n_hereditary_excl <- nrow(hereditary_exclusion)
n_base_cohort     <- nrow(base_cohort)
n_prior_mal_excl  <- nrow(prior_malignancy_ids)
n_after_prior_mal <- nrow(base_cohort |> anti_join(prior_malignancy_ids, by = "ID"))

flowchart <- tibble(
  STEP = c(
    "1. NHIRD 全人口",
    "2. 年齡符合（2016-2023 間曾為 0-18 歲）",
    "3. 排除：遺傳性腫瘤疾患",
    "4. 排除：2016 年前既有惡性腫瘤確診",
    "5. 最終 base cohort（三個 Study 共用）"
  ),
  N_REMAINING = c(
    n_pers_info,
    n_age_eligible,
    n_base_cohort,
    n_after_prior_mal,
    n_after_prior_mal
  ),
  N_EXCLUDED_THIS_STEP = c(
    NA_integer_,
    n_pers_info - n_age_eligible,
    n_age_eligible - n_base_cohort,
    n_prior_mal_excl,
    NA_integer_
  ),
  NOTE = c(
    "pers_info 全筆",
    "出生年 1997–2023；NHIRD 民國 89–112 年均有紀錄",
    "ICD-10：Q850/1/8/9、Z1501/2/9；查詢範圍 2000–2023",
    "FIRST_MALIGNANCY_DATE < 2016-01-01；查詢範圍 2000–2023",
    "進入 Study 1 / 2 / 3 的共同起點"
  )
)

flowchart |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "cohort-flowchart.csv"))

cat("\n========== Cohort Attrition ==========\n")
print(flowchart, n = Inf)

## 3.5 Helper: Attach Income Quartile & Urban Status to a Cohort --------------
# 給定一個 cohort（含 ID 與 INDEX_DATE），取 INDEX_DATE 所在年份的
# INCOME_Q 與 URBAN_STATUS 作為該人的共變項。
#
# 邏輯：
#   (1) 依 INDEX_DATE 算 INDEX_YEAR
#   (2) left_join(enrol_income_urban, by = c(ID, INDEX_YEAR = YEAR))
#   (3) 若該人該年無家戶承保紀錄（極少數情況：父母皆無工作或在保身分外）
#       → INCOME_Q / URBAN_STATUS = NA，分析時作 missing indicator
#
# C-1 修正後 enrol_income_urban 已使用家戶主投保者 proxy，覆蓋率應大幅提升。

attach_income_urban <- function(cohort, enrol_df = enrol_income_urban) {
  cohort |>
    mutate(INDEX_YEAR = as.integer(year(INDEX_DATE))) |>
    left_join(
      enrol_df |> select(ID, YEAR, INCOME_Q, URBAN_STATUS),
      by = c("ID", "INDEX_YEAR" = "YEAR")
    )
}

# 4. Study 1: Population-Level Cohort -----------------------------------------
# 暴露：研究期間有任何 CT 掃描（0-18 歲）
# 對照：研究期間無 CT 掃描
# 結果：首次 CT 後 >2 年發生惡性腫瘤
# 控制：年齡、性別、薪資四分位、六都/非六都

## 4.1 Build Study 1 Cohort --------------------------------------------------

log_progress("===== 4.1 Study 1 cohort =====")

# F-8：base_cohort 已從 pers_info 衍生，已包含 ID_S 與 BIRTH_DATE。
# 原版本在這裡又做 left_join(pers_info |> ... select(ID, ID_S))，會導致
# 兩邊都有 ID_S → dplyr 自動重新命名為 ID_S.x / ID_S.y → 後續 mutate
# 用 ID_S 找不到欄位而報錯。修正方式：刪除這個重複的 join。
study1_cohort <- base_cohort |>
  anti_join(prior_malignancy_ids, by = "ID") |>
  left_join(ct_summary, by = "ID") |>
  mutate(
    EXPOSED       = !is.na(FIRST_CT_DATE),
    INDEX_DATE    = if_else(EXPOSED, FIRST_CT_DATE, STUDY_START),
    AGE_AT_INDEX  = calc_age_years(BIRTH_DATE, INDEX_DATE),
    # I-5：分開保留兩個對齊年齡，描述統計時兩組可比
    AGE_AT_2016   = calc_age_years(BIRTH_DATE, STUDY_START),
    FIRST_CT_AGE  = if_else(EXPOSED, calc_age_years(BIRTH_DATE, FIRST_CT_DATE), NA_integer_),
    SEX           = if_else(ID_S == "1", "Male", "Female")
  ) |>
  filter(AGE_AT_INDEX <= MAX_AGE_AT_INDEX) |>
  left_join(malignancy_dx, by = "ID") |>
  mutate(
    ELIGIBLE_MALIGNANCY_DATE = if_else(
      !is.na(FIRST_MALIGNANCY_DATE) &
        FIRST_MALIGNANCY_DATE > INDEX_DATE + years(LATENCY_YEARS),
      FIRST_MALIGNANCY_DATE,
      NA_Date_
    ),
    OUTCOME     = !is.na(ELIGIBLE_MALIGNANCY_DATE),
    CENSOR_DATE = pmin(coalesce(ELIGIBLE_MALIGNANCY_DATE, STUDY_END), STUDY_END),
    TIME_MONTHS = as.numeric(
      interval(INDEX_DATE + years(LATENCY_YEARS), CENSOR_DATE) / months(1)
    )
  ) |>
  filter(TIME_MONTHS >= 0) |>
  attach_income_urban()

## 4.2 Study 1 Analysis -------------------------------------------------------

### 4.2.1 Descriptive Statistics ----------------------------------------------
# I-5 修正：原 MEDIAN_AGE 用 AGE_AT_INDEX，但暴露組 INDEX = FIRST_CT_DATE、
# 對照組 INDEX = STUDY_START，兩組年齡不可比。改成：
#   - MEDIAN_AGE_AT_2016：兩組共同基準（2016-01-01）的年齡，可比
#   - MEDIAN_FIRST_CT_AGE：僅暴露組有效（NA 表示對照組）

study1_descriptive <- study1_cohort |>
  group_by(EXPOSED) |>
  summarise(
    N                    = n(),
    N_MALE               = sum(SEX == "Male"),
    N_FEMALE             = sum(SEX == "Female"),
    MEDIAN_AGE_AT_2016   = median(AGE_AT_2016, na.rm = TRUE),
    MEDIAN_FIRST_CT_AGE  = median(FIRST_CT_AGE, na.rm = TRUE),
    N_MALIGNANCY         = sum(OUTCOME),
    INCIDENCE_RATE_K_PM  = sum(OUTCOME) / sum(TIME_MONTHS) * 1000,
    N_METRO              = sum(URBAN_STATUS == "Metro", na.rm = TRUE),
    N_NONMETRO           = sum(URBAN_STATUS == "Non-Metro", na.rm = TRUE),
    N_NA_URBAN           = sum(is.na(URBAN_STATUS)),
    .groups              = "drop"
  )

study1_descriptive |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-descriptive.csv"))

### 4.2.2 Risk Ratio (Crude) -------------------------------------------------

study1_rr_overall <- study1_cohort |>
  group_by(EXPOSED) |>
  summarise(
    N_EVENT = sum(OUTCOME),
    N_TOTAL = n(),
    RISK    = N_EVENT / N_TOTAL,
    .groups = "drop"
  ) |>
  summarise(
    RR    = RISK[EXPOSED] / RISK[!EXPOSED],
    RR_CI = paste0(
      "(95% CI: ",
      round(exp(log(RR) - 1.96 * sqrt(
        1/N_EVENT[EXPOSED] + 1/N_EVENT[!EXPOSED] -
          1/N_TOTAL[EXPOSED] - 1/N_TOTAL[!EXPOSED])), 2),
      " - ",
      round(exp(log(RR) + 1.96 * sqrt(
        1/N_EVENT[EXPOSED] + 1/N_EVENT[!EXPOSED] -
          1/N_TOTAL[EXPOSED] - 1/N_TOTAL[!EXPOSED])), 2),
      ")"
    )
  )

# M-5：改用直接常數，避免 by = character() 的 cross join
N_TOTAL_TRUE  <- sum(study1_cohort$EXPOSED)
N_TOTAL_FALSE <- sum(!study1_cohort$EXPOSED)

study1_rr_by_type <- study1_cohort |>
  filter(OUTCOME) |>
  count(MALIGNANCY_TYPE, EXPOSED) |>
  pivot_wider(
    names_from   = EXPOSED,
    values_from  = n,
    names_prefix = "N_EVENT_",
    values_fill  = 0L
  ) |>
  mutate(
    N_TOTAL_TRUE  = N_TOTAL_TRUE,
    N_TOTAL_FALSE = N_TOTAL_FALSE,
    RR = (N_EVENT_TRUE  / N_TOTAL_TRUE) /
      (N_EVENT_FALSE / N_TOTAL_FALSE)
  )

study1_rr_by_type |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-rr-by-malignancy-type.csv"))

### 4.2.3 Linear Probability Model（加回 INCOME_Q、URBAN_STATUS）------------

study1_lm <- lm(
  OUTCOME ~ EXPOSED + SEX + AGE_AT_INDEX + INCOME_Q + URBAN_STATUS,
  data = study1_cohort
)

tidy(study1_lm, conf.int = TRUE) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-linear-regression.csv"))

### 4.2.4 Logistic Regression（加回 INCOME_Q、URBAN_STATUS）-----------------

study1_logit <- glm(
  OUTCOME ~ EXPOSED + SEX + AGE_AT_INDEX + INCOME_Q + URBAN_STATUS,
  data   = study1_cohort,
  family = binomial(link = "logit")
)

tidy(study1_logit, conf.int = TRUE, exponentiate = TRUE) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-logistic-regression.csv"))

### 4.2.5 Cox Proportional Hazards (副分析, C-2) ----------------------------
# 主分析仍以研究計畫指定的 LM / Logistic 為準（§4.2.3 / §4.2.4）。
# 此處新增 Cox 比例風險模型作為副分析，將 person-time 納入考量，
# 處理暴露/對照組 follow-up 長度不對稱的問題（暴露組從 FIRST_CT_DATE
# 起算、對照組從 STUDY_START 起算）。

study1_cox <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS,
  data = study1_cohort
)

tidy(study1_cox, conf.int = TRUE, exponentiate = TRUE) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-cox-regression.csv"))

# 模型診斷：比例風險假設檢定（Schoenfeld residuals）
study1_cox_zph <- cox.zph(study1_cox)
study1_cox_zph_table <- tibble(
  TERM    = rownames(study1_cox_zph$table),
  CHISQ   = study1_cox_zph$table[, "chisq"],
  DF      = study1_cox_zph$table[, "df"],
  P_VALUE = study1_cox_zph$table[, "p"]
)
study1_cox_zph_table |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-cox-zph.csv"))

### 4.2.6 Stratified Analyses (加回 income / urban 分層) ---------------------

study1_strat_sex    <- stratify_rr(study1_cohort, SEX)
study1_strat_agegrp <- stratify_rr(
  study1_cohort |>
    mutate(AGE_GROUP = cut(AGE_AT_INDEX,
                           breaks = c(0, 5, 10, 15, 18),
                           include.lowest = TRUE)),
  AGE_GROUP
)
study1_strat_cttype <- stratify_rr(
  study1_cohort |> filter(EXPOSED),
  FIRST_CT_TYPE
)
study1_strat_income <- stratify_rr(study1_cohort, INCOME_Q)
study1_strat_urban  <- stratify_rr(study1_cohort, URBAN_STATUS)

list(
  sex      = study1_strat_sex,
  age_grp  = study1_strat_agegrp,
  ct_type  = study1_strat_cttype,
  income_q = study1_strat_income,
  urban    = study1_strat_urban
) |>
  imap(~ write_csv(.x, file.path(
    OUTPUT_TABLE_PATH, paste0("study1-strat-", .y, ".csv")
  )))

# 5. Study 2: Sibling-Matched Cohort ------------------------------------------
# 暴露：研究期間有 CT（0-18 歲）
# 對照：同胞手足在相同期間無 CT
# 結果：首次 CT 後 >2 年發生惡性腫瘤
# 分層：性別一致性、首次 CT 年齡組、薪資四分位、六都/非六都

## 5.1 Identify Sibling Pairs -------------------------------------------------

log_progress("===== 5.1 Sibling pairs =====")

sib_by_father <- enrol_relation_ext |>
  filter(!is.na(EFF_F)) |>
  select(ID, EFF_F) |>
  inner_join(
    enrol_relation_ext |> filter(!is.na(EFF_F)) |> select(ID_2 = ID, EFF_F),
    by = "EFF_F"
  ) |>
  filter(ID < ID_2) |>
  select(ID_1 = ID, ID_2)

sib_by_mother <- enrol_relation_ext |>
  filter(!is.na(EFF_M)) |>
  select(ID, EFF_M) |>
  inner_join(
    enrol_relation_ext |> filter(!is.na(EFF_M)) |> select(ID_2 = ID, EFF_M),
    by = "EFF_M"
  ) |>
  filter(ID < ID_2) |>
  select(ID_1 = ID, ID_2)

sibling_pairs <- bind_rows(sib_by_father, sib_by_mother) |>
  distinct(ID_1, ID_2) |>
  semi_join(base_cohort, by = c("ID_1" = "ID")) |>
  semi_join(base_cohort, by = c("ID_2" = "ID"))

# I-1：對齊研究計畫 "sex concordance/discordance between siblings"
# 在 pair 層級標記 Both_Male / Both_Female / Mixed
sibling_pairs <- sibling_pairs |>
  left_join(
    pers_info |> select(ID, SEX_1 = ID_S),
    by = c("ID_1" = "ID")
  ) |>
  left_join(
    pers_info |> select(ID, SEX_2 = ID_S),
    by = c("ID_2" = "ID")
  ) |>
  mutate(
    PAIR_SEX = case_when(
      SEX_1 == "1" & SEX_2 == "1" ~ "Both_Male",
      SEX_1 == "2" & SEX_2 == "2" ~ "Both_Female",
      !is.na(SEX_1) & !is.na(SEX_2) ~ "Mixed",
      TRUE                          ~ NA_character_
    )
  ) |>
  select(ID_1, ID_2, PAIR_SEX)

write_rds(sibling_pairs, file.path(INTERMEDIATE_PATH, "sibling-pairs.rds"))
log_progress(sprintf("已儲存 sibling-pairs.rds：%s 對",
                     format(nrow(sibling_pairs), big.mark = ",")))

## 5.2 Build Study 2 Cohort --------------------------------------------------

# I-7：select 加入 BIRTH_DATE 與 AGE_AT_2016 / FIRST_CT_AGE
study2_exposed <- study1_cohort |>
  filter(EXPOSED) |>
  select(ID, BIRTH_DATE, INDEX_DATE, FIRST_CT_DATE,
         AGE_AT_INDEX, AGE_AT_2016, FIRST_CT_AGE, SEX,
         OUTCOME, ELIGIBLE_MALIGNANCY_DATE, TIME_MONTHS,
         FIRST_CT_TYPE, N_CT_TOTAL,
         INCOME_Q, URBAN_STATUS, MALIGNANCY_TYPE)

# sibling_pairs 以 ID_1 < ID_2 去重，不保證 ID_1 為暴露者，故雙向展開
# 雙向展開時帶上 PAIR_SEX 與配對者 ID（FAMILY_ID = pair）
study2_control <- bind_rows(
  sibling_pairs |>
    transmute(ID_EXPOSED = ID_1, ID_CONTROL = ID_2, PAIR_SEX,
              PAIR_ID = paste(pmin(ID_1, ID_2), pmax(ID_1, ID_2), sep = "_")),
  sibling_pairs |>
    transmute(ID_EXPOSED = ID_2, ID_CONTROL = ID_1, PAIR_SEX,
              PAIR_ID = paste(pmin(ID_1, ID_2), pmax(ID_1, ID_2), sep = "_"))
) |>
  inner_join(
    study2_exposed |> select(ID_EXPOSED = ID, INDEX_DATE),
    by = "ID_EXPOSED"
  ) |>
  select(-ID_EXPOSED) |>
  anti_join(ct_summary, by = c("ID_CONTROL" = "ID")) |>
  anti_join(prior_malignancy_ids, by = c("ID_CONTROL" = "ID")) |>
  left_join(
    pers_info |>
      mutate(BIRTH_DATE = resolve_birth_date(ID_BIRTHYM)) |>
      select(ID, BIRTH_DATE, ID_S),
    by = c("ID_CONTROL" = "ID")
  ) |>
  mutate(
    AGE_AT_INDEX = calc_age_years(BIRTH_DATE, INDEX_DATE),
    AGE_AT_2016  = calc_age_years(BIRTH_DATE, STUDY_START),
    FIRST_CT_AGE = NA_integer_,   # control 無 CT
    SEX          = if_else(ID_S == "1", "Male", "Female"),
    EXPOSED      = FALSE
  ) |>
  # C-3：必須在 INDEX_DATE 那天介於 0–18 歲（已出生且未滿 19）
  filter(AGE_AT_INDEX >= 0, AGE_AT_INDEX <= MAX_AGE_AT_INDEX) |>
  left_join(malignancy_dx |> rename(ID_CONTROL = ID), by = "ID_CONTROL") |>
  mutate(
    ELIGIBLE_MALIGNANCY_DATE = if_else(
      !is.na(FIRST_MALIGNANCY_DATE) &
        FIRST_MALIGNANCY_DATE > INDEX_DATE + years(LATENCY_YEARS),
      FIRST_MALIGNANCY_DATE,
      NA_Date_
    ),
    OUTCOME     = !is.na(ELIGIBLE_MALIGNANCY_DATE),
    CENSOR_DATE = pmin(coalesce(ELIGIBLE_MALIGNANCY_DATE, STUDY_END), STUDY_END),
    TIME_MONTHS = as.numeric(
      interval(INDEX_DATE + years(LATENCY_YEARS), CENSOR_DATE) / months(1)
    )
  ) |>
  filter(TIME_MONTHS >= 0) |>
  rename(ID = ID_CONTROL) |>
  attach_income_urban()

# 暴露組對齊：bind 前先把 PAIR_SEX / PAIR_ID 帶到 exposed 那邊
# （以同人作為 exposed，可能對應多個 pair；保留所有對應，每對一列）
study2_exposed_paired <- bind_rows(
  sibling_pairs |>
    filter(ID_1 %in% study2_exposed$ID) |>
    transmute(ID = ID_1, PAIR_SEX,
              PAIR_ID = paste(pmin(ID_1, ID_2), pmax(ID_1, ID_2), sep = "_")),
  sibling_pairs |>
    filter(ID_2 %in% study2_exposed$ID) |>
    transmute(ID = ID_2, PAIR_SEX,
              PAIR_ID = paste(pmin(ID_1, ID_2), pmax(ID_1, ID_2), sep = "_"))
) |>
  inner_join(study2_exposed, by = "ID") |>
  mutate(EXPOSED = TRUE)

study2_cohort <- bind_rows(
  study2_exposed_paired,
  study2_control
) |>
  # I-1：SEX_CONCORDANCE 改為 pair-level（直接沿用 PAIR_SEX）
  mutate(SEX_CONCORDANCE = PAIR_SEX)

## 5.3 Study 2 Analysis -------------------------------------------------------

### 5.3.1 Risk Ratio ----------------------------------------------------------

study2_rr_overall <- study2_cohort |>
  group_by(EXPOSED) |>
  summarise(
    N_EVENT = sum(OUTCOME),
    N_TOTAL = n(),
    RISK    = N_EVENT / N_TOTAL,
    .groups = "drop"
  ) |>
  summarise(RR = RISK[EXPOSED] / RISK[!EXPOSED])

study2_rr_by_type <- study2_cohort |>
  filter(OUTCOME) |>
  group_by(MALIGNANCY_TYPE, EXPOSED) |>
  summarise(N = n(), .groups = "drop") |>
  pivot_wider(names_from = EXPOSED, values_from = N,
              names_prefix = "N_", values_fill = 0L) |>
  mutate(RR = N_TRUE / N_FALSE)

study2_rr_by_type |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-rr-by-malignancy-type.csv"))

### 5.3.2 Cox Proportional Hazards (副分析, C-2) ----------------------------
# 主分析仍為 RR（研究計畫指定）；此處加 Cox 將 person-time 納入考量。
# 注意：本版本未處理 within-pair correlation（C-4 暫不處理），故 SE 可能略低估。

study2_cox <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS,
  data = study2_cohort
)

tidy(study2_cox, conf.int = TRUE, exponentiate = TRUE) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-cox-regression.csv"))

study2_cox_zph <- cox.zph(study2_cox)
study2_cox_zph_table <- tibble(
  TERM    = rownames(study2_cox_zph$table),
  CHISQ   = study2_cox_zph$table[, "chisq"],
  DF      = study2_cox_zph$table[, "df"],
  P_VALUE = study2_cox_zph$table[, "p"]
)
study2_cox_zph_table |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-cox-zph.csv"))

### 5.3.3 Stratified Analyses (加回 income / urban 分層) --------------------

study2_strat_agegrp <- stratify_rr(
  study2_cohort |>
    mutate(AGE_GROUP = cut(AGE_AT_INDEX,
                           breaks = c(0, 5, 10, 15, 18),
                           include.lowest = TRUE)),
  AGE_GROUP
)
study2_strat_sex_concordance <- stratify_rr(study2_cohort, SEX_CONCORDANCE)
study2_strat_income <- stratify_rr(study2_cohort, INCOME_Q)
study2_strat_urban  <- stratify_rr(study2_cohort, URBAN_STATUS)

list(
  age_grp         = study2_strat_agegrp,
  sex_concordance = study2_strat_sex_concordance,
  income_q        = study2_strat_income,
  urban           = study2_strat_urban
) |>
  imap(~ write_csv(.x, file.path(
    OUTPUT_TABLE_PATH, paste0("study2-strat-", .y, ".csv")
  )))

# 6. Study 3: Appendicitis Cohort ---------------------------------------------
# 族群：0-18 歲，2016-2023 年因闌尾炎接受闌尾切除術
#       【主分析】OPDTO 醫令確認（Group A） — appendicitis_appendectomy
#       【敏感性】聯集（A + B；含 IPDTE 住院 proxy） — appendicitis_appendectomy_sens
# 暴露：闌尾切除術【前後 7 天】內有非頭部型 CT（腹部 CT proxy）
# 對照：無 CT
# 結果：手術後 >2 年發生惡性腫瘤
# 控制：年齡、性別、薪資四分位、六都/非六都
#
# C-5 修正：
#   (1) 暴露窗從 [-3, +1] 天擴大為 [-7, +7] 天，與 §2.5 Group A
#       闌尾切除術 ±7 天的標準對齊；臨床上術中／術後 CT 也可被觀測
#   (2) 主分析限縮 Group A（OPDTO 醫令確認）。Group B（IPDTE proxy）
#       因本專案無 IPDTO，住院期間 CT 醫令無法觀測，幾乎一定被誤標
#       EXPOSED = FALSE，故移到敏感性分析

EXPOSURE_WINDOW_DAYS <- 7L

## 6.1 Build Study 3 Cohort --------------------------------------------------

log_progress("===== 6.1 Study 3 cohort =====")

# 共用 helper：給定一個 appendicitis 表，回傳 study3 cohort
build_study3_cohort <- function(appendicitis_df, label) {
  base <- appendicitis_df |>
    semi_join(base_cohort, by = "ID") |>
    anti_join(prior_malignancy_ids, by = "ID") |>
    filter(APPENDIX_DATE >= STUDY_START, APPENDIX_DATE <= STUDY_END) |>
    left_join(
      pers_info |>
        mutate(BIRTH_DATE = resolve_birth_date(ID_BIRTHYM)) |>
        select(ID, BIRTH_DATE, ID_S),
      by = "ID"
    ) |>
    mutate(
      AGE_AT_APPENDIX = calc_age_years(BIRTH_DATE, APPENDIX_DATE),
      # I-5：對 Study 3 而言，AGE_AT_APPENDIX 對所有人都是 INDEX 那天的年齡
      #       本身就可比；因此不另設 AGE_AT_2016。但保留 FIRST_CT_AGE 供描述。
      SEX             = if_else(ID_S == "1", "Male", "Female")
    ) |>
    filter(AGE_AT_APPENDIX <= MAX_AGE_AT_INDEX, AGE_AT_APPENDIX >= 0)
  
  # C-5：暴露判斷 — 術前後 ±7 天內有非頭部型 CT（腹部 CT proxy）
  ct_flag <- ct_records |>
    filter(ORDER_CODE %in% abdomen_ct_proxy_codes) |>
    inner_join(base |> select(ID, APPENDIX_DATE), by = "ID") |>
    filter(
      FUNC_DATE >= APPENDIX_DATE - days(EXPOSURE_WINDOW_DAYS),
      FUNC_DATE <= APPENDIX_DATE + days(EXPOSURE_WINDOW_DAYS)
    ) |>
    distinct(ID) |>
    mutate(EXPOSED = TRUE)
  
  # 同時記錄首次腹部 CT 的日期（描述用）
  first_ab_ct <- ct_records |>
    filter(ORDER_CODE %in% abdomen_ct_proxy_codes) |>
    inner_join(base |> select(ID, APPENDIX_DATE), by = "ID") |>
    filter(
      FUNC_DATE >= APPENDIX_DATE - days(EXPOSURE_WINDOW_DAYS),
      FUNC_DATE <= APPENDIX_DATE + days(EXPOSURE_WINDOW_DAYS)
    ) |>
    group_by(ID) |>
    arrange(FUNC_DATE, ORDER_CODE, .by_group = TRUE) |>
    slice(1) |>
    ungroup() |>
    select(ID, FIRST_AB_CT_DATE = FUNC_DATE)
  
  out <- base |>
    left_join(ct_flag, by = "ID") |>
    replace_na(list(EXPOSED = FALSE)) |>
    left_join(first_ab_ct, by = "ID") |>
    mutate(
      FIRST_CT_AGE = if_else(EXPOSED,
                             calc_age_years(BIRTH_DATE, FIRST_AB_CT_DATE),
                             NA_integer_)
    ) |>
    left_join(malignancy_dx, by = "ID") |>
    mutate(
      INDEX_DATE               = APPENDIX_DATE,
      ELIGIBLE_MALIGNANCY_DATE = if_else(
        !is.na(FIRST_MALIGNANCY_DATE) &
          FIRST_MALIGNANCY_DATE > INDEX_DATE + years(LATENCY_YEARS),
        FIRST_MALIGNANCY_DATE,
        NA_Date_
      ),
      OUTCOME     = !is.na(ELIGIBLE_MALIGNANCY_DATE),
      CENSOR_DATE = pmin(coalesce(ELIGIBLE_MALIGNANCY_DATE, STUDY_END), STUDY_END),
      TIME_MONTHS = as.numeric(
        interval(INDEX_DATE + years(LATENCY_YEARS), CENSOR_DATE) / months(1)
      )
    ) |>
    filter(TIME_MONTHS >= 0) |>
    attach_income_urban() |>
    mutate(COHORT_LABEL = label)
  
  out
}

# 主分析（Group A only）
study3_cohort <- build_study3_cohort(appendicitis_appendectomy, "main_groupA")

# 敏感性分析（Group A + Group B 聯集）
study3_cohort_sens <- build_study3_cohort(appendicitis_appendectomy_sens,
                                          "sensitivity_AB")

log_progress(sprintf("Study 3 主分析 cohort（Group A）：%s 人，EXPOSED = %s",
                     format(nrow(study3_cohort), big.mark = ","),
                     format(sum(study3_cohort$EXPOSED), big.mark = ",")))
log_progress(sprintf("Study 3 敏感性 cohort（A+B 聯集）：%s 人，EXPOSED = %s",
                     format(nrow(study3_cohort_sens), big.mark = ","),
                     format(sum(study3_cohort_sens$EXPOSED), big.mark = ",")))

## 6.2 Study 3 Analysis -------------------------------------------------------

### 6.2.1 Descriptive Statistics ----------------------------------------------
# I-5：Study 3 中 AGE_AT_APPENDIX 已是兩組共同基準（INDEX = APPENDIX_DATE
# 對暴露/對照組都成立），可比；額外呈現 FIRST_CT_AGE（暴露組才有意義）。

study3_descriptive <- study3_cohort |>
  group_by(EXPOSED) |>
  summarise(
    N                   = n(),
    N_MALE              = sum(SEX == "Male"),
    N_FEMALE            = sum(SEX == "Female"),
    MEDIAN_AGE_APPENDIX = median(AGE_AT_APPENDIX, na.rm = TRUE),
    MEDIAN_FIRST_CT_AGE = median(FIRST_CT_AGE, na.rm = TRUE),
    N_MALIGNANCY        = sum(OUTCOME),
    INCIDENCE_RATE_K_PM = sum(OUTCOME) / sum(TIME_MONTHS) * 1000,
    N_METRO             = sum(URBAN_STATUS == "Metro", na.rm = TRUE),
    N_NONMETRO          = sum(URBAN_STATUS == "Non-Metro", na.rm = TRUE),
    N_NA_URBAN          = sum(is.na(URBAN_STATUS)),
    .groups             = "drop"
  )

study3_descriptive |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-descriptive.csv"))

### 6.2.2 Risk Ratio ----------------------------------------------------------

study3_rr_overall <- study3_cohort |>
  group_by(EXPOSED) |>
  summarise(
    N_EVENT = sum(OUTCOME),
    N_TOTAL = n(),
    RISK    = N_EVENT / N_TOTAL,
    .groups = "drop"
  ) |>
  summarise(RR = RISK[EXPOSED] / RISK[!EXPOSED])

study3_rr_by_type <- study3_cohort |>
  filter(OUTCOME) |>
  group_by(MALIGNANCY_TYPE, EXPOSED) |>
  summarise(N = n(), .groups = "drop") |>
  pivot_wider(names_from = EXPOSED, values_from = N,
              names_prefix = "N_", values_fill = 0L) |>
  mutate(RR = N_TRUE / N_FALSE)

study3_rr_by_type |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-rr-by-malignancy-type.csv"))

### 6.2.3 Linear Regression（加回 INCOME_Q、URBAN_STATUS）-------------------

study3_lm <- lm(
  OUTCOME ~ EXPOSED + SEX + AGE_AT_APPENDIX + INCOME_Q + URBAN_STATUS,
  data = study3_cohort
)

tidy(study3_lm, conf.int = TRUE) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-linear-regression.csv"))

### 6.2.4 Logistic Regression（加回 INCOME_Q、URBAN_STATUS）------------------

study3_logit <- glm(
  OUTCOME ~ EXPOSED + SEX + AGE_AT_APPENDIX + INCOME_Q + URBAN_STATUS,
  data   = study3_cohort,
  family = binomial(link = "logit")
)

tidy(study3_logit, conf.int = TRUE, exponentiate = TRUE) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-logistic-regression.csv"))

### 6.2.5 Cox Proportional Hazards (副分析, C-2) ----------------------------

study3_cox <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_APPENDIX +
    INCOME_Q + URBAN_STATUS,
  data = study3_cohort
)

tidy(study3_cox, conf.int = TRUE, exponentiate = TRUE) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-cox-regression.csv"))

study3_cox_zph <- cox.zph(study3_cox)
study3_cox_zph_table <- tibble(
  TERM    = rownames(study3_cox_zph$table),
  CHISQ   = study3_cox_zph$table[, "chisq"],
  DF      = study3_cox_zph$table[, "df"],
  P_VALUE = study3_cox_zph$table[, "p"]
)
study3_cox_zph_table |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-cox-zph.csv"))

### 6.2.6 Stratified Analyses (加回 income / urban 分層) --------------------

study3_strat_sex    <- stratify_rr(study3_cohort, SEX)
study3_strat_agegrp <- stratify_rr(
  study3_cohort |>
    mutate(AGE_GROUP = cut(AGE_AT_APPENDIX,
                           breaks = c(0, 5, 10, 15, 18),
                           include.lowest = TRUE)),
  AGE_GROUP
)
study3_strat_income <- stratify_rr(study3_cohort, INCOME_Q)
study3_strat_urban  <- stratify_rr(study3_cohort, URBAN_STATUS)

list(
  sex      = study3_strat_sex,
  age_grp  = study3_strat_agegrp,
  income_q = study3_strat_income,
  urban    = study3_strat_urban
) |>
  imap(~ write_csv(.x, file.path(
    OUTPUT_TABLE_PATH, paste0("study3-strat-", .y, ".csv")
  )))

### 6.2.7 Sensitivity: 含 IPDTE 住院 proxy 的整體 RR (C-5) ------------------

study3_rr_sens_overall <- study3_cohort_sens |>
  group_by(EXPOSED) |>
  summarise(
    N_EVENT = sum(OUTCOME),
    N_TOTAL = n(),
    RISK    = N_EVENT / N_TOTAL,
    .groups = "drop"
  ) |>
  summarise(RR_SENS_AB = RISK[EXPOSED] / RISK[!EXPOSED])

study3_rr_sens_overall |>
  mutate(NOTE = "Sensitivity: A + B union; CT exposure may be underestimated for Group B due to missing IPDTO") |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-rr-sensitivity-AB.csv"))

# 不同 APPX_SOURCE 的描述（暴露率比對）
study3_appx_source_breakdown <- study3_cohort_sens |>
  left_join(
    appendicitis_appendectomy_sens |> select(ID, APPX_SOURCE),
    by = "ID"
  ) |>
  group_by(APPX_SOURCE, EXPOSED) |>
  summarise(N = n(), .groups = "drop") |>
  pivot_wider(names_from = EXPOSED, values_from = N,
              names_prefix = "EXPOSED_", values_fill = 0L) |>
  mutate(
    EXPOSURE_PCT = EXPOSED_TRUE / (EXPOSED_TRUE + EXPOSED_FALSE) * 100
  )

study3_appx_source_breakdown |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study3-appx-source-exposure-pct.csv"))

# 7. Summary Tables and Figures -----------------------------------------------

## 7.1 Combined RR Summary Table -----------------------------------------------

rr_summary <- bind_rows(
  study1_rr_overall |> mutate(Study = "Study 1 (Population Cohort)"),
  study2_rr_overall |> mutate(Study = "Study 2 (Sibling-Matched)"),
  study3_rr_overall |> mutate(Study = "Study 3 (Appendicitis Cohort)")
)

rr_summary |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "rr-summary-all-studies.csv"))

## 7.2 CT Type Distribution (Descriptive) -------------------------------------

ct_bodypart_summary <- ct_records |>
  filter(FUNC_DATE >= STUDY_START, FUNC_DATE <= STUDY_END) |>
  semi_join(base_cohort, by = "ID") |>
  count(CT_TYPE) |>
  arrange(desc(n)) |>
  mutate(PCT = n / sum(n) * 100)

ct_bodypart_summary |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "ct-type-distribution.csv"))

## 7.3 Figure: CT Exposure by Body Part and Age Group -------------------------

ct_age_bodypart_fig <- ct_records |>
  filter(FUNC_DATE >= STUDY_START, FUNC_DATE <= STUDY_END) |>
  semi_join(base_cohort, by = "ID") |>
  left_join(
    pers_info |>
      mutate(BIRTH_DATE = resolve_birth_date(ID_BIRTHYM)) |>
      select(ID, BIRTH_DATE),
    by = "ID"
  ) |>
  mutate(
    AGE_AT_CT = calc_age_years(BIRTH_DATE, FUNC_DATE),
    AGE_GROUP = cut(AGE_AT_CT,
                    breaks = c(0, 2, 5, 10, 15, 18),
                    include.lowest = TRUE,
                    labels = c("0-2", "3-5", "6-10", "11-15", "16-18"))
  ) |>
  count(CT_TYPE, AGE_GROUP) |>
  ggplot(aes(x = AGE_GROUP, y = n, fill = CT_TYPE)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "CT Scan Frequency by CT Type and Age Group (0-18 years, 2016-2023)",
    subtitle = "Taiwan NHIRD｜Head-type CT = 33067-33069B；Body-type CT = 33070-33072B",
    x        = "Age Group (years)",
    y        = "Number of CT Scans",
    fill     = "CT Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(OUTPUT_FIGURE_PATH, "ct-bodypart-age-distribution.pdf"),
  plot     = ct_age_bodypart_fig,
  width    = 10,
  height   = 6
)

## 7.4 Figure: Kaplan-Meier (Study 1) ----------------------------------------

km_study1 <- survfit(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED,
  data = study1_cohort
)

km_study1_plot <- km_study1 |>
  tidy() |>
  mutate(
    CIF   = 1 - estimate,
    Group = if_else(strata == "EXPOSED=TRUE", "CT Exposed", "CT Unexposed")
  ) |>
  ggplot(aes(x = time, y = CIF, color = Group)) +
  geom_step() +
  geom_ribbon(
    aes(ymin = 1 - conf.high, ymax = 1 - conf.low, fill = Group),
    alpha = 0.15, color = NA
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01)) +
  labs(
    title    = "Cumulative Incidence of Malignancy by CT Exposure Status",
    subtitle = "Study 1: Population Cohort, latency >2 years (Taiwan NHIRD 2016-2023)",
    x        = "Months from Index Date + 2 Years",
    y        = "Cumulative Incidence of Malignancy",
    color    = "Group",
    fill     = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(OUTPUT_FIGURE_PATH, "km-malignancy-study1.pdf"),
  plot     = km_study1_plot,
  width    = 8,
  height   = 5
)

# 8. Cohort Counts Summary ----------------------------------------------------
# 匯出一份統整的人數表：從建檔階段、base cohort 篩選、
# 到三個 study cohort 以及各分層類別，每一層的 N_PERSONS（distinct ID）
# 與 N_RECORDS（紀錄筆數，若適用）。

log_progress("===== 8. 輸出 cohort-counts.csv =====")

safe_n_distinct <- function(x) {
  if (is.null(x) || length(x) == 0) 0L else dplyr::n_distinct(x, na.rm = TRUE)
}

count_row <- function(category, subcategory, n_persons, n_records = NA_integer_, note = NA_character_) {
  tibble(
    CATEGORY     = category,
    SUBCATEGORY  = subcategory,
    N_PERSONS    = as.integer(n_persons),
    N_RECORDS    = as.integer(n_records),
    NOTE         = note
  )
}

## 8.1 建檔階段人數 -----------------------------------------------------------

counts_build <- bind_rows(
  count_row("1. Build / Intermediate",
            "pers-info (全母體)",
            safe_n_distinct(pers_info$ID),
            nrow(pers_info),
            "2.1 pers_info"),
  count_row("1. Build / Intermediate",
            "eligible_ids (年齡 0-18 於 2016-2023)",
            length(eligible_ids),
            length(eligible_ids),
            "2.1 研究年齡範圍預篩"),
  count_row("1. Build / Intermediate",
            "ct-records-all.rds (門診 CT 事件)",
            safe_n_distinct(ct_records_raw$ID),
            nrow(ct_records_raw),
            "2.2 OPDTO × OPDTE"),
  count_row("1. Build / Intermediate",
            "malignancy-dx.rds (首次惡性腫瘤診斷)",
            safe_n_distinct(malignancy_raw$ID),
            nrow(malignancy_raw),
            "2.3 OPDTE + IPDTE; 2000-2023"),
  count_row("1. Build / Intermediate",
            "hereditary-records-full.rds (遺傳疾患紀錄)",
            safe_n_distinct(hereditary_records_raw$ID),
            nrow(hereditary_records_raw),
            "2.4 OPDTE + IPDTE; 2000-2023"),
  count_row("1. Build / Intermediate",
            "hereditary-exclusion-ids.rds (遺傳疾患排除清單)",
            safe_n_distinct(hereditary_exclusion_raw$ID),
            nrow(hereditary_exclusion_raw),
            "2.4 每人首次"),
  count_row("1. Build / Intermediate",
            "appendicitis-appendectomy.rds: 組A (OP 醫令確認, 主分析)",
            safe_n_distinct(appx_group_A$ID),
            nrow(appx_group_A),
            "2.5 OPDTE 診斷 + OPDTO 74002B/74004B ±7 天"),
  count_row("1. Build / Intermediate",
            "appendicitis-appendectomy.rds: 組B (IP 住院 proxy, 敏感性)",
            safe_n_distinct(appx_group_B$ID),
            nrow(appx_group_B),
            "2.5 IPDTE 闌尾炎診斷 (K35/K36/K37 三碼前綴)"),
  count_row("1. Build / Intermediate",
            "appendicitis-appendectomy-main.rds (僅組 A，主分析)",
            safe_n_distinct(appendicitis_appendectomy_main$ID),
            nrow(appendicitis_appendectomy_main),
            "C-5：主分析使用"),
  count_row("1. Build / Intermediate",
            "appendicitis-appendectomy.rds: 聯集 A+B (敏感性)",
            safe_n_distinct(appendicitis_appendectomy_full$ID),
            nrow(appendicitis_appendectomy_full),
            "2.5 每人最早一次"),
  count_row("1. Build / Intermediate",
            "sibling-pairs.rds (手足配對)",
            safe_n_distinct(c(sibling_pairs$ID_1, sibling_pairs$ID_2)),
            nrow(sibling_pairs),
            "5.1 N_RECORDS = 手足對數"),
  count_row("1. Build / Intermediate",
            "enrol-income-urban.rds (人-年薪資/縣市)",
            safe_n_distinct(enrol_income_urban_raw$ID),
            nrow(enrol_income_urban_raw),
            "2.7 N_RECORDS = 人-年筆數")
)

## 8.2 Base Cohort 篩選流程（以 distinct ID 統計）---------------------------

n_age_eligible    <- length(eligible_ids)
n_base_cohort     <- nrow(base_cohort)
n_prior_mal_excl  <- nrow(prior_malignancy_ids)
n_after_prior_mal <- nrow(base_cohort |> anti_join(prior_malignancy_ids, by = "ID"))

counts_attrition <- bind_rows(
  count_row("2. Base Cohort Attrition",
            "Step 1: NHIRD 全人口",
            safe_n_distinct(pers_info$ID),
            NA_integer_,
            "pers_info"),
  count_row("2. Base Cohort Attrition",
            "Step 2: 年齡符合 (2016-2023 間曾 0-18 歲)",
            n_age_eligible,
            NA_integer_,
            NA_character_),
  count_row("2. Base Cohort Attrition",
            "Step 3: 排除遺傳性腫瘤疾患 (剩餘人數)",
            n_base_cohort,
            NA_integer_,
            sprintf("排除 %s 人", format(n_age_eligible - n_base_cohort, big.mark = ","))),
  count_row("2. Base Cohort Attrition",
            "Step 4: 排除 2016 年前既有惡性腫瘤 (最終 base)",
            n_after_prior_mal,
            NA_integer_,
            sprintf("排除 %s 人", format(n_prior_mal_excl, big.mark = ",")))
)

## 8.3 三個 Study 的 Exposure 分佈 -------------------------------------------

counts_study_exposure <- bind_rows(
  # Study 1
  study1_cohort |>
    count(EXPOSED) |>
    mutate(
      CATEGORY    = "3. Study 1 (Population)",
      SUBCATEGORY = paste0("EXPOSED = ", EXPOSED),
      N_PERSONS   = as.integer(n),
      N_RECORDS   = as.integer(n),
      NOTE        = "CT 暴露 vs 未暴露"
    ) |>
    select(CATEGORY, SUBCATEGORY, N_PERSONS, N_RECORDS, NOTE),
  # Study 2
  study2_cohort |>
    count(EXPOSED) |>
    mutate(
      CATEGORY    = "3. Study 2 (Sibling-Matched)",
      SUBCATEGORY = paste0("EXPOSED = ", EXPOSED),
      N_PERSONS   = as.integer(n),
      N_RECORDS   = as.integer(n),
      NOTE        = "暴露者 vs 未暴露手足"
    ) |>
    select(CATEGORY, SUBCATEGORY, N_PERSONS, N_RECORDS, NOTE),
  # Study 3
  study3_cohort |>
    count(EXPOSED) |>
    mutate(
      CATEGORY    = "3. Study 3 (Appendicitis)",
      SUBCATEGORY = paste0("EXPOSED = ", EXPOSED),
      N_PERSONS   = as.integer(n),
      N_RECORDS   = as.integer(n),
      NOTE        = "腹部 CT vs 無 CT (闌尾炎族群內)"
    ) |>
    select(CATEGORY, SUBCATEGORY, N_PERSONS, N_RECORDS, NOTE)
)

## 8.4 三個 Study 的分層人數 --------------------------------------------------

# Helper：依單一分層變項算 EXPOSED × Stratum 表
stratum_counts <- function(cohort_df, cohort_name, strat_col) {
  col_sym <- sym(strat_col)
  cohort_df |>
    mutate(STRAT = as.character(!!col_sym)) |>
    count(EXPOSED, STRAT) |>
    mutate(
      CATEGORY    = cohort_name,
      SUBCATEGORY = sprintf("%s = %s | EXPOSED = %s", strat_col, STRAT, EXPOSED),
      N_PERSONS   = as.integer(n),
      N_RECORDS   = as.integer(n),
      NOTE        = NA_character_
    ) |>
    select(CATEGORY, SUBCATEGORY, N_PERSONS, N_RECORDS, NOTE)
}

# 事先備好年齡分組欄位
study1_cohort_counts <- study1_cohort |>
  mutate(AGE_GROUP = as.character(cut(AGE_AT_INDEX,
                                      breaks = c(0, 5, 10, 15, 18),
                                      include.lowest = TRUE)))
study2_cohort_counts <- study2_cohort |>
  mutate(AGE_GROUP = as.character(cut(AGE_AT_INDEX,
                                      breaks = c(0, 5, 10, 15, 18),
                                      include.lowest = TRUE)))
study3_cohort_counts <- study3_cohort |>
  mutate(AGE_GROUP = as.character(cut(AGE_AT_APPENDIX,
                                      breaks = c(0, 5, 10, 15, 18),
                                      include.lowest = TRUE)))

counts_strata <- bind_rows(
  # Study 1
  stratum_counts(study1_cohort_counts, "4. Study 1 × Sex",          "SEX"),
  stratum_counts(study1_cohort_counts, "4. Study 1 × Age Group",    "AGE_GROUP"),
  stratum_counts(study1_cohort_counts, "4. Study 1 × Income Q",     "INCOME_Q"),
  stratum_counts(study1_cohort_counts, "4. Study 1 × Urban Status", "URBAN_STATUS"),
  # Study 1 的 CT 類型（僅對 EXPOSED 有意義）
  study1_cohort |>
    filter(EXPOSED) |>
    count(FIRST_CT_TYPE) |>
    transmute(
      CATEGORY    = "4. Study 1 × First CT Type",
      SUBCATEGORY = sprintf("FIRST_CT_TYPE = %s (EXPOSED only)", FIRST_CT_TYPE),
      N_PERSONS   = as.integer(n),
      N_RECORDS   = as.integer(n),
      NOTE        = NA_character_
    ),
  # Study 2
  stratum_counts(study2_cohort_counts, "4. Study 2 × Age Group",       "AGE_GROUP"),
  stratum_counts(study2_cohort_counts, "4. Study 2 × Sex Concordance", "SEX_CONCORDANCE"),
  stratum_counts(study2_cohort_counts, "4. Study 2 × Income Q",        "INCOME_Q"),
  stratum_counts(study2_cohort_counts, "4. Study 2 × Urban Status",    "URBAN_STATUS"),
  # Study 3
  stratum_counts(study3_cohort_counts, "4. Study 3 × Sex",          "SEX"),
  stratum_counts(study3_cohort_counts, "4. Study 3 × Age Group",    "AGE_GROUP"),
  stratum_counts(study3_cohort_counts, "4. Study 3 × Income Q",     "INCOME_Q"),
  stratum_counts(study3_cohort_counts, "4. Study 3 × Urban Status", "URBAN_STATUS"),
  # Study 3 × 闌尾炎來源（敏感性 cohort 才有 A vs B 對比；主分析全為 A）
  study3_cohort_sens |>
    left_join(
      appendicitis_appendectomy_sens |> select(ID, APPX_SOURCE),
      by = "ID"
    ) |>
    count(EXPOSED, APPX_SOURCE) |>
    transmute(
      CATEGORY    = "4. Study 3 (Sensitivity) × Appendicitis Source",
      SUBCATEGORY = sprintf("APPX_SOURCE = %s | EXPOSED = %s", APPX_SOURCE, EXPOSED),
      N_PERSONS   = as.integer(n),
      N_RECORDS   = as.integer(n),
      NOTE        = "op_op_medorder = OP 診斷+OPDTO 醫令 (主分析); ip_diagnosis_proxy = IPDTE 住院 proxy (敏感性)"
    )
)

## 8.5 結果事件（發生惡性腫瘤）人數 -------------------------------------------

counts_outcome <- bind_rows(
  study1_cohort |>
    group_by(EXPOSED) |>
    summarise(n_event = sum(OUTCOME), n_total = n(), .groups = "drop") |>
    transmute(
      CATEGORY    = "5. Study 1 Outcome (Malignancy)",
      SUBCATEGORY = sprintf("EXPOSED = %s", EXPOSED),
      N_PERSONS   = as.integer(n_event),
      N_RECORDS   = as.integer(n_total),
      NOTE        = "N_PERSONS = 發生惡性腫瘤; N_RECORDS = cohort 總人數"
    ),
  study2_cohort |>
    group_by(EXPOSED) |>
    summarise(n_event = sum(OUTCOME), n_total = n(), .groups = "drop") |>
    transmute(
      CATEGORY    = "5. Study 2 Outcome (Malignancy)",
      SUBCATEGORY = sprintf("EXPOSED = %s", EXPOSED),
      N_PERSONS   = as.integer(n_event),
      N_RECORDS   = as.integer(n_total),
      NOTE        = "N_PERSONS = 發生惡性腫瘤; N_RECORDS = cohort 總人數"
    ),
  study3_cohort |>
    group_by(EXPOSED) |>
    summarise(n_event = sum(OUTCOME), n_total = n(), .groups = "drop") |>
    transmute(
      CATEGORY    = "5. Study 3 Outcome (Malignancy)",
      SUBCATEGORY = sprintf("EXPOSED = %s", EXPOSED),
      N_PERSONS   = as.integer(n_event),
      N_RECORDS   = as.integer(n_total),
      NOTE        = "N_PERSONS = 發生惡性腫瘤; N_RECORDS = cohort 總人數"
    )
)

## 8.6 合併輸出 ---------------------------------------------------------------

cohort_counts_all <- bind_rows(
  counts_build,
  counts_attrition,
  counts_study_exposure,
  counts_strata,
  counts_outcome
)

cohort_counts_all |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "cohort-counts.csv"))

log_progress(sprintf("已儲存 cohort-counts.csv：%s 列", nrow(cohort_counts_all)))

cat("\n========== Cohort Counts Summary ==========\n")
print(cohort_counts_all, n = Inf)

log_progress("===== 全部完成 =====")