# 0. Description --------------------------------------------------------------
# Author: [Author Name]
# Last Modified: 2026-05-10 (v03g)
# Goal: Examine the association between CT scan exposure in childhood and
# subsequent malignancy diagnosis using Taiwan's National Health
# Insurance Research Database (NHIRD), corresponding to three study
# designs:
#         Study 1 - Population-level cohort (all children aged 0-18)
#         Study 2 - Sibling/twin matched cohort
#         Study 3 - Appendicitis-restricted cohort
#
# ═══════════════════════════════════════════════════════════════════════
# 本版本（v03g）變更：依 2026-05-08 meeting comments 全面增補。
# ═══════════════════════════════════════════════════════════════════════
# 【G-1】§2.2 ct-records-all.rds 新增 CT_SETTING 欄位（ER vs Clinic）
#        - OPDTE 的 case_type 為健保就醫類別代碼：
#            "02" = 急診 → ER
#            "01"/"03"/其他 = 一般門診 + 連續處方箋等 → Clinic
#        - 在 collect 後依 case_type 衍生 CT_SETTING ∈ {"ER", "Clinic"}
#        - 同時保留 case_type 原值供日後驗證
#        - ⚠ 此項變更需重建 ct-records-all.rds
#
# 【G-2】§3.2 ct_summary 加入 FIRST_CT_SETTING（同日多筆時 ER 優先）
#
# 【G-3】三個 study 的 cohort（§4.1 / §5.2 / §6.1）統一新增衍生欄位：
#        - AGE_AT_DIAGNOSIS：發生事件者的診斷時年齡（NA 若無事件）
#        - FIRST_CT_AGE_GROUP：暴露組 first CT 時的年齡分組（NA 若未暴露）
#        - AGE_AT_DIAGNOSIS_GROUP：發生事件者的診斷時年齡分組
#
# 【G-4】§4.2 / §5.3 / §6.2 描述統計全面重構：
#        - 新增 make_descriptive() helper，三個 study 共用同一張 schema
#        - Study 2 補上原本缺失的描述統計表（study2-descriptive.csv）
#        - 三個 study 的描述統計欄位統一：sex / age_index / age_first_ct /
#          age_at_diagnosis / income / urban / ct_setting / outcome
#
# 【G-5】三個 study 新增分層分析：
#        - by FIRST_CT_SETTING（ER vs Clinic），僅暴露組（study1/3）
#          或全 cohort（study2 因 EXPOSED 那邊沒有 unexposed 對照）
#        - by AGE_AT_DIAGNOSIS_GROUP（呼應 PDF 研究計畫所列）
#        - by FIRST_CT_AGE_GROUP（呼應 PDF "stratified by age at CT exposure"）
#
# 【G-6】Sensitivity analyses 全面建立：
#        - run_with_latency(yrs)：可變 latency wrapper（主分析 2y、敏感性 5y）
#        - CT setting sensitivity：限 Clinic CT 為暴露的 sub-analysis
#        - Multiple-CT dose-response：N_CT_GROUP ∈ {1, 2-3, ≥4}
#        - Study 2 within-pair clustering：Cox + cluster(PAIR_ID)
#        - 全部結果寫入 OUTPUT_TABLE_PATH/sensitivity/ 子目錄
#
# 【G-7】§5.3.1 study2_rr_by_type 公式從 case-ratio 改為正確 risk-ratio
#        （原版 N_TRUE / N_FALSE 是案例比，正確為 (N_TRUE/N_TOT_TRUE) /
#        (N_FALSE/N_TOT_FALSE)）
#
# ═══════════════════════════════════════════════════════════════════════
# 第二批變更（2026-05-10 後）：
# ═══════════════════════════════════════════════════════════════════════
# 【H-1】Study 1 新增 Time-varying exposure 平行分析（保留原版）
#        - 對齊 Mathews 2013 BMJ / Smoll 2023 AJNR / Pearce 2012 Lancet 標準
#        - 模式 (i)：lag period 內算未暴露 person-years
#          每人在 [entry, first_CT + lag) 期間 exposed=0
#                 [first_CT + lag, exit)  期間 exposed=1
#          無 CT 者整段都是未暴露
#        - 用 survival::tmerge 把 study1_cohort 切成 long format
#        - 平行於原版 fixed-INDEX_DATE Cox，輸出 study1-cox-tv.csv
#        - 副版：baseline N_CT_GROUP（1 / 2-3 / >=4）做 dose-response，
#          沿用原 cohort、不時變
#
# 【H-2】三個 study 新增 type-specific latency sensitivity
#        - run_with_latency(yrs, type) wrapper：給定 latency 年數重跑 cohort
#        - 主分析仍用統一 LATENCY_YEARS = 2L（呼應 Mathews）
#        - 敏感性 1：solid tumor 用 5y latency、hematologic 用 2y latency
#          （對齊 Pearce 2012：「leukemia 至少 2y、solid 至少 5y 才會發生」）
#        - 敏感性 2：所有 type 統一用 5y latency
#        - 結果寫入 OUTPUT_SENS_PATH/latency-{config}/
#
# 【H-3】Study 2 Cox model 加 strata(PAIR_ID) 平行版本
#        - 原版 §5.3.2（無 strata）保留：屬於 "extended cohort" 設計，
#          只是把 sibling 拉進當 control，沒做 within-family confounding 控制
#        - 新增 §5.3.4 加 strata(PAIR_ID)：真正的 sibling-matched design
#          對齊 Lichtenstein NEJM 2000 / D'Onofrio sibling comparison 文獻
#        - 兩版都輸出，命名為 study2-cox-regression.csv（無 strata）
#          與 study2-cox-strata.csv（有 strata）
#        - 注意：strata(PAIR_ID) 後僅 within-pair 有 outcome 變異的 pair
#          貢獻資訊，N 會大幅縮水，HR 可能不穩定（這是 sibling design
#          固有代價，非 bug）
#
# 【H-4】Study window 從 2016-2023 擴展為 2000-2023（第二階段）
#        - STUDY_START 仍維持 2016-01-01 作為「主分析入組起點」
#        - 新增 EXTENDED_STUDY_START = 2000-01-01 作為敏感性分析窗
#        - 出生年範圍從 1997-2023 擴展為 1982-2023（讓 1982 出生者
#          在 2000-01-01 為 18 歲，仍符合 0-18 歲入組）
#        - 主分析 cohort 命名 base_cohort（沿用 2016-2023 設計）
#        - 敏感性 cohort 命名 base_cohort_extended（2000-2023）
#        - 寫入 OUTPUT_SENS_PATH/extended-window/
#
# 【H-5】ICD-9 惡性腫瘤分 11 類對應（為 H-4 服務）
#        - 03f 對 CODE_VERSION == "ICD-9" 統一歸 "icd9_unclassified"
#        - 本版新增 classify_malignancy_icd9() helper：
#            ICD-9 140-149 → head_and_neck
#            150-159 → intestinal
#            160-165 → chest_and_lung
#            170    → bone
#            171    → cnt
#            174-184 → breast_and_female
#            185-189 → urinary_and_fertile
#            191-192 → brain_and_cns
#            193-194 → endocrine
#            195-199 → other_malignant
#            200-208 → lymphoma_and_leukemia
#        - 與 ICD-10 11 類同名，可直接對照分析
#        - ⚠ 此 mapping 為初版，臨床端若需更精細請提出
#
# 【H-6】§3.3 Prior malignancy exclusion 改為個人化 entry date
#        - 原版以 STUDY_START（2016-01-01）為界排除「2016 前」的惡性腫瘤
#        - H-4 擴大 window 後，每個人應以「自己的 study entry date」為界
#          entry_date = max(birth_date, EXTENDED_STUDY_START)
#        - 1985 出生者的 entry = 2000-01-01（他們在 1990-1999 期間發生
#          的惡性腫瘤無法觀察，這是已知 limitation）
#        - 主分析（2016-2023 window）邏輯不變，只在 H-4 敏感性 cohort 套用
#
# 【H-7】§2.7 ENROL/income 重建範圍延伸到 2000-2023（H-4 的相依檔）
#        - 原版只跑 2016-2023 共 96 個月，無法支援 H-4 extended cohort
#          中 1982-1996 出生者在 2000-2015 期間的 income / urban 賦值
#        - 本版啟用 build 範圍 2000-2023 共 288 個月（耗時長，但必要）
#        - 注意：AMT_CUTOFF = 2000 點未依年度通膨調整；
#          INCOME_Q 為「該年內相對排名」，跨年比較需謹慎
#        - 加 pre-flight check：列出 2000-2015 各月檔存在性 + 第一個讀
#          到的檔印 names() 確認欄位結構一致
#
# 【H-8】Prior malignancy 排除改用 effective start date
#        - 對暴露者：effective_start_date = max(entry_date, first_CT_date)
#          → 排除「first CT 之前」發生的 malignancy
#          理由：CT 暴露之前已有的 cancer 不該歸因為「兒童期 CT 暴露之後」
#          的事件
#        - 對未暴露者：effective_start_date = entry_date
#        - 主分析 cohort 因 entry = 2016-01-01、CT 也在 2016 之後，
#          通常 first_CT > entry，所以這個改動會比舊邏輯多排除一些人
#        - H-4 extended cohort：1982-1996 出生者 entry = 2000-01-01，
#          first_CT 可能在他們 25 歲時，需排除 25 歲前發生的 cancer
#
# 【H-9】撤回「19 歲時 censor」的提案 — 文獻不這樣做
#        Mathews 2013 / Pearce 2012 / EPI-CT 2023 / Korea 2025 的標準做法：
#        「兒童期 CT 暴露」是 inclusion 條件，但 follow-up 不在 18 歲
#        就 censor，而是延續到 study_end / 死亡 / 第一次 cancer 之最早。
#        - Mathews：1985 掃 CT 的 5 歲小孩追蹤到 2007 年（27 歲）
#        - Pearce：「主分析不限制成年 follow-up；敏感性分析才限到
#                   28 歲（brain）/ 25 歲（leukemia），但結果差異很小」
#        - EPI-CT：51% 的事件發生在 20+ 歲
#        理由：輻射致癌潛伏期跨成年（實體瘤 5-10y、有時 20y+），
#        若 18 歲就 censor，輻射效應根本還沒有時間發生
#
#        H-4 extended cohort 的 1982-1996 出生者，他們的 follow-up 大半
#        在成年期，這符合文獻標準做法，不是 bug。
#
# 【H-10】AGE_AT_DIAGNOSIS_GROUP 在 H-4 extended cohort 擴大分組
#        - 主分析（2016-2023）維持 c(0, 5, 10, 15, 20, Inf)
#          → labels = 0-5 / 6-10 / 11-15 / 16-20 / 21+
#          這個 cohort 最大年齡 27 歲（1997 出生 → 2023 年），
#          21+ 一格夠用
#        - H-4 extended cohort：擴成 c(0, 5, 10, 15, 20, 30, 40, Inf)
#          → labels = 0-5 / 6-10 / 11-15 / 16-20 / 21-30 / 31-40 / 41+
#          1982 出生者在 2023 年最大 41 歲，需要更細分組

# ═══════════════════════════════════════════════════════════════════════
# RDS 重建清單（依 H-1 ~ H-10 變更整理）
# ═══════════════════════════════════════════════════════════════════════
# 需重建（取消對應 §2.x 的註解，跑一次後重新註解回去）：
#   ✅ ct-records-all.rds      — G-1（CT_SETTING）+ H-4（範圍 2000-2023）
#   ✅ enrol-income-urban.rds  — H-7（範圍 2000-2023）
#   ✅ malignancy-dx.rds       — 建議重建以套用 H-5（ICD-9 分類），
#                                  但程式碼裡有 in-memory 補修邏輯，
#                                  若不重建也可運作
#
# 不需重建（邏輯沒變或本來就實際執行）：
#   ⚪ hereditary-records-full.rds
#   ⚪ hereditary-exclusion-ids.rds
#   ⚪ appendicitis-appendectomy.rds（§2.5 本來就實際執行）
#   ⚪ appendicitis-appendectomy-main.rds（§2.5 同上）
#   ⚪ enrol-relation-ext.rds（§2.6 本來就實際執行）
#   ⚪ sibling-pairs.rds（§5.1 本來就實際執行）
# ═══════════════════════════════════════════════════════════════════════
#
# v03f (2026-04-27) 變更：對照 00-basic_rate.R 修正 arrow lazy query 問題。
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
OUTPUT_SENS_PATH    <- file.path(OUTPUT_TABLE_PATH, "sensitivity")

# 確保資料夾存在
dir.create(INTERMEDIATE_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TABLE_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIGURE_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_SENS_PATH, recursive = TRUE, showWarnings = FALSE)

# Study window
STUDY_START      <- ymd(20160101)
STUDY_END        <- ymd(20231231)
MAX_AGE_AT_INDEX <- 18L
LATENCY_YEARS    <- 2L

# H-4：第二階段擴展窗口（敏感性分析）
EXTENDED_STUDY_START <- ymd(20000101)

# H-2：type-specific latency（敏感性分析用）
# 對齊 Pearce 2012 Lancet：白血病 ≥2y、實體瘤 ≥5y 才會發生
LATENCY_HEMATOLOGIC <- 2L
LATENCY_SOLID       <- 5L

# H-1：time-varying analysis 的 lag period（同 LATENCY_YEARS 主版）
TV_LAG_YEARS <- 2L

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

## 1.4.0a Malignancy Type Classification (ICD-10 + ICD-9, H-5) ----------------
# ICD-10 11 類分類（已在 §2.3 build malignancy 時用 case_when 寫死，這裡集中
# 為 helper，方便 H-4 敏感性分析重複使用）
classify_malignancy_icd10 <- function(icd3) {
  # icd3 為 3 字元 character，如 "C16"、"C50"
  icd_num <- suppressWarnings(as.integer(substr(icd3, 2, 3)))
  dplyr::case_when(
    is.na(icd_num)         ~ NA_character_,
    icd_num <= 14          ~ "head_and_neck",
    icd_num <= 26          ~ "intestinal",
    icd_num <= 39          ~ "chest_and_lung",
    icd_num <= 41          ~ "bone",
    icd_num <= 49          ~ "cnt",
    icd_num <= 58          ~ "breast_and_female",
    icd_num <= 68          ~ "urinary_and_fertile",
    icd_num <= 72          ~ "brain_and_cns",
    icd_num <= 75          ~ "endocrine",
    icd_num <= 80          ~ "other_malignant",
    icd_num <= 96          ~ "lymphoma_and_leukemia",
    TRUE                   ~ "other_malignant"
  )
}

# H-5：ICD-9 11 類分類（與 ICD-10 同名，可直接對照分析）
# 本 mapping 為初版，臨床端若需更精細的對應請提出
#   ICD-9 140-149 = 唇/口腔/咽 → head_and_neck
#   ICD-9 150-159 = 消化道（食道/胃/腸/肝/胰 等）→ intestinal
#   ICD-9 160-165 = 呼吸/胸內 → chest_and_lung
#   ICD-9 170    = 骨/關節軟骨 → bone
#   ICD-9 171    = 結締/其他軟組織 → cnt
#   ICD-9 174-175 = 乳房（女男）→ breast_and_female
#   ICD-9 179-184 = 子宮/卵巢等女性生殖 → breast_and_female（合併）
#   ICD-9 185-189 = 攝護腺/睪丸/腎/膀胱/泌尿 → urinary_and_fertile
#   ICD-9 191-192 = 腦/CNS → brain_and_cns
#   ICD-9 193-194 = 甲狀腺/其他內分泌 → endocrine
#   ICD-9 195-199 = 其他惡性 → other_malignant
#   ICD-9 200-208 = 淋巴/造血 → lymphoma_and_leukemia
classify_malignancy_icd9 <- function(icd3) {
  # icd3 為 3 字元 character，如 "150"、"174"
  icd_num <- suppressWarnings(as.integer(icd3))
  dplyr::case_when(
    is.na(icd_num)                                      ~ NA_character_,
    icd_num >= 140 & icd_num <= 149                     ~ "head_and_neck",
    icd_num >= 150 & icd_num <= 159                     ~ "intestinal",
    icd_num >= 160 & icd_num <= 165                     ~ "chest_and_lung",
    icd_num == 170                                      ~ "bone",
    icd_num == 171                                      ~ "cnt",
    icd_num >= 174 & icd_num <= 184                     ~ "breast_and_female",
    icd_num >= 185 & icd_num <= 189                     ~ "urinary_and_fertile",
    icd_num >= 191 & icd_num <= 192                     ~ "brain_and_cns",
    icd_num >= 193 & icd_num <= 194                     ~ "endocrine",
    icd_num >= 195 & icd_num <= 199                     ~ "other_malignant",
    icd_num >= 200 & icd_num <= 208                     ~ "lymphoma_and_leukemia",
    TRUE                                                 ~ "other_malignant"
  )
}

# H-2：判定某個 11 類是否屬於 hematologic（type-specific latency 用）
# 主要是 lymphoma_and_leukemia；其餘皆視為 solid
HEMATOLOGIC_TYPES <- c("lymphoma_and_leukemia")
is_hematologic <- function(malignancy_type) {
  malignancy_type %in% HEMATOLOGIC_TYPES
}

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

# G-5 強化版 stratify：同時輸出 N_EVENT × N_TOTAL × RR(95% CI)
# 用於對外輸出時需要可驗證的 N（HWDC 規範要 cell N >= 3 才能釋出）
stratify_rr_full <- function(data, strat_var, label = "stratum") {
  cell <- data |>
    mutate(
      .strat   = as.character({{ strat_var }}),
      .strat   = if_else(is.na(.strat), "Missing", .strat),
      .exposed = if_else(EXPOSED, "TRUE", "FALSE")
    ) |>
    group_by(.strat, .exposed) |>
    summarise(
      N_EVENT = sum(OUTCOME),
      N_TOTAL = dplyr::n(),
      .groups = "drop"
    )
  wide <- cell |>
    pivot_wider(
      names_from = .exposed,
      values_from = c(N_EVENT, N_TOTAL),
      values_fill = 0L
    )
  # 確保 N_EVENT_TRUE / N_EVENT_FALSE / N_TOTAL_TRUE / N_TOTAL_FALSE 都存在
  for (col in c("N_EVENT_TRUE", "N_EVENT_FALSE", "N_TOTAL_TRUE", "N_TOTAL_FALSE")) {
    if (!col %in% names(wide)) wide[[col]] <- 0L
  }
  wide |>
    rename(STRATUM = .strat) |>
    mutate(
      RISK_E   = if_else(N_TOTAL_TRUE  > 0, N_EVENT_TRUE  / N_TOTAL_TRUE,  NA_real_),
      RISK_U   = if_else(N_TOTAL_FALSE > 0, N_EVENT_FALSE / N_TOTAL_FALSE, NA_real_),
      RR       = if_else(N_EVENT_FALSE > 0 & N_EVENT_TRUE > 0,
                         RISK_E / RISK_U, NA_real_),
      LOG_SE   = if_else(N_EVENT_FALSE > 0 & N_EVENT_TRUE > 0,
                         sqrt(1/N_EVENT_TRUE + 1/N_EVENT_FALSE
                              - 1/N_TOTAL_TRUE - 1/N_TOTAL_FALSE),
                         NA_real_),
      RR_LCI   = exp(log(RR) - 1.96 * LOG_SE),
      RR_UCI   = exp(log(RR) + 1.96 * LOG_SE),
      STRATUM_VAR = label
    ) |>
    select(STRATUM_VAR, STRATUM, N_EVENT_TRUE, N_TOTAL_TRUE,
           N_EVENT_FALSE, N_TOTAL_FALSE,
           RISK_E, RISK_U, RR, RR_LCI, RR_UCI)
}

# G-1 CT_SETTING 衍生 helper：把 case_type 對到 ER vs Clinic
# 健保就醫類別代碼（OPDTE.case_type）：
#   "02" = 急診 → ER
#   "01" = 一般門診（含特約門診）→ Clinic
#   "03" = 連續處方箋 → Clinic
#   其他罕見代碼一律歸 Clinic（保守處理；若實際資料有意外值會在 build log 印出）
classify_ct_setting <- function(case_type) {
  ct <- as.character(case_type)
  ct <- trimws(ct)
  dplyr::case_when(
    is.na(ct) | ct == "" ~ NA_character_,
    ct == "02"           ~ "ER",
    TRUE                 ~ "Clinic"
  )
}

# G-3 年齡分組 helper：給定一個年齡 vector 與一組 break，回傳 character factor
# - first CT 與 age at index 用 [0,5] (6,10] (10,15] (15,18]
# - age at diagnosis 用 [0,5] (5,10] (10,15] (15,Inf) 因為診斷可能晚於 18 歲
make_age_group <- function(age, breaks = c(0, 5, 10, 15, 18),
                           labels = NULL, right = TRUE) {
  if (is.null(labels)) {
    labels <- paste0(head(breaks, -1) + as.integer(right & head(breaks, -1) > 0),
                     "-",
                     tail(breaks, -1))
  }
  out <- cut(age, breaks = breaks, labels = labels,
             include.lowest = TRUE, right = right)
  as.character(out)
}

# G-3 multiple-CT 分組（暴露組 N_CT_TOTAL 的次數分箱）
# 0 → "0"（理論上 cohort 內不會出現，因為 0 等於 unexposed）
# 1 → "1"；2-3 → "2-3"；>=4 → ">=4"
make_n_ct_group <- function(n_ct) {
  dplyr::case_when(
    is.na(n_ct) | n_ct == 0L ~ NA_character_,
    n_ct == 1L               ~ "1",
    n_ct <= 3L               ~ "2-3",
    TRUE                     ~ ">=4"
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

# 主分析 eligible_ids：研究期間 2016-2023，0-18 歲，出生年 1997-2023
eligible_ids <- pers_info_with_birth |>
  filter(
    !is.na(BIRTH_DATE),
    BIRTH_DATE <= STUDY_END,
    BIRTH_DATE + years(MAX_AGE_AT_INDEX + 1L) > STUDY_START
  ) |>
  pull(ID)

log_progress(sprintf("eligible_ids 筆數 = %s（2016-2023 主分析）",
                     format(length(eligible_ids), big.mark = ",")))
if (length(eligible_ids) == 0L) {
  stop("eligible_ids 為 0 — 請檢查 pers_info$ID_BIRTHYM 的內容")
}

# H-4：擴展 eligible_ids（2000-2023 敏感性分析用）
# 出生年範圍 1982-2023，讓 1982 出生者在 2000-01-01 為 18 歲
eligible_ids_extended <- pers_info_with_birth |>
  filter(
    !is.na(BIRTH_DATE),
    BIRTH_DATE <= STUDY_END,
    BIRTH_DATE + years(MAX_AGE_AT_INDEX + 1L) > EXTENDED_STUDY_START
  ) |>
  pull(ID)

log_progress(sprintf("eligible_ids_extended 筆數 = %s（2000-2023 H-4 敏感性）",
                     format(length(eligible_ids_extended), big.mark = ",")))

## 2.2 CT Scan Records (OPDTO only; IPDTO 不可用) -----------------------------
# 【建檔步驟:對齊 03i 風格已註解】
# ⚠ 資料來源:
#   - 門診:OPDTO(含 drug_no)join OPDTE(含 id, func_date, case_type)
#   - 本專案無 IPDTO,不抓住院 CT 醫令
#   - OPDTE 月份格式:H_NHI_OPDTE{roc_yr}{mm}_10.parquet
# ⚠ 二次跑分析時直接 read_rds 即可;要重建檔請取消下方註解區塊
#
# 【G-1】新增 CT_SETTING 欄位（ER vs Clinic）
#   OPDTE.case_type 為健保就醫類別代碼：
#     "02" = 急診 → ER
#     "01"/"03"/其他 = 一般門診/連續處方箋等 → Clinic
#   保留 CASE_TYPE 原值供日後驗證
#
# 【H-4】CT 紀錄 build 範圍從 2016-2023 擴展為 2000-2023
#   主分析仍只用 2016 起的 CT，但讓 ct-records-all.rds 包含全期間，
#   方便 H-4 敏感性 cohort 直接讀同一個檔。
#   id 篩選用 eligible_ids_extended（涵蓋 1982-2023 出生者）。

# log_progress("===== 建檔 2.2:門診 CT 紀錄(OPDTO join OPDTE, 2000-2023)=====")
# 
# CT_BUILD_YEARS <- 2000:2023
# N_OP <- length(CT_BUILD_YEARS) * 12L
# ct_op_list <- vector("list", N_OP)
# idx <- 1L
# for (y in CT_BUILD_YEARS) {
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
#     # G-1：select 加上 case_type
#     dt_e <- open_dataset(file_e) |>
#       rename_with(tolower) |>
#       select(fee_ym, appl_date, appl_type, case_type, seq_no, hosp_id,
#              id, func_date) |>
#       filter(id %in% eligible_ids_extended) |>
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
#     ),
#     # G-1：CT_SETTING 在 collect 後依 case_type 衍生
#     #     case_type 在 join_keys 中已經 normalise（toupper + trimws）
#     ct_setting = classify_ct_setting(case_type)
#   ) |>
#   select(id, func_date, order_code, source, ct_type, ct_setting, case_type) |>
#   rename(ID = id, FUNC_DATE = func_date, ORDER_CODE = order_code,
#          SOURCE = source, CT_TYPE = ct_type, CT_SETTING = ct_setting,
#          CASE_TYPE = case_type)
# 
# # G-1 診斷：確認 case_type 分布是否符合預期（"02" 急診 vs 其他）
# log_progress("===== §2.2 case_type 分布診斷 =====")
# print(ct_records_raw |> count(CASE_TYPE) |> arrange(desc(n)))
# log_progress("===== §2.2 CT_SETTING 分布診斷 =====")
# print(ct_records_raw |> count(CT_SETTING) |> arrange(desc(n)))
# log_progress("===== §2.2 CT_TYPE × CT_SETTING 交叉診斷 =====")
# print(ct_records_raw |> count(CT_TYPE, CT_SETTING))
# 
# log_progress(sprintf("CT 紀錄合併完成（2000-2023）:%s 筆",
#                      format(nrow(ct_records_raw), big.mark = ",")))
# 
# # H-4 診斷：CT 紀錄按年分布（看 2000-2015 的 CT 量是否合理）
# log_progress("===== §2.2 各年 ORDER_CODE 數量診斷（H-4） =====")
# print(ct_records_raw |> mutate(YEAR = year(FUNC_DATE)) |>
#         count(YEAR, ORDER_CODE) |>
#         pivot_wider(names_from = ORDER_CODE, values_from = n,
#                     values_fill = 0L),
#       n = Inf)
# 
# write_rds(ct_records_raw, file.path(INTERMEDIATE_PATH, "ct-records-all.rds"))
# log_progress("已儲存 ct-records-all.rds")

# 改回 read_rds 模式（不重建）
ct_records_raw <- read_rds(file.path(INTERMEDIATE_PATH, "ct-records-all.rds"))
log_progress(sprintf("讀入 ct-records-all.rds：%s 筆",
                     format(nrow(ct_records_raw), big.mark = ",")))

ct_records <- ct_records_raw

# 安全網：若程式上方未執行 build（例如使用者手動把 §2.2 註解掉重跑分析），
# 而舊版 rds 沒有 CT_SETTING 欄位，至少不要讓下游分析掛掉
if (!"CT_SETTING" %in% names(ct_records)) {
  warning("ct-records-all.rds 沒有 CT_SETTING 欄位（舊版 v03f 之前建立）。",
          "請取消 §2.2 註解重建以啟用 ER/Clinic 分層分析。",
          "目前先以 NA 填入。")
  ct_records$CT_SETTING <- NA_character_
  ct_records$CASE_TYPE  <- NA_character_
}

log_progress(sprintf("ct_records 就緒:%s 筆", format(nrow(ct_records), big.mark = ",")))

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
#       filter(id %in% eligible_ids_extended) |>
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
#     filter(id %in% eligible_ids_extended) |>
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
#     # H-5：ICD-9 與 ICD-10 都用對應 helper 分 11 類
#     #     原版 03f：CODE_VERSION == "ICD-9" → "icd9_unclassified"（無法做 by-type）
#     #     03g 起：用 classify_malignancy_icd9() / icd10() 各自分類
#     MALIGNANCY_TYPE = if_else(
#       CODE_VERSION == "ICD-9",
#       classify_malignancy_icd9(ICD3),
#       classify_malignancy_icd10(ICD3)
#     )
#   )
# 
# log_progress(sprintf("惡性腫瘤首次診斷:%s 人(含 ICD-9 段 %s 人,ICD-10 段 %s 人)",
#                      format(nrow(malignancy_raw), big.mark = ","),
#                      format(sum(malignancy_raw$CODE_VERSION == "ICD-9"), big.mark = ","),
#                      format(sum(malignancy_raw$CODE_VERSION == "ICD-10"), big.mark = ",")))
# # H-5 診斷：印出 ICD-9 段的 11 類分布，確認分類是否合理
# log_progress("===== §2.3 ICD-9 段 MALIGNANCY_TYPE 分布診斷（H-5） =====")
# print(malignancy_raw |> filter(CODE_VERSION == "ICD-9") |> count(MALIGNANCY_TYPE))
# log_progress("===== §2.3 ICD-10 段 MALIGNANCY_TYPE 分布診斷 =====")
# print(malignancy_raw |> filter(CODE_VERSION == "ICD-10") |> count(MALIGNANCY_TYPE))
# write_rds(malignancy_raw, file.path(INTERMEDIATE_PATH, "malignancy-dx.rds"))
# log_progress("已儲存 malignancy-dx.rds")

# 改回 read_rds 模式（不重建）
# ⚠ 注意：舊 rds 若仍以 "icd9_unclassified" 標記 ICD-9 紀錄，
#   下方 safety net 會用 H-5 helper 即時補修分類
malignancy_raw <- read_rds(file.path(INTERMEDIATE_PATH, "malignancy-dx.rds"))
log_progress(sprintf("讀入 malignancy-dx.rds：%s 人",
                     format(nrow(malignancy_raw), big.mark = ",")))

malignancy_dx <- malignancy_raw

# Safety net：若 §2.3 build 區塊被使用者手動註解、改用舊 rds，
# 而舊 rds 的 ICD-9 紀錄都標 "icd9_unclassified"，這裡即時補修一次
n_unclass_before <- sum(malignancy_dx$MALIGNANCY_TYPE == "icd9_unclassified",
                        na.rm = TRUE)
if (n_unclass_before > 0L) {
  log_progress(sprintf("發現 %s 筆 'icd9_unclassified' 紀錄，套用 H-5 helper 即時補修",
                       format(n_unclass_before, big.mark = ",")))
  malignancy_dx <- malignancy_dx |>
    mutate(
      MALIGNANCY_TYPE = if_else(
        MALIGNANCY_TYPE == "icd9_unclassified" & CODE_VERSION == "ICD-9",
        classify_malignancy_icd9(ICD3),
        MALIGNANCY_TYPE
      )
    )
  log_progress("===== H-5 補修後 ICD-9 段 MALIGNANCY_TYPE 分布 =====")
  print(malignancy_dx |> filter(CODE_VERSION == "ICD-9") |> count(MALIGNANCY_TYPE))
}
log_progress(sprintf("malignancy_dx 就緒:%s 人", format(nrow(malignancy_dx), big.mark = ",")))

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
# 【建檔步驟：已註解（前次跑完已寫出 rds），預設僅 read_rds】
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
# ⚠ 要重建檔請取消本區塊原 build code（已從本檔移除以縮短行數）；
#    歷史版本見 git log 中的 v01b 之前版本，或對照 v03f §2.5。
#
# 改回 read_rds 模式（不重建）
appendicitis_appendectomy_full <- read_rds(
  file.path(INTERMEDIATE_PATH, "appendicitis-appendectomy.rds")
)
appendicitis_appendectomy_main <- read_rds(
  file.path(INTERMEDIATE_PATH, "appendicitis-appendectomy-main.rds")
)
log_progress(sprintf("讀入 appendicitis-appendectomy.rds：%s 人（A∪B 聯集）",
                     format(nrow(appendicitis_appendectomy_full), big.mark = ",")))
log_progress(sprintf("讀入 appendicitis-appendectomy-main.rds：%s 人（僅 Group A，主分析）",
                     format(nrow(appendicitis_appendectomy_main), big.mark = ",")))

# 後續 §6 主分析使用 main；敏感性分析使用 full
appendicitis_appendectomy      <- appendicitis_appendectomy_main
appendicitis_appendectomy_sens <- appendicitis_appendectomy_full
## 2.6 Sibling Relationship File (Study 2) ------------------------------------
# 【建檔步驟：已註解（前次跑完已寫出 rds），預設僅 read_rds】
# ⚠ 要重建檔：取消下方註解的 build 區塊
#
# log_progress("===== 建檔 2.6:手足配對關係 =====")
# enrol_relation_ext <- open_dataset(
#   file.path(PROCESSED_DATA_PATH, "enrol-relation-extend.parquet")
# ) |>
#   select(ID, F_ID, M_ID, Guess_F, Guess_M) |>
#   collect() |>
#   mutate(
#     EFF_F = coalesce(as.character(F_ID), as.character(Guess_F)),
#     EFF_M = coalesce(as.character(M_ID), as.character(Guess_M))
#   ) |>
#   select(ID, EFF_F, EFF_M)
# write_rds(enrol_relation_ext, file.path(INTERMEDIATE_PATH, "enrol-relation-ext.rds"))
# log_progress("已儲存 enrol-relation-ext.rds")

# 改回 read_rds 模式（不重建）
enrol_relation_ext <- read_rds(file.path(INTERMEDIATE_PATH, "enrol-relation-ext.rds"))
log_progress(sprintf("讀入 enrol-relation-ext.rds：%s 筆",
                     format(nrow(enrol_relation_ext), big.mark = ",")))
## 2.7 Annual Enrollment: Urban Status & Income Quartile (記憶體優化版) ----
# 【建檔步驟：實際執行（這次跑的目的就是要重建此檔）】
#
# 改寫重點（相對 v01b 原版）：
#   ① Step 2 的 id1_amt / id_status / id_roc filter 全部推回 Arrow 階段
#   ② Step 1 與 Step 2 共用同一個 open_dataset 物件
#   ③ 年內就 bind_rows + gc，避免堆 288 個 list 元素到結尾才 bind
#   ④ Step 2 用 relevant_id1 把 Arrow query 進一步縮小
#   ⑤ 拆成 2000-2015 / 2016-2023 兩個 chunk，分別 write_rds 後合併
#
# ⚠ 官方手冊 H_NHI_ENROL 真實欄位（小寫對齊 00-basic_rate.R）：
#     id, id1, id1_amt, id1_city, id_status, id_roc, prem_ym
# ⚠ ENROL 月檔無 _10 sample group 後綴，find_parquet 必須傳 suffix = ""
# ⚠ AMT_CUTOFF = 2000 點未依年度通膨調整；INCOME_Q 為「該年內相對排名」，
#   跨年比較需謹慎（在報告 limitations 寫清楚）

log_progress("===== 建檔 2.7：承保檔 income/urban（記憶體優化版 + 分 chunk） =====")

# H-7 範圍：2000-2023（為 H-4 extended cohort 服務）
# 拆成兩段：
#   chunk_history = 2000-2015（歷史回查段，僅 H-4 敏感性用）
#   chunk_main    = 2016-2023（主分析段；對應原 03f rds 範圍）
ENROL_CHUNKS <- list(
  history = 2000:2015,
  main    = 2016:2023
)

# Step 1（小孩 → 主投保者 id1）和 Step 2（主投保者本人薪資/縣市）的欄位集合
# 小寫對齊手冊真實欄位名（F-6）
child_select_cols <- c("id", "id1", "prem_ym", "id_roc")
ins_select_cols   <- c("id", "id1", "id1_amt", "id1_city", "prem_ym",
                       "id_status", "id_roc")

# Pre-flight：每年的存在月數摘要（替代原版 288 行 detail）
log_progress("===== §2.7 pre-flight：ENROL 各年月檔存在數 =====")
preflight_summary <- list()
first_file_names_printed <- FALSE
for (y in unlist(ENROL_CHUNKS)) {
  n_exist <- 0L
  for (m in 1:12) {
    file_en <- find_parquet("ENROL", y, m, suffix = "")
    if (!is.na(file_en)) {
      n_exist <- n_exist + 1L
      if (!first_file_names_printed) {
        log_progress(sprintf("第一個找到的 ENROL 檔：%d-%02d", y, m))
        log_progress(sprintf("  欄位名：%s",
                             paste(names(open_dataset(file_en) |>
                                           rename_with(tolower)),
                                   collapse = ", ")))
        first_file_names_printed <- TRUE
      }
    }
  }
  preflight_summary[[as.character(y)]] <- n_exist
}
preflight_df <- tibble(
  year     = as.integer(names(preflight_summary)),
  n_months = as.integer(unlist(preflight_summary))
)
print(preflight_df, n = Inf)
log_progress(sprintf("ENROL 月檔總存在數：%d / %d",
                     sum(preflight_df$n_months), nrow(preflight_df) * 12L))

# ──────────────────────────────────────────────────────────────────────
# 月內處理函式：給定一個月 (y, m)，回傳 list(child = ..., ins = ...)
#   - child：該月 eligible 小孩 → 主投保者 id1 對應
#   - ins  ：該月主投保者本人薪資與縣市（已 Arrow 端篩 amt/status/roc）
# 兩個物件都已縮到「該月有意義」的最小集合
# ──────────────────────────────────────────────────────────────────────
process_one_month <- function(y, m) {
  file_en <- find_parquet("ENROL", y, m, suffix = "")
  if (is.na(file_en)) return(list(child = NULL, ins = NULL))
  
  ds <- open_dataset(file_en) |> rename_with(tolower)
  
  # Step 1：小孩 → 主投保者 id1
  #   Arrow 階段：select + filter(id %in% eligible_ids_extended) + filter(id_roc == "0")
  #   注意：id_roc == "0" 是字串常數比較，arrow 完全支援，可以下推
  child_raw <- tryCatch(
    ds |>
      select(any_of(child_select_cols)) |>
      filter(id %in% eligible_ids_extended,
             id_roc == "0") |>
      collect(),
    error = function(e) {
      log_progress(sprintf("  ⚠ Step1 collect 失敗 (%d-%02d)：%s",
                           y, m, conditionMessage(e)))
      tibble()
    }
  )
  
  if (nrow(child_raw) == 0L) {
    return(list(child = NULL, ins = NULL))
  }
  
  child_clean <- child_raw |>
    mutate(
      id      = as.character(id),
      id1     = as.character(id1),
      prem_ym = as.character(prem_ym),
      year    = suppressWarnings(as.integer(substr(prem_ym, 1, 4)))
    ) |>
    filter(year == y) |>
    distinct(id, id1, year, prem_ym)
  
  if (nrow(child_clean) == 0L) {
    return(list(child = NULL, ins = NULL))
  }
  
  # Step 2：主投保者本人薪資與縣市
  #   Arrow 階段：把所有常數比較的 filter 都推下去
  #   - id == id1（主投保者本人）
  #   - id %in% relevant_id1（縮到 Step 1 出現過的家戶主）
  #   - !is.na(id1_amt) & id1_amt >= AMT_CUTOFF
  #   - id_status %in% c("1","2","3")
  #   - id_roc == "0"
  relevant_id1 <- unique(child_clean$id1)
  
  ins_raw <- tryCatch(
    ds |>
      select(any_of(ins_select_cols)) |>
      filter(
        id == id1,
        id %in% relevant_id1,
        !is.na(id1_amt),
        id1_amt   >= AMT_CUTOFF,
        id_status %in% c("1", "2", "3"),
        id_roc    == "0"
      ) |>
      collect(),
    error = function(e) {
      log_progress(sprintf("  ⚠ Step2 collect 失敗 (%d-%02d)：%s",
                           y, m, conditionMessage(e)))
      tibble()
    }
  )
  
  if (nrow(ins_raw) == 0L) {
    # Step1 有結果但 Step2 沒有：保留 child（後續 inner_join 自然會丟掉這個月）
    return(list(child = child_clean, ins = NULL))
  }
  
  ins_clean <- ins_raw |>
    mutate(
      id        = as.character(id),
      id_status = as.character(id_status),
      id_roc    = as.character(id_roc),
      id1_city  = as.character(id1_city),
      prem_ym   = as.character(prem_ym),
      year      = suppressWarnings(as.integer(substr(prem_ym, 1, 4)))
    ) |>
    filter(year == y) |>
    transmute(id1 = id, id1_amt, id1_city, prem_ym, year)
  
  list(child = child_clean, ins = ins_clean)
}

# ──────────────────────────────────────────────────────────────────────
# Chunk 處理函式：處理一個年度區段（如 2000:2015 或 2016:2023）
#   回傳：一個 list(child = combined_tbl, ins = combined_tbl)
# 每年結束就 bind_rows + gc，避免堆過多月份 list
# ──────────────────────────────────────────────────────────────────────
process_chunk <- function(years_in_chunk, chunk_label) {
  log_progress(sprintf(">>> 處理 chunk: %s (years %d-%d)",
                       chunk_label, min(years_in_chunk), max(years_in_chunk)))
  
  child_by_year <- vector("list", length(years_in_chunk))
  ins_by_year   <- vector("list", length(years_in_chunk))
  
  for (yi in seq_along(years_in_chunk)) {
    y <- years_in_chunk[yi]
    log_progress(sprintf("  [Enrol] year %d", y))
    
    child_month <- vector("list", 12)
    ins_month   <- vector("list", 12)
    
    for (m in 1:12) {
      res <- process_one_month(y, m)
      if (!is.null(res$child)) child_month[[m]] <- res$child
      if (!is.null(res$ins))   ins_month[[m]]   <- res$ins
    }
    
    # 年內合併並縮小
    child_year <- bind_rows(child_month) |>
      distinct(id, id1, year, prem_ym)
    ins_year   <- bind_rows(ins_month)
    
    child_by_year[[yi]] <- child_year
    ins_by_year[[yi]]   <- ins_year
    
    log_progress(sprintf("    年內統計 %d：child %s 筆 / ins %s 筆",
                         y,
                         format(nrow(child_year), big.mark = ","),
                         format(nrow(ins_year),   big.mark = ",")))
    
    rm(child_month, ins_month, child_year, ins_year)
    gc(verbose = FALSE)
  }
  
  list(
    child = bind_rows(child_by_year),
    ins   = bind_rows(ins_by_year)
  )
}

# 跑兩個 chunk
chunk_results <- list()
for (cname in names(ENROL_CHUNKS)) {
  chunk_results[[cname]] <- process_chunk(ENROL_CHUNKS[[cname]], cname)
  log_progress(sprintf("Chunk %s 完成：child %s 筆 / ins %s 筆",
                       cname,
                       format(nrow(chunk_results[[cname]]$child), big.mark = ","),
                       format(nrow(chunk_results[[cname]]$ins),   big.mark = ",")))
  gc(verbose = FALSE)
}

# 合併兩個 chunk
child_ins_monthly <- bind_rows(
  chunk_results$history$child,
  chunk_results$main$child
)
ins_amt_monthly <- bind_rows(
  chunk_results$history$ins,
  chunk_results$main$ins
)

# 釋放 chunk_results 記憶體
rm(chunk_results); gc(verbose = FALSE)

log_progress(sprintf("child_ins_monthly：%s 筆（小孩-月-主投保者）",
                     format(nrow(child_ins_monthly), big.mark = ",")))
log_progress(sprintf("ins_amt_monthly  ：%s 筆（主投保者-月-薪資）",
                     format(nrow(ins_amt_monthly), big.mark = ",")))

# 型別正規化（保險起見）
child_ins_monthly <- child_ins_monthly |>
  mutate(
    id1     = as.character(id1),
    prem_ym = as.character(prem_ym),
    year    = as.integer(year)
  )
ins_amt_monthly <- ins_amt_monthly |>
  mutate(
    id1     = as.character(id1),
    prem_ym = as.character(prem_ym),
    year    = as.integer(year)
  )

# ──────────────────────────────────────────────────────────────────────
# 小孩 ↔ 家戶主 月度 inner_join → 推導年度薪資 / 主要縣市 / 四分位
# ──────────────────────────────────────────────────────────────────────
household_monthly <- child_ins_monthly |>
  inner_join(ins_amt_monthly, by = c("id1", "year", "prem_ym"))

log_progress(sprintf("household_monthly：%s 筆（小孩-月，含家戶薪資/縣市）",
                     format(nrow(household_monthly), big.mark = ",")))

# 釋放兩個原始 monthly（合併後不再需要）
rm(child_ins_monthly, ins_amt_monthly); gc(verbose = FALSE)

annual_amt <- household_monthly |>
  group_by(id, year) |>
  summarise(
    TOTAL_AMT = sum(id1_amt, na.rm = TRUE),
    N_MONTHS  = dplyr::n(),
    .groups   = "drop"
  ) |>
  rename(ID = id, YEAR = year)

annual_city <- household_monthly |>
  filter(!is.na(id1_city), id1_city != "") |>
  group_by(id, year, id1_city) |>
  summarise(N_MONTHS_CITY = dplyr::n(), .groups = "drop") |>
  group_by(id, year) |>
  arrange(desc(N_MONTHS_CITY), id1_city, .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  select(ID = id, YEAR = year, MAIN_CITY = id1_city)

# household_monthly 用完了
rm(household_monthly); gc(verbose = FALSE)

# 年度四分位（同年內相對排名）
annual_amt_quartile <- annual_amt |>
  group_by(YEAR) |>
  mutate(
    INCOME_Q = cut(
      TOTAL_AMT,
      breaks         = quantile(TOTAL_AMT, probs = c(0, 0.25, 0.5, 0.75, 1),
                                na.rm = TRUE, type = 7),
      labels         = c("Q1", "Q2", "Q3", "Q4"),
      include.lowest = TRUE
    )
  ) |>
  ungroup()

enrol_income_urban_raw <- annual_amt_quartile |>
  left_join(annual_city, by = c("ID", "YEAR")) |>
  mutate(URBAN_STATUS = classify_urban(MAIN_CITY)) |>
  select(ID, YEAR, TOTAL_AMT, N_MONTHS, INCOME_Q, MAIN_CITY, URBAN_STATUS)

log_progress(sprintf("enrol_income_urban_raw：%s 筆（人-年，跨 2000-2023）",
                     format(nrow(enrol_income_urban_raw), big.mark = ",")))

# 各年覆蓋人數診斷
log_progress("===== §2.7 各年覆蓋人數診斷 =====")
print(enrol_income_urban_raw |>
        group_by(YEAR) |>
        summarise(N_persons     = dplyr::n_distinct(ID),
                  N_with_income = sum(!is.na(INCOME_Q)),
                  N_with_urban  = sum(!is.na(URBAN_STATUS)),
                  .groups = "drop"),
      n = Inf)

n_eligible_with_income <- enrol_income_urban_raw |>
  filter(!is.na(INCOME_Q)) |>
  distinct(ID) |>
  nrow()
log_progress(sprintf("eligible_ids_extended 中取得 INCOME_Q 的人數：%s / %s (%.1f%%)",
                     format(n_eligible_with_income, big.mark = ","),
                     format(length(eligible_ids_extended), big.mark = ","),
                     100 * n_eligible_with_income / length(eligible_ids_extended)))

write_rds(enrol_income_urban_raw,
          file.path(INTERMEDIATE_PATH, "enrol-income-urban.rds"))
log_progress("已儲存 enrol-income-urban.rds")

enrol_income_urban <- enrol_income_urban_raw
log_progress(sprintf("enrol_income_urban 就緒：%s 筆（人-年）",
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
# G-2：新增 FIRST_CT_SETTING（同日多筆時 ER 優先）
# H-1：保留 N_CT_TOTAL，後續做 dose-response（baseline N_CT_GROUP）副版

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
    # FIRST_CT_TYPE：取同人最早一筆（同日按 ORDER_CODE 字典序）
    ct_records |>
      filter(FUNC_DATE >= STUDY_START, FUNC_DATE <= STUDY_END) |>
      group_by(ID) |>
      arrange(FUNC_DATE, ORDER_CODE, .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      select(ID, FIRST_CT_TYPE = CT_TYPE),
    by = "ID"
  ) |>
  left_join(
    # G-2：FIRST_CT_SETTING — 取同人最早 CT 日期那筆的 setting
    #     同日多筆時 ER 優先（ER 在 arrange 時排在前 — 因為 "Clinic" > "ER" 字典序，
    #     需要明確 desc 或用 case_when 處理；此處用 factor level 強制 ER 優先）
    ct_records |>
      filter(FUNC_DATE >= STUDY_START, FUNC_DATE <= STUDY_END) |>
      mutate(SETTING_PRIO = if_else(CT_SETTING == "ER", 1L, 2L)) |>
      group_by(ID) |>
      arrange(FUNC_DATE, SETTING_PRIO, ORDER_CODE, .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      select(ID, FIRST_CT_SETTING = CT_SETTING),
    by = "ID"
  ) |>
  mutate(
    # H-1 副版：dose-response 用 N_CT_GROUP（1 / 2-3 / >=4）
    N_CT_GROUP = make_n_ct_group(N_CT_TOTAL)
  )

## 3.3 Prior Malignancy Exclusion Flag ----------------------------------------
# H-8：改用 effective start date 為界
#   暴露者：effective_start = max(BIRTH_DATE, STUDY_START, FIRST_CT_DATE)
#           → 排除「first CT 之前」的 malignancy
#   未暴露者：effective_start = max(BIRTH_DATE, STUDY_START)
#           → 排除「STUDY_START 之前」的 malignancy
#
# 理由：CT 暴露之前已有 cancer 不該歸因為「兒童期 CT 暴露之後」事件。
# 主分析 cohort 因 entry = 2016-01-01、CT 也在 2016 之後，
# 通常 first_CT > entry → 暴露者排除範圍會比舊邏輯多一些 lead-time

prior_malignancy_ids <- base_cohort |>
  left_join(ct_summary |> select(ID, FIRST_CT_DATE), by = "ID") |>
  mutate(
    EFFECTIVE_START = pmax(
      BIRTH_DATE,
      STUDY_START,
      coalesce(FIRST_CT_DATE, STUDY_START),
      na.rm = TRUE
    )
  ) |>
  inner_join(malignancy_dx |> select(ID, FIRST_MALIGNANCY_DATE),
             by = "ID") |>
  filter(FIRST_MALIGNANCY_DATE < EFFECTIVE_START) |>
  distinct(ID)

log_progress(sprintf("§3.3 prior malignancy 排除人數（H-8 effective start date）：%s",
                     format(nrow(prior_malignancy_ids), big.mark = ",")))

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

# G-4：統一描述統計 helper（三個 study 共用）
# 輸入：cohort（已 attach income/urban、已有 OUTCOME / TIME_MONTHS / SEX 等欄位）
# 輸出：long-format tibble，每個變項一個 row block；按 EXPOSED 欄展開
#
# 欄位 schema：
#   VARIABLE       變項名稱（例：N、SEX、AGE_AT_INDEX、INCOME_Q）
#   CATEGORY       次類別（例：Male、Female、Q1、Metro）
#   EXPOSED_TRUE   暴露組數值
#   EXPOSED_FALSE  未暴露組數值
#   TOTAL          全 cohort 數值
#   STAT_TYPE      "n" / "n_pct" / "median_iqr" / "mean_sd"
#
# 參數 age_var：個人年齡變項名稱（study1/2 = "AGE_AT_INDEX"，study3 = "AGE_AT_APPENDIX"）
make_descriptive <- function(cohort, age_var = "AGE_AT_INDEX",
                             age_at_2016_var = "AGE_AT_2016") {
  
  total_n <- nrow(cohort)
  exp_n   <- sum(cohort$EXPOSED, na.rm = TRUE)
  unx_n   <- sum(!cohort$EXPOSED, na.rm = TRUE)
  
  # --- block 1: N ---
  blk_n <- tibble(
    VARIABLE = "N", CATEGORY = "",
    EXPOSED_TRUE  = as.character(exp_n),
    EXPOSED_FALSE = as.character(unx_n),
    TOTAL         = as.character(total_n),
    STAT_TYPE     = "n"
  )
  
  # 用於 % 顯示
  fmt_n_pct <- function(n_e, n_u, n_t) {
    list(
      e = sprintf("%s (%.1f%%)", format(n_e, big.mark = ","),
                  100 * n_e / max(exp_n, 1)),
      u = sprintf("%s (%.1f%%)", format(n_u, big.mark = ","),
                  100 * n_u / max(unx_n, 1)),
      t = sprintf("%s (%.1f%%)", format(n_t, big.mark = ","),
                  100 * n_t / max(total_n, 1))
    )
  }
  
  # 通用：依某 categorical 欄位算 EXPOSED × CATEGORY 表
  block_categorical <- function(var_name, df = cohort) {
    if (!var_name %in% names(df)) return(NULL)
    grp <- df |>
      mutate(.cat = as.character(.data[[var_name]])) |>
      mutate(.cat = if_else(is.na(.cat), "Missing", .cat)) |>
      group_by(.cat, EXPOSED) |>
      summarise(n = dplyr::n(), .groups = "drop") |>
      pivot_wider(names_from = EXPOSED, values_from = n,
                  names_prefix = "n_", values_fill = 0L)
    
    # 確保有 n_TRUE / n_FALSE 兩欄
    if (!"n_TRUE"  %in% names(grp)) grp$n_TRUE  <- 0L
    if (!"n_FALSE" %in% names(grp)) grp$n_FALSE <- 0L
    
    grp |>
      mutate(n_total = n_TRUE + n_FALSE) |>
      arrange(.cat) |>
      mutate(
        VARIABLE = var_name,
        CATEGORY = .cat,
        EXPOSED_TRUE  = sprintf("%s (%.1f%%)", format(n_TRUE, big.mark = ","),
                                100 * n_TRUE / max(exp_n, 1)),
        EXPOSED_FALSE = sprintf("%s (%.1f%%)", format(n_FALSE, big.mark = ","),
                                100 * n_FALSE / max(unx_n, 1)),
        TOTAL         = sprintf("%s (%.1f%%)", format(n_total, big.mark = ","),
                                100 * n_total / max(total_n, 1)),
        STAT_TYPE     = "n_pct"
      ) |>
      select(VARIABLE, CATEGORY, EXPOSED_TRUE, EXPOSED_FALSE, TOTAL, STAT_TYPE)
  }
  
  # 通用：依某 numeric 欄位算 median (IQR)
  block_numeric_median <- function(var_name, label = var_name, df = cohort) {
    if (!var_name %in% names(df)) return(NULL)
    f <- function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return("NA")
      sprintf("%.1f (%.1f-%.1f)",
              stats::median(x),
              stats::quantile(x, 0.25),
              stats::quantile(x, 0.75))
    }
    tibble(
      VARIABLE = label, CATEGORY = "median (IQR)",
      EXPOSED_TRUE  = f(df[[var_name]][df$EXPOSED]),
      EXPOSED_FALSE = f(df[[var_name]][!df$EXPOSED]),
      TOTAL         = f(df[[var_name]]),
      STAT_TYPE     = "median_iqr"
    )
  }
  
  # block 2-9：各變項
  blk_sex      <- block_categorical("SEX")
  blk_ageidx   <- block_numeric_median(age_var,
                                       label = paste0(age_var, " (years)"))
  blk_age2016  <- if (age_at_2016_var %in% names(cohort))
    block_numeric_median(age_at_2016_var,
                         label = paste0(age_at_2016_var, " (years)"))
  else NULL
  blk_first_ct_age <- block_numeric_median("FIRST_CT_AGE",
                                           label = "Age at first CT (years, exposed only)")
  blk_first_ct_age_grp <- block_categorical("FIRST_CT_AGE_GROUP")
  blk_age_dx   <- block_numeric_median("AGE_AT_DIAGNOSIS",
                                       label = "Age at malignancy diagnosis (years, events only)")
  blk_age_dx_grp <- block_categorical("AGE_AT_DIAGNOSIS_GROUP")
  blk_income   <- block_categorical("INCOME_Q")
  blk_urban    <- block_categorical("URBAN_STATUS")
  blk_ct_set   <- block_categorical("FIRST_CT_SETTING")
  blk_ct_type  <- block_categorical("FIRST_CT_TYPE")
  blk_n_ct_grp <- block_categorical("N_CT_GROUP")
  
  # block 10: outcome 與 follow-up
  n_event_e <- sum(cohort$OUTCOME &  cohort$EXPOSED, na.rm = TRUE)
  n_event_u <- sum(cohort$OUTCOME & !cohort$EXPOSED, na.rm = TRUE)
  blk_event <- tibble(
    VARIABLE = "Malignancy events", CATEGORY = "",
    EXPOSED_TRUE  = sprintf("%s (%.2f%%)", format(n_event_e, big.mark = ","),
                            100 * n_event_e / max(exp_n, 1)),
    EXPOSED_FALSE = sprintf("%s (%.2f%%)", format(n_event_u, big.mark = ","),
                            100 * n_event_u / max(unx_n, 1)),
    TOTAL         = sprintf("%s (%.2f%%)",
                            format(n_event_e + n_event_u, big.mark = ","),
                            100 * (n_event_e + n_event_u) / max(total_n, 1)),
    STAT_TYPE     = "n_pct"
  )
  blk_pyears <- if ("TIME_MONTHS" %in% names(cohort)) {
    py_e <- sum(cohort$TIME_MONTHS[cohort$EXPOSED],  na.rm = TRUE) / 12
    py_u <- sum(cohort$TIME_MONTHS[!cohort$EXPOSED], na.rm = TRUE) / 12
    tibble(
      VARIABLE = "Person-years", CATEGORY = "",
      EXPOSED_TRUE  = format(round(py_e), big.mark = ","),
      EXPOSED_FALSE = format(round(py_u), big.mark = ","),
      TOTAL         = format(round(py_e + py_u), big.mark = ","),
      STAT_TYPE     = "n"
    )
  } else NULL
  
  bind_rows(
    blk_n,
    blk_sex,
    blk_ageidx,
    blk_age2016,
    blk_first_ct_age,
    blk_first_ct_age_grp,
    blk_age_dx,
    blk_age_dx_grp,
    blk_income,
    blk_urban,
    blk_ct_set,
    blk_ct_type,
    blk_n_ct_grp,
    blk_event,
    blk_pyears
  )
}

# H-1 helper：把固定 INDEX_DATE 的 cohort 切成 time-varying long format
# 採 Mathews 2013 模式 (i)：lag period 內算未暴露 person-years
#
# 輸入 cohort 必要欄位：
#   ID, BIRTH_DATE, EXPOSED, FIRST_CT_DATE, OUTCOME, ELIGIBLE_MALIGNANCY_DATE,
#   SEX, INCOME_Q, URBAN_STATUS, FIRST_CT_SETTING, FIRST_CT_TYPE,
#   AGE_AT_2016, FIRST_CT_AGE, MALIGNANCY_TYPE
#
# 輸出 long-format tibble：每人 1-2 列
#   未暴露者 1 列：[ENTRY, EXIT)
#   暴露者 2 列：
#     列1 = [ENTRY, FIRST_CT_DATE + lag) 期間 EXPOSED_TV = 0
#     列2 = [FIRST_CT_DATE + lag, EXIT)  期間 EXPOSED_TV = 1
#
# 參數：
#   study_start    研究起點（用於計算 entry date）
#   study_end      研究終點（censoring 點）
#   lag_days       lag period in days（預設 TV_LAG_YEARS * 365.25）
#
# 注意：
#   - 事件發生時間：若 OUTCOME=TRUE 且 ELIGIBLE_MALIGNANCY_DATE < EXIT，
#     EXIT 改設為 ELIGIBLE_MALIGNANCY_DATE
#   - 若 first_CT + lag >= EXIT，那個人在 lag 內就 censor 了，整段都算未暴露
#     （這就是 Mathews 模式 (i) 的精神）
#   - 時間單位：days from study_start
build_tv_cohort <- function(cohort,
                            study_start = STUDY_START,
                            study_end   = STUDY_END,
                            lag_days    = round(TV_LAG_YEARS * 365.25)) {
  
  # 把每個人的 entry / exit / event_date / transfer_date（若暴露）算好
  prep <- cohort |>
    mutate(
      ENTRY_DATE    = pmax(study_start, BIRTH_DATE),
      EVENT_DATE    = if_else(OUTCOME, ELIGIBLE_MALIGNANCY_DATE, NA_Date_),
      EXIT_DATE     = pmin(coalesce(EVENT_DATE, study_end), study_end),
      TRANSFER_DATE = if_else(EXPOSED, FIRST_CT_DATE + days(lag_days), NA_Date_),
      # 整數天數（自 study_start 起算）
      tstart_d  = as.integer(ENTRY_DATE - study_start),
      texit_d   = as.integer(EXIT_DATE  - study_start),
      ttrans_d  = as.integer(TRANSFER_DATE - study_start)
    ) |>
    filter(texit_d > tstart_d)  # 排除 entry >= exit 者
  
  # 未暴露者：1 列
  unexp <- prep |>
    filter(!EXPOSED) |>
    transmute(
      ID, BIRTH_DATE, SEX, INCOME_Q, URBAN_STATUS,
      FIRST_CT_SETTING, FIRST_CT_TYPE, FIRST_CT_AGE,
      MALIGNANCY_TYPE,
      AGE_AT_2016,
      tstart      = tstart_d,
      tstop       = texit_d,
      EXPOSED_TV  = 0L,
      OUTCOME_TV  = as.integer(OUTCOME)
    )
  
  # 暴露者：依 transfer_date 與 exit_date 的相對位置決定切幾段
  exp_data <- prep |> filter(EXPOSED)
  
  # 情況 A：transfer_date >= exit_date
  #   → 整段算未暴露（事件如果在 lag 內發生，OUTCOME 就在這列上）
  exp_lag_only <- exp_data |>
    filter(ttrans_d >= texit_d) |>
    transmute(
      ID, BIRTH_DATE, SEX, INCOME_Q, URBAN_STATUS,
      FIRST_CT_SETTING, FIRST_CT_TYPE, FIRST_CT_AGE,
      MALIGNANCY_TYPE,
      AGE_AT_2016,
      tstart      = tstart_d,
      tstop       = texit_d,
      EXPOSED_TV  = 0L,
      OUTCOME_TV  = as.integer(OUTCOME)
    )
  
  # 情況 B：transfer_date < exit_date
  #   → 切兩段：[entry, transfer) unexposed (event=0)；[transfer, exit) exposed (event=OUTCOME)
  exp_split <- exp_data |> filter(ttrans_d < texit_d)
  
  exp_unx_part <- exp_split |>
    transmute(
      ID, BIRTH_DATE, SEX, INCOME_Q, URBAN_STATUS,
      FIRST_CT_SETTING, FIRST_CT_TYPE, FIRST_CT_AGE,
      MALIGNANCY_TYPE,
      AGE_AT_2016,
      tstart      = tstart_d,
      # transfer date 之前算未暴露
      tstop       = pmax(ttrans_d, tstart_d + 1L),  # 確保 tstop > tstart
      EXPOSED_TV  = 0L,
      OUTCOME_TV  = 0L  # event 必然發生在 transfer 之後（因情況 B）
    )
  exp_exp_part <- exp_split |>
    transmute(
      ID, BIRTH_DATE, SEX, INCOME_Q, URBAN_STATUS,
      FIRST_CT_SETTING, FIRST_CT_TYPE, FIRST_CT_AGE,
      MALIGNANCY_TYPE,
      AGE_AT_2016,
      tstart      = pmax(ttrans_d, tstart_d + 1L),  # 與 unexp_part 對接
      tstop       = texit_d,
      EXPOSED_TV  = 1L,
      OUTCOME_TV  = as.integer(OUTCOME)
    )
  
  out <- bind_rows(unexp, exp_lag_only, exp_unx_part, exp_exp_part) |>
    filter(tstop > tstart) |>
    arrange(ID, tstart)
  
  out
}

# 4. Study 1: Population-Level Cohort -----------------------------------------
# 暴露：研究期間有任何 CT 掃描（0-18 歲）
# 對照：研究期間無 CT 掃描
# 結果：首次 CT 後 >2 年發生惡性腫瘤
# 控制：年齡、性別、薪資四分位、六都/非六都

## 4.1 Build Study 1 Cohort --------------------------------------------------

log_progress("===== 4.1 Study 1 cohort =====")

# H-2：把 study1 cohort 構建抽成 helper，方便不同 latency 的敏感性分析重複使用
# 也支援 type-specific latency（給每個人依其 MALIGNANCY_TYPE 套用不同 latency）
#
# 參數：
#   base_df        基礎 cohort（已排除遺傳/已知惡性者）
#   ct_summary_df  ct_summary 表（含 FIRST_CT_DATE / FIRST_CT_TYPE / FIRST_CT_SETTING / N_CT_GROUP）
#   malignancy_df  惡性腫瘤首次診斷表（含 MALIGNANCY_TYPE）
#   prior_mal_ids  先期已診斷惡性腫瘤者（要排除）
#   study_start    研究起點（主分析 = STUDY_START；H-4 = EXTENDED_STUDY_START）
#   study_end      研究終點
#   latency_years  number 或 NULL；
#                  若為 number → 所有人套同一 latency
#                  若為 NULL  → 採 type-specific：
#                    is_hematologic(MALIGNANCY_TYPE) → LATENCY_HEMATOLOGIC
#                    其餘                            → LATENCY_SOLID
#                    無事件者                         → LATENCY_HEMATOLOGIC（保守，不犧牲未發生事件者的 person-time）
build_study1_cohort <- function(base_df,
                                ct_summary_df,
                                malignancy_df,
                                prior_mal_ids,
                                study_start    = STUDY_START,
                                study_end      = STUDY_END,
                                latency_years  = LATENCY_YEARS) {
  
  cohort <- base_df |>
    anti_join(prior_mal_ids, by = "ID") |>
    left_join(ct_summary_df, by = "ID") |>
    mutate(
      EXPOSED       = !is.na(FIRST_CT_DATE),
      INDEX_DATE    = if_else(EXPOSED, FIRST_CT_DATE, study_start),
      AGE_AT_INDEX  = calc_age_years(BIRTH_DATE, INDEX_DATE),
      # I-5：兩組共同基準年齡（study_start 那天的年齡），可比較
      AGE_AT_2016   = calc_age_years(BIRTH_DATE, study_start),
      FIRST_CT_AGE  = if_else(EXPOSED, calc_age_years(BIRTH_DATE, FIRST_CT_DATE), NA_integer_),
      SEX           = if_else(ID_S == "1", "Male", "Female"),
      # G-3：分組欄位（暴露組才有 FIRST_CT_AGE_GROUP，對照組為 NA）
      FIRST_CT_AGE_GROUP = if_else(
        EXPOSED & !is.na(FIRST_CT_AGE),
        make_age_group(FIRST_CT_AGE, breaks = c(0, 5, 10, 15, 18)),
        NA_character_
      )
    ) |>
    filter(AGE_AT_INDEX <= MAX_AGE_AT_INDEX, AGE_AT_INDEX >= 0) |>
    left_join(malignancy_df, by = "ID")
  
  # H-2：依 latency_years 參數計算 ELIGIBLE_MALIGNANCY_DATE
  if (is.null(latency_years)) {
    # type-specific：依 MALIGNANCY_TYPE 套用 LATENCY_HEMATOLOGIC 或 LATENCY_SOLID
    cohort <- cohort |>
      mutate(
        APPLIED_LATENCY = if_else(
          !is.na(MALIGNANCY_TYPE) & is_hematologic(MALIGNANCY_TYPE),
          LATENCY_HEMATOLOGIC,
          LATENCY_SOLID
        )
      )
  } else {
    cohort <- cohort |>
      mutate(APPLIED_LATENCY = as.integer(latency_years))
  }
  
  cohort |>
    mutate(
      ELIGIBLE_MALIGNANCY_DATE = if_else(
        !is.na(FIRST_MALIGNANCY_DATE) &
          FIRST_MALIGNANCY_DATE > INDEX_DATE + years(APPLIED_LATENCY),
        FIRST_MALIGNANCY_DATE,
        NA_Date_
      ),
      OUTCOME = !is.na(ELIGIBLE_MALIGNANCY_DATE),
      # G-3：AGE_AT_DIAGNOSIS — 發生事件者的診斷時年齡
      AGE_AT_DIAGNOSIS = if_else(
        OUTCOME,
        calc_age_years(BIRTH_DATE, ELIGIBLE_MALIGNANCY_DATE),
        NA_integer_
      ),
      AGE_AT_DIAGNOSIS_GROUP = if_else(
        OUTCOME & !is.na(AGE_AT_DIAGNOSIS),
        make_age_group(AGE_AT_DIAGNOSIS,
                       breaks = c(0, 5, 10, 15, 20, Inf),
                       labels = c("0-5", "6-10", "11-15", "16-20", "21+")),
        NA_character_
      ),
      CENSOR_DATE = pmin(coalesce(ELIGIBLE_MALIGNANCY_DATE, study_end), study_end),
      TIME_MONTHS = as.numeric(
        interval(INDEX_DATE + years(APPLIED_LATENCY), CENSOR_DATE) / months(1)
      )
    ) |>
    filter(TIME_MONTHS >= 0) |>
    attach_income_urban()
}

# 主分析 cohort（latency = LATENCY_YEARS = 2L）
study1_cohort <- build_study1_cohort(
  base_df        = base_cohort,
  ct_summary_df  = ct_summary,
  malignancy_df  = malignancy_dx,
  prior_mal_ids  = prior_malignancy_ids,
  study_start    = STUDY_START,
  study_end      = STUDY_END,
  latency_years  = LATENCY_YEARS
)

log_progress(sprintf("study1_cohort（主分析，latency=%dy）：%s 人；EXPOSED=%s；OUTCOME=%s",
                     LATENCY_YEARS,
                     format(nrow(study1_cohort), big.mark = ","),
                     format(sum(study1_cohort$EXPOSED), big.mark = ","),
                     format(sum(study1_cohort$OUTCOME), big.mark = ",")))

## 4.2 Study 1 Analysis -------------------------------------------------------

### 4.2.1 Descriptive Statistics ----------------------------------------------
# G-4：用統一 make_descriptive() helper，三個 study 共用同一張 schema
# 同時保留舊版 group_by(EXPOSED) 風格的精簡表（study1-descriptive-summary.csv）
# 供向後相容

study1_descriptive_full <- make_descriptive(
  study1_cohort,
  age_var         = "AGE_AT_INDEX",
  age_at_2016_var = "AGE_AT_2016"
)
study1_descriptive_full |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-descriptive.csv"))

# 舊版 summary 表保留（向後相容）
study1_descriptive_summary <- study1_cohort |>
  group_by(EXPOSED) |>
  summarise(
    N                    = n(),
    N_MALE               = sum(SEX == "Male"),
    N_FEMALE             = sum(SEX == "Female"),
    MEDIAN_AGE_AT_2016   = median(AGE_AT_2016, na.rm = TRUE),
    MEDIAN_FIRST_CT_AGE  = median(FIRST_CT_AGE, na.rm = TRUE),
    # G-3：新增 AGE_AT_DIAGNOSIS 描述
    MEDIAN_AGE_AT_DIAGNOSIS = median(AGE_AT_DIAGNOSIS, na.rm = TRUE),
    N_MALIGNANCY         = sum(OUTCOME),
    INCIDENCE_RATE_K_PM  = sum(OUTCOME) / sum(TIME_MONTHS) * 1000,
    N_METRO              = sum(URBAN_STATUS == "Metro", na.rm = TRUE),
    N_NONMETRO           = sum(URBAN_STATUS == "Non-Metro", na.rm = TRUE),
    N_NA_URBAN           = sum(is.na(URBAN_STATUS)),
    .groups              = "drop"
  )
study1_descriptive_summary |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-descriptive-summary.csv"))

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
# 此處新增 Cox 比例風險模型作為副分析，將 person-time 納入考量。
# ⚠ 此 Cox 用「固定 INDEX_DATE」設計（暴露組從 FIRST_CT_DATE 起算、
#   對照組從 STUDY_START 起算），有 index date 不對稱偏誤。
#   §4.2.7 提供 time-varying exposure 平行版本（H-1，對齊文獻標準）。

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

### 4.2.6 Stratified Analyses ------------------------------------------------
# G-5：原有 sex / age_grp / ct_type / income / urban 分層保留，新增：
#   - FIRST_CT_SETTING（ER vs Clinic）— 僅暴露組
#   - AGE_AT_DIAGNOSIS_GROUP（呼應研究計畫所列 "stratified by age at diagnosis"）
#   - FIRST_CT_AGE_GROUP（呼應研究計畫所列 "stratified by age at CT exposure"）

study1_strat_sex    <- stratify_rr_full(study1_cohort, SEX,
                                        label = "SEX")
study1_strat_agegrp <- stratify_rr_full(
  study1_cohort |>
    mutate(AGE_GROUP = cut(AGE_AT_INDEX,
                           breaks = c(0, 5, 10, 15, 18),
                           include.lowest = TRUE)),
  AGE_GROUP, label = "AGE_GROUP"
)
study1_strat_cttype <- stratify_rr_full(
  study1_cohort |> filter(EXPOSED),
  FIRST_CT_TYPE, label = "FIRST_CT_TYPE (exposed only)"
)
study1_strat_income <- stratify_rr_full(study1_cohort, INCOME_Q,
                                        label = "INCOME_Q")
study1_strat_urban  <- stratify_rr_full(study1_cohort, URBAN_STATUS,
                                        label = "URBAN_STATUS")

# G-5 新增：CT setting (ER vs Clinic) — 僅暴露組
study1_strat_ctsetting <- stratify_rr_full(
  study1_cohort |> filter(EXPOSED),
  FIRST_CT_SETTING, label = "FIRST_CT_SETTING (exposed only)"
)

# G-5 新增：FIRST_CT_AGE_GROUP — 暴露組內部依 first CT 年齡比較
# 注意：這個分層的對照基線是「整個未暴露組」（不是同年齡的未暴露組）
#   每個 first_CT_age_group 對到「整個 unexposed cohort」算 RR
#   要做「同年齡組內 RR」需另外設計（暫不做）
study1_strat_firstctage <- stratify_rr_full(
  study1_cohort |>
    mutate(FIRST_CT_AGE_GROUP_OR_NONE = if_else(EXPOSED,
                                                FIRST_CT_AGE_GROUP,
                                                "Unexposed")),
  FIRST_CT_AGE_GROUP_OR_NONE,
  label = "FIRST_CT_AGE_GROUP"
)

# G-5 新增：AGE_AT_DIAGNOSIS_GROUP（事件發生者的診斷年齡分布）
# 注意：這個是「在發生事件者中，依診斷年齡分組的 EXPOSED vs UNEXPOSED 比例」
#       不是 RR；用 cross-tab 呈現
study1_dx_age_breakdown <- study1_cohort |>
  filter(OUTCOME) |>
  count(AGE_AT_DIAGNOSIS_GROUP, EXPOSED) |>
  pivot_wider(names_from = EXPOSED, values_from = n,
              names_prefix = "N_", values_fill = 0L) |>
  mutate(
    PCT_EXPOSED = N_TRUE / (N_TRUE + N_FALSE) * 100
  )
study1_dx_age_breakdown |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-strat-age-at-diagnosis.csv"))

list(
  sex          = study1_strat_sex,
  age_grp      = study1_strat_agegrp,
  ct_type      = study1_strat_cttype,
  income_q     = study1_strat_income,
  urban        = study1_strat_urban,
  ct_setting   = study1_strat_ctsetting,
  firstctage   = study1_strat_firstctage
) |>
  imap(~ write_csv(.x, file.path(
    OUTPUT_TABLE_PATH, paste0("study1-strat-", .y, ".csv")
  )))

### 4.2.7 Time-Varying Cox (H-1, Mathews 2013 模式) -------------------------
# 對齊 Mathews 2013 BMJ / Smoll 2023 AJNR / Pearce 2012 Lancet 標準：
# 每人從 study entry 開始貢獻 unexposed person-years，曾做 CT 者
# 在 first_CT + lag 後轉為 exposed。lag 期內仍算未暴露。
#
# 這個分析平行於 §4.2.5（保留固定 INDEX_DATE 版本以利對照）

log_progress("===== 4.2.7 Study 1 Time-Varying Cox（H-1） =====")

study1_tv <- build_tv_cohort(
  study1_cohort,
  study_start = STUDY_START,
  study_end   = STUDY_END,
  lag_days    = round(TV_LAG_YEARS * 365.25)
)

log_progress(sprintf("study1_tv：%s 列（包含 unexposed 與 exposed person-time 切段）",
                     format(nrow(study1_tv), big.mark = ",")))

# attach INCOME_Q / URBAN_STATUS — 用 ENTRY 那年（study_start 那年）
# 由於 build_tv_cohort 不帶 INCOME_Q / URBAN_STATUS，這裡再 join 一次
study1_tv <- study1_tv |>
  left_join(
    enrol_income_urban |>
      filter(YEAR == year(STUDY_START)) |>
      select(ID, INCOME_Q, URBAN_STATUS),
    by = "ID"
  )

study1_cox_tv <- coxph(
  Surv(tstart, tstop, OUTCOME_TV) ~ EXPOSED_TV + SEX + AGE_AT_2016 +
    INCOME_Q + URBAN_STATUS,
  data = study1_tv
)

tidy(study1_cox_tv, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(MODEL = "time_varying_lag2y") |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-cox-tv.csv"))

# PH assumption check
study1_cox_tv_zph <- cox.zph(study1_cox_tv)
tibble(
  TERM    = rownames(study1_cox_tv_zph$table),
  CHISQ   = study1_cox_tv_zph$table[, "chisq"],
  DF      = study1_cox_tv_zph$table[, "df"],
  P_VALUE = study1_cox_tv_zph$table[, "p"]
) |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-cox-tv-zph.csv"))

### 4.2.8 Dose-Response: N_CT_GROUP (H-1 副版, 對齊 Mathews) ---------------
# Mathews 2013：以「總 CT 次數」為 baseline covariate（不時變）
# 每多一次 IRR +0.16；我們改成 group 1 / 2-3 / >=4 看 trend

log_progress("===== 4.2.8 Study 1 N_CT_GROUP dose-response（H-1 副版） =====")

study1_dose <- study1_cohort |>
  mutate(
    N_CT_GROUP_F = factor(
      if_else(EXPOSED, N_CT_GROUP, "0"),
      levels = c("0", "1", "2-3", ">=4")
    )
  )

study1_cox_dose <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ N_CT_GROUP_F + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS,
  data = study1_dose
)
tidy(study1_cox_dose, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(MODEL = "n_ct_group_dose_response") |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-cox-dose-response.csv"))

# 同時印出 N_CT_GROUP 分布（暴露組內）
study1_n_ct_distrib <- study1_cohort |>
  filter(EXPOSED) |>
  count(N_CT_GROUP) |>
  mutate(PCT = n / sum(n) * 100)
study1_n_ct_distrib |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-n-ct-distribution.csv"))

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
# G-3：select 加入 AGE_AT_DIAGNOSIS / AGE_AT_DIAGNOSIS_GROUP / FIRST_CT_AGE_GROUP
#      / FIRST_CT_SETTING / N_CT_GROUP
study2_exposed <- study1_cohort |>
  filter(EXPOSED) |>
  select(ID, BIRTH_DATE, INDEX_DATE, FIRST_CT_DATE,
         AGE_AT_INDEX, AGE_AT_2016, FIRST_CT_AGE,
         FIRST_CT_AGE_GROUP, AGE_AT_DIAGNOSIS, AGE_AT_DIAGNOSIS_GROUP,
         SEX,
         OUTCOME, ELIGIBLE_MALIGNANCY_DATE, TIME_MONTHS,
         FIRST_CT_TYPE, FIRST_CT_SETTING, N_CT_TOTAL, N_CT_GROUP,
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
    # G-3：control 沒有 first CT，FIRST_CT_AGE_GROUP / SETTING / TYPE 都 NA
    FIRST_CT_AGE_GROUP = NA_character_,
    FIRST_CT_SETTING   = NA_character_,
    FIRST_CT_TYPE      = NA_character_,
    N_CT_GROUP         = NA_character_,
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
    # G-3：AGE_AT_DIAGNOSIS / GROUP
    AGE_AT_DIAGNOSIS = if_else(
      OUTCOME,
      calc_age_years(BIRTH_DATE, ELIGIBLE_MALIGNANCY_DATE),
      NA_integer_
    ),
    AGE_AT_DIAGNOSIS_GROUP = if_else(
      OUTCOME & !is.na(AGE_AT_DIAGNOSIS),
      make_age_group(AGE_AT_DIAGNOSIS,
                     breaks = c(0, 5, 10, 15, 20, Inf),
                     labels = c("0-5", "6-10", "11-15", "16-20", "21+")),
      NA_character_
    ),
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

### 5.3.0 Descriptive Statistics (G-4 補) -----------------------------------
# 03f 缺失 Study 2 的描述統計表，本版補上，schema 與 Study 1 對齊

study2_descriptive_full <- make_descriptive(
  study2_cohort,
  age_var         = "AGE_AT_INDEX",
  age_at_2016_var = "AGE_AT_2016"
)
study2_descriptive_full |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-descriptive.csv"))

# pair-level 統計（暴露/未暴露 sibling 各幾人、PAIR_SEX 分布等）
study2_pair_summary <- study2_cohort |>
  group_by(EXPOSED) |>
  summarise(
    N                       = n(),
    N_PAIRS                 = dplyr::n_distinct(PAIR_ID),
    N_MALE                  = sum(SEX == "Male"),
    N_FEMALE                = sum(SEX == "Female"),
    MEDIAN_AGE_AT_2016      = median(AGE_AT_2016, na.rm = TRUE),
    MEDIAN_FIRST_CT_AGE     = median(FIRST_CT_AGE, na.rm = TRUE),
    MEDIAN_AGE_AT_DIAGNOSIS = median(AGE_AT_DIAGNOSIS, na.rm = TRUE),
    N_MALIGNANCY            = sum(OUTCOME),
    INCIDENCE_RATE_K_PM     = sum(OUTCOME) / max(sum(TIME_MONTHS), 1) * 1000,
    .groups                 = "drop"
  )
study2_pair_summary |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-descriptive-summary.csv"))

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

# G-7：原版公式 RR = N_TRUE / N_FALSE 是「案例比」不是 RR
#       正確：RR = (N_EVENT_TRUE / N_TOTAL_TRUE) / (N_EVENT_FALSE / N_TOTAL_FALSE)
N2_TOTAL_TRUE  <- sum(study2_cohort$EXPOSED)
N2_TOTAL_FALSE <- sum(!study2_cohort$EXPOSED)

study2_rr_by_type <- study2_cohort |>
  filter(OUTCOME) |>
  count(MALIGNANCY_TYPE, EXPOSED) |>
  pivot_wider(names_from = EXPOSED, values_from = n,
              names_prefix = "N_EVENT_", values_fill = 0L) |>
  mutate(
    N_TOTAL_TRUE  = N2_TOTAL_TRUE,
    N_TOTAL_FALSE = N2_TOTAL_FALSE,
    # G-7：正確 RR 公式
    RR = (N_EVENT_TRUE  / N_TOTAL_TRUE) /
      (N_EVENT_FALSE / N_TOTAL_FALSE)
  )

study2_rr_by_type |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-rr-by-malignancy-type.csv"))

### 5.3.2 Cox Proportional Hazards (副分析, C-2) ----------------------------
# ⚠ 此版本「沒有」處理 within-pair correlation 也沒有 strata(PAIR_ID)，
#   屬於 "extended cohort" 設計：只是把 sibling 拉進當 control，沒做
#   within-family confounding 控制。SE 也可能略低估。
# §5.3.4 提供 strata(PAIR_ID) 平行版本（H-3，真正的 sibling-matched design）

study2_cox <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS,
  data = study2_cohort
)

tidy(study2_cox, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(MODEL = "extended_cohort_no_strata") |>
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

### 5.3.3 Stratified Analyses ------------------------------------------------
# G-5 新增：FIRST_CT_SETTING、AGE_AT_DIAGNOSIS_GROUP、FIRST_CT_AGE_GROUP
# 注意 Study 2 的 SEX 分層改用 SEX_CONCORDANCE（pair-level）

study2_strat_agegrp <- stratify_rr_full(
  study2_cohort |>
    mutate(AGE_GROUP = cut(AGE_AT_INDEX,
                           breaks = c(0, 5, 10, 15, 18),
                           include.lowest = TRUE)),
  AGE_GROUP, label = "AGE_GROUP"
)
study2_strat_sex_concordance <- stratify_rr_full(study2_cohort, SEX_CONCORDANCE,
                                                 label = "SEX_CONCORDANCE")
study2_strat_income <- stratify_rr_full(study2_cohort, INCOME_Q,
                                        label = "INCOME_Q")
study2_strat_urban  <- stratify_rr_full(study2_cohort, URBAN_STATUS,
                                        label = "URBAN_STATUS")

# G-5 新增：FIRST_CT_SETTING — 把 control 標 "Unexposed"，這樣可以
#   分別計算 ER vs Clinic vs Unexposed 的 RR
study2_strat_ctsetting <- stratify_rr_full(
  study2_cohort |>
    mutate(FIRST_CT_SETTING_OR_NONE = if_else(EXPOSED, FIRST_CT_SETTING, "Unexposed")),
  FIRST_CT_SETTING_OR_NONE, label = "FIRST_CT_SETTING"
)

# G-5 新增：FIRST_CT_AGE_GROUP（暴露組依首次 CT 年齡分組，control 標 "Unexposed"）
study2_strat_firstctage <- stratify_rr_full(
  study2_cohort |>
    mutate(FIRST_CT_AGE_GROUP_OR_NONE = if_else(EXPOSED, FIRST_CT_AGE_GROUP, "Unexposed")),
  FIRST_CT_AGE_GROUP_OR_NONE, label = "FIRST_CT_AGE_GROUP"
)

# G-5 新增：AGE_AT_DIAGNOSIS_GROUP — 在事件發生者中的暴露/未暴露分布
study2_dx_age_breakdown <- study2_cohort |>
  filter(OUTCOME) |>
  count(AGE_AT_DIAGNOSIS_GROUP, EXPOSED) |>
  pivot_wider(names_from = EXPOSED, values_from = n,
              names_prefix = "N_", values_fill = 0L) |>
  mutate(PCT_EXPOSED = N_TRUE / (N_TRUE + N_FALSE) * 100)
study2_dx_age_breakdown |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-strat-age-at-diagnosis.csv"))

list(
  age_grp         = study2_strat_agegrp,
  sex_concordance = study2_strat_sex_concordance,
  income_q        = study2_strat_income,
  urban           = study2_strat_urban,
  ct_setting      = study2_strat_ctsetting,
  firstctage      = study2_strat_firstctage
) |>
  imap(~ write_csv(.x, file.path(
    OUTPUT_TABLE_PATH, paste0("study2-strat-", .y, ".csv")
  )))

### 5.3.4 Cox with strata(PAIR_ID) — 真正的 sibling-matched design (H-3) ----
# 對齊 Lichtenstein NEJM 2000 / D'Onofrio sibling comparison 文獻
# 每個 pair 為一個 stratum，自動處理 within-family unobserved confounding
# ⚠ 注意：strata 後僅 within-pair 有 outcome 變異的 pair 貢獻資訊
#   N 會大幅縮水，HR 可能不穩定（這是 sibling design 固有代價，非 bug）

log_progress("===== 5.3.4 Study 2 Cox with strata(PAIR_ID)（H-3） =====")

# Sibling pair 內 outcome 變異統計（informative pair 數量）
sibling_pair_informativeness <- study2_cohort |>
  group_by(PAIR_ID) |>
  summarise(
    n_in_pair    = dplyr::n(),
    n_event      = sum(OUTCOME),
    n_exposed    = sum(EXPOSED),
    informative  = n_event > 0 & n_event < n_in_pair,
    .groups      = "drop"
  )
n_informative_pairs <- sum(sibling_pair_informativeness$informative)
log_progress(sprintf("Study 2 informative pairs（pair 內既有 event 也有 non-event）：%s / %s",
                     format(n_informative_pairs, big.mark = ","),
                     format(nrow(sibling_pair_informativeness), big.mark = ",")))

# 若 informative pair 太少，stratified Cox 可能 fail；try-catch 包裝
study2_cox_strata <- tryCatch(
  coxph(
    Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
      INCOME_Q + URBAN_STATUS + strata(PAIR_ID),
    data = study2_cohort
  ),
  error = function(e) {
    log_progress(sprintf("Study 2 strata(PAIR_ID) Cox 失敗：%s", e$message))
    NULL
  }
)

if (!is.null(study2_cox_strata)) {
  tidy(study2_cox_strata, conf.int = TRUE, exponentiate = TRUE) |>
    mutate(MODEL = "sibling_matched_strata_pair") |>
    write_csv(file.path(OUTPUT_TABLE_PATH, "study2-cox-strata.csv"))
} else {
  tibble(
    MODEL = "sibling_matched_strata_pair",
    NOTE  = "Failed to fit; likely too few informative pairs"
  ) |>
    write_csv(file.path(OUTPUT_TABLE_PATH, "study2-cox-strata.csv"))
}

# 同時跑一個輕量版：cluster(PAIR_ID) 修正 SE（不處理 family confounding）
study2_cox_cluster <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS + cluster(PAIR_ID),
  data = study2_cohort
)
tidy(study2_cox_cluster, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(MODEL = "extended_cohort_cluster_pair") |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-cox-cluster.csv"))

# 6. Study 3: Appendicitis Cohort ---------------------------------------------
# ⚠ Study 3 整段先暫停執行（user request 2026-05-10）：
#    NHIRD 目前無 IPDTO，住院期間 CT 醫令無法觀測，導致 §2.5 闌尾切除術
#    cohort 與 §6.x 主分析 cohort 都太小、無法給出穩定估計。
#    暫時用 if (FALSE) 包裝整段，方便日後資料補齊後一次啟用。
#
# 注意：§6.5 / §6.6 / §6.7 / §6.8 是 Study 1 的 sensitivity，不屬於 Study 3，
#       仍保留執行。
if (FALSE) {
  
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
    
    # 同時記錄首次腹部 CT 的日期、CT_TYPE、CT_SETTING（描述用）
    # G-1：取第一次窗內腹部 CT 的 setting（同日多筆 ER 優先）
    first_ab_ct <- ct_records |>
      filter(ORDER_CODE %in% abdomen_ct_proxy_codes) |>
      inner_join(base |> select(ID, APPENDIX_DATE), by = "ID") |>
      filter(
        FUNC_DATE >= APPENDIX_DATE - days(EXPOSURE_WINDOW_DAYS),
        FUNC_DATE <= APPENDIX_DATE + days(EXPOSURE_WINDOW_DAYS)
      ) |>
      mutate(SETTING_PRIO = if_else(CT_SETTING == "ER", 1L, 2L)) |>
      group_by(ID) |>
      arrange(FUNC_DATE, SETTING_PRIO, ORDER_CODE, .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      select(ID,
             FIRST_AB_CT_DATE    = FUNC_DATE,
             FIRST_CT_TYPE       = CT_TYPE,
             FIRST_CT_SETTING    = CT_SETTING)
    
    # 該 cohort 暴露者 N_CT 次數（窗內）
    n_ct_in_window <- ct_records |>
      filter(ORDER_CODE %in% abdomen_ct_proxy_codes) |>
      inner_join(base |> select(ID, APPENDIX_DATE), by = "ID") |>
      filter(
        FUNC_DATE >= APPENDIX_DATE - days(EXPOSURE_WINDOW_DAYS),
        FUNC_DATE <= APPENDIX_DATE + days(EXPOSURE_WINDOW_DAYS)
      ) |>
      count(ID, name = "N_CT_TOTAL")
    
    out <- base |>
      left_join(ct_flag, by = "ID") |>
      replace_na(list(EXPOSED = FALSE)) |>
      left_join(first_ab_ct, by = "ID") |>
      left_join(n_ct_in_window, by = "ID") |>
      mutate(
        FIRST_CT_AGE = if_else(EXPOSED,
                               calc_age_years(BIRTH_DATE, FIRST_AB_CT_DATE),
                               NA_integer_),
        # G-3：FIRST_CT_AGE_GROUP（暴露組才有意義）
        FIRST_CT_AGE_GROUP = if_else(
          EXPOSED & !is.na(FIRST_CT_AGE),
          make_age_group(FIRST_CT_AGE, breaks = c(0, 5, 10, 15, 18)),
          NA_character_
        ),
        # G-3：N_CT_GROUP
        N_CT_GROUP = if_else(EXPOSED, make_n_ct_group(N_CT_TOTAL), NA_character_)
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
        # G-3：AGE_AT_DIAGNOSIS / GROUP
        AGE_AT_DIAGNOSIS = if_else(
          OUTCOME,
          calc_age_years(BIRTH_DATE, ELIGIBLE_MALIGNANCY_DATE),
          NA_integer_
        ),
        AGE_AT_DIAGNOSIS_GROUP = if_else(
          OUTCOME & !is.na(AGE_AT_DIAGNOSIS),
          make_age_group(AGE_AT_DIAGNOSIS,
                         breaks = c(0, 5, 10, 15, 20, Inf),
                         labels = c("0-5", "6-10", "11-15", "16-20", "21+")),
          NA_character_
        ),
        CENSOR_DATE = pmin(coalesce(ELIGIBLE_MALIGNANCY_DATE, STUDY_END), STUDY_END),
        TIME_MONTHS = as.numeric(
          interval(INDEX_DATE + years(LATENCY_YEARS), CENSOR_DATE) / months(1)
        ),
        # I-5 補：Study 3 中 AGE_AT_2016 對描述統計也有意義（出生年距離 STUDY_START 的年齡）
        AGE_AT_2016 = calc_age_years(BIRTH_DATE, STUDY_START)
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
  # G-4：用統一 make_descriptive() helper
  
  study3_descriptive_full <- make_descriptive(
    study3_cohort,
    age_var         = "AGE_AT_APPENDIX",
    age_at_2016_var = "AGE_AT_2016"
  )
  study3_descriptive_full |>
    write_csv(file.path(OUTPUT_TABLE_PATH, "study3-descriptive.csv"))
  
  # 舊版精簡 summary 保留
  study3_descriptive_summary <- study3_cohort |>
    group_by(EXPOSED) |>
    summarise(
      N                       = n(),
      N_MALE                  = sum(SEX == "Male"),
      N_FEMALE                = sum(SEX == "Female"),
      MEDIAN_AGE_APPENDIX     = median(AGE_AT_APPENDIX, na.rm = TRUE),
      MEDIAN_FIRST_CT_AGE     = median(FIRST_CT_AGE, na.rm = TRUE),
      MEDIAN_AGE_AT_DIAGNOSIS = median(AGE_AT_DIAGNOSIS, na.rm = TRUE),
      N_MALIGNANCY            = sum(OUTCOME),
      INCIDENCE_RATE_K_PM     = sum(OUTCOME) / max(sum(TIME_MONTHS), 1) * 1000,
      N_METRO                 = sum(URBAN_STATUS == "Metro", na.rm = TRUE),
      N_NONMETRO              = sum(URBAN_STATUS == "Non-Metro", na.rm = TRUE),
      N_NA_URBAN              = sum(is.na(URBAN_STATUS)),
      .groups                 = "drop"
    )
  study3_descriptive_summary |>
    write_csv(file.path(OUTPUT_TABLE_PATH, "study3-descriptive-summary.csv"))
  
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
  
  # G-7：與 Study 2 同樣修正 RR by type 公式
  N3_TOTAL_TRUE  <- sum(study3_cohort$EXPOSED)
  N3_TOTAL_FALSE <- sum(!study3_cohort$EXPOSED)
  
  study3_rr_by_type <- study3_cohort |>
    filter(OUTCOME) |>
    count(MALIGNANCY_TYPE, EXPOSED) |>
    pivot_wider(names_from = EXPOSED, values_from = n,
                names_prefix = "N_EVENT_", values_fill = 0L) |>
    mutate(
      N_TOTAL_TRUE  = N3_TOTAL_TRUE,
      N_TOTAL_FALSE = N3_TOTAL_FALSE,
      RR = (N_EVENT_TRUE  / N_TOTAL_TRUE) /
        (N_EVENT_FALSE / N_TOTAL_FALSE)
    )
  
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
  
  ### 6.2.6 Stratified Analyses ------------------------------------------------
  # G-5 新增：FIRST_CT_SETTING、AGE_AT_DIAGNOSIS_GROUP、FIRST_CT_AGE_GROUP
  
  study3_strat_sex    <- stratify_rr_full(study3_cohort, SEX,
                                          label = "SEX")
  study3_strat_agegrp <- stratify_rr_full(
    study3_cohort |>
      mutate(AGE_GROUP = cut(AGE_AT_APPENDIX,
                             breaks = c(0, 5, 10, 15, 18),
                             include.lowest = TRUE)),
    AGE_GROUP, label = "AGE_GROUP"
  )
  study3_strat_income <- stratify_rr_full(study3_cohort, INCOME_Q,
                                          label = "INCOME_Q")
  study3_strat_urban  <- stratify_rr_full(study3_cohort, URBAN_STATUS,
                                          label = "URBAN_STATUS")
  
  # G-5 新增：CT setting (ER vs Clinic) — 暴露組才有意義；
  # 把 control 標 "Unexposed" 一起做 RR
  study3_strat_ctsetting <- stratify_rr_full(
    study3_cohort |>
      mutate(FIRST_CT_SETTING_OR_NONE = if_else(EXPOSED, FIRST_CT_SETTING,
                                                "Unexposed")),
    FIRST_CT_SETTING_OR_NONE, label = "FIRST_CT_SETTING"
  )
  
  # G-5 新增：FIRST_CT_AGE_GROUP（暴露組依首次 CT 年齡分組，control 標 "Unexposed"）
  study3_strat_firstctage <- stratify_rr_full(
    study3_cohort |>
      mutate(FIRST_CT_AGE_GROUP_OR_NONE = if_else(EXPOSED, FIRST_CT_AGE_GROUP,
                                                  "Unexposed")),
    FIRST_CT_AGE_GROUP_OR_NONE, label = "FIRST_CT_AGE_GROUP"
  )
  
  # G-5 新增：AGE_AT_DIAGNOSIS_GROUP — 在事件發生者中的暴露/未暴露分布
  study3_dx_age_breakdown <- study3_cohort |>
    filter(OUTCOME) |>
    count(AGE_AT_DIAGNOSIS_GROUP, EXPOSED) |>
    pivot_wider(names_from = EXPOSED, values_from = n,
                names_prefix = "N_", values_fill = 0L) |>
    mutate(PCT_EXPOSED = N_TRUE / (N_TRUE + N_FALSE) * 100)
  study3_dx_age_breakdown |>
    write_csv(file.path(OUTPUT_TABLE_PATH, "study3-strat-age-at-diagnosis.csv"))
  
  list(
    sex          = study3_strat_sex,
    age_grp      = study3_strat_agegrp,
    income_q     = study3_strat_income,
    urban        = study3_strat_urban,
    ct_setting   = study3_strat_ctsetting,
    firstctage   = study3_strat_firstctage
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
  
}  # end of if (FALSE) — Study 3 暫停執行區塊

# ════════════════════════════════════════════════════════════════════════════
# 6.5–6.7 SENSITIVITY ANALYSES (H-2 / H-4)
# ════════════════════════════════════════════════════════════════════════════
# 主分析（§4-§6）使用 LATENCY_YEARS = 2L、STUDY_START = 2016-01-01。
# 這個 section 跑：
#   §6.5  H-2a: solid 5y / hematologic 2y type-specific latency
#   §6.6  H-2b: 統一 5y latency（所有 cancer 都 5y）
#   §6.7  H-4 : Extended window 2000-2023（與 H-2a 結合：type-specific latency）
#
# 所有結果寫入 OUTPUT_SENS_PATH (= outputs/tables/sensitivity/) 子目錄。
#
# ⚠ 5y latency 在 2016-2023 (8 年) window 下會大幅截短 person-time：
#   2018+ 做 CT 的人實際可觀察 < 3 年。Power 對 solid tumor 訊號會偏弱。
#   H-4 結合 extended window 是為了補 follow-up 時間，看是否能恢復 power。
# ════════════════════════════════════════════════════════════════════════════

log_progress("===== 6.5 Sensitivity: type-specific latency (solid 5y / heme 2y) =====")

## 6.5 H-2a: type-specific latency ------------------------------------------

study1_cohort_typespec <- build_study1_cohort(
  base_df        = base_cohort,
  ct_summary_df  = ct_summary,
  malignancy_df  = malignancy_dx,
  prior_mal_ids  = prior_malignancy_ids,
  study_start    = STUDY_START,
  study_end      = STUDY_END,
  latency_years  = NULL  # → type-specific
)

log_progress(sprintf("study1_cohort_typespec：%s 人；EXPOSED=%s；OUTCOME=%s",
                     format(nrow(study1_cohort_typespec), big.mark = ","),
                     format(sum(study1_cohort_typespec$EXPOSED), big.mark = ","),
                     format(sum(study1_cohort_typespec$OUTCOME), big.mark = ",")))

# Cox（與主分析同設計，只差 latency）
study1_cox_typespec <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS,
  data = study1_cohort_typespec
)
tidy(study1_cox_typespec, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(MODEL = "study1_cox_solid5y_heme2y") |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-cox-typespec-latency.csv"))

# RR by type — heme 用 2y、solid 用 5y 各自比較
N1TS_TOTAL_TRUE  <- sum(study1_cohort_typespec$EXPOSED)
N1TS_TOTAL_FALSE <- sum(!study1_cohort_typespec$EXPOSED)

study1_rr_by_type_typespec <- study1_cohort_typespec |>
  filter(OUTCOME) |>
  count(MALIGNANCY_TYPE, EXPOSED) |>
  pivot_wider(names_from = EXPOSED, values_from = n,
              names_prefix = "N_EVENT_", values_fill = 0L) |>
  mutate(
    N_TOTAL_TRUE  = N1TS_TOTAL_TRUE,
    N_TOTAL_FALSE = N1TS_TOTAL_FALSE,
    APPLIED_LATENCY_YEARS = if_else(is_hematologic(MALIGNANCY_TYPE),
                                    LATENCY_HEMATOLOGIC, LATENCY_SOLID),
    RR = (N_EVENT_TRUE  / N_TOTAL_TRUE) /
      (N_EVENT_FALSE / N_TOTAL_FALSE)
  )
study1_rr_by_type_typespec |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-rr-by-type-typespec-latency.csv"))

## 6.6 H-2b: 統一 5y latency ------------------------------------------------

log_progress("===== 6.6 Sensitivity: 統一 5y latency =====")

study1_cohort_5y <- build_study1_cohort(
  base_df        = base_cohort,
  ct_summary_df  = ct_summary,
  malignancy_df  = malignancy_dx,
  prior_mal_ids  = prior_malignancy_ids,
  study_start    = STUDY_START,
  study_end      = STUDY_END,
  latency_years  = 5L
)

log_progress(sprintf("study1_cohort_5y：%s 人；EXPOSED=%s；OUTCOME=%s",
                     format(nrow(study1_cohort_5y), big.mark = ","),
                     format(sum(study1_cohort_5y$EXPOSED), big.mark = ","),
                     format(sum(study1_cohort_5y$OUTCOME), big.mark = ",")))

study1_cox_5y <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS,
  data = study1_cohort_5y
)
tidy(study1_cox_5y, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(MODEL = "study1_cox_uniform_5y") |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-cox-5y-latency.csv"))

study1_rr_overall_5y <- study1_cohort_5y |>
  group_by(EXPOSED) |>
  summarise(
    N_EVENT = sum(OUTCOME),
    N_TOTAL = n(),
    RISK    = N_EVENT / N_TOTAL,
    .groups = "drop"
  ) |>
  summarise(
    RR = RISK[EXPOSED] / RISK[!EXPOSED],
    LATENCY_YEARS = 5L
  )
study1_rr_overall_5y |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-rr-overall-5y-latency.csv"))

## 6.7 H-4: Extended window 2000-2023 + type-specific latency ---------------
# ⚠ 此 section 需要 base_cohort_extended（含 1982-2023 出生者）
#    與 prior_malignancy_ids_extended（個人化 entry date 排除 prior cancer）

log_progress("===== 6.7 Sensitivity: Extended window 2000-2023 (H-4) =====")

# H-4 / H-6：擴展 base cohort 到 1982-2023 出生者
base_cohort_extended <- pers_info |>
  anti_join(hereditary_exclusion, by = "ID") |>
  mutate(
    BIRTH_DATE = resolve_birth_date(ID_BIRTHYM),
    AGE_2000   = calc_age_years(BIRTH_DATE, EXTENDED_STUDY_START),
    AGE_2023   = calc_age_years(BIRTH_DATE, STUDY_END),
    # H-6：個人化 entry date — 取 max(birth_date, EXTENDED_STUDY_START)
    ENTRY_DATE = pmax(BIRTH_DATE, EXTENDED_STUDY_START)
  ) |>
  filter(
    !is.na(BIRTH_DATE),
    BIRTH_DATE <= STUDY_END,
    BIRTH_DATE + years(MAX_AGE_AT_INDEX + 1L) > EXTENDED_STUDY_START
  )

log_progress(sprintf("base_cohort_extended：%s 人",
                     format(nrow(base_cohort_extended), big.mark = ",")))

# H-8：用 effective start date 排除 prior malignancy
#   暴露者：effective_start = max(BIRTH_DATE, EXTENDED_STUDY_START, FIRST_CT_DATE)
#   未暴露者：effective_start = max(BIRTH_DATE, EXTENDED_STUDY_START)
# 注意：1982-1990 出生者進入研究時可能 10-18 歲，1990s 期間若有 malignancy
#      （NHIRD 之前看不到）會被誤分類為「無 prior」 — 這是已知 limitation
#
# ⚠ 此處需要 ct_summary_extended，所以放在 ct_summary_extended 之後計算
# 為避免 forward reference，先建 ct_summary_extended，再算 prior

# H-4：用 EXTENDED_STUDY_START 重建 ct_summary（涵蓋 2000-2023 全部 CT）
ct_summary_extended <- ct_records |>
  filter(FUNC_DATE >= EXTENDED_STUDY_START, FUNC_DATE <= STUDY_END) |>
  group_by(ID) |>
  summarise(
    FIRST_CT_DATE = min(FUNC_DATE),
    N_CT_TOTAL    = dplyr::n(),
    CT_TYPES      = paste(sort(unique(CT_TYPE)), collapse = "; "),
    .groups = "drop"
  ) |>
  left_join(
    ct_records |>
      filter(FUNC_DATE >= EXTENDED_STUDY_START, FUNC_DATE <= STUDY_END) |>
      group_by(ID) |>
      arrange(FUNC_DATE, ORDER_CODE, .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      select(ID, FIRST_CT_TYPE = CT_TYPE),
    by = "ID"
  ) |>
  left_join(
    ct_records |>
      filter(FUNC_DATE >= EXTENDED_STUDY_START, FUNC_DATE <= STUDY_END) |>
      mutate(SETTING_PRIO = if_else(CT_SETTING == "ER", 1L, 2L)) |>
      group_by(ID) |>
      arrange(FUNC_DATE, SETTING_PRIO, ORDER_CODE, .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      select(ID, FIRST_CT_SETTING = CT_SETTING),
    by = "ID"
  ) |>
  mutate(N_CT_GROUP = make_n_ct_group(N_CT_TOTAL))

# H-8：H-4 extended cohort 的 prior malignancy 排除
prior_malignancy_ids_extended <- base_cohort_extended |>
  left_join(ct_summary_extended |> select(ID, FIRST_CT_DATE), by = "ID") |>
  mutate(
    EFFECTIVE_START = pmax(
      BIRTH_DATE,
      EXTENDED_STUDY_START,
      coalesce(FIRST_CT_DATE, EXTENDED_STUDY_START),
      na.rm = TRUE
    )
  ) |>
  inner_join(malignancy_dx |> select(ID, FIRST_MALIGNANCY_DATE),
             by = "ID") |>
  filter(FIRST_MALIGNANCY_DATE < EFFECTIVE_START) |>
  distinct(ID)

log_progress(sprintf("prior_malignancy_ids_extended（H-8 effective start date）：%s 人",
                     format(nrow(prior_malignancy_ids_extended), big.mark = ",")))

# Build extended cohort（用 type-specific latency）
study1_cohort_extended <- build_study1_cohort(
  base_df        = base_cohort_extended,
  ct_summary_df  = ct_summary_extended,
  malignancy_df  = malignancy_dx,
  prior_mal_ids  = prior_malignancy_ids_extended,
  study_start    = EXTENDED_STUDY_START,
  study_end      = STUDY_END,
  latency_years  = NULL  # type-specific
)

# H-10：H-4 extended cohort 的 AGE_AT_DIAGNOSIS_GROUP 重新分組（擴大到 41+）
# build_study1_cohort 內部用主分析的分組 c(0, 5, 10, 15, 20, Inf)，
# 對 H-4 中 1982 出生者最大 41 歲不夠用。這裡覆寫一次。
study1_cohort_extended <- study1_cohort_extended |>
  mutate(
    AGE_AT_DIAGNOSIS_GROUP = if_else(
      OUTCOME & !is.na(AGE_AT_DIAGNOSIS),
      make_age_group(AGE_AT_DIAGNOSIS,
                     breaks = c(0, 5, 10, 15, 20, 30, 40, Inf),
                     labels = c("0-5", "6-10", "11-15", "16-20",
                                "21-30", "31-40", "41+")),
      NA_character_
    )
  )

log_progress(sprintf("study1_cohort_extended：%s 人；EXPOSED=%s；OUTCOME=%s",
                     format(nrow(study1_cohort_extended), big.mark = ","),
                     format(sum(study1_cohort_extended$EXPOSED), big.mark = ","),
                     format(sum(study1_cohort_extended$OUTCOME), big.mark = ",")))

# H-10 診斷：印出 H-4 cohort 的 AGE_AT_DIAGNOSIS_GROUP 分布
log_progress("===== §6.7 H-4 cohort AGE_AT_DIAGNOSIS_GROUP 分布 =====")
print(study1_cohort_extended |>
        filter(OUTCOME) |>
        count(AGE_AT_DIAGNOSIS_GROUP, EXPOSED) |>
        pivot_wider(names_from = EXPOSED, values_from = n,
                    names_prefix = "N_", values_fill = 0L))

# H-10 診斷：印出 H-4 cohort 的 AGE_AT_DIAGNOSIS_GROUP 分布
log_progress("===== §6.7 H-4 cohort AGE_AT_DIAGNOSIS_GROUP 分布 =====")
print(study1_cohort_extended |>
        filter(OUTCOME) |>
        count(AGE_AT_DIAGNOSIS_GROUP, EXPOSED) |>
        pivot_wider(names_from = EXPOSED, values_from = n,
                    names_prefix = "N_", values_fill = 0L))

# G-4：H-4 cohort 的描述統計表（Table 1 等價）
study1_extended_descriptive <- make_descriptive(
  study1_cohort_extended,
  age_var         = "AGE_AT_INDEX",
  age_at_2016_var = "AGE_AT_2016"
)
study1_extended_descriptive |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-extended-descriptive.csv"))

study1_cox_extended <- coxph(
  Surv(TIME_MONTHS, OUTCOME) ~ EXPOSED + SEX + AGE_AT_INDEX +
    INCOME_Q + URBAN_STATUS,
  data = study1_cohort_extended
)
tidy(study1_cox_extended, conf.int = TRUE, exponentiate = TRUE) |>
  mutate(MODEL = "study1_cox_extended_2000_typespec") |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-cox-extended-2000.csv"))

# RR by type
N1ext_T <- sum(study1_cohort_extended$EXPOSED)
N1ext_F <- sum(!study1_cohort_extended$EXPOSED)
study1_rr_by_type_extended <- study1_cohort_extended |>
  filter(OUTCOME) |>
  count(MALIGNANCY_TYPE, EXPOSED) |>
  pivot_wider(names_from = EXPOSED, values_from = n,
              names_prefix = "N_EVENT_", values_fill = 0L) |>
  mutate(
    N_TOTAL_TRUE  = N1ext_T,
    N_TOTAL_FALSE = N1ext_F,
    APPLIED_LATENCY_YEARS = if_else(is_hematologic(MALIGNANCY_TYPE),
                                    LATENCY_HEMATOLOGIC, LATENCY_SOLID),
    RR = (N_EVENT_TRUE  / N_TOTAL_TRUE) /
      (N_EVENT_FALSE / N_TOTAL_FALSE)
  )
study1_rr_by_type_extended |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-rr-by-type-extended-2000.csv"))

# Time-varying Cox 也跑一個 extended 版本（H-1 + H-4 結合）
log_progress("===== 6.7b Time-varying Cox on extended cohort =====")
study1_tv_extended <- build_tv_cohort(
  study1_cohort_extended,
  study_start = EXTENDED_STUDY_START,
  study_end   = STUDY_END,
  lag_days    = round(TV_LAG_YEARS * 365.25)
) |>
  left_join(
    enrol_income_urban |>
      filter(YEAR == year(EXTENDED_STUDY_START)) |>
      select(ID, INCOME_Q, URBAN_STATUS),
    by = "ID"
  )

study1_cox_tv_extended <- tryCatch(
  coxph(
    Surv(tstart, tstop, OUTCOME_TV) ~ EXPOSED_TV + SEX + AGE_AT_2016 +
      INCOME_Q + URBAN_STATUS,
    data = study1_tv_extended
  ),
  error = function(e) {
    log_progress(sprintf("Time-varying Cox on extended cohort failed: %s",
                         e$message))
    NULL
  }
)

if (!is.null(study1_cox_tv_extended)) {
  tidy(study1_cox_tv_extended, conf.int = TRUE, exponentiate = TRUE) |>
    mutate(MODEL = "study1_cox_tv_extended_2000") |>
    write_csv(file.path(OUTPUT_SENS_PATH, "study1-cox-tv-extended-2000.csv"))
}

## 6.8 Sensitivity Summary Table -------------------------------------------
# 把所有 sensitivity 的 EXPOSED HR 拉出來放成一張表，方便比較

extract_exposed_hr <- function(cox_model, label) {
  if (is.null(cox_model)) {
    return(tibble(
      MODEL = label, HR = NA_real_, LCI = NA_real_, UCI = NA_real_,
      P_VALUE = NA_real_
    ))
  }
  td <- tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)
  exp_row <- td |>
    filter(term %in% c("EXPOSEDTRUE", "EXPOSED_TV", "EXPOSEDTRUE")) |>
    slice(1)
  if (nrow(exp_row) == 0) {
    return(tibble(
      MODEL = label, HR = NA_real_, LCI = NA_real_, UCI = NA_real_,
      P_VALUE = NA_real_
    ))
  }
  tibble(
    MODEL   = label,
    HR      = exp_row$estimate,
    LCI     = exp_row$conf.low,
    UCI     = exp_row$conf.high,
    P_VALUE = exp_row$p.value
  )
}

sensitivity_summary <- bind_rows(
  extract_exposed_hr(study1_cox,           "main_2y_fixed_index"),
  extract_exposed_hr(study1_cox_tv,        "main_2y_time_varying"),
  extract_exposed_hr(study1_cox_typespec,  "sens_typespec_solid5y_heme2y"),
  extract_exposed_hr(study1_cox_5y,        "sens_uniform_5y"),
  extract_exposed_hr(study1_cox_extended,  "sens_extended_2000_typespec"),
  extract_exposed_hr(study1_cox_tv_extended, "sens_extended_2000_time_varying")
)
sensitivity_summary |>
  write_csv(file.path(OUTPUT_SENS_PATH, "study1-sensitivity-summary.csv"))

cat("\n========== Sensitivity Summary (Study 1 EXPOSED HR) ==========\n")
print(sensitivity_summary)

# 7. Summary Tables and Figures -----------------------------------------------

## 7.1 Combined RR Summary Table -----------------------------------------------

rr_summary <- bind_rows(
  study1_rr_overall |> mutate(Study = "Study 1 (Population Cohort)"),
  study2_rr_overall |> mutate(Study = "Study 2 (Sibling-Matched)")
  # Study 3 暫停執行，rr_summary 暫不含
)

rr_summary |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "rr-summary-all-studies.csv"))

## 7.1b Publication-Quality Table 1: Baseline Characteristics -----------------
# 對齊 NEJM / JAMA / BMJ 期刊規範的 Table 1 格式
# 欄位：Characteristic | Overall | Exposed (CT) | Unexposed | SMD | P-value
#
# 統計方法：
#   - 連續變項：median (IQR) → Wilcoxon rank-sum；mean (SD) → t-test
#     （兒童 CT cohort 多用 median，因 age 與 follow-up time 偏態）
#   - 類別變項：n (%) → χ² test（cell N < 5 用 Fisher exact）
#   - SMD (Standardized Mean Difference)：
#     連續 = (mean_e - mean_u) / sqrt((SD_e^2 + SD_u^2) / 2)
#     類別 = sqrt(sum((p_e - p_u)^2 / (p_avg(1-p_avg))))
#     ⚠ SMD > 0.1 通常視為「不平衡」（Austin 2009, J Clin Epidemiol）
#
# 參考文獻 Table 1 結構：
#   - Mathews 2013 BMJ Table 1（澳洲 CT cohort）
#   - Krille 2015 (EPI-CT 設計論文) Table 1
#   - Korea 2025 NHIS Table 1

# Helper: 連續變項 row（median + IQR）
table1_row_continuous <- function(cohort, var_name, label, transform_fun = identity) {
  df <- cohort |>
    mutate(.v = transform_fun(.data[[var_name]])) |>
    filter(!is.na(.v))
  if (nrow(df) == 0) return(NULL)
  
  overall <- df$.v
  exp_v   <- df$.v[df$EXPOSED]
  unx_v   <- df$.v[!df$EXPOSED]
  
  # SMD for continuous
  mean_e <- mean(exp_v, na.rm = TRUE)
  mean_u <- mean(unx_v, na.rm = TRUE)
  sd_e   <- sd(exp_v,   na.rm = TRUE)
  sd_u   <- sd(unx_v,   na.rm = TRUE)
  smd    <- if (is.na(sd_e) || is.na(sd_u) || (sd_e^2 + sd_u^2) == 0) {
    NA_real_
  } else {
    (mean_e - mean_u) / sqrt((sd_e^2 + sd_u^2) / 2)
  }
  
  # P-value: Wilcoxon rank-sum
  pval <- tryCatch(
    {
      if (length(exp_v) > 0 && length(unx_v) > 0) {
        suppressWarnings(wilcox.test(exp_v, unx_v)$p.value)
      } else NA_real_
    },
    error = function(e) NA_real_
  )
  
  fmt <- function(x) sprintf("%.1f [%.1f, %.1f]",
                             stats::median(x),
                             stats::quantile(x, 0.25),
                             stats::quantile(x, 0.75))
  
  tibble(
    Characteristic = label,
    Category       = "Median [IQR]",
    Overall        = fmt(overall),
    Exposed        = fmt(exp_v),
    Unexposed      = fmt(unx_v),
    SMD            = if (is.na(smd)) "—" else sprintf("%.3f", abs(smd)),
    P_value        = if (is.na(pval)) "—" else
      if (pval < 0.001) "<0.001" else sprintf("%.3f", pval)
  )
}

# Helper: 類別變項 row（n (%) 每個 level 一列；第一列附 SMD/P-value）
table1_rows_categorical <- function(cohort, var_name, label) {
  if (!var_name %in% names(cohort)) return(NULL)
  
  df <- cohort |>
    mutate(.v = as.character(.data[[var_name]])) |>
    mutate(.v = if_else(is.na(.v) | .v == "", "Missing", .v))
  
  total_n <- nrow(df)
  exp_n   <- sum(df$EXPOSED, na.rm = TRUE)
  unx_n   <- sum(!df$EXPOSED, na.rm = TRUE)
  if (total_n == 0) return(NULL)
  
  levs <- sort(unique(df$.v))
  
  # SMD for categorical
  smd_cat <- tryCatch({
    p_e <- table(df$.v[df$EXPOSED])  / max(exp_n, 1)
    p_u <- table(df$.v[!df$EXPOSED]) / max(unx_n, 1)
    all_levs <- union(names(p_e), names(p_u))
    p_e_v <- as.numeric(p_e[all_levs]); p_e_v[is.na(p_e_v)] <- 0
    p_u_v <- as.numeric(p_u[all_levs]); p_u_v[is.na(p_u_v)] <- 0
    if (length(all_levs) == 1) {
      # binary
      p_e1 <- p_e_v[1]; p_u1 <- p_u_v[1]
      p_avg <- (p_e1 + p_u1) / 2
      if (p_avg %in% c(0, 1)) NA_real_ else
        (p_e1 - p_u1) / sqrt(p_avg * (1 - p_avg))
    } else {
      # multinomial SMD（Yang & Dalton 2012）
      diff <- p_e_v - p_u_v
      avg  <- (p_e_v + p_u_v) / 2
      # 用對角覆變矩陣近似
      sqrt(sum(diff^2 / pmax(avg * (1 - avg), 1e-10)))
    }
  }, error = function(e) NA_real_)
  
  # P-value: χ² test on 2 × K table
  tab <- df |>
    count(EXPOSED, .v) |>
    pivot_wider(names_from = .v, values_from = n, values_fill = 0L) |>
    select(-EXPOSED) |>
    as.matrix()
  pval <- tryCatch(
    suppressWarnings(chisq.test(tab)$p.value),
    error = function(e) NA_real_
  )
  
  # 各 level 一列
  rows <- lapply(seq_along(levs), function(i) {
    lv  <- levs[i]
    n_t <- sum(df$.v == lv)
    n_e <- sum(df$.v == lv & df$EXPOSED)
    n_u <- sum(df$.v == lv & !df$EXPOSED)
    tibble(
      Characteristic = if (i == 1) label else "",
      Category       = lv,
      Overall        = sprintf("%s (%.1f%%)", format(n_t, big.mark = ","),
                               100 * n_t / max(total_n, 1)),
      Exposed        = sprintf("%s (%.1f%%)", format(n_e, big.mark = ","),
                               100 * n_e / max(exp_n, 1)),
      Unexposed      = sprintf("%s (%.1f%%)", format(n_u, big.mark = ","),
                               100 * n_u / max(unx_n, 1)),
      SMD            = if (i == 1) {
        if (is.na(smd_cat)) "—" else sprintf("%.3f", abs(smd_cat))
      } else "",
      P_value        = if (i == 1) {
        if (is.na(pval)) "—" else
          if (pval < 0.001) "<0.001" else sprintf("%.3f", pval)
      } else ""
    )
  })
  bind_rows(rows)
}

# 共用：建一個 study 的 Table 1
build_table1 <- function(cohort, study_label, age_var = "AGE_AT_INDEX") {
  total_n <- nrow(cohort)
  exp_n   <- sum(cohort$EXPOSED, na.rm = TRUE)
  unx_n   <- sum(!cohort$EXPOSED, na.rm = TRUE)
  
  # 第一列：N（總人數）
  header <- tibble(
    Characteristic = "Total participants",
    Category       = "",
    Overall        = format(total_n, big.mark = ","),
    Exposed        = format(exp_n,   big.mark = ","),
    Unexposed      = format(unx_n,   big.mark = ","),
    SMD            = "",
    P_value        = ""
  )
  
  # Person-years
  py <- if ("TIME_MONTHS" %in% names(cohort)) {
    py_t <- sum(cohort$TIME_MONTHS, na.rm = TRUE) / 12
    py_e <- sum(cohort$TIME_MONTHS[cohort$EXPOSED],  na.rm = TRUE) / 12
    py_u <- sum(cohort$TIME_MONTHS[!cohort$EXPOSED], na.rm = TRUE) / 12
    tibble(
      Characteristic = "Person-years of follow-up",
      Category       = "",
      Overall        = format(round(py_t), big.mark = ","),
      Exposed        = format(round(py_e), big.mark = ","),
      Unexposed      = format(round(py_u), big.mark = ","),
      SMD            = "",
      P_value        = ""
    )
  } else NULL
  
  # Median follow-up
  fu_median <- if ("TIME_MONTHS" %in% names(cohort)) {
    table1_row_continuous(cohort, "TIME_MONTHS",
                          "Follow-up time (months)")
  } else NULL
  
  # Sex
  sex_block <- table1_rows_categorical(cohort, "SEX", "Sex")
  
  # Age at index
  age_block <- table1_row_continuous(
    cohort, age_var, sprintf("Age at index (years, %s)", age_var)
  )
  
  # Age at 2016（兩組共同基準）
  age2016_block <- if ("AGE_AT_2016" %in% names(cohort)) {
    table1_row_continuous(cohort, "AGE_AT_2016",
                          "Age at study start (2016-01-01) (years)")
  } else NULL
  
  # First CT age（僅暴露組有意義 — SMD/P-value 不計算）
  first_ct_age_block <- if ("FIRST_CT_AGE" %in% names(cohort) &&
                            sum(!is.na(cohort$FIRST_CT_AGE)) > 0) {
    exp_first_ct <- cohort$FIRST_CT_AGE[cohort$EXPOSED]
    tibble(
      Characteristic = "Age at first CT (years, exposed only)",
      Category       = "Median [IQR]",
      Overall        = sprintf("%.1f [%.1f, %.1f]",
                               stats::median(exp_first_ct, na.rm = TRUE),
                               stats::quantile(exp_first_ct, 0.25, na.rm = TRUE),
                               stats::quantile(exp_first_ct, 0.75, na.rm = TRUE)),
      Exposed        = sprintf("%.1f [%.1f, %.1f]",
                               stats::median(exp_first_ct, na.rm = TRUE),
                               stats::quantile(exp_first_ct, 0.25, na.rm = TRUE),
                               stats::quantile(exp_first_ct, 0.75, na.rm = TRUE)),
      Unexposed      = "—",
      SMD            = "—",
      P_value        = "—"
    )
  } else NULL
  
  # First CT age group（類別）
  first_ct_age_grp_block <- if ("FIRST_CT_AGE_GROUP" %in% names(cohort)) {
    sub <- cohort |>
      mutate(FIRST_CT_AGE_GROUP = if_else(EXPOSED, FIRST_CT_AGE_GROUP, NA_character_))
    table1_rows_categorical(sub |> filter(EXPOSED),
                            "FIRST_CT_AGE_GROUP",
                            "Age at first CT group (exposed only)")
  } else NULL
  
  # Age at diagnosis（僅事件者有意義）
  age_dx_block <- if ("AGE_AT_DIAGNOSIS" %in% names(cohort) &&
                      sum(!is.na(cohort$AGE_AT_DIAGNOSIS)) > 0) {
    table1_row_continuous(cohort |> filter(OUTCOME),
                          "AGE_AT_DIAGNOSIS",
                          "Age at malignancy diagnosis (years, events only)")
  } else NULL
  
  age_dx_grp_block <- if ("AGE_AT_DIAGNOSIS_GROUP" %in% names(cohort)) {
    table1_rows_categorical(
      cohort |> filter(OUTCOME),
      "AGE_AT_DIAGNOSIS_GROUP",
      "Age at malignancy diagnosis group (events only)"
    )
  } else NULL
  
  # Income quartile
  income_block <- table1_rows_categorical(cohort, "INCOME_Q",
                                          "Income quartile")
  
  # Urban status
  urban_block <- table1_rows_categorical(cohort, "URBAN_STATUS",
                                         "Urban status")
  
  # CT setting (ER vs Clinic)
  ct_setting_block <- if ("FIRST_CT_SETTING" %in% names(cohort)) {
    sub <- cohort |> filter(EXPOSED)
    rows <- table1_rows_categorical(sub, "FIRST_CT_SETTING",
                                    "CT setting at first scan (exposed only)")
    # 替換 SMD / P-value（因為只有暴露組）
    if (!is.null(rows) && nrow(rows) > 0L) {
      rows$SMD     <- ""
      rows$P_value <- ""
    }
    rows
  } else NULL
  
  # First CT type
  ct_type_block <- if ("FIRST_CT_TYPE" %in% names(cohort)) {
    sub <- cohort |> filter(EXPOSED)
    rows <- table1_rows_categorical(sub, "FIRST_CT_TYPE",
                                    "First CT type (exposed only)")
    if (!is.null(rows) && nrow(rows) > 0L) {
      rows$SMD     <- ""
      rows$P_value <- ""
    }
    rows
  } else NULL
  
  # Number of CT scans group
  n_ct_block <- if ("N_CT_GROUP" %in% names(cohort)) {
    sub <- cohort |> filter(EXPOSED)
    rows <- table1_rows_categorical(sub, "N_CT_GROUP",
                                    "Number of CT scans (exposed only)")
    if (!is.null(rows) && nrow(rows) > 0L) {
      rows$SMD     <- ""
      rows$P_value <- ""
    }
    rows
  } else NULL
  
  # Outcome: malignancy
  n_event_e <- sum(cohort$OUTCOME &  cohort$EXPOSED, na.rm = TRUE)
  n_event_u <- sum(cohort$OUTCOME & !cohort$EXPOSED, na.rm = TRUE)
  outcome_row <- tibble(
    Characteristic = "Malignancy events",
    Category       = "",
    Overall        = sprintf("%s (%.2f%%)",
                             format(n_event_e + n_event_u, big.mark = ","),
                             100 * (n_event_e + n_event_u) / max(total_n, 1)),
    Exposed        = sprintf("%s (%.2f%%)",
                             format(n_event_e, big.mark = ","),
                             100 * n_event_e / max(exp_n, 1)),
    Unexposed      = sprintf("%s (%.2f%%)",
                             format(n_event_u, big.mark = ","),
                             100 * n_event_u / max(unx_n, 1)),
    SMD            = "",
    P_value        = ""
  )
  
  bind_rows(
    header,
    py,
    fu_median,
    sex_block,
    age_block,
    age2016_block,
    first_ct_age_block,
    first_ct_age_grp_block,
    age_dx_block,
    age_dx_grp_block,
    income_block,
    urban_block,
    ct_setting_block,
    ct_type_block,
    n_ct_block,
    outcome_row
  ) |>
    mutate(Study = study_label) |>
    select(Study, everything())
}

log_progress("===== 7.1b 建立 publication-style Table 1 =====")
study1_table1 <- build_table1(study1_cohort, "Study 1 (Population)")
study2_table1 <- build_table1(study2_cohort, "Study 2 (Sibling-Matched)")

# 各 study 各自存一份
study1_table1 |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study1-table1.csv"))
study2_table1 |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "study2-table1.csv"))

# 合併版（兩個 study 並排）
table1_combined <- bind_rows(study1_table1, study2_table1)
table1_combined |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "table1-baseline-characteristics.csv"))

log_progress("已輸出 Table 1：study1-table1.csv / study2-table1.csv / table1-baseline-characteristics.csv")

## 7.2 CT Type / Setting Distribution (Descriptive) --------------------------
# G-1：原 ct-type-distribution.csv 保留，新增 ct-type-by-setting.csv 雙向交叉表
# 【2026-05-10 修正】CT 分布只保留 0-18 歲時做的 CT 紀錄
#   （研究計畫的「兒童 CT」定義，>=19 歲的紀錄不算暴露）

# 共用前置：算出每筆 CT 在掃描日的年齡，過濾 0-18 歲
ct_records_pediatric <- ct_records |>
  filter(FUNC_DATE >= STUDY_START, FUNC_DATE <= STUDY_END) |>
  semi_join(base_cohort, by = "ID") |>
  left_join(
    pers_info |>
      mutate(BIRTH_DATE = resolve_birth_date(ID_BIRTHYM)) |>
      select(ID, BIRTH_DATE),
    by = "ID"
  ) |>
  mutate(AGE_AT_CT = calc_age_years(BIRTH_DATE, FUNC_DATE)) |>
  filter(AGE_AT_CT >= 0, AGE_AT_CT <= MAX_AGE_AT_INDEX)

log_progress(sprintf("§7.2 ct_records_pediatric (0-18 歲)：%s 筆",
                     format(nrow(ct_records_pediatric), big.mark = ",")))

ct_bodypart_summary <- ct_records_pediatric |>
  count(CT_TYPE) |>
  arrange(desc(n)) |>
  mutate(PCT = n / sum(n) * 100)

ct_bodypart_summary |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "ct-type-distribution.csv"))

# G-1 新增：CT_TYPE × CT_SETTING 交叉表
ct_type_setting_xtab <- ct_records_pediatric |>
  count(CT_TYPE, CT_SETTING) |>
  pivot_wider(
    names_from   = CT_SETTING,
    values_from  = n,
    names_prefix = "n_",
    values_fill  = 0L
  ) |>
  mutate(TOTAL = rowSums(across(starts_with("n_")), na.rm = TRUE))

ct_type_setting_xtab |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "ct-type-by-setting.csv"))

# G-1 額外診斷：ER vs Clinic 各 CT 代碼的比例
ct_setting_by_code <- ct_records_pediatric |>
  count(ORDER_CODE, CT_SETTING) |>
  pivot_wider(
    names_from   = CT_SETTING,
    values_from  = n,
    names_prefix = "n_",
    values_fill  = 0L
  )
ct_setting_by_code |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "ct-code-by-setting.csv"))

## 7.3 Figure: CT Exposure by Body Part and Age Group -------------------------
# 沿用 §7.2 的 ct_records_pediatric（0-18 歲）

ct_age_bodypart_fig <- ct_records_pediatric |>
  mutate(
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
    select(CATEGORY, SUBCATEGORY, N_PERSONS, N_RECORDS, NOTE)
  # Study 3 暫停執行，不納入 cohort-counts
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
# Study 3 cohort 暫停執行，跳過

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
  stratum_counts(study2_cohort_counts, "4. Study 2 × Urban Status",    "URBAN_STATUS")
  # Study 3 strata 暫停執行
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
    )
  # Study 3 outcome 暫停執行
)

## 8.5b Sensitivity Cohorts 人數（H-2 / H-4） -------------------------------

counts_sensitivity <- bind_rows(
  tibble(
    CATEGORY    = "6. Sensitivity Cohorts",
    SUBCATEGORY = "Study 1 main (latency 2y, fixed index)",
    N_PERSONS   = nrow(study1_cohort),
    N_RECORDS   = nrow(study1_cohort),
    NOTE        = sprintf("OUTCOME=%s; EXPOSED=%s",
                          format(sum(study1_cohort$OUTCOME), big.mark = ","),
                          format(sum(study1_cohort$EXPOSED), big.mark = ","))
  ),
  tibble(
    CATEGORY    = "6. Sensitivity Cohorts",
    SUBCATEGORY = "Study 1 typespec (heme 2y / solid 5y)",
    N_PERSONS   = nrow(study1_cohort_typespec),
    N_RECORDS   = nrow(study1_cohort_typespec),
    NOTE        = sprintf("OUTCOME=%s",
                          format(sum(study1_cohort_typespec$OUTCOME), big.mark = ","))
  ),
  tibble(
    CATEGORY    = "6. Sensitivity Cohorts",
    SUBCATEGORY = "Study 1 uniform 5y latency",
    N_PERSONS   = nrow(study1_cohort_5y),
    N_RECORDS   = nrow(study1_cohort_5y),
    NOTE        = sprintf("OUTCOME=%s",
                          format(sum(study1_cohort_5y$OUTCOME), big.mark = ","))
  ),
  tibble(
    CATEGORY    = "6. Sensitivity Cohorts",
    SUBCATEGORY = "Study 1 extended 2000-2023 (typespec)",
    N_PERSONS   = nrow(study1_cohort_extended),
    N_RECORDS   = nrow(study1_cohort_extended),
    NOTE        = sprintf("OUTCOME=%s",
                          format(sum(study1_cohort_extended$OUTCOME), big.mark = ","))
  ),
  tibble(
    CATEGORY    = "6. Sensitivity Cohorts",
    SUBCATEGORY = "Study 1 time-varying long-format rows",
    N_PERSONS   = safe_n_distinct(study1_tv$ID),
    N_RECORDS   = nrow(study1_tv),
    NOTE        = "Long-format split rows; N_PERSONS = distinct IDs"
  )
)

## 8.6 合併輸出 ---------------------------------------------------------------

cohort_counts_all <- bind_rows(
  counts_build,
  counts_attrition,
  counts_study_exposure,
  counts_strata,
  counts_outcome,
  counts_sensitivity
)

cohort_counts_all |>
  write_csv(file.path(OUTPUT_TABLE_PATH, "cohort-counts.csv"))

log_progress(sprintf("已儲存 cohort-counts.csv：%s 列", nrow(cohort_counts_all)))

cat("\n========== Cohort Counts Summary ==========\n")
print(cohort_counts_all, n = Inf)

log_progress("===== 全部完成 =====")