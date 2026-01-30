# ==============================================================================
# PROJECT METADATA
# ==============================================================================
# Data Source: Worldscope - Fundamentals Annual
# Raw data:    https://docs.google.com/spreadsheets/d/15aiHK6ctJdRPhf0NuW3pzQUEkRSxFRbgNtVE_3URcbA/edit?gid=1099793521#gid=1099793521
# Clean data:  https://docs.google.com/spreadsheets/d/1NZ63R4YSf6_f2EEIrgCiJepIDCCipQ5wExsLoVZuu5U/edit?gid=1690865623#gid=1690865623

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, fixest, ggplot2, modelsummary, psych, tidyr, gridExtra) 

# ==============================================================================
# 1. DATA IMPORT
# ==============================================================================
df_raw <- read.csv('/Users/ruslansmirnov/Downloads/G7 (excl Canada) - j5ttstwj1qq5hdcx.csv', 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)

# ==============================================================================
# 2. DATA CLEANING & WINSORIZATION
# ==============================================================================

df_clean <- df_raw %>%
  rename(
    company_id   = `Code`,                      
    year_val     = `Year_`,
    nation       = `NATION`,
    mkt_cap_raw  = `MARKET CAPITALIZATION (U.S.$)`,
    equity_raw   = `COMMON EQUITY (U.S.$)`,
    assets_raw   = `TOTAL ASSETS (U.S.$)`,
    leverage_raw = `TOTAL DEBT % COMMON EQUITY`,
    roe_raw      = `RETURN ON EQUITY - TOTAL (%)`,
    growth_raw   = `BOOK VALUE PER SHARE - 1 YR ANNUAL GROWTH`,
    div_raw      = `COMMON DIVIDENDS (CASH)`,
    sic          = `SIC CODE 1`
  ) %>%
  mutate(
    # Clean numeric columns
    across(c(mkt_cap_raw, equity_raw, roe_raw, assets_raw, growth_raw, leverage_raw, div_raw), 
           ~ as.numeric(gsub(",", ".", gsub("[^0-9,.-]", "", .)))),
    
    # Fix Dividends: Assume NA equals 0
    div_raw = ifelse(is.na(div_raw), 0, div_raw)
  ) %>%
  
  # Filter missing core data
  filter(!is.na(mkt_cap_raw), !is.na(equity_raw), !is.na(assets_raw)) %>%
  filter(equity_raw > 0, mkt_cap_raw > 0) %>% 
  filter(!is.na(roe_raw), !is.na(growth_raw), !is.na(leverage_raw), !is.na(sic)) %>%
  
  select(company_id, year_val, nation, mkt_cap_raw, equity_raw, roe_raw, 
         assets_raw, growth_raw, leverage_raw, div_raw, sic) %>%

  # --- WINSORIZATION STEP ---
  mutate(
    roe_w      = psych::winsor(roe_raw, trim = 0.01, na.rm = TRUE),
    growth_w   = psych::winsor(growth_raw, trim = 0.01, na.rm = TRUE),
    leverage_w = psych::winsor(leverage_raw, trim = 0.01, na.rm = TRUE),
    mkt_cap_w  = psych::winsor(mkt_cap_raw, trim = 0.01, na.rm = TRUE),
    equity_w   = psych::winsor(equity_raw, trim = 0.01, na.rm = TRUE),
    assets_w   = psych::winsor(assets_raw, trim = 0.01, na.rm = TRUE)
  ) %>%

  mutate(
    PB = mkt_cap_w / equity_w,
        PB_w = psych::winsor(PB, trim = 0.01, na.rm = TRUE),
    log_PB = log(PB_w),
    log_mkt_cap = log(mkt_cap_w), 
    log_equity  = log(equity_w),
    is_common_law = ifelse(nation %in% c("UNITED STATES", "UNITED KINGDOM", "USA", "GBR"), 1, 0),
    Legal_System  = ifelse(is_common_law == 1, "Common Law", "Code Law"),
    is_dividend_payer = ifelse(div_raw > 0, 1, 0),
    Size     = log(assets_w),                        
    Industry = as.factor(substr(sprintf("%04d", sic), 1, 2)), 
    Year_f   = as.factor(year_val)
  )

# ==============================================================================
# 3. CALCULATING WEIGHTS
# ==============================================================================
df_weighted <- df_clean %>%
  group_by(nation) %>%
  mutate(count_country = n()) %>%
  ungroup() %>%
  mutate(Weight_Var = 1 / count_country)

cat("\n--- DATA DISTRIBUTION (FULL SAMPLE) ---\n")
print(table(df_weighted$nation)) 

# ==============================================================================
# 4. MAIN MODELS (Using Winsorized Variables)
# ==============================================================================

models_PB <- list(
  "OLS (Unweighted)" = feols(log_PB ~ is_common_law + roe_w + growth_w + Size + leverage_w | Industry + Year_f, 
                             data = df_weighted),
  "WLS (Weighted)"   = feols(log_PB ~ is_common_law + roe_w + growth_w + Size + leverage_w | Industry + Year_f, 
                             data = df_weighted, 
                             weights = ~Weight_Var)
)

cat("\n--- H1: P/B RATIO RESULTS ---\n")
msummary(models_PB, stars = TRUE, title = "H1: P/B Ratio (Winsorized & Weighted)")

models_PB <- list(
  "OLS (Unweighted)" = feols(log_PB ~ is_common_law + roe_w + growth_w + Size + leverage_w | Industry + Year_f, 
                             data = df_weighted),
  "WLS (Weighted)"   = feols(log_PB ~ is_common_law + roe_w + growth_w + Size + leverage_w | Industry + Year_f, 
                             data = df_weighted, 
                             weights = ~Weight_Var)
)

cat("\n--- H1: P/B RATIO RESULTS ---\n")


# ==============================================================================
# 5. ADDITIONAL HYPOTHESES (H2 & H3)
# ==============================================================================

# H2: Market Cap
model_h2_wls <- feols(log_mkt_cap ~ is_common_law + Size + leverage_w | Industry + Year_f, 
                      data = df_weighted, 
                      weights = ~Weight_Var)

# H3: Book Value
model_h3_wls <- feols(log_equity ~ is_common_law + Size + leverage_w | Industry + Year_f, 
                      data = df_weighted, 
                      weights = ~Weight_Var)

models_extra <- list(
  "H2: Mkt Cap (WLS)"  = model_h2_wls,
  "H3: Book Val (WLS)" = model_h3_wls
)

cat("\n--- H2 & H3: COMPONENTS RESULTS ---\n")
msummary(models_extra, stars = TRUE, title = "Testing H2 & H3 (Winsorized)")

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================

# Plot 1: H1 (P/B Ratio)
p_pb <- ggplot(df_weighted, aes(x = Legal_System, y = log_PB, fill = Legal_System)) +
  geom_boxplot(outlier.alpha = 0.1, varwidth = TRUE) + 
  scale_fill_manual(values = c("Code Law" = "#e74c3c", "Common Law" = "#2c3e50")) +
  labs(title = "H1: P/B Ratio", subtitle = "Winsorized at 1%", y = "Log(P/B)", x = "") +
  theme_minimal() + theme(legend.position = "none", plot.title = element_text(face = "bold"))

# Plot 2: H2 (Market Cap)
p_cap <- ggplot(df_weighted, aes(x = Legal_System, y = log_mkt_cap, fill = Legal_System)) +
  geom_boxplot(outlier.alpha = 0.1, varwidth = TRUE) +
  scale_fill_manual(values = c("Code Law" = "#e74c3c", "Common Law" = "#2c3e50")) +
  labs(title = "H2: Market Cap", y = "Log(Market Cap)", x = "") +
  theme_minimal() + theme(legend.position = "none")

# Plot 3: H3 (Book Value)
p_bv <- ggplot(df_weighted, aes(x = Legal_System, y = log_equity, fill = Legal_System)) +
  geom_boxplot(outlier.alpha = 0.1, varwidth = TRUE) +
  scale_fill_manual(values = c("Code Law" = "#e74c3c", "Common Law" = "#2c3e50")) +
  labs(title = "H3: Book Value", y = "Log(Book Value)", x = "") +
  theme_minimal() + theme(legend.position = "none")

# Layout
grid.arrange(p_pb, ncol = 1, top = "Hypothesis 1: P/B Ratio Analysis")
grid.arrange(p_cap, p_bv, ncol = 2, top = "Hypotheses 2 & 3: Components Analysis")

# ==============================================================================
# 7. DESCRIPTIVE STATISTICS 
# ==============================================================================

df_table <- df_weighted %>%
  select(
    `P/B Ratio` = PB_w,
    `Market Capitalization (Log)` = log_mkt_cap,
    `Book Value of Equity (Log)` = log_equity,
    `Total Assets (Log)` = Size,
    `Return on Equity (ROE)` = roe_w,
    `Financial Leverage` = leverage_w,
    `Annual Growth` = growth_w,
    Legal_System
  )

datasummary(
  All(df_table) ~ Legal_System * (Mean + SD + Median),
  data = df_table,
  title = "Table 1: Descriptive Statistics by Legal Origin (Winsorized)",
  fmt = 2, 
  notes = list(
    "Notes: The sample consists of publicly traded firms from G7 countries (excluding Canada).",
    "Continuous variables are winsorized at the 1% level (trim=0.01) to mitigate outliers.",
    "P/B Ratio is Market Cap divided by Common Equity.",
    "Code Law includes France, Germany, Italy, Japan. Common Law includes USA and UK."
  )
)

if (!require("writexl")) install.packages("writexl")
library(writexl)

## OBSERVATIONS BY YEAR AND COUNTRY
obs_table <- df_weighted %>%
  count(nation, year_val) %>%
  pivot_wider(names_from = year_val, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(nation)

final_obs_table <- obs_table %>%
  bind_rows(summarise(., across(where(is.numeric), sum), nation = "Total"))

## COUNTRY WEIGHTS AND LEGAL SYSTEM SUMMARY
weights_table <- df_weighted %>%
  group_by(nation) %>%
  summarise(
    Legal_System = first(Legal_System),
    Observations = n(),
    Weight_Value = first(Weight_Var),
    Contribution_Pct = (Weight_Value * Observations) / sum(df_weighted$Weight_Var) * 100
  ) %>%
  arrange(Legal_System, desc(Observations))

cat("\n--- TABLE C: OBSERVATIONS BY YEAR ---\n")
print(final_obs_table)

cat("\n--- COUNTRY WEIGHTS AND LEGAL SYSTEM SUMMARY ---\n")
print(weights_table)
