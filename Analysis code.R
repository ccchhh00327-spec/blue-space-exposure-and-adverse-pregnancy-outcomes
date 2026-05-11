#Analytical codes directly associated with cohort participant information and codes for macro-scale health impact analyses are available at: 
#GitHub: https://github.com/ccchhh00327-spec/blue-space-exposure-and-adverse-pregnancy-outcomes.git
#logistic regression/Cox regression analyses
library(survival)
run_logistic_beta_se <- function(data, outcome_vars, main_predictors, covariates = NULL) {
  results_list <- list()
  
  for (yvar in outcome_vars) {
    for (main_pred in main_predictors) {
      
      # ---------- crude model（无协变量） ----------
      formula_crude <- as.formula(paste(yvar, "~", main_pred))
      model_crude <- glm(formula_crude, data = data, family = binomial)
      
      tidy_crude <- broom::tidy(model_crude) %>%
        filter(term == main_pred) %>%
        mutate(
          outcome = yvar,
          predictor = main_pred,
          model = "crude"
        ) %>%
        select(model, outcome, predictor, estimate, std.error, p.value) %>%
        rename(beta = estimate, se = std.error)
      
      results_list[[paste(yvar, main_pred, "crude", sep = "_")]] <- tidy_crude
      
      
      # ---------- adjusted model（有协变量） ----------
      if (!is.null(covariates)) {
        formula_adj <- as.formula(
          paste(yvar, "~", main_pred, "+", paste(covariates, collapse = " + "))
        )
      } else {
        formula_adj <- as.formula(paste(yvar, "~", main_pred))
      }
      
      model_adj <- glm(formula_adj, data = data, family = binomial)
      
      tidy_adj <- broom::tidy(model_adj) %>%
        filter(term == main_pred) %>%
        mutate(
          outcome = yvar,
          predictor = main_pred,
          model = "adjusted"
        ) %>%
        select(model, outcome, predictor, estimate, std.error, p.value) %>%
        rename(beta = estimate, se = std.error)
      
      results_list[[paste(yvar, main_pred, "adjusted", sep = "_")]] <- tidy_adj
      
    }
  }
  
  final_res <- bind_rows(results_list)
  return(final_res)
}
main_predictors <- c("MNDWI_buffer500_prepreg1yr", "MNDWI_buffer500_prepregT3", "MNDWI_buffer500_prepregP6","MNDWI_buffer500_prepregT2","MNDWI_buffer500_prepregT1","MNDWI_buffer500_pregT1","MNDWI_buffer500_pregT2","MNDWI_buffer500_pregG6","MNDWI_buffer500_pregfull","MNDWI_buffer500_p12")
outcome_vars <- c("matern_cvd","matern_gest_htn","matern_preecl","ap_anx","ap_deps","ap_ment_dis","pp_deps","perin_ment_dis")
covariates <- c("age", "edu_band", "incom", "ethn_bin","temp_pregfull","precip_pregfull")
res <- run_logistic_beta_se(
  data = water_data,
  outcome_vars = outcome_vars,
  main_predictors = main_predictors,
  covariates = covariates
)

run_cox_beta_se <- function(data, outcome_vars, main_predictors, covariates = NULL) {
  results_list <- list()
  
  for (yvar in outcome_vars) {
    for (main_pred in main_predictors) {
      
      # ---------- crude model ----------
      formula_crude <- as.formula(
        paste("Surv(gest_age,", yvar, ") ~", main_pred)
      )
      
      model_crude <- coxph(formula_crude, data = data)
      
      tidy_crude <- broom::tidy(model_crude) %>%
        filter(term == main_pred) %>%
        mutate(
          outcome = yvar,
          predictor = main_pred,
          model = "crude"
        ) %>%
        select(model, outcome, predictor, estimate, std.error, p.value) %>%
        rename(beta = estimate, se = std.error)
      
      results_list[[paste(yvar, main_pred, "crude", sep = "_")]] <- tidy_crude
      
      
      # ---------- adjusted model ----------
      if (!is.null(covariates)) {
        formula_adj <- as.formula(
          paste("Surv(gest_age,", yvar, ") ~", 
                main_pred, "+", paste(covariates, collapse = " + "))
        )
      } else {
        formula_adj <- as.formula(
          paste("Surv(gest_age,", yvar, ") ~", main_pred)
        )
      }
      
      model_adj <- coxph(formula_adj, data = data)
      
      tidy_adj <- broom::tidy(model_adj) %>%
        filter(term == main_pred) %>%
        mutate(
          outcome = yvar,
          predictor = main_pred,
          model = "adjusted"
        ) %>%
        select(model, outcome, predictor, estimate, std.error, p.value) %>%
        rename(beta = estimate, se = std.error)
      
      results_list[[paste(yvar, main_pred, "adjusted", sep = "_")]] <- tidy_adj
      
    }
  }
  
  final_res <- bind_rows(results_list)
  return(final_res)
}
main_predictors <- c("MNDWI_buffer500_prepreg1yr", "MNDWI_buffer500_prepregT3", "MNDWI_buffer500_prepregP6","MNDWI_buffer500_prepregT2","MNDWI_buffer500_prepregT1","MNDWI_buffer500_pregT1","MNDWI_buffer500_pregT2","MNDWI_buffer500_pregG6","MNDWI_buffer500_pregfull","MNDWI_buffer500_p12")
outcome_vars <- c("ptb","lbw","sga","stl_bth")
covariates <- c("age", "edu_band", "incom", "ethn_bin","temp_pregfull","precip_pregfull")
res_cox <- run_cox_beta_se(
  data = df,
  outcome_vars = outcome_vars,
  main_predictors = main_predictors,
  covariates = covariates
)

#forest plot
beta<-res%>%mutate(OR=exp(0.01*beta),lower=exp(0.01*(beta-1.96*se)),upper=exp(0.01*(beta+1.96*se)))
df_matern_cvd <- beta %>%
  filter(outcome == 'matern_cvd')
df_matern_cvd <- df_matern_cvd %>%
  mutate(
    predictor = factor(predictor),             
    model = factor(model, levels = c("adjusted", "crude")),  # crude vs adjusted
    shape_group = model,
    color_group = predictor
  )
df_matern_cvd$predictor <- factor(df_matern_cvd$predictor, levels = rev(c(
  "MNDWI_buffer500_prepreg1yr",
  "MNDWI_buffer500_prepregT3",
  "MNDWI_buffer500_prepregP6",
  "MNDWI_buffer500_prepregT2",
  "MNDWI_buffer500_prepregT1",
  "MNDWI_buffer500_pregT1",
  "MNDWI_buffer500_pregT2",
  "MNDWI_buffer500_pregG6",
  "MNDWI_buffer500_pregfull",
  "MNDWI_buffer500_p12"
)))
library(ggplot2)
predictor_levels <- levels(df_matern_cvd$predictor)

custom_colors <- c(
  "MNDWI_buffer500_prepreg1yr" = "#999EA2", 
  "MNDWI_buffer500_prepregT3"  = "#A39183",
  "MNDWI_buffer500_prepregP6"  = "#CCD2CC", 
  "MNDWI_buffer500_prepregT2"  = "#DBD2C9",
  "MNDWI_buffer500_prepregT1"  = "#976666", 
  "MNDWI_buffer500_pregT1"     = "#C09D9B", 
  "MNDWI_buffer500_pregT2"     = "#BEBEBE", 
  "MNDWI_buffer500_pregG6"     = "#7A848D", 
  "MNDWI_buffer500_pregfull"   = "#A9B7AA", 
  "MNDWI_buffer500_p12"        = "#D4C3AA"  
)
library(colorspace)
custom_colors_darker <- darken(custom_colors, amount = 0.2)

ggplot(df_matern_cvd, aes(
  x = predictor,
  y = OR,
  ymin = lower,
  ymax = upper,
  color = color_group,
  shape = shape_group
)) +
  geom_pointrange(
    position = position_dodge(width = 0.6),
    size = 1.0
  ) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  coord_flip() +
  scale_shape_manual(
    values = c("crude" = 1, "adjusted" = 16),
    name = "Model"
  ) +
  scale_color_manual(
    values = custom_colors_darker, 
    name = "Gestational period"
  ) +
  guides(
    color = guide_legend(order = 1),   # Gestational period 在上
    shape = guide_legend(order = 2)    # Model 在下
  ) +
  labs(
    x = "",
    y = "Odds Ratio (95% CI)") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )+scale_y_continuous(breaks = c(0.990, 0.995, 1.000))+scale_x_discrete(
    labels = c(
      "MNDWI_buffer500_prepreg1yr" = "Preconception 0-12 months", 
      "MNDWI_buffer500_prepregT3"  = "Preconception 6-12 months",
      "MNDWI_buffer500_prepregP6"  = "Preconception 0-6 months", 
      "MNDWI_buffer500_prepregT2"  = "Preconception 3-6 months",
      "MNDWI_buffer500_prepregT1"  = "Preconception 0-3 months", 
      "MNDWI_buffer500_pregT1"     = "First Trimester", 
      "MNDWI_buffer500_pregT2"     = "Second Trimester", 
      "MNDWI_buffer500_pregG6"     = "Conception 0-6 months", 
      "MNDWI_buffer500_pregfull"   = "Full Pregnancy", 
      "MNDWI_buffer500_p12"        = "Periconception 12 months"  
    ))


#DLNM+Cox
library(dlnm)
#######对出生数据库和暴露数据库进行匹配
df<-merge(monthly_mndwi_list_buffer500,zebra[,c("seq","ptb","lbw","sga","stl_bth","age","edu_band","ethn_bin","incom","temp_pregfull","precip_pregfull")],by="seq")

df<-df[order(df$seq),] 
df$id <- as.integer(factor(df$seq, levels = unique(df$seq)))  # 重新连续编号
View(df)
#######确定终止时间，即出现结局的时间
ftime <- sort(unique(df$preg_month[df$lbw==1])) 
ftime<-unique(floor(ftime)) 
list(ftime)
#将原始数据中的每个区间在给定的点上进行切割，其中gestation表示生存时间, lbw表示结局
birthspl <-survSplit(Surv(preg_month,lbw)~., df[,c("id","seq","preg_month","lbw","age","edu_band","ethn_bin","incom","temp_pregfull","precip_pregfull")], cut=ftime, start="gesst",end="gesexit",event="lbw")
#######对数据集进行排序，针对孕周 创建每周的风险集
birthspl <- birthspl[order(birthspl$id,birthspl$gesexit),] 
birthspl$riskset <- as.numeric(factor(birthspl$gesexit,levels=ftime))
##建立EXPOSURE PROFILE暴露框架，产生矩阵
data_mndwi <- df[order(df$id),]
data_mndwi <- df[,c(4:26)]
fexpfull <- function(exp) {
  expfull <- as.numeric(exp)
  names(expfull) <-1:23
  return(expfull)}
rdata_mndwi <- t(apply(data_mndwi,1,function(x) fexpfull(x[1:23])))
###建立暴露历史矩阵 matrixofexposurehistories
rdata_mndwi<- do.call(rbind,lapply(seq(nrow(birthspl)),
                                   function(i) exphist(rdata_mndwi[birthspl$id[i],],(birthspl$gesexit[i]),c(0,22))))
####建立交叉基函数，定义滞后效应中样条函数的节点
cb <-crossbasis(rdata_mndwi,lag=c(0,22),argvar=list(fun="ns",df=4),arglag=list(fun="ns",df=3))
###birthspl<-merge(birthspl,zebra[,c("seq","age","edu_band","ethn_bin","incom","temp_pregfull","precip_pregfull")],by="seq",all.x=T)
##构建模型
model1<- coxph(Surv(gesst,gesexit,lbw)~cb+age+incom+edu_band+ethn_bin+temp_pregfull+precip_pregfull, birthspl)
##预测模型
pred1 <- crosspred(cb,model1,cen=0,cumul=T, by=0.005)
########可视化展示
par(mfrow = c(1, 3))
plot(pred1, var = 0.01, col = "red",
xlim = c(0, 15),
ylim = c(0.85,1.05),
ylab = "Hazard Ratio (95%CI)",
main = "MNDWI=0.01")
plot(pred1, var = 0.02, col = "red",
xlim = c(0, 15),
ylim = c(0.85,1.05),
main = "MNDWI=0.02")
plot(pred1, var = 0.05, col = "red",
xlim = c(0, 15),
ylim = c(0.85,1.05),
main = "MNDWI=0.05")


#rcs+kernel density
dens <- density(water_data$MNDWI_buffer500_p12, na.rm = TRUE)
dens_df <- data.frame(
  x = dens$x,
  density = dens$y
)
hr_max<-1
scale_factor <- (hr_max * 0.50) / max(dens_df$density)
dens_df$density_scaled <- dens_df$density * scale_factor
library(rms)
dd<-datadist(water_data)
dd$limits$MNDWI_buffer500_p12[2]<-min(water_data$MNDWI_buffer500_p12,na.rm = T)
options(datadist = 'dd')
fit<-cph(Surv(gest_age,lbw)~rcs(MNDWI_buffer500_p12,3)+age+ethn_bin+edu_band+incom+temp_pregfull+precip_pregfull,data=water_data,x=T,y=T)
pred<-rms::Predict(fit,MNDWI_buffer500_p12,fun=exp,type="predictions",ref.zero = T,conf.int = 0.95)

p <- ggplot() +geom_area(
  data = dens_df,aes(x = x, y = density_scaled),
  fill = "#AAB5C7",
  alpha = 0.2
) +geom_ribbon(
  data = pred,
  aes(x = MNDWI_buffer500_p12, ymin = lower, ymax = upper),
  fill = "#2C5C7A",
  alpha = 0.2
) +
  geom_line(
    data = pred,
    aes(x = MNDWI_buffer500_p12, y = yhat),
    linewidth = 1,
    color = "#2C5C7A"
  ) +geom_boxplot(
    data = water_data,
    aes(x = MNDWI_buffer500_p12, y = -0.02),
    width = 0.05,
    fill = "#AAB5C7",
    color = "#AAB5C7",
    alpha = 0.8,
    outlier.shape = NA
  ) +
  scale_y_continuous(
    name = "HR (95% CI)",
    limits = c(-0.2, hr_max),
    sec.axis = sec_axis(~ . / scale_factor, name = "Kernel Density")
  ) +coord_cartesian(xlim = c(0, 0.2)) +
  theme_classic() +theme(axis.title.y.right = element_text(color = "black"),axis.text.y.right = element_text(color = "black")) 
plot(p)


