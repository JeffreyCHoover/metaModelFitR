library(MetaModelFitR)
library(tidyverse)
library(here)
library(glue)
library(metafor)
library(readxl)
library(runjags)
library(rjags)

options(scipen = 999)

HDL <- read_excel("HDL.xlsx")
LDL <- read_excel("LDL.xlsx")

hdl_meta <- rma(data = HDL, yi = HDL$d, vi = HDL$V, method = "FE")
trimfill(hdl_meta)

tiff(filename = "hdl_funnel.tiff", width = 4, height = 4, units = "in",
     res = 300)
funnel(hdl_meta)
dev.off()

#mc <- monte_carlo_ma(meta = hdl_meta, fixed = TRUE, fileName = "hdl_monte")

ppmc <- ppmc_ma(fileName = "hdl_ppmc", meta = hdl_meta, fixed = TRUE)

#mc
ppmc


ldl_meta <- rma(data = LDL, yi = LDL$d, vi = LDL$V)
trimfill(ldl_meta)
trimfill(ldl_meta, side = "right")
regtest(ldl_meta)

tiff(filename = "ldl_funnel.tiff", width = 4, height = 4, units = "in",
     res = 300)
funnel(ldl_meta)
dev.off()

ppmc_ldl <- ppmc_ma(fileName = "ldl_meta", meta = ldl_meta, fixed = FALSE)

ppmc_ldl


test1 <- read_excel("contrived.xlsx", sheet = "balanced") %>%
  select(id, d, V)

test2 <- read_excel("contrived.xlsx", sheet = "outlier") %>%
  select(id, d, V)

test3 <- read_excel("contrived.xlsx", sheet = "subgroups") %>%
  select(id, d, V)

test4 <- read_excel("contrived.xlsx", sheet = "subgroups") %>%
  select(id, d, V) %>%
  filter(id < 7)

test5 <- read_excel("contrived.xlsx", sheet = "subgroups") %>%
  select(id, d, V) %>%
  filter(id > 6)

test1_meta <- rma(yi = d, vi = V, data = test1, method = "REML")
regtest(test1_meta, predictor = "vi")
trimfill(test1_meta, side = "left")
trimfill(test1_meta, side = "right")
funnel(test1_meta)
test2_meta <- rma(yi = d, vi = V, data = test2, method = "REML")
regtest(test2_meta, predictor = "vi")
trimfill(test2_meta, side = "left")
trimfill(test2_meta, side = "right")
funnel(test2_meta)
test3_meta <- rma(yi = d, vi = V, data = test3, method = "REML")
regtest(test3_meta, predictor = "vi")
trimfill(test3_meta, side = "left")
trimfill(test3_meta, side = "right")
tiff(filename = "test3-funnel.tiff", height = 4, width = 4, units = "in", res = 300)
funnel(test3_meta)
dev.off()
test4_meta <- rma(yi = d, vi = V, data = test4, method = "REML")
test5_meta <- rma(yi = d, vi = V, data = test5, method = "REML")
regtest(test5_meta, predictor = "vi")
trimfill(test5_meta, side = "left")
trimfill(test5_meta, side = "right")
funnel(test5_meta)

test1_ppmc <- ppmc_ma(fileName = "test1_ppmc", meta = test1_meta, fixed = FALSE)
test1_ppmc

#test1_mc <- monte_carlo_ma(meta = test1_meta, fixed = FALSE,
#                           fileName = "test1_monte")
#test1_mc

test2_ppmc <- ppmc_ma(fileName = "test2_ppmc", meta = test2_meta, fixed = FALSE)
test2_ppmc

#test2_mc <- monte_carlo_ma(meta = test2_meta, fixed = FALSE,
#                           fileName = "test2_monte")
#test2_mc

test3_ppmc <- ppmc_ma(fileName = "test3_ppmc", meta = test3_meta, fixed = FALSE)
test3_ppmc
#test3_mc <- monte_carlo_ma(meta = test3_meta, fixed = FALSE,
#                           fileName = "test3_monte")
#test3_mc

test4_ppmc <- ppmc_ma(fileName = "test4_ppmc", meta = test4_meta, fixed = FALSE)
test4_ppmc
#test4_mc <- monte_carlo_ma(meta = test4_meta, fixed = FALSE,
#                           fileName = "test4_monte")
#test4_mc

test5_ppmc <- ppmc_ma(fileName = "test5_ppmc", meta = test5_meta, fixed = FALSE)
test5_ppmc
#test5_mc <- monte_carlo_ma(meta = test5_meta, fixed = FALSE,
#                           fileName = "test5_monte")
#test5_mc

data_print <- c(test1_ppmc, test2_ppmc, test3_ppmc, test4_ppmc, test5_ppmc,
                ppmc)

saveRDS(data_print, file = "manuscript-output")
