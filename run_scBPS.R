
scRNA <- readRDS("scRNA_annotated.Rds")

source("scBPS_module.R")
run_scBPS_pipeline(
  scrna_obj = scRNA,
  celltype_col = "cellType_1",
  gwas = "Dutch207", # 默认是荷兰207个肠道菌群数据, 其他暂未处理
  output_dir = "./",
  Only_Plot = FALSE
)
