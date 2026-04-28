library(ieugwasr)
library(VariantAnnotation)
library(biomaRt)
library(SeuratDisk)

run_scBPS_pipeline <- function(
    scrna_obj,
    celltype_col = "cellType_1",
    gwas = "Dutch207",
    output_dir = "./",
    Only_Plot = FALSE) {
  setwd(output_dir)
  dirs <- list.dirs(".", recursive = FALSE)
  matched <- dirs[grepl("scBPS", dirs)]

  if (length(matched) == 0) {
    dir.create("scBPS")
    setwd("scBPS")
  } else {
    setwd(matched[1])
  }

  if (!dir.exists("output")) dir.create("output")

  BPS_AUC_path <- "./output/BPS_AUC.txt"
  BPS_AUC_P_path <- "./output/pvalue_AUC.txt"
  if (Only_Plot) {
    if (!file.exists(BPS_AUC_path)) stop("BPS评分文件不存在: ", BPS_AUC_path)
    if (!file.exists(BPS_AUC_P_path)) stop("BPS评分文件不存在: ", BPS_AUC_P_path)
    # BPS评分绘图，筛选关键细胞-关键菌群
    source("./code_module/BPS_plot_func.R")
    BPS_plot(
      BPS_AUC_path = BPS_AUC_path,
      BPS_AUC_P_path = BPS_AUC_P_path
    )
    return(invisible(TRUE))
  } else {
    # 抽样转为scBPS需要格式
    if (file.exists("output/sc_adata.h5ad")) {
      message("sc_adata.h5ad already exists, skipping Seurat to h5ad conversion.")
    } else {
      message("Processing Seurat object and converting to h5ad format...")
      source("./code_module/stratified_sample_seurat.R")
      scRNA <- stratified_sample_seurat(scrna_obj)

      SaveH5Seurat(scRNA, filename = "output/sc_adata.h5Seurat", overwrite = TRUE)
      Convert("output/sc_adata.h5Seurat", dest = "h5ad", overwrite = TRUE)

      cell_anno <- data.frame(
        cell_id = colnames(scRNA),
        cell_annotation = scRNA[[celltype_col]]
      )
      write.table(
        cell_anno, "output/cell.annotation.txt",
        row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
      )
    }
    # 运行scBPS
    system2("bash",
      args = c(
        "./code_module/BPS_score.sh",
        gwas
      )
    )

    # BPS评分绘图，筛选关键细胞-关键菌群
    source("./code_module/BPS_plot_func.R")
    BPS_plot(
      BPS_AUC_path = BPS_AUC_path,
      BPS_AUC_P_path = BPS_AUC_P_path
    )

    invisible(TRUE)
  }
}
