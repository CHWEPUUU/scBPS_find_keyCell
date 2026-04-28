stratified_sample_seurat <- function(scrna_obj, target_cells = 20000, 
                                    celltype_col = NULL) {
  set.seed(9527)
  
  meta <- scrna_obj@meta.data
  if (is.null(celltype_col)) {
    if ("cellType_1" %in% colnames(meta)) {
      celltype_col <- "cellType_1"
    } else {
      stop("未指定 celltype_col，且未找到默认的 cellType_1 列")
    }
  }

  # 如果细胞数不超过 target_cells，直接返回
  if (ncol(scrna_obj) <= 30000) {
    message("细胞数未超过 30000，无需抽样")
    return(scrna_obj)
  }

  meta$cell_id <- rownames(meta)

  # 2. 计算各 celltype 数量与比例
  strat_n <- as.data.frame(table(meta[[celltype_col]]), stringsAsFactors = FALSE)
  colnames(strat_n) <- c("strata", "n")
  strat_n$prop <- strat_n$n / sum(strat_n$n)

  # 3. 按比例分配样本数（floor）
  strat_n$alloc <- floor(target_cells * strat_n$prop)

  # 4. 最大余数法补齐到 target_cells
  need <- target_cells - sum(strat_n$alloc)
  if (need > 0) {
    strat_n$frac <- target_cells * strat_n$prop - strat_n$alloc
    idx <- order(strat_n$frac, decreasing = TRUE)[seq_len(need)]
    strat_n$alloc[idx] <- strat_n$alloc[idx] + 1
  }

  # 5. 不超过原始数量（保护稀有细胞）
  strat_n$alloc <- pmin(strat_n$alloc, strat_n$n)

  # 6. 分层抽样
  sampled_cells <- unlist(lapply(seq_len(nrow(strat_n)), function(i) {
    s <- strat_n$strata[i]
    n_sample <- strat_n$alloc[i]
    if (n_sample <= 0) return(character(0))
    ids <- meta$cell_id[meta[[celltype_col]] == s]
    sample(ids, n_sample)
  }))

  # 7. 返回抽样后的 Seurat 对象
  scRNA_sub <- subset(scrna_obj, cells = sampled_cells)

  message("已按 ", celltype_col, " 分层抽样：",
          ncol(scrna_obj), " -> ", ncol(scRNA_sub), " cells")

  return(scRNA_sub)
}
