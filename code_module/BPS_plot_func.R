library(AUCell)
library(ggplot2)
library(reshape2)
library(dplyr)

BPS_plot <- function(BPS_AUC_path, BPS_AUC_P_path) {
     BPS_AUC <- read.table(
          BPS_AUC_path,
          header = TRUE, row.names = 1, sep = "\t",
          check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE,
          quote = ""
     )
     BPS_AUC_P <- read.table(
          BPS_AUC_P_path,
          header = TRUE, row.names = 1, sep = "\t",
          check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE,
          quote = ""
     )

     # 仅保留AUC列（去掉missing、ncells等非AUC信息）
     auc_cols <- setdiff(colnames(BPS_AUC), c("missing", "ncells"))
     BPS_AUC <- BPS_AUC[, auc_cols, drop = FALSE]

     auc_mat <- matrix(as.numeric(as.matrix(t(BPS_AUC))), nrow = 1)
     rownames(auc_mat) <- "BPS_AUC"

     cells_assignment <- AUCell_exploreThresholds(auc_mat, plotHist = FALSE, assign = TRUE)
     cutoff <- suppressWarnings(as.numeric(cells_assignment$BPS_AUC$aucThr$selected))

     # 出图
     pdf("1.BPS_AUC_histogram.pdf", width = 6, height = 4)
     hist_obj <- hist(
          as.vector(auc_mat),
          breaks = 207,
          col = "#adc7d9",
          border = NA,
          main = "",
          xlab = expression(BPS[AUC] ~ " histogram")
     )
     abline(v = cutoff, col = "red", lwd = 3)
     text(
          cutoff,
          max(hist_obj$counts),
          labels = paste0("Cutoff = ", round(cutoff, 3)),
          pos = 4,
          col = "#000000"
     )
     dev.off()

     # 1. 计算 FDR（按 trait 方向校正）
     BPS_FDR <- apply(BPS_AUC_P, 2, p.adjust, method = "fdr")

     # 2. 合并 AUC 和 FDR
     BPS_AUC$celltype <- rownames(BPS_AUC)
     AUC_long <- reshape2::melt(BPS_AUC,
          id.vars = "celltype",
          variable.name = "trait", value.name = "AUC"
     )
     FDR_long <- reshape2::melt(BPS_FDR, varnames = c("celltype", "trait"), value.name = "FDR")

     merged <- merge(AUC_long, FDR_long, by = c("celltype", "trait"))


     df <- merged
     # df <- df[is.finite(df$AUC) & is.finite(df$FDR) & df$FDR > 0, , drop = FALSE]
     # 计算 -log10(FDR)
     df$neglogFDR <- -log10(df$FDR)

     # 分类
     df$class <- with(df, ifelse(FDR < 0.01 & AUC > 0.5, "Stringent",
          ifelse(FDR < 0.01 & AUC > cutoff & AUC < 0.5, "Moderate",
               ifelse(FDR > 0.01 & FDR < 0.05 & AUC > cutoff & AUC < 0.5, "Lenient",
                    "Nonsig"
               )
          )
     ))

     N_sig <- sum(df$class %in% c("Stringent", "Moderate", "Lenient"))
     n_stringent <- sum(df$class == "Stringent")
     n_moderate <- sum(df$class == "Moderate")
     n_lenient <- sum(df$class == "Lenient")

     # 颜色
     cols <- c(
          "Stringent" = "#681c21",
          "Moderate"  = "#427497",
          "Lenient"   = "#619581",
          "Nonsig"    = "grey80"
     )

     # 统计数量
     table(df$class)

     # 绘图
     p_scatter <- ggplot(df, aes(x = AUC, y = neglogFDR, color = class)) +
          geom_point(alpha = 0.8, size = 2) +
          scale_color_manual(values = cols) +
          geom_vline(xintercept = 0.5, linetype = "dashed", color = "#722c2c") +
          geom_vline(xintercept = cutoff, linetype = "dashed", color = "#722c2c") +
          geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#722c2c") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#722c2c") +
          annotate("text",
               x = min(df$AUC, na.rm = TRUE) + 0.20,
               y = max(df$neglogFDR, na.rm = TRUE) * 1.06,
               label = paste0("Moderate (N=", format(n_moderate, big.mark = ","), ")"),
               hjust = 0, vjust = 1, color = cols["Moderate"], size = 4, fontface = "italic"
          ) +
          annotate("text",
               x = min(df$AUC, na.rm = TRUE) + 0.55,
               y = max(df$neglogFDR, na.rm = TRUE) * 1.06,
               label = paste0("Stringent (N=", format(n_stringent, big.mark = ","), ")"),
               hjust = 0, vjust = 1, color = cols["Stringent"], size = 4, fontface = "italic"
          ) +
          annotate("text",
               x = min(df$AUC, na.rm = TRUE) + 0.45,
               y = max(df$neglogFDR, na.rm = TRUE) * 0.50,
               label = paste0("Lenient (N=", format(n_lenient, big.mark = ","), ")"),
               hjust = 0, vjust = 1, color = cols["Lenient"], size = 4, fontface = "italic"
          ) +
          labs(
               x = expression(BPS[AUC] ~ " scores"),
               y = expression(-log[10](FDR)),
               color = "Class"
          ) +
          theme_bw(base_size = 16) +
          theme(
               panel.grid = element_blank(), # 去掉所有背景线条
               panel.border = element_blank(), # 去掉外边框
               axis.line = element_line(color = "black") # 只保留坐标轴
          )


     p_scatter <- p_scatter + theme(legend.position = "none")
     ggsave("2.BPS_AUC_scatter.pdf", p_scatter, width = 5, height = 5)

     bar_df <- df %>%
          filter(class %in% c("Stringent", "Moderate", "Lenient")) %>%
          group_by(celltype, class) %>%
          summarise(n = n(), .groups = "drop") %>%
          ungroup()
     # 计算每个 celltype 的总数
     order_df <- bar_df %>%
          group_by(celltype) %>%
          summarise(total = sum(n), .groups = "drop") %>%
          arrange(desc(total))
     # 按总数排序 celltype
     bar_df$celltype <- factor(bar_df$celltype, levels = order_df$celltype)
     bar_df$class <- factor(bar_df$class, levels = c("Stringent", "Moderate", "Lenient"))

     p_bar <- ggplot(bar_df, aes(x = celltype, y = n, fill = class)) +
          geom_bar(stat = "identity", width = 0.6) + # bar 更宽，距离更近
          scale_fill_manual(values = cols) +
          labs(x = "", y = "Number", fill = "") +
          theme_bw(base_size = 16) +
          theme(
               panel.grid = element_blank(), # 去掉所有背景线条
               panel.border = element_blank(), # 去掉外边框
               axis.line = element_line(color = "black"), # 只保留坐标轴
               axis.text.x = element_text(angle = 30, hjust = 1),
               legend.position = c(0.55, 1.2),
               legend.justification = c(0, 1),
               legend.background = element_rect(fill = alpha("white", 0), color = "NA"),
               legend.key = element_rect(fill = alpha("white", 0), color = NA)
          ) +
          scale_y_continuous(expand = c(0, 0))

     width <- 2.5 + 0.25 * (length(unique(bar_df$celltype)))
     ggsave("3.BPS_AUC_bar.pdf", p_bar, width = width, height = 4)

     # 数量最多的细胞为关键细胞
     key_cell <- bar_df %>%
          group_by(celltype) %>%
          summarise(total = sum(n), .groups = "drop") %>%
          arrange(desc(total)) %>%
          slice_head(n = 5) %>%
          pull(celltype) %>%
          head(1)
     write.table(key_cell, "key_cell.txt", row.names = FALSE)

     df_key <- df %>%
          filter(celltype == key_cell) %>%
          arrange(desc(AUC))
     write.table(df_key, "key_cell_BPS_AUC.txt", row.names = FALSE)

     ## 复现环形图：关键细胞中不同 class 的数量分布
     pie_df <- df_key %>%
          filter(class %in% c("Stringent", "Moderate", "Lenient", "Nonsig")) %>%
          dplyr::count(class, name = "n")
          
     pie_df$class <- factor(pie_df$class, levels = c("Stringent", "Moderate", "Lenient", "Nonsig"))
     pie_df <- pie_df[order(pie_df$class), ]
     pie_df$label <- paste0("n = ", pie_df$n)

     total_n <- sum(pie_df$n)

     p_pie <- ggplot(pie_df, aes(x = 2, y = n, fill = class)) +
          geom_col(color = "black", width = 1) +
          coord_polar(theta = "y") +
          xlim(0.5, 2.5) +
          scale_fill_manual(values = cols) +
          geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "white", size = 5) +
          annotate("text",
               x = 0.55, y = 0, label = paste0("Total: ", total_n),
               color = "red", size = 6, fontface = "bold"
          ) +
          labs(fill = "") +
          theme_void(base_size = 16) +
          theme(legend.position = "right")
     ggsave("4.BPS_AUC_pie.pdf", p_pie, width = 6, height = 4)
}
