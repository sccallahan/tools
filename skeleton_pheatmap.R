map_hkg <- pheatmap(data_norm, cluster_rows = T, cluster_cols = T, 
                        show_rownames = F, show_colnames = F, fontsize_row = 18, 
                        legend = T, filename = NA, fontsize_col = 10, scale = "row", 
                        fontface="bold", clustering_distance_rows = "correlation", 
                        clustering_distance_cols = "euclidean", clustering_method = "ward.D2" ,
                        annotation = annot1, annotation_colors = anno_color,
                        color = viridis(n = 5))
