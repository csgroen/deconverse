
# Spatial ------------
#' Plot spatial correlation of cell type proportions
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay with `deconverse::deconvolute` results
#' @param method a method to plot. If NULL and only one method ran, that method will be plotted.
#' @param level a level of annotation to plot. If NULL, the finest level will be plotted.
#'
#' @return a ggplot of spatial correlation between cell types
#'
#' @import tidyverse
#' @importFrom corrr correlate autoplot
#' @export

plt_spatial_correlation <- function(spatial_obj, method = NULL, level = NULL) {
    res_df <- .get_spatial_results(spatial_obj, method, level)

    # Plot ---
    plt <- res_df %>%
        select(starts_with("frac_")) %>%
        rename_with(~ str_remove(., "frac_")) %>%
        correlate() %>%
        autoplot(triangular = "lower") +
        geom_text(aes(label = round(r, digits = 3)), size = 3) +
        scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)) +
        labs(fill = "Corr. coef.") +
        guides(fill = guide_colorbar()) +
        theme(legend.title = element_text(size = 10),
              axis.text = element_text(color = "black"))

    plt
}

#' Plot spatial correlation of cell type proportions
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay with `deconverse::deconvolute` results
#' @param method a method to plot. If NULL and only one method ran, that method will be plotted.
#' @param level a level of annotation to plot. If NULL, the finest level will be plotted.
#' @param min_prop minimal proportion of a cell type to be considered as occuring in a given spot
#'
#' @returns a ggnetwork plot of co-occurance of cell types in spots
#'
#' @import tidyverse
#' @importFrom network network set.vertex.attribute get.vertex.attribute
#' @importFrom ggnetwork ggnetwork geom_edges geom_nodes geom_nodetext_repel theme_blank
#' @export
plt_colocalization <-  function(spatial_obj, method = NULL, level = NULL, min_prop = 0.01) {
    res_df <- .get_spatial_results(spatial_obj, method, level)

    pop_occurance <- res_df %>%
        select(starts_with("frac_")) %>%
        rename_with(~ str_remove(., "frac_")) %>%
        mutate(across(.fns = ~ if_else(. < min_prop, 0, 1), .cols = everything()))

    spatial_nodes <- tibble(population = colnames(pop_occurance),
                            n_occurances = colSums(pop_occurance))

    spatial_edges <- combn(colnames(pop_occurance), 2) %>%
        t() %>%
        as.data.frame() %>%
        tibble() %>%
        rename(pop1 = V1, pop2 = V2) %>%
        rowwise() %>%
        mutate(co_occ_events = map2_dbl(pop1, pop2, .get_cooc, pop_occurance)) %>%
        filter(co_occ_events > 0)


    n <- network(spatial_edges)
    vertex_names <- get.vertex.attribute(n, "vertex.names")
    set.vertex.attribute(n,
                         "total_occurances",
                         pull(spatial_nodes, n_occurances, population)[vertex_names])

    set.seed(0)
    plt <- ggplot(ggnetwork(n), aes(x = x, y = y, xend =  xend, yend = yend)) +
        geom_edges(aes(alpha = co_occ_events, lwd = co_occ_events), alpha = 0.2, curvature = 0.1) +
        geom_nodes(aes(fill = total_occurances, size = total_occurances), pch = 21, color = "black") +
        geom_nodetext_repel(aes(label = vertex.names)) +
        theme_blank() +
        guides(size = "none") +
        labs(fill = "Total occurances", lwd = "Co-occurance events")

    return(list(plot = plt, co_occurance_table = spatial_edges))
}

#' Plot correlation between methods on spatial data
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay with `deconverse::deconvolute` results
#' @param level a level of annotation to plot. If NULL, the finest level will be plotted.
#' @param nrow number of rows of the plot facets
#' @param cor_fill_limits a vector with two numbers, representing the lower and
#' upper limits of the color scale for the correlation
#' @param cor_text_size a number, representing the size of the text of the
#' correlation in the plot tiles.
#'
#' @returns a list containing the plot and table for the comparisons
#'
#' @import tidyverse
#' @export
plt_method_correlation <- function(spatial_obj, level = 1, nrow = 1,
                                   cor_fill_limits = c(0,1), cor_text_size = 3) {

    # Get results
    assert(class(spatial_obj) == "Seurat")
    assert(!is.null(spatial_obj@tools[["deconverse"]]))
    res <- spatial_obj@tools$deconverse
    methods <- setdiff(names(res), c("major_population", "populations"))
    assert(length(methods) > 1)
    res <- res[methods]
    if(is.null(level)) {
        level <- length(res[[1]])
    }
    res_tb <- lapply(res, function(r) { r[[level]] }) %>% bind_rows()

    message("Undetected populations will be removed")
    cor_tb <- res_tb %>%
        pivot_wider(id_cols = c("sample"), names_from = "method", values_from = starts_with("frac")) %>%
        rename_with(~ str_remove(., "frac_")) %>%
        select(where(~ sum(. != 0, na.rm = TRUE) > 0)) %>%
        corrr::correlate() %>%
        corrr::stretch() %>%
        mutate(population1 = str_remove(x, "_.*$"),
               population2 = str_remove(y, "_.*$"),
               method1 = str_remove(x, "^.*_"),
               method2 =  str_remove(y, "^.*_")) %>%
        filter(population1 == population2,  method1 != method2, !is.na(r)) %>%
        select(population = population1, method1, method2, r) %>%
        mutate(methods = paste0(method1, "|", method2))

    plt <- cor_tb %>%
        mutate(cor = if_else(r < mean(cor_fill_limits), "low", "high")) %>%
        ggplot(aes(method1, population, fill = r)) +
        facet_wrap(~ method2, scales = "free", nrow = nrow) +
        geom_tile() +
        geom_text(aes(label = signif(r, 2), color = cor), size = cor_text_size) +
        scale_color_manual(values = c("low" = "white", "high" = "black")) +
        scale_fill_viridis_c(limits = cor_fill_limits, oob = scales::squish) +
        labs(x = "Method", y = "Cell population", fill = "Cor") +
        guides(color = "none") +
        theme_sparse3()

    return(list(plot = plt, table = cor_tb))
}

#' Plot percentage of concordance of major population between methods on spatial data
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay with `deconverse::deconvolute` results
#'
#' @returns a list containing the plot and table for the comparisons
#'
#' @import tidyverse
#' @export
plt_method_concordance <- function(spatial_obj) {
    assert(class(spatial_obj) == "Seurat")
    assert(!is.null(spatial_obj@tools[["deconverse"]]))

    # Get results
    major_pop <- spatial_obj@tools$deconverse$major_population
    res <- spatial_obj@meta.data[,major_pop] %>%
        rownames_to_column("umi") %>%
        tibble() %>%
        rename_with(~ str_remove(., "_major_population"))
    res_tb <-     res %>%
        pivot_longer(cols = -umi, names_to = "method", values_to = "major_population")
    methods <- unique(res_tb$method)

    method_concordance_df <- combn(methods, 2) %>%
        t() %>%
        as.data.frame()

    concord <- sapply(1:nrow(method_concordance_df), function(i) {
        method1 <- method_concordance_df[i,1]
        method2 <- method_concordance_df[i,2]

        concord <- res_tb %>%
            filter(method %in% c(method1, method2)) %>%
            group_by(umi) %>%
            summarize(concordance = length(unique(major_population)) == 1) %>%
            pull(concordance) %>% sum()/nrow(res)
    })
    method_concordance_df <- method_concordance_df %>%
        mutate(concordance = concord) %>%
        rename(method1 = V1, method2 = V2)


    plt <- method_concordance_df %>%
        mutate(concord_disc = if_else(concordance > 0.5, "high", "low"),
               method1 = fct_reorder(method1, concordance),
               method2 = fct_reorder(method2, concordance)) %>%
        ggplot(aes(method1, method2, fill = concordance)) +
        geom_tile() +
        geom_text(aes(label = paste0(signif(concordance,2)*100, "%"), color = concord_disc)) +
        scale_fill_viridis_c(option = "A", limits = c(0, 1), labels = scales::percent_format()) +
        scale_color_manual(values = c("low" = "white", "high" = "black")) +
        labs(x = "Methods", y = "Methods", fill = "Major population\nconcordance") +
        guides(color = "none") +
        theme_sparse3()

    return(list(plot = plt, table = method_concordance_df))
}





# Helpers ------------------
.get_spatial_results <- function(spatial_obj, method, level) {
    assert(class(spatial_obj) == "Seurat")
    assert(!is.null(spatial_obj@tools$deconverse))

    res <- spatial_obj@tools$deconverse
    methods <- subset(names(res), names(res) %in% deconvolution_methods())

    if(is.null(method) & length(methods) == 1) {
        method <- methods
        message("-- Showing results for `", method, "` method")
    } else if(!(method %in% methods)) {
        stop("`method` provided is not valid. If more than one method has been computed, `method` must be specified.")
    }
    res <- res[[method]]

    if(is.null(level)) {
        level <- length(res)
        if(level > 1) message("-- Showing results for annotation level ", level)
    } else {
        assert(level <= length(res) | level %in% paste0("l", 1:length(res)))
    }
    res_df <- res[[level]]
    return(res_df)
}

.get_cooc <- function(pop1, pop2, pop_occurance) {
    pop_occurance %>%
        select(pop1 = all_of(pop1), pop2 = all_of(pop2)) %>%
        filter(pop1 == 1, pop2 == 1) %>%
        count() %>%
        pull(n)
}
