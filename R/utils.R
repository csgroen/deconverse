#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

# Utility
left_df_join <- function(x, y, sep = ".") {
    #-- Fixed column names
    fixed_dfs <- .rename_columns(x,y,sep)
    x <- fixed_dfs[[1]]; y <- fixed_dfs[[2]]; rm(fixed_dfs)
    #-- Get y not in x
    missing_rownames <- rownames(x)[! rownames(x) %in% rownames(y)]
    #-- Fill missing rownames in y for cbind
    ydf_fill <- data.frame(
        matrix(NA,
               nrow = length(missing_rownames),
               ncol = ncol(y),
               dimnames = list(missing_rownames, colnames(y)))
    )
    new_y <- rbind(y, ydf_fill)
    #-- Merge
    merged_df <- cbind(x, new_y[rownames(x),])
    if(ncol(new_y) == 1) {
        colnames(merged_df)[ncol(merged_df)] <- colnames(new_y)
    }
    return(merged_df)
}

.rename_columns <- function(x, y, sep) {
    repeated_cols <- intersect(colnames(x), colnames(y))
    if(length(repeated_cols) > 0) {
        x <- ._fix_colnames(x, repeated_cols, "x", sep)
        y <- ._fix_colnames(y, repeated_cols, "y", sep)
    }
    return(list(x,y))
}

._fix_colnames <- function(df, repeated_cols, name, sep) {
    new_cols <- ifelse(colnames(df) %in% repeated_cols,
                       paste(colnames(df), name, sep = sep),
                       colnames(df))
    colnames(df) <- new_cols
    return(df)
}

