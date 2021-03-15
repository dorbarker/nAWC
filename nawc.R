suppressPackageStartupMessages(
    {
        usePackage <- function(p) {
            if (!is.element(p, installed.packages()[,1])) 
                install.packages(p, dep = TRUE)
            library(p, character.only = TRUE)
        }
        
        usePackage("ggplot2")
        usePackage("magrittr")
        usePackage("reshape2")
        usePackage("optparse")
        
        source("wallace.R")
    }
)


arguments <- function() {
    list(
        make_option(c("-i", "--input"),
                    help = "Input table of clusters over thresholds",
                    metavar = "FILE"),
        
        make_option(c("-o", "--outdir"),
                    help = "Output directory - will be created if it doesn't already exist",
                    metavar = "DIR"),
        
        make_option(c("-d", "--delimiter"),
                    default = "\t",
                    help = "Delimiter character for input [TAB]",
                    metavar = "CHAR")
    ) %>% 
    OptionParser(option_list = .) %>% 
    parse_args()
}

shannon <- compiler::cmpfun(function(clusters, base=2) {
    # Calculates the Shannon Entropy for a group of clusters
    # uses base 2, returning bits by default
    
    p_i <- function(i) {
        sum(clusters == i) / length(clusters)
    }
    
    N <- unique(clusters)
    
    -sum(sapply(N, function(i) p_i(i) * log(p_i(i), base = base)))
})

neighbour_awc <- compiler::cmpfun(function(clusts, i) {
    # Calculates the Adjusted Wallace Coefficent of the
    # ith column versus the i-1th column
    
    cur <- clusts[ , i]
    
    if ((i - 1) == 0 || length(unique(cur)) == 1) {
        result <- NA
        
    } else {
        
        suppressWarnings({
            result <- adj_wallace(cur, clusts[, i-1]) %>%
            use_series("Adjusted_Wallace_A_vs_B")
        })
    }
    
    result
})

singleton_proportion <- function(x) {
    # Determines what proportion of clusters
    # have only a single member genome
    
    singletons <- sum(table(x) == 1)
    
    p_singleton <- singletons / length(unique(x))
    
    p_singleton
    
}

stats_df <- function(clusters) {
    # Calculates various statistics for clusters,
    # and binds them together in a data.frame
    
    thresholds <-
        clusters %>%
        colnames %>%
        as.integer
    
    entropy <- sapply(clusters, shannon)
    
    nawc <-
        clusters %>%
        seq_along %>%
        sapply(neighbour_awc, clusts=clusters)
    
    n_clusts <- sapply(clusters, function(x) length(unique(x)))
    
    # p_singletons <- sapply(clusters, singleton_proportion)
    
    df <- data.frame("Threshold" = thresholds,
                     "Neighbour AWC" = nawc,
                     "Shannon (bits)" = entropy,
                     # "Number of Clusters" = n_clusts,
                     # "Proportion of Singleton Clusters" = p_singletons,
                     check.names = FALSE)
    df
}

process <- function(input_file, outdir, delimiter) {

    if (! dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }
   
    
    clusters <-
        input_file %>% 
        read.table(sep = delimiter,
                   row.names = 1, 
                   header = TRUE, 
                   check.names = FALSE)
    
    cluster_stats <- stats_df(clusters)
    
    generate_plot(cluster_stats, outdir)
    
    write.table(cluster_stats, 
                file = file.path(outdir, "stats.tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
}

generate_plot <- function(cluster_stats, outdir) {

    m <- melt(cluster_stats, id.vars = "Threshold")
    
    suppressWarnings({
        # because NAs raise warnings in plots,
        # but we're okay with the NAs being there
        
        p <- ggplot(m, aes(x = Threshold, y = value)) +
            geom_step() +
            facet_grid(variable ~ ., scales = "free_y", switch = "y") +
            labs(x = "goeBURST Threshold", y = "")
        
        ggsave(file.path(outdir, "stability_plot.png"),
               device = "png",
               plot = p,
               width = 16,
               height = 9,
               units = "in",
               dpi = 320)
    }) 
}

main <- function() {
    
    args <- arguments()
    
    process(args$input, args$outdir, args$delimiter)
    

}
main()
