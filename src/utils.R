rccToDat = function(fileName) {
  library(dplyr); library(tidyr);
  lines <- data.frame(values = readLines(fileName))
  dat <- suppressWarnings(separate(data = lines, col = values,
                                   sep = ",", into = c("CodeClass", "Name", "Accession",
                                                       "Count")))
  
  ind <- grep("<[A-Z]", dat$CodeClass)
  attr <- rep(NA, nrow(dat))
  for (i in 1:length(ind)) attr[ind[i]:nrow(dat)] <- grep("<[A-Z]",
                                                          dat$CodeClass, value = TRUE)[i]
  dat <- dat %>% mutate(CodeClass = paste(CodeClass, gsub(" ",
                                                          "", chartr("<>", "  ", attr)), sep = "_"), fileName = fileName)
  dat <- dat[-grep("<", dat$CodeClass), ]
  dat <- dat[!is.na(dat$Name), ]
  
  ## split flow cell data (properties) and biological (gene)
  ## data
  techDat <- dat[1:(grep("CodeClass", dat$CodeClass) - 1),
  ] %>% dplyr::select(-c(Accession:Count)) %>% spread(CodeClass,
                                                      Name)
  bioDat <- dat[(grep("CodeClass", dat$CodeClass) + 1):nrow(dat),
  ]
  
  ## combine techDat and bioDat
  Dat <- full_join(techDat, bioDat, by = "fileName")
  
  return(Dat)
}



#' customeTheme function for ggplot
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param sizeStripFont font of size of facet labels
#' @param xAngle angle of x-axis labels
#' @param hjust horizontal justification 0-left, 0.5-center, 1-right
#' @param vjust vertical justification 0-low, 0.5-middle, 1-high
#' @param xSize font size of x-axis label
#' @param ySize font size of y-axis label
#' @param xAxisSize font size of x-axis label title
#' @param yAxisSize fotn size of y-axis label title
#' @export
customTheme = function(sizeStripFont, xAngle, hjust, vjust, xSize,
                       ySize, xAxisSize, yAxisSize) {
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 1), strip.text.x = element_text(size = sizeStripFont),
        strip.text.y = element_text(size = sizeStripFont), axis.text.x = element_text(angle = xAngle,
                                                                                      hjust = hjust, vjust = vjust, size = xSize, color = "black"),
        axis.text.y = element_text(size = ySize, color = "black"),
        axis.title.x = element_text(size = xAxisSize, color = "black"),
        axis.title.y = element_text(size = yAxisSize, color = "black"),
        panel.background = element_rect(fill = "white", color = "black"))
}



#' Calculate desriptive statistics for each group and compare with lm() or lme()
#'
#' takes in data and determines the summary statistics (Mean & SD) and also compares the levels of the groups (binary) variable
#' @param demo is an (nxp) dataset
#' @param groups specifies the column name of a binary variable in demo
#' @param variables is a vector of column names to be compared between the 2 groups
#' @param paired boolean (T/F) - default = FALSE (repeated measures or not)
#' @param pairing is a column name that containing the pairing information
#' @export
descriptiveStat = function(demo, groups, variables, paired = FALSE, pairing = NULL){
  library(dplyr)
  library(tidyr)
  library(broom)
  
  if(all(paired)){
    X <- demo[, c(variables, groups, pairing), drop = FALSE]
    colnames(X) <- c(variables, "Group", "Pairing")
    lvls <- levels(X$Group)
    meanSD <- X %>% gather(Variable, Value, -c(Group, Pairing)) %>% dplyr::group_by(Variable,
                                                                                    Group) %>% dplyr::summarise(MEAN = mean(Value, na.rm = TRUE),
                                                                                                                SD = sd(Value, na.rm = TRUE))
    
    pval0 <- X %>% gather(Variable, Value, -c(Group, Pairing)) %>% dplyr::group_by(Variable) %>%
      nest() %>% dplyr::mutate(model = purrr::map(data,
                                                  ~tryCatch(lme(Value ~ Group, random = ~ 1 | Pairing, data = .), error = function(e) NA)
      ))
    pval <- do.call(rbind, lapply(pval0$model, function(i){
      tryCatch(summary(i)$tTable[2,], error = function(e) NA)
    })) %>%
      data.frame %>% mutate(Variable = variables, term = paste("Group", lvls[2]),
                            BH.FDR = p.adjust(p.value, "BH"))
  } else {
    X <- demo[, c(variables, groups), drop = FALSE]
    colnames(X) <- c(variables, "Group")
    lvls <- levels(X$Group)
    meanSD <- X %>% 
      gather(Variable, Value, -Group) %>% 
      dplyr::group_by(Variable, Group) %>%
      dplyr::summarise(MEAN = mean(Value, na.rm = TRUE),
                       SD = sd(Value, na.rm = TRUE))
    
    pval <- X %>% 
      tidyr::gather(Variable, Value, -Group) %>% 
      dplyr::group_by(Variable) %>%
      tidyr::nest() %>% 
      dplyr::mutate(model = purrr::map(data, ~ .x %>% 
                                         summarise(model = list(broom::tidy(lm(Value ~ Group)))))) %>%
      dplyr::select(-data) %>% 
      tidyr::unnest(model) %>%
      tidyr::unnest(model) %>%
      group_by(Variable) %>% dplyr::slice(2)
    pval$BH.FDR <- p.adjust(pval$p.value, "BH")
  }
  
  return(list(meanSD = meanSD, pval = pval))
}


compVar= function (demo, eset, variables, ncomp = 10)
{
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  pcaX <- prcomp(eset, scale. = TRUE, center = TRUE)
  pval <- do.call(rbind, lapply(variables, function(i) {
    apply(pcaX$x, 2, function(j) {
      predictor <- demo[, i]
      if (class(predictor) == "factor") {
        if (nlevels(predictor) == 2) {
          coef(summary(lm(as.numeric(j) ~ predictor)))[2,
                                                       "Pr(>|t|)"]
        }
        else {
          anova(lm(as.numeric(j) ~ predictor))["predictor",
                                               "Pr(>F)"]
        }
      }
      else {
        coef(summary(lm(as.numeric(j) ~ predictor)))[2,
                                                     "Pr(>|t|)"]
      }
    })
  }))
  rownames(pval) <- variables
  colnames(pval) <- paste(colnames(pval), paste0(round(100 *
                                                         (pcaX$sdev^2/sum(pcaX$sdev^2)), 1), "%"), sep = "-")
  pval <- pval[, 1:ncomp]
  pvalheatmap <- pval
  pvalheatmap[pvalheatmap < 0.01] <- 0.01
  pvalheatmap[pvalheatmap > 0.1] <- 1
  pvalheatmap[pvalheatmap > 0.01 & pvalheatmap < 0.05] <- 0.05
  pvalheatmap[pvalheatmap > 0.05 & pvalheatmap < 0.1] <- 0.1
  pvalheatmap[pvalheatmap == "0.01"] <- "p < 0.01"
  pvalheatmap[pvalheatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalheatmap[pvalheatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalheatmap[pvalheatmap == "1"] <- "p > 0.10"
  p <- pvalheatmap %>% as.data.frame %>% mutate(Variable = rownames(.)) %>%
    gather(Threshold, Value, -Variable) %>% mutate(Threshold = factor(Threshold,
                                                                      levels = unique(Threshold))) %>% mutate(Variable = factor(as.character(Variable),
                                                                                                                                levels = variables)) %>%
    mutate(Value = factor(Value, levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10"))) %>%
    ggplot(aes(Threshold, Variable)) +
    geom_tile(aes(fill = Value), colour = "white") + scale_fill_manual(values = rev(brewer.pal(n = 8,
                                                                                               name = "Blues")[c(2, 4, 6, 8)])) + customTheme(sizeStripFont = 10,
                                                                                                                                              xAngle = 40, hjust = 1, vjust = 1, xSize = 10, ySize = 10,
                                                                                                                                              xAxisSize = 10, yAxisSize = 10) + xlab("") + ylab("")
  return(list(pval = pval, pvalheatmap = pvalheatmap, p = p))
}


library(enrichR)
plot_pathway_network = function(genes=c(), pathway_db, fdr_cutoff=0.01, rm_pathways=NULL, title=""){
  # Core wrapping function
  ## https://stackoverflow.com/questions/20241065/r-barplot-wrapping-long-text-labels
  wrap.it <- function(x, len)
  { 
    sapply(x, function(y) paste(strwrap(y, len), 
                                collapse = "\n"), 
           USE.NAMES = FALSE)
  }
  
  
  # Call this function with a list or vector
  wrap.labels <- function(x, len)
  {
    if (is.list(x))
    {
      lapply(x, wrap.it, len)
    } else {
      wrap.it(x, len)
    }
  }
  
  ## run enrichment analysis
  enrichedAll <- enrichr(unique(genes), pathway_db)
  sig_gsets <- do.call(rbind, enrichedAll) %>%
    dplyr::mutate(database = rep(names(enrichedAll), sapply(enrichedAll, nrow))) %>%
    dplyr::filter(Adjusted.P.value < fdr_cutoff)
  
  ## rm certain pathways
  sig_gsets <- subset(sig_gsets, !(Term %in% rm_pathways))
  sig_gsets <- sig_gsets[order(sapply(strsplit(sig_gsets$Genes, ";"), length), decreasing = TRUE), ]
  
  ## commute number of common gene members
  mat <- matrix(0, nr = nrow(sig_gsets), nc = nrow(sig_gsets))
  rownames(mat) <- colnames(mat) <- sig_gsets$Term
  for(term1 in sig_gsets$Term){
    for(term2 in sig_gsets$Term){
      genes1 <- strsplit(subset(sig_gsets, Term == term1)$Genes, ";")[[1]]
      genes2 <- strsplit(subset(sig_gsets, Term == term2)$Genes, ";")[[1]]
      mat[term1, term2] <- length(intersect(genes1, genes2))
    }
  }
  rownames(mat) <- colnames(mat) <- wrap.labels(trimws(sapply(strsplit(rownames(mat), "WP"), function(i) i[[1]])), 30)
  
  ## plot network
  network <- simplify(graph_from_adjacency_matrix(mat, mode = "undirected", 
                                                  diag = FALSE, weighted = TRUE))
  
  # plot it
  nodeSize <- 20*sapply(strsplit(sig_gsets$Genes, ";"), length)/10
  E(network)$width <- E(network)$weight
  plot(network, 
       vertex.size = nodeSize, 
       layout = layout_in_circle(network), 
       edge.weight=100, main=title)
  
  return(list(sig_gsets=sig_gsets))
}


## join points to centroids
StatCentSeg <- ggplot2::ggproto("StatCentSeg", Stat,
                                compute_group = function(data, scales, params,
                                                         cfun=median) {
                                  data$xend <- cfun(data$x)
                                  data$yend <- cfun(data$y)
                                  return(data)
                                },
                                required_aes = c("x", "y")
)
stat_centseg <- function(mapping = NULL, data = NULL, geom = "segment",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, cfun=median, ...) {
  layer(
    stat = StatCentSeg, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, cfun = cfun, ...)
  )
}



zip_nPure = function(.x, .fields = NULL, .simplify = FALSE) {
  if (length(.x) == 0)
    return(list())
  if (is.null(.fields)) {
    if (is.null(names(.x[[1]]))) {
      .fields <- seq_along(.x[[1]])
    } else {
      .fields <- stats::setNames(names(.x[[1]]), names(.x[[1]]))
    }
  } else {
    if (is.character(.fields) && is.null(names(.fields))) {
      names(.fields) <- .fields
    }
  }
  out <- lapply(.fields, function(i) lapply(.x, .subset2, i))
  if (.simplify)
    out <- lapply(out, simplify_if_possible)
  out
}