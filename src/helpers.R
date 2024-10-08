#'
#'
#'
#'
#'

vega_20_scanpy = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2",  
                   "#bcbd22","#17becf","#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5",
                   "#c49c94","#f7b6d2","#dbdb8d","#9edae5","#ad494a","#8c6d31")  

zeileis_26 = c("#023fa5","#7d87b9","#bec1d4","#d6bcc0","#bb7784","#8e063b","#4a6fe3",
               "#8595e1","#b5bbe3","#e6afb9","#e07b91","#d33f6a","#11c638","#8dd593",
               "#c6dec7","#ead3c6","#f0b98d","#ef9708","#0fcfc0","#9cded6","#d5eae7",
               "#f3e1eb","#f6c4e1","#f79cd4","#7f7f7f","#c7c7c7","#1CE6FF","#336600")

godsnot_64 = c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059",
               "#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF",
               "#997D87","#5A0007","#809693","#FEFFE6","#1B4400","#4FC601","#3B5DFF",
               "#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92",
               "#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299",
               "#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED",
               "#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062",
               "#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578",
               "#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C",
               "#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F",
               "#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757",
               "#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C",
               "#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F",
               "#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55",
               "#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C")

my_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  # axis.title.x=element_blank(),
  # axis.title.y=element_blank(),
  # legend.position="none",
  panel.background=element_blank(),
  # panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

my_theme_nolegened <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  # axis.title.x=element_blank(),
  # axis.title.y=element_blank(),
  legend.position="none",
  panel.background=element_blank(),
  # panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

my_blank_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  legend.position="none",
  panel.background=element_blank(),
  # panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)





#'
#'
#'

my_plot <- function( data, column, embedding = "tsne" ) { 
  
  ggplot(data, aes_string(x = paste0("X_", embedding, "1"), 
                          y = paste0("X_", embedding, "2"), 
                          color = column)) +
    geom_point(shape = 46) + 
    # my_theme_nolegened + 
    scale_color_gradient(low="lightgrey", high="darkgreen") +
    labs(title = column)
  
}



#'
#'
#'
#'


subset_reshape <- function( aprobs, ch, gmap, K, marker ) {
  
  chr_aprobs <- aprobs[[ch]]
  
  out <- lapply(rownames(chr_aprobs), function(r) {
    
    tmp <- round(t(chr_aprobs[r, , ]), 2)
    tmp_order <- t(apply(tmp, 1, order))
    allele1 <- LETTERS[tmp_order[,8]]
    allele2 <- LETTERS[tmp_order[,7]]
    color1 <- CCcolors[tmp_order[,8]]
    color2 <- CCcolors[tmp_order[,7]]
    homozygous <- apply(tmp, 1, function(x) any(x > 0.8))
    geno <- paste0(allele1, allele2)
    geno[homozygous] <- paste0(allele1[homozygous], allele1[homozygous])
    geno <- sapply(lapply(strsplit(geno, ""), sort), paste0, collapse="")
    
    tmp %>%
      data.frame() %>%
      rownames_to_column(var = "marker") %>%
      mutate(allele1 = allele1,
             allele2 = allele2,
             color1 = color1,
             color2 = color2,
             homozygous = homozygous,
             geno = geno,
             id = r) %>%
      select(id, marker, allele1, allele2, color1, color2, homozygous, geno)
    
  }) %>% 
    do.call(rbind, .) %>%
    left_join( gmap )
  
  
  if ( is.null(marker) ) {
    
    Kchr <- K[[ch]]
    Kd <- as.dist(Kchr)
    Khc <- hclust(Kd)
    
    out <- out %>%
      mutate(id = factor(id, levels = Khc$labels[Khc$order]),
             marker = factor(marker, levels = gmap$marker[gmap$chr == ch]))
    
  }else{
    id_order <- out$id[out$marker == marker][order(out$color1[out$marker == marker])]
    
    out <- out %>%
      mutate(id = factor(id, levels = id_order),
             marker = factor(marker, levels = gmap$marker[gmap$chr == ch]))
  }
  
  out
  
}




#'
#'
#'
#'
#'

markers_window <- function( gmap, marker, window = 1e6 ) {
  
  marker_pos <- gmap[gmap$marker == marker, "pos"]
  marker_chr <- gmap[gmap$marker == marker, "chr"]
  
  gmap %>% 
    filter( between(pos, marker_pos - window, marker_pos + window), 
            chr == marker_chr ) %$%
    marker
  
}




#'
#'
#'
#'
#'


plot_founder <- function( df, markers ) {
  
  df %>%
    filter(marker %in% markers) %>%
    ggplot(., aes(x = marker, y = id, fill = I(color1))) + 
    geom_tile()
  
}

