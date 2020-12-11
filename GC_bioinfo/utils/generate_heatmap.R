#! /usr/bin/Rscript

library(RColorBrewer)
library(graphics)

args <- commandArgs(trailingOnly=TRUE)

input_matrix <- args[1]
color_palette <- args[2]
output_filename <- args[3]
gamma <- as.double(args[4])
min_value <- as.double(args[5])
max_value <- as.double(args[6])
center_value <- as.double(args[7])

matrix <- as.matrix(read.csv(input_matrix,header=FALSE,sep='\t'))

# Need to rotate the matrix for the image function
rotate <- function(x) t(apply(x, 2, rev))
data <- rotate(matrix)

dimensions <- dim(data)
rows <- dimensions[1]
cols <- dimensions[2]

if (color_palette == "gray") {
  my_palette <- gray.colors(256, start=0, end=1, gamma=gamma, rev=TRUE)
} else if (color_palette == "red/blue"){
  my_palette <- colorRampPalette(c("blue", "white", "red"))(255)
} else if (color_palette == "cyan/magenta"){
  my_palette <- colorRampPalette(c("magenta", "white", "cyan"))(255)
}

tiff(file=output_filename, width=rows, height=cols, units="px", res=600)

par(mar=c(0,0,0,0))

if (color_palette == "gray") {
  image(data, xaxt= "n", yaxt= "n", bty="n", col=my_palette)
} else {
  # Make the color scale from min_value to the max_value centered on a value (center value will be white for red/blue heatmaps)
  breaks <- c(seq(min_value, center_value, length.out=129), seq(center_value + (max_value/255), max_value, length.out=127))
  image(data, xaxt= "n", yaxt= "n", bty="n", col=my_palette, breaks=breaks)
}

