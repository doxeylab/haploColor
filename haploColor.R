#------------------------------------------------------------------
#haplocolor
#Feb 15, 2020
#Author: Andrew Doxey, doxey.uwaterloo.ca

#Summary: Given a  protein or nucleotide alignment, this algorithm will produce an image of uniquely haplotype/recombination blocks 
#Input file: an alignment file (fasta)

#Algorithm:
#Paint all residues of first sequence color 1
#Where other sequences match the first, paint them color 1
#Identify sequence most dissimilar to first
#Paint it color 2 where it has not yet been colored
#Where other sequences match the second, paint them color 2
#Identify sequence most dissimilar to second
#etc.
#repeat above until all sequences have been completely colored 

#------------------------------------------------------------------


# required packages 
# use install.packages() if not already installed
library(pheatmap)
library(plotrix)
library(seqinr)


# read in sequence alignment
# alignment should only include variable positions
aln <- read.alignment("example.aln",format="fasta")

# convert alignment to matrix format
tb <- as.matrix(aln)


#uncomment following 10 lines if you want to automatically remove invariant sites in alignment
#keep = vector(length = ncol(tb))
#for (i in 1:ncol(tb))
#{	if (length(unique(tb[,1])) > 1)
#	{	keep[i] = 1
#	}
#	else
#	{	keep[i] = 0
#	}
#}
#tb <- tb[,which(keep == 1)]


# initialize a second matrix to keep track of haplotype 'colors'
colorMatrix = matrix(ncol=ncol(tb),nrow=nrow(tb))

# fill matrix with zeroes
colorMatrix[] <- 0

# the first sequence in the matrix (tb[1,]) starts as the reference
# assign a "1" to all residues in the alignment matching the reference 
for (i in 1:nrow(tb))
{
	if (length(which(tb[i,] == tb[1,])) > 0)
	{
		colorMatrix[i,tb[i,] == tb[1,]] = 1
	}
}

#uncomment if you want to visualize this initial coloring
#pheatmap(colorMatrix,cluster_cols=F)


#repeat loop
#this will repeat the process above until there are no more colors to assign

continue = 1
count = 1
while (continue == 1)
{	count = count + 1

	zeroes <- vector(length = nrow(colorMatrix))

	# first we must count the number of 0s in each row (sequence)
	# the seq with the most 0s is most dissimilar from the last reference sequence
	for (i in 1:nrow(colorMatrix))
	{
		zeroes[i] <- length(which(colorMatrix[i,] == 0))
	}

	if (sum(zeroes) > 0)  # checks to make sure there are still sequences left with 0s; otherwise algorithm is done
	{

		#picks the next reference sequence. If there are ties, it will take the first of the group.
		nextSequence <- which(zeroes == max(zeroes))[1]


		for (i in 1:nrow(colorMatrix))  # for each sequence in alignment
		{	
			for (k in which(colorMatrix[nextSequence,] == 0))  # for each position k that has not yet been assigned a color
			{	if (tb[i,k] == tb[nextSequence,k])
				{  colorMatrix[i,k] = count   # assign the 'color' (count) if the kth residue matches the kth residue in the reference
				}
			}
		}
		#uncomment if you want to see an image for each iteration of the algorithm
		#pheatmap(colorMatrix,cluster_cols=F,cluster_rows=F)
	} else
	{ continue = 0
	}

}

#now, colors have been assigned to all residues
#the number of colors is stored in count

#below, we will create a plot of the alignment including all the colors
#but we set max colors to 30 colors b/c above that, it becomes hard to distinguish between colors
#or fewer if there are fewer than 30
colorMax = min(count,30)


colors <- colorMatrix
colors[] <- "gray"  #gray will be the default color
#the 30 color palette is included below
colPalette <- c("#201923", "#ffffff", "#fcff5d", "#7dfc00", "#0ec434", "#228c68", "#8ad8e8", "#235b54", "#29bdab", "#3998f5", "#37294f", "#277da7", "#3750db", "#f22020", "#991919", "#ffcba5", "#e68f66", "#c56133", "#96341c", "#632819", "#ffc413", "#f47a22", "#2f2aa0", "#b732cc", "#772b9d", "#f07cab", "#d30b94", "#edeff3", "#c3a5b4", "#946aa2", "#5d4c86")
for (i in 1:colorMax)
{
	colors[colorMatrix == i] <- colPalette[i]
}

#generate the colored matrix
color2D.matplot(colorMatrix,cellcolors=colors,border=NA)
