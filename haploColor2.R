
matr = alignment2matrix(aln)
tb = matr

keep = vector(length = ncol(tb))
for (i in 1:ncol(tb))
{	if (length(unique(tb[,i])) > 1)
	{	keep[i] = 1
	}
	else
	{	keep[i] = 0
	}
}
tb <- tb[,which(keep == 1)]

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
		#zeroes[remove] = 0
		nextSequence <- which(zeroes == max(zeroes))[1]
		colorMatrix[nextSequence,which(colorMatrix[nextSequence,] == 0)] = count
		for (i in 1:nrow(colorMatrix))  # for each sequence in alignment
		{	
			for (k in which(colorMatrix[nextSequence,] == count))  # for each position k that has not yet been assigned a color
			{	if (tb[i,k] == tb[nextSequence,k])
				{  	if (colorMatrix[i,k] == 0)
					{	colorMatrix[i,k] = count   # assign the 'color' (count) if the kth residue matches the kth residue in the reference
					}
				}
			}
		}
		#uncomment if you want to see an image for each iteration of the algorithm
	} 
	else
	{ continue = 0
	}
	if (count > 16)
	{	continue = 0
	}
}

residues = which(keep == 1)

domains = residues
domains[] = 0
domains[match(intersect(residues,1:543),residues)] = 1
domains[match(intersect(residues,544:829),residues)] = 2
domains[match(intersect(residues,830:1833),residues)] = 3
domains[match(intersect(residues,1833:length(tcdbvector)),residues)] = 4

domaincols = list(domains = c("0" = "white","1" = rgb(220/255,83/255,0),"2" = rgb(0,189/255,126/255),"3" = rgb(252/255,0,125/255),"4" = rgb(200/255,164/255,29/255)))

foo = as.data.frame(domains)
domain_annotation = HeatmapAnnotation(df = foo, col=domaincols,which="column", show_legend = TRUE)


ha = HeatmapAnnotation(df = as.data.frame(clusters), col=cols,which="row", show_legend = TRUE)
Heatmap(colorMatrix,row_split=numClusters,cluster_columns=F,cluster_rows = tree,col=list("0" = "gray", "1" = "black", "2" = "springgreen4","3" = "goldenrod2", "4" = "dark blue", "5" = "red", "6" = "light blue", "7" = "slateblue4", "8" = "purple", "9" = "green", "10" = "pink", "11" = "dodgerblue1","12" = "cyan", "13" = "firebrick", "14" = "brown", "15" = "peachpuff", "16" = "forestgreen", "17" = "red","18" = "steelblue2"),left_annotation = ha,bottom_annotation=domain_annotation)



