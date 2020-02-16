# haploColor
Color a sequence alignment (nucleotide or protein) to visualize haplotype/recombination blocks

Given an initial alignment of variable sites (e.g., example.aln), haploColor.R will process it to facilitate visualization of 'haplotype' structures and recombination blocks.

### Color assignment algorithm:

1. Assign first sequence as reference.
2. Paint all residues of reference a unique color C.
3. Where other sequences match the reference, paint them color C
4. Identify sequence most dissimilar to the reference, and assign it as the new reference
5. Repeat steps 2-3 until all sequences are completely colored.

### Block assignment algorithm (in progress!):

This is a greedy algorithm that still has some issues.

For each sequence: 
  For its most common to least common colors:
    Assign that color to a block defined by its min to max position if its density is > a threshold
    


![](https://github.com/doxeylab/haploColor/raw/master/haploColor.png "haploColor visualization")



