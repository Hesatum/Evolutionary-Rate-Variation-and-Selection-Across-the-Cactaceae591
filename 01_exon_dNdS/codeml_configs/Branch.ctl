      seqfile = LOCUS.phy
     treefile = tree.nwk
      outfile = LOCUS_Branch.txt

        noisy = 1              * 0-9: output detail
      verbose = 1              * More or less detailed report in outfile
      seqtype = 1              * Codon data
        ndata = 1              * Number of data sets or loci
        icode = 0              * Universal genetic code
    cleandata = 1              * Remove sites with ambiguity data?

        model = 2         * Models for omega varying across lineages
      NSsites = 0          * Models for omega varying across sites
    CodonFreq = 7        * Codon frequencies (F3x4)
      estFreq = 0              * Use observed frequencies
        clock = 0              * No molecular clock
    fix_omega = 0         * Estimate omega
        omega = 0.5        * Initial omega value
