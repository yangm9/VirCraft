$protein = $ARGV[0];
unless (-e "$ARGV[1]/homolog.fasta") {
    open OUT, ">", "$ARGV[1]/homolog.fasta" or die "Can not create file homolog.fasta, $!\n";
    open IN, $protein or die "Can not open file $protein, $!\n";
    $_ = <IN>;
    s/\s.*\s?$/\n/; s/[^\w\n>]/_/g;
    print OUT;
    while (<IN>) {
        if (m/^>/) {
            s/\s.*\s?$/\n/; s/[^\w\n>]/_/g;
            print OUT "\n$_";
        }
        else {
            s/\s+?$//g;
        $_ = uc($_);
            print OUT;
        }
    }
    close IN;
    close OUT;
}

