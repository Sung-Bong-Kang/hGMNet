#!/usr/bin/perl -w
# rsID P : from stdin
# ID gene-list(comma-separated) : as the 1st command-line argument, db18_20000.txt.gz from the GSA-SNP package
# k-th best : as the 2nd argument (default, 2nd)
#
my $genelistfile = shift;
my $genelist = loadGenes( $genelistfile );

my $kth = shift || 2;

my $geneP;
my $nin = 0;
my $nun = 0;
while( <> ) {
	next unless /^rs\d+/;
	chomp; chop if substr($_,-1,1) eq "\r";
	my ($rs, $p) = split /\s/;
	my $id = substr($rs,2);
	unless( defined $genelist->{$id} ) {
		$nun++;
		next;
	}
	foreach my $g ( @{$genelist->{$id}} ) {
		push @{$geneP->{$g}}, [$p,$rs];
	}
	$nin++;
}

print STDERR $nin, " SNPs looked up, ", $nun, " unmapped\n";
print STDERR scalar( keys %$geneP ), " genes\n";

my @res = ();
foreach my $g ( keys %$geneP ) {
	my @sorted = sort { $a->[0] <=> $b->[0] } @{$geneP->{$g}};
	my $k = $kth - 1;
	$k = $#sorted if scalar(@sorted) < $kth;	# if less than k-th terms, the last one is filled in
	push @res, [ $g, @{$sorted[$k]} ];
}

foreach my $s ( sort { $a->[1] <=> $b->[1] } @res ) {
	print join("\t", @$s), "\n";
}

sub loadGenes {
	my $gzfile = shift;
	open(G, "gzip -dc $gzfile |") || die "$!\n";
	my $genelist;
	while( <G> ) {
		next unless /^\d+/;
		chomp; chop if substr($_,-1,1) eq "\r";
		my ($id, $genes) = split /\t/;
		my @genes = split /,/, $genes;
		$genelist->{$id} = \@genes;
	}
	close(G);
	print STDERR scalar( keys %$genelist ), " SNPs loaded into DB\n";
	return $genelist;
}
