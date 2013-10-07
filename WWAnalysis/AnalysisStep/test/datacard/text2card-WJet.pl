use strict;
use warnings;
my %collect = ();
foreach my $file (glob("WJetsYieldsData_ForTable_[01]j.txt")) {
    my ($j) = ($file =~ /WJetsYieldsData_ForTable(_.j)\./);
    open FILE, $file;
    my $header = <FILE>;
    $_ = <FILE>; s/^\s+//;
    my @columns = split(/\s+/, $_);
    while (<FILE>) {
        my @items = split(/\s+/, $_);
        my $mass = shift(@items);
        foreach my $i (0 .. $#columns) {
            $collect{$columns[$i].$j} .=  sprintf('%s  %s  %s', $mass, $items[$i*3], $items[$i*3+2])."\n";
        }
    }
    close FILE;
}
foreach (sort(keys(%collect))) {
    my $text = $collect{$_};
    s/mu/m/g; s/el/e/g;
    open CARD, "> WJet_$_.txt";
    print CARD $text;
    close CARD;
    print "> WJet_$_.txt\n";
}
