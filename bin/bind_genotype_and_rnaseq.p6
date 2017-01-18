#!/usr/bin/env perl6
#

sub MAIN($sampleIDs-path, $bam-dir, *@timepoints) {
  my @file-paths = $bam-dir.IO.dir(test => / '.' bam $/).map(*.Str);
  my @sampleIDs = $sampleIDs-path.IO.lines;

  say "sample\ttimepoint\tbam";

  for @sampleIDs -> $sample-id {
    for @timepoints -> $timepoint {
      my $bam-path = select-bam $sample-id, @file-paths, $timepoint;
      say "$sample-id\t$timepoint\t{$bam-path.IO.abspath}";
    }
  }
}

sub select-bam($sampleID, @file-paths, $stimulation-time) {
  my %names;
  for @file-paths -> $path {
    given $path.IO.basename {
      # Bam file name is a (gluten specific) T-Cell Clone
      when /(TCC '-'? \d+) .* t(\d+)/ {
	%names.push($0 => $path.IO.abspath) if $1 == $stimulation-time;
      }
    }
  }
  for %names.kv -> $k, $v is rw {
    if $v.elems > 1 {
      # if a single TCC<ID> has multiple technical replicates,
      # return the BAM file which is first after sorting path names
      # (so semantically random, but technically deterministic)
      $v = @($v.sort)[0];
    }
  }
  return %names{$sampleID};
}
