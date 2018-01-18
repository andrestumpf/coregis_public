#!/usr/bin/env perl
# ./gen_xml_mask.pl 200 500 mask.tif mask.xml

use warnings;
use strict;

my $sizex=$ARGV[0];
my $sizey=$ARGV[1];
my $mask_filename=$ARGV[2];
my $output_name=$ARGV[3];

my $content = &get_maks_xml_content($mask_filename, $sizex, $sizey);

# write the xml file
open(OUT, ">$output_name");
print OUT "$content";
close(OUT);


sub get_maks_xml_content()
{
	my ($local_name, $sizex, $sizey) = @_;
	my $text=<< "EOF";
<?xml version="1.0" ?>
<FileOriMnt>
     <NameFileMnt>./$local_name</NameFileMnt>
     <NombrePixels>$sizex $sizey</NombrePixels>
     <OriginePlani>0 0</OriginePlani>
     <ResolutionPlani>1 1</ResolutionPlani>
     <OrigineAlti>0</OrigineAlti>
     <ResolutionAlti>1</ResolutionAlti>
     <Geometrie>eGeomMNTFaisceauIm1PrCh_Px1D</Geometrie>
</FileOriMnt>
EOF
	return $text;

}