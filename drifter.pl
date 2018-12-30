#!/usr/bin/env perl

# Programme: drifter v.2.2
# A programme to simulate changing allele frequencies through time
# Institution: Arizona Research Laboratories, University of Arizona
# Date: 2005
# Author: Murray Cox
# Developed on perl v5.8.7 build for MSWin32-x86-multi-thread

# Confirmed still working with perl v5.18.2
# December 2018

# Copyright (C) 2005 Murray P. Cox

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Time::localtime;
use Getopt::Long;

# Disable filehandle buffering on STDOUT
my $stdout_filehandle = select(STDOUT);
$| = 1;
select($stdout_filehandle);

# Global parameters
my $NUMBER_OF_GENERATIONS = undef;
my $NUMBER_OF_ITERATIONS = undef;
my $MAXIMUM_ITERATIONS_ALLOWED = 10000;

# Global variables
my (@start_time,@end_time);
my ($number_of_populations, $number_of_different_alleles);
my ($population_name, $population_size, $sample_size);
my (@file_lines, @data, @local_data);
my (@allele_names, @local_names, @local_frequencies);
my (@boundary_array, @count_array);

# Set flag options
my ($all_generations, $vertical_output_on, $horizontal_output_on, $timefile_output_on);
my ($user_defined_seed, $generations_on, $iterations_on, $sampling_only_on, $verbose_on);

GetOptions(
	"a" => \$all_generations,				# If boolean true, list all generations
	"v" => \$vertical_output_on,		# If boolean true, vertical output on
	"h" => \$horizontal_output_on,	# If boolean true, horizontal output on
	"t" => \$timefile_output_on,		# If boolean true, timefile output on	
	"s" => \$sampling_only_on,			# If boolean true, sampling output only on
	"b" => \$verbose_on,						# If boolean true, verbose comments on
	"g=s" => \$generations_on,			# If defined, holds number of generations
	"i=s" => \$iterations_on,				# If defined, holds number of iterations
	"d=s" => \$user_defined_seed		# If defined, holds user defined random number seed
	); 

# Global usage information
my $usage = "Error: Correct usage is: drifter <input_file> <-flags>
                         -b   verbose on-screen information
                         -h   generate horizontal format file output
                         -v   generate vertical format file output
                         -a   output generation values (one population only)
                         -t   generate timefile output
                         -s   enable sampling only

                         -g N integer value for number of generations
                         -i N integer value for number of iterations
                         -d N integer value for random number seed\n\n";

# Explicitly calling random number seed
if($user_defined_seed) {
	srand($user_defined_seed);
}
else {
	srand();
}

# Subroutine functions
sub READ_DATA() {

	# Check for drifter file format
	if ($file_lines[0] !~ m/#drifter/i) {
		die "Whoa! Is this really a drifter file?\n";
	}

	# Read in data
	foreach my $lines (@file_lines) {
		
		if ($lines =~ m/poplabel = \[(\w+)\]/i) {
			push(@data,$1);
			$number_of_populations++;
		}
		elsif ($lines =~ m/Np = \[(\d+)\]/i) { push(@data,$1); }
		elsif ($lines =~ m/Ns = \[(\d+)\]/i) { push(@data,$1); }
		elsif ($lines =~ m/\[(\w+)\]\[(\d+\.\d+)\]/i) {
			push(@data,$1);
			push(@data,$2);
			push(@allele_names,$1);
		}
		elsif ($lines =~ m/;/i) { push(@data,"END"); }
	}
	
	# Get sorted list of unique allele names
	my %seen = ();
	my @unique = grep { ! $seen{$_} ++ } @allele_names;
	@allele_names = sort(@unique);
	$number_of_different_alleles = scalar(@allele_names);
	
	# Cleanup
	undef(@file_lines);
}

sub READ_LOCAL_DATA() {
	
	# Collect population data from @local_data array
	$population_name = shift(@local_data);
	$population_size = shift(@local_data);
	$sample_size = shift(@local_data);
	
	# Shift data onto @local_names and @local_frequencies arrays
	@local_names = ();
	@local_frequencies = ();
	
	until ($local_data[0] =~ m/END/i) {
		my $entry_one = shift(@local_data);
		my $entry_two = shift(@local_data);
		push(@local_names,$entry_one);
		push(@local_frequencies,$entry_two);
		}
	shift(@local_data);	
}

sub CREATE_BOUNDARY_ARRAY() {

	# Create the boundary array for the first iteration
	@boundary_array = ();
	my $local_number_of_alleles = scalar(@local_names);
	my $last_frequency_value = 0;
	for (my $a = 0; $a < $local_number_of_alleles; $a++)
		{
		my $summation = $local_frequencies[$a] + $last_frequency_value;
		push(@boundary_array,$summation);	
		$last_frequency_value = $summation;
		}
}

sub RANDOMLY_POPULATE_FREQUENCIES() {
	
	# Randomly generate frequencies for the t+1 generation
	@count_array = ();
	my $boundary_array_size = scalar(@boundary_array);
	
	# Set all @count_array entries to zero
	for (my $n = 0; $n < $boundary_array_size; $n++)
		{ $count_array[$n] = 0; }
	
	# Populate @count_array entries
	for (my $z = 0; $z < $population_size; $z++)
		{
		# Add 0.0000001 because rand(x) never returns x
		my $iterated_freq = rand($boundary_array[-1] + 0.0000001);
		
		COUNTARRAY:for (my $y = 0; $y < $boundary_array_size; $y++) {
		
			while ($iterated_freq < $boundary_array[$y])
				{
				$count_array[$y]++;
				last COUNTARRAY;
				}
			}	
		}
}

sub CALCULATE_FREQUENCIES() {

	@local_frequencies = ();
	my $count_total = 0;
	foreach my $array_count (@count_array) {
		$count_total += $array_count;
	}

	foreach my $increment_total (@count_array) 
		{
		my $temporary_frequency = $increment_total / $count_total;
		push(@local_frequencies,$temporary_frequency);
		
		# Output frequencies for each generation
		if ($all_generations) {
			my $g_output = sprintf("%.5f", $temporary_frequency);		
			print AOUTPUT "$g_output\t";
			}
		}
}

# Generate a random permutation of @population in place
sub FISHER_YATES_SHUFFLE() {

	my $array = shift;
	
	for (my $i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}

sub GENERATE_SAMPLE_FREQUENCIES() {
	
	# Set up subroutine variables
	my @population = ();
	my $individuals;
	my @population_count_array;
	
	# Populate array
	my $r = 0;
	while ($local_names[$r]) {
		$individuals = int(( $population_size * $local_frequencies[$r] ) + 0.5 );
			for (my $q = 0; $q < $individuals; $q++)
				{ push(@population,$local_names[$r]); }
		$r++;
	}
	
	# Shuffle the array
	&FISHER_YATES_SHUFFLE( \@population );
	
	# Determine number of alleles
	my $local_number_of_alleles = scalar(@local_names);
	
	# Check that $sample_size is not larger than the @population array
	my $pop_array_size = scalar(@population);
	my $local_sample_size = $sample_size;
	if ($local_sample_size > $pop_array_size)
		{ $local_sample_size = $pop_array_size; }

	# Establish zero entries is @population_count_array
	for (my $array_setup = 0; $array_setup < $local_number_of_alleles; $array_setup++ ) {
		$population_count_array[$array_setup] = 0;
		}
		
	# Pick the sample number from the array and generate frequencies 	
	for (my $sampling_counter = 0; $sampling_counter < $local_sample_size; $sampling_counter++) {
		
		for (my $s = 0; $s < $local_number_of_alleles; $s++) {
			
			if ($population[$sampling_counter] =~ $local_names[$s]) {
				$population_count_array[$s]++;
			}
		}
	}

	# Calculate frequencies from the @population_count_array
	@local_frequencies = ();
	for (my $z = 0; $z < $local_number_of_alleles; $z++) {
		my $this_allele_frequency = $population_count_array[$z] / $sample_size;
		push(@local_frequencies,$this_allele_frequency);
	}

}

sub PRINT_FREQUENCIES() 	{
	
	# Setting up population names in output files
	if ($vertical_output_on) {
		print VOUTPUT "$population_name\t";
	}
	if ($horizontal_output_on) {
		print HOUTPUT "$population_name\t";
	}
	if ($all_generations) {
		print AOUTPUT "$population_name\t";
	}

	my $array_size = scalar(@local_names);
	
	foreach my $name (@allele_names) {
	
		my $match_found = 0;
		for (my $e = 0; $e < $array_size; $e++) {
		if ($name =~ $local_names[$e]) {
			my $allele_frequency = sprintf("%.5f", $local_frequencies[$e]);
			if ($vertical_output_on) {
				print VOUTPUT "$allele_frequency\t";
			}
			if ($horizontal_output_on) {
				print HOUTPUT "$allele_frequency\t";
			}
			if ($all_generations) {
				print AOUTPUT "$allele_frequency\t";
			}
			$match_found = 1;
			}
		}	
		
		if ($match_found == 0) {
			if ($vertical_output_on) {
				print VOUTPUT "0\t";
			}
			if ($horizontal_output_on) {
				print HOUTPUT "0\t";
			}
			if ($all_generations) {
				print AOUTPUT "0\t";
			}
		}
	}

	if ($vertical_output_on) {
		print VOUTPUT "\n";
	}
	
}

# Open input nexus file
my $input_file = shift;
if (!defined $input_file || ! -e $input_file || -d $input_file) {
	die $usage;
}
open(INPUT,"<$input_file");

# Generate file names
my $input_file_name;
if($input_file =~ m/(\w+)\.\w+/i) {
	$input_file_name = $1;
}

my $vertical_output_name = $input_file_name . "_vertical.out";
my $horizontal_output_name = $input_file_name . "_horizontal.out";
my $allgen_output_name = $input_file_name . "_allgen.out";
my $timefile_output_name = $input_file_name . "_timefile.out";

# Tests for vertical output file request
if($vertical_output_on || (!$vertical_output_on && !$horizontal_output_on)) {
	if (-e $vertical_output_name) {
		die "\nDang it! I can't make $vertical_output_name\.\nDon't use filenames more than once...!\n";
	}
	open(VOUTPUT,">$vertical_output_name") or die "Aaaggh! Can't create $vertical_output_name\....\n";
	$vertical_output_on = "true";
}

# Tests for horizontal output file request
if($horizontal_output_on) {
	if (-e $horizontal_output_name) {
		die "\nDang it! I can't make $horizontal_output_name\.\nDon't use filenames more than once...!\n\n";
	}
	open(HOUTPUT,">$horizontal_output_name") or die "Aaaggh! Can't create $horizontal_output_name\....\n";
}

# Tests for all_generations output file request
if ( $all_generations ) {
	if (-e $allgen_output_name) {
		die "\nDang it! I can't make $allgen_output_name\.\nDon't use filenames more than once...!\n\n";
	}
	open(AOUTPUT,">$allgen_output_name") or die "Aaaggh! Can't create $horizontal_output_name\....\n";
}

# Tests for timefile output file request
if($timefile_output_on) {
	if (-e $timefile_output_name) {
		die "\nDang it! I can't make $timefile_output_name\.\n";
	}
	open(TOUTPUT,">>$timefile_output_name") or die "Aaaggh! Can't create $timefile_output_name\....\n";
}

#Beginning interaction with user
if ($verbose_on) {
	print STDOUT "\nWelcome to drifter v.2.2, dude....\n";
}

#Reading in input file and closing
@file_lines = <INPUT>;
foreach my $file_string (@file_lines) {
	chomp($file_string);
}
close INPUT or die "\n Having trouble with the gate: can't close $input_file\n\n";

#Requesting number of generations
if ($generations_on) {
	$NUMBER_OF_GENERATIONS = $generations_on;
}
elsif ($sampling_only_on) {
	$NUMBER_OF_GENERATIONS = 0;
}
else {
	print STDOUT "And for how many generations? [Default: 0]  ";
	chomp($NUMBER_OF_GENERATIONS = <STDIN>);
	if (!$NUMBER_OF_GENERATIONS) {
		$NUMBER_OF_GENERATIONS = 0;
	}
}

#Requesting number of iterations
if ($iterations_on) {
	$NUMBER_OF_ITERATIONS = $iterations_on;
}
else {
	print STDOUT "And iterate for how many times? [Default: 0]  ";
	chomp($NUMBER_OF_ITERATIONS = <STDIN>);
	if (!$NUMBER_OF_ITERATIONS) {
		$NUMBER_OF_ITERATIONS = 0;
	}
	if ($NUMBER_OF_ITERATIONS > $MAXIMUM_ITERATIONS_ALLOWED) {
		die "You want more than $MAXIMUM_ITERATIONS_ALLOWED iterations? Are you insane?!\n";
	}
}

# Read in data file
READ_DATA();

# Die if all_generations option enabled and more than one population
if ( $all_generations ) {
	if ( $number_of_populations > 1) {
		print STDOUT "\nThe -a flag is illegal with more than one population\n";
		exit();
	}
}

if ($verbose_on) {	
	print STDOUT "\nLook's like we've got $number_of_populations populations and $number_of_different_alleles different alleles, bud.\n";
	print STDOUT "And we're drifting for $NUMBER_OF_GENERATIONS generations, right?\n\n";
	print STDOUT "So, let's get down to those $NUMBER_OF_ITERATIONS iterations then....\n\n";
}

# Setting up VOUTPUT file
if ($vertical_output_on) {
	print VOUTPUT "Number of iterations: $NUMBER_OF_ITERATIONS\n";
	print VOUTPUT "Number of generations: $NUMBER_OF_GENERATIONS\n\n";
	foreach my $a_names (@allele_names) {
		print VOUTPUT "\t$a_names";
	}
	print VOUTPUT "\n";
}

# Setting up HOUTPUT file
if ($horizontal_output_on) {
	for ( my $h = 0; $h < $number_of_populations; $h++ )
		{
		foreach my $name_of_allele (@allele_names) {
			print HOUTPUT "\t$name_of_allele";
		}
		print HOUTPUT "\t";
		}
	print HOUTPUT "\n";
}

# Setting up AOUTPUT file
if ($all_generations) {
	for ( my $h = 0; $h < $number_of_populations; $h++ )
		{
		foreach my $name_of_allele (@allele_names) {
			print AOUTPUT "\t$name_of_allele";
		}
		print AOUTPUT "\t";
		}
	print AOUTPUT "\n";
}

# Print beginning time point
if ($timefile_output_on) {
	$Time::tm = localtime;
	(my $MDAY, my $HOUR, my $MINUTE, my $SECOND)
		= ($Time::tm->mday, $Time::tm->hour, $Time::tm->min, $Time::tm->sec);
	print TOUTPUT "Operations started -- day $MDAY: hr $HOUR: min $MINUTE: sec $SECOND\n";
	@start_time = ($MDAY, $HOUR, $MINUTE, $SECOND);
}

#Begin requested iterations
for(my $iterations_counter = 0; $iterations_counter < $NUMBER_OF_ITERATIONS; $iterations_counter++) {
		
	#creates a local copy of @data array
	@local_data = @data;

	#Begin iterations of populations
	for(my $population_counter = 0; $population_counter < $number_of_populations; $population_counter++) {
		
		&READ_LOCAL_DATA();
		&CREATE_BOUNDARY_ARRAY();
		
		#Begin number of generations iterations
		for(my $generation_counter = 0; $generation_counter < $NUMBER_OF_GENERATIONS; $generation_counter++) {

				# Generation counter incremented for all_generations file
				my $local_g_counter = $generation_counter + 1;
				if ($all_generations) {
					print AOUTPUT "G$local_g_counter\t";
				}
			
				# Creates individual data for a single population
				&RANDOMLY_POPULATE_FREQUENCIES();
							
				# Calculates frequency in t+1 generation
				&CALCULATE_FREQUENCIES();
			
				# Generates subsequent boundary arrays
				&CREATE_BOUNDARY_ARRAY();
			
        } #End of generations iterations loop
        	
      #Creates data for final population
			&GENERATE_SAMPLE_FREQUENCIES();
			&PRINT_FREQUENCIES();
     	
         	if ($verbose_on) {
         		print STDOUT ".";
         	}
         	
		} #End of populations iterations loop

		if ($vertical_output_on) {
			print VOUTPUT ";\n";
		}
		if ($horizontal_output_on) {
			print HOUTPUT "\n";
		}
		if ($all_generations) {
			print AOUTPUT "\n";
		}
		if ($verbose_on) {
			print STDOUT "|";
		}
		
	} #End of iterations loop

# Print end time point
if ($timefile_output_on) {
	$Time::tm = localtime;
	(my $MDAY, my $HOUR, my $MINUTE, my $SECOND)
		= ($Time::tm->mday, $Time::tm->hour, $Time::tm->min, $Time::tm->sec);
	print TOUTPUT "Operations finished -- day $MDAY: hr $HOUR: min $MINUTE: sec $SECOND\n";
	@end_time = ($MDAY, $HOUR, $MINUTE, $SECOND);
	my $number_of_seconds = ($end_time[3] - $start_time[3]);
	$number_of_seconds += 60 * ($end_time[2] - $start_time[2]);
	$number_of_seconds += 60 * 60 * ($end_time[1] - $start_time[1]);
	$number_of_seconds += 60 * 60 * 60 * ($end_time[0] - $start_time[0]);
	print TOUTPUT "Operations took $number_of_seconds seconds\n";
}

# Close filehandles
if ($vertical_output_on) {
	close VOUTPUT or die "\nDang! Can't close $vertical_output_name\n\n";
}
if ($horizontal_output_on) {
	close HOUTPUT or die "\nCan't close $horizontal_output_name\n\n";
}
if ($all_generations) {
	close AOUTPUT or die "\nCan't close $allgen_output_name\n\n";
}

# Confirm correct programme termination
if ($verbose_on) {
	print STDOUT "\n\nThat's the end of this here rodeo, folks. Thank you all for coming!\n";
}
