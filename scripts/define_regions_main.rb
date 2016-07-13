#!/usr/bin/env ruby
require 'optparse'
require_relative 'define_regions_functions'
require 'fileutils'

$round_parm = 3
$prec_parm = 0.001
$cM = 0.1
$bp = 0

input = ""
output = Dir.pwd
sInputPvals = ""

# Parser for the arguments
# -i input
# -m cM (normally 0.1(default) or 0.2)
# -r bp region (normally 50kb)

OptionParser.new do |opts|
  opts.banner = "Usage: your_app [options]"
  opts.on('-i [ARG]', '--input_i [ARG]', "Specify the input_i") do |v|
    input = v
  end

  opts.on('-o [ARG]', '--output_o [ARG]', "Specify the output_o") do |v|
    output = v
  end

  opts.on('-m [ARG]', '--cmregion_n [ARG]', "Specify the cM region m") do |v|
    $cM = v.to_f
  end

  opts.on('-r [ARG]', '--bpregion_r [ARG]', "Specify the bp region r") do |v|
    $bp = v.to_i
  end
    
  opts.on('--version', 'Display the version') do 
    puts "VERSION"
    exit
  end
  opts.on('-h', '--help', 'Display this help') do 
    puts opts
    exit
  end
end.parse!

if (input == "")  #or keys_snp.size! 
  puts 'Error please check arguemnts of the input files'
  abort
end

if $bp!=0 then $cM=0 end
puts "input: " + input
puts "cM: " + $cM.to_s
puts "bp region: " + $bp.to_s

# read in the indexSNPs
cont_snps = read_file(input,true)
if cont_snps.nil?
    abort
end

a = Time.now

if $bp==0
  # calculates the extended region for every indexSNP
  data,snps_not_in_range, num_cols = calculation(cont_snps)
  output_data = [data,snps_not_in_range]
  #puts output_data[1]
  puts "SNPs not in range: " + snps_not_in_range.size.to_s
  cm_output = true
else
  output_data = cont_snps
  num_cols = cont_snps[1].split("\t").size
  cm_output = false
end
b = Time.now
puts "Compilation time: " + ((b-a)/60).to_s
puts "Number of columns in input file :" + num_cols.to_s
write_output(output_data,cm_output,output)

#if $bp==0
#  path_bound_red = "output/regions/region_boundaries_"+ $cM.to_s.gsub(".","") + "cm.txt"
#else
#  path_bound_red = "output/regions/region_boundaries_"+ $bp.to_s + "bp.txt"
#end

if $bp==0
  path_bound_red = output + $cM.to_s.gsub(".","") + "cm.txt"
else
  path_bound_red = output + $bp.to_s + "bp.txt"
end
