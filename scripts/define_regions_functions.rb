#!/usr/bin/ruby
def read_file(arg,skip)
    begin
	if skip then i =1 else i=0 end
        content = File.readlines(arg)
        return content[i,content.size] #need task
    rescue
        puts 'Cannot open ' + arg
    end       
end

def calc_cmregions(rounding,bp,cMcount,num_bp,inter,calc_pos,cMbp)   
    if (rounding > $prec_parm )
       # puts "case 1"
        bp = bp + num_bp
        cMcount = cMcount + inter
        snp_pos = calc_pos
        return false, bp, cMcount, snp_pos

    #if reached, start searching for 0.2cM
    elsif (rounding == $prec_parm or rounding == 0)
      #  puts "case 2"
        bp = bp + num_bp
        cMcount = cMcount + inter
        return true, bp, cMcount
       
    #if overreached calculate using rest/cMbp
    else
       # puts "case 3"
        num_bp = (($cM-cMcount)/cMbp).round              
        cMcount = cMcount + num_bp*cMbp 
        bp = bp + num_bp
        return true, bp, cMcount       
     end    
end

def bimodal_calculation(i1,snp_pos1,cont_chr)
    regFWD = []
    regBWD = []
    snp_pos = snp_pos1
    cMcountFwd =0
    bpfwd = 0

    #forward
    # for each position
    (i1..cont_chr.size-1).each_with_index do |x,e|
        line_chr = cont_chr[x]
        keys_chr = line_chr.split(" ")
        if cont_chr.size-1 == e+i1
          regFWD = true, bpfwd, cMcountFwd
          break
        end

        ending = cont_chr[i1+e+1].split(" ")[1].to_i-1
        cMbp = keys_chr[2].to_f/(10**6)
        num_bp = ending - snp_pos      
        inter = num_bp*cMbp    
        rounding = (($cM-(inter+cMcountFwd))*100).round/100.0
        status, bpfwd, cMcountFwd, *rest = calc_cmregions(rounding,bpfwd,cMcountFwd,num_bp,inter,ending,cMbp)
        if status
            cMcountFwd = (cMcountFwd*10000).round/10000.0
            regFWD = true, bpfwd, cMcountFwd    
            break  
        else
            snp_pos  = rest[0]
        end
        if x == cont_chr.size-1
            cMcountFwd = (cMcountFwd*10000).round/10000.0
            regFWD =  false, bpfwd, cMcountFwd
            break
        end        
    end
    #backward
    snp_pos = snp_pos1
    y = 0
    cMcountBwd = 0
    bpbwd = 0
     i1.downto(0) do |x|
        line_chr = cont_chr[x]
        keys_chr = line_chr.split(" ")
        start = keys_chr[1].to_i
        cMbp = keys_chr[2].to_f/(10**6)
        num_bp = snp_pos - start      
        inter = num_bp*cMbp    #cM for given num_bp
        rounding = (($cM-(inter+cMcountBwd))*100).round/100.0
        status, bpbwd, cMcountBwd, *rest = calc_cmregions(rounding,bpbwd,cMcountBwd,num_bp,inter,start,cMbp)
        if status
            cMcountBwd = (cMcountBwd*10000).round/10000.0
            regBWD = true, bpbwd, cMcountBwd   
            break 
        else
            snp_pos  = rest[0]
            snp_pos = snp_pos -1
        end
        if x == 0
            cMcountBwd = (cMcountBwd*10000).round/10000.0
            regBWD = false, bpbwd, cMcountBwd
            break
        end
    end  
   return regFWD, regBWD 
end 

def calculation(cont_snps)
  snps_total = cont_snps.size
  puts "number of snps: " + snps_total.to_s
  puts "estimated time: " + ((snps_total/6)).to_s + "sec"
  cur_chr = -1
  cont_chr = ""
  #rates_path = "/Users/katja/Desktop/Manchester/Detection/Bayesian/recombination_rates/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr"
  # HapMap recombination rates
  rates_path = "source_data/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr"
  output_data = []
  not_in_range_snps = []
  num_cols = 0
  
  # go over all indexSNPS
  cont_snps.each_with_index do |line_snp,y|
      keys_snp = line_snp.split("\t") 
      num_cols = keys_snp.size
      # file should contain ID, CHROM, START-POS (optional GENE)
      if (num_cols != 3) and (num_cols != 4)  #or keys_snp.size!
        keys_snp = line_snp.split(" ")
        num_cols = keys_snp.size

        if (num_cols != 3) and (num_cols != 4)  #or keys_snp.size! 
          puts 'Check if file with SNPs is tab or space-delimited and contains 3 or 4 columns'
          puts 'Number of columns: ' + num_cols.to_s
          abort
        end
      end

      # progress counter
      per20 = (snps_total*2/10).round
      per50 = (snps_total*5/10).round
      per80 = (snps_total*8/10).round
      per10 = (snps_total*1/10).round
      if per50 == y
          puts "50% of SNPs are worked out"
      elsif per20 == y 
          puts "20% of SNPs are worked out"
      elsif per80 ==y
         puts "80% of SNPs are worked out"
      elsif per10 ==y
         puts "10% of SNPs are worked out"            
      end

      # read in the rates for each individual chromosome
      snp_name = keys_snp[0]
      snp_pos = keys_snp[2].to_i      
      if keys_snp[1]!= cur_chr
          cur_chr = keys_snp[1]
          chr_file = rates_path+cur_chr+".txt"
          cont_chr = read_file(chr_file,true)
      end
      inrange = false
      #can be optmized as SNPs are sorted accodring to chr and pos
      size_chr = cont_chr.size-1  # remove header

      # go through all regions of the recombination rate file associated with the specific chromosome
      # columns are Chromosome  Position(bp)  Rate(cM/Mb) Map(cM)
      cont_chr.each_with_index do |line_chr,i|
          keys_chr = line_chr.split(" ")        
          start = keys_chr[1].to_i
          if size_chr == i
            break
          end
          ending = (cont_chr[i+1].split(" "))[1].to_i-1
          # look if indexSNP is withing regions of the recombination file 
          if snp_pos.between?(start,ending)
              inrange = true
              regFWD, regBWD = bimodal_calculation(i,snp_pos,cont_chr)  

              output_data.push ([snp_name, keys_snp[1], snp_pos, regFWD, regBWD ]) 
              break
          else
              next
          end
      end
      if not inrange        
          not_in_range_snps.push([snp_name,keys_snp[1], snp_pos]) 
      end 
  end
 
  return output_data,not_in_range_snps,num_cols
end

def write_output(output_data,cmoutput)
	#output cM regions
	t = "\t"
	if cmoutput 
	#	puts(output_data)
		out, snps_not_range = output_data
		cM_m = $cM.to_s.gsub(".", "")
		# path = "outputs/regions_"+ $cM.to_s + "cm.txt"
		# path_bound = "outputs/boundaries_"+ $cM.to_s + "cm.txt"
		# path_bound_red = "outputs/red_boundaries_"+ $cM.to_s + "cm.txt"
		path = "output/regions/supplementary/regions_"+ cM_m + "cm.txt"
		path_bound = "output/regions/supplementary/boundaries_"+ cM_m + "cm.txt"
		path_bound_red = "output/regions/red_boundaries_"+ cM_m + "cm.txt"
		dirs = [path,path_bound,path_bound_red]
		dirs.each do |dir|
			dirname = File.dirname(dir)
			unless File.directory?(dirname)
				FileUtils.mkdir_p(dirname)
			end
		end	

		o_file = File.open(path,"w")
		b_file = File.open(path_bound,"w")
		rb_file = File.open(path_bound_red, "w")
		
		header = "SNP" +t + "chr" +  t + "Pos" + t + "reached" +t+ "FWD_region(0.1cM)" + t +  "FWD0.1cM" + t + "reached" + t + "BWD_region(0.1cM)" + t + "BWD0.1cM"  + "\n"
		o_file.write(header)
		header = "SNP" +t + "chr" +  t + "Pos"  + t +  "up_bound" + t + "low_bound" + "\n"
		b_file.write(header)
		out.each do |d|
			#puts(d)
			d.flatten!
			snp, chr, pos, reached, fwd_reg01,  fwd01,  reached2, bwd_reg01,  bwd01 = d
			# puts(d)			
			o_file.write(snp + t + chr + t + pos.to_s + t + reached.to_s + t + fwd_reg01.to_s + t + fwd01.to_s + t + reached2.to_s + t + bwd_reg01.to_s + t + bwd01.to_s  + "\n")
			upper = pos + fwd_reg01
			lower = pos - bwd_reg01
			if lower < 0 then lower = 0 end
			b_file.write(snp + t + chr + t + pos.to_s + t + upper.to_s + t + lower.to_s + "\n" )
			rb_file.write(snp + t +"chr"+chr+":"+lower.to_s+"-"+upper.to_s + "\n")
		end

		o_file.close()

		if !snps_not_range.empty?
			path = "output/regions/supplementary/snps_not_in_range.txt"
			#path = "outputs/snps_not_in_range.txt"
			dirname = File.dirname(path)
			unless File.directory?(dirname)
			FileUtils.mkdir_p(dirname)
			end

			o_file = File.open(path,"w")
			snps_not_range.each do |d|
				snp, chr, pos, p_value, maf = d
				o_file.write(snp + t + chr + t + pos.to_s + "\n")
				lower = pos - 50000
				if lower <0 then lower =0 end
				b_file.write(snp + t + chr + t + pos.to_s + t + (pos+50000).to_s + t +lower.to_s+ "\n")
				rb_file.write(snp + t + "chr" +chr+ t + lower.to_s+"-"+(pos+50000).to_s+ "\n")
			end
			o_file.close()
		end
		b_file.close()
		rb_file.close()

	elsif !cmoutput 
		puts("BP")
		path_bound_red = "output/regions/red_boundaries_"+ $bp.to_s + "bp.txt"
		dirname = File.dirname(path_bound_red)
		unless File.directory?(dirname)
			FileUtils.mkdir_p(dirname)
		end

		rb_file = File.open(path_bound_red, "w")

		output_data.each_with_index do |line_snp,y|
     		keys_snp = line_snp.split("\t")
     		snp = keys_snp[0] 
     		chr = keys_snp[1]
     		pos = keys_snp[2].to_i
     		lower = pos - $bp
		if lower < 0  then lower = 0 end
     		upper = pos + $bp
     		if keys_snp.size != 3 #or keys_snp.size!
         	 	puts 'Check if file with SNPs is tab-delimited and contains 3 columns'
         i		abort
      		end
      		rb_file.write(snp + t +"chr"+chr+":"+lower.to_s+"-"+upper.to_s + "\n")
      	end	
        rb_file.close()
	end
end

def retr_region_snps(path,num_cols,sInputPvals)
    file = read_file(path,false)

    file_pvalues = read_file(sInputPvals,false)

    file.each_with_index do |line,i|
	puts i
        keys = line.split("\t")
        indexSNP = keys[0]
	puts indexSNP
        reg = keys[1].strip
        patt_reg_snps = 'output/regions/pri_reg_' + i.to_s + '.txt'
        dirname = File.dirname(patt_reg_snps)
    	unless File.directory?(dirname)
    	FileUtils.mkdir_p(dirname)
    	end

        writer = open(patt_reg_snps, 'w') 

        if i == 0
            writer.write("#SNPID \t CHROM \t POS \t A1 \t A2 \t A1_UNAFF \t PVAL \t indexSNP \t GENE")
        else 
            # retrieve chromosome, start and end position of the indexSNP
            # example: reg = chr1:113843881-114538820
            parts = reg.split("\:")
            chromosome = parts[0]
            start = parts[1].split("\-")[0]
            ende = parts[1].split("\-")[1]
	   
	    puts indexSNP 
	    puts chromosome
	    puts start
	    puts ende 

	
            # go through all lines of the file that contains the pvalue and MAF of all SNPS
            file_pvalues.each_with_index do |line_pvlaue,j|
                tokens = line_pvlaue.split("\t")
                num_cols = tokens.size

                if (num_cols < 7) 
                    tokens = line_pvlaue.split(" ")
                    num_cols = tokens.size

                    if (num_cols < 7) 
                        tokens = line_pvlaue.split("\,")
                        num_cols = tokens.size

                        if (num_cols < 7) 
                            puts 'Check if file with pvalues for SNPs is tab, comma or space-delimited and contains 7 columns'
                            puts 'Number of columns: ' + num_cols.to_s
                            abort
                        end
                    end
                end

                # look if position of SNP is between start and end and the chromosome is the same
                # if so write line to file
                if ( ("chr" + tokens[1]) == chromosome && tokens[2] >= start && tokens[2] <= ende)

                    writer.write("\n" + tokens[0] + "\t" + chromosome + "\t" + tokens[2] + 
                        "\t" + tokens[3] + "\t" + tokens[4] + "\t" + tokens[5] + 
                        "\t" + tokens[6] + "\t"  + indexSNP + "\t" + "--")

                    # use this for JIA data
                    #writer.write("\n" + tokens[0] + "\t" + chromosome + "\t" + tokens[2] + 
                    #    "\t" + tokens[3] + "\t" + tokens[4] + "\t" + tokens[5] + 
                    #    "\t" + tokens[8] + "\t"  + indexSNP + "\t" + "--")

                end
            end
       end

       writer.close

    end


    # concatenating all output files (regarding SNPs)
    befehl = 'rm -f output/regions/regional_snps.txt'
    system befehl
    befehl = 'cat output/regions/pri_reg*.txt > output/regions/regional_snps.txt'
    system befehl
    #befehl = 'rm -f output/regions/pri_reg*.txt'
    #system befehl

end
