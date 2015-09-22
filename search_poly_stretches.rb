#!/mnt/prostlocal/programs/ruby/2.2.2/bin/ruby


## TODO: we search poly-R or Y stretches with
## - minimum length (bspw 12)
## - einer maximalen rate von R respektive Y Interruptions (bspw 10%)

## TODO add param for 'at least x % of G in seed'
## reverse complement search?

require 'nyaplot'
require 'bio'

class SearchPolyStretches

  def initialize(fastq, min_stretch_length, mismatches, min_guanine_cytosine_rate, just_plot)

    $out_dir = "../poly_out/#{File.basename(fastq).split('.')[0]}/"
    Dir.mkdir($out_dir) unless Dir.exist?($out_dir)

    samples_h = {File.basename(fastq).split('.')[0] => 'TY1-03112014`'}

    if just_plot

      plot_overview

      ####################
      ## sortmerna
      ####################
      sortmerna("../sortmerna", samples_h)

    else

      params = "min-#{min_stretch_length}_mis-#{mismatches}_gc-#{min_guanine_cytosine_rate}"
      $out_dir = "../poly_out/#{File.basename(fastq).split('.')[0]}/#{params}/"
      Dir.mkdir($out_dir) unless Dir.exist?($out_dir)

      $min_stretch_length = min_stretch_length

      puts "\n\n#######################\n Input:\t #{fastq}\n\n"
      $fastq = File.open(fastq, 'r')

      # build all regular expressions based on given params
      $regex_purine_a = []
      $regex_pyrimidine_a = []
      build_regex(min_stretch_length, mismatches)

      # run through the read file and check for each read if at least one regexp matches
      line_count = 0
      read_id = nil
      $candidates = {:purin => {}, :pyrimidin => {}}
      $read_count = 0
      $fastq.each do |fastq_line|
        line_count += 1
        read_id = fastq_line if line_count == 1
        if line_count == 2
          $read_count += 1
          read_seq = fastq_line.chomp.upcase
          read_fasta = Bio::FastaFormat.new(">#{read_id}\n#{read_seq}")
          get_candidate(read_fasta, $regex_purine_a, %w(A G), min_guanine_cytosine_rate)
          get_candidate(read_fasta, $regex_pyrimidine_a, %w(T C), min_guanine_cytosine_rate)
          puts "processed #{$read_count} reads..." if $read_count.modulo(100000) == 0
        end
        line_count = 0 if line_count == 4
      end

      puts "Found #{$candidates[:purin].size} poly-purine and #{$candidates[:pyrimidin].size} poly-pyrimidine candidate stretches with a minimum length of #{min_stretch_length} nt and #{mismatches} allowed mismatches with a G/C content of at least #{min_guanine_cytosine_rate*100}% in #{$read_count} inspected reads."

      $length_distribution = { :purin => {}, :pyrimidin => {} }
      $length_distribution.each do |type, hash|
        50.times do |i|
          hash[i] = 0 unless i < min_stretch_length
        end
      end

      # write csv
      write_csv($candidates)
      # write fasta & fill length_counter
      write_fasta($candidates)

      # count based on poly-stretch length and plot
      #plot_length_distribution

      # write log file
      log = File.open("#{$out_dir}/log.txt",'w')
      log << "##{fastq}\n\n##{params}\n\n"
      log << "input_reads:\t#{$read_count}\n"
      log << "candidate_purin_reads:\t#{$candidates[:purin].size}\n"
      log << "candidate_pyrimidin_reads:\t#{$candidates[:pyrimidin].size}\n"
      log << "purin_rate:\t#{($candidates[:purin].size.to_f / $read_count * 100).round(2)}%\n"
      log << "pyrimidin_rate:\t#{($candidates[:pyrimidin].size.to_f / $read_count * 100).round(2)}%\n"
      log.close

    end

  end

  #################
  ## 1) for each sample build a pie chart with the amount of 5s, 5.8s, 18s and 28s rRNA
  ## 2) for all six samples build a bar plot with the amount of rRNA and 'other' reads
  def sortmerna(input_dir, samples_h)

    samples, labels = [], []
    rrna_5s_a, rrna_5_8s_a, rrna_18s_a, rrna_28s_a, rrna_all = [], [], [], [], []

    if File.exist?("#{input_dir}/#{samples_h.keys[0]}_aligned.log")
      puts "\t\tRun SortMeRna Statistics..."

      samples_h.each do |sample, label|
      samples.push(sample)
      labels.push(label)
      read = false
      all = 0
        File.open("#{input_dir}/#{sample}_aligned.log",'r').each do |line|
          read = true if line.include?('By database:')
          if read
            if line.include?('rfam-5s-database-id98.fasta')
              value = line.split("\t")[2].chomp.sub('%','').to_f
              rrna_5s_a.push(value)
              all += value
            end
            if line.include?('rfam-5.8s-database-id98.fasta')
              value = line.split("\t")[2].chomp.sub('%','').to_f
              rrna_5_8s_a.push(value)
              all += value
            end
            if line.include?('silva-euk-18s-id95.fasta')
              value = line.split("\t")[2].chomp.sub('%','').to_f
              rrna_18s_a.push(value)
              all += value
            end
            if line.include?('silva-euk-28s-id98.fasta')
              value = line.split("\t")[2].chomp.sub('%','').to_f
              rrna_28s_a.push(value)
              all += value
            end
        end
      end
      rrna_all.push(all.round(2))
    end
    df = Nyaplot::DataFrame.new({:sample => samples, :label => labels, :rrna_5s => rrna_5s_a, :rrna_5_8s => rrna_5_8s_a, :rrna_18s => rrna_18s_a, :rrna_28s => rrna_28s_a, :all_rrna => rrna_all })
      df = Nyaplot::DataFrame.new({:label => %w(5s 5.8s 18s 28s all), :rrna => rrna_5s_a+rrna_5_8s_a+rrna_18s_a+rrna_28s_a+rrna_all })

    colors = Nyaplot::Colors.qual
    frame = Nyaplot::Frame.new

    #[:all_rrna, :rrna_5s, :rrna_5_8s, :rrna_18s, :rrna_28s].each do |rrna|
      plot = Nyaplot::Plot.new
      plot.configure do
        x_label('rRNA type')
        y_label("% rRNA")
        yrange([0,100])
        legend(true)
      end
      #bar = plot.add_with_df(df, :bar, :label, rrna) # x-> column :label, y-> column :rrna
      bar = plot.add_with_df(df, :bar, :label, :rrna)
      bar.color(colors)
      frame.add(plot)

    frame.export_html("#{$out_dir}/sortmerna.html")
    frame_html = File.open("#{$out_dir}/sortmerna.html",'a')
    frame_html << "\n\n<p>\n#{df_to_html_table(df.to_html)}\n</p>\n\n"
    frame_html.close
    end
  end

  def df_to_html_table(html)
    refactored_html = html
    refactored_html.sub!('<table>','<table><thead>')
    refactored_html.sub!('</tr>','</tr></thead><tbody>')
    refactored_html.sub!('</table>','</tbody></table>')
    refactored_html
  end

  def get_candidate(read, regex_a, type, rate)

    #puts "read: #{read.definition}"

    read_seq = read.seq

    stretch = (0..0)
    regex_a.each do |regexp|
      match = read_seq.match(regexp)
      if match
        #puts "\t#{regexp.to_s}\t#{match[0]}\t#{match.begin(0)}\t#{match.end(0)}"
        match_begin = match.begin(0)
        match_end = match.end(0)
        if stretch.size < (match_end-match_begin)
          stretch = (match_begin..match_end)
        end
      end
    end
    if stretch.size > 1
      extend_stretch_seed(read, stretch, type, rate)
    end
  end

  def extend_stretch_seed(read, stretch, type, rate)
    #puts "\tSEED\t#{read.seq[stretch.begin, (stretch.end-stretch.begin)]}"
    nt_a = read.seq.scan(/./)
    seed_begin = stretch.begin; seed_end = stretch.end

    extension_5_prime = 0
    extension_3_prime = 0

    # extend 5'
    nt_a[0..seed_begin].reverse.each do |nt|
      if type.include?(nt)
        extension_5_prime += 1
      else
        break
      end
    end

    # extend 3'
    nt_a[seed_end..read.length].each do |nt|
      if type.include?(nt)
        extension_3_prime += 1
      else
        break
      end
    end

    #puts "\t5' Extension: #{extension_5_prime}"
    #puts "\t3' Extension: #{extension_3_prime}"

    candidate_seq = read.seq[seed_begin-extension_5_prime+1, seed_end-seed_begin+extension_5_prime+extension_3_prime-1]

    # check rate of guanine/cytosine
    if type.include?('A') && type.include?('G')
      candidate_rate = (candidate_seq.count('G').to_f / candidate_seq.size).round(2)
    else
      candidate_rate = (candidate_seq.count('C').to_f / candidate_seq.size).round(2)
    end

    if candidate_rate >= rate
      if type.include?('A') && type.include?('G')
        $candidates[:purin][read.definition] = candidate_seq
      else
        $candidates[:pyrimidin][read.definition] = candidate_seq
      end
    end

  end

  def write_csv(candidates)
    candidates.each do |type, candidate_h|
      out = File.open("#{$out_dir}/#{type}.csv",'w')
      out << "#Read\tSeq\n"
      candidate_h.each do |read_id, candidate_seq|
        out << "#{read_id}\t#{candidate_seq}\n"
      end
      out.close
    end
  end

  def write_fasta(candidates)
    candidates.each do |type, candidate_h|
      out = File.open("#{$out_dir}/#{type}.fa",'w')
      candidate_h.each do |read_id, candidate_seq|
        length = candidate_seq.length
        if $length_distribution[type][length]
          $length_distribution[type][length] += 1
        else
          $length_distribution[type][length] = 1
        end
        out << ">#{read_id}\n#{candidate_seq}\n"
      end
      out.close
    end
  end

  def build_regex(min_length, mismatches)
    purine = '[AG]'
    pyrimidine = '[TC]'

    [purine, pyrimidine].each do |set|
      start = 1
      (min_length-2).times do
        regexp = ''
        regexp << set << "{#{start}}"
        regexp << ".{#{mismatches}}"
        regexp << set << "{#{min_length-start-1}}"
        start += 1
        $regex_purine_a.push(Regexp.new(regexp)) if set == '[AG]'
        $regex_pyrimidine_a.push(Regexp.new(regexp)) if set == '[TC]'
      end
    end
  end

  def plot_overview

    # collect data
    length_a = []; count_a = []; percent_a = []; type_a = []; mismatches_a = []; gc_a = []; input_reads = []
    Dir.glob("#{$out_dir}/min*").each do |dir|
      params = File.basename(dir).split('_')
      length = params[0].split('-')[1].to_i
      mis = params[1].split('-')[1]
      gc = params[2].split('-')[1]
      if File.exist?("#{dir}/log.txt")
        2.times do
          length_a.push(length)
          mismatches_a.push(mis)
          gc_a.push(gc)
        end
        type_a.push('purin')
        type_a.push('pyrimidin')
        File.open("#{dir}/log.txt",'r').each do |line|
          2.times do input_reads.push(line.split("\t")[1].chomp) end if line.start_with?('input_reads')
          count_a.push(line.split("\t")[1].to_i) if line.start_with?('candidate_purin_reads')
          count_a.push(pyrimidin_reads = line.split("\t")[1].to_i) if line.start_with?('candidate_pyrimidin_reads')
          percent_a.push(line.split("\t")[1].sub('%','').to_f) if line.start_with?('purin_rate')
          percent_a.push(line.split("\t")[1].sub('%','').to_f) if line.start_with?('pyrimidin_rate')
        end
      end
    end

    df = Nyaplot::DataFrame.new({:length => length_a, :count => count_a, :mismatches => mismatches_a, :type => type_a, :percent => percent_a, :gc => gc_a, :input_reads => input_reads})

    plot = Nyaplot::Plot.new
    plot.x_label('Minimum length [nt]')
    plot.y_label('# poly-R/Y [reads]')
    sc = plot.add_with_df(df, :scatter, :length, :count)
    sc.tooltip_contents([:input_reads, :type, :mismatches, :gc, :percent])

    colors = Nyaplot::Colors.qual
    sc.color(colors)
    sc.fill_by(:type)
    sc.shape_by(:mismatches)

    plot.legend(true)

    frame_overview = Nyaplot::Frame.new
    frame_overview.add(plot)
    frame_overview.export_html("#{$out_dir}/overview.html")

    ## add a header
    overivew_file = File.open("#{$out_dir}/overview.html",'r')
    tmp_file = File.open("#{$out_dir}/overview_tmp.html",'w')
    overivew_file.each do |line|
      if line.start_with?('<body>')
        tmp_file << "\n\n<body><br><h1>#{File.basename($out_dir)} Statistic Overview</h1>\n"
      end
      tmp_file << line.sub('<body>','')
    end
    overivew_file.close; tmp_file.close; `cp #{tmp_file.path} #{overivew_file.path}`


    %w(purin pyrimidin).each do |this_type|
      [0,1,2].each do |this_mis|

        create_plot = false

        length_a = []; count_a = []; percent_a = []; type_a = []; mismatches_a = []; gc_a = []; input_reads = []
        Dir.glob("#{$out_dir}/min*").each do |dir|
          params = File.basename(dir).split('_')
          length = params[0].split('-')[1].to_i
          mis = params[1].split('-')[1]
          gc = params[2].split('-')[1]

          if File.exist?("#{dir}/log.txt") && mis.to_s == this_mis.to_s
            create_plot = true
            length_a.push(length)
            mismatches_a.push(mis)
            gc_a.push(gc)
            File.open("#{dir}/log.txt",'r').each do |line|
              input_reads.push(line.split("\t")[1].chomp) if line.start_with?('input_reads')
              count_a.push(line.split("\t")[1].to_i) if line.start_with?("candidate_#{this_type}_reads")
              percent_a.push(line.split("\t")[1].sub('%','').to_f) if line.start_with?("#{this_type}_rate")
            end
          end
        end

        if create_plot

          # sort length and count arrays
          length_h = {}; c = 0
          length_a.each do |length|
            length_h[length] = c
            c += 1
          end

          length_a_sorted = []; count_a_sorted = []; percent_a_sorted = []

          length_a.sort!
          length_a.each do |length|
            pos = length_h[length]
            count_a_sorted.push(count_a[pos])
            percent_a_sorted.push(percent_a[pos])
            length_a_sorted.push(length.to_s)
          end

          #####################
          ## ABSOLUTE READ COUNT PLOTS
          if false
          header = "Absolute read counts of #{File.basename($out_dir)}: poly-#{this_type}, allowed mismatches=#{this_mis}, G+C>=#{gc_a[0]}, input reads=#{input_reads[0]}"
          plot2 = Nyaplot::Plot.new
          plot2.x_label('Minimum length [nt]')
          if this_type == 'purin'
            plot2.y_label('# poly-R [reads]')
          else
            plot2.y_label('# poly-Y [reads]')
          end
          bar = plot2.add(:bar, length_a_sorted, count_a_sorted)
          bar.color(colors)
          plot2.legend(true)


          plot2.export_html('tmp.html')
          html_file = File.open('tmp.html','r')
          html = ''; read = false
          html_file.each do |l|
            read = false if l.start_with?('</body>')
            read = true if l.start_with?('<body>')
            if read
              html << l.sub('<body>','')
            end
          end
          html_file.close; File.delete(html_file)

          overivew_file = File.open("#{$out_dir}/overview.html",'r')
          tmp_file = File.open("#{$out_dir}/overview_tmp.html",'w')
          overivew_file.each do |line|
            if line.start_with?('</body>')
              tmp_file << "\n\n<br><br><br><h2>#{header}</h2>\n"
              tmp_file << html
            end
            tmp_file << line
          end
          overivew_file.close; tmp_file.close; `cp #{tmp_file.path} #{overivew_file.path}`
          end

          #####################
          ## PERCENT READ COUNT PLOTS
          header = "Percentage of read counts of #{File.basename($out_dir)}: poly-#{this_type}, allowed mismatches=#{this_mis}, G+C>=#{gc_a[0]}, input reads=#{input_reads[0]}"
          plot2 = Nyaplot::Plot.new
          plot2.x_label('Minimum length [nt]')
          if this_type == 'purin'
            plot2.y_label('# poly-R [reads]')
          else
            plot2.y_label('# poly-Y [reads]')
          end
          bar = plot2.add(:bar, length_a_sorted, percent_a_sorted)
          bar.color(colors)
          plot2.legend(true)


          plot2.export_html('tmp.html')
          html_file = File.open('tmp.html','r')
          html = ''; read = false
          html_file.each do |l|
            read = false if l.start_with?('</body>')
            read = true if l.start_with?('<body>')
            if read
              html << l.sub('<body>','')
            end
          end
          html_file.close; File.delete(html_file)

          overivew_file = File.open("#{$out_dir}/overview.html",'r')
          tmp_file = File.open("#{$out_dir}/overview_tmp.html",'w')
          overivew_file.each do |line|
            if line.start_with?('</body>')
              tmp_file << "\n\n<br><br><br><h2>#{header}</h2>\n"
              tmp_file << html
            end
            tmp_file << line
          end
          overivew_file.close; tmp_file.close; `cp #{tmp_file.path} #{overivew_file.path}`

        end
      end
    end



  end


end

fastq = '../Bierhoff-TY1-03112014_S1_L001_R1_001.cutadapt.trimmed.internal.fastq'
#fastq = '../sortmerna/Bierhoff-TY1-03112014_S1_L001_R1_001_other.cutadapt.trimmed.internal.fastq'
#fastq = '/mnt/manja-geborgt/reads/reads/cellprogramming_fli_sarmistha/shSCR_iPSc_1_TGACCA_L006_R1_001.fastq'

min_guanine_cytosine_rate = 0.4

mismatches = 0
(12..25).each do |min_stretch_length|
  #SearchPolyStretches.new(fastq, min_stretch_length, mismatches, min_guanine_cytosine_rate, false)
end

mismatches = 1
(12..25).each do |min_stretch_length|
  #SearchPolyStretches.new(fastq, min_stretch_length, mismatches, min_guanine_cytosine_rate, false)
end

mismatches = 2
(12..25).each do |min_stretch_length|
  #SearchPolyStretches.new(fastq, min_stretch_length, mismatches, min_guanine_cytosine_rate, false)
end

SearchPolyStretches.new(fastq, nil, nil, nil, true)

