require 'bio'

# To filter the different alignments and find the orthologues I will select the hits with a evalue smaller than 10e-6 and an overlap >= 50%
# I took this parameters from the paper:
# Gabriel Moreno-Hagelsieb, Kristen Latimer, Choosing BLAST options for better detection of orthologs as reciprocal best hits, Bioinformatics, Volume 24, Issue 3, 1 February 2008, Pages 319–324, https://doi.org/10.1093/bioinformatics/btm585

def store_fastafile(file)
  
  fastafile = Bio::FastaFormat.open(file)
  return fastafile
end

def create_blastdatabase(file1, type1, file2, type2)
  
  system("rm -r Databases")
  system("mkdir Databases")
  system("makeblastdb -in #{file1} -dbtype '#{type1}' -out ./Databases/#{file1}")
  system("makeblastdb -in #{file2} -dbtype '#{type2}' -out ./Databases/#{file2}")
end

def find_orthologues(file1, file2, fastafile1, fastafile2)

best_hit_fasta1 = {}
best_hit_fasta2 = {}
common_hits = []

fastafile1.each do |gene|
  gene_id = gene.definition.split("|")[0]
  query = gene.seq.to_s
  seq = Bio::Sequence.auto(gene.seq.to_s)
  if seq.guess.to_s.split("::")[2] == "AA"
    factory1 = Bio::Blast.local('tblastn', "./Databases/#{file1}")
    report1 = factory1.query(query)
    report1.each do |hit|
      if hit.evalue < 10e-6 and hit.overlap/seq.length.to_f >= 0.5
        best_hit_fasta1[gene_id] = hit.definition.split("|")[0]
      end 
    end 
  elsif seq.guess.to_s.split("::")[2] == "NA"
    factory1 = Bio::Blast.local('blastx', "./Databases/#{file2}")
    report1 = factory1.query(query)
    report1.each do |hit|
      if hit.evalue < 10e-6 and hit.overlap/seq.translate.length.to_f >= 0.5
        best_hit_fasta1[gene_id] = hit.definition.split("|")[0]
      end 
    end
  end
end 
fastafile2.each do |gene|
  gene_id = gene.definition.split("|")[0]
  query = gene.seq.to_s
  seq = Bio::Sequence.auto(gene.seq.to_s)
  if seq.guess.to_s.split("::")[2] == "AA"
    factory2 = Bio::Blast.local('tblastn', "./Databases/#{file1}")
    report2 = factory2.query(query)
    report2.each do |hit|
      if hit.evalue < 10e-6 and hit.overlap/seq.length.to_f >= 0.5
        best_hit_fasta2[gene_id] = hit.definition.split("|")[0]
      end 
    end 
  elsif seq.guess.to_s.split("::")[2] == "NA"
    factory2 = Bio::Blast.local('blastx', "./Databases/#{file2}")
    report2 = factory2.query(query)
    report2.each do |hit|
      if hit.evalue < 10e-6 and hit.overlap/seq.translate.length.to_f >= 0.5
        best_hit_fasta2[gene_id] = hit.definition.split("|")[0]
      end 
    end
  end
end

best_hit_fasta1.each_key do |key|
  best_hit_fasta2.each_key do |key1|
  if best_hit_fasta1[key] == key1 and best_hit_fasta2[key1] == key
    orthologues = [key, key1]
    common_hits.append(orthologues)
  end 
  end 
end

if common_hits.length != 0
  puts ("#{common_hits.length} orthologues found")
else
  puts ("No orthologues")
end

return common_hits

end

def write_report(file_name, orthologues_list)
  out = File.open(file_name, "w")
  out.puts ("Pairs of orthologues found for the species Arabidopsis and Pombe")
  out.puts ("There are a total #{orthologues_list.length} pairs of orthologues which you can find below")
  orthologues_list.each do |pair|
    out.puts ("#{pair[0]}\t-\t#{pair[1]}")
  end  
  out.close
  puts ("The pairs of orthologues found have been written in the file #{file_name}")
end

fastafile1 = store_fastafile("Arabidopsis.fa")
fastafile2 = store_fastafile("pombe.fa")
create_blastdatabase("Arabidopsis.fa", "nucl", "pombe.fa", "prot")
common_hits = find_orthologues("Arabidopsis.fa", "pombe.fa", fastafile1, fastafile2)
write_report("orthologues_pairs.txt", common_hits)

