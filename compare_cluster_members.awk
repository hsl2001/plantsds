BEGIN {
  FS = "[ \t]+"
  coverage = coverage > 0 ? coverage : 1
  read_file += 0
  bins[1] = "2"
  bins[2] = "3"
  bins[3] = "4"
  bins[4] = "5"
  bins[5] = "6-10"
  bins[6] = "11-20"
  bins[7] = "21-50"
  bins[8] = "51+"
}

FNR == 1 {
  file_count = ARGIND
}

!/^#/ {
  members[ARGIND, $4]++
}

function bin(size) {
  return size <= 5 ? size : size <= 10 ? "6-10" :
         size <= 20 ? "11-20" : size <= 50 ? "21-50" : "51+"
}

END {
  if (!read_file)
    read_file = file_count == 1 ? 1 : 2
  for (key in members) {
    split(key, part, SUBSEP)
    file = part[1]
    copies = file == read_file ? int(members[key] / coverage + 0.5) : members[key]
    histogram[file, bin(copies)]++
    clusters[file]++
  }

  printf "#cluster_members"
  for (file = 1; file <= file_count; file++)
    printf "\t%s_count\t%s_pct", ARGV[file], ARGV[file]
  printf "\n"
    printf "#read_file\t%s\n#estimated_haploid_coverage\t%.3f\n",
      ARGV[read_file], coverage

  for (bin_index = 1; bin_index <= length(bins); bin_index++) {
    printf "%s", bins[bin_index]
    for (file = 1; file <= file_count; file++) {
      count = histogram[file, bins[bin_index]]
      printf "\t%d\t%.3f", count, 100 * count / clusters[file]
    }
    printf "\n"
  }
}