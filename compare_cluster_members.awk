BEGIN {
  FS = "[ \t]+"
  coverage = coverage > 0 ? coverage : 1
  read_file += 0
}

FNR == 1 {
  file_count = ARGIND
}

!/^#/ {
  members[ARGIND, $4]++
}

END {
  if (!read_file)
    read_file = file_count == 1 ? 1 : 2
  for (key in members) {
    split(key, part, SUBSEP)
    file = part[1]
    copies = file == read_file ? int(members[key] / coverage + 0.5) : members[key]
    histogram[file, copies]++
    clusters[file]++
    if (copies > max_copies)
      max_copies = copies
  }

  printf "#cluster_members"
  for (file = 1; file <= file_count; file++)
    printf "\t%s_count\t%s_pct", ARGV[file], ARGV[file]
  printf "\n"
    printf "#read_file\t%s\n#estimated_haploid_coverage\t%.3f\n",
      ARGV[read_file], coverage

  for (copies = 1; copies <= max_copies; copies++) {
    printf "%d", copies
    for (file = 1; file <= file_count; file++) {
      count = histogram[file, copies]
      printf "\t%d\t%.3f", count, 100 * count / clusters[file]
    }
    printf "\n"
  }
}