#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/kseq.h"
#include "segtrace.h"

/* Reader initialization */
KSEQ_INIT(gzFile, gzread)

// ==============================================================
// SECTION 1: ENTRY POINT & CLI PARSING
// ==============================================================

static void print_usage(void) {
  printf("Segtrace: Segmental Duplication Tracer\n\n"
         "Usage: segtrace [options] fasta1 [fasta2 ...]\n\n"
         "Options:\n"
         "  -k: kmer size (default: 17)\n"
         "  -s: scale factor (default: 16)\n"
         "  -w: window size in bp (default: 1024)\n"
         "  -t: step size in bp (default: 0 [auto: 33%% of window size])\n"
         "  -b: minimum valid bases per window (default: 0 [auto: 25%% of "
         "window size])\n"
         "  -c: minimum copies per genome/file to report (default: 2 "
         "[1=duplication map, >=3=polyploid])\n"
         "  -m: filter soft-masked bases (treat lowercase a/c/g/t as invalid)\n"
         "  -o: output file prefix (default: segtrace)\n"
         "  -p: number of threads (default: 8)\n"
         "  -h, --help: show this help message\n\n");
}

int main(int argc, char **argv) {
  if (argc < 2) {
    print_usage();
    return 1;
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      print_usage();
      return 0;
    }
  }

  uint32_t kmer_size = 17;
  uint64_t scale = 16;
  size_t window_size = 1024, step_size = 0, min_bases = 0;
  uint32_t min_copies = 2;
  const char *out_prefix = "segtrace";
  int n_threads = 8, filter_masked = 0;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:w:t:b:c:o:p:mh", 0)) >= 0) {
    if (c == 'h') {
      print_usage();
      return 0;
    } else if (c == 'k')
      kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      scale = (uint64_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'w')
      window_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 't')
      step_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'b')
      min_bases = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'c')
      min_copies = (uint32_t)atoi(opt.arg);
    else if (c == 'o')
      out_prefix = opt.arg;
    else if (c == 'p') {
      n_threads = atoi(opt.arg);
      if (n_threads < 1)
        n_threads = 1;
    }
    else if (c == 'm')
      filter_masked = 1;
    else
      return 1;
  }
  if (min_copies < 1)
    min_copies = 1;
  if (step_size == 0)
    step_size = window_size / 3;
  if (min_bases == 0)
    min_bases = window_size / 4;
  if (opt.ind == argc) {
    fprintf(stderr, "[ERROR] Input FASTA files are required.\n");
    return 1;
  }

  int num_files = argc - opt.ind;
  char **files = &argv[opt.ind];

  Segtrace r = {.hash_window = kmer_size, .hash_seed = 42};
  memset(r.base_lookup, -1, sizeof(r.base_lookup));
  for (int8_t code = 0; code < 4; code++) {
    uint8_t base = (uint8_t)"ACGT"[code];
    r.base_lookup[base] = code;
    if (!filter_masked)
      r.base_lookup[base + ('a' - 'A')] = code;
  }

  GlobalWindows gw =
      extract_all_windows(files, num_files, &r, scale, window_size,
                          step_size, min_bases, n_threads);

  fprintf(stderr,
          "[segtrace] Discovering candidates and computing distances...\n");
  CandidateGraph graph =
      discover_and_compute(gw.all_hashes, gw.coords, gw.num_sketches,
                           window_size, step_size, n_threads, r.hash_window);

  free(gw.all_hashes);

  SegtraceDupRegion *dup_regions = NULL;
  size_t n_dup_regions = 0;
  build_duplicate_loci(&graph, gw.num_sketches, gw.coords, gw.seq_lens,
                       step_size, window_size, &dup_regions, &n_dup_regions);
  cluster_duplicate_loci(&graph, gw.coords, dup_regions, n_dup_regions);
  free_candidate_graph(&graph);

  free(gw.coords);

  size_t n_filtered =
      filter_regions_by_copy_count(dup_regions, n_dup_regions, min_copies);

  write_dup_bed(out_prefix, dup_regions, n_filtered, gw.seq_lens,
                window_size < MIN_SD_LEN ? window_size : MIN_SD_LEN);

  free(dup_regions);
  for (size_t i = 0; i < gw.num_seqs; i++) {
    free(gw.seq_lens[i].genome);
    free(gw.seq_lens[i].seq);
  }
  free(gw.seq_lens);
  return 0;
}

// ==============================================================
// SECTION 2: ROLLING NTHASH & WINDOW EXTRACTION
// ==============================================================

static inline uint64_t rol64(uint64_t v, unsigned int n) {
  n &= 63;
  return (v << n) | (v >> ((64 - n) & 63));
}

static inline uint64_t ror64(uint64_t v, unsigned int n) {
  n &= 63;
  return (v >> n) | (v << ((64 - n) & 63));
}

static const uint64_t NTHASH_H[4] = {
    0x3c8bf4f53c8bf4f5ULL, // A
    0x04c903a704c903a7ULL, // C
    0x2b8104c92b8104c9ULL, // G
    0x2e0600d3fd09e083ULL  // T
};

static inline void extract_hash_direct(const Segtrace *r, uint32_t *out_hashes,
                                       size_t *out_size, uint32_t threshold,
                                       const uint8_t *seq, size_t len) {
  uint32_t k = r->hash_window;
  if (len < k) {
    *out_size = 0;
    return;
  }

  uint64_t f_hash = 0, r_hash = 0;
  size_t valid_len = 0;
  size_t count = 0;

  for (size_t i = 0; i < len; i++) {
    int8_t b = r->base_lookup[seq[i]];
    if (b < 0) {
      valid_len = 0;
      f_hash = 0;
      r_hash = 0;
      continue;
    }

    if (valid_len < k) {
      int8_t b_rc = b ^ 3;
      f_hash ^= rol64(NTHASH_H[b], k - 1 - (uint32_t)valid_len);
      r_hash ^= rol64(NTHASH_H[b_rc], (uint32_t)valid_len);
      valid_len++;
    } else {
      int8_t b_out = r->base_lookup[seq[i - k]];
      f_hash = rol64(f_hash, 1) ^ rol64(NTHASH_H[b_out], k) ^ NTHASH_H[b];
      r_hash = ror64(r_hash, 1) ^ ror64(NTHASH_H[b_out ^ 3], 1) ^
               rol64(NTHASH_H[b ^ 3], k - 1);
    }

    if (valid_len >= k) {
      uint64_t canonical = (f_hash < r_hash) ? f_hash : r_hash;
      uint32_t h = mix_hash(canonical, r->hash_seed);
      if (h < threshold && count < 2048) {
        out_hashes[count++] = h;
      }
    }
  }

  if (count > 1) {
    qsort(out_hashes, count, sizeof(uint32_t), compare_uint32);
    size_t u = 0;
    for (size_t i = 0; i < count; i++) {
      if (i == 0 || out_hashes[i] != out_hashes[i - 1])
        out_hashes[u++] = out_hashes[i];
    }
    count = u;
  }
  *out_size = count;
}

static void seq_chunk_worker(void *data, long i, int tid) {
  (void)tid;
  SeqChunkJob *job = &((SeqChunkJob *)data)[i];
  uint32_t current_window_idx =
      (uint32_t)(job->chunk_start_idx / job->step_size);

  uint32_t local_hashes[2048];

  size_t idx = job->chunk_start_idx;
  size_t valid_bases = 0;
  for (size_t j = 0; j < job->window_size; j++) {
    if (job->r->base_lookup[job->seq_ptr[idx + j]] >= 0)
      valid_bases++;
  }

  for (; idx + job->window_size <= job->chunk_end_idx;
       idx += job->step_size, current_window_idx++) {

    size_t sketch_size = 0;
    if (valid_bases >= job->min_bases) {
      extract_hash_direct(job->r, local_hashes, &sketch_size, job->threshold,
                          job->seq_ptr + idx, job->window_size);
    }

    DA_RESERVE(job->coords, job->cap_coords, job->num_coords + 1);
    WindowCoord *wc = &job->coords[job->num_coords++];
    wc->seq_id = job->seq_id;
    wc->window_idx = current_window_idx;
    wc->sketch_size = (uint16_t)sketch_size;

    size_t h_idx = job->num_hashes;
    if (sketch_size > 0) {
      DA_RESERVE(job->hashes, job->cap_hashes, job->num_hashes + sketch_size);
      memcpy(job->hashes + h_idx, local_hashes, sketch_size * sizeof(uint32_t));
      job->num_hashes += sketch_size;
    }
    wc->sketch_offset = (uint32_t)h_idx;

    if (idx + job->step_size + job->window_size <= job->chunk_end_idx) {
      for (size_t k = 0; k < job->step_size; k++) {
        if (job->r->base_lookup[job->seq_ptr[idx + k]] >= 0)
          valid_bases--;
        if (job->r->base_lookup[job->seq_ptr[idx + job->window_size + k]] >=
            0)
          valid_bases++;
      }
    }
  }
}

GlobalWindows extract_all_windows(char **files, int num_files,
                                  const Segtrace *r, uint64_t scale,
                                  size_t window_size, size_t step_size,
                                  size_t min_bases, int n_threads) {
  fprintf(stderr, "[segtrace] Extracting windows across genomes...\n");
  GlobalWindows gw = {0};
  uint32_t threshold = (uint32_t)(UINT32_MAX / scale);

  for (int f = 0; f < num_files; f++) {
    char bname[256];
    get_basename(files[f], bname, sizeof(bname));

    gzFile fp = gzopen(files[f], "r");
    if (!fp) {
      fprintf(stderr, "[WARNING] Failed to open FASTA file: %s\n", files[f]);
      continue;
    }
    kseq_t *ks = kseq_init(fp);
    if (!ks) {
      gzclose(fp);
      continue;
    }

    while (kseq_read(ks) >= 0) {
      size_t len = ks->seq.l;
      if (len < window_size)
        continue;

      DA_RESERVE(gw.seq_lens, gw.cap_seqs, gw.num_seqs + 1);
      gw.seq_lens[gw.num_seqs].genome = strdup(bname);
      gw.seq_lens[gw.num_seqs].seq = strdup(ks->name.s);
      gw.seq_lens[gw.num_seqs].file_id = (uint32_t)f;
      uint32_t seq_id = (uint32_t)gw.num_seqs++;

      uint8_t *seq_ptr = (uint8_t *)ks->seq.s;

      size_t est_windows =
          (len >= window_size) ? ((len - window_size) / step_size + 1) : 0;
      DA_RESERVE(gw.coords, gw.cap_sketches,
                 gw.num_sketches + est_windows + 16);
      DA_RESERVE(gw.all_hashes, gw.cap_all_hashes,
                 gw.num_all_hashes + (est_windows + 16) * 96);

      size_t chunk_size = len / (n_threads * 4);
      if (chunk_size < 100000)
        chunk_size = 100000;
      chunk_size = ((chunk_size + step_size - 1) / step_size) * step_size;

      size_t cap_jobs = 16, num_jobs = 0;
      SeqChunkJob *jobs = malloc(cap_jobs * sizeof(SeqChunkJob));

      for (size_t c_start = 0; c_start <= len - window_size;
           c_start += chunk_size) {
        size_t c_end = c_start + chunk_size + window_size - step_size;
        if (c_end > len)
          c_end = len;

        DA_RESERVE(jobs, cap_jobs, num_jobs + 1);
        jobs[num_jobs++] = (SeqChunkJob){.r = r,
                                         .threshold = threshold,
                                         .window_size = window_size,
                                         .step_size = step_size,
                                         .min_bases = min_bases,
                                         .seq_id = seq_id,
                                         .seq_ptr = seq_ptr,
                                         .chunk_start_idx = c_start,
                                         .chunk_end_idx = c_end};
      }

      kt_for(n_threads, seq_chunk_worker, jobs, num_jobs);

      for (size_t j = 0; j < num_jobs; j++) {
        SeqChunkJob *job = &jobs[j];
        if (job->num_coords > 0) {
          DA_RESERVE(gw.all_hashes, gw.cap_all_hashes,
                     gw.num_all_hashes + job->num_hashes);
          DA_RESERVE(gw.coords, gw.cap_sketches,
                     gw.num_sketches + job->num_coords);
          size_t base_h_offset = gw.num_all_hashes;
          memcpy(gw.all_hashes + base_h_offset, job->hashes,
                 job->num_hashes * sizeof(uint32_t));
          gw.num_all_hashes += job->num_hashes;

          size_t base_c_offset = gw.num_sketches;
          for (size_t k = 0; k < job->num_coords; k++) {
            WindowCoord wc = job->coords[k];
            wc.sketch_offset += (uint32_t)base_h_offset;
            gw.coords[base_c_offset + k] = wc;
          }
          gw.num_sketches += job->num_coords;
        }
        free(job->hashes);
        free(job->coords);
      }
      free(jobs);
    }
    kseq_destroy(ks);
    gzclose(fp);
  }
  return gw;
}

// ==============================================================
// SECTION 3: CANDIDATE DISCOVERY & DISTANCE COMPUTATION
// ==============================================================

static inline size_t calculate_sketch_dist(const uint32_t *a, size_t n_a,
                                           const uint32_t *b, size_t n_b) {
  size_t i = 0, j = 0, shared = 0;
  while (i < n_a && j < n_b) {
    uint32_t va = a[i], vb = b[j];
    shared += (va == vb);
    i += (va <= vb);
    j += (va >= vb);
  }
  return shared;
}

static inline size_t calculate_window_dist(const uint32_t *all_hashes,
                                           const WindowCoord *wa,
                                           const WindowCoord *wb) {
  return calculate_sketch_dist(all_hashes + wa->sketch_offset, wa->sketch_size,
                               all_hashes + wb->sketch_offset, wb->sketch_size);
}

static inline int windows_overlap(const DiscoverComputeData *w, uint32_t wa,
                                  uint32_t wb) {
  if (w->coords[wa].seq_id != w->coords[wb].seq_id)
    return 0;
  size_t window_distance =
      (size_t)ABS_DIFF(w->coords[wa].window_idx, w->coords[wb].window_idx);
  return window_distance * w->step_size < w->window_size;
}

static inline size_t max_sketch_size(const DiscoverComputeData *w, uint32_t wa,
                                     uint32_t wb) {
  return w->coords[wa].sketch_size > w->coords[wb].sketch_size
             ? w->coords[wa].sketch_size
             : w->coords[wb].sketch_size;
}

static inline size_t required_shared(const DiscoverComputeData *w, uint32_t wa,
                                     uint32_t wb) {
  size_t max_size = max_sketch_size(w, wa, wb);
  size_t min_shared = (size_t)ceil((double)max_size * w->p_kmer);
  return min_shared < 3 ? 3 : min_shared;
}

static inline int matching_window_pair(const DiscoverComputeData *w,
                                       long long wa, long long wb,
                                       uint32_t seq_a, uint32_t seq_b) {
  if (wa < 0 || wa >= (long long)w->n_windows || wb < 0 ||
      wb >= (long long)w->n_windows || w->coords[wa].seq_id != seq_a ||
      w->coords[wb].seq_id != seq_b ||
      windows_overlap(w, (uint32_t)wa, (uint32_t)wb))
    return 0;

  size_t min_shared = required_shared(w, (uint32_t)wa, (uint32_t)wb);
  return calculate_window_dist(w->all_hashes, &w->coords[wa],
                               &w->coords[wb]) >= min_shared;
}

static inline int check_collinear_neighbor(const DiscoverComputeData *w,
                                           uint32_t wa, uint32_t wb) {
  uint32_t seq_a = w->coords[wa].seq_id;
  uint32_t seq_b = w->coords[wb].seq_id;

  const int dir_a[] = {1, -1, 1, -1};
  const int dir_b[] = {1, -1, -1, 1};

  // Pass 1: Exact diagonals
  for (int d = 0; d < 4; d++) {
    int da = dir_a[d], db = dir_b[d];
    for (int step = 1; step <= MAX_COLLINEAR_LOOKAHEAD; step++) {
      long long next_a = (long long)wa + da * step;
      long long next_b = (long long)wb + db * step;
      if (matching_window_pair(w, next_a, next_b, seq_a, seq_b))
        return 1;
    }
  }

  // Pass 2: Gapped / Indel lookahead
  for (int d = 0; d < 4; d++) {
    int da = dir_a[d], db = dir_b[d];
    for (int step_a = 1; step_a <= MAX_COLLINEAR_LOOKAHEAD; step_a++) {
      for (int step_b = 1; step_b <= MAX_COLLINEAR_LOOKAHEAD; step_b++) {
        if (step_a == step_b)
          continue;
        long long next_a = (long long)wa + da * step_a;
        long long next_b = (long long)wb + db * step_b;
        if (matching_window_pair(w, next_a, next_b, seq_a, seq_b))
          return 1;
      }
    }
  }

  return 0;
}

static void discover_compute_worker(void *data, long idx, int tid) {
  DiscoverComputeData *w_data = (DiscoverComputeData *)data;
  long p = (long)w_data->batch_start + idx;
  PartitionBucket *b = &w_data->buckets[p];
  if (b->size == 0)
    return;

  qsort(b->entries, b->size, sizeof(HashWindowEntry), compare_hash_entry);

  size_t i = 0;
  while (i < b->size) {
    size_t j = i + 1;
    while (j < b->size && b->entries[j].hash == b->entries[i].hash)
      j++;
    size_t run_len = j - i;

    if (run_len >= 2 && run_len <= MAX_KMER_FREQ) {
      for (size_t a = i; a < j; a++) {
        size_t b_max =
            a + 1 + MAX_PAIR_COMPARISONS < j ? a + 1 + MAX_PAIR_COMPARISONS : j;
        for (size_t b_idx = a + 1; b_idx < b_max; b_idx++) {
          uint32_t wa = b->entries[a].window_id,
                   wb = b->entries[b_idx].window_id;
          if (windows_overlap(w_data, wa, wb))
            continue;

          uint64_t pk = encode_pair(wa, wb);
          if (bloom_test_and_set(w_data->t_bloom[tid], pk))
            continue;

          size_t min_shared = required_shared(w_data, wa, wb);

          size_t shared = calculate_window_dist(
              w_data->all_hashes, &w_data->coords[wa], &w_data->coords[wb]);
          if (shared < min_shared ||
              !check_collinear_neighbor(w_data, wa, wb))
            continue;

          size_t max_size = max_sketch_size(w_data, wa, wb);
          uint32_t score =
              (uint32_t)((shared * UINT8_MAX + max_size / 2) / max_size);
          DA_PUSH(w_data->t_pairs[tid], w_data->t_n_pairs[tid],
                  w_data->t_cap_pairs[tid],
                  ((CandidatePair){
                      wa | ((score & UINT32_C(0x0f))
                            << CANDIDATE_SCORE_SHIFT),
                      wb | ((score >> 4) << CANDIDATE_SCORE_SHIFT)}));
        }
      }
    }
    i = j;
  }
}

CandidateGraph discover_and_compute(const uint32_t *all_hashes,
                                    const WindowCoord *coords,
                                    size_t n_windows, size_t window_size,
                                    size_t step_size, int n_threads,
                                    uint32_t kmer_size) {
  if (n_windows > (size_t)CANDIDATE_WINDOW_MASK + 1) {
    fprintf(stderr, "[ERROR] Too many windows for candidate encoding\n");
    exit(1);
  }

  DiscoverComputeData w = {
      .all_hashes = all_hashes,
      .coords = coords,
      .n_windows = n_windows,
      .window_size = window_size,
      .step_size = step_size,
      .p_kmer = pow(MIN_IDENTITY, (double)kmer_size),
      .buckets = calloc(NUM_PARTITIONS, sizeof(PartitionBucket)),
      .t_bloom = malloc(n_threads * sizeof(uint8_t *)),
      .t_pairs = calloc(n_threads, sizeof(CandidatePair *)),
      .t_n_pairs = calloc(n_threads, sizeof(size_t)),
      .t_cap_pairs = calloc(n_threads, sizeof(size_t))};

  if (!w.buckets || !w.t_bloom || !w.t_pairs) {
    fprintf(stderr, "[ERROR] Memory allocation failed\n");
    exit(1);
  }
  for (int t = 0; t < n_threads; t++)
    w.t_bloom[t] = calloc(BLOOM_SIZE_BYTES, 1);

  uint32_t part_size = (uint32_t)(UINT32_MAX / NUM_PARTITIONS);
  uint16_t *win_curr_pos = calloc(n_windows, sizeof(uint16_t));

  for (size_t batch_start = 0; batch_start < NUM_PARTITIONS;
       batch_start += BATCH_PARTITIONS) {
    size_t batch_end = batch_start + BATCH_PARTITIONS;
    if (batch_end > NUM_PARTITIONS)
      batch_end = NUM_PARTITIONS;

    w.batch_start = batch_start;
    size_t batch_count = batch_end - batch_start;
    uint32_t max_hash = (batch_end < NUM_PARTITIONS)
                            ? (uint32_t)(batch_end * part_size)
                            : UINT32_MAX;

    for (size_t win = 0; win < n_windows; win++) {
      uint32_t off = coords[win].sketch_offset;
      uint16_t sz = coords[win].sketch_size;
      uint16_t pos = win_curr_pos[win];

      while (pos < sz) {
        uint32_t val = all_hashes[off + pos];
        if (val >= max_hash)
          break;
        size_t p = (size_t)(val / part_size);
        if (p >= NUM_PARTITIONS)
          p = NUM_PARTITIONS - 1;
        DA_PUSH(w.buckets[p].entries, w.buckets[p].size, w.buckets[p].cap,
                ((HashWindowEntry){val, (uint32_t)win}));
        pos++;
      }
      win_curr_pos[win] = pos;
    }

    kt_for(n_threads, discover_compute_worker, &w, (long)batch_count);

    for (size_t p = batch_start; p < batch_end; p++) {
      free(w.buckets[p].entries);
    }
  }

  for (int t = 0; t < n_threads; t++)
    free(w.t_bloom[t]);
  free(w.buckets);
  free(w.t_bloom);
  free(w.t_cap_pairs);
  free(win_curr_pos);

  return (CandidateGraph){
      .pairs = w.t_pairs, .counts = w.t_n_pairs, .n_threads = n_threads};
}

// ==============================================================
// SECTION 4: CLUSTERING, LOCUS MERGING & COPY FILTERING
// ==============================================================

static inline uint32_t candidate_window(uint32_t encoded) {
  return encoded & CANDIDATE_WINDOW_MASK;
}

static inline uint32_t candidate_score(CandidatePair pair) {
  return (pair.a >> CANDIDATE_SCORE_SHIFT) |
         ((pair.b >> CANDIDATE_SCORE_SHIFT) << 4);
}

void build_duplicate_loci(const CandidateGraph *graph, size_t num_windows,
                          WindowCoord *coords, const GenomeSeqLen *seq_lens,
                          size_t step_size, size_t window_size,
                          SegtraceDupRegion **out_regions,
                          size_t *out_n_regions) {
  for (size_t i = 0; i < num_windows; i++) {
    coords[i].sketch_offset = UINT32_MAX;
    coords[i].sketch_size = 0;
  }
  for (int t = 0; t < graph->n_threads; t++) {
    for (size_t i = 0; i < graph->counts[t]; i++) {
      CandidatePair pair = graph->pairs[t][i];
      coords[candidate_window(pair.a)].sketch_size = 1;
      coords[candidate_window(pair.b)].sketch_size = 1;
    }
  }

  SegtraceDupRegion *regions = NULL;
  size_t n_regions = 0, cap_regions = 0;
  uint32_t previous_seq = UINT32_MAX;
  uint32_t previous_window = 0;

  for (size_t i = 0; i < num_windows; i++) {
    if (coords[i].sketch_size == 0)
      continue;

    uint32_t seq_id = coords[i].seq_id;
    uint32_t window_idx = coords[i].window_idx;
    size_t start = (size_t)window_idx * step_size;
    size_t end = start + window_size;

    if (n_regions > 0 && seq_id == previous_seq &&
        window_idx - previous_window <= MAX_COLLINEAR_LOOKAHEAD) {
      regions[n_regions - 1].end = end;
    } else {
      DA_PUSH(regions, n_regions, cap_regions,
              ((SegtraceDupRegion){.seq_id = seq_id,
                                   .file_id = seq_lens[seq_id].file_id,
                                   .start = start,
                                   .end = end,
                                   .cluster_id = 0,
                                   .partner_id = UINT32_MAX}));
    }
    coords[i].sketch_offset = (uint32_t)(n_regions - 1);
    previous_seq = seq_id;
    previous_window = window_idx;
  }

  *out_regions = regions;
  *out_n_regions = n_regions;
}

static int compare_dup_region_by_cluster_file(const void *a, const void *b) {
  const SegtraceDupRegion *ra = (const SegtraceDupRegion *)a,
                          *rb = (const SegtraceDupRegion *)b;
  if (ra->cluster_id != rb->cluster_id)
    return CMP(ra->cluster_id, rb->cluster_id);
  if (ra->file_id != rb->file_id)
    return CMP(ra->file_id, rb->file_id);
  if (ra->seq_id != rb->seq_id)
    return CMP(ra->seq_id, rb->seq_id);
  if (ra->start != rb->start)
    return CMP(ra->start, rb->start);
  return CMP(ra->end, rb->end);
}

void cluster_duplicate_loci(const CandidateGraph *graph,
                            const WindowCoord *coords,
                            SegtraceDupRegion *regions, size_t n_regions) {
  if (n_regions == 0)
    return;

  for (int t = 0; t < graph->n_threads; t++) {
    for (size_t i = 0; i < graph->counts[t]; i++) {
      CandidatePair pair = graph->pairs[t][i];
      uint32_t region_a = coords[candidate_window(pair.a)].sketch_offset;
      uint32_t region_b = coords[candidate_window(pair.b)].sketch_offset;
      if (region_a == UINT32_MAX || region_b == UINT32_MAX)
        continue;
      if (region_a == region_b) {
        if (regions[region_a].partner_id == UINT32_MAX)
          regions[region_a].partner_id = INTERNAL_DUPLICATION_ID;
        continue;
      }

      uint32_t score = candidate_score(pair);
      if (score > regions[region_a].cluster_id ||
          (score == regions[region_a].cluster_id &&
           region_b < regions[region_a].partner_id)) {
        regions[region_a].cluster_id = score;
        regions[region_a].partner_id = region_b;
      }
      if (score > regions[region_b].cluster_id ||
          (score == regions[region_b].cluster_id &&
           region_a < regions[region_b].partner_id)) {
        regions[region_b].cluster_id = score;
        regions[region_b].partner_id = region_a;
      }
    }
  }

  UnionFind uf;
  init_unionfind(&uf, n_regions);
  for (size_t i = 0; i < n_regions; i++) {
    if (regions[i].partner_id < n_regions)
      union_unionfind(&uf, (uint32_t)i, regions[i].partner_id);
  }

  uint32_t *cluster_map = calloc(n_regions, sizeof(uint32_t));
  uint32_t next_cluster_id = 1;
  for (size_t i = 0; i < n_regions; i++) {
    uint32_t root = find_unionfind(&uf, (uint32_t)i);
    if (cluster_map[root] == 0)
      cluster_map[root] = next_cluster_id++;
    regions[i].cluster_id = cluster_map[root];
  }
  free(cluster_map);
  free_unionfind(&uf);

  qsort(regions, n_regions, sizeof(SegtraceDupRegion),
        compare_dup_region_by_cluster_file);
}

void free_candidate_graph(CandidateGraph *graph) {
  for (int t = 0; t < graph->n_threads; t++)
    free(graph->pairs[t]);
  free(graph->pairs);
  free(graph->counts);
}

size_t filter_regions_by_copy_count(SegtraceDupRegion *regions, size_t n,
                                    uint32_t min_copies) {
  if (n == 0 || min_copies <= 1)
    return n;

  size_t out_count = 0;
  size_t i = 0;
  while (i < n) {
    size_t j = i + 1;
    while (j < n && regions[j].cluster_id == regions[i].cluster_id &&
           regions[j].file_id == regions[i].file_id) {
      j++;
    }
    size_t copy_count = j - i;
    if (copy_count >= min_copies ||
        (min_copies == 2 && copy_count == 1 &&
         regions[i].partner_id == INTERNAL_DUPLICATION_ID)) {
      for (size_t k = i; k < j; k++) {
        regions[out_count++] = regions[k];
      }
    }
    i = j;
  }
  return out_count;
}

void write_dup_bed(const char *out_prefix, const SegtraceDupRegion *dup_regions,
                   size_t n_merged, const GenomeSeqLen *seq_lens,
                   size_t min_sd_len) {
  if (n_merged == 0)
    return;

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
  FILE *out_bed = fopen(path_buf, "w");
  if (!out_bed) {
    fprintf(stderr, "[ERROR] Cannot open output file: %s\n", path_buf);
    return;
  }

  fprintf(out_bed, "#chrom\tstart\tend\tcluster_id\n");
  for (size_t k = 0; k < n_merged; k++) {
    if (dup_regions[k].end - dup_regions[k].start >= min_sd_len) {
      uint32_t seq_i = dup_regions[k].seq_id;
      fprintf(out_bed, "%s-%s\t%zu\t%zu\t%u\n", seq_lens[seq_i].genome,
              seq_lens[seq_i].seq, dup_regions[k].start, dup_regions[k].end,
              dup_regions[k].cluster_id);
    }
  }
  fclose(out_bed);
}

// ==============================================================
// SECTION 5: CORE ALGORITHMS & UTILITIES
// ==============================================================

void init_unionfind(UnionFind *uf, size_t n) {
  uf->parent = malloc(n * sizeof(uint32_t));
  uf->rank = calloc(n, sizeof(uint8_t));
  for (size_t i = 0; i < n; i++)
    uf->parent[i] = (uint32_t)i;
}

uint32_t find_unionfind(UnionFind *uf, uint32_t x) {
  while (uf->parent[x] != x) {
    uint32_t next = uf->parent[x];
    uf->parent[x] = uf->parent[next];
    x = next;
  }
  return x;
}

void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b) {
  uint32_t root_a = find_unionfind(uf, a), root_b = find_unionfind(uf, b);
  if (root_a != root_b) {
    if (uf->rank[root_a] < uf->rank[root_b])
      uf->parent[root_a] = root_b;
    else if (uf->rank[root_a] > uf->rank[root_b])
      uf->parent[root_b] = root_a;
    else {
      uf->parent[root_b] = root_a;
      uf->rank[root_a]++;
    }
  }
}

void free_unionfind(UnionFind *uf) {
  free(uf->parent);
  free(uf->rank);
}

void get_basename(const char *filename, char *basename, size_t size) {
  const char *last_slash = strrchr(filename, '/');
  const char *name = last_slash ? last_slash + 1 : filename;
  strncpy(basename, name, size - 1);
  basename[size - 1] = '\0';

  char *dot = strrchr(basename, '.');
  if (dot && (strcmp(dot, ".gz") == 0 || strcmp(dot, ".bgz") == 0)) {
    *dot = '\0';
  }
  dot = strrchr(basename, '.');
  if (dot && (strcmp(dot, ".fa") == 0 || strcmp(dot, ".fna") == 0 ||
              strcmp(dot, ".fasta") == 0 || strcmp(dot, ".fastq") == 0 ||
              strcmp(dot, ".fq") == 0)) {
    *dot = '\0';
  }
}

inline uint32_t mix_hash(uint64_t hash_value, uint64_t seed) {
  hash_value ^= seed;
  hash_value ^= hash_value >> 33;
  hash_value *= MIX_CONST1;
  hash_value ^= hash_value >> 33;
  hash_value *= MIX_CONST2;
  hash_value ^= hash_value >> 33;
  return (uint32_t)hash_value;
}

uint64_t encode_pair(uint32_t a, uint32_t b) {
  return a < b ? ((uint64_t)a << 32) | b : ((uint64_t)b << 32) | a;
}

static inline uint64_t splitmix64(uint64_t x) {
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  return x ^ (x >> 31);
}

int bloom_test_and_set(uint8_t *bloom, uint64_t key) {
  uint64_t h = splitmix64(key);
  uint32_t h1 = (uint32_t)h & BLOOM_MASK;
  uint32_t h2 = (uint32_t)(h >> 32) & BLOOM_MASK;
  int was_set =
      ((bloom[h1 >> 3] >> (h1 & 7)) & 1) & ((bloom[h2 >> 3] >> (h2 & 7)) & 1);
  bloom[h1 >> 3] |= (uint8_t)(1 << (h1 & 7));
  bloom[h2 >> 3] |= (uint8_t)(1 << (h2 & 7));
  return was_set;
}

int compare_uint32(const void *a, const void *b) {
  uint32_t va = *(const uint32_t *)a, vb = *(const uint32_t *)b;
  return (va > vb) - (va < vb);
}

int compare_hash_entry(const void *a, const void *b) {
  const HashWindowEntry *ea = (const HashWindowEntry *)a,
                        *eb = (const HashWindowEntry *)b;
  return ea->hash != eb->hash ? CMP(ea->hash, eb->hash)
                              : CMP(ea->window_id, eb->window_id);
}
