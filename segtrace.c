#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/khash.h"
#include "klib/kseq.h"
#include "segtrace.h"

/* Reader initialization */
KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(genome_map, uint32_t)

/* 2-bit nucleotide encoding: A = 00, C = 01, G = 10, T = 11. Soft-masked
 * (lowercase) = -1 */
const int8_t BASE_LOOKUP[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1};

const int8_t BASE_LOOKUP_NO_MASK[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1};

// ==============================================================
// SECTION 1: ENTRY POINT & CLI PARSING
// ==============================================================

void print_usage(void) {
  printf("Segtrace: Segmental Duplication Tracer\n\n"
         "Usage: segtrace [options] fasta1 [fasta2 ...]\n\n"
         "Options:\n"
         "  -k: kmer size (default: 17)\n"
         "  -s: scale factor (default: 16)\n"
         "  -w: window size in bp (default: 1024)\n"
         "  -t: step size in bp (default: 0 [auto: 33%% of window size])\n"
         "  -b: minimum valid bases per window (default: 0 [auto: 25%% of "
         "window size])\n"
         "  -m: do not filter soft-masked bases (treat lowercase a/c/g/t as "
         "valid)\n"
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
    if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      print_usage();
      return 0;
    }
  }

  uint32_t def_kmer_size = 17;
  uint64_t def_scale = 16, def_hash_seed = 42;
  size_t window_size = 1024, step_size = 0, min_bases = 0, flank_size = 256;
  const char *out_prefix = "segtrace";
  int n_threads = 8, filter_masked = 1;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:e:w:t:b:d:o:p:D:mh", 0)) >= 0) {
    if (c == 'h') {
      print_usage();
      return 0;
    } else if (c == 'k')
      def_kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      def_scale = (uint64_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'w')
      window_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 't')
      step_size = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'b')
      min_bases = (size_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'o')
      out_prefix = opt.arg;
    else if (c == 'p')
      n_threads = atoi(opt.arg) < 1 ? 1 : atoi(opt.arg);
    else if (c == 'm')
      filter_masked = 0;
    else
      return 1;
  }
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

  Segtrace r;
  init_segtrace(&r, def_kmer_size, filter_masked);
  r.hash_seed = def_hash_seed;

  uint32_t *all_hashes = NULL;
  WindowCoord *coords = NULL;
  size_t num_sketches = 0, num_seqs = 0;
  GenomeSeqLen *seq_lens = NULL;

  StreamWorkerData *workers =
      extract_all_windows(files, num_files, &r, def_scale, window_size,
                          step_size, min_bases, n_threads);
  merge_global_data(workers, num_files, out_prefix, &all_hashes, &coords,
                    &num_sketches, &seq_lens, &num_seqs);
  free(workers);

  UnionFind uf;
  init_unionfind(&uf, num_sketches);

  fprintf(stderr,
          "[segtrace] Discovering candidates and computing distances...\n");
  discover_and_compute(all_hashes, coords, num_sketches, window_size, step_size,
                       n_threads, r.hash_window, &uf);

  SegtraceDupRegion *dup_regions = NULL;
  size_t n_dup_regions = 0;
  build_duplicate_regions(&uf, num_sketches, coords, step_size, window_size,
                          &dup_regions, &n_dup_regions);

  free_unionfind(&uf);
  free(all_hashes);
  all_hashes = NULL;
  free(coords);
  coords = NULL;

  size_t n_merged = merge_dup_regions(dup_regions, n_dup_regions, window_size);

  fprintf(stderr,
          "[INFO] Extracting flanking sequences for sub-clustering...\n");
  extract_flankings(files, num_files, &r, def_scale, dup_regions, n_merged,
                    n_threads, flank_size, seq_lens, num_seqs);

  fprintf(stderr, "[INFO] Sub-clustering based on flanking similarities...\n");
  perform_subclustering(dup_regions, n_merged, n_threads, r.hash_window);

  write_dup_bed(out_prefix, dup_regions, n_merged, seq_lens);

  for (size_t i = 0; i < n_merged; i++) {
    free(dup_regions[i].flank_sketch.hashes);
  }
  free(dup_regions);
  for (size_t i = 0; i < num_seqs; i++) {
    free(seq_lens[i].genome);
    free(seq_lens[i].seq);
  }
  free(seq_lens);
  return 0;
}

// ==============================================================
// SECTION 2: WINDOW EXTRACTION & FASTA STREAMING
// ==============================================================

static void seq_chunk_worker(void *data, long i, int tid) {
  (void)tid;
  SeqChunkJob *job = &((SeqChunkJob *)data)[i];
  uint32_t current_window_idx =
      (uint32_t)(job->chunk_start_idx / job->step_size);

  for (size_t idx = job->chunk_start_idx;
       idx + job->window_size <= job->chunk_end_idx &&
       idx + job->window_size <= job->seq_len;
       idx += job->step_size, current_window_idx++) {
    size_t valid_bases = 0;
    for (size_t j = 0; j < job->window_size; j++) {
      if (job->base_lookup[job->seq_ptr[idx + j]] >= 0)
        valid_bases++;
    }
    size_t sketch_size = 0;
    uint32_t *hashes = NULL;
    if (valid_bases >= job->min_bases) {
      HashPool pool;
      init_hash_pool(&pool, UINT32_MAX / job->scale);
      extract_hash(job->r, &pool, job->seq_ptr + idx, job->window_size);
      finalize_hash_pool(&pool, &hashes, &sketch_size);
    }

    DA_RESERVE(job->coords, job->cap_coords, job->num_coords + 1);
    WindowCoord *wc = &job->coords[job->num_coords++];
    wc->seq_id = job->seq_id;
    wc->sketch_offset = 0;
    wc->window_idx = current_window_idx;
    wc->sketch_size = (uint16_t)sketch_size;
    wc->flags = 0;

    size_t h_idx = job->num_hashes;
    if (sketch_size > 0) {
      DA_RESERVE(job->hashes, job->cap_hashes, job->num_hashes + sketch_size);
      memcpy(job->hashes + h_idx, hashes, sketch_size * sizeof(uint32_t));
      job->num_hashes += sketch_size;
    }

    wc->sketch_offset = (uint32_t)h_idx;
    free(hashes);
  }
}

StreamWorkerData *extract_all_windows(char **files, int num_files,
                                      const Segtrace *r, uint64_t scale,
                                      size_t window_size, size_t step_size,
                                      size_t min_bases, int n_threads) {
  fprintf(stderr, "[segtrace] Extracting windows across genomes...\n");
  StreamWorkerData *workers = calloc(1, sizeof(StreamWorkerData));
  if (!workers) {
    fprintf(stderr, "[ERROR] Memory allocation failed\n");
    exit(1);
  }

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

      DA_RESERVE(workers[0].seq_lens, workers[0].cap_seqs,
                 workers[0].num_seqs + 1);
      workers[0].seq_lens[workers[0].num_seqs].genome = strdup(bname);
      workers[0].seq_lens[workers[0].num_seqs].seq = strdup(ks->name.s);
      uint32_t seq_id = (uint32_t)workers[0].num_seqs++;

      uint8_t *seq_copy = malloc(len + 1);
      memcpy(seq_copy, ks->seq.s, len + 1);

      size_t chunk_size = len / (n_threads * 4);
      if (chunk_size < 100000)
        chunk_size = 100000;
      chunk_size = ((chunk_size + step_size - 1) / step_size) * step_size;

      size_t cap_jobs = 16, num_jobs = 0;
      SeqChunkJob *jobs = malloc(cap_jobs * sizeof(SeqChunkJob));

      for (size_t c_start = 0; c_start < len; c_start += chunk_size) {
        size_t c_end = c_start + chunk_size + window_size - step_size;
        if (c_start + window_size > len)
          break;

        DA_RESERVE(jobs, cap_jobs, num_jobs + 1);
        jobs[num_jobs++] = (SeqChunkJob){.r = r,
                                         .base_lookup = r->base_lookup,
                                         .scale = scale,
                                         .window_size = window_size,
                                         .step_size = step_size,
                                         .min_bases = min_bases,
                                         .seq_id = seq_id,
                                         .seq_ptr = seq_copy,
                                         .seq_len = len,
                                         .chunk_start_idx = c_start,
                                         .chunk_end_idx = c_end,
                                         .hashes = NULL,
                                         .num_hashes = 0,
                                         .cap_hashes = 0,
                                         .coords = NULL,
                                         .num_coords = 0,
                                         .cap_coords = 0};
      }

      kt_for(n_threads, seq_chunk_worker, jobs, num_jobs);

      for (size_t j = 0; j < num_jobs; j++) {
        SeqChunkJob *job = &jobs[j];
        if (job->num_coords > 0) {
          size_t base_h_offset = workers[0].num_all_hashes;
          DA_RESERVE(workers[0].all_hashes, workers[0].cap_all_hashes,
                     workers[0].num_all_hashes + job->num_hashes);
          memcpy(workers[0].all_hashes + base_h_offset, job->hashes,
                 job->num_hashes * sizeof(uint32_t));
          workers[0].num_all_hashes += job->num_hashes;

          size_t base_c_offset = workers[0].num_sketches;
          DA_RESERVE(workers[0].coords, workers[0].cap_sketches,
                     workers[0].num_sketches + job->num_coords);

          for (size_t k = 0; k < job->num_coords; k++) {
            WindowCoord wc = job->coords[k];
            wc.sketch_offset += (uint32_t)base_h_offset;
            workers[0].coords[base_c_offset + k] = wc;
          }
          workers[0].num_sketches += job->num_coords;
        }
        free(job->hashes);
        free(job->coords);
      }
      free(jobs);
      free(seq_copy);
    }
    kseq_destroy(ks);
    gzclose(fp);
  }
  return workers;
}

void merge_global_data(StreamWorkerData *workers, int num_files,
                       const char *out_prefix, uint32_t **out_all_hashes,
                       WindowCoord **out_coords, size_t *out_num_sketches,
                       GenomeSeqLen **out_seq_lens, size_t *out_num_seqs) {
  (void)out_prefix;
  (void)num_files;
  size_t total_sketches = workers[0].num_sketches;
  size_t total_seqs = workers[0].num_seqs;

  *out_all_hashes = workers[0].all_hashes;
  *out_coords = workers[0].coords;
  *out_num_sketches = total_sketches;
  *out_seq_lens = workers[0].seq_lens;
  *out_num_seqs = total_seqs;

  workers[0].all_hashes = NULL;
  workers[0].coords = NULL;
  workers[0].seq_lens = NULL;
}

// ==============================================================
// SECTION 3: CANDIDATE DISCOVERY & DISTANCE COMPUTATION
// ==============================================================

void discover_and_compute(const uint32_t *all_hashes, WindowCoord *coords,
                          size_t n_windows, size_t window_size,
                          size_t step_size, int n_threads, uint32_t kmer_size,
                          UnionFind *uf) {
  DiscoverComputeData w = {.all_hashes = all_hashes,
                           .coords = coords,
                           .n_windows = n_windows,
                           .window_size = window_size,
                           .step_size = step_size,
                           .kmer_size = kmer_size,
                           .uf = uf,
                           .buckets =
                               calloc(NUM_PARTITIONS, sizeof(PartitionBucket)),
                           .t_bloom = malloc(n_threads * sizeof(uint8_t *))};
  if (!w.buckets || !w.t_bloom) {
    fprintf(stderr, "[ERROR] Memory allocation failed\n");
    exit(1);
  }
  for (int t = 0; t < n_threads; t++)
    w.t_bloom[t] = calloc(BLOOM_SIZE_BYTES, 1);

  uint32_t part_size = UINT32_MAX / NUM_PARTITIONS;

  for (size_t batch_start = 0; batch_start < NUM_PARTITIONS;
       batch_start += BATCH_PARTITIONS) {
    size_t batch_end = batch_start + BATCH_PARTITIONS;
    if (batch_end > NUM_PARTITIONS)
      batch_end = NUM_PARTITIONS;

    w.batch_start = batch_start;
    size_t batch_count = batch_end - batch_start;

    for (size_t p = batch_start; p < batch_end; p++) {
      w.buckets[p].size = 0;
      w.buckets[p].cap = 0;
      w.buckets[p].entries = NULL;
    }

    for (size_t win = 0; win < n_windows; win++) {
      const uint32_t *h = all_hashes + coords[win].sketch_offset;
      size_t sz = coords[win].sketch_size;
      for (size_t k = 0; k < sz; k++) {
        uint32_t val = h[k];
        size_t p = (size_t)(val / part_size);
        if (p >= NUM_PARTITIONS)
          p = NUM_PARTITIONS - 1;
        if (p >= batch_start && p < batch_end) {
          DA_PUSH(w.buckets[p].entries, w.buckets[p].size, w.buckets[p].cap,
                  ((HashWindowEntry){val, (uint32_t)win}));
        }
      }
    }

    kt_for(n_threads, discover_compute_worker, &w, (long)batch_count);

    for (size_t p = batch_start; p < batch_end; p++) {
      free(w.buckets[p].entries);
      w.buckets[p].entries = NULL;
      w.buckets[p].size = 0;
      w.buckets[p].cap = 0;
    }
  }

  for (int t = 0; t < n_threads; t++)
    free(w.t_bloom[t]);
  free(w.buckets);
  free(w.t_bloom);
}

SegtraceDistResult calculate_window_dist(const uint32_t *all_hashes,
                                         const WindowCoord *wa,
                                         const WindowCoord *wb) {
  SegtraceSketch sa = {.sketch_size = wa->sketch_size,
                       .hashes = (uint32_t *)(all_hashes + wa->sketch_offset)};
  SegtraceSketch sb = {.sketch_size = wb->sketch_size,
                       .hashes = (uint32_t *)(all_hashes + wb->sketch_offset)};
  return calculate_segtrace_dist(&sa, &sb);
}

static inline int check_collinear_neighbor(DiscoverComputeData *w, uint32_t wa,
                                           uint32_t wb, size_t min_shared) {
  const int dir_a[] = {1, -1, 1, -1};
  const int dir_b[] = {1, -1, -1, 1};

  uint32_t seq_a = w->coords[wa].seq_id;
  uint32_t seq_b = w->coords[wb].seq_id;

  for (int d = 0; d < 4; d++) {
    for (int step_a = 1; step_a <= MAX_COLLINEAR_LOOOKAHEAD; step_a++) {
      for (int step_b = 1; step_b <= MAX_COLLINEAR_LOOOKAHEAD; step_b++) {

        long long next_a = (long long)wa + (long long)dir_a[d] * step_a;
        long long next_b = (long long)wb + (long long)dir_b[d] * step_b;

        if (next_a >= 0 && next_a < (long long)w->n_windows && next_b >= 0 &&
            next_b < (long long)w->n_windows &&
            w->coords[next_a].seq_id == seq_a &&
            w->coords[next_b].seq_id == seq_b) {

          if (w->coords[next_a].sketch_size > 0 &&
              w->coords[next_b].sketch_size > 0) {
            SegtraceDistResult res = calculate_window_dist(
                w->all_hashes, &w->coords[next_a], &w->coords[next_b]);
            if (res.shared_hashes >= min_shared)
              return 1;
          }
        }
      }
    }
  }
  return 0;
}

void discover_compute_worker(void *data, long idx, int tid) {
  DiscoverComputeData *w_data = (DiscoverComputeData *)data;
  long p = (long)w_data->batch_start + idx;
  PartitionBucket *b = &w_data->buckets[p];
  if (b->size == 0)
    return;

  qsort(b->entries, b->size, sizeof(HashWindowEntry), compare_hash_entry);
  memset(w_data->t_bloom[tid], 0, BLOOM_SIZE_BYTES);

  double p_kmer = pow(0.90, (double)w_data->kmer_size);
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
          size_t start_wa =
              (size_t)w_data->coords[wa].window_idx * w_data->step_size;
          size_t start_wb =
              (size_t)w_data->coords[wb].window_idx * w_data->step_size;

          if (w_data->coords[wa].seq_id == w_data->coords[wb].seq_id &&
              ABS_DIFF(start_wa, start_wb) < w_data->window_size)
            continue;

          uint64_t pk = encode_pair(wa, wb);
          if (bloom_test_and_set(w_data->t_bloom[tid], pk, BLOOM_MASK))
            continue;

          size_t min_sz =
              w_data->coords[wa].sketch_size < w_data->coords[wb].sketch_size
                  ? w_data->coords[wa].sketch_size
                  : w_data->coords[wb].sketch_size;
          size_t min_shared = (size_t)ceil((double)min_sz * p_kmer) * 2;
          if (min_shared < 2)
            min_shared = 2;

          SegtraceDistResult d = calculate_window_dist(
              w_data->all_hashes, &w_data->coords[wa], &w_data->coords[wb]);
          if (d.shared_hashes >= min_shared) {
            if (check_collinear_neighbor(w_data, wa, wb, min_shared)) {
              union_unionfind(w_data->uf, wa, wb);
            }
          }
        }
      }
    }
    i = j;
  }
}

// ==============================================================
// SECTION 4: CLUSTERING, LOCUS MERGING & FLANKING SUBCLUSTERING
// ==============================================================

void build_duplicate_regions(UnionFind *uf, size_t num_sketches,
                             WindowCoord *coords, size_t step_size,
                             size_t window_size,
                             SegtraceDupRegion **out_regions,
                             size_t *out_n_regions) {

  uint32_t *comp_size = calloc(num_sketches, sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++) {
    comp_size[find_unionfind(uf, (uint32_t)i)]++;
  }

  uint32_t *cluster_map = calloc(num_sketches, sizeof(uint32_t));
  uint32_t next_cluster_id = 1;
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t root = find_unionfind(uf, (uint32_t)i);
    if (comp_size[root] >= 2 && cluster_map[root] == 0) {
      cluster_map[root] = next_cluster_id++;
    }
  }

  size_t n_dup_regions = 0, cap_dup_regions = 0;
  SegtraceDupRegion *dup_regions = NULL;

  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t root_i = find_unionfind(uf, (uint32_t)i);
    uint32_t cid = cluster_map[root_i];
    if (cid == 0) {
      continue;
    }

    uint32_t seq_i = coords[i].seq_id;
    size_t start = (size_t)coords[i].window_idx * step_size;
    size_t end = start + window_size;

    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((SegtraceDupRegion){.seq_id = seq_i,
                                 .start = start,
                                 .end = end,
                                 .cluster_id = cid,
                                 .subcluster_id = 0,
                                 .flank_sketch = {0},
                                 .window_idx = coords[i].window_idx}));
  }

  free(comp_size);
  free(cluster_map);

  *out_regions = dup_regions;
  *out_n_regions = n_dup_regions;
}

static int compare_dup_region_by_pos(const void *a, const void *b) {
  const SegtraceDupRegion *ra = (const SegtraceDupRegion *)a,
                          *rb = (const SegtraceDupRegion *)b;
  if (ra->seq_id != rb->seq_id)
    return CMP(ra->seq_id, rb->seq_id);
  if (ra->start != rb->start)
    return CMP(ra->start, rb->start);
  return CMP(ra->end, rb->end);
}

static int compare_dup_region_by_cluster(const void *a, const void *b) {
  const SegtraceDupRegion *ra = (const SegtraceDupRegion *)a,
                          *rb = (const SegtraceDupRegion *)b;
  if (ra->cluster_id != rb->cluster_id)
    return CMP(ra->cluster_id, rb->cluster_id);
  if (ra->seq_id != rb->seq_id)
    return CMP(ra->seq_id, rb->seq_id);
  if (ra->start != rb->start)
    return CMP(ra->start, rb->start);
  return CMP(ra->end, rb->end);
}

size_t merge_dup_regions(SegtraceDupRegion *regions, size_t n,
                         size_t window_size) {
  if (n <= 1)
    return n;
  qsort(regions, n, sizeof(SegtraceDupRegion), compare_dup_region_by_cluster);

  size_t out = 0;
  for (size_t i = 1; i < n; i++) {
    if (regions[i].cluster_id == regions[out].cluster_id &&
        regions[i].seq_id == regions[out].seq_id &&
        regions[i].start <= regions[out].end + MERGE_COEFF * window_size) {
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
      if (regions[i].window_idx > regions[out].window_idx)
        regions[out].window_idx = regions[i].window_idx;
    } else {
      out++;
      if (out != i)
        regions[out] = regions[i];
    }
  }
  return out + 1;
}

void extract_flankings(char **files, int num_files, const Segtrace *r,
                       uint64_t scale, SegtraceDupRegion *regions,
                       size_t n_regions, int n_threads, size_t flank_size,
                       const GenomeSeqLen *seq_lens, size_t num_seqs) {
  qsort(regions, n_regions, sizeof(SegtraceDupRegion),
        compare_dup_region_by_pos);
  FlankingWorkerData w = {files,     r,          scale,    regions,
                          n_regions, flank_size, seq_lens, num_seqs};
  kt_for(n_threads, extract_flankings_worker, &w, num_files);
}

static void find_seq_id_range(const SegtraceDupRegion *regions, size_t n,
                              uint32_t seq_id, size_t *first, size_t *last) {
  size_t low = 0, high = n;
  while (low < high) {
    size_t mid = low + (high - low) / 2;
    if (regions[mid].seq_id < seq_id)
      low = mid + 1;
    else
      high = mid;
  }
  *first = low;

  low = 0;
  high = n;
  while (low < high) {
    size_t mid = low + (high - low) / 2;
    if (regions[mid].seq_id <= seq_id)
      low = mid + 1;
    else
      high = mid;
  }
  *last = low;
}

void extract_flankings_worker(void *data, long f, int tid) {
  (void)tid;
  FlankingWorkerData *w = (FlankingWorkerData *)data;

  char bname[256];
  get_basename(w->files[f], bname, sizeof(bname));
  gzFile fp = gzopen(w->files[f], "r");
  if (!fp)
    return;
  kseq_t *ks = kseq_init(fp);
  if (!ks) {
    gzclose(fp);
    return;
  }

  while (kseq_read(ks) >= 0) {
    int found_seq = 0;
    uint32_t target_seq_id = 0;
    for (size_t s = 0; s < w->num_seqs; s++) {
      if (strcmp(w->seq_lens[s].genome, bname) == 0 &&
          strcmp(w->seq_lens[s].seq, ks->name.s) == 0) {
        target_seq_id = (uint32_t)s;
        found_seq = 1;
        break;
      }
    }
    if (!found_seq)
      continue;

    size_t first, last;
    find_seq_id_range(w->regions, w->n_regions, target_seq_id, &first, &last);

    for (size_t i = first; i < last; i++) {
      size_t start = w->regions[i].start, end = w->regions[i].end,
             flank_size = w->flank_size;
      size_t seq_len = ks->seq.l;
      size_t start_clamped = start > seq_len ? seq_len : start;
      size_t end_clamped = end > seq_len ? seq_len : end;
      if (end_clamped < start_clamped)
        end_clamped = start_clamped;

      size_t left_start =
          start_clamped > flank_size ? start_clamped - flank_size : 0;
      size_t right_end = end_clamped + flank_size > seq_len
                             ? seq_len
                             : end_clamped + flank_size;
      size_t left_len = start_clamped - left_start;
      size_t right_len = right_end > end_clamped ? right_end - end_clamped : 0;

      free(w->regions[i].flank_sketch.hashes);
      w->regions[i].flank_sketch.hashes = NULL;
      w->regions[i].flank_sketch.sketch_size = 0;

      if (left_len + right_len == 0)
        continue;

      HashPool pool;
      init_hash_pool(&pool, UINT32_MAX / w->scale);
      if (left_len > 0)
        extract_hash(w->r, &pool, (const uint8_t *)(ks->seq.s + left_start),
                     left_len);
      if (right_len > 0)
        extract_hash(w->r, &pool, (const uint8_t *)(ks->seq.s + end_clamped),
                     right_len);
      finalize_hash_pool(&pool, &w->regions[i].flank_sketch.hashes,
                         &w->regions[i].flank_sketch.sketch_size);
    }
  }
  kseq_destroy(ks);
  gzclose(fp);
}

void perform_subclustering(SegtraceDupRegion *regions, size_t n_merged,
                           int n_threads, uint32_t kmer_size) {
  if (n_merged <= 1)
    return;

  qsort(regions, n_merged, sizeof(SegtraceDupRegion),
        compare_dup_region_by_cluster);

  ClusterSpan *spans = NULL;
  size_t n_spans = 0, cap_spans = 0;

  size_t i = 0;
  while (i < n_merged) {
    size_t j = i + 1;
    while (j < n_merged && regions[i].cluster_id == regions[j].cluster_id)
      j++;
    DA_PUSH(spans, n_spans, cap_spans, ((ClusterSpan){i, j - i}));
    i = j;
  }

  UnionFind sub_uf;
  init_unionfind(&sub_uf, n_merged);

  SubclusterData w = {.regions = regions,
                      .n_merged = n_merged,
                      .kmer_size = kmer_size,
                      .spans = spans,
                      .t_pairs = calloc(n_threads, sizeof(SubclusterPair *)),
                      .t_n_pairs = calloc(n_threads, sizeof(size_t)),
                      .t_cap_pairs = calloc(n_threads, sizeof(size_t)),
                      .t_bloom = calloc(n_threads, sizeof(uint8_t *))};
  for (int t = 0; t < n_threads; t++)
    w.t_bloom[t] = calloc(SUBCLUSTER_BLOOM_SIZE_BYTES, 1);

  kt_for(n_threads, process_subcluster, &w, n_spans);

  for (int t = 0; t < n_threads; t++) {
    for (size_t k = 0; k < w.t_n_pairs[t]; k++) {
      union_unionfind(&sub_uf, w.t_pairs[t][k].i, w.t_pairs[t][k].j);
    }
    if (w.t_pairs[t])
      free(w.t_pairs[t]);
    free(w.t_bloom[t]);
  }
  free(w.t_bloom);
  free(w.t_pairs);
  free(w.t_n_pairs);
  free(w.t_cap_pairs);
  free(spans);

  uint32_t *mapping = calloc(n_merged, sizeof(uint32_t));
  uint32_t current_id = 1;

  for (size_t k = 0; k < n_merged; k++) {
    uint32_t p = find_unionfind(&sub_uf, (uint32_t)k);
    if (mapping[p] == 0)
      mapping[p] = current_id++;
    regions[k].subcluster_id = mapping[p];
  }
  free(mapping);
  free_unionfind(&sub_uf);
}

static inline void check_and_eval_flank_pair(SubclusterData *w, int tid,
                                             size_t ra, size_t rb,
                                             double p_kmer) {
  if (w->regions[ra].flank_sketch.sketch_size == 0 ||
      w->regions[rb].flank_sketch.sketch_size == 0)
    return;

  size_t min_sz = w->regions[ra].flank_sketch.sketch_size <
                          w->regions[rb].flank_sketch.sketch_size
                      ? w->regions[ra].flank_sketch.sketch_size
                      : w->regions[rb].flank_sketch.sketch_size;
  size_t min_shared = (size_t)ceil((double)min_sz * p_kmer) * 2;
  if (min_shared < 2)
    min_shared = 2;

  SegtraceDistResult d = calculate_segtrace_dist(&w->regions[ra].flank_sketch,
                                                 &w->regions[rb].flank_sketch);
  if (d.shared_hashes >= min_shared) {
    DA_PUSH(w->t_pairs[tid], w->t_n_pairs[tid], w->t_cap_pairs[tid],
            ((SubclusterPair){(uint32_t)ra, (uint32_t)rb}));
  }
}

void process_subcluster(void *data, long s, int tid) {
  SubclusterData *w = (SubclusterData *)data;
  size_t start = w->spans[s].start, count = w->spans[s].count;
  if (count <= 1)
    return;

  double p_kmer = pow(0.90, (double)w->kmer_size);

  size_t total_flank_hashes = 0;
  for (size_t a = 0; a < count; a++)
    total_flank_hashes += w->regions[start + a].flank_sketch.sketch_size;
  if (total_flank_hashes == 0)
    return;

  HashWindowEntry *entries =
      malloc(total_flank_hashes * sizeof(HashWindowEntry));
  if (!entries)
    return;

  size_t n_entries = 0;
  for (size_t a = 0; a < count; a++) {
    size_t ra = start + a;
    const SegtraceSketch *sk = &w->regions[ra].flank_sketch;
    for (size_t k = 0; k < sk->sketch_size; k++) {
      entries[n_entries++] = (HashWindowEntry){sk->hashes[k], (uint32_t)a};
    }
  }

  qsort(entries, n_entries, sizeof(HashWindowEntry), compare_hash_entry);
  uint8_t *bloom = w->t_bloom[tid];
  memset(bloom, 0, SUBCLUSTER_BLOOM_SIZE_BYTES);

  size_t i = 0;
  while (i < n_entries) {
    size_t j = i + 1;
    while (j < n_entries && entries[j].hash == entries[i].hash)
      j++;
    size_t run_len = j - i;
    if (run_len >= 2) {
      for (size_t a = i; a < j; a++) {
        size_t b_max =
            a + 1 + MAX_PAIR_COMPARISONS < j ? a + 1 + MAX_PAIR_COMPARISONS : j;
        for (size_t b = a + 1; b < b_max; b++) {
          uint32_t la = entries[a].window_id, lb = entries[b].window_id;
          if (la == lb)
            continue;
          uint64_t pk = encode_pair(la, lb);
          if (bloom_test_and_set(bloom, pk, SUBCLUSTER_BLOOM_MASK))
            continue;
          check_and_eval_flank_pair(w, tid, start + la, start + lb, p_kmer);
        }
      }
    }
    i = j;
  }
  free(entries);
}

// ==============================================================
// SECTION 5: REPORTING & FILE OUTPUT WRITERS
// ==============================================================

void write_dup_bed(const char *out_prefix, SegtraceDupRegion *dup_regions,
                   size_t n_merged, const GenomeSeqLen *seq_lens) {
  if (n_merged == 0)
    return;
  qsort(dup_regions, n_merged, sizeof(SegtraceDupRegion),
        compare_dup_region_by_cluster);

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
  FILE *out_bed = fopen(path_buf, "w");
  if (!out_bed) {
    fprintf(stderr, "[ERROR] Cannot open output file: %s\n", path_buf);
    return;
  }

  fprintf(out_bed, "#chrom\tstart\tend\tcluster_id\tsubcluster_id\n");
  for (size_t k = 0; k < n_merged; k++) {
    if (dup_regions[k].end - dup_regions[k].start >= MIN_SD_LEN) {
      uint32_t seq_i = dup_regions[k].seq_id;
      fprintf(out_bed, "%s-%s\t%zu\t%zu\t%u\t%u\n", seq_lens[seq_i].genome,
              seq_lens[seq_i].seq, dup_regions[k].start, dup_regions[k].end,
              dup_regions[k].cluster_id, dup_regions[k].subcluster_id);
    }
  }
  fclose(out_bed);
}

// ==============================================================
// SECTION 6: CORE ALGORITHMS & UTILITIES
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

void init_segtrace(Segtrace *r, size_t hash_window, int filter_masked) {
  r->hash_window = (uint32_t)hash_window;
  r->filter_masked = filter_masked;
  r->base_lookup = filter_masked ? BASE_LOOKUP : BASE_LOOKUP_NO_MASK;
}

__attribute__((hot)) void extract_hash(const Segtrace *r, HashPool *pool,
                                       const uint8_t *seq, size_t len) {
  uint32_t k = r->hash_window;
  if (len < k)
    return;

  uint64_t f_hash = 0, r_hash = 0;
  size_t valid_len = 0;

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
      if (b_out < 0 || b_out > 3) {
        valid_len = 0;
        f_hash = 0;
        r_hash = 0;
        continue;
      }
      int8_t b_in = b;
      f_hash = rol64(f_hash, 1) ^ rol64(NTHASH_H[b_out], k) ^ NTHASH_H[b_in];
      r_hash = ror64(r_hash, 1) ^ ror64(NTHASH_H[b_out ^ 3], 1) ^
               rol64(NTHASH_H[b_in ^ 3], k - 1);
    }

    if (valid_len >= k) {
      uint64_t canonical = (f_hash < r_hash) ? f_hash : r_hash;
      insert_hash_pool(pool, mix_hash(canonical, r->hash_seed));
    }
  }
}

void init_hash_pool(HashPool *pool, uint32_t threshold) {
  pool->size = 0;
  pool->cap = 16;
  pool->hash_threshold = threshold;
  pool->hashes = malloc(pool->cap * sizeof(uint32_t));
}

void insert_hash_pool(HashPool *pool, uint32_t h) {
  if (h >= pool->hash_threshold)
    return;
  DA_PUSH(pool->hashes, pool->size, pool->cap, h);
}

void finalize_hash_pool(HashPool *pool, uint32_t **out_hashes,
                        size_t *out_size) {
  if (pool->size == 0) {
    if (pool->hashes)
      free(pool->hashes);
    *out_hashes = NULL;
    *out_size = 0;
    return;
  }
  qsort(pool->hashes, pool->size, sizeof(uint32_t), compare_uint32);
  size_t unique_count = 0;
  for (size_t i = 0; i < pool->size; i++) {
    if (i == 0 || pool->hashes[i] != pool->hashes[i - 1])
      pool->hashes[unique_count++] = pool->hashes[i];
  }
  *out_hashes = realloc(pool->hashes, unique_count * sizeof(uint32_t));
  *out_size = unique_count;
}

SegtraceDistResult calculate_segtrace_dist(const SegtraceSketch *ref,
                                           const SegtraceSketch *query) {
  SegtraceDistResult res = {0};
  if (!ref || !query || ref->sketch_size == 0 || query->sketch_size == 0)
    return res;

  size_t i = 0, j = 0, shared = 0;
  while (i < ref->sketch_size && j < query->sketch_size) {
    if (ref->hashes[i] == query->hashes[j]) {
      shared++;
      i++;
      j++;
    } else if (ref->hashes[i] < query->hashes[j]) {
      i++;
    } else {
      j++;
    }
  }
  res.shared_hashes = shared;
  return res;
}

void init_unionfind(UnionFind *uf, size_t n) {
  uf->n = n;
  uf->parent = malloc(n * sizeof(uint32_t));
  uf->rank = calloc(n, sizeof(uint8_t));
  pthread_mutex_init(&uf->lock, NULL);
  for (size_t i = 0; i < n; i++)
    uf->parent[i] = (uint32_t)i;
}

uint32_t find_unionfind(UnionFind *uf, uint32_t x) {
  /* Iterative path-splitting to avoid stack overflow on deep chains */
  while (uf->parent[x] != x) {
    uint32_t next = uf->parent[x];
    uf->parent[x] = uf->parent[next]; /* path splitting */
    x = next;
  }
  return x;
}

void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b) {
  pthread_mutex_lock(&uf->lock);
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
  pthread_mutex_unlock(&uf->lock);
}

void free_unionfind(UnionFind *uf) {
  if (uf->parent)
    free(uf->parent);
  if (uf->rank)
    free(uf->rank);
  pthread_mutex_destroy(&uf->lock);
  uf->parent = NULL;
  uf->rank = NULL;
  uf->n = 0;
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

int bloom_test_and_set(uint8_t *bloom, uint64_t key, uint32_t mask) {
  uint64_t h = splitmix64(key);
  uint32_t h1 = (uint32_t)h & mask, h2 = (uint32_t)(h >> 32) & mask;
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
