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

/* 2-bit nucleotide encoding: A = 00, C = 01, G = 10, T = 11 */
const int8_t BASE_LOOKUP[256] = {
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
         "  -k: kmer size (default: 21)\n"
         "  -s: scale factor (default: 16)\n"
         "  -e: hash seed (default: 42)\n"
         "  -w: window size in bp (default: 1024)\n"
         "  -t: step size in bp (default: 0 [auto: 50%% of window size])\n"
         "  -b: minimum valid bases per window (default: 1000)\n"
         "  -d: maximum distance to consider as copy (default: 0.10)\n"
         "  -D: sub-cluster distance threshold (default: 0.3)\n"
         "  -f: flanking size in bp for sub-clustering (default: 1000)\n"
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

  uint32_t def_kmer_size = 25;
  uint64_t def_scale = 16, def_hash_seed = 42;
  size_t window_size = 1024, step_size = 0, min_bases = 1000, flank_size = 1000;
  const char *out_prefix = "segtrace";
  int n_threads = 8;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:e:w:t:b:d:o:p:D:f:h", 0)) >= 0) {
    if (c == 'h') {
      print_usage();
      return 0;
    } else if (c == 'k')
      def_kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      def_scale = (uint64_t)strtoull(opt.arg, NULL, 10);
    else if (c == 'e')
      def_hash_seed = strtoull(opt.arg, NULL, 0);
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
    else if (c == 'f')
      flank_size = (size_t)strtoull(opt.arg, NULL, 10);
    else
      return 1;
  }

  if (step_size == 0)
    step_size = window_size / 3;
  if (opt.ind == argc) {
    fprintf(stderr, "[ERROR] Input FASTA files are required.\n");
    return 1;
  }

  int num_files = argc - opt.ind;
  char **files = &argv[opt.ind];

  Segtrace r;
  init_segtrace(&r, def_kmer_size);
  r.hash_seed = def_hash_seed;

  uint64_t *all_hashes = NULL;
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
  discover_and_compute(all_hashes, coords, num_sketches, window_size, n_threads,
                       r.hash_window, &uf);

  SegtraceDupRegion *dup_regions = NULL;
  size_t n_dup_regions = 0;
  build_duplicate_regions(&uf, num_sketches, window_size, num_files, files,
                          seq_lens, coords, &dup_regions, &n_dup_regions);

  free_unionfind(&uf);
  free(all_hashes);
  all_hashes = NULL;
  free(coords);
  coords = NULL;

  size_t n_merged = merge_dup_regions(dup_regions, n_dup_regions);

  fprintf(stderr,
          "[INFO] Extracting flanking sequences for sub-clustering...\n");
  extract_flankings(files, num_files, &r, def_scale, dup_regions, n_merged,
                    n_threads, flank_size);

  fprintf(stderr, "[INFO] Sub-clustering based on flanking similarities...\n");
  perform_subclustering(dup_regions, n_merged, n_threads, r.hash_window);

  write_dup_bed(out_prefix, dup_regions, n_merged);
  write_dup_bedpe(out_prefix, dup_regions, n_merged);

  for (size_t i = 0; i < n_merged; i++) {
    free(dup_regions[i].chrom);
    free(dup_regions[i].cluster_id);
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
       idx + job->window_size <= job->chunk_end_idx;
       idx += job->step_size, current_window_idx++) {
    size_t valid_bases = 0;
    for (size_t j = 0; j < job->window_size; j++) {
      if (BASE_LOOKUP[job->seq_ptr[idx + j]] >= 0)
        valid_bases++;
    }
    if (valid_bases < job->min_bases)
      continue;

    HashPool pool;
    init_hash_pool(&pool, UINT64_MAX / job->scale);
    extract_hash(job->r, &pool, job->seq_ptr + idx, job->window_size);

    uint64_t *hashes = NULL;
    size_t sketch_size = 0;
    finalize_hash_pool(&pool, &hashes, &sketch_size);

    if (sketch_size > 0) {
      DA_RESERVE(job->coords, job->cap_coords, job->num_coords + 1);
      WindowCoord *wc = &job->coords[job->num_coords++];
      wc->seq_id = job->seq_id;
      wc->start = idx;
      wc->end = idx + job->window_size;
      wc->window_idx = current_window_idx;

      size_t h_idx = job->num_hashes;
      DA_RESERVE(job->hashes, job->cap_hashes, job->num_hashes + sketch_size);
      memcpy(job->hashes + h_idx, hashes, sketch_size * sizeof(uint64_t));
      job->num_hashes += sketch_size;

      wc->sketch_offset = h_idx;
      wc->sketch_size = (uint32_t)sketch_size;
      free(hashes);
    }
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

  size_t cap_jobs = 16, num_jobs = 0;
  SeqChunkJob *jobs = malloc(cap_jobs * sizeof(SeqChunkJob));

  for (int f = 0; f < num_files; f++) {
    char bname[256];
    get_basename(files[f], bname, sizeof(bname));

    gzFile fp = gzopen(files[f], "r");
    if (!fp)
      continue;
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

      for (size_t c_start = 0; c_start < len; c_start += chunk_size) {
        size_t c_end = c_start + chunk_size;
        if (c_end > len)
          c_end = len;
        if (c_end - c_start < window_size)
          break;

        DA_RESERVE(jobs, cap_jobs, num_jobs + 1);
        jobs[num_jobs++] = (SeqChunkJob){.r = r,
                                         .scale = scale,
                                         .window_size = window_size,
                                         .step_size = step_size,
                                         .min_bases = min_bases,
                                         .seq_id = seq_id,
                                         .seq_ptr = seq_copy,
                                         .seq_len = len,
                                         .chunk_start_idx = c_start,
                                         .chunk_end_idx = c_end};
      }
    }
    kseq_destroy(ks);
    gzclose(fp);
  }

  kt_for(n_threads, seq_chunk_worker, jobs, num_jobs);

  for (size_t j = 0; j < num_jobs; j++) {
    SeqChunkJob *job = &jobs[j];
    if (job->num_coords > 0) {
      size_t base_h_offset = workers[0].num_all_hashes;
      DA_RESERVE(workers[0].all_hashes, workers[0].cap_all_hashes,
                 workers[0].num_all_hashes + job->num_hashes);
      memcpy(workers[0].all_hashes + base_h_offset, job->hashes,
             job->num_hashes * sizeof(uint64_t));
      workers[0].num_all_hashes += job->num_hashes;

      size_t base_c_offset = workers[0].num_sketches;
      DA_RESERVE(workers[0].coords, workers[0].cap_sketches,
                 workers[0].num_sketches + job->num_coords);

      for (size_t k = 0; k < job->num_coords; k++) {
        WindowCoord wc = job->coords[k];
        wc.sketch_offset += base_h_offset;
        workers[0].coords[base_c_offset + k] = wc;
      }
      workers[0].num_sketches += job->num_coords;
    }
    free(job->hashes);
    free(job->coords);
  }

  uint8_t **freed_seqs = NULL;
  size_t n_freed = 0, cap_freed = 0;
  for (size_t j = 0; j < num_jobs; j++) {
    uint8_t *p = (uint8_t *)jobs[j].seq_ptr;
    int already = 0;
    for (size_t k = 0; k < n_freed; k++) {
      if (freed_seqs[k] == p) {
        already = 1;
        break;
      }
    }
    if (!already) {
      DA_RESERVE(freed_seqs, cap_freed, n_freed + 1);
      freed_seqs[n_freed++] = p;
      free(p);
    }
  }
  free(freed_seqs);
  free(jobs);
  return workers;
}

void merge_global_data(StreamWorkerData *workers, int num_files,
                       const char *out_prefix, uint64_t **out_all_hashes,
                       WindowCoord **out_coords, size_t *out_num_sketches,
                       GenomeSeqLen **out_seq_lens, size_t *out_num_seqs) {
  size_t total_hashes = 0, total_sketches = 0, total_seqs = 0;
  for (int i = 0; i < num_files; i++) {
    total_hashes += workers[i].num_all_hashes;
    total_sketches += workers[i].num_sketches;
    total_seqs += workers[i].num_seqs;
  }

  uint64_t *all_hashes =
      total_hashes ? malloc(total_hashes * sizeof(uint64_t)) : NULL;
  WindowCoord *coords =
      total_sketches ? malloc(total_sketches * sizeof(WindowCoord)) : NULL;
  GenomeSeqLen *seq_lens =
      total_seqs ? malloc(total_seqs * sizeof(GenomeSeqLen)) : NULL;
  size_t g_hash_offset = 0, g_sketch_offset = 0, g_seq_offset = 0;

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.window.bed", out_prefix);
  FILE *bed_fp = fopen(path_buf, "w");

  for (int i = 0; i < num_files; i++) {
    StreamWorkerData *w = &workers[i];
    if (w->num_all_hashes > 0)
      memcpy(all_hashes + g_hash_offset, w->all_hashes,
             w->num_all_hashes * sizeof(uint64_t));
    if (w->num_seqs > 0)
      memcpy(seq_lens + g_seq_offset, w->seq_lens,
             w->num_seqs * sizeof(GenomeSeqLen));

    for (size_t j = 0; j < w->num_sketches; j++) {
      WindowCoord c = w->coords[j];
      c.sketch_offset += g_hash_offset;
      c.seq_id += g_seq_offset;
      coords[g_sketch_offset + j] = c;

      if (bed_fp) {
        char chr_name[512];
        snprintf(chr_name, sizeof(chr_name), "%s-%s", seq_lens[c.seq_id].genome,
                 seq_lens[c.seq_id].seq);
        fprintf(bed_fp, "%s\t%zu\t%zu\t%s_%zu_%zu\n", chr_name, c.start, c.end,
                chr_name, c.start, c.end);
      }
    }

    g_hash_offset += w->num_all_hashes;
    g_sketch_offset += w->num_sketches;
    g_seq_offset += w->num_seqs;
    free(w->all_hashes);
    free(w->coords);
    free(w->seq_lens);
  }
  if (bed_fp)
    fclose(bed_fp);

  *out_all_hashes = all_hashes;
  *out_coords = coords;
  *out_num_sketches = total_sketches;
  *out_seq_lens = seq_lens;
  *out_num_seqs = total_seqs;
}

// ==============================================================
// SECTION 3: CANDIDATE DISCOVERY & DISTANCE COMPUTATION
// ==============================================================

void discover_and_compute(const uint64_t *all_hashes, WindowCoord *coords,
                          size_t n_windows, size_t window_size, int n_threads,
                          uint32_t kmer_size, UnionFind *uf) {
  DiscoverComputeData w = {.all_hashes = all_hashes,
                           .coords = coords,
                           .n_windows = n_windows,
                           .window_size = window_size,
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

  uint64_t part_size = UINT64_MAX / NUM_PARTITIONS;
  for (size_t win = 0; win < n_windows; win++) {
    const uint64_t *h = all_hashes + coords[win].sketch_offset;
    size_t sz = coords[win].sketch_size;
    for (size_t k = 0; k < sz; k++) {
      uint64_t val = h[k];
      size_t p = (size_t)(val / part_size);
      if (p >= NUM_PARTITIONS)
        p = NUM_PARTITIONS - 1;
      DA_PUSH(w.buckets[p].entries, w.buckets[p].size, w.buckets[p].cap,
              ((HashWindowEntry){val, (uint32_t)win}));
    }
  }

  kt_for(n_threads, discover_compute_worker, &w, NUM_PARTITIONS);

  for (int t = 0; t < n_threads; t++)
    free(w.t_bloom[t]);
  for (size_t p = 0; p < NUM_PARTITIONS; p++)
    free(w.buckets[p].entries);
  free(w.buckets);
  free(w.t_bloom);
}

SegtraceDistResult calculate_window_dist(const uint64_t *all_hashes,
                                         const WindowCoord *wa,
                                         const WindowCoord *wb,
                                         uint32_t kmer_size) {
  SegtraceSketch sa = {.sketch_size = wa->sketch_size,
                       .hashes = (uint64_t *)(all_hashes + wa->sketch_offset)};
  SegtraceSketch sb = {.sketch_size = wb->sketch_size,
                       .hashes = (uint64_t *)(all_hashes + wb->sketch_offset)};
  return calculate_segtrace_dist(&sa, &sb, kmer_size);
}

static inline int check_collinear_neighbor(DiscoverComputeData *w, uint32_t wa,
                                           uint32_t wb, size_t min_shared) {
  // 1. Forward collinear neighbor: (wa + 1, wb + 1)
  if (wa + 1 < w->n_windows && wb + 1 < w->n_windows &&
      w->coords[wa + 1].seq_id == w->coords[wa].seq_id &&
      w->coords[wb + 1].seq_id == w->coords[wb].seq_id) {
    if (calculate_window_dist(w->all_hashes, &w->coords[wa + 1],
                              &w->coords[wb + 1], w->kmer_size)
            .shared_hashes >= min_shared)
      return 1;
  }
  // 2. Inverted collinear neighbor: (wa + 1, wb - 1)
  if (wa + 1 < w->n_windows && wb > 0 &&
      w->coords[wa + 1].seq_id == w->coords[wa].seq_id &&
      w->coords[wb - 1].seq_id == w->coords[wb].seq_id) {
    if (calculate_window_dist(w->all_hashes, &w->coords[wa + 1],
                              &w->coords[wb - 1], w->kmer_size)
            .shared_hashes >= min_shared)
      return 1;
  }
  // 3. Backward collinear neighbor: (wa - 1, wb - 1)
  if (wa > 0 && wb > 0 && w->coords[wa - 1].seq_id == w->coords[wa].seq_id &&
      w->coords[wb - 1].seq_id == w->coords[wb].seq_id) {
    if (calculate_window_dist(w->all_hashes, &w->coords[wa - 1],
                              &w->coords[wb - 1], w->kmer_size)
            .shared_hashes >= min_shared)
      return 1;
  }
  return 0;
}

void discover_compute_worker(void *data, long p, int tid) {
  DiscoverComputeData *w_data = (DiscoverComputeData *)data;
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

    if (run_len >= 2 && run_len <= MAX_RUN_LEN) {
      for (size_t a = i; a < j; a++) {
        for (size_t b_idx = a + 1; b_idx < j; b_idx++) {
          uint32_t wa = b->entries[a].window_id,
                   wb = b->entries[b_idx].window_id;
          if (w_data->coords[wa].seq_id == w_data->coords[wb].seq_id &&
              ABS_DIFF(w_data->coords[wa].start, w_data->coords[wb].start) <
                  w_data->window_size)
            continue;

          uint64_t pk = encode_pair(wa, wb);
          if (bloom_test_and_set(w_data->t_bloom[tid], pk, BLOOM_MASK))
            continue;

          size_t min_sz =
              w_data->coords[wa].sketch_size < w_data->coords[wb].sketch_size
                  ? w_data->coords[wa].sketch_size
                  : w_data->coords[wb].sketch_size;
          size_t min_shared = (size_t)floor((double)min_sz * p_kmer);
          if (min_shared < 1)
            min_shared = 1;

          SegtraceDistResult d =
              calculate_window_dist(w_data->all_hashes, &w_data->coords[wa],
                                    &w_data->coords[wb], w_data->kmer_size);
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
                             size_t window_size, int num_files, char **files,
                             GenomeSeqLen *seq_lens, WindowCoord *coords,
                             SegtraceDupRegion **out_regions,
                             size_t *out_n_regions) {
  (void)num_files;
  (void)files;

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

  size_t i = 0;
  while (i < num_sketches) {
    uint32_t root_i = find_unionfind(uf, (uint32_t)i);
    if (cluster_map[root_i] == 0) {
      i++;
      continue;
    }

    uint32_t seq_i = coords[i].seq_id;
    size_t min_start = coords[i].start;
    size_t max_end = coords[i].end;
    uint32_t min_cid = cluster_map[root_i];
    uint32_t win_idx = coords[i].window_idx;

    size_t j = i + 1;
    while (j < num_sketches && coords[j].seq_id == seq_i &&
           coords[j].start <= max_end + window_size) {
      uint32_t root_j = find_unionfind(uf, (uint32_t)j);
      if (cluster_map[root_j] != 0) {
        if (coords[j].end > max_end)
          max_end = coords[j].end;
        if (cluster_map[root_j] < min_cid)
          min_cid = cluster_map[root_j];
      }
      j++;
    }

    char label[32];
    snprintf(label, sizeof(label), "%u", min_cid);

    char chrom_name[512];
    snprintf(chrom_name, sizeof(chrom_name), "%s-%s", seq_lens[seq_i].genome,
             seq_lens[seq_i].seq);

    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((SegtraceDupRegion){.chrom = strdup(chrom_name),
                                 .start = min_start,
                                 .end = max_end,
                                 .cluster_id = strdup(label),
                                 .copy_count = comp_size[root_i],
                                 .subcluster_id = 0,
                                 .flank_sketch = {0},
                                 .window_idx = win_idx}));
    i = j;
  }

  free(comp_size);
  free(cluster_map);

  *out_regions = dup_regions;
  *out_n_regions = n_dup_regions;
}

static int compare_dup_region_by_pos(const void *a, const void *b) {
  const SegtraceDupRegion *ra = (const SegtraceDupRegion *)a,
                          *rb = (const SegtraceDupRegion *)b;
  int c_chr = strcmp(ra->chrom, rb->chrom);
  if (c_chr != 0)
    return c_chr;
  if (ra->start != rb->start)
    return CMP(ra->start, rb->start);
  return CMP(ra->end, rb->end);
}

static int compare_dup_region_by_cluster(const void *a, const void *b) {
  const SegtraceDupRegion *ra = (const SegtraceDupRegion *)a,
                          *rb = (const SegtraceDupRegion *)b;
  uint32_t ca = (uint32_t)strtoul(ra->cluster_id, NULL, 10);
  uint32_t cb = (uint32_t)strtoul(rb->cluster_id, NULL, 10);
  if (ca != cb)
    return CMP(ca, cb);
  int c_chr = strcmp(ra->chrom, rb->chrom);
  if (c_chr != 0)
    return c_chr;
  if (ra->start != rb->start)
    return CMP(ra->start, rb->start);
  return CMP(ra->end, rb->end);
}

size_t merge_dup_regions(SegtraceDupRegion *regions, size_t n) {
  if (n <= 1)
    return n;
  qsort(regions, n, sizeof(SegtraceDupRegion), compare_dup_region_by_cluster);

  size_t out = 0;
  for (size_t i = 1; i < n; i++) {
    if (strcmp(regions[i].cluster_id, regions[out].cluster_id) == 0 &&
        strcmp(regions[i].chrom, regions[out].chrom) == 0 &&
        regions[i].start <= regions[out].end) {
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
      if (regions[i].window_idx > regions[out].window_idx)
        regions[out].window_idx = regions[i].window_idx;
      free(regions[i].cluster_id);
      free(regions[i].chrom);
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
                       size_t n_regions, int n_threads, size_t flank_size) {
  qsort(regions, n_regions, sizeof(SegtraceDupRegion),
        compare_dup_region_by_pos);
  FlankingWorkerData w = {files, r, scale, regions, n_regions, flank_size};
  kt_for(n_threads, extract_flankings_worker, &w, num_files);
}

static void find_chrom_range(const SegtraceDupRegion *regions, size_t n,
                             const char *chr_name, size_t *first,
                             size_t *last) {
  size_t low = 0, high = n;
  while (low < high) {
    size_t mid = low + (high - low) / 2;
    if (strcmp(regions[mid].chrom, chr_name) < 0)
      low = mid + 1;
    else
      high = mid;
  }
  *first = low;
  high = n;
  while (low < high) {
    size_t mid = low + (high - low) / 2;
    if (strcmp(regions[mid].chrom, chr_name) <= 0)
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
    char chr_name[512];
    snprintf(chr_name, sizeof(chr_name), "%s-%s", bname, ks->name.s);

    size_t first, last;
    find_chrom_range(w->regions, w->n_regions, chr_name, &first, &last);

    for (size_t i = first; i < last; i++) {
      size_t start = w->regions[i].start, end = w->regions[i].end,
             flank_size = w->flank_size;
      size_t left_start = start > flank_size ? start - flank_size : 0;
      size_t right_end =
          end + flank_size > ks->seq.l ? ks->seq.l : end + flank_size;
      size_t left_len = start - left_start, right_len = right_end - end;

      free(w->regions[i].flank_sketch.hashes);
      w->regions[i].flank_sketch.hashes = NULL;
      w->regions[i].flank_sketch.sketch_size = 0;

      uint8_t *flank_seq = malloc(left_len + right_len);
      if (left_len > 0)
        memcpy(flank_seq, ks->seq.s + left_start, left_len);
      if (right_len > 0)
        memcpy(flank_seq + left_len, ks->seq.s + end, right_len);

      HashPool pool;
      init_hash_pool(&pool, UINT64_MAX / w->scale);
      extract_hash(w->r, &pool, flank_seq, left_len + right_len);
      finalize_hash_pool(&pool, &w->regions[i].flank_sketch.hashes,
                         &w->regions[i].flank_sketch.sketch_size);
      free(flank_seq);
    }
  }
  kseq_destroy(ks);
  gzclose(fp);
}

typedef struct {
  uint64_t hash;
  uint32_t local_idx;
} FlankHashEntry;

static int compare_flank_hash_entry(const void *a, const void *b) {
  const FlankHashEntry *ea = (const FlankHashEntry *)a,
                       *eb = (const FlankHashEntry *)b;
  return ea->hash != eb->hash ? CMP(ea->hash, eb->hash)
                              : CMP(ea->local_idx, eb->local_idx);
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
    while (j < n_merged &&
           strcmp(regions[i].cluster_id, regions[j].cluster_id) == 0)
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
  size_t min_shared = (size_t)floor((double)min_sz * p_kmer);
  if (min_shared < 1)
    min_shared = 1;

  SegtraceDistResult d = calculate_segtrace_dist(
      &w->regions[ra].flank_sketch, &w->regions[rb].flank_sketch, w->kmer_size);
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

  double p_kmer = pow(0.80, (double)w->kmer_size);

  if (count <= 64) {
    for (size_t a = 0; a < count; a++) {
      for (size_t b = a + 1; b < count; b++) {
        check_and_eval_flank_pair(w, tid, start + a, start + b, p_kmer);
      }
    }
  } else {
    size_t total_flank_hashes = 0;
    for (size_t a = 0; a < count; a++)
      total_flank_hashes += w->regions[start + a].flank_sketch.sketch_size;
    if (total_flank_hashes == 0)
      return;

    FlankHashEntry *entries =
        malloc(total_flank_hashes * sizeof(FlankHashEntry));
    if (!entries)
      return;

    size_t n_entries = 0;
    for (size_t a = 0; a < count; a++) {
      size_t ra = start + a;
      const SegtraceSketch *sk = &w->regions[ra].flank_sketch;
      for (size_t k = 0; k < sk->sketch_size; k++) {
        entries[n_entries++] = (FlankHashEntry){sk->hashes[k], (uint32_t)a};
      }
    }

    qsort(entries, n_entries, sizeof(FlankHashEntry), compare_flank_hash_entry);
    uint8_t *bloom = w->t_bloom[tid];
    memset(bloom, 0, SUBCLUSTER_BLOOM_SIZE_BYTES);

    size_t i = 0;
    while (i < n_entries) {
      size_t j = i + 1;
      while (j < n_entries && entries[j].hash == entries[i].hash)
        j++;
      size_t run_len = j - i;
      if (run_len >= 2 && run_len <= MAX_RUN_LEN) {
        for (size_t a = i; a < j; a++) {
          for (size_t b = a + 1; b < j; b++) {
            uint32_t la = entries[a].local_idx, lb = entries[b].local_idx;
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
}

// ==============================================================
// SECTION 5: REPORTING & FILE OUTPUT WRITERS
// ==============================================================

void write_dup_bed(const char *out_prefix, SegtraceDupRegion *dup_regions,
                   size_t n_merged) {
  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bed", out_prefix);
  FILE *out_bed = fopen(path_buf, "w");
  if (!out_bed) {
    fprintf(stderr, "[ERROR] Cannot open output file: %s\n", path_buf);
    return;
  }

  fprintf(out_bed, "#chrom\tstart\tend\tcluster_id\tsubcluster_id\n");
  for (size_t i = 0; i < n_merged; i++) {
    fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\n", dup_regions[i].chrom,
            dup_regions[i].start, dup_regions[i].end, dup_regions[i].cluster_id,
            dup_regions[i].subcluster_id);
  }
  fclose(out_bed);
}

void write_dup_bedpe(const char *out_prefix, SegtraceDupRegion *dup_regions,
                     size_t n_merged) {
  if (n_merged == 0)
    return;
  SegtraceDupRegion *regions_copy =
      malloc(n_merged * sizeof(SegtraceDupRegion));
  if (!regions_copy)
    return;
  memcpy(regions_copy, dup_regions, n_merged * sizeof(SegtraceDupRegion));
  qsort(regions_copy, n_merged, sizeof(SegtraceDupRegion),
        compare_dup_region_by_cluster);

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bedpe", out_prefix);
  FILE *out_bedpe = fopen(path_buf, "w");
  if (!out_bedpe) {
    free(regions_copy);
    return;
  }

  size_t i = 0;
  while (i < n_merged) {
    size_t j = i + 1;
    while (j < n_merged &&
           strcmp(regions_copy[i].cluster_id, regions_copy[j].cluster_id) == 0)
      j++;

    size_t cluster_size = j - i;
    if (cluster_size > 1 && cluster_size <= 100) {
      for (size_t a = i; a < j; a++) {
        for (size_t b = a + 1; b < j; b++) {
          const char *c1 = regions_copy[a].chrom, *c2 = regions_copy[b].chrom;
          size_t s1 = regions_copy[a].start, e1 = regions_copy[a].end;
          size_t s2 = regions_copy[b].start, e2 = regions_copy[b].end;
          int cmp = strcmp(c1, c2);
          int swap =
              (cmp > 0) || (cmp == 0 && (s1 > s2 || (s1 == s2 && e1 > e2)));
          if (swap)
            fprintf(out_bedpe, "%s\t%zu\t%zu\t%s\t%zu\t%zu\n", c2, s2, e2, c1,
                    s1, e1);
          else
            fprintf(out_bedpe, "%s\t%zu\t%zu\t%s\t%zu\t%zu\n", c1, s1, e1, c2,
                    s2, e2);
        }
      }
    }
    i = j;
  }
  fclose(out_bedpe);
  free(regions_copy);
}

// ==============================================================
// SECTION 6: CORE ALGORITHMS & UTILITIES
// ==============================================================

void init_segtrace(Segtrace *r, size_t hash_window) {
  size_t k = hash_window < 32 ? hash_window : 32;
  uint32_t kmer_bits = 2 * (uint32_t)k;
  r->hash_window = k;
  r->remover_mask =
      (kmer_bits > 2) ? (((uint64_t)1 << (kmer_bits - 2)) - 1) : 0;
  r->kmer_bits = kmer_bits;
  r->rc_shift = (kmer_bits > 0) ? (kmer_bits - 2) : 0;
}

__attribute__((hot)) void extract_hash(const Segtrace *r, HashPool *pool,
                                       const uint8_t *seq, size_t len) {
  if (len < r->hash_window)
    return;
  uint64_t kmer = 0, kmer_rc = 0;
  size_t l = 0;
  for (size_t i = 0; i < len; i++) {
    int8_t c = BASE_LOOKUP[seq[i]];
    if (c < 0) {
      l = 0;
      kmer = 0;
      kmer_rc = 0;
      continue;
    }
    kmer = ((kmer & r->remover_mask) << 2) | (uint64_t)c;
    kmer_rc = (kmer_rc >> 2) | (((uint64_t)(c ^ 3)) << r->rc_shift);
    l++;
    if (l >= r->hash_window) {
      uint64_t canonical = (kmer < kmer_rc) ? kmer : kmer_rc;
      insert_hash_pool(pool, mix_hash(canonical, r->hash_seed));
    }
  }
}

void init_hash_pool(HashPool *pool, uint64_t threshold) {
  pool->size = 0;
  pool->cap = 16;
  pool->hash_threshold = threshold;
  pool->hashes = malloc(pool->cap * sizeof(uint64_t));
}

void insert_hash_pool(HashPool *pool, uint64_t h) {
  if (h >= pool->hash_threshold)
    return;
  DA_PUSH(pool->hashes, pool->size, pool->cap, h);
}

void finalize_hash_pool(HashPool *pool, uint64_t **out_hashes,
                        size_t *out_size) {
  if (pool->size == 0) {
    if (pool->hashes)
      free(pool->hashes);
    *out_hashes = NULL;
    *out_size = 0;
    return;
  }
  qsort(pool->hashes, pool->size, sizeof(uint64_t), compare_uint64);
  size_t unique_count = 0;
  for (size_t i = 0; i < pool->size; i++) {
    if (i == 0 || pool->hashes[i] != pool->hashes[i - 1])
      pool->hashes[unique_count++] = pool->hashes[i];
  }
  *out_hashes = realloc(pool->hashes, unique_count * sizeof(uint64_t));
  *out_size = unique_count;
}

SegtraceDistResult calculate_segtrace_dist(const SegtraceSketch *ref,
                                           const SegtraceSketch *query,
                                           uint32_t kmer_size) {
  SegtraceDistResult res = {0.0, 1.0, 0};
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
  size_t min_sz = ref->sketch_size < query->sketch_size ? ref->sketch_size
                                                        : query->sketch_size;
  res.containment = (double)shared / (double)min_sz;
  res.distance = 1.0 - pow(res.containment, 1.0 / (double)kmer_size);
  return res;
}

void init_unionfind(UnionFind *uf, size_t n) {
  uf->n = n;
  uf->parent = malloc(n * sizeof(uint32_t));
  uf->rank = calloc(n, sizeof(uint8_t));
  for (size_t i = 0; i < n; i++)
    uf->parent[i] = (uint32_t)i;
}

uint32_t find_unionfind(UnionFind *uf, uint32_t x) {
  if (uf->parent[x] != x)
    uf->parent[x] = find_unionfind(uf, uf->parent[x]);
  return uf->parent[x];
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
  if (uf->parent)
    free(uf->parent);
  if (uf->rank)
    free(uf->rank);
  uf->parent = NULL;
  uf->rank = NULL;
  uf->n = 0;
}

void get_basename(const char *filename, char *basename, size_t size) {
  const char *last_slash = strrchr(filename, '/');
  const char *name = last_slash ? last_slash + 1 : filename;
  const char *dot = strrchr(name, '.');
  if (dot && dot != name) {
    size_t len = dot - name;
    if (len >= size)
      len = size - 1;
    strncpy(basename, name, len);
    basename[len] = '\0';
  } else {
    strncpy(basename, name, size - 1);
    basename[size - 1] = '\0';
  }
}

inline uint64_t mix_hash(uint64_t hash_value, uint64_t seed) {
  hash_value ^= seed;
  hash_value ^= hash_value >> 33;
  hash_value *= MIX_CONST1;
  hash_value ^= hash_value >> 33;
  hash_value *= MIX_CONST2;
  hash_value ^= hash_value >> 33;
  return hash_value;
}

size_t lower_bound_u64(const uint64_t *arr, size_t n, uint64_t target) {
  size_t lo = 0, hi = n;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (arr[mid] < target)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
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

int compare_uint64(const void *a, const void *b) {
  uint64_t va = *(const uint64_t *)a, vb = *(const uint64_t *)b;
  return (va > vb) - (va < vb);
}

int compare_hash_entry(const void *a, const void *b) {
  const HashWindowEntry *ea = (const HashWindowEntry *)a,
                        *eb = (const HashWindowEntry *)b;
  return (ea->hash > eb->hash) - (ea->hash < eb->hash);
}
