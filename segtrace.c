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
         "  -k: kmer size (default: 15)\n"
         "  -s: scale factor (default: 8)\n"
         "  -e: hash seed (default: 42)\n"
         "  -w: window size in bp (default: 1024)\n"
         "  -t: step size in bp (default: 0 [auto: 33%% of window size])\n"
         "  -b: minimum valid bases per window (default: 0 [auto: 25%% of "
         "window size])\n"
         "  -m: not filtering soft-masked bases (treat lowercase a/c/g/t as "
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

  uint32_t def_kmer_size = 15;
  uint64_t def_scale = 8, def_hash_seed = 42;
  size_t window_size = 1024, step_size = 0, min_bases = 0, flank_size = 2048;
  const char *out_prefix = "segtrace";
  int n_threads = 8, filter_masked = 1;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:e:w:t:b:d:o:p:D:f:mh", 0)) >=
         0) {
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
  build_duplicate_regions(&uf, num_sketches, seq_lens, coords, &dup_regions,
                          &n_dup_regions);

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

  size_t chunk_len = job->chunk_end_idx - job->chunk_start_idx;
  PosHashPool chunk_pool = {NULL, 0, 0};
  extract_hash_chunk(job->r, &chunk_pool, job->seq_ptr + job->chunk_start_idx,
                     chunk_len, UINT64_MAX / job->scale);

  uint32_t current_window_idx =
      (uint32_t)(job->chunk_start_idx / job->step_size);
  size_t p_start = 0, p_end = 0;

  for (size_t idx = job->chunk_start_idx;
       idx + job->window_size <= job->chunk_end_idx;
       idx += job->step_size, current_window_idx++) {

    size_t w_start = idx - job->chunk_start_idx;
    size_t w_end = w_start + job->window_size;
    size_t min_i = w_start + job->r->hash_window - 1;
    size_t max_i = w_end;

    while (p_start < chunk_pool.size && chunk_pool.entries[p_start].pos < min_i)
      p_start++;
    while (p_end < chunk_pool.size && chunk_pool.entries[p_end].pos < max_i)
      p_end++;

    size_t valid_bases = 0;
    for (size_t j = 0; j < job->window_size; j++) {
      if (job->base_lookup[job->seq_ptr[idx + j]] >= 0)
        valid_bases++;
    }

    size_t sketch_size = 0;
    uint64_t *hashes = NULL;
    if (valid_bases >= job->min_bases) {
      size_t count = p_end - p_start;
      if (count > 0) {
        uint64_t *tmp_hashes = malloc(count * sizeof(uint64_t));
        for (size_t k = 0; k < count; k++)
          tmp_hashes[k] = chunk_pool.entries[p_start + k].hash;
        qsort(tmp_hashes, count, sizeof(uint64_t), compare_uint64);
        size_t unique_count = 0;
        for (size_t k = 0; k < count; k++) {
          if (k == 0 || tmp_hashes[k] != tmp_hashes[k - 1])
            tmp_hashes[unique_count++] = tmp_hashes[k];
        }
        hashes = realloc(tmp_hashes, unique_count * sizeof(uint64_t));
        sketch_size = unique_count;
      }
    }

    DA_RESERVE(job->coords, job->cap_coords, job->num_coords + 1);
    WindowCoord *wc = &job->coords[job->num_coords++];
    wc->seq_id = job->seq_id;
    wc->start = idx;
    wc->end = idx + job->window_size;
    wc->window_idx = current_window_idx;

    size_t h_idx = job->num_hashes;
    if (sketch_size > 0) {
      DA_RESERVE(job->hashes, job->cap_hashes, job->num_hashes + sketch_size);
      memcpy(job->hashes + h_idx, hashes, sketch_size * sizeof(uint64_t));
      job->num_hashes += sketch_size;
    }

    wc->sketch_offset = h_idx;
    wc->sketch_size = (uint32_t)sketch_size;
    if (hashes)
      free(hashes);
  }
  if (chunk_pool.entries)
    free(chunk_pool.entries);
}

static void process_and_gather_jobs(StreamWorkerData *worker, SeqChunkJob *jobs,
                                    size_t num_jobs, int n_threads) {
  if (num_jobs == 0)
    return;
  kt_for(n_threads, seq_chunk_worker, jobs, num_jobs);

  for (size_t j = 0; j < num_jobs; j++) {
    SeqChunkJob *job = &jobs[j];
    if (job->num_coords > 0) {
      size_t base_h_offset = worker->num_all_hashes;
      DA_RESERVE(worker->all_hashes, worker->cap_all_hashes,
                 worker->num_all_hashes + job->num_hashes);
      memcpy(worker->all_hashes + base_h_offset, job->hashes,
             job->num_hashes * sizeof(uint64_t));
      worker->num_all_hashes += job->num_hashes;

      size_t base_c_offset = worker->num_sketches;
      DA_RESERVE(worker->coords, worker->cap_sketches,
                 worker->num_sketches + job->num_coords);

      for (size_t k = 0; k < job->num_coords; k++) {
        WindowCoord wc = job->coords[k];
        wc.sketch_offset += base_h_offset;
        worker->coords[base_c_offset + k] = wc;
      }
      worker->num_sketches += job->num_coords;
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
  size_t batch_seq_len = 0;

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
      if (!seq_copy) {
        fprintf(stderr, "[ERROR] Memory allocation failed for seq_copy\n");
        exit(1);
      }
      memcpy(seq_copy, ks->seq.s, len + 1);
      batch_seq_len += len;

      size_t chunk_size = len / (n_threads * 4);
      if (chunk_size < 100000)
        chunk_size = 100000;
      chunk_size = ((chunk_size + step_size - 1) / step_size) * step_size;

      for (size_t c_start = 0; c_start < len; c_start += chunk_size) {
        size_t c_end = c_start + chunk_size + window_size - step_size;
        if (c_end > len)
          c_end = len;
        if (c_end - c_start < window_size)
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
                                         .chunk_end_idx = c_end};
      }

      if (num_jobs >= (size_t)(n_threads * 4) || batch_seq_len > 250000000) {
        process_and_gather_jobs(&workers[0], jobs, num_jobs, n_threads);
        num_jobs = 0;
        batch_seq_len = 0;
      }
    }
    kseq_destroy(ks);
    gzclose(fp);
  }

  process_and_gather_jobs(&workers[0], jobs, num_jobs, n_threads);
  free(jobs);
  return workers;
}

void merge_global_data(StreamWorkerData *workers, int num_files,
                       const char *out_prefix, uint64_t **out_all_hashes,
                       WindowCoord **out_coords, size_t *out_num_sketches,
                       GenomeSeqLen **out_seq_lens, size_t *out_num_seqs) {
  (void)out_prefix;
  (void)num_files;

  *out_all_hashes = workers[0].all_hashes;
  *out_coords = workers[0].coords;
  *out_num_sketches = workers[0].num_sketches;
  *out_seq_lens = workers[0].seq_lens;
  *out_num_seqs = workers[0].num_seqs;
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

      PartitionBucket *b = &w.buckets[p];
      size_t rem = b->size % BUCKET_BLOCK_SIZE;
      if (rem == 0) {
        HashBlock *nb = malloc(sizeof(HashBlock));
        nb->next = NULL;
        if (!b->head) {
          b->head = b->tail = nb;
        } else {
          b->tail->next = nb;
          b->tail = nb;
        }
      }
      b->tail->entries[rem] = (HashWindowEntry){val, (uint32_t)win};
      b->size++;
    }
  }

  kt_for(n_threads, discover_compute_worker, &w, NUM_PARTITIONS);

  for (int t = 0; t < n_threads; t++)
    free(w.t_bloom[t]);
  for (size_t p = 0; p < NUM_PARTITIONS; p++) {
    HashBlock *curr = w.buckets[p].head;
    while (curr) {
      HashBlock *next = curr->next;
      free(curr);
      curr = next;
    }
  }
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
  const int dir_a[] = {1, -1, 1, -1};
  const int dir_b[] = {1, -1, -1, 1};

  for (int d = 0; d < 4; d++) {
    for (int step_a = 1; step_a <= MAX_COLLINEAR_LOOOKAHEAD; step_a++) {
      for (int step_b = 1; step_b <= MAX_COLLINEAR_LOOOKAHEAD; step_b++) {
        long long next_a = (long long)wa + dir_a[d] * step_a;
        long long next_b = (long long)wb + dir_b[d] * step_b;

        if (next_a >= 0 && next_a < (long long)w->n_windows && next_b >= 0 &&
            next_b < (long long)w->n_windows &&
            w->coords[next_a].seq_id == w->coords[wa].seq_id &&
            w->coords[next_b].seq_id == w->coords[wb].seq_id) {

          if (w->coords[next_a].sketch_size > 0 &&
              w->coords[next_b].sketch_size > 0) {
            if (calculate_window_dist(w->all_hashes, &w->coords[next_a],
                                      &w->coords[next_b], w->kmer_size)
                    .shared_hashes >= min_shared)
              return 1;
          }
        }
      }
    }
  }
  return 0;
}

void discover_compute_worker(void *data, long p, int tid) {
  DiscoverComputeData *w_data = (DiscoverComputeData *)data;
  PartitionBucket *b = &w_data->buckets[p];
  if (b->size == 0)
    return;

  HashWindowEntry *flat_entries = malloc(b->size * sizeof(HashWindowEntry));
  size_t idx = 0;
  HashBlock *curr = b->head;
  while (curr) {
    size_t count = (b->size - idx > BUCKET_BLOCK_SIZE) ? BUCKET_BLOCK_SIZE
                                                       : (b->size - idx);
    memcpy(flat_entries + idx, curr->entries, count * sizeof(HashWindowEntry));
    idx += count;
    curr = curr->next;
  }

  qsort(flat_entries, b->size, sizeof(HashWindowEntry), compare_hash_entry);
  memset(w_data->t_bloom[tid], 0, BLOOM_SIZE_BYTES);

  double p_kmer = pow(0.90, (double)w_data->kmer_size);
  size_t i = 0;
  while (i < b->size) {
    size_t j = i + 1;
    while (j < b->size && flat_entries[j].hash == flat_entries[i].hash)
      j++;
    size_t run_len = j - i;

    if (run_len >= 2 && run_len <= MAX_KMER_FREQ) {
      for (size_t a = i; a < j; a++) {
        size_t b_max =
            a + 1 + MAX_PAIR_COMPARISONS < j ? a + 1 + MAX_PAIR_COMPARISONS : j;
        for (size_t b_idx = a + 1; b_idx < b_max; b_idx++) {
          uint32_t wa = flat_entries[a].window_id,
                   wb = flat_entries[b_idx].window_id;
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
          size_t min_shared = (size_t)ceil((double)min_sz * p_kmer) + 1;
          if (min_shared < 2)
            min_shared = 2;

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
  free(flat_entries);
}

// ==============================================================
// SECTION 4: CLUSTERING, LOCUS MERGING & FLANKING SUBCLUSTERING
// ==============================================================

void build_duplicate_regions(UnionFind *uf, size_t num_sketches,
                             GenomeSeqLen *seq_lens, WindowCoord *coords,
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

    char label[32];
    snprintf(label, sizeof(label), "%u", cid);

    char chrom_name[512];
    uint32_t seq_i = coords[i].seq_id;
    snprintf(chrom_name, sizeof(chrom_name), "%s-%s", seq_lens[seq_i].genome,
             seq_lens[seq_i].seq);

    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((SegtraceDupRegion){.chrom = strdup(chrom_name),
                                 .start = coords[i].start,
                                 .end = coords[i].end,
                                 .cluster_id = strdup(label),
                                 .copy_count = comp_size[root_i],
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

      if (left_len + right_len == 0)
        continue;

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
  size_t min_shared = (size_t)ceil((double)min_sz * p_kmer) + 1;
  if (min_shared < 2)
    min_shared = 2;

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
                   size_t n_merged) {
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
  size_t i = 0;
  while (i < n_merged) {
    size_t j = i + 1;
    while (j < n_merged &&
           strcmp(dup_regions[i].cluster_id, dup_regions[j].cluster_id) == 0)
      j++;
    size_t cluster_size = j - i;
    if (cluster_size >= 2) {
      int valid_regions = 0;
      for (size_t k = i; k < j; k++) {
        if (dup_regions[k].end - dup_regions[k].start >= MIN_SD_LEN)
          valid_regions++;
      }
      if (valid_regions >= 2) {
        for (size_t k = i; k < j; k++) {
          if (dup_regions[k].end - dup_regions[k].start >= MIN_SD_LEN) {
            fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\n", dup_regions[k].chrom,
                    dup_regions[k].start, dup_regions[k].end,
                    dup_regions[k].cluster_id, dup_regions[k].subcluster_id);
          }
        }
      }
    }
    i = j;
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

__attribute__((hot)) void extract_hash_chunk(const Segtrace *r,
                                             PosHashPool *pool,
                                             const uint8_t *seq, size_t len,
                                             uint64_t threshold) {
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
      int8_t b_in = b;
      f_hash = rol64(f_hash, 1) ^ rol64(NTHASH_H[b_out], k) ^ NTHASH_H[b_in];
      r_hash = ror64(r_hash, 1) ^ ror64(NTHASH_H[b_out ^ 3], 1) ^
               rol64(NTHASH_H[b_in ^ 3], k - 1);
    }

    if (valid_len >= k) {
      uint64_t canonical = (f_hash < r_hash) ? f_hash : r_hash;
      uint64_t h = mix_hash(canonical, r->hash_seed);
      if (h < threshold) {
        DA_PUSH(pool->entries, pool->size, pool->cap,
                ((PosHash){(uint32_t)i, h}));
      }
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
  return ea->hash != eb->hash ? CMP(ea->hash, eb->hash)
                              : CMP(ea->window_id, eb->window_id);
}
