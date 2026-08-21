#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/kseq.h"
#include "segtrace.h"

/* Reader initialization */
KSEQ_INIT(gzFile, gzread)

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
         "  -c: minimum copies per genome/file to report (default: 2 "
         "[1=duplication map, >=3=polyploid])\n"
         "  -m: filter soft-masked bases (treat lowercase a/c/g/t as valid)\n"
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
      def_kmer_size = (uint32_t)atoi(opt.arg);
    else if (c == 's')
      def_scale = (uint64_t)strtoull(opt.arg, NULL, 10);
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
    else if (c == 'p')
      n_threads = atoi(opt.arg) < 1 ? 1 : atoi(opt.arg);
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

  Segtrace r;
  init_segtrace(&r, def_kmer_size, filter_masked);
  r.hash_seed = def_hash_seed;

  GlobalWindows gw =
      extract_all_windows(files, num_files, &r, def_scale, window_size,
                          step_size, min_bases, n_threads);

  UnionFind uf;
  init_unionfind(&uf, gw.num_sketches);

  fprintf(stderr,
          "[segtrace] Discovering candidates and computing distances...\n");
  discover_and_compute(gw.all_hashes, gw.coords, gw.num_sketches, window_size,
                       step_size, n_threads, r.hash_window, &uf);

  SegtraceDupRegion *dup_regions = NULL;
  size_t n_dup_regions = 0;
  build_duplicate_regions(&uf, gw.num_sketches, gw.coords, gw.seq_lens,
                          step_size, window_size, &dup_regions, &n_dup_regions);

  free_unionfind(&uf);
  free(gw.all_hashes);
  gw.all_hashes = NULL;
  free(gw.coords);
  gw.coords = NULL;

  size_t n_merged = merge_dup_regions(dup_regions, n_dup_regions, window_size);
  size_t n_filtered =
      filter_regions_by_copy_count(dup_regions, n_merged, min_copies);

  size_t min_sd_len = window_size < MIN_SD_LEN ? window_size : MIN_SD_LEN;
  write_dup_bed(out_prefix, dup_regions, n_filtered, gw.seq_lens, min_sd_len);

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
      if (b_out < 0 || b_out > 3) {
        valid_len = 0;
        f_hash = 0;
        r_hash = 0;
        continue;
      }
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
  if (idx + job->window_size > job->chunk_end_idx ||
      idx + job->window_size > job->seq_len)
    return;

  size_t valid_bases = 0;
  for (size_t j = 0; j < job->window_size; j++) {
    if (job->base_lookup[job->seq_ptr[idx + j]] >= 0)
      valid_bases++;
  }

  for (; idx + job->window_size <= job->chunk_end_idx &&
         idx + job->window_size <= job->seq_len;
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

    size_t next_idx = idx + job->step_size;
    if (next_idx + job->window_size <= job->chunk_end_idx &&
        next_idx + job->window_size <= job->seq_len) {
      for (size_t k = 0; k < job->step_size; k++) {
        if (job->base_lookup[job->seq_ptr[idx + k]] >= 0)
          valid_bases--;
        if (job->base_lookup[job->seq_ptr[idx + job->window_size + k]] >= 0)
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

      for (size_t c_start = 0; c_start < len; c_start += chunk_size) {
        size_t c_end = c_start + chunk_size + window_size - step_size;
        if (c_start + window_size > len)
          break;

        DA_RESERVE(jobs, cap_jobs, num_jobs + 1);
        jobs[num_jobs++] = (SeqChunkJob){.r = r,
                                         .base_lookup = r->base_lookup,
                                         .threshold = threshold,
                                         .window_size = window_size,
                                         .step_size = step_size,
                                         .min_bases = min_bases,
                                         .seq_id = seq_id,
                                         .seq_ptr = seq_ptr,
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

static inline size_t calculate_sketch_dist_fast(const uint32_t *a, size_t n_a,
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
  if (wa->sketch_size == 0 || wb->sketch_size == 0)
    return 0;
  return calculate_sketch_dist_fast(
      all_hashes + wa->sketch_offset, wa->sketch_size,
      all_hashes + wb->sketch_offset, wb->sketch_size);
}

static inline int check_collinear_neighbor(DiscoverComputeData *w, uint32_t wa,
                                           uint32_t wb, size_t min_shared) {
  uint32_t seq_a = w->coords[wa].seq_id;
  uint32_t seq_b = w->coords[wb].seq_id;
  long long n_win = (long long)w->n_windows;

  const int dir_a[] = {1, -1, 1, -1};
  const int dir_b[] = {1, -1, -1, 1};

  // Pass 1: Exact diagonals
  for (int d = 0; d < 4; d++) {
    int da = dir_a[d], db = dir_b[d];
    for (int step = 1; step <= MAX_COLLINEAR_LOOKAHEAD; step++) {
      long long next_a = (long long)wa + da * step;
      long long next_b = (long long)wb + db * step;
      if (next_a >= 0 && next_a < n_win && next_b >= 0 && next_b < n_win &&
          w->coords[next_a].seq_id == seq_a &&
          w->coords[next_b].seq_id == seq_b) {
        if (w->coords[next_a].sketch_size > 0 &&
            w->coords[next_b].sketch_size > 0) {
          size_t shared = calculate_window_dist(
              w->all_hashes, &w->coords[next_a], &w->coords[next_b]);
          if (shared >= min_shared)
            return 1;
        }
      }
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
        if (next_a >= 0 && next_a < n_win && next_b >= 0 && next_b < n_win &&
            w->coords[next_a].seq_id == seq_a &&
            w->coords[next_b].seq_id == seq_b) {
          if (w->coords[next_a].sketch_size > 0 &&
              w->coords[next_b].sketch_size > 0) {
            size_t shared = calculate_window_dist(
                w->all_hashes, &w->coords[next_a], &w->coords[next_b]);
            if (shared >= min_shared)
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

  double p_kmer = w_data->p_kmer;
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
          size_t min_shared = (size_t)ceil((double)min_sz * p_kmer);
          if (min_shared < 2)
            min_shared = 2;

          size_t shared = calculate_window_dist(
              w_data->all_hashes, &w_data->coords[wa], &w_data->coords[wb]);
          if (shared >= min_shared) {
            if (check_collinear_neighbor(w_data, wa, wb, min_shared)) {
              DA_PUSH(w_data->t_pairs[tid], w_data->t_n_pairs[tid],
                      w_data->t_cap_pairs[tid], ((CandidatePair){wa, wb}));
            }
          }
        }
      }
    }
    i = j;
  }
}

void discover_and_compute(const uint32_t *all_hashes, const WindowCoord *coords,
                          size_t n_windows, size_t window_size,
                          size_t step_size, int n_threads, uint32_t kmer_size,
                          UnionFind *uf) {
  DiscoverComputeData w = {
      .all_hashes = all_hashes,
      .coords = coords,
      .n_windows = n_windows,
      .window_size = window_size,
      .step_size = step_size,
      .kmer_size = kmer_size,
      .p_kmer = pow(0.80, (double)kmer_size),
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

    for (size_t p = batch_start; p < batch_end; p++) {
      w.buckets[p].size = 0;
      w.buckets[p].cap = 0;
      w.buckets[p].entries = NULL;
    }

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

    for (int t = 0; t < n_threads; t++) {
      for (size_t k = 0; k < w.t_n_pairs[t]; k++) {
        union_unionfind(uf, w.t_pairs[t][k].a, w.t_pairs[t][k].b);
      }
      w.t_n_pairs[t] = 0;
    }

    for (size_t p = batch_start; p < batch_end; p++) {
      free(w.buckets[p].entries);
      w.buckets[p].entries = NULL;
      w.buckets[p].size = 0;
      w.buckets[p].cap = 0;
    }
  }

  free(win_curr_pos);
  for (int t = 0; t < n_threads; t++) {
    free(w.t_bloom[t]);
    free(w.t_pairs[t]);
  }
  free(w.buckets);
  free(w.t_bloom);
  free(w.t_pairs);
  free(w.t_n_pairs);
  free(w.t_cap_pairs);
}

// ==============================================================
// SECTION 4: CLUSTERING, LOCUS MERGING & COPY FILTERING
// ==============================================================

void build_duplicate_regions(UnionFind *uf, size_t num_sketches,
                             const WindowCoord *coords,
                             const GenomeSeqLen *seq_lens, size_t step_size,
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
    if (cid == 0)
      continue;

    uint32_t seq_i = coords[i].seq_id;
    size_t start = (size_t)coords[i].window_idx * step_size;
    size_t end = start + window_size;

    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((SegtraceDupRegion){.seq_id = seq_i,
                                 .file_id = seq_lens[seq_i].file_id,
                                 .start = start,
                                 .end = end,
                                 .cluster_id = cid}));
  }

  free(comp_size);
  free(cluster_map);

  *out_regions = dup_regions;
  *out_n_regions = n_dup_regions;
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

size_t merge_dup_regions(SegtraceDupRegion *regions, size_t n,
                         size_t window_size) {
  if (n <= 1)
    return n;
  qsort(regions, n, sizeof(SegtraceDupRegion),
        compare_dup_region_by_cluster_file);

  size_t out = 0;
  for (size_t i = 1; i < n; i++) {
    if (regions[i].cluster_id == regions[out].cluster_id &&
        regions[i].file_id == regions[out].file_id &&
        regions[i].seq_id == regions[out].seq_id &&
        regions[i].start <= regions[out].end + MERGE_COEFF * window_size) {
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
    } else {
      out++;
      if (out != i)
        regions[out] = regions[i];
    }
  }
  return out + 1;
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
    if (copy_count >= min_copies) {
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

void init_segtrace(Segtrace *r, size_t hash_window, int filter_masked) {
  r->hash_window = (uint32_t)hash_window;
  r->filter_masked = filter_masked;
  r->base_lookup = filter_masked ? BASE_LOOKUP : BASE_LOOKUP_NO_MASK;
}

void init_unionfind(UnionFind *uf, size_t n) {
  uf->n = n;
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
