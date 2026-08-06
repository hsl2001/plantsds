#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/khash.h"
#include "klib/kseq.h"
#include "segtrace.h"

/* Reader initiation */
KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(genome_map, uint32_t)

/* Segtrace encoding: A = 00, C = 01, G = 10, T = 11 */
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
// 1. ENTRY POINT & CLI
// ==============================================================

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

  // ==============================================================
  // CLI defaults
  // ==============================================================

  uint32_t def_kmer_size = 21;
  uint64_t def_scale = 16;
  uint64_t def_hash_seed = 42;
  size_t window_size = 1024;
  size_t step_size = 0; /* 0 = auto (window/2) */
  size_t min_bases = 1000;
  double max_dist = 0.15;
  const char *out_prefix = "segtrace";
  int n_threads = 8;
  uint32_t adjacency_threshold = 2;
  double subcluster_dist = 0.2; /* -1.0: auto */
  size_t flank_size = 2000;

  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:e:w:t:b:d:o:p:a:D:f:h", 0)) >=
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
    else if (c == 'd')
      max_dist = atof(opt.arg);
    else if (c == 'o')
      out_prefix = opt.arg;
    else if (c == 'p')
      n_threads = atoi(opt.arg) < 1 ? 1 : atoi(opt.arg);
    else if (c == 'a')
      adjacency_threshold = (uint32_t)atoi(opt.arg);
    else if (c == 'D')
      subcluster_dist = atof(opt.arg);
    else if (c == 'f')
      flank_size = (size_t)strtoull(opt.arg, NULL, 10);
    else
      return 1;
  }

  if (subcluster_dist < 0.0)
    subcluster_dist = max_dist;

  if (step_size == 0)
    step_size = window_size * 4 / 5;

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
  size_t num_sketches = 0;
  GenomeSeqLen *seq_lens = NULL;
  size_t num_seqs = 0;

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
  discover_and_compute(all_hashes, coords, num_sketches, window_size, max_dist,
                       n_threads, r.hash_window, &uf);

  SegtraceDupRegion *dup_regions = NULL;
  size_t n_dup_regions = 0;
  build_duplicate_regions(&uf, num_sketches, num_files, files, seq_lens, coords,
                          &dup_regions, &n_dup_regions);

  size_t n_merged =
      merge_dup_regions(dup_regions, n_dup_regions, adjacency_threshold);

  fprintf(stderr,
          "[INFO] Extracting flanking sequences for sub-clustering...\n");
  extract_flankings(files, num_files, &r, def_scale, dup_regions, n_merged,
                    n_threads, flank_size);

  fprintf(stderr, "[INFO] Sub-clustering based on flanking similarities...\n");
  perform_subclustering(dup_regions, n_merged, subcluster_dist, n_threads,
                        r.hash_window);

  write_dup_bed(out_prefix, dup_regions, n_merged);
  write_dup_bedpe(out_prefix, dup_regions, n_merged);

  for (size_t i = 0; i < n_merged; i++) {
    free(dup_regions[i].chrom);
    free(dup_regions[i].cluster_id);
    free(dup_regions[i].flank_sketch.hashes);
  }
  free(dup_regions);
  free_unionfind(&uf);
  free(all_hashes);
  free(coords);
  for (size_t i = 0; i < num_seqs; i++) {
    free(seq_lens[i].genome);
    free(seq_lens[i].seq);
  }
  free(seq_lens);
  return 0;
}

void print_usage(void) {
  printf("Segtrace: Segmental Duplication Tracer\n\n"
         "Usage: segtrace [options] fasta1 [fasta2 ...]\n\n"
         "Options:\n"
         "  -k: kmer size (default: 21)\n"
         "  -s: scale factor (default: 10)\n"
         "  -e: hash seed (default: 42)\n"
         "  -w: window size in bp (default: 1000)\n"
         "  -t: step size in bp (default: 0 [auto: 80%% of window size])\n"
         "  -b: minimum valid bases per window (default: 1000)\n"
         "  -d: maximum distance to consider as copy (default: 0.15)\n"
         "  -D: sub-cluster distance threshold (default: 0.2)\n"
         "  -f: flanking size in bp for sub-clustering (default: 2000)\n"
         "  -a: adjacency threshold for merging regions (default: 2)\n"
         "  -o: output file prefix (default: segtrace)\n"
         "  -p: number of threads (default: 8)\n"
         "  -h, --help: show this help message\n"
         "\n");
}

// ==============================================================
// 3. WINDOW EXTRACTION
// ==============================================================

StreamWorkerData *extract_all_windows(char **files, int num_files,
                                      const Segtrace *r, uint64_t scale,
                                      size_t window_size, size_t step_size,
                                      size_t min_bases, int n_threads) {
  fprintf(stderr, "[segtrace] Extracting windows across pangenome...\n");
  StreamWorkerData *workers = calloc(num_files, sizeof(StreamWorkerData));
  if (!workers) {
    fprintf(stderr, "[ERROR] Failed to allocate memory for workers\n");
    exit(1);
  }
  for (int i = 0; i < num_files; i++) {
    workers[i].filename = files[i];
    get_basename(files[i], workers[i].bname, sizeof(workers[i].bname));
    workers[i].r = r;
    workers[i].scale = scale;
    workers[i].window_size = window_size;
    workers[i].step_size = step_size;
    workers[i].min_bases = min_bases;
  }
  kt_for(n_threads, stream_pangenome_worker, workers, num_files);
  return workers;
}

void stream_pangenome_worker(void *data, long i, int tid) {
  (void)tid;
  StreamWorkerData *w = &((StreamWorkerData *)data)[i];

  gzFile fp = gzopen(w->filename, "r");
  if (!fp)
    return;
  kseq_t *ks = kseq_init(fp);
  if (!ks) {
    gzclose(fp);
    return;
  }

  while (kseq_read(ks) >= 0) {
    size_t len = ks->seq.l;

    DA_RESERVE(w->seq_lens, w->cap_seqs, w->num_seqs + 1);
    w->seq_lens[w->num_seqs].genome = strdup(w->bname);
    w->seq_lens[w->num_seqs].seq = strdup(ks->name.s);
    w->num_seqs++;

    /* Pre-allocate estimated windows to reduce realloc overhead */
    size_t est_windows =
        len >= w->window_size ? (len - w->window_size) / w->step_size + 1 : 0;
    DA_RESERVE(w->coords, w->cap_sketches, w->num_sketches + est_windows);

    uint32_t current_window_idx = 0;
    for (size_t idx = 0; idx + w->window_size <= len;
         idx += w->step_size, current_window_idx++) {
      size_t valid_bases = 0;
      for (size_t j = 0; j < w->window_size; j++) {
        if (BASE_LOOKUP[(uint8_t)ks->seq.s[idx + j]] >= 0)
          valid_bases++;
      }
      if (valid_bases < w->min_bases)
        continue;

      DA_RESERVE(w->coords, w->cap_sketches, w->num_sketches + 1);

      WindowCoord *wc = &w->coords[w->num_sketches];

      wc->seq_id = (uint32_t)(w->num_seqs - 1);
      wc->start = idx;
      wc->end = idx + w->window_size;
      wc->window_idx = current_window_idx;

      /* Extract sketch, write to disk, free immediately */
      HashPool pool;
      init_hash_pool(&pool, UINT64_MAX / w->scale);
      extract_hash(w->r, &pool, (const uint8_t *)ks->seq.s + idx,
                   w->window_size);

      uint64_t *hashes = NULL;
      size_t sketch_size = 0;
      finalize_hash_pool(&pool, &hashes, &sketch_size);

      if (sketch_size > 0) {
        size_t hash_idx = w->num_all_hashes;
        if (w->cap_all_hashes == 0) {
          w->cap_all_hashes = 524288;
          w->all_hashes = malloc(w->cap_all_hashes * sizeof(uint64_t));
        }
        DA_RESERVE(w->all_hashes, w->cap_all_hashes,
                   w->num_all_hashes + sketch_size);
        memcpy(w->all_hashes + hash_idx, hashes,
               sketch_size * sizeof(uint64_t));
        w->num_all_hashes += sketch_size;

        wc->sketch_offset = hash_idx;
        wc->sketch_size = (uint32_t)sketch_size;
        free(hashes);
        w->num_sketches++;
      }
      /* finalize_hash_pool handles cleanup when sketch_size == 0 */
    }
  }
  kseq_destroy(ks);
  gzclose(fp);
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

  size_t g_hash_offset = 0;
  size_t g_sketch_offset = 0;
  size_t g_seq_offset = 0;

  char path_buf[PATH_MAX];
  FILE *bed_fp;
  snprintf(path_buf, sizeof(path_buf), "%s.window.bed", out_prefix);
  bed_fp = fopen(path_buf, "w");

  for (int i = 0; i < num_files; i++) {
    StreamWorkerData *w = &workers[i];

    if (w->num_all_hashes > 0) {
      memcpy(all_hashes + g_hash_offset, w->all_hashes,
             w->num_all_hashes * sizeof(uint64_t));
    }

    if (w->num_seqs > 0) {
      memcpy(seq_lens + g_seq_offset, w->seq_lens,
             w->num_seqs * sizeof(GenomeSeqLen));
    }

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
// 4. DISCOVERY & DISTANCE CALCULATION
// ==============================================================

/* Fused discovery + distance computation + UF union.
 * Eliminates the intermediate cand_pairs array entirely. */
void discover_and_compute(const uint64_t *all_hashes, WindowCoord *coords,
                          size_t n_windows, size_t window_size, double max_dist,
                          int n_threads, uint32_t kmer_size, UnionFind *uf) {
  DiscoverComputeData w;
  w.all_hashes = all_hashes;
  w.coords = coords;
  w.n_windows = n_windows;
  w.window_size = window_size;
  w.max_dist = max_dist;
  w.kmer_size = kmer_size;
  w.t_edges = calloc(n_threads, sizeof(SegtraceDupEdge *));
  w.t_n_edges = calloc(n_threads, sizeof(size_t));
  w.t_cap_edges = calloc(n_threads, sizeof(size_t));
  w.t_bloom = malloc(n_threads * sizeof(uint8_t *));
  if (!w.t_edges || !w.t_n_edges || !w.t_cap_edges || !w.t_bloom) {
    fprintf(stderr, "[ERROR] Failed to allocate memory for discovery data\n");
    exit(1);
  }
  for (int t = 0; t < n_threads; t++)
    w.t_bloom[t] = calloc(BLOOM_SIZE_BYTES, 1);

  kt_for(n_threads, discover_compute_worker, &w, NUM_PARTITIONS);

  /* Union all discovered edges into UF */
  size_t total_edges = 0;
  for (int t = 0; t < n_threads; t++) {
    for (size_t k = 0; k < w.t_n_edges[t]; k++) {
      union_unionfind(uf, w.t_edges[t][k].win_a, w.t_edges[t][k].win_b);
    }
    total_edges += w.t_n_edges[t];
    free(w.t_edges[t]);
    free(w.t_bloom[t]);
  }
  free(w.t_edges);
  free(w.t_n_edges);
  free(w.t_cap_edges);
  free(w.t_bloom);

  fprintf(stderr, "[segtrace] Total edges after distance filter: %zu\n",
          total_edges);
}

void discover_compute_worker(void *data, long p, int tid) {
  DiscoverComputeData *w_data = (DiscoverComputeData *)data;

  /* Partition boundaries: divide [0, UINT64_MAX] into NUM_PARTITIONS */
  uint64_t part_size = UINT64_MAX / NUM_PARTITIONS;
  uint64_t lo = part_size * (uint64_t)p;
  uint64_t hi = (p == NUM_PARTITIONS - 1) ? UINT64_MAX
                                          : part_size * (uint64_t)(p + 1) - 1;

  /* Pre-count entries for exact allocation (avoids realloc overhead) */
  size_t total_in_partition = 0;
  for (size_t w = 0; w < w_data->n_windows; w++) {
    const uint64_t *h = w_data->all_hashes + w_data->coords[w].sketch_offset;
    size_t sz = w_data->coords[w].sketch_size;
    size_t s = lower_bound_u64(h, sz, lo);
    size_t e = (hi == UINT64_MAX) ? sz : lower_bound_u64(h, sz, hi + 1);
    total_in_partition += e - s;
  }

  if (total_in_partition == 0)
    return;

  /* Allocate exact size for entries */
  HashWindowEntry *entries =
      malloc(total_in_partition * sizeof(HashWindowEntry));
  size_t n_entries = 0;

  for (size_t w = 0; w < w_data->n_windows; w++) {
    const uint64_t *h = w_data->all_hashes + w_data->coords[w].sketch_offset;
    size_t sz = w_data->coords[w].sketch_size;
    size_t s = lower_bound_u64(h, sz, lo);
    size_t e = (hi == UINT64_MAX) ? sz : lower_bound_u64(h, sz, hi + 1);

    for (size_t k = s; k < e; k++)
      entries[n_entries++] = (HashWindowEntry){h[k], (uint32_t)w};
  }

  /* Sort by hash value within this partition */
  qsort(entries, n_entries, sizeof(HashWindowEntry), compare_hash_entry);

  /* Reset bloom filter for this partition */
  memset(w_data->t_bloom[tid], 0, BLOOM_SIZE_BYTES);

  /* Scan runs of identical hashes -> compute distances immediately */
  size_t i = 0;
  while (i < n_entries) {
    size_t j = i + 1;
    while (j < n_entries && entries[j].hash == entries[i].hash)
      j++;
    size_t run_len = j - i;

    if (run_len >= 2 && run_len <= MAX_RUN_LEN) {
      for (size_t a = i; a < j; a++) {
        for (size_t b = a + 1; b < j; b++) {
          uint32_t wa = entries[a].window_id;
          uint32_t wb = entries[b].window_id;
          /* Skip overlapping windows on same chromosome */
          if (w_data->coords[wa].seq_id == w_data->coords[wb].seq_id &&
              ABS_DIFF(w_data->coords[wa].start, w_data->coords[wb].start) <
                  w_data->window_size)
            continue;

          /* Bloom-filter approximate dedup (cheap, replaces khash) */
          uint64_t pk = encode_pair(wa, wb);
          if (bloom_test_and_set(w_data->t_bloom[tid], pk))
            continue;

          /* Compute distance immediately (fused Phase 2) */
          SegtraceDistResult d =
              calculate_window_dist(w_data->all_hashes, &w_data->coords[wa],
                                    &w_data->coords[wb], w_data->kmer_size);
          if (d.distance < w_data->max_dist) {
            DA_PUSH(w_data->t_edges[tid], w_data->t_n_edges[tid],
                    w_data->t_cap_edges[tid],
                    ((SegtraceDupEdge){wa, wb, d.distance}));
          }
        }
      }
    }
    i = j;
  }
  free(entries);
}

/* Helper: compute distance between two windows using in-memory hashes */
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

// ==============================================================
void build_duplicate_regions(UnionFind *uf, size_t num_sketches, int num_files,
                             char **files, GenomeSeqLen *seq_lens,
                             WindowCoord *coords,
                             SegtraceDupRegion **out_regions,
                             size_t *out_n_regions) {
  /* O(1) genome -> file_id lookup via hash map (replaces O(N*F) loop) */
  khash_t(genome_map) *gmap = kh_init(genome_map);
  for (int f = 0; f < num_files; f++) {
    char bname[256];
    get_basename(files[f], bname, sizeof(bname));
    int ret;
    khiter_t k = kh_put(genome_map, gmap, strdup(bname), &ret);
    kh_val(gmap, k) = (uint32_t)f;
  }

  uint32_t *genome_id = malloc(num_sketches * sizeof(uint32_t));
  for (size_t i = 0; i < num_sketches; i++) {
    const char *gname = seq_lens[coords[i].seq_id].genome;
    khiter_t k = kh_get(genome_map, gmap, gname);
    genome_id[i] = (k != kh_end(gmap)) ? kh_val(gmap, k) : 0;
  }

  for (khiter_t k = kh_begin(gmap); k != kh_end(gmap); ++k) {
    if (kh_exist(gmap, k))
      free((char *)kh_key(gmap, k));
  }
  kh_destroy(genome_map, gmap);
  free(genome_id);

  /* Compact family IDs: num_sketches -> n_families (typically 100x smaller) */
  uint32_t *fam_id = calloc(num_sketches, sizeof(uint32_t));
  uint32_t n_families = 0;
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = find_unionfind(uf, (uint32_t)i);
    if (fam_id[fam] == 0)
      fam_id[fam] = ++n_families;
  }

  uint8_t *final_is_sd = calloc(n_families + 1, sizeof(uint8_t));
  char **hub_label = calloc(n_families + 1, sizeof(char *));
  uint32_t next_cluster_id = 1;
  uint32_t *comp_size = calloc(n_families + 1, sizeof(uint32_t));

  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = find_unionfind(uf, (uint32_t)i);
    uint32_t fid = fam_id[fam];
    comp_size[fid]++;
  }

  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = find_unionfind(uf, (uint32_t)i);
    uint32_t fid = fam_id[fam];
    if (comp_size[fid] >= 2) {
      final_is_sd[fid] = 1;
      if (!hub_label[fid]) {
        char buf[64];
        snprintf(buf, sizeof(buf), "%u", next_cluster_id++);
        hub_label[fid] = strdup(buf);
      }
    }
  }

  size_t n_dup_regions = 0, cap_dup_regions = 0;
  SegtraceDupRegion *dup_regions = NULL;
  for (size_t i = 0; i < num_sketches; i++) {
    uint32_t fam = find_unionfind(uf, (uint32_t)i);
    uint32_t fid = fam_id[fam];
    if (!final_is_sd[fid])
      continue;

    const char *label = hub_label[fid] ? hub_label[fid] : "unknown";

    char chrom_name[512];
    snprintf(chrom_name, sizeof(chrom_name), "%s-%s",
             seq_lens[coords[i].seq_id].genome, seq_lens[coords[i].seq_id].seq);

    DA_PUSH(dup_regions, n_dup_regions, cap_dup_regions,
            ((SegtraceDupRegion){.chrom = strdup(chrom_name),
                                 .start = coords[i].start,
                                 .end = coords[i].end,
                                 .cluster_id = strdup(label),
                                 .copy_count = comp_size[fid],
                                 .subcluster_id = 0,
                                 .flank_sketch = {0},
                                 .window_idx = coords[i].window_idx}));
  }
  free(final_is_sd);
  free(comp_size);
  free(fam_id);

  for (uint32_t i = 0; i <= n_families; i++) {
    if (hub_label[i])
      free(hub_label[i]);
  }
  free(hub_label);

  *out_regions = dup_regions;
  *out_n_regions = n_dup_regions;
}

/* Merge adjacent/overlapping regions in the same SD family.
 * Returns the new count of merged regions. */
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

/* Merge adjacent/overlapping regions in the same chromosome.
 * Returns the new count of merged regions. */
size_t merge_dup_regions(SegtraceDupRegion *regions, size_t n,
                         uint32_t adjacency_threshold) {
  if (n <= 1)
    return n;

  /* 1. Sort by position (chrom, start, end) */
  qsort(regions, n, sizeof(SegtraceDupRegion), compare_dup_region_by_pos);

  /* 2. Merge all overlapping/adjacent regions on the same chromosome & keep min
   * cluster_id */
  size_t out = 0;
  for (size_t i = 1; i < n; i++) {
    if (strcmp(regions[i].chrom, regions[out].chrom) == 0 &&
        (regions[i].window_idx <=
             regions[out].window_idx + adjacency_threshold ||
         regions[i].start <= regions[out].end)) {
      /* Expand boundary */
      if (regions[i].end > regions[out].end)
        regions[out].end = regions[i].end;
      if (regions[i].window_idx > regions[out].window_idx)
        regions[out].window_idx = regions[i].window_idx;

      /* Pick minimum cluster_id in the overlapping bundle */
      uint32_t c_out = (uint32_t)strtoul(regions[out].cluster_id, NULL, 10);
      uint32_t c_i = (uint32_t)strtoul(regions[i].cluster_id, NULL, 10);
      if (c_i < c_out) {
        free(regions[out].cluster_id);
        char buf[64];
        snprintf(buf, sizeof(buf), "%u", c_i);
        regions[out].cluster_id = strdup(buf);
      }

      free(regions[i].cluster_id);
      free(regions[i].chrom);
    } else {
      out++;
      if (out != i)
        regions[out] = regions[i];
    }
  }

  /* Return all merged regions without dropping any */
  return out + 1;
}

void extract_flankings(char **files, int num_files, const Segtrace *r,
                       uint64_t scale, SegtraceDupRegion *regions,
                       size_t n_regions, int n_threads, size_t flank_size) {
  FlankingWorkerData w = {files, r, scale, regions, n_regions, flank_size};
  kt_for(n_threads, extract_flankings_worker, &w, num_files);
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

    for (size_t i = 0; i < w->n_regions; i++) {
      if (strcmp(w->regions[i].chrom, chr_name) == 0) {
        size_t start = w->regions[i].start;
        size_t end = w->regions[i].end;
        size_t flank_size = w->flank_size;
        size_t left_start = start > flank_size ? start - flank_size : 0;
        size_t right_end =
            end + flank_size > ks->seq.l ? ks->seq.l : end + flank_size;

        size_t left_len = start - left_start;
        size_t right_len = right_end - end;

        /* Free previous flanking sketch if being overwritten */
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
  }
  kseq_destroy(ks);
  gzclose(fp);
}

static int compare_dup_region_by_cluster(const void *a, const void *b);
static inline uint64_t splitmix64(uint64_t x);

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
                           double max_dist, int n_threads, uint32_t kmer_size) {
  if (n_merged <= 1)
    return;

  /* Group regions into contiguous spans of the same cluster_id */
  qsort(regions, n_merged, sizeof(SegtraceDupRegion),
        compare_dup_region_by_cluster);

  ClusterSpan *spans = NULL;
  size_t n_spans = 0, cap_spans = 0;

  size_t i = 0;
  while (i < n_merged) {
    size_t j = i + 1;
    while (j < n_merged &&
           strcmp(regions[i].cluster_id, regions[j].cluster_id) == 0) {
      j++;
    }
    DA_PUSH(spans, n_spans, cap_spans, ((ClusterSpan){i, j - i}));
    i = j;
  }

  UnionFind sub_uf;
  init_unionfind(&sub_uf, n_merged);

  SubclusterData w;
  w.regions = regions;
  w.max_dist = max_dist;
  w.n_merged = n_merged;
  w.kmer_size = kmer_size;
  w.spans = spans;
  w.t_pairs = calloc(n_threads, sizeof(SubclusterPair *));
  w.t_n_pairs = calloc(n_threads, sizeof(size_t));
  w.t_cap_pairs = calloc(n_threads, sizeof(size_t));

  kt_for(n_threads, process_subcluster, &w, n_spans);

  for (int t = 0; t < n_threads; t++) {
    for (size_t k = 0; k < w.t_n_pairs[t]; k++) {
      union_unionfind(&sub_uf, w.t_pairs[t][k].i, w.t_pairs[t][k].j);
    }
    if (w.t_pairs[t])
      free(w.t_pairs[t]);
  }
  free(w.t_pairs);
  free(w.t_n_pairs);
  free(w.t_cap_pairs);
  free(spans);

  /* Map union-find parents to subcluster_id */
  uint32_t *mapping = calloc(n_merged, sizeof(uint32_t));
  uint32_t current_id = 1;

  for (size_t k = 0; k < n_merged; k++) {
    uint32_t p = find_unionfind(&sub_uf, (uint32_t)k);
    if (mapping[p] == 0) {
      mapping[p] = current_id++;
    }
    regions[k].subcluster_id = mapping[p];
  }
  free(mapping);
  free_unionfind(&sub_uf);
}

void process_subcluster(void *data, long s, int tid) {
  SubclusterData *w = (SubclusterData *)data;
  size_t start = w->spans[s].start;
  size_t count = w->spans[s].count;

  if (count <= 1)
    return;

  if (count <= 64) {
    /* Fast all-pairs direct comparison for small clusters */
    for (size_t a = 0; a < count; a++) {
      size_t ra = start + a;
      if (w->regions[ra].flank_sketch.sketch_size == 0)
        continue;
      for (size_t b = a + 1; b < count; b++) {
        size_t rb = start + b;
        if (w->regions[rb].flank_sketch.sketch_size == 0)
          continue;

        SegtraceDistResult d =
            calculate_segtrace_dist(&w->regions[ra].flank_sketch,
                                    &w->regions[rb].flank_sketch, w->kmer_size);
        if (d.distance < w->max_dist) {
          DA_PUSH(w->t_pairs[tid], w->t_n_pairs[tid], w->t_cap_pairs[tid],
                  ((SubclusterPair){(uint32_t)ra, (uint32_t)rb}));
        }
      }
    }
  } else {
    /* Inverted k-mer index candidate extraction for large clusters */
    size_t total_flank_hashes = 0;
    for (size_t a = 0; a < count; a++) {
      total_flank_hashes += w->regions[start + a].flank_sketch.sketch_size;
    }

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

    /* Small 1MB Bloom filter to deduplicate candidate pairs within this cluster
     */
    uint32_t bloom_bits = (1 << 23);
    uint32_t mask = bloom_bits - 1;
    uint8_t *bloom = calloc(bloom_bits / 8, 1);
    if (!bloom) {
      free(entries);
      return;
    }

    size_t i = 0;
    while (i < n_entries) {
      size_t j = i + 1;
      while (j < n_entries && entries[j].hash == entries[i].hash)
        j++;
      size_t run_len = j - i;
      if (run_len >= 2 && run_len <= 1000) {
        for (size_t a = i; a < j; a++) {
          for (size_t b = a + 1; b < j; b++) {
            uint32_t la = entries[a].local_idx;
            uint32_t lb = entries[b].local_idx;
            if (la == lb)
              continue;

            uint64_t pk = encode_pair(la, lb);
            uint64_t h = splitmix64(pk);
            uint32_t h1 = (uint32_t)h & mask;
            uint32_t h2 = (uint32_t)(h >> 32) & mask;
            int was_set = ((bloom[h1 >> 3] >> (h1 & 7)) & 1) &
                          ((bloom[h2 >> 3] >> (h2 & 7)) & 1);
            bloom[h1 >> 3] |= (uint8_t)(1 << (h1 & 7));
            bloom[h2 >> 3] |= (uint8_t)(1 << (h2 & 7));
            if (was_set)
              continue;

            size_t ra = start + la;
            size_t rb = start + lb;
            SegtraceDistResult d = calculate_segtrace_dist(
                &w->regions[ra].flank_sketch, &w->regions[rb].flank_sketch,
                w->kmer_size);
            if (d.distance < w->max_dist) {
              DA_PUSH(w->t_pairs[tid], w->t_n_pairs[tid], w->t_cap_pairs[tid],
                      ((SubclusterPair){(uint32_t)ra, (uint32_t)rb}));
            }
          }
        }
      }
      i = j;
    }
    free(bloom);
    free(entries);
  }
}

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

  uint32_t max_subcluster = 0;
  for (size_t i = 0; i < n_merged; i++) {
    fprintf(out_bed, "%s\t%zu\t%zu\t%s\t%u\n", dup_regions[i].chrom,
            dup_regions[i].start, dup_regions[i].end, dup_regions[i].cluster_id,
            dup_regions[i].subcluster_id);

    if (dup_regions[i].subcluster_id > max_subcluster)
      max_subcluster = dup_regions[i].subcluster_id;
  }
  fclose(out_bed);
}

typedef struct {
  const char *c1;
  size_t s1, e1;
  const char *c2;
  size_t s2, e2;
} BedpePair;

static int compare_dup_region_by_cluster(const void *a, const void *b) {
  const SegtraceDupRegion *ra = (const SegtraceDupRegion *)a,
                          *rb = (const SegtraceDupRegion *)b;
  int c_id = strcmp(ra->cluster_id, rb->cluster_id);
  if (c_id != 0)
    return c_id;
  int c_chr = strcmp(ra->chrom, rb->chrom);
  if (c_chr != 0)
    return c_chr;
  if (ra->start != rb->start)
    return CMP(ra->start, rb->start);
  return CMP(ra->end, rb->end);
}

static int compare_bedpe_pair(const void *a, const void *b) {
  const BedpePair *pa = (const BedpePair *)a, *pb = (const BedpePair *)b;
  int c;
  c = strcmp(pa->c1, pb->c1);
  if (c != 0)
    return c;
  if (pa->s1 != pb->s1)
    return CMP(pa->s1, pb->s1);
  if (pa->e1 != pb->e1)
    return CMP(pa->e1, pb->e1);
  c = strcmp(pa->c2, pb->c2);
  if (c != 0)
    return c;
  if (pa->s2 != pb->s2)
    return CMP(pa->s2, pb->s2);
  return CMP(pa->e2, pb->e2);
}

void write_dup_bedpe(const char *out_prefix, SegtraceDupRegion *dup_regions,
                     size_t n_merged) {
  if (n_merged == 0)
    return;

  SegtraceDupRegion *regions_copy =
      malloc(n_merged * sizeof(SegtraceDupRegion));
  memcpy(regions_copy, dup_regions, n_merged * sizeof(SegtraceDupRegion));
  qsort(regions_copy, n_merged, sizeof(SegtraceDupRegion),
        compare_dup_region_by_cluster);

  BedpePair *pairs = NULL;
  size_t n_pairs = 0, cap_pairs = 0;

  size_t i = 0;
  while (i < n_merged) {
    size_t j = i + 1;
    while (j < n_merged && strcmp(regions_copy[i].cluster_id,
                                  regions_copy[j].cluster_id) == 0) {
      j++;
    }

    for (size_t a = i; a < j; a++) {
      for (size_t b = a + 1; b < j; b++) {
        const char *c1 = regions_copy[a].chrom;
        size_t s1 = regions_copy[a].start;
        size_t e1 = regions_copy[a].end;

        const char *c2 = regions_copy[b].chrom;
        size_t s2 = regions_copy[b].start;
        size_t e2 = regions_copy[b].end;

        int swap = 0;
        int cmp_chrom = strcmp(c1, c2);
        if (cmp_chrom > 0)
          swap = 1;
        else if (cmp_chrom == 0) {
          if (s1 > s2)
            swap = 1;
          else if (s1 == s2 && e1 > e2)
            swap = 1;
        }

        if (swap) {
          DA_PUSH(pairs, n_pairs, cap_pairs,
                  ((BedpePair){c2, s2, e2, c1, s1, e1}));
        } else {
          DA_PUSH(pairs, n_pairs, cap_pairs,
                  ((BedpePair){c1, s1, e1, c2, s2, e2}));
        }
      }
    }
    i = j;
  }
  free(regions_copy);

  if (n_pairs > 0) {
    qsort(pairs, n_pairs, sizeof(BedpePair), compare_bedpe_pair);
  }

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.dup.bedpe", out_prefix);
  FILE *out_bedpe = fopen(path_buf, "w");
  if (!out_bedpe) {
    fprintf(stderr, "[ERROR] Cannot open output file: %s\n", path_buf);
    free(pairs);
    return;
  }

  for (size_t p = 0; p < n_pairs; p++) {
    fprintf(out_bedpe, "%s\t%zu\t%zu\t%s\t%zu\t%zu\n", pairs[p].c1, pairs[p].s1,
            pairs[p].e1, pairs[p].c2, pairs[p].s2, pairs[p].e2);
  }
  fclose(out_bedpe);
  free(pairs);
}

// ==============================================================
// 6. CORE ALGORITHMS
// ==============================================================

void init_segtrace(Segtrace *r, size_t hash_window) {
  size_t k = hash_window < 32 ? hash_window : 32; /* 2*32=64 bits */
  uint32_t kmer_bits = 2 * (uint32_t)k;

  uint64_t remover_mask =
      (kmer_bits > 2) ? (((uint64_t)1 << (kmer_bits - 2)) - 1) : 0;
  /* remover mask to forget previous base */

  r->hash_window = k;
  r->remover_mask = remover_mask;
  r->kmer_bits = kmer_bits;
  r->rc_shift =
      (kmer_bits > 0) ? (kmer_bits - 2) : 0; /* reverse_complement shift */
}

/* Extract segtrace hash and insert in HashPool */
__attribute__((hot)) void extract_hash(const Segtrace *r, HashPool *pool,
                                       const uint8_t *seq, size_t len) {
  size_t K = r->hash_window;
  size_t valid = 0;
  uint64_t fwd = 0;
  uint64_t rev = 0;
  uint64_t mask = r->remover_mask;
  uint32_t rc_shift = r->rc_shift;

  for (size_t idx = 0; idx < len; idx++) {
    int8_t lv = BASE_LOOKUP[seq[idx]];
    if (lv < 0) {
      fwd = 0;
      rev = 0;
      valid = 0;
      continue;
    }

    fwd = ((fwd & mask) << 2) | (uint8_t)lv;
    rev = (rev >> 2) | ((uint64_t)(lv ^ 3) << rc_shift);

    if (valid < K)
      valid++;
    if (valid < K)
      continue;

    uint64_t canon = fwd < rev ? fwd : rev;
    uint64_t h = mix_hash(canon, r->hash_seed);

    if (h < pool->hash_threshold)
      insert_hash_pool(pool, h);
  }
}

void init_hash_pool(HashPool *pool, uint64_t threshold) {
  *pool = (HashPool){.hash_threshold = threshold};
}

void insert_hash_pool(HashPool *pool, uint64_t h) {
  if (h >= pool->hash_threshold)
    return;
  DA_PUSH(pool->hashes, pool->size, pool->cap, h);
}

void finalize_hash_pool(HashPool *pool, uint64_t **out_hashes,
                        size_t *out_size) {
  size_t n = pool->size;
  if (n) {
    qsort(pool->hashes, n, sizeof(uint64_t), compare_uint64);
    /* Deduplicate in-place */
    size_t u = 0;
    for (size_t i = 0; i < n; i++) {
      if (u == 0 || pool->hashes[i] != pool->hashes[u - 1])
        pool->hashes[u++] = pool->hashes[i];
    }
    n = u;
    uint64_t *nh = realloc(pool->hashes, n * sizeof(uint64_t));
    if (nh)
      pool->hashes = nh;
  } else {
    /* pool had capacity but no hashes survived — free the buffer */
    free(pool->hashes);
  }
  *out_size = n;
  *out_hashes = n ? pool->hashes : NULL;
  pool->hashes = NULL;
  pool->size = pool->cap = 0;
}

/* Calculate distance between two sketch sets */
SegtraceDistResult calculate_segtrace_dist(const SegtraceSketch *ref,
                                           const SegtraceSketch *query,
                                           uint32_t kmer_size) {
  SegtraceDistResult res = {0.0, 1.0, 0};

  size_t shared = 0, i = 0, j = 0;

  /* Set operations with two pointer algorithm; ref and query are sorted. */
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

  /* Calculate distance with containment method */
  res.shared_hashes = shared;
  if (ref->sketch_size > 0 && query->sketch_size > 0) {
    res.containment =
        0.5 * shared * (1.0 / ref->sketch_size + 1.0 / query->sketch_size);
    res.distance = 1.0 - pow(res.containment, 1.0 / (double)kmer_size);
  }
  return res;
}

void init_unionfind(UnionFind *uf, size_t n) {
  uf->n = n;
  uf->parent = (uint32_t *)malloc(n * sizeof(uint32_t));
  uf->rank = (uint8_t *)calloc(n, sizeof(uint8_t));
  if (!uf->parent || !uf->rank) {
    fprintf(stderr, "[ERROR] Failed to allocate UnionFind for %zu elements\n",
            n);
    exit(1);
  }
  for (size_t i = 0; i < n; i++)
    uf->parent[i] = (uint32_t)i;
}

uint32_t find_unionfind(UnionFind *uf, uint32_t x) {
  while (uf->parent[x] != x) {
    uf->parent[x] = uf->parent[uf->parent[x]]; /* path splitting */
    x = uf->parent[x];
  }
  return x;
}

void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b) {
  a = find_unionfind(uf, a);
  b = find_unionfind(uf, b);
  if (a == b)
    return;
  if (uf->rank[a] < uf->rank[b])
    SWAP(uint32_t, a, b);
  uf->parent[b] = a;
  if (uf->rank[a] == uf->rank[b])
    uf->rank[a]++;
}

void free_unionfind(UnionFind *uf) {
  free(uf->parent);
  free(uf->rank);
}

// ==============================================================
// 7. UTILITIES
// ==============================================================

void get_basename(const char *filename, char *basename, size_t size) {
  const char *slash = strrchr(filename, '/');
  const char *base = slash ? slash + 1 : filename;
  strncpy(basename, base, size - 1);
  basename[size - 1] = '\0';
  char *dot = strchr(basename, '.');
  if (dot)
    *dot = '\0';
}

uint64_t mix_hash(uint64_t hash_value, uint64_t seed) {
  __uint128_t p =
      ((__uint128_t)hash_value ^ MIX_CONST1) * ((__uint128_t)seed ^ MIX_CONST2);
  return (uint64_t)(p ^ (p >> 64));
}

/* Binary search: first index where arr[i] >= target */
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

/* Candidate pair stored as single uint64: (min << 32) | max */
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

/* Wrapper for macro to use in `qsort`*/
int compare_uint64(const void *a, const void *b) {
  return CMP(*(const uint64_t *)a, *(const uint64_t *)b);
}

int compare_hash_entry(const void *a, const void *b) {
  const HashWindowEntry *ea = (const HashWindowEntry *)a,
                        *eb = (const HashWindowEntry *)b;
  return ea->hash != eb->hash ? CMP(ea->hash, eb->hash)
                              : CMP(ea->window_id, eb->window_id);
}

int compare_dup_region(const void *a, const void *b) {
  const SegtraceDupRegion *ra = (const SegtraceDupRegion *)a,
                          *rb = (const SegtraceDupRegion *)b;
  int c_id = strcmp(ra->cluster_id, rb->cluster_id);
  if (c_id != 0)
    return c_id;
  int c_chr = strcmp(ra->chrom, rb->chrom);
  if (c_chr != 0)
    return c_chr;
  if (ra->start != rb->start)
    return CMP(ra->start, rb->start);
  return CMP(ra->end, rb->end);
}
