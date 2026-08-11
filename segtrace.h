#ifndef SEGTRACE_H
#define SEGTRACE_H

#include <pthread.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#ifdef __cplusplus
extern "C" {
#endif

// ==============================================================
// CONSTANTS & MACROS
// ==============================================================
#define MIX_CONST1 0xa0761d6478bd642fULL
#define MIX_CONST2 0xe7037ed1a0b428dbULL

/* Generic dynamic-array capacity reserve */
#define DA_RESERVE(arr, cap, req_cap)                                          \
  do {                                                                         \
    if ((req_cap) > (cap)) {                                                   \
      (cap) = (cap) ? (cap) : 16;                                              \
      while ((cap) < (req_cap))                                                \
        (cap) *= 2;                                                            \
      void *tmp_da_ptr = realloc((arr), (cap) * sizeof(*(arr)));               \
      if (!tmp_da_ptr) {                                                       \
        fprintf(stderr, "[ERROR] Out of memory in DA_RESERVE\n");              \
        exit(1);                                                               \
      }                                                                        \
      (arr) = tmp_da_ptr;                                                      \
    }                                                                          \
  } while (0)

/* Generic dynamic-array push */
#define DA_PUSH(arr, n, cap, val)                                              \
  do {                                                                         \
    DA_RESERVE(arr, cap, (n) + 1);                                             \
    (arr)[(n)++] = (val);                                                      \
  } while (0)

#define CMP(a, b) (((a) > (b)) - ((a) < (b)))

#define ABS_DIFF(a, b) ((a) > (b) ? (a) - (b) : (b) - (a))

#define NUM_PARTITIONS 512
#define MAX_KMER_FREQ 32
#define MAX_PAIR_COMPARISONS 64
#define MAX_COLLINEAR_LOOOKAHEAD 8
#define MIN_SD_LEN 1000
#define MERGE_COEFF 10

#define BLOOM_SIZE_BITS (1 << 24)
#define BLOOM_SIZE_BYTES (BLOOM_SIZE_BITS / 8)
#define BLOOM_MASK (BLOOM_SIZE_BITS - 1)

#define SUBCLUSTER_BLOOM_SIZE_BITS (BLOOM_SIZE_BITS / 32)
#define SUBCLUSTER_BLOOM_SIZE_BYTES (SUBCLUSTER_BLOOM_SIZE_BITS / 8)
#define SUBCLUSTER_BLOOM_MASK (SUBCLUSTER_BLOOM_SIZE_BITS - 1)

// ==============================================================
// CORE DATA STRUCTURES
// ==============================================================
typedef struct {
  uint32_t hash_window;
  uint64_t hash_seed;
  int filter_masked;
  const int8_t *base_lookup;
} Segtrace;

typedef struct {
  size_t sketch_size;
  uint64_t *hashes;
} SegtraceSketch;

typedef struct {
  size_t shared_hashes;
} SegtraceDistResult;

typedef struct {
  uint32_t *parent;
  uint8_t *rank;
  size_t n;
  pthread_mutex_t lock;
} UnionFind;

typedef struct {
  char *chrom;
  size_t start;
  size_t end;
  char *cluster_id;
  uint32_t subcluster_id;
  SegtraceSketch flank_sketch;
  uint32_t window_idx;
} SegtraceDupRegion;

// ==============================================================
// INTERNAL DATA STRUCTURES
// ==============================================================
typedef struct {
  size_t size;
  size_t cap;
  uint64_t hash_threshold;
  uint64_t *hashes;
} HashPool;

typedef struct {
  uint32_t seq_id;
  size_t start;
  size_t end;
  size_t sketch_offset;
  uint32_t sketch_size;
  uint32_t window_idx;
} WindowCoord;

typedef struct {
  char *genome;
  char *seq;
} GenomeSeqLen;

typedef struct {
  const Segtrace *r;
  uint64_t scale;
  size_t window_size;
  size_t step_size;
  size_t min_bases;
  uint64_t *all_hashes;
  size_t num_all_hashes;
  size_t cap_all_hashes;
  WindowCoord *coords;
  size_t num_sketches;
  size_t cap_sketches;
  GenomeSeqLen *seq_lens;
  size_t num_seqs;
  size_t cap_seqs;
} StreamWorkerData;

typedef struct {
  char **files;
  const Segtrace *r;
  uint64_t scale;
  SegtraceDupRegion *regions;
  size_t n_regions;
  size_t flank_size;
} FlankingWorkerData;

typedef struct {
  uint32_t i, j;
} SubclusterPair;

typedef struct {
  size_t start;
  size_t count;
} ClusterSpan;

typedef struct {
  SegtraceDupRegion *regions;
  size_t n_merged;
  SubclusterPair **t_pairs;
  size_t *t_n_pairs;
  size_t *t_cap_pairs;
  uint32_t kmer_size;
  ClusterSpan *spans;
  uint8_t **t_bloom;
} SubclusterData;

typedef struct {
  const Segtrace *r;
  const int8_t *base_lookup;
  uint64_t scale;
  size_t window_size;
  size_t step_size;
  size_t min_bases;
  uint32_t seq_id;
  const uint8_t *seq_ptr;
  size_t seq_len;
  size_t chunk_start_idx;
  size_t chunk_end_idx;
  uint64_t *hashes;
  size_t num_hashes;
  size_t cap_hashes;
  WindowCoord *coords;
  size_t num_coords;
  size_t cap_coords;
} SeqChunkJob;

typedef struct {
  uint64_t hash;
  uint32_t window_id;
} HashWindowEntry;

typedef struct {
  HashWindowEntry *entries;
  size_t size;
  size_t cap;
} PartitionBucket;

typedef struct {
  const uint64_t *all_hashes;
  WindowCoord *coords;
  size_t n_windows;
  size_t window_size;
  uint32_t kmer_size;
  UnionFind *uf;
  PartitionBucket *buckets;
  uint8_t **t_bloom;
} DiscoverComputeData;

extern const int8_t BASE_LOOKUP[256];

// ==============================================================
// FUNCTION DECLARATIONS
// ==============================================================
void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);

// 1. ENTRY POINT & CLI
void print_usage(void);

// 3. WINDOW EXTRACTION
StreamWorkerData *extract_all_windows(char **files, int num_files,
                                      const Segtrace *r, uint64_t scale,
                                      size_t window_size, size_t step_size,
                                      size_t min_bases, int n_threads);
void merge_global_data(StreamWorkerData *workers, int num_files,
                       const char *out_prefix, uint64_t **out_all_hashes,
                       WindowCoord **out_coords, size_t *out_num_sketches,
                       GenomeSeqLen **out_seq_lens, size_t *out_num_seqs);

// 4. DISCOVERY & DISTANCE CALCULATION
void discover_and_compute(const uint64_t *all_hashes, WindowCoord *coords,
                          size_t n_windows, size_t window_size, int n_threads,
                          uint32_t kmer_size, UnionFind *uf);
void discover_compute_worker(void *data, long p, int tid);
SegtraceDistResult calculate_window_dist(const uint64_t *all_hashes,
                                         const WindowCoord *wa,
                                         const WindowCoord *wb,
                                         uint32_t kmer_size);

// 5. REGION CLUSTERING & OUTPUT
void build_duplicate_regions(UnionFind *uf, size_t num_sketches,
                             GenomeSeqLen *seq_lens, WindowCoord *coords,
                             SegtraceDupRegion **out_regions,
                             size_t *out_n_regions);
size_t merge_dup_regions(SegtraceDupRegion *regions, size_t n,
                         size_t window_size);
void extract_flankings(char **files, int num_files, const Segtrace *r,
                       uint64_t scale, SegtraceDupRegion *regions,
                       size_t n_regions, int n_threads, size_t flank_size);
void extract_flankings_worker(void *data, long f, int tid);
void perform_subclustering(SegtraceDupRegion *regions, size_t n_merged,
                           int n_threads, uint32_t kmer_size);
void process_subcluster(void *data, long i, int tid);
void write_dup_bed(const char *out_prefix, SegtraceDupRegion *dup_regions,
                   size_t n_merged);

// 6. CORE ALGORITHMS
void init_segtrace(Segtrace *r, size_t hash_window, int filter_masked);
void extract_hash(const Segtrace *r, HashPool *pool, const uint8_t *seq,
                  size_t len);
void init_hash_pool(HashPool *pool, uint64_t threshold);
void insert_hash_pool(HashPool *pool, uint64_t h);
void finalize_hash_pool(HashPool *pool, uint64_t **out_hashes,
                        size_t *out_size);
SegtraceDistResult calculate_segtrace_dist(const SegtraceSketch *ref,
                                           const SegtraceSketch *query);
void init_unionfind(UnionFind *uf, size_t n);
uint32_t find_unionfind(UnionFind *uf, uint32_t x);
void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b);
void free_unionfind(UnionFind *uf);

// 7. UTILITIES
void get_basename(const char *filename, char *basename, size_t size);
uint64_t mix_hash(uint64_t hash_value, uint64_t seed);
uint64_t encode_pair(uint32_t a, uint32_t b);
int bloom_test_and_set(uint8_t *bloom, uint64_t key, uint32_t mask);
int compare_uint64(const void *a, const void *b);
int compare_hash_entry(const void *a, const void *b);

#ifdef __cplusplus
}
#endif

#endif
