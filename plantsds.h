#ifndef PLANTSDS_H
#define PLANTSDS_H

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
#define MIX_CONST1 0xff51afd7ed558ccdULL
#define MIX_CONST2 0xc4ceb9fe1a85ec53ULL

/* Generic dynamic-array capacity reserve */
#define DA_RESERVE(arr, cap, req_cap)                                          \
  do {                                                                         \
    if ((req_cap) > (cap)) {                                                   \
      (cap) = (cap) ? (cap) : 16;                                              \
      while ((cap) < (req_cap))                                                \
        (cap) *= 2;                                                            \
      (arr) = realloc((arr), (cap) * sizeof(*(arr)));                          \
    }                                                                          \
  } while (0)

/* Generic dynamic-array push */
#define DA_PUSH(arr, n, cap, val)                                              \
  do {                                                                         \
    DA_RESERVE(arr, cap, (n) + 1);                                             \
    (arr)[(n)++] = (val);                                                      \
  } while (0)

#define CMP(a, b) (((a) > (b)) - ((a) < (b)))
#define SWAP(type, a, b)                                                       \
  do {                                                                         \
    type _t = (a);                                                             \
    (a) = (b);                                                                 \
    (b) = _t;                                                                  \
  } while (0)

#define ABS_DIFF(a, b) ((a) > (b) ? (a) - (b) : (b) - (a))

#define NUM_PARTITIONS 256
#define MAX_RUN_LEN 100

#define BLOOM_SIZE_BITS (1 << 22)
#define BLOOM_SIZE_BYTES (BLOOM_SIZE_BITS / 8)
#define BLOOM_MASK (BLOOM_SIZE_BITS - 1)

// ==============================================================
// CORE DATA STRUCTURES
// ==============================================================
typedef struct {
  uint32_t hash_window;
  __uint128_t remover_mask;
  uint32_t kmer_bits;
  uint32_t rc_shift;
  uint64_t hash_seed;
} Plantsds;

typedef struct {
  size_t sketch_size;
  uint64_t *hashes;
} PlantsdsSketch;

typedef struct {
  double containment;
  double distance;
  size_t shared_hashes;
} PlantsdsDistResult;

typedef struct {
  uint32_t win_a;
  uint32_t win_b;
  double distance;
} PlantsdsDupEdge;

typedef struct {
  uint32_t *parent;
  uint8_t *rank;
  size_t n;
} UnionFind;

typedef struct {
  char *chrom;
  size_t start;
  size_t end;
  char *cluster_id;
  uint32_t copy_count;
  uint32_t subcluster_id;
  PlantsdsSketch flank_sketch;
  uint32_t window_idx;
} PlantsdsDupRegion;

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
  const char *filename;
  char bname[256];
  const Plantsds *r;
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
  const Plantsds *r;
  uint64_t scale;
  PlantsdsDupRegion *regions;
  size_t n_regions;
  size_t flank_size;
} FlankingWorkerData;

typedef struct {
  uint32_t i, j;
} SubclusterPair;

typedef struct {
  PlantsdsDupRegion *regions;
  double max_dist;
  size_t n_merged;
  SubclusterPair **t_pairs;
  size_t *t_n_pairs;
  size_t *t_cap_pairs;
  uint32_t kmer_size;
} SubclusterData;

typedef struct {
  uint64_t hash;
  uint32_t window_id;
} HashWindowEntry;

typedef struct {
  const uint64_t *all_hashes;
  WindowCoord *coords;
  size_t n_windows;
  size_t window_size;
  double max_dist;
  uint32_t kmer_size;
  PlantsdsDupEdge **t_edges;
  size_t *t_n_edges;
  size_t *t_cap_edges;
  uint8_t **t_bloom;
} DiscoverComputeData;

extern const int8_t BASE_LOOKUP[256];

// ==============================================================
// FUNCTION DECLARATIONS
// ==============================================================
void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);

// 1. ENTRY POINT & CLI
int run_dup(int argc, char **argv);
void print_usage(void);

// 2. PIPELINE ORCHESTRATION
int run_pangenome(int num_files, char **files, size_t flank_size,
                  const Plantsds *r, uint64_t scale, size_t window_size,
                  size_t step_size, size_t min_bases, double max_dist,
                  int min_copy, int max_copy, const char *out_prefix,
                  int n_threads, uint32_t adjacency_threshold,
                  double subcluster_dist);

// 3. WINDOW EXTRACTION
StreamWorkerData *extract_all_windows(char **files, int num_files,
                                      const Plantsds *r, uint64_t scale,
                                      size_t window_size, size_t step_size,
                                      size_t min_bases, int n_threads);
void stream_pangenome_worker(void *data, long i, int tid);
void merge_global_data(StreamWorkerData *workers, int num_files,
                       const char *out_prefix, uint64_t **out_all_hashes,
                       WindowCoord **out_coords, size_t *out_num_sketches,
                       GenomeSeqLen **out_seq_lens, size_t *out_num_seqs);

// 4. DISCOVERY & DISTANCE CALCULATION
void discover_and_compute(const uint64_t *all_hashes, WindowCoord *coords,
                          size_t n_windows, size_t window_size, double max_dist,
                          int n_threads, uint32_t kmer_size, UnionFind *uf);
void discover_compute_worker(void *data, long p, int tid);
PlantsdsDistResult calculate_window_dist(const uint64_t *all_hashes,
                                         const WindowCoord *wa,
                                         const WindowCoord *wb,
                                         uint32_t kmer_size);

// 5. REGION CLUSTERING & OUTPUT
void build_duplicate_regions(UnionFind *uf, size_t num_sketches, int num_files,
                             char **files, GenomeSeqLen *seq_lens,
                             WindowCoord *coords, int min_copy, int max_copy,
                             PlantsdsDupRegion **out_regions,
                             size_t *out_n_regions);
size_t merge_dup_regions(PlantsdsDupRegion *regions, size_t n,
                         uint32_t adjacency_threshold);
void extract_flankings(char **files, int num_files, const Plantsds *r,
                       uint64_t scale, PlantsdsDupRegion *regions,
                       size_t n_regions, size_t flank_size, int n_threads);
void extract_flankings_worker(void *data, long f, int tid);
void perform_subclustering(PlantsdsDupRegion *regions, size_t n_merged,
                           double max_dist, int n_threads, uint32_t kmer_size);
void process_subcluster(void *data, long i, int tid);
void write_dup_bed(const char *out_prefix, PlantsdsDupRegion *dup_regions,
                   size_t n_merged);

// 6. CORE ALGORITHMS
void init_plantsds(Plantsds *r, size_t hash_window);
void extract_hash(const Plantsds *r, HashPool *pool, const uint8_t *seq,
                  size_t len);
void init_hash_pool(HashPool *pool, uint64_t threshold);
void insert_hash_pool(HashPool *pool, uint64_t h);
void finalize_hash_pool(HashPool *pool, uint64_t **out_hashes,
                        size_t *out_size);
PlantsdsDistResult calculate_plantsds_dist(const PlantsdsSketch *ref,
                                           const PlantsdsSketch *query,
                                           uint32_t kmer_size);
void init_unionfind(UnionFind *uf, size_t n);
uint32_t find_unionfind(UnionFind *uf, uint32_t x);
void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b);
void free_unionfind(UnionFind *uf);

// 7. UTILITIES
void get_basename(const char *filename, char *basename, size_t size);
uint64_t mix_hash(__uint128_t hash_value, uint64_t seed);
uint64_t reverse_bits64(uint64_t n);
__uint128_t reverse_bits128(__uint128_t n);
size_t lower_bound_u64(const uint64_t *arr, size_t n, uint64_t target);
uint64_t encode_pair(uint32_t a, uint32_t b);
int bloom_test_and_set(uint8_t *bloom, uint64_t key);
int compare_uint64(const void *a, const void *b);
int compare_hash_entry(const void *a, const void *b);
int compare_dup_region(const void *a, const void *b);

#ifdef __cplusplus
}
#endif

#endif
