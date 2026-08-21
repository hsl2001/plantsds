#ifndef SEGTRACE_H
#define SEGTRACE_H

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

#define NUM_PARTITIONS 1024
#define BATCH_PARTITIONS 128
#define MAX_KMER_FREQ 256
#define MAX_PAIR_COMPARISONS 16
#define MAX_COLLINEAR_LOOKAHEAD 8
#define MIN_SD_LEN 1000
#define MERGE_COEFF 2

#define BLOOM_SIZE_BITS (1 << 24)
#define BLOOM_SIZE_BYTES (BLOOM_SIZE_BITS / 8)
#define BLOOM_MASK (BLOOM_SIZE_BITS - 1)

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
  uint32_t *parent;
  uint8_t *rank;
  size_t n;
} UnionFind;

typedef struct {
  uint32_t seq_id;
  uint32_t file_id;
  size_t start;
  size_t end;
  uint32_t cluster_id;
} SegtraceDupRegion;

typedef struct {
  uint32_t seq_id;
  uint32_t sketch_offset;
  uint32_t window_idx;
  uint16_t sketch_size;
} WindowCoord;

typedef struct {
  char *genome;
  char *seq;
  uint32_t file_id;
} GenomeSeqLen;

typedef struct {
  uint32_t a;
  uint32_t b;
} CandidatePair;

typedef struct {
  uint32_t hash;
  uint32_t window_id;
} HashWindowEntry;

typedef struct {
  HashWindowEntry *entries;
  size_t size;
  size_t cap;
} PartitionBucket;

typedef struct {
  uint32_t *all_hashes;
  size_t num_all_hashes;
  size_t cap_all_hashes;
  WindowCoord *coords;
  size_t num_sketches;
  size_t cap_sketches;
  GenomeSeqLen *seq_lens;
  size_t num_seqs;
  size_t cap_seqs;
} GlobalWindows;

typedef struct {
  const Segtrace *r;
  const int8_t *base_lookup;
  uint32_t threshold;
  size_t window_size;
  size_t step_size;
  size_t min_bases;
  uint32_t seq_id;
  const uint8_t *seq_ptr;
  size_t seq_len;
  size_t chunk_start_idx;
  size_t chunk_end_idx;
  uint32_t *hashes;
  size_t num_hashes;
  size_t cap_hashes;
  WindowCoord *coords;
  size_t num_coords;
  size_t cap_coords;
} SeqChunkJob;

typedef struct {
  const uint32_t *all_hashes;
  const WindowCoord *coords;
  size_t n_windows;
  size_t window_size;
  size_t step_size;
  uint32_t kmer_size;
  double p_kmer;
  PartitionBucket *buckets;
  uint8_t **t_bloom;
  CandidatePair **t_pairs;
  size_t *t_n_pairs;
  size_t *t_cap_pairs;
  size_t batch_start;
} DiscoverComputeData;

extern const int8_t BASE_LOOKUP[256];
extern const int8_t BASE_LOOKUP_NO_MASK[256];

// ==============================================================
// FUNCTION DECLARATIONS
// ==============================================================
void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);

// 1. CLI & UTILITIES
void print_usage(void);
void get_basename(const char *filename, char *basename, size_t size);
uint32_t mix_hash(uint64_t hash_value, uint64_t seed);
uint64_t encode_pair(uint32_t a, uint32_t b);
int bloom_test_and_set(uint8_t *bloom, uint64_t key, uint32_t mask);
int compare_uint32(const void *a, const void *b);
int compare_hash_entry(const void *a, const void *b);

// 2. CORE & SKETCHING
void init_segtrace(Segtrace *r, size_t hash_window, int filter_masked);
void init_unionfind(UnionFind *uf, size_t n);
uint32_t find_unionfind(UnionFind *uf, uint32_t x);
void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b);
void free_unionfind(UnionFind *uf);

// 3. WINDOW EXTRACTION
GlobalWindows extract_all_windows(char **files, int num_files,
                                  const Segtrace *r, uint64_t scale,
                                  size_t window_size, size_t step_size,
                                  size_t min_bases, int n_threads);

// 4. DISCOVERY & DISTANCE COMPUTATION
void discover_and_compute(const uint32_t *all_hashes, const WindowCoord *coords,
                          size_t n_windows, size_t window_size,
                          size_t step_size, int n_threads, uint32_t kmer_size,
                          UnionFind *uf);
void discover_compute_worker(void *data, long idx, int tid);

// 5. REGION CLUSTERING, COPY FILTERING & OUTPUT
void build_duplicate_regions(UnionFind *uf, size_t num_sketches,
                             const WindowCoord *coords,
                             const GenomeSeqLen *seq_lens, size_t step_size,
                             size_t window_size,
                             SegtraceDupRegion **out_regions,
                             size_t *out_n_regions);
size_t merge_dup_regions(SegtraceDupRegion *regions, size_t n,
                         size_t window_size);
size_t filter_regions_by_copy_count(SegtraceDupRegion *regions, size_t n,
                                    uint32_t min_copies);
void write_dup_bed(const char *out_prefix, const SegtraceDupRegion *dup_regions,
                   size_t n_merged, const GenomeSeqLen *seq_lens,
                   size_t min_sd_len);

#ifdef __cplusplus
}
#endif

#endif
