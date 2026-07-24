#ifndef PLANTSDS_H
#define PLANTSDS_H

#include <stddef.h>
#include <stdint.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
#ifdef __cplusplus
extern "C" {
#endif

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

/* (HashWindowEntry is defined internally in plantsds.c) */

/* Edge in the duplication graph */
typedef struct {
  uint32_t win_a;
  uint32_t win_b;
  double distance;
} PlantsdsDupEdge;

typedef struct {
  uint32_t *parent;
  uint32_t *rank;
  size_t n;
} UnionFind;

/* Merged SD region */
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

void init_plantsds(Plantsds *r, size_t hash_window);
PlantsdsDistResult calculate_plantsds_dist(const PlantsdsSketch *ref,
                                           const PlantsdsSketch *query,
                                           uint32_t kmer_size);

/* Union-Find operations */
void init_unionfind(UnionFind *uf, size_t n);
uint32_t find_unionfind(UnionFind *uf, uint32_t x);
void union_unionfind(UnionFind *uf, uint32_t a, uint32_t b);
void free_unionfind(UnionFind *uf);

#ifdef __cplusplus
}
#endif

#endif
