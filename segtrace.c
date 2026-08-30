#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "klib/ketopt.h"
#include "klib/kseq.h"
#include "segtrace.h"

/* read 모드(-r) 플래그. extract_all_windows가 윈도우 대신 read 전체를
 * 스케치하도록 전환한다 */
int g_segtrace_read_mode = 0;

/* kseq 리더 인스턴스화: gzread 기반이라 일반/gzip FASTA 모두 읽을 수 있다 */
KSEQ_INIT(gzFile, gzread)

// ==============================================================
// SECTION 1: ENTRY POINT & CLI PARSING
// ==============================================================

/* 명령행 사용법과 각 옵션의 기본값을 출력한다 */
static void print_usage(void) {
  printf("Segtrace: Segmental Tracer\n\n"
         "Usage: segtrace [options] fasta1 [fasta2 ...]\n\n"
         "Options:\n"
         "  -k: kmer size (default: 19)\n"
         "  -s: scale factor (default: 16)\n"
         "  -w: window size in bp (default: 1024)\n"
         "  -t: step size in bp (default: 0 [auto: 33%% of window size])\n"
         "  -b: minimum valid bases per window (default: 0 [auto: 25%% of "
         "window size])\n"
         "  -c: minimum copies per genome/file to report (default: 1)\n"
         "  -r: read mode - each FASTA/FASTQ record is one segment; window "
         "and step are set automatically (window = whole record)\n"
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
  /* --help는 위치와 상관없이 동작하도록 전체 인자를 먼저 훑는다 */
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      print_usage();
      return 0;
    }
  }

  /* 파라미터 기본값:
   * kmer_size  = ntHash k-mer 길이
   * scale      = 스케치 축소율 (해시값 하위 1/scale만 샘플링)
   * window/step/min_bases = 윈도우 크기, 이동 간격, 최소 유효 염기 수
   *   (0이면 윈도우 크기에서 자동 유도)
   * min_copies = 파일(유전체)당 보고할 최소 복제 수 */
  uint32_t kmer_size = 19;
  uint64_t scale = 16;
  size_t window_size = 1024, step_size = 0, min_bases = 0;
  uint32_t min_copies = 1;
  const char *out_prefix = "segtrace";
  int n_threads = 8, filter_masked = 0, read_mode = 0;

  /* 단일 대시 옵션 파싱 (ketopt: getopt의 경량 대체) */
  ketopt_t opt = KETOPT_INIT;
  int c;
  while ((c = ketopt(&opt, argc, argv, 1, "k:s:w:t:b:c:o:p:rmh", 0)) >= 0) {
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
    } else if (c == 'm')
      filter_masked = 1;
    else if (c == 'r')
      read_mode = 1;
    else
      return 1;
  }
  g_segtrace_read_mode = read_mode;
  /* read 모드에서는 윈도우/스텝을 사용자가 신경 쓸 필요가 없다:
   * window_size = 0이 "레코드 전체가 윈도우 하나"라는 뜻의 내부 신호이고,
   * step_size는 read 경계가 곧 윈도우 경계이므로 쓰이지 않는다.
   * -w/-t를 함께 줘도 read 모드에서는 무시한다 */
  if (read_mode) {
    window_size = 0;
    step_size = 0;
  }
  /* 자동 파라미터 결정:
   * step_size는 윈도우의 1/3 (인접 윈도우가 3배 중첩),
   * min_bases는 윈도우의 1/4 (N이 75% 이상인 윈도우는 스케치하지 않음).
   * read 모드에서는 step이 쓰이지 않으므로 건너뛴다 */
  if (!read_mode && step_size == 0)
    step_size = window_size / 3;
  if (min_bases == 0)
    min_bases = read_mode ? 1 : window_size / 4;
  if (opt.ind == argc) {
    fprintf(stderr, "[ERROR] Input FASTA files are required.\n");
    return 1;
  }

  int num_files = argc - opt.ind;
  char **files = &argv[opt.ind];

  /* 염기 -> 2bit 코드(A=0,C=1,G=2,T=3) 룩업 테이블 구성.
   * -1은 N 등 유효하지 않은 염기. -m 옵션이 없으면 소문자(soft-masked)도
   * 유효 염기로 취급한다. hash_seed=42는 스케치 재현성을 위한 고정 시드. */
  Segtrace r = {.hash_window = kmer_size, .hash_seed = 42};
  memset(r.base_lookup, -1, sizeof(r.base_lookup));
  for (int8_t code = 0; code < 4; code++) {
    uint8_t base = (uint8_t)"ACGT"[code];
    r.base_lookup[base] = code;
    if (!filter_masked)
      r.base_lookup[base + ('a' - 'A')] = code;
  }

  void *thread_pool = n_threads > 1 ? kt_forpool_init(n_threads) : NULL;

  /* [1단계] 모든 FASTA를 윈도우(read 모드에서는 read) 단위로 나누고 각
   * 단위의 스케치(해시 집합) 추출 */
  GlobalWindows gw =
      extract_all_windows(files, num_files, &r, scale, window_size,
                          step_size, min_bases, n_threads, thread_pool);

  if (read_mode) {
    /* read 모드: 동일 스케치를 가진 read를 그룹화해 그룹 크기(관측 개수)를
     * 구하고, 크기-1 그룹의 봉우리로 배수체(haploid) 커버리지를 추정한다.
     * 세그먼트 복제수 = round(그룹 크기 * scale / C1):
     * 스케치당 기대 개수는 ~lambda/s 이므로 실제 깊이 lambda = 크기 * scale */
    ReadGroupStat *groups = NULL;
    size_t n_groups = group_identical_reads(gw.all_hashes, gw.coords,
                                            gw.num_sketches, n_threads,
                                            thread_pool, &groups);
    if (thread_pool)
      kt_forpool_destroy(thread_pool);

    double hap_cov =
        estimate_haploid_coverage(groups, n_groups, scale, window_size);
    fprintf(stderr,
            "[segtrace] Read mode: %zu windows, %zu distinct segments, "
            "haploid coverage ~ %.2fx\n",
            gw.num_sketches, n_groups, hap_cov);
    write_read_seg_bed(out_prefix, gw.coords, gw.seq_lens, groups, n_groups,
                       scale, window_size, hap_cov, min_copies);

    free(groups);
    free(gw.all_hashes);
    free(gw.coords);
    for (size_t i = 0; i < gw.num_seqs; i++) {
      free(gw.seq_lens[i].genome);
      free(gw.seq_lens[i].seq);
    }
    free(gw.seq_lens);
    return 0;
  }

  /* [2단계] 해시를 공유하는 윈도우 쌍을 후보로 찾고 유사도(공유 해시 수) 계산 */
  fprintf(stderr,
          "[segtrace] Discovering candidates and computing distances...\n");
  CandidateGraph graph =
      discover_and_compute(gw.all_hashes, gw.coords, gw.num_sketches,
                 window_size, step_size, n_threads, r.hash_window,
                           thread_pool);
  if (thread_pool)
    kt_forpool_destroy(thread_pool);

  free(gw.all_hashes); /* 이후 단계에서는 원본 해시 배열이 필요 없음 */

  /* [3단계] 후보에 포함된 윈도우를 같은 서열 내에서 연속 구간(locus)으로 병합 */
  SegtraceDupRegion *dup_regions = NULL;
  size_t n_dup_regions = 0;
  build_duplicate_loci(&graph, gw.num_sketches, gw.coords, gw.seq_lens,
                       step_size, window_size, &dup_regions, &n_dup_regions);
  /* [4단계] 유사 구간 쌍을 union-find로 묶어 클러스터(복제 패밀리) id 부여 */
  cluster_duplicate_loci(&graph, gw.coords, dup_regions, n_dup_regions);
  free_candidate_graph(&graph);

  free(gw.coords);

  /* [5단계] 유전체당 복제 수(min_copies) 미만인 클러스터 그룹 제거 */
  size_t n_filtered =
      filter_regions_by_copy_count(dup_regions, n_dup_regions, min_copies);

  /* [6단계] BED 형식으로 출력 (최소 SD 길이 미만 구간은 제외) */
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

/* 64비트 순환 좌회전. (64 - n) & 63 트릭으로 n=0일 때 64비트 시프트라는
 * 정의되지 않은 동작(UB)을 회피한다 */
static inline uint64_t rol64(uint64_t v, unsigned int n) {
  n &= 63;
  return (v << n) | (v >> ((64 - n) & 63));
}

/* 64비트 순환 우회전 (역상보 롤링 해시 갱신용) */
static inline uint64_t ror64(uint64_t v, unsigned int n) {
  n &= 63;
  return (v >> n) | (v << ((64 - n) & 63));
}

/* ntHash 염기별 64비트 시드 테이블. 각 염기는 이 값을 k-mer 내 위치에 따라
 * 순환회전시켜 XOR 하는 방식으로 해시에 반영된다 */
static const uint64_t NTHASH_H[4] = {
    0x3c8bf4f53c8bf4f5ULL, // A
    0x04c903a704c903a7ULL, // C
    0x2b8104c92b8104c9ULL, // G
    0x2e0600d3fd09e083ULL  // T
};

/* 선택된 [0, threshold) 해시를 32비트 전체 범위에 단조 확장한다.
 * 동일성 및 정렬 순서는 유지하면서 파티션 계산을 곱셈+shift로 단순화한다. */
static inline uint32_t normalize_sampled_hash(uint32_t hash,
                                              uint32_t threshold) {
  return (uint32_t)(((uint64_t)hash << 32) / threshold);
}

/* 윈도우 하나의 스케치(해시 집합) 추출.
 * ntHash 롤링 해시로 각 k-mer의 정방향/역상보 해시를 염기당 O(1)에 갱신하고,
 * 정규(canonical) 해시가 threshold 미만인 k-mer만 샘플링한다.
 * (threshold = UINT32_MAX/scale이므로 전체의 1/scale. minimizer처럼 서열 내용과
 *  무관하게 해시값 기준으로 선택하므로 두 서열에서 같은 k-mer가 선택됨이 보장됨)
 * 채택된 해시는 정렬+중복 제거 후 out_hashes에 저장한다. */
static inline void extract_hash_direct(const Segtrace *r, uint32_t *out_hashes,
                                       size_t *out_size, uint32_t threshold,
                                       const uint8_t *seq, size_t len) {
  uint32_t k = r->hash_window;
  if (len < k) {
    *out_size = 0;
    return;
  }

  uint64_t f_hash = 0, r_hash = 0;
  size_t valid_len = 0; /* 현재 위치까지 연속으로 유효한 염기 수 */
  size_t count = 0;

  for (size_t i = 0; i < len; i++) {
    int8_t b = r->base_lookup[seq[i]];
    if (b < 0) {
      /* N 등 무효 염기: 롤링 상태를 리셋하고 k-mer 누적을 다시 시작 */
      valid_len = 0;
      f_hash = 0;
      r_hash = 0;
      continue;
    }

    if (valid_len < k) {
      /* 첫 k-mer 구축 단계: 들어온 염기를 위치에 맞게 회전시켜 XOR.
       * b ^ 3은 상보 염기(A<->T, C<->G) */
      int8_t b_rc = b ^ 3;
      f_hash ^= rol64(NTHASH_H[b], k - 1 - (uint32_t)valid_len);
      r_hash ^= rol64(NTHASH_H[b_rc], (uint32_t)valid_len);
      valid_len++;
    } else {
      /* 롤링 갱신: 윈도우에서 빠지는 염기(i-k)의 기여를 제거하고 새로
       * 들어오는 염기(i)의 기여를 추가. 전체를 다시 해싱하지 않아 O(1) */
      int8_t b_out = r->base_lookup[seq[i - k]];
      f_hash = rol64(f_hash, 1) ^ rol64(NTHASH_H[b_out], k) ^ NTHASH_H[b];
      r_hash = ror64(r_hash, 1) ^ ror64(NTHASH_H[b_out ^ 3], 1) ^
               rol64(NTHASH_H[b ^ 3], k - 1);
    }

    if (valid_len >= k) {
      /* 정방향/역상보 중 작은 값을 정규 해시로 사용해 가닥 방향과 무관하게
       * 동일 k-mer에 동일 해시를 부여. mix_hash로 엔트로피를 32비트 전체에
       * 퍼뜨린 뒤 threshold 미만이면 스케치에 채택 */
      uint64_t canonical = (f_hash < r_hash) ? f_hash : r_hash;
      uint32_t h = mix_hash(canonical, r->hash_seed);
      if (h < threshold && count < MAX_SKETCH_SIZE) {
        out_hashes[count++] = normalize_sampled_hash(h, threshold);
      }
    }
  }

  /* 정렬 + 중복 제거: 이후 두 스케치의 교집합을 merge-join으로 세기 위한 전처리 */
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

/* 스레드 풀 작업: 염색체의 한 chunk 구간을 담당해 윈도우를 step_size씩
 * 이동하며 스케치를 만든다. 스레드 간 충돌을 피하기 위해 결과는 job 내부의
 * 로컬 배열에 쌓고, 종료 후 메인 스레드가 전역 배열로 병합한다. */
static void seq_chunk_worker(void *data, long i, int tid) {
  (void)tid;
  SeqChunkJob *job = &((SeqChunkJob *)data)[i];
  /* chunk 시작 오프셋으로부터 전역 윈도우 인덱스 복원
   * (chunk 크기가 step_size의 배수이므로 인덱스가 어긋나지 않음) */
  uint32_t current_window_idx =
      (uint32_t)(job->chunk_start_idx / job->step_size);

  uint32_t local_hashes[MAX_SKETCH_SIZE];

  size_t idx = job->chunk_start_idx;
  size_t valid_bases = 0;
  /* 첫 윈도우의 유효 염기 수만 직접 계산; 이후에는 아래에서 슬라이딩 갱신 */
  for (size_t j = 0; j < job->window_size; j++) {
    if (job->r->base_lookup[job->seq_ptr[idx + j]] >= 0)
      valid_bases++;
  }

  for (; idx + job->window_size <= job->chunk_end_idx;
       idx += job->step_size, current_window_idx++) {

    /* 유효 염기가 부족한(N-rich) 윈도우는 스케치하지 않고 빈 상태로 기록 */
    size_t sketch_size = 0;
    if (valid_bases >= job->min_bases) {
      extract_hash_direct(job->r, local_hashes, &sketch_size, job->threshold,
                          job->seq_ptr + idx, job->window_size);
    }

    /* 윈도우 메타데이터(서열 id, 윈도우 인덱스, 스케치 위치/크기) 기록 */
    DA_RESERVE(job->coords, job->cap_coords, job->num_coords + 1);
    WindowCoord *wc = &job->coords[job->num_coords++];
    wc->seq_id = job->seq_id;
    wc->window_idx = current_window_idx;
    wc->sketch_size = (uint16_t)sketch_size;
    /* read 모드(step_size >= window_size)에서는 윈도우 길이 = read 길이 */
    wc->read_len = (uint32_t)job->window_size;
    wc->sample_idx = 0;

    size_t h_idx = job->num_hashes;
    if (sketch_size > 0) {
      DA_RESERVE(job->hashes, job->cap_hashes, job->num_hashes + sketch_size);
      memcpy(job->hashes + h_idx, local_hashes, sketch_size * sizeof(uint32_t));
      job->num_hashes += sketch_size;
    }
    wc->sketch_offset = (uint64_t)h_idx;

    /* 다음 윈도우를 위해 유효 염기 수를 슬라이딩 갱신:
     * 앞에서 빠지는 step_size개를 빼고 뒤에서 들어오는 step_size개를 더함 */
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

/* 모든 입력 FASTA를 순회하며 전역 윈도우 테이블(gw.coords)과 전역 스케치
 * 해시 테이블(gw.all_hashes)을 구축한다.
 * threshold: 해시가 UINT32_MAX/scale 미만인 k-mer만 채택 (1/scale 샘플링) */
GlobalWindows extract_all_windows(char **files, int num_files,
                                  const Segtrace *r, uint64_t scale,
                                  size_t window_size, size_t step_size,
                                  size_t min_bases, int n_threads,
                                  void *thread_pool) {
  fprintf(stderr, "[segtrace] Extracting windows across genomes...\n");
  GlobalWindows gw = {0};
  size_t num_all_hashes = 0, cap_all_hashes = 0;
  size_t cap_sketches = 0, cap_seqs = 0;
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
      /* read 모드에서는 어떤 길이의 레코드도 하나의 스케치로 취급한다 */
      if (window_size && len < window_size)
        continue;

      DA_RESERVE(gw.seq_lens, cap_seqs, gw.num_seqs + 1);
      gw.seq_lens[gw.num_seqs].genome = strdup(bname);
      gw.seq_lens[gw.num_seqs].seq = strdup(ks->name.s);
      gw.seq_lens[gw.num_seqs].file_id = (uint32_t)f;
      uint32_t seq_id = (uint32_t)gw.num_seqs++;

      uint8_t *seq_ptr = (uint8_t *)ks->seq.s;

      /* window_size == 0 (read 모드): 서열 전체가 윈도우 하나 */
      size_t eff_window = window_size ? window_size : len;
      size_t eff_step = window_size ? step_size : len;
      size_t seq_windows =
          window_size ? (len - window_size) / step_size + 1 : 1;
      DA_RESERVE(gw.coords, cap_sketches, gw.num_sketches + seq_windows);
      DA_RESERVE(gw.all_hashes, cap_all_hashes,
             num_all_hashes + seq_windows * 96);

      /* 병렬화를 위해 염색체를 chunk로 분할.
       * 스레드 수의 4배로 쪼개 부하 균형을 맞추고 최소 100kb를 보장하며,
       * step_size의 배수로 맞춰 chunk 경계에서도 윈도우 인덱스가 어긋나지 않게 함 */
      size_t chunk_size = len / ((size_t)n_threads * 4);
      if (chunk_size < 100000)
        chunk_size = 100000;
      chunk_size = ((chunk_size + eff_step - 1) / eff_step) * eff_step;

      size_t cap_jobs = 16, num_jobs = 0;
      SeqChunkJob *jobs = malloc(cap_jobs * sizeof(SeqChunkJob));

      /* read 모드(step >= window)에서는 read 전체가 하나의 윈도우이므로
       * 마지막 위치 len - window도 포함해야 read가 누락되지 않는다 */
      for (size_t c_start = 0;; c_start += chunk_size) {
        /* chunk 끝을 window_size - step_size만큼 연장해 경계에 걸친 윈도우가
         * 누락되지 않도록 한다 (인접 chunk와 겹침) */
        size_t c_end = c_start + chunk_size + eff_window - eff_step;
        if (c_end > len)
          c_end = len;

        DA_RESERVE(jobs, cap_jobs, num_jobs + 1);
        jobs[num_jobs++] = (SeqChunkJob){.r = r,
                                         .threshold = threshold,
                                         .scale = scale,
                                         .window_size = eff_window,
                                         .step_size = eff_step,
                                         .min_bases = min_bases,
                                         .seq_id = seq_id,
                                         .seq_ptr = seq_ptr,
                                         .chunk_start_idx = c_start,
                                         .chunk_end_idx = c_end};
        if (c_end >= len)
          break;
      }

      /* chunk들을 스레드 풀로 병렬 스케치 */
      kt_forpool(thread_pool, seq_chunk_worker, jobs, (long)num_jobs);

      /* 각 chunk의 로컬 결과를 전역 배열로 병합.
       * sketch_offset은 로컬 기준이므로 전역 해시 배열 기준으로 보정한다 */
      for (size_t j = 0; j < num_jobs; j++) {
        SeqChunkJob *job = &jobs[j];
        if (job->num_coords > 0) {
          DA_RESERVE(gw.all_hashes, cap_all_hashes,
               num_all_hashes + job->num_hashes);
          DA_RESERVE(gw.coords, cap_sketches,
               gw.num_sketches + job->num_coords);
          size_t base_h_offset = num_all_hashes;
          memcpy(gw.all_hashes + base_h_offset, job->hashes,
                 job->num_hashes * sizeof(uint32_t));
          num_all_hashes += job->num_hashes;

          size_t base_c_offset = gw.num_sketches;
          for (size_t k = 0; k < job->num_coords; k++) {
            WindowCoord wc = job->coords[k];
            wc.sketch_offset += (uint64_t)base_h_offset;
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

  /* read 모드: 동일 스케치를 가진 read를 READ_SAMPLE_SLOTS개 샘플 슬롯으로
   * 나눠 그룹 크기 해상도를 높인다. 첫 샘플 해시의 상위 4비트로 슬롯을
   * 정하면 같은 스케치는 항상 같은 슬롯에 들어가고, 슬롯당 read 수의
   * 기대치는 전체 깊이의 1/READ_SAMPLE_SLOTS가 된다. */
  if (g_segtrace_read_mode) {
    for (size_t i = 0; i < gw.num_sketches; i++) {
      WindowCoord *wc = &gw.coords[i];
      if (wc->sketch_size == 0)
        continue;
      uint32_t first = gw.all_hashes[wc->sketch_offset];
      wc->sample_idx = first >> READ_SLOT_SHIFT;
    }
  }
  return gw;
}

// ==============================================================
// SECTION 3: CANDIDATE DISCOVERY & DISTANCE COMPUTATION
// ==============================================================

/* 두 정렬된 스케치의 교집합 크기(공유 해시 개수)를 merge-join으로 계산.
 * 이름은 dist지만 실제로는 유사도(클수록 유사)를 반환한다 */
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

/* 윈도우 메타데이터에서 스케치 위치를 찾아 교집합 크기를 계산하는 래퍼 */
static inline size_t calculate_window_dist(const uint32_t *all_hashes,
                                           const WindowCoord *wa,
                                           const WindowCoord *wb) {
  return calculate_sketch_dist(all_hashes + wa->sketch_offset, wa->sketch_size,
                               all_hashes + wb->sketch_offset, wb->sketch_size);
}

/* 두 윈도우가 같은 서열 위에서 좌표가 겹치는지 검사.
 * 겹치는 윈도우는 서열이 본질적으로 같으므로 비교 대상에서 제외한다 */
static inline int windows_overlap(const DiscoverComputeData *w, uint32_t wa,
                                  uint32_t wb) {
  if (w->coords[wa].seq_id != w->coords[wb].seq_id)
    return 0;
  size_t window_distance =
      (size_t)ABS_DIFF(w->coords[wa].window_idx, w->coords[wb].window_idx);
  return window_distance * w->step_size < w->window_size;
}

/* 두 윈도우 중 더 큰 스케치 크기 반환 (유사도 정규화의 분모) */
static inline size_t max_sketch_size(const DiscoverComputeData *w, uint32_t wa,
                                     uint32_t wb) {
  return w->coords[wa].sketch_size > w->coords[wb].sketch_size
             ? w->coords[wa].sketch_size
             : w->coords[wb].sketch_size;
}

/* 매치로 인정하기 위한 최소 공유 해시 수.
 * p_kmer = identity^k = "k-mer 하나가 두 서열 간에 보존될 확률"이므로
 * 기대 공유 스케치 수는 대략 (스케치 크기) x p_kmer.
 * 하한 3개를 둬 스케치가 작은 윈도우의 노이즈 매치를 걸러낸다 */
static inline size_t required_shared(const DiscoverComputeData *w, uint32_t wa,
                                     uint32_t wb) {
  size_t max_size = max_sketch_size(w, wa, wb);
  size_t min_shared = (size_t)ceil((double)max_size * w->p_kmer);
  return min_shared < 3 ? 3 : min_shared;
}

/* collinear 탐색용 헬퍼: (wa, wb)가 배열 범위 안이고, 기대하는 서열 쌍
 * (seq_a, seq_b)과 일치하며, 겹치지 않고, 유사도 기준을 통과하는지 검사 */
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

/* 후보 쌍 주변에 '연쇄적인(collinear)' 유사 윈도우가 더 있는지 검사.
 * 진짜 segmental duplication은 여러 윈도우에 걸쳐 대각선 상에 연속으로
 * 나타나므로, 고립된 단일 윈도우 매치(우연한 반복 서열 등)를 걸러내는 역할.
 * dir_a/dir_b 조합으로 정방향/역방향, 양쪽 진행 방향의 대각선 4가지를 검사 */
static inline int check_collinear_neighbor(const DiscoverComputeData *w,
                                           uint32_t wa, uint32_t wb) {
  uint32_t seq_a = w->coords[wa].seq_id;
  uint32_t seq_b = w->coords[wb].seq_id;

  const int dir_a[] = {1, -1, 1, -1};
  const int dir_b[] = {1, -1, -1, 1};

  // Pass 1: 정확한 대각선 (indel 없이 양쪽이 같은 step으로 진행하는 경우)
  for (int d = 0; d < 4; d++) {
    int da = dir_a[d], db = dir_b[d];
    for (int step = 1; step <= MAX_COLLINEAR_LOOKAHEAD; step++) {
      long long next_a = (long long)wa + da * step;
      long long next_b = (long long)wb + db * step;
      if (matching_window_pair(w, next_a, next_b, seq_a, seq_b))
        return 1;
    }
  }

  // Pass 2: indel이 낀 대각선 (양쪽의 진행 step이 다른 경우까지 허용)
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

typedef struct {
  const uint32_t *all_hashes;
  const WindowCoord *coords;
  uint16_t *win_curr_pos;
  size_t n_windows;
  size_t n_blocks;
  size_t batch_start;
  size_t batch_count;
  size_t *block_offsets;
  HashWindowEntry *entries;
} BucketBuildData;

static inline size_t hash_partition(uint32_t hash) {
  return (size_t)(((uint64_t)hash * NUM_PARTITIONS) >> 32);
}

/* 각 block이 담당하는 윈도우에서 이번 배치의 파티션별 엔트리 수를 센다. */
static void count_bucket_entries(void *data, long idx, int tid) {
  (void)tid;
  BucketBuildData *build = (BucketBuildData *)data;
  size_t block = (size_t)idx;
  size_t win_start = build->n_windows * block / build->n_blocks;
  size_t win_end = build->n_windows * (block + 1) / build->n_blocks;
  size_t *counts = build->block_offsets + block * build->batch_count;

  for (size_t win = win_start; win < win_end; win++) {
    uint64_t off = build->coords[win].sketch_offset;
    uint16_t size = build->coords[win].sketch_size;
    uint16_t pos = build->win_curr_pos[win];
    while (pos < size) {
      uint32_t hash = build->all_hashes[off + pos];
      size_t partition = hash_partition(hash);
      if (partition >= build->batch_start + build->batch_count)
        break;
      counts[partition - build->batch_start]++;
      pos++;
    }
  }
}

/* prefix sum으로 미리 배정된 block별 구간에 엔트리를 락 없이 기록한다. */
static void scatter_bucket_entries(void *data, long idx, int tid) {
  (void)tid;
  BucketBuildData *build = (BucketBuildData *)data;
  size_t block = (size_t)idx;
  size_t win_start = build->n_windows * block / build->n_blocks;
  size_t win_end = build->n_windows * (block + 1) / build->n_blocks;
  size_t *offsets = build->block_offsets + block * build->batch_count;

  for (size_t win = win_start; win < win_end; win++) {
    uint64_t off = build->coords[win].sketch_offset;
    uint16_t size = build->coords[win].sketch_size;
    uint16_t pos = build->win_curr_pos[win];
    while (pos < size) {
      uint32_t hash = build->all_hashes[off + pos];
      size_t partition = hash_partition(hash);
      if (partition >= build->batch_start + build->batch_count)
        break;
      size_t local_partition = partition - build->batch_start;
      build->entries[offsets[local_partition]++] =
          (HashWindowEntry){hash, (uint32_t)win};
      pos++;
    }
    build->win_curr_pos[win] = pos;
  }
}

/* 파티션(해시 구간) 하나를 담당하는 스레드 작업.
 * 버킷을 (hash, window_id)로 정렬하면 같은 해시를 가진 윈도우들이 연속된
 * run을 이루고, run 내부의 윈도우 쌍이 곧 "해시를 공유하는 후보 쌍"이다. */
static void discover_compute_worker(void *data, long idx, int tid) {
  DiscoverComputeData *w_data = (DiscoverComputeData *)data;
  long p = (long)w_data->batch_start + idx;
  PartitionBucket *b = &w_data->buckets[p];
  if (b->size == 0)
    return;

  qsort(b->entries, b->size, sizeof(HashWindowEntry), compare_hash_entry);

  /* 같은 해시 값을 가진 run을 하나씩 순회 */
  size_t i = 0;
  while (i < b->size) {
    size_t j = i + 1;
    while (j < b->size && b->entries[j].hash == b->entries[i].hash)
      j++;
    size_t run_len = j - i;

    /* run이 너무 크면(초고빈도 k-mer, 저복잡도/반복 서열) 비교 폭발 방지를
     * 위해 건너뛴다. 2 이상이어야 공유 쌍이 존재 */
    if (run_len >= 2 && run_len <= MAX_KMER_FREQ) {
      /* 각 윈도우당 인접한 MAX_PAIR_COMPARISONS개의 이웃만 비교해
       * run 내 비교 횟수를 선형으로 제한 (버스트 방지) */
      for (size_t a = i; a < j; a++) {
        size_t b_max =
            a + 1 + MAX_PAIR_COMPARISONS < j ? a + 1 + MAX_PAIR_COMPARISONS : j;
        for (size_t b_idx = a + 1; b_idx < b_max; b_idx++) {
          uint32_t wa = b->entries[a].window_id,
                   wb = b->entries[b_idx].window_id;
          /* 같은 서열의 겹치는 윈도우끼리는 비교하지 않음 */
          if (windows_overlap(w_data, wa, wb))
            continue;

          /* 같은 쌍이 다른 해시 run이나 스레드에서 반복 발견될 수 있으므로
           * 공유 cache-blocked 블룸필터로 중복 제거 */
          uint64_t pk = encode_pair(wa, wb);
          if (bloom_test_and_set(w_data->bloom, pk))
            continue;

          size_t min_shared = required_shared(w_data, wa, wb);

          /* 전체 스케치 교집합 크기로 유사도 산정. 기준 미달이거나
           * collinear 이웃이 없으면(고립 매치면) 탈락 */
          size_t shared = calculate_window_dist(
              w_data->all_hashes, &w_data->coords[wa], &w_data->coords[wb]);
          if (shared < min_shared ||
              !check_collinear_neighbor(w_data, wa, wb))
            continue;

              /* score = 공유 비율을 0~255로 정규화 (반올림 포함). */
          size_t max_size = max_sketch_size(w_data, wa, wb);
          uint32_t score =
              (uint32_t)((shared * UINT8_MAX + max_size / 2) / max_size);
          DA_PUSH(w_data->t_pairs[tid], w_data->t_n_pairs[tid],
                  w_data->t_cap_pairs[tid],
                  ((CandidatePair){wa, wb, (uint8_t)score}));
        }
      }
    }
    i = j;
  }
}

/* 후보 쌍 탐색의 메인 드라이버.
 * 32비트 해시 공간을 NUM_PARTITIONS개로 나누고 BATCH_PARTITIONS개씩 묶어
 * 배치 처리해 메모리 사용량을 제한한다 (외부 정렬과 유사한 발상).
 * 같은 해시는 항상 같은 파티션에 들어가므로 분할이 결과에 영향을 주지 않는다. */
CandidateGraph discover_and_compute(const uint32_t *all_hashes,
                                    const WindowCoord *coords,
                                    size_t n_windows, size_t window_size,
                                    size_t step_size, int n_threads,
                                    uint32_t kmer_size,
                                    void *thread_pool) {
  if (n_windows > (size_t)CANDIDATE_WINDOW_MASK + 1) {
    /* 쌍 인코딩이 윈도우 id에 32비트만 쓰므로 상한 검사 */
    fprintf(stderr, "[ERROR] Too many windows for candidate encoding\n");
    exit(1);
  }

  DiscoverComputeData w = {
      .all_hashes = all_hashes,
      .coords = coords,
      .n_windows = n_windows,
      .window_size = window_size,
      .step_size = step_size,
      /* p_kmer = identity^k: k-mer 하나가 두 서열 간에 보존될 확률.
       * 공유 스케치 수의 기대치는 대략 (스케치 크기) x p_kmer */
      .p_kmer = pow(MIN_IDENTITY, (double)kmer_size),
      .buckets = calloc(NUM_PARTITIONS, sizeof(PartitionBucket)),
      .bloom = calloc(BLOOM_NUM_WORDS, sizeof(uint64_t)),
      .t_pairs = calloc(n_threads, sizeof(CandidatePair *)),
      .t_n_pairs = calloc(n_threads, sizeof(size_t)),
      .t_cap_pairs = calloc(n_threads, sizeof(size_t))};

  if (!w.buckets || !w.bloom || !w.t_pairs || !w.t_n_pairs ||
      !w.t_cap_pairs) {
    fprintf(stderr, "[ERROR] Memory allocation failed\n");
    exit(1);
  }

  /* 윈도우별로 "이번 배치까지 몇 번째 스케치 해시를 흘려보냈는지" 기록.
   * 스케치가 정렬되어 있으므로 배치마다 처음부터 다시 훑지 않아도 된다 */
  uint16_t *win_curr_pos = calloc(n_windows, sizeof(uint16_t));
  if (n_windows && !win_curr_pos) {
    fprintf(stderr, "[ERROR] Memory allocation failed\n");
    exit(1);
  }
  size_t n_blocks = (size_t)n_threads * 8;
  if (n_blocks > n_windows)
    n_blocks = n_windows;
  size_t *block_offsets =
      n_blocks ? calloc(n_blocks * BATCH_PARTITIONS, sizeof(size_t)) : NULL;
  if (n_blocks && !block_offsets) {
    fprintf(stderr, "[ERROR] Memory allocation failed\n");
    exit(1);
  }

  for (size_t batch_start = 0; batch_start < NUM_PARTITIONS;
       batch_start += BATCH_PARTITIONS) {
    size_t batch_end = batch_start + BATCH_PARTITIONS;
    if (batch_end > NUM_PARTITIONS)
      batch_end = NUM_PARTITIONS;

    w.batch_start = batch_start;
    size_t batch_count = batch_end - batch_start;
    /* block별 count와 prefix sum으로 정확한 크기의 연속 버킷을 만든 뒤,
     * 각 block의 전용 구간에 병렬로 엔트리를 기록한다. */
    if (n_blocks)
      memset(block_offsets, 0,
             n_blocks * batch_count * sizeof(*block_offsets));
    BucketBuildData build = {.all_hashes = all_hashes,
                             .coords = coords,
                             .win_curr_pos = win_curr_pos,
                             .n_windows = n_windows,
                             .n_blocks = n_blocks,
                             .batch_start = batch_start,
                             .batch_count = batch_count,
                             .block_offsets = block_offsets};
    if (n_blocks)
      kt_forpool(thread_pool, count_bucket_entries, &build, (long)n_blocks);

    size_t total_entries = 0;
    for (size_t local = 0; local < batch_count; local++) {
      PartitionBucket *bucket = &w.buckets[batch_start + local];
      bucket->size = 0;
      for (size_t block = 0; block < n_blocks; block++)
        bucket->size += block_offsets[block * batch_count + local];
      total_entries += bucket->size;
    }

    HashWindowEntry *batch_entries =
        total_entries ? malloc(total_entries * sizeof(*batch_entries)) : NULL;
    if (total_entries && !batch_entries) {
      fprintf(stderr, "[ERROR] Memory allocation failed\n");
      exit(1);
    }
    size_t partition_offset = 0;
    for (size_t local = 0; local < batch_count; local++) {
      PartitionBucket *bucket = &w.buckets[batch_start + local];
      bucket->entries = bucket->size ? batch_entries + partition_offset : NULL;
      size_t block_offset = partition_offset;
      for (size_t block = 0; block < n_blocks; block++) {
        size_t offset_idx = block * batch_count + local;
        size_t count = block_offsets[offset_idx];
        block_offsets[offset_idx] = block_offset;
        block_offset += count;
      }
      partition_offset += bucket->size;
    }
    build.entries = batch_entries;
    if (n_blocks)
      kt_forpool(thread_pool, scatter_bucket_entries, &build, (long)n_blocks);

    /* 배치 내 파티션들을 병렬 처리 */
    kt_forpool(thread_pool, discover_compute_worker, &w, (long)batch_count);

    /* 배치가 끝난 연속 버킷 메모리는 즉시 해제 */
    free(batch_entries);
    for (size_t p = batch_start; p < batch_end; p++) {
      w.buckets[p] = (PartitionBucket){0};
    }
  }

  free(w.buckets);
  free(w.bloom);
  free(w.t_cap_pairs);
  free(win_curr_pos);
  free(block_offsets);

  return (CandidateGraph){
      .pairs = w.t_pairs, .counts = w.t_n_pairs, .n_threads = n_threads};
}

// ==============================================================
// SECTION 3b: READ-MODE EXACT GROUPING & COVERAGE-BASED COPY NUMBER
// ==============================================================

/* 스케치(정렬된 해시 배열)의 FNV-1a 지문.
 * 스케치는 정규화된 32비트 해시의 정렬 배열이므로 같은 스케치는 항상 같은
 * 지문을 가지고, 다른 스케치가 같은 지문을 가질 확률은 ~2^-64이다.
 * 지문 비교만으로 그룹화하므로 해시 항목 배열을 추가로 만들지 않아
 * 메모리가 O(윈도우 수)에 머문다. */
static uint64_t sketch_fingerprint(const uint32_t *h, size_t n) {
  uint64_t fp = UINT64_C(1469598103934665603);
  for (size_t i = 0; i < n; i++) {
    uint32_t v = h[i];
    for (int b = 0; b < 4; b++) {
      fp ^= (v >> (b * 8)) & 0xff;
      fp *= UINT64_C(1099511628211);
    }
  }
  return fp;
}

static int compare_read_group(const void *a, const void *b) {
  const ReadGroupStat *ga = (const ReadGroupStat *)a,
                      *gb = (const ReadGroupStat *)b;
  if (ga->sketch_size != gb->sketch_size)
    return CMP(ga->sketch_size, gb->sketch_size);
  if (ga->sample_idx != gb->sample_idx)
    return CMP(ga->sample_idx, gb->sample_idx);
  return CMP(ga->fingerprint, gb->fingerprint);
}

/* 슬롯을 무시하고 스케치 자체만 비교 (슬롯 합산용 2단계 정렬) */
static int compare_read_group_noslot(const void *a, const void *b) {
  const ReadGroupStat *ga = (const ReadGroupStat *)a,
                      *gb = (const ReadGroupStat *)b;
  if (ga->sketch_size != gb->sketch_size)
    return CMP(ga->sketch_size, gb->sketch_size);
  return CMP(ga->fingerprint, gb->fingerprint);
}

/* 스케치가 완전히 같은 read를 그룹화해 그룹 크기(관측 개수)를 구한다.
 * (sketch_size, fingerprint)로 정렬하면 동일 스케치가 연속 run을 이루고,
 * run 크기가 곧 그 세그먼트를 커버하는 read 수다. 빈 스케치는 제외한다.
 * k-mer 스페이스가 충분히 크므로 (4^k >> 윈도우 수) 서열이 다른 두 read가
 * 완전히 같은 스케치를 가질 확률은 무시할 수 있다. */
size_t group_identical_reads(const uint32_t *all_hashes,
                             const WindowCoord *coords, size_t n_windows,
                             int n_threads, void *thread_pool,
                             ReadGroupStat **out_groups) {
  fprintf(stderr, "[segtrace] Grouping reads by identical sketches...\n");
  (void)n_threads;
  (void)thread_pool;

  ReadGroupStat *groups = NULL;
  size_t n_groups = 0, cap_groups = 0;
  for (size_t i = 0; i < n_windows; i++) {
    if (coords[i].sketch_size == 0)
      continue;
    DA_PUSH(groups, n_groups, cap_groups,
            ((ReadGroupStat){(uint32_t)i, 1,
                             sketch_fingerprint(
                                 all_hashes + coords[i].sketch_offset,
                                 coords[i].sketch_size),
                             coords[i].read_len, coords[i].sketch_size,
                             coords[i].sample_idx}));
  }
  if (n_groups == 0) {
    *out_groups = groups;
    return 0;
  }

  qsort(groups, n_groups, sizeof(ReadGroupStat), compare_read_group);

  /* 1단계: 동일 (size, sample_idx, fingerprint) run을 하나의 슬롯 그룹으로
   * 압축. 슬롯 그룹의 count가 곧 그 슬롯에서 관측된 read 수다 */
  size_t out = 0;
  size_t i = 0;
  while (i < n_groups) {
    size_t j = i + 1;
    while (j < n_groups && groups[j].sketch_size == groups[i].sketch_size &&
           groups[j].sample_idx == groups[i].sample_idx &&
           groups[j].fingerprint == groups[i].fingerprint)
      j++;
    groups[out] = groups[i];
    groups[out].count = (uint32_t)(j - i);
    out++;
    i = j;
  }
  n_groups = out;

  /* 2단계: 슬롯 키를 제외하고 (size, fingerprint)만으로 다시 정렬해
   * 같은 스케치의 16개 슬롯을 하나의 세그먼트로 합친다 */
  qsort(groups, n_groups, sizeof(ReadGroupStat), compare_read_group_noslot);
  out = 0;
  i = 0;
  while (i < n_groups) {
    size_t j = i + 1;
    uint64_t total = groups[i].count;
    while (j < n_groups && groups[j].sketch_size == groups[i].sketch_size &&
           groups[j].fingerprint == groups[i].fingerprint) {
      total += groups[j].count;
      j++;
    }
    groups[out] = groups[i];
    /* 슬롯 합산이 uint32를 넘는 극단적인 심층 데이터에서는 포화시킨다 */
    groups[out].count =
        total > UINT32_MAX ? UINT32_MAX : (uint32_t)total;
    out++;
    i = j;
  }
  *out_groups = groups;
  return out;
}

/* 그룹 크기(관측 read 수) 히스토그램의 첫 번째(최소 크기) 봉우리를 찾아
 * 배수체(haploid) 커버리지로 삼는다.
 * 수학적 근거: 각 read는 첫 샘플 해시의 상위 4비트가 가리키는 정확히
 * 하나의 슬롯에 속하므로, 슬롯들은 한 세그먼트의 read를 서로 겹치지 않게
 * 분할한다. 따라서 슬롯당 관측수는 전체 read 수의 불편 추정량이고,
 * 슬롯당 관측수 히스토그램의 최빈값이 곧 1-copy(haploid) 세그먼트의
 * 총 read 수 C1이다 (genomescope가 k-mer 빈도 히스토그램의 첫 봉우리를
 * 쓰는 것과 같은 원리). 유전체에서 가장 흔한 세그먼트가 1-copy라는
 * 가정이 성립할 때만 유효하다. 추정이 불가능하면(그룹 없음) 0을
 * 돌려준다. */
double estimate_haploid_coverage(const ReadGroupStat *groups, size_t n_groups,
                                 uint64_t scale, size_t window_size) {
  /* scale과 window_size는 현재 추정기가 쓰지 않지만, 향후 스케치 크기 기반
   * 추정으로 확장할 때 시그니처를 유지하기 위해 남겨둔다 */
  (void)scale;
  (void)window_size;
  if (n_groups == 0)
    return 0.0;

  uint32_t max_count = 0;
  for (size_t i = 0; i < n_groups; i++)
    if (groups[i].count > max_count)
      max_count = groups[i].count;

  size_t *hist = calloc((size_t)max_count + 1, sizeof(size_t));
  if (!hist) {
    fprintf(stderr, "[ERROR] Memory allocation failed\n");
    exit(1);
  }
  for (size_t i = 0; i < n_groups; i++)
    hist[groups[i].count]++;

  /* 최빈 count. count=1(한 번만 관측된 read)은 대부분 시퀀싱 오류
   * k-mer/리드이므로 genomescope 관례에 따라 후보에서 제외하고
   * count>=2부터 최빈값을 찾는다. 모든 그룹이 count=1이면 mode=1이다 */
  size_t mode = 1;
  for (size_t i = 2; i <= max_count; i++)
    if (hist[i] > hist[mode])
      mode = i;
  free(hist);

  /* 슬롯당 count가 곧 전체 read 수의 추정치다 (각 read가 정확히 하나의
   * 슬롯에만 속하므로 슬롯당 count 평균은 전체 count의 불편 추정량).
   * 별도 환산 없이 mode 자체가 haploid 세그먼트의 read 수다. */
  return (double)mode;
}

/* read 모드 출력: 각 그룹(고유 세그먼트)을 대표 read 이름으로 BED에 쓴다.
 * 그룹은 (스케치, 슬롯) 단위로 나뉘어 있으므로, 같은 스케치의 슬롯들을
 * 합쳐 전체 read 수를 복원한 뒤 복제수 = round(총 read 수 / C1)를 계산한다.
 * C1은 estimate_haploid_coverage가 준 슬롯 환산 haploid 깊이.
 * copy < min_copies인 그룹은 제외한다 (기본 -c 1이면 모두 출력). */
void write_read_seg_bed(const char *out_prefix, const WindowCoord *coords,
                        const GenomeSeqLen *seq_lens,
                        const ReadGroupStat *groups, size_t n_groups,
                        uint64_t scale, size_t window_size,
                        double haploid_coverage, uint32_t min_copies) {
  if (n_groups == 0)
    return;
  if (min_copies < 1)
    min_copies = 1;

  char path_buf[PATH_MAX];
  snprintf(path_buf, sizeof(path_buf), "%s.seg.bed", out_prefix);
  FILE *out = fopen(path_buf, "w");
  if (!out) {
    fprintf(stderr, "[ERROR] Cannot open output file: %s\n", path_buf);
    return;
  }

  /* scale/window_size는 현재 출력에 쓰이지 않지만, 추정기와 시그니처를
   * 맞추기 위해 유지한다 */
  (void)scale;
  (void)window_size;
  fprintf(out, "#chrom\tstart\tend\tread_count\test_depth\test_copies\n");

  /* group_identical_reads가 이미 슬롯을 합산해 그룹당 총 read 수를
   * 넣어두었으므로, 여기서는 그룹당 한 행만 출력하면 된다 */
  for (size_t i = 0; i < n_groups; i++) {
    const ReadGroupStat *g = &groups[i];
    double depth = (double)g->count;
    uint32_t copies = haploid_coverage > 0.0
                          ? (uint32_t)(depth / haploid_coverage + 0.5)
                          : 1;
    if (copies < min_copies)
      continue;
    uint32_t seq_i = coords[g->window_id].seq_id;
    fprintf(out, "%s-%s\t0\t%u\t%u\t%.1f\t%u\n", seq_lens[seq_i].genome,
            seq_lens[seq_i].seq, (uint32_t)g->read_len, g->count, depth,
            copies);
  }
  fclose(out);
}




/* 인코딩된 값에서 윈도우 id(하위 28비트) 추출 */
static inline uint32_t candidate_window(uint32_t encoded) {
  return encoded;
}

static inline uint32_t candidate_score(CandidatePair pair) {
  return pair.score;
}

/* 후보 쌍에 한 번이라도 등장한 윈도우들을 같은 서열 내 인접 윈도우끼리
 * 하나의 연속 구간(locus)으로 병합한다.
 * coords[].sketch_offset/size 필드를 재활용해:
 *   sketch_size   = 후보 포함 여부 플래그
 *   sketch_offset = 해당 윈도우가 속한 region 인덱스 */
void build_duplicate_loci(const CandidateGraph *graph, size_t num_windows,
                          WindowCoord *coords, const GenomeSeqLen *seq_lens,
                          size_t step_size, size_t window_size,
                          SegtraceDupRegion **out_regions,
                          size_t *out_n_regions) {
  /* 1) 마킹 초기화: 모든 윈도우를 "후보 아님" 상태로 */
  for (size_t i = 0; i < num_windows; i++) {
    coords[i].sketch_offset = UINT64_MAX;
    coords[i].sketch_size = 0;
  }
  for (int t = 0; t < graph->n_threads; t++) {
    for (size_t i = 0; i < graph->counts[t]; i++) {
      CandidatePair pair = graph->pairs[t][i];
      coords[candidate_window(pair.a)].sketch_size = 1;
      coords[candidate_window(pair.b)].sketch_size = 1;
    }
  }

  /* 2) 마킹된 윈도우를 순서대로 스캔하며, 같은 서열에서 윈도우 인덱스 차이가
   * MAX_COLLINEAR_LOOKAHEAD 이내면 같은 구간으로 병합.
   * coords가 서열/인덱스 순으로 생성되므로 단순 선형 스캔으로 충분하다 */
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

    int merge = n_regions > 0 && seq_id == previous_seq &&
          window_idx - previous_window <= MAX_COLLINEAR_LOOKAHEAD;

    if (merge) {
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

/* 정렬 기준: 클러스터 id -> 파일 -> 서열 -> 좌표 순.
 * 이 정렬 덕분에 (클러스터, 파일) 그룹이 연속 배치되어 copy count 세기가 쉬워짐 */
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

/* 구간 간 클러스터링:
 * 1패스 - 각 구간에 대해 score가 가장 높은 파트너 1개만 기록
 *         (같은 locus 내부의 self-match는 별도 카피가 아닌 잡음이므로
 *         best-partner 경쟁에서 완전히 제외한다),
 * 2패스 - best-partner 간선들을 union-find로 연결해 클러스터 id 부여,
 * 마지막으로 클러스터/파일/좌표 순으로 정렬한다. */
void cluster_duplicate_loci(const CandidateGraph *graph,
                            const WindowCoord *coords,
                            SegtraceDupRegion *regions, size_t n_regions) {
  if (n_regions == 0)
    return;

  /* 1패스: best partner 선정 (cluster_id 필드를 임시로 best score 저장에 사용) */
  for (int t = 0; t < graph->n_threads; t++) {
    for (size_t i = 0; i < graph->counts[t]; i++) {
      CandidatePair pair = graph->pairs[t][i];
      uint64_t region_a = coords[candidate_window(pair.a)].sketch_offset;
      uint64_t region_b = coords[candidate_window(pair.b)].sketch_offset;
      if (region_a == UINT64_MAX || region_b == UINT64_MAX ||
          region_a == region_b)
        continue;

      uint32_t score = candidate_score(pair);
      if (score > regions[region_a].cluster_id ||
          (score == regions[region_a].cluster_id &&
           region_b < regions[region_a].partner_id)) {
        regions[region_a].cluster_id = score;
        regions[region_a].partner_id = (uint32_t)region_b;
      }
      if (score > regions[region_b].cluster_id ||
          (score == regions[region_b].cluster_id &&
           region_a < regions[region_b].partner_id)) {
        regions[region_b].cluster_id = score;
        regions[region_b].partner_id = (uint32_t)region_a;
      }
    }
  }

  /* 2패스: best-partner 관계로 union-find를 수행해 연결 요소마다 클러스터 부여 */
  UnionFind uf;
  init_unionfind(&uf, n_regions);
  for (size_t i = 0; i < n_regions; i++) {
    if (regions[i].partner_id < n_regions)
      union_unionfind(&uf, (uint32_t)i, regions[i].partner_id);
  }

  /* root -> 1부터 시작하는 연속된 클러스터 번호로 재매핑 */
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

/* 스레드별로 쌓인 후보 쌍 배열들을 모두 해제 */
void free_candidate_graph(CandidateGraph *graph) {
  for (int t = 0; t < graph->n_threads; t++)
    free(graph->pairs[t]);
  free(graph->pairs);
  free(graph->counts);
}

/* 클러스터에 한 개 이상의 파일에서 복제 수가 min_copies 이상인 구간이
 * 있으면, 각 (클러스터, 파일) 그룹 중 파일별 복제 수가 min_copies 이상인
 * 그룹만 남긴다. 한 유전체 안에만 존재하는 클러스터도 출력 대상이다.
 * regions가 cluster/file 순으로 정렬된 상태이므로 각 그룹은 연속 구간이고,
 * in-place로 압축한 뒤 살아남은 개수를 반환한다.
 * 탈락한 클러스터로 인해 번호가 듬성듬성해지지 않도록, 살아남은 클러스터만
 * 등장 순서대로 1부터 다시 번호를 매긴다. */
size_t filter_regions_by_copy_count(SegtraceDupRegion *regions, size_t n,
                                    uint32_t min_copies) {
  if (n == 0)
    return 0;
  if (min_copies < 1)
    min_copies = 1;

  size_t out_count = 0;
  uint32_t next_cluster_id = 1;
  size_t ci = 0;
  while (ci < n) {
    size_t cj = ci + 1;
    while (cj < n && regions[cj].cluster_id == regions[ci].cluster_id)
      cj++;

    size_t cluster_out_start = out_count;
    size_t i = ci;
    while (i < cj) {
      size_t j = i + 1;
      while (j < cj && regions[j].file_id == regions[i].file_id)
        j++;
      if (j - i >= min_copies) {
        for (size_t k = i; k < j; k++) {
          regions[out_count++] = regions[k];
        }
      }
      i = j;
    }

    if (out_count > cluster_out_start) {
      uint32_t out_cluster_id = next_cluster_id++;
      for (size_t k = cluster_out_start; k < out_count; k++)
        regions[k].cluster_id = out_cluster_id;
    }
    ci = cj;
  }
  return out_count;
}

/* 최종 구간을 "<prefix>.dup.bed"에 기록.
 * chrom 컬럼은 "파일명-서열명" 형태, 4번째 컬럼은 클러스터 id.
 * min_sd_len 미만의 짧은 구간은 출력하지 않는다 */
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

/* Union-Find (Disjoint Set) 자료구조: path halving + union-by-rank 사용 */
void init_unionfind(UnionFind *uf, size_t n) {
  uf->parent = malloc(n * sizeof(uint32_t));
  uf->rank = calloc(n, sizeof(uint8_t));
  for (size_t i = 0; i < n; i++)
    uf->parent[i] = (uint32_t)i;
}

/* 루트를 찾아가며 경로 상의 노드를 조부모에 연결 (path halving) */
uint32_t find_unionfind(UnionFind *uf, uint32_t x) {
  while (uf->parent[x] != x) {
    uint32_t next = uf->parent[x];
    uf->parent[x] = uf->parent[next];
    x = next;
  }
  return x;
}

/* rank가 낮은 트리를 높은 트리 아래에 붙이는 union-by-rank */
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

/* 경로와 알려진 확장자(.fa/.fasta/.gz 등)를 제거한 순수 파일명 추출.
 * 출력 BED의 genome 라벨로 사용된다 */
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

/* wyhash 스타일의 64->32비트 최종 믹서.
 * ntHash 값의 상위 비트 엔트로피를 32비트 전체에 퍼뜨려
 * threshold 기반 샘플링이 편향되지 않게 한다 */
inline uint32_t mix_hash(uint64_t hash_value, uint64_t seed) {
  hash_value ^= seed;
  hash_value ^= hash_value >> 33;
  hash_value *= MIX_CONST1;
  hash_value ^= hash_value >> 33;
  hash_value *= MIX_CONST2;
  hash_value ^= hash_value >> 33;
  return (uint32_t)hash_value;
}

/* 순서 무관하게 (a,b) 쌍을 하나의 64비트 키로 인코딩 (블룸필터 키용) */
uint64_t encode_pair(uint32_t a, uint32_t b) {
  return a < b ? ((uint64_t)a << 32) | b : ((uint64_t)b << 32) | a;
}

/* splitmix64 최종 해시: 인접한 키 값들이 블룸필터에 고르게 퍼지도록 함 */
static inline uint64_t splitmix64(uint64_t x) {
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  return x ^ (x >> 31);
}

/* 공유 cache-blocked 블룸필터 조회 겸 삽입. 한 키의 세 비트를 같은 64비트
 * word에 배치해 메모리 접근을 한 cache line으로 제한하고, atomic OR로
 * 여러 worker가 락 없이 공유한다. */
int bloom_test_and_set(uint64_t *bloom, uint64_t key) {
  uint64_t h = splitmix64(key);
  uint32_t word_idx = (uint32_t)h & (BLOOM_NUM_WORDS - 1);
  uint64_t bits = (UINT64_C(1) << ((h >> 22) & 63)) |
                  (UINT64_C(1) << ((h >> 36) & 63)) |
                  (UINT64_C(1) << ((h >> 50) & 63));
  uint64_t old =
      __atomic_fetch_or(&bloom[word_idx], bits, __ATOMIC_RELAXED);
  return (old & bits) == bits;
}

/* qsort 비교 함수: uint32 오름차순 */
int compare_uint32(const void *a, const void *b) {
  uint32_t va = *(const uint32_t *)a, vb = *(const uint32_t *)b;
  return (va > vb) - (va < vb);
}

/* qsort 비교 함수: 해시 오름차순, 같은 해시면 윈도우 id 오름차순 */
int compare_hash_entry(const void *a, const void *b) {
  const HashWindowEntry *ea = (const HashWindowEntry *)a,
                        *eb = (const HashWindowEntry *)b;
  return ea->hash != eb->hash ? CMP(ea->hash, eb->hash)
                              : CMP(ea->window_id, eb->window_id);
}
