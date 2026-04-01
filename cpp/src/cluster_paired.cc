/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.

*/

#include "cluster_paired.h"

#include "align_simd.h"
#include "attributes.h"
#include "dbindex_paired.h"
#include "mask.h"
#include "minheap.h"
#include "otutable.h"
#include "searchcore_paired.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/xpthread.hpp"
#include "vsearch.h"

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <pthread.h>
#include <utility>
#include <vector>

namespace {

auto constexpr invalid_target_seqno_paired =
    std::numeric_limits<unsigned int>::max();

static int tophits_paired = 0;
static int seqcount_paired = 0;
static int clusters_paired = 0;
static int longest_end_paired = 0;
static int64_t nucleotide_count_paired = 0;

static std::vector<record_paired_s> *records_paired = nullptr;
static std::vector<unsigned int> target_lengths_paired;
static std::vector<unsigned int> target_seqnos_r1_paired;
static std::vector<unsigned int> target_seqnos_r2_paired;

struct clusterinfo_s_paired {
  int seqno = 0;
  int clusterno = 0;
  char *cigar_r1 = nullptr;
  char *cigar_r2 = nullptr;
  int target_seqno = -1;
  int strand = 0;
  double id = 0.0;
  bool perfect = false;
};

static std::vector<clusterinfo_s_paired> clusterinfo_paired;
static std::vector<int> centroid_seqnos_paired;
static std::vector<int64_t> const *cluster_abundance_for_sort_paired = nullptr;
static std::FILE *fp_uc_paired = nullptr;

static pthread_attr_t attr_paired;

struct thread_info_s_paired {
  pthread_t thread;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  int work = 0;
  int query_first = 0;
  int query_count = 0;
};

static thread_info_s_paired *ti_paired = nullptr;
static searchinfo_s_paired *si_plus_paired = nullptr;
static searchinfo_s_paired *si_minus_paired = nullptr;

inline auto compare_byclusterno_paired(void const *a, void const *b) -> int {
  auto const *lhs = static_cast<clusterinfo_s_paired const *>(a);
  auto const *rhs = static_cast<clusterinfo_s_paired const *>(b);

  if (lhs->clusterno < rhs->clusterno) {
    return -1;
  }
  if (lhs->clusterno > rhs->clusterno) {
    return +1;
  }
  if (lhs->seqno < rhs->seqno) {
    return -1;
  }
  if (lhs->seqno > rhs->seqno) {
    return +1;
  }
  return 0;
}

inline auto compare_byclusterabundance_paired(void const *a, void const *b)
    -> int {
  auto const *lhs = static_cast<clusterinfo_s_paired const *>(a);
  auto const *rhs = static_cast<clusterinfo_s_paired const *>(b);

  auto const lhs_abundance =
      (*cluster_abundance_for_sort_paired)[static_cast<std::size_t>(
          lhs->clusterno)];
  auto const rhs_abundance =
      (*cluster_abundance_for_sort_paired)[static_cast<std::size_t>(
          rhs->clusterno)];

  if (lhs_abundance > rhs_abundance) {
    return -1;
  }
  if (lhs_abundance < rhs_abundance) {
    return +1;
  }
  if (lhs->clusterno < rhs->clusterno) {
    return -1;
  }
  if (lhs->clusterno > rhs->clusterno) {
    return +1;
  }
  if (lhs->seqno < rhs->seqno) {
    return -1;
  }
  if (lhs->seqno > rhs->seqno) {
    return +1;
  }
  return 0;
}

auto cluster_query_core_paired(searchinfo_s_paired *si) -> void {
  auto const &record = (*records_paired)[si->query_no];
  si->query_head_len = static_cast<int>(record.header.size());
  si->query_head = const_cast<char *>(record.header.c_str());
  si->qsize = record.abundance;

  if (si->strand != 0) {
    si->qseqlen_r1 = static_cast<int>(record.qsequence_r2.size());
    si->qseqlen_r2 = static_cast<int>(record.qsequence_r1.size());
    reverse_complement(si->qsequence_r1, record.qsequence_r2.c_str(),
                       si->qseqlen_r1);
    reverse_complement(si->qsequence_r2, record.qsequence_r1.c_str(),
                       si->qseqlen_r2);
  } else {
    si->qseqlen_r1 = static_cast<int>(record.qsequence_r1.size());
    si->qseqlen_r2 = static_cast<int>(record.qsequence_r2.size());
    std::strcpy(si->qsequence_r1, record.qsequence_r1.c_str());
    std::strcpy(si->qsequence_r2, record.qsequence_r2.c_str());
  }

  search_onequery_paired(si, opt_qmask);
}

auto cluster_worker_paired(int64_t const t) -> void {
  for (int q = 0; q < ti_paired[t].query_count; q++) {
    cluster_query_core_paired(si_plus_paired + ti_paired[t].query_first + q);
    if (opt_strand > 1) {
      cluster_query_core_paired(si_minus_paired + ti_paired[t].query_first + q);
    }
  }
}

auto threads_worker_paired(void *vp) -> void * {
  auto const t = (int64_t)vp;
  auto *tip = ti_paired + t;
  xpthread_mutex_lock(&tip->mutex);
  while (tip->work >= 0) {
    if (tip->work == 0) {
      xpthread_cond_wait(&tip->cond, &tip->mutex);
    }
    if (tip->work > 0) {
      cluster_worker_paired(t);
      tip->work = 0;
      xpthread_cond_signal(&tip->cond);
    }
  }
  xpthread_mutex_unlock(&tip->mutex);
  return nullptr;
}

auto threads_wakeup_paired(int const queries) -> void {
  int const threads = queries > opt_threads ? opt_threads : queries;
  int queries_rest = queries;
  int threads_rest = threads;
  int query_next = 0;

  for (int t = 0; t < threads; t++) {
    auto *tip = ti_paired + t;

    tip->query_first = query_next;
    tip->query_count = (queries_rest + threads_rest - 1) / threads_rest;
    queries_rest -= tip->query_count;
    query_next += tip->query_count;
    --threads_rest;

    xpthread_mutex_lock(&tip->mutex);
    tip->work = 1;
    xpthread_cond_signal(&tip->cond);
    xpthread_mutex_unlock(&tip->mutex);
  }

  for (int t = 0; t < threads; t++) {
    auto *tip = ti_paired + t;
    xpthread_mutex_lock(&tip->mutex);
    while (tip->work > 0) {
      xpthread_cond_wait(&tip->cond, &tip->mutex);
    }
    xpthread_mutex_unlock(&tip->mutex);
  }
}

auto threads_init_paired() -> void {
  xpthread_attr_init(&attr_paired);
  xpthread_attr_setdetachstate(&attr_paired, PTHREAD_CREATE_JOINABLE);

  ti_paired = static_cast<thread_info_s_paired *>(
      xmalloc(opt_threads * sizeof(thread_info_s_paired)));

  for (int t = 0; t < opt_threads; t++) {
    auto *tip = ti_paired + t;
    tip->work = 0;
    xpthread_mutex_init(&tip->mutex, nullptr);
    xpthread_cond_init(&tip->cond, nullptr);
    xpthread_create(&tip->thread, &attr_paired, threads_worker_paired,
                    (void *)(int64_t)t);
  }
}

auto threads_exit_paired() -> void {
  for (int t = 0; t < opt_threads; t++) {
    auto *tip = ti_paired + t;

    xpthread_mutex_lock(&tip->mutex);
    tip->work = -1;
    xpthread_cond_signal(&tip->cond);
    xpthread_mutex_unlock(&tip->mutex);

    xpthread_join(tip->thread, nullptr);
    xpthread_cond_destroy(&tip->cond);
    xpthread_mutex_destroy(&tip->mutex);
  }

  if (ti_paired != nullptr) {
    xfree(ti_paired);
    ti_paired = nullptr;
  }
  xpthread_attr_destroy(&attr_paired);
}

auto cluster_query_init_paired(searchinfo_s_paired *si) -> void {
  si->query_no = 0;
  si->strand = 0;
  si->qsize = 1;
  si->query_head_len = 0;
  si->query_head = nullptr;
  si->qseqlen_r1 = 0;
  si->qseqlen_r2 = 0;
  si->kmersamplecount_r1 = 0;
  si->kmersamplecount_r2 = 0;
  si->kmersample_r1 = nullptr;
  si->kmersample_r2 = nullptr;
  si->hit_count = 0;
  si->accepts = 0;
  si->rejects = 0;
  si->finalized = 0;
  si->lma_r1 = nullptr;
  si->lma_r2 = nullptr;
  si->seq_alloc = longest_end_paired + 1;
  si->qsequence_r1 = static_cast<char *>(xmalloc(si->seq_alloc));
  si->qsequence_r2 = static_cast<char *>(xmalloc(si->seq_alloc));
  si->kmers =
      static_cast<count_t *>(xmalloc((seqcount_paired * sizeof(count_t)) + 32));
  si->hits = static_cast<hit_paired_s *>(
      xmalloc(sizeof(hit_paired_s) * static_cast<std::size_t>(tophits_paired)));

  si->uh_r1 = unique_init();
  si->uh_r2 = unique_init();
  si->m = minheap_init(tophits_paired);
  si->s_r1 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);
  si->s_r2 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);

  si->target_seqnos_r1 = &target_seqnos_r1_paired;
  si->target_seqnos_r2 = &target_seqnos_r2_paired;
}

auto cluster_query_exit_paired(searchinfo_s_paired *si) -> void {
  if (si->s_r1 != nullptr) {
    search16_exit(si->s_r1);
    si->s_r1 = nullptr;
  }
  if (si->s_r2 != nullptr) {
    search16_exit(si->s_r2);
    si->s_r2 = nullptr;
  }
  if (si->uh_r1 != nullptr) {
    unique_exit(si->uh_r1);
    si->uh_r1 = nullptr;
  }
  if (si->uh_r2 != nullptr) {
    unique_exit(si->uh_r2);
    si->uh_r2 = nullptr;
  }
  if (si->m != nullptr) {
    minheap_exit(si->m);
    si->m = nullptr;
  }
  if (si->qsequence_r1 != nullptr) {
    xfree(si->qsequence_r1);
    si->qsequence_r1 = nullptr;
  }
  if (si->qsequence_r2 != nullptr) {
    xfree(si->qsequence_r2);
    si->qsequence_r2 = nullptr;
  }
  if (si->hits != nullptr) {
    xfree(si->hits);
    si->hits = nullptr;
  }
  if (si->kmers != nullptr) {
    xfree(si->kmers);
    si->kmers = nullptr;
  }
  si->query_head_len = 0;
  si->query_head = nullptr;
}

auto relabel_otu_paired(int const clusterno, record_paired_s const &record)
    -> std::string {
  if (opt_relabel != nullptr) {
    return std::string{opt_relabel} + std::to_string(clusterno + 1);
  }
  if (opt_relabel_self) {
    return record.qsequence_r1;
  }
  if (opt_relabel_sha1) {
    std::vector<char> label(len_hex_dig_sha1);
    get_hex_seq_digest_sha1(label.data(), record.qsequence_r1.c_str(),
                            static_cast<int>(record.qsequence_r1.size()));
    return label.data();
  }
  if (opt_relabel_md5) {
    std::vector<char> label(len_hex_dig_md5);
    get_hex_seq_digest_md5(label.data(), record.qsequence_r1.c_str(),
                           static_cast<int>(record.qsequence_r1.size()));
    return label.data();
  }
  return record.header;
}

auto cluster_core_results_hit_paired(hit_paired_s *best, int const clusterno,
                                     char *query_head, int const query_head_len,
                                     int const qseqlen, int64_t const qsize)
    -> void {
  if ((opt_otutabout != nullptr) or (opt_mothur_shared_out != nullptr) or
      (opt_biomout != nullptr)) {
      if ((opt_relabel != nullptr) or opt_relabel_self or opt_relabel_sha1 or
          opt_relabel_md5) {
      auto const label =
          relabel_otu_paired(clusterno, (*records_paired)[best->target]);
      otutable_add(query_head, label.c_str(), qsize);
    } else {
      otutable_add(query_head,
                   (*records_paired)[best->target].header.c_str(), qsize);
    }
  }

  if (fp_uc_paired != nullptr) {
    auto const perfect = (best->mismatches == 0) and
                         (best->internal_indels == 0) and
                         (best->internal_gaps == 0);
    std::fprintf(fp_uc_paired, "H\t%d\t%d\t%.1f\t%c\t0\t0\t%s\t", clusterno,
                 qseqlen, best->id, (best->strand != 0) ? '-' : '+',
                 perfect ? "=" : "*");
    header_fprint_strip(fp_uc_paired, query_head, query_head_len, opt_xsize,
                        opt_xee, opt_xlength);
    std::fprintf(fp_uc_paired, "\t");
    header_fprint_strip(
        fp_uc_paired, (*records_paired)[best->target].header.c_str(),
        static_cast<int64_t>((*records_paired)[best->target].header.size()),
        opt_xsize, opt_xee, opt_xlength);
    std::fprintf(fp_uc_paired, "\n");
  }
}

auto cluster_core_results_nohit_paired(int const clusterno,
                                       char *query_head,
                                       int const query_head_len,
                                       int const qseqlen,
                                       record_paired_s const &record,
                                       int64_t const qsize) -> void {
  if ((opt_otutabout != nullptr) or (opt_mothur_shared_out != nullptr) or
      (opt_biomout != nullptr)) {
    if ((opt_relabel != nullptr) or opt_relabel_self or opt_relabel_sha1 or
        opt_relabel_md5) {
      auto const label = relabel_otu_paired(clusterno, record);
      otutable_add(query_head, label.c_str(), qsize);
    } else {
      otutable_add(query_head, query_head, qsize);
    }
  }

  if (fp_uc_paired != nullptr) {
    std::fprintf(fp_uc_paired, "S\t%d\t%d\t*\t*\t*\t*\t*\t", clusterno,
                 qseqlen);
    header_fprint_strip(fp_uc_paired, query_head, query_head_len, opt_xsize,
                        opt_xee, opt_xlength);
    std::fprintf(fp_uc_paired, "\t*\n");
  }
}

auto clear_hits_paired(searchinfo_s_paired *si) -> void {
  for (int i = 0; i < si->hit_count; ++i) {
    if (si->hits[i].aligned) {
      if (si->hits[i].r1.nwalignment != nullptr) {
        xfree(si->hits[i].r1.nwalignment);
        si->hits[i].r1.nwalignment = nullptr;
      }
      if (si->hits[i].r2.nwalignment != nullptr) {
        xfree(si->hits[i].r2.nwalignment);
        si->hits[i].r2.nwalignment = nullptr;
      }
    }
    si->hits[i] = hit_paired_s{};
  }
  si->hit_count = 0;
}

auto cluster_core_parallel_paired() -> void {
  threads_init_paired();

  constexpr static int queries_per_thread = 1;
  int const max_queries = queries_per_thread * opt_threads;

  si_plus_paired = static_cast<searchinfo_s_paired *>(
      xmalloc(max_queries * sizeof(searchinfo_s_paired)));

  if (opt_strand > 1) {
    si_minus_paired = static_cast<searchinfo_s_paired *>(
        xmalloc(max_queries * sizeof(searchinfo_s_paired)));
  }

  for (int i = 0; i < max_queries; i++) {
    cluster_query_init_paired(si_plus_paired + i);
    si_plus_paired[i].strand = 0;
    if (opt_strand > 1) {
      cluster_query_init_paired(si_minus_paired + i);
      si_minus_paired[i].strand = 1;
    }
  }

  std::vector<int> extra_list(static_cast<std::size_t>(max_queries));

  struct Scoring scoring;
  scoring.match = opt_match;
  scoring.mismatch = opt_mismatch;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_open_query_left = opt_gap_open_query_left;
  scoring.gap_open_target_left = opt_gap_open_target_left;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_open_target_interior = opt_gap_open_target_interior;
  scoring.gap_open_query_right = opt_gap_open_query_right;
  scoring.gap_open_target_right = opt_gap_open_target_right;
  scoring.gap_extension_query_left = opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opt_gap_extension_target_left;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_extension_target_interior = opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opt_gap_extension_target_right;

  LinearMemoryAligner lma_r1(scoring);
  LinearMemoryAligner lma_r2(scoring);

  auto lastlength = std::numeric_limits<int>::max();
  int seqno = 0;
  int64_t sum_nucleotides = 0;

  progress_init("Clustering", static_cast<uint64_t>(nucleotide_count_paired));

  while (seqno < seqcount_paired) {
    int queries = 0;

    for (int i = 0; i < max_queries; i++) {
      if (seqno < seqcount_paired) {
        auto const length =
            static_cast<int>(target_lengths_paired[static_cast<std::size_t>(
                seqno)]);

#if 1
        if ((opt_cluster_smallmem != nullptr) and (opt_usersort == 0) and
            (length > lastlength)) {
          fatal("Sequences not sorted by length and --usersort not specified.");
        }
#endif

        lastlength = length;

        si_plus_paired[i].query_no = seqno;
        si_plus_paired[i].strand = 0;

        if (opt_strand > 1) {
          si_minus_paired[i].query_no = seqno;
          si_minus_paired[i].strand = 1;
        }

        ++queries;
        ++seqno;
      }
    }

    threads_wakeup_paired(queries);

    int extra_count = 0;

    for (int i = 0; i < queries; i++) {
      auto *si_p = si_plus_paired + i;
      auto *si_m = opt_strand > 1 ? si_minus_paired + i : nullptr;

      for (int s = 0; s < opt_strand; s++) {
        auto *si = (s != 0) ? si_m : si_p;

        int added = 0;

        if (extra_count != 0) {
          for (int j = 0; j < extra_count; j++) {
            auto const *sic = si_plus_paired + extra_list[j];
            auto const shared =
                unique_count_shared(*si->uh_r1, static_cast<int>(opt_wordlength),
                                    static_cast<int>(sic->kmersamplecount_r1),
                                    sic->kmersample_r1) +
                unique_count_shared(*si->uh_r2, static_cast<int>(opt_wordlength),
                                    static_cast<int>(sic->kmersamplecount_r2),
                                    sic->kmersample_r2);

            if (search_enough_kmers_paired(*si, shared)) {
              auto const seqno_r1 = target_seqnos_r1_paired[sic->query_no];
              auto const seqno_r2 = target_seqnos_r2_paired[sic->query_no];
              auto const length = static_cast<unsigned int>(
                  db_getsequencelen(seqno_r1) + db_getsequencelen(seqno_r2));

              int x = si->hit_count;
              while (
                  (x > 0) and
                  ((si->hits[x - 1].count < shared) or
                   ((si->hits[x - 1].count == shared) and
                    (static_cast<unsigned int>(
                         db_getsequencelen(
                             target_seqnos_r1_paired[si->hits[x - 1].target]) +
                         db_getsequencelen(
                             target_seqnos_r2_paired[si->hits[x - 1].target])) >
                     length)))) {
                --x;
              }

              if (x < opt_maxaccepts + opt_maxrejects - 1) {
                if (si->hit_count >= opt_maxaccepts + opt_maxrejects - 1) {
                  if (si->hits[si->hit_count - 1].aligned) {
                    xfree(si->hits[si->hit_count - 1].r1.nwalignment);
                    xfree(si->hits[si->hit_count - 1].r2.nwalignment);
                  }
                  si->hits[si->hit_count - 1] = hit_paired_s{};
                  --si->hit_count;
                }

                for (int z = si->hit_count; z > x; z--) {
                  si->hits[z] = si->hits[z - 1];
                }

                auto *hit = si->hits + x;
                ++si->hit_count;

                *hit = hit_paired_s{};
                hit->target = sic->query_no;
                hit->target_seqno_r1 = seqno_r1;
                hit->target_seqno_r2 = seqno_r2;
                hit->strand = si->strand;
                hit->count = shared;

                ++added;
              }
            }
          }
        }

        if (added != 0) {
          si->rejects = 0;
          si->accepts = 0;

          for (int t = 0; t < si->hit_count; t++) {
            si->hits[t].accepted = false;
            si->hits[t].rejected = false;
          }

          for (int t = 0;
               (si->accepts < opt_maxaccepts) and
               (si->rejects < opt_maxrejects) and (t < si->hit_count);
               ++t) {
            auto *hit = si->hits + t;

            if (not hit->aligned) {
              auto const target = hit->target;
              if (search_acceptable_unaligned_paired(*si, target)) {
                auto const target_r1 = target_seqnos_r1_paired[target];
                auto const target_r2 = target_seqnos_r2_paired[target];

                unsigned int nwtarget_r1 = target_r1;
                unsigned int nwtarget_r2 = target_r2;

                int64_t nwscore_r1 = 0;
                int64_t nwscore_r2 = 0;
                int64_t nwalignmentlength_r1 = 0;
                int64_t nwalignmentlength_r2 = 0;
                int64_t nwmatches_r1 = 0;
                int64_t nwmatches_r2 = 0;
                int64_t nwmismatches_r1 = 0;
                int64_t nwmismatches_r2 = 0;
                int64_t nwgaps_r1 = 0;
                int64_t nwgaps_r2 = 0;

                char *nwcigar_r1 = nullptr;
                char *nwcigar_r2 = nullptr;

                CELL snwscore_r1 = 0;
                CELL snwscore_r2 = 0;
                unsigned short snwalignmentlength_r1 = 0;
                unsigned short snwalignmentlength_r2 = 0;
                unsigned short snwmatches_r1 = 0;
                unsigned short snwmatches_r2 = 0;
                unsigned short snwmismatches_r1 = 0;
                unsigned short snwmismatches_r2 = 0;
                unsigned short snwgaps_r1 = 0;
                unsigned short snwgaps_r2 = 0;

                search16(si->s_r1, 1, &nwtarget_r1, &snwscore_r1,
                         &snwalignmentlength_r1, &snwmatches_r1,
                         &snwmismatches_r1, &snwgaps_r1, &nwcigar_r1);
                search16(si->s_r2, 1, &nwtarget_r2, &snwscore_r2,
                         &snwalignmentlength_r2, &snwmatches_r2,
                         &snwmismatches_r2, &snwgaps_r2, &nwcigar_r2);

                auto const tseqlen_r1 = db_getsequencelen(target_r1);
                auto const tseqlen_r2 = db_getsequencelen(target_r2);

                if (snwscore_r1 == std::numeric_limits<short>::max()) {
                  auto *tseq = db_getsequence(target_r1);

                  if (nwcigar_r1 != nullptr) {
                    xfree(nwcigar_r1);
                  }

                  nwcigar_r1 = xstrdup(lma_r1.align(
                      si->qsequence_r1, tseq, si->qseqlen_r1, tseqlen_r1));

                  lma_r1.alignstats(nwcigar_r1, si->qsequence_r1, tseq,
                                    &nwscore_r1, &nwalignmentlength_r1,
                                    &nwmatches_r1, &nwmismatches_r1,
                                    &nwgaps_r1);
                } else {
                  nwscore_r1 = snwscore_r1;
                  nwalignmentlength_r1 = snwalignmentlength_r1;
                  nwmatches_r1 = snwmatches_r1;
                  nwmismatches_r1 = snwmismatches_r1;
                  nwgaps_r1 = snwgaps_r1;
                }

                if (snwscore_r2 == std::numeric_limits<short>::max()) {
                  auto *tseq = db_getsequence(target_r2);

                  if (nwcigar_r2 != nullptr) {
                    xfree(nwcigar_r2);
                  }

                  nwcigar_r2 = xstrdup(lma_r2.align(
                      si->qsequence_r2, tseq, si->qseqlen_r2, tseqlen_r2));

                  lma_r2.alignstats(nwcigar_r2, si->qsequence_r2, tseq,
                                    &nwscore_r2, &nwalignmentlength_r2,
                                    &nwmatches_r2, &nwmismatches_r2,
                                    &nwgaps_r2);
                } else {
                  nwscore_r2 = snwscore_r2;
                  nwalignmentlength_r2 = snwalignmentlength_r2;
                  nwmatches_r2 = snwmatches_r2;
                  nwmismatches_r2 = snwmismatches_r2;
                  nwgaps_r2 = snwgaps_r2;
                }

                if (nwcigar_r1 == nullptr) {
                  nwcigar_r1 = xstrdup("");
                }
                if (nwcigar_r2 == nullptr) {
                  nwcigar_r2 = xstrdup("");
                }

                hit->aligned = true;
                hit->r1.nwalignment = nwcigar_r1;
                hit->r1.nwscore = static_cast<int>(nwscore_r1);
                hit->r1.nwdiff =
                    static_cast<int>(nwalignmentlength_r1 - nwmatches_r1);
                hit->r1.matches = static_cast<int>(nwmatches_r1);
                hit->r1.mismatches = static_cast<int>(nwmismatches_r1);
                hit->r1.nwgaps = static_cast<int>(nwgaps_r1);
                hit->r1.nwindels = static_cast<int>(
                    nwalignmentlength_r1 - nwmatches_r1 - nwmismatches_r1);
                hit->r1.nwalignmentlength =
                    static_cast<int>(nwalignmentlength_r1);
                hit->r1.nwid = (hit->r1.nwalignmentlength != 0)
                                   ? (100.0 * (hit->r1.nwalignmentlength -
                                               hit->r1.nwdiff) /
                                      hit->r1.nwalignmentlength)
                                   : 0.0;
                hit->r1.shortest =
                    std::min(si->qseqlen_r1, static_cast<int>(tseqlen_r1));
                hit->r1.longest =
                    std::max(si->qseqlen_r1, static_cast<int>(tseqlen_r1));
                align_trim(&hit->r1);

                hit->r2.nwalignment = nwcigar_r2;
                hit->r2.nwscore = static_cast<int>(nwscore_r2);
                hit->r2.nwdiff =
                    static_cast<int>(nwalignmentlength_r2 - nwmatches_r2);
                hit->r2.matches = static_cast<int>(nwmatches_r2);
                hit->r2.mismatches = static_cast<int>(nwmismatches_r2);
                hit->r2.nwgaps = static_cast<int>(nwgaps_r2);
                hit->r2.nwindels = static_cast<int>(
                    nwalignmentlength_r2 - nwmatches_r2 - nwmismatches_r2);
                hit->r2.nwalignmentlength =
                    static_cast<int>(nwalignmentlength_r2);
                hit->r2.nwid = (hit->r2.nwalignmentlength != 0)
                                   ? (100.0 * (hit->r2.nwalignmentlength -
                                               hit->r2.nwdiff) /
                                      hit->r2.nwalignmentlength)
                                   : 0.0;
                hit->r2.shortest =
                    std::min(si->qseqlen_r2, static_cast<int>(tseqlen_r2));
                hit->r2.longest =
                    std::max(si->qseqlen_r2, static_cast<int>(tseqlen_r2));
                align_trim(&hit->r2);

                auto const qseqlen = si->qseqlen_r1 + si->qseqlen_r2;
                auto const tseqlen =
                    static_cast<int>(tseqlen_r1 + tseqlen_r2);

                hit->matches = hit->r1.matches + hit->r2.matches;
                hit->mismatches = hit->r1.mismatches + hit->r2.mismatches;
                hit->nwgaps = hit->r1.nwgaps + hit->r2.nwgaps;
                hit->nwindels = hit->r1.nwindels + hit->r2.nwindels;
                hit->nwalignmentlength =
                    hit->r1.nwalignmentlength + hit->r2.nwalignmentlength;
                hit->internal_alignmentlength =
                    hit->r1.internal_alignmentlength +
                    hit->r2.internal_alignmentlength;
                hit->internal_gaps = hit->r1.internal_gaps + hit->r2.internal_gaps;
                hit->internal_indels =
                    hit->r1.internal_indels + hit->r2.internal_indels;
                hit->shortest = std::min(qseqlen, tseqlen);
                hit->longest = std::max(qseqlen, tseqlen);
                hit->nwdiff = hit->nwalignmentlength - hit->matches;
                hit->nwscore = static_cast<int>(nwscore_r1 + nwscore_r2);
                if (hit->nwalignmentlength != 0) {
                  hit->nwid = 100.0 * (hit->nwalignmentlength - hit->nwdiff) /
                              hit->nwalignmentlength;
                } else {
                  hit->nwid = 0.0;
                }

                hit->mismatches_total = hit->mismatches;
                hit->nwgaps_total = hit->nwgaps;
                hit->nwalignment_cols_total = hit->nwalignmentlength;
                hit->internal_alignment_cols_total =
                    hit->internal_alignmentlength;
                hit->internal_gaps_total = hit->internal_gaps;
                hit->internal_indels_total = hit->internal_indels;
              } else {
                hit->rejected = true;
                ++si->rejects;
              }
            }

            if (not hit->rejected) {
              if (search_acceptable_aligned_paired(*si, hit)) {
                ++si->accepts;
              } else {
                ++si->rejects;
              }
            }
          }

          int const old_hit_count = si->hit_count;
          int new_hit_count = si->hit_count;
          for (int t = si->hit_count - 1; t >= 0; t--) {
            auto const *hit = si->hits + t;
            if (not hit->accepted and not hit->rejected) {
              new_hit_count = t;
            }
          }
          for (int t = new_hit_count; t < old_hit_count; ++t) {
            if (si->hits[t].aligned) {
              xfree(si->hits[t].r1.nwalignment);
              xfree(si->hits[t].r2.nwalignment);
            }
            si->hits[t] = hit_paired_s{};
          }
          si->hit_count = new_hit_count;
        }
      }

      hit_paired_s *best = nullptr;
      if (opt_sizeorder) {
        best = search_findbest2_bysize_paired(si_p, si_m);
      } else {
        best = search_findbest2_byid_paired(si_p, si_m);
      }

      int const myseqno = si_p->query_no;
      auto const qseqlen = static_cast<int>(
          target_lengths_paired[static_cast<std::size_t>(myseqno)]);

      if (best != nullptr) {
        auto const target = best->target;
        auto const perfect = (best->mismatches == 0) and
                             (best->internal_indels == 0) and
                             (best->internal_gaps == 0);

        cluster_core_results_hit_paired(best,
                                        clusterinfo_paired[target].clusterno,
                                        si_p->query_head,
                                        si_p->query_head_len, qseqlen,
                                        si_p->qsize);

        clusterinfo_paired[myseqno].seqno = myseqno;
        clusterinfo_paired[myseqno].clusterno =
            clusterinfo_paired[target].clusterno;
        clusterinfo_paired[myseqno].cigar_r1 = best->r1.nwalignment;
        clusterinfo_paired[myseqno].cigar_r2 = best->r2.nwalignment;
        clusterinfo_paired[myseqno].target_seqno = target;
        clusterinfo_paired[myseqno].strand = best->strand;
        clusterinfo_paired[myseqno].id = best->id;
        clusterinfo_paired[myseqno].perfect = perfect;
        best->r1.nwalignment = nullptr;
        best->r2.nwalignment = nullptr;
      } else {
        extra_list[extra_count] = i;
        ++extra_count;

        clusterinfo_paired[myseqno].seqno = myseqno;
        clusterinfo_paired[myseqno].clusterno = clusters_paired;
        clusterinfo_paired[myseqno].cigar_r1 = nullptr;
        clusterinfo_paired[myseqno].cigar_r2 = nullptr;
        clusterinfo_paired[myseqno].target_seqno = -1;
        clusterinfo_paired[myseqno].strand = 0;
        clusterinfo_paired[myseqno].id = 0.0;
        clusterinfo_paired[myseqno].perfect = false;
        centroid_seqnos_paired.push_back(myseqno);

        if (target_seqnos_r1_paired[myseqno] == invalid_target_seqno_paired) {
          auto const &record = (*records_paired)[myseqno];
          auto const left_seqno =
              static_cast<unsigned int>(db_getsequencecount());

          db_add(false, record.header.c_str(), record.qsequence_r1.c_str(),
                 nullptr, record.header.size(), record.qsequence_r1.size(),
                 std::max<int64_t>(record.abundance, 1));
          db_add(false, record.header_r2.c_str(), record.qsequence_r2.c_str(),
                 nullptr, record.header_r2.size(), record.qsequence_r2.size(),
                 std::max<int64_t>(record.abundance, 1));

          target_seqnos_r1_paired[myseqno] = left_seqno;
          target_seqnos_r2_paired[myseqno] = left_seqno + 1U;
        }

        dbindex_addsequence_paired(static_cast<unsigned int>(myseqno), opt_qmask);

        cluster_core_results_nohit_paired(clusters_paired, si_p->query_head,
                                          si_p->query_head_len, qseqlen,
                                          (*records_paired)[myseqno],
                                          si_p->qsize);
        ++clusters_paired;
      }

      clear_hits_paired(si_p);
      if (si_m != nullptr) {
        clear_hits_paired(si_m);
      }

      sum_nucleotides += static_cast<int64_t>(target_lengths_paired[myseqno]);
    }

    progress_update(static_cast<uint64_t>(sum_nucleotides));
  }
  progress_done();

  for (int i = 0; i < max_queries; i++) {
    cluster_query_exit_paired(si_plus_paired + i);
    if (opt_strand > 1) {
      cluster_query_exit_paired(si_minus_paired + i);
    }
  }

  xfree(si_plus_paired);
  si_plus_paired = nullptr;
  if (si_minus_paired != nullptr) {
    xfree(si_minus_paired);
    si_minus_paired = nullptr;
  }

  threads_exit_paired();
}

auto cluster_core_serial_paired() -> void {
  std::array<searchinfo_s_paired, 1> si_p{{}};
  std::array<searchinfo_s_paired, 1> si_m{{}};

  cluster_query_init_paired(si_p.data());
  if (opt_strand > 1) {
    cluster_query_init_paired(si_m.data());
  }

  auto lastlength = std::numeric_limits<int>::max();
  progress_init("Clustering", static_cast<uint64_t>(seqcount_paired));
  for (int seqno = 0; seqno < seqcount_paired; seqno++) {
    auto const length =
        static_cast<int>(target_lengths_paired[static_cast<std::size_t>(seqno)]);

#if 1
    if ((opt_cluster_smallmem != nullptr) and (opt_usersort == 0) and
        (length > lastlength)) {
      fatal("Sequences not sorted by length and --usersort not specified.");
    }
#endif

    lastlength = length;

    si_p[0].query_no = seqno;
    si_p[0].strand = 0;
    cluster_query_core_paired(si_p.data());

    if (opt_strand > 1) {
      si_m[0].query_no = seqno;
      si_m[0].strand = 1;
      cluster_query_core_paired(si_m.data());
    }

    hit_paired_s *best = nullptr;
    if (opt_sizeorder) {
      best = search_findbest2_bysize_paired(si_p.data(), si_m.data());
    } else {
      best = search_findbest2_byid_paired(si_p.data(), si_m.data());
    }

    auto const qseqlen = static_cast<int>(
        target_lengths_paired[static_cast<std::size_t>(seqno)]);

    if (best != nullptr) {
      auto const target = best->target;
      auto const perfect = (best->mismatches == 0) and
                           (best->internal_indels == 0) and
                           (best->internal_gaps == 0);

      cluster_core_results_hit_paired(best,
                                      clusterinfo_paired[target].clusterno,
                                      si_p[0].query_head,
                                      si_p[0].query_head_len, qseqlen,
                                      si_p[0].qsize);

      clusterinfo_paired[seqno].seqno = seqno;
      clusterinfo_paired[seqno].clusterno =
          clusterinfo_paired[target].clusterno;
      clusterinfo_paired[seqno].cigar_r1 = best->r1.nwalignment;
      clusterinfo_paired[seqno].cigar_r2 = best->r2.nwalignment;
      clusterinfo_paired[seqno].target_seqno = target;
      clusterinfo_paired[seqno].strand = best->strand;
      clusterinfo_paired[seqno].id = best->id;
      clusterinfo_paired[seqno].perfect = perfect;
      best->r1.nwalignment = nullptr;
      best->r2.nwalignment = nullptr;
    } else {
      clusterinfo_paired[seqno].seqno = seqno;
      clusterinfo_paired[seqno].clusterno = clusters_paired;
      clusterinfo_paired[seqno].cigar_r1 = nullptr;
      clusterinfo_paired[seqno].cigar_r2 = nullptr;
      clusterinfo_paired[seqno].target_seqno = -1;
      clusterinfo_paired[seqno].strand = 0;
      clusterinfo_paired[seqno].id = 0.0;
      clusterinfo_paired[seqno].perfect = false;
      centroid_seqnos_paired.push_back(seqno);

      if (target_seqnos_r1_paired[seqno] == invalid_target_seqno_paired) {
        auto const &record = (*records_paired)[seqno];
        auto const left_seqno = static_cast<unsigned int>(db_getsequencecount());

        db_add(false, record.header.c_str(), record.qsequence_r1.c_str(),
               nullptr, record.header.size(), record.qsequence_r1.size(),
               std::max<int64_t>(record.abundance, 1));
        db_add(false, record.header_r2.c_str(), record.qsequence_r2.c_str(),
               nullptr, record.header_r2.size(), record.qsequence_r2.size(),
               std::max<int64_t>(record.abundance, 1));

        target_seqnos_r1_paired[seqno] = left_seqno;
        target_seqnos_r2_paired[seqno] = left_seqno + 1U;
      }

      dbindex_addsequence_paired(static_cast<unsigned int>(seqno), opt_qmask);

      cluster_core_results_nohit_paired(clusters_paired, si_p[0].query_head,
                                        si_p[0].query_head_len, qseqlen,
                                        (*records_paired)[seqno], si_p[0].qsize);
      ++clusters_paired;
    }

    clear_hits_paired(si_p.data());
    if (opt_strand > 1) {
      clear_hits_paired(si_m.data());
    }

    progress_update(static_cast<uint64_t>(seqno));
  }
  progress_done();

  cluster_query_exit_paired(si_p.data());
  if (opt_strand > 1) {
    cluster_query_exit_paired(si_m.data());
  }
}

auto cluster_paired(struct Parameters const &parameters, char *dbname,
                    char *cmdline, char *progheader) -> void {
  (void)dbname;
  (void)cmdline;
  (void)progheader;

  auto *left_h = fastx_open(parameters.opt_cluster_unoise);
  auto *right_h =
      parameters.opt_interleaved ? nullptr : fastx_open(parameters.opt_reverse);

  std::vector<record_paired_s> records;

  while (fastx_next(left_h, not opt_notrunclabels,
                    chrmap_no_change_vector.data())) {
    auto record = record_paired_s{};
    record.header = fastx_get_header(left_h);
    record.abundance = fastx_get_abundance(left_h);
    record.first_seen = static_cast<int64_t>(records.size());

    auto const left_len =
        static_cast<std::size_t>(fastx_get_sequence_length(left_h));
    record.qsequence_r1.assign(fastx_get_sequence(left_h), left_len);

    if (parameters.opt_interleaved) {
      if (!fastx_next(left_h, not opt_notrunclabels,
                      chrmap_no_change_vector.data())) {
        fatal("Odd number of records in interleaved paired FASTX input %s; "
              "expected left/right entries",
              parameters.opt_cluster_unoise);
      }
      record.header_r2 = fastx_get_header(left_h);
      if (paired_header_key_paired(record.header) !=
          paired_header_key_paired(record.header_r2)) {
        auto const message = std::string{"Paired FASTX headers differ ("} +
                             record.header + " vs " + record.header_r2 + ")";
        fatal(message.c_str());
      }
      auto const right_len =
          static_cast<std::size_t>(fastx_get_sequence_length(left_h));
      record.qsequence_r2.assign(fastx_get_sequence(left_h), right_len);
    } else {
      if (!fastx_next(right_h, not opt_notrunclabels,
                      chrmap_no_change_vector.data())) {
        fatal(
            "More forward records than reverse records in paired FASTX input");
      }
      record.header_r2 = fastx_get_header(right_h);
      if (paired_header_key_paired(record.header) !=
          paired_header_key_paired(record.header_r2)) {
        auto const message = std::string{"Paired FASTX headers differ ("} +
                             record.header + " vs " + record.header_r2 + ")";
        fatal(message.c_str());
      }
      auto const right_len =
          static_cast<std::size_t>(fastx_get_sequence_length(right_h));
      record.qsequence_r2.assign(fastx_get_sequence(right_h), right_len);
    }

    if ((opt_qmask == MASK_DUST) or
        ((opt_qmask == MASK_SOFT) and (opt_hardmask != 0))) {
      std::vector<char> left_buf(record.qsequence_r1.begin(),
                                 record.qsequence_r1.end());
      left_buf.push_back('\0');
      std::vector<char> right_buf(record.qsequence_r2.begin(),
                                  record.qsequence_r2.end());
      right_buf.push_back('\0');

      if (opt_qmask == MASK_DUST) {
        dust(left_buf.data(), static_cast<int>(record.qsequence_r1.size()));
        dust(right_buf.data(), static_cast<int>(record.qsequence_r2.size()));
      } else {
        hardmask(left_buf.data(), static_cast<int>(record.qsequence_r1.size()));
        hardmask(right_buf.data(),
                 static_cast<int>(record.qsequence_r2.size()));
      }

      record.qsequence_r1.assign(left_buf.data(), record.qsequence_r1.size());
      record.qsequence_r2.assign(right_buf.data(), record.qsequence_r2.size());
    }

    records.push_back(record);
  }

  if ((not parameters.opt_interleaved) and
      fastx_next(right_h, not opt_notrunclabels,
                 chrmap_no_change_vector.data())) {
    fatal("More reverse records than forward records in paired FASTX input");
  }

  if (right_h != nullptr) {
    fastx_close(right_h);
  }
  fastx_close(left_h);

  if (records.empty()) {
    fatal("Input for paired cluster_unoise is empty");
  }

  std::vector<record_paired_s> filtered_records;
  filtered_records.reserve(records.size());

  int64_t filtered_small = 0;
  int64_t filtered_large = 0;
  for (auto const &record : records) {
    if (record.abundance < opt_minsize) {
      ++filtered_small;
      continue;
    }
    if (record.abundance > opt_maxsize) {
      ++filtered_large;
      continue;
    }
    filtered_records.push_back(record);
  }

  if (filtered_records.empty()) {
    fatal("No TAV records pass --minsize/--maxsize filtering");
  }

  std::sort(filtered_records.begin(), filtered_records.end(),
            [](record_paired_s const &lhs, record_paired_s const &rhs) {
              if (lhs.abundance != rhs.abundance) {
                return lhs.abundance > rhs.abundance;
              }
              if (lhs.header != rhs.header) {
                return lhs.header < rhs.header;
              }
              if (lhs.qsequence_r1 != rhs.qsequence_r1) {
                return lhs.qsequence_r1 < rhs.qsequence_r1;
              }
              return lhs.qsequence_r2 < rhs.qsequence_r2;
            });

  records_paired = &filtered_records;
  seqcount_paired = static_cast<int>(filtered_records.size());
  clusters_paired = 0;
  longest_end_paired = 0;
  nucleotide_count_paired = 0;

  target_lengths_paired.assign(static_cast<std::size_t>(seqcount_paired), 0U);
  target_seqnos_r1_paired.assign(static_cast<std::size_t>(seqcount_paired),
                                 invalid_target_seqno_paired);
  target_seqnos_r2_paired.assign(static_cast<std::size_t>(seqcount_paired),
                                 invalid_target_seqno_paired);
  clusterinfo_paired.assign(static_cast<std::size_t>(seqcount_paired), {});
  centroid_seqnos_paired.clear();

  for (int seqno = 0; seqno < seqcount_paired; ++seqno) {
    auto const &record = filtered_records[static_cast<std::size_t>(seqno)];
    target_lengths_paired[static_cast<std::size_t>(seqno)] =
        static_cast<unsigned int>(record.qsequence_r1.size() +
                                  record.qsequence_r2.size());
    longest_end_paired =
        std::max(longest_end_paired,
                 static_cast<int>(std::max(record.qsequence_r1.size(),
                                           record.qsequence_r2.size())));
    nucleotide_count_paired += static_cast<int64_t>(record.qsequence_r1.size() +
                                                    record.qsequence_r2.size());
  }

  if ((opt_maxrejects == 0) or (opt_maxrejects > seqcount_paired)) {
    opt_maxrejects = seqcount_paired;
  }
  if ((opt_maxaccepts == 0) or (opt_maxaccepts > seqcount_paired)) {
    opt_maxaccepts = seqcount_paired;
  }

  tophits_paired = opt_maxrejects + opt_maxaccepts + MAXDELAYED;
  tophits_paired = std::min(tophits_paired, seqcount_paired);

  if (opt_log != nullptr) {
    uint64_t const slots = 1ULL
                           << (static_cast<uint64_t>(opt_wordlength) << 1ULL);
    std::fprintf(fp_log, "\n");
    std::fprintf(fp_log, "      Alphabet  nt\n");
    std::fprintf(fp_log, "    Word width  %" PRId64 "\n", opt_wordlength);
    std::fprintf(fp_log, "     Word ones  %" PRId64 "\n", opt_wordlength);
    std::fprintf(fp_log, "        Spaced  No\n");
    std::fprintf(fp_log, "        Hashed  No\n");
    std::fprintf(fp_log, "         Coded  No\n");
    std::fprintf(fp_log, "       Stepped  No\n");
    std::fprintf(fp_log, "         Slots  %" PRIu64 " (%.1fk)\n", slots,
                 slots / 1000.0);
    std::fprintf(fp_log, "       DBAccel  100%%\n");
    std::fprintf(fp_log, "\n");
  }

  show_rusage();

  db_free();
  db_setinfo(false, 0, 0, 0, std::numeric_limits<uint64_t>::max(), 0);
  dbindex_prepare_paired(records_paired, 1, opt_qmask);

  otutable_init();

  if (parameters.opt_uc != nullptr) {
    fp_uc_paired = fopen_output(parameters.opt_uc);
    if (fp_uc_paired == nullptr) {
      fatal("Unable to open uc file for writing");
    }
  }

  if (opt_threads == 1) {
    cluster_core_serial_paired();
  } else {
    cluster_core_parallel_paired();
  }

  std::vector<int64_t> cluster_abundance(
      static_cast<std::size_t>(clusters_paired), 0);
  for (int seqno = 0; seqno < seqcount_paired; ++seqno) {
    auto const clusterno =
        clusterinfo_paired[static_cast<std::size_t>(seqno)].clusterno;
    cluster_abundance[static_cast<std::size_t>(clusterno)] +=
        opt_sizein ? filtered_records[static_cast<std::size_t>(seqno)].abundance
                   : 1;
  }

  auto const minmax_elements =
      std::minmax_element(cluster_abundance.cbegin(), cluster_abundance.cend());
  auto const abundance_min =
      cluster_abundance.empty() ? 0 : *std::get<0>(minmax_elements);
  auto const abundance_max =
      cluster_abundance.empty() ? 0 : *std::get<1>(minmax_elements);
  int const singletons = std::count(cluster_abundance.cbegin(),
                                    cluster_abundance.cend(), int64_t{1});

  std::vector<clusterinfo_s_paired> clusterinfo_sorted = clusterinfo_paired;
  progress_init("Sorting clusters", static_cast<uint64_t>(clusters_paired));
  cluster_abundance_for_sort_paired = &cluster_abundance;
  if (opt_clusterout_sort != 0) {
    std::qsort(clusterinfo_sorted.data(),
               static_cast<std::size_t>(seqcount_paired),
               sizeof(clusterinfo_s_paired), compare_byclusterabundance_paired);
  } else {
    std::qsort(clusterinfo_sorted.data(),
               static_cast<std::size_t>(seqcount_paired),
               sizeof(clusterinfo_s_paired), compare_byclusterno_paired);
  }
  cluster_abundance_for_sort_paired = nullptr;
  progress_done();

  std::FILE *fp_centroids_left = nullptr;
  std::FILE *fp_fastaout_left = nullptr;
  std::FILE *fp_centroids_right = nullptr;
  bool const write_paired_centroids =
      (opt_centroids != nullptr) or (parameters.opt_fastaout != nullptr);
  bool share_left_outputs = false;

  if (opt_centroids != nullptr) {
    fp_centroids_left = fopen_output(opt_centroids);
    if (fp_centroids_left == nullptr) {
      fatal("Unable to open centroids file for writing");
    }
  }
  if (parameters.opt_fastaout != nullptr) {
    if ((opt_centroids != nullptr) and
        (std::strcmp(parameters.opt_fastaout, opt_centroids) == 0)) {
      fp_fastaout_left = fp_centroids_left;
      share_left_outputs = true;
    } else {
      fp_fastaout_left = fopen_output(parameters.opt_fastaout);
      if (fp_fastaout_left == nullptr) {
        fatal("Unable to open left-anchor FASTA output file for writing");
      }
    }
  }
  if (write_paired_centroids) {
    fp_centroids_right = fopen_output(parameters.opt_fastaout_rev);
    if (fp_centroids_right == nullptr) {
      fatal("Unable to open right-anchor FASTA output file for writing");
    }
  }

  progress_init("Writing clusters", static_cast<uint64_t>(seqcount_paired));

  std::FILE *fp_clusters = nullptr;
  std::vector<char> fn_clusters;
  static constexpr auto space_for_cluster_id = 25;
  if (opt_clusters != nullptr) {
    fn_clusters.resize(std::strlen(opt_clusters) + space_for_cluster_id);
  }

  std::vector<record_paired_s> centroid_records;
  centroid_records.reserve(static_cast<std::size_t>(clusters_paired));
  int lastcluster = -1;
  for (int i = 0; i < seqcount_paired; ++i) {
    auto const &info = clusterinfo_sorted[static_cast<std::size_t>(i)];
    auto const seqno = info.seqno;
    auto const clusterno = info.clusterno;
    auto const &record = filtered_records[static_cast<std::size_t>(seqno)];

    if (clusterno != lastcluster) {
      if ((opt_clusters != nullptr) and (lastcluster != -1)) {
        std::fclose(fp_clusters);
      }

      if (opt_clusters != nullptr) {
        std::snprintf(fn_clusters.data(), fn_clusters.size(), "%s%d",
                      opt_clusters, clusterno);
        fp_clusters = fopen_output(fn_clusters.data());
        if (fp_clusters == nullptr) {
          fatal("Unable to open clusters file for writing");
        }
      }

      auto centroid_record = filtered_records[static_cast<std::size_t>(
          centroid_seqnos_paired[static_cast<std::size_t>(clusterno)])];
      centroid_record.abundance =
          cluster_abundance[static_cast<std::size_t>(clusterno)];
      centroid_records.push_back(centroid_record);

      auto const clusterid = (opt_clusterout_id != 0) ? clusterno : -1;

      if (fp_centroids_left != nullptr) {
        fasta_print_general(
            fp_centroids_left, nullptr, centroid_record.qsequence_r1.c_str(),
            static_cast<int>(centroid_record.qsequence_r1.size()),
            centroid_record.header.c_str(),
            static_cast<int>(centroid_record.header.size()),
            static_cast<unsigned int>(centroid_record.abundance), clusterno + 1,
            -1.0, -1, clusterid, nullptr, 0.0);
      }
      if ((fp_fastaout_left != nullptr) and
          (fp_fastaout_left != fp_centroids_left)) {
        fasta_print_general(
            fp_fastaout_left, nullptr, centroid_record.qsequence_r1.c_str(),
            static_cast<int>(centroid_record.qsequence_r1.size()),
            centroid_record.header.c_str(),
            static_cast<int>(centroid_record.header.size()),
            static_cast<unsigned int>(centroid_record.abundance), clusterno + 1,
            -1.0, -1, clusterid, nullptr, 0.0);
      }
      if (fp_centroids_right != nullptr) {
        fasta_print_general(
            fp_centroids_right, nullptr, centroid_record.qsequence_r2.c_str(),
            static_cast<int>(centroid_record.qsequence_r2.size()),
            centroid_record.header_r2.c_str(),
            static_cast<int>(centroid_record.header_r2.size()),
            static_cast<unsigned int>(centroid_record.abundance), clusterno + 1,
            -1.0, -1, clusterid, nullptr, 0.0);
      }

      if (fp_uc_paired != nullptr) {
        std::fprintf(fp_uc_paired, "C\t%d\t%" PRId64 "\t*\t*\t*\t*\t*\t",
                     clusterno,
                     cluster_abundance[static_cast<std::size_t>(clusterno)]);
        header_fprint_strip(fp_uc_paired, centroid_record.header.c_str(),
                            static_cast<int64_t>(centroid_record.header.size()),
                            opt_xsize, opt_xee, opt_xlength);
        std::fprintf(fp_uc_paired, "\t*\n");
      }

      lastcluster = clusterno;
    }

    if (opt_clusters != nullptr) {
      fasta_print_general(fp_clusters, nullptr, record.qsequence_r1.c_str(),
                          static_cast<int>(record.qsequence_r1.size()),
                          record.header.c_str(),
                          static_cast<int>(record.header.size()),
                          static_cast<unsigned int>(record.abundance), 0, -1.0,
                          -1, -1, nullptr, 0.0);
      fasta_print_general(fp_clusters, nullptr, record.qsequence_r2.c_str(),
                          static_cast<int>(record.qsequence_r2.size()),
                          record.header_r2.c_str(),
                          static_cast<int>(record.header_r2.size()),
                          static_cast<unsigned int>(record.abundance), 0, -1.0,
                          -1, -1, nullptr, 0.0);
    }

    progress_update(static_cast<uint64_t>(i));
  }

  if ((opt_clusters != nullptr) and (lastcluster != -1)) {
    std::fclose(fp_clusters);
  }
  progress_done();

  if (parameters.opt_tabbedout != nullptr) {
    auto *fp = fopen_output(parameters.opt_tabbedout);
    if (fp == nullptr) {
      fatal("Unable to open TAV catalog output file for writing");
    }
    std::fprintf(fp, "tav_id\tabundance\theader\tleft_anchor\tright_anchor\n");
    for (auto const &record : centroid_records) {
      std::fprintf(fp, "%s\t%" PRId64 "\t%s\t%s\t%s\n", record.header.c_str(),
                   record.abundance, record.header.c_str(),
                   record.qsequence_r1.c_str(), record.qsequence_r2.c_str());
    }
    std::fclose(fp);
  }

  if (fp_centroids_right != nullptr) {
    std::fclose(fp_centroids_right);
  }
  if (fp_fastaout_left != nullptr and (not share_left_outputs)) {
    std::fclose(fp_fastaout_left);
  }
  if (fp_centroids_left != nullptr) {
    std::fclose(fp_centroids_left);
  }

  if (clusters_paired < 1) {
    if (not opt_quiet) {
      std::fprintf(stderr, "Clusters: 0\n");
      std::fprintf(stderr, "Singletons: 0\n");
    }
    if (opt_log != nullptr) {
      std::fprintf(fp_log, "Clusters: 0\n");
      std::fprintf(fp_log, "Singletons: 0\n");
    }
  } else {
    if (not opt_quiet) {
      std::fprintf(stderr,
                   "Clusters: %d Size min %" PRId64 ", max %" PRId64
                   ", avg %.1f\n",
                   clusters_paired, abundance_min, abundance_max,
                   1.0 * seqcount_paired / clusters_paired);
      std::fprintf(stderr,
                   "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
                   singletons, 100.0 * singletons / seqcount_paired,
                   100.0 * singletons / clusters_paired);
      std::fprintf(stderr,
                   "Paired cluster_unoise: %zu input pairs (%" PRId64
                   " filtered by minsize, %" PRId64
                   " filtered by maxsize) -> %zu centroids\n",
                   records.size(), filtered_small, filtered_large,
                   centroid_records.size());
    }

    if (opt_log != nullptr) {
      std::fprintf(fp_log,
                   "Clusters: %d Size min %" PRId64 ", max %" PRId64
                   ", avg %.1f\n",
                   clusters_paired, abundance_min, abundance_max,
                   1.0 * seqcount_paired / clusters_paired);
      std::fprintf(fp_log,
                   "Singletons: %d, %.1f%% of seqs, %.1f%% of clusters\n",
                   singletons, 100.0 * singletons / seqcount_paired,
                   100.0 * singletons / clusters_paired);
      std::fprintf(fp_log,
                   "Paired cluster_unoise: %zu input pairs (%" PRId64
                   " filtered by minsize, %" PRId64
                   " filtered by maxsize) -> %zu centroids\n",
                   records.size(), filtered_small, filtered_large,
                   centroid_records.size());
      std::fprintf(fp_log, "\n");
    }
  }

  if (opt_biomout != nullptr) {
    auto *fp_biom = fopen_output(opt_biomout);
    if (fp_biom == nullptr) {
      fatal("Unable to open OTU table (biom 1.0 format) output file for "
            "writing");
    }
    otutable_print_biomout(fp_biom);
    std::fclose(fp_biom);
  }

  if (opt_otutabout != nullptr) {
    auto *fp_otu = fopen_output(opt_otutabout);
    if (fp_otu == nullptr) {
      fatal("Unable to open OTU table (text format) output file for writing");
    }
    otutable_print_otutabout(fp_otu);
    std::fclose(fp_otu);
  }

  if (opt_mothur_shared_out != nullptr) {
    auto *fp_mothur = fopen_output(opt_mothur_shared_out);
    if (fp_mothur == nullptr) {
      fatal("Unable to open OTU table (mothur format) output file for "
            "writing");
    }
    otutable_print_mothur_shared_out(fp_mothur);
    std::fclose(fp_mothur);
  }

  otutable_done();

  if (not opt_quiet and (clusters_paired < 1)) {
    std::fprintf(stderr,
                 "Paired cluster_unoise: %zu input pairs (%" PRId64
                 " filtered by minsize, %" PRId64
                 " filtered by maxsize) -> 0 centroids\n",
                 records.size(), filtered_small, filtered_large);
  } else if ((opt_log != nullptr) and (clusters_paired < 1)) {
    std::fprintf(fp_log,
                 "Paired cluster_unoise: %zu input pairs (%" PRId64
                 " filtered by minsize, %" PRId64
                 " filtered by maxsize) -> 0 centroids\n\n",
                 records.size(), filtered_small, filtered_large);
  }

  if (fp_uc_paired != nullptr) {
    std::fclose(fp_uc_paired);
    fp_uc_paired = nullptr;
  }

  for (auto &clusterinfo : clusterinfo_paired) {
    if (clusterinfo.cigar_r1 != nullptr) {
      xfree(clusterinfo.cigar_r1);
      clusterinfo.cigar_r1 = nullptr;
    }
    if (clusterinfo.cigar_r2 != nullptr) {
      xfree(clusterinfo.cigar_r2);
      clusterinfo.cigar_r2 = nullptr;
    }
  }

  dbindex_free_paired();
  db_free();
  records_paired = nullptr;
  show_rusage();
}

} // namespace

auto cluster_unoise_paired(struct Parameters const &parameters, char *cmdline,
                           char *progheader) -> void {
  cluster_paired(parameters, parameters.opt_cluster_unoise, cmdline,
                 progheader);
}
