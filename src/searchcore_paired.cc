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

#include "searchcore_paired.h"
#include "align_simd.h"
#include "dbindex_paired.h"
#include "linmemalign.h"
#include "minheap.h"
#include "unique.h"
#include "utils/seqcmp.hpp"
#include "vsearch.h"

#include "utils/span.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>

namespace {

auto make_hits_span_paired(searchinfo_s_paired const *search_info)
    -> Span<hit_paired_s> {
  assert(search_info != nullptr);
  assert(search_info->hit_count >= 0);
  auto const length = static_cast<std::size_t>(search_info->hit_count);
  return Span<hit_paired_s>{search_info->hits, length};
}

auto count_number_of_hits_to_keep_paired(searchinfo_s_paired const *search_info)
    -> std::size_t {
  if (search_info == nullptr) {
    return std::size_t{0};
  }
  auto const hits = make_hits_span_paired(search_info);
  return static_cast<std::size_t>(std::count_if(
      hits.cbegin(), hits.cend(), [](hit_paired_s const &hit) -> bool {
        return hit.accepted or hit.weak;
      }));
}

auto copy_over_hits_to_be_kept_paired(std::vector<hit_paired_s> &hits,
                                      searchinfo_s_paired const *search_info)
    -> void {
  if (search_info == nullptr) {
    return;
  }
  for (auto const &hit : make_hits_span_paired(search_info)) {
    if (hit.accepted or hit.weak) {
      hits.emplace_back(hit);
    }
  }
}

auto free_rejected_alignments_paired(searchinfo_s_paired *search_info)
    -> void {
  if (search_info == nullptr) {
    return;
  }
  for (auto &hit : make_hits_span_paired(search_info)) {
    if (not(hit.accepted or hit.weak) and hit.aligned) {
      if (hit.r1.nwalignment != nullptr) {
        xfree(hit.r1.nwalignment);
        hit.r1.nwalignment = nullptr;
      }
      if (hit.r2.nwalignment != nullptr) {
        xfree(hit.r2.nwalignment);
        hit.r2.nwalignment = nullptr;
      }
    }
  }
}

} // end of anonymous namespace

/* per thread data */

inline auto hit_compare_byid_typed_paired(hit_paired_s const *lhs,
                                          hit_paired_s const *rhs) -> int {
  /*
    Order:
    accepted, then rejected (weak)
    high id, then low id
    early target, then late target
  */

  if (lhs->rejected < rhs->rejected) {
    return -1;
  }
  if (lhs->rejected > rhs->rejected) {
    return +1;
  }
  if (lhs->aligned > rhs->aligned) {
    return -1;
  }
  if (lhs->aligned < rhs->aligned) {
    return +1;
  }
  if (lhs->aligned == 0) {
    return 0;
  }
  if (lhs->id > rhs->id) {
    return -1;
  }
  if (lhs->id < rhs->id) {
    return +1;
  }
  if (lhs->target < rhs->target) {
    return -1;
  }
  if (lhs->target > rhs->target) {
    return +1;
  }
  return 0;
}

inline auto hit_compare_bysize_typed_paired(hit_paired_s const *lhs,
                                            hit_paired_s const *rhs) -> int {
  // high abundance, then low abundance
  // high id, then low id
  // early target, then late target

  if (lhs->rejected < rhs->rejected) {
    return -1;
  }
  if (lhs->rejected > rhs->rejected) {
    return +1;
  }
  if (lhs->rejected == 1) {
    return 0;
  }

  if (lhs->aligned > rhs->aligned) {
    return -1;
  }
  if (lhs->aligned < rhs->aligned) {
    return +1;
  }
  if (lhs->aligned == 0) {
    return 0;
  }

  auto const lhs_abundance = db_getabundance(lhs->target_seqno_r1);
  auto const rhs_abundance = db_getabundance(rhs->target_seqno_r1);
  if (lhs_abundance > rhs_abundance) {
    return -1;
  }
  if (lhs_abundance < rhs_abundance) {
    return +1;
  }

  if (lhs->id > rhs->id) {
    return -1;
  }
  if (lhs->id < rhs->id) {
    return +1;
  }

  if (lhs->target < rhs->target) {
    return -1;
  }
  if (lhs->target > rhs->target) {
    return +1;
  }
  return 0;
}

auto hit_compare_byid_paired(void const *lhs, void const *rhs) -> int {
  return hit_compare_byid_typed_paired((hit_paired_s const *)lhs,
                                       (hit_paired_s const *)rhs);
}

auto hit_compare_bysize_paired(void const *lhs, void const *rhs) -> int {
  return hit_compare_bysize_typed_paired((hit_paired_s const *)lhs,
                                         (hit_paired_s const *)rhs);
}

auto search_enough_kmers_paired(searchinfo_s_paired const &searchinfo,
                                unsigned int const count) -> bool {
  auto const kmersamplecount =
      searchinfo.kmersamplecount_r1 + searchinfo.kmersamplecount_r2;
  return (count >= opt_minwordmatches) or (count >= kmersamplecount);
}

auto search_topscores_paired(searchinfo_s_paired *searchinfo) -> void {
  /*
    Count the kmer hits in each database sequence and
    make a sorted list of a given number (th)
    of the database sequences with the highest number of matching kmers.
    These are stored in the min heap array.
  */

  /* count kmer hits in the database sequences */
  auto const indexed_count = static_cast<int>(dbindex_getcount_paired());
  assert(searchinfo->target_seqnos_r1 != nullptr);
  assert(searchinfo->target_seqnos_r2 != nullptr);

  /* zero counts */
  std::memset(searchinfo->kmers, 0, indexed_count * sizeof(count_t));

  minheap_clear(searchinfo->m);

  for (auto i = 0U; i < searchinfo->kmersamplecount_r1; i++) {
    auto const kmer = searchinfo->kmersample_r1[i];
    auto *bitmap = dbindex_getbitmap_paired(dbindex_r1_key_paired(kmer));

    if (bitmap != nullptr) {
#ifdef __x86_64__
      if (ssse3_present != 0) {
        increment_counters_from_bitmap_ssse3(searchinfo->kmers, bitmap,
                                             indexed_count);
      } else {
        increment_counters_from_bitmap_sse2(searchinfo->kmers, bitmap,
                                            indexed_count);
      }
#else
      increment_counters_from_bitmap(searchinfo->kmers, bitmap, indexed_count);
#endif
    } else {
      auto *list = dbindex_getmatchlist_paired(dbindex_r1_key_paired(kmer));
      auto const count =
          dbindex_getmatchcount_paired(dbindex_r1_key_paired(kmer));
      for (auto j = 0U; j < count; j++) {
        searchinfo->kmers[list[j]]++;
      }
    }
  }

  for (auto i = 0U; i < searchinfo->kmersamplecount_r2; i++) {
    auto const kmer = searchinfo->kmersample_r2[i];
    auto *bitmap = dbindex_getbitmap_paired(dbindex_r2_key_paired(kmer));

    if (bitmap != nullptr) {
#ifdef __x86_64__
      if (ssse3_present != 0) {
        increment_counters_from_bitmap_ssse3(searchinfo->kmers, bitmap,
                                             indexed_count);
      } else {
        increment_counters_from_bitmap_sse2(searchinfo->kmers, bitmap,
                                            indexed_count);
      }
#else
      increment_counters_from_bitmap(searchinfo->kmers, bitmap, indexed_count);
#endif
    } else {
      auto *list = dbindex_getmatchlist_paired(dbindex_r2_key_paired(kmer));
      auto const count =
          dbindex_getmatchcount_paired(dbindex_r2_key_paired(kmer));
      for (auto j = 0U; j < count; j++) {
        searchinfo->kmers[list[j]]++;
      }
    }
  }

  auto const minmatches =
      std::min(static_cast<unsigned int>(opt_minwordmatches),
               searchinfo->kmersamplecount_r1 + searchinfo->kmersamplecount_r2);

  for (auto i = 0; i < indexed_count; i++) {
    auto const count = searchinfo->kmers[i];
    if (count >= minmatches) {
      auto const seqno = dbindex_getmapping_paired(i);
      auto const target_index = static_cast<std::size_t>(seqno);
      auto const target_seqno_r1 =
          (*searchinfo->target_seqnos_r1)[target_index];
      auto const target_seqno_r2 =
          (*searchinfo->target_seqnos_r2)[target_index];
      elem_t novel;
      novel.count = count;
      novel.seqno = seqno;
      novel.length =
          static_cast<unsigned int>(db_getsequencelen(target_seqno_r1) +
                                    db_getsequencelen(target_seqno_r2));

      minheap_add(searchinfo->m, &novel);
    }
  }

  minheap_sort(searchinfo->m);
}

auto search_acceptable_unaligned_paired(searchinfo_s_paired const &searchinfo,
                                        int const target) -> bool {
  /* consider whether a hit satisfies accepted criteria before alignment */

  // true: needs further consideration
  // false: reject

  auto const target_seqno_r1 =
      (*searchinfo.target_seqnos_r1)[static_cast<std::size_t>(target)];
  auto const target_seqno_r2 =
      (*searchinfo.target_seqnos_r2)[static_cast<std::size_t>(target)];
  auto const *dlabel = db_getheader(target_seqno_r1);
  auto const *dseq_r1 = db_getsequence(target_seqno_r1);
  auto const *dseq_r2 = db_getsequence(target_seqno_r2);
  int64_t const dseqlen_r1 = db_getsequencelen(target_seqno_r1);
  int64_t const dseqlen_r2 = db_getsequencelen(target_seqno_r2);

  struct search_unaligned_numeric_filters_s numeric_filters {};
  numeric_filters.qsize = searchinfo.qsize;
  numeric_filters.tsize = db_getabundance(target_seqno_r1);
  numeric_filters.qlen =
      static_cast<double>(searchinfo.qseqlen_r1 + searchinfo.qseqlen_r2);
  numeric_filters.tlen = static_cast<double>(dseqlen_r1 + dseqlen_r2);

  if (not search_unaligned_numeric_filters_pass(numeric_filters)) {
    return false;
  }

  return (
      /* idprefix */
      ((searchinfo.qseqlen_r1 >= opt_idprefix) and
       (dseqlen_r1 >= opt_idprefix) and
       (seqcmp(searchinfo.qsequence_r1, dseq_r1, opt_idprefix) == 0)) and
      /* idsuffix */
      ((searchinfo.qseqlen_r2 >= opt_idsuffix) and
       (dseqlen_r2 >= opt_idsuffix) and
       (seqcmp(searchinfo.qsequence_r2 + searchinfo.qseqlen_r2 - opt_idsuffix,
               dseq_r2 + dseqlen_r2 - opt_idsuffix, opt_idsuffix) == 0)) and
      /* self */
      ((opt_self == 0) or (std::strcmp(searchinfo.query_head, dlabel) != 0)) and
      /* selfid */
      ((opt_selfid == 0) or (searchinfo.qseqlen_r1 != dseqlen_r1) or
       (searchinfo.qseqlen_r2 != dseqlen_r2) or
       (seqcmp(searchinfo.qsequence_r1, dseq_r1, searchinfo.qseqlen_r1) != 0) or
       (seqcmp(searchinfo.qsequence_r2, dseq_r2, searchinfo.qseqlen_r2) != 0)));
}

auto search_acceptable_aligned_paired(searchinfo_s_paired const &searchinfo,
                                      hit_paired_s *hit) -> bool {
  auto const target_index = static_cast<std::size_t>(hit->target);
  auto const target_seqno_r1 = (*searchinfo.target_seqnos_r1)[target_index];
  auto const target_seqno_r2 = (*searchinfo.target_seqnos_r2)[target_index];

  struct search_aligned_filter_input_s filter_input {};
  filter_input.mismatches = hit->mismatches;
  filter_input.nwgaps = hit->nwgaps;
  filter_input.nwalignmentlength = hit->nwalignmentlength;
  filter_input.internal_alignmentlength = hit->internal_alignmentlength;
  filter_input.internal_gaps = hit->internal_gaps;
  filter_input.internal_indels = hit->internal_indels;
  filter_input.matches = hit->matches;
  filter_input.query_len = searchinfo.qseqlen_r1 + searchinfo.qseqlen_r2;
  filter_input.target_len = static_cast<int>(db_getsequencelen(target_seqno_r1) +
                                             db_getsequencelen(target_seqno_r2));
  filter_input.trim_left_total = hit->r1.trim_q_left + hit->r1.trim_t_left +
                                 hit->r2.trim_q_left + hit->r2.trim_t_left;
  filter_input.trim_right_total = hit->r1.trim_q_right + hit->r1.trim_t_right +
                                  hit->r2.trim_q_right + hit->r2.trim_t_right;

  auto const metrics = search_aligned_compute_identity_metrics(filter_input);
  hit->shortest = metrics.shortest;
  hit->longest = metrics.longest;
  hit->id0 = metrics.id0;
  hit->id1 = metrics.id1;
  hit->id2 = metrics.id2;
  hit->id3 = metrics.id3;
  hit->id4 = metrics.id4;
  hit->id = metrics.id;
  hit->mid = metrics.mid;

  if (not search_aligned_threshold_filters_pass(filter_input, metrics)) {
    /* rejected */
    hit->rejected = true;
    hit->weak = false;
    return false;
  }

  if (opt_cluster_unoise != nullptr) {
    auto target_abundance = static_cast<int64_t>(db_getabundance(target_seqno_r1));
    if (target_abundance < 1) {
      target_abundance = 1;
    }
    auto const skew = 1.0 * searchinfo.qsize / target_abundance;
    auto const beta = 1.0 / std::pow(2, (1.0 * opt_unoise_alpha * hit->mismatches) + 1);

    if (skew <= beta or hit->mismatches == 0) {
      /* accepted */
      hit->accepted = true;
      hit->weak = false;
      return true;
    }
    /* rejected, but weak hit */
    hit->rejected = true;
    hit->weak = true;
    return false;
  }

  if (hit->id >= 100.0 * opt_id) {
    /* accepted */
    hit->accepted = true;
    hit->weak = false;
    return true;
  }
  /* rejected, but weak hit */
  hit->rejected = true;
  hit->weak = true;
  return false;
}

auto align_delayed_paired(searchinfo_s_paired *searchinfo) -> void {
  /* compute global alignment */

  std::array<unsigned int, MAXDELAYED> target_list_r1{{}};
  std::array<unsigned int, MAXDELAYED> target_list_r2{{}};
  std::array<CELL, MAXDELAYED> nwscore_list_r1{{}};
  std::array<CELL, MAXDELAYED> nwscore_list_r2{{}};
  std::array<unsigned short, MAXDELAYED> nwalignmentlength_list_r1{{}};
  std::array<unsigned short, MAXDELAYED> nwalignmentlength_list_r2{{}};
  std::array<unsigned short, MAXDELAYED> nwmatches_list_r1{{}};
  std::array<unsigned short, MAXDELAYED> nwmatches_list_r2{{}};
  std::array<unsigned short, MAXDELAYED> nwmismatches_list_r1{{}};
  std::array<unsigned short, MAXDELAYED> nwmismatches_list_r2{{}};
  std::array<unsigned short, MAXDELAYED> nwgaps_list_r1{{}};
  std::array<unsigned short, MAXDELAYED> nwgaps_list_r2{{}};
  std::array<char *, MAXDELAYED> nwcigar_list_r1{{}};
  std::array<char *, MAXDELAYED> nwcigar_list_r2{{}};

  int target_count = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++) {
    hit_paired_s const *hit = searchinfo->hits + x;
    if (not hit->rejected) {
      auto const target = static_cast<std::size_t>(hit->target);
      target_list_r1[target_count] = (*searchinfo->target_seqnos_r1)[target];
      target_list_r2[target_count] = (*searchinfo->target_seqnos_r2)[target];
      target_count++;
    }
  }

  if (target_count != 0) {
    search16(searchinfo->s_r1, target_count, target_list_r1.data(),
             nwscore_list_r1.data(), nwalignmentlength_list_r1.data(),
             nwmatches_list_r1.data(), nwmismatches_list_r1.data(),
             nwgaps_list_r1.data(), nwcigar_list_r1.data());
    search16(searchinfo->s_r2, target_count, target_list_r2.data(),
             nwscore_list_r2.data(), nwalignmentlength_list_r2.data(),
             nwmatches_list_r2.data(), nwmismatches_list_r2.data(),
             nwgaps_list_r2.data(), nwcigar_list_r2.data());
  }

  int i = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++) {
    /* maxrejects or maxaccepts reached - ignore remaining hits */
    if ((searchinfo->rejects < opt_maxrejects) and
        (searchinfo->accepts < opt_maxaccepts)) {
      hit_paired_s *hit = searchinfo->hits + x;

      if (hit->rejected) {
        searchinfo->rejects++;
      } else {
        auto const target_r1 = target_list_r1[i];
        auto const target_r2 = target_list_r2[i];

        int64_t nwscore_r1 = nwscore_list_r1[i];
        int64_t nwscore_r2 = nwscore_list_r2[i];
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

        if (nwscore_r1 == std::numeric_limits<short>::max()) {
          char *dseq = db_getsequence(target_r1);
          int64_t const dseqlen = db_getsequencelen(target_r1);

          if (nwcigar_list_r1[i] != nullptr) {
            xfree(nwcigar_list_r1[i]);
          }

          nwcigar_r1 = xstrdup(searchinfo->lma_r1->align(
              searchinfo->qsequence_r1, dseq, searchinfo->qseqlen_r1, dseqlen));

          searchinfo->lma_r1->alignstats(nwcigar_r1, searchinfo->qsequence_r1,
                                         dseq, &nwscore_r1,
                                         &nwalignmentlength_r1, &nwmatches_r1,
                                         &nwmismatches_r1, &nwgaps_r1);
        } else {
          nwalignmentlength_r1 = nwalignmentlength_list_r1[i];
          nwmatches_r1 = nwmatches_list_r1[i];
          nwmismatches_r1 = nwmismatches_list_r1[i];
          nwgaps_r1 = nwgaps_list_r1[i];
          nwcigar_r1 = nwcigar_list_r1[i];
        }
        if (nwcigar_r1 == nullptr) {
          nwcigar_r1 = xstrdup("");
        }

        if (nwscore_r2 == std::numeric_limits<short>::max()) {
          char *dseq = db_getsequence(target_r2);
          int64_t const dseqlen = db_getsequencelen(target_r2);

          if (nwcigar_list_r2[i] != nullptr) {
            xfree(nwcigar_list_r2[i]);
          }

          nwcigar_r2 = xstrdup(searchinfo->lma_r2->align(
              searchinfo->qsequence_r2, dseq, searchinfo->qseqlen_r2, dseqlen));

          searchinfo->lma_r2->alignstats(nwcigar_r2, searchinfo->qsequence_r2,
                                         dseq, &nwscore_r2,
                                         &nwalignmentlength_r2, &nwmatches_r2,
                                         &nwmismatches_r2, &nwgaps_r2);
        } else {
          nwalignmentlength_r2 = nwalignmentlength_list_r2[i];
          nwmatches_r2 = nwmatches_list_r2[i];
          nwmismatches_r2 = nwmismatches_list_r2[i];
          nwgaps_r2 = nwgaps_list_r2[i];
          nwcigar_r2 = nwcigar_list_r2[i];
        }
        if (nwcigar_r2 == nullptr) {
          nwcigar_r2 = xstrdup("");
        }

        auto const query_len = searchinfo->qseqlen_r1 + searchinfo->qseqlen_r2;
        auto const target_len = static_cast<int>(db_getsequencelen(target_r1) +
                                                 db_getsequencelen(target_r2));

        hit->aligned = true;
        hit->r1.nwalignment = nwcigar_r1;
        hit->r1.nwscore = static_cast<int>(nwscore_r1);
        hit->r1.nwdiff = static_cast<int>(nwalignmentlength_r1 - nwmatches_r1);
        hit->r1.nwgaps = nwgaps_r1;
        hit->r1.nwindels = static_cast<int>(nwalignmentlength_r1 - nwmatches_r1 -
                                            nwmismatches_r1);
        hit->r1.nwalignmentlength = nwalignmentlength_r1;
        hit->r1.matches = static_cast<int>(nwmatches_r1);
        hit->r1.mismatches = static_cast<int>(nwmismatches_r1);
        hit->r1.nwid = (hit->r1.nwalignmentlength != 0)
                           ? (100.0 * (hit->r1.nwalignmentlength -
                                       hit->r1.nwdiff) /
                              hit->r1.nwalignmentlength)
                           : 0.0;
        hit->r1.shortest = std::min(searchinfo->qseqlen_r1,
                                    static_cast<int>(db_getsequencelen(target_r1)));
        hit->r1.longest = std::max(searchinfo->qseqlen_r1,
                                   static_cast<int>(db_getsequencelen(target_r1)));
        align_trim(&hit->r1);

        hit->r2.nwalignment = nwcigar_r2;
        hit->r2.nwscore = static_cast<int>(nwscore_r2);
        hit->r2.nwdiff = static_cast<int>(nwalignmentlength_r2 - nwmatches_r2);
        hit->r2.nwgaps = nwgaps_r2;
        hit->r2.nwindels = static_cast<int>(nwalignmentlength_r2 - nwmatches_r2 -
                                            nwmismatches_r2);
        hit->r2.nwalignmentlength = nwalignmentlength_r2;
        hit->r2.matches = static_cast<int>(nwmatches_r2);
        hit->r2.mismatches = static_cast<int>(nwmismatches_r2);
        hit->r2.nwid = (hit->r2.nwalignmentlength != 0)
                           ? (100.0 * (hit->r2.nwalignmentlength -
                                       hit->r2.nwdiff) /
                              hit->r2.nwalignmentlength)
                           : 0.0;
        hit->r2.shortest = std::min(searchinfo->qseqlen_r2,
                                    static_cast<int>(db_getsequencelen(target_r2)));
        hit->r2.longest = std::max(searchinfo->qseqlen_r2,
                                   static_cast<int>(db_getsequencelen(target_r2)));
        align_trim(&hit->r2);

        hit->matches = hit->r1.matches + hit->r2.matches;
        hit->mismatches = hit->r1.mismatches + hit->r2.mismatches;
        hit->nwgaps = hit->r1.nwgaps + hit->r2.nwgaps;
        hit->nwindels = hit->r1.nwindels + hit->r2.nwindels;
        hit->nwalignmentlength =
            hit->r1.nwalignmentlength + hit->r2.nwalignmentlength;
        hit->internal_alignmentlength =
            hit->r1.internal_alignmentlength + hit->r2.internal_alignmentlength;
        hit->internal_gaps = hit->r1.internal_gaps + hit->r2.internal_gaps;
        hit->internal_indels = hit->r1.internal_indels + hit->r2.internal_indels;
        hit->shortest = std::min(query_len, target_len);
        hit->longest = std::max(query_len, target_len);
        hit->nwdiff = hit->nwalignmentlength - hit->matches;
        hit->nwscore = static_cast<int>(nwscore_r1 + nwscore_r2);
        if (hit->nwalignmentlength != 0) {
          hit->nwid =
              100.0 * (hit->nwalignmentlength - hit->nwdiff) / hit->nwalignmentlength;
        } else {
          hit->nwid = 0.0;
        }

        hit->mismatches_total = hit->mismatches;
        hit->nwgaps_total = hit->nwgaps;
        hit->nwalignment_cols_total = hit->nwalignmentlength;
        hit->internal_alignment_cols_total = hit->internal_alignmentlength;
        hit->internal_gaps_total = hit->internal_gaps;
        hit->internal_indels_total = hit->internal_indels;

        /* test accept/reject criteria after alignment */
        if (search_acceptable_aligned_paired(*searchinfo, hit)) {
          searchinfo->accepts++;
        } else {
          searchinfo->rejects++;
        }

        ++i;
      }
    }
  }

  /* free ignored alignments */
  while (i < target_count) {
    xfree(nwcigar_list_r1[i]);
    xfree(nwcigar_list_r2[i]);
    ++i;
  }

  searchinfo->finalized = searchinfo->hit_count;
}

auto search_onequery_paired(searchinfo_s_paired *searchinfo, int seqmask)
    -> void {
  searchinfo->hit_count = 0;

  search16_qprep(searchinfo->s_r1, searchinfo->qsequence_r1,
                 searchinfo->qseqlen_r1);
  search16_qprep(searchinfo->s_r2, searchinfo->qsequence_r2,
                 searchinfo->qseqlen_r2);

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

  searchinfo->lma_r1 = new LinearMemoryAligner(scoring);
  searchinfo->lma_r2 = new LinearMemoryAligner(scoring);

  /* extract unique kmer samples from query*/
  unique_count(searchinfo->uh_r1, static_cast<unsigned int>(opt_wordlength),
               searchinfo->qseqlen_r1, searchinfo->qsequence_r1,
               &searchinfo->kmersamplecount_r1, &searchinfo->kmersample_r1,
               seqmask);
  unique_count(searchinfo->uh_r2, static_cast<unsigned int>(opt_wordlength),
               searchinfo->qseqlen_r2, searchinfo->qsequence_r2,
               &searchinfo->kmersamplecount_r2, &searchinfo->kmersample_r2,
               seqmask);

  /* find database sequences with the most kmer hits */
  search_topscores_paired(searchinfo);

  /* analyse targets with the highest number of kmer hits */
  searchinfo->accepts = 0;
  searchinfo->rejects = 0;
  searchinfo->finalized = 0;

  int delayed = 0;

  while ((searchinfo->finalized + delayed <
          opt_maxaccepts + opt_maxrejects - 1) and
         (searchinfo->rejects < opt_maxrejects) and
         (searchinfo->accepts < opt_maxaccepts) and
         (not minheap_isempty(searchinfo->m))) {
    elem_t const e = minheap_poplast(searchinfo->m);

    hit_paired_s *hit = searchinfo->hits + searchinfo->hit_count;
    *hit = hit_paired_s{};

    hit->target = static_cast<int>(e.seqno);
    hit->target_seqno_r1 =
        (*searchinfo->target_seqnos_r1)[static_cast<std::size_t>(e.seqno)];
    hit->target_seqno_r2 =
        (*searchinfo->target_seqnos_r2)[static_cast<std::size_t>(e.seqno)];
    hit->count = e.count;
    hit->strand = searchinfo->strand;
    hit->rejected = false;
    hit->accepted = false;
    hit->aligned = false;
    hit->weak = false;

    /* Test some accept/reject criteria before alignment */
    if (search_acceptable_unaligned_paired(*searchinfo, e.seqno)) {
      ++delayed;
    } else {
      hit->rejected = true;
    }

    searchinfo->hit_count++;

    if (delayed == MAXDELAYED) {
      align_delayed_paired(searchinfo);
      delayed = 0;
    }
  }
  if (delayed > 0) {
    align_delayed_paired(searchinfo);
  }

  delete searchinfo->lma_r1;
  delete searchinfo->lma_r2;
}

auto search_findbest2_byid_paired(searchinfo_s_paired *si_p,
                                  searchinfo_s_paired *si_m) -> hit_paired_s * {
  hit_paired_s *best = nullptr;

  for (int i = 0; i < si_p->hit_count; i++) {
    if ((best == nullptr) or
        (hit_compare_byid_typed_paired(si_p->hits + i, best) < 0)) {
      best = si_p->hits + i;
    }
  }

  if (opt_strand > 1) {
    for (int i = 0; i < si_m->hit_count; i++) {
      if ((best == nullptr) or
          (hit_compare_byid_typed_paired(si_m->hits + i, best) < 0)) {
        best = si_m->hits + i;
      }
    }
  }

  if ((best != nullptr) and not best->accepted) {
    best = nullptr;
  }

  return best;
}

auto search_findbest2_bysize_paired(searchinfo_s_paired *si_p,
                                    searchinfo_s_paired *si_m)
    -> hit_paired_s * {
  hit_paired_s *best = nullptr;

  for (int i = 0; i < si_p->hit_count; i++) {
    if ((best == nullptr) or
        (hit_compare_bysize_typed_paired(si_p->hits + i, best) < 0)) {
      best = si_p->hits + i;
    }
  }

  if (opt_strand > 1) {
    for (int i = 0; i < si_m->hit_count; i++) {
      if ((best == nullptr) or
          (hit_compare_bysize_typed_paired(si_m->hits + i, best) < 0)) {
        best = si_m->hits + i;
      }
    }
  }

  if ((best != nullptr) and not best->accepted) {
    best = nullptr;
  }

  return best;
}

auto search_joinhits_paired(searchinfo_s_paired *si_plus,
                            searchinfo_s_paired *si_minus,
                            std::vector<hit_paired_s> &hits) -> void {
  /* join and sort accepted and weak hits from both strands */
  /* free the remaining alignments */

  auto const counter = count_number_of_hits_to_keep_paired(si_plus) +
                       count_number_of_hits_to_keep_paired(si_minus);

  /* allocate new array of hits */
  hits.reserve(counter);

  copy_over_hits_to_be_kept_paired(hits, si_plus);
  copy_over_hits_to_be_kept_paired(hits, si_minus);

  free_rejected_alignments_paired(si_plus);
  free_rejected_alignments_paired(si_minus);

  /* last, sort the hits */
  std::qsort(hits.data(), counter, sizeof(hit_paired_s),
             hit_compare_byid_paired);
}
