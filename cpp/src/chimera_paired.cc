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

#include "chimera.h"
#include "chimera_paired.h"
#include "align_simd.h"
#include "attributes.h"
#include "dbindex_paired.h"
#include "linmemalign.h"
#include "mask.h"
#include "minheap.h"
#include "unique.h"
#include "utils/cigar.hpp"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/span.hpp"
#include "utils/xpthread.hpp"
#include "vsearch.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <limits>
#include <pthread.h>
#include <string>
#include <vector>

namespace {

constexpr auto parts_paired = 4;
constexpr auto maxparts_paired = 100;
constexpr auto window_paired = 32;
constexpr auto few_paired = 4;
constexpr auto maxcandidates_paired = few_paired * maxparts_paired;
constexpr auto rejects_paired = 16;
constexpr auto chimera_id_paired = 0.55;

static int tophits_paired = 0;
static int records_count_paired = 0;
static pthread_attr_t attr_paired;
static pthread_t *pthread_paired = nullptr;

static pthread_mutex_t mutex_input_paired;
static pthread_mutex_t mutex_output_paired;
static unsigned int seqno_paired = 0;
static uint64_t progress_paired = 0;
static int chimera_count_paired = 0;
static int nonchimera_count_paired = 0;
static int borderline_count_paired = 0;
static int total_count_paired = 0;
static int64_t chimera_abundance_paired = 0;
static int64_t nonchimera_abundance_paired = 0;
static int64_t borderline_abundance_paired = 0;
static int64_t total_abundance_paired = 0;
static int longest_end_paired = 0;

static std::FILE *fp_chimeras_left_paired = nullptr;
static std::FILE *fp_chimeras_right_paired = nullptr;
static std::FILE *fp_nonchimeras_left_paired = nullptr;
static std::FILE *fp_nonchimeras_right_paired = nullptr;
static std::FILE *fp_chimeras_tsv_paired = nullptr;
static std::FILE *fp_nonchimeras_tsv_paired = nullptr;
static std::FILE *fp_uchimealns_paired = nullptr;
static std::FILE *fp_uchimeout_paired = nullptr;
static std::FILE *fp_borderline_paired = nullptr;
static std::FILE *fp_tabbedout_paired = nullptr;

static std::vector<record_paired_s> records_paired;
static std::vector<unsigned int> target_seqnos_r1_paired;
static std::vector<unsigned int> target_seqnos_r2_paired;

struct chimera_info_s_paired {
  int query_alloc_r1 = 0;
  int query_alloc_r2 = 0;
  int head_alloc = 0;

  unsigned int query_no = 0;
  std::vector<char> query_head;
  int query_head_len = 0;
  int64_t query_size = 0;
  std::vector<char> query_seq_r1;
  int query_len_r1 = 0;
  std::vector<char> query_seq_r2;
  int query_len_r2 = 0;

  std::array<searchinfo_s_paired, parts_paired> si{{}};

  std::array<unsigned int, maxcandidates_paired> cand_list{{}};
  std::array<unsigned int, maxcandidates_paired> cand_seqnos_r1{{}};
  std::array<unsigned int, maxcandidates_paired> cand_seqnos_r2{{}};
  int cand_count = 0;

  struct s16info_s *s_r1 = nullptr;
  struct s16info_s *s_r2 = nullptr;

  std::array<CELL, maxcandidates_paired> snwscore_r1{{}};
  std::array<unsigned short, maxcandidates_paired> snwalignmentlength_r1{{}};
  std::array<unsigned short, maxcandidates_paired> snwmatches_r1{{}};
  std::array<unsigned short, maxcandidates_paired> snwmismatches_r1{{}};
  std::array<unsigned short, maxcandidates_paired> snwgaps_r1{{}};
  std::array<int64_t, maxcandidates_paired> nwscore_r1{{}};
  std::array<int64_t, maxcandidates_paired> nwalignmentlength_r1{{}};
  std::array<int64_t, maxcandidates_paired> nwmatches_r1{{}};
  std::array<int64_t, maxcandidates_paired> nwmismatches_r1{{}};
  std::array<int64_t, maxcandidates_paired> nwgaps_r1{{}};
  std::array<char *, maxcandidates_paired> nwcigar_r1{{}};

  std::array<CELL, maxcandidates_paired> snwscore_r2{{}};
  std::array<unsigned short, maxcandidates_paired> snwalignmentlength_r2{{}};
  std::array<unsigned short, maxcandidates_paired> snwmatches_r2{{}};
  std::array<unsigned short, maxcandidates_paired> snwmismatches_r2{{}};
  std::array<unsigned short, maxcandidates_paired> snwgaps_r2{{}};
  std::array<int64_t, maxcandidates_paired> nwscore_r2{{}};
  std::array<int64_t, maxcandidates_paired> nwalignmentlength_r2{{}};
  std::array<int64_t, maxcandidates_paired> nwmatches_r2{{}};
  std::array<int64_t, maxcandidates_paired> nwmismatches_r2{{}};
  std::array<int64_t, maxcandidates_paired> nwgaps_r2{{}};
  std::array<char *, maxcandidates_paired> nwcigar_r2{{}};

  int match_size = 0;
  std::vector<int> match;

  int parents_found = 0;
  std::array<int, maxparents> best_parents{{}};

  double best_h = 0.0;
  bool report_has_parents = false;
  int report_parent_a = -1;
  int report_parent_b = -1;
  int report_best_one = 0;
  int report_best_two = 0;
  int report_delta = 0;
  char const *report_breakpoint_class = "NONE";
};

static chimera_info_s_paired *cia_paired = nullptr;

} // namespace

auto realloc_arrays_paired(chimera_info_s_paired *ci) -> void {
  auto const maxhlen = std::max(ci->query_head_len, 1);
  if (maxhlen > ci->head_alloc) {
    ci->query_head.resize(static_cast<std::size_t>(maxhlen + 1));
    ci->head_alloc = maxhlen;
  }

  auto const maxqlen_r1 = std::max(ci->query_len_r1, 1);
  if (maxqlen_r1 > ci->query_alloc_r1) {
    ci->query_seq_r1.resize(static_cast<std::size_t>(maxqlen_r1 + 1));
    ci->query_alloc_r1 = maxqlen_r1;
  }

  auto const maxqlen_r2 = std::max(ci->query_len_r2, 1);
  if (maxqlen_r2 > ci->query_alloc_r2) {
    ci->query_seq_r2.resize(static_cast<std::size_t>(maxqlen_r2 + 1));
    ci->query_alloc_r2 = maxqlen_r2;
  }

  auto const total_len = std::max(ci->query_len_r1 + ci->query_len_r2, 1);
  if (total_len > ci->match_size) {
    ci->match.resize(static_cast<std::size_t>(maxcandidates_paired * total_len));
    ci->match_size = total_len;
  }

  auto const maxpartlen =
      std::max((total_len + parts_paired - 1) / parts_paired, 100);

  for (auto &search_info : ci->si) {
    if ((maxpartlen + 1) > search_info.seq_alloc) {
      search_info.seq_alloc = maxpartlen + 1;
      search_info.qsequence_r1 = static_cast<char *>(
          xrealloc(search_info.qsequence_r1,
                   static_cast<std::size_t>(search_info.seq_alloc)));
      search_info.qsequence_r2 = static_cast<char *>(
          xrealloc(search_info.qsequence_r2,
                   static_cast<std::size_t>(search_info.seq_alloc)));
    }
  }
}

auto reset_matches_paired(chimera_info_s_paired *ci) -> void {
  std::fill(ci->match.begin(), ci->match.end(), 0);
}

auto query_init_paired(searchinfo_s_paired *search_info) -> void {
  search_info->query_no = 0;
  search_info->strand = 0;
  search_info->qsize = 1;
  search_info->query_head_len = 0;
  search_info->query_head = nullptr;
  search_info->qseqlen_r1 = 0;
  search_info->qseqlen_r2 = 0;
  search_info->seq_alloc = 1;
  search_info->qsequence_r1 = static_cast<char *>(xmalloc(1));
  search_info->qsequence_r2 = static_cast<char *>(xmalloc(1));
  search_info->qsequence_r1[0] = '\0';
  search_info->qsequence_r2[0] = '\0';
  search_info->kmers = static_cast<count_t *>(
      xmalloc((static_cast<std::size_t>(records_count_paired) * sizeof(count_t)) +
              32));
  search_info->hits = static_cast<hit_paired_s *>(
      xmalloc(sizeof(hit_paired_s) * static_cast<std::size_t>(tophits_paired)));
  search_info->hit_count = 0;
  search_info->accepts = 0;
  search_info->rejects = 0;
  search_info->finalized = 0;
  search_info->uh_r1 = unique_init();
  search_info->uh_r2 = unique_init();
  search_info->m = minheap_init(tophits_paired);
  search_info->s_r1 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);
  search_info->s_r2 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);
  search_info->target_seqnos_r1 = &target_seqnos_r1_paired;
  search_info->target_seqnos_r2 = &target_seqnos_r2_paired;
}

auto query_exit_paired(searchinfo_s_paired *search_info) -> void {
  if (search_info->s_r1 != nullptr) {
    search16_exit(search_info->s_r1);
    search_info->s_r1 = nullptr;
  }
  if (search_info->s_r2 != nullptr) {
    search16_exit(search_info->s_r2);
    search_info->s_r2 = nullptr;
  }
  if (search_info->uh_r1 != nullptr) {
    unique_exit(search_info->uh_r1);
    search_info->uh_r1 = nullptr;
  }
  if (search_info->uh_r2 != nullptr) {
    unique_exit(search_info->uh_r2);
    search_info->uh_r2 = nullptr;
  }
  if (search_info->m != nullptr) {
    minheap_exit(search_info->m);
    search_info->m = nullptr;
  }
  if (search_info->qsequence_r1 != nullptr) {
    xfree(search_info->qsequence_r1);
    search_info->qsequence_r1 = nullptr;
  }
  if (search_info->qsequence_r2 != nullptr) {
    xfree(search_info->qsequence_r2);
    search_info->qsequence_r2 = nullptr;
  }
  if (search_info->kmers != nullptr) {
    xfree(search_info->kmers);
    search_info->kmers = nullptr;
  }
  if (search_info->hits != nullptr) {
    xfree(search_info->hits);
    search_info->hits = nullptr;
  }
}

auto partition_query_paired(chimera_info_s_paired *ci) -> void {
  auto const total_len = ci->query_len_r1 + ci->query_len_r2;
  auto rest = total_len;
  auto offset = 0;

  for (auto i = 0; i < parts_paired; ++i) {
    auto const part_len = (rest + (parts_paired - i - 1)) / (parts_paired - i);
    auto const start = offset;
    auto const end = offset + part_len;
    auto const start_r1 = std::min(start, ci->query_len_r1);
    auto const end_r1 = std::min(end, ci->query_len_r1);
    auto const len_r1 = std::max(end_r1 - start_r1, 0);
    auto const start_r2 = std::max(start - ci->query_len_r1, 0);
    auto const end_r2 =
        std::max(std::min(end - ci->query_len_r1, ci->query_len_r2), 0);
    auto const len_r2 = std::max(end_r2 - start_r2, 0);

    auto &search_info = ci->si[static_cast<std::size_t>(i)];
    search_info.query_no = static_cast<int>(ci->query_no);
    search_info.strand = 0;
    search_info.qsize = ci->query_size;
    search_info.query_head_len = ci->query_head_len;
    search_info.query_head = ci->query_head.data();
    search_info.qseqlen_r1 = len_r1;
    search_info.qseqlen_r2 = len_r2;

    if (len_r1 > 0) {
      std::copy(ci->query_seq_r1.data() + start_r1, ci->query_seq_r1.data() + end_r1,
                search_info.qsequence_r1);
    }
    search_info.qsequence_r1[len_r1] = '\0';

    if (len_r2 > 0) {
      std::copy(ci->query_seq_r2.data() + start_r2, ci->query_seq_r2.data() + end_r2,
                search_info.qsequence_r2);
    }
    search_info.qsequence_r2[len_r2] = '\0';

    rest -= part_len;
    offset += part_len;
  }
}

auto find_matches_one_end_paired(chimera_info_s_paired *ci, int cand_index,
                                 int offset, char const *qseq, int qlen,
                                 unsigned int target_seqno, char *cigar_start)
    -> void {
  if ((cigar_start == nullptr) or (qlen <= 0)) {
    return;
  }

  auto const *tseq = db_getsequence(target_seqno);
  auto qpos = 0;
  auto tpos = 0;
  auto const total_len = ci->query_len_r1 + ci->query_len_r2;
  auto const cigar_length = std::strlen(cigar_start);
  auto const cigar_pairs =
      parse_cigar_string(Span<char>{cigar_start, cigar_length});

  for (auto const &a_pair : cigar_pairs) {
    auto const operation = a_pair.first;
    auto const runlength = static_cast<int>(a_pair.second);
    switch (operation) {
    case Operation::match:
      for (auto j = 0; j < runlength; ++j) {
        if ((qpos < qlen) and
            ((map_4bit(qseq[qpos]) & map_4bit(tseq[tpos])) != 0U)) {
          ci->match[(cand_index * total_len) + offset + qpos] = 1;
        }
        ++qpos;
        ++tpos;
      }
      break;

    case Operation::insertion:
      tpos += runlength;
      break;

    case Operation::deletion:
      qpos += runlength;
      break;
    }
  }
}

auto find_matches_paired(chimera_info_s_paired *ci) -> void {
  for (auto i = 0; i < ci->cand_count; ++i) {
    find_matches_one_end_paired(ci, i, 0, ci->query_seq_r1.data(), ci->query_len_r1,
                                ci->cand_seqnos_r1[static_cast<std::size_t>(i)],
                                ci->nwcigar_r1[static_cast<std::size_t>(i)]);
    find_matches_one_end_paired(ci, i, ci->query_len_r1, ci->query_seq_r2.data(),
                                ci->query_len_r2,
                                ci->cand_seqnos_r2[static_cast<std::size_t>(i)],
                                ci->nwcigar_r2[static_cast<std::size_t>(i)]);
  }
}

auto find_best_parents_paired(chimera_info_s_paired *ci) -> int {
  reset_matches_paired(ci);
  find_matches_paired(ci);

  std::array<int, 2> best_parent_cand{{-1, -1}};
  auto const found = select_best_two_parents_from_match_matrix(
      ci->match, ci->cand_count, ci->query_len_r1 + ci->query_len_r2,
      window_paired, best_parent_cand);

  for (auto f = 0U; f < best_parent_cand.size(); ++f) {
    ci->best_parents[f] = best_parent_cand[f];
  }
  ci->parents_found = found ? 2 : 0;
  return static_cast<int>(found);
}

auto eval_parents_paired(chimera_info_s_paired *ci) -> Status {
  auto status = Status::no_alignment;
  ci->report_has_parents = false;

  auto const cand_a = ci->best_parents[0];
  auto const cand_b = ci->best_parents[1];
  if ((cand_a < 0) or (cand_b < 0)) {
    ci->best_h = 0.0;
    return status;
  }

  auto build_end_alignments = [](char const *query, int query_len,
                                 unsigned int target_seqno_a,
                                 unsigned int target_seqno_b, char *cigar_a,
                                 char *cigar_b, std::string &q_aln,
                                 std::string &a_aln,
                                 std::string &b_aln) -> void {
    auto const *target_a = db_getsequence(target_seqno_a);
    auto const *target_b = db_getsequence(target_seqno_b);
    auto const cigar_a_length = static_cast<int>(std::strlen(cigar_a));
    auto const cigar_b_length = static_cast<int>(std::strlen(cigar_b));
    auto const cigar_a_operations = parse_cigar_string(
        Span<char>{cigar_a, static_cast<std::size_t>(cigar_a_length)});
    auto const cigar_b_operations = parse_cigar_string(
        Span<char>{cigar_b, static_cast<std::size_t>(cigar_b_length)});

    std::vector<int> maxi(static_cast<std::size_t>(query_len + 1), 0);

    auto fill_maxi = [&](std::vector<std::pair<Operation, long long>> const &ops)
        -> void {
      auto qpos = 0;
      for (auto const &operation_and_run : ops) {
        auto const operation = operation_and_run.first;
        auto const run = static_cast<int>(operation_and_run.second);
        if (operation == Operation::insertion) {
          maxi[static_cast<std::size_t>(qpos)] =
              std::max(maxi[static_cast<std::size_t>(qpos)], run);
        } else {
          qpos += run;
        }
      }
    };

    fill_maxi(cigar_a_operations);
    fill_maxi(cigar_b_operations);

    q_aln.clear();
    q_aln.reserve(static_cast<std::size_t>(query_len));
    for (auto i = 0; i < query_len; ++i) {
      q_aln.append(static_cast<std::size_t>(maxi[static_cast<std::size_t>(i)]), '-');
      q_aln.push_back(map_uppercase(query[i]));
    }
    q_aln.append(static_cast<std::size_t>(maxi[static_cast<std::size_t>(query_len)]),
                 '-');

    auto build_parent_alignment = [&](char const *target,
                                      std::vector<std::pair<Operation, long long>> const &ops,
                                      std::string &out) -> void {
      out.clear();
      out.reserve(q_aln.size());

      auto qpos = 0;
      auto tpos = 0;
      auto inserted = false;

      for (auto const &operation_and_run : ops) {
        auto const operation = operation_and_run.first;
        auto const run = static_cast<int>(operation_and_run.second);

        if (operation == Operation::insertion) {
          auto const ins = maxi[static_cast<std::size_t>(qpos)];
          for (auto i = 0; i < ins; ++i) {
            if ((i < run) and (target[tpos] != '\0')) {
              out.push_back(map_uppercase(target[tpos]));
              ++tpos;
            } else {
              out.push_back('-');
            }
          }
          inserted = true;
        } else {
          for (auto i = 0; i < run; ++i) {
            if (not inserted) {
              out.append(static_cast<std::size_t>(maxi[static_cast<std::size_t>(qpos)]), '-');
            }
            if (operation == Operation::match) {
              out.push_back(map_uppercase(target[tpos]));
              ++tpos;
            } else {
              out.push_back('-');
            }
            ++qpos;
            inserted = false;
          }
        }
      }

      if (not inserted) {
        out.append(static_cast<std::size_t>(maxi[static_cast<std::size_t>(qpos)]), '-');
      }
    };

    build_parent_alignment(target_a, cigar_a_operations, a_aln);
    build_parent_alignment(target_b, cigar_b_operations, b_aln);
  };

  std::string left_qaln;
  std::string left_aaln;
  std::string left_baln;
  build_end_alignments(ci->query_seq_r1.data(), ci->query_len_r1,
                       ci->cand_seqnos_r1[static_cast<std::size_t>(cand_a)],
                       ci->cand_seqnos_r1[static_cast<std::size_t>(cand_b)],
                       ci->nwcigar_r1[static_cast<std::size_t>(cand_a)],
                       ci->nwcigar_r1[static_cast<std::size_t>(cand_b)],
                       left_qaln, left_aaln, left_baln);

  std::string right_qaln;
  std::string right_aaln;
  std::string right_baln;
  build_end_alignments(ci->query_seq_r2.data(), ci->query_len_r2,
                       ci->cand_seqnos_r2[static_cast<std::size_t>(cand_a)],
                       ci->cand_seqnos_r2[static_cast<std::size_t>(cand_b)],
                       ci->nwcigar_r2[static_cast<std::size_t>(cand_a)],
                       ci->nwcigar_r2[static_cast<std::size_t>(cand_b)],
                       right_qaln, right_aaln, right_baln);

  auto const left_aln_len = static_cast<int>(left_qaln.size());
  std::string q_aln = left_qaln + right_qaln;
  std::string a_aln = left_aaln + right_aaln;
  std::string b_aln = left_baln + right_baln;

  auto const alnlen = static_cast<int>(q_aln.size());
  if (alnlen == 0) {
    ci->best_h = 0.0;
    return status;
  }

  std::string diffs(static_cast<std::size_t>(alnlen), ' ');
  std::string votes(static_cast<std::size_t>(alnlen), ' ');
  std::string model(static_cast<std::size_t>(alnlen), ' ');
  std::vector<bool> ignore(static_cast<std::size_t>(alnlen), false);

  for (auto i = 0; i < alnlen; ++i) {
    auto const qsym = map_4bit(q_aln[static_cast<std::size_t>(i)]);
    auto const asym = map_4bit(a_aln[static_cast<std::size_t>(i)]);
    auto const bsym = map_4bit(b_aln[static_cast<std::size_t>(i)]);

    if ((qsym == 0U) or (asym == 0U) or (bsym == 0U)) {
      ignore[static_cast<std::size_t>(i)] = true;
      if (i > 0) {
        ignore[static_cast<std::size_t>(i - 1)] = true;
      }
      if (i < alnlen - 1) {
        ignore[static_cast<std::size_t>(i + 1)] = true;
      }
    }

    if (is_ambiguous_4bit(qsym) or is_ambiguous_4bit(asym) or
        is_ambiguous_4bit(bsym)) {
      ignore[static_cast<std::size_t>(i)] = true;
    }

    if ((asym != 0U) and (asym != qsym)) {
      a_aln[static_cast<std::size_t>(i)] =
          static_cast<char>(std::tolower(a_aln[static_cast<std::size_t>(i)]));
    }
    if ((bsym != 0U) and (bsym != qsym)) {
      b_aln[static_cast<std::size_t>(i)] =
          static_cast<char>(std::tolower(b_aln[static_cast<std::size_t>(i)]));
    }

    if ((qsym != 0U) and (asym != 0U) and (bsym != 0U)) {
      if (asym == bsym) {
        diffs[static_cast<std::size_t>(i)] = (qsym == asym) ? ' ' : 'N';
      } else if (qsym == asym) {
        diffs[static_cast<std::size_t>(i)] = 'A';
      } else if (qsym == bsym) {
        diffs[static_cast<std::size_t>(i)] = 'B';
      } else {
        diffs[static_cast<std::size_t>(i)] = '?';
      }
    }
  }

  auto sumA = 0;
  auto sumB = 0;
  auto sumN = 0;
  for (auto i = 0; i < alnlen; ++i) {
    if (ignore[static_cast<std::size_t>(i)]) {
      continue;
    }
    auto const diff = diffs[static_cast<std::size_t>(i)];
    if (diff == 'A') {
      ++sumA;
    } else if (diff == 'B') {
      ++sumB;
    } else if (diff != ' ') {
      ++sumN;
    }
  }

  auto left_n = 0;
  auto left_a = 0;
  auto left_y = 0;
  auto right_n = sumA;
  auto right_a = sumN;
  auto right_y = sumB;
  auto best_h = -1.0;
  auto best_i = -1;
  auto best_is_reverse = false;
  auto best_left_y = 0;
  auto best_left_n = 0;
  auto best_left_a = 0;
  auto best_right_y = 0;
  auto best_right_n = 0;
  auto best_right_a = 0;

  std::vector<int> scan_index;
  scan_index.reserve(static_cast<std::size_t>(alnlen));
  for (auto i = 0; i < left_aln_len; ++i) {
    scan_index.push_back(i);
  }
  for (auto i = alnlen - 1; i >= left_aln_len; --i) {
    scan_index.push_back(i);
  }

  std::vector<int> alignment_to_scan_index(static_cast<std::size_t>(alnlen), -1);
  for (auto scan_i = 0; scan_i < static_cast<int>(scan_index.size()); ++scan_i) {
    alignment_to_scan_index[static_cast<std::size_t>(
        scan_index[static_cast<std::size_t>(scan_i)])] = scan_i;
  }

  for (auto scan_i = 0; scan_i < static_cast<int>(scan_index.size()); ++scan_i) {
    auto const i = scan_index[static_cast<std::size_t>(scan_i)];
    if (ignore[static_cast<std::size_t>(i)]) {
      continue;
    }
    auto const diff = diffs[static_cast<std::size_t>(i)];
    if (diff == ' ') {
      continue;
    }

    if (diff == 'A') {
      ++left_y;
      --right_n;
    } else if (diff == 'B') {
      ++left_n;
      --right_y;
    } else {
      ++left_a;
      --right_a;
    }

    if ((left_y > left_n) and (right_y > right_n)) {
      auto const left_h = static_cast<double>(left_y) /
                          ((opt_xn * (static_cast<double>(left_n) + opt_dn)) +
                           static_cast<double>(left_a));
      auto const right_h = static_cast<double>(right_y) /
                           ((opt_xn * (static_cast<double>(right_n) + opt_dn)) +
                            static_cast<double>(right_a));
      auto const h = left_h * right_h;
      if (h > best_h) {
        best_h = h;
        best_i = scan_i;
        best_is_reverse = false;
        best_left_y = left_y;
        best_left_n = left_n;
        best_left_a = left_a;
        best_right_y = right_y;
        best_right_n = right_n;
        best_right_a = right_a;
      }
    } else if ((left_n > left_y) and (right_n > right_y)) {
      auto const left_h = static_cast<double>(left_n) /
                          ((opt_xn * (static_cast<double>(left_y) + opt_dn)) +
                           static_cast<double>(left_a));
      auto const right_h = static_cast<double>(right_n) /
                           ((opt_xn * (static_cast<double>(right_y) + opt_dn)) +
                            static_cast<double>(right_a));
      auto const h = left_h * right_h;
      if (h > best_h) {
        best_h = h;
        best_i = scan_i;
        best_is_reverse = true;
        best_left_y = left_n;
        best_left_n = left_y;
        best_left_a = left_a;
        best_right_y = right_n;
        best_right_n = right_y;
        best_right_a = right_a;
      }
    }
  }

  ci->best_h = best_h > 0.0 ? best_h : 0.0;
  if (best_h < 0.0) {
    return status;
  }

  status = Status::low_score;

  if (best_is_reverse) {
    for (auto i = 0; i < alnlen; ++i) {
      auto &diff = diffs[static_cast<std::size_t>(i)];
      if (diff == 'A') {
        diff = 'B';
      } else if (diff == 'B') {
        diff = 'A';
      }
    }
  }

  for (auto i = 0; i < alnlen; ++i) {
    auto const scan_i = alignment_to_scan_index[static_cast<std::size_t>(i)];
    auto const m = (scan_i <= best_i) ? 'A' : 'B';
    model[static_cast<std::size_t>(i)] = m;

    auto vote = ' ';
    if (not ignore[static_cast<std::size_t>(i)]) {
      auto const diff = diffs[static_cast<std::size_t>(i)];
      if ((diff == 'A') or (diff == 'B')) {
        vote = (diff == m) ? '+' : '!';
      } else if ((diff == 'N') or (diff == '?')) {
        vote = '0';
      }
    }
    votes[static_cast<std::size_t>(i)] = vote;
    if (vote == '!') {
      diffs[static_cast<std::size_t>(i)] =
          static_cast<char>(std::tolower(diffs[static_cast<std::size_t>(i)]));
    }
  }

  for (auto scan_i = best_i + 1; scan_i < static_cast<int>(scan_index.size());
       ++scan_i) {
    auto const i = scan_index[static_cast<std::size_t>(scan_i)];
    auto const diff = diffs[static_cast<std::size_t>(i)];
    if ((diff == ' ') or (diff == 'A')) {
      model[static_cast<std::size_t>(i)] = 'x';
    } else {
      break;
    }
  }

  auto const index_a = best_is_reverse ? 1U : 0U;
  auto const index_b = best_is_reverse ? 0U : 1U;
  auto const &aln_a = best_is_reverse ? b_aln : a_aln;
  auto const &aln_b = best_is_reverse ? a_aln : b_aln;

  auto match_QA = 0;
  auto match_QB = 0;
  auto match_AB = 0;
  auto match_QM = 0;
  auto cols = 0;
  for (auto i = 0; i < alnlen; ++i) {
    if (ignore[static_cast<std::size_t>(i)]) {
      continue;
    }
    ++cols;
    auto const qsym = map_4bit(q_aln[static_cast<std::size_t>(i)]);
    auto const asym = map_4bit(aln_a[static_cast<std::size_t>(i)]);
    auto const bsym = map_4bit(aln_b[static_cast<std::size_t>(i)]);
    auto const scan_i = alignment_to_scan_index[static_cast<std::size_t>(i)];
    auto const msym = (scan_i <= best_i) ? asym : bsym;

    if (qsym == asym) {
      ++match_QA;
    }
    if (qsym == bsym) {
      ++match_QB;
    }
    if (asym == bsym) {
      ++match_AB;
    }
    if (qsym == msym) {
      ++match_QM;
    }
  }

  if (cols <= 0) {
    ci->best_h = 0.0;
    return Status::no_alignment;
  }

  auto const pair_target_a =
      ci->cand_list[static_cast<std::size_t>(ci->best_parents[index_a])];
  auto const pair_target_b =
      ci->cand_list[static_cast<std::size_t>(ci->best_parents[index_b])];
  auto const target_seqno_a_r1 =
      ci->cand_seqnos_r1[static_cast<std::size_t>(ci->best_parents[index_a])];
  auto const target_seqno_a_r2 =
      ci->cand_seqnos_r2[static_cast<std::size_t>(ci->best_parents[index_a])];
  auto const target_seqno_b_r1 =
      ci->cand_seqnos_r1[static_cast<std::size_t>(ci->best_parents[index_b])];
  auto const target_seqno_b_r2 =
      ci->cand_seqnos_r2[static_cast<std::size_t>(ci->best_parents[index_b])];

  auto const QA = 100.0 * match_QA / cols;
  auto const QB = 100.0 * match_QB / cols;
  auto const AB = 100.0 * match_AB / cols;
  auto const QT = std::max(QA, QB);
  auto const QM = 100.0 * match_QM / cols;
  auto const divdiff = QM - QT;
  auto const sumL = best_left_n + best_left_a + best_left_y;
  auto const sumR = best_right_n + best_right_a + best_right_y;

  if ((match_QM == cols) and (QT < 100.0)) {
    status = Status::chimeric;
  }

  ci->report_has_parents = true;
  ci->report_parent_a = static_cast<int>(pair_target_a);
  ci->report_parent_b = static_cast<int>(pair_target_b);
  ci->report_best_one = std::max(match_QA, match_QB);
  ci->report_best_two = match_QM;
  ci->report_delta = ci->report_best_two - ci->report_best_one;
  if (best_i < 0) {
    ci->report_breakpoint_class = "NONE";
  } else if (left_aln_len == 0) {
    ci->report_breakpoint_class = "RIGHT_BREAK";
  } else if (best_i < (left_aln_len - 1)) {
    ci->report_breakpoint_class = "LEFT_BREAK";
  } else if (best_i == (left_aln_len - 1)) {
    ci->report_breakpoint_class = "MIDDLE_BREAK";
  } else {
    ci->report_breakpoint_class = "RIGHT_BREAK";
  }

  xpthread_mutex_lock(&mutex_output_paired);

  if ((fp_uchimealns_paired != nullptr) and (status == Status::chimeric)) {
    auto write_wrapped = [](std::FILE *fp, char const *label,
                            std::string const &seq) -> void {
      auto width = opt_alignwidth;
      if (width <= 0) {
        width = static_cast<int>(seq.size());
      }
      if (width <= 0) {
        std::fprintf(fp, "%-8s\n", label);
        return;
      }

      auto first = true;
      for (std::size_t start = 0; start < seq.size();
           start += static_cast<std::size_t>(width)) {
        auto const chunk =
            std::min(static_cast<std::size_t>(width), seq.size() - start);
        std::fprintf(fp, "%-8s %.*s\n", first ? label : "",
                     static_cast<int>(chunk), seq.data() + start);
        first = false;
      }
    };

    std::fprintf(fp_uchimealns_paired, "\n");
    std::fprintf(fp_uchimealns_paired,
                 "------------------------------------------------------------"
                 "------------\n");
    std::fprintf(fp_uchimealns_paired, "Query   (%5d nt) ",
                 ci->query_len_r1 + ci->query_len_r2);
    header_fprint_strip(fp_uchimealns_paired, ci->query_head.data(),
                        ci->query_head_len, opt_xsize, opt_xee, opt_xlength);
    std::fprintf(fp_uchimealns_paired, "\nParentA (%5" PRIu64 " nt) ",
                 db_getsequencelen(target_seqno_a_r1) +
                     db_getsequencelen(target_seqno_a_r2));
    header_fprint_strip(fp_uchimealns_paired, db_getheader(target_seqno_a_r1),
                        db_getheaderlen(target_seqno_a_r1), opt_xsize,
                        opt_xee, opt_xlength);
    std::fprintf(fp_uchimealns_paired, "\nParentB (%5" PRIu64 " nt) ",
                 db_getsequencelen(target_seqno_b_r1) +
                     db_getsequencelen(target_seqno_b_r2));
    header_fprint_strip(fp_uchimealns_paired, db_getheader(target_seqno_b_r1),
                        db_getheaderlen(target_seqno_b_r1), opt_xsize,
                        opt_xee, opt_xlength);
    std::fprintf(fp_uchimealns_paired, "\nClass %s, Score %.4f\n\n",
                 ci->report_breakpoint_class, ci->best_h);

    auto const left_q = q_aln.substr(0U, static_cast<std::size_t>(left_aln_len));
    auto const left_a = a_aln.substr(0U, static_cast<std::size_t>(left_aln_len));
    auto const left_b = b_aln.substr(0U, static_cast<std::size_t>(left_aln_len));
    auto const left_d = diffs.substr(0U, static_cast<std::size_t>(left_aln_len));
    auto const left_v = votes.substr(0U, static_cast<std::size_t>(left_aln_len));
    auto const left_m = model.substr(0U, static_cast<std::size_t>(left_aln_len));
    auto const right_q = q_aln.substr(static_cast<std::size_t>(left_aln_len));
    auto const right_a = a_aln.substr(static_cast<std::size_t>(left_aln_len));
    auto const right_b = b_aln.substr(static_cast<std::size_t>(left_aln_len));
    auto const right_d = diffs.substr(static_cast<std::size_t>(left_aln_len));
    auto const right_v = votes.substr(static_cast<std::size_t>(left_aln_len));
    auto const right_m = model.substr(static_cast<std::size_t>(left_aln_len));

    write_wrapped(fp_uchimealns_paired, "A_LEFT", left_a);
    write_wrapped(fp_uchimealns_paired, "Q_LEFT", left_q);
    write_wrapped(fp_uchimealns_paired, "B_LEFT", left_b);
    write_wrapped(fp_uchimealns_paired, "DIFF_L", left_d);
    write_wrapped(fp_uchimealns_paired, "VOTE_L", left_v);
    write_wrapped(fp_uchimealns_paired, "MODL_L", left_m);
    std::fprintf(fp_uchimealns_paired, "\n");
    write_wrapped(fp_uchimealns_paired, "A_RIGHT", right_a);
    write_wrapped(fp_uchimealns_paired, "Q_RIGHT", right_q);
    write_wrapped(fp_uchimealns_paired, "B_RIGHT", right_b);
    write_wrapped(fp_uchimealns_paired, "DIFF_R", right_d);
    write_wrapped(fp_uchimealns_paired, "VOTE_R", right_v);
    write_wrapped(fp_uchimealns_paired, "MODL_R", right_m);
    std::fprintf(fp_uchimealns_paired, "\n");

    std::fprintf(fp_uchimealns_paired,
                 "Ids. QA %.1f%%, QB %.1f%%, AB %.1f%%, QModel %.1f%%, Div. %+.1f%%\n",
                 QA, QB, AB, QM, divdiff);
    std::fprintf(fp_uchimealns_paired,
                 "Diffs Left %d: N %d, A %d, Y %d (%.1f%%); Right %d: N %d, A %d, Y %d (%.1f%%), Score %.4f\n",
                 sumL, best_left_n, best_left_a, best_left_y,
                 (sumL > 0) ? (100.0 * best_left_y / sumL) : 0.0, sumR,
                 best_right_n, best_right_a, best_right_y,
                 (sumR > 0) ? (100.0 * best_right_y / sumR) : 0.0, ci->best_h);
  }

  if (fp_uchimeout_paired != nullptr) {
    std::fprintf(fp_uchimeout_paired, "%.4f\t", ci->best_h);
    header_fprint_strip(fp_uchimeout_paired, ci->query_head.data(),
                        ci->query_head_len, opt_xsize, opt_xee, opt_xlength);
    std::fprintf(fp_uchimeout_paired, "\t");
    header_fprint_strip(fp_uchimeout_paired, db_getheader(target_seqno_a_r1),
                        db_getheaderlen(target_seqno_a_r1), opt_xsize,
                        opt_xee, opt_xlength);
    std::fprintf(fp_uchimeout_paired, "\t");
    header_fprint_strip(fp_uchimeout_paired, db_getheader(target_seqno_b_r1),
                        db_getheaderlen(target_seqno_b_r1), opt_xsize,
                        opt_xee, opt_xlength);
    std::fprintf(fp_uchimeout_paired, "\t");

    if (opt_uchimeout5 == 0) {
      if (QA >= QB) {
        header_fprint_strip(fp_uchimeout_paired, db_getheader(target_seqno_a_r1),
                            db_getheaderlen(target_seqno_a_r1), opt_xsize,
                            opt_xee, opt_xlength);
      } else {
        header_fprint_strip(fp_uchimeout_paired, db_getheader(target_seqno_b_r1),
                            db_getheaderlen(target_seqno_b_r1), opt_xsize,
                            opt_xee, opt_xlength);
      }
      std::fprintf(fp_uchimeout_paired, "\t");
    }

    std::fprintf(fp_uchimeout_paired,
                 "%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t"
                 "%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%c\n",
                 QM, QA, QB, AB, QT, best_left_y, best_left_n, best_left_a,
                 best_right_y, best_right_n, best_right_a, divdiff,
                 status == Status::chimeric ? 'Y' : 'N');
  }

  xpthread_mutex_unlock(&mutex_output_paired);
  return status;
}

auto chimera_thread_init_paired(chimera_info_s_paired *ci) -> void {
  for (auto &search_info : ci->si) {
    query_init_paired(&search_info);
  }

  ci->s_r1 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);
  ci->s_r2 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);
}

auto chimera_thread_exit_paired(chimera_info_s_paired *ci) -> void {
  if (ci->s_r1 != nullptr) {
    search16_exit(ci->s_r1);
    ci->s_r1 = nullptr;
  }
  if (ci->s_r2 != nullptr) {
    search16_exit(ci->s_r2);
    ci->s_r2 = nullptr;
  }

  for (auto &search_info : ci->si) {
    query_exit_paired(&search_info);
  }
}

auto chimera_thread_core_paired(chimera_info_s_paired *ci) -> uint64_t {
  chimera_thread_init_paired(ci);

  std::array<hit_paired_s, maxcandidates_paired> allhits_list{};

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

  while (true) {
    auto current_seqno = 0U;

    xpthread_mutex_lock(&mutex_input_paired);
    if (seqno_paired < records_paired.size()) {
      current_seqno = seqno_paired;
      auto const &record = records_paired[static_cast<std::size_t>(current_seqno)];
      ci->query_no = current_seqno;
      ci->query_head_len = static_cast<int>(record.header.size());
      ci->query_len_r1 = static_cast<int>(record.qsequence_r1.size());
      ci->query_len_r2 = static_cast<int>(record.qsequence_r2.size());
      ci->query_size = record.abundance;
      realloc_arrays_paired(ci);
      std::memcpy(ci->query_head.data(), record.header.c_str(),
                  static_cast<std::size_t>(ci->query_head_len));
      ci->query_head[static_cast<std::size_t>(ci->query_head_len)] = '\0';
      std::memcpy(ci->query_seq_r1.data(), record.qsequence_r1.c_str(),
                  static_cast<std::size_t>(ci->query_len_r1));
      ci->query_seq_r1[static_cast<std::size_t>(ci->query_len_r1)] = '\0';
      std::memcpy(ci->query_seq_r2.data(), record.qsequence_r2.c_str(),
                  static_cast<std::size_t>(ci->query_len_r2));
      ci->query_seq_r2[static_cast<std::size_t>(ci->query_len_r2)] = '\0';
    } else {
      xpthread_mutex_unlock(&mutex_input_paired);
      break;
    }
    xpthread_mutex_unlock(&mutex_input_paired);

    auto status = Status::no_parents;
    ci->cand_count = 0;
    ci->parents_found = 0;
    ci->best_h = 0.0;
    ci->report_has_parents = false;
    ci->report_parent_a = -1;
    ci->report_parent_b = -1;
    ci->report_best_one = ci->query_len_r1 + ci->query_len_r2;
    ci->report_best_two = ci->report_best_one;
    ci->report_delta = 0;
    ci->report_breakpoint_class = "NONE";

    partition_query_paired(ci);
    auto allhits_count = 0;
    if ((ci->query_len_r1 + ci->query_len_r2) >= parts_paired) {
      for (auto i = 0; i < parts_paired; ++i) {
        auto &search_info = ci->si[static_cast<std::size_t>(i)];
        unsigned int kmers_r1 = 0;
        unsigned int const *kmersample_r1 = nullptr;
        unique_count(search_info.uh_r1, opt_wordlength, search_info.qseqlen_r1,
                     search_info.qsequence_r1, &kmers_r1, &kmersample_r1,
                     opt_qmask);

        unsigned int kmers_r2 = 0;
        unsigned int const *kmersample_r2 = nullptr;
        unique_count(search_info.uh_r2, opt_wordlength, search_info.qseqlen_r2,
                     search_info.qsequence_r2, &kmers_r2, &kmersample_r2,
                     opt_qmask);

        auto const kmers_total = kmers_r1 + kmers_r2;
        if (kmers_total == 0U) {
          continue;
        }

        auto const indexed_count = static_cast<int>(dbindex_getcount_paired());
        std::memset(search_info.kmers, 0,
                    static_cast<std::size_t>(indexed_count) * sizeof(count_t));
        minheap_clear(search_info.m);

        for (auto j = 0U; j < kmers_r1; ++j) {
          auto const kmer = dbindex_r1_key_paired(kmersample_r1[j]);
          auto *bitmap = dbindex_getbitmap_paired(kmer);
          if (bitmap != nullptr) {
#ifdef __x86_64__
            if (ssse3_present != 0) {
              increment_counters_from_bitmap_ssse3(search_info.kmers, bitmap,
                                                   indexed_count);
            } else {
              increment_counters_from_bitmap_sse2(search_info.kmers, bitmap,
                                                  indexed_count);
            }
#else
            increment_counters_from_bitmap(search_info.kmers, bitmap,
                                           indexed_count);
#endif
          } else {
            auto *list = dbindex_getmatchlist_paired(kmer);
            auto const count = dbindex_getmatchcount_paired(kmer);
            for (auto k = 0U; k < count; ++k) {
              search_info.kmers[list[k]]++;
            }
          }
        }

        for (auto j = 0U; j < kmers_r2; ++j) {
          auto const kmer = dbindex_r2_key_paired(kmersample_r2[j]);
          auto *bitmap = dbindex_getbitmap_paired(kmer);
          if (bitmap != nullptr) {
#ifdef __x86_64__
            if (ssse3_present != 0) {
              increment_counters_from_bitmap_ssse3(search_info.kmers, bitmap,
                                                   indexed_count);
            } else {
              increment_counters_from_bitmap_sse2(search_info.kmers, bitmap,
                                                  indexed_count);
            }
#else
            increment_counters_from_bitmap(search_info.kmers, bitmap,
                                           indexed_count);
#endif
          } else {
            auto *list = dbindex_getmatchlist_paired(kmer);
            auto const count = dbindex_getmatchcount_paired(kmer);
            for (auto k = 0U; k < count; ++k) {
              search_info.kmers[list[k]]++;
            }
          }
        }

        auto const minmatches =
            std::min(static_cast<unsigned int>(opt_minwordmatches), kmers_total);
        for (auto j = 0; j < indexed_count; ++j) {
          auto const count = search_info.kmers[j];
          if (count < minmatches) {
            continue;
          }

          auto const target = dbindex_getmapping_paired(static_cast<unsigned int>(j));
          auto const &parent = records_paired[static_cast<std::size_t>(target)];
          if (parent.abundance <
              static_cast<int64_t>(opt_abskew * static_cast<double>(ci->query_size))) {
            continue;
          }

          elem_t novel;
          novel.count = count;
          novel.seqno = target;
          novel.length = static_cast<unsigned int>(
              db_getsequencelen(target_seqnos_r1_paired[static_cast<std::size_t>(target)]) +
              db_getsequencelen(target_seqnos_r2_paired[static_cast<std::size_t>(target)]));
          minheap_add(search_info.m, &novel);
        }

        minheap_sort(search_info.m);
        while ((not minheap_isempty(search_info.m)) and
               (allhits_count < maxcandidates_paired)) {
          auto const e = minheap_poplast(search_info.m);
          allhits_list[static_cast<std::size_t>(allhits_count)] = hit_paired_s{};
          allhits_list[static_cast<std::size_t>(allhits_count)].target =
              static_cast<int>(e.seqno);
          allhits_list[static_cast<std::size_t>(allhits_count)].accepted = true;
          ++allhits_count;
        }
      }
    }
    for (auto i = 0; i < allhits_count; ++i) {
      auto const target = allhits_list[static_cast<std::size_t>(i)].target;
      auto duplicate = false;
      for (auto k = 0; k < ci->cand_count; ++k) {
        if (ci->cand_list[static_cast<std::size_t>(k)] ==
            static_cast<unsigned int>(target)) {
          duplicate = true;
          break;
        }
      }
      if ((not duplicate) and (ci->cand_count < maxcandidates_paired)) {
        ci->cand_list[static_cast<std::size_t>(ci->cand_count)] =
            static_cast<unsigned int>(target);
        ci->cand_seqnos_r1[static_cast<std::size_t>(ci->cand_count)] =
            target_seqnos_r1_paired[static_cast<std::size_t>(target)];
        ci->cand_seqnos_r2[static_cast<std::size_t>(ci->cand_count)] =
            target_seqnos_r2_paired[static_cast<std::size_t>(target)];
        ++ci->cand_count;
      }

      if (allhits_list[static_cast<std::size_t>(i)].r1.nwalignment != nullptr) {
        xfree(allhits_list[static_cast<std::size_t>(i)].r1.nwalignment);
        allhits_list[static_cast<std::size_t>(i)].r1.nwalignment = nullptr;
      }
      if (allhits_list[static_cast<std::size_t>(i)].r2.nwalignment != nullptr) {
        xfree(allhits_list[static_cast<std::size_t>(i)].r2.nwalignment);
        allhits_list[static_cast<std::size_t>(i)].r2.nwalignment = nullptr;
      }
    }

    if (ci->cand_count > 0) {
      search16_qprep(ci->s_r1, ci->query_seq_r1.data(), ci->query_len_r1);
      search16_qprep(ci->s_r2, ci->query_seq_r2.data(), ci->query_len_r2);

      search16(ci->s_r1, ci->cand_count, ci->cand_seqnos_r1.data(),
               ci->snwscore_r1.data(), ci->snwalignmentlength_r1.data(),
               ci->snwmatches_r1.data(), ci->snwmismatches_r1.data(),
               ci->snwgaps_r1.data(), ci->nwcigar_r1.data());
      search16(ci->s_r2, ci->cand_count, ci->cand_seqnos_r2.data(),
               ci->snwscore_r2.data(), ci->snwalignmentlength_r2.data(),
               ci->snwmatches_r2.data(), ci->snwmismatches_r2.data(),
               ci->snwgaps_r2.data(), ci->nwcigar_r2.data());

      for (auto i = 0; i < ci->cand_count; ++i) {
        auto const target_r1 = ci->cand_seqnos_r1[static_cast<std::size_t>(i)];
        auto const target_r2 = ci->cand_seqnos_r2[static_cast<std::size_t>(i)];

        if (ci->snwscore_r1[static_cast<std::size_t>(i)] ==
            std::numeric_limits<short>::max()) {
          auto *tseq = const_cast<char *>(db_getsequence(target_r1));
          auto const tseqlen = static_cast<int64_t>(db_getsequencelen(target_r1));
          if (ci->nwcigar_r1[static_cast<std::size_t>(i)] != nullptr) {
            xfree(ci->nwcigar_r1[static_cast<std::size_t>(i)]);
          }
          auto *nwcigar = xstrdup(lma_r1.align(ci->query_seq_r1.data(), tseq,
                                               ci->query_len_r1, tseqlen));
          ci->nwcigar_r1[static_cast<std::size_t>(i)] = nwcigar;
          lma_r1.alignstats(nwcigar, ci->query_seq_r1.data(), tseq,
                            &ci->nwscore_r1[static_cast<std::size_t>(i)],
                            &ci->nwalignmentlength_r1[static_cast<std::size_t>(i)],
                            &ci->nwmatches_r1[static_cast<std::size_t>(i)],
                            &ci->nwmismatches_r1[static_cast<std::size_t>(i)],
                            &ci->nwgaps_r1[static_cast<std::size_t>(i)]);
        } else {
          ci->nwscore_r1[static_cast<std::size_t>(i)] =
              ci->snwscore_r1[static_cast<std::size_t>(i)];
          ci->nwalignmentlength_r1[static_cast<std::size_t>(i)] =
              ci->snwalignmentlength_r1[static_cast<std::size_t>(i)];
          ci->nwmatches_r1[static_cast<std::size_t>(i)] =
              ci->snwmatches_r1[static_cast<std::size_t>(i)];
          ci->nwmismatches_r1[static_cast<std::size_t>(i)] =
              ci->snwmismatches_r1[static_cast<std::size_t>(i)];
          ci->nwgaps_r1[static_cast<std::size_t>(i)] =
              ci->snwgaps_r1[static_cast<std::size_t>(i)];
        }

        if (ci->snwscore_r2[static_cast<std::size_t>(i)] ==
            std::numeric_limits<short>::max()) {
          auto *tseq = const_cast<char *>(db_getsequence(target_r2));
          auto const tseqlen = static_cast<int64_t>(db_getsequencelen(target_r2));
          if (ci->nwcigar_r2[static_cast<std::size_t>(i)] != nullptr) {
            xfree(ci->nwcigar_r2[static_cast<std::size_t>(i)]);
          }
          auto *nwcigar = xstrdup(lma_r2.align(ci->query_seq_r2.data(), tseq,
                                               ci->query_len_r2, tseqlen));
          ci->nwcigar_r2[static_cast<std::size_t>(i)] = nwcigar;
          lma_r2.alignstats(nwcigar, ci->query_seq_r2.data(), tseq,
                            &ci->nwscore_r2[static_cast<std::size_t>(i)],
                            &ci->nwalignmentlength_r2[static_cast<std::size_t>(i)],
                            &ci->nwmatches_r2[static_cast<std::size_t>(i)],
                            &ci->nwmismatches_r2[static_cast<std::size_t>(i)],
                            &ci->nwgaps_r2[static_cast<std::size_t>(i)]);
        } else {
          ci->nwscore_r2[static_cast<std::size_t>(i)] =
              ci->snwscore_r2[static_cast<std::size_t>(i)];
          ci->nwalignmentlength_r2[static_cast<std::size_t>(i)] =
              ci->snwalignmentlength_r2[static_cast<std::size_t>(i)];
          ci->nwmatches_r2[static_cast<std::size_t>(i)] =
              ci->snwmatches_r2[static_cast<std::size_t>(i)];
          ci->nwmismatches_r2[static_cast<std::size_t>(i)] =
              ci->snwmismatches_r2[static_cast<std::size_t>(i)];
          ci->nwgaps_r2[static_cast<std::size_t>(i)] =
              ci->snwgaps_r2[static_cast<std::size_t>(i)];
        }
      }

      if (find_best_parents_paired(ci) != 0) {
        status = eval_parents_paired(ci);
      }
    }

    xpthread_mutex_lock(&mutex_output_paired);

    ++total_count_paired;
    total_abundance_paired += ci->query_size;

    auto const query_total_len = ci->query_len_r1 + ci->query_len_r2;
    auto const *parent_a_header = ci->query_head.data();
    auto const *parent_b_header = ci->query_head.data();
    if (ci->report_has_parents) {
      parent_a_header = db_getheader(target_seqnos_r1_paired[static_cast<std::size_t>(ci->report_parent_a)]);
      parent_b_header = db_getheader(target_seqnos_r1_paired[static_cast<std::size_t>(ci->report_parent_b)]);
    }

    if (fp_tabbedout_paired != nullptr) {
      std::fprintf(fp_tabbedout_paired, "%s\t%s\t%s\t%s\t%d\t%d\t%d\n",
                   ci->query_head.data(), parent_a_header, parent_b_header,
                   ci->report_breakpoint_class, ci->report_best_one,
                   ci->report_best_two, ci->report_delta);
    }

    if (status == Status::chimeric) {
      ++chimera_count_paired;
      chimera_abundance_paired += ci->query_size;

      if (fp_chimeras_left_paired != nullptr) {
        fasta_print_general(fp_chimeras_left_paired, nullptr,
                            ci->query_seq_r1.data(), ci->query_len_r1,
                            ci->query_head.data(), ci->query_head_len,
                            static_cast<unsigned int>(ci->query_size),
                            chimera_count_paired, -1.0, -1, -1,
                            opt_fasta_score ? "uchime_denovo" : nullptr,
                            ci->best_h);
      }
      if (fp_chimeras_right_paired != nullptr) {
        fasta_print_general(fp_chimeras_right_paired, nullptr,
                            ci->query_seq_r2.data(), ci->query_len_r2,
                            records_paired[static_cast<std::size_t>(current_seqno)]
                                .header_r2.c_str(),
                            static_cast<int>(records_paired[static_cast<std::size_t>(
                                                 current_seqno)]
                                                 .header_r2.size()),
                            static_cast<unsigned int>(ci->query_size),
                            chimera_count_paired, -1.0, -1, -1,
                            opt_fasta_score ? "uchime_denovo" : nullptr,
                            ci->best_h);
      }
      if (fp_chimeras_tsv_paired != nullptr) {
        std::fprintf(fp_chimeras_tsv_paired, "%s\t%" PRId64 "\t%s\t%s\t%s\n",
                     ci->query_head.data(), ci->query_size, ci->query_head.data(),
                     ci->query_seq_r1.data(), ci->query_seq_r2.data());
      }
    }

    if (status == Status::suspicious) {
      ++borderline_count_paired;
      borderline_abundance_paired += ci->query_size;
    }

    if (status < Status::suspicious) {
      ++nonchimera_count_paired;
      nonchimera_abundance_paired += ci->query_size;

      if ((status < Status::low_score) and (fp_uchimeout_paired != nullptr)) {
        std::fprintf(fp_uchimeout_paired, "%.4f\t", ci->best_h);
        header_fprint_strip(fp_uchimeout_paired, ci->query_head.data(),
                            ci->query_head_len, opt_xsize, opt_xee,
                            opt_xlength);
        if (opt_uchimeout5 != 0) {
          std::fprintf(fp_uchimeout_paired,
                       "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
        } else {
          std::fprintf(fp_uchimeout_paired,
                       "\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
        }
      }

      if (fp_nonchimeras_left_paired != nullptr) {
        fasta_print_general(fp_nonchimeras_left_paired, nullptr,
                            ci->query_seq_r1.data(), ci->query_len_r1,
                            ci->query_head.data(), ci->query_head_len,
                            static_cast<unsigned int>(ci->query_size),
                            nonchimera_count_paired, -1.0, -1, -1,
                            opt_fasta_score ? "uchime_denovo" : nullptr,
                            ci->best_h);
      }
      if (fp_nonchimeras_right_paired != nullptr) {
        fasta_print_general(fp_nonchimeras_right_paired, nullptr,
                            ci->query_seq_r2.data(), ci->query_len_r2,
                            records_paired[static_cast<std::size_t>(current_seqno)]
                                .header_r2.c_str(),
                            static_cast<int>(records_paired[static_cast<std::size_t>(
                                                 current_seqno)]
                                                 .header_r2.size()),
                            static_cast<unsigned int>(ci->query_size),
                            nonchimera_count_paired, -1.0, -1, -1,
                            opt_fasta_score ? "uchime_denovo" : nullptr,
                            ci->best_h);
      }
      if (fp_nonchimeras_tsv_paired != nullptr) {
        std::fprintf(fp_nonchimeras_tsv_paired,
                     "%s\t%" PRId64 "\t%s\t%s\t%s\n", ci->query_head.data(),
                     ci->query_size, ci->query_head.data(), ci->query_seq_r1.data(),
                     ci->query_seq_r2.data());
      }

      dbindex_addsequence_paired(current_seqno, opt_qmask);
    }

    for (auto i = 0; i < ci->cand_count; ++i) {
      if (ci->nwcigar_r1[static_cast<std::size_t>(i)] != nullptr) {
        xfree(ci->nwcigar_r1[static_cast<std::size_t>(i)]);
        ci->nwcigar_r1[static_cast<std::size_t>(i)] = nullptr;
      }
      if (ci->nwcigar_r2[static_cast<std::size_t>(i)] != nullptr) {
        xfree(ci->nwcigar_r2[static_cast<std::size_t>(i)]);
        ci->nwcigar_r2[static_cast<std::size_t>(i)] = nullptr;
      }
    }

    progress_paired += static_cast<uint64_t>(query_total_len);
    progress_update(progress_paired);
    ++seqno_paired;

    xpthread_mutex_unlock(&mutex_output_paired);
  }

  chimera_thread_exit_paired(ci);
  return 0;
}

auto chimera_thread_worker_paired(void *vp) -> void * {
  return reinterpret_cast<void *>(
      chimera_thread_core_paired(cia_paired + reinterpret_cast<int64_t>(vp)));
}

auto chimera_threads_run_paired() -> void {
  xpthread_attr_init(&attr_paired);
  xpthread_attr_setdetachstate(&attr_paired, PTHREAD_CREATE_JOINABLE);

  for (int64_t t = 0; t < opt_threads; ++t) {
    xpthread_create(pthread_paired + t, &attr_paired, chimera_thread_worker_paired,
                    reinterpret_cast<void *>(t));
  }

  for (int64_t t = 0; t < opt_threads; ++t) {
    xpthread_join(pthread_paired[t], nullptr);
  }

  xpthread_attr_destroy(&attr_paired);
}

auto open_chimera_file_paired(std::FILE **output_stream, char const *name)
    -> void {
  if (name != nullptr) {
    *output_stream = fopen_output(name);
    if (*output_stream == nullptr) {
      fatal("Unable to open file %s for writing", name);
    }
  } else {
    *output_stream = nullptr;
  }
}

auto close_chimera_file_paired(std::FILE *output_stream) -> void {
  if (output_stream != nullptr) {
    std::fclose(output_stream);
  }
}

auto uchime3_denovo_paired(struct Parameters const &parameters) -> void {
  open_chimera_file_paired(&fp_chimeras_left_paired, opt_chimeras);
  open_chimera_file_paired(&fp_chimeras_right_paired, opt_chimeras_r2);
  open_chimera_file_paired(&fp_nonchimeras_left_paired, opt_nonchimeras);
  open_chimera_file_paired(&fp_nonchimeras_right_paired, opt_nonchimeras_r2);
  open_chimera_file_paired(&fp_chimeras_tsv_paired, opt_chimeras_tsv);
  open_chimera_file_paired(&fp_nonchimeras_tsv_paired, opt_nonchimeras_tsv);
  open_chimera_file_paired(&fp_uchimealns_paired, opt_uchimealns);
  open_chimera_file_paired(&fp_uchimeout_paired, opt_uchimeout);
  open_chimera_file_paired(&fp_borderline_paired, opt_borderline);
  open_chimera_file_paired(&fp_tabbedout_paired, parameters.opt_tabbedout);

  if (fp_chimeras_tsv_paired != nullptr) {
    std::fprintf(fp_chimeras_tsv_paired,
                 "tav_id\tabundance\theader\tleft_anchor\tright_anchor\n");
  }
  if (fp_nonchimeras_tsv_paired != nullptr) {
    std::fprintf(fp_nonchimeras_tsv_paired,
                 "tav_id\tabundance\theader\tleft_anchor\tright_anchor\n");
  }
  if (fp_tabbedout_paired != nullptr) {
    std::fprintf(fp_tabbedout_paired,
                 "query_tav\tparent_a\tparent_b\tbreakpoint_class\tbest_one_score\tbest_two_score\tdelta\n");
  }

  opt_maxaccepts = few_paired;
  opt_maxrejects = rejects_paired;
  opt_id = chimera_id_paired;
  opt_self = 1;
  opt_selfid = 1;
  opt_threads = 1;
  opt_maxsizeratio = 1.0 / opt_abskew;

  tophits_paired = opt_maxaccepts + opt_maxrejects;
  seqno_paired = 0;
  progress_paired = 0;
  chimera_count_paired = 0;
  nonchimera_count_paired = 0;
  borderline_count_paired = 0;
  total_count_paired = 0;
  chimera_abundance_paired = 0;
  nonchimera_abundance_paired = 0;
  borderline_abundance_paired = 0;
  total_abundance_paired = 0;
  longest_end_paired = 0;
  records_paired.clear();
  target_seqnos_r1_paired.clear();
  target_seqnos_r2_paired.clear();

  auto append_record = [&](std::string const &left_header,
                           std::string const &right_header,
                           int64_t left_abundance,
                           std::string const &left_sequence,
                           std::string const &right_sequence) -> void {
    if (paired_header_key_paired(left_header) !=
        paired_header_key_paired(right_header)) {
      auto const message = std::string{"Paired FASTX headers differ ("} +
                           left_header + " vs " + right_header + ")";
      fatal(message.c_str());
    }
    record_paired_s record;
    record.header = left_header;
    record.header_r2 = right_header;
    record.qsequence_r1 = left_sequence;
    record.qsequence_r2 = right_sequence;

    if ((opt_qmask == MASK_DUST) or
        ((opt_qmask == MASK_SOFT) and (opt_hardmask != 0))) {
      std::vector<char> left_buffer(record.qsequence_r1.begin(),
                                    record.qsequence_r1.end());
      std::vector<char> right_buffer(record.qsequence_r2.begin(),
                                     record.qsequence_r2.end());
      left_buffer.push_back('\0');
      right_buffer.push_back('\0');

      if (opt_qmask == MASK_DUST) {
        dust(left_buffer.data(), static_cast<int>(record.qsequence_r1.size()));
        dust(right_buffer.data(), static_cast<int>(record.qsequence_r2.size()));
      } else {
        hardmask(left_buffer.data(), static_cast<int>(record.qsequence_r1.size()));
        hardmask(right_buffer.data(), static_cast<int>(record.qsequence_r2.size()));
      }

      record.qsequence_r1.assign(left_buffer.data(), record.qsequence_r1.size());
      record.qsequence_r2.assign(right_buffer.data(), record.qsequence_r2.size());
    }

    record.abundance = left_abundance;
    record.first_seen = static_cast<int64_t>(records_paired.size());
    records_paired.push_back(std::move(record));
  };

  if (parameters.opt_interleaved) {
    auto *fastx_h = fastx_open(parameters.opt_uchime3_denovo);
    if (fastx_h == nullptr) {
      fatal("Unrecognized interleaved paired FASTX input: %s",
            parameters.opt_uchime3_denovo);
    }

    while (fastx_next(fastx_h, not opt_notrunclabels,
                      chrmap_no_change_vector.data())) {
      auto const left_header = std::string{fastx_get_header(fastx_h)};
      auto const left_abundance = fastx_get_abundance(fastx_h);
      auto const left_sequence =
          std::string{fastx_get_sequence(fastx_h),
                      static_cast<std::size_t>(fastx_get_sequence_length(fastx_h))};

      if (not fastx_next(fastx_h, not opt_notrunclabels,
                         chrmap_no_change_vector.data())) {
        fatal("Odd number of records in interleaved paired FASTX input %s; expected left/right entries",
              parameters.opt_uchime3_denovo);
      }

      append_record(left_header, std::string{fastx_get_header(fastx_h)},
                    left_abundance, left_sequence,
                    std::string{fastx_get_sequence(fastx_h),
                                static_cast<std::size_t>(
                                    fastx_get_sequence_length(fastx_h))});
    }

    fastx_close(fastx_h);
  } else {
    auto *left_h = fastx_open(parameters.opt_uchime3_denovo);
    if (left_h == nullptr) {
      fatal("Unrecognized left FASTX input: %s", parameters.opt_uchime3_denovo);
    }
    auto *right_h = fastx_open(parameters.opt_reverse);
    if (right_h == nullptr) {
      fastx_close(left_h);
      fatal("Unrecognized right FASTX input: %s", parameters.opt_reverse);
    }

    while (fastx_next(left_h, not opt_notrunclabels,
                      chrmap_no_change_vector.data())) {
      if (not fastx_next(right_h, not opt_notrunclabels,
                         chrmap_no_change_vector.data())) {
        fatal("More forward records than reverse records in paired FASTX input");
      }

      append_record(
          std::string{fastx_get_header(left_h)},
          std::string{fastx_get_header(right_h)},
          fastx_get_abundance(left_h),
          std::string{fastx_get_sequence(left_h),
                      static_cast<std::size_t>(fastx_get_sequence_length(left_h))},
          std::string{fastx_get_sequence(right_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(right_h))});
    }

    if (fastx_next(right_h, not opt_notrunclabels,
                   chrmap_no_change_vector.data())) {
      fatal("More reverse records than forward records in paired FASTX input");
    }

    fastx_close(right_h);
    fastx_close(left_h);
  }

  if (records_paired.empty()) {
    fatal("Input for paired uchime3_denovo is empty");
  }

  std::sort(records_paired.begin(), records_paired.end(),
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

  db_free();
  db_setinfo(false, 0, 0, 0, std::numeric_limits<uint64_t>::max(), 0);

  auto progress_total = uint64_t{0};
  for (auto const &record : records_paired) {
    auto const left_seqno = static_cast<unsigned int>(db_getsequencecount());
    db_add(false, record.header.c_str(), record.qsequence_r1.c_str(), nullptr,
           record.header.size(), record.qsequence_r1.size(), record.abundance);
    auto const right_seqno = static_cast<unsigned int>(db_getsequencecount());
    db_add(false, record.header_r2.c_str(), record.qsequence_r2.c_str(),
           nullptr, record.header_r2.size(), record.qsequence_r2.size(),
           record.abundance);
    target_seqnos_r1_paired.push_back(left_seqno);
    target_seqnos_r2_paired.push_back(right_seqno);
    longest_end_paired = std::max(
        longest_end_paired,
        std::max(static_cast<int>(record.qsequence_r1.size()),
                 static_cast<int>(record.qsequence_r2.size())));
  }
  progress_total = db_getnucleotidecount();

  records_count_paired = static_cast<int>(records_paired.size());
  dbindex_prepare_paired(&records_paired, 1, opt_qmask);

  if (parameters.opt_log != nullptr) {
    std::fprintf(fp_log, "%8.2f  xn\n", opt_xn);
    std::fprintf(fp_log, "%8.2f  dn\n", opt_dn);
    std::fprintf(fp_log, "%8.2f  xa\n", 1.0);
    std::fprintf(fp_log, "%8.2f  id\n", opt_id);
    std::fprintf(fp_log, "%8d  maxp\n\n", 2);
  }

  std::vector<pthread_t> pthread_v(static_cast<std::size_t>(opt_threads));
  pthread_paired = pthread_v.data();
  std::vector<chimera_info_s_paired> cia_v(static_cast<std::size_t>(opt_threads));
  cia_paired = cia_v.data();

  xpthread_mutex_init(&mutex_input_paired, nullptr);
  xpthread_mutex_init(&mutex_output_paired, nullptr);

  progress_init("Detecting chimeras", progress_total);
  chimera_threads_run_paired();
  progress_done();

  xpthread_mutex_destroy(&mutex_output_paired);
  xpthread_mutex_destroy(&mutex_input_paired);

  if (not parameters.opt_quiet) {
    if (total_count_paired > 0) {
      std::fprintf(stderr,
                   "Found %d (%.1f%%) chimeras and %d (%.1f%%) non-chimeras in %d unique paired sequences.\n",
                   chimera_count_paired,
                   100.0 * chimera_count_paired / total_count_paired,
                   nonchimera_count_paired,
                   100.0 * nonchimera_count_paired / total_count_paired,
                   total_count_paired);
    } else {
      std::fprintf(stderr,
                   "Found %d chimeras and %d non-chimeras in %d unique paired sequences.\n",
                   chimera_count_paired, nonchimera_count_paired,
                   total_count_paired);
    }

    if (total_abundance_paired > 0) {
      std::fprintf(stderr,
                   "Taking abundance information into account, this corresponds to\n%" PRId64 " (%.1f%%) chimeras and %" PRId64 " (%.1f%%) non-chimeras in %" PRId64 " total paired sequences.\n",
                   chimera_abundance_paired,
                   100.0 * chimera_abundance_paired / total_abundance_paired,
                   nonchimera_abundance_paired,
                   100.0 * nonchimera_abundance_paired / total_abundance_paired,
                   total_abundance_paired);
    } else {
      std::fprintf(stderr,
                   "Taking abundance information into account, this corresponds to\n%" PRId64 " chimeras and %" PRId64 " non-chimeras in %" PRId64 " total paired sequences.\n",
                   chimera_abundance_paired, nonchimera_abundance_paired,
                   total_abundance_paired);
    }
  }

  dbindex_free_paired();
  db_free();

  close_chimera_file_paired(fp_tabbedout_paired);
  close_chimera_file_paired(fp_borderline_paired);
  close_chimera_file_paired(fp_uchimeout_paired);
  close_chimera_file_paired(fp_uchimealns_paired);
  close_chimera_file_paired(fp_nonchimeras_tsv_paired);
  close_chimera_file_paired(fp_chimeras_tsv_paired);
  close_chimera_file_paired(fp_nonchimeras_right_paired);
  close_chimera_file_paired(fp_nonchimeras_left_paired);
  close_chimera_file_paired(fp_chimeras_right_paired);
  close_chimera_file_paired(fp_chimeras_left_paired);

  show_rusage();
}
