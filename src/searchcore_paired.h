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

#ifndef SEARCHCORE_PAIRED_H
#define SEARCHCORE_PAIRED_H

#include "cpu.h"

#include <cstdint>
#include <string>
#include <vector>

struct uhandle_s;
struct minheap_s;
struct s16info_s;
class LinearMemoryAligner;

struct record_paired_s {
  std::string header;
  std::string qsequence_r1;
  std::string qsequence_r2;
  int64_t abundance = 0;
  int64_t first_seen = 0;
};

struct align_stats_paired_s {
  int matches = 0;
  int mismatches = 0;
  int nwgaps = 0;
  int nwindels = 0;
  int nwalignmentlength = 0;
  int internal_alignmentlength = 0;
  int internal_gaps = 0;
  int internal_indels = 0;
  int trim_q_left = 0;
  int trim_q_right = 0;
  int trim_t_left = 0;
  int trim_t_right = 0;
  int trim_aln_left = 0;
  int trim_aln_right = 0;
  double id = 0.0;
  std::string nwalignment;
};

struct hit_paired_s {
  int target = -1;
  int strand = 0;
  unsigned int count = 0;
  bool accepted = false;
  bool rejected = false;
  bool aligned = false;
  bool weak = false;

  align_stats_paired_s r1;
  align_stats_paired_s r2;

  int mismatches_total = 0;
  int nwgaps_total = 0;
  int nwalignment_cols_total = 0;
  int internal_alignment_cols_total = 0;
  int internal_gaps_total = 0;
  int internal_indels_total = 0;
  double id = 0.0;
  double mid = 0.0;
};

struct searchinfo_s_paired {
  int query_no = 0;
  int strand = 0;
  int64_t qsize = 0;

  std::string query_head;

  int qseqlen_r1 = 0;
  int qseqlen_r2 = 0;
  std::vector<char> qsequence_r1_v;
  std::vector<char> qsequence_r2_v;
  char *qsequence_r1 = nullptr;
  char *qsequence_r2 = nullptr;

  unsigned int kmersamplecount_r1 = 0;
  unsigned int kmersamplecount_r2 = 0;
  unsigned int const *kmersample_r1 = nullptr;
  unsigned int const *kmersample_r2 = nullptr;

  std::vector<count_t> kmers_v;
  count_t *kmers = nullptr;

  std::vector<unsigned int> const *target_seqnos_r1 = nullptr;
  std::vector<unsigned int> const *target_seqnos_r2 = nullptr;

  std::vector<hit_paired_s> hits_v;
  hit_paired_s *hits = nullptr;
  int hit_count = 0;

  int accepts = 0;
  int rejects = 0;
  int finalized = 0;

  struct uhandle_s *uh_r1 = nullptr;
  struct uhandle_s *uh_r2 = nullptr;
  struct s16info_s *s_r1 = nullptr;
  struct s16info_s *s_r2 = nullptr;
  LinearMemoryAligner *lma_r1 = nullptr;
  LinearMemoryAligner *lma_r2 = nullptr;
  struct minheap_s *m = nullptr;
};

auto search_topscores_paired(searchinfo_s_paired *searchinfo) -> void;

auto search_onequery_paired(searchinfo_s_paired *searchinfo, int seqmask)
    -> void;

auto search_findbest2_byid_paired(searchinfo_s_paired *si_p,
                                  searchinfo_s_paired *si_m)
    -> hit_paired_s *;

auto search_findbest2_bysize_paired(searchinfo_s_paired *si_p,
                                    searchinfo_s_paired *si_m)
    -> hit_paired_s *;

auto search_acceptable_unaligned_paired(searchinfo_s_paired const &searchinfo,
                                        int target) -> bool;

auto search_acceptable_aligned_paired(searchinfo_s_paired const &searchinfo,
                                      hit_paired_s *hit) -> bool;

auto search_joinhits_paired(searchinfo_s_paired *si_plus,
                            searchinfo_s_paired *si_minus,
                            std::vector<hit_paired_s> &hits) -> void;

auto search_enough_kmers_paired(searchinfo_s_paired const &searchinfo,
                                unsigned int count) -> bool;

#endif
