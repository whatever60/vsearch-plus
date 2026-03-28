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

#include "searchcore.h"

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

struct hit_paired_s {
  int target = -1;
  unsigned int target_seqno_r1 = 0;
  unsigned int target_seqno_r2 = 0;
  int strand = 0;
  unsigned int count = 0;
  bool accepted = false;
  bool rejected = false;
  bool aligned = false;
  bool weak = false;

  int nwscore = 0;
  int nwdiff = 0;
  int matches = 0;
  int mismatches = 0;
  int nwgaps = 0;
  int nwindels = 0;
  int nwalignmentlength = 0;
  int internal_alignmentlength = 0;
  int internal_gaps = 0;
  int internal_indels = 0;
  int shortest = 0;
  int longest = 0;
  double nwid = 0.0;
  double id0 = 0.0;
  double id1 = 0.0;
  double id2 = 0.0;
  double id3 = 0.0;
  double id4 = 0.0;
  double id = 0.0;
  double mid = 0.0;

  struct hit r1{};
  struct hit r2{};

  int mismatches_total = 0;
  int nwgaps_total = 0;
  int nwalignment_cols_total = 0;
  int internal_alignment_cols_total = 0;
  int internal_gaps_total = 0;
  int internal_indels_total = 0;
};

struct searchinfo_s_paired {
  int query_no = 0;
  int strand = 0;
  int64_t qsize = 0;
  int query_head_len = 0;
  int query_head_alloc = 0;
  char *query_head = nullptr;

  int qseqlen_r1 = 0;
  int qseqlen_r2 = 0;
  int seq_alloc = 0;
  char *qsequence_r1 = nullptr;
  char *qsequence_r2 = nullptr;

  unsigned int kmersamplecount_r1 = 0;
  unsigned int kmersamplecount_r2 = 0;
  unsigned int const *kmersample_r1 = nullptr;
  unsigned int const *kmersample_r2 = nullptr;

  count_t *kmers = nullptr;

  std::vector<unsigned int> const *target_seqnos_r1 = nullptr;
  std::vector<unsigned int> const *target_seqnos_r2 = nullptr;

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
                                  searchinfo_s_paired *si_m) -> hit_paired_s *;

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
