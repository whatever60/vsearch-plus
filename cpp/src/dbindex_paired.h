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

#ifndef DBINDEX_PAIRED_H
#define DBINDEX_PAIRED_H

#include "bitmap.h"
#include "searchcore_paired.h"

#include <cstdint>
#include <vector>

struct uhandle_s;

extern unsigned int *kmercount_paired;
extern uint64_t *kmerhash_paired;
extern unsigned int *kmerindex_paired;
extern struct bitmap_s **kmerbitmap_paired;
extern unsigned int *dbindex_map_paired;
extern unsigned int dbindex_count_paired;
extern unsigned int kmerhashsize_paired;
extern uint64_t kmerindexsize_paired;
extern uhandle_s *dbindex_uh_paired;

auto dbindex_r1_key_paired(unsigned int kmer) -> unsigned int;
auto dbindex_r2_key_paired(unsigned int kmer) -> unsigned int;

auto dbindex_prepare_paired(std::vector<record_paired_s> const *records,
                            int use_bitmap, int seqmask) -> void;
auto dbindex_addallsequences_paired(int seqmask) -> void;
auto dbindex_addsequence_paired(unsigned int seqno, int seqmask) -> void;
auto dbindex_free_paired() -> void;

auto dbindex_getbitmap_paired(unsigned int kmer) -> unsigned char *;
auto dbindex_getmatchcount_paired(unsigned int kmer) -> unsigned int;
auto dbindex_getmatchlist_paired(unsigned int kmer) -> unsigned int *;
auto dbindex_getmapping_paired(unsigned int index) -> unsigned int;
auto dbindex_getcount_paired() -> unsigned int;

#endif
