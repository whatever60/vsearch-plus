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

#include "dbindex_paired.h"

#include "bitmap.h"
#include "unique.h"
#include "vsearch.h"

#include <cassert>
#include <cstring>
#include <iterator>

unsigned int *kmercount_paired;
uint64_t *kmerhash_paired;
unsigned int *kmerindex_paired;
struct bitmap_s **kmerbitmap_paired;
unsigned int *dbindex_map_paired;
unsigned int kmerhashsize_paired;
uint64_t kmerindexsize_paired;
unsigned int dbindex_count_paired;
uhandle_s *dbindex_uh_paired;

namespace {

constexpr unsigned int bitmap_threshold_paired = 8;

static unsigned int bitmap_mincount_paired;
static std::vector<record_paired_s> const *records_for_dbindex_paired = nullptr;

} // namespace

auto dbindex_r1_key_paired(unsigned int const kmer) -> unsigned int {
  return kmer;
}

auto dbindex_r2_key_paired(unsigned int const kmer) -> unsigned int {
  return kmer + (1U << (2 * opt_wordlength));
}

auto dbindex_getbitmap_paired(unsigned int const kmer) -> unsigned char * {
  auto *a_bitmap_s = *std::next(kmerbitmap_paired, kmer);
  if (a_bitmap_s != nullptr) {
    return a_bitmap_s->bitmap;
  }
  return nullptr;
}

auto dbindex_getmatchcount_paired(unsigned int const kmer) -> unsigned int {
  return *std::next(kmercount_paired, kmer);
}

auto dbindex_getmatchlist_paired(unsigned int const kmer) -> unsigned int * {
  return std::next(kmerindex_paired, *std::next(kmerhash_paired, kmer));
}

auto dbindex_getmapping_paired(unsigned int const index) -> unsigned int {
  return *std::next(dbindex_map_paired, index);
}

auto dbindex_getcount_paired() -> unsigned int { return dbindex_count_paired; }

auto dbindex_addsequence_paired(unsigned int const seqno, int const seqmask)
    -> void {
#if 0
  std::printf("Adding paired seqno %d as index element no %d\n", seqno,
              dbindex_count_paired);
#endif

  assert(records_for_dbindex_paired != nullptr);
  auto const &record = (*records_for_dbindex_paired)[seqno];

  unsigned int uniquecount_r1 = 0;
  unsigned int const *uniquelist_r1 = nullptr;
  unique_count(dbindex_uh_paired, opt_wordlength, record.qsequence_r1.size(),
               record.qsequence_r1.c_str(), &uniquecount_r1, &uniquelist_r1,
               seqmask);

  unsigned int uniquecount_r2 = 0;
  unsigned int const *uniquelist_r2 = nullptr;
  unique_count(dbindex_uh_paired, opt_wordlength, record.qsequence_r2.size(),
               record.qsequence_r2.c_str(), &uniquecount_r2, &uniquelist_r2,
               seqmask);

  dbindex_map_paired[dbindex_count_paired] = seqno;

  for (auto i = 0U; i < uniquecount_r1; i++) {
    auto const kmer = dbindex_r1_key_paired(uniquelist_r1[i]);
    if (kmerbitmap_paired[kmer] != nullptr) {
      ++kmercount_paired[kmer];
      bitmap_set(kmerbitmap_paired[kmer], dbindex_count_paired);
    } else {
      kmerindex_paired[kmerhash_paired[kmer] + kmercount_paired[kmer]] =
          dbindex_count_paired;
      ++kmercount_paired[kmer];
    }
  }

  for (auto i = 0U; i < uniquecount_r2; i++) {
    auto const kmer = dbindex_r2_key_paired(uniquelist_r2[i]);
    if (kmerbitmap_paired[kmer] != nullptr) {
      ++kmercount_paired[kmer];
      bitmap_set(kmerbitmap_paired[kmer], dbindex_count_paired);
    } else {
      kmerindex_paired[kmerhash_paired[kmer] + kmercount_paired[kmer]] =
          dbindex_count_paired;
      ++kmercount_paired[kmer];
    }
  }

  ++dbindex_count_paired;
}

auto dbindex_addallsequences_paired(int const seqmask) -> void {
  assert(records_for_dbindex_paired != nullptr);
  auto const seqcount =
      static_cast<unsigned int>(records_for_dbindex_paired->size());
  progress_init("Creating k-mer index", seqcount);
  for (auto seqno = 0U; seqno < seqcount; seqno++) {
    dbindex_addsequence_paired(seqno, seqmask);
    progress_update(seqno);
  }
  progress_done();
}

auto dbindex_prepare_paired(std::vector<record_paired_s> const *records,
                            int const use_bitmap, int const seqmask) -> void {
  records_for_dbindex_paired = records;
  dbindex_uh_paired = unique_init();

  assert(records_for_dbindex_paired != nullptr);
  auto const seqcount =
      static_cast<unsigned int>(records_for_dbindex_paired->size());
  auto const single_end_kmerhashsize = 1U << (2 * opt_wordlength);
  kmerhashsize_paired = 2U * single_end_kmerhashsize;

  kmercount_paired =
      (unsigned int *)xmalloc(kmerhashsize_paired * sizeof(unsigned int));
  std::memset(kmercount_paired, 0, kmerhashsize_paired * sizeof(unsigned int));

  progress_init("Counting k-mers", seqcount);
  for (auto seqno = 0U; seqno < seqcount; seqno++) {
    auto const &record = (*records_for_dbindex_paired)[seqno];

    unsigned int uniquecount_r1 = 0;
    unsigned int const *uniquelist_r1 = nullptr;
    unique_count(dbindex_uh_paired, opt_wordlength, record.qsequence_r1.size(),
                 record.qsequence_r1.c_str(), &uniquecount_r1, &uniquelist_r1,
                 seqmask);
    for (auto i = 0U; i < uniquecount_r1; i++) {
      ++kmercount_paired[dbindex_r1_key_paired(uniquelist_r1[i])];
    }

    unsigned int uniquecount_r2 = 0;
    unsigned int const *uniquelist_r2 = nullptr;
    unique_count(dbindex_uh_paired, opt_wordlength, record.qsequence_r2.size(),
                 record.qsequence_r2.c_str(), &uniquecount_r2, &uniquelist_r2,
                 seqmask);
    for (auto i = 0U; i < uniquecount_r2; i++) {
      ++kmercount_paired[dbindex_r2_key_paired(uniquelist_r2[i])];
    }

    progress_update(seqno);
  }
  progress_done();

  if (use_bitmap != 0) {
    bitmap_mincount_paired = seqcount / bitmap_threshold_paired;
  } else {
    bitmap_mincount_paired = seqcount + 1;
  }

  kmerbitmap_paired = (struct bitmap_s **)xmalloc(kmerhashsize_paired *
                                                  sizeof(struct bitmap_s *));
  std::memset(kmerbitmap_paired, 0,
              kmerhashsize_paired * sizeof(struct bitmap_s *));

  kmerhash_paired =
      (uint64_t *)xmalloc((kmerhashsize_paired + 1) * sizeof(uint64_t));
  uint64_t sum = 0;
  for (auto i = 0U; i < kmerhashsize_paired; i++) {
    kmerhash_paired[i] = sum;
    if (kmercount_paired[i] >= bitmap_mincount_paired) {
      kmerbitmap_paired[i] = bitmap_init(seqcount + 127);
      bitmap_reset_all(kmerbitmap_paired[i]);
    } else {
      sum += kmercount_paired[i];
    }
  }
  kmerindexsize_paired = sum;
  kmerhash_paired[kmerhashsize_paired] = sum;

  std::memset(kmercount_paired, 0, kmerhashsize_paired * sizeof(unsigned int));

  kmerindex_paired =
      (unsigned int *)xmalloc(kmerindexsize_paired * sizeof(unsigned int));
  dbindex_map_paired = (unsigned int *)xmalloc(seqcount * sizeof(unsigned int));

  dbindex_count_paired = 0;

  show_rusage();
}

auto dbindex_free_paired() -> void {
  xfree(kmerhash_paired);
  xfree(kmerindex_paired);
  xfree(kmercount_paired);
  xfree(dbindex_map_paired);

  for (auto kmer = 0U; kmer < kmerhashsize_paired; kmer++) {
    if (kmerbitmap_paired[kmer] != nullptr) {
      bitmap_free(kmerbitmap_paired[kmer]);
    }
  }
  xfree(kmerbitmap_paired);
  unique_exit(dbindex_uh_paired);
  records_for_dbindex_paired = nullptr;
}
