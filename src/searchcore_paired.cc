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

#include "vsearch.h"
#include "searchcore_paired.h"

#include "minheap.h"
#include "unique.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>

namespace {

auto count_number_of_hits_to_keep_paired(searchinfo_s_paired const * search_info) -> std::size_t
{
  if (search_info == nullptr)
    {
      return std::size_t{0};
    }

  auto kept = std::size_t{0};
  for (auto const & hit : search_info->hits)
    {
      if (hit.accepted or hit.weak)
        {
          ++kept;
        }
    }
  return kept;
}


auto copy_over_hits_to_be_kept_paired(std::vector<hit_paired_s> & hits,
                                      searchinfo_s_paired const * search_info) -> void
{
  if (search_info == nullptr)
    {
      return;
    }

  for (auto const & hit : search_info->hits)
    {
      if (hit.accepted or hit.weak)
        {
          hits.emplace_back(hit);
        }
    }
}


auto free_rejected_alignments_paired(searchinfo_s_paired const * search_info) -> void
{
  (void) search_info;
  /*
    stock searchcore frees rejected CIGAR allocations here.
    paired hit storage currently keeps alignment strings in std::string,
    so there is no xfree-equivalent needed in this helper.
  */
}

}  // namespace

auto hit_compare_byid_paired(void const * lhs_ptr, void const * rhs_ptr) -> int
{
  auto const * lhs = static_cast<hit_paired_s const *>(lhs_ptr);
  auto const * rhs = static_cast<hit_paired_s const *>(rhs_ptr);

  if (lhs->rejected < rhs->rejected)
    {
      return -1;
    }
  if (lhs->rejected > rhs->rejected)
    {
      return +1;
    }

  if (lhs->id > rhs->id)
    {
      return -1;
    }
  if (lhs->id < rhs->id)
    {
      return +1;
    }

  if (lhs->target < rhs->target)
    {
      return -1;
    }
  if (lhs->target > rhs->target)
    {
      return +1;
    }

  return 0;
}



auto search_topscores_paired(searchinfo_s_paired * searchinfo,
                             int wordlength) -> void
{
  (void) wordlength;
  auto const indexed_count = static_cast<int>(searchinfo->target_count);

  std::memset(searchinfo->kmers, 0, static_cast<std::size_t>(indexed_count) * sizeof(unsigned short));

  minheap_clear(searchinfo->m);

  for (auto i = 0U; i < searchinfo->kmersamplecount_r1; ++i)
    {
      auto const kmer = searchinfo->kmersample_r1[i];
      auto const it = searchinfo->target_kmer_index->r1_postings.find(kmer);
      if (it == searchinfo->target_kmer_index->r1_postings.end())
        {
          continue;
        }
      for (auto const target_index : it->second)
        {
          ++searchinfo->kmers[static_cast<std::size_t>(target_index)];
        }
    }

  for (auto i = 0U; i < searchinfo->kmersamplecount_r2; ++i)
    {
      auto const kmer = searchinfo->kmersample_r2[i];
      auto const it = searchinfo->target_kmer_index->r2_postings.find(kmer);
      if (it == searchinfo->target_kmer_index->r2_postings.end())
        {
          continue;
        }
      for (auto const target_index : it->second)
        {
          ++searchinfo->kmers[static_cast<std::size_t>(target_index)];
        }
    }

  auto const kmersamplecount_total = searchinfo->kmersamplecount_r1 + searchinfo->kmersamplecount_r2;
  auto const minmatches = std::min(static_cast<unsigned int>(opt_minwordmatches), kmersamplecount_total);

  for (auto i = 0; i < indexed_count; ++i)
    {
      auto const count = searchinfo->kmers[static_cast<std::size_t>(i)];
      if (count >= minmatches)
        {
          elem_t novel;
          novel.count = count;
          novel.seqno = static_cast<unsigned int>(i);
          novel.length = 0U;
          if ((searchinfo->target_lengths != nullptr) and (static_cast<std::size_t>(i) < searchinfo->target_lengths->size()))
            {
              novel.length = (*searchinfo->target_lengths)[static_cast<std::size_t>(i)];
            }

          minheap_add(searchinfo->m, & novel);
        }
    }

  minheap_sort(searchinfo->m);
}


auto search_onequery_paired(searchinfo_s_paired * searchinfo,
                            int seqmask) -> void
{
  searchinfo->hit_count = 0;
  searchinfo->hits.clear();

  unique_count(searchinfo->uh,
               static_cast<unsigned int>(opt_wordlength),
               searchinfo->qseqlen_r1,
               searchinfo->qsequence_r1,
               & searchinfo->kmersamplecount_r1,
               & searchinfo->kmersample_r1,
               seqmask);

  unique_count(searchinfo->uh,
               static_cast<unsigned int>(opt_wordlength),
               searchinfo->qseqlen_r2,
               searchinfo->qsequence_r2,
               & searchinfo->kmersamplecount_r2,
               & searchinfo->kmersample_r2,
               seqmask);

  search_topscores_paired(searchinfo, static_cast<int>(opt_wordlength));

  searchinfo->accepts = 0;
  searchinfo->rejects = 0;
  searchinfo->finalized = 0;

  auto delayed = 0;

  while ((searchinfo->finalized + delayed < opt_maxaccepts + opt_maxrejects - 1) and
         (searchinfo->rejects < opt_maxrejects) and
         (searchinfo->accepts < opt_maxaccepts) and
         (not minheap_isempty(searchinfo->m)))
    {
      elem_t const e = minheap_poplast(searchinfo->m);

      hit_paired_s hit;
      hit.target = static_cast<int>(e.seqno);
      hit.count = e.count;
      hit.strand = searchinfo->strand;
      hit.rejected = false;
      hit.accepted = false;
      hit.aligned = false;
      hit.weak = false;

      ++delayed;
      searchinfo->hits.push_back(hit);
      ++searchinfo->hit_count;

      if (delayed == static_cast<int>(MAXDELAYED))
        {
          searchinfo->finalized = searchinfo->hit_count;
          delayed = 0;
        }
    }

  if (delayed > 0)
    {
      searchinfo->finalized = searchinfo->hit_count;
    }
}


auto search_joinhits_paired(searchinfo_s_paired * si_plus,
                            searchinfo_s_paired * si_minus,
                            std::vector<hit_paired_s> & hits) -> void
{
  auto const counter = count_number_of_hits_to_keep_paired(si_plus) +
                       count_number_of_hits_to_keep_paired(si_minus);

  hits.reserve(counter);

  copy_over_hits_to_be_kept_paired(hits, si_plus);
  copy_over_hits_to_be_kept_paired(hits, si_minus);

  free_rejected_alignments_paired(si_plus);
  free_rejected_alignments_paired(si_minus);

  if (counter > 0U)
    {
      std::qsort(hits.data(), counter, sizeof(hit_paired_s), hit_compare_byid_paired);
    }
}
