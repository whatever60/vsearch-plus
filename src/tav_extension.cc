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
#include "attributes.h"
#include "minheap.h"
#include "otutable.h"
#include "tav_extension.h"
#include "mask.h"
#include "showalign.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include <algorithm>
#include <cmath>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct TavRecord
{
  std::string tav_id;
  std::string header;
  std::string left;
  std::string right;
  int64_t abundance = 0;
  int64_t first_seen = 0;
};

struct TavAssignment
{
  int index = -1;
  int d_left = 0;
  int d_right = 0;
  int d_total = 0;
  double id_left = 0.0;
  double id_right = 0.0;
  double id_total = 0.0;
};

struct TavAlignStats
{
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

struct TavHit
{
  int centroid_index = -1;
  unsigned int kmer_score = 0;
  bool rejected = false;
  bool accepted = false;
  bool weak = false;
  bool aligned = false;
  TavAlignStats left;
  TavAlignStats right;
  int mismatches_total = 0;
  int nwgaps_total = 0;
  int nwalignment_cols_total = 0;
  int internal_alignment_cols_total = 0;
  int internal_gaps_total = 0;
  int internal_indels_total = 0;
  double id = 0.0;
  double mid = 0.0;
};

struct TavKmerIndex
{
  std::unordered_map<unsigned int, std::vector<int>> left_postings;
  std::unordered_map<unsigned int, std::vector<int>> right_postings;
};

struct TavLoadResult
{
  std::vector<TavRecord> records;
  int64_t filtered_short = 0;
  int64_t filtered_long = 0;
};

auto align_one_end_stock_style(std::string const & query,
                               std::string const & target,
                               LinearMemoryAligner & lma) -> TavAlignStats
{
  TavAlignStats stats;
  std::vector<char> query_buf(query.begin(), query.end());
  query_buf.push_back('\0');
  std::vector<char> target_buf(target.begin(), target.end());
  target_buf.push_back('\0');

  auto * cigar = xstrdup(lma.align(query_buf.data(),
                                   target_buf.data(),
                                   static_cast<int64_t>(query.size()),
                                   static_cast<int64_t>(target.size())));

  int64_t nwscore = 0;
  int64_t nwalignmentlength = 0;
  int64_t nwmatches = 0;
  int64_t nwmismatches = 0;
  int64_t nwgaps = 0;
  lma.alignstats(cigar,
                 query_buf.data(),
                 target_buf.data(),
                 & nwscore,
                 & nwalignmentlength,
                 & nwmatches,
                 & nwmismatches,
                 & nwgaps);

  struct hit h {};
  h.nwalignment = cigar;
  h.nwscore = static_cast<int>(nwscore);
  h.nwalignmentlength = static_cast<int>(nwalignmentlength);
  h.nwgaps = static_cast<int>(nwgaps);
  h.nwdiff = h.nwalignmentlength - static_cast<int>(nwmatches);
  h.nwindels = h.nwalignmentlength - static_cast<int>(nwmatches) - static_cast<int>(nwmismatches);
  h.matches = h.nwalignmentlength - h.nwdiff;
  h.mismatches = h.nwdiff - h.nwindels;
  h.shortest = std::min(static_cast<int>(query.size()), static_cast<int>(target.size()));
  h.longest = std::max(static_cast<int>(query.size()), static_cast<int>(target.size()));

  align_trim(&h);

  stats.matches = h.matches;
  stats.mismatches = h.mismatches;
  stats.nwgaps = h.nwgaps;
  stats.nwindels = h.nwindels;
  stats.nwalignmentlength = h.nwalignmentlength;
  stats.internal_alignmentlength = h.internal_alignmentlength;
  stats.internal_gaps = h.internal_gaps;
  stats.internal_indels = h.internal_indels;
  stats.trim_q_left = h.trim_q_left;
  stats.trim_q_right = h.trim_q_right;
  stats.trim_t_left = h.trim_t_left;
  stats.trim_t_right = h.trim_t_right;
  stats.trim_aln_left = h.trim_aln_left;
  stats.trim_aln_right = h.trim_aln_right;
  stats.id = h.id;
  stats.nwalignment = cigar;

  xfree(cigar);

  return stats;
}


auto paired_unaligned_filters_pass(TavRecord const & query,
                                   TavRecord const & target) -> bool
{
  auto const qsize = query.abundance;
  auto const tsize = target.abundance;

  auto const qlen = static_cast<double>(query.left.size() + query.right.size());
  auto const tlen = static_cast<double>(target.left.size() + target.right.size());

  if (qsize > opt_maxqsize)
    {
      return false;
    }
  if (tsize < opt_mintsize)
    {
      return false;
    }
  if (static_cast<double>(qsize) < (opt_minsizeratio * static_cast<double>(tsize)))
    {
      return false;
    }
  if (static_cast<double>(qsize) > (opt_maxsizeratio * static_cast<double>(tsize)))
    {
      return false;
    }
  if (qlen < (opt_minqt * tlen))
    {
      return false;
    }
  if (qlen > (opt_maxqt * tlen))
    {
      return false;
    }

  auto const shorter = std::min(qlen, tlen);
  auto const longer = std::max(qlen, tlen);
  if (shorter < (opt_minsl * longer))
    {
      return false;
    }
  if (shorter > (opt_maxsl * longer))
    {
      return false;
    }

  if (opt_idprefix > 0)
    {
      auto const p = static_cast<std::size_t>(opt_idprefix);
      if ((query.left.size() < p) or (target.left.size() < p) or
          (query.right.size() < p) or (target.right.size() < p))
        {
          return false;
        }
      if ((query.left.compare(0, p, target.left, 0, p) != 0) or
          (query.right.compare(0, p, target.right, 0, p) != 0))
        {
          return false;
        }
    }

  if (opt_idsuffix > 0)
    {
      auto const s = static_cast<std::size_t>(opt_idsuffix);
      if ((query.left.size() < s) or (target.left.size() < s) or
          (query.right.size() < s) or (target.right.size() < s))
        {
          return false;
        }
      if ((query.left.compare(query.left.size() - s, s, target.left, target.left.size() - s, s) != 0) or
          (query.right.compare(query.right.size() - s, s, target.right, target.right.size() - s, s) != 0))
        {
          return false;
        }
    }

  if ((opt_self != 0) and (query.header == target.header))
    {
      return false;
    }
  if ((opt_selfid != 0) and
      (query.left == target.left) and
      (query.right == target.right))
    {
      return false;
    }

  return true;
}


auto paired_aligned_filters_pass(TavRecord const & query,
                                 TavRecord const & target,
                                 TavHit & hit,
                                 LinearMemoryAligner & lma) -> bool
{
  hit.left = align_one_end_stock_style(query.left, target.left, lma);
  hit.right = align_one_end_stock_style(query.right, target.right, lma);
  hit.aligned = true;

  hit.mismatches_total = hit.left.mismatches + hit.right.mismatches;
  hit.nwgaps_total = hit.left.nwgaps + hit.right.nwgaps;
  hit.nwalignment_cols_total = hit.left.nwalignmentlength + hit.right.nwalignmentlength;
  hit.internal_alignment_cols_total = hit.left.internal_alignmentlength + hit.right.internal_alignmentlength;
  hit.internal_gaps_total = hit.left.internal_gaps + hit.right.internal_gaps;
  hit.internal_indels_total = hit.left.internal_indels + hit.right.internal_indels;

  auto const matches_total = hit.left.matches + hit.right.matches;
  auto const ungapped_cols_total = matches_total + hit.mismatches_total;
  auto const query_len = static_cast<int>(query.left.size() + query.right.size());
  auto const target_len = static_cast<int>(target.left.size() + target.right.size());
  auto const shortest_total = std::min(query_len, target_len);
  auto const longest_total = std::max(query_len, target_len);

  auto const id0 = (shortest_total > 0)
    ? (100.0 * static_cast<double>(matches_total) / static_cast<double>(shortest_total))
    : 0.0;
  auto const id1 = (hit.nwalignment_cols_total > 0)
    ? (100.0 * static_cast<double>(matches_total) / static_cast<double>(hit.nwalignment_cols_total))
    : 0.0;
  auto const id2 = (hit.internal_alignment_cols_total > 0)
    ? (100.0 * static_cast<double>(matches_total) / static_cast<double>(hit.internal_alignment_cols_total))
    : 0.0;
  auto const id3 = (longest_total > 0)
    ? std::max(0.0,
               100.0 * (1.0 - (static_cast<double>(hit.mismatches_total + hit.nwgaps_total) /
                               static_cast<double>(longest_total))))
    : 0.0;
  auto const id4 = id1;

  switch (opt_iddef)
    {
    case 0:
      hit.id = id0;
      break;
    case 1:
      hit.id = id1;
      break;
    case 2:
      hit.id = id2;
      break;
    case 3:
      hit.id = id3;
      break;
    case 4:
      hit.id = id4;
      break;
    default:
      hit.id = id2;
      break;
    }

  if (ungapped_cols_total > 0)
    {
      hit.mid = 100.0 * static_cast<double>(matches_total) /
                static_cast<double>(ungapped_cols_total);
    }

  if (hit.id < (100.0 * opt_weak_id))
    {
      return false;
    }
  if (hit.mismatches_total > opt_maxsubs)
    {
      return false;
    }
  if (hit.internal_gaps_total > opt_maxgaps)
    {
      return false;
    }
  if ((opt_leftjust != 0) and
      ((hit.left.trim_q_left + hit.left.trim_t_left +
        hit.right.trim_q_left + hit.right.trim_t_left) > 0))
    {
      return false;
    }
  if ((opt_rightjust != 0) and
      ((hit.left.trim_q_right + hit.left.trim_t_right +
        hit.right.trim_q_right + hit.right.trim_t_right) > 0))
    {
      return false;
    }
  if (hit.internal_alignment_cols_total < opt_mincols)
    {
      return false;
    }
  if (static_cast<double>(ungapped_cols_total) < (opt_query_cov * static_cast<double>(query_len)))
    {
      return false;
    }
  if (static_cast<double>(ungapped_cols_total) < (opt_target_cov * static_cast<double>(target_len)))
    {
      return false;
    }
  if (hit.id > (100.0 * opt_maxid))
    {
      return false;
    }
  if (hit.mid < opt_mid)
    {
      return false;
    }
  if ((hit.mismatches_total + hit.internal_indels_total) > opt_maxdiffs)
    {
      return false;
    }

  return true;
}


auto trim_header_to_id(char const * header) -> std::string
{
  std::string h = header;
  auto const p = h.find_first_of(" \t");
  if (p != std::string::npos)
    {
      h.resize(p);
    }
  return h;
}


auto get_anchor_len(struct Parameters const & parameters,
                    int64_t len1,
                    int64_t len2) -> int64_t
{
  (void) parameters;
  auto const min_len = std::min(len1, len2);
  if (opt_fastq_trunclen > 0)
    {
      return std::min(min_len, opt_fastq_trunclen);
    }
  return min_len;
}


auto load_paired_records_from_fastx(char const * left_filename,
                                    char const * right_filename) -> TavLoadResult
{
  auto * left_h = fastx_open(left_filename);
  auto * right_h = fastx_open(right_filename);

  TavLoadResult result;

  while (fastx_next(left_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
    {
      if (!fastx_next(right_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
        {
          fatal("More forward records than reverse records in paired cluster_unoise input");
        }

      auto const left_len = static_cast<int64_t>(fastx_get_sequence_length(left_h));
      auto const right_len = static_cast<int64_t>(fastx_get_sequence_length(right_h));

      auto const left_short = left_len < opt_minseqlength;
      auto const right_short = right_len < opt_minseqlength;
      auto const left_long = left_len > opt_maxseqlength;
      auto const right_long = right_len > opt_maxseqlength;
      auto const drop_short = (opt_filter == 0) ? (left_short or right_short) : (left_short and right_short);
      auto const drop_long = (opt_filter == 0) ? (left_long or right_long) : (left_long and right_long);
      if (drop_short)
        {
          ++result.filtered_short;
          continue;
        }
      if (drop_long)
        {
          ++result.filtered_long;
          continue;
        }

      auto const left_seq = fastx_get_sequence(left_h);
      auto const right_seq = fastx_get_sequence(right_h);

      TavRecord r;
      r.tav_id = "";
      r.header = fastx_get_header(left_h);
      r.left.assign(left_seq, static_cast<std::size_t>(left_len));
      r.right.assign(right_seq, static_cast<std::size_t>(right_len));

      if ((opt_qmask == MASK_DUST) or ((opt_qmask == MASK_SOFT) and (opt_hardmask != 0)))
        {
          std::vector<char> left_buf(r.left.size() + 1U);
          std::copy(r.left.begin(), r.left.end(), left_buf.begin());
          left_buf[r.left.size()] = '\0';
          std::vector<char> right_buf(r.right.size() + 1U);
          std::copy(r.right.begin(), r.right.end(), right_buf.begin());
          right_buf[r.right.size()] = '\0';

          if (opt_qmask == MASK_DUST)
            {
              dust(left_buf.data(), static_cast<int>(r.left.size()));
              dust(right_buf.data(), static_cast<int>(r.right.size()));
            }
          else
            {
              hardmask(left_buf.data(), static_cast<int>(r.left.size()));
              hardmask(right_buf.data(), static_cast<int>(r.right.size()));
            }

          r.left.assign(left_buf.data(), r.left.size());
          r.right.assign(right_buf.data(), r.right.size());
        }

      r.abundance = fastx_get_abundance(left_h);
      r.first_seen = static_cast<int64_t>(result.records.size());
      result.records.push_back(r);
    }

  if (fastx_next(right_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
    {
      fatal("More reverse records than forward records in paired cluster_unoise input");
    }

  fastx_close(right_h);
  fastx_close(left_h);

  return result;
}


auto write_catalog(std::FILE * fp, std::vector<TavRecord> const & records) -> void
{
  std::fprintf(fp, "tav_id\tabundance\theader\tleft_anchor\tright_anchor\n");
  for (auto const & r : records)
    {
      std::fprintf(fp, "%s\t%" PRId64 "\t%s\t%s\t%s\n",
                   r.tav_id.c_str(),
                   r.abundance,
                   r.header.c_str(),
                   r.left.c_str(),
                   r.right.c_str());
    }
}


auto write_paired_fasta_interleaved(std::FILE * fp,
                                    std::vector<TavRecord> const & records,
                                    std::vector<double> const & scores,
                                    char const * score_name) -> void
{
  int ordinal = 1;
  for (std::size_t i = 0; i < records.size(); ++i)
    {
      auto const & r = records[i];
      auto const score = (i < scores.size()) ? scores[i] : 0.0;
      auto left_header = r.header + "/1";
      auto right_header = r.header + "/2";

      fasta_print_general(fp,
                          nullptr,
                          r.left.c_str(),
                          static_cast<int>(r.left.size()),
                          left_header.c_str(),
                          static_cast<int>(left_header.size()),
                          static_cast<unsigned int>(r.abundance),
                          ordinal,
                          -1.0,
                          -1,
                          -1,
                          opt_fasta_score ? score_name : nullptr,
                          score);
      ++ordinal;

      fasta_print_general(fp,
                          nullptr,
                          r.right.c_str(),
                          static_cast<int>(r.right.size()),
                          right_header.c_str(),
                          static_cast<int>(right_header.size()),
                          static_cast<unsigned int>(r.abundance),
                          ordinal,
                          -1.0,
                          -1,
                          -1,
                          opt_fasta_score ? score_name : nullptr,
                          score);
      ++ordinal;
    }
}


auto write_left_right_fasta(struct Parameters const & parameters,
                            std::vector<TavRecord> const & records) -> void
{
  std::FILE * fp_left = nullptr;
  std::FILE * fp_right = nullptr;

  if (parameters.opt_fastaout != nullptr)
    {
      fp_left = fopen_output(parameters.opt_fastaout);
      if (fp_left == nullptr)
        {
          fatal("Unable to open left-anchor FASTA output file for writing");
        }
    }

  if (parameters.opt_fastaout_rev != nullptr)
    {
      fp_right = fopen_output(parameters.opt_fastaout_rev);
      if (fp_right == nullptr)
        {
          fatal("Unable to open right-anchor FASTA output file for writing");
        }
    }

  int64_t index = 1;
  for (auto const & r : records)
    {
      if (fp_left != nullptr)
        {
          fasta_print_general(fp_left,
                              nullptr,
                              r.left.c_str(),
                              static_cast<int64_t>(r.left.size()),
                              r.tav_id.c_str(),
                              static_cast<int64_t>(r.tav_id.size()),
                              r.abundance,
                              index,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }
      if (fp_right != nullptr)
        {
          fasta_print_general(fp_right,
                              nullptr,
                              r.right.c_str(),
                              static_cast<int64_t>(r.right.size()),
                              r.tav_id.c_str(),
                              static_cast<int64_t>(r.tav_id.size()),
                              r.abundance,
                              index,
                              -1.0,
                              -1,
                              -1,
                              nullptr,
                              0.0);
        }
      ++index;
    }

  if (fp_left != nullptr)
    {
      std::fclose(fp_left);
    }
  if (fp_right != nullptr)
    {
      std::fclose(fp_right);
    }
}

} // namespace


auto tav_is_catalog_file(char const * filename) -> bool
{
  std::ifstream in(filename);
  if (!in)
    {
      return false;
    }
  std::string line;
  if (!std::getline(in, line))
    {
      return false;
    }
  if (line.rfind("tav_id\t", 0) == 0)
    {
      return true;
    }
  return false;
}


auto tav_fastx_uniques(struct Parameters const & parameters) -> void
{
  if (parameters.opt_reverse == nullptr)
    {
      fatal("TAV fastx_uniques requires paired reads with --reverse");
    }
  if ((parameters.opt_tabbedout == nullptr) and
      (parameters.opt_fastaout == nullptr) and
      (parameters.opt_fastaout_rev == nullptr))
    {
      fatal("TAV fastx_uniques requires output via --tabbedout and/or --fastaout/--fastaout_rev");
    }

  auto * fwd = fastx_open(parameters.opt_fastx_uniques);
  auto * rev = fastx_open(parameters.opt_reverse);

  std::unordered_map<std::string, std::size_t> pair_to_index;
  std::vector<TavRecord> records;
  int64_t read_ordinal = 0;

  int64_t skipped_short = 0;

  while (fastx_next(fwd, not opt_notrunclabels, chrmap_no_change_vector.data()))
    {
      if (!fastx_next(rev, not opt_notrunclabels, chrmap_no_change_vector.data()))
        {
          fatal("More forward reads than reverse reads");
        }

      auto const f_len = static_cast<int64_t>(fastx_get_sequence_length(fwd));
      auto const r_len = static_cast<int64_t>(fastx_get_sequence_length(rev));
      auto const anchor_len = get_anchor_len(parameters, f_len, r_len);
      if ((anchor_len <= 0) or (f_len < anchor_len) or (r_len < anchor_len))
        {
          ++skipped_short;
          continue;
        }

      std::string left(fastx_get_sequence(fwd), static_cast<std::size_t>(anchor_len));

      std::vector<char> revcomp(static_cast<std::size_t>(r_len) + 1U);
      reverse_complement(revcomp.data(), fastx_get_sequence(rev), r_len);
      std::string right(revcomp.data(), static_cast<std::size_t>(anchor_len));

      auto const key = left + '\t' + right;
      auto const abundance = fastx_get_abundance(fwd);
      auto pos = pair_to_index.find(key);
      if (pos == pair_to_index.end())
        {
          TavRecord r;
          r.tav_id = "";
          r.header = fastx_get_header(fwd);
          r.left = left;
          r.right = right;
          r.abundance = abundance;
          r.first_seen = read_ordinal;
          pair_to_index.emplace(key, records.size());
          records.push_back(r);
        }
      else
        {
          records[pos->second].abundance += abundance;
        }

      ++read_ordinal;
    }

  if (fastx_next(rev, not opt_notrunclabels, chrmap_no_change_vector.data()))
    {
      fatal("More reverse reads than forward reads");
    }

  fastx_close(rev);
  fastx_close(fwd);

  std::sort(records.begin(), records.end(),
            [](TavRecord const & lhs, TavRecord const & rhs) {
              if (lhs.abundance != rhs.abundance)
                {
                  return lhs.abundance > rhs.abundance;
                }
              if (lhs.header != rhs.header)
                {
                  return lhs.header < rhs.header;
                }
              return lhs.first_seen < rhs.first_seen;
            });

  for (std::size_t i = 0; i < records.size(); ++i)
    {
      records[i].tav_id = records[i].header;
    }

  if (parameters.opt_tabbedout != nullptr)
    {
      auto * fp = fopen_output(parameters.opt_tabbedout);
      if (fp == nullptr)
        {
          fatal("Unable to open TAV catalog output file for writing");
        }
      write_catalog(fp, records);
      std::fclose(fp);
    }

  write_left_right_fasta(parameters, records);

  if (!opt_quiet)
    {
      std::fprintf(stderr,
                   "TAV fastx_uniques: %zu unique paired anchors, %" PRId64 " pairs skipped for anchor length constraints\n",
                   records.size(),
                   skipped_short);
    }
}


auto tav_cluster_unoise(struct Parameters const & parameters) -> void
{
  if (parameters.opt_reverse == nullptr)
    {
      fatal("Paired cluster_unoise requires --reverse (catalog mode is not supported)");
    }

  auto load_result = load_paired_records_from_fastx(parameters.opt_cluster_unoise,
                                                    parameters.opt_reverse);
  auto & input_records = load_result.records;
  auto const loaded_pairs = input_records.size() + static_cast<std::size_t>(load_result.filtered_short + load_result.filtered_long);

  if (input_records.empty())
    {
      fatal("Input for paired cluster_unoise is empty");
    }

  std::vector<TavRecord> records;
  records.reserve(input_records.size());

  int64_t filtered_small = 0;
  int64_t filtered_large = 0;
  for (auto const & r : input_records)
    {
      if (r.abundance < opt_minsize)
        {
          ++filtered_small;
          continue;
        }
      if (r.abundance > opt_maxsize)
        {
          ++filtered_large;
          continue;
        }
      records.push_back(r);
    }

  if (records.empty())
    {
      fatal("No TAV records pass --minsize/--maxsize filtering");
    }

  std::sort(records.begin(), records.end(),
            [](TavRecord const & lhs, TavRecord const & rhs) {
              if (lhs.abundance != rhs.abundance)
                {
                  return lhs.abundance > rhs.abundance;
                }
              if (lhs.header != rhs.header)
                {
                  return lhs.header < rhs.header;
                }
              if (lhs.left != rhs.left)
                {
                  return lhs.left < rhs.left;
                }
              return lhs.right < rhs.right;
            });

  std::vector<TavRecord> centroids;
  TavKmerIndex kmer_index;
  auto * kmer_unique_handle = unique_init();
  auto const wordlength = static_cast<int>(opt_wordlength);
  struct Scoring scoring;
  scoring.match = opt_match;
  scoring.mismatch = opt_mismatch;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_open_query_left = opt_gap_open_query_left;
  scoring.gap_open_target_left = opt_gap_open_target_left;
  scoring.gap_open_target_interior = opt_gap_open_target_interior;
  scoring.gap_open_query_right = opt_gap_open_query_right;
  scoring.gap_open_target_right = opt_gap_open_target_right;
  scoring.gap_extension_query_left = opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opt_gap_extension_target_left;
  scoring.gap_extension_target_interior = opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opt_gap_extension_target_right;
  LinearMemoryAligner lma(scoring);

  int64_t total_candidates = 0;
  int64_t total_aligned = 0;
  int64_t total_accepted = 0;
  int64_t total_weak = 0;
  int64_t total_rejected = 0;

  for (auto const & q : records)
    {
      std::vector<TavHit> hits;

      if (!centroids.empty())
        {
          std::vector<unsigned int> centroid_kmer_scores(centroids.size(), 0U);

          unsigned int qk_left = 0;
          unsigned int const * qk_left_list = nullptr;
          unique_count(kmer_unique_handle,
                       wordlength,
                       static_cast<int>(q.left.size()),
                       q.left.c_str(),
                       &qk_left,
                       &qk_left_list,
                       opt_qmask);
          for (auto i = 0U; i < qk_left; ++i)
            {
              auto const left_it = kmer_index.left_postings.find(qk_left_list[i]);
              if (left_it == kmer_index.left_postings.end())
                {
                  continue;
                }
              for (auto const centroid_index : left_it->second)
                {
                  ++centroid_kmer_scores[static_cast<std::size_t>(centroid_index)];
                }
            }

          unsigned int qk_right = 0;
          unsigned int const * qk_right_list = nullptr;
          unique_count(kmer_unique_handle,
                       wordlength,
                       static_cast<int>(q.right.size()),
                       q.right.c_str(),
                       &qk_right,
                       &qk_right_list,
                       opt_qmask);
          for (auto i = 0U; i < qk_right; ++i)
            {
              auto const right_it = kmer_index.right_postings.find(qk_right_list[i]);
              if (right_it == kmer_index.right_postings.end())
                {
                  continue;
                }
              for (auto const centroid_index : right_it->second)
                {
                  ++centroid_kmer_scores[static_cast<std::size_t>(centroid_index)];
                }
            }

          auto const qk_total = qk_left + qk_right;
          auto const minmatches = std::min(static_cast<unsigned int>(opt_minwordmatches), qk_total);
          auto const maxaccepts_effective =
            (opt_maxaccepts == 0)
            ? static_cast<int64_t>(centroids.size())
            : std::min<int64_t>(opt_maxaccepts, static_cast<int64_t>(centroids.size()));
          auto const maxrejects_effective =
            (opt_maxrejects == 0)
            ? static_cast<int64_t>(centroids.size())
            : std::min<int64_t>(opt_maxrejects, static_cast<int64_t>(centroids.size()));

          std::vector<std::pair<int, unsigned int>> candidates;
          auto const tophits = std::max<int64_t>(
            1,
            std::min<int64_t>(static_cast<int64_t>(centroids.size()),
                              maxrejects_effective + maxaccepts_effective + static_cast<int64_t>(MAXDELAYED)));
          auto * candidate_heap = minheap_init(static_cast<int>(tophits));
          for (std::size_t i = 0; i < centroids.size(); ++i)
            {
              auto const score = centroid_kmer_scores[i];
              if (score >= minmatches)
                {
                  elem_t novel;
                  novel.count = score;
                  novel.seqno = static_cast<unsigned int>(i);
                  novel.length = static_cast<unsigned int>(centroids[i].left.size() + centroids[i].right.size());
                  minheap_add(candidate_heap, &novel);
                }
            }

          minheap_sort(candidate_heap);
          candidates.reserve(static_cast<std::size_t>(candidate_heap->count));
          while (!minheap_isempty(candidate_heap))
            {
              auto const e = minheap_poplast(candidate_heap);
              candidates.emplace_back(static_cast<int>(e.seqno), e.count);
            }
          minheap_exit(candidate_heap);

          total_candidates += static_cast<int64_t>(candidates.size());

          auto accepts = 0;
          auto rejects = 0;
          std::vector<std::size_t> delayed;
          delayed.reserve(MAXDELAYED);

          auto process_delayed = [&]() -> void {
            for (auto const hit_index : delayed)
              {
                if ((accepts >= maxaccepts_effective) or (rejects >= maxrejects_effective))
                  {
                    break;
                  }

                auto & hit = hits[hit_index];
                auto const & centroid = centroids[hit.centroid_index];
                ++total_aligned;

                if (!paired_aligned_filters_pass(q, centroid, hit, lma))
                  {
                    hit.rejected = true;
                    hit.weak = false;
                    ++rejects;
                    ++total_rejected;
                    continue;
                  }

                auto const skew = static_cast<double>(q.abundance) /
                                  static_cast<double>(centroid.abundance);
                auto const beta = std::pow(2.0,
                                           -1.0 - (opt_unoise_alpha * static_cast<double>(hit.mismatches_total)));
                if ((skew <= beta) or (hit.mismatches_total == 0))
                  {
                    hit.accepted = true;
                    hit.weak = false;
                    ++accepts;
                    ++total_accepted;
                  }
                else
                  {
                    hit.rejected = true;
                    hit.weak = true;
                    ++rejects;
                    ++total_rejected;
                    ++total_weak;
                  }
              }
            delayed.clear();
          };

          for (auto const & c : candidates)
            {
              if ((accepts >= maxaccepts_effective) or (rejects >= maxrejects_effective))
                {
                  break;
                }

              TavHit hit;
              hit.centroid_index = c.first;
              hit.kmer_score = c.second;
              hits.push_back(hit);

              if (paired_unaligned_filters_pass(q, centroids[c.first]))
                {
                  delayed.push_back(hits.size() - 1U);
                  if (delayed.size() == MAXDELAYED)
                    {
                      process_delayed();
                    }
                }
              else
                {
                  hits.back().rejected = true;
                  hits.back().weak = false;
                  ++rejects;
                  ++total_rejected;
                }
            }

          if (!delayed.empty())
            {
              process_delayed();
            }
        }

      auto best_index = -1;
      double best_id = -1.0;
      for (auto const & hit : hits)
        {
          if (!hit.accepted)
            {
              continue;
            }
          auto const idx = hit.centroid_index;
          if ((best_index < 0) or
              (centroids[idx].abundance > centroids[best_index].abundance) or
              ((centroids[idx].abundance == centroids[best_index].abundance) and (hit.id > best_id)) or
              ((centroids[idx].abundance == centroids[best_index].abundance) and (hit.id == best_id) and (idx < best_index)))
            {
              best_index = idx;
              best_id = hit.id;
            }
        }

      if (best_index >= 0)
        {
          centroids[best_index].abundance += q.abundance;
        }
      else
        {
          centroids.push_back(q);
          auto const centroid_index = static_cast<int>(centroids.size() - 1U);

          unsigned int left_unique_count = 0;
          unsigned int const * left_unique_list = nullptr;
          unique_count(kmer_unique_handle,
                       wordlength,
                       static_cast<int>(centroids.back().left.size()),
                       centroids.back().left.c_str(),
                       &left_unique_count,
                       &left_unique_list,
                       opt_qmask);
          for (auto i = 0U; i < left_unique_count; ++i)
            {
              kmer_index.left_postings[left_unique_list[i]].push_back(centroid_index);
            }

          unsigned int right_unique_count = 0;
          unsigned int const * right_unique_list = nullptr;
          unique_count(kmer_unique_handle,
                       wordlength,
                       static_cast<int>(centroids.back().right.size()),
                       centroids.back().right.c_str(),
                       &right_unique_count,
                       &right_unique_list,
                       opt_qmask);
          for (auto i = 0U; i < right_unique_count; ++i)
            {
              kmer_index.right_postings[right_unique_list[i]].push_back(centroid_index);
            }
        }
    }

  unique_exit(kmer_unique_handle);

  for (std::size_t i = 0; i < centroids.size(); ++i)
    {
      centroids[i].tav_id = centroids[i].header;
    }

  if (opt_centroids != nullptr)
    {
      if (parameters.opt_fastaout_rev == nullptr)
        {
          fatal("Paired cluster_unoise with --centroids also requires --fastaout_rev for paired centroid FASTA output");
        }

      auto centroid_params = parameters;
      centroid_params.opt_fastaout = opt_centroids;
      write_left_right_fasta(centroid_params, centroids);
    }

  if (parameters.opt_tabbedout != nullptr)
    {
      auto * fp = fopen_output(parameters.opt_tabbedout);
      if (fp == nullptr)
        {
          fatal("Unable to open denoised TAV tabbed output file for writing");
        }
      write_catalog(fp, centroids);
      std::fclose(fp);
    }

  if (parameters.opt_fastaout != nullptr)
    {
      write_left_right_fasta(parameters, centroids);
    }

  if (!opt_quiet)
    {
      auto const filter_mode_label = (opt_filter == 0) ? "any" : "both";
      std::fprintf(stderr,
                   "TAV cluster_unoise: %zu input pairs (%" PRId64 " filtered short, %" PRId64 " filtered long, filter=%s; %" PRId64 " filtered by minsize, %" PRId64 " filtered by maxsize) -> %zu centroids; candidates=%" PRId64 ", aligned=%" PRId64 ", accepted=%" PRId64 ", weak=%" PRId64 ", rejected=%" PRId64 "\n",
                   loaded_pairs,
                   load_result.filtered_short,
                   load_result.filtered_long,
                   filter_mode_label,
                   filtered_small,
                   filtered_large,
                   centroids.size(),
                   total_candidates,
                   total_aligned,
                   total_accepted,
                   total_weak,
                   total_rejected);
    }
}


auto tav_uchime3_denovo(struct Parameters const & parameters) -> void
{
  if (parameters.opt_reverse == nullptr)
    {
      fatal("Paired uchime3_denovo requires --reverse (catalog mode is not supported)");
    }

  auto load_result = load_paired_records_from_fastx(parameters.opt_uchime3_denovo,
                                                    parameters.opt_reverse);
  auto records = std::move(load_result.records);
  if (records.empty())
    {
      fatal("Input for paired uchime3_denovo is empty");
    }

  for (auto & r : records)
    {
      r.tav_id = r.header;
    }

  std::sort(records.begin(), records.end(),
            [](TavRecord const & lhs, TavRecord const & rhs) {
              if (lhs.abundance != rhs.abundance)
                {
                  return lhs.abundance > rhs.abundance;
                }
              if (lhs.header != rhs.header)
                {
                  return lhs.header < rhs.header;
                }
              if (lhs.left != rhs.left)
                {
                  return lhs.left < rhs.left;
                }
              return lhs.right < rhs.right;
            });

  std::vector<TavRecord> nonchim;
  std::vector<TavRecord> chim;
  std::vector<double> nonchim_scores;
  std::vector<double> chim_scores;
  std::vector<std::size_t> parent_pool_indices;
  TavKmerIndex parent_kmer_index;
  parent_pool_indices.reserve(records.size());
  auto * parent_unique_handle = unique_init();
  auto const chimera_wordlength = static_cast<int>(opt_wordlength);
  static constexpr int chimera_window = 32;
  static constexpr int chimera_maxaccepts = 4;
  static constexpr int chimera_maxrejects = 16;
  static constexpr int chimera_tophits = chimera_maxaccepts + chimera_maxrejects;

  struct Scoring scoring;
  scoring.match = opt_match;
  scoring.mismatch = opt_mismatch;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_open_query_left = opt_gap_open_query_left;
  scoring.gap_open_target_left = opt_gap_open_target_left;
  scoring.gap_open_target_interior = opt_gap_open_target_interior;
  scoring.gap_open_query_right = opt_gap_open_query_right;
  scoring.gap_open_target_right = opt_gap_open_target_right;
  scoring.gap_extension_query_left = opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opt_gap_extension_target_left;
  scoring.gap_extension_target_interior = opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opt_gap_extension_target_right;
  LinearMemoryAligner lma(scoring);

  auto print_stripped_header = [](std::FILE * fp, std::string const & header) -> void {
    header_fprint_strip(fp,
                        header.c_str(),
                        static_cast<int64_t>(header.size()),
                        opt_xsize,
                        opt_xee,
                        opt_xlength);
  };

  std::FILE * fp_report = nullptr;
  if (parameters.opt_tabbedout != nullptr)
    {
      fp_report = fopen_output(parameters.opt_tabbedout);
      if (fp_report == nullptr)
        {
          fatal("Unable to open TAV chimera report output file for writing");
        }
      std::fprintf(fp_report,
                   "query_tav\tparent_a\tparent_b\tbreakpoint_class\tbest_one_score\tbest_two_score\tdelta\n");
    }

  std::FILE * fp_uchimeout_paired = nullptr;
  if (opt_uchimeout != nullptr)
    {
      fp_uchimeout_paired = fopen_output(opt_uchimeout);
      if (fp_uchimeout_paired == nullptr)
        {
          fatal("Unable to open uchimeout output file for writing");
        }
    }

  std::FILE * fp_uchimealns_paired = nullptr;
  if (opt_uchimealns != nullptr)
    {
      fp_uchimealns_paired = fopen_output(opt_uchimealns);
      if (fp_uchimealns_paired == nullptr)
        {
          fatal("Unable to open uchimealns output file for writing");
        }
    }

  std::FILE * fp_borderline_paired = nullptr;
  if (opt_borderline != nullptr)
    {
      fp_borderline_paired = fopen_output(opt_borderline);
      if (fp_borderline_paired == nullptr)
        {
          fatal("Unable to open borderline output file for writing");
        }
    }

  auto write_uchimeout_no_parents = [&](TavRecord const & q, double const best_h) -> void {
    if (fp_uchimeout_paired == nullptr)
      {
        return;
      }

    std::fprintf(fp_uchimeout_paired, "%.4f\t", std::max(0.0, best_h));
    print_stripped_header(fp_uchimeout_paired, q.header);
    if (opt_uchimeout5 != 0)
      {
        std::fprintf(fp_uchimeout_paired,
                     "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
      }
    else
      {
        std::fprintf(fp_uchimeout_paired,
                     "\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n");
      }
  };

  auto add_parent_to_index = [&](std::size_t const parent_index) -> void {
    unsigned int left_unique_count = 0;
    unsigned int const * left_unique_list = nullptr;
    unique_count(parent_unique_handle,
                 chimera_wordlength,
                 static_cast<int>(records[parent_index].left.size()),
                 records[parent_index].left.c_str(),
                 &left_unique_count,
                 &left_unique_list,
                 opt_qmask);
    for (auto i = 0U; i < left_unique_count; ++i)
      {
        parent_kmer_index.left_postings[left_unique_list[i]].push_back(static_cast<int>(parent_index));
      }

    unsigned int right_unique_count = 0;
    unsigned int const * right_unique_list = nullptr;
    unique_count(parent_unique_handle,
                 chimera_wordlength,
                 static_cast<int>(records[parent_index].right.size()),
                 records[parent_index].right.c_str(),
                 &right_unique_count,
                 &right_unique_list,
                 opt_qmask);
    for (auto i = 0U; i < right_unique_count; ++i)
      {
        parent_kmer_index.right_postings[right_unique_list[i]].push_back(static_cast<int>(parent_index));
      }
  };

  auto compute_match_vector = [&](std::string const & query,
                                  std::string const & target,
                                  std::string const & cigar) -> std::vector<unsigned char> {
    std::vector<unsigned char> match(query.size(), 0U);
    auto qpos = std::size_t{0};
    auto tpos = std::size_t{0};

    char const * p = cigar.c_str();
    while (*p != '\0')
      {
        int run = 1;
        int scanlength = 0;
        if (std::sscanf(p, "%d%n", &run, &scanlength) == 1)
          {
            p += scanlength;
          }
        if (*p == '\0')
          {
            break;
          }
        char const op = *p;
        ++p;

        if (op == 'M')
          {
            for (int j = 0; j < run; ++j)
              {
                if ((qpos < query.size()) and (tpos < target.size()) and
                    ((map_4bit(query[qpos]) & map_4bit(target[tpos])) != 0U))
                  {
                    match[qpos] = 1U;
                  }
                if (qpos < query.size())
                  {
                    ++qpos;
                  }
                if (tpos < target.size())
                  {
                    ++tpos;
                  }
              }
          }
        else if (op == 'I')
          {
            tpos += static_cast<std::size_t>(run);
          }
        else if (op == 'D')
          {
            qpos += static_cast<std::size_t>(run);
          }
      }
    return match;
  };

  auto align_cigar_string = [&](std::string const & query,
                                std::string const & target) -> std::string {
    std::vector<char> qbuf(query.begin(), query.end());
    qbuf.push_back('\0');
    std::vector<char> tbuf(target.begin(), target.end());
    tbuf.push_back('\0');
    auto * cigar = xstrdup(lma.align(qbuf.data(),
                                     tbuf.data(),
                                     static_cast<int64_t>(query.size()),
                                     static_cast<int64_t>(target.size())));
    std::string out(cigar);
    xfree(cigar);
    return out;
  };

  struct TavParentCandidate
  {
    std::size_t record_index = 0;
    unsigned int kmer_score = 0U;
    std::string cigar_left;
    std::string cigar_right;
    std::vector<unsigned char> match_axis;
  };

  auto select_best_parent_pair = [&](std::vector<TavParentCandidate> const & candidates,
                                     std::size_t const query_len_total) -> std::pair<int, int> {
    if ((candidates.size() < 2U) or (query_len_total < static_cast<std::size_t>(chimera_window)))
      {
        return {-1, -1};
      }

    auto const cand_count = candidates.size();
    std::vector<std::vector<unsigned char>> match(cand_count);
    std::vector<std::vector<int>> smooth(cand_count, std::vector<int>(query_len_total, 0));
    std::vector<int> maxsmooth(query_len_total, 0);
    for (std::size_t i = 0; i < cand_count; ++i)
      {
        match[i] = candidates[i].match_axis;
      }

    std::vector<bool> cand_selected(cand_count, false);
    std::array<int, 2> best_parent_idx {{-1, -1}};

    for (auto round = 0; round < 2; ++round)
      {
        if (round > 0)
          {
            auto const prev = best_parent_idx[static_cast<std::size_t>(round - 1)];
            if (prev < 0)
              {
                break;
              }
            for (std::size_t qpos = static_cast<std::size_t>(chimera_window - 1); qpos < query_len_total; ++qpos)
              {
                if ((smooth[static_cast<std::size_t>(prev)][qpos] == maxsmooth[qpos]) and (maxsmooth[qpos] > 0))
                  {
                    auto const start = qpos + 1U - static_cast<std::size_t>(chimera_window);
                    for (std::size_t i = start; i <= qpos; ++i)
                      {
                        for (std::size_t c = 0; c < cand_count; ++c)
                          {
                            match[c][i] = 0U;
                          }
                      }
                  }
              }
          }

        std::fill(maxsmooth.begin(), maxsmooth.end(), 0);
        for (std::size_t c = 0; c < cand_count; ++c)
          {
            if (cand_selected[c])
              {
                continue;
              }

            auto sum = 0;
            for (std::size_t qpos = 0; qpos < query_len_total; ++qpos)
              {
                sum += match[c][qpos];
                if (qpos >= static_cast<std::size_t>(chimera_window))
                  {
                    sum -= match[c][qpos - static_cast<std::size_t>(chimera_window)];
                  }
                if (qpos >= static_cast<std::size_t>(chimera_window - 1))
                  {
                    smooth[c][qpos] = sum;
                    maxsmooth[qpos] = std::max(maxsmooth[qpos], sum);
                  }
              }
          }

        std::vector<int> wins(cand_count, 0);
        for (std::size_t qpos = static_cast<std::size_t>(chimera_window - 1); qpos < query_len_total; ++qpos)
          {
            if (maxsmooth[qpos] == 0)
              {
                continue;
              }
            for (std::size_t c = 0; c < cand_count; ++c)
              {
                if (cand_selected[c])
                  {
                    continue;
                  }
                if (smooth[c][qpos] == maxsmooth[qpos])
                  {
                    ++wins[c];
                  }
              }
          }

        auto maxwins = 0;
        auto best = -1;
        for (std::size_t c = 0; c < cand_count; ++c)
          {
            if (cand_selected[c])
              {
                continue;
              }
            if (wins[c] > maxwins)
              {
                maxwins = wins[c];
                best = static_cast<int>(c);
              }
          }

        if (best < 0)
          {
            break;
          }
        best_parent_idx[static_cast<std::size_t>(round)] = best;
        cand_selected[static_cast<std::size_t>(best)] = true;
      }

    return {best_parent_idx[0], best_parent_idx[1]};
  };

  struct TavPairEval
  {
    bool valid = false;
    bool reverse_orientation = false;
    int best_i = -1;
    double best_h = 0.0;
    int cols = 0;
    int match_QA = 0;
    int match_QB = 0;
    int match_AB = 0;
    int match_QM = 0;
    double QA = 0.0;
    double QB = 0.0;
    double AB = 0.0;
    double QT = 0.0;
    double QM = 0.0;
    double divdiff = 0.0;
    int left_y = 0;
    int left_n = 0;
    int left_a = 0;
    int right_y = 0;
    int right_n = 0;
    int right_a = 0;
    std::size_t left_aln_len = 0U;
    std::string q_aln;
    std::string a_aln;
    std::string b_aln;
    std::string diffs;
    std::string votes;
    std::string model;
  };

  auto evaluate_parent_pair = [&](TavRecord const & q,
                                  TavRecord const & parent_a,
                                  TavRecord const & parent_b,
                                  std::string const & cigar_left_a,
                                  std::string const & cigar_right_a,
                                  std::string const & cigar_left_b,
                                  std::string const & cigar_right_b) -> TavPairEval {
    TavPairEval eval;

    auto build_end_alignments = [](std::string const & query,
                                   std::string const & target_a,
                                   std::string const & target_b,
                                   std::string const & cigar_a,
                                   std::string const & cigar_b,
                                   std::string & q_aln,
                                   std::string & a_aln,
                                   std::string & b_aln) -> void {
      std::vector<int> maxi(query.size() + 1U, 0);

      auto fill_maxi = [&](std::string const & cigar) -> void {
        auto qpos = std::size_t{0};
        char const * p = cigar.c_str();
        while (*p != '\0')
          {
            int run = 1;
            int scanlength = 0;
            if (std::sscanf(p, "%d%n", &run, &scanlength) == 1)
              {
                p += scanlength;
              }
            if (*p == '\0')
              {
                break;
              }
            auto const op = *p;
            ++p;

            if (op == 'I')
              {
                if (qpos < maxi.size())
                  {
                    maxi[qpos] = std::max(maxi[qpos], run);
                  }
              }
            else if ((op == 'M') or (op == 'D'))
              {
                qpos += static_cast<std::size_t>(run);
              }
          }
      };

      fill_maxi(cigar_a);
      fill_maxi(cigar_b);

      q_aln.clear();
      q_aln.reserve(query.size());
      for (std::size_t i = 0; i < query.size(); ++i)
        {
          q_aln.append(static_cast<std::size_t>(maxi[i]), '-');
          q_aln.push_back(map_uppercase(query[i]));
        }
      q_aln.append(static_cast<std::size_t>(maxi[query.size()]), '-');

      auto build_parent_alignment = [&](std::string const & target,
                                        std::string const & cigar,
                                        std::string & out) -> void {
        out.clear();
        out.reserve(q_aln.size());

        auto qpos = std::size_t{0};
        auto tpos = std::size_t{0};
        auto inserted = false;

        char const * p = cigar.c_str();
        while (*p != '\0')
          {
            int run = 1;
            int scanlength = 0;
            if (std::sscanf(p, "%d%n", &run, &scanlength) == 1)
              {
                p += scanlength;
              }
            if (*p == '\0')
              {
                break;
              }
            auto const op = *p;
            ++p;

            if (op == 'I')
              {
                auto const ins = (qpos < maxi.size()) ? maxi[qpos] : 0;
                for (auto i = 0; i < ins; ++i)
                  {
                    if ((i < run) and (tpos < target.size()))
                      {
                        out.push_back(map_uppercase(target[tpos]));
                        ++tpos;
                      }
                    else
                      {
                        out.push_back('-');
                      }
                  }
                inserted = true;
              }
            else if ((op == 'M') or (op == 'D'))
              {
                for (int i = 0; i < run; ++i)
                  {
                    if (not inserted)
                      {
                        auto const ins = (qpos < maxi.size()) ? maxi[qpos] : 0;
                        out.append(static_cast<std::size_t>(ins), '-');
                      }

                    if (op == 'M')
                      {
                        if (tpos < target.size())
                          {
                            out.push_back(map_uppercase(target[tpos]));
                            ++tpos;
                          }
                        else
                          {
                            out.push_back('N');
                          }
                      }
                    else
                      {
                        out.push_back('-');
                      }

                    ++qpos;
                    inserted = false;
                  }
              }
          }

        if (not inserted)
          {
            auto const ins = (qpos < maxi.size()) ? maxi[qpos] : 0;
            out.append(static_cast<std::size_t>(ins), '-');
          }
      };

      build_parent_alignment(target_a, cigar_a, a_aln);
      build_parent_alignment(target_b, cigar_b, b_aln);
    };

    std::string left_qaln;
    std::string left_aaln;
    std::string left_baln;
    build_end_alignments(q.left,
                         parent_a.left,
                         parent_b.left,
                         cigar_left_a,
                         cigar_left_b,
                         left_qaln,
                         left_aaln,
                         left_baln);

    std::string right_qaln;
    std::string right_aaln;
    std::string right_baln;
    build_end_alignments(q.right,
                         parent_a.right,
                         parent_b.right,
                         cigar_right_a,
                         cigar_right_b,
                         right_qaln,
                         right_aaln,
                         right_baln);

    eval.left_aln_len = left_qaln.size();
    eval.q_aln = left_qaln + right_qaln;
    eval.a_aln = left_aaln + right_aaln;
    eval.b_aln = left_baln + right_baln;

    auto const alnlen = static_cast<int>(eval.q_aln.size());
    if (alnlen == 0)
      {
        return eval;
      }

    eval.diffs.assign(static_cast<std::size_t>(alnlen), ' ');
    std::vector<bool> ignore(static_cast<std::size_t>(alnlen), false);

    for (auto i = 0; i < alnlen; ++i)
      {
        auto const qsym = map_4bit(eval.q_aln[static_cast<std::size_t>(i)]);
        auto const p1sym = map_4bit(eval.a_aln[static_cast<std::size_t>(i)]);
        auto const p2sym = map_4bit(eval.b_aln[static_cast<std::size_t>(i)]);

        if ((qsym == 0U) or (p1sym == 0U) or (p2sym == 0U))
          {
            ignore[static_cast<std::size_t>(i)] = true;
            if (i > 0)
              {
                ignore[static_cast<std::size_t>(i - 1)] = true;
              }
            if (i < (alnlen - 1))
              {
                ignore[static_cast<std::size_t>(i + 1)] = true;
              }
          }

        if (is_ambiguous_4bit(qsym) or is_ambiguous_4bit(p1sym) or is_ambiguous_4bit(p2sym))
          {
            ignore[static_cast<std::size_t>(i)] = true;
          }

        if ((p1sym != 0U) and (p1sym != qsym))
          {
            eval.a_aln[static_cast<std::size_t>(i)] = std::tolower(eval.a_aln[static_cast<std::size_t>(i)]);
          }
        if ((p2sym != 0U) and (p2sym != qsym))
          {
            eval.b_aln[static_cast<std::size_t>(i)] = std::tolower(eval.b_aln[static_cast<std::size_t>(i)]);
          }

        char diff = ' ';
        if ((qsym != 0U) and (p1sym != 0U) and (p2sym != 0U))
          {
            if (p1sym == p2sym)
              {
                diff = (qsym == p1sym) ? ' ' : 'N';
              }
            else
              {
                if (qsym == p1sym)
                  {
                    diff = 'A';
                  }
                else if (qsym == p2sym)
                  {
                    diff = 'B';
                  }
                else
                  {
                    diff = '?';
                  }
              }
          }
        eval.diffs[static_cast<std::size_t>(i)] = diff;
      }

    auto sumA = 0;
    auto sumB = 0;
    auto sumN = 0;
    for (auto i = 0; i < alnlen; ++i)
      {
        if (ignore[static_cast<std::size_t>(i)])
          {
            continue;
          }
        auto const diff = eval.diffs[static_cast<std::size_t>(i)];
        if (diff == 'A')
          {
            ++sumA;
          }
        else if (diff == 'B')
          {
            ++sumB;
          }
        else if (diff != ' ')
          {
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

    for (auto i = 0; i < alnlen; ++i)
      {
        if (ignore[static_cast<std::size_t>(i)])
          {
            continue;
          }
        auto const diff = eval.diffs[static_cast<std::size_t>(i)];
        if (diff == ' ')
          {
            continue;
          }

        if (diff == 'A')
          {
            ++left_y;
            --right_n;
          }
        else if (diff == 'B')
          {
            ++left_n;
            --right_y;
          }
        else
          {
            ++left_a;
            --right_a;
          }

        if ((left_y > left_n) and (right_y > right_n))
          {
            auto const left_h = static_cast<double>(left_y) /
                                ((opt_xn * (static_cast<double>(left_n) + opt_dn)) + static_cast<double>(left_a));
            auto const right_h = static_cast<double>(right_y) /
                                 ((opt_xn * (static_cast<double>(right_n) + opt_dn)) + static_cast<double>(right_a));
            auto const h = left_h * right_h;
            if (h > best_h)
              {
                best_h = h;
                best_i = i;
                best_is_reverse = false;
                best_left_y = left_y;
                best_left_n = left_n;
                best_left_a = left_a;
                best_right_y = right_y;
                best_right_n = right_n;
                best_right_a = right_a;
              }
          }
        else if ((left_n > left_y) and (right_n > right_y))
          {
            auto const left_h = static_cast<double>(left_n) /
                                ((opt_xn * (static_cast<double>(left_y) + opt_dn)) + static_cast<double>(left_a));
            auto const right_h = static_cast<double>(right_n) /
                                 ((opt_xn * (static_cast<double>(right_y) + opt_dn)) + static_cast<double>(right_a));
            auto const h = left_h * right_h;
            if (h > best_h)
              {
                best_h = h;
                best_i = i;
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

    if (best_h < 0.0)
      {
        return eval;
      }

    if (best_is_reverse)
      {
        for (auto i = 0; i < alnlen; ++i)
          {
            auto & diff = eval.diffs[static_cast<std::size_t>(i)];
            if (diff == 'A')
              {
                diff = 'B';
              }
            else if (diff == 'B')
              {
                diff = 'A';
              }
          }
      }

    eval.votes.assign(static_cast<std::size_t>(alnlen), ' ');
    eval.model.assign(static_cast<std::size_t>(alnlen), ' ');
    for (auto i = 0; i < alnlen; ++i)
      {
        auto const m = (i <= best_i) ? 'A' : 'B';
        eval.model[static_cast<std::size_t>(i)] = m;
        char vote = ' ';
        if (not ignore[static_cast<std::size_t>(i)])
          {
            auto const diff = eval.diffs[static_cast<std::size_t>(i)];
            if ((diff == 'A') or (diff == 'B'))
              {
                vote = (diff == m) ? '+' : '!';
              }
            else if ((diff == 'N') or (diff == '?'))
              {
                vote = '0';
              }
          }
        eval.votes[static_cast<std::size_t>(i)] = vote;
        if (vote == '!')
          {
            eval.diffs[static_cast<std::size_t>(i)] =
              static_cast<char>(std::tolower(eval.diffs[static_cast<std::size_t>(i)]));
          }
      }

    for (auto i = best_i + 1; i < alnlen; ++i)
      {
        auto const diff = eval.diffs[static_cast<std::size_t>(i)];
        if ((diff == ' ') or (diff == 'A'))
          {
            eval.model[static_cast<std::size_t>(i)] = 'x';
          }
        else
          {
            break;
          }
      }

    auto const & aln_a = best_is_reverse ? eval.b_aln : eval.a_aln;
    auto const & aln_b = best_is_reverse ? eval.a_aln : eval.b_aln;
    auto match_QA = 0;
    auto match_QB = 0;
    auto match_AB = 0;
    auto match_QM = 0;
    auto cols = 0;
    for (auto i = 0; i < alnlen; ++i)
      {
        if (ignore[static_cast<std::size_t>(i)])
          {
            continue;
          }
        ++cols;
        auto const qsym = map_4bit(eval.q_aln[static_cast<std::size_t>(i)]);
        auto const asym = map_4bit(aln_a[static_cast<std::size_t>(i)]);
        auto const bsym = map_4bit(aln_b[static_cast<std::size_t>(i)]);
        auto const msym = (i <= best_i) ? asym : bsym;
        if (qsym == asym)
          {
            ++match_QA;
          }
        if (qsym == bsym)
          {
            ++match_QB;
          }
        if (asym == bsym)
          {
            ++match_AB;
          }
        if (qsym == msym)
          {
            ++match_QM;
          }
      }

    if (cols <= 0)
      {
        return eval;
      }

    eval.valid = true;
    eval.reverse_orientation = best_is_reverse;
    eval.best_i = best_i;
    eval.best_h = best_h;
    eval.cols = cols;
    eval.match_QA = match_QA;
    eval.match_QB = match_QB;
    eval.match_AB = match_AB;
    eval.match_QM = match_QM;
    eval.QA = 100.0 * static_cast<double>(match_QA) / static_cast<double>(cols);
    eval.QB = 100.0 * static_cast<double>(match_QB) / static_cast<double>(cols);
    eval.AB = 100.0 * static_cast<double>(match_AB) / static_cast<double>(cols);
    eval.QT = std::max(eval.QA, eval.QB);
    eval.QM = 100.0 * static_cast<double>(match_QM) / static_cast<double>(cols);
    eval.divdiff = eval.QM - eval.QT;
    eval.left_y = best_left_y;
    eval.left_n = best_left_n;
    eval.left_a = best_left_a;
    eval.right_y = best_right_y;
    eval.right_n = best_right_n;
    eval.right_a = best_right_a;
    return eval;
  };

  for (std::size_t qi = 0; qi < records.size(); ++qi)
    {
      auto const & q = records[qi];
      auto best_one = static_cast<int>(q.left.size() + q.right.size());
      auto best_two = best_one;
      auto delta = 0;
      auto best_class = std::string{"NONE"};
      auto best_h_out = 0.0;
      auto best_left_y = 0;
      auto best_left_n = 0;
      auto best_left_a = 0;
      auto best_right_y = 0;
      auto best_right_n = 0;
      auto best_right_a = 0;
      auto best_ai = -1;
      auto best_bi = -1;
      TavPairEval best_eval;
      auto has_parent_pair = false;

      std::vector<unsigned int> parent_kmer_scores(records.size(), 0U);
      unsigned int qk_left = 0;
      unsigned int const * qk_left_list = nullptr;
      unique_count(parent_unique_handle,
                   chimera_wordlength,
                   static_cast<int>(q.left.size()),
                   q.left.c_str(),
                   &qk_left,
                   &qk_left_list,
                   opt_qmask);
      for (auto i = 0U; i < qk_left; ++i)
        {
          auto const left_it = parent_kmer_index.left_postings.find(qk_left_list[i]);
          if (left_it == parent_kmer_index.left_postings.end())
            {
              continue;
            }
          for (auto const parent_idx : left_it->second)
            {
              ++parent_kmer_scores[static_cast<std::size_t>(parent_idx)];
            }
        }

      unsigned int qk_right = 0;
      unsigned int const * qk_right_list = nullptr;
      unique_count(parent_unique_handle,
                   chimera_wordlength,
                   static_cast<int>(q.right.size()),
                   q.right.c_str(),
                   &qk_right,
                   &qk_right_list,
                   opt_qmask);
      for (auto i = 0U; i < qk_right; ++i)
        {
          auto const right_it = parent_kmer_index.right_postings.find(qk_right_list[i]);
          if (right_it == parent_kmer_index.right_postings.end())
            {
              continue;
            }
          for (auto const parent_idx : right_it->second)
            {
              ++parent_kmer_scores[static_cast<std::size_t>(parent_idx)];
            }
        }

      auto const qk_total = qk_left + qk_right;
      auto const minmatches = std::min(static_cast<unsigned int>(opt_minwordmatches), qk_total);
      std::vector<std::size_t> candidate_indices;
      auto * candidate_heap = minheap_init(chimera_tophits);
      for (auto const parent_idx : parent_pool_indices)
        {
          auto const & p = records[parent_idx];
          if (p.abundance < static_cast<int64_t>(opt_abskew * q.abundance))
            {
              continue;
            }
          auto const score = parent_kmer_scores[parent_idx];
          if (score >= minmatches)
            {
              elem_t novel;
              novel.count = score;
              novel.seqno = static_cast<unsigned int>(parent_idx);
              novel.length = static_cast<unsigned int>(p.left.size() + p.right.size());
              minheap_add(candidate_heap, &novel);
            }
        }
      minheap_sort(candidate_heap);
      while (!minheap_isempty(candidate_heap))
        {
          auto const e = minheap_poplast(candidate_heap);
          candidate_indices.push_back(static_cast<std::size_t>(e.seqno));
        }
      minheap_exit(candidate_heap);

      std::vector<TavParentCandidate> candidates;
      candidates.reserve(candidate_indices.size());
      for (auto const parent_idx : candidate_indices)
        {
          TavParentCandidate c;
          c.record_index = parent_idx;
          c.kmer_score = parent_kmer_scores[parent_idx];

          auto const & p = records[parent_idx];
          c.cigar_left = align_cigar_string(q.left, p.left);
          c.cigar_right = align_cigar_string(q.right, p.right);

          auto left_match = compute_match_vector(q.left, p.left, c.cigar_left);
          auto right_match = compute_match_vector(q.right, p.right, c.cigar_right);
          c.match_axis.reserve(left_match.size() + right_match.size());
          c.match_axis.insert(c.match_axis.end(), left_match.begin(), left_match.end());
          c.match_axis.insert(c.match_axis.end(), right_match.begin(), right_match.end());
          candidates.push_back(std::move(c));
        }

      auto const best_parent_pair = select_best_parent_pair(candidates, q.left.size() + q.right.size());
      if ((best_parent_pair.first >= 0) and (best_parent_pair.second >= 0))
        {
          auto const & cand_a = candidates[static_cast<std::size_t>(best_parent_pair.first)];
          auto const & cand_b = candidates[static_cast<std::size_t>(best_parent_pair.second)];
          auto const & parent_a = records[cand_a.record_index];
          auto const & parent_b = records[cand_b.record_index];
          best_eval = evaluate_parent_pair(q,
                                           parent_a,
                                           parent_b,
                                           cand_a.cigar_left,
                                           cand_a.cigar_right,
                                           cand_b.cigar_left,
                                           cand_b.cigar_right);
          if (best_eval.valid)
            {
              has_parent_pair = true;
              best_ai = static_cast<int>(cand_a.record_index);
              best_bi = static_cast<int>(cand_b.record_index);
              if (best_eval.reverse_orientation)
                {
                  std::swap(best_ai, best_bi);
                }
              best_h_out = std::max(0.0, best_eval.best_h);
              best_left_y = best_eval.left_y;
              best_left_n = best_eval.left_n;
              best_left_a = best_eval.left_a;
              best_right_y = best_eval.right_y;
              best_right_n = best_eval.right_n;
              best_right_a = best_eval.right_a;
              best_one = std::max(best_eval.match_QA, best_eval.match_QB);
              best_two = best_eval.match_QM;
              delta = best_two - best_one;
              if (best_eval.best_i < 0)
                {
                  best_class = "NONE";
                }
              else if (best_eval.left_aln_len == 0U)
                {
                  best_class = "RIGHT_BREAK";
                }
              else
                {
                  auto const boundary = static_cast<int>(best_eval.left_aln_len) - 1;
                  if (best_eval.best_i < boundary)
                    {
                      best_class = "LEFT_BREAK";
                    }
                  else if (best_eval.best_i == boundary)
                    {
                      best_class = "MIDDLE_BREAK";
                    }
                  else
                    {
                      best_class = "RIGHT_BREAK";
                    }
                }
            }
        }

      auto best_a = q.tav_id;
      auto best_b = q.tav_id;
      auto is_chim = false;
      if (has_parent_pair and (best_ai >= 0) and (best_bi >= 0))
        {
          best_a = records[static_cast<std::size_t>(best_ai)].tav_id;
          best_b = records[static_cast<std::size_t>(best_bi)].tav_id;
          is_chim = (best_eval.match_QM == best_eval.cols) and (best_eval.QT < 100.0);
        }

      if (is_chim)
        {
          chim.push_back(q);
          chim_scores.push_back(best_h_out);
        }
      else
        {
          nonchim.push_back(q);
          nonchim_scores.push_back(best_h_out);
          parent_pool_indices.push_back(qi);
          add_parent_to_index(qi);
        }

      if (fp_report != nullptr)
        {
          std::fprintf(fp_report,
                       "%s\t%s\t%s\t%s\t%d\t%d\t%d\n",
                       q.tav_id.c_str(),
                       best_a.c_str(),
                       best_b.c_str(),
                       best_class.c_str(),
                       best_one,
                       best_two,
                       delta);
        }

      if (fp_uchimeout_paired != nullptr)
        {
          if ((not has_parent_pair) or (best_ai < 0) or (best_bi < 0) or (not best_eval.valid))
            {
              write_uchimeout_no_parents(q, best_h_out);
            }
          else
            {
              auto const & parent_a = records[static_cast<std::size_t>(best_ai)];
              auto const & parent_b = records[static_cast<std::size_t>(best_bi)];

              std::fprintf(fp_uchimeout_paired, "%.4f\t", best_h_out);
              print_stripped_header(fp_uchimeout_paired, q.header);
              std::fprintf(fp_uchimeout_paired, "\t");
              print_stripped_header(fp_uchimeout_paired, parent_a.header);
              std::fprintf(fp_uchimeout_paired, "\t");
              print_stripped_header(fp_uchimeout_paired, parent_b.header);
              std::fprintf(fp_uchimeout_paired, "\t");
              if (opt_uchimeout5 == 0)
                {
                  if (best_eval.QA >= best_eval.QB)
                    {
                      print_stripped_header(fp_uchimeout_paired, parent_a.header);
                    }
                  else
                    {
                      print_stripped_header(fp_uchimeout_paired, parent_b.header);
                    }
                  std::fprintf(fp_uchimeout_paired, "\t");
                }

              std::fprintf(fp_uchimeout_paired,
                           "%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t"
                           "%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%c\n",
                           best_eval.QM,
                           best_eval.QA,
                           best_eval.QB,
                           best_eval.AB,
                           best_eval.QT,
                           best_left_y,
                           best_left_n,
                           best_left_a,
                           best_right_y,
                           best_right_n,
                           best_right_a,
                           best_eval.divdiff,
                           is_chim ? 'Y' : 'N');
            }
        }

      if ((fp_uchimealns_paired != nullptr) and is_chim and (best_ai >= 0) and (best_bi >= 0) and best_eval.valid)
        {
          auto const & parent_a = records[static_cast<std::size_t>(best_ai)];
          auto const & parent_b = records[static_cast<std::size_t>(best_bi)];

          auto write_wrapped = [&](char const * label, std::string const & seq) -> void {
            auto width = opt_alignwidth;
            if (width <= 0)
              {
                width = static_cast<int>(seq.size());
              }
            if (width <= 0)
              {
                std::fprintf(fp_uchimealns_paired, "%-8s\n", label);
                return;
              }
            bool first = true;
            for (std::size_t start = 0; start < seq.size(); start += static_cast<std::size_t>(width))
              {
                auto const chunk = std::min(static_cast<std::size_t>(width), seq.size() - start);
                std::fprintf(fp_uchimealns_paired,
                             "%-8s %.*s\n",
                             first ? label : "",
                             static_cast<int>(chunk),
                             seq.data() + start);
                first = false;
              }
          };

          std::fprintf(fp_uchimealns_paired, "\n");
          std::fprintf(fp_uchimealns_paired,
                       "------------------------------------------------------------------------\n");
          std::fprintf(fp_uchimealns_paired, "Query   (%5zu nt) ", q.left.size() + q.right.size());
          print_stripped_header(fp_uchimealns_paired, q.header);
          std::fprintf(fp_uchimealns_paired, "\nParentA (%5zu nt) ", parent_a.left.size() + parent_a.right.size());
          print_stripped_header(fp_uchimealns_paired, parent_a.header);
          std::fprintf(fp_uchimealns_paired, "\nParentB (%5zu nt) ", parent_b.left.size() + parent_b.right.size());
          print_stripped_header(fp_uchimealns_paired, parent_b.header);
          std::fprintf(fp_uchimealns_paired, "\nClass %s, Score %.4f\n\n", best_class.c_str(), best_h_out);

          auto const left_q = best_eval.q_aln.substr(0U, best_eval.left_aln_len);
          auto const left_a = best_eval.a_aln.substr(0U, best_eval.left_aln_len);
          auto const left_b = best_eval.b_aln.substr(0U, best_eval.left_aln_len);
          auto const left_d = best_eval.diffs.substr(0U, best_eval.left_aln_len);
          auto const left_v = best_eval.votes.substr(0U, best_eval.left_aln_len);
          auto const left_m = best_eval.model.substr(0U, best_eval.left_aln_len);
          auto const right_q = best_eval.q_aln.substr(best_eval.left_aln_len);
          auto const right_a = best_eval.a_aln.substr(best_eval.left_aln_len);
          auto const right_b = best_eval.b_aln.substr(best_eval.left_aln_len);
          auto const right_d = best_eval.diffs.substr(best_eval.left_aln_len);
          auto const right_v = best_eval.votes.substr(best_eval.left_aln_len);
          auto const right_m = best_eval.model.substr(best_eval.left_aln_len);

          write_wrapped("A_LEFT", left_a);
          write_wrapped("Q_LEFT", left_q);
          write_wrapped("B_LEFT", left_b);
          write_wrapped("DIFF_L", left_d);
          write_wrapped("VOTE_L", left_v);
          write_wrapped("MODL_L", left_m);
          std::fprintf(fp_uchimealns_paired, "\n");
          write_wrapped("A_RIGHT", right_a);
          write_wrapped("Q_RIGHT", right_q);
          write_wrapped("B_RIGHT", right_b);
          write_wrapped("DIFF_R", right_d);
          write_wrapped("VOTE_R", right_v);
          write_wrapped("MODL_R", right_m);
          std::fprintf(fp_uchimealns_paired, "\n");

          std::fprintf(fp_uchimealns_paired,
                       "Ids. QA %.1f%%, QB %.1f%%, AB %.1f%%, QModel %.1f%%, Div. %+.1f%%\n",
                       best_eval.QA,
                       best_eval.QB,
                       best_eval.AB,
                       best_eval.QM,
                       best_eval.divdiff);
          std::fprintf(fp_uchimealns_paired,
                       "Diffs Left: N %d, A %d, Y %d; Right: N %d, A %d, Y %d\n",
                       best_left_n,
                       best_left_a,
                       best_left_y,
                       best_right_n,
                       best_right_a,
                       best_right_y);
        }
    }

  if (fp_report != nullptr)
    {
      std::fclose(fp_report);
    }
  if (fp_uchimeout_paired != nullptr)
    {
      std::fclose(fp_uchimeout_paired);
    }
  if (fp_uchimealns_paired != nullptr)
    {
      std::fclose(fp_uchimealns_paired);
    }
  if (fp_borderline_paired != nullptr)
    {
      std::fclose(fp_borderline_paired);
    }
  unique_exit(parent_unique_handle);

  if (opt_nonchimeras != nullptr)
    {
      auto * fp = fopen_output(opt_nonchimeras);
      if (fp == nullptr)
        {
          fatal("Unable to open nonchimeras output file for writing");
        }
      write_paired_fasta_interleaved(fp, nonchim, nonchim_scores, "uchime_denovo");
      std::fclose(fp);
    }
  if (opt_chimeras != nullptr)
    {
      auto * fp = fopen_output(opt_chimeras);
      if (fp == nullptr)
        {
          fatal("Unable to open chimeras output file for writing");
        }
      write_paired_fasta_interleaved(fp, chim, chim_scores, "uchime_denovo");
      std::fclose(fp);
    }
  if (opt_nonchimeras_tsv != nullptr)
    {
      auto * fp = fopen_output(opt_nonchimeras_tsv);
      if (fp == nullptr)
        {
          fatal("Unable to open nonchimeras TSV output file for writing");
        }
      write_catalog(fp, nonchim);
      std::fclose(fp);
    }
  if (opt_chimeras_tsv != nullptr)
    {
      auto * fp = fopen_output(opt_chimeras_tsv);
      if (fp == nullptr)
        {
          fatal("Unable to open chimeras TSV output file for writing");
        }
      write_catalog(fp, chim);
      std::fclose(fp);
    }

  if (!opt_quiet)
    {
      std::fprintf(stderr,
                   "TAV uchime3_denovo: %zu input pairs (%" PRId64 " filtered short, %" PRId64 " filtered long, filter=%s) -> %zu nonchimeric, %zu chimeric\n",
                   records.size() + static_cast<std::size_t>(load_result.filtered_short + load_result.filtered_long),
                   load_result.filtered_short,
                   load_result.filtered_long,
                   (opt_filter == 0) ? "any" : "both",
                   nonchim.size(),
                   chim.size());
    }
}


auto tav_usearch_global(struct Parameters const & parameters, char * cmdline, char * progheader) -> void
{
  if (parameters.opt_reverse == nullptr)
    {
      fatal("TAV usearch_global requires paired query reads with --reverse");
    }
  if (parameters.opt_db == nullptr)
    {
      fatal("Database filename not specified with --db");
    }

  auto mask_sequence = [](std::string & seq, int const mask_mode) -> void {
    if (seq.empty())
      {
        return;
      }
    if ((mask_mode != MASK_DUST) and (not ((mask_mode == MASK_SOFT) and (opt_hardmask != 0))))
      {
        return;
      }

    std::vector<char> buf(seq.size() + 1U);
    std::copy(seq.begin(), seq.end(), buf.begin());
    buf[seq.size()] = '\0';

    if (mask_mode == MASK_DUST)
      {
        dust(buf.data(), static_cast<int>(seq.size()));
      }
    else
      {
        hardmask(buf.data(), static_cast<int>(seq.size()));
      }

    seq.assign(buf.data(), seq.size());
  };

  auto trim_pair_suffix = [](std::string const & header) -> std::string {
    auto id = trim_header_to_id(header.c_str());
    if ((id.size() > 2U) and (id[id.size() - 2U] == '/') and ((id.back() == '1') or (id.back() == '2')))
      {
        id.resize(id.size() - 2U);
      }
    return id;
  };

  auto append_db_record = [&](std::vector<TavRecord> & out,
                              std::string const & left_header,
                              std::string const & right_header,
                              std::string const & left_seq,
                              std::string const & right_seq,
                              int64_t const left_abundance,
                              int64_t const right_abundance) -> void {
    auto const left_len = static_cast<int64_t>(left_seq.size());
    auto const right_len = static_cast<int64_t>(right_seq.size());
    auto const anchor_len = get_anchor_len(parameters, left_len, right_len);
    if (anchor_len <= 0)
      {
        return;
      }

    auto const left_short = left_len < opt_minseqlength;
    auto const right_short = right_len < opt_minseqlength;
    auto const left_long = left_len > opt_maxseqlength;
    auto const right_long = right_len > opt_maxseqlength;
    auto const drop_short = (opt_filter == 0) ? (left_short or right_short) : (left_short and right_short);
    auto const drop_long = (opt_filter == 0) ? (left_long or right_long) : (left_long and right_long);
    if (drop_short or drop_long)
      {
        return;
      }

    TavRecord r;
    r.left = left_seq.substr(0U, static_cast<std::size_t>(anchor_len));
    r.right = right_seq.substr(0U, static_cast<std::size_t>(anchor_len));
    mask_sequence(r.left, opt_dbmask);
    mask_sequence(r.right, opt_dbmask);
    r.abundance = left_abundance;
    r.header = trim_pair_suffix(left_header);
    auto const right_base = trim_pair_suffix(right_header);
    if ((right_base != r.header) and (not opt_quiet))
      {
        std::fprintf(stderr,
                     "Warning: paired database headers differ (%s vs %s); using left header base\n",
                     left_header.c_str(),
                     right_header.c_str());
      }
    if ((right_abundance != left_abundance) and (not opt_quiet))
      {
        std::fprintf(stderr,
                     "Warning: paired database abundances differ (%" PRId64 " vs %" PRId64 "); using left abundance\n",
                     left_abundance,
                     right_abundance);
      }
    r.tav_id = r.header;
    r.first_seen = static_cast<int64_t>(out.size());
    out.push_back(std::move(r));
  };

  auto load_paired_db_from_interleaved_fastx = [&](char const * filename) -> std::vector<TavRecord> {
    std::vector<TavRecord> out;
    auto * db_h = fastx_open(filename);
    if (db_h == nullptr)
      {
        fatal("Unrecognized file type for paired database (not proper FASTA or FASTQ format): %s", filename);
      }

    while (fastx_next(db_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
      {
        auto const left_header = std::string{fastx_get_header(db_h)};
        auto const left_len = static_cast<int64_t>(fastx_get_sequence_length(db_h));
        auto const left_abundance = fastx_get_abundance(db_h);
        auto const left = std::string{fastx_get_sequence(db_h), static_cast<std::size_t>(left_len)};

        if (!fastx_next(db_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
          {
            fatal("Odd number of records in paired FASTX database %s; expected interleaved left/right entries", filename);
          }

        auto const right_header = std::string{fastx_get_header(db_h)};
        auto const right_len = static_cast<int64_t>(fastx_get_sequence_length(db_h));
        auto const right_abundance = fastx_get_abundance(db_h);
        auto const right = std::string{fastx_get_sequence(db_h), static_cast<std::size_t>(right_len)};

        append_db_record(out,
                         left_header,
                         right_header,
                         left,
                         right,
                         left_abundance,
                         right_abundance);
      }

    fastx_close(db_h);
    return out;
  };

  auto load_paired_db_from_split_fastx = [&](char const * left_filename,
                                              char const * right_filename) -> std::vector<TavRecord> {
    std::vector<TavRecord> out;

    auto * db_left_h = fastx_open(left_filename);
    if (db_left_h == nullptr)
      {
        fatal("Unrecognized file type for paired left database (not proper FASTA or FASTQ format): %s", left_filename);
      }
    auto * db_right_h = fastx_open(right_filename);
    if (db_right_h == nullptr)
      {
        fastx_close(db_left_h);
        fatal("Unrecognized file type for paired right database (not proper FASTA or FASTQ format): %s", right_filename);
      }

    while (fastx_next(db_left_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
      {
        if (!fastx_next(db_right_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
          {
            fastx_close(db_right_h);
            fastx_close(db_left_h);
            fatal("More left database records than right database records");
          }

        auto const left_header = std::string{fastx_get_header(db_left_h)};
        auto const right_header = std::string{fastx_get_header(db_right_h)};
        auto const left_len = static_cast<int64_t>(fastx_get_sequence_length(db_left_h));
        auto const right_len = static_cast<int64_t>(fastx_get_sequence_length(db_right_h));
        auto const left_abundance = fastx_get_abundance(db_left_h);
        auto const right_abundance = fastx_get_abundance(db_right_h);
        auto const left = std::string{fastx_get_sequence(db_left_h), static_cast<std::size_t>(left_len)};
        auto const right = std::string{fastx_get_sequence(db_right_h), static_cast<std::size_t>(right_len)};

        append_db_record(out,
                         left_header,
                         right_header,
                         left,
                         right,
                         left_abundance,
                         right_abundance);
      }

    if (fastx_next(db_right_h, not opt_notrunclabels, chrmap_no_change_vector.data()))
      {
        fastx_close(db_right_h);
        fastx_close(db_left_h);
        fatal("More right database records than left database records");
      }

    fastx_close(db_right_h);
    fastx_close(db_left_h);
    return out;
  };

  if (tav_is_catalog_file(parameters.opt_db) or
      ((parameters.opt_db2 != nullptr) and tav_is_catalog_file(parameters.opt_db2)))
    {
      fatal("Paired usearch_global catalog input is not supported; provide paired FASTA/FASTQ database records");
    }

  auto db = (parameters.opt_db2 != nullptr)
    ? load_paired_db_from_split_fastx(parameters.opt_db, parameters.opt_db2)
    : load_paired_db_from_interleaved_fastx(parameters.opt_db);

  if (db.empty())
    {
      if (parameters.opt_db2 != nullptr)
        {
          fatal("No paired database records loaded from the provided --db/--db2 files");
        }
      fatal("No paired database records loaded from %s", parameters.opt_db);
    }

  for (auto & r : db)
    {
      if (r.header.empty())
        {
          r.header = r.tav_id;
        }
      if (r.tav_id.empty())
        {
          r.tav_id = trim_header_to_id(r.header.c_str());
        }
      if (r.abundance < 1)
        {
          r.abundance = 1;
        }
    }

  if (db.empty())
    {
      fatal("Paired usearch_global database is empty");
    }

  auto write_fasta_record = [](std::FILE * fp,
                               std::string const & header,
                               std::string const & sequence,
                               int64_t const abundance,
                               int & ordinal) -> void {
    if (fp == nullptr)
      {
        return;
      }

    auto const bounded_abundance = std::max<int64_t>(abundance, 1);
    auto const fasta_abundance =
      static_cast<unsigned int>(std::min<int64_t>(bounded_abundance, std::numeric_limits<unsigned int>::max()));

    fasta_print_general(fp,
                        nullptr,
                        sequence.c_str(),
                        static_cast<int>(sequence.size()),
                        header.c_str(),
                        static_cast<int>(header.size()),
                        fasta_abundance,
                        ordinal,
                        -1.0,
                        -1,
                        -1,
                        nullptr,
                        0.0);
    ++ordinal;
  };

  auto write_query_pair_interleaved = [&](std::FILE * fp,
                                          std::string const & header_base,
                                          std::string const & left,
                                          std::string const & right,
                                          int64_t const abundance,
                                          int & ordinal) -> void {
    auto const left_header = header_base + "/1";
    auto const right_header = header_base + "/2";
    write_fasta_record(fp, left_header, left, abundance, ordinal);
    write_fasta_record(fp, right_header, right, abundance, ordinal);
  };

  auto write_query_pair_split = [&](std::FILE * fp_left,
                                    std::FILE * fp_right,
                                    std::string const & header_base,
                                    std::string const & left,
                                    std::string const & right,
                                    int64_t const abundance,
                                    int & left_ordinal,
                                    int & right_ordinal) -> void {
    auto const left_header = header_base + "/1";
    auto const right_header = header_base + "/2";
    write_fasta_record(fp_left, left_header, left, abundance, left_ordinal);
    write_fasta_record(fp_right, right_header, right, abundance, right_ordinal);
  };

  auto write_db_pair_interleaved = [&](std::FILE * fp,
                                       TavRecord const & rec,
                                       int64_t const abundance,
                                       int & ordinal) -> void {
    write_query_pair_interleaved(fp, rec.header, rec.left, rec.right, abundance, ordinal);
  };

  auto write_db_pair_split = [&](std::FILE * fp_left,
                                 std::FILE * fp_right,
                                 TavRecord const & rec,
                                 int64_t const abundance,
                                 int & left_ordinal,
                                 int & right_ordinal) -> void {
    write_query_pair_split(fp_left, fp_right, rec.header, rec.left, rec.right, abundance, left_ordinal, right_ordinal);
  };

  std::FILE * fp_alnout = nullptr;
  if (opt_alnout != nullptr)
    {
      fp_alnout = fopen_output(opt_alnout);
      if (fp_alnout == nullptr)
        {
          fatal("Unable to open alnout output file for writing");
        }
      std::fprintf(fp_alnout, "%s\n", cmdline);
      std::fprintf(fp_alnout, "%s\n", progheader);
    }

  std::FILE * fp_userout = nullptr;
  if (opt_userout != nullptr)
    {
      fp_userout = fopen_output(opt_userout);
      if (fp_userout == nullptr)
        {
          fatal("Unable to open paired userout output file for writing");
        }
    }

  std::FILE * fp_blast6 = nullptr;
  if (opt_blast6out != nullptr)
    {
      fp_blast6 = fopen_output(opt_blast6out);
      if (fp_blast6 == nullptr)
        {
          fatal("Unable to open blast6-like output file for writing");
        }
    }

  std::FILE * fp_uc = nullptr;
  if (parameters.opt_uc != nullptr)
    {
      fp_uc = fopen_output(parameters.opt_uc);
      if (fp_uc == nullptr)
        {
          fatal("Unable to open uc output file for writing");
        }
    }

  std::FILE * fp_matched = nullptr;
  if (opt_matched != nullptr)
    {
      fp_matched = fopen_output(opt_matched);
      if (fp_matched == nullptr)
        {
          fatal("Unable to open matched output file for writing");
        }
    }

  std::FILE * fp_matched2 = nullptr;
  if (opt_matched2 != nullptr)
    {
      fp_matched2 = fopen_output(opt_matched2);
      if (fp_matched2 == nullptr)
        {
          fatal("Unable to open matched2 output file for writing");
        }
    }

  std::FILE * fp_notmatched = nullptr;
  if (opt_notmatched != nullptr)
    {
      fp_notmatched = fopen_output(opt_notmatched);
      if (fp_notmatched == nullptr)
        {
          fatal("Unable to open notmatched output file for writing");
        }
    }

  std::FILE * fp_notmatched2 = nullptr;
  if (opt_notmatched2 != nullptr)
    {
      fp_notmatched2 = fopen_output(opt_notmatched2);
      if (fp_notmatched2 == nullptr)
        {
          fatal("Unable to open notmatched2 output file for writing");
        }
    }

  std::FILE * fp_dbmatched = nullptr;
  if (parameters.opt_dbmatched != nullptr)
    {
      fp_dbmatched = fopen_output(parameters.opt_dbmatched);
      if (fp_dbmatched == nullptr)
        {
          fatal("Unable to open dbmatched output file for writing");
        }
    }

  std::FILE * fp_dbmatched2 = nullptr;
  if (parameters.opt_dbmatched2 != nullptr)
    {
      fp_dbmatched2 = fopen_output(parameters.opt_dbmatched2);
      if (fp_dbmatched2 == nullptr)
        {
          fatal("Unable to open dbmatched2 output file for writing");
        }
    }

  std::FILE * fp_dbnotmatched = nullptr;
  if (parameters.opt_dbnotmatched != nullptr)
    {
      fp_dbnotmatched = fopen_output(parameters.opt_dbnotmatched);
      if (fp_dbnotmatched == nullptr)
        {
          fatal("Unable to open dbnotmatched output file for writing");
        }
    }

  std::FILE * fp_dbnotmatched2 = nullptr;
  if (parameters.opt_dbnotmatched2 != nullptr)
    {
      fp_dbnotmatched2 = fopen_output(parameters.opt_dbnotmatched2);
      if (fp_dbnotmatched2 == nullptr)
        {
          fatal("Unable to open dbnotmatched2 output file for writing");
        }
    }

  auto const any_otu_table_output = (opt_otutabout != nullptr) or
                                    (opt_biomout != nullptr) or
                                    (opt_mothur_shared_out != nullptr);
  if (any_otu_table_output)
    {
      otutable_init();
    }

  TavKmerIndex db_kmer_index;
  auto * kmer_unique_handle = unique_init();
  auto const wordlength = static_cast<int>(opt_wordlength);
  for (std::size_t i = 0; i < db.size(); ++i)
    {
      unsigned int left_unique_count = 0;
      unsigned int const * left_unique_list = nullptr;
      unique_count(kmer_unique_handle,
                   wordlength,
                   static_cast<int>(db[i].left.size()),
                   db[i].left.c_str(),
                   &left_unique_count,
                   &left_unique_list,
                   opt_dbmask);
      for (auto k = 0U; k < left_unique_count; ++k)
        {
          db_kmer_index.left_postings[left_unique_list[k]].push_back(static_cast<int>(i));
        }

      unsigned int right_unique_count = 0;
      unsigned int const * right_unique_list = nullptr;
      unique_count(kmer_unique_handle,
                   wordlength,
                   static_cast<int>(db[i].right.size()),
                   db[i].right.c_str(),
                   &right_unique_count,
                   &right_unique_list,
                   opt_dbmask);
      for (auto k = 0U; k < right_unique_count; ++k)
        {
          db_kmer_index.right_postings[right_unique_list[k]].push_back(static_cast<int>(i));
        }
    }

  struct Scoring scoring;
  scoring.match = opt_match;
  scoring.mismatch = opt_mismatch;
  scoring.gap_open_query_interior = opt_gap_open_query_interior;
  scoring.gap_extension_query_interior = opt_gap_extension_query_interior;
  scoring.gap_open_query_left = opt_gap_open_query_left;
  scoring.gap_open_target_left = opt_gap_open_target_left;
  scoring.gap_open_target_interior = opt_gap_open_target_interior;
  scoring.gap_open_query_right = opt_gap_open_query_right;
  scoring.gap_open_target_right = opt_gap_open_target_right;
  scoring.gap_extension_query_left = opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opt_gap_extension_target_left;
  scoring.gap_extension_target_interior = opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opt_gap_extension_target_right;
  LinearMemoryAligner lma(scoring);

  std::vector<uint64_t> dbmatched_abundance(db.size(), 0U);

  auto * qf = fastx_open(parameters.opt_usearch_global);
  auto * qr = fastx_open(parameters.opt_reverse);
  if (qf == nullptr)
    {
      fatal("Unrecognized query file type (not proper FASTA or FASTQ format): %s", parameters.opt_usearch_global);
    }
  if (qr == nullptr)
    {
      fatal("Unrecognized reverse query file type (not proper FASTA or FASTQ format): %s", parameters.opt_reverse);
    }

  int64_t queries = 0;
  uint64_t queries_abundance = 0;
  int64_t assigned = 0;
  uint64_t assigned_abundance = 0;
  int64_t filtered_short = 0;
  int64_t filtered_long = 0;

  int matched_ordinal = 1;
  int matched2_ordinal = 1;
  int notmatched_ordinal = 1;
  int notmatched2_ordinal = 1;
  int dbmatched_ordinal = 1;
  int dbmatched2_ordinal = 1;
  int dbnotmatched_ordinal = 1;
  int dbnotmatched2_ordinal = 1;

  while (fastx_next(qf, not opt_notrunclabels, chrmap_no_change_vector.data()))
    {
      if (!fastx_next(qr, not opt_notrunclabels, chrmap_no_change_vector.data()))
        {
          fatal("More forward queries than reverse queries");
        }

      ++queries;
      auto const query_abundance = std::max<int64_t>(static_cast<int64_t>(fastx_get_abundance(qf)), 1);
      queries_abundance += static_cast<uint64_t>(query_abundance);

      auto const qf_len = static_cast<int64_t>(fastx_get_sequence_length(qf));
      auto const qr_len = static_cast<int64_t>(fastx_get_sequence_length(qr));
      auto const anchor_len = std::max<int64_t>(0, get_anchor_len(parameters, qf_len, qr_len));

      auto const query_header = std::string{fastx_get_header(qf)};
      std::string q_left(fastx_get_sequence(qf), static_cast<std::size_t>(anchor_len));
      std::vector<char> q_rc(static_cast<std::size_t>(qr_len) + 1U);
      reverse_complement(q_rc.data(), fastx_get_sequence(qr), qr_len);
      std::string q_right(q_rc.data(), static_cast<std::size_t>(anchor_len));
      mask_sequence(q_left, parameters.opt_qmask);
      mask_sequence(q_right, parameters.opt_qmask);

      TavRecord qrec;
      qrec.header = query_header;
      qrec.tav_id = trim_header_to_id(query_header.c_str());
      qrec.left = q_left;
      qrec.right = q_right;
      qrec.abundance = query_abundance;

      std::vector<TavHit> hits;

      std::vector<unsigned int> db_kmer_scores(db.size(), 0U);
      unsigned int qk_left = 0;
      unsigned int const * qk_left_list = nullptr;
      unique_count(kmer_unique_handle,
                   wordlength,
                   static_cast<int>(qrec.left.size()),
                   qrec.left.c_str(),
                   &qk_left,
                   &qk_left_list,
                   parameters.opt_qmask);
      for (auto i = 0U; i < qk_left; ++i)
        {
          auto const left_it = db_kmer_index.left_postings.find(qk_left_list[i]);
          if (left_it == db_kmer_index.left_postings.end())
            {
              continue;
            }
          for (auto const db_index : left_it->second)
            {
              ++db_kmer_scores[static_cast<std::size_t>(db_index)];
            }
        }

      unsigned int qk_right = 0;
      unsigned int const * qk_right_list = nullptr;
      unique_count(kmer_unique_handle,
                   wordlength,
                   static_cast<int>(qrec.right.size()),
                   qrec.right.c_str(),
                   &qk_right,
                   &qk_right_list,
                   parameters.opt_qmask);
      for (auto i = 0U; i < qk_right; ++i)
        {
          auto const right_it = db_kmer_index.right_postings.find(qk_right_list[i]);
          if (right_it == db_kmer_index.right_postings.end())
            {
              continue;
            }
          for (auto const db_index : right_it->second)
            {
              ++db_kmer_scores[static_cast<std::size_t>(db_index)];
            }
        }

      auto const qk_total = qk_left + qk_right;
      auto const minmatches = std::min(static_cast<unsigned int>(opt_minwordmatches), qk_total);
      auto const maxaccepts_effective =
        (opt_maxaccepts == 0)
        ? static_cast<int64_t>(db.size())
        : std::min<int64_t>(opt_maxaccepts, static_cast<int64_t>(db.size()));
      auto const maxrejects_effective =
        (opt_maxrejects == 0)
        ? static_cast<int64_t>(db.size())
        : std::min<int64_t>(opt_maxrejects, static_cast<int64_t>(db.size()));
      auto const maxhits_considered = std::max<int64_t>(1, maxrejects_effective + maxaccepts_effective - 1);
      auto const tophits = std::max<int64_t>(
        1,
        std::min<int64_t>(static_cast<int64_t>(db.size()),
                          maxrejects_effective + maxaccepts_effective + static_cast<int64_t>(MAXDELAYED)));

      std::vector<std::pair<int, unsigned int>> candidates;
      auto * candidate_heap = minheap_init(static_cast<int>(tophits));
      for (std::size_t i = 0; i < db.size(); ++i)
        {
          auto const score = db_kmer_scores[i];
          if (score >= minmatches)
            {
              elem_t novel;
              novel.count = score;
              novel.seqno = static_cast<unsigned int>(i);
              novel.length = static_cast<unsigned int>(db[i].left.size() + db[i].right.size());
              minheap_add(candidate_heap, &novel);
            }
        }
      minheap_sort(candidate_heap);
      candidates.reserve(static_cast<std::size_t>(candidate_heap->count));
      while (!minheap_isempty(candidate_heap))
        {
          auto const e = minheap_poplast(candidate_heap);
          candidates.emplace_back(static_cast<int>(e.seqno), e.count);
        }
      minheap_exit(candidate_heap);

      auto accepts = 0;
      auto rejects = 0;
      std::vector<std::size_t> pending;
      pending.reserve(MAXDELAYED);

      auto process_pending = [&]() -> void {
        for (auto const hit_index : pending)
          {
            if ((accepts >= maxaccepts_effective) or (rejects >= maxrejects_effective))
              {
                break;
              }

            auto & hit = hits[hit_index];
            if (hit.rejected)
              {
                ++rejects;
                continue;
              }

            auto const & candidate = db[static_cast<std::size_t>(hit.centroid_index)];

            if (!paired_aligned_filters_pass(qrec, candidate, hit, lma))
              {
                hit.rejected = true;
                hit.weak = false;
                ++rejects;
                continue;
              }

            if (hit.id >= (100.0 * opt_id))
              {
                hit.accepted = true;
                hit.weak = false;
                ++accepts;
              }
            else
              {
                hit.rejected = true;
                hit.weak = true;
                ++rejects;
              }
          }
        pending.clear();
      };

      for (auto const & c : candidates)
        {
          if ((accepts >= maxaccepts_effective) or
              (rejects >= maxrejects_effective) or
              (static_cast<int64_t>(hits.size()) >= maxhits_considered))
            {
              break;
            }

          TavHit hit;
          hit.centroid_index = c.first;
          hit.kmer_score = c.second;
          hits.push_back(hit);

          if (paired_unaligned_filters_pass(qrec, db[static_cast<std::size_t>(c.first)]))
            {
              // Defer alignment to keep stock-like delayed batching.
            }
          else
            {
              hits.back().rejected = true;
              hits.back().weak = false;
            }

          pending.push_back(hits.size() - 1U);
          if (pending.size() == MAXDELAYED)
            {
              process_pending();
            }
        }

      if (!pending.empty())
        {
          process_pending();
        }

      std::vector<TavHit const *> report_hits;
      report_hits.reserve(hits.size());
      for (auto const & hit : hits)
        {
          if (hit.accepted or hit.weak)
            {
              report_hits.push_back(&hit);
            }
        }

      std::sort(report_hits.begin(),
                report_hits.end(),
                [](TavHit const * lhs, TavHit const * rhs) -> bool {
                  if (lhs->rejected != rhs->rejected)
                    {
                      return lhs->rejected < rhs->rejected;
                    }
                  if (lhs->aligned != rhs->aligned)
                    {
                      return lhs->aligned > rhs->aligned;
                    }
                  if (lhs->aligned and (lhs->id != rhs->id))
                    {
                      return lhs->id > rhs->id;
                    }
                  return lhs->centroid_index < rhs->centroid_index;
                });

      auto const query_len = static_cast<int64_t>(qrec.left.size() + qrec.right.size());
      auto const toreport = std::min<int64_t>(opt_maxhits, static_cast<int64_t>(report_hits.size()));

      if (!report_hits.empty())
        {
          ++assigned;
          assigned_abundance += static_cast<uint64_t>(query_abundance);

          auto const dbmatch_weight = static_cast<uint64_t>(opt_sizein ? query_abundance : 1);
          for (auto const * hp : report_hits)
            {
              dbmatched_abundance[static_cast<std::size_t>(hp->centroid_index)] += dbmatch_weight;
            }

          auto const * best_hit = report_hits.front();
          auto const & best_target = db[static_cast<std::size_t>(best_hit->centroid_index)];
          if (any_otu_table_output)
            {
              otutable_add(query_header.c_str(), best_target.header.c_str(), qrec.abundance);
            }

          auto const top_hit_id = best_hit->id;
          for (int64_t t = 0; t < toreport; ++t)
            {
              auto const * hp = report_hits[static_cast<std::size_t>(t)];
              if ((opt_top_hits_only != 0) and (hp->id < top_hit_id))
                {
                  break;
                }

              auto const & target = db[static_cast<std::size_t>(hp->centroid_index)];
              auto const target_len = static_cast<int64_t>(target.left.size() + target.right.size());
              auto const d_left = hp->left.mismatches + hp->left.internal_indels;
              auto const d_right = hp->right.mismatches + hp->right.internal_indels;
              auto const d_total = d_left + d_right;

              if (fp_userout != nullptr)
                {
                  std::fprintf(fp_userout,
                               "%s\t%s\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\n",
                               query_header.c_str(),
                               target.header.c_str(),
                               hp->left.id / 100.0,
                               hp->right.id / 100.0,
                               hp->id / 100.0,
                               d_left,
                               d_right,
                               d_total);
                }

              if (fp_blast6 != nullptr)
                {
                  std::fprintf(fp_blast6,
                               "%s\t%s\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%" PRId64 "\t%d\t%d\n",
                               query_header.c_str(),
                               target.header.c_str(),
                               hp->id,
                               hp->internal_alignment_cols_total,
                               hp->mismatches_total,
                               hp->internal_gaps_total,
                               1,
                               static_cast<int>(query_len),
                               1,
                               target_len,
                               -1,
                               0);
                }

              if ((fp_uc != nullptr) and ((t == 0) or (opt_uc_allhits != 0)))
                {
                  auto const is_perfect = (d_total == 0);
                  std::fprintf(fp_uc,
                               "H\t%d\t%" PRId64 "\t%.1f\t+\t0\t0\t%s\t",
                               hp->centroid_index,
                               query_len,
                               hp->id,
                               is_perfect ? "=" : "*");
                  header_fprint_strip(fp_uc,
                                      query_header.c_str(),
                                      static_cast<int64_t>(query_header.size()),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                  std::fprintf(fp_uc, "\t");
                  header_fprint_strip(fp_uc,
                                      target.header.c_str(),
                                      static_cast<int64_t>(target.header.size()),
                                      opt_xsize,
                                      opt_xee,
                                      opt_xlength);
                  std::fprintf(fp_uc, "\n");
                }
            }

          if (fp_alnout != nullptr)
            {
              std::fprintf(fp_alnout, "\n");
              std::fprintf(fp_alnout, "Query >%s\n", query_header.c_str());
              std::fprintf(fp_alnout, " %%Id   TLen  Target\n");

              auto const top_hit_id = report_hits.front()->id;
              for (int64_t t = 0; t < toreport; ++t)
                {
                  auto const * hp = report_hits[static_cast<std::size_t>(t)];
                  if ((opt_top_hits_only != 0) and (hp->id < top_hit_id))
                    {
                      break;
                    }
                  auto const & target = db[static_cast<std::size_t>(hp->centroid_index)];
                  auto const target_len = static_cast<int64_t>(target.left.size() + target.right.size());
                  std::fprintf(fp_alnout,
                               "%3.0f%% %6" PRId64 "  %s\n",
                               hp->id,
                               target_len,
                               target.header.c_str());
                }

              for (int64_t t = 0; t < toreport; ++t)
                {
                  auto const * hp = report_hits[static_cast<std::size_t>(t)];
                  if ((opt_top_hits_only != 0) and (hp->id < top_hit_id))
                    {
                      break;
                    }

                  auto const & target = db[static_cast<std::size_t>(hp->centroid_index)];
                  auto const target_len = static_cast<int64_t>(target.left.size() + target.right.size());
                  auto const qlenlen = std::snprintf(nullptr, 0, "%" PRId64, query_len);
                  auto const tlenlen = std::snprintf(nullptr, 0, "%" PRId64, target_len);
                  auto const numwidth = std::max(qlenlen, tlenlen);
                  auto const rowlen = (opt_rowlen == 0)
                    ? static_cast<int>(query_len + target_len)
                    : static_cast<int>(opt_rowlen);

                  std::fprintf(fp_alnout, "\n");
                  std::fprintf(fp_alnout, " Query %*" PRId64 "nt >%s\n", numwidth, query_len, query_header.c_str());
                  std::fprintf(fp_alnout, "Target %*" PRId64 "nt >%s\n", numwidth, target_len, target.header.c_str());
                  std::fprintf(fp_alnout,
                               " Pair stats: id=%.1f%%, left=%.1f%%, right=%.1f%%\n",
                               hp->id,
                               hp->left.id,
                               hp->right.id);

                  auto const left_cigar_len = static_cast<int64_t>(hp->left.nwalignment.size())
                    - hp->left.trim_aln_left
                    - hp->left.trim_aln_right;
                  auto const right_cigar_len = static_cast<int64_t>(hp->right.nwalignment.size())
                    - hp->right.trim_aln_left
                    - hp->right.trim_aln_right;

                  std::fprintf(fp_alnout, " Left end (R1)\n");
                  align_show(fp_alnout,
                             qrec.left.c_str(),
                             static_cast<int64_t>(qrec.left.size()),
                             hp->left.trim_q_left,
                             "Qry1",
                             target.left.c_str(),
                             static_cast<int64_t>(target.left.size()),
                             hp->left.trim_t_left,
                             "Tgt1",
                             hp->left.nwalignment.c_str() + hp->left.trim_aln_left,
                             left_cigar_len,
                             numwidth,
                             4,
                             rowlen,
                             0);

                  std::fprintf(fp_alnout, " Right end (R2)\n");
                  align_show(fp_alnout,
                             qrec.right.c_str(),
                             static_cast<int64_t>(qrec.right.size()),
                             hp->right.trim_q_left,
                             "Qry2",
                             target.right.c_str(),
                             static_cast<int64_t>(target.right.size()),
                             hp->right.trim_t_left,
                             "Tgt2",
                             hp->right.nwalignment.c_str() + hp->right.trim_aln_left,
                             right_cigar_len,
                             numwidth,
                             4,
                             rowlen,
                             0);

                  auto const matches_total = hp->left.matches + hp->right.matches;
                  auto const indel_pct =
                    (hp->internal_alignment_cols_total > 0)
                    ? (100.0 * static_cast<double>(hp->internal_indels_total) /
                       static_cast<double>(hp->internal_alignment_cols_total))
                    : 0.0;
                  std::fprintf(fp_alnout,
                               "\n%d cols, %d ids (%3.1f%%), %d gaps (%3.1f%%)\n",
                               hp->internal_alignment_cols_total,
                               matches_total,
                               hp->id,
                               hp->internal_indels_total,
                               indel_pct);
                }
            }

          if (fp_matched != nullptr)
            {
              if (fp_matched2 != nullptr)
                {
                  write_query_pair_split(fp_matched,
                                         fp_matched2,
                                         query_header,
                                         qrec.left,
                                         qrec.right,
                                         qrec.abundance,
                                         matched_ordinal,
                                         matched2_ordinal);
                }
              else
                {
                  write_query_pair_interleaved(fp_matched,
                                               query_header,
                                               qrec.left,
                                               qrec.right,
                                               qrec.abundance,
                                               matched_ordinal);
                }
            }
        }
      else
        {
          if ((fp_alnout != nullptr) and (opt_output_no_hits != 0))
            {
              std::fprintf(fp_alnout, "\n");
              std::fprintf(fp_alnout, "Query >%s\n", query_header.c_str());
              std::fprintf(fp_alnout, "No hits\n");
            }

          if ((fp_userout != nullptr) and (opt_output_no_hits != 0))
            {
              std::fprintf(fp_userout,
                           "%s\t*\t0.000000\t0.000000\t0.000000\t0\t0\t0\n",
                           query_header.c_str());
            }

          if ((fp_blast6 != nullptr) and (opt_output_no_hits != 0))
            {
              std::fprintf(fp_blast6,
                           "%s\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n",
                           query_header.c_str());
            }

          if (fp_uc != nullptr)
            {
              std::fprintf(fp_uc, "N\t*\t*\t*\t.\t*\t*\t*\t");
              header_fprint_strip(fp_uc,
                                  query_header.c_str(),
                                  static_cast<int64_t>(query_header.size()),
                                  opt_xsize,
                                  opt_xee,
                                  opt_xlength);
              std::fprintf(fp_uc, "\t*\n");
            }

          if (fp_notmatched != nullptr)
            {
              if (fp_notmatched2 != nullptr)
                {
                  write_query_pair_split(fp_notmatched,
                                         fp_notmatched2,
                                         query_header,
                                         qrec.left,
                                         qrec.right,
                                         qrec.abundance,
                                         notmatched_ordinal,
                                         notmatched2_ordinal);
                }
              else
                {
                  write_query_pair_interleaved(fp_notmatched,
                                               query_header,
                                               qrec.left,
                                               qrec.right,
                                               qrec.abundance,
                                               notmatched_ordinal);
                }
            }

          if (any_otu_table_output)
            {
              if (opt_unknown_name != nullptr)
                {
                  otutable_add(query_header.c_str(), opt_unknown_name, qrec.abundance);
                }
              else
                {
                  otutable_add(query_header.c_str(), nullptr, qrec.abundance);
                }
            }
        }
    }

  if (fastx_next(qr, not opt_notrunclabels, chrmap_no_change_vector.data()))
    {
      fatal("More reverse queries than forward queries");
    }

  fastx_close(qr);
  fastx_close(qf);

  if (any_otu_table_output)
    {
      for (auto const & target : db)
        {
          otutable_add(nullptr, target.header.c_str(), 0);
        }
      if (opt_unknown_name != nullptr)
        {
          otutable_add(nullptr, opt_unknown_name, 0);
        }

      if (opt_biomout != nullptr)
        {
          auto * fp_biom = fopen_output(opt_biomout);
          if (fp_biom == nullptr)
            {
              fatal("Unable to open OTU table (biom 1.0 format) output file for writing");
            }
          otutable_print_biomout(fp_biom);
          std::fclose(fp_biom);
        }

      if (opt_otutabout != nullptr)
        {
          auto * fp_otu = fopen_output(opt_otutabout);
          if (fp_otu == nullptr)
            {
              fatal("Unable to open OTU table (classic format) output file for writing");
            }
          otutable_print_otutabout(fp_otu);
          std::fclose(fp_otu);
        }

      if (opt_mothur_shared_out != nullptr)
        {
          auto * fp_mothur = fopen_output(opt_mothur_shared_out);
          if (fp_mothur == nullptr)
            {
              fatal("Unable to open OTU table (mothur format) output file for writing");
            }
          otutable_print_mothur_shared_out(fp_mothur);
          std::fclose(fp_mothur);
        }

      otutable_done();
    }

  if ((fp_dbmatched != nullptr) or (fp_dbmatched2 != nullptr) or
      (fp_dbnotmatched != nullptr) or (fp_dbnotmatched2 != nullptr))
    {
      for (std::size_t i = 0; i < db.size(); ++i)
        {
          if ((dbmatched_abundance[i] > 0U) and (fp_dbmatched != nullptr))
            {
              if (fp_dbmatched2 != nullptr)
                {
                  write_db_pair_split(fp_dbmatched,
                                      fp_dbmatched2,
                                      db[i],
                                      static_cast<int64_t>(dbmatched_abundance[i]),
                                      dbmatched_ordinal,
                                      dbmatched2_ordinal);
                }
              else
                {
                  write_db_pair_interleaved(fp_dbmatched,
                                            db[i],
                                            static_cast<int64_t>(dbmatched_abundance[i]),
                                            dbmatched_ordinal);
                }
            }
          if ((dbmatched_abundance[i] == 0U) and (fp_dbnotmatched != nullptr))
            {
              if (fp_dbnotmatched2 != nullptr)
                {
                  write_db_pair_split(fp_dbnotmatched,
                                      fp_dbnotmatched2,
                                      db[i],
                                      db[i].abundance,
                                      dbnotmatched_ordinal,
                                      dbnotmatched2_ordinal);
                }
              else
                {
                  write_db_pair_interleaved(fp_dbnotmatched,
                                            db[i],
                                            db[i].abundance,
                                            dbnotmatched_ordinal);
                }
            }
        }
    }

  if (fp_alnout != nullptr)
    {
      std::fclose(fp_alnout);
    }
  if (fp_userout != nullptr)
    {
      std::fclose(fp_userout);
    }
  if (fp_blast6 != nullptr)
    {
      std::fclose(fp_blast6);
    }
  if (fp_uc != nullptr)
    {
      std::fclose(fp_uc);
    }
  if (fp_matched != nullptr)
    {
      std::fclose(fp_matched);
    }
  if (fp_matched2 != nullptr)
    {
      std::fclose(fp_matched2);
    }
  if (fp_notmatched != nullptr)
    {
      std::fclose(fp_notmatched);
    }
  if (fp_notmatched2 != nullptr)
    {
      std::fclose(fp_notmatched2);
    }
  if (fp_dbmatched != nullptr)
    {
      std::fclose(fp_dbmatched);
    }
  if (fp_dbmatched2 != nullptr)
    {
      std::fclose(fp_dbmatched2);
    }
  if (fp_dbnotmatched != nullptr)
    {
      std::fclose(fp_dbnotmatched);
    }
  if (fp_dbnotmatched2 != nullptr)
    {
      std::fclose(fp_dbnotmatched2);
    }

  unique_exit(kmer_unique_handle);

  if (!opt_quiet)
    {
      std::fprintf(stderr,
                   "TAV usearch_global: %" PRId64 " query pairs (%" PRId64 " filtered short, %" PRId64 " filtered long) -> %" PRId64 " matched (%.2f%%)\n",
                   queries,
                   filtered_short,
                   filtered_long,
                   assigned,
                   (queries > 0) ? (100.0 * static_cast<double>(assigned) / static_cast<double>(queries)) : 0.0);
      if ((opt_sizein != 0) and (queries_abundance > 0U))
        {
          std::fprintf(stderr,
                       "TAV usearch_global abundance-weighted matches: %" PRIu64 " of %" PRIu64 " (%.2f%%)\n",
                       assigned_abundance,
                       queries_abundance,
                       100.0 * static_cast<double>(assigned_abundance) / static_cast<double>(queries_abundance));
        }
    }
}
