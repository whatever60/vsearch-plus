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

#include "derep_paired.h"

#include "searchcore_paired.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "vsearch.h"

#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

auto derep_paired(struct Parameters const &parameters) -> void {
  if (parameters.opt_interleaved and (parameters.opt_reverse != nullptr)) {
    fatal("Paired fastx_uniques input must be either split (R1 + R2 "
          "positional input) or --interleaved, not both");
  }
  if ((not parameters.opt_interleaved) and (parameters.opt_reverse == nullptr)) {
    fatal("Paired fastx_uniques requires paired input: provide R2 as "
          "positional input or set --interleaved");
  }
  if ((parameters.opt_tabbedout == nullptr) and
      (parameters.opt_fastaout == nullptr) and
      (parameters.opt_fastaout_rev == nullptr)) {
    fatal("Paired fastx_uniques requires output via --tabbedout and/or "
          "--fastaout/--fastaout_rev");
  }

  std::unordered_map<std::string, std::size_t> pair_to_index;
  std::vector<record_paired_s> records;
  int64_t read_ordinal = 0;
  int64_t skipped_short = 0;

  if (parameters.opt_interleaved) {
    auto *fastx_h = fastx_open(parameters.opt_fastx_uniques);
    if (fastx_h == nullptr) {
      fatal("Unrecognized interleaved paired FASTX input: %s",
            parameters.opt_fastx_uniques);
    }

    while (fastx_next(fastx_h, not opt_notrunclabels,
                      chrmap_no_change_vector.data())) {
      auto const left_header = std::string{fastx_get_header(fastx_h)};
      auto const left_abundance = fastx_get_abundance(fastx_h);
      auto const left_sequence =
          std::string{fastx_get_sequence(fastx_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(fastx_h))};

      if (not fastx_next(fastx_h, not opt_notrunclabels,
                         chrmap_no_change_vector.data())) {
        fatal("Odd number of records in interleaved paired FASTX input %s; "
              "expected left/right entries",
              parameters.opt_fastx_uniques);
      }

      auto const right_header = std::string{fastx_get_header(fastx_h)};
      if (paired_header_key_paired(left_header) !=
          paired_header_key_paired(right_header)) {
        auto const message = std::string{"Paired FASTX headers differ ("} +
                             left_header + " vs " + right_header + ")";
        fatal(message.c_str());
      }
      auto const right_sequence =
          std::string{fastx_get_sequence(fastx_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(fastx_h))};
      auto anchor_len = std::min<int64_t>(
          static_cast<int64_t>(left_sequence.size()),
          static_cast<int64_t>(right_sequence.size()));
      if (opt_fastq_trunclen > 0) {
        anchor_len = std::min<int64_t>(anchor_len, opt_fastq_trunclen);
      }
      if ((anchor_len <= 0) or
          (static_cast<int64_t>(left_sequence.size()) < anchor_len) or
          (static_cast<int64_t>(right_sequence.size()) < anchor_len)) {
        ++skipped_short;
      } else {
        auto const left_anchor =
            left_sequence.substr(0U, static_cast<std::size_t>(anchor_len));
        auto const right_anchor =
            right_sequence.substr(0U, static_cast<std::size_t>(anchor_len));
        auto const key = left_anchor + '\t' + right_anchor;
        auto const pos = pair_to_index.find(key);
        if (pos == pair_to_index.end()) {
          record_paired_s record;
          record.header = left_header;
          record.header_r2 = right_header;
          record.qsequence_r1 = left_anchor;
          record.qsequence_r2 = right_anchor;
          record.abundance = left_abundance;
          record.first_seen = read_ordinal;
          pair_to_index.emplace(key, records.size());
          records.push_back(std::move(record));
        } else {
          records[pos->second].abundance += left_abundance;
        }
      }
      ++read_ordinal;
    }

    fastx_close(fastx_h);
  } else {
    auto *left_h = fastx_open(parameters.opt_fastx_uniques);
    if (left_h == nullptr) {
      fatal("Unrecognized left FASTX input: %s", parameters.opt_fastx_uniques);
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

      auto const left_header = std::string{fastx_get_header(left_h)};
      auto const right_header = std::string{fastx_get_header(right_h)};
      if (paired_header_key_paired(left_header) !=
          paired_header_key_paired(right_header)) {
        auto const message = std::string{"Paired FASTX headers differ ("} +
                             left_header + " vs " + right_header + ")";
        fatal(message.c_str());
      }
      auto const left_abundance = fastx_get_abundance(left_h);
      auto const left_sequence =
          std::string{fastx_get_sequence(left_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(left_h))};
      auto const right_sequence =
          std::string{fastx_get_sequence(right_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(right_h))};
      auto anchor_len = std::min<int64_t>(
          static_cast<int64_t>(left_sequence.size()),
          static_cast<int64_t>(right_sequence.size()));
      if (opt_fastq_trunclen > 0) {
        anchor_len = std::min<int64_t>(anchor_len, opt_fastq_trunclen);
      }
      if ((anchor_len <= 0) or
          (static_cast<int64_t>(left_sequence.size()) < anchor_len) or
          (static_cast<int64_t>(right_sequence.size()) < anchor_len)) {
        ++skipped_short;
      } else {
        auto const left_anchor =
            left_sequence.substr(0U, static_cast<std::size_t>(anchor_len));
        auto const right_anchor =
            right_sequence.substr(0U, static_cast<std::size_t>(anchor_len));
        auto const key = left_anchor + '\t' + right_anchor;
        auto const pos = pair_to_index.find(key);
        if (pos == pair_to_index.end()) {
          record_paired_s record;
          record.header = left_header;
          record.header_r2 = right_header;
          record.qsequence_r1 = left_anchor;
          record.qsequence_r2 = right_anchor;
          record.abundance = left_abundance;
          record.first_seen = read_ordinal;
          pair_to_index.emplace(key, records.size());
          records.push_back(std::move(record));
        } else {
          records[pos->second].abundance += left_abundance;
        }
      }
      ++read_ordinal;
    }

    if (fastx_next(right_h, not opt_notrunclabels,
                   chrmap_no_change_vector.data())) {
      fatal("More reverse records than forward records in paired FASTX input");
    }

    fastx_close(right_h);
    fastx_close(left_h);
  }

  std::sort(records.begin(), records.end(),
            [](record_paired_s const &lhs, record_paired_s const &rhs) {
              if (lhs.abundance != rhs.abundance) {
                return lhs.abundance > rhs.abundance;
              }
              if (lhs.header != rhs.header) {
                return lhs.header < rhs.header;
              }
              return lhs.first_seen < rhs.first_seen;
            });

  if (parameters.opt_tabbedout != nullptr) {
    auto *fp = fopen_output(parameters.opt_tabbedout);
    if (fp == nullptr) {
      fatal("Unable to open paired fastx_uniques catalog output file for "
            "writing");
    }
    std::fprintf(fp, "tav_id\tabundance\theader\tleft_anchor\tright_anchor\n");
    for (auto const &record : records) {
      std::fprintf(fp, "%s\t%" PRId64 "\t%s\t%s\t%s\n",
                   record.header.c_str(), record.abundance,
                   record.header.c_str(), record.qsequence_r1.c_str(),
                   record.qsequence_r2.c_str());
    }
    std::fclose(fp);
  }

  std::FILE *fp_left = nullptr;
  std::FILE *fp_right = nullptr;
  if (parameters.opt_fastaout != nullptr) {
    fp_left = fopen_output(parameters.opt_fastaout);
    if (fp_left == nullptr) {
      fatal("Unable to open paired fastx_uniques left FASTA output file for "
            "writing");
    }
  }
  if (parameters.opt_fastaout_rev != nullptr) {
    fp_right = fopen_output(parameters.opt_fastaout_rev);
    if (fp_right == nullptr) {
      if (fp_left != nullptr) {
        std::fclose(fp_left);
      }
      fatal("Unable to open paired fastx_uniques right FASTA output file for "
            "writing");
    }
  }

  int64_t ordinal = 1;
  for (auto const &record : records) {
    if (fp_left != nullptr) {
      fasta_print_general(fp_left, nullptr, record.qsequence_r1.c_str(),
                          static_cast<int>(record.qsequence_r1.size()),
                          record.header.c_str(),
                          static_cast<int>(record.header.size()),
                          static_cast<unsigned int>(record.abundance), ordinal,
                          -1.0, -1, -1, nullptr, 0.0);
    }
    if (fp_right != nullptr) {
      fasta_print_general(fp_right, nullptr, record.qsequence_r2.c_str(),
                          static_cast<int>(record.qsequence_r2.size()),
                          record.header_r2.c_str(),
                          static_cast<int>(record.header_r2.size()),
                          static_cast<unsigned int>(record.abundance), ordinal,
                          -1.0, -1, -1, nullptr, 0.0);
    }
    ++ordinal;
  }

  if (fp_left != nullptr) {
    std::fclose(fp_left);
  }
  if (fp_right != nullptr) {
    std::fclose(fp_right);
  }

  if (not opt_quiet) {
    std::fprintf(stderr,
                 "Paired fastx_uniques: %zu unique paired anchors, %" PRId64
                 " pairs skipped for anchor length constraints\n",
                 records.size(), skipped_short);
  }
}
