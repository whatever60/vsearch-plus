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


  GNU General Public License version 3

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The BSD 2-Clause License

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "vsearch.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include <algorithm>  // std::min, std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow, std::signbit
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose
#include <cstdlib>  // std::exit, EXIT_FAILURE
#include <limits>
#include <string>


inline auto fastq_get_qual(char const quality_symbol) -> int
{
  int const quality_score = quality_symbol - opt_fastq_ascii;

  if (quality_score < opt_fastq_qmin)
    {
      std::fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
              PRId64 ")\n",
              quality_score, opt_fastq_qmin);
      if (fp_log != nullptr)
        {
          std::fprintf(stderr,
                  "\n\nFatal error: FASTQ quality value (%d) below qmin (%"
                  PRId64 ")\n",
                  quality_score, opt_fastq_qmin);
        }
      std::exit(EXIT_FAILURE);
    }
  else if (quality_score > opt_fastq_qmax)
    {
      std::fprintf(stderr,
              "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
              PRId64 ")\n",
              quality_score, opt_fastq_qmax);
      std::fprintf(stderr,
              "By default, quality values range from 0 to 41.\n"
              "To allow higher quality values, "
              "please use the option --fastq_qmax %d\n", quality_score);
      if (fp_log != nullptr)
        {
          std::fprintf(fp_log,
                  "\n\nFatal error: FASTQ quality value (%d) above qmax (%"
                  PRId64 ")\n",
                  quality_score, opt_fastq_qmax);
          std::fprintf(fp_log,
                  "By default, quality values range from 0 to 41.\n"
                  "To allow higher quality values, "
                  "please use the option --fastq_qmax %d\n", quality_score);
        }
      std::exit(EXIT_FAILURE);
    }
  return quality_score;
}


struct analysis_res
{
  bool discarded = false;
  bool truncated = false;
  int start = 0;
  int length = 0;
  double ee = -1.0;
};


inline auto ee_thresholds_pass(double const ee, int const length) -> bool
{
  if (ee > opt_fastq_maxee)
    {
      return false;
    }
  if ((length > 0) and ((ee / length) > opt_fastq_maxee_rate))
    {
      return false;
    }
  return true;
}


auto analyse(fastx_handle input_handle, bool const apply_ee_filters) -> struct analysis_res
{
  auto const fastq_trunclen = static_cast<int>(opt_fastq_trunclen);
  auto const fastq_trunclen_keep = static_cast<int>(opt_fastq_trunclen_keep);
  struct analysis_res res;
  res.length = static_cast<int>(fastx_get_sequence_length(input_handle));
  auto const old_length = res.length;

  /* strip left (5') end */
  if (opt_fastq_stripleft < res.length)
    {
      res.start += opt_fastq_stripleft;
      res.length -= opt_fastq_stripleft;
    }
  else
    {
      res.start = res.length;
      res.length = 0;
    }

  /* strip right (3') end */
  if (opt_fastq_stripright < res.length)
    {
      res.length -= opt_fastq_stripright;
    }
  else
    {
      res.length = 0;
    }

  /* truncate trailing (3') part */
  if (opt_fastq_trunclen >= 0)
    {
      res.length = std::min(res.length, fastq_trunclen);
    }

  /* truncate trailing (3') part, but keep if short */
  if (opt_fastq_trunclen_keep >= 0)
    {
      res.length = std::min(res.length, fastq_trunclen_keep);
    }

  if (input_handle->is_fastq)
    {
      /* truncate by quality and expected errors (ee) */
      res.ee = 0.0;
      static constexpr auto base = 10.0;
      auto const * quality_symbols = fastx_get_quality(input_handle) + res.start;
      for (auto i = 0; i < res.length; ++i)
        {
          auto const quality_score = fastq_get_qual(quality_symbols[i]);
          auto const expected_error = std::pow(base, -quality_score / base);
          res.ee += expected_error;

          if ((quality_score <= opt_fastq_truncqual) or
              (res.ee > opt_fastq_truncee) or
              (res.ee > opt_fastq_truncee_rate * (i + 1)))
            {
              res.ee -= expected_error;
              res.length = i;
              break;
            }

          if (quality_score < opt_fastq_minqual)
            {
              res.discarded = true;
            }
        }

      /* filter by expected errors (ee) */
      if (apply_ee_filters and (not ee_thresholds_pass(res.ee, res.length)))
        {
          res.discarded = true;
        }
    }

  /* filter by length */
  if ((opt_fastq_trunclen >= 0) and (res.length < opt_fastq_trunclen))
    {
      res.discarded = true;
    }
  if (res.length < opt_fastq_minlen)
    {
      res.discarded = true;
    }
  if (res.length > opt_fastq_maxlen)
    {
      res.discarded = true;
    }

  /* filter by n's */  // refactoring: std::count_if();
  int64_t ncount = 0;
  auto const * nucleotides = fastx_get_sequence(input_handle) + res.start;
  for (auto i = 0; i < res.length; ++i)
    {
      auto const nucleotide = nucleotides[i];
      if ((nucleotide == 'N') or (nucleotide == 'n'))
        {
          ++ncount;
        }
    }
  if (ncount > opt_fastq_maxns)
    {
      res.discarded = true;
    }

  /* filter by abundance */
  auto const abundance = fastx_get_abundance(input_handle);
  if (abundance < opt_minsize)
    {
      res.discarded = true;
    }
  if (abundance > opt_maxsize)
    {
      res.discarded = true;
    }

  res.truncated = res.length < old_length;

  return res;
}


auto filter(bool const fastq_only, char * filename) -> void
{
  static constexpr auto dbl_max = std::numeric_limits<double>::max();  // refactoring: redundant?
  static constexpr auto long_min = std::numeric_limits<long>::min();

  if ((opt_fastqout == nullptr) and (opt_fastaout == nullptr) and
      (opt_fastqout_discarded == nullptr) and (opt_fastaout_discarded == nullptr) and
      (opt_fastqout_rev == nullptr) and (opt_fastaout_rev == nullptr) and
      (opt_fastqout_discarded_rev == nullptr) and (opt_fastaout_discarded_rev == nullptr))
    {
      fatal("No output files specified");
    }

  auto * forward_handle = fastx_open(filename);
  fastx_handle reverse_handle = nullptr;  // refactoring: direct initialization

  if (forward_handle == nullptr)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  if (not (forward_handle->is_fastq or forward_handle->is_empty))
    {
      if (fastq_only)
        {
          fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");
      }
      else if (opt_eeout or (opt_fastq_ascii != 33) or opt_fastq_eeout or
               (opt_fastq_maxee < dbl_max) or
               (opt_fastq_maxee_rate < dbl_max) or (opt_fastqout != nullptr) or
               (opt_fastq_qmax < 41) or (opt_fastq_qmin > 0) or
               (opt_fastq_truncee < dbl_max) or
               (opt_fastq_truncee_rate < dbl_max) or
               (opt_fastq_truncqual < long_min) or
               (opt_fastq_minqual > 0) or
               (opt_fastqout_discarded != nullptr) or
               (opt_fastqout_discarded_rev != nullptr) or
               (opt_fastqout_rev != nullptr))
        {
          fatal("The following options are not accepted with the fastx_filter command when the input is a FASTA file, because quality scores are not available: eeout, fastq_ascii, fastq_eeout, fastq_maxee, fastq_maxee_rate, fastq_minqual, fastq_out, fastq_qmax, fastq_qmin, fastq_truncee, fastq_truncee_rate, fastq_truncqual,  fastqout_discarded, fastqout_discarded_rev, fastqout_rev");
        }
    }

  auto const filesize = fastx_get_size(forward_handle);

  if (opt_reverse != nullptr)
    {
      reverse_handle = fastx_open(opt_reverse);

      if (reverse_handle == nullptr)
        {
          fatal("Unrecognized file type (not proper FASTA or FASTQ format) for reverse reads");
        }

      if (forward_handle->is_fastq != reverse_handle->is_fastq)
        {
          fatal("The forward and reverse input sequence must in the same format, either FASTA or FASTQ");
        }
    }

  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_fastaout_discarded = nullptr;
  std::FILE * fp_fastqout_discarded = nullptr;

  std::FILE * fp_fastaout_rev = nullptr;
  std::FILE * fp_fastqout_rev = nullptr;
  std::FILE * fp_fastaout_discarded_rev = nullptr;
  std::FILE * fp_fastqout_discarded_rev = nullptr;

  if (opt_fastaout != nullptr)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (fp_fastaout == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout != nullptr)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (fp_fastqout == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (opt_fastaout_discarded != nullptr)
    {
      fp_fastaout_discarded = fopen_output(opt_fastaout_discarded);
      if (fp_fastaout_discarded == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout_discarded != nullptr)
    {
      fp_fastqout_discarded = fopen_output(opt_fastqout_discarded);
      if (fp_fastqout_discarded == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (reverse_handle != nullptr)
    {
      if (opt_fastaout_rev != nullptr)
        {
          fp_fastaout_rev = fopen_output(opt_fastaout_rev);
          if (fp_fastaout_rev == nullptr)
            {
              fatal("Unable to open FASTA output file for writing");
            }
        }

      if (opt_fastqout_rev != nullptr)
        {
          fp_fastqout_rev = fopen_output(opt_fastqout_rev);
          if (fp_fastqout_rev == nullptr)
            {
              fatal("Unable to open FASTQ output file for writing");
            }
        }

      if (opt_fastaout_discarded_rev != nullptr)
        {
          fp_fastaout_discarded_rev = fopen_output(opt_fastaout_discarded_rev);
          if (fp_fastaout_discarded_rev == nullptr)
            {
              fatal("Unable to open FASTA output file for writing");
            }
        }

      if (opt_fastqout_discarded_rev != nullptr)
        {
          fp_fastqout_discarded_rev = fopen_output(opt_fastqout_discarded_rev);
          if (fp_fastqout_discarded_rev == nullptr)
            {
              fatal("Unable to open FASTQ output file for writing");
            }
        }
    }

  progress_init("Reading input file", filesize);

  int64_t kept = 0;
  int64_t discarded = 0;
  int64_t truncated = 0;

  while (fastx_next(forward_handle, false, chrmap_no_change_vector.data()))
    {
      if ((reverse_handle != nullptr) and not fastx_next(reverse_handle, false, chrmap_no_change_vector.data()))
        {
          fatal("More forward reads than reverse reads");
        }

      struct analysis_res res1;
      res1.ee = 0.0;
      struct analysis_res res2;

      res1 = analyse(forward_handle, true);
      if (reverse_handle != nullptr)
        {
          res2 = analyse(reverse_handle, true);
        }

      if (res1.discarded or res2.discarded)
        {
          /* discard the sequence(s) */

          ++discarded;

          if (opt_fastaout_discarded != nullptr)
            {
              fasta_print_general(fp_fastaout_discarded,
                                  nullptr,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_abundance(forward_handle),
                                  discarded,
                                  res1.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout_discarded != nullptr)
            {
              fastq_print_general(fp_fastqout_discarded,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_quality(forward_handle) + res1.start,
                                  fastx_get_abundance(forward_handle),
                                  discarded,
                                  res1.ee);
            }

          if (reverse_handle != nullptr)
            {
              if (opt_fastaout_discarded_rev != nullptr)
                {
                  fasta_print_general(fp_fastaout_discarded_rev,
                                      nullptr,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_abundance(reverse_handle),
                                      discarded,
                                      res2.ee,
                                      -1,
                                      -1,
                                      nullptr,
                                      0.0);
                }

              if (opt_fastqout_discarded_rev != nullptr)
                {
                  fastq_print_general(fp_fastqout_discarded_rev,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_quality(reverse_handle) + res2.start,
                                      fastx_get_abundance(reverse_handle),
                                      discarded,
                                      res2.ee);
                }
            }
        }
      else
        {
          /* keep the sequence(s) */

          ++kept;

          if (res1.truncated or res2.truncated)
            {
              ++truncated;
            }

          if (opt_fastaout != nullptr)
            {
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_abundance(forward_handle),
                                  kept,
                                  res1.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout != nullptr)
            {
              fastq_print_general(fp_fastqout,
                                  fastx_get_sequence(forward_handle) + res1.start,
                                  res1.length,
                                  fastx_get_header(forward_handle),
                                  fastx_get_header_length(forward_handle),
                                  fastx_get_quality(forward_handle) + res1.start,
                                  fastx_get_abundance(forward_handle),
                                  kept,
                                  res1.ee);
            }

          if (reverse_handle != nullptr)
            {
              if (opt_fastaout_rev != nullptr)
                {
                  fasta_print_general(fp_fastaout_rev,
                                      nullptr,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_abundance(reverse_handle),
                                      kept,
                                      res2.ee,
                                      -1,
                                      -1,
                                      nullptr,
                                      0.0);
                }

              if (opt_fastqout_rev != nullptr)
                {
                  fastq_print_general(fp_fastqout_rev,
                                      fastx_get_sequence(reverse_handle) + res2.start,
                                      res2.length,
                                      fastx_get_header(reverse_handle),
                                      fastx_get_header_length(reverse_handle),
                                      fastx_get_quality(reverse_handle) + res2.start,
                                      fastx_get_abundance(reverse_handle),
                                      kept,
                                      res2.ee);
                }
            }
        }

      progress_update(fastx_get_position(forward_handle));
    }

  progress_done();

  if ((reverse_handle != nullptr) and fastx_next(reverse_handle, false, chrmap_no_change_vector.data()))
    {
      fatal("More reverse reads than forward reads");
    }

  if (not opt_quiet)
    {
      std::fprintf(stderr,
              "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
              kept,
              truncated,
              discarded);
    }

  if (opt_log != nullptr)
    {
      std::fprintf(fp_log,
              "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
              kept,
              truncated,
              discarded);
    }

  if (reverse_handle != nullptr)
    {
      if (opt_fastaout_rev != nullptr)
        {
          std::fclose(fp_fastaout_rev);
        }

      if (opt_fastqout_rev != nullptr)
        {
          std::fclose(fp_fastqout_rev);
        }

      if (opt_fastaout_discarded_rev != nullptr)
        {
          std::fclose(fp_fastaout_discarded_rev);
        }

      if (opt_fastqout_discarded_rev != nullptr)
        {
          std::fclose(fp_fastqout_discarded_rev);
        }

      fastx_close(reverse_handle);
    }

  if (opt_fastaout != nullptr)
    {
      std::fclose(fp_fastaout);
    }

  if (opt_fastqout != nullptr)
    {
      std::fclose(fp_fastqout);
    }

  if (opt_fastaout_discarded != nullptr)
    {
      std::fclose(fp_fastaout_discarded);
    }

  if (opt_fastqout_discarded != nullptr)
    {
      std::fclose(fp_fastqout_discarded);
    }

  fastx_close(forward_handle);
}


auto filter_paired_ext_fastq(char * filename, bool const interleaved) -> void
{
  if ((opt_fastqout == nullptr) and
      (opt_fastqout2 == nullptr) and
      (opt_fastqout_discarded == nullptr) and
      (opt_fastqout_discarded2 == nullptr) and
      (opt_fastaout == nullptr) and
      (opt_fastaout2 == nullptr) and
      (opt_fastaout_discarded == nullptr) and
      (opt_fastaout_discarded2 == nullptr))
    {
      fatal("No output files specified");
    }

  auto * forward_handle = fastx_open(filename);
  if (forward_handle == nullptr)
    {
      fatal("Unrecognized file type (not proper FASTA or FASTQ format)");
    }

  if (not (forward_handle->is_fastq or forward_handle->is_empty))
    {
      fatal("FASTA input files not allowed with fastq_filter, consider using fastx_filter command instead");
    }

  fastx_handle reverse_handle = nullptr;
  if (not interleaved)
    {
      if (opt_reverse == nullptr)
        {
          fatal("No reverse reads file specified for paired fastq_filter extension mode");
        }

      reverse_handle = fastx_open(opt_reverse);
      if (reverse_handle == nullptr)
        {
          fatal("Unrecognized file type (not proper FASTA or FASTQ format) for reverse reads");
        }

      if (forward_handle->is_fastq != reverse_handle->is_fastq)
        {
          fatal("The forward and reverse input sequence must in the same format, either FASTA or FASTQ");
        }
    }

  std::FILE * fp_fastaout = nullptr;
  std::FILE * fp_fastaout2 = nullptr;
  std::FILE * fp_fastqout = nullptr;
  std::FILE * fp_fastqout2 = nullptr;
  std::FILE * fp_fastaout_discarded = nullptr;
  std::FILE * fp_fastaout_discarded2 = nullptr;
  std::FILE * fp_fastqout_discarded = nullptr;
  std::FILE * fp_fastqout_discarded2 = nullptr;

  if (opt_fastaout != nullptr)
    {
      fp_fastaout = fopen_output(opt_fastaout);
      if (fp_fastaout == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastaout2 != nullptr)
    {
      fp_fastaout2 = fopen_output(opt_fastaout2);
      if (fp_fastaout2 == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout != nullptr)
    {
      fp_fastqout = fopen_output(opt_fastqout);
      if (fp_fastqout == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (opt_fastqout2 != nullptr)
    {
      fp_fastqout2 = fopen_output(opt_fastqout2);
      if (fp_fastqout2 == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (opt_fastaout_discarded != nullptr)
    {
      fp_fastaout_discarded = fopen_output(opt_fastaout_discarded);
      if (fp_fastaout_discarded == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastaout_discarded2 != nullptr)
    {
      fp_fastaout_discarded2 = fopen_output(opt_fastaout_discarded2);
      if (fp_fastaout_discarded2 == nullptr)
        {
          fatal("Unable to open FASTA output file for writing");
        }
    }

  if (opt_fastqout_discarded != nullptr)
    {
      fp_fastqout_discarded = fopen_output(opt_fastqout_discarded);
      if (fp_fastqout_discarded == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  if (opt_fastqout_discarded2 != nullptr)
    {
      fp_fastqout_discarded2 = fopen_output(opt_fastqout_discarded2);
      if (fp_fastqout_discarded2 == nullptr)
        {
          fatal("Unable to open FASTQ output file for writing");
        }
    }

  struct read_snapshot
  {
    std::string header;
    std::string sequence;
    std::string quality;
    int64_t abundance = 0;
  };

  auto const snapshot_read = [](fastx_handle handle) -> struct read_snapshot
    {
      struct read_snapshot read;

      auto const header_length = fastx_get_header_length(handle);
      read.header.assign(fastx_get_header(handle), header_length);

      auto const sequence_length = fastx_get_sequence_length(handle);
      read.sequence.assign(fastx_get_sequence(handle), sequence_length);
      read.abundance = fastx_get_abundance(handle);

      if (handle->is_fastq)
        {
          read.quality.assign(fastx_get_quality(handle), sequence_length);
        }

      return read;
    };

  progress_init("Reading input file", fastx_get_size(forward_handle));

  int64_t kept = 0;
  int64_t discarded = 0;
  int64_t truncated = 0;

  while (fastx_next(forward_handle, false, chrmap_no_change_vector.data()))
    {
      auto const left_read = snapshot_read(forward_handle);
      auto const left_result = analyse(forward_handle, false);

      struct read_snapshot right_read;
      struct analysis_res right_result;

      if (interleaved)
        {
          if (not fastx_next(forward_handle, false, chrmap_no_change_vector.data()))
            {
              fatal("Odd number of records in interleaved paired FASTQ input %s; expected R1/R2 entries", filename);
            }
        }
      else if (not fastx_next(reverse_handle, false, chrmap_no_change_vector.data()))
        {
          fatal("More forward reads than reverse reads");
        }

      auto * const right_source = interleaved ? forward_handle : reverse_handle;
      right_read = snapshot_read(right_source);
      right_result = analyse(right_source, false);

      bool pair_discarded = left_result.discarded or right_result.discarded;
      auto const pair_ee = left_result.ee + right_result.ee;
      auto const pair_length = left_result.length + right_result.length;

      if ((not pair_discarded) and (not ee_thresholds_pass(pair_ee, pair_length)))
        {
          pair_discarded = true;
        }

      if (pair_discarded)
        {
          ++discarded;

          if (opt_fastaout_discarded != nullptr)
            {
              fasta_print_general(fp_fastaout_discarded,
                                  nullptr,
                                  left_read.sequence.data() + left_result.start,
                                  left_result.length,
                                  left_read.header.data(),
                                  static_cast<int>(left_read.header.size()),
                                  left_read.abundance,
                                  discarded,
                                  left_result.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout_discarded != nullptr)
            {
              fastq_print_general(fp_fastqout_discarded,
                                  left_read.sequence.data() + left_result.start,
                                  left_result.length,
                                  left_read.header.data(),
                                  static_cast<int>(left_read.header.size()),
                                  left_read.quality.data() + left_result.start,
                                  left_read.abundance,
                                  discarded,
                                  left_result.ee);
            }

          if (opt_fastaout_discarded2 != nullptr)
            {
              fasta_print_general(fp_fastaout_discarded2,
                                  nullptr,
                                  right_read.sequence.data() + right_result.start,
                                  right_result.length,
                                  right_read.header.data(),
                                  static_cast<int>(right_read.header.size()),
                                  right_read.abundance,
                                  discarded,
                                  right_result.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout_discarded2 != nullptr)
            {
              fastq_print_general(fp_fastqout_discarded2,
                                  right_read.sequence.data() + right_result.start,
                                  right_result.length,
                                  right_read.header.data(),
                                  static_cast<int>(right_read.header.size()),
                                  right_read.quality.data() + right_result.start,
                                  right_read.abundance,
                                  discarded,
                                  right_result.ee);
            }
        }
      else
        {
          ++kept;

          if (left_result.truncated or right_result.truncated)
            {
              ++truncated;
            }

          if (opt_fastaout != nullptr)
            {
              fasta_print_general(fp_fastaout,
                                  nullptr,
                                  left_read.sequence.data() + left_result.start,
                                  left_result.length,
                                  left_read.header.data(),
                                  static_cast<int>(left_read.header.size()),
                                  left_read.abundance,
                                  kept,
                                  left_result.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout != nullptr)
            {
              fastq_print_general(fp_fastqout,
                                  left_read.sequence.data() + left_result.start,
                                  left_result.length,
                                  left_read.header.data(),
                                  static_cast<int>(left_read.header.size()),
                                  left_read.quality.data() + left_result.start,
                                  left_read.abundance,
                                  kept,
                                  left_result.ee);
            }

          if (opt_fastaout2 != nullptr)
            {
              fasta_print_general(fp_fastaout2,
                                  nullptr,
                                  right_read.sequence.data() + right_result.start,
                                  right_result.length,
                                  right_read.header.data(),
                                  static_cast<int>(right_read.header.size()),
                                  right_read.abundance,
                                  kept,
                                  right_result.ee,
                                  -1,
                                  -1,
                                  nullptr,
                                  0.0);
            }

          if (opt_fastqout2 != nullptr)
            {
              fastq_print_general(fp_fastqout2,
                                  right_read.sequence.data() + right_result.start,
                                  right_result.length,
                                  right_read.header.data(),
                                  static_cast<int>(right_read.header.size()),
                                  right_read.quality.data() + right_result.start,
                                  right_read.abundance,
                                  kept,
                                  right_result.ee);
            }
        }

      progress_update(fastx_get_position(forward_handle));
    }

  progress_done();

  if ((not interleaved) and
      (reverse_handle != nullptr) and
      fastx_next(reverse_handle, false, chrmap_no_change_vector.data()))
    {
      fatal("More reverse reads than forward reads");
    }

  if (not opt_quiet)
    {
      std::fprintf(stderr,
              "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
              kept,
              truncated,
              discarded);
    }

  if (opt_log != nullptr)
    {
      std::fprintf(fp_log,
              "%" PRId64 " sequences kept (of which %" PRId64 " truncated), %" PRId64 " sequences discarded.\n",
              kept,
              truncated,
              discarded);
    }

  if (opt_fastaout != nullptr)
    {
      std::fclose(fp_fastaout);
    }

  if (opt_fastaout2 != nullptr)
    {
      std::fclose(fp_fastaout2);
    }

  if (opt_fastqout != nullptr)
    {
      std::fclose(fp_fastqout);
    }

  if (opt_fastqout2 != nullptr)
    {
      std::fclose(fp_fastqout2);
    }

  if (opt_fastaout_discarded != nullptr)
    {
      std::fclose(fp_fastaout_discarded);
    }

  if (opt_fastaout_discarded2 != nullptr)
    {
      std::fclose(fp_fastaout_discarded2);
    }

  if (opt_fastqout_discarded != nullptr)
    {
      std::fclose(fp_fastqout_discarded);
    }

  if (opt_fastqout_discarded2 != nullptr)
    {
      std::fclose(fp_fastqout_discarded2);
    }

  if (reverse_handle != nullptr)
    {
      fastx_close(reverse_handle);
    }

  fastx_close(forward_handle);
}


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  auto check_parameters(struct Parameters const & parameters) -> void {
    auto const is_negative = std::signbit(parameters.opt_fastq_truncee_rate);
    if (is_negative) {
      fatal("--fastq_truncee_rate cannot be negative");
    }

    if (parameters.opt_fastq_minqual < 0) {
      fatal("--fastq_minqual cannot be negative");
    }
  }

}  // end of anonymous namespace


auto fastq_filter(struct Parameters const & parameters) -> void
{
  check_parameters(parameters);
  filter(true, parameters.opt_fastq_filter);
}


auto fastq_filter_paired_ext(struct Parameters const & parameters) -> void
{
  check_parameters(parameters);
  filter_paired_ext_fastq(parameters.opt_fastq_filter, parameters.opt_interleaved);
}


auto fastx_filter(struct Parameters const & parameters) -> void
{
  check_parameters(parameters);
  filter(false, parameters.opt_fastx_filter);
}
