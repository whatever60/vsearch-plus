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

#include "align_simd.h"
#include "attributes.h"
#include "dbindex_paired.h"
#include "mask.h"
#include "minheap.h"
#include "otutable.h"
#include "searchcore_paired.h"
#include "showalign.h"
#include "udb.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/xpthread.hpp"
#include "vsearch.h"

#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <pthread.h>
#include <string>
#include <vector>

namespace {

auto get_anchor_len_paired(int64_t const len1, int64_t const len2) -> int64_t {
  auto const min_len = std::min(len1, len2);
  if (opt_fastq_trunclen > 0) {
    return std::min<int64_t>(min_len, opt_fastq_trunclen);
  }
  return min_len;
}

auto is_catalog_file_paired(char const *filename) -> bool {
  auto *fp = std::fopen(filename, "r");
  if (fp == nullptr) {
    return false;
  }

  char line[16] = {0};
  auto const read_ok = (std::fgets(line, sizeof(line), fp) != nullptr);
  std::fclose(fp);

  return read_ok and (std::strncmp(line, "tav_id\t", 7) == 0);
}

auto write_fasta_record_paired(FILE *fp, char const *header,
                               char const *sequence, int const seqlen,
                               int64_t const abundance, int &ordinal) -> void {
  if (fp == nullptr) {
    return;
  }

  fasta_print_general(fp, nullptr, sequence, seqlen, header,
                      static_cast<int>(std::strlen(header)), abundance,
                      ordinal, -1.0, -1, -1, nullptr, 0.0);
  ++ordinal;
}

auto write_query_pair_split_paired(FILE *fp_left, FILE *fp_right,
                                   char const *header_r1,
                                   char const *header_r2,
                                   char const *qsequence_r1,
                                   int const qseqlen_r1,
                                   char const *qsequence_r2,
                                   int const qseqlen_r2,
                                   int64_t const abundance,
                                   int &ordinal) -> void {
  write_fasta_record_paired(fp_left, header_r1, qsequence_r1, qseqlen_r1,
                            abundance, ordinal);
  write_fasta_record_paired(fp_right, header_r2, qsequence_r2, qseqlen_r2,
                            abundance, ordinal);
}

auto write_db_pair_split_paired(FILE *fp_left, FILE *fp_right,
                                record_paired_s const &record,
                                int64_t const abundance, int &ordinal)
    -> void {
  write_query_pair_split_paired(fp_left, fp_right, record.header.c_str(),
                                record.header_r2.c_str(),
                                record.qsequence_r1.c_str(),
                                static_cast<int>(record.qsequence_r1.size()),
                                record.qsequence_r2.c_str(),
                                static_cast<int>(record.qsequence_r2.size()),
                                abundance, ordinal);
}

static searchinfo_s_paired *si_plus;
static searchinfo_s_paired *si_minus;
static pthread_t *pthread;

static int tophits;
static int seqcount;
static pthread_attr_t attr;
static fastx_handle query_fastx_h_left;
static fastx_handle query_fastx_h_right;
static char const *query_filename_left;
static std::vector<record_paired_s> db_records;
static std::vector<unsigned int> target_seqnos_r1;
static std::vector<unsigned int> target_seqnos_r2;

static pthread_mutex_t mutex_input;
static pthread_mutex_t mutex_output;
static int qmatches;
static uint64_t qmatches_abundance;
static int queries;
static uint64_t queries_abundance;
static uint64_t *dbmatched;
static FILE *fp_alnout = nullptr;
static FILE *fp_userout = nullptr;
static FILE *fp_blast6out = nullptr;
static FILE *fp_uc = nullptr;
static FILE *fp_matched = nullptr;
static FILE *fp_matched2 = nullptr;
static FILE *fp_notmatched = nullptr;
static FILE *fp_notmatched2 = nullptr;
static FILE *fp_dbmatched = nullptr;
static FILE *fp_dbmatched2 = nullptr;
static FILE *fp_dbnotmatched = nullptr;
static FILE *fp_dbnotmatched2 = nullptr;
static FILE *fp_otutabout = nullptr;
static FILE *fp_mothur_shared_out = nullptr;
static FILE *fp_biomout = nullptr;

static int count_matched = 0;
static int count_notmatched = 0;

} // namespace

auto search_output_results_paired(std::vector<hit_paired_s> const &hits,
                                  char const *query_head,
                                  int const query_head_len,
                                  char const *query_head_r2,
                                  int const query_head_len_r2,
                                  int const qseqlen_r1,
                                  char const *qsequence_r1,
                                  int const qseqlen_r2,
                                  char const *qsequence_r2,
                                  int64_t const qsize) -> void {
  (void)query_head_len_r2;
  xpthread_mutex_lock(&mutex_output);

  auto const toreport =
      std::min(opt_maxhits, static_cast<int64_t>(hits.size()));
  auto const query_len = static_cast<int64_t>(qseqlen_r1 + qseqlen_r2);

  if (toreport != 0) {
    if ((opt_otutabout != nullptr) or (opt_mothur_shared_out != nullptr) or
        (opt_biomout != nullptr)) {
      otutable_add(query_head, db_getheader(hits[0].target_seqno_r1), qsize);
    }

    auto const top_hit_id = hits[0].id;

    for (auto t = 0; t < toreport; t++) {
      auto const *hp = &hits[static_cast<std::size_t>(t)];

      if ((opt_top_hits_only != 0) and (hp->id < top_hit_id)) {
        break;
      }

      auto const target_head = db_getheader(hp->target_seqno_r1);
      auto const target_head_len = db_getheaderlen(hp->target_seqno_r1);
      auto const target_len = static_cast<int64_t>(
          db_getsequencelen(hp->target_seqno_r1) +
          db_getsequencelen(hp->target_seqno_r2));
      auto const d_left = hp->r1.mismatches + hp->r1.internal_indels;
      auto const d_right = hp->r2.mismatches + hp->r2.internal_indels;
      auto const d_total = d_left + d_right;

      if (fp_userout != nullptr) {
        std::fprintf(fp_userout, "%s\t%s\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\n",
                     query_head, target_head, hp->r1.id / 100.0,
                     hp->r2.id / 100.0, hp->id / 100.0, d_left, d_right,
                     d_total);
      }

      if (fp_blast6out != nullptr) {
        std::fprintf(fp_blast6out,
                     "%s\t%s\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%" PRId64
                     "\t%d\t%d\n",
                     query_head, target_head, hp->id,
                     hp->internal_alignment_cols_total, hp->mismatches_total,
                     hp->internal_gaps_total, 1, static_cast<int>(query_len),
                     1, target_len, -1, 0);
      }

      if ((fp_uc != nullptr) and ((t == 0) or (opt_uc_allhits != 0))) {
        std::fprintf(fp_uc, "H\t%d\t%" PRId64 "\t%.1f\t+\t0\t0\t%s\t",
                     hp->target, query_len, hp->id,
                     (d_total == 0) ? "=" : "*");
        header_fprint_strip(fp_uc, query_head, query_head_len, opt_xsize,
                            opt_xee, opt_xlength);
        std::fprintf(fp_uc, "\t");
        header_fprint_strip(fp_uc, target_head, target_head_len, opt_xsize,
                            opt_xee, opt_xlength);
        std::fprintf(fp_uc, "\n");
      }
    }

    if (fp_alnout != nullptr) {
      std::fprintf(fp_alnout, "\n");
      std::fprintf(fp_alnout, "Query >%s\n", query_head);
      std::fprintf(fp_alnout, " %%Id   TLen  Target\n");

      auto const top_hit_id = hits[0].id;
      for (auto t = 0; t < toreport; t++) {
        auto const *hp = &hits[static_cast<std::size_t>(t)];
        if ((opt_top_hits_only != 0) and (hp->id < top_hit_id)) {
          break;
        }
        auto const target_len = static_cast<int64_t>(
            db_getsequencelen(hp->target_seqno_r1) +
            db_getsequencelen(hp->target_seqno_r2));
        std::fprintf(fp_alnout, "%3.0f%% %6" PRId64 "  %s\n", hp->id,
                     target_len, db_getheader(hp->target_seqno_r1));
      }

      for (auto t = 0; t < toreport; t++) {
        auto const *hp = &hits[static_cast<std::size_t>(t)];
        if ((opt_top_hits_only != 0) and (hp->id < top_hit_id)) {
          break;
        }

        auto const target_r1 = db_getsequence(hp->target_seqno_r1);
        auto const target_r2 = db_getsequence(hp->target_seqno_r2);
        auto const target_len = static_cast<int64_t>(
            db_getsequencelen(hp->target_seqno_r1) +
            db_getsequencelen(hp->target_seqno_r2));
        auto const qlenlen = std::snprintf(nullptr, 0, "%" PRId64, query_len);
        auto const tlenlen =
            std::snprintf(nullptr, 0, "%" PRId64, target_len);
        auto const numwidth = std::max(qlenlen, tlenlen);
        auto const rowlen = (opt_rowlen == 0)
                                ? static_cast<int>(query_len + target_len)
                                : static_cast<int>(opt_rowlen);
        auto const left_cigar_len =
            static_cast<int64_t>(std::strlen(hp->r1.nwalignment)) -
            hp->r1.trim_aln_left - hp->r1.trim_aln_right;
        auto const right_cigar_len =
            static_cast<int64_t>(std::strlen(hp->r2.nwalignment)) -
            hp->r2.trim_aln_left - hp->r2.trim_aln_right;
        auto const matches_total = hp->r1.matches + hp->r2.matches;
        auto const indel_pct =
            (hp->internal_alignment_cols_total > 0)
                ? (100.0 * static_cast<double>(hp->internal_indels_total) /
                   static_cast<double>(hp->internal_alignment_cols_total))
                : 0.0;

        std::fprintf(fp_alnout, "\n");
        std::fprintf(fp_alnout, " Query %*" PRId64 "nt >%s\n", numwidth,
                     query_len, query_head);
        std::fprintf(fp_alnout, "Target %*" PRId64 "nt >%s\n", numwidth,
                     target_len, db_getheader(hp->target_seqno_r1));
        std::fprintf(fp_alnout,
                     " Pair stats: id=%.1f%%, left=%.1f%%, right=%.1f%%\n",
                     hp->id, hp->r1.id, hp->r2.id);

        std::fprintf(fp_alnout, " Left end (R1)\n");
        align_show(fp_alnout, qsequence_r1, qseqlen_r1, hp->r1.trim_q_left,
                   "Qry1", target_r1, db_getsequencelen(hp->target_seqno_r1),
                   hp->r1.trim_t_left, "Tgt1",
                   hp->r1.nwalignment + hp->r1.trim_aln_left, left_cigar_len,
                   numwidth, 4, rowlen, 0);

        std::fprintf(fp_alnout, " Right end (R2)\n");
        align_show(fp_alnout, qsequence_r2, qseqlen_r2, hp->r2.trim_q_left,
                   "Qry2", target_r2, db_getsequencelen(hp->target_seqno_r2),
                   hp->r2.trim_t_left, "Tgt2",
                   hp->r2.nwalignment + hp->r2.trim_aln_left, right_cigar_len,
                   numwidth, 4, rowlen, 0);

        std::fprintf(fp_alnout,
                     "\n%d cols, %d ids (%3.1f%%), %d gaps (%3.1f%%)\n",
                     hp->internal_alignment_cols_total, matches_total, hp->id,
                     hp->internal_indels_total, indel_pct);
      }
    }
  } else {
    if ((opt_otutabout != nullptr) or (opt_mothur_shared_out != nullptr) or
        (opt_biomout != nullptr)) {
      if (opt_unknown_name != nullptr) {
        otutable_add(query_head, opt_unknown_name, qsize);
      } else {
        otutable_add(query_head, nullptr, qsize);
      }
    }

    if ((fp_alnout != nullptr) and (opt_output_no_hits != 0)) {
      std::fprintf(fp_alnout, "\n");
      std::fprintf(fp_alnout, "Query >%s\n", query_head);
      std::fprintf(fp_alnout, "No hits\n");
    }

    if ((fp_userout != nullptr) and (opt_output_no_hits != 0)) {
      std::fprintf(fp_userout,
                   "%s\t*\t0.000000\t0.000000\t0.000000\t0\t0\t0\n",
                   query_head);
    }

    if ((fp_blast6out != nullptr) and (opt_output_no_hits != 0)) {
      std::fprintf(fp_blast6out,
                   "%s\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n",
                   query_head);
    }

    if (fp_uc != nullptr) {
      std::fprintf(fp_uc, "N\t*\t*\t*\t.\t*\t*\t*\t");
      header_fprint_strip(fp_uc, query_head, query_head_len, opt_xsize,
                          opt_xee, opt_xlength);
      std::fprintf(fp_uc, "\t*\n");
    }
  }

  if (!hits.empty()) {
    ++count_matched;
    if (opt_matched != nullptr) {
      write_query_pair_split_paired(fp_matched, fp_matched2, query_head,
                                    query_head_r2,
                                    qsequence_r1, qseqlen_r1, qsequence_r2,
                                    qseqlen_r2, qsize, count_matched);
    }
  } else {
    ++count_notmatched;
    if (opt_notmatched != nullptr) {
      write_query_pair_split_paired(fp_notmatched, fp_notmatched2, query_head,
                                    query_head_r2,
                                    qsequence_r1, qseqlen_r1, qsequence_r2,
                                    qseqlen_r2, qsize, count_notmatched);
    }
  }

  auto const dbmatch_weight = static_cast<uint64_t>(opt_sizein ? qsize : 1);
  for (auto const &hit : hits) {
    if (hit.accepted or hit.weak) {
      dbmatched[hit.target] += dbmatch_weight;
    }
  }

  xpthread_mutex_unlock(&mutex_output);
}

auto search_query_paired(int64_t const t) -> int {
  for (int s = 0; s < opt_strand; s++) {
    auto *si = (s != 0) ? si_minus + t : si_plus + t;

    if (opt_qmask == MASK_DUST) {
      dust(si->qsequence_r1, si->qseqlen_r1);
      dust(si->qsequence_r2, si->qseqlen_r2);
    } else if ((opt_qmask == MASK_SOFT) and (opt_hardmask != 0)) {
      hardmask(si->qsequence_r1, si->qseqlen_r1);
      hardmask(si->qsequence_r2, si->qseqlen_r2);
    }

    search_onequery_paired(si, opt_qmask);
  }

  std::vector<hit_paired_s> hits;
  search_joinhits_paired(si_plus + t, opt_strand > 1 ? si_minus + t : nullptr,
                         hits);

  search_output_results_paired(hits, si_plus[t].query_head,
                               si_plus[t].query_head_len,
                               si_plus[t].query_head_r2,
                               si_plus[t].query_head_len_r2,
                               si_plus[t].qseqlen_r1, si_plus[t].qsequence_r1,
                               si_plus[t].qseqlen_r2, si_plus[t].qsequence_r2,
                               si_plus[t].qsize);

  for (auto const &hit : hits) {
    if (hit.aligned) {
      xfree(hit.r1.nwalignment);
      xfree(hit.r2.nwalignment);
    }
  }

  return static_cast<int>(hits.size());
}

auto search_thread_run_paired(int64_t const t) -> void {
  while (true) {
    xpthread_mutex_lock(&mutex_input);

    if (fastx_next(query_fastx_h_left, (opt_notrunclabels == 0),
                   chrmap_no_change_vector.data())) {
      std::string query_head{fastx_get_header(query_fastx_h_left)};
      auto const query_no = fastx_get_seqno(query_fastx_h_left);
      auto const qsize =
          std::max<int64_t>(fastx_get_abundance(query_fastx_h_left), 1);
      auto const qseqlen_r1 = fastx_get_sequence_length(query_fastx_h_left);
      std::string qsequence_r1{fastx_get_sequence(query_fastx_h_left),
                               static_cast<std::size_t>(qseqlen_r1)};

      std::string query_head_r2;
      std::string qsequence_r2;
      int qseqlen_r2 = 0;
      if (query_fastx_h_right == nullptr) {
        if (not fastx_next(query_fastx_h_left, (opt_notrunclabels == 0),
                           chrmap_no_change_vector.data())) {
          fatal("Odd number of records in interleaved paired FASTX input %s; "
                "expected left/right entries",
                query_filename_left);
        }
        query_head_r2 = fastx_get_header(query_fastx_h_left);
        qseqlen_r2 = fastx_get_sequence_length(query_fastx_h_left);
        qsequence_r2.assign(fastx_get_sequence(query_fastx_h_left),
                            static_cast<std::size_t>(qseqlen_r2));
      } else {
        if (not fastx_next(query_fastx_h_right, (opt_notrunclabels == 0),
                           chrmap_no_change_vector.data())) {
          fatal("More forward queries than reverse queries");
        }
        query_head_r2 = fastx_get_header(query_fastx_h_right);
        qseqlen_r2 = fastx_get_sequence_length(query_fastx_h_right);
        qsequence_r2.assign(fastx_get_sequence(query_fastx_h_right),
                            static_cast<std::size_t>(qseqlen_r2));
      }
      if (paired_header_key_paired(query_head) !=
          paired_header_key_paired(query_head_r2)) {
        auto const message = std::string{"Paired query headers differ ("} +
                             query_head + " vs " + query_head_r2 + ")";
        fatal(message.c_str());
      }

      auto const anchor_len =
          static_cast<int>(std::max<int64_t>(0, get_anchor_len_paired(
                                                   qseqlen_r1, qseqlen_r2)));

      for (int s = 0; s < opt_strand; s++) {
        auto *si = (s != 0) ? si_minus + t : si_plus + t;
        si->query_head_len = static_cast<int>(query_head.size());
        si->query_head_len_r2 = static_cast<int>(query_head_r2.size());
        si->qseqlen_r1 = anchor_len;
        si->qseqlen_r2 = anchor_len;
        si->query_no = query_no;
        si->qsize = qsize;
        si->strand = s;

        if (si->query_head_len + 1 > si->query_head_alloc) {
          si->query_head_alloc = si->query_head_len + 2001;
          si->query_head = static_cast<char *>(
              xrealloc(si->query_head, static_cast<std::size_t>(si->query_head_alloc)));
        }
        if (si->query_head_len_r2 + 1 > si->query_head_alloc_r2) {
          si->query_head_alloc_r2 = si->query_head_len_r2 + 2001;
          si->query_head_r2 = static_cast<char *>(
              xrealloc(si->query_head_r2,
                       static_cast<std::size_t>(si->query_head_alloc_r2)));
        }

        if (anchor_len + 1 > si->seq_alloc) {
          si->seq_alloc = anchor_len + 2001;
          si->qsequence_r1 = static_cast<char *>(
              xrealloc(si->qsequence_r1, static_cast<std::size_t>(si->seq_alloc)));
          si->qsequence_r2 = static_cast<char *>(
              xrealloc(si->qsequence_r2, static_cast<std::size_t>(si->seq_alloc)));
        }
      }

      std::strcpy(si_plus[t].query_head, query_head.c_str());
      std::strcpy(si_plus[t].query_head_r2, query_head_r2.c_str());
      std::memcpy(si_plus[t].qsequence_r1, qsequence_r1.data(),
                  static_cast<std::size_t>(anchor_len));
      si_plus[t].qsequence_r1[anchor_len] = '\0';
      std::memcpy(si_plus[t].qsequence_r2, qsequence_r2.data(),
                  static_cast<std::size_t>(anchor_len));
      si_plus[t].qsequence_r2[anchor_len] = '\0';

      auto const progress = fastx_get_position(query_fastx_h_left);

      xpthread_mutex_unlock(&mutex_input);

      if (opt_strand > 1) {
        std::strcpy(si_minus[t].query_head, si_plus[t].query_head);
        std::strcpy(si_minus[t].query_head_r2, si_plus[t].query_head_r2);
        reverse_complement(si_minus[t].qsequence_r1, si_plus[t].qsequence_r1,
                           si_plus[t].qseqlen_r1);
        reverse_complement(si_minus[t].qsequence_r2, si_plus[t].qsequence_r2,
                           si_plus[t].qseqlen_r2);
      }

      auto const match = search_query_paired(t);

      xpthread_mutex_lock(&mutex_output);
      ++queries;
      queries_abundance += qsize;
      if (match != 0) {
        ++qmatches;
        qmatches_abundance += qsize;
      }
      progress_update(progress);
      xpthread_mutex_unlock(&mutex_output);
    } else {
      xpthread_mutex_unlock(&mutex_input);
      break;
    }
  }
}

auto search_thread_init_paired(searchinfo_s_paired *si) -> void {
  si->uh_r1 = unique_init();
  si->uh_r2 = unique_init();
  si->kmers = static_cast<count_t *>(
      xmalloc((static_cast<std::size_t>(seqcount) * sizeof(count_t)) + 32));
  si->m = minheap_init(tophits);
  si->hits = static_cast<hit_paired_s *>(
      xmalloc(sizeof(hit_paired_s) * tophits * opt_strand));
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = nullptr;
  si->query_head_alloc_r2 = 0;
  si->query_head_r2 = nullptr;
  si->seq_alloc = 0;
  si->qsequence_r1 = nullptr;
  si->qsequence_r2 = nullptr;
  si->target_seqnos_r1 = &target_seqnos_r1;
  si->target_seqnos_r2 = &target_seqnos_r2;
  si->s_r1 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);
  si->s_r2 = search16_init(
      opt_match, opt_mismatch, opt_gap_open_query_left,
      opt_gap_open_target_left, opt_gap_open_query_interior,
      opt_gap_open_target_interior, opt_gap_open_query_right,
      opt_gap_open_target_right, opt_gap_extension_query_left,
      opt_gap_extension_target_left, opt_gap_extension_query_interior,
      opt_gap_extension_target_interior, opt_gap_extension_query_right,
      opt_gap_extension_target_right);
}

auto search_thread_exit_paired(searchinfo_s_paired *si) -> void {
  search16_exit(si->s_r1);
  search16_exit(si->s_r2);
  unique_exit(si->uh_r1);
  unique_exit(si->uh_r2);
  xfree(si->hits);
  minheap_exit(si->m);
  xfree(si->kmers);
  if (si->query_head != nullptr) {
    xfree(si->query_head);
  }
  if (si->query_head_r2 != nullptr) {
    xfree(si->query_head_r2);
  }
  if (si->qsequence_r1 != nullptr) {
    xfree(si->qsequence_r1);
  }
  if (si->qsequence_r2 != nullptr) {
    xfree(si->qsequence_r2);
  }
}

auto search_thread_worker_paired(void *vp) -> void * {
  auto const t = (int64_t)vp;
  search_thread_run_paired(t);
  return nullptr;
}

auto search_thread_worker_run_paired() -> void {
  xpthread_attr_init(&attr);
  xpthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (int t = 0; t < opt_threads; t++) {
    search_thread_init_paired(si_plus + t);
    if (si_minus != nullptr) {
      search_thread_init_paired(si_minus + t);
    }
    xpthread_create(pthread + t, &attr, search_thread_worker_paired,
                    (void *)(int64_t)t);
  }

  for (int t = 0; t < opt_threads; t++) {
    xpthread_join(pthread[t], nullptr);
    search_thread_exit_paired(si_plus + t);
    if (si_minus != nullptr) {
      search_thread_exit_paired(si_minus + t);
    }
  }

  xpthread_attr_destroy(&attr);
}

auto search_prep_paired(char *cmdline, char *progheader) -> void {
  if (opt_alnout != nullptr) {
    fp_alnout = fopen_output(opt_alnout);
    if (fp_alnout == nullptr) {
      fatal("Unable to open alignment output file for writing");
    }
    std::fprintf(fp_alnout, "%s\n", cmdline);
    std::fprintf(fp_alnout, "%s\n", progheader);
  }

  if (opt_userout != nullptr) {
    fp_userout = fopen_output(opt_userout);
    if (fp_userout == nullptr) {
      fatal("Unable to open paired userout output file for writing");
    }
  }

  if (opt_blast6out != nullptr) {
    fp_blast6out = fopen_output(opt_blast6out);
    if (fp_blast6out == nullptr) {
      fatal("Unable to open blast6-like output file for writing");
    }
  }

  if (opt_uc != nullptr) {
    fp_uc = fopen_output(opt_uc);
    if (fp_uc == nullptr) {
      fatal("Unable to open uc output file for writing");
    }
  }

  if (opt_matched != nullptr) {
    fp_matched = fopen_output(opt_matched);
    if (fp_matched == nullptr) {
      fatal("Unable to open matched output file for writing");
    }
  }

  if (opt_matched2 != nullptr) {
    fp_matched2 = fopen_output(opt_matched2);
    if (fp_matched2 == nullptr) {
      fatal("Unable to open matched2 output file for writing");
    }
  }

  if (opt_notmatched != nullptr) {
    fp_notmatched = fopen_output(opt_notmatched);
    if (fp_notmatched == nullptr) {
      fatal("Unable to open notmatched output file for writing");
    }
  }

  if (opt_notmatched2 != nullptr) {
    fp_notmatched2 = fopen_output(opt_notmatched2);
    if (fp_notmatched2 == nullptr) {
      fatal("Unable to open notmatched2 output file for writing");
    }
  }

  if (opt_otutabout != nullptr) {
    fp_otutabout = fopen_output(opt_otutabout);
    if (fp_otutabout == nullptr) {
      fatal("Unable to open OTU table (text format) output file for writing");
    }
  }

  if (opt_mothur_shared_out != nullptr) {
    fp_mothur_shared_out = fopen_output(opt_mothur_shared_out);
    if (fp_mothur_shared_out == nullptr) {
      fatal("Unable to open OTU table (mothur format) output file for writing");
    }
  }

  if (opt_biomout != nullptr) {
    fp_biomout = fopen_output(opt_biomout);
    if (fp_biomout == nullptr) {
      fatal("Unable to open OTU table (biom 1.0 format) output file for writing");
    }
  }

  if (udb_detect_isudb(opt_db) or
      ((opt_db2 != nullptr) and udb_detect_isudb(opt_db2))) {
    fatal("Paired usearch_global UDB input is not supported; provide paired FASTA/FASTQ database records");
  }

  if (is_catalog_file_paired(opt_db) or
      ((opt_db2 != nullptr) and is_catalog_file_paired(opt_db2))) {
    fatal("Paired usearch_global catalog input is not supported; provide paired FASTA/FASTQ database records");
  }

  db_records.clear();
  target_seqnos_r1.clear();
  target_seqnos_r2.clear();

  auto append_record = [&](std::string const &left_header,
                           std::string const &right_header,
                           std::string const &left_sequence,
                           std::string const &right_sequence,
                           int64_t const left_abundance,
                           int64_t const right_abundance) -> void {
    auto const left_len = static_cast<int64_t>(left_sequence.size());
    auto const right_len = static_cast<int64_t>(right_sequence.size());
    auto const anchor_len = get_anchor_len_paired(left_len, right_len);
    if (anchor_len <= 0) {
      return;
    }

    auto const left_short = left_len < opt_minseqlength;
    auto const right_short = right_len < opt_minseqlength;
    auto const left_long = left_len > opt_maxseqlength;
    auto const right_long = right_len > opt_maxseqlength;
    auto const drop_short =
        (opt_filter == 0) ? (left_short or right_short)
                          : (left_short and right_short);
    auto const drop_long =
        (opt_filter == 0) ? (left_long or right_long)
                          : (left_long and right_long);
    if (drop_short or drop_long) {
      return;
    }

    if (paired_header_key_paired(left_header) !=
        paired_header_key_paired(right_header)) {
      auto const message = std::string{"Paired database headers differ ("} +
                           left_header + " vs " + right_header + ")";
      fatal(message.c_str());
    }

    record_paired_s record;
    record.header = left_header;
    record.header_r2 = right_header;
    if ((right_abundance != left_abundance) and (not opt_quiet)) {
      std::fprintf(stderr,
                   "Warning: paired database abundances differ (%" PRId64
                   " vs %" PRId64 "); using left abundance\n",
                   left_abundance, right_abundance);
    }

    record.qsequence_r1.assign(left_sequence.data(),
                               static_cast<std::size_t>(anchor_len));
    record.qsequence_r2.assign(right_sequence.data(),
                               static_cast<std::size_t>(anchor_len));
    if ((opt_dbmask == MASK_DUST) or
        ((opt_dbmask == MASK_SOFT) and (opt_hardmask != 0))) {
      std::vector<char> left_buffer(record.qsequence_r1.begin(),
                                    record.qsequence_r1.end());
      std::vector<char> right_buffer(record.qsequence_r2.begin(),
                                     record.qsequence_r2.end());
      left_buffer.push_back('\0');
      right_buffer.push_back('\0');

      if (opt_dbmask == MASK_DUST) {
        dust(left_buffer.data(), static_cast<int>(record.qsequence_r1.size()));
        dust(right_buffer.data(), static_cast<int>(record.qsequence_r2.size()));
      } else {
        hardmask(left_buffer.data(),
                 static_cast<int>(record.qsequence_r1.size()));
        hardmask(right_buffer.data(),
                 static_cast<int>(record.qsequence_r2.size()));
      }

      record.qsequence_r1.assign(left_buffer.data(),
                                 record.qsequence_r1.size());
      record.qsequence_r2.assign(right_buffer.data(),
                                 record.qsequence_r2.size());
    }
    record.abundance = std::max<int64_t>(left_abundance, 1);
    record.first_seen = static_cast<int64_t>(db_records.size());
    db_records.push_back(std::move(record));
  };

  if (opt_db2 != nullptr) {
    auto *db_left_h = fastx_open(opt_db);
    if (db_left_h == nullptr) {
      fatal("Unrecognized file type for paired left database (not proper FASTA or FASTQ format): %s",
            opt_db);
    }
    auto *db_right_h = fastx_open(opt_db2);
    if (db_right_h == nullptr) {
      fastx_close(db_left_h);
      fatal("Unrecognized file type for paired right database (not proper FASTA or FASTQ format): %s",
            opt_db2);
    }

    while (fastx_next(db_left_h, (opt_notrunclabels == 0),
                      chrmap_no_change_vector.data())) {
      if (not fastx_next(db_right_h, (opt_notrunclabels == 0),
                         chrmap_no_change_vector.data())) {
        fatal("More left database records than right database records");
      }

      append_record(
          std::string{fastx_get_header(db_left_h)},
          std::string{fastx_get_header(db_right_h)},
          std::string{fastx_get_sequence(db_left_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(db_left_h))},
          std::string{fastx_get_sequence(db_right_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(db_right_h))},
          fastx_get_abundance(db_left_h), fastx_get_abundance(db_right_h));
    }

    if (fastx_next(db_right_h, (opt_notrunclabels == 0),
                   chrmap_no_change_vector.data())) {
      fatal("More right database records than left database records");
    }

    fastx_close(db_right_h);
    fastx_close(db_left_h);
  } else {
    auto *db_h = fastx_open(opt_db);
    if (db_h == nullptr) {
      fatal("Unrecognized file type for paired database (not proper FASTA or FASTQ format): %s",
            opt_db);
    }

    while (fastx_next(db_h, (opt_notrunclabels == 0),
                      chrmap_no_change_vector.data())) {
      auto const left_header = std::string{fastx_get_header(db_h)};
      auto const left_sequence = std::string{
          fastx_get_sequence(db_h),
          static_cast<std::size_t>(fastx_get_sequence_length(db_h))};
      auto const left_abundance = fastx_get_abundance(db_h);

      if (not fastx_next(db_h, (opt_notrunclabels == 0),
                         chrmap_no_change_vector.data())) {
        fatal("Odd number of records in paired FASTX database %s; expected interleaved left/right entries",
              opt_db);
      }

      append_record(
          left_header, std::string{fastx_get_header(db_h)}, left_sequence,
          std::string{fastx_get_sequence(db_h),
                      static_cast<std::size_t>(
                          fastx_get_sequence_length(db_h))},
          left_abundance, fastx_get_abundance(db_h));
    }

    fastx_close(db_h);
  }

  if (db_records.empty()) {
    if (opt_db2 != nullptr) {
      fatal("No paired database records loaded from the provided --db/--db2 files");
    }
    fatal("No paired database records loaded from %s", opt_db);
  }

  db_free();
  db_setinfo(false, 0, 0, 0, std::numeric_limits<uint64_t>::max(), 0);

  for (auto const &record : db_records) {
    auto const left_seqno = static_cast<unsigned int>(db_getsequencecount());
    db_add(false, record.header.c_str(), record.qsequence_r1.c_str(), nullptr,
           record.header.size(), record.qsequence_r1.size(), record.abundance);
    auto const right_seqno = static_cast<unsigned int>(db_getsequencecount());
    db_add(false, record.header_r2.c_str(), record.qsequence_r2.c_str(),
           nullptr, record.header_r2.size(), record.qsequence_r2.size(),
           record.abundance);
    target_seqnos_r1.push_back(left_seqno);
    target_seqnos_r2.push_back(right_seqno);
  }

  show_rusage();

  seqcount = static_cast<int>(db_records.size());
  dbindex_prepare_paired(&db_records, 1, opt_dbmask);
  dbindex_addallsequences_paired(opt_dbmask);

  if ((opt_maxrejects == 0) or (opt_maxrejects > seqcount)) {
    opt_maxrejects = seqcount;
  }
  if ((opt_maxaccepts == 0) or (opt_maxaccepts > seqcount)) {
    opt_maxaccepts = seqcount;
  }

  tophits = opt_maxrejects + opt_maxaccepts + MAXDELAYED;
  tophits = std::min(tophits, seqcount);
}

auto search_done_paired() -> void {
  dbindex_free_paired();
  db_free();

  if (fp_matched != nullptr) {
    std::fclose(fp_matched);
  }
  if (fp_matched2 != nullptr) {
    std::fclose(fp_matched2);
  }
  if (fp_notmatched != nullptr) {
    std::fclose(fp_notmatched);
  }
  if (fp_notmatched2 != nullptr) {
    std::fclose(fp_notmatched2);
  }
  if (fp_uc != nullptr) {
    std::fclose(fp_uc);
  }
  if (fp_blast6out != nullptr) {
    std::fclose(fp_blast6out);
  }
  if (fp_userout != nullptr) {
    std::fclose(fp_userout);
    clean_up();
  }
  if (fp_alnout != nullptr) {
    std::fclose(fp_alnout);
  }

  show_rusage();
}

auto usearch_global_paired(struct Parameters const &parameters, char *cmdline,
                           char *progheader) -> void {
  search_prep_paired(cmdline, progheader);

  if (parameters.opt_dbmatched != nullptr) {
    fp_dbmatched = fopen_output(parameters.opt_dbmatched);
    if (fp_dbmatched == nullptr) {
      fatal("Unable to open dbmatched output file for writing");
    }
  }

  if (opt_dbmatched2 != nullptr) {
    fp_dbmatched2 = fopen_output(opt_dbmatched2);
    if (fp_dbmatched2 == nullptr) {
      fatal("Unable to open dbmatched2 output file for writing");
    }
  }

  if (parameters.opt_dbnotmatched != nullptr) {
    fp_dbnotmatched = fopen_output(parameters.opt_dbnotmatched);
    if (fp_dbnotmatched == nullptr) {
      fatal("Unable to open dbnotmatched output file for writing");
    }
  }

  if (opt_dbnotmatched2 != nullptr) {
    fp_dbnotmatched2 = fopen_output(opt_dbnotmatched2);
    if (fp_dbnotmatched2 == nullptr) {
      fatal("Unable to open dbnotmatched2 output file for writing");
    }
  }

  dbmatched = static_cast<uint64_t *>(xmalloc(static_cast<std::size_t>(seqcount) *
                                              sizeof(uint64_t)));
  std::memset(dbmatched, 0, static_cast<std::size_t>(seqcount) * sizeof(uint64_t));

  otutable_init();

  qmatches = 0;
  qmatches_abundance = 0;
  queries = 0;
  queries_abundance = 0;
  count_matched = 0;
  count_notmatched = 0;
  query_filename_left = parameters.opt_usearch_global;

  query_fastx_h_left = fastx_open(parameters.opt_usearch_global);
  if (query_fastx_h_left == nullptr) {
    fatal("Unrecognized query file type (not proper FASTA or FASTQ format): %s",
          parameters.opt_usearch_global);
  }

  if (parameters.opt_interleaved) {
    query_fastx_h_right = nullptr;
  } else {
    query_fastx_h_right = fastx_open(parameters.opt_reverse);
    if (query_fastx_h_right == nullptr) {
      fatal("Unrecognized split R2 query file type (not proper FASTA or FASTQ format): %s",
            parameters.opt_reverse);
    }
  }

  si_plus = static_cast<searchinfo_s_paired *>(
      xmalloc(static_cast<std::size_t>(opt_threads) * sizeof(searchinfo_s_paired)));
  if (opt_strand > 1) {
    si_minus = static_cast<searchinfo_s_paired *>(xmalloc(
        static_cast<std::size_t>(opt_threads) * sizeof(searchinfo_s_paired)));
  } else {
    si_minus = nullptr;
  }
  pthread = static_cast<pthread_t *>(
      xmalloc(static_cast<std::size_t>(opt_threads) * sizeof(pthread_t)));

  xpthread_mutex_init(&mutex_input, nullptr);
  xpthread_mutex_init(&mutex_output, nullptr);

  progress_init("Searching", fastx_get_size(query_fastx_h_left));
  search_thread_worker_run_paired();
  progress_done();

  xpthread_mutex_destroy(&mutex_output);
  xpthread_mutex_destroy(&mutex_input);

  xfree(pthread);
  xfree(si_plus);
  if (si_minus != nullptr) {
    xfree(si_minus);
  }

  if ((query_fastx_h_right != nullptr) and
      fastx_next(query_fastx_h_right, (opt_notrunclabels == 0),
                 chrmap_no_change_vector.data())) {
    fatal("More reverse queries than forward queries");
  }

  if (query_fastx_h_right != nullptr) {
    fastx_close(query_fastx_h_right);
  }
  fastx_close(query_fastx_h_left);

  if (!opt_quiet) {
    std::fprintf(stderr, "Matching unique query pairs: %d of %d", qmatches,
                 queries);
    if (queries > 0) {
      std::fprintf(stderr, " (%.2f%%)", 100.0 * qmatches / queries);
    }
    std::fprintf(stderr, "\n");
    if (opt_sizein) {
      std::fprintf(stderr,
                   "Matching total query pairs: %" PRIu64 " of %" PRIu64,
                   qmatches_abundance, queries_abundance);
      if (queries_abundance > 0) {
        std::fprintf(stderr, " (%.2f%%)",
                     100.0 * qmatches_abundance / queries_abundance);
      }
      std::fprintf(stderr, "\n");
    }
  }

  if (opt_log != nullptr) {
    std::fprintf(fp_log, "Matching unique query pairs: %d of %d", qmatches,
                 queries);
    if (queries > 0) {
      std::fprintf(fp_log, " (%.2f%%)", 100.0 * qmatches / queries);
    }
    std::fprintf(fp_log, "\n");
    if (opt_sizein) {
      std::fprintf(fp_log,
                   "Matching total query pairs: %" PRIu64 " of %" PRIu64,
                   qmatches_abundance, queries_abundance);
      if (queries_abundance > 0) {
        std::fprintf(fp_log, " (%.2f%%)",
                     100.0 * qmatches_abundance / queries_abundance);
      }
      std::fprintf(fp_log, "\n");
    }
  }

  if ((opt_otutabout != nullptr) or (opt_mothur_shared_out != nullptr) or
      (opt_biomout != nullptr)) {
    for (auto i = 0; i < seqcount; i++) {
      if (dbmatched[i] == 0U) {
        otutable_add(nullptr, db_records[static_cast<std::size_t>(i)].header.c_str(),
                     0);
      }
    }
    if (opt_unknown_name != nullptr) {
      otutable_add(nullptr, opt_unknown_name, 0);
    }
  }

  if (opt_biomout != nullptr) {
    otutable_print_biomout(fp_biomout);
    std::fclose(fp_biomout);
  }
  if (opt_otutabout != nullptr) {
    otutable_print_otutabout(fp_otutabout);
    std::fclose(fp_otutabout);
  }
  if (opt_mothur_shared_out != nullptr) {
    otutable_print_mothur_shared_out(fp_mothur_shared_out);
    std::fclose(fp_mothur_shared_out);
  }

  otutable_done();

  int count_dbmatched = 0;
  int count_dbnotmatched = 0;
  if ((fp_dbmatched != nullptr) or (fp_dbnotmatched != nullptr)) {
    for (auto i = 0; i < seqcount; i++) {
      auto const &record = db_records[static_cast<std::size_t>(i)];
      if (dbmatched[i] != 0U) {
        ++count_dbmatched;
        if (fp_dbmatched != nullptr) {
          write_db_pair_split_paired(fp_dbmatched, fp_dbmatched2, record,
                                     static_cast<int64_t>(dbmatched[i]),
                                     count_dbmatched);
        }
      } else {
        ++count_dbnotmatched;
        if (fp_dbnotmatched != nullptr) {
          write_db_pair_split_paired(fp_dbnotmatched, fp_dbnotmatched2, record,
                                     record.abundance, count_dbnotmatched);
        }
      }
    }
  }

  xfree(dbmatched);

  if (fp_dbmatched != nullptr) {
    std::fclose(fp_dbmatched);
  }
  if (fp_dbmatched2 != nullptr) {
    std::fclose(fp_dbmatched2);
  }
  if (fp_dbnotmatched != nullptr) {
    std::fclose(fp_dbnotmatched);
  }
  if (fp_dbnotmatched2 != nullptr) {
    std::fclose(fp_dbnotmatched2);
  }

  search_done_paired();
}
