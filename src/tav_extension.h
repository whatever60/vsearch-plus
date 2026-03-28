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

#ifndef TAV_EXTENSION_H
#define TAV_EXTENSION_H

/* Legacy paired-end extension surface.
   Kept only as historical reference; no active CLI path or build target
   depends on these declarations anymore. */

struct Parameters;

auto tav_is_catalog_file(char const *filename) -> bool;

auto tav_fastx_uniques(struct Parameters const &parameters) -> void;

auto tav_cluster_unoise(struct Parameters const &parameters) -> void;

auto tav_uchime3_denovo(struct Parameters const &parameters) -> void;

auto tav_usearch_global(struct Parameters const &parameters, char *cmdline,
                        char *progheader) -> void;

#endif
