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

#include "cluster_paired.h"

#include "tav_extension.h"

/*
  Temporary compatibility shim.
  We route paired cluster_unoise calls through this entrypoint first,
  then replace internals incrementally with stock-callstack-parity code.
*/
auto cluster_unoise_paired(struct Parameters const &parameters) -> void {
  tav_cluster_unoise(parameters);
}
