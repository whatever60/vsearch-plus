# Paired Orientation Unification (R2 Kept Native)

This note documents the paired-orientation policy for the native paired pipeline.

## Policy

- External input/output orientation is now unified: R1 and R2 stay in their original read orientation.
- We no longer reverse-complement R2 as a persistent data transformation in the native paired VSEARCH path.
- Internal reverse handling is allowed only when required by a specific algorithm.

## Command-level behavior

For the native paired path (`fastx_uniques -> cluster_unoise -> uchime3_denovo -> usearch_global`) and paired RDP classification:

1. `fastx_uniques` (native paired mode)
- R2 is not reverse-complemented at ingestion.
- Right anchors are written in native R2 orientation.

2. `cluster_unoise` (native paired mode)
- Uses paired records from `fastx_uniques` in native orientation.
- No persistent R2 RC transform.

3. `uchime3_denovo` (native paired mode)
- Uses paired records in native orientation.
- Breakpoint scan logic is orientation-aware for the right segment:
  - left segment scanned from start to end
  - right segment scanned from end to start
  This preserves the intended combined-axis chimera model while keeping stored R2 unchanged.

4. `usearch_global` (native paired mode)
- Query and database right-end anchors are compared in native R2 orientation.
- No query-side forced R2 RC at ingestion.

5. RDP paired classifier (Java extension)
- Query-side orientation normalization remains internal to stock-RDP logic:
  - `trainingInfo.isSeqReversed(...)`
  - `ClassifierSequence.getReversedSeq()`
- This is an internal scoring detail; it does not change external FASTA/FASTQ orientation contracts.

## Developer guidance

- If you add new paired steps, keep persisted/output R2 in native orientation unless there is a hard compatibility reason not to.
- If a method needs orientation normalization for scoring, do it inside that method and keep it local to computation.
- Document any internal RC usage explicitly in parity notes.
