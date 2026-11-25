# Assembly findings vocabulary

Vocabulary for the HTML report Main Finding and Detail columns.

## Plain words first

- An oligo is a short DNA fragment in the pool.
- COMMEC is the screening tool that looks for sequences of concern.
- OliGraph does not screen. It only tries to join short fragments into longer assembled pieces called contigs.
- COMMEC then screens those contigs.
- The report shows one row for each original oligo.
- Status is the final combined call after step 5.
- Main Finding says which evidence path produced a concern.
- Detailed findings name the workflow steps behind that label.

## Main findings summary

| Main Finding | Meaning |
| --- | --- |
| COMMEC direct | COMMEC found a concern in this short fragment itself, before any assembly. This is workflow step 2A. |
| COMMEC on OliGraph contig | This fragment alone did not show a concern. COMMEC found a concern after OliGraph joined it into a longer contig. Assembly is step 2B. Contig screening is step 4B. |
| COMMEC direct and OliGraph contig | COMMEC found a concern both in the original fragment and in the longer assembled piece. |
| No concern | COMMEC did not find a concern in the original fragment or in any assembled contig that includes it. |
| Incomplete | Screening could not finish for a practical reason, such as a fragment that is too short or too long, no contig, exclusion from assembly, or a failed screen. |

## Detailed findings

| What happened | Steps | Detailed finding |
| --- | --- | --- |
| Concern only on the oligo | 2A | 2A: COMMEC found a concern in this oligo |
| Concern only on a contig that includes this oligo | 2B + 4B | 4B: COMMEC found a concern in an OliGraph contig |
| Concern on both | 2A + 4B | 2A+4B: COMMEC found a concern in the oligo and its OliGraph contig |
| No concern on either | 2A + 4B | 2A+4B: COMMEC found no concern in the oligo or its OliGraph contig |
| Oligo too short; contig clear | 2A skip + 4B clear | 2A skipped (too short). 4B: no concern in OliGraph contig |
| Oligo and contig too short | 2A skip + 4B skip | 2A and 4B skipped (too short) |
| Oligo too long; contig clear | 2A skip + 4B clear | 2A skipped (too long). 4B: no concern in OliGraph contig |
| OliGraph made no contig for it | 2B miss | 2B: OliGraph built no contig for this oligo |
| Not eligible to assemble | 2B excluded | 2B: not eligible for OliGraph assembly |
| Contig screen did not finish | 4B failed | 4B: COMMEC did not finish screening the OliGraph contig |
| Direct screen failed | 2A failed | 2A: COMMEC failed while screening this oligo |

## Workflow steps used in Detail

| Step | Meaning |
| --- | --- |
| 2A | COMMEC screens each original oligo |
| 2B | OliGraph assembles eligible oligos into contigs |
| 4B | COMMEC screens those contigs and attributes concerns to source oligos |
| 5 | Combine 2A and 4B for each original oligo |

Step 6 on the flowchart is a pool classification rule. The current product report does not impose a pool verdict. It shows oligo counts only.
