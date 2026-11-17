CDK9 AND CDK7 PAPER — COMPETITION ASSAY FOLD CHANGE DATA
=========================================================
Prepared by: Jon Acosta  |  Last updated: 2026-09-13
For figure assembly by Sam Gould

OVERVIEW
--------
Each .txt file in this folder contains fold-change (FC) values from A549
competition assays. These are the numbers that go directly into GraphPad
Prism / figure panels. FC = (BFP% at T1 treated) / (BFP% at T1 untreated),
normalized to the T0 baseline within each replicate.

All values are means ± SEM across 3 biological replicates unless noted.

SOURCE DATA
-----------
Raw BFP% values and per-replicate FC calculations live in the Excel/XLS files
in the parent folder (Cdk9 and Cdk7 Paper Raw Data). These .txt files are
clean exports drawn from those source files — verified match.


FILE-BY-FILE GUIDE
==================

--- CDK9 ---

Cdk9 FC.txt
  Variants : NT (non-targeting), L156F, C106Y
  Drug     : Control vs KI-Cdk9d-32 (CDK9 PROTAC)
  Source   : Replicate 3/CDK Variants Comp Assays.xlsx  →  Round 3-JA tab

Cdk9 FC New Variant.txt
  Variants : ABE NT, CBE NT, L156F, I157V
  Drug     : Control vs KI-Cdk9d-32
  Source   : CDK9 new variant/CDK9 New Variant Comp Assay.xlsx  →  Round 3-JA tab

CBE Cdk9 FC Normalized New Variant + Drugs.txt
  Variants : CBE NT, L156F
  Drugs    : Control, KI-Cdk9d-32 (PROTAC), KB-0472 (KB), Atuveciclib (ATU), Voruciclib (VORU)
  Source   : A549 Competition Assay .../96 Well - U bottom/CBE Cdk9 + Drugs.xls

Cdk9 FC Normalized New Variant + Drugs.txt
  Variants : ABE NT, I157V
  Drugs    : Control, KI-Cdk9d-32, KB-0472, Atuveciclib, Voruciclib
  Source   : A549 Competition Assay .../96 Well - U bottom/ABE Cdk9 + Drugs.xls

--- CDK7 ---

Cdk7 FC.txt
  Variants : NT, L18F
  Drug     : Control vs SY-5609 (CDK7 inhibitor)
  Source   : Replicate 3/CDK Variants Comp Assays.xlsx  →  Round 3-JA tab

Cdk7 FC Normalized New Variant.txt
  Variants : ABE NT, V90A
  Drug     : Control vs SY-5609
  Source   : A549 Competition Assay .../96 Well - U bottom/ABE Cdk7 V90A SY-5609.xls

--- CDK12 ---

Cdk12 FC.txt
  Variants : NT, G731K, I733V
  Drug     : Control vs HQ461 (CDK12/13 inhibitor)
  Source   : Replicate 3/CDK Variants Comp Assays.xlsx  →  Round 3-JA tab

Cdk12 FC ZLY.txt
  Variants : ABE NT, CBE NT, G731K, I733V
  Drug     : Control vs ZLY-025
  Source   : ⚠ Source Excel not in this folder — data to be added by Jon

--- CDK13 ---

Cdk13 FC.txt
  Variants : NT, I711V
  Drug     : Control vs HQ461
  Source   : Replicate 3/CDK Variants Comp Assays.xlsx  →  Round 3-JA tab

Cdk13 FC ZLY.txt
  Variants : ABE NT, CBE NT, I711V
  Drug     : Control vs ZLY-025
  Source   : ⚠ Source Excel not in this folder — data to be added by Jon


COLUMN FORMAT (all files)
-------------------------
  Variant | Drug | Rep1 | Rep2 | Rep3 | Mean | SEM

DRUG ABBREVIATIONS
------------------
  PROTAC  = KI-Cdk9d-32   (selective CDK9 degrader)
  KB      = KB-0472
  ATU     = Atuveciclib
  VORU    = Voruciclib
  HQ461   = CDK12/13 inhibitor
  ZLY     = ZLY-025        (CDK12/13 inhibitor)
  SY      = SY-5609        (CDK7 inhibitor)
