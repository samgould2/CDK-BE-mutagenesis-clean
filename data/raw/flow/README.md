# Flow cytometry

The exports behind the titering, MOI, selection-check and sort panels, one
directory per acquisition, named for its date and what was run.

`.fcs` files are as acquired. The `.wsp` files are the FlowJo workspaces that
define the gating.

**The workspaces' file paths were rewritten.** FlowJo stores an absolute URI
for every sample, which on the acquiring machine began with a home directory
and, in one case, a lab member's own Dropbox folder. Each URI now starts at
the acquisition directory instead -- `file:2024_10_3_A549_titer/RFP%20and%20GFP/
Specimen_001_1.fcs` -- so it names the sample without naming anyone, and under
`data/raw/flow/` it resolves against this directory. Nothing else in the
workspaces was touched, and all eight still parse. No code reads them: the
panels are drawn from the exported tables.
