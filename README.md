# IBD

In this repository, I'm working on solving Inflammatory Bowel Diseases (IBD).
I know - very modest...
I do realize that the problem is very complex and requires several breakthroughs in immunology, microbiology, drug design, diagnostics, and probably many other fields.
But I still can't wrap my head around the fact that currently patients are treated with a trial-and-error approach, and that there is no way to predict whether a patient will respond to a given therapy.
That's why, we'll start with this problem.
**The result should be a diagnostic test used by clinicians to assign patients to the right therapy.**

I've found an absolute ****load of data just lying around the Internet.
They're all from different studies, different labs, different countries, different years, different technologies, different everything, so it should be fun to try to make sense of it all.
But so far, I've identified nearly 1200 (sic!) samples of expression data from mucosal biopsies of IBD patients.
Some of them with response information so I think it's reeeeally worth a shot!

## Current state of affairs
We have 3 platforms covering 90% of the samples and metadata processing for all these datasets.
The samples are normalized between datasets, so they can be analyzed together.
Tested the first, very simple model for infliximab prediction with ~0.75 AUC ROC.
It selected 5 important features ('DCBLD1', 'IL13RA2', 'CSGALNACT2', 'WNK2', 'SNAPC1'), out of which IL13RA2 and WNK2 have been identified in a previous study which uses parts of this data, so this confirms to some extent that our platform works!

Next steps:
- further modeling,
- semi-supervised learning,
- going back and improving the whole pipeline.
