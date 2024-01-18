# IBD

## Intro
In this repository, I'm working on figuring out Inflammatory Bowel Diseases (IBD).
I know - very modest...
I do realize that the problem is very complex and requires several breakthroughs in immunology, microbiology, drug design, diagnostics, and probably many other fields.
But I still can't wrap my head around the fact that currently patients are treated with a trial-and-error approach, and that there is no way to predict whether a patient will respond to a given therapy. That's why, we'll start with this problem.
**The result should be a diagnostic test used by clinicians to assign patients to the right therapy.**

I've found an absolute ****load of data just lying around the Internet.
They're all from different studies, different labs, different countries, different years, different technologies, different everything, so it should be fun to try to make sense of it all.
But so far, I've collected nearly 1200 (sic!) samples of expression data from mucosal biopsies of IBD patients.
Some of them with response information so I think it's reeeeally worth a shot!

## Setup

After downloading the repo, first install the dependencies.
```
pip install -r requirements
```

Then, you can install the package (```pip install .```), or run it directly from the repo.
The package is runnable and currently serves one main purpose - to download and process the data.
You can run it like this:
```
python -m ibd run-pipeline -d GSE11223,GSE75214,GSE6731
```
This will download the data from the 3 datasets, process it, and save the results in the ```db``` directory.

You can also run the pipeline for all datasets in the database (currently in ```data/DAQ.csv``` file) by omitting the ```-d``` flag.
```
python -m ibd run-pipeline
`````

In order for it to work, there needs to be:
- a class to process the expression data for a given platform (currently supported platforms are: ```GPL1708```, ```GPL6244```, ```GPL570```, and ```GPL17996```, which cover 90% of the samples),
- a class to process the metadata, which is unfortunately dataset specific.

The package also contains a mock for runnig an experiment, but that's for later.
For now, I'll conduct all of the experiments in notebooks, and go back to this only after improving the whole data processing pipeline, which is a loooong way to go.

## Current state of affairs
I have implemented a comprehensive modeling pipeline in the ```experiment_playground.ipynb``` notebook.
The output model predicts Infliximab response in patients with CD and UC with 0.78 cross-validation AUC.
It identified 5 most importang genes: DCBLD1, IL13RA2, CSGALNACT2, WNK2, SNAPC1.
Some of those don't seem to make any sense, but some really, really do (e.g., IL13RA2 and WNK2).
Furthermore, the correlation analysis shows that most of these genes actually don't have any other genes correlated above 0.8, with the exception of CSGALNACT2.
This gene in and of itself doesn't seem to make sense at all, but it's correlated with a bunch of other genes that really do make sense: IL1R1, FGF7, LY96, MMP19.

Next steps:
- checking how the model performs on healthy controls,
- semi-supervised learning,
- going back and improving the whole pipeline.
