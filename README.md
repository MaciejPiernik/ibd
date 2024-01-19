# IBD

## Intro
The goal of this project if to figure out Inflammatory Bowel Diseases (IBD).
I do realize that the problem is very complex and requires several breakthroughs in immunology, microbiology, drug design, diagnostics, and probably many other fields.
But I still can't wrap my head around the fact that currently patients are treated with a trial-and-error approach.
That's why, we're starting with the problem of drug response prediction.

## Current status
So far, I've found ~2500 samples of expression data from mucosal biopsies of IBD patients, spanning across 32 datasets and 13 different platforms (see ```data/DAQ.csv```).
All of the datasets are from the Gene Expression Omnibus (GEO).
To analyze this data together, I've implemented a data processing pipeline for 3 platforms and several datasets, allowing me to play around with ~1k samples.
The expressions are 
I have also implemented a comprehensive modeling pipeline in the ```ibd/notebooks/experiment_playground.ipynb``` notebook.
The output model predicts Infliximab response in patients with CD and UC with 0.78 cross-validation AUC.
It identified 5 most importang genes: DCBLD1, IL13RA2, CSGALNACT2, WNK2, SNAPC1.
Some of those don't seem to make any sense, but some really, really do (e.g., IL13RA2 and WNK2).
Furthermore, the correlation analysis shows that most of these genes actually don't have any other genes correlated above 0.8, with the exception of CSGALNACT2.
This gene in and of itself doesn't seem to make sense at all, but it's correlated with a bunch of other genes that do: IL1R1, FGF7, LY96, MMP19.

I checked how the model performs on healthy controls, and it indicates that the vast majority of them are predicted to be responders.
This made me think if drug response prediction is the same as predicting if a person should receive a treatment or not, and I have to say, I'm not so sure anymore.
Of course it is also very plausible that the model is simply wrong.
For example, it doesn't consider the severity of the disease, so maybe it simply reflects it (which probably is a good predictor of response to any treatment).
I also implemented a diagnostic model (0.88 CV AUC) to see if maybe some of the genes will be shared with the response model, but I didn't find any.

## Roadmap

### Next steps
- Semi-supervised model for Infliximab response prediction
- Adding next platforms to gain access to more data
- Going back and improving the whole pipeline
- Adding an actual database
- Collecting more data

### [Further On Up The Road](https://www.youtube.com/watch?v=h5aVK70P88k)
- Processing raw microarray data instead of processed expressions
- Productization: dockerization, orchestration, automation, cloud processing, ...
- Automating the process of building predictive models
- Foundation model for gene expression
- Application for non-coding users

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
