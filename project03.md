Freie Universität Berlin, Robert-Koch Institute

Practical course: SARS-2 Bioinformatics & Data Science

Max von Kleist, Martin Hölzer

# Project 3

## Development of a risk score to analyze SARS-CoV-2 genomes

**Deadline**: 04.10.2024, 13:00

*The project should be worked out in groups of two or three students. Students should document their progress and present their work on Friday; 04.10.2024; 13:30am to the lecturers and students. The talk should be about 30-40 min, and allowing for 10-15min of questions.*

**a) The presentation should be uploaded via whiteboard as `Project3.pdf`, no later than the above stated deadline.**

An exemplary structure of the talk: 
*	Title including group composition (possibly a link to the code repository)
*	Short background (what is this good for?) & description of task
*	Methods: A step-by-step explanation of what has been done; troubleshooting
*	Results: A summary of the outcomes and analysis
*	Discussion of results; aspects that were unclear; how the workflow could be improved 
*	Outlook

**b) Codes for conducting the project tasks should be zipped into `Project3.zip` and uploaded either via whiteboard or made accessible via GitHub or GitLab.** In the latter case, a link to the repo should be provided on the slides and the repo must be accessible to the lecturers.

**Detailed Tasks:**

1) Download SARS-CoV-2 sequence set from [https://osf.io/vd7gj](https://osf.io/vd7gj) 
2) Align to the Index SARS-CoV-2 sequence from the Wuhan outbreak (NCBI: [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2))
  * Find mutant positions (e.g. CIGAR string or via some tool)
  * Find mutant positions in the Spike gene
3) Perform basic statistics
  * Number of ambiguous nucleotides per sequence
  * Number of mutations (SNPs, INDELs)
  * Project and cluster sequences (for example: PCA, kmeans)
4) Download SARS-CoV-2 antibody escape scores provided by the _Bloom Lab_ (paper: [https://www.science.org/doi/10.1126/science.abf9302](https://www.science.org/doi/10.1126/science.abf9302), online resource: [https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps](https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps))
5) Develop a "risk" score, e.g. based on the Spike sequences and the Bloom scores 
  * Among your example data, what's the most dangerous SARS-CoV-2 sample in your eyes?

