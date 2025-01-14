# Sample manifest structure

The only user-specified input to the workflow is the sample manifest. In the `example_manifest.csv` file shown above, there are four fields:

-   **finalreport_file**: a full path the GigaMUGA genotype file obtained from Neogen/TransnetYX.

-   **project_id:** A descriptive name given to the project to which each genotype file and covariate file (below) correspond. Nextflow concatenates all processed intensities and R/qtl2 processed data at the project level.

-   **covar_file:** A covariate file following convention described in the R/qtl2 documentation. This file will be referenced to create the cross object for each set of genotypes. The first three fields of this file must be `id, sex, gen`.

    The sample names provided in the `id` field must match the sample names in the provided **finalreport_file**. The pipeline is not built to parse every iteration of how one file's sample name can be manipulated to match another.

-   **cross_type:** The type of mouse cross of the experiment. For now, the only cross type than can be handled is 'do'. Other common crosses, in particular cross types enabling analyses in CC strains, CC-RIX mice, BXD strains, will be added soon.
