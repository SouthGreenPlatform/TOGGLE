# Known Issues

## Release 0.3.+

* *BWA ALN* driving to a 'EOF error' on some FASTQ:
 * you are probably using a *BWA* version between 0.7.2 and 0.7.6. Please use either a 0.6x version or a 0.7.8+ version.
* Error with options incorporating an equal symbol in their option value.
 * If you have an option such as **-l hostname=myNode**, it will block if you provide the option basically as **-l=hostname=MyNode** in the software config file.
 * You can overrride this behaviour using **hostanme==MyNode** to correct it (see *Manual.md*).
* If adaptors are not informed in the software.config for cutadapt (-b options), it will drive to an error in *Cutadapt* step
* Cleaner and compressor will drive to an error if a step is too early compressed , i.e. before the pipeline needs it.
* *picardToolsValidateSamFile* will always stop the pipeline on an error if the SAM is not perfectly correct. thus, if you want to use this tool, use it as the last one in the pipeline, that will drive whatever to an error.
