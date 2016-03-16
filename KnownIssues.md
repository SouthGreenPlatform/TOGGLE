# Known Issues

## Release 0.3.0

* BWA ALN driving to a 'EOF error' on some FASTQ:
 * you are probably using a BWA version between 0.7.2 and 0.7.6. Please use either a 0.6x version or a 0.7.8+ version.
* Error with options incorporating an equal symbol in their option value.
 * If you have an option such as **-l hostname=myNode**, it will block. A dev version is available at https://github.com/SouthGreenPlatform/TOGGLE-DEV for that.
* If adaptors are not informed in the software.config for cutadapt (-b options), it will drive to an error in Cutadapt step
* Cleaner and compressor will drive to an error if a step is too early compressed , ie before the pipeline needs it.
