# xtandem-tsv-converter
Java tool to convert X!Tandem result files to a simple tab-delimited text format.

If the source MGF file is available as well, the spectrum's title will be added to the output file. Additionally, if the MGF file contains an annotated sequence (using the "SEQ=" or the "SEQUENCE=" field) this sequence will be available as well.

## Usage
```
java -jar xtandem-tsv-converter.jar [X!Tandem result file] [FDR] [(optional) source MGF file]
```
