# ROOT to PDF file conversion

## Output

One single pdf file with first all TH1 histograms, and then all TH2 histograms.
### FILE SUBJECT TO CHANGE, READ THIS BEFORE USING IT

## Usage

If you want a .pdf file with the same name as .root file, in the same folder as that file, use

>./ROOTtoPDF.sh [ROOT_file_you_want_to_convert.root]

If you want a .pdf file with different  name somewhere else, use

>./ROOTtoPDF.sh [ROOT_file_you_want_to_convert.root] [path_to_new_pdf_file]

.pdf extension at the end is not necessary, but it might still be written