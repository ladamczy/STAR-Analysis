# ROOT to PDF file conversion

## Output

One single pdf file with TH1 and TH2 histograms in order of saving in .root file.
#### SCRIPT SUBJECT TO CHANGE, READ THIS BEFORE USING IT

## Usage

### Preparation

1. Make a symbolic link with
   >ln -s path_to_ROOTtoPDF.sh_in_shared_folder

   in place you wish to use it from
   
2. use it as listed below (paths are supposed to be relative to symbolic link position)
3. (optional) you can declare your config in Styles.h file

### Commands

If you want a .pdf file with the same name as .root file, in the same folder as that file, use

>./ROOTtoPDF.sh [ROOT_file_you_want_to_convert.root]

If you want a .pdf file with different name somewhere else, use

>./ROOTtoPDF.sh [ROOT_file_you_want_to_convert.root] [path_to_new_pdf_file]

.pdf extension at the end is not necessary, but it might still be written

If you want special config, you can do it this way:

>./ROOTtoPDF.sh -m [main config number] -o [1D histogram number] -t [2D histogram number] [ROOT_file_you_want_to_convert.root]

or

>./ROOTtoPDF.sh -m [main config number] -o [1D histogram number] -t [2D histogram number] [ROOT_file_you_want_to_convert.root] [path_to_new_pdf_file]

If you do not specify a config number, the default one will be used.