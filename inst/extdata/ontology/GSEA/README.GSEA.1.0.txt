Gene Set Enrichment Analysis (GSEA-P-R) Readme page -- Broad Institute -- July 2005

These are the instructions to run the R version of the GSEA program (GSEA-P-R.ZIP). Notice that there is a more user friendly version of GSEA-P written in Java and that has a web interface (www.broad.mit.edu/gsea). If you want to run GSEA and you are not a programmer or a computational biologist that version may be a better choice. The R version is intended for more computational experienced biologists, bioinformaticians or computational biologists. 

The GSEA-P-R program described here reflects the version of the methodology described and used in the Subramanian and Tamayo et al 2005 paper. For details about the method and the content of the output please see Supporting Information for that paper.

You need to install R release 2.0 or later.

- Copy the GSEA-P-R.ZIP file to your computer. 
- Unzip the file GSEA-P-R.ZIP using the option to create subdirectories.
  This should create the following files and subdirectories:

GSEA program and functions in R (all the GSEA code is conatined there):

GSEA/GSEA-P-R/GSEA.1.0.R        

Directory with input datasets, gct and cls files:
  GSEA/GSEA-P-R/Datasets/        
                         Gender.gct
                         Gender.cls
                         Leukemia.gct
                         Leukemia.cls
                         Lung_Boston.gct
                         Lung_Boston.cls
                         Lung_Michigan.gct
                         Lung_Michigan.cls
                         Lung_Stanford.gct
                         Lung_Stanford.cls
                         Lung_Bost_maxed_common_Mich_Bost.gct
                         Lung_Mich_maxed_common_Mich_Bost.gct
                         P53.gct
                         P53.cls

Directory with gene set databases, gmt files:
  GSEA/GSEA-P-R/GeneSetDatabases/
                                 C1.gmt
                                 C2.gmt
                                 C3.gmt
                                 C4.gmt
                                 Lung_Boston_poor_outcome.gmt
                                 Lung_Michigan_poor_outcome.gmt

Directories with results of running the examples described in the paper:

  GSEA/GSEA-P-R/Gender_C1/
                          Gender_C2
                          Leukemia_C1
                          Lung_Boston_C2
                          Lung_Stanford_C2 
                          Lung_Michigan_C2
                          Lung_Boston_outcome 
                          Lung_Michigan_outcome
                          P53_C2

The top 20 high scoring gene sets are reported in table S2 (Supporting Information).

One page R scripts to run the examples described in the paper:

  GSEA/GSEA-P-R/
                Run.Gender_C1.R
                Run.Gender_C2.R
                Run.Leukemia_C1.R
                Run.Lung_Boston_C2.R
                Run.Lung_Stanford_C2.R
                Run.Lung_Michigan_C2.R
                Run.Lung_Boston_outcome.R
                Run.Lung_Michigan_outcome.R
                Run.P53_C2.R

To run, for example, the Leukemia dataset with the C1 gene set database go to the file GSEA/GSEA-P-R/Run.Leukemia_C1.R and change the file pathnames to reflect the location of the GSEA directory in your machine. For example if you expanded the ZIP file under your directory "C:/my_directory" you need to change the line: 

GSEA.program.location <- "d:/CGP2005/GSEA/GSEA-P-R/GSEA.1.0.R"  

To:

GSEA.program.location <- "c:my_directory/GSEA/GSEA-P-R/GSEA.1.0.R"

 And the same change to each pathname in that file: you need to replace "d:/CGP2005" with "C"/my_directory".

 You may also want to change the line:

doc.string            = "Leukemia_C1",

To:

doc.string            = "my_run_of_Leukemia_C1",

or any other prefix label you want to give your results. This way you won't overwrite the original results that come in those directories and can use them for comparison with the results of you own run. 

After the pathnames have been changed to reflect the location of the directories in your machine to run GSEA program just open the R GUI and paste the content of the Run.<example>.R files on it.  Fro example to run the Leukemia vs. C1 example use the contents of the file "Run.Leukemia_C1.R" The program is self-contained and should run and produce the results under the directory "C:my_directory/GSEA/GSEA-P-R/Leukemia_C1". These files are set up with the parameters used in the examples of the paper (e.g. to produce detailed results for the significant and top 20 gene sets). You may want to start using these parameters and change them only when needed and when you get mnore experience with the program. For details of what are the effect of changing some of the parameters see the Supporting Information document.

If you want to run a completely new dataset the easiest way is:

i) Create a new directory: e.g. GSEA/GSEA-P-R/my_dataset, where you can store the inputs and outputs of running GSEA on those files. 
ii) Convert manually your files to *.gct (expression dataset) and *.cls (phenotype labels)
iii) Use Run.Leukemia_C1.R as a template to make a new script to run your data.
iv) Change the relevant pathnames to point to your input files in directory my_dataset. Change the doc.string to an approprote prefix name for your files.
v) Cut and paste the contents of this new script file in the R GUI to run it. The results will be stored in my_directory.

The GSEA-P-R program reads input files in *.gct, *.cls and *.gmt formats. As you can see from the examples's files these are simple tab separated ASCII files. If your datasets are not in this format you can use a text editor to convert them. If you start with a tab separated ASCII file tipically the conversion would consist in  modifying the header lines on top of the file.

If you have questions or problems running or using the program please  send them to gsea@broad.mit.edu. Also lets us know if you find GSEA a useful tool in your work.




